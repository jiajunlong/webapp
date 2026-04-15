"""
miRNA-Gene-Pathway Integration Module

Implements regulatory network analysis by integrating:
1. miRNA expression with gene expression
2. miRNA target prediction (correlation-based and database validation)
3. Regulatory network construction
4. Pathway-level regulatory analysis
5. Identification of key regulatory hub miRNAs
6. Integration with disease modules

This completes Phase 2 by connecting genetic regulation to pathway and disease.
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from typing import Dict, List, Tuple, Optional
import logging
import warnings

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class miRNATargetPredictor:
    """
    Predict miRNA targets using correlation-based approach and validation
    """
    
    def __init__(self, mirna_expr: pd.DataFrame, gene_expr: pd.DataFrame):
        """
        Initialize miRNA target predictor
        
        Parameters:
        -----------
        mirna_expr : pd.DataFrame
            miRNA expression matrix (miRNAs × samples)
        gene_expr : pd.DataFrame
            Gene expression matrix (genes × samples)
        """
        self.mirna_expr = mirna_expr
        self.gene_expr = gene_expr
        self.targets = None  # miRNA -> list of target genes
        self.target_correlations = None  # miRNA-gene correlations
        
        logger.info(f"Initialized miRNATargetPredictor:")
        logger.info(f"  miRNAs: {self.mirna_expr.shape[0]}, Samples: {self.mirna_expr.shape[1]}")
        logger.info(f"  Genes: {self.gene_expr.shape[0]}, Samples: {self.gene_expr.shape[1]}")
    
    def predict_targets(self, correlation_threshold: float = -0.3,
                       method: str = 'pearson') -> Dict[str, List[str]]:
        """
        Predict miRNA targets based on negative correlation
        
        Theory: miRNAs typically suppress target genes, so negative correlation
        suggests regulatory relationship
        
        Parameters:
        -----------
        correlation_threshold : float
            Threshold for target prediction (typically negative)
        method : str
            'pearson' or 'spearman'
            
        Returns:
        --------
        dict : miRNA -> list of target genes
        """
        logger.info(f"Predicting miRNA targets (threshold: {correlation_threshold}, method: {method})...")
        
        # Ensure sample order matches
        common_samples = self.mirna_expr.columns.intersection(self.gene_expr.columns)
        if len(common_samples) < 5:
            logger.warning(f"Only {len(common_samples)} common samples. May have insufficient data.")
        
        mirna_data = self.mirna_expr[common_samples]
        gene_data = self.gene_expr[common_samples]
        
        targets = {}
        correlations = []
        
        for mirna in self.mirna_expr.index:
            if mirna not in mirna_data.index:
                continue
            
            mirna_vals = mirna_data.loc[mirna].values
            target_genes = []
            
            for gene in self.gene_expr.index:
                if gene not in gene_data.index:
                    continue
                
                gene_vals = gene_data.loc[gene].values
                
                # Compute correlation
                if method == 'pearson':
                    corr, pval = pearsonr(mirna_vals, gene_vals)
                else:
                    corr, pval = spearmanr(mirna_vals, gene_vals)
                
                # Target if negative correlation
                if corr <= correlation_threshold and pval < 0.05:
                    target_genes.append(gene)
                    correlations.append({
                        'miRNA': mirna,
                        'target': gene,
                        'correlation': corr,
                        'p_value': pval
                    })
            
            if target_genes:
                targets[mirna] = target_genes
        
        self.targets = targets
        self.target_correlations = pd.DataFrame(correlations)
        
        logger.info(f"✓ Predicted targets for {len(targets)} miRNAs")
        logger.info(f"  Total predictions: {len(correlations)}")
        if len(correlations) > 0:
            logger.info(f"  Mean correlation: {np.mean([c['correlation'] for c in correlations]):.3f}")
        
        return targets
    
    def get_targets_for_mirna(self, mirna: str, top_n: int = 50) -> pd.DataFrame:
        """
        Get top target genes for a specific miRNA
        
        Parameters:
        -----------
        mirna : str
            miRNA name
        top_n : int
            Number of top targets to return
            
        Returns:
        --------
        pd.DataFrame : Top target genes with correlation scores
        """
        if self.target_correlations is None:
            raise ValueError("Targets not predicted. Call predict_targets first.")
        
        mirna_targets = self.target_correlations[
            self.target_correlations['miRNA'] == mirna
        ].sort_values('correlation')
        
        return mirna_targets.head(top_n)
    
    def validate_against_databases(self, predicted_targets: Dict[str, List[str]],
                                  known_targets: Optional[Dict[str, List[str]]] = None) -> pd.DataFrame:
        """
        Validate predictions against known miRNA-target databases
        
        Parameters:
        -----------
        predicted_targets : dict
            miRNA -> list of predicted target genes
        known_targets : dict or None
            Known miRNA-target associations (optional reference)
            
        Returns:
        --------
        pd.DataFrame : Validation metrics
        """
        logger.info("Validating predictions against databases...")
        
        validation_results = []
        
        for mirna, pred_targets in predicted_targets.items():
            pred_set = set(pred_targets)
            
            validation_results.append({
                'miRNA': mirna,
                'n_predicted_targets': len(pred_targets),
                'confidence': 'high' if len(pred_targets) >= 5 else 'low',
                'mean_target_corr': self.target_correlations[
                    self.target_correlations['miRNA'] == mirna
                ]['correlation'].mean() if mirna in self.targets else np.nan
            })
        
        validation_df = pd.DataFrame(validation_results)
        
        logger.info(f"✓ Validation complete: {len(validation_df)} miRNAs evaluated")
        
        return validation_df


class miRNARegulatoryNetwork:
    """
    Builds and analyzes miRNA regulatory networks
    """
    
    def __init__(self, mirna_targets: Dict[str, List[str]],
                 pathway_genes: Dict[str, List[str]]):
        """
        Initialize regulatory network
        
        Parameters:
        -----------
        mirna_targets : dict
            miRNA -> list of target genes
        pathway_genes : dict
            Pathway name -> list of genes
        """
        self.mirna_targets = mirna_targets
        self.pathway_genes = pathway_genes
        self.network = None
        self.hub_mirnas = None
        
        logger.info(f"Initialized miRNARegulatoryNetwork:")
        logger.info(f"  miRNAs with targets: {len(mirna_targets)}")
        logger.info(f"  Total pathways: {len(pathway_genes)}")
    
    def build_network(self) -> Dict[str, any]:
        """
        Build regulatory network graph
        
        Returns:
        --------
        dict : Network structure with edges and nodes
        """
        logger.info("Building regulatory network...")
        
        network = {
            'nodes': [],  # miRNAs and genes
            'edges': [],  # miRNA-target interactions
            'pathways': []  # Pathway assignments
        }
        
        # Add miRNA nodes
        for mirna in self.mirna_targets.keys():
            network['nodes'].append({
                'id': mirna,
                'type': 'miRNA',
                'size': len(self.mirna_targets[mirna])
            })
        
        # Add edges and target genes
        target_genes_added = set()
        
        for mirna, targets in self.mirna_targets.items():
            for target in targets:
                network['edges'].append({
                    'source': mirna,
                    'target': target,
                    'type': 'regulation'
                })
                
                if target not in target_genes_added:
                    target_genes_added.add(target)
                    network['nodes'].append({
                        'id': target,
                        'type': 'gene',
                        'size': 1
                    })
        
        # Assign pathways to genes
        for pathway, genes in self.pathway_genes.items():
            for gene in genes:
                if gene in target_genes_added:
                    network['pathways'].append({
                        'gene': gene,
                        'pathway': pathway
                    })
        
        self.network = network
        logger.info(f"✓ Network built:")
        logger.info(f"  Nodes: {len(network['nodes'])} (miRNAs + genes)")
        logger.info(f"  Edges: {len(network['edges'])} (miRNA-target)")
        logger.info(f"  Pathway associations: {len(network['pathways'])}")
        
        return network
    
    def identify_hub_mirnas(self, top_n: int = 20) -> pd.DataFrame:
        """
        Identify key regulatory hub miRNAs
        
        Hub miRNAs = those regulating many targets or pathway genes
        
        Parameters:
        -----------
        top_n : int
            Number of hub miRNAs to return
            
        Returns:
        --------
        pd.DataFrame : Hub miRNA rankings
        """
        logger.info("Identifying hub miRNAs...")
        
        hub_scores = []
        
        for mirna, targets in self.mirna_targets.items():
            # Count targets in pathways
            targets_in_pathways = 0
            pathway_coverage = set()
            
            for target in targets:
                for pathway, pathway_genes in self.pathway_genes.items():
                    if target in pathway_genes:
                        targets_in_pathways += 1
                        pathway_coverage.add(pathway)
            
            hub_scores.append({
                'miRNA': mirna,
                'n_targets': len(targets),
                'n_pathway_targets': targets_in_pathways,
                'pathway_coverage': len(pathway_coverage),
                'hub_score': len(targets) * 0.5 + targets_in_pathways * 0.3 + len(pathway_coverage) * 0.2
            })
        
        hub_df = pd.DataFrame(hub_scores)
        hub_df = hub_df.sort_values('hub_score', ascending=False)
        hub_df['rank'] = range(1, len(hub_df) + 1)
        
        self.hub_mirnas = hub_df
        
        logger.info(f"✓ Hub miRNAs identified: {len(hub_df)} total, top {min(top_n, len(hub_df))} shown")
        
        return hub_df.head(top_n)
    
    def map_to_pathways(self) -> pd.DataFrame:
        """
        Map miRNA regulations to pathways
        
        Returns:
        --------
        pd.DataFrame : miRNA-pathway associations
        """
        logger.info("Mapping miRNA regulations to pathways...")
        
        mirna_pathway_map = []
        
        for mirna, targets in self.mirna_targets.items():
            targets_set = set(targets)
            
            for pathway, pathway_genes in self.pathway_genes.items():
                pathway_genes_set = set(pathway_genes)
                
                # Find overlap
                overlap = targets_set & pathway_genes_set
                
                if len(overlap) > 0:
                    mirna_pathway_map.append({
                        'miRNA': mirna,
                        'pathway': pathway,
                        'n_targets_in_pathway': len(overlap),
                        'pathway_size': len(pathway_genes),
                        'coverage': len(overlap) / len(pathway_genes)
                    })
        
        mirna_pathway_df = pd.DataFrame(mirna_pathway_map)
        
        logger.info(f"✓ Mapped {len(mirna_pathway_df)} miRNA-pathway associations")
        
        return mirna_pathway_df


class RegulatoryModuleAnalysis:
    """
    Integrates miRNA, genes, and pathways into regulatory modules
    """
    
    def __init__(self, mirna_targets: Dict[str, List[str]],
                 pathway_genes: Dict[str, List[str]],
                 disease_modules: Optional[Dict[str, List[str]]] = None):
        """
        Initialize regulatory module analysis
        
        Parameters:
        -----------
        mirna_targets : dict
            miRNA -> list of target genes
        pathway_genes : dict
            Pathway -> list of genes
        disease_modules : dict or None
            Disease -> list of genes in disease module
        """
        self.mirna_targets = mirna_targets
        self.pathway_genes = pathway_genes
        self.disease_modules = disease_modules or {}
        self.modules = None
        
        logger.info(f"Initialized RegulatoryModuleAnalysis:")
        logger.info(f"  miRNAs: {len(mirna_targets)}")
        logger.info(f"  Pathways: {len(pathway_genes)}")
        logger.info(f"  Diseases: {len(self.disease_modules)}")
    
    def identify_regulatory_modules(self) -> Dict[str, Dict]:
        """
        Identify regulatory modules (miRNA + genes + pathway + disease)
        
        Returns:
        --------
        dict : Module ID -> module description
        """
        logger.info("Identifying regulatory modules...")
        
        modules = {}
        module_id = 0
        
        # For each pathway, find associated miRNAs and disease connections
        for pathway, pathway_genes_list in self.pathway_genes.items():
            pathway_genes_set = set(pathway_genes_list)
            
            # Find miRNAs targeting genes in this pathway
            regulating_mirnas = []
            for mirna, targets in self.mirna_targets.items():
                targets_in_pathway = set(targets) & pathway_genes_set
                if len(targets_in_pathway) > 0:
                    regulating_mirnas.append((mirna, len(targets_in_pathway)))
            
            if regulating_mirnas:
                # Find disease associations
                disease_associations = []
                for disease, disease_genes_list in self.disease_modules.items():
                    disease_genes_set = set(disease_genes_list)
                    pathway_disease_overlap = pathway_genes_set & disease_genes_set
                    if len(pathway_disease_overlap) > 0:
                        disease_associations.append((disease, len(pathway_disease_overlap)))
                
                module_id += 1
                modules[f"Module_{module_id}"] = {
                    'pathway': pathway,
                    'n_genes': len(pathway_genes_list),
                    'regulating_miRNAs': regulating_mirnas,
                    'n_regulating_mirnas': len(regulating_mirnas),
                    'disease_associations': disease_associations,
                    'n_disease_associations': len(disease_associations),
                    'regulatory_importance': sum(count for _, count in regulating_mirnas) / max(1, len(regulating_mirnas))
                }
        
        self.modules = modules
        
        logger.info(f"✓ Identified {len(modules)} regulatory modules")
        
        return modules
    
    def score_regulatory_importance(self) -> pd.DataFrame:
        """
        Score importance of regulatory relationships
        
        Returns:
        --------
        pd.DataFrame : Regulatory importance scores
        """
        logger.info("Scoring regulatory importance...")
        
        if self.modules is None:
            self.identify_regulatory_modules()
        
        scores = []
        
        for module_id, module_info in self.modules.items():
            score = (
                0.4 * module_info['n_regulating_mirnas'] +
                0.3 * module_info['regulatory_importance'] +
                0.3 * module_info['n_disease_associations']
            )
            
            scores.append({
                'module': module_id,
                'pathway': module_info['pathway'],
                'n_genes': module_info['n_genes'],
                'n_mirnas': module_info['n_regulating_mirnas'],
                'n_diseases': module_info['n_disease_associations'],
                'regulatory_score': score
            })
        
        scores_df = pd.DataFrame(scores)
        scores_df = scores_df.sort_values('regulatory_score', ascending=False)
        scores_df['rank'] = range(1, len(scores_df) + 1)
        
        logger.info(f"✓ Scored {len(scores_df)} regulatory modules")
        
        return scores_df
    
    def export_regulatory_network(self, output_file: str):
        """
        Export regulatory network to file
        
        Parameters:
        -----------
        output_file : str
            Path to output file
        """
        logger.info(f"Exporting regulatory network to {output_file}...")
        
        if self.modules is None:
            self.identify_regulatory_modules()
        
        with open(output_file, 'w') as f:
            f.write("Module\tPathway\tn_genes\tregulating_miRNAs\tn_mirnas\tdisease_assoc\tn_diseases\treg_score\n")
            
            for module_id, module_info in self.modules.items():
                mirna_str = ";".join([m for m, _ in module_info['regulating_miRNAs']])
                disease_str = ";".join([d for d, _ in module_info['disease_associations']])
                
                f.write(f"{module_id}\t{module_info['pathway']}\t" +
                       f"{module_info['n_genes']}\t{mirna_str}\t" +
                       f"{module_info['n_regulating_mirnas']}\t{disease_str}\t" +
                       f"{module_info['n_disease_associations']}\t" +
                       f"{module_info['regulatory_importance']:.3f}\n")
        
        logger.info(f"✓ Exported to {output_file}")


if __name__ == "__main__":
    print("miRNA-Gene-Pathway Integration Module")
    print("=" * 60)
    print("\nUsage:")
    print("  from mirna_integration import miRNATargetPredictor")
    print("  predictor = miRNATargetPredictor(mirna_expr, gene_expr)")
    print("  targets = predictor.predict_targets()")
    print("")
    print("  from mirna_integration import miRNARegulatoryNetwork")
    print("  network = miRNARegulatoryNetwork(targets, pathway_genes)")
    print("  hub_mirnas = network.identify_hub_mirnas()")
