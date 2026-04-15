"""
Pathway Activity Scoring Module
Implements multiple pathway scoring methods from gene expression data.
- GSVA-like scoring (mean z-score of pathway genes per sample)
- Mean expression scoring (baseline)
- Supports TCGA and other expression matrices
"""

import pandas as pd
import numpy as np
from scipy.stats import zscore
from typing import Dict, List, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PathwayActivityScorer:
    """Score pathway activity from gene expression data"""
    
    def __init__(self, pathway_genes: Dict[str, List[str]]):
        """
        Initialize pathway activity scorer
        
        Parameters:
        -----------
        pathway_genes : dict
            Mapping of pathway name -> list of gene symbols
        """
        self.pathway_genes = pathway_genes
        self.pathway_activity = None
        self.gene_expr = None
        self.mapped_pathways = {}  # Track mapping success
        
    def load_expression_data(self, expr_file: str) -> Tuple[int, int]:
        """
        Load gene expression matrix from file
        
        Parameters:
        -----------
        expr_file : str
            Path to CSV/TSV expression file (genes × samples)
            First column should be gene names/IDs
            
        Returns:
        --------
        tuple : (n_genes, n_samples) dimensions
        """
        logger.info(f"Loading expression data from {expr_file}")
        self.gene_expr = pd.read_csv(expr_file, index_col=0)
        
        # Validate data
        if self.gene_expr.shape[0] == 0 or self.gene_expr.shape[1] == 0:
            raise ValueError("Expression matrix is empty")
        
        logger.info(f"✓ Loaded expression: {self.gene_expr.shape} (genes × samples)")
        logger.info(f"  Gene names: {self.gene_expr.index[:5].tolist()}...")
        logger.info(f"  Samples: {self.gene_expr.columns[:5].tolist()}...")
        
        return self.gene_expr.shape
    
    def _map_pathway_genes(self) -> Dict[str, List[str]]:
        """
        Map pathway genes to expression matrix indices
        
        Returns:
        --------
        dict : pathway -> list of genes in expression matrix
        """
        if self.gene_expr is None:
            raise ValueError("Expression data not loaded. Call load_expression_data first.")
        
        expr_genes = set(self.gene_expr.index)
        mapped = {}
        
        for pathway, genes in self.pathway_genes.items():
            genes_in_expr = [g for g in genes if g in expr_genes]
            if len(genes_in_expr) >= 2:  # Require at least 2 genes
                mapped[pathway] = genes_in_expr
        
        n_mapped = len(mapped)
        n_total = len(self.pathway_genes)
        pct = 100 * n_mapped / n_total if n_total > 0 else 0
        
        logger.info(f"Pathway mapping: {n_mapped}/{n_total} pathways ({pct:.1f}%) " +
                   f"have ≥2 genes in expression matrix")
        
        self.mapped_pathways = mapped
        return mapped
    
    def score_gsva(self) -> pd.DataFrame:
        """
        Score pathways using GSVA-like approach
        
        Method: For each pathway and sample, compute z-scores of genes,
        then take the mean z-score as pathway activity.
        
        This is a simplified version of GSVA that is computationally fast.
        
        Returns:
        --------
        pd.DataFrame : pathway activity matrix (pathways × samples)
            - Rows: pathway names
            - Columns: sample names
            - Values: pathway activity scores (z-score normalized)
        """
        logger.info("Computing pathway activity scores (GSVA method)")
        
        mapped = self._map_pathway_genes()
        
        # Initialize output matrix
        activity = pd.DataFrame(
            np.nan,
            index=mapped.keys(),
            columns=self.gene_expr.columns,
            dtype=np.float32
        )
        
        # Score each pathway
        for pathway_idx, (pathway, genes) in enumerate(mapped.items()):
            if (pathway_idx + 1) % 50 == 0:
                logger.info(f"  Scoring pathway {pathway_idx + 1}/{len(mapped)}")
            
            # Get expression values for this pathway's genes
            pathway_expr = self.gene_expr.loc[genes]
            
            # Z-score normalize genes across samples
            z_scores = pathway_expr.apply(zscore, axis=0, nan_policy='omit')
            
            # Mean z-score per sample = pathway activity
            activity.loc[pathway] = z_scores.mean(axis=0)
        
        self.pathway_activity = activity
        logger.info(f"✓ Computed pathway activity: {activity.shape}")
        logger.info(f"  Activity range: [{activity.min().min():.2f}, {activity.max().max():.2f}]")
        
        return self.pathway_activity
    
    def score_mean(self) -> pd.DataFrame:
        """
        Simple mean expression method
        
        Returns average expression level of pathway genes per sample.
        Useful as baseline / validation method.
        
        Returns:
        --------
        pd.DataFrame : pathway activity matrix (pathways × samples)
        """
        logger.info("Computing pathway activity scores (Mean method)")
        
        mapped = self._map_pathway_genes()
        
        # Initialize output matrix
        activity = pd.DataFrame(
            np.nan,
            index=mapped.keys(),
            columns=self.gene_expr.columns,
            dtype=np.float32
        )
        
        # Score each pathway
        for pathway, genes in mapped.items():
            pathway_expr = self.gene_expr.loc[genes]
            activity.loc[pathway] = pathway_expr.mean(axis=0)
        
        self.pathway_activity = activity
        logger.info(f"✓ Computed pathway activity: {activity.shape}")
        
        return self.pathway_activity
    
    def get_pathway_statistics(self) -> pd.DataFrame:
        """
        Get summary statistics for each pathway
        
        Returns:
        --------
        pd.DataFrame : pathway statistics
            Columns: pathway, mean_activity, std_activity, min, max, n_genes
        """
        if self.pathway_activity is None:
            raise ValueError("Pathway activity not computed. Call score_gsva or score_mean first.")
        
        mapped = self.mapped_pathways
        stats = []
        
        for pathway in self.pathway_activity.index:
            activity_vals = self.pathway_activity.loc[pathway]
            stats.append({
                'pathway': pathway,
                'mean_activity': activity_vals.mean(),
                'std_activity': activity_vals.std(),
                'min_activity': activity_vals.min(),
                'max_activity': activity_vals.max(),
                'n_genes': len(mapped.get(pathway, []))
            })
        
        return pd.DataFrame(stats)
    
    def save_pathway_activity(self, output_file: str):
        """Save pathway activity matrix to file"""
        if self.pathway_activity is None:
            raise ValueError("Pathway activity not computed")
        
        self.pathway_activity.to_csv(output_file)
        logger.info(f"✓ Saved pathway activity to {output_file}")


if __name__ == "__main__":
    # Test with TCGA data
    import sys
    
    logger.info("Testing PathwayActivityScorer with TCGA-COAD data")
    
    # Load pathway data
    logger.info("Loading pathway genes...")
    pathway_df = pd.read_csv("data/pathway(基因名映射版).tsv", sep='\t')
    pathway_genes = {}
    
    for _, row in pathway_df.iterrows():
        pathway_name = row['Pathway_Name']
        genes_str = row['Gene']
        
        if pd.isna(genes_str):
            continue
        
        # Parse genes (may be comma or semicolon separated)
        genes = [g.strip() for g in str(genes_str).replace(';', ',').split(',')]
        genes = [g for g in genes if g and g != 'NA']
        
        if genes:
            pathway_genes[pathway_name] = genes
    
    logger.info(f"Loaded {len(pathway_genes)} pathways")
    
    # Initialize scorer
    scorer = PathwayActivityScorer(pathway_genes)
    
    # Load expression data
    scorer.load_expression_data('TCGA-COAD/filtered_hiseq_data.csv')
    
    # Score pathways
    logger.info("\nScoring pathways with GSVA method...")
    activity = scorer.score_gsva()
    
    # Get statistics
    stats = scorer.get_pathway_statistics()
    logger.info("\nPathway Statistics (first 10):")
    logger.info(stats.head(10).to_string())
    
    # Save results
    scorer.save_pathway_activity('pathway_activity_matrix.csv')
    stats.to_csv('pathway_statistics.csv', index=False)
    
    logger.info("\n✅ Pathway activity scoring complete!")
