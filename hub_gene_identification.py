"""
Hub Gene Identification Module

Identifies the most important genes within each pathway by combining:
1. Network centrality (degree + betweenness)
2. Expression variance in TCGA samples
3. Pathway-specific metrics (number of pathway connections)

This bridges the Gene → Pathway level in the multi-scale framework.
"""

import pandas as pd
import numpy as np
import networkx as nx
from typing import Dict, List, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class HubGeneIdentifier:
    """
    Identify hub genes within pathways.
    
    Hub genes are those that:
    - Have high degree centrality (many connections within pathway)
    - Have high betweenness centrality (bridge between pathway regions)
    - Show high expression variance (variable across patients)
    - Have high connectivity within their pathway
    """
    
    def __init__(self, pathway_genes: Dict[str, List[str]], 
                 gene_network: Optional[nx.Graph] = None):
        """
        Initialize hub gene identifier.
        
        Parameters:
        -----------
        pathway_genes : dict
            Mapping of pathway name -> list of gene symbols
        gene_network : nx.Graph or None
            Network graph with genes as nodes. If None, will use pathway membership
            as the graph structure (simple pathway graph).
        """
        self.pathway_genes = pathway_genes
        self.gene_network = gene_network
        self.hub_scores = {}  # pathway -> DataFrame of hub scores
        self.expr_data = None
        
        logger.info(f"Initialized HubGeneIdentifier with {len(pathway_genes)} pathways")
        if gene_network is not None:
            logger.info(f"  Gene network: {gene_network.number_of_nodes()} nodes, " +
                       f"{gene_network.number_of_edges()} edges")
    
    def load_expression_data(self, expr_file: str):
        """
        Load gene expression data to calculate expression variance.
        
        Parameters:
        -----------
        expr_file : str
            Path to CSV/TSV expression file (genes × samples)
        """
        logger.info(f"Loading expression data from {expr_file}")
        self.expr_data = pd.read_csv(expr_file, index_col=0)
        logger.info(f"✓ Loaded expression: {self.expr_data.shape} (genes × samples)")
    
    def _build_pathway_subgraph(self, pathway: str, genes: List[str]) -> nx.Graph:
        """
        Build subgraph containing only genes in the pathway.
        
        If a full gene_network was provided, extract the subgraph.
        Otherwise, create a complete graph (all genes connected in pathway).
        
        Parameters:
        -----------
        pathway : str
            Pathway name
        genes : List[str]
            List of gene symbols in pathway
            
        Returns:
        --------
        nx.Graph : Subgraph of pathway genes
        """
        genes_set = set(genes)
        
        if self.gene_network is not None:
            # Extract subgraph from full network
            # Include only nodes in the pathway
            available_genes = [g for g in genes if g in self.gene_network]
            if len(available_genes) > 0:
                subgraph = self.gene_network.subgraph(available_genes).copy()
            else:
                # No genes in network; create complete graph
                subgraph = nx.complete_graph(len(genes))
                # Map node indices to gene names
                mapping = {i: genes[i] for i in range(len(genes))}
                subgraph = nx.relabel_nodes(subgraph, mapping)
        else:
            # No network provided; use pathway membership as complete graph
            # (all pathway members are potentially connected)
            subgraph = nx.complete_graph(len(genes))
            mapping = {i: genes[i] for i in range(len(genes))}
            subgraph = nx.relabel_nodes(subgraph, mapping)
        
        return subgraph
    
    def calculate_hub_score(self, pathway: str, 
                           normalize: bool = True) -> pd.DataFrame:
        """
        Calculate hub score for genes in a pathway.
        
        Hub score combines:
        1. Degree centrality (within pathway subgraph)
        2. Betweenness centrality (bridges between genes)
        3. Expression variance (biological variability)
        4. Pathway connectivity (within-pathway connections)
        
        Parameters:
        -----------
        pathway : str
            Pathway name
        normalize : bool
            If True, normalize scores to [0, 1] range
            
        Returns:
        --------
        pd.DataFrame
            Columns: gene, degree, betweenness, expr_var, pathway_connectivity,
                     hub_score, rank
        """
        if pathway not in self.pathway_genes:
            raise ValueError(f"Pathway '{pathway}' not in pathway_genes")
        
        genes = self.pathway_genes[pathway]
        
        # Build pathway subgraph
        subgraph = self._build_pathway_subgraph(pathway, genes)
        
        scores = []
        
        for gene in genes:
            metrics = {'gene': gene}
            
            # 1. Degree centrality (number of connections)
            if gene in subgraph:
                degree = subgraph.degree(gene)
            else:
                degree = 0
            metrics['degree'] = degree
            
            # 2. Betweenness centrality (bridges between other genes)
            if gene in subgraph and subgraph.number_of_nodes() > 2:
                try:
                    # Calculate betweenness for the whole subgraph
                    betweenness_dict = nx.betweenness_centrality(subgraph)
                    betweenness = betweenness_dict.get(gene, 0.0)
                except:
                    betweenness = 0.0
            else:
                betweenness = 0.0
            metrics['betweenness'] = betweenness
            
            # 3. Expression variance (from TCGA data)
            expr_var = 0.0
            if self.expr_data is not None and gene in self.expr_data.index:
                expr_values = self.expr_data.loc[gene]
                expr_var = float(expr_values.var())
            metrics['expr_variance'] = expr_var
            
            # 4. Pathway connectivity (count of pathway neighbors)
            pathway_connectivity = subgraph.degree(gene) if gene in subgraph else 0
            metrics['pathway_connectivity'] = pathway_connectivity
            
            # 5. Hub score (weighted combination)
            # Normalize components first if requested
            if normalize and subgraph.number_of_nodes() > 1:
                max_degree = max([subgraph.degree(n) for n in subgraph.nodes()]) or 1
                norm_degree = degree / max_degree
            else:
                norm_degree = degree
            
            norm_betweenness = betweenness  # Already normalized by nx.betweenness_centrality
            
            # Normalize expression variance
            if self.expr_data is not None:
                all_variances = self.expr_data.var(axis=1)
                max_var = all_variances.max()
                norm_expr_var = expr_var / max_var if max_var > 0 else 0
            else:
                norm_expr_var = 0
            
            # Composite hub score: weighted average
            # Weights: degree=40%, betweenness=30%, expression variance=30%
            hub_score = (0.4 * norm_degree + 
                        0.3 * norm_betweenness + 
                        0.3 * norm_expr_var)
            
            metrics['hub_score'] = hub_score
            scores.append(metrics)
        
        # Create DataFrame and rank
        df = pd.DataFrame(scores)
        df = df.sort_values('hub_score', ascending=False)
        df['rank'] = range(1, len(df) + 1)
        
        # Store in cache
        self.hub_scores[pathway] = df
        
        return df.reset_index(drop=True)
    
    def get_top_hub_genes(self, pathway: str, top_n: int = 10) -> pd.DataFrame:
        """
        Get top N hub genes in a pathway.
        
        Parameters:
        -----------
        pathway : str
            Pathway name
        top_n : int
            Number of top genes to return
            
        Returns:
        --------
        pd.DataFrame
            Top N genes with hub scores
        """
        if pathway not in self.hub_scores:
            self.calculate_hub_score(pathway)
        
        df = self.hub_scores[pathway]
        return df.head(top_n)
    
    def calculate_all_hub_genes(self, pathways: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
        """
        Calculate hub genes for multiple pathways.
        
        Parameters:
        -----------
        pathways : List[str] or None
            List of pathway names to calculate. If None, calculate for all pathways.
            
        Returns:
        --------
        dict : pathway -> DataFrame of hub scores
        """
        if pathways is None:
            pathways = list(self.pathway_genes.keys())
        
        logger.info(f"Calculating hub genes for {len(pathways)} pathways...")
        
        for i, pathway in enumerate(pathways):
            if (i + 1) % max(1, len(pathways) // 10) == 0:
                logger.info(f"  Processed {i + 1}/{len(pathways)} pathways...")
            
            try:
                self.calculate_hub_score(pathway)
            except Exception as e:
                logger.warning(f"  Error processing '{pathway}': {e}")
        
        logger.info(f"✓ Calculated hub genes for {len(self.hub_scores)} pathways")
        return self.hub_scores
    
    def export_hub_genes_summary(self, output_file: str, top_n: int = 10):
        """
        Export hub genes for all pathways to a CSV file.
        
        Parameters:
        -----------
        output_file : str
            Path to output CSV file
        top_n : int
            Number of top genes per pathway to export
        """
        logger.info(f"Exporting hub genes (top {top_n} per pathway) to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("pathway,rank,gene,hub_score,degree,betweenness,expr_variance,pathway_connectivity\n")
            
            for pathway in sorted(self.hub_scores.keys()):
                df = self.hub_scores[pathway].head(top_n)
                for _, row in df.iterrows():
                    f.write(f"{pathway},{int(row['rank'])},{row['gene']}," +
                           f"{row['hub_score']:.6f},{int(row['degree'])}," +
                           f"{row['betweenness']:.6f},{row['expr_variance']:.6f}," +
                           f"{int(row['pathway_connectivity'])}\n")
        
        logger.info(f"✓ Exported to {output_file}")
    
    def get_hub_gene_summary(self, top_n_pathways: int = 20, 
                            top_n_genes: int = 5) -> pd.DataFrame:
        """
        Get a summary table of top hub genes from top pathways.
        
        Parameters:
        -----------
        top_n_pathways : int
            Number of top pathways to include
        top_n_genes : int
            Number of top genes per pathway
            
        Returns:
        --------
        pd.DataFrame
            Summary table with columns: pathway, gene_rank, gene, hub_score
        """
        summary = []
        
        # Get pathways with highest mean hub score
        pathway_scores = {}
        for pathway, df in self.hub_scores.items():
            pathway_scores[pathway] = df['hub_score'].mean()
        
        top_pathways = sorted(pathway_scores.items(), 
                             key=lambda x: x[1], 
                             reverse=True)[:top_n_pathways]
        
        for pathway, _ in top_pathways:
            if pathway in self.hub_scores:
                df = self.hub_scores[pathway].head(top_n_genes)
                for _, row in df.iterrows():
                    summary.append({
                        'pathway': pathway,
                        'gene_rank': int(row['rank']),
                        'gene': row['gene'],
                        'hub_score': row['hub_score'],
                        'degree': int(row['degree'])
                    })
        
        return pd.DataFrame(summary)


# ============================================================================
# Utility functions for batch processing
# ============================================================================

def identify_hub_genes_in_pathway(pathway: str, 
                                  genes: List[str],
                                  expr_file: str,
                                  gene_network: Optional[nx.Graph] = None) -> pd.DataFrame:
    """
    Convenient function to identify hub genes in a single pathway.
    
    Parameters:
    -----------
    pathway : str
        Pathway name
    genes : List[str]
        Genes in the pathway
    expr_file : str
        Path to expression data
    gene_network : nx.Graph or None
        Gene network (optional)
        
    Returns:
    --------
    pd.DataFrame
        Hub gene scores for the pathway
    """
    pathway_dict = {pathway: genes}
    identifier = HubGeneIdentifier(pathway_dict, gene_network)
    identifier.load_expression_data(expr_file)
    return identifier.calculate_hub_score(pathway)


if __name__ == "__main__":
    # Example usage
    print("Hub Gene Identification Module")
    print("=" * 60)
    print("\nUsage:")
    print("  from hub_gene_identification import HubGeneIdentifier")
    print("  identifier = HubGeneIdentifier(pathway_genes_dict, gene_network)")
    print("  identifier.load_expression_data('data/TCGA-COAD/filtered_hiseq_data.csv')")
    print("  hub_scores = identifier.calculate_all_hub_genes()")
    print("  top_genes = identifier.get_top_hub_genes('Glycolysis / Gluconeogenesis', top_n=10)")
