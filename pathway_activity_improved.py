"""
Improved Pathway Activity Scoring with log transformation
"""

import pandas as pd
import numpy as np
from scipy.stats import zscore
from typing import Dict, List
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PathwayActivityScorerImproved:
    """Score pathway activity with log transformation for better variance"""
    
    def __init__(self, pathway_genes: Dict[str, List[str]]):
        self.pathway_genes = pathway_genes
        self.pathway_activity = None
        self.gene_expr = None
        self.mapped_pathways = {}
        
    def load_expression_data(self, expr_file: str):
        """Load and log-transform expression data"""
        logger.info(f"Loading expression data from {expr_file}")
        self.gene_expr = pd.read_csv(expr_file, index_col=0)
        
        # Log transform (add pseudocount to avoid log(0))
        self.gene_expr = np.log2(self.gene_expr + 1)
        
        logger.info(f"✓ Loaded and log2-transformed: {self.gene_expr.shape}")
        return self.gene_expr.shape
    
    def _map_pathway_genes(self) -> Dict[str, List[str]]:
        """Map pathway genes to expression matrix"""
        if self.gene_expr is None:
            raise ValueError("Expression data not loaded")
        
        expr_genes = set(self.gene_expr.index)
        mapped = {}
        
        for pathway, genes in self.pathway_genes.items():
            genes_in_expr = [g for g in genes if g in expr_genes]
            if len(genes_in_expr) >= 2:
                mapped[pathway] = genes_in_expr
        
        logger.info(f"Mapped {len(mapped)}/{len(self.pathway_genes)} pathways")
        self.mapped_pathways = mapped
        return mapped
    
    def score_gsva(self) -> pd.DataFrame:
        """Score pathways using GSVA approach on log-transformed data"""
        logger.info("Computing pathway activity scores (GSVA on log-transformed data)")
        
        mapped = self._map_pathway_genes()
        
        # Initialize with float64 to avoid dtype issues
        activity = pd.DataFrame(
            np.nan,
            index=mapped.keys(),
            columns=self.gene_expr.columns,
            dtype=np.float64
        )
        
        # Score each pathway
        for pathway_idx, (pathway, genes) in enumerate(mapped.items()):
            if (pathway_idx + 1) % 50 == 0:
                logger.info(f"  Scoring pathway {pathway_idx + 1}/{len(mapped)}")
            
            # Get expression values
            pathway_expr = self.gene_expr.loc[genes]
            
            # Z-score normalize
            z_scores = pathway_expr.apply(zscore, axis=0, nan_policy='omit')
            
            # Mean z-score
            activity.loc[pathway] = z_scores.mean(axis=0)
        
        self.pathway_activity = activity
        logger.info(f"✓ Computed pathway activity: {activity.shape}")
        logger.info(f"  Activity range: [{activity.min().min():.3f}, {activity.max().max():.3f}]")
        logger.info(f"  Mean activity: {activity.values.mean():.3f} ± {activity.values.std():.3f}")
        
        return self.pathway_activity


if __name__ == "__main__":
    # Load pathways
    logger.info("Loading pathways...")
    pathway_df = pd.read_csv("data/pathway(基因名映射版).tsv", sep='\t')
    pathway_genes = {}
    
    for _, row in pathway_df.iterrows():
        pathway_name = row['Pathway_Name']
        genes_str = row['Gene']
        
        if pd.isna(genes_str):
            continue
        
        genes = [g.strip() for g in str(genes_str).replace(';', ',').split(',')]
        genes = [g for g in genes if g and g != 'NA']
        
        if genes:
            pathway_genes[pathway_name] = genes
    
    # Score pathways
    scorer = PathwayActivityScorerImproved(pathway_genes)
    scorer.load_expression_data('TCGA-COAD/filtered_hiseq_data.csv')
    activity = scorer.score_gsva()
    
    # Save
    activity.to_csv('pathway_activity_log_matrix.csv')
    logger.info("✓ Saved pathway activity to pathway_activity_log_matrix.csv")
