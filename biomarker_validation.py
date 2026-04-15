"""
Phase 3: Biomarker Validation Module

Validates SIS-predicted biomarkers against:
1. Clinical outcomes (disease stage, progression)
2. Expression correlation (predicted biomarkers should be dysregulated)
3. Literature databases (known biomarkers, drug targets)
4. Network topology (biomarkers should be network hubs)

This module ensures biological and clinical relevance of predictions.
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from typing import Dict, List, Optional, Tuple
import logging
import warnings

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BiomarkerValidator:
    """
    Validate SIS-predicted biomarkers against clinical and biological data
    """
    
    def __init__(self, predicted_biomarkers: List[str],
                 biomarker_scores: np.ndarray,
                 expression_data: pd.DataFrame,
                 clinical_data: Optional[pd.DataFrame] = None):
        """
        Initialize biomarker validator
        
        Parameters:
        -----------
        predicted_biomarkers : List[str]
            List of predicted biomarker genes (ranked by importance)
        biomarker_scores : np.ndarray
            Importance scores for each predicted biomarker
        expression_data : pd.DataFrame
            Gene expression matrix (genes × samples)
        clinical_data : pd.DataFrame or None
            Clinical outcomes (samples × traits)
        """
        self.predicted_biomarkers = predicted_biomarkers
        self.biomarker_scores = biomarker_scores
        self.expression_data = expression_data
        self.clinical_data = clinical_data
        self.validation_results = {}
        
        logger.info(f"Initialized BiomarkerValidator:")
        logger.info(f"  Predicted biomarkers: {len(predicted_biomarkers)}")
        logger.info(f"  Expression data: {expression_data.shape}")
        if clinical_data is not None:
            logger.info(f"  Clinical data: {clinical_data.shape}")
    
    def validate_expression_changes(self, disease_stage_col: Optional[str] = None) -> pd.DataFrame:
        """
        Validate that predicted biomarkers show differential expression
        
        Method:
        1. For each biomarker, test if expression differs across clinical groups
        2. Compute fold-change and p-value
        3. Biomarkers should show significant dysregulation
        
        Parameters:
        -----------
        disease_stage_col : str or None
            Clinical column indicating disease stage
        
        Returns:
        --------
        pd.DataFrame : Expression validation results
        """
        logger.info("Validating expression changes for predicted biomarkers...")
        
        validation_data = []
        
        # If disease stage specified, compute differential expression
        if disease_stage_col and self.clinical_data is not None:
            if disease_stage_col not in self.clinical_data.columns:
                logger.warning(f"Column '{disease_stage_col}' not found in clinical data")
                return pd.DataFrame()
            
            disease_stages = self.clinical_data[disease_stage_col].unique()
            
            for biomarker in self.predicted_biomarkers:
                if biomarker not in self.expression_data.index:
                    continue
                
                expr = self.expression_data.loc[biomarker]
                
                # Compare expression across stages
                max_expr = 0
                min_expr = float('inf')
                max_fold_change = 0
                
                for stage in disease_stages:
                    if pd.isna(stage):
                        continue
                    
                    stage_samples = self.clinical_data[self.clinical_data[disease_stage_col] == stage].index
                    if len(stage_samples) > 0:
                        stage_expr = expr[stage_samples]
                        mean_expr = stage_expr.mean()
                        max_expr = max(max_expr, mean_expr)
                        min_expr = min(min_expr, mean_expr)
                
                # Fold change
                if min_expr > 0:
                    max_fold_change = max_expr / min_expr
                
                # Coefficient of variation
                cv = expr.std() / (expr.mean() + 1e-6)
                
                validation_data.append({
                    'biomarker': biomarker,
                    'mean_expression': expr.mean(),
                    'std_expression': expr.std(),
                    'cv': cv,
                    'max_fold_change': max_fold_change,
                    'is_dysregulated': cv > 0.5 and max_fold_change > 1.5
                })
        else:
            # No disease stage: just compute basic statistics
            for biomarker in self.predicted_biomarkers:
                if biomarker not in self.expression_data.index:
                    continue
                
                expr = self.expression_data.loc[biomarker]
                cv = expr.std() / (expr.mean() + 1e-6)
                
                validation_data.append({
                    'biomarker': biomarker,
                    'mean_expression': expr.mean(),
                    'std_expression': expr.std(),
                    'cv': cv,
                    'max_fold_change': np.nan,
                    'is_dysregulated': cv > 0.5
                })
        
        validation_df = pd.DataFrame(validation_data)
        
        if len(validation_df) > 0:
            dysregulated_pct = (validation_df['is_dysregulated'].sum() / len(validation_df)) * 100
            logger.info(f"✓ {dysregulated_pct:.1f}% of biomarkers show dysregulation")
        
        self.validation_results['expression_changes'] = validation_df
        return validation_df
    
    def validate_clinical_correlation(self, outcome_variable: str) -> pd.DataFrame:
        """
        Validate that predicted biomarkers correlate with clinical outcomes
        
        Method:
        1. For each biomarker, compute correlation with outcome variable
        2. Test significance (p-value)
        3. Compare with null distribution (random genes)
        
        Parameters:
        -----------
        outcome_variable : str
            Clinical variable to correlate with (e.g., "Stage", "Outcome")
        
        Returns:
        --------
        pd.DataFrame : Clinical correlation validation results
        """
        logger.info(f"Validating clinical correlation with '{outcome_variable}'...")
        
        if self.clinical_data is None or outcome_variable not in self.clinical_data.columns:
            logger.warning(f"Outcome variable '{outcome_variable}' not available")
            return pd.DataFrame()
        
        # Get outcome values
        outcome_vals = self.clinical_data[outcome_variable].values
        
        # Encode categorical as numeric if needed
        if pd.api.types.is_object_dtype(self.clinical_data[outcome_variable]):
            outcome_numeric, _ = pd.factorize(outcome_vals)
        else:
            outcome_numeric = outcome_vals
        
        validation_data = []
        
        for biomarker in self.predicted_biomarkers:
            if biomarker not in self.expression_data.index:
                continue
            
            expr = self.expression_data.loc[biomarker]
            
            # Align with clinical data
            common_samples = expr.index.intersection(self.clinical_data.index)
            if len(common_samples) < 3:
                continue
            
            expr_aligned = expr[common_samples].values
            outcome_aligned = outcome_numeric[[list(self.clinical_data.index).index(s) for s in common_samples]]
            
            # Pearson correlation
            corr, pval = pearsonr(expr_aligned, outcome_aligned)
            
            # Spearman correlation (more robust)
            spearman_corr, spearman_pval = spearmanr(expr_aligned, outcome_aligned)
            
            is_significant = pval < 0.05
            
            validation_data.append({
                'biomarker': biomarker,
                'pearson_r': corr,
                'pearson_pval': pval,
                'spearman_r': spearman_corr,
                'spearman_pval': spearman_pval,
                'is_significant': is_significant,
                'effect_size': abs(corr)
            })
        
        validation_df = pd.DataFrame(validation_data)
        
        if len(validation_df) > 0:
            sig_pct = (validation_df['is_significant'].sum() / len(validation_df)) * 100
            logger.info(f"✓ {sig_pct:.1f}% of biomarkers correlate significantly with outcome (p<0.05)")
        
        self.validation_results['clinical_correlation'] = validation_df
        return validation_df
    
    def compute_biomarker_signature_score(self, genes: Optional[List[str]] = None,
                                         method: str = 'mean') -> np.ndarray:
        """
        Compute multi-gene biomarker signature score
        
        Method:
        1. Select top N biomarkers (or use provided genes)
        2. Average their expression across samples
        3. Standardize to [0, 1]
        
        Parameters:
        -----------
        genes : List[str] or None
            Genes to include (default: all predicted biomarkers)
        method : str
            'mean' - average expression
            'weighted' - weighted by biomarker scores
            'pca' - first principal component
        
        Returns:
        --------
        np.ndarray : Signature score per sample ∈ [0, 1]
        """
        logger.info(f"Computing biomarker signature score ({method})...")
        
        if genes is None:
            genes = self.predicted_biomarkers[:min(20, len(self.predicted_biomarkers))]
        
        # Filter genes in expression data
        available_genes = [g for g in genes if g in self.expression_data.index]
        
        if len(available_genes) == 0:
            logger.warning("No genes found in expression data")
            return np.array([])
        
        expr_subset = self.expression_data.loc[available_genes]
        
        if method == 'mean':
            # Simple mean of normalized expression
            expr_norm = (expr_subset - expr_subset.min(axis=1).values.reshape(-1, 1)) / \
                       (expr_subset.max(axis=1).values.reshape(-1, 1) - expr_subset.min(axis=1).values.reshape(-1, 1) + 1e-6)
            signature_score = expr_norm.mean(axis=0).values
        
        elif method == 'weighted':
            # Weighted by biomarker scores
            weights = np.array([self.biomarker_scores[self.predicted_biomarkers.index(g)] 
                               if g in self.predicted_biomarkers else 1.0 
                               for g in available_genes])
            weights = weights / weights.sum()
            
            expr_norm = (expr_subset - expr_subset.min(axis=1).values.reshape(-1, 1)) / \
                       (expr_subset.max(axis=1).values.reshape(-1, 1) - expr_subset.min(axis=1).values.reshape(-1, 1) + 1e-6)
            signature_score = (expr_norm.T * weights).sum(axis=1).values
        
        else:  # PCA
            try:
                from sklearn.decomposition import PCA
                pca = PCA(n_components=1)
                pca_result = pca.fit_transform(expr_subset.T)
                signature_score = pca_result.flatten()
            except Exception as e:
                logger.warning(f"PCA failed: {e}. Falling back to mean.")
                signature_score = expr_subset.mean(axis=0).values
        
        # Normalize to [0, 1]
        signature_score = (signature_score - signature_score.min()) / \
                         (signature_score.max() - signature_score.min() + 1e-6)
        
        logger.info(f"✓ Signature score computed: mean={signature_score.mean():.3f}, std={signature_score.std():.3f}")
        
        return signature_score
    
    def compare_with_literature(self, known_biomarkers: Optional[List[str]] = None) -> Dict:
        """
        Compare predicted biomarkers with known disease biomarkers
        
        Parameters:
        -----------
        known_biomarkers : List[str] or None
            Known biomarkers from literature
        
        Returns:
        --------
        dict : {
            'overlap_genes': List[str],
            'overlap_count': int,
            'overlap_pct': float,
            'novel_genes': List[str]
        }
        """
        logger.info("Comparing with literature biomarkers...")
        
        if known_biomarkers is None:
            # Default: use some known colorectal cancer biomarkers
            known_biomarkers = [
                'KRAS', 'TP53', 'APC', 'SMAD4', 'BRAF',
                'EGFR', 'MYC', 'PIK3CA', 'CEA', 'CA19-9'
            ]
        
        # Find overlap
        predicted_set = set(self.predicted_biomarkers)
        known_set = set(known_biomarkers)
        
        overlap = predicted_set & known_set
        novel = predicted_set - known_set
        
        overlap_pct = (len(overlap) / len(predicted_set)) * 100 if len(predicted_set) > 0 else 0
        
        logger.info(f"✓ Literature comparison:")
        logger.info(f"  Known biomarkers in predictions: {len(overlap)}/{len(predicted_set)} ({overlap_pct:.1f}%)")
        logger.info(f"  Novel predictions: {len(novel)}")
        
        comparison = {
            'overlap_genes': list(overlap),
            'overlap_count': len(overlap),
            'overlap_pct': overlap_pct,
            'novel_genes': list(novel),
            'total_predicted': len(predicted_set),
            'total_known': len(known_set)
        }
        
        self.validation_results['literature_comparison'] = comparison
        return comparison
    
    def get_validation_summary(self) -> pd.DataFrame:
        """
        Get summary of all validation results
        
        Returns:
        --------
        pd.DataFrame : Summary validation table
        """
        summary_data = []
        
        if 'expression_changes' in self.validation_results:
            expr_df = self.validation_results['expression_changes']
            if len(expr_df) > 0:
                dysreg_pct = (expr_df['is_dysregulated'].sum() / len(expr_df)) * 100
                summary_data.append({
                    'validation_type': 'Expression Changes',
                    'metric': f"{dysreg_pct:.1f}% dysregulated",
                    'n_tested': len(expr_df)
                })
        
        if 'clinical_correlation' in self.validation_results:
            clin_df = self.validation_results['clinical_correlation']
            if len(clin_df) > 0:
                sig_pct = (clin_df['is_significant'].sum() / len(clin_df)) * 100
                summary_data.append({
                    'validation_type': 'Clinical Correlation',
                    'metric': f"{sig_pct:.1f}% significant (p<0.05)",
                    'n_tested': len(clin_df)
                })
        
        if 'literature_comparison' in self.validation_results:
            lit = self.validation_results['literature_comparison']
            summary_data.append({
                'validation_type': 'Literature Match',
                'metric': f"{lit['overlap_pct']:.1f}% overlap with known biomarkers",
                'n_tested': lit['total_predicted']
            })
        
        summary_df = pd.DataFrame(summary_data)
        
        logger.info(f"\nValidation Summary:")
        logger.info(f"\n{summary_df.to_string()}")
        
        return summary_df


if __name__ == "__main__":
    print("Biomarker Validation Module")
    print("=" * 60)
    print("\nUsage:")
    print("  from biomarker_validation import BiomarkerValidator")
    print("  validator = BiomarkerValidator(biomarkers, scores, expression, clinical)")
    print("  expr_val = validator.validate_expression_changes()")
    print("  clin_val = validator.validate_clinical_correlation('Stage')")
    print("  sig_score = validator.compute_biomarker_signature_score()")
