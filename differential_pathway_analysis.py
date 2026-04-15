"""
Differential Pathway Analysis Module
Compare pathway activity across patient clinical groups (age, sex, disease stage)
using statistical tests with multiple testing correction.
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, f_oneway, mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
from typing import Dict, List, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DifferentialPathwayAnalysis:
    """Compare pathway activity across clinical groups"""
    
    def __init__(self, pathway_activity: pd.DataFrame, 
                 clinical_data: pd.DataFrame):
        """
        Initialize differential pathway analysis
        
        Parameters:
        -----------
        pathway_activity : pd.DataFrame
            Pathway activity matrix (pathways × samples)
            Rows: pathway names
            Columns: sample IDs (must match clinical_data index)
        clinical_data : pd.DataFrame
            Clinical metadata (samples × variables)
            Index: sample IDs
            Columns: clinical variables (age_group, gender, disease_stage, etc.)
        """
        # Ensure index alignment
        common_samples = pathway_activity.columns.intersection(clinical_data.index)
        if len(common_samples) < len(pathway_activity.columns) * 0.9:
            logger.warning(f"Only {len(common_samples)}/{len(pathway_activity.columns)} " +
                          "samples matched between pathway activity and clinical data")
        
        self.pathway_activity = pathway_activity[common_samples]
        self.clinical_data = clinical_data.loc[common_samples]
        self.results = {}
        
        logger.info(f"Initialized with {len(self.pathway_activity)} pathways, " +
                   f"{len(self.pathway_activity.columns)} samples")
    
    def compare_by_group(self, clinical_var: str, method='auto') -> pd.DataFrame:
        """
        Compare pathway activity between patient groups
        
        Parameters:
        -----------
        clinical_var : str
            Name of clinical variable to group by (e.g., 'age_group', 'gender', 'disease_stage')
        method : str
            Statistical test: 'auto' (infers from group count), 'ttest', 'anova', 'kruskal'
            - 'auto': t-test for 2 groups, ANOVA for 3+ groups (parametric)
            - 'ttest': t-test (2 groups only)
            - 'anova': One-way ANOVA (2+ groups, parametric)
            - 'kruskal': Kruskal-Wallis (non-parametric, use if data non-normal)
        
        Returns:
        --------
        pd.DataFrame : Results table with columns
            - pathway: pathway name
            - statistic: test statistic
            - pvalue: unadjusted p-value
            - padj: FDR-corrected p-value
            - significant: boolean (padj < 0.05)
            - group: clinical variable tested
            - n_groups: number of groups compared
            - group_labels: comma-separated group names
        """
        
        logger.info(f"Comparing pathway activity by '{clinical_var}'")
        
        # Get unique groups
        groups = sorted(self.clinical_data[clinical_var].unique())
        n_groups = len(groups)
        
        logger.info(f"  Found {n_groups} groups: {groups}")
        
        # Infer method if 'auto'
        if method == 'auto':
            if n_groups == 2:
                method = 'ttest'
            else:
                method = 'anova'
        
        logger.info(f"  Using method: {method}")
        
        results = []
        
        # Test each pathway
        for pathway_idx, pathway in enumerate(self.pathway_activity.index):
            if (pathway_idx + 1) % 50 == 0:
                logger.info(f"  Testing pathway {pathway_idx + 1}/{len(self.pathway_activity)}")
            
            pathway_vals = self.pathway_activity.loc[pathway]
            
            # Get values per group (aligned with clinical data index)
            group_vals = []
            group_sizes = []
            
            for group in groups:
                mask = self.clinical_data[clinical_var] == group
                vals = pathway_vals[mask].dropna().values
                group_vals.append(vals)
                group_sizes.append(len(vals))
            
            # Perform statistical test
            stat, pval = np.nan, np.nan
            
            if method == 'ttest':
                if n_groups != 2:
                    logger.warning(f"t-test requires 2 groups, got {n_groups}. Skipping.")
                    continue
                if len(group_vals[0]) > 0 and len(group_vals[1]) > 0:
                    stat, pval = ttest_ind(group_vals[0], group_vals[1])
            
            elif method == 'anova':
                if n_groups < 2:
                    logger.warning(f"ANOVA requires ≥2 groups")
                    continue
                valid_groups = [v for v in group_vals if len(v) > 0]
                if len(valid_groups) >= 2:
                    stat, pval = f_oneway(*valid_groups)
            
            elif method == 'kruskal':
                if n_groups < 2:
                    logger.warning(f"Kruskal-Wallis requires ≥2 groups")
                    continue
                valid_groups = [v for v in group_vals if len(v) > 0]
                if len(valid_groups) >= 2:
                    stat, pval = kruskal(*valid_groups)
            
            results.append({
                'pathway': pathway,
                'statistic': stat,
                'pvalue': pval,
                'group': clinical_var,
                'n_groups': n_groups,
                'group_labels': ','.join([str(g) for g in groups]),
                'group_sizes': ','.join([str(s) for s in group_sizes])
            })
        
        # Convert to DataFrame
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction (FDR-Benjamini-Hochberg)
        logger.info(f"  Applying FDR correction ({len(results_df)} tests)")
        
        valid_pvals = ~results_df['pvalue'].isna()
        reject = np.zeros(len(results_df), dtype=bool)
        qval = np.ones(len(results_df))
        
        if valid_pvals.sum() > 0:
            reject[valid_pvals], qval[valid_pvals], _, _ = multipletests(
                results_df.loc[valid_pvals, 'pvalue'],
                method='fdr_bh'
            )
        
        results_df['padj'] = qval
        results_df['significant'] = reject
        
        # Sort by p-value
        results_df = results_df.sort_values('pvalue')
        
        n_sig = results_df['significant'].sum()
        logger.info(f"  Found {n_sig} significant pathways (FDR < 0.05)")
        
        self.results[clinical_var] = results_df
        return results_df
    
    def get_significant_pathways(self, clinical_var: str, padj_threshold=0.05) -> List[str]:
        """Get list of significantly different pathways"""
        if clinical_var not in self.results:
            raise ValueError(f"No results for {clinical_var}. Run compare_by_group first.")
        
        sig_df = self.results[clinical_var]
        sig_pathways = sig_df[sig_df['padj'] < padj_threshold]['pathway'].tolist()
        
        return sig_pathways
    
    def save_results(self, clinical_var: str, output_prefix: str):
        """Save results to CSV"""
        if clinical_var not in self.results:
            raise ValueError(f"No results for {clinical_var}")
        
        results_df = self.results[clinical_var]
        output_file = f"{output_prefix}_{clinical_var}_results.csv"
        results_df.to_csv(output_file, index=False)
        logger.info(f"✓ Saved results to {output_file}")
    
    def compare_all_clinical_vars(self, clinical_vars: List[str]) -> Dict[str, pd.DataFrame]:
        """Run differential analysis for multiple clinical variables"""
        results_dict = {}
        
        for var in clinical_vars:
            if var in self.clinical_data.columns:
                logger.info(f"\nAnalyzing '{var}'...")
                results_dict[var] = self.compare_by_group(var, method='auto')
            else:
                logger.warning(f"Variable '{var}' not found in clinical data")
        
        return results_dict
    
    def summary_table(self) -> pd.DataFrame:
        """Create summary of all analyses"""
        summary_rows = []
        
        for clinical_var, results_df in self.results.items():
            n_sig = results_df['significant'].sum()
            top_pathway = results_df.iloc[0]['pathway'] if len(results_df) > 0 else "N/A"
            top_pval = results_df.iloc[0]['pvalue'] if len(results_df) > 0 else np.nan
            
            summary_rows.append({
                'clinical_variable': clinical_var,
                'n_pathways_tested': len(results_df),
                'n_significant_pathways': n_sig,
                'pct_significant': 100 * n_sig / len(results_df) if len(results_df) > 0 else 0,
                'top_pathway': top_pathway,
                'top_pvalue': top_pval
            })
        
        return pd.DataFrame(summary_rows)


if __name__ == "__main__":
    import sys
    
    logger.info("Testing DifferentialPathwayAnalysis with TCGA-COAD data")
    
    # Load pathway activity
    logger.info("Loading pathway activity...")
    pathway_activity = pd.read_csv('pathway_activity_matrix.csv', index_col=0)
    logger.info(f"Loaded: {pathway_activity.shape}")
    
    # Load clinical data
    logger.info("Loading clinical data...")
    clinical_df = pd.read_csv('TCGA-COAD/filtered_clinical.csv', index_col='case_submitter_id')
    
    # Add age groups
    clinical_df['age_group'] = pd.cut(
        clinical_df['age_at_index'],
        bins=[0, 50, 70, 120],
        labels=['0-50', '50-70', '70+']
    )
    
    logger.info(f"Loaded: {clinical_df.shape}")
    
    # Initialize analyzer
    analyzer = DifferentialPathwayAnalysis(pathway_activity, clinical_df)
    
    # Compare by clinical variables
    logger.info("\n--- Testing differential analysis ---")
    
    # Test age groups
    logger.info("\nAnalyzing by age group...")
    age_results = analyzer.compare_by_group('age_group', method='anova')
    logger.info(f"Age group results (top 5):\n{age_results.head(5)[['pathway', 'pvalue', 'padj', 'significant']].to_string()}")
    
    # Test gender
    logger.info("\nAnalyzing by gender...")
    gender_results = analyzer.compare_by_group('gender', method='ttest')
    logger.info(f"Gender results (top 5):\n{gender_results.head(5)[['pathway', 'pvalue', 'padj', 'significant']].to_string()}")
    
    # Get summary
    logger.info("\nAnalysis Summary:")
    summary = analyzer.summary_table()
    logger.info(summary.to_string())
    
    # Save results
    analyzer.save_results('age_group', 'differential_pathway')
    analyzer.save_results('gender', 'differential_pathway')
    
    logger.info("\n✅ Differential pathway analysis complete!")
