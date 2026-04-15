import pandas as pd
import numpy as np
import logging
from differential_pathway_analysis import DifferentialPathwayAnalysis

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("Testing DifferentialPathwayAnalysis - FINAL TEST")

# Load pathway activity
logger.info("Loading pathway activity...")
pathway_activity = pd.read_csv('pathway_activity_matrix.csv', index_col=0)
logger.info(f"Pathway activity shape: {pathway_activity.shape}")

# Load clinical data
logger.info("Loading clinical data...")
clinical_df = pd.read_csv('TCGA-COAD/filtered_clinical.csv')
clinical_df = clinical_df.set_index('case_submitter_id')

# Create mapping from expression sample IDs to case submitter IDs
expr_samples = list(pathway_activity.columns)
expr_to_case = {}
for expr_id in expr_samples:
    case_id = '-'.join(expr_id.split('-')[:3])
    expr_to_case[expr_id] = case_id

# Filter pathway_activity to only include samples in clinical data
common_cases = [expr_to_case[e] for e in expr_samples if expr_to_case[e] in clinical_df.index]
common_expr_samples = [e for e in expr_samples if expr_to_case[e] in common_cases]

# Reindex pathway activity and clinical data with same samples
pathway_activity_matched = pathway_activity[common_expr_samples].copy()
pathway_activity_matched.columns = [expr_to_case[c] for c in common_expr_samples]

clinical_df_matched = clinical_df.loc[pathway_activity_matched.columns]

logger.info(f"Matched samples: {len(pathway_activity_matched.columns)}")
logger.info(f"Pathway activity shape: {pathway_activity_matched.shape}")
logger.info(f"Clinical data shape: {clinical_df_matched.shape}")

# Initialize analyzer
analyzer = DifferentialPathwayAnalysis(pathway_activity_matched, clinical_df_matched)

# Compare by gender
logger.info("\n" + "="*70)
logger.info("Testing by GENDER")
logger.info("="*70)
gender_results = analyzer.compare_by_group('gender', method='auto')
logger.info(f"Total pathways tested: {len(gender_results)}")
logger.info(f"Significant pathways (FDR < 0.05): {gender_results['significant'].sum()}")

if len(gender_results) > 0:
    logger.info("\nTop 10 by p-value:")
    print(gender_results[['pathway', 'pvalue', 'padj', 'significant']].head(10).to_string(index=False))

# Compare by disease stage (if enough variation)
if 'ajcc_pathologic_stage' in clinical_df_matched.columns:
    logger.info("\n" + "="*70)
    logger.info("Testing by DISEASE STAGE")
    logger.info("="*70)
    stage_results = analyzer.compare_by_group('ajcc_pathologic_stage', method='auto')
    logger.info(f"Total pathways tested: {len(stage_results)}")
    logger.info(f"Significant pathways (FDR < 0.05): {stage_results['significant'].sum()}")
    
    if len(stage_results) > 0:
        logger.info("\nTop 10 by p-value:")
        print(stage_results[['pathway', 'pvalue', 'padj', 'significant']].head(10).to_string(index=False))

# Summary
logger.info("\n" + "="*70)
logger.info("ANALYSIS SUMMARY")
logger.info("="*70)
summary = analyzer.summary_table()
print(summary.to_string(index=False))

# Save results
analyzer.save_results('gender', 'phase1_diff_pathway')
if 'ajcc_pathologic_stage' in clinical_df_matched.columns:
    analyzer.save_results('ajcc_pathologic_stage', 'phase1_diff_pathway')

logger.info("\n✅ Differential pathway analysis SUCCESSFUL!")
