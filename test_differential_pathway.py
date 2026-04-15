import pandas as pd
import numpy as np
import logging
from differential_pathway_analysis import DifferentialPathwayAnalysis

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("Testing DifferentialPathwayAnalysis with TCGA-COAD data")

# Load pathway activity
logger.info("Loading pathway activity...")
pathway_activity = pd.read_csv('pathway_activity_matrix.csv', index_col=0)
logger.info(f"Loaded pathway activity: {pathway_activity.shape}")

# Load clinical data - use full clinical.tsv and filter to 255 samples
logger.info("Loading clinical data...")
clinical_full = pd.read_csv('TCGA-COAD/clinical.tsv', sep='\t')
clinical_filtered = pd.read_csv('TCGA-COAD/filtered_clinical.csv')

# Get sample IDs that are in pathway activity
samples_in_activity = set(pathway_activity.columns)
samples_in_filtered = set(clinical_filtered['case_submitter_id'])
overlapping_samples = samples_in_activity & samples_in_filtered

logger.info(f"Overlapping samples: {len(overlapping_samples)}")

# Filter to overlapping samples
clinical_df = clinical_filtered[clinical_filtered['case_submitter_id'].isin(overlapping_samples)].copy()
clinical_df = clinical_df.set_index('case_submitter_id')

# Get age from full clinical data
clinical_full_map = {}
for _, row in clinical_full.iterrows():
    if pd.notna(row.get('case_submitter_id')) and pd.notna(row.get('age_at_index')):
        clinical_full_map[row['case_submitter_id']] = row['age_at_index']

# Add age to filtered data
clinical_df['age_at_index'] = clinical_df.index.map(clinical_full_map)

# Add age groups (only for samples with age data)
clinical_df_with_age = clinical_df.dropna(subset=['age_at_index'])
clinical_df_with_age['age_group'] = pd.cut(
    clinical_df_with_age['age_at_index'],
    bins=[0, 50, 70, 150],
    labels=['0-50', '50-70', '70+']
)

logger.info(f"Clinical data shape: {clinical_df.shape}")
logger.info(f"Clinical data with age: {clinical_df_with_age.shape}")
logger.info(f"Available clinical variables: {clinical_df.columns.tolist()}")

# Initialize analyzer
analyzer = DifferentialPathwayAnalysis(pathway_activity, clinical_df_with_age)

# Compare by clinical variables
logger.info("\n--- Differential Pathway Analysis ---")

# Test gender
if 'gender' in clinical_df_with_age.columns:
    logger.info("\nAnalyzing by gender...")
    gender_results = analyzer.compare_by_group('gender', method='auto')
    logger.info(f"Gender analysis - Significant pathways: {gender_results['significant'].sum()}")
    if gender_results['significant'].sum() > 0:
        logger.info("Top significant pathways (by p-value):")
        logger.info(gender_results[gender_results['significant']][['pathway', 'pvalue', 'padj']].head(3).to_string())

# Test age groups
if 'age_group' in clinical_df_with_age.columns:
    logger.info("\nAnalyzing by age group...")
    age_results = analyzer.compare_by_group('age_group', method='auto')
    logger.info(f"Age group analysis - Significant pathways: {age_results['significant'].sum()}")
    if age_results['significant'].sum() > 0:
        logger.info("Top significant pathways (by p-value):")
        logger.info(age_results[age_results['significant']][['pathway', 'pvalue', 'padj']].head(3).to_string())

# Test disease stage
if 'ajcc_pathologic_stage' in clinical_df_with_age.columns:
    logger.info("\nAnalyzing by disease stage...")
    stage_results = analyzer.compare_by_group('ajcc_pathologic_stage', method='auto')
    logger.info(f"Disease stage analysis - Significant pathways: {stage_results['significant'].sum()}")
    if stage_results['significant'].sum() > 0:
        logger.info("Top significant pathways (by p-value):")
        logger.info(stage_results[stage_results['significant']][['pathway', 'pvalue', 'padj']].head(3).to_string())

# Summary
logger.info("\n--- Analysis Summary ---")
summary = analyzer.summary_table()
logger.info(summary.to_string())

# Save results
analyzer.save_results('gender', 'differential_pathway')
if 'age_group' in clinical_df_with_age.columns:
    analyzer.save_results('age_group', 'differential_pathway')
if 'ajcc_pathologic_stage' in clinical_df_with_age.columns:
    analyzer.save_results('ajcc_pathologic_stage', 'differential_pathway')

logger.info("\n✅ Differential pathway analysis complete!")
