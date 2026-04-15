"""
Full test of HubGeneIdentifier on all 347 pathways.
Outputs hub gene summary for all pathways.
"""

import pandas as pd
import sys
import logging
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

sys.path.insert(0, '/Users/jaber/Downloads/webapp')

from hub_gene_identification import HubGeneIdentifier

# Load pathway genes
pathway_file = '/Users/jaber/Downloads/webapp/data/pathway(基因名映射版).tsv'
pathway_df = pd.read_csv(pathway_file, sep='\t')
pathway_dict = {}
for _, row in pathway_df.iterrows():
    pathway_name = row['Pathway_Name']
    gene_str = row['Gene']
    if pd.notna(gene_str):
        genes = [g.strip() for g in gene_str.replace(';', ',').split(',')]
        genes = [g for g in genes if g and g != 'NA']
        if len(genes) > 0:
            pathway_dict[pathway_name] = genes

logger.info(f"Loaded {len(pathway_dict)} pathways")

# Initialize and calculate
identifier = HubGeneIdentifier(pathway_dict, gene_network=None)  # No network this time
identifier.load_expression_data('/Users/jaber/Downloads/webapp/TCGA-COAD/filtered_hiseq_data.csv')

start = time.time()
all_hub_scores = identifier.calculate_all_hub_genes()
elapsed = time.time() - start

logger.info(f"✓ Calculated hub genes for {len(all_hub_scores)} pathways in {elapsed:.1f} seconds")

# Export summary
output_file = '/Users/jaber/Downloads/webapp/phase1_hub_genes_all_pathways.csv'
identifier.export_hub_genes_summary(output_file, top_n=10)

# Check file
df = pd.read_csv(output_file)
logger.info(f"✓ Exported {len(df)} records (top 10 genes × {len(all_hub_scores)} pathways)")
logger.info(f"\nSample records:")
print(df.head(15).to_string(index=False))

logger.info(f"\n✓ SUCCESS - Full hub gene calculation completed")
