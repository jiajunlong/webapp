"""
Test script for pathway_visualizations module

Tests all visualization functions with real TCGA-COAD data.
"""

import pandas as pd
import numpy as np
import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

sys.path.insert(0, '/Users/jaber/Downloads/webapp')

from pathway_visualizations import (
    plot_pathway_activity_heatmap,
    plot_pathway_violin,
    plot_hub_genes_bar,
    plot_differential_pathways,
)
from hub_gene_identification import HubGeneIdentifier


def test_visualizations():
    """Test all visualization functions."""
    print("\n" + "="*70)
    print("PATHWAY VISUALIZATION TESTS")
    print("="*70)
    
    # Load data
    logger.info("Loading test data...")
    
    # Pathway activity
    pathway_activity = pd.read_csv(
        '/Users/jaber/Downloads/webapp/pathway_activity_matrix.csv',
        index_col=0
    )
    logger.info(f"✓ Loaded pathway activity: {pathway_activity.shape}")
    
    # Clinical data
    clinical_data = pd.read_csv(
        '/Users/jaber/Downloads/webapp/TCGA-COAD/filtered_clinical.csv',
        index_col=0
    )
    logger.info(f"✓ Loaded clinical data: {clinical_data.shape}")
    
    # Differential results
    diff_results = pd.read_csv(
        '/Users/jaber/Downloads/webapp/phase1_diff_pathway_gender_results.csv'
    )
    logger.info(f"✓ Loaded differential results: {diff_results.shape}")
    
    # Hub genes
    pathway_file = '/Users/jaber/Downloads/webapp/data/pathway(基因名映射版).tsv'
    pathway_df = pd.read_csv(pathway_file, sep='\t')
    pathway_dict = {}
    for _, row in pathway_df.iterrows():
        name = row['Pathway_Name']
        gene_str = row['Gene']
        if pd.notna(gene_str):
            genes = [g.strip() for g in gene_str.replace(';', ',').split(',')]
            genes = [g for g in genes if g and g != 'NA']
            if len(genes) > 0:
                pathway_dict[name] = genes
    
    identifier = HubGeneIdentifier(pathway_dict)
    identifier.load_expression_data('/Users/jaber/Downloads/webapp/TCGA-COAD/filtered_hiseq_data.csv')
    
    # Calculate hub genes for one pathway (for speed)
    test_pathway = 'Glycolysis / Gluconeogenesis - Homo sapiens (human)'
    hub_genes = identifier.calculate_hub_score(test_pathway)
    logger.info(f"✓ Calculated hub genes for {test_pathway}")
    
    # =====================================================================
    # TEST 1: Heatmap
    # =====================================================================
    print("\n" + "-"*70)
    print("TEST 1: Pathway Activity Heatmap")
    print("-"*70)
    
    try:
        fig = plot_pathway_activity_heatmap(pathway_activity, clinical_data, 'gender')
        logger.info("✓ Created heatmap successfully")
        
        # Save to HTML
        output_file = '/Users/jaber/Downloads/webapp/test_heatmap.html'
        fig.write_html(output_file)
        logger.info(f"✓ Saved to {output_file}")
    except Exception as e:
        logger.error(f"✗ Heatmap test failed: {e}")
    
    # =====================================================================
    # TEST 2: Violin plot
    # =====================================================================
    print("\n" + "-"*70)
    print("TEST 2: Pathway Violin Plot")
    print("-"*70)
    
    try:
        fig = plot_pathway_violin(pathway_activity, clinical_data,
                                 test_pathway, 'gender')
        logger.info("✓ Created violin plot successfully")
        
        output_file = '/Users/jaber/Downloads/webapp/test_violin.html'
        fig.write_html(output_file)
        logger.info(f"✓ Saved to {output_file}")
    except Exception as e:
        logger.error(f"✗ Violin plot test failed: {e}")
    
    # =====================================================================
    # TEST 3: Hub genes bar chart
    # =====================================================================
    print("\n" + "-"*70)
    print("TEST 3: Hub Genes Bar Chart")
    print("-"*70)
    
    try:
        fig = plot_hub_genes_bar(hub_genes, top_n=15)
        logger.info("✓ Created bar chart successfully")
        
        output_file = '/Users/jaber/Downloads/webapp/test_hub_bar.html'
        fig.write_html(output_file)
        logger.info(f"✓ Saved to {output_file}")
    except Exception as e:
        logger.error(f"✗ Bar chart test failed: {e}")
    
    # =====================================================================
    # TEST 4: Differential pathways
    # =====================================================================
    print("\n" + "-"*70)
    print("TEST 4: Differential Pathways Plot")
    print("-"*70)
    
    try:
        fig = plot_differential_pathways(diff_results, top_n=20)
        logger.info("✓ Created differential pathway plot successfully")
        
        output_file = '/Users/jaber/Downloads/webapp/test_diff_pathways.html'
        fig.write_html(output_file)
        logger.info(f"✓ Saved to {output_file}")
    except Exception as e:
        logger.error(f"✗ Differential plot test failed: {e}")
    
    print("\n" + "="*70)
    print("✓ ALL VISUALIZATION TESTS COMPLETED")
    print("="*70 + "\n")


if __name__ == "__main__":
    try:
        test_visualizations()
    except Exception as e:
        logger.error(f"Test suite failed: {e}", exc_info=True)
        sys.exit(1)
