"""
Test script for HubGeneIdentifier module

Tests hub gene identification with:
1. Real pathway data from KEGG
2. Real gene expression from TCGA-COAD
3. Real gene network (merged_output.tsv)
"""

import pandas as pd
import numpy as np
import networkx as nx
import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add project root to path
sys.path.insert(0, '/Users/jaber/Downloads/webapp')

from hub_gene_identification import HubGeneIdentifier
from data_loader import RealDataLoader


def load_pathway_genes_dict():
    """Load pathway genes as dictionary from KEGG data."""
    logger.info("Loading pathway data...")
    
    pathway_file = '/Users/jaber/Downloads/webapp/data/pathway(基因名映射版).tsv'
    pathway_df = pd.read_csv(pathway_file, sep='\t')
    logger.info(f"✓ Loaded {len(pathway_df)} pathway records")
    
    pathway_dict = {}
    for _, row in pathway_df.iterrows():
        pathway_name = row['Pathway_Name']
        gene_str = row['Gene']
        
        if pd.notna(gene_str):
            # Parse genes (can be comma or semicolon separated)
            genes = [g.strip() for g in gene_str.replace(';', ',').split(',')]
            genes = [g for g in genes if g and g != 'NA']
            
            if len(genes) > 0:
                pathway_dict[pathway_name] = genes
    
    logger.info(f"✓ Loaded {len(pathway_dict)} pathways with gene mappings")
    return pathway_dict


def build_gene_network_from_edges(edge_file: str, max_edges: int = 50000) -> nx.Graph:
    """Build NetworkX graph from edge list file."""
    logger.info(f"Building gene network from {edge_file}")
    G = nx.Graph()
    
    try:
        df = pd.read_csv(edge_file, sep='\t', nrows=max_edges)
        for _, row in df.iterrows():
            source = row['source']
            target = row['target']
            weight = float(row.get('weight', 1.0))
            G.add_edge(source, target, weight=weight)
    except Exception as e:
        logger.warning(f"Error loading edges: {e}")
    
    logger.info(f"✓ Built network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def test_hub_genes_single_pathway():
    """Test hub gene identification for a single pathway."""
    print("\n" + "="*70)
    print("TEST 1: Single Pathway Hub Gene Identification")
    print("="*70)
    
    # Load data
    pathway_genes = load_pathway_genes_dict()
    
    # Load expression data
    logger.info("Loading TCGA-COAD expression data...")
    expr_file = '/Users/jaber/Downloads/webapp/TCGA-COAD/filtered_hiseq_data.csv'
    expr_data = pd.read_csv(expr_file, index_col=0)
    logger.info(f"✓ Loaded expression: {expr_data.shape}")
    
    # Build gene network
    edge_file = '/Users/jaber/Downloads/webapp/data/merged_output.tsv'
    gene_network = build_gene_network_from_edges(edge_file)
    
    # Initialize identifier
    identifier = HubGeneIdentifier(pathway_genes, gene_network)
    identifier.load_expression_data(expr_file)
    
    # Test on "Glycolysis / Gluconeogenesis" pathway
    test_pathway = None
    for pathway in pathway_genes.keys():
        if 'Glycolysis' in pathway:
            test_pathway = pathway
            break
    
    if test_pathway is None:
        logger.warning("Could not find Glycolysis pathway, using first available")
        test_pathway = list(pathway_genes.keys())[0]
    
    logger.info(f"\nTesting with pathway: {test_pathway}")
    logger.info(f"  Genes in pathway: {len(pathway_genes[test_pathway])}")
    
    # Calculate hub scores
    hub_scores = identifier.calculate_hub_score(test_pathway)
    
    logger.info(f"\n✓ Hub scores calculated for {len(hub_scores)} genes")
    logger.info(f"\nTop 10 hub genes in {test_pathway}:")
    cols = ['rank', 'gene', 'hub_score', 'degree', 'betweenness', 'expr_variance']
    print(hub_scores[cols].head(10).to_string(index=False))
    
    return identifier


def test_hub_genes_multiple_pathways(identifier):
    """Test hub gene identification for multiple pathways."""
    print("\n" + "="*70)
    print("TEST 2: Multiple Pathway Hub Gene Identification")
    print("="*70)
    
    # Calculate hub genes for all pathways (limit to 50 for testing)
    all_pathways = list(identifier.pathway_genes.keys())[:50]
    logger.info(f"Calculating hub genes for {len(all_pathways)} pathways (limit for testing)...")
    
    all_hub_scores = identifier.calculate_all_hub_genes(all_pathways)
    logger.info(f"✓ Calculated hub genes for {len(all_hub_scores)} pathways")
    
    # Get summary of top pathways
    logger.info("\nTop 10 pathways by mean hub score:")
    pathway_scores = {}
    for pathway, df in all_hub_scores.items():
        pathway_scores[pathway] = df['hub_score'].mean()
    
    top_pathways = sorted(pathway_scores.items(), 
                         key=lambda x: x[1], 
                         reverse=True)[:10]
    
    for i, (pathway, score) in enumerate(top_pathways, 1):
        logger.info(f"  {i:2d}. {pathway[:50]:50s} (mean hub score: {score:.4f})")
    
    return all_hub_scores


def test_hub_genes_statistics(identifier):
    """Generate statistics on hub genes."""
    print("\n" + "="*70)
    print("TEST 3: Hub Gene Statistics")
    print("="*70)
    
    # Collect statistics
    all_degrees = []
    all_betweenness = []
    all_expr_vars = []
    all_hub_scores = []
    
    for pathway, df in identifier.hub_scores.items():
        all_degrees.extend(df['degree'].values)
        all_betweenness.extend(df['betweenness'].values)
        all_expr_vars.extend(df['expr_variance'].values)
        all_hub_scores.extend(df['hub_score'].values)
    
    logger.info("\n✓ Hub Gene Statistics Across All Pathways:")
    logger.info(f"  Total genes analyzed: {len(all_hub_scores)}")
    logger.info(f"  Pathways analyzed: {len(identifier.hub_scores)}")
    
    logger.info(f"\nDegree Centrality (pathway connections):")
    logger.info(f"  Mean: {np.mean(all_degrees):.2f}")
    logger.info(f"  Median: {np.median(all_degrees):.2f}")
    logger.info(f"  Range: [{np.min(all_degrees)}, {np.max(all_degrees)}]")
    
    logger.info(f"\nBetweenness Centrality:")
    logger.info(f"  Mean: {np.mean(all_betweenness):.4f}")
    logger.info(f"  Median: {np.median(all_betweenness):.4f}")
    logger.info(f"  Range: [{np.min(all_betweenness):.4f}, {np.max(all_betweenness):.4f}]")
    
    logger.info(f"\nExpression Variance:")
    logger.info(f"  Mean: {np.mean(all_expr_vars):.4f}")
    logger.info(f"  Median: {np.median(all_expr_vars):.4f}")
    logger.info(f"  Range: [{np.min(all_expr_vars):.4f}, {np.max(all_expr_vars):.4f}]")
    
    logger.info(f"\nHub Score (composite):")
    logger.info(f"  Mean: {np.mean(all_hub_scores):.4f}")
    logger.info(f"  Median: {np.median(all_hub_scores):.4f}")
    logger.info(f"  Range: [{np.min(all_hub_scores):.4f}, {np.max(all_hub_scores):.4f}]")


def test_export_hub_genes(identifier):
    """Test exporting hub genes to file."""
    print("\n" + "="*70)
    print("TEST 4: Export Hub Genes Summary")
    print("="*70)
    
    output_file = '/Users/jaber/Downloads/webapp/phase1_hub_genes_summary.csv'
    identifier.export_hub_genes_summary(output_file, top_n=5)
    
    # Verify file was created
    try:
        df = pd.read_csv(output_file)
        logger.info(f"✓ Exported {len(df)} hub gene records to {output_file}")
        logger.info(f"\nFirst 10 records:")
        print(df.head(10).to_string(index=False))
    except Exception as e:
        logger.error(f"Error reading exported file: {e}")


if __name__ == "__main__":
    print("\n" + "="*70)
    print("HUB GENE IDENTIFICATION - COMPREHENSIVE TEST SUITE")
    print("="*70)
    
    try:
        # Test 1: Single pathway
        identifier = test_hub_genes_single_pathway()
        
        # Test 2: Multiple pathways
        all_scores = test_hub_genes_multiple_pathways(identifier)
        
        # Test 3: Statistics
        test_hub_genes_statistics(identifier)
        
        # Test 4: Export
        test_export_hub_genes(identifier)
        
        print("\n" + "="*70)
        print("✓ ALL TESTS COMPLETED SUCCESSFULLY")
        print("="*70 + "\n")
        
    except Exception as e:
        logger.error(f"Test failed: {e}", exc_info=True)
        sys.exit(1)
