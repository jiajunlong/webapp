"""
Unit Tests for Disease Module Detection (Phase 2)

Tests for:
- DiseaseNetworkBuilder: Gene-disease associations, PPI network, subnetwork extraction
- CommunityDetector: Community detection algorithms, metrics
- ModuleSeparationMetrics: Disease separation, comorbidity prediction
"""

import unittest
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import sys
import tempfile
import os

sys.path.insert(0, '/Users/jaber/Downloads/webapp')

from disease_module_detection import (
    DiseaseNetworkBuilder,
    CommunityDetector,
    ModuleSeparationMetrics
)


class TestDiseaseNetworkBuilder(unittest.TestCase):
    """Test DiseaseNetworkBuilder functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.builder = DiseaseNetworkBuilder('data/gene_disease.tsv')
    
    def test_initialization(self):
        """Test DiseaseNetworkBuilder initialization"""
        self.assertIsNotNone(self.builder)
        self.assertEqual(self.builder.ppi_network, None)
        print("✓ test_initialization passed")
    
    def test_load_gene_disease_associations(self):
        """Test loading gene-disease associations"""
        gene_disease_map = self.builder.load_gene_disease_associations()
        
        self.assertGreater(len(gene_disease_map), 0)
        self.assertGreater(len(self.builder.disease_metadata), 0)
        
        # Check structure
        for disease, genes in list(gene_disease_map.items())[:5]:
            self.assertIsInstance(genes, set)
            self.assertGreater(len(genes), 0)
        
        print(f"✓ test_load_gene_disease_associations passed ({len(gene_disease_map)} diseases)")
    
    def test_load_ppi_network(self):
        """Test PPI network loading"""
        # Must load gene-disease associations first
        self.builder.load_gene_disease_associations()
        ppi = self.builder.load_ppi_network('simple')
        
        self.assertIsNotNone(ppi)
        self.assertGreater(ppi.number_of_nodes(), 0)
        self.assertGreater(ppi.number_of_edges(), 0)
        self.assertEqual(self.builder.ppi_source, 'simple')
        
        print(f"✓ test_load_ppi_network passed ({ppi.number_of_nodes()} nodes)")
    
    def test_build_disease_subnetwork(self):
        """Test disease subnetwork extraction"""
        self.builder.load_gene_disease_associations()
        self.builder.load_ppi_network('simple')
        
        # Get a disease
        disease = list(self.builder.gene_disease_map.keys())[0]
        subnetwork = self.builder.build_disease_subnetwork(disease)
        
        self.assertIsNotNone(subnetwork)
        self.assertGreater(subnetwork.number_of_nodes(), 0)
        self.assertLessEqual(subnetwork.number_of_nodes(), 
                            len(self.builder.gene_disease_map[disease]))
        
        print(f"✓ test_build_disease_subnetwork passed")
    
    def test_build_all_disease_subnetworks(self):
        """Test building subnetworks for all diseases"""
        self.builder.load_gene_disease_associations()
        self.builder.load_ppi_network('simple')
        
        # Use limited diseases for speed
        limited_diseases = list(self.builder.gene_disease_map.keys())[:50]
        
        for disease in limited_diseases:
            self.builder.build_disease_subnetwork(disease)
        
        self.assertGreater(len(self.builder.disease_subnetworks), 0)
        
        print(f"✓ test_build_all_disease_subnetworks passed ({len(self.builder.disease_subnetworks)} diseases)")
    
    def test_get_disease_connectivity_stats(self):
        """Test network statistics computation"""
        self.builder.load_gene_disease_associations()
        self.builder.load_ppi_network('simple')
        
        disease = list(self.builder.gene_disease_map.keys())[0]
        self.builder.build_disease_subnetwork(disease)
        
        stats = self.builder.get_disease_connectivity_stats(disease)
        
        self.assertIn('n_nodes', stats)
        self.assertIn('density', stats)
        self.assertGreaterEqual(stats['density'], 0)
        self.assertLessEqual(stats['density'], 1)
        
        print(f"✓ test_get_disease_connectivity_stats passed")


class TestCommunityDetector(unittest.TestCase):
    """Test CommunityDetector functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.detector = CommunityDetector()
        
        # Create a simple test network with clear communities
        self.test_network = nx.Graph()
        
        # Community 1
        self.test_network.add_edges_from([
            ('A', 'B'), ('B', 'C'), ('C', 'A'),
            ('A', 'D'), ('B', 'D')
        ])
        
        # Community 2
        self.test_network.add_edges_from([
            ('X', 'Y'), ('Y', 'Z'), ('Z', 'X'),
            ('X', 'W'), ('Y', 'W')
        ])
        
        # Few connections between communities
        self.test_network.add_edge('D', 'X')
    
    def test_initialization(self):
        """Test CommunityDetector initialization"""
        self.assertIsNotNone(self.detector)
        self.assertEqual(len(self.detector.communities), 0)
        print("✓ test_initialization passed")
    
    def test_detect_communities_greedy(self):
        """Test greedy community detection"""
        communities = self.detector.detect_communities_greedy(self.test_network)
        
        self.assertGreater(len(communities), 0)
        
        # All nodes should be in some community
        all_nodes = set()
        for comm in communities:
            all_nodes.update(comm)
        
        self.assertEqual(all_nodes, set(self.test_network.nodes()))
        
        print(f"✓ test_detect_communities_greedy passed ({len(communities)} communities)")
    
    def test_detect_communities_label_propagation(self):
        """Test label propagation community detection"""
        communities = self.detector.detect_communities_label_propagation(self.test_network)
        
        self.assertGreater(len(communities), 0)
        
        print(f"✓ test_detect_communities_label_propagation passed ({len(communities)} communities)")
    
    def test_compute_community_metrics(self):
        """Test community metrics computation"""
        communities = self.detector.detect_communities_greedy(self.test_network)
        metrics = self.detector.compute_community_metrics(self.test_network, communities)
        
        self.assertGreater(len(metrics), 0)
        self.assertIn('n_nodes', metrics.columns)
        self.assertIn('density', metrics.columns)
        
        # All densities should be valid
        for density in metrics['density']:
            self.assertGreaterEqual(density, 0)
            self.assertLessEqual(density, 1)
        
        print(f"✓ test_compute_community_metrics passed")


class TestModuleSeparationMetrics(unittest.TestCase):
    """Test ModuleSeparationMetrics functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        # Create two simple networks
        self.ppi = nx.Graph()
        self.ppi.add_edges_from([
            (1, 2), (2, 3), (3, 4), (4, 5),
            (6, 7), (7, 8), (8, 9), (9, 10),
            (5, 6)  # Connection between networks
        ])
        
        disease_subnetworks = {
            'Disease_A': self.ppi.subgraph([1, 2, 3, 4, 5]).copy(),
            'Disease_B': self.ppi.subgraph([6, 7, 8, 9, 10]).copy(),
        }
        
        self.calculator = ModuleSeparationMetrics(self.ppi, disease_subnetworks)
    
    def test_initialization(self):
        """Test ModuleSeparationMetrics initialization"""
        self.assertIsNotNone(self.calculator)
        self.assertEqual(self.calculator.ppi_network.number_of_nodes(), 10)
        print("✓ test_initialization passed")
    
    def test_compute_network_separation(self):
        """Test network separation computation"""
        separation = self.calculator.compute_network_separation('Disease_A', 'Disease_B')
        
        self.assertIsNotNone(separation)
        self.assertGreater(separation, 0)
        self.assertLess(separation, np.inf)
        
        print(f"✓ test_compute_network_separation passed (separation={separation:.3f})")
    
    def test_compute_all_disease_pairs(self):
        """Test pairwise disease separation computation"""
        df = self.calculator.compute_all_disease_pairs()
        
        self.assertGreater(len(df), 0)
        self.assertIn('disease1', df.columns)
        self.assertIn('disease2', df.columns)
        self.assertIn('separation', df.columns)
        
        print(f"✓ test_compute_all_disease_pairs passed ({len(df)} pairs)")
    
    def test_compute_comorbidity_scores(self):
        """Test comorbidity score computation"""
        separation_df = self.calculator.compute_all_disease_pairs()
        comorbidity_df = self.calculator.compute_comorbidity_scores(separation_df)
        
        self.assertGreater(len(comorbidity_df), 0)
        self.assertIn('comorbidity', comorbidity_df.columns)
        
        # All comorbidity scores should be in [0, 1]
        for score in comorbidity_df['comorbidity']:
            self.assertGreaterEqual(score, 0)
            self.assertLessEqual(score, 1)
        
        print(f"✓ test_compute_comorbidity_scores passed")


class TestIntegration(unittest.TestCase):
    """Integration tests for Phase 2"""
    
    def test_full_pipeline(self):
        """Test full disease module detection pipeline"""
        # Initialize
        builder = DiseaseNetworkBuilder('data/gene_disease.tsv')
        builder.load_gene_disease_associations()
        builder.load_ppi_network('simple')
        
        # Build subnetworks for first 20 diseases
        diseases = list(builder.gene_disease_map.keys())[:20]
        for disease in diseases:
            builder.build_disease_subnetwork(disease)
        
        # Compute statistics
        stats = builder.compute_network_statistics_summary()
        self.assertGreater(len(stats), 0)
        
        # Run community detection on a sample
        detector = CommunityDetector()
        sample_disease = diseases[5]
        sample_network = builder.disease_subnetworks.get(sample_disease)
        
        if sample_network and sample_network.number_of_nodes() > 2:
            communities = detector.detect_communities_greedy(sample_network)
            self.assertGreater(len(communities), 0)
        
        print("✓ test_full_pipeline passed")


# ============================================================================
# Test Suite
# ============================================================================

def run_tests():
    """Run all tests and print summary"""
    print("=" * 70)
    print("PHASE 2: DISEASE MODULE DETECTION - TEST SUITE")
    print("=" * 70)
    
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all tests
    suite.addTests(loader.loadTestsFromTestCase(TestDiseaseNetworkBuilder))
    suite.addTests(loader.loadTestsFromTestCase(TestCommunityDetector))
    suite.addTests(loader.loadTestsFromTestCase(TestModuleSeparationMetrics))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegration))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {100 * (result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun:.1f}%")
    print("=" * 70)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
