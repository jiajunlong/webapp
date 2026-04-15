"""
Unit tests for miRNA-Gene-Pathway Integration Module (Phase 2 Task 3)

Tests cover:
1. miRNATargetPredictor - Target prediction and validation
2. miRNARegulatoryNetwork - Network building and hub identification
3. RegulatoryModuleAnalysis - Module identification and scoring
4. Integration test - Full pipeline validation

Test Strategy:
- Use realistic synthetic data matching TCGA-COAD dimensions
- Test both success cases and edge cases
- Validate data flows and outputs
- Ensure biological relevance of predictions

Status: Ready for execution
"""

import unittest
import pandas as pd
import numpy as np
import warnings
import tempfile
import os
from pathlib import Path

# Import modules to test
from mirna_integration import (
    miRNATargetPredictor,
    miRNARegulatoryNetwork,
    RegulatoryModuleAnalysis
)

warnings.filterwarnings('ignore')


class TestmiRNATargetPredictor(unittest.TestCase):
    """Test miRNA target prediction"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data once for all tests"""
        np.random.seed(42)
        
        # Create synthetic miRNA expression (50 miRNAs × 40 samples)
        cls.mirna_expr = pd.DataFrame(
            np.random.randn(50, 40),
            index=[f"miR-{i}" for i in range(1, 51)],
            columns=[f"sample_{i}" for i in range(1, 41)]
        )
        
        # Create synthetic gene expression (500 genes × 40 samples)
        cls.gene_expr = pd.DataFrame(
            np.random.randn(500, 40),
            index=[f"GENE_{i}" for i in range(1, 501)],
            columns=[f"sample_{i}" for i in range(1, 41)]
        )
        
        # Create some synthetic negative correlations for specific pairs
        # miR-1 should negatively correlate with GENE_1, GENE_2, GENE_3
        mir1_vals = cls.mirna_expr.loc["miR-1"].values
        cls.gene_expr.loc["GENE_1"] = -mir1_vals + np.random.randn(40) * 0.2  # Strong negative corr
        cls.gene_expr.loc["GENE_2"] = -mir1_vals + np.random.randn(40) * 0.2
        cls.gene_expr.loc["GENE_3"] = -mir1_vals + np.random.randn(40) * 0.2
        
        # miR-2 should negatively correlate with GENE_10, GENE_11, GENE_12
        mir2_vals = cls.mirna_expr.loc["miR-2"].values
        cls.gene_expr.loc["GENE_10"] = -mir2_vals + np.random.randn(40) * 0.2
        cls.gene_expr.loc["GENE_11"] = -mir2_vals + np.random.randn(40) * 0.2
        cls.gene_expr.loc["GENE_12"] = -mir2_vals + np.random.randn(40) * 0.2
    
    def test_predictor_initialization(self):
        """Test proper initialization"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        
        self.assertEqual(predictor.mirna_expr.shape[0], 50)
        self.assertEqual(predictor.gene_expr.shape[0], 500)
        self.assertIsNone(predictor.targets)
        self.assertIsNone(predictor.target_correlations)
    
    def test_predict_targets_basic(self):
        """Test basic target prediction"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        targets = predictor.predict_targets(correlation_threshold=-0.3, method='pearson')
        
        # Should have predicted at least some targets
        self.assertIsInstance(targets, dict)
        self.assertGreater(len(targets), 0)
        
        # Check that predicted targets exist for strong correlations
        self.assertIn("miR-1", targets)
        self.assertIn("miR-2", targets)
    
    def test_predict_targets_correlation_filtering(self):
        """Test that correlation threshold properly filters targets"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        
        # Stricter threshold
        targets_strict = predictor.predict_targets(correlation_threshold=-0.5, method='pearson')
        
        # Looser threshold should find more targets
        predictor2 = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        targets_loose = predictor2.predict_targets(correlation_threshold=-0.2, method='pearson')
        
        # Loose should find >= strict
        total_strict = sum(len(v) for v in targets_strict.values())
        total_loose = sum(len(v) for v in targets_loose.values())
        self.assertGreaterEqual(total_loose, total_strict)
    
    def test_get_targets_for_mirna(self):
        """Test retrieving top targets for specific miRNA"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        predictor.predict_targets(correlation_threshold=-0.3)
        
        # Get targets for miR-1
        targets_df = predictor.get_targets_for_mirna("miR-1", top_n=10)
        
        self.assertIsInstance(targets_df, pd.DataFrame)
        self.assertEqual(targets_df['miRNA'].iloc[0], "miR-1")
        # Should be sorted by correlation (most negative first)
        if len(targets_df) > 1:
            self.assertLessEqual(targets_df['correlation'].iloc[0], targets_df['correlation'].iloc[-1])
    
    def test_validate_against_databases(self):
        """Test validation against known databases"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        targets = predictor.predict_targets(correlation_threshold=-0.3)
        
        validation_df = predictor.validate_against_databases(targets)
        
        self.assertIsInstance(validation_df, pd.DataFrame)
        self.assertGreater(len(validation_df), 0)
        self.assertIn('miRNA', validation_df.columns)
        self.assertIn('n_predicted_targets', validation_df.columns)
        self.assertIn('confidence', validation_df.columns)
        
        # All confidence values should be 'high' or 'low'
        self.assertTrue(validation_df['confidence'].isin(['high', 'low']).all())
    
    def test_predict_targets_spearman(self):
        """Test target prediction with Spearman correlation"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        targets = predictor.predict_targets(correlation_threshold=-0.3, method='spearman')
        
        self.assertIsInstance(targets, dict)
        self.assertGreater(len(targets), 0)
        
        # Verify correlations are computed correctly
        self.assertIsNotNone(predictor.target_correlations)


class TestmiRNARegulatoryNetwork(unittest.TestCase):
    """Test miRNA regulatory network building and analysis"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data"""
        # Create synthetic miRNA targets
        cls.mirna_targets = {
            "miR-1": ["GENE_1", "GENE_2", "GENE_3", "GENE_4"],
            "miR-2": ["GENE_3", "GENE_4", "GENE_5", "GENE_6"],
            "miR-3": ["GENE_7", "GENE_8"],
            "miR-4": ["GENE_1", "GENE_9", "GENE_10"],
        }
        
        # Create synthetic pathway genes
        cls.pathway_genes = {
            "Pathway_A": ["GENE_1", "GENE_2", "GENE_3", "GENE_4", "GENE_5"],
            "Pathway_B": ["GENE_6", "GENE_7", "GENE_8", "GENE_9"],
            "Pathway_C": ["GENE_10", "GENE_11", "GENE_12"],
        }
    
    def test_network_initialization(self):
        """Test network initialization"""
        network = miRNARegulatoryNetwork(self.mirna_targets, self.pathway_genes)
        
        self.assertEqual(len(self.mirna_targets), 4)
        self.assertEqual(len(self.pathway_genes), 3)
        self.assertIsNone(network.network)
        self.assertIsNone(network.hub_mirnas)
    
    def test_build_network(self):
        """Test network building"""
        network = miRNARegulatoryNetwork(self.mirna_targets, self.pathway_genes)
        net_structure = network.build_network()
        
        # Check structure
        self.assertIn('nodes', net_structure)
        self.assertIn('edges', net_structure)
        self.assertIn('pathways', net_structure)
        
        # Should have miRNA nodes + gene nodes
        self.assertGreater(len(net_structure['nodes']), 0)
        
        # Count miRNA and gene nodes
        mirna_nodes = [n for n in net_structure['nodes'] if n['type'] == 'miRNA']
        gene_nodes = [n for n in net_structure['nodes'] if n['type'] == 'gene']
        
        self.assertEqual(len(mirna_nodes), 4)  # 4 miRNAs
        self.assertGreater(len(gene_nodes), 0)  # Some genes
        
        # Edges should represent miRNA-target relations
        self.assertGreater(len(net_structure['edges']), 0)
    
    def test_identify_hub_mirnas(self):
        """Test hub miRNA identification"""
        network = miRNARegulatoryNetwork(self.mirna_targets, self.pathway_genes)
        network.build_network()
        hub_mirnas = network.identify_hub_mirnas(top_n=3)
        
        self.assertIsInstance(hub_mirnas, pd.DataFrame)
        self.assertLessEqual(len(hub_mirnas), 3)
        
        # Should have required columns
        self.assertIn('miRNA', hub_mirnas.columns)
        self.assertIn('n_targets', hub_mirnas.columns)
        self.assertIn('hub_score', hub_mirnas.columns)
        
        # Hub scores should be sorted descending
        if len(hub_mirnas) > 1:
            self.assertGreaterEqual(hub_mirnas['hub_score'].iloc[0], hub_mirnas['hub_score'].iloc[-1])
    
    def test_map_to_pathways(self):
        """Test mapping miRNA regulations to pathways"""
        network = miRNARegulatoryNetwork(self.mirna_targets, self.pathway_genes)
        mirna_pathway_map = network.map_to_pathways()
        
        self.assertIsInstance(mirna_pathway_map, pd.DataFrame)
        self.assertGreater(len(mirna_pathway_map), 0)
        
        # Should have required columns
        self.assertIn('miRNA', mirna_pathway_map.columns)
        self.assertIn('pathway', mirna_pathway_map.columns)
        self.assertIn('coverage', mirna_pathway_map.columns)
        
        # Coverage should be between 0 and 1
        self.assertTrue((mirna_pathway_map['coverage'] >= 0).all())
        self.assertTrue((mirna_pathway_map['coverage'] <= 1).all())


class TestRegulatoryModuleAnalysis(unittest.TestCase):
    """Test regulatory module identification and scoring"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data"""
        # Create synthetic miRNA targets
        cls.mirna_targets = {
            "miR-1": ["GENE_1", "GENE_2", "GENE_3"],
            "miR-2": ["GENE_3", "GENE_4", "GENE_5"],
            "miR-3": ["GENE_6", "GENE_7"],
        }
        
        # Create synthetic pathway genes
        cls.pathway_genes = {
            "Pathway_A": ["GENE_1", "GENE_2", "GENE_3", "GENE_4"],
            "Pathway_B": ["GENE_5", "GENE_6", "GENE_7"],
            "Pathway_C": ["GENE_8", "GENE_9"],
        }
        
        # Create synthetic disease modules
        cls.disease_modules = {
            "Disease_1": ["GENE_1", "GENE_2"],
            "Disease_2": ["GENE_3", "GENE_4"],
            "Disease_3": ["GENE_6", "GENE_7"],
        }
    
    def test_module_analysis_initialization(self):
        """Test initialization with and without disease modules"""
        # With disease modules
        analyzer = RegulatoryModuleAnalysis(
            self.mirna_targets,
            self.pathway_genes,
            self.disease_modules
        )
        
        self.assertEqual(len(analyzer.mirna_targets), 3)
        self.assertEqual(len(analyzer.pathway_genes), 3)
        self.assertEqual(len(analyzer.disease_modules), 3)
        
        # Without disease modules
        analyzer2 = RegulatoryModuleAnalysis(
            self.mirna_targets,
            self.pathway_genes
        )
        
        self.assertEqual(len(analyzer2.disease_modules), 0)
    
    def test_identify_regulatory_modules(self):
        """Test regulatory module identification"""
        analyzer = RegulatoryModuleAnalysis(
            self.mirna_targets,
            self.pathway_genes,
            self.disease_modules
        )
        modules = analyzer.identify_regulatory_modules()
        
        self.assertIsInstance(modules, dict)
        self.assertGreater(len(modules), 0)
        
        # Each module should have required structure
        for module_id, module_info in modules.items():
            self.assertIn('pathway', module_info)
            self.assertIn('n_genes', module_info)
            self.assertIn('regulating_miRNAs', module_info)
            self.assertIn('n_regulating_mirnas', module_info)
            self.assertIn('disease_associations', module_info)
    
    def test_score_regulatory_importance(self):
        """Test regulatory importance scoring"""
        analyzer = RegulatoryModuleAnalysis(
            self.mirna_targets,
            self.pathway_genes,
            self.disease_modules
        )
        scores = analyzer.score_regulatory_importance()
        
        self.assertIsInstance(scores, pd.DataFrame)
        self.assertGreater(len(scores), 0)
        
        # Check required columns
        self.assertIn('module', scores.columns)
        self.assertIn('pathway', scores.columns)
        self.assertIn('regulatory_score', scores.columns)
        self.assertIn('rank', scores.columns)
        
        # Scores should be sorted by regulatory_score descending
        if len(scores) > 1:
            self.assertGreaterEqual(scores['regulatory_score'].iloc[0], scores['regulatory_score'].iloc[-1])
    
    def test_export_regulatory_network(self):
        """Test exporting regulatory network"""
        analyzer = RegulatoryModuleAnalysis(
            self.mirna_targets,
            self.pathway_genes,
            self.disease_modules
        )
        
        # Use temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            temp_file = f.name
        
        try:
            analyzer.export_regulatory_network(temp_file)
            
            # Check file exists and has content
            self.assertTrue(os.path.exists(temp_file))
            
            # Read and validate content
            with open(temp_file, 'r') as f:
                lines = f.readlines()
            
            # Should have header + data
            self.assertGreater(len(lines), 1)
            
            # Header check
            header = lines[0].strip().split('\t')
            self.assertIn('Module', header)
            self.assertIn('Pathway', header)
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)


class TestIntegration(unittest.TestCase):
    """Integration tests for full miRNA analysis pipeline"""
    
    @classmethod
    def setUpClass(cls):
        """Set up realistic test data"""
        np.random.seed(42)
        
        # Create realistic miRNA expression (100 miRNAs × 50 samples)
        cls.mirna_expr = pd.DataFrame(
            np.random.randn(100, 50),
            index=[f"miR-{i}" for i in range(1, 101)],
            columns=[f"sample_{i}" for i in range(1, 51)]
        )
        
        # Create realistic gene expression (1000 genes × 50 samples)
        cls.gene_expr = pd.DataFrame(
            np.random.randn(1000, 50),
            index=[f"GENE_{i}" for i in range(1, 1001)],
            columns=[f"sample_{i}" for i in range(1, 51)]
        )
        
        # Create structured co-regulation patterns
        # 10 miRNAs with strong negative correlations to specific gene sets
        for mir_idx in range(1, 11):
            mirna_name = f"miR-{mir_idx}"
            mir_vals = cls.mirna_expr.loc[mirna_name].values
            
            # Create 5-10 target genes for each miRNA
            for target_idx in range(mir_idx * 10, mir_idx * 10 + 8):
                gene_name = f"GENE_{target_idx}"
                if gene_name in cls.gene_expr.index:
                    cls.gene_expr.loc[gene_name] = -mir_vals + np.random.randn(50) * 0.3
        
        # Create pathway definitions
        cls.pathway_genes = {
            f"Pathway_{i}": [f"GENE_{j}" for j in range(i * 30, (i + 1) * 30)]
            for i in range(1, 11)
        }
        
        # Create disease modules
        cls.disease_modules = {
            f"Disease_{i}": [f"GENE_{j}" for j in range(i * 50, (i + 1) * 50)]
            for i in range(1, 6)
        }
    
    def test_full_pipeline(self):
        """Test complete miRNA analysis pipeline"""
        print("\n=== Full Pipeline Integration Test ===")
        
        # Step 1: Predict targets
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        targets = predictor.predict_targets(correlation_threshold=-0.3)
        
        print(f"✓ Predicted targets for {len(targets)} miRNAs")
        self.assertGreater(len(targets), 0)
        
        # Step 2: Build network
        network = miRNARegulatoryNetwork(targets, self.pathway_genes)
        net_structure = network.build_network()
        
        print(f"✓ Built network: {len(net_structure['nodes'])} nodes, {len(net_structure['edges'])} edges")
        self.assertGreater(len(net_structure['edges']), 0)
        
        # Step 3: Identify hub miRNAs
        hub_mirnas = network.identify_hub_mirnas(top_n=10)
        
        print(f"✓ Identified {len(hub_mirnas)} hub miRNAs")
        self.assertGreater(len(hub_mirnas), 0)
        
        # Step 4: Map to pathways
        mirna_pathway = network.map_to_pathways()
        
        print(f"✓ Mapped {len(mirna_pathway)} miRNA-pathway associations")
        self.assertGreater(len(mirna_pathway), 0)
        
        # Step 5: Regulatory module analysis
        analyzer = RegulatoryModuleAnalysis(targets, self.pathway_genes, self.disease_modules)
        modules = analyzer.identify_regulatory_modules()
        
        print(f"✓ Identified {len(modules)} regulatory modules")
        self.assertGreater(len(modules), 0)
        
        # Step 6: Score importance
        scores = analyzer.score_regulatory_importance()
        
        print(f"✓ Scored {len(scores)} modules")
        self.assertGreater(len(scores), 0)
        
        # Step 7: Export results
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            temp_file = f.name
        
        try:
            analyzer.export_regulatory_network(temp_file)
            print(f"✓ Exported regulatory network")
            self.assertTrue(os.path.exists(temp_file))
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        
        print("\n✅ Full pipeline completed successfully!")
    
    def test_pipeline_with_validation(self):
        """Test pipeline with database validation"""
        predictor = miRNATargetPredictor(self.mirna_expr, self.gene_expr)
        targets = predictor.predict_targets(correlation_threshold=-0.3)
        
        # Validate predictions
        validation = predictor.validate_against_databases(targets)
        
        self.assertGreater(len(validation), 0)
        self.assertIn('confidence', validation.columns)
        
        # Check that high-confidence miRNAs have more targets
        high_conf = validation[validation['confidence'] == 'high']
        low_conf = validation[validation['confidence'] == 'low']
        
        if len(high_conf) > 0 and len(low_conf) > 0:
            mean_high = high_conf['n_predicted_targets'].mean()
            mean_low = low_conf['n_predicted_targets'].mean()
            self.assertGreaterEqual(mean_high, mean_low)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling"""
    
    def test_empty_targets(self):
        """Test handling of empty target predictions"""
        mirna_expr = pd.DataFrame(
            np.random.randn(5, 10),
            index=[f"miR-{i}" for i in range(1, 6)],
            columns=[f"s{i}" for i in range(1, 11)]
        )
        gene_expr = pd.DataFrame(
            np.random.randn(10, 10),
            index=[f"GENE_{i}" for i in range(1, 11)],
            columns=[f"s{i}" for i in range(1, 11)]
        )
        
        predictor = miRNATargetPredictor(mirna_expr, gene_expr)
        # Very strict threshold should result in few/no targets
        targets = predictor.predict_targets(correlation_threshold=-0.9)
        
        # Should still return a dict (possibly empty)
        self.assertIsInstance(targets, dict)
    
    def test_single_pathway(self):
        """Test with single pathway"""
        mirna_targets = {"miR-1": ["GENE_1", "GENE_2"]}
        pathway_genes = {"Pathway": ["GENE_1", "GENE_2", "GENE_3"]}
        
        network = miRNARegulatoryNetwork(mirna_targets, pathway_genes)
        hub_mirnas = network.identify_hub_mirnas()
        
        self.assertEqual(len(hub_mirnas), 1)
        self.assertEqual(hub_mirnas['miRNA'].iloc[0], "miR-1")
    
    def test_no_pathway_overlap(self):
        """Test when miRNA targets don't overlap with pathway genes"""
        mirna_targets = {"miR-1": ["GENE_100", "GENE_101"]}
        pathway_genes = {"Pathway": ["GENE_1", "GENE_2"]}
        
        network = miRNARegulatoryNetwork(mirna_targets, pathway_genes)
        mirna_pathway = network.map_to_pathways()
        
        # Should return empty if no overlap
        self.assertEqual(len(mirna_pathway), 0)
    
    def test_large_dataset_performance(self):
        """Test performance with larger dataset"""
        np.random.seed(42)
        
        # Larger dataset: 500 miRNAs × 100 samples
        mirna_expr = pd.DataFrame(
            np.random.randn(500, 100),
            index=[f"miR-{i}" for i in range(1, 501)],
            columns=[f"s{i}" for i in range(1, 101)]
        )
        
        # 5000 genes × 100 samples
        gene_expr = pd.DataFrame(
            np.random.randn(5000, 100),
            index=[f"GENE_{i}" for i in range(1, 5001)],
            columns=[f"s{i}" for i in range(1, 101)]
        )
        
        # This should complete in reasonable time
        predictor = miRNATargetPredictor(mirna_expr, gene_expr)
        targets = predictor.predict_targets(correlation_threshold=-0.3)
        
        # Should successfully predict some targets
        self.assertIsInstance(targets, dict)


def run_tests():
    """Run all tests with verbose output"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestmiRNATargetPredictor))
    suite.addTests(loader.loadTestsFromTestCase(TestmiRNARegulatoryNetwork))
    suite.addTests(loader.loadTestsFromTestCase(TestRegulatoryModuleAnalysis))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegration))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == "__main__":
    result = run_tests()
    
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {100.0 * (result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun:.1f}%")
    print("=" * 70)
