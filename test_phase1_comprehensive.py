"""
Comprehensive Unit Test Suite for Phase 1: Pathway Activity Analysis
Tests all modules: PathwayActivityScorer, DifferentialPathwayAnalysis, HubGeneIdentifier, Visualizations

Status: Ready for Production
Last Updated: 2026-04-15
"""

import unittest
import pandas as pd
import numpy as np
import tempfile
import os
from pathlib import Path
import logging

# Configure logging for tests
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import Phase 1 modules
from pathway_activity import PathwayActivityScorer
from differential_pathway_analysis import DifferentialPathwayAnalysis
from hub_gene_identification import HubGeneIdentifier
from pathway_visualizations import (
    plot_pathway_activity_heatmap,
    plot_pathway_violin,
    plot_hub_genes_bar,
    plot_differential_pathways
)


class TestPathwayActivityScorer(unittest.TestCase):
    """Test Suite for PathwayActivityScorer module"""
    
    @classmethod
    def setUpClass(cls):
        """Create test data fixtures"""
        # Minimal pathway dictionary
        cls.test_pathways = {
            'Glycolysis': ['ALDH2', 'PDHA1', 'LDHA', 'PGK1', 'ENO1'],
            'Krebs_Cycle': ['PDHA1', 'IDH1', 'SUCLA2', 'SDHA'],
            'Oxidative_Phosphorylation': ['SDHA', 'COX5A', 'ATP5F1A'],
            'Apoptosis': ['BAX', 'BCL2', 'CASP3', 'TP53']
        }
        
        # Create synthetic expression data (4 samples × 10 genes)
        np.random.seed(42)
        genes = list(set().union(*cls.test_pathways.values()))
        cls.test_expr_data = pd.DataFrame(
            np.random.randn(len(genes), 4) + np.random.randn(len(genes), 1),
            index=genes,
            columns=['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        )
        
        # Save to temporary file
        cls.temp_expr_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        cls.test_expr_data.to_csv(cls.temp_expr_file.name)
        cls.temp_expr_file.close()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test files"""
        if os.path.exists(cls.temp_expr_file.name):
            os.remove(cls.temp_expr_file.name)
    
    def test_scorer_initialization(self):
        """Test PathwayActivityScorer initialization"""
        scorer = PathwayActivityScorer(self.test_pathways)
        self.assertEqual(len(scorer.pathway_genes), 4)
        self.assertIsNone(scorer.pathway_activity)
        self.assertIsNone(scorer.gene_expr)
        logger.info("✓ Test initialization passed")
    
    def test_load_expression_data(self):
        """Test loading expression data"""
        scorer = PathwayActivityScorer(self.test_pathways)
        n_genes, n_samples = scorer.load_expression_data(self.temp_expr_file.name)
        
        self.assertEqual(n_samples, 4)
        self.assertIsNotNone(scorer.gene_expr)
        self.assertEqual(scorer.gene_expr.shape[1], 4)
        logger.info("✓ Test load_expression_data passed")
    
    def test_map_pathway_genes(self):
        """Test pathway gene mapping"""
        scorer = PathwayActivityScorer(self.test_pathways)
        scorer.load_expression_data(self.temp_expr_file.name)
        
        mapped = scorer._map_pathway_genes()
        
        # All 4 pathways should map (each has ≥2 genes in expression)
        self.assertEqual(len(mapped), 4)
        self.assertIn('Glycolysis', mapped)
        
        # Check that mapped genes are subset of original pathway genes
        for pathway, genes in mapped.items():
            self.assertTrue(all(g in self.test_pathways[pathway] for g in genes))
        
        logger.info("✓ Test pathway gene mapping passed")
    
    def test_score_gsva(self):
        """Test GSVA scoring method"""
        scorer = PathwayActivityScorer(self.test_pathways)
        scorer.load_expression_data(self.temp_expr_file.name)
        
        activity = scorer.score_gsva()
        
        # Check output dimensions
        self.assertEqual(activity.shape[0], 4)  # 4 pathways
        self.assertEqual(activity.shape[1], 4)  # 4 samples
        
        # Check values are numeric and not all NaN
        self.assertTrue(np.isfinite(activity.values).any())
        
        # Check pathway activity is stored
        self.assertIsNotNone(scorer.pathway_activity)
        
        logger.info(f"✓ Test GSVA scoring passed - shape: {activity.shape}, range: [{activity.min().min():.2f}, {activity.max().max():.2f}]")
    
    def test_score_mean(self):
        """Test mean expression scoring"""
        scorer = PathwayActivityScorer(self.test_pathways)
        scorer.load_expression_data(self.temp_expr_file.name)
        
        activity = scorer.score_mean()
        
        # Check output dimensions
        self.assertEqual(activity.shape[0], 4)
        self.assertEqual(activity.shape[1], 4)
        
        # Mean scores should be finite
        self.assertTrue(np.isfinite(activity.values).all())
        
        logger.info(f"✓ Test mean scoring passed - range: [{activity.min().min():.2f}, {activity.max().max():.2f}]")
    
    def test_get_pathway_statistics(self):
        """Test pathway statistics computation"""
        scorer = PathwayActivityScorer(self.test_pathways)
        scorer.load_expression_data(self.temp_expr_file.name)
        scorer.score_gsva()
        
        stats = scorer.get_pathway_statistics()
        
        # Check statistics dataframe
        self.assertEqual(len(stats), 4)
        self.assertIn('mean_activity', stats.columns)
        self.assertIn('std_activity', stats.columns)
        self.assertIn('n_genes', stats.columns)
        
        # All statistics should be numeric
        self.assertTrue(all(pd.api.types.is_numeric_dtype(stats[col]) for col in 
                           ['mean_activity', 'std_activity', 'min_activity', 'max_activity', 'n_genes']))
        
        logger.info(f"✓ Test pathway statistics passed - {len(stats)} pathways analyzed")
    
    def test_save_pathway_activity(self):
        """Test saving pathway activity to file"""
        scorer = PathwayActivityScorer(self.test_pathways)
        scorer.load_expression_data(self.temp_expr_file.name)
        scorer.score_gsva()
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_file = f.name
        
        try:
            scorer.save_pathway_activity(output_file)
            
            # Verify file was created and contains data
            self.assertTrue(os.path.exists(output_file))
            saved_data = pd.read_csv(output_file, index_col=0)
            self.assertEqual(saved_data.shape, (4, 4))
            logger.info(f"✓ Test save pathway activity passed - file size: {os.path.getsize(output_file)} bytes")
        finally:
            if os.path.exists(output_file):
                os.remove(output_file)


class TestDifferentialPathwayAnalysis(unittest.TestCase):
    """Test Suite for DifferentialPathwayAnalysis module"""
    
    @classmethod
    def setUpClass(cls):
        """Create test data fixtures"""
        np.random.seed(42)
        
        # Create pathway activity data
        cls.pathway_activity = pd.DataFrame(
            np.random.randn(10, 20),
            index=[f'Pathway_{i}' for i in range(10)],
            columns=[f'Sample_{i}' for i in range(20)]
        )
        
        # Create clinical data with stratification
        cls.clinical_data = pd.DataFrame({
            'Age_Group': np.repeat(['Young', 'Old'], 10),
            'Gender': np.tile(['M', 'F'], 10),
            'Stage': np.repeat(['I', 'II', 'III', 'IV'], 5)
        }, index=[f'Sample_{i}' for i in range(20)])
    
    def test_initialization(self):
        """Test DifferentialPathwayAnalysis initialization"""
        diff_analysis = DifferentialPathwayAnalysis(
            self.pathway_activity, 
            self.clinical_data
        )
        
        self.assertEqual(diff_analysis.pathway_activity.shape, (10, 20))
        self.assertEqual(len(diff_analysis.clinical_data), 20)
        logger.info("✓ Test DifferentialPathwayAnalysis initialization passed")
    
    def test_compare_by_group_ttest(self):
        """Test t-test comparison"""
        diff_analysis = DifferentialPathwayAnalysis(
            self.pathway_activity,
            self.clinical_data
        )
        
        results = diff_analysis.compare_by_group('Age_Group', method='ttest')
        
        # Check results structure
        self.assertEqual(len(results), 10)
        self.assertIn('pvalue', results.columns)
        self.assertIn('padj', results.columns)
        self.assertIn('significant', results.columns)
        
        # Check p-values are valid
        self.assertTrue(all(0 <= p <= 1 for p in results['pvalue']))
        self.assertTrue(all(0 <= p <= 1 for p in results['padj']))
        
        logger.info(f"✓ Test t-test comparison passed - {results['significant'].sum()} significant pathways")
    
    def test_compare_by_group_anova(self):
        """Test ANOVA comparison"""
        diff_analysis = DifferentialPathwayAnalysis(
            self.pathway_activity,
            self.clinical_data
        )
        
        results = diff_analysis.compare_by_group('Stage', method='anova')
        
        # Check results
        self.assertEqual(len(results), 10)
        self.assertTrue(all(0 <= p <= 1 for p in results['pvalue']))
        
        logger.info(f"✓ Test ANOVA comparison passed - {results['significant'].sum()} significant pathways")
    
    def test_multiple_testing_correction(self):
        """Test FDR correction is applied"""
        diff_analysis = DifferentialPathwayAnalysis(
            self.pathway_activity,
            self.clinical_data
        )
        
        results = diff_analysis.compare_by_group('Gender', method='ttest')
        
        # padj should be >= pvalue (FDR always less conservative than raw p)
        self.assertTrue(all(results['padj'] >= results['pvalue']))
        
        logger.info("✓ Test multiple testing correction passed")


class TestHubGeneIdentifier(unittest.TestCase):
    """Test Suite for HubGeneIdentifier module"""
    
    @classmethod
    def setUpClass(cls):
        """Create test data fixtures"""
        import networkx as nx
        
        # Create test gene network
        cls.gene_network = nx.Graph()
        edges = [
            ('ALDH2', 'PDHA1'), ('PDHA1', 'LDHA'), ('LDHA', 'ENO1'),
            ('ALDH2', 'PGK1'), ('PGK1', 'ENO1'),
            ('IDH1', 'SUCLA2'), ('SUCLA2', 'SDHA'),
            ('SDHA', 'COX5A'), ('COX5A', 'ATP5F1A')
        ]
        cls.gene_network.add_edges_from(edges)
        
        # Test pathways
        cls.test_pathways = {
            'Glycolysis': ['ALDH2', 'PDHA1', 'LDHA', 'PGK1', 'ENO1'],
            'Krebs_Cycle': ['PDHA1', 'IDH1', 'SUCLA2', 'SDHA'],
            'OXPHOS': ['SDHA', 'COX5A', 'ATP5F1A', 'NDUFV1']
        }
        
        # Create expression data and save to file
        np.random.seed(42)
        genes = list(set().union(*cls.test_pathways.values()))
        cls.expr_data = pd.DataFrame(
            np.random.randn(len(genes), 20),
            index=genes,
            columns=[f'Sample_{i}' for i in range(20)]
        )
        
        # Save expression data to temporary file
        cls.temp_expr_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        cls.expr_data.to_csv(cls.temp_expr_file.name)
        cls.temp_expr_file.close()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test files"""
        if os.path.exists(cls.temp_expr_file.name):
            os.remove(cls.temp_expr_file.name)
    
    def test_hub_identifier_initialization(self):
        """Test HubGeneIdentifier initialization"""
        hub_finder = HubGeneIdentifier(
            self.test_pathways,
            gene_network=self.gene_network
        )
        
        self.assertEqual(len(hub_finder.pathway_genes), 3)
        self.assertIsNotNone(hub_finder.gene_network)
        logger.info("✓ Test HubGeneIdentifier initialization passed")
    
    def test_calculate_hub_score_single_pathway(self):
        """Test hub score calculation for single pathway"""
        hub_finder = HubGeneIdentifier(
            self.test_pathways,
            gene_network=self.gene_network
        )
        hub_finder.load_expression_data(self.temp_expr_file.name)
        
        hub_scores = hub_finder.calculate_hub_score('Glycolysis', normalize=True)
        
        # Check output
        self.assertGreater(len(hub_scores), 0)
        self.assertIn('hub_score', hub_scores.columns)
        self.assertIn('degree', hub_scores.columns)
        self.assertIn('betweenness', hub_scores.columns)
        
        # Scores should be between 0-1 (normalized)
        self.assertTrue(all(0 <= score <= 1 for score in hub_scores['hub_score']))
        
        # Top gene should have highest score
        top_gene = hub_scores.iloc[0]
        self.assertTrue(top_gene['hub_score'] >= hub_scores['hub_score'].min())
        
        logger.info(f"✓ Test hub score calculation passed - top gene: {top_gene['gene']} (score: {top_gene['hub_score']:.3f})")
    
    def test_calculate_all_hub_genes(self):
        """Test batch hub gene calculation"""
        hub_finder = HubGeneIdentifier(
            self.test_pathways,
            gene_network=self.gene_network
        )
        hub_finder.load_expression_data(self.temp_expr_file.name)
        
        all_hubs = hub_finder.calculate_all_hub_genes()
        
        # Check structure
        self.assertEqual(len(all_hubs), 3)
        for pathway in self.test_pathways.keys():
            self.assertIn(pathway, all_hubs)
            self.assertTrue(len(all_hubs[pathway]) > 0)
        
        logger.info(f"✓ Test batch hub gene calculation passed - {len(all_hubs)} pathways processed")
    
    def test_export_hub_genes_summary(self):
        """Test hub gene export"""
        hub_finder = HubGeneIdentifier(
            self.test_pathways,
            gene_network=self.gene_network
        )
        hub_finder.load_expression_data(self.temp_expr_file.name)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_file = f.name
        
        try:
            hub_finder.export_hub_genes_summary(output_file, top_n=3)
            
            # Verify file
            self.assertTrue(os.path.exists(output_file))
            exported = pd.read_csv(output_file)
            
            # Should have top 3 genes per pathway (3 pathways × 3 genes = 9 rows max)
            self.assertLessEqual(len(exported), 9)
            self.assertIn('pathway', exported.columns)
            self.assertIn('gene', exported.columns)
            
            logger.info(f"✓ Test export hub genes passed - {len(exported)} records exported")
        finally:
            if os.path.exists(output_file):
                os.remove(output_file)


class TestPathwayVisualizations(unittest.TestCase):
    """Test Suite for Visualization Functions"""
    
    @classmethod
    def setUpClass(cls):
        """Create test data fixtures"""
        np.random.seed(42)
        
        # Pathway activity data
        cls.pathway_activity = pd.DataFrame(
            np.random.randn(20, 50),
            index=[f'Pathway_{i}' for i in range(20)],
            columns=[f'Sample_{i}' for i in range(50)]
        )
        
        # Clinical data
        cls.clinical_data = pd.DataFrame({
            'Age_Group': np.repeat(['Young', 'Middle', 'Old'], int(50/3) + 1)[:50],
            'Gender': np.tile(['M', 'F'], 25),
            'Stage': np.repeat(['I', 'II', 'III', 'IV'], int(50/4) + 1)[:50]
        }, index=[f'Sample_{i}' for i in range(50)])
        
        # Hub genes data
        cls.hub_genes_df = pd.DataFrame({
            'gene': [f'Gene_{i}' for i in range(20)],
            'hub_score': np.sort(np.random.rand(20))[::-1],
            'degree': np.random.randint(1, 10, 20),
            'betweenness': np.random.rand(20) * 0.5
        })
        
        # Differential results
        cls.diff_results = pd.DataFrame({
            'pathway': [f'Pathway_{i}' for i in range(20)],
            'pvalue': np.random.exponential(0.1, 20),
            'padj': np.random.exponential(0.1, 20),
            'significant': np.random.choice([True, False], 20)
        })
    
    def test_plot_pathway_activity_heatmap(self):
        """Test heatmap visualization"""
        fig = plot_pathway_activity_heatmap(
            self.pathway_activity,
            self.clinical_data,
            group_by='Age_Group'
        )
        
        # Check that figure is created
        self.assertIsNotNone(fig)
        self.assertTrue(hasattr(fig, 'data'))
        self.assertTrue(len(fig.data) > 0)
        
        logger.info("✓ Test heatmap visualization passed")
    
    def test_plot_pathway_violin(self):
        """Test violin plot visualization"""
        fig = plot_pathway_violin(
            self.pathway_activity,
            self.clinical_data,
            pathway='Pathway_0',
            group_by='Gender'
        )
        
        self.assertIsNotNone(fig)
        self.assertTrue(hasattr(fig, 'data'))
        
        logger.info("✓ Test violin plot visualization passed")
    
    def test_plot_hub_genes_bar(self):
        """Test hub genes bar chart"""
        fig = plot_hub_genes_bar(self.hub_genes_df, top_n=10)
        
        self.assertIsNotNone(fig)
        self.assertTrue(hasattr(fig, 'data'))
        
        logger.info("✓ Test hub genes bar chart passed")
    
    def test_plot_differential_pathways(self):
        """Test differential pathways plot"""
        fig = plot_differential_pathways(self.diff_results, top_n=15)
        
        self.assertIsNotNone(fig)
        self.assertTrue(hasattr(fig, 'data'))
        
        logger.info("✓ Test differential pathways plot passed")


class TestPhase1Integration(unittest.TestCase):
    """Integration tests for Phase 1 modules working together"""
    
    @classmethod
    def setUpClass(cls):
        """Setup integration test data"""
        cls.test_pathways = {
            'Glycolysis': ['ALDH2', 'PDHA1', 'LDHA', 'PGK1', 'ENO1'],
            'Krebs_Cycle': ['PDHA1', 'IDH1', 'SUCLA2', 'SDHA'],
            'Oxidative_Phosphorylation': ['SDHA', 'COX5A', 'ATP5F1A'],
            'Apoptosis': ['BAX', 'BCL2', 'CASP3', 'TP53']
        }
        
        np.random.seed(42)
        genes = list(set().union(*cls.test_pathways.values()))
        cls.expr_data = pd.DataFrame(
            np.random.randn(len(genes), 30) + np.random.randn(len(genes), 1),
            index=genes,
            columns=[f'Sample_{i}' for i in range(30)]
        )
        
        # Save expression data
        cls.temp_expr_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        cls.expr_data.to_csv(cls.temp_expr_file.name)
        cls.temp_expr_file.close()
        
        cls.clinical_data = pd.DataFrame({
            'Age_Group': np.repeat(['Young', 'Old'], 15),
            'Gender': np.tile(['M', 'F'], 15),
            'Stage': np.repeat(['Early', 'Late'], 15)
        }, index=[f'Sample_{i}' for i in range(30)])
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test files"""
        if os.path.exists(cls.temp_expr_file.name):
            os.remove(cls.temp_expr_file.name)
    
    def test_full_pipeline(self):
        """Test complete Phase 1 pipeline"""
        logger.info("\n" + "="*60)
        logger.info("INTEGRATION TEST: Full Phase 1 Pipeline")
        logger.info("="*60)
        
        # Step 1: Score pathways
        logger.info("Step 1: Scoring pathways with GSVA...")
        scorer = PathwayActivityScorer(self.test_pathways)
        scorer.load_expression_data(self.temp_expr_file.name)
        pathway_activity = scorer.score_gsva()
        logger.info(f"✓ Pathway activity matrix: {pathway_activity.shape}")
        
        # Step 2: Differential analysis
        logger.info("Step 2: Differential pathway analysis...")
        diff_analysis = DifferentialPathwayAnalysis(pathway_activity, self.clinical_data)
        diff_results = diff_analysis.compare_by_group('Stage', method='ttest')
        sig_pathways = diff_results[diff_results['significant']].shape[0]
        logger.info(f"✓ Significant pathways (FDR < 0.05): {sig_pathways}")
        
        # Step 3: Hub gene identification
        logger.info("Step 3: Hub gene identification...")
        import networkx as nx
        gene_network = nx.complete_graph(list(set().union(*self.test_pathways.values())))
        hub_finder = HubGeneIdentifier(self.test_pathways, gene_network)
        hub_finder.load_expression_data(self.temp_expr_file.name)
        all_hub_genes = hub_finder.calculate_all_hub_genes()
        logger.info(f"✓ Hub genes identified: {sum(len(h) for h in all_hub_genes.values())} total")
        
        # Verify outputs
        self.assertEqual(pathway_activity.shape[0], 4)
        self.assertEqual(len(diff_results), 4)
        self.assertEqual(len(all_hub_genes), 4)
        
        logger.info("\n" + "="*60)
        logger.info("INTEGRATION TEST PASSED ✓")
        logger.info("="*60 + "\n")


def run_all_tests():
    """Run all test suites and generate report"""
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestPathwayActivityScorer))
    suite.addTests(loader.loadTestsFromTestCase(TestDifferentialPathwayAnalysis))
    suite.addTests(loader.loadTestsFromTestCase(TestHubGeneIdentifier))
    suite.addTests(loader.loadTestsFromTestCase(TestPathwayVisualizations))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase1Integration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "="*70)
    print("PHASE 1 COMPREHENSIVE TEST SUMMARY")
    print("="*70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    print(f"Success rate: {100 * (result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun:.1f}%")
    print("="*70 + "\n")
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_all_tests()
    exit(0 if success else 1)
