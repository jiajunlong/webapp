"""
Integration Tests for Phase 3: Cross-Scale SIS Biomarker Discovery

Tests end-to-end Phase 1 → Phase 2 → Phase 3 pipeline:
- Parameter extraction from disease modules
- SIS network propagation
- Biomarker validation
- Cross-scale consistency

Test Coverage: 8 integration tests
Status: Ready for execution
"""

import unittest
import numpy as np
import pandas as pd
import networkx as nx
from typing import Dict, List
import logging

# Import Phase 3 modules
from parameter_extraction import ParameterExtractor
from sis_network_propagation import SISNetworkPropagation
from biomarker_validation import BiomarkerValidator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestPhase3EndToEndPipeline(unittest.TestCase):
    """End-to-end tests for complete Phase 1-2-3 pipeline"""
    
    def setUp(self):
        """
        Create simulated Phase 1-2 data:
        - Gene expression (Phase 1)
        - Disease modules (Phase 2)
        - Clinical data
        """
        np.random.seed(42)
        
        # ==================== Phase 1 Simulation ====================
        # Gene expression: 100 genes × 120 samples (80 healthy, 40 disease)
        self.n_genes_total = 100
        self.n_genes_in_module = 25
        self.n_healthy = 80
        self.n_disease = 40
        self.n_samples = self.n_healthy + self.n_disease
        
        # Simulate baseline expression with HIGHER baseline noise (SD=3.0) for stronger dysregulation signal
        healthy_expr = np.random.normal(5, 3.0, (self.n_genes_total, self.n_healthy))
        
        # Disease expression: dysregulated disease module
        disease_expr = np.random.normal(5, 3.0, (self.n_genes_total, self.n_disease))
        
        # Make disease module genes dysregulated with strong signal
        disease_module_genes = list(range(self.n_genes_in_module))
        for i in disease_module_genes:
            disease_expr[i, :] += np.random.normal(4, 1.0, self.n_disease)
        
        expr_data = np.hstack([healthy_expr, disease_expr])
        
        self.expr_data = pd.DataFrame(
            expr_data,
            index=[f"Gene_{i}" for i in range(self.n_genes_total)],
            columns=[f"Sample_{i}" for i in range(self.n_samples)]
        )
        
        # ==================== Phase 2 Simulation ====================
        # Disease module: 25 genes with PPI network
        self.disease_module_name = "TestDisease_Module"
        self.disease_module_genes = [f"Gene_{i}" for i in range(self.n_genes_in_module)]
        
        # Create PPI network (realistic scale-free topology)
        self.ppi_network = nx.barabasi_albert_graph(self.n_genes_in_module, 3)
        self.ppi_network = nx.relabel_nodes(
            self.ppi_network,
            {i: self.disease_module_genes[i] for i in range(self.n_genes_in_module)}
        )
        
        # Store disease modules mapping
        self.disease_modules = {
            self.disease_module_name: self.disease_module_genes
        }
        
        # Create adjacency matrix from PPI network
        self.adjacency_matrices = {
            self.disease_module_name: nx.to_numpy_array(self.ppi_network)
        }
        
        # ==================== Clinical Data ====================
        self.clinical_data = pd.DataFrame({
            'disease_status': ['healthy'] * self.n_healthy + ['disease'] * self.n_disease,
            'survival_days': np.random.exponential(500, self.n_samples) + 100,
            'age': np.random.randint(30, 80, self.n_samples),
            'stage': np.random.choice(['I', 'II', 'III'], self.n_samples)
        }, index=[f"Sample_{i}" for i in range(self.n_samples)])
        
        # Correlate survival with disease module expression
        for gene in self.disease_module_genes[:5]:
            idx = int(gene.split('_')[1])
            self.clinical_data['survival_days'] += self.expr_data.loc[gene].values * 80
    
    def test_parameter_extraction_pipeline(self):
        """
        Test Phase 2 → Phase 3 bridge: Extract parameters from disease modules
        """
        logger.info("Testing parameter extraction from disease modules...")
        
        extractor = ParameterExtractor(
            disease_modules=self.disease_modules,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data,
            adjacency_matrices=self.adjacency_matrices
        )
        
        # Extract parameters for disease module (per-module extraction)
        module_params = extractor.extract_all_parameters(
            module_name=self.disease_module_name,
            disease_group='disease_status'
        )
        
        # Verify parameter extraction
        self.assertIn('beta', module_params)
        self.assertIn('gamma', module_params)
        self.assertIn('initial_infection', module_params)
        self.assertIn('adjacency', module_params)
        
        # Verify parameter ranges
        self.assertGreaterEqual(module_params['beta'], 0)
        self.assertLessEqual(module_params['beta'], 1)
        
        self.assertGreater(module_params['gamma'], 0)
        self.assertLessEqual(module_params['gamma'], 1)
        
        self.assertEqual(len(module_params['initial_infection']), 
                        len(self.disease_module_genes))
        
        logger.info(f"✓ Extracted parameters: β={module_params['beta']:.3f}, "
                   f"γ={module_params['gamma']:.3f}")
    
    def test_sis_propagation_pipeline(self):
        """
        Test Phase 3 core: Run SIS dynamics on extracted parameters
        """
        logger.info("Testing SIS network propagation...")
        
        # Step 1: Extract parameters
        extractor = ParameterExtractor(
            disease_modules=self.disease_modules,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data,
            adjacency_matrices=self.adjacency_matrices
        )
        
        module_params = extractor.extract_all_parameters(
            module_name=self.disease_module_name,
            disease_group='disease_status'
        )
        
        # Step 2: Run SIS dynamics with CORRECT constructor signature
        # Parameters: adjacency_matrix, gene_names, parameters (dict)
        propagation = SISNetworkPropagation(
            adjacency_matrix=module_params['adjacency'],
            gene_names=self.disease_module_genes,
            parameters={
                'beta': module_params['beta'],
                'gamma': module_params['gamma'],
                'initial_infection': module_params['initial_infection']
            }
        )
        
        propagation.run_dynamics(n_steps=500, n_runs=50)
        
        # Verify SIS output
        persistence_scores = propagation.persistence_scores
        self.assertIsNotNone(persistence_scores)
        self.assertEqual(len(persistence_scores), len(self.disease_module_genes))
        self.assertTrue((persistence_scores >= 0).all())
        self.assertTrue((persistence_scores <= 1).all())
        
        logger.info(f"✓ SIS dynamics completed, persistence scores computed")
    
    def test_biomarker_validation_pipeline(self):
        """
        Test Phase 3 validation: Validate SIS-predicted biomarkers
        """
        logger.info("Testing biomarker validation...")
        
        # Step 1: Get biomarkers from SIS
        extractor = ParameterExtractor(
            disease_modules=self.disease_modules,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data,
            adjacency_matrices=self.adjacency_matrices
        )
        
        module_params = extractor.extract_all_parameters(
            module_name=self.disease_module_name,
            disease_group='disease_status'
        )
        
        propagation = SISNetworkPropagation(
            adjacency_matrix=module_params['adjacency'],
            gene_names=self.disease_module_genes,
            parameters={
                'beta': module_params['beta'],
                'gamma': module_params['gamma'],
                'initial_infection': module_params['initial_infection']
            }
        )
        
        propagation.run_dynamics(n_steps=500, n_runs=50)
        
        # Get top biomarkers by persistence
        persistence_df = pd.DataFrame({
            'gene': self.disease_module_genes,
            'persistence': propagation.persistence_scores
        }).sort_values('persistence', ascending=False)
        
        biomarkers = persistence_df.head(10)['gene'].tolist()
        
        # Step 2: Validate biomarkers with CORRECT constructor signature
        # Parameters: predicted_biomarkers, biomarker_scores, expression_data, clinical_data
        biomarker_scores = persistence_df.set_index('gene')['persistence'].to_dict()
        scores_array = np.array([biomarker_scores[b] for b in biomarkers])
        
        validator = BiomarkerValidator(
            predicted_biomarkers=biomarkers,
            biomarker_scores=scores_array,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data
        )
        
        # Expression validation
        expr_val = validator.validate_expression_changes(
            disease_stage_col='disease_status'
        )
        
        self.assertEqual(len(expr_val), len(biomarkers))
        self.assertIn('cv', expr_val.columns)
        self.assertIn('is_dysregulated', expr_val.columns)
        
        # Clinical correlation
        clinical_val = validator.validate_clinical_correlation(
            outcome_variable='survival_days'
        )
        
        self.assertEqual(len(clinical_val), len(biomarkers))
        self.assertIn('pearson_r', clinical_val.columns)
        self.assertIn('is_significant', clinical_val.columns)
        
        logger.info(f"✓ Validated {len(biomarkers)} biomarkers")
        logger.info(f"  - Dysregulated: {expr_val['is_dysregulated'].sum()}/{len(biomarkers)}")
        logger.info(f"  - Significant: {clinical_val['is_significant'].sum()}/{len(biomarkers)}")
    
    def test_cross_scale_consistency(self):
        """
        Test consistency across scales: gene → module → disease
        
        Verifies that biomarkers identified from disease modules:
        1. Belong to the disease module
        2. Show dysregulation in disease samples
        3. Correlate with clinical outcomes
        4. Have realistic parameter values
        """
        logger.info("Testing cross-scale consistency...")
        
        # Extract parameters
        extractor = ParameterExtractor(
            disease_modules=self.disease_modules,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data,
            adjacency_matrices=self.adjacency_matrices
        )
        
        module_params = extractor.extract_all_parameters(
            module_name=self.disease_module_name,
            disease_group='disease_status'
        )
        
        # Run SIS dynamics
        propagation = SISNetworkPropagation(
            adjacency_matrix=module_params['adjacency'],
            gene_names=self.disease_module_genes,
            parameters={
                'beta': module_params['beta'],
                'gamma': module_params['gamma'],
                'initial_infection': module_params['initial_infection']
            }
        )
        
        propagation.run_dynamics(n_steps=500, n_runs=50)
        
        persistence_df = pd.DataFrame({
            'gene': self.disease_module_genes,
            'persistence': propagation.persistence_scores
        }).sort_values('persistence', ascending=False)
        
        biomarkers = persistence_df.head(10)['gene'].tolist()
        
        # Check 1: Biomarkers are in disease module
        for biomarker in biomarkers:
            self.assertIn(biomarker, self.disease_module_genes,
                         f"{biomarker} not in disease module")
        
        logger.info(f"✓ All {len(biomarkers)} biomarkers verified in disease module")
        
        # Check 2: Biomarkers show dysregulation
        biomarker_scores = persistence_df.set_index('gene')['persistence'].to_dict()
        scores_array = np.array([biomarker_scores[b] for b in biomarkers])
        
        validator = BiomarkerValidator(
            predicted_biomarkers=biomarkers,
            biomarker_scores=scores_array,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data
        )
        
        expr_val = validator.validate_expression_changes(
            disease_stage_col='disease_status'
        )
        
        dysregulated_pct = (expr_val['is_dysregulated'].sum() / len(biomarkers)) * 100
        self.assertGreater(dysregulated_pct, 40,
                          f"Only {dysregulated_pct}% of biomarkers show dysregulation")
        
        logger.info(f"✓ {dysregulated_pct:.1f}% of biomarkers dysregulated")
        
        # Check 3: Parameter values reasonable for disease
        self.assertGreater(module_params['beta'], 0.1,
                          "Transmission rate too low for active disease module")
        self.assertLess(module_params['gamma'], 0.8,
                       "Recovery rate too high for persistent disease")
        
        logger.info(f"✓ Parameter values consistent with disease dynamics")
    
    def test_multiple_disease_modules(self):
        """
        Test pipeline with multiple disease modules
        """
        logger.info("Testing pipeline with multiple disease modules...")
        
        # Create second disease module
        module2_genes = [f"Gene_{i}" for i in range(self.n_genes_in_module, 
                                                     2*self.n_genes_in_module)]
        
        # Dysregulate module 2
        for i in range(self.n_genes_in_module, 2*self.n_genes_in_module):
            self.expr_data.iloc[i, self.n_healthy:] += np.random.normal(3.5, 0.8, 
                                                                         self.n_disease)
        
        disease_modules_multi = {
            "Disease_Module_1": self.disease_module_genes,
            "Disease_Module_2": module2_genes
        }
        
        # Create adjacency matrices for both
        module2_network = nx.complete_graph(len(module2_genes))
        module2_network = nx.relabel_nodes(
            module2_network,
            {i: module2_genes[i] for i in range(len(module2_genes))}
        )
        
        adjacency_multi = {
            "Disease_Module_1": self.adjacency_matrices[self.disease_module_name],
            "Disease_Module_2": nx.to_numpy_array(module2_network)
        }
        
        # Extract parameters for both modules
        extractor = ParameterExtractor(
            disease_modules=disease_modules_multi,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data,
            adjacency_matrices=adjacency_multi
        )
        
        all_params = {}
        for module_name in disease_modules_multi.keys():
            module_params = extractor.extract_all_parameters(
                module_name=module_name,
                disease_group='disease_status'
            )
            all_params[module_name] = module_params
        
        self.assertEqual(len(all_params), 2)
        for module_name in disease_modules_multi.keys():
            self.assertIn(module_name, all_params)
        
        logger.info(f"✓ Successfully processed {len(all_params)} disease modules")


class TestPhase3ParameterSensitivity(unittest.TestCase):
    """Test sensitivity of biomarker discovery to parameter variations"""
    
    def setUp(self):
        """Create test network"""
        np.random.seed(42)
        
        self.n_genes = 20
        self.n_samples = 80
        
        # Gene expression
        expr_data = np.random.normal(5, 1, (self.n_genes, self.n_samples))
        
        self.expr_data = pd.DataFrame(
            expr_data,
            index=[f"Gene_{i}" for i in range(self.n_genes)],
            columns=[f"Sample_{i}" for i in range(self.n_samples)]
        )
        
        # Clinical data
        self.clinical_data = pd.DataFrame({
            'disease': ['healthy']*40 + ['disease']*40
        }, index=[f"Sample_{i}" for i in range(self.n_samples)])
        
        # Network
        self.network = nx.complete_graph(self.n_genes)
        self.network = nx.relabel_nodes(
            self.network,
            {i: f"Gene_{i}" for i in range(self.n_genes)}
        )
        
        # Disease module
        self.genes = [f"Gene_{i}" for i in range(self.n_genes)]
        
        # Create adjacency matrix
        self.adjacency = nx.to_numpy_array(self.network)
    
    def test_transmission_rate_effect(self):
        """Test effect of transmission rate on biomarker discovery"""
        logger.info("Testing transmission rate sensitivity...")
        
        initial_infection = np.random.uniform(0, 0.5, self.n_genes)
        
        results_by_beta = {}
        
        for beta in [0.1, 0.3, 0.5, 0.7]:
            propagation = SISNetworkPropagation(
                adjacency_matrix=self.adjacency,
                gene_names=self.genes,
                parameters={
                    'beta': beta,
                    'gamma': 0.3,
                    'initial_infection': initial_infection
                }
            )
            
            propagation.run_dynamics(n_steps=300, n_runs=30)
            persistence = propagation.persistence_scores
            # Count genes with high persistence
            high_persist = (persistence > 0.5).sum()
            results_by_beta[beta] = high_persist
        
        # Higher β should result in more persistent genes
        self.assertLessEqual(results_by_beta[0.1], results_by_beta[0.7],
                            "Higher transmission rate should increase persistence")
        
        logger.info(f"✓ High persistence genes by β: {results_by_beta}")
    
    def test_recovery_rate_effect(self):
        """Test effect of recovery rate on biomarker discovery"""
        logger.info("Testing recovery rate sensitivity...")
        
        initial_infection = np.random.uniform(0, 0.5, self.n_genes)
        
        results_by_gamma = {}
        
        for gamma in [0.1, 0.3, 0.5, 0.7]:
            propagation = SISNetworkPropagation(
                adjacency_matrix=self.adjacency,
                gene_names=self.genes,
                parameters={
                    'beta': 0.4,
                    'gamma': gamma,
                    'initial_infection': initial_infection
                }
            )
            
            propagation.run_dynamics(n_steps=300, n_runs=30)
            persistence = propagation.persistence_scores
            high_persist = (persistence > 0.5).sum()
            results_by_gamma[gamma] = high_persist
        
        # Higher γ (faster recovery) should result in fewer persistent genes
        self.assertGreaterEqual(results_by_gamma[0.1], results_by_gamma[0.7],
                               "Higher recovery rate should decrease persistence")
        
        logger.info(f"✓ High persistence genes by γ: {results_by_gamma}")


class TestPhase3ValidationConsistency(unittest.TestCase):
    """Test consistency of validation metrics across different biomarker sets"""
    
    def setUp(self):
        """Create data with known structure"""
        np.random.seed(42)
        
        n_genes = 50
        n_samples = 100
        
        # Expression with clear dysregulation for first 10 genes
        expr_data = np.random.normal(5, 1, (n_genes, n_samples))
        expr_data[:10, 50:] += np.random.normal(2.5, 0.5, (10, 50))
        
        self.expr_data = pd.DataFrame(
            expr_data,
            index=[f"Gene_{i}" for i in range(n_genes)],
            columns=[f"Sample_{i}" for i in range(n_samples)]
        )
        
        self.clinical_data = pd.DataFrame({
            'disease': ['healthy']*50 + ['disease']*50,
            'outcome': np.random.exponential(400, n_samples)
        }, index=[f"Sample_{i}" for i in range(n_samples)])
        
        # Correlate outcome with dysregulated genes
        for i in range(10):
            self.clinical_data['outcome'] += self.expr_data.iloc[i, :].values * 50
    
    def test_biomarker_set_comparison(self):
        """Test validation consistency for different biomarker sets"""
        logger.info("Testing biomarker validation consistency...")
        
        # True biomarkers (first 10 dysregulated genes)
        true_biomarkers = [f"Gene_{i}" for i in range(10)]
        
        # Random biomarkers for comparison
        random_biomarkers = [f"Gene_{i}" for i in range(20, 30)]
        
        # Validate true biomarkers
        true_scores = np.linspace(0.9, 0.7, len(true_biomarkers))
        validator_true = BiomarkerValidator(
            predicted_biomarkers=true_biomarkers,
            biomarker_scores=true_scores,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data
        )
        
        expr_val_true = validator_true.validate_expression_changes(
            disease_stage_col='disease'
        )
        clinical_val_true = validator_true.validate_clinical_correlation(
            outcome_variable='outcome'
        )
        
        # Validate random biomarkers
        random_scores = np.linspace(0.5, 0.3, len(random_biomarkers))
        validator_random = BiomarkerValidator(
            predicted_biomarkers=random_biomarkers,
            biomarker_scores=random_scores,
            expression_data=self.expr_data,
            clinical_data=self.clinical_data
        )
        
        expr_val_random = validator_random.validate_expression_changes(
            disease_stage_col='disease'
        )
        clinical_val_random = validator_random.validate_clinical_correlation(
            outcome_variable='outcome'
        )
        
        # True biomarkers should have more dysregulation
        true_dysreg = expr_val_true['is_dysregulated'].sum()
        random_dysreg = expr_val_random['is_dysregulated'].sum()
        
        self.assertGreaterEqual(true_dysreg, random_dysreg,
                               "True biomarkers should show more dysregulation")
        
        logger.info(f"✓ True biomarkers: {true_dysreg}/10 dysregulated")
        logger.info(f"✓ Random biomarkers: {random_dysreg}/10 dysregulated")


if __name__ == '__main__':
    # Run tests with verbose output
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestPhase3EndToEndPipeline))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase3ParameterSensitivity))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase3ValidationConsistency))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "="*70)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print("="*70)
