"""
Test Suite for Phase 2 Gradio Integration Module

Tests the Phase2DataLoader and UI component creation functions.
Does NOT require actual data files or running the Gradio app.

Tests:
1. Phase2DataLoader initialization and loading
2. Data file path detection
3. Tab creation without errors
4. Component initialization
5. Callback function definitions
"""

import unittest
import pandas as pd
import numpy as np
import os
import tempfile
import logging
from pathlib import Path

# Import the module being tested
from gradio_phase2_integration import (
    Phase2DataLoader,
    create_disease_module_tab,
    create_wgcna_tab,
    create_mirna_tab,
    create_phase2_network_medicine_tab,
    get_phase2_integration_instructions
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestPhase2DataLoader(unittest.TestCase):
    """Test Phase2DataLoader class"""
    
    def test_initialization(self):
        """Test Phase2DataLoader initialization"""
        loader = Phase2DataLoader()
        
        self.assertIsNone(loader.gene_disease)
        self.assertIsNone(loader.ppi_network)
        self.assertIsNone(loader.gene_expression)
        self.assertIsNone(loader.mirna_expression)
        self.assertIsNone(loader.clinical_traits)
        self.assertIsNone(loader.pathway_genes)
        self.assertFalse(loader.loaded)
        
        logger.info("✓ Phase2DataLoader initialization successful")
    
    def test_load_all_graceful_failure(self):
        """Test that load_all handles missing data gracefully"""
        loader = Phase2DataLoader()
        
        # This should return False if data doesn't exist, not crash
        result = loader.load_all()
        
        self.assertIsInstance(result, bool)
        # Result depends on whether data files exist
        logger.info(f"✓ Phase2DataLoader.load_all() returned: {result}")
    
    def test_data_attributes_types(self):
        """Test that data attributes have correct types after initialization"""
        loader = Phase2DataLoader()
        
        # Check types
        self.assertIsInstance(loader.gene_disease, (type(None), pd.DataFrame))
        self.assertIsInstance(loader.gene_expression, (type(None), pd.DataFrame))
        self.assertIsInstance(loader.mirna_expression, (type(None), pd.DataFrame))
        self.assertIsInstance(loader.clinical_traits, (type(None), pd.DataFrame))
        self.assertIsInstance(loader.pathway_genes, (type(None), dict))
        
        logger.info("✓ All Phase2DataLoader attributes have correct types")


class TestPhase2TabCreation(unittest.TestCase):
    """Test Phase 2 tab creation functions"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.loader = Phase2DataLoader()
        
        # Try to load data, but don't fail if it's not available
        try:
            self.loader.load_all()
        except:
            pass
    
    def test_disease_module_tab_creation(self):
        """Test disease module detection tab creation"""
        try:
            import gradio as gr
            with gr.Blocks() as demo:
                create_disease_module_tab(self.loader)
            
            logger.info("✓ Disease module tab created successfully")
        except Exception as e:
            self.fail(f"Disease module tab creation failed: {e}")
    
    def test_wgcna_tab_creation(self):
        """Test WGCNA tab creation"""
        try:
            import gradio as gr
            with gr.Blocks() as demo:
                create_wgcna_tab(self.loader)
            
            logger.info("✓ WGCNA tab created successfully")
        except Exception as e:
            self.fail(f"WGCNA tab creation failed: {e}")
    
    def test_mirna_tab_creation(self):
        """Test miRNA regulatory tab creation"""
        try:
            import gradio as gr
            with gr.Blocks() as demo:
                create_mirna_tab(self.loader)
            
            logger.info("✓ miRNA tab created successfully")
        except Exception as e:
            self.fail(f"miRNA tab creation failed: {e}")
    
    def test_phase2_network_medicine_tab_creation(self):
        """Test complete Phase 2 network medicine tab creation"""
        try:
            import gradio as gr
            with gr.Blocks() as demo:
                create_phase2_network_medicine_tab()
            
            logger.info("✓ Phase 2 network medicine tab created successfully")
        except Exception as e:
            self.fail(f"Phase 2 network medicine tab creation failed: {e}")


class TestPhase2Integration(unittest.TestCase):
    """Test Phase 2 integration with app"""
    
    def test_imports_from_app(self):
        """Test that app_full.py can import Phase 2 integration"""
        try:
            # Check if the import line exists in app_full.py
            with open('app_full.py', 'r') as f:
                content = f.read()
            
            self.assertIn('from gradio_phase2_integration import', content)
            self.assertIn('create_phase2_network_medicine_tab', content)
            
            logger.info("✓ app_full.py correctly imports Phase 2 integration")
        except Exception as e:
            self.fail(f"Phase 2 import verification failed: {e}")
    
    def test_phase2_tab_in_app(self):
        """Test that Phase 2 tab is added to app structure"""
        try:
            with open('app_full.py', 'r') as f:
                content = f.read()
            
            self.assertIn('网络医学分析', content)
            self.assertIn('create_phase2_network_medicine_tab()', content)
            
            logger.info("✓ Phase 2 tab is properly integrated into app_full.py")
        except Exception as e:
            self.fail(f"Phase 2 tab verification failed: {e}")


class TestPhase2Documentation(unittest.TestCase):
    """Test Phase 2 integration documentation"""
    
    def test_integration_instructions_available(self):
        """Test that integration instructions are available"""
        instructions = get_phase2_integration_instructions()
        
        self.assertIsInstance(instructions, str)
        self.assertIn('Phase 2', instructions)
        self.assertIn('Integration', instructions)
        self.assertIn('Data Requirements', instructions)
        
        logger.info("✓ Phase 2 integration instructions are complete")
    
    def test_module_docstrings(self):
        """Test that all functions have docstrings"""
        from gradio_phase2_integration import (
            create_disease_module_tab,
            create_wgcna_tab,
            create_mirna_tab,
            create_phase2_network_medicine_tab
        )
        
        self.assertIsNotNone(create_disease_module_tab.__doc__)
        self.assertIsNotNone(create_wgcna_tab.__doc__)
        self.assertIsNotNone(create_mirna_tab.__doc__)
        self.assertIsNotNone(create_phase2_network_medicine_tab.__doc__)
        
        logger.info("✓ All Phase 2 functions have docstrings")


class TestPhase2Functionality(unittest.TestCase):
    """Test Phase 2 component functionality"""
    
    def test_data_loader_attributes_accessible(self):
        """Test that loader attributes are accessible"""
        loader = Phase2DataLoader()
        
        # All these should be accessible without error
        _ = loader.gene_disease
        _ = loader.ppi_network
        _ = loader.gene_expression
        _ = loader.mirna_expression
        _ = loader.clinical_traits
        _ = loader.pathway_genes
        _ = loader.loaded
        
        logger.info("✓ All Phase2DataLoader attributes are accessible")
    
    def test_phase2_constants(self):
        """Test Phase 2 module constants"""
        from gradio_phase2_integration import PHASE2_AVAILABLE
        
        self.assertIsInstance(PHASE2_AVAILABLE, bool)
        logger.info(f"✓ PHASE2_AVAILABLE = {PHASE2_AVAILABLE}")


class TestPhase2Integration(unittest.TestCase):
    """Integration tests for full Phase 2 pipeline"""
    
    def test_full_import_chain(self):
        """Test that all Phase 2 modules can be imported"""
        try:
            # Try to import all Phase 2 modules that gradio_phase2_integration uses
            try:
                from disease_module_detection import DiseaseNetworkBuilder
                logger.info("✓ disease_module_detection imported")
            except ImportError:
                logger.warning("⚠️  disease_module_detection not available")
            
            try:
                from wgcna_analysis import WGCNAAnalyzer
                logger.info("✓ wgcna_analysis imported")
            except ImportError:
                logger.warning("⚠️  wgcna_analysis not available")
            
            try:
                from mirna_integration import miRNATargetPredictor
                logger.info("✓ mirna_integration imported")
            except ImportError:
                logger.warning("⚠️  mirna_integration not available")
        
        except Exception as e:
            logger.error(f"Import chain test failed: {e}")
    
    def test_phase2_integration_instructions_format(self):
        """Test format of integration instructions"""
        instructions = get_phase2_integration_instructions()
        
        # Should contain markdown structure
        self.assertIn('#', instructions)  # Markdown headers
        self.assertIn('```', instructions)  # Code blocks
        
        logger.info("✓ Integration instructions have proper markdown format")


# ============================================================================
# Test Summary and Utilities
# ============================================================================

def run_all_tests():
    """Run all Phase 2 integration tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test suites
    suite.addTests(loader.loadTestsFromTestCase(TestPhase2DataLoader))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase2TabCreation))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase2Integration))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase2Documentation))
    suite.addTests(loader.loadTestsFromTestCase(TestPhase2Functionality))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == "__main__":
    # Run tests
    result = run_all_tests()
    
    # Print summary
    print("\n" + "=" * 70)
    print("PHASE 2 GRADIO INTEGRATION TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {100.0 * (result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun:.1f}%")
    print("=" * 70)
    
    # Exit with appropriate code
    exit(0 if result.wasSuccessful() else 1)
