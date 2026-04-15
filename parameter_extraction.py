"""
Phase 3: Parameter Extraction Module

Extracts SIS model parameters (β, γ, initial infection, network structure)
from Phase 2 disease modules and Phase 1 pathway data.

Components:
1. Transmission rate (β) extraction from network connectivity and expression variance
2. Recovery rate (γ) extraction from pathway redundancy
3. Initial infection state from differential expression
4. Network structure from co-expression and PPI data

This module bridges molecular-scale data to population-scale epidemic dynamics.
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from typing import Dict, List, Tuple, Optional
import logging
import warnings

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ParameterExtractor:
    """
    Extract SIS model parameters from disease modules and networks
    
    Parameters are estimated from:
    - Network topology (for β)
    - Pathway redundancy (for γ)
    - Differential expression (for I₀)
    - Co-expression/PPI data (for adjacency matrix)
    """
    
    def __init__(self, disease_modules: Dict[str, List[str]],
                 expression_data: pd.DataFrame,
                 clinical_data: Optional[pd.DataFrame] = None,
                 adjacency_matrices: Optional[Dict[str, np.ndarray]] = None):
        """
        Initialize parameter extractor
        
        Parameters:
        -----------
        disease_modules : dict
            Module name → list of gene names
        expression_data : pd.DataFrame
            Gene expression matrix (genes × samples)
        clinical_data : pd.DataFrame or None
            Clinical traits (samples × traits) for disease vs control
        adjacency_matrices : dict or None
            Pre-computed adjacency matrices per module {module_name: adjacency_matrix}
        """
        self.disease_modules = disease_modules
        self.expression_data = expression_data
        self.clinical_data = clinical_data
        self.adjacency_matrices = adjacency_matrices or {}
        self.parameters = {}  # Store extracted parameters
        
        logger.info(f"Initialized ParameterExtractor:")
        logger.info(f"  Modules: {len(disease_modules)}")
        logger.info(f"  Expression: {expression_data.shape[0]} genes × {expression_data.shape[1]} samples")
        if clinical_data is not None:
            logger.info(f"  Clinical: {clinical_data.shape}")
    
    def extract_transmission_rate(self, module_name: str,
                                 connectivity_weight: float = 0.6,
                                 variance_weight: float = 0.4) -> float:
        """
        Extract transmission rate (β) from module structure
        
        Theory: Higher connected modules with more variable expression
        transmit dysregulation more effectively
        
        Formula:
        --------
        β = min(1, α * connectivity * expression_variability)
        
        where:
          connectivity = avg_degree / max_degree
          expression_variability = std(expr) / mean(expr)  [coefficient of variation]
          α ∈ [0.3, 0.7] (calibration, default 0.5)
        
        Parameters:
        -----------
        module_name : str
            Name of disease module
        connectivity_weight : float
            Weight for connectivity component (default 0.6)
        variance_weight : float
            Weight for variance component (default 0.4)
        
        Returns:
        --------
        float : β ∈ [0, 1]
        """
        logger.info(f"Extracting β for module '{module_name}'...")
        
        if module_name not in self.disease_modules:
            logger.warning(f"Module {module_name} not found. Using default β=0.3")
            return 0.3
        
        genes = self.disease_modules[module_name]
        
        # Get expressions for genes in module
        module_genes = [g for g in genes if g in self.expression_data.index]
        if len(module_genes) == 0:
            logger.warning(f"No genes found in expression data for {module_name}")
            return 0.3
        
        expr_subset = self.expression_data.loc[module_genes]
        
        # Compute connectivity
        if module_name in self.adjacency_matrices:
            A = self.adjacency_matrices[module_name]
            degrees = A.sum(axis=1)
            avg_degree = degrees.mean()
            max_degree = degrees.max()
            if max_degree > 0:
                connectivity = avg_degree / max_degree
            else:
                connectivity = 0.5
        else:
            # Default: assume modular network with moderate connectivity
            connectivity = 0.5
        
        # Compute expression variability (coefficient of variation)
        mean_expr = expr_subset.mean(axis=1).mean()
        std_expr = expr_subset.values.flatten().std()
        
        if mean_expr > 0:
            cv = std_expr / mean_expr
        else:
            cv = 0.5
        
        # Normalize CV to [0, 1]
        cv_normalized = min(1.0, cv / 2.0)  # Assume max reasonable CV is ~2
        
        # Combine components
        combined = connectivity_weight * connectivity + variance_weight * cv_normalized
        alpha = 0.5  # Calibration parameter
        beta = min(1.0, alpha * combined)
        
        logger.info(f"  β calculation:")
        logger.info(f"    Connectivity: {connectivity:.3f}")
        logger.info(f"    Expression variability (CV): {cv:.3f}")
        logger.info(f"    Combined: {combined:.3f}")
        logger.info(f"    β = {beta:.3f}")
        
        return beta
    
    def extract_recovery_rate(self, module_name: str,
                             redundancy_factor: Optional[float] = None) -> float:
        """
        Extract recovery rate (γ) from module robustness
        
        Theory: Modules with high pathway redundancy recover slowly
        (disease is more persistent)
        
        Formula:
        --------
        γ = 1 / (1 + redundancy_factor)
        
        where:
          redundancy_factor = n_alternative_paths / n_genes
          high redundancy → low γ (difficult to recover)
          low redundancy → high γ (easy to recover)
        
        Parameters:
        -----------
        module_name : str
            Name of disease module
        redundancy_factor : float or None
            Pre-computed redundancy (if None, estimated from network)
        
        Returns:
        --------
        float : γ ∈ [0, 1]
        """
        logger.info(f"Extracting γ for module '{module_name}'...")
        
        if module_name not in self.disease_modules:
            logger.warning(f"Module {module_name} not found. Using default γ=0.3")
            return 0.3
        
        if redundancy_factor is None:
            # Estimate from network structure
            if module_name in self.adjacency_matrices:
                A = self.adjacency_matrices[module_name]
                
                # Count alternative paths: sum of 2-hop connections
                n_genes = A.shape[0]
                two_hop = np.dot(A, A)  # A^2 matrix
                n_paths = two_hop.sum() / 2  # Divide by 2 for symmetry
                
                redundancy_factor = n_paths / max(1, n_genes)
            else:
                # Default: moderate redundancy
                redundancy_factor = 2.0
        
        # Ensure redundancy_factor is positive
        redundancy_factor = max(0.1, redundancy_factor)
        
        # Compute γ
        gamma = 1.0 / (1.0 + redundancy_factor)
        gamma = min(1.0, max(0.01, gamma))  # Clamp to [0.01, 1]
        
        logger.info(f"  γ calculation:")
        logger.info(f"    Redundancy factor: {redundancy_factor:.3f}")
        logger.info(f"    γ = {gamma:.3f}")
        
        return gamma
    
    def extract_initial_infection(self, module_name: str,
                                 disease_group: Optional[str] = None,
                                 fc_threshold: float = 0.5,
                                 pval_threshold: float = 0.05) -> np.ndarray:
        """
        Extract initial infection state from differential expression
        
        Theory: Dysregulated genes start "infected"; normal genes start susceptible
        
        Formula:
        --------
        For each gene g:
          log2FC[g] = mean(expr_disease[g]) - mean(expr_normal[g])
          p-value[g] = t-test between disease and normal groups
          
          I₀[g] = {
            min(1, max(0, |log2FC[g]| / fc_threshold))   if p-value < pval_threshold
            0                                             otherwise
          }
        
        Parameters:
        -----------
        module_name : str
            Name of disease module
        disease_group : str or None
            Clinical variable indicating disease state (e.g., "Stage")
        fc_threshold : float
            Fold-change threshold for scaling (default 0.5 log2 units)
        pval_threshold : float
            P-value threshold for significance (default 0.05)
        
        Returns:
        --------
        np.ndarray : Initial infection vector I₀ ∈ [0,1]^n
        """
        logger.info(f"Extracting I₀ for module '{module_name}'...")
        
        if module_name not in self.disease_modules:
            logger.warning(f"Module {module_name} not found. Using default I₀")
            return np.array([])
        
        genes = self.disease_modules[module_name]
        module_genes = [g for g in genes if g in self.expression_data.index]
        
        if len(module_genes) == 0:
            logger.warning(f"No genes in expression data for {module_name}")
            return np.array([])
        
        # Get expression for module genes
        expr_subset = self.expression_data.loc[module_genes]
        
        I0 = np.zeros(len(module_genes))
        
        # If disease group specified, compute differential expression
        if disease_group is not None and self.clinical_data is not None:
            if disease_group not in self.clinical_data.columns:
                logger.warning(f"Disease group '{disease_group}' not in clinical data")
            else:
                clinical_vals = self.clinical_data[disease_group]
                
                # Get disease vs non-disease samples
                disease_samples = clinical_vals[clinical_vals != 0].index
                healthy_samples = clinical_vals[clinical_vals == 0].index
                
                if len(disease_samples) > 0 and len(healthy_samples) > 0:
                    disease_expr = expr_subset[disease_samples]
                    healthy_expr = expr_subset[healthy_samples]
                    
                    for i, gene in enumerate(module_genes):
                        if gene in disease_expr.index and gene in healthy_expr.index:
                            disease_vals = disease_expr.loc[gene].values
                            healthy_vals = healthy_expr.loc[gene].values
                            
                            # T-test
                            if len(disease_vals) > 1 and len(healthy_vals) > 1:
                                _, pval = ttest_ind(disease_vals, healthy_vals)
                                
                                if pval < pval_threshold:
                                    # Compute log2 fold change
                                    mean_disease = disease_vals.mean()
                                    mean_healthy = healthy_vals.mean()
                                    
                                    if mean_healthy > 0:
                                        log2fc = np.log2(mean_disease / mean_healthy)
                                    else:
                                        log2fc = np.log2(mean_disease + 1)
                                    
                                    # Scale to [0, 1]
                                    I0[i] = min(1.0, max(0.0, abs(log2fc) / fc_threshold))
        else:
            # No disease group: use expression variance as proxy for dysregulation
            mean_expr = expr_subset.mean(axis=1)
            std_expr = expr_subset.std(axis=1)
            
            # Coefficient of variation
            cv = std_expr / (mean_expr + 1e-6)
            cv_norm = cv / cv.max() if cv.max() > 0 else cv
            
            I0 = np.array(cv_norm)
        
        logger.info(f"  I₀ extraction:")
        logger.info(f"    Mean I₀: {I0.mean():.3f}")
        logger.info(f"    Std I₀: {I0.std():.3f}")
        logger.info(f"    High infection (I₀ > 0.5): {(I0 > 0.5).sum()} genes")
        
        return I0
    
    def extract_network_structure(self, module_name: str,
                                 correlation_threshold: float = 0.3) -> np.ndarray:
        """
        Build adjacency matrix from co-expression or PPI data
        
        If pre-computed adjacency is available, returns it.
        Otherwise, computes from expression correlation.
        
        Parameters:
        -----------
        module_name : str
            Name of disease module
        correlation_threshold : float
            Minimum correlation for edge (default 0.3)
        
        Returns:
        --------
        np.ndarray : Weighted adjacency matrix (n × n)
        """
        logger.info(f"Extracting network structure for module '{module_name}'...")
        
        if module_name not in self.disease_modules:
            logger.warning(f"Module {module_name} not found.")
            return np.array([])
        
        # Check if pre-computed
        if module_name in self.adjacency_matrices:
            logger.info(f"  Using pre-computed adjacency matrix")
            return self.adjacency_matrices[module_name]
        
        genes = self.disease_modules[module_name]
        module_genes = [g for g in genes if g in self.expression_data.index]
        
        if len(module_genes) < 2:
            logger.warning(f"Module has < 2 genes. Using empty adjacency.")
            return np.zeros((1, 1))
        
        # Compute correlation matrix
        expr_subset = self.expression_data.loc[module_genes]
        corr_matrix = expr_subset.T.corr(method='pearson').values
        
        # Threshold and weight
        n = len(module_genes)
        adjacency = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                if abs(corr_matrix[i, j]) >= correlation_threshold:
                    # Weight by correlation strength
                    adjacency[i, j] = abs(corr_matrix[i, j])
                    adjacency[j, i] = abs(corr_matrix[i, j])
        
        logger.info(f"  Network structure:")
        logger.info(f"    Genes: {len(module_genes)}")
        logger.info(f"    Edges (|r| > {correlation_threshold}): {(adjacency > 0).sum() // 2}")
        logger.info(f"    Mean correlation: {corr_matrix[np.triu_indices_from(corr_matrix, k=1)].mean():.3f}")
        
        return adjacency
    
    def extract_all_parameters(self, module_name: str,
                              disease_group: Optional[str] = None) -> Dict[str, any]:
        """
        Extract all SIS parameters for a module
        
        Parameters:
        -----------
        module_name : str
            Name of disease module
        disease_group : str or None
            Clinical variable for disease state
        
        Returns:
        --------
        dict : {
            'module': str,
            'beta': float,
            'gamma': float,
            'initial_infection': np.ndarray,
            'adjacency': np.ndarray,
            'n_genes': int,
            'n_edges': int
        }
        """
        logger.info(f"\nExtracting all parameters for '{module_name}'...")
        
        beta = self.extract_transmission_rate(module_name)
        gamma = self.extract_recovery_rate(module_name)
        I0 = self.extract_initial_infection(module_name, disease_group)
        A = self.extract_network_structure(module_name)
        
        params = {
            'module': module_name,
            'beta': beta,
            'gamma': gamma,
            'initial_infection': I0,
            'adjacency': A,
            'n_genes': len(self.disease_modules.get(module_name, [])),
            'n_edges': (A > 0).sum() // 2 if A.size > 0 else 0
        }
        
        self.parameters[module_name] = params
        
        logger.info(f"✓ Parameters extracted for '{module_name}'")
        logger.info(f"  β={beta:.3f}, γ={gamma:.3f}, I₀_mean={I0.mean():.3f}")
        
        return params
    
    def get_parameter_summary(self) -> pd.DataFrame:
        """
        Get summary of all extracted parameters
        
        Returns:
        --------
        pd.DataFrame : Parameter summary table
        """
        summary_data = []
        
        for module_name, params in self.parameters.items():
            summary_data.append({
                'module': module_name,
                'n_genes': params['n_genes'],
                'n_edges': params['n_edges'],
                'beta': params['beta'],
                'gamma': params['gamma'],
                'I0_mean': params['initial_infection'].mean() if len(params['initial_infection']) > 0 else 0,
                'I0_max': params['initial_infection'].max() if len(params['initial_infection']) > 0 else 0
            })
        
        summary_df = pd.DataFrame(summary_data)
        logger.info(f"\nParameter Summary:")
        logger.info(f"\n{summary_df.to_string()}")
        
        return summary_df


if __name__ == "__main__":
    print("Parameter Extraction Module for Phase 3: SIS-as-Biomarker-Discovery")
    print("=" * 70)
    print("\nUsage:")
    print("  from parameter_extraction import ParameterExtractor")
    print("  extractor = ParameterExtractor(disease_modules, expression_data)")
    print("  params = extractor.extract_all_parameters(module_name)")
    print("  summary = extractor.get_parameter_summary()")
