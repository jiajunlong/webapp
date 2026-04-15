"""
Phase 3: SIS Network Propagation Engine

Runs stochastic SIS (Susceptible-Infected-Susceptible) dynamics on disease networks
to identify persistent biomarkers.

Key Concept:
- Nodes: Genes in a disease module network
- "Infection": Dysregulation (abnormal expression/activity)
- Transmission: Dysregulation spreads through network via regulatory/co-expression edges
- Recovery: Genes normalize their expression
- Biomarkers: Genes that remain dysregulated (high persistence)

Algorithm:
1. Initialize with I₀ (initial dysregulation from expression data)
2. For each timestep: transmission step (infected nodes infect neighbors with probability β)
3. For each timestep: recovery step (infected nodes recover with probability γ)
4. Track which genes remain "infected" over time
5. Genes with high persistence are robust biomarkers
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import Dict, List, Tuple, Optional
import logging
import warnings

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SISNetworkPropagation:
    """
    Run stochastic SIS dynamics on disease networks for biomarker discovery
    """
    
    def __init__(self, adjacency_matrix: np.ndarray,
                 gene_names: List[str],
                 parameters: Dict[str, float]):
        """
        Initialize SIS propagation engine
        
        Parameters:
        -----------
        adjacency_matrix : np.ndarray
            Weighted adjacency matrix (n × n)
        gene_names : List[str]
            Names of genes (nodes) in network
        parameters : dict
            {
                'beta': float,  # Transmission rate
                'gamma': float,  # Recovery rate
                'initial_infection': np.ndarray  # I₀ vector
            }
        """
        self.A = adjacency_matrix
        self.gene_names = gene_names
        self.beta = parameters.get('beta', 0.3)
        self.gamma = parameters.get('gamma', 0.3)
        self.I0 = parameters.get('initial_infection', None)
        
        self.n_nodes = len(gene_names)
        
        # Initialize I₀ if not provided
        if self.I0 is None:
            self.I0 = np.zeros(self.n_nodes)
            self.I0[0] = 1.0  # Start with first node infected
        
        # Validate dimensions
        if len(self.I0) != self.n_nodes:
            logger.warning(f"I₀ length {len(self.I0)} != n_nodes {self.n_nodes}")
            self.I0 = np.resize(self.I0, self.n_nodes)
        
        self.infection_history = None
        self.persistence_scores = None
        self.biomarkers = None
        
        logger.info(f"Initialized SISNetworkPropagation:")
        logger.info(f"  Network: {self.n_nodes} nodes, {(self.A > 0).sum() // 2} edges")
        logger.info(f"  Parameters: β={self.beta:.3f}, γ={self.gamma:.3f}")
        logger.info(f"  Initial infection: {self.I0.mean():.3f} (mean)")
    
    def run_dynamics(self, n_steps: int = 1000, n_runs: int = 100,
                    random_seed: Optional[int] = None) -> Dict:
        """
        Run stochastic SIS dynamics for multiple runs
        
        Algorithm:
        ----------
        For each run:
          I = I₀  # Initialize infection vector
          For each timestep:
            1. Transmission: For each infected node i, infect neighbors j with prob β*A[i,j]
            2. Recovery: For each infected node i, recover with prob γ
            3. Record infection state I
        
        Parameters:
        -----------
        n_steps : int
            Number of simulation steps (default 1000)
        n_runs : int
            Number of stochastic realizations (default 100)
        random_seed : int or None
            Random seed for reproducibility
        
        Returns:
        --------
        dict : {
            'infection_history': np.ndarray (n_runs × n_steps × n_nodes),
            'persistence_scores': np.ndarray (n_nodes,),
            'biomarkers': list of top genes,
            'statistics': dict with summary stats
        }
        """
        logger.info(f"Running SIS dynamics: {n_runs} runs × {n_steps} steps...")
        
        if random_seed is not None:
            np.random.seed(random_seed)
        
        # Initialize storage
        infection_history = np.zeros((n_runs, n_steps, self.n_nodes))
        
        for run in range(n_runs):
            if (run + 1) % max(1, n_runs // 5) == 0:
                logger.info(f"  Progress: {run + 1}/{n_runs} runs completed")
            
            I = self.I0.copy().astype(float)
            
            for t in range(n_steps):
                infection_history[run, t, :] = I
                
                # Transmission step
                I_new = I.copy()
                for i in range(self.n_nodes):
                    if I[i] > 0:  # Node i is infected
                        for j in range(self.n_nodes):
                            if i != j and self.A[i, j] > 0 and I[j] == 0:  # j is susceptible, edge exists
                                # Transmission probability
                                p_transmit = self.beta * self.A[i, j] * I[i]
                                if np.random.rand() < p_transmit:
                                    I_new[j] = 1.0
                
                # Recovery step
                for i in range(self.n_nodes):
                    if I_new[i] > 0:  # Node is infected
                        if np.random.rand() < self.gamma:
                            I_new[i] = 0.0  # Recovery
                
                I = I_new
        
        self.infection_history = infection_history
        logger.info(f"✓ Dynamics completed")
        
        # Compute persistence scores
        self._compute_persistence_scores()
        
        # Identify biomarkers
        self._identify_biomarkers()
        
        # Compile results
        results = {
            'infection_history': infection_history,
            'persistence_scores': self.persistence_scores,
            'biomarkers': self.biomarkers,
            'statistics': self._compute_statistics(infection_history)
        }
        
        return results
    
    def _compute_persistence_scores(self) -> np.ndarray:
        """
        Compute persistence score for each gene
        
        Persistence = (# timesteps infected) / (n_runs × n_steps)
        
        Returns:
        --------
        np.ndarray : Persistence scores ∈ [0, 1]
        """
        if self.infection_history is None:
            raise ValueError("Run dynamics first")
        
        n_runs, n_steps, n_nodes = self.infection_history.shape
        persistence = np.zeros(n_nodes)
        
        for i in range(n_nodes):
            # Count timesteps where node i is infected
            n_infected = (self.infection_history[:, :, i] > 0).sum()
            persistence[i] = n_infected / (n_runs * n_steps)
        
        self.persistence_scores = persistence
        
        logger.info(f"Computed persistence scores:")
        logger.info(f"  Mean: {persistence.mean():.3f}")
        logger.info(f"  Std: {persistence.std():.3f}")
        logger.info(f"  Max: {persistence.max():.3f}")
        
        return persistence
    
    def _identify_biomarkers(self, percentile: float = 75.0,
                            min_persistence: float = 0.3) -> List[str]:
        """
        Identify genes as biomarkers based on persistence
        
        Logic:
        1. Genes with persistence > percentile (e.g., 75th) are biomarkers
        2. AND persistence > min_persistence threshold
        3. Ranked by persistence score
        
        Parameters:
        -----------
        percentile : float
            Percentile threshold (default 75)
        min_persistence : float
            Minimum absolute persistence threshold (default 0.3)
        
        Returns:
        --------
        List[str] : Biomarker gene names, ranked by persistence
        """
        if self.persistence_scores is None:
            raise ValueError("Run dynamics and compute persistence first")
        
        # Compute threshold
        threshold = np.percentile(self.persistence_scores, percentile)
        threshold = max(threshold, min_persistence)
        
        # Find biomarkers
        biomarker_indices = np.where(self.persistence_scores >= threshold)[0]
        biomarker_genes = [self.gene_names[i] for i in biomarker_indices]
        
        # Sort by persistence (descending)
        biomarker_persistence = self.persistence_scores[biomarker_indices]
        sorted_idx = np.argsort(-biomarker_persistence)
        biomarker_genes = [biomarker_genes[i] for i in sorted_idx]
        
        self.biomarkers = biomarker_genes
        
        logger.info(f"Identified {len(biomarker_genes)} biomarkers:")
        logger.info(f"  Persistence threshold: {threshold:.3f}")
        logger.info(f"  Top 5 biomarkers:")
        for i, gene in enumerate(biomarker_genes[:5]):
            idx = self.gene_names.index(gene)
            logger.info(f"    {i+1}. {gene}: {self.persistence_scores[idx]:.3f}")
        
        return biomarker_genes
    
    def _compute_statistics(self, infection_history: np.ndarray) -> Dict:
        """
        Compute summary statistics of SIS dynamics
        
        Returns:
        --------
        dict : Statistics including mean/max/min infected, extinction time, etc.
        """
        n_runs, n_steps, n_nodes = infection_history.shape
        
        # Number of infected nodes over time (averaged across runs)
        n_infected_per_step = (infection_history > 0).sum(axis=2).mean(axis=0)
        
        # Final state statistics
        final_infected = (infection_history[:, -1, :] > 0).sum(axis=1)
        
        stats = {
            'mean_infected_per_step': n_infected_per_step.mean(),
            'max_infected_per_step': n_infected_per_step.max(),
            'min_infected_per_step': n_infected_per_step.min(),
            'mean_final_infected': final_infected.mean(),
            'extinction_rate': (final_infected == 0).sum() / n_runs,
            'mean_persistence': self.persistence_scores.mean(),
            'max_persistence': self.persistence_scores.max(),
            'n_biomarkers': len(self.biomarkers) if self.biomarkers else 0
        }
        
        return stats
    
    def get_biomarker_table(self, top_n: int = 20) -> pd.DataFrame:
        """
        Get biomarkers as a ranked DataFrame
        
        Parameters:
        -----------
        top_n : int
            Number of top biomarkers to return
        
        Returns:
        --------
        pd.DataFrame : Biomarker ranking table
        """
        if self.biomarkers is None or self.persistence_scores is None:
            raise ValueError("Run dynamics first")
        
        biomarker_data = []
        for i, gene in enumerate(self.biomarkers[:top_n]):
            gene_idx = self.gene_names.index(gene)
            
            # Compute network metrics
            gene_degree = (self.A[gene_idx, :] > 0).sum()
            gene_strength = self.A[gene_idx, :].sum()
            
            biomarker_data.append({
                'rank': i + 1,
                'gene': gene,
                'persistence': self.persistence_scores[gene_idx],
                'degree': gene_degree,
                'strength': gene_strength,
                'I0': self.I0[gene_idx]
            })
        
        biomarker_df = pd.DataFrame(biomarker_data)
        return biomarker_df
    
    def get_infection_dynamics(self) -> Dict[str, np.ndarray]:
        """
        Get infection dynamics for visualization
        
        Returns:
        --------
        dict : {
            'mean_infection': np.ndarray (n_steps,),  # Mean infection fraction
            'infection_std': np.ndarray (n_steps,),  # Std of infection fraction
            'n_infected_per_step': np.ndarray (n_steps,)  # Mean # infected genes
        }
        """
        if self.infection_history is None:
            raise ValueError("Run dynamics first")
        
        n_runs, n_steps, n_nodes = self.infection_history.shape
        
        # Infection fraction per step
        infection_fraction = (self.infection_history > 0).sum(axis=2) / n_nodes
        mean_fraction = infection_fraction.mean(axis=0)
        std_fraction = infection_fraction.std(axis=0)
        n_infected = infection_fraction.mean(axis=0) * n_nodes
        
        return {
            'mean_infection': mean_fraction,
            'infection_std': std_fraction,
            'n_infected_per_step': n_infected
        }
    
    def compute_weighted_biomarker_score(self, centrality_weights: Optional[np.ndarray] = None) -> pd.DataFrame:
        """
        Compute weighted biomarker score combining persistence + network centrality
        
        Score = 0.6 × persistence + 0.4 × normalized_centrality
        
        Parameters:
        -----------
        centrality_weights : np.ndarray or None
            Pre-computed centrality scores (e.g., betweenness centrality)
        
        Returns:
        --------
        pd.DataFrame : Genes ranked by weighted score
        """
        if self.persistence_scores is None:
            raise ValueError("Run dynamics first")
        
        # Compute centrality if not provided
        if centrality_weights is None:
            # Use degree as simple centrality metric
            centrality_weights = (self.A > 0).sum(axis=1)
        
        # Normalize both to [0, 1]
        persistence_norm = (self.persistence_scores - self.persistence_scores.min()) / \
                          (self.persistence_scores.max() - self.persistence_scores.min() + 1e-6)
        centrality_norm = (centrality_weights - centrality_weights.min()) / \
                         (centrality_weights.max() - centrality_weights.min() + 1e-6)
        
        # Compute weighted score
        weighted_score = 0.6 * persistence_norm + 0.4 * centrality_norm
        
        # Create results table
        results_data = []
        for i in range(len(self.gene_names)):
            results_data.append({
                'gene': self.gene_names[i],
                'persistence': self.persistence_scores[i],
                'centrality': centrality_weights[i],
                'weighted_score': weighted_score[i],
                'rank': 0
            })
        
        results_df = pd.DataFrame(results_data)
        results_df = results_df.sort_values('weighted_score', ascending=False)
        results_df['rank'] = range(1, len(results_df) + 1)
        
        return results_df


if __name__ == "__main__":
    print("SIS Network Propagation Engine")
    print("=" * 60)
    print("\nUsage:")
    print("  from sis_network_propagation import SISNetworkPropagation")
    print("  propagator = SISNetworkPropagation(adjacency, gene_names, parameters)")
    print("  results = propagator.run_dynamics(n_runs=100, n_steps=1000)")
    print("  biomarkers_df = propagator.get_biomarker_table()")
