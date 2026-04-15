"""
Disease Module Detection Module (Phase 2)

Implements network medicine approaches to detect disease modules and predict
comorbidities by analyzing disease-disease and disease-pathway relationships
through protein interaction networks.

Features:
- Load gene-disease associations from knowledge bases
- Map genes to PPI networks (STRING/BioGRID)
- Detect disease modules using community detection
- Compute module separation and comorbidity scores
- Integrate with Phase 1 pathway activity scores

Status: In Development
"""

import pandas as pd
import numpy as np
import networkx as nx
from typing import Dict, List, Tuple, Optional, Set
import logging
from pathlib import Path
import os
from collections import defaultdict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DiseaseNetworkBuilder:
    """
    Build and manage disease networks from gene-disease associations and PPI.
    
    This class:
    1. Loads gene-disease associations from knowledge bases
    2. Maps genes to a protein interaction network
    3. Extracts disease-specific subnetworks
    4. Prepares data for community detection
    """
    
    def __init__(self, gene_disease_file: Optional[str] = None):
        """
        Initialize disease network builder.
        
        Parameters:
        -----------
        gene_disease_file : str or None
            Path to gene-disease association file (TSV format)
            Required columns: disease_name, gene_symbol
            If None, will attempt to find 'data/gene_disease.tsv'
        """
        self.gene_disease_file = gene_disease_file or "data/gene_disease.tsv"
        self.gene_disease_map = {}  # disease -> set of genes
        self.disease_metadata = {}  # disease -> metadata dict
        self.ppi_network = None
        self.ppi_source = None
        self.disease_subnetworks = {}  # disease -> subnetwork
        
        logger.info("Initialized DiseaseNetworkBuilder")
        
    def load_gene_disease_associations(self) -> Dict[str, Set[str]]:
        """
        Load gene-disease associations from file.
        
        Returns:
        --------
        dict : disease_name -> set of gene symbols
        """
        logger.info(f"Loading gene-disease associations from {self.gene_disease_file}")
        
        if not os.path.exists(self.gene_disease_file):
            logger.error(f"File not found: {self.gene_disease_file}")
            return {}
        
        try:
            df = pd.read_csv(self.gene_disease_file, sep='\t', low_memory=False)
            logger.info(f"  Loaded {len(df)} rows")
            
            # Group by disease
            for disease_name, group in df.groupby('disease_name'):
                # Get unique genes for this disease
                genes = set()
                
                # Try gene_symbol column
                if 'gene_symbol' in group.columns:
                    gene_symbols = group['gene_symbol'].dropna()
                    genes.update([str(g).strip() for g in gene_symbols if g and str(g) != 'nan'])
                
                # Also try disease_gene column (alternative name)
                if 'disease_gene' in group.columns:
                    gene_names = group['disease_gene'].dropna()
                    genes.update([str(g).strip() for g in gene_names if g and str(g) != 'nan'])
                
                if genes:
                    self.gene_disease_map[disease_name] = genes
                    
                    # Store metadata
                    self.disease_metadata[disease_name] = {
                        'n_genes': len(genes),
                        'category': group['disease_category'].iloc[0] if 'disease_category' in group.columns else 'Unknown'
                    }
            
            logger.info(f"✓ Loaded {len(self.gene_disease_map)} unique diseases")
            logger.info(f"  Total gene-disease associations: {sum(len(g) for g in self.gene_disease_map.values())}")
            
            return self.gene_disease_map
            
        except Exception as e:
            logger.error(f"Error loading gene-disease associations: {e}")
            return {}
    
    def load_ppi_network(self, source: str = 'simple', network_file: Optional[str] = None) -> nx.Graph:
        """
        Load protein interaction network.
        
        Parameters:
        -----------
        source : str
            Network source: 'string' (download), 'biogrid', or 'simple'
            'simple' creates a minimal example network from existing genes
        network_file : str or None
            If provided, load network from this file instead
            
        Returns:
        --------
        nx.Graph : Protein interaction network
        """
        logger.info(f"Loading PPI network from source: {source}")
        
        if network_file and os.path.exists(network_file):
            logger.info(f"Loading from file: {network_file}")
            try:
                df = pd.read_csv(network_file, sep='\t', header=None, names=['protein1', 'protein2'])
                G = nx.Graph()
                for _, row in df.iterrows():
                    G.add_edge(str(row['protein1']), str(row['protein2']))
                self.ppi_network = G
                self.ppi_source = 'file'
                logger.info(f"✓ Loaded PPI network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
                return G
            except Exception as e:
                logger.warning(f"Error loading from file: {e}, will create simple network")
        
        if source == 'simple':
            # Create a simple network from existing genes for testing
            logger.info("Creating simple PPI network from existing genes")
            G = nx.Graph()
            
            # Collect all genes
            all_genes = set()
            for genes in self.gene_disease_map.values():
                all_genes.update(genes)
            
            logger.info(f"  Creating network with {len(all_genes)} genes")
            
            # Add genes as nodes
            G.add_nodes_from(all_genes)
            
            # Create some edges based on simple rules (for demonstration)
            # In real scenario, would load from STRING or BioGRID
            all_genes_list = list(all_genes)
            np.random.seed(42)
            
            # Add edges randomly but with some structure
            n_edges = min(len(all_genes_list) * 2, 1000)  # ~2 edges per gene on average
            for _ in range(n_edges):
                if len(all_genes_list) > 1:
                    g1, g2 = np.random.choice(all_genes_list, 2, replace=False)
                    G.add_edge(g1, g2)
            
            self.ppi_network = G
            self.ppi_source = 'simple'
            logger.info(f"✓ Created simple PPI network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
            return G
        
        else:
            logger.warning(f"Source '{source}' not implemented, creating simple network")
            return self.load_ppi_network('simple')
    
    def build_disease_subnetwork(self, disease: str) -> Optional[nx.Graph]:
        """
        Extract disease-specific subnetwork from PPI.
        
        Parameters:
        -----------
        disease : str
            Disease name (must be in gene_disease_map)
            
        Returns:
        --------
        nx.Graph or None : Subnetwork containing genes associated with disease
        """
        if disease not in self.gene_disease_map:
            logger.warning(f"Disease '{disease}' not found in gene-disease map")
            return None
        
        if self.ppi_network is None:
            logger.error("PPI network not loaded. Call load_ppi_network first.")
            return None
        
        disease_genes = self.gene_disease_map[disease]
        
        # Get genes that are in both disease set and PPI network
        genes_in_network = [g for g in disease_genes if g in self.ppi_network]
        
        logger.info(f"Building subnetwork for {disease}")
        logger.info(f"  Disease genes: {len(disease_genes)}, in PPI: {len(genes_in_network)}")
        
        if len(genes_in_network) < 2:
            logger.warning(f"  Only {len(genes_in_network)} genes in PPI network, cannot build subnetwork")
            return None
        
        # Extract induced subgraph
        subgraph = self.ppi_network.subgraph(genes_in_network).copy()
        
        logger.info(f"  Subnetwork: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges")
        
        self.disease_subnetworks[disease] = subgraph
        return subgraph
    
    def build_all_disease_subnetworks(self) -> Dict[str, nx.Graph]:
        """
        Build subnetworks for all diseases.
        
        Returns:
        --------
        dict : disease -> subnetwork
        """
        logger.info(f"Building subnetworks for {len(self.gene_disease_map)} diseases...")
        
        for i, disease in enumerate(self.gene_disease_map.keys()):
            if (i + 1) % max(1, len(self.gene_disease_map) // 10) == 0:
                logger.info(f"  Processed {i + 1}/{len(self.gene_disease_map)} diseases")
            
            self.build_disease_subnetwork(disease)
        
        logger.info(f"✓ Built subnetworks for {len(self.disease_subnetworks)} diseases")
        return self.disease_subnetworks
    
    def get_disease_connectivity_stats(self, disease: str) -> Dict[str, float]:
        """
        Compute network statistics for a disease subnetwork.
        
        Parameters:
        -----------
        disease : str
            Disease name
            
        Returns:
        --------
        dict : Statistics including nodes, edges, density, average clustering
        """
        if disease not in self.disease_subnetworks:
            self.build_disease_subnetwork(disease)
        
        subgraph = self.disease_subnetworks.get(disease)
        
        if subgraph is None or subgraph.number_of_nodes() < 2:
            return {}
        
        stats = {
            'n_nodes': subgraph.number_of_nodes(),
            'n_edges': subgraph.number_of_edges(),
            'density': nx.density(subgraph),
            'avg_degree': 2 * subgraph.number_of_edges() / max(1, subgraph.number_of_nodes()),
            'avg_clustering': nx.average_clustering(subgraph) if subgraph.number_of_nodes() > 2 else 0,
            'is_connected': nx.is_connected(subgraph)
        }
        
        return stats
    
    def compute_network_statistics_summary(self) -> pd.DataFrame:
        """
        Compute statistics for all disease subnetworks.
        
        Returns:
        --------
        pd.DataFrame : Statistics for each disease
        """
        stats_list = []
        
        for disease in sorted(self.disease_subnetworks.keys()):
            stats = self.get_disease_connectivity_stats(disease)
            if stats:
                stats['disease'] = disease
                stats['category'] = self.disease_metadata.get(disease, {}).get('category', 'Unknown')
                stats_list.append(stats)
        
        return pd.DataFrame(stats_list)
    
    def export_disease_network_summary(self, output_file: str):
        """
        Export disease network summary to CSV.
        
        Parameters:
        -----------
        output_file : str
            Path to output CSV file
        """
        logger.info(f"Exporting disease network summary to {output_file}")
        
        df = self.compute_network_statistics_summary()
        df.to_csv(output_file, index=False)
        
        logger.info(f"✓ Exported {len(df)} disease networks")


class CommunityDetector:
    """
    Detect communities (disease modules) in networks using multiple algorithms.
    
    Methods:
    - Louvain: Fast, scalable community detection
    - Leiden: Improved version of Louvain
    - Label Propagation: Fast alternative
    """
    
    def __init__(self):
        """Initialize community detector"""
        self.communities = {}  # network_id -> list of communities
        self.modularity_scores = {}
        
        logger.info("Initialized CommunityDetector")
    
    def detect_communities_louvain(self, network: nx.Graph, seed: int = 42) -> List[Set[str]]:
        """
        Detect communities using Louvain algorithm.
        
        Parameters:
        -----------
        network : nx.Graph
            Input network
        seed : int
            Random seed for reproducibility
            
        Returns:
        --------
        list : List of communities (each community is a set of nodes)
        """
        try:
            import community as community_louvain
        except ImportError:
            logger.warning("python-louvain not installed, using label propagation instead")
            return self.detect_communities_label_propagation(network)
        
        logger.info("Detecting communities with Louvain algorithm")
        
        # Detect communities
        partition = community_louvain.best_partition(network)
        
        # Convert partition to communities
        communities = defaultdict(set)
        for node, comm_id in partition.items():
            communities[comm_id].add(node)
        
        communities = list(communities.values())
        modularity = community_louvain.modularity(partition, network)
        
        logger.info(f"✓ Detected {len(communities)} communities, modularity: {modularity:.3f}")
        
        return communities
    
    def detect_communities_label_propagation(self, network: nx.Graph) -> List[Set[str]]:
        """
        Detect communities using label propagation.
        
        Parameters:
        -----------
        network : nx.Graph
            Input network
            
        Returns:
        --------
        list : List of communities (each community is a set of nodes)
        """
        logger.info("Detecting communities with label propagation")
        
        # Use NetworkX label propagation
        try:
            from networkx.algorithms import community
            communities = list(community.label_propagation_communities(network))
            logger.info(f"✓ Detected {len(communities)} communities")
            return communities
        except Exception as e:
            logger.warning(f"Label propagation failed: {e}, using greedy modularity")
            return self.detect_communities_greedy(network)
    
    def detect_communities_greedy(self, network: nx.Graph) -> List[Set[str]]:
        """
        Detect communities using greedy modularity optimization.
        
        Parameters:
        -----------
        network : nx.Graph
            Input network
            
        Returns:
        --------
        list : List of communities (each community is a set of nodes)
        """
        logger.info("Detecting communities with greedy modularity optimization")
        
        from networkx.algorithms import community
        communities = list(community.greedy_modularity_communities(network))
        
        logger.info(f"✓ Detected {len(communities)} communities")
        return communities
    
    def compute_community_metrics(self, network: nx.Graph, communities: List[Set[str]]) -> pd.DataFrame:
        """
        Compute metrics for detected communities.
        
        Parameters:
        -----------
        network : nx.Graph
            Input network
        communities : list
            List of communities (each is a set of nodes)
            
        Returns:
        --------
        pd.DataFrame : Community metrics
        """
        metrics = []
        
        for comm_id, community_nodes in enumerate(communities):
            subgraph = network.subgraph(community_nodes)
            
            # Compute metrics
            metric = {
                'community_id': comm_id,
                'n_nodes': len(community_nodes),
                'n_edges': subgraph.number_of_edges(),
                'density': nx.density(subgraph) if len(community_nodes) > 1 else 0,
                'avg_degree': 2 * subgraph.number_of_edges() / max(1, len(community_nodes)),
                'avg_clustering': nx.average_clustering(subgraph) if len(community_nodes) > 2 else 0
            }
            
            metrics.append(metric)
        
        return pd.DataFrame(metrics)


if __name__ == "__main__":
    print("Disease Module Detection Module")
    print("=" * 60)
    print("\nUsage:")
    print("  from disease_module_detection import DiseaseNetworkBuilder")
    print("  builder = DiseaseNetworkBuilder('data/gene_disease.tsv')")
    print("  builder.load_gene_disease_associations()")
    print("  builder.load_ppi_network('simple')")
    print("  builder.build_all_disease_subnetworks()")


class ModuleSeparationMetrics:
    """
    Compute module separation and comorbidity metrics between diseases.
    
    Based on Menche et al. (2015) "Uncovering disease-disease relationships
    through the incomplete interactome" (Science).
    
    Key metric: Network separation (s_d1_d2) measures how well separated
    two disease modules are in the network.
    """
    
    def __init__(self, ppi_network: nx.Graph, disease_subnetworks: Dict[str, nx.Graph]):
        """
        Initialize module separation calculator.
        
        Parameters:
        -----------
        ppi_network : nx.Graph
            Full PPI network
        disease_subnetworks : dict
            Mapping of disease -> disease subnetwork
        """
        self.ppi_network = ppi_network
        self.disease_subnetworks = disease_subnetworks
        self.separation_matrix = None
        self.comorbidity_scores = None
        
        logger.info(f"Initialized ModuleSeparationMetrics with {len(disease_subnetworks)} diseases")
    
    def compute_network_separation(self, disease1: str, disease2: str) -> Optional[float]:
        """
        Compute network separation between two disease modules.
        
        Parameters:
        -----------
        disease1, disease2 : str
            Disease names
            
        Returns:
        --------
        float or None : Network separation score (lower = more connected = higher comorbidity)
        """
        if disease1 not in self.disease_subnetworks or disease2 not in self.disease_subnetworks:
            return None
        
        g1 = self.disease_subnetworks[disease1]
        g2 = self.disease_subnetworks[disease2]
        
        if g1.number_of_nodes() < 2 or g2.number_of_nodes() < 2:
            return None
        
        # Compute average shortest path between modules
        # Use -1 for infinite distances (disconnected nodes)
        
        separation_scores = []
        
        for n1 in g1.nodes():
            for n2 in g2.nodes():
                try:
                    if nx.has_path(self.ppi_network, n1, n2):
                        dist = nx.shortest_path_length(self.ppi_network, n1, n2)
                        separation_scores.append(dist)
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    separation_scores.append(np.inf)
        
        if not separation_scores:
            return np.inf
        
        # Use median separation (robust to outliers)
        separation = np.median([s for s in separation_scores if s != np.inf])
        
        return float(separation) if separation != np.inf else np.inf
    
    def compute_all_disease_pairs(self) -> pd.DataFrame:
        """
        Compute separation for all disease pairs.
        
        Returns:
        --------
        pd.DataFrame : Pairwise disease separations
        """
        diseases = list(self.disease_subnetworks.keys())
        n_diseases = len(diseases)
        
        logger.info(f"Computing separation for {n_diseases} diseases ({n_diseases*(n_diseases-1)//2} pairs)")
        
        pairs = []
        for i, d1 in enumerate(diseases):
            if (i + 1) % max(1, n_diseases // 10) == 0:
                logger.info(f"  Processed {i + 1}/{n_diseases} diseases")
            
            for d2 in diseases[i+1:]:
                separation = self.compute_network_separation(d1, d2)
                if separation is not None and separation != np.inf:
                    pairs.append({
                        'disease1': d1,
                        'disease2': d2,
                        'separation': separation
                    })
        
        df = pd.DataFrame(pairs)
        
        if len(df) > 0:
            logger.info(f"✓ Computed separations for {len(df)} disease pairs")
            logger.info(f"  Separation range: [{df['separation'].min():.3f}, {df['separation'].max():.3f}]")
            self.separation_matrix = df
            return df
        else:
            logger.warning("No disease pairs with finite separation")
            return pd.DataFrame()
    
    def compute_comorbidity_scores(self, separation_matrix: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Convert network separation to comorbidity scores.
        
        Lower separation → higher comorbidity (more likely to co-occur)
        
        Parameters:
        -----------
        separation_matrix : pd.DataFrame or None
            Pairwise separations. If None, compute from all pairs.
            
        Returns:
        --------
        pd.DataFrame : Comorbidity scores (inverse of separation)
        """
        if separation_matrix is None:
            if self.separation_matrix is None:
                self.compute_all_disease_pairs()
            separation_matrix = self.separation_matrix
        
        if len(separation_matrix) == 0:
            return pd.DataFrame()
        
        # Normalize separation to [0, 1] comorbidity score
        min_sep = separation_matrix['separation'].min()
        max_sep = separation_matrix['separation'].max()
        
        if max_sep == min_sep:
            separation_matrix['comorbidity'] = 0.5
        else:
            # Invert: high separation → low comorbidity
            separation_matrix['comorbidity'] = 1 - (
                (separation_matrix['separation'] - min_sep) / (max_sep - min_sep)
            )
        
        self.comorbidity_scores = separation_matrix.sort_values('comorbidity', ascending=False)
        
        logger.info(f"✓ Computed comorbidity scores for {len(self.comorbidity_scores)} pairs")
        
        return self.comorbidity_scores
    
    def get_comorbidities_for_disease(self, disease: str, top_n: int = 20) -> pd.DataFrame:
        """
        Get predicted comorbidities for a specific disease.
        
        Parameters:
        -----------
        disease : str
            Disease name
        top_n : int
            Number of top comorbidities to return
            
        Returns:
        --------
        pd.DataFrame : Top comorbid diseases
        """
        if self.comorbidity_scores is None or len(self.comorbidity_scores) == 0:
            self.compute_comorbidity_scores()
        
        if len(self.comorbidity_scores) == 0:
            return pd.DataFrame()
        
        # Find all pairs involving this disease
        mask = (self.comorbidity_scores['disease1'] == disease) | \
               (self.comorbidity_scores['disease2'] == disease)
        
        result = self.comorbidity_scores[mask].copy()
        
        # Ensure disease1 is always the query disease for consistency
        mask_swap = result['disease2'] == disease
        result.loc[mask_swap, ['disease1', 'disease2']] = \
            result.loc[mask_swap, ['disease2', 'disease1']].values
        
        return result[['disease1', 'disease2', 'separation', 'comorbidity']].head(top_n)
    
    def predict_comorbidities(self, threshold: float = 0.5) -> pd.DataFrame:
        """
        Predict comorbidities above a threshold.
        
        Parameters:
        -----------
        threshold : float
            Comorbidity threshold for prediction
            
        Returns:
        --------
        pd.DataFrame : Predicted comorbidities
        """
        if self.comorbidity_scores is None or len(self.comorbidity_scores) == 0:
            self.compute_comorbidity_scores()
        
        if len(self.comorbidity_scores) == 0:
            return pd.DataFrame()
        
        predicted = self.comorbidity_scores[self.comorbidity_scores['comorbidity'] >= threshold]
        
        logger.info(f"✓ Predicted {len(predicted)} comorbidities with score >= {threshold}")
        
        return predicted.sort_values('comorbidity', ascending=False)


# Export convenience function
def analyze_disease_modules(gene_disease_file: str = "data/gene_disease.tsv",
                           n_diseases: Optional[int] = None) -> Tuple[DiseaseNetworkBuilder, CommunityDetector]:
    """
    Convenience function to set up and analyze disease modules.
    
    Parameters:
    -----------
    gene_disease_file : str
        Path to gene-disease association file
    n_diseases : int or None
        Limit analysis to top N diseases (for quick testing)
        
    Returns:
    --------
    tuple : (DiseaseNetworkBuilder, CommunityDetector)
    """
    builder = DiseaseNetworkBuilder(gene_disease_file)
    builder.load_gene_disease_associations()
    builder.load_ppi_network('simple')
    
    diseases = list(builder.gene_disease_map.keys())
    if n_diseases:
        diseases = diseases[:n_diseases]
    
    logger.info(f"Analyzing {len(diseases)} diseases")
    
    for disease in diseases:
        builder.build_disease_subnetwork(disease)
    
    return builder, CommunityDetector()
