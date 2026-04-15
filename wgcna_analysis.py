"""
WGCNA (Weighted Gene Co-expression Network Analysis) Module

Implements weighted gene co-expression network construction and analysis:
1. Soft power threshold optimization
2. Weighted correlation network construction
3. Dynamic tree cutting for module identification
4. Module eigengene computation
5. Module-trait correlation analysis
6. Hub gene identification within co-expression modules

This bridges Phase 1 pathway activity to disease modules from Phase 2 Task 1.
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.preprocessing import StandardScaler
from typing import Dict, List, Tuple, Optional
import logging
import warnings

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class WGCNAAnalyzer:
    """
    Weighted Gene Co-expression Network Analysis
    
    Implements the core WGCNA methodology:
    - Constructs weighted correlation network from gene expression
    - Identifies co-expression modules using hierarchical clustering
    - Computes module eigengenes (principal component per module)
    - Correlates modules with clinical traits
    """
    
    def __init__(self, expr_data: pd.DataFrame):
        """
        Initialize WGCNA analyzer
        
        Parameters:
        -----------
        expr_data : pd.DataFrame
            Gene expression matrix (genes × samples)
            Rows: gene names, Columns: sample names
            Values: log2-transformed expression counts
        """
        self.expr_data = expr_data
        self.n_genes, self.n_samples = expr_data.shape
        self.soft_power = None
        self.adjacency = None
        self.tom = None  # Topological Overlap Matrix
        self.modules = None  # gene -> module assignment
        self.eigengenes = None  # Module eigengenes (module × sample)
        self.hub_genes = {}  # module -> hub genes
        self.module_colors = None
        
        logger.info(f"Initialized WGCNAAnalyzer: {self.n_genes} genes × {self.n_samples} samples")
    
    def select_soft_power(self, power_range: range = range(1, 21),
                         correlation_type: str = 'pearson',
                         target_r2: float = 0.9) -> Tuple[int, np.ndarray]:
        """
        Select soft power threshold for network construction
        
        Uses scale-free topology fitting to find power that best fits
        scale-free network model (R² > target_r2).
        
        Parameters:
        -----------
        power_range : range
            Range of powers to test (default: 1-20)
        correlation_type : str
            'pearson' or 'spearman'
        target_r2 : float
            Target R² for scale-free fit (default: 0.9)
            
        Returns:
        --------
        tuple : (optimal_power, r2_values)
        """
        logger.info(f"Selecting soft power (testing powers {power_range.start}-{power_range.stop-1})...")
        
        # Compute correlation matrix
        if correlation_type == 'pearson':
            corr_matrix = self.expr_data.T.corr(method='pearson').values
        else:
            corr_matrix = self.expr_data.T.corr(method='spearman').values
        
        # Make sure correlations are in [0, 1]
        corr_matrix = np.abs(corr_matrix)
        np.fill_diagonal(corr_matrix, 0)
        
        r2_values = []
        
        for power in power_range:
            # Raise correlation to power to get adjacency
            adjacency = corr_matrix ** power
            
            # Compute connectivity (degree)
            connectivity = adjacency.sum(axis=1)
            
            # Scale-free fitting: log(connectivity) vs log(frequency)
            # Bin connectivities
            k_min = connectivity.min()
            k_max = connectivity.max()
            bins = np.logspace(np.log10(k_min + 0.001), np.log10(k_max + 1), 50)
            
            pk, _ = np.histogram(connectivity, bins=bins)
            pk = pk / (pk.sum() + 1e-10)  # Normalize
            
            # Fit log-log relationship
            bin_centers = (bins[:-1] + bins[1:]) / 2
            
            # Remove zero counts
            valid_idx = pk > 0
            if valid_idx.sum() < 2:
                r2_values.append(0.0)
                continue
            
            x = np.log10(bin_centers[valid_idx])
            y = np.log10(pk[valid_idx])
            
            # Linear fit
            coeffs = np.polyfit(x, y, 1)
            y_fit = np.polyval(coeffs, x)
            
            # Compute R²
            ss_res = ((y - y_fit) ** 2).sum()
            ss_tot = ((y - y.mean()) ** 2).sum()
            r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            r2_values.append(r2)
            logger.info(f"  Power {power}: R² = {r2:.3f}")
        
        r2_values = np.array(r2_values)
        
        # Select power with R² closest to target
        if (r2_values >= target_r2).any():
            optimal_power = power_range[np.where(r2_values >= target_r2)[0][0]]
            logger.info(f"✓ Selected soft power: {optimal_power} (R² = {r2_values[optimal_power-power_range.start]:.3f})")
        else:
            # Use power with highest R²
            optimal_power = power_range[np.argmax(r2_values)]
            logger.info(f"⚠️  No power reached R² > {target_r2}. Using power {optimal_power} (R² = {max(r2_values):.3f})")
        
        self.soft_power = optimal_power
        return optimal_power, r2_values
    
    def build_network(self, soft_power: Optional[int] = None,
                     correlation_type: str = 'pearson') -> np.ndarray:
        """
        Build weighted gene correlation network
        
        Parameters:
        -----------
        soft_power : int or None
            Power for adjacency computation. If None, will select automatically.
        correlation_type : str
            'pearson' or 'spearman'
            
        Returns:
        --------
        np.ndarray : Adjacency matrix (genes × genes)
        """
        if soft_power is None:
            self.select_soft_power(correlation_type=correlation_type)
            soft_power = self.soft_power
        else:
            self.soft_power = soft_power
        
        logger.info(f"Building network (soft power = {soft_power})...")
        
        # Compute correlation matrix
        if correlation_type == 'pearson':
            corr_matrix = self.expr_data.T.corr(method='pearson').values
        else:
            corr_matrix = self.expr_data.T.corr(method='spearman').values
        
        # Make sure correlations are in [0, 1]
        corr_matrix = np.abs(corr_matrix)
        np.fill_diagonal(corr_matrix, 0)
        
        # Compute adjacency
        adjacency = corr_matrix ** soft_power
        self.adjacency = adjacency
        
        logger.info(f"✓ Network built: {adjacency.shape[0]} nodes, " +
                   f"Mean adjacency = {adjacency.mean():.4f}")
        
        return adjacency
    
    def compute_tom(self, adjacency: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute Topological Overlap Matrix (TOM)
        
        TOM measures how many neighbors two genes share,
        providing more robust distance for clustering.
        
        Parameters:
        -----------
        adjacency : np.ndarray or None
            Adjacency matrix. If None, uses stored adjacency.
            
        Returns:
        --------
        np.ndarray : TOM matrix (genes × genes)
        """
        if adjacency is None:
            if self.adjacency is None:
                raise ValueError("Adjacency matrix not computed. Call build_network first.")
            adjacency = self.adjacency
        
        logger.info("Computing Topological Overlap Matrix...")
        
        n = adjacency.shape[0]
        tom = np.zeros_like(adjacency)
        
        # TOM[i,j] = (w[i,j] + sum_k(w[i,k]*w[k,j])) / (min(k[i], k[j]) + 1 - w[i,j])
        k = adjacency.sum(axis=1)  # Connectivity of each gene
        
        for i in range(n):
            if (i + 1) % max(1, n // 10) == 0:
                logger.info(f"  Processed {i + 1}/{n} genes...")
            
            for j in range(i, n):
                overlap = np.dot(adjacency[i, :], adjacency[j, :])
                denom = min(k[i], k[j]) + 1 - adjacency[i, j]
                tom[i, j] = (adjacency[i, j] + overlap) / denom if denom > 0 else 0
                tom[j, i] = tom[i, j]
        
        self.tom = tom
        logger.info(f"✓ TOM computed: {tom.shape}")
        
        return tom
    
    def identify_modules(self, distance_threshold: float = 0.25,
                        min_module_size: int = 30,
                        merge_threshold: float = 0.85) -> Dict[str, np.ndarray]:
        """
        Identify co-expression modules using hierarchical clustering
        and dynamic tree cutting
        
        Parameters:
        -----------
        distance_threshold : float
            Height threshold for cutting dendrogram (default: 0.25)
        min_module_size : int
            Minimum genes per module (default: 30)
        merge_threshold : float
            Correlation threshold for merging similar modules (default: 0.85)
            
        Returns:
        --------
        dict : module_color -> list of gene indices
        """
        logger.info("Identifying co-expression modules...")
        
        if self.tom is None:
            self.compute_tom()
        
        # Convert TOM to distance
        distance = 1 - self.tom
        np.fill_diagonal(distance, 0)
        
        # Hierarchical clustering
        logger.info("  Performing hierarchical clustering...")
        # Use only upper triangle for linkage
        from scipy.spatial.distance import squareform
        distance_condensed = squareform(distance, checks=False)
        linkage_matrix = linkage(distance_condensed, method='average')
        
        # Dynamic tree cutting using fcluster (fixed API)
        logger.info("  Applying dynamic tree cutting...")
        # fcluster requires criterion parameter
        clusters = fcluster(linkage_matrix, t=distance_threshold, criterion='distance')
        
        # Filter small modules and assign colors
        modules = {}
        color_map = {}
        color_palette = plt_to_color_palette(max(clusters) + 2)
        
        color_idx = 0
        for cluster_id in np.unique(clusters):
            genes_in_cluster = np.where(clusters == cluster_id)[0]
            
            if len(genes_in_cluster) >= min_module_size:
                color = color_palette[color_idx % len(color_palette)]
                modules[color] = genes_in_cluster
                color_map[cluster_id] = color
                color_idx += 1
        
        logger.info(f"  Found {len(modules)} modules (min size: {min_module_size})")
        
        # Merge highly correlated modules
        logger.info("  Merging similar modules...")
        modules = self._merge_modules(modules, merge_threshold)
        
        logger.info(f"✓ Final modules: {len(modules)}")
        
        # Create module assignment
        module_assignment = {}
        gene_names = self.expr_data.index.tolist()
        
        for color, gene_indices in modules.items():
            for idx in gene_indices:
                module_assignment[gene_names[idx]] = color
        
        self.modules = module_assignment
        self.module_colors = list(modules.keys())
        
        return modules
    
    def _merge_modules(self, modules: Dict[str, np.ndarray],
                      threshold: float = 0.85) -> Dict[str, np.ndarray]:
        """Merge modules with correlation > threshold"""
        
        # Compute module eigengenes
        eigengenes = {}
        for color, gene_indices in modules.items():
            # Principal component of module expression
            X = self.expr_data.iloc[gene_indices].T.values
            if X.shape[0] < X.shape[1]:
                # More genes than samples; use mean instead
                eigengenes[color] = X.mean(axis=1)
            else:
                # Compute PC1
                from sklearn.decomposition import PCA
                pca = PCA(n_components=1)
                eigengenes[color] = pca.fit_transform(X).flatten()
        
        # Compute pairwise correlations
        merged = False
        colors = list(modules.keys())
        
        for i in range(len(colors)):
            for j in range(i + 1, len(colors)):
                color_i = colors[i]
                color_j = colors[j]
                
                if color_i not in modules or color_j not in modules:
                    continue
                
                # Compute correlation
                corr = np.abs(np.corrcoef(eigengenes[color_i], eigengenes[color_j])[0, 1])
                
                if corr > threshold:
                    # Merge j into i
                    modules[color_i] = np.concatenate([
                        modules[color_i],
                        modules[color_j]
                    ])
                    del modules[color_j]
                    merged = True
        
        if merged:
            return self._merge_modules(modules, threshold)
        
        return modules
    
    def compute_eigengenes(self) -> pd.DataFrame:
        """
        Compute module eigengenes (principal component per module)
        
        Returns:
        --------
        pd.DataFrame : Eigengenes matrix (modules × samples)
        """
        logger.info("Computing module eigengenes...")
        
        if self.modules is None:
            raise ValueError("Modules not identified. Call identify_modules first.")
        
        from sklearn.decomposition import PCA
        
        eigengenes = {}
        module_genes = {}
        
        # Group genes by module
        for gene, module_color in self.modules.items():
            if module_color not in module_genes:
                module_genes[module_color] = []
            module_genes[module_color].append(gene)
        
        # Compute eigengene for each module
        for color, genes in module_genes.items():
            expr_subset = self.expr_data.loc[genes].T  # samples × genes
            
            if expr_subset.shape[1] >= 2:
                pca = PCA(n_components=1)
                eigengene = pca.fit_transform(expr_subset.values).flatten()
                eigengenes[color] = eigengene
            else:
                # Single gene; just use its expression
                eigengenes[color] = expr_subset.iloc[:, 0].values
        
        self.eigengenes = pd.DataFrame(
            eigengenes,
            index=self.expr_data.columns,
            columns=list(eigengenes.keys())
        )
        
        logger.info(f"✓ Computed eigengenes: {self.eigengenes.shape}")
        
        return self.eigengenes
    
    def correlate_with_traits(self, traits_data: pd.DataFrame,
                             min_correlation: float = 0.3) -> pd.DataFrame:
        """
        Correlate module eigengenes with clinical traits
        
        Parameters:
        -----------
        traits_data : pd.DataFrame
            Clinical traits (samples × traits)
            Rows must match expression sample order
        min_correlation : float
            Minimum correlation to report
            
        Returns:
        --------
        pd.DataFrame : Module-trait correlations with p-values
        """
        logger.info("Correlating module eigengenes with traits...")
        
        if self.eigengenes is None:
            self.compute_eigengenes()
        
        # Ensure trait data is numeric
        traits_numeric = pd.DataFrame()
        for col in traits_data.columns:
            if pd.api.types.is_numeric_dtype(traits_data[col]):
                traits_numeric[col] = traits_data[col]
            else:
                # Encode categorical as numeric
                traits_numeric[col] = pd.factorize(traits_data[col])[0]
        
        results = []
        
        for module_color in self.eigengenes.columns:
            module_eigengene = self.eigengenes[module_color].values
            
            for trait in traits_numeric.columns:
                trait_values = traits_numeric[trait].values
                
                # Compute correlation and p-value
                corr, pval = pearsonr(module_eigengene, trait_values)
                
                if abs(corr) >= min_correlation:
                    results.append({
                        'module': module_color,
                        'trait': trait,
                        'correlation': corr,
                        'p_value': pval,
                        'significant': pval < 0.05
                    })
        
        correlation_df = pd.DataFrame(results)
        
        if len(correlation_df) > 0:
            # FDR correction
            from scipy.stats import rankdata
            n_tests = len(correlation_df)
            correlation_df['padj'] = correlation_df['p_value'] * n_tests  # Bonferroni
            correlation_df['padj'] = np.minimum(correlation_df['padj'], 1.0)
        
        logger.info(f"✓ Found {len(correlation_df)} module-trait correlations (|r| >= {min_correlation})")
        
        return correlation_df
    
    def identify_hub_genes(self, module_color: str, top_n: int = 10) -> pd.DataFrame:
        """
        Identify hub genes in a module
        
        Hub genes = highest connectivity within module
        
        Parameters:
        -----------
        module_color : str
            Module color
        top_n : int
            Number of top hub genes to return
            
        Returns:
        --------
        pd.DataFrame : Hub genes with connectivity scores
        """
        if self.modules is None or module_color not in [m for m in self.modules.values()]:
            raise ValueError(f"Module {module_color} not found")
        
        # Get genes in module
        gene_names = self.expr_data.index.tolist()
        module_genes_idx = [i for i, (gene, color) in enumerate([(gene_names[i], self.modules.get(gene_names[i])) for i in range(len(gene_names))]) if color == module_color]
        module_genes = [gene_names[i] for i in module_genes_idx]
        
        if self.adjacency is None:
            raise ValueError("Adjacency matrix not computed")
        
        # Compute connectivity within module
        hub_scores = []
        for gene_idx in module_genes_idx:
            # Sum of adjacencies to other genes in same module
            connectivity = self.adjacency[gene_idx, module_genes_idx].sum()
            
            hub_scores.append({
                'gene': gene_names[gene_idx],
                'connectivity': connectivity,
                'rank': 0
            })
        
        hub_df = pd.DataFrame(hub_scores)
        hub_df = hub_df.sort_values('connectivity', ascending=False)
        hub_df['rank'] = range(1, len(hub_df) + 1)
        
        return hub_df.head(top_n)


class ModuleTraitCorrelation:
    """
    Analyzes correlation between co-expression modules and clinical traits
    """
    
    def __init__(self, eigengenes: pd.DataFrame, traits_data: pd.DataFrame):
        """
        Initialize module-trait analysis
        
        Parameters:
        -----------
        eigengenes : pd.DataFrame
            Module eigengenes (samples × modules)
        traits_data : pd.DataFrame
            Clinical traits (samples × traits)
        """
        self.eigengenes = eigengenes
        self.traits_data = traits_data
        self.correlations = None
    
    def compute_correlations(self, method: str = 'pearson') -> pd.DataFrame:
        """
        Compute module-trait correlations
        
        Parameters:
        -----------
        method : str
            'pearson' or 'spearman'
            
        Returns:
        --------
        pd.DataFrame : Correlation matrix (modules × traits)
        """
        logger.info(f"Computing module-trait correlations ({method})...")
        
        # Convert traits to numeric
        traits_numeric = pd.DataFrame()
        for col in self.traits_data.columns:
            if pd.api.types.is_numeric_dtype(self.traits_data[col]):
                traits_numeric[col] = self.traits_data[col]
            else:
                traits_numeric[col] = pd.factorize(self.traits_data[col])[0]
        
        if method == 'pearson':
            correlations = self.eigengenes.corrwith(traits_numeric, method='pearson')
        else:
            correlations = self.eigengenes.corrwith(traits_numeric, method='spearman')
        
        self.correlations = correlations
        
        logger.info(f"✓ Computed {len(correlations)} correlations")
        
        return correlations


class ModuleComparison:
    """
    Compares WGCNA co-expression modules with disease modules
    """
    
    def __init__(self, wgcna_modules: Dict[str, List[str]],
                 disease_modules: Dict[str, List[str]]):
        """
        Initialize module comparison
        
        Parameters:
        -----------
        wgcna_modules : dict
            WGCNA module_color -> list of genes
        disease_modules : dict
            Disease name -> list of genes in disease module
        """
        self.wgcna_modules = wgcna_modules
        self.disease_modules = disease_modules
    
    def compute_overlap(self) -> pd.DataFrame:
        """
        Compute overlap between WGCNA and disease modules
        
        Returns:
        --------
        pd.DataFrame : Overlap metrics (WGCNA module × disease)
        """
        logger.info("Computing module overlap...")
        
        overlaps = []
        
        for wgcna_color, wgcna_genes in self.wgcna_modules.items():
            wgcna_set = set(wgcna_genes)
            
            for disease, disease_genes in self.disease_modules.items():
                disease_set = set(disease_genes)
                
                # Jaccard similarity
                intersection = len(wgcna_set & disease_set)
                union = len(wgcna_set | disease_set)
                jaccard = intersection / union if union > 0 else 0
                
                overlaps.append({
                    'wgcna_module': wgcna_color,
                    'disease': disease,
                    'overlap_count': intersection,
                    'jaccard_similarity': jaccard,
                    'wgcna_module_size': len(wgcna_genes),
                    'disease_module_size': len(disease_genes)
                })
        
        overlap_df = pd.DataFrame(overlaps)
        logger.info(f"✓ Computed {len(overlap_df)} overlaps")
        
        return overlap_df


# Utility functions
def plt_to_color_palette(n_colors: int) -> List[str]:
    """Generate a palette of distinct colors"""
    colors = [
        'red', 'blue', 'green', 'yellow', 'purple', 'orange',
        'brown', 'pink', 'gray', 'olive', 'cyan', 'magenta',
        'teal', 'navy', 'maroon', 'lime', 'aqua', 'salmon'
    ]
    return (colors * ((n_colors // len(colors)) + 1))[:n_colors]


if __name__ == "__main__":
    print("WGCNA Analysis Module")
    print("=" * 60)
    print("\nUsage:")
    print("  from wgcna_analysis import WGCNAAnalyzer")
    print("  analyzer = WGCNAAnalyzer(expr_data)")
    print("  analyzer.select_soft_power()")
    print("  analyzer.build_network()")
    print("  analyzer.identify_modules()")
    print("  eigengenes = analyzer.compute_eigengenes()")
