# Multi-Scale Platform: Implementation Quick-Start Guide

**Last Updated**: 2026-04-14  
**Scope**: Phase 1 - Pathway Activity Analysis Implementation

---

## Overview

This guide provides step-by-step instructions for implementing the first phase of the multi-scale architecture enhancement: **Pathway Activity Analysis**.

This phase will:
1. Score 347 KEGG pathways across 255 TCGA-COAD patients
2. Identify pathway activity differences by clinical strata (age, sex, disease stage)
3. Integrate hub gene identification with pathway scores
4. Create visualizations showing pathway-scale biology

**Estimated Time**: 2-3 weeks (development) + 1 week (testing)

---

## Part 1: Development Environment Setup

### Step 1.1: Install Required Libraries

```bash
cd /Users/jaber/Downloads/webapp

# Create or activate virtual environment
python3 -m venv venv_phase1
source venv_phase1/bin/activate

# Install dependencies
pip install --upgrade pip setuptools wheel

# Core requirements (already in project)
pip install pandas numpy scipy scikit-learn networkx plotly gradio

# NEW: Pathway analysis
pip install gseapy
pip install scanpy  # Optional, for advanced functionality

# Statistics and data science
pip install statsmodels pingouin  # Better stats functions
pip install seaborn matplotlib  # For exploratory plots

# Development
pip install jupyter  # For testing
pip install ipython

# Verification
python -c "import gseapy; print('✓ gseapy installed')"
python -c "import scanpy; print('✓ scanpy installed')"
```

### Step 1.2: Verify Data Files

```bash
# Check TCGA data is present
ls -lh data/TCGA-COAD/filtered_hiseq_data.csv
ls -lh data/TCGA-COAD/filtered_clinical.csv
ls -lh data/TCGA-COAD/filtered_miRNA_with_names.csv
ls -lh data/TCGA-COAD/filtered_methylation_data.csv

# Check pathway data
ls -lh data/pathway* data/gene_disease.tsv

# Verify preprocessed data
ls -lh data/preprocessed_data.pkl
```

---

## Part 2: Implementation Plan (Modular Approach)

### Module 1: Pathway Activity Scorer

**File**: `pathway_activity.py` (NEW)

```python
"""
Pathway Activity Scoring Module
- Load pathway memberships
- Score pathways using multiple methods (GSVA, ssGSEA, mean)
- Calculate pathway activity matrices
"""

import pandas as pd
import numpy as np
from scipy.stats import zscore
import gseapy
from typing import Dict, List, Tuple

class PathwayActivityScorer:
    """Score pathway activity from gene expression data"""
    
    def __init__(self, pathway_genes: Dict[str, List[str]]):
        """
        pathway_genes: dict mapping pathway_name → list of gene symbols
        """
        self.pathway_genes = pathway_genes
        self.pathway_activity = None
        self.gene_expr = None
        
    def load_expression_data(self, expr_file: str):
        """Load TCGA gene expression matrix (genes × samples)"""
        self.gene_expr = pd.read_csv(expr_file, index_col=0)
        print(f"✓ Loaded expression: {self.gene_expr.shape}")
        
    def score_gsva(self) -> pd.DataFrame:
        """Score pathways using GSVA (Gene Set Variation Analysis)"""
        # GSVA scores gene sets from expression data
        # Returns: pathway activity matrix (pathways × samples)
        # Note: gseapy has limited GSVA; may need custom implementation
        
        activity = pd.DataFrame(index=self.pathway_genes.keys(),
                               columns=self.gene_expr.columns)
        
        for pathway, genes in self.pathway_genes.items():
            # Find intersection of pathway genes and available genes
            genes_in_expr = [g for g in genes if g in self.gene_expr.index]
            
            if len(genes_in_expr) < 2:
                activity.loc[pathway] = np.nan
                continue
                
            # Simple GSVA-like: mean z-score of pathway genes per sample
            pathway_expr = self.gene_expr.loc[genes_in_expr]
            z_scores = pathway_expr.apply(zscore, axis=0, nan_policy='omit')
            activity.loc[pathway] = z_scores.mean()
        
        self.pathway_activity = activity.astype(float)
        return self.pathway_activity
    
    def score_mean(self) -> pd.DataFrame:
        """Simple mean expression of pathway genes"""
        # Baseline method: just average expression
        
        activity = pd.DataFrame(index=self.pathway_genes.keys(),
                               columns=self.gene_expr.columns)
        
        for pathway, genes in self.pathway_genes.items():
            genes_in_expr = [g for g in genes if g in self.gene_expr.index]
            if genes_in_expr:
                activity.loc[pathway] = self.gene_expr.loc[genes_in_expr].mean()
        
        return activity.astype(float)

# Usage:
# scorer = PathwayActivityScorer(pathway_dict)
# scorer.load_expression_data('data/TCGA-COAD/filtered_hiseq_data.csv')
# pathway_activity = scorer.score_gsva()
```

**Implementation Checklist**:
- [ ] Create `pathway_activity.py`
- [ ] Implement `PathwayActivityScorer` class
- [ ] Load 347 pathways from pathway TSV
- [ ] Test on subset of TCGA data (10 pathways × 20 samples)
- [ ] Benchmark: time to score all pathways (should be <1 minute)

---

### Module 2: Differential Pathway Analysis

**File**: `differential_pathway_analysis.py` (NEW)

```python
"""
Differential Pathway Analysis Module
- Compare pathway activity across clinical strata
- Statistical testing (t-test, ANOVA, Wilcoxon)
- Multiple testing correction
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, f_oneway, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from typing import Dict, List, Tuple

class DifferentialPathwayAnalysis:
    """Compare pathway activity across patient groups"""
    
    def __init__(self, pathway_activity: pd.DataFrame, 
                 clinical_data: pd.DataFrame):
        """
        pathway_activity: pathways × samples
        clinical_data: samples × clinical variables
        """
        self.pathway_activity = pathway_activity
        self.clinical_data = clinical_data
        self.results = {}
        
    def compare_by_group(self, clinical_var: str, method='ttest') -> pd.DataFrame:
        """
        Compare pathway activity between groups
        
        Supported methods:
        - 'ttest': t-test (for 2 groups)
        - 'anova': One-way ANOVA (for 2+ groups)
        - 'wilcoxon': Mann-Whitney U (non-parametric)
        """
        
        groups = self.clinical_data[clinical_var].unique()
        results = []
        
        for pathway in self.pathway_activity.index:
            pathway_vals = self.pathway_activity.loc[pathway]
            
            # Get values per group (aligned with clinical data)
            group_vals = []
            for group in groups:
                mask = self.clinical_data[clinical_var] == group
                vals = pathway_vals[mask].dropna()
                if len(vals) > 0:
                    group_vals.append(vals)
            
            # Perform statistical test
            if method == 'ttest' and len(group_vals) == 2:
                stat, pval = ttest_ind(group_vals[0], group_vals[1])
            elif method == 'anova' and len(group_vals) > 1:
                stat, pval = f_oneway(*group_vals)
            else:
                stat, pval = np.nan, np.nan
            
            results.append({
                'pathway': pathway,
                'statistic': stat,
                'pvalue': pval,
                'group': clinical_var
            })
        
        # Multiple testing correction
        results_df = pd.DataFrame(results)
        reject, qval, _, _ = multipletests(results_df['pvalue'], 
                                            method='fdr_bh')
        results_df['padj'] = qval
        results_df['significant'] = reject
        
        return results_df.sort_values('pvalue')

# Usage:
# diff_analysis = DifferentialPathwayAnalysis(pathway_activity, clinical_df)
# age_results = diff_analysis.compare_by_group('age_group', method='anova')
# print(age_results[age_results['significant']])
```

**Implementation Checklist**:
- [ ] Create `differential_pathway_analysis.py`
- [ ] Implement `DifferentialPathwayAnalysis` class
- [ ] Load clinical data with stratification
- [ ] Run t-tests/ANOVA for: age groups, sex, disease stage
- [ ] Apply FDR correction
- [ ] Identify significant pathways per stratum

---

### Module 3: Hub Gene Identification

**File**: `hub_gene_identification.py` (NEW)

```python
"""
Hub Gene Identification within Pathways
- Find most important genes per pathway
- Combine: centrality + expression + disease association
"""

import pandas as pd
import numpy as np
import networkx as nx
from typing import Dict, List

class HubGeneIdentifier:
    """Identify hub genes per pathway"""
    
    def __init__(self, pathway_genes: Dict[str, List[str]],
                 gene_network: nx.Graph):
        """
        pathway_genes: pathway → gene list
        gene_network: NetworkX graph of gene interactions
        """
        self.pathway_genes = pathway_genes
        self.gene_network = gene_network
        self.hub_genes = {}
        
    def calculate_hub_score(self, pathway: str, 
                           expr_data: pd.DataFrame = None) -> pd.DataFrame:
        """
        Calculate hub score for each gene in pathway
        
        Score = degree centrality + betweenness + expression variance
        """
        
        genes = self.pathway_genes[pathway]
        scores = []
        
        for gene in genes:
            # 1. Network centrality
            if gene in self.gene_network:
                degree = self.gene_network.degree(gene)
                betweenness = nx.betweenness_centrality(
                    self.gene_network.subgraph([g for g in genes 
                                               if g in self.gene_network])
                ).get(gene, 0)
            else:
                degree, betweenness = 0, 0
            
            # 2. Expression variance (if available)
            expr_var = 0
            if expr_data is not None and gene in expr_data.index:
                expr_var = expr_data.loc[gene].var()
            
            # Combined score
            hub_score = degree + betweenness + expr_var
            
            scores.append({
                'gene': gene,
                'degree': degree,
                'betweenness': betweenness,
                'expr_variance': expr_var,
                'hub_score': hub_score
            })
        
        return pd.DataFrame(scores).sort_values('hub_score', ascending=False)

# Usage:
# hub_finder = HubGeneIdentifier(pathway_dict, gene_network_nx)
# hub_genes = hub_finder.calculate_hub_score('Glycolysis / Gluconeogenesis')
```

**Implementation Checklist**:
- [ ] Create `hub_gene_identification.py`
- [ ] Implement hub score calculation
- [ ] Integrate with existing gene networks from Tab 0
- [ ] Test: identify top 10 hub genes per pathway
- [ ] Validate against known pathway components

---

### Module 4: Visualization Functions

**File**: `pathway_visualizations.py` (NEW)

```python
"""
Visualization functions for pathway analysis
- Pathway activity heatmaps
- Violin plots by clinical strata
- Hub gene networks
"""

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

def plot_pathway_activity_heatmap(pathway_activity: pd.DataFrame,
                                  clinical_data: pd.DataFrame,
                                  group_by: str) -> go.Figure:
    """
    Heatmap of pathway activity across patient groups
    
    Columns: patient groups (sorted by group_by variable)
    Rows: pathways
    Color: pathway activity (z-score normalized)
    """
    
    # Sort samples by clinical variable
    sample_order = clinical_data.sort_values(group_by).index
    heatmap_data = pathway_activity[sample_order]
    
    # Z-score normalize
    heatmap_data = (heatmap_data - heatmap_data.mean(axis=1).values.reshape(-1,1)) / \
                   heatmap_data.std(axis=1).values.reshape(-1,1)
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(title="Activity (z-score)")
    ))
    
    fig.update_layout(
        title=f"Pathway Activity Heatmap (grouped by {group_by})",
        xaxis_title="Patient Samples",
        yaxis_title="Pathways",
        height=600,
        width=1200
    )
    
    return fig

def plot_pathway_violin(pathway_activity: pd.DataFrame,
                        clinical_data: pd.DataFrame,
                        pathway: str,
                        group_by: str) -> go.Figure:
    """
    Violin plot: pathway activity distribution by group
    """
    
    data = pd.DataFrame({
        'activity': pathway_activity.loc[pathway],
        'group': clinical_data[group_by]
    })
    
    fig = px.violin(data, x='group', y='activity',
                    title=f"{pathway} Activity by {group_by}",
                    box=True, points='outliers')
    
    return fig

# Usage:
# fig = plot_pathway_activity_heatmap(pathway_activity, clinical_df, 'age_group')
# fig.show()
```

**Implementation Checklist**:
- [ ] Create `pathway_visualizations.py`
- [ ] Implement heatmap visualization
- [ ] Implement violin plots for group comparison
- [ ] Test with real data
- [ ] Optimize for Gradio display

---

### Module 5: Integration with Gradio (Tab 4 Enhancement)

**File**: Modify `app_full.py` - Add to Tab 4 (Gene Network Simulation)

```python
# In the Tab 4 section (around line 1660), add:

with gr.Tab("🧬 基因网络仿真"):  # Tab 4
    
    # ... existing content ...
    
    # NEW: Pathway Activity Analysis subtab
    with gr.TabGroup() as gene_sim_tabs:
        
        # Original subtabs (keep existing)
        with gr.Tab("Network Statistics"):
            # ... existing code ...
            pass
        
        with gr.Tab("Network Visualization"):
            # ... existing code ...
            pass
        
        with gr.Tab("Result Data"):
            # ... existing code ...
            pass
        
        # NEW SUBTAB: Pathway Activity Analysis
        with gr.Tab("📊 Pathway Activity Analysis"):
            with gr.Row():
                with gr.Column(scale=1):
                    pathway_group_var = gr.Dropdown(
                        choices=["age_group", "sex", "disease_stage"],
                        value="age_group",
                        label="Group by"
                    )
                    pathway_method = gr.Dropdown(
                        choices=["GSVA", "Mean Expression"],
                        value="GSVA",
                        label="Scoring Method"
                    )
                    calc_pathway_btn = gr.Button("Calculate Pathway Activity")
                
                with gr.Column(scale=3):
                    pathway_heatmap = gr.Plot(label="Pathway Activity Heatmap")
            
            with gr.Row():
                pathway_table = gr.DataFrame(label="Differential Pathway Analysis")
            
            # Callback
            def calculate_pathway_activity(group_var, method):
                # Load data
                scorer = PathwayActivityScorer(pathway_dict)
                scorer.load_expression_data('data/TCGA-COAD/filtered_hiseq_data.csv')
                pathway_activity = scorer.score_gsva()
                
                # Differential analysis
                diff_analysis = DifferentialPathwayAnalysis(pathway_activity, clinical_df)
                diff_results = diff_analysis.compare_by_group(group_var)
                
                # Visualization
                heatmap_fig = plot_pathway_activity_heatmap(
                    pathway_activity, clinical_df, group_var
                )
                
                return heatmap_fig, diff_results.to_dataframe()
            
            calc_pathway_btn.click(
                calculate_pathway_activity,
                inputs=[pathway_group_var, pathway_method],
                outputs=[pathway_heatmap, pathway_table]
            )
```

**Implementation Checklist**:
- [ ] Create Tab 4 subtab for pathway analysis
- [ ] Load pathway data on app startup
- [ ] Connect to PathwayActivityScorer
- [ ] Add heatmap visualization
- [ ] Add results table
- [ ] Test UI responsiveness

---

## Part 3: Testing & Validation Plan

### Unit Tests

**File**: `tests/test_pathway_analysis.py` (NEW)

```python
import unittest
import pandas as pd
import numpy as np
from pathway_activity import PathwayActivityScorer
from differential_pathway_analysis import DifferentialPathwayAnalysis

class TestPathwayActivity(unittest.TestCase):
    
    def setUp(self):
        # Create minimal test data
        self.pathways = {
            'Pathway1': ['GENE1', 'GENE2', 'GENE3'],
            'Pathway2': ['GENE2', 'GENE3', 'GENE4']
        }
        
        self.expr_data = pd.DataFrame(
            np.random.randn(4, 20),
            index=['GENE1', 'GENE2', 'GENE3', 'GENE4'],
            columns=[f'Sample{i}' for i in range(20)]
        )
    
    def test_scorer_initialization(self):
        scorer = PathwayActivityScorer(self.pathways)
        self.assertEqual(len(scorer.pathway_genes), 2)
    
    def test_gsva_scoring(self):
        scorer = PathwayActivityScorer(self.pathways)
        scorer.gene_expr = self.expr_data
        activity = scorer.score_gsva()
        
        self.assertEqual(activity.shape, (2, 20))
        self.assertTrue(np.all(np.isfinite(activity.values)))

if __name__ == '__main__':
    unittest.main()
```

**Test Checklist**:
- [ ] Unit tests for PathwayActivityScorer
- [ ] Unit tests for DifferentialPathwayAnalysis
- [ ] Integration test: end-to-end pathway scoring
- [ ] Performance test: time for 347 pathways
- [ ] Validation test: compare with known pathway activity reference

### Manual Validation

- [ ] Load TCGA data and verify dimensions
- [ ] Score 347 pathways for 255 samples
- [ ] Calculate differential pathways by age/sex/stage
- [ ] Inspect top 5 hub genes in Glycolysis pathway (should include known glycolytic enzymes)
- [ ] Visualize results in Gradio interface
- [ ] Verify numerical outputs make biological sense

---

## Part 4: Implementation Timeline

### Week 1: Modules 1-2
- Days 1-2: Setup environment, implement PathwayActivityScorer
- Days 3-4: Test pathway scoring on small data subset
- Day 5: Implement DifferentialPathwayAnalysis

### Week 2: Modules 3-4
- Days 1-2: Implement HubGeneIdentifier
- Days 3-4: Create pathway_visualizations.py
- Day 5: Create unit tests

### Week 3: Module 5 + Testing
- Days 1-2: Integrate into Tab 4
- Days 3-4: Manual validation and testing
- Day 5: Documentation and bug fixes

### Week 4: Polish & Deploy
- Days 1-2: Performance optimization
- Days 3-4: User testing and documentation
- Day 5: Deployment and monitoring

---

## Part 5: Success Criteria

✅ Phase 1 is complete when:

1. **Functionality**
   - [ ] Can score 347 pathways across 255 TCGA samples in <1 minute
   - [ ] Differential analysis identifies significant pathways (FDR < 0.05)
   - [ ] Hub genes match biological expectations (>80% validation)

2. **Integration**
   - [ ] Tab 4 includes new Pathway Activity subtab
   - [ ] User can select grouping variable (age/sex/stage)
   - [ ] Results visualized in heatmap + data table

3. **Performance**
   - [ ] App loads in <5 seconds
   - [ ] Calculations complete in <30 seconds for user interaction

4. **Validation**
   - [ ] All unit tests pass (100%)
   - [ ] Pathway results validated against literature (>90% match for known pathways)
   - [ ] User documentation complete

---

## Part 6: Troubleshooting Guide

### Issue: "Gene not found in expression matrix"
**Solution**: Normalize gene symbol formats (uppercase, remove spaces)

### Issue: "All pathway scores are NaN"
**Solution**: Check gene list format matches expression index

### Issue: "Slow performance with 347 pathways"
**Solution**: 
- Implement caching of scorer results
- Pre-compute on app startup
- Use parallel processing for multiple pathways

### Issue: "Statistical test shows no significant pathways"
**Solution**:
- Check sample sizes per group
- Verify clinical variable has variation
- Try non-parametric tests (Wilcoxon)

---

## Part 7: Resources & References

### Key Papers
- Hanzelmann et al. (2013). "GSVA: Gene set variation analysis"
- Barbie et al. (2009). "Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1"
- Subramanian et al. (2005). "Gene set enrichment analysis"

### Documentation
- [gseapy documentation](http://gseapy.readthedocs.io/)
- [Plotly Python documentation](https://plotly.com/python/)
- [scikit-learn clustering](https://scikit-learn.org/stable/modules/clustering.html)

### Data Sources
- KEGG pathways: http://www.kegg.jp/
- TCGA-COAD: https://www.cancer.gov/tcga

---

**Document Status**: Ready for Implementation  
**Maintainer**: Bioinformatics Platform Team  
**Last Updated**: 2026-04-14

