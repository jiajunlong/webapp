# Next Steps: Implementation Roadmap

**Current Status:** Framework complete, documentation ready, roadmap defined  
**Target:** Implement 6 literature-validated cross-scale analysis methods  
**Timeline:** 10-22 days (staggered implementation)

---

## Implementation Priority Order

### Phase 1: Gene→Pathway Bridge (Week 1)

#### **Task 1: Implement ssGSEA** (1-2 days)
**What:** Single-Sample Gene Set Enrichment Analysis for pathway activity scoring  
**Why:** Creates the bridge from individual genes → pathway activity patterns  
**How:**
```python
# Install: pip install gseapy
import gseapy
# Input: 14,521 × 257 expression matrix + 359 pathways
# Output: 359 × 257 pathway activity matrix
# Add to: Tab 7 (cascade analysis)
```
**Reference:** Barbie et al., Nature (2009)  
**Success:** Pathway activity heatmap visualization in Tab 7

#### **Task 2: Differential Pathway Analysis** (1-2 days)
**What:** Test which pathways are differentially active across clinical groups  
**Why:** Connects pathway activity to disease stage, age, gender  
**How:**
```python
# Use ssGSEA output + clinical.tsv
# Tests: Wilcoxon rank-sum or Kruskal-Wallis across groups
# Output: significant pathways with p-values, effect sizes
# Visualization: heatmaps stratified by clinical groups
```
**Reference:** Standard TCGA analysis workflow  
**Success:** Identify disease-stage-specific pathways

---

### Phase 2: Tissue-Level Discovery (Week 2-3)

#### **Task 3: miRNA-Gene Correlation** (2-3 days)
**What:** Correlate 619 miRNAs with 14,521 genes to identify regulatory relationships  
**Why:** Adds regulatory layer; you already have the data!  
**How:**
```python
# Compute: Pearson correlations (619 × 14,521 matrix)
# Filter: |r| > 0.6, p < 0.01
# Map: correlations to pathways
# Visualization: bipartite network (miRNA ← gene ← pathway)
```
**Reference:** Integration of complementary data types  
**Success:** Identify top regulatory miRNA-gene-pathway triplets

#### **Task 4: WGCNA Co-expression Modules** (3-5 days)
**What:** Discover gene co-expression modules and correlate with clinical traits  
**Why:** Identifies tissue-level biological structures (the "systems")  
**How:**
```python
# Input: 14,521 × 257 expression + clinical traits
# Method: Weighted correlation network → hierarchical clustering
# Output: 10-30 gene modules + module-trait correlations
# Visualization: dendrograms, trait correlation heatmap
```
**Reference:** Langfelder & Horvath, BMC Bioinformatics (2008) — ~15,000 citations  
**Success:** Discover clinically-relevant gene modules

---

### Phase 3: Disease-Layer Integration (Week 4-5)

#### **Task 5: Disease Module Detection** (3-4 days)
**What:** Map disease-associated genes onto protein-protein interaction network  
**Why:** Identifies disease modules and inter-disease similarity  
**How:**
```python
# 1. Download PPI from STRING DB (confidence > 0.4)
# 2. Map 2,503 disease-gene associations onto PPI
# 3. Find connected components (disease modules)
# 4. Compute disease-disease similarity
# Visualization: disease module network
```
**Reference:** Barabási et al., Nature Reviews Genetics (2011); Menche et al., Science (2015)  
**Success:** Show disease-disease relationships via module overlap

---

### Phase 4: Novel Cross-Scale Method

#### **Task 6: SIS-as-Biomarker-Discovery** (2-3 days)
**What:** Use SIS epidemic dynamics on gene networks to identify biomarkers  
**Why:** Novel application: epidemiology meets systems biology  
**How:**
```python
# 1. Use MRNetB-inferred gene regulatory network
# 2. Seed with differentially expressed genes (high "infection")
# 3. Run SIS dynamics for 100-1000 steps
# 4. Genes with persistent "infection" = biomarkers
# 5. Validate: check if predicted biomarkers correlate with outcomes
```
**Reference:** Cowen et al., Nature Reviews Genetics (2017) — network propagation analogy  
**Success:** Identify novel biomarkers via SIS mechanism

---

## Implementation Sequence

```
Start → Task 1 → Task 2 → Task 3 + Task 4 (parallel) → Task 5 → Task 6 → End

        ↓ssGSEA    ↓Differential        ↓miRNA        ↓WGCNA
       Gene→      Pathway→         Gene→Pathway  Tissue→
       Pathway    Clinical         Regulatory    Modules
       
                                                ↓Disease Modules
                                                Tissue→Disease
                                                
                                                ↓SIS
                                                Cross-scale
```

---

## Quick Reference: What to Code

### Task 1: ssGSEA
```python
# File: tcga_coad_simulator.py or new ssgsea_module.py
import gseapy

def compute_ssgsea_pathways(expression_matrix, pathway_genes):
    """
    Input: genes × samples matrix, dict of pathways
    Output: pathways × samples activity matrix
    """
    # Use gseapy.ssgsea() implementation
    # Return DataFrame with pathway activity scores
```

### Task 2: Differential Analysis
```python
# File: cross_scale_engine.py
from scipy.stats import wilcoxon, kruskal

def differential_pathway_activity(pathway_activity, clinical_groups):
    """
    Input: pathway activity matrix, clinical groupings
    Output: p-values, effect sizes for each pathway
    """
    # Test each pathway across clinical groups
    # Correct for multiple testing (FDR)
    # Return sorted results
```

### Task 3: miRNA Correlation
```python
# File: new mirna_analysis.py
from scipy.stats import pearsonr

def compute_mirna_gene_correlation(mirna_expr, gene_expr):
    """
    Input: miRNA × samples, gene × samples
    Output: filtered correlation pairs
    """
    # Compute all correlations
    # Filter by |r| > 0.6, p < 0.01
    # Map to pathways
```

### Task 4: WGCNA
```python
# File: new wgcna_analysis.py
# Option 1: Pure Python using correlation + hierarchical clustering
# Option 2: Use gseapy or scikit-learn's hierarchical clustering

def discover_coexpression_modules(expression_matrix, clinical_traits):
    """
    Input: genes × samples, clinical data
    Output: gene module assignments, module-trait correlations
    """
    # Build weighted correlation network
    # Hierarchical clustering
    # Module eigengene calculation
```

---

## File Organization

Suggested organization for new code:

```
webapp/
├── app_full.py (main UI)
├── cross_scale_engine.py (orchestration)
├── model_library.py (model catalog)
├── ssgsea_module.py (Task 1)
├── mirna_analysis.py (Task 3)
├── wgcna_analysis.py (Task 4)
├── disease_modules.py (Task 5)
├── sis_biomarker_discovery.py (Task 6)
├── social_network_sim.py (existing)
├── tcga_coad_simulator.py (existing)
└── [data and documentation]
```

---

## Testing Strategy

For each task:
1. **Test with small dataset** (subset of genes/samples)
2. **Verify output shape** (expected dimensions)
3. **Check visualization** (plots render correctly)
4. **Add to Tab 7** incrementally
5. **Get human feedback** before next task

---

## Documentation to Update

After each implementation:
1. Update `PROJECT_STATUS.md` with completion
2. Add method documentation to code
3. Update Tab 7 descriptions
4. Keep `QUICK_REFERENCE.txt` current

---

## Dependencies to Install

```bash
# Core (already installed)
pip install gradio plotly networkx pandas numpy

# For new tasks
pip install gseapy scipy scikit-learn  # Task 1-3, 4, 5
# For WGCNA: either gseapy or native implementation
# For disease modules: STRING DB download (free)
```

---

## Success Criteria (Overall)

- [ ] ssGSEA produces 359×257 pathway activity matrix
- [ ] Differential analysis identifies stage-specific pathways
- [ ] miRNA correlations map to regulatory patterns
- [ ] WGCNA discovers 10-30 modules with clinical trait correlation
- [ ] Disease modules show inter-disease similarity
- [ ] SIS biomarkers predict clinical outcomes
- [ ] All visualizations in Tab 7
- [ ] Documentation updated
- [ ] No new compilation errors

---

## Estimated Timeline

| Task | Effort | Difficulty | Cumulative |
|------|--------|------------|-----------|
| 1. ssGSEA | 1-2d | Low | 1-2d |
| 2. Differential | 1-2d | Low | 2-4d |
| 3. miRNA | 2-3d | Medium | 4-7d |
| 4. WGCNA | 3-5d | Medium | 7-12d |
| 5. Disease modules | 3-4d | Medium | 10-16d |
| 6. SIS biomarkers | 2-3d | Medium | 12-19d |
| **Total** | **12-19d** | **Avg: Medium** | **12-19d** |

**Notes:**
- Tasks 3 & 4 can run in parallel (reduced to 10-16 total)
- Each task adds functionality; platform is usable after Task 2
- Tasks are independent; can be reordered if needed

---

## Start Here

1. Read `PROJECT_STATUS.md` — Detailed methods & references
2. Pick **Task 1 (ssGSEA)** — Lowest barrier to entry
3. Check gseapy documentation for `ssgsea()` function
4. Add method to `tcga_coad_simulator.py`
5. Wire into Tab 7 cascade analysis
6. Test with subset of data
7. Move to Task 2

---

**Good luck! The foundation is solid. Now it's about implementing well-established methods.** 🚀

