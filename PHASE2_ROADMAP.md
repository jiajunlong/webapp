# Phase 2 Roadmap: Network Medicine & Disease Module Detection

**Status**: Ready for Implementation  
**Timeline**: 2-3 weeks  
**Estimated Effort**: 80-120 hours  
**Prerequisites**: Phase 1 Complete ✅

---

## Overview

Phase 2 bridges the **Cell/Tissue → Disease/Population** scales by implementing network medicine approaches. This phase connects expression-based pathway activity (Phase 1) to gene-disease associations through protein interaction networks and identifies disease modules.

---

## Phase 2 Architecture

```
Phase 1 Output (Pathway Activity)
    ↓
    ├─→ Disease Module Detection
    │   ├─ Map genes to PPI network (STRING/BioGRID)
    │   ├─ Community detection (Louvain/Leiden)
    │   └─ Module separation metrics
    │
    ├─→ WGCNA Co-expression Modules
    │   ├─ Cluster correlated genes
    │   ├─ Module-trait correlations
    │   └─ Identify key module genes
    │
    └─→ miRNA-Gene-Pathway Integration
        ├─ miRNA-target prediction
        ├─ Regulatory network inference
        └─ Map to pathways/modules
            ↓
        Phase 2 Output (Disease Modules + Regulations)
            ↓
        → Phase 3 (Parameter Propagation to SIS)
```

---

## Phase 2 Task Breakdown

### Task 1: Disease Module Detection (Primary)

**File**: `disease_module_detection.py` (NEW - ~400 lines)

**Components**:
1. `DiseaseNetworkBuilder` class
   - Load gene-disease associations from gene_disease.tsv
   - Download/cache PPI network (STRING or BioGRID)
   - Map genes to PPI network nodes
   - Build disease-specific subnetworks

2. `CommunityDetector` class
   - Implement Louvain community detection (via `python-louvain` or `leidenalg`)
   - Extract communities (disease modules)
   - Rank modules by internal cohesion/external separation

3. `ModuleSeparationMetrics` class
   - Compute module separation (Barabási framework)
   - Calculate Jaccard similarity between disease modules
   - Identify overlap/hierarchy relationships
   - Predict comorbidity based on module proximity

**Inputs**:
- `gene_disease.tsv` (existing) - gene-disease associations
- PPI network (STRING, ~4M interactions) - download on first use
- Gene expression (Phase 1 pathway scores)

**Outputs**:
- Disease modules (list of genes per disease)
- Module separation matrix (disease × disease)
- Module activity scores (per patient)

**Key Methods**:
```python
# Example usage:
builder = DiseaseNetworkBuilder("data/gene_disease.tsv")
builder.load_ppi_network("string")  # Downloads if needed

# Detect modules for a disease
modules = builder.detect_disease_modules("Colorectal Neoplasms")
separation = builder.compute_module_separation()

# Predict comorbidity
comorbidities = builder.predict_comorbidities(threshold=0.8)
```

**Literature Basis**:
- Menche et al. (2015) - "Uncovering disease-disease relationships through the incomplete interactome", Science
- Barabási et al. (2011) - "Network medicine: a network-based approach to human disease", Nature Rev Genetics

---

### Task 2: WGCNA Co-expression Modules (Secondary)

**File**: `wgcna_analysis.py` (NEW - ~350 lines)

**Components**:
1. `WGCNAAnalyzer` class
   - Compute soft power threshold for network construction
   - Build weighted correlation network
   - Identify co-expression modules (dynamic tree cutting)
   - Compute module eigengenes

2. `ModuleTraitCorrelation` class
   - Correlate module eigengenes with clinical traits
   - Statistical testing and FDR correction
   - Identify clinically relevant modules

3. `ModuleComparison` class
   - Compare WGCNA modules with disease modules
   - Compute overlap metrics
   - Identify common core genes

**Inputs**:
- Gene expression (14,520 × 255)
- Clinical traits (age, stage, gender, disease status)

**Outputs**:
- Module assignments (gene → module ID)
- Eigengene matrix (module × sample)
- Module-trait correlation matrix
- Hub genes per module

**Key Methods**:
```python
# Example usage:
analyzer = WGCNAAnalyzer(expr_data)
modules = analyzer.build_network(soft_power=8)
eigengenes = analyzer.compute_eigengenes()

trait_corr = WGCNAAnalyzer.correlate_with_traits(eigengenes, clinical_data)
hub_genes = analyzer.identify_hub_genes_per_module()
```

**Literature Basis**:
- Langfelder & Horvath (2008) - "WGCNA: an R package for weighted correlation network analysis"

---

### Task 3: miRNA-Gene-Pathway Integration (Tertiary)

**File**: `mirna_integration.py` (NEW - ~300 lines)

**Components**:
1. `miRNATargetPredictor` class
   - Correlate miRNA with gene expression
   - Identify predicted targets (correlation-based)
   - Validate against miRDB/TargetScan databases

2. `miRNARegulatoryNetwork` class
   - Build miRNA → gene regulatory network
   - Map regulatory edges to pathways
   - Identify key regulatory hub miRNAs

3. `RegulatoryModuleAnalysis` class
   - Integrate miRNA into disease modules
   - Identify regulatory modules (miRNA + genes + pathway)
   - Score regulatory importance

**Inputs**:
- miRNA expression (619 × 255)
- Gene expression (14,520 × 255)
- Pathway mappings (Phase 1)
- Gene-disease associations

**Outputs**:
- miRNA-target correlation network
- Regulatory hub miRNAs
- Regulatory modules
- miRNA-pathway associations

---

### Task 4: Integration into Gradio App

**File**: Modify `app_full.py` - Add to Tab 5 (Model Library)

**New Subtab Structure**:
```
Tab 5: 📚 模型库 (Model Library)
├── [existing content]
└── [NEW] 🔗 网络医学分析 (Network Medicine)
    ├── Disease Module Detection
    │   ├── Select disease
    │   ├── View module genes
    │   ├── Comorbidity predictions
    │   └── Module separation heatmap
    │
    ├── WGCNA Co-expression
    │   ├── View modules
    │   ├── Module-trait correlation heatmap
    │   ├── Hub genes per module
    │   └── Module overlap with disease modules
    │
    └── miRNA Regulation
        ├── Select miRNA
        ├── Target genes
        ├── Regulatory pathway mapping
        └── Regulatory module analysis
```

**UI Components**:
- Disease selector dropdown
- Module visualization (force-directed network graph)
- Separation matrix heatmap
- Hub gene tables
- Comorbidity prediction results

---

## Implementation Priorities

### Priority 1 (Weeks 1-1.5): Core Disease Module Detection
- [ ] Implement `DiseaseNetworkBuilder`
- [ ] Integrate STRING PPI network (auto-download)
- [ ] Implement community detection (Louvain)
- [ ] Compute module separation metrics
- [ ] Create basic Gradio interface

**Dependencies**: networkx, python-louvain, requests (for PPI download)

### Priority 2 (Week 2): WGCNA Integration
- [ ] Implement `WGCNAAnalyzer`
- [ ] Soft power threshold optimization
- [ ] Module-trait correlations
- [ ] Add to Gradio Tab 5

**Dependencies**: scikit-learn (for hierarchical clustering)

### Priority 3 (Week 2.5): miRNA Integration & Polish
- [ ] Implement `miRNATargetPredictor`
- [ ] Build regulatory networks
- [ ] Cross-validate with databases
- [ ] Finalize UI and documentation

**Dependencies**: pandas, scipy

---

## Expected Features on Completion

### Disease Module Detection Capabilities
- ✅ Map 4,322 genes to ~18,000 PPI edges from STRING
- ✅ Detect 50-200 disease modules per disease
- ✅ Compute module separation for all 2,503 diseases
- ✅ Predict comorbidities with 75-85% biological relevance
- ✅ Visualize as force-directed network graphs
- ✅ Compare with Phase 1 pathway activity

### WGCNA Capabilities
- ✅ Identify 30-60 co-expression modules from TCGA-COAD
- ✅ Correlate modules with clinical traits
- ✅ Identify trait-specific and core modules
- ✅ Extract hub genes (top 10-50 per module)
- ✅ Compare with disease modules (overlap analysis)

### miRNA Integration
- ✅ Link 619 miRNAs to 4,322 genes
- ✅ Build tissue-specific regulatory networks
- ✅ Validate predictions against 3 databases
- ✅ Identify hub regulatory miRNAs
- ✅ Map regulations to pathways

---

## Testing Strategy

### Unit Tests (25-30 tests)
- `test_disease_network_builder.py` (8 tests)
  - PPI network loading
  - Gene mapping
  - Subnetwork extraction
  - Disease module detection

- `test_wgcna_analyzer.py` (8 tests)
  - Soft power optimization
  - Network construction
  - Eigengene computation
  - Trait correlation

- `test_mirna_integration.py` (6 tests)
  - Target prediction
  - Regulatory network construction
  - Database validation

- `test_module_metrics.py` (4 tests)
  - Module separation computation
  - Comorbidity prediction
  - Overlap metrics

### Integration Tests (10 tests)
- End-to-end disease module detection
- WGCNA module validation
- Cross-scale integration with Phase 1
- Gradio UI functionality

### Validation Tests (5 tests)
- Known disease-disease associations (literature validation)
- Module genes match known disease pathways
- WGCNA modules correlate with clinical phenotypes

**Target**: 40-45 tests, 100% pass rate

---

## Data Requirements

### Input Data (Already Available)
- ✅ Gene expression: 14,520 genes × 255 samples
- ✅ Gene-disease associations: 4,322 genes × 2,503 diseases
- ✅ miRNA expression: 619 × 255 samples
- ✅ Clinical traits: age, stage, gender

### Required External Data (Auto-Download)
- PPI network (STRING v12): ~18,000 edges, 4,000 nodes (20 MB)
- miRNA target databases: miRDB, TargetScan (optional, reference only)

### Output Data Size
- Disease modules: ~10-20 MB (per disease)
- WGCNA modules: ~5 MB
- Module separation matrix: ~100 MB (2,503 × 2,503)
- Regulatory networks: ~50 MB

---

## Performance Targets

| Component | Expected Time | Hardware |
|-----------|---------------|----------|
| PPI network loading | <5 sec | SSD |
| Community detection | ~2 min | CPU |
| Module separation (all diseases) | ~5-10 min | CPU (parallel) |
| WGCNA (full dataset) | ~10-15 min | CPU (parallel) |
| miRNA target prediction | ~3-5 min | CPU |
| **Full pipeline** | **~30-35 min** | CPU (8 cores) |

---

## Deliverables

### Code (3 new modules)
- [ ] `disease_module_detection.py` - ~400 lines
- [ ] `wgcna_analysis.py` - ~350 lines
- [ ] `mirna_integration.py` - ~300 lines
- [ ] Update `app_full.py` - Add Tab 5 subtab

### Documentation
- [ ] PHASE2_IMPLEMENTATION_GUIDE.md
- [ ] PHASE2_TEST_REPORT.md
- [ ] API documentation for all 3 modules
- [ ] User guide for Gradio interface

### Tests
- [ ] 40-45 unit tests
- [ ] 10 integration tests
- [ ] Validation against literature

### Visualizations
- [ ] Disease module force-directed networks
- [ ] Module separation heatmaps
- [ ] WGCNA module-trait correlation heatmaps
- [ ] miRNA regulatory network diagrams

---

## Known Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|-----------|
| PPI network incomplete | Module detection may miss genes | Use multiple networks (BioGRID, Reactome) |
| Community detection is stochastic | Results may vary | Set random seed, run multiple iterations |
| WGCNA slow with 14k genes | Analysis may take 20+ min | Pre-filter genes, use approximate algorithms |
| miRNA targets hard to validate | Predictions may be noisy | Use correlation thresholds, validate subset |

---

## Phase 2 → Phase 3 Bridge

Phase 2 outputs feed directly into Phase 3 (Cross-Scale Parameter Propagation):

```
Phase 1 Output: Pathway activity scores
    ↓
Phase 2 Output: Disease modules, module activity
    ↓
Phase 3 Input: Map patient-specific pathway scores → 
               disease module predictions →
               SIS epidemic parameters (β, γ per module)
```

This enables truly multi-scale analysis where individual gene effects propagate through pathways to population-level disease dynamics.

---

## Success Criteria

Phase 2 is complete when:

1. **Functionality** ✅
   - [ ] Disease modules detected for 90%+ of diseases
   - [ ] Module separation computed for all 2,503 × 2,503 disease pairs
   - [ ] WGCNA identifies 40-60 modules with clinical significance
   - [ ] miRNA-gene links established for 80%+ of genes

2. **Validation** ✅
   - [ ] Disease modules validated against known pathways (80%+ match)
   - [ ] WGCNA modules correlate with clinical traits (p < 0.05)
   - [ ] Comorbidity predictions match literature (75%+ validation)

3. **Performance** ✅
   - [ ] Full pipeline runs in <1 hour on CPU
   - [ ] Interactive Gradio interface responds in <10 sec
   - [ ] Scales to 3000+ diseases

4. **Integration** ✅
   - [ ] Seamless Tab 5 integration in Gradio app
   - [ ] All Phase 1 + Phase 2 results accessible
   - [ ] No breaking changes to existing functionality

5. **Documentation** ✅
   - [ ] 100% code documentation
   - [ ] 40-45 unit tests passing
   - [ ] User guide and troubleshooting
   - [ ] API documentation

---

## Timeline Overview

```
Week 1:     Core disease module detection + PPI integration
Week 2:     WGCNA implementation + module comparison
Week 2.5:   miRNA integration + Gradio interface
Week 3:     Testing, validation, documentation + Phase 3 planning
```

**Target Completion**: End of Week 3 (2026-05-03)

---

**Status**: Ready for Implementation  
**Next Step**: Task #11 (Disease Module Detection) - Begin Week of 2026-04-22
