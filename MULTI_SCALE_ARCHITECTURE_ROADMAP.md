# Multi-Scale Bioinformatics Platform: Architecture Roadmap

**Date**: 2026-04-14  
**Status**: Research & Planning Phase  
**Objective**: Reorganize platform around scientifically-grounded multi-scale framework

---

## Executive Summary

The existing cross-scale gene & social network simulation platform contains substantial data across multiple biological scales but lacks systematic integration. This roadmap identifies:

1. **Current state**: Inventory of all data sources and existing analyses
2. **Scale coverage gaps**: Where data/methods are weak
3. **Literature-based methods**: Published algorithms for cross-scale analysis
4. **Implementation priorities**: Sequenced roadmap for platform enhancement

---

## Part 1: Current Multi-Scale Architecture

### Data Landscape by Scale

```
SCALE 1: GENE/PATHWAY LEVEL (MOLECULAR SCALE)
├─ Data: 4,322 genes, 2,503 diseases, 347 pathways
├─ Source: KEGG gene_disease.tsv + pathway files
├─ Status: ✅ Comprehensive knowledge base
└─ Analysis: Gene networks, pathway membership, drug associations

SCALE 2: PATIENT-LEVEL EXPRESSION (CELLULAR/TISSUE SCALE)  
├─ Data: 14,520 genes, 619 miRNAs, 14,520 CpG sites
├─ Samples: 255 TCGA-COAD tissue samples with clinical metadata
├─ Status: ✅ Rich omics data, partial clinical integration
└─ Analysis: MRNetB network inference, stratified by age/sex/stage

SCALE 3: POPULATION DYNAMICS (POPULATION SCALE)
├─ Data: Synthetic social networks (50-200 nodes, 2-10 communities)
├─ Status: ✅ SIS epidemic simulator implemented
└─ Analysis: Transmission dynamics, community effects

SCALE 4: CROSS-SCALE INTEGRATION (EMERGING)
├─ Data: Hub genes → TCGA networks → SIS parameters
├─ Status: ⚠️ Conceptual framework exists, limited implementation
└─ Analysis: Multi-Scale Linkage Analysis tab (Tab 6) - needs development
```

### Existing Computational Models

| Model | Scale | Algorithm | Implementation | Status |
|-------|-------|-----------|-----------------|--------|
| Gene Interaction Network | Molecular | Graph theory, degree centrality | ✅ Complete | Tab 0 |
| Gene Regulatory Network | Molecular | Directed graphs, regulatory cascade | ✅ Complete | Tab 0 |
| IS Coefficient | Molecular | Influence score f(coverage, overlap, density) | ✅ Complete | Tab 1 |
| MRNetB | Cellular | Mutual information network inference | ✅ Complete | Tab 4 |
| Network Flow Entropy (MFE) | Cellular | Information-theoretic entropy | 📋 Listed | Tab 5 |
| SIS Dynamics | Population | Epidemic transmission on networks | ✅ Complete | Tab 3 |

---

## Part 2: Identified Gaps & Opportunities

### Gap 1: Gene → Pathway Cross-Scale Analysis
**Current State**: 
- Gene-pathway membership is static lookup tables
- No quantitative "gene importance" within pathways

**Needed Methods**:
- **Gene Set Enrichment Analysis (GSEA)** - identify significantly enriched gene sets
- **Pathway activity scoring** - quantify pathway activation from expression data
  - ssGSEA (single sample GSEA)
  - GSVA (Gene Set Variation Analysis)
  - PLAGE (Pathway Level Analysis of Gene Expression)
- **Hub gene identification** - topological + expression-based importance

**Data Availability**: ✅ Have all needed data
- Gene symbols ✓
- Pathway memberships ✓
- Expression data ✓
- Clinical metadata ✓

---

### Gap 2: Pathway → Cell/Tissue Scale
**Current State**:
- TCGA patient stratification (age, sex, stage) exists
- No pathway-level activity analysis across strata

**Needed Methods**:
- **Pathway activity scoring across patients** - ssGSEA/GSVA on expression matrices
- **Differential pathway analysis** - pathway activity changes by clinical group
  - LIMMA: Linear models for expression analysis across groups
  - edgeR/DESeq2: For count-based pathway scoring
- **Pathway-methylation integration** - epigenetic regulation of pathways
- **miRNA-pathway targeting** - identify pathways regulated by differentially expressed miRNAs

**Data Availability**: ✅ Complete
- Expression matrices ✓
- miRNA expression ✓
- Methylation data ✓
- Clinical stratification ✓

---

### Gap 3: Cell/Tissue → Disease/Population Level
**Current State**:
- Patient networks inferred from TCGA (MRNetB)
- SIS model runs on synthetic social networks
- **No connection between tissue-level insights and population dynamics**

**Needed Methods**:
- **Network medicine approach** (Barabási et al.)
  - Disease modules: Identify connected gene subnetworks associated with disease
  - Degree-based prioritization: Hub genes in disease modules as drug targets
- **Patient stratification networks** - Network-based patient clustering
  - Use expression similarity to stratify patients into phenotype groups
  - Map to disease progression/outcome
- **Disease signature propagation** - From molecular networks to clinical phenotypes
  - Disease-specific biomarker networks
  - Prognostic signatures combining pathway activity + clinical variables

**Data Availability**: ✅ Available but underutilized
- TCGA networks from MRNetB ✓
- Clinical outcomes (disease stage) ✓
- Expression signatures ✓
- Missing: Explicit survival/outcome data (could use stage as proxy)

---

### Gap 4: Population-Level Disease Network
**Current State**:
- SIS model runs on social networks (not disease-related)
- No integration with disease-disease networks

**Needed Methods**:
- **Disease-disease network** (co-occurrence of genes/pathways)
  - Build from gene_disease.tsv: diseases sharing genes/pathways
  - Network medicine: Disease comorbidity networks
- **Transmission on disease networks** - SIS on realistic disease-comorbidity structure
  - Simulate disease "spread" through shared pathways
  - Population-level disease cascade models

**Data Availability**: ✅ Available
- Gene-disease associations ✓
- Pathway-disease associations ✓
- Pathway-pathway network ✓

---

### Gap 5: Visualization & Integration
**Current State**:
- Network visualizations (force-directed graphs)
- Statistical tables and bar charts
- No cross-scale flow visualization

**Needed Visualizations**:
- **Alluvial/Sankey diagrams**: Gene → Pathway → Disease → Population flow
- **Multi-layer network views**: Show same genes across scales
- **Heatmaps with hierarchical clustering**: Pathway activity across patient strata
- **Network comparison plots**: Side-by-side comparison of networks at different scales
- **Interactive multi-scale dashboard**: Brush/filter at one scale, see effects on others

---

## Part 3: Recommended Implementation Roadmap

### Phase 1: Pathway Activity Analysis (2-3 weeks)
**Goal**: Connect expression data to pathway-level biology

**Tasks**:
1. Implement ssGSEA pathway activity scoring
   - Input: Expression matrix + pathway gene lists
   - Output: Pathway activity scores per patient
   - Libraries: Recommended use `scanpy` (Python) or `immunedeconv` packages

2. Calculate differential pathway activity by clinical strata
   - Age groups: [0-50, 50-70, 70+]
   - Sex: [M, F]
   - Disease stage: [IA, IB, IIA, IIB, IIIA, IIIB, IIIC, IV]
   - Statistics: T-tests with multiple testing correction (FDR)

3. Integrate with Gene → Pathway analysis
   - Hub genes identified from expression + network topology
   - Genes with highest loading in pathway activity signature

**Output**: New Tab or Tab 4 enhancement
- Pathway activity heatmap by patient strata
- Top pathways per stratum
- Hub gene identification per pathway

**Data Used**:
- `filtered_hiseq_data.csv` (14,520 × 255)
- `pathway(最终版).tsv` gene memberships
- `filtered_clinical.csv` strata labels

---

### Phase 2: Network Medicine & Disease Modules (2-3 weeks)
**Goal**: Identify disease-relevant gene modules connecting TCGA to KEGG diseases

**Tasks**:
1. Build gene-disease network from KEGG
   - Nodes: Genes from 2,503 diseases
   - Edges: Gene co-occurrence in disease pathways
   - Weights: Number of shared pathways/diseases

2. Identify disease modules (densely connected subnetworks)
   - Algorithm: Community detection (Louvain or InfoMAP)
   - Filter: Modules enriched for disease genes + TCGA expression changes

3. Integrate TCGA patient networks with disease modules
   - Project TCGA MRNetB networks onto disease module structure
   - Identify which disease modules are dysregulated by patient strata

4. Patient stratification by pathway activity + network topology
   - Combine pathway scores + network centrality
   - K-means or hierarchical clustering of patients
   - Map to disease stage/outcomes

**Output**: New "Disease Module Analysis" Tab
- Disease module network visualization
- Module membership (genes + patient association)
- Module activity change across patient strata
- Biomarker genes per module

**Data Used**:
- `gene_disease.tsv` (all 2,503 diseases)
- TCGA MRNetB networks (Tab 4 output)
- Pathway-disease mappings

---

### Phase 3: Cross-Scale Parameter Propagation (2-3 weeks)
**Goal**: Link molecular insights through to population-scale disease dynamics

**Tasks**:
1. Extract disease-specific SIS parameters from molecular scale
   - Hub gene connectivity (from disease modules) → transmission rate β
   - Pathway redundancy/robustness → recovery rate γ  
   - Network clustering coefficient → community structure

2. Build disease-disease networks for SIS model
   - Nodes: Diseases from KEGG (or top disease modules)
   - Edges: Gene/pathway co-occurrence
   - Run SIS on disease-disease network instead of social network

3. Implement multi-disease population dynamics
   - Parameter cascade: Gene network → Disease modules → SIS parameters
   - Show "disease spread" through shared pathways

4. Clinical outcome integration
   - Disease stage as disease progression measure
   - Validate parameter predictions with TCGA outcomes

**Output**: Enhance Tab 6 "Multi-Scale Linkage Analysis"
- Cascade parameter flow: Gene → Module → Disease → Population
- Disease progression network simulation
- Prediction of disease co-occurrence / comorbidity

**Data Used**:
- Hub genes from disease modules (Phase 2)
- Network topology metrics
- `gene_disease.tsv` disease relationships
- SIS simulator framework

---

### Phase 4: Cross-Scale Visualization & Dashboard (2-3 weeks)
**Goal**: Create unified multi-scale visualization interface

**Tasks**:
1. Implement Sankey/Alluvial diagram generator
   - Input: Gene selection
   - Flow: Gene → Pathway → Disease → Population
   - Show pathway participation, disease association, population prevalence

2. Multi-layer network visualization
   - Layer 1: KEGG gene networks
   - Layer 2: TCGA expression networks  
   - Layer 3: Disease modules
   - Interactive navigation: Zoom/pan across scales

3. Interactive heatmap dashboard
   - Rows: Pathways or genes or diseases
   - Columns: Patient strata (age/sex/stage groups)
   - Values: Pathway activity / expression / centrality
   - Brush/filter interaction: Select in one layer, highlight in others

4. Comparison and export features
   - Side-by-side network comparison across scales
   - Export cross-scale analysis results (SVG/PDF)

**Output**: Enhanced multi-scale dashboard
- Sankey diagram of selected gene's involvement
- Interactive network explorer
- Stratified pathway activity heatmaps
- Cross-scale summary statistics

**Data Used**: All scales (aggregated from Phases 1-3)

---

## Part 4: Literature & Method References

### Key Methods to Implement

#### Pathway Activity Scoring
**Methods**:
- **ssGSEA** (Barbie et al., 2009)
  - Single-sample Gene Set Enrichment Analysis
  - Outputs: Continuous pathway activity score per sample
  - Use case: Score all TCGA samples for pathway activity

- **GSVA** (Hänzelmann et al., 2013)
  - Gene Set Variation Analysis
  - Non-parametric method for activity scoring
  - Handles small gene sets and batch effects

- **PLAGE** (Tomfohr et al., 2005)
  - Pathway Level Analysis of Gene Expression
  - Uses PCA to derive pathway activity scores

**Recommended**: Use GSVA for robustness and flexibility

#### Network Medicine
**Key Papers**:
- Barabási et al. (2011) "Network medicine: a network-based approach to human disease"
- Menche et al. (2015) "Disease networks: uncovering disease-disease relationships"
- Key concept: Disease modules as connected components of protein-protein interaction networks

**Methods**:
- Community detection: Louvain, InfoMAP, or spectral clustering
- Module disease association: Hypergeometric enrichment
- Hub gene identification: Degree + betweenness centrality + expression changes

#### Differential Pathway Analysis
**Tools/Methods**:
- **GSEA** (Subramanian et al., 2005)
  - Rank genes by t-statistic or fold-change
  - Test for enrichment in gene sets

- **LIMMA** (Ritchie et al., 2015)
  - Linear models for microarray/RNA-seq
  - Model clinical covariates (age, sex, stage)
  - Direct pathway activity comparison

#### Network Inference
**Already Implemented**:
- **MRNetB** (Meyer et al., 2007)
  - Mutual Information-based network inference
  - Used in Tab 4 - no changes needed

#### Multi-Scale Visualization
**Methods**:
- **Alluvial diagrams** - show flow across hierarchical levels
- **Multi-layer networks** - represent scale-specific relationships
- **Interactive brushing** - coordinated views across scales
- **Graph layouts**: Hierarchical layout for gene→pathway→disease→population

---

## Part 5: Implementation Checklist

### Priority 1 (High Impact, Feasible): Pathway Activity Analysis
- [ ] Implement GSVA pathway activity scoring function
- [ ] Apply to all TCGA samples
- [ ] Calculate differential pathway activity by clinical strata
- [ ] Create visualization: pathway activity heatmap
- [ ] Add to Tab 4 or create new Tab

### Priority 2 (High Value, Moderate Effort): Gene-Disease Network & Modules
- [ ] Build gene-disease co-occurrence network from KEGG
- [ ] Implement community detection
- [ ] Identify disease-associated gene modules
- [ ] Overlay TCGA networks onto modules
- [ ] Create visualization: module networks with patient stratification

### Priority 3 (Integration & Validation): Cross-Scale Parameter Propagation
- [ ] Develop parameter extraction: hub gene metrics → SIS parameters
- [ ] Build disease-disease network
- [ ] Implement disease-level SIS dynamics
- [ ] Validate with TCGA disease stage data
- [ ] Enhance Tab 6 with cascade visualization

### Priority 4 (Polish & Documentation): Unified Dashboard
- [ ] Implement Sankey/Alluvial visualizations
- [ ] Create interactive multi-layer network explorer
- [ ] Build cross-scale comparison tools
- [ ] Documentation and user guide
- [ ] Performance optimization

---

## Part 6: Technical Dependencies & Libraries

### Required Python Libraries

```python
# Pathway analysis
gsva              # Gene Set Variation Analysis (or equivalent)
gseapy            # GSEA implementation
scikit-learn      # For clustering, PCA, metrics

# Network analysis
networkx          # Network construction & algorithms
igraph            # Fast community detection
# (or) python-louvain  # Louvain algorithm

# Statistics
scipy             # Statistical tests
pandas            # Data manipulation
numpy             # Numerical operations

# Visualization
plotly            # Interactive plots (already used)
matplotlib/seaborn  # Static plots
plotly-orca       # Export to SVG/PNG
# (or) altair / vega-lite  # Alternative declarative viz

# Already in use
gradio            # Web framework
networkx          # Network visualization

# Data processing
scanpy            # Bioinformatics data structures
anndata           # HDF5-backed data containers
```

### Existing Framework

✅ Already has:
- Gradio web framework (Tab structure ready)
- Plotly for interactive visualization
- Data loading infrastructure
- Clinical metadata integration
- Expression data matrices loaded
- Network visualization components

---

## Part 7: Success Metrics

### Upon Completion:

1. **Pathway Activity Layer** ✅
   - Can score 347 pathways across 255 TCGA samples
   - Identify pathway differences by age/sex/stage
   - Integrate with Tab 4 for patient stratification

2. **Disease Module Layer** ✅
   - Identify 10-50 disease modules from 2,503 diseases
   - Map genes to modules
   - Show module activity in TCGA samples

3. **Cross-Scale Integration** ✅
   - Parameter propagation: Gene → Disease module → SIS model
   - Disease-disease network with realistic SIS parameters
   - Multi-disease cascade simulation

4. **User Experience** ✅
   - Unified multi-scale visualization dashboard
   - Interactive navigation across scales
   - Clear biological interpretation of cross-scale relationships

---

## Part 8: Next Steps

### Immediate (This Week)
1. ✅ Complete data inventory (DONE - see PROJECT_DATA_INVENTORY.md)
2. ✅ Identify architectural gaps (DONE - this document)
3. [ ] Research and validate literature methods for each gap
4. [ ] Create detailed technical specifications for Phase 1

### Short-term (Next 2 Weeks)
1. Prioritize which phases to tackle first
2. Set up development environment with required libraries
3. Implement Phase 1 (Pathway Activity Analysis)
4. Create initial pathway activity visualizations

### Medium-term (Weeks 3-8)
1. Implement Phases 2-3 (Network Medicine + Parameter Propagation)
2. Integrate with existing tabs
3. Begin cross-scale visualization (Phase 4)

### Long-term (Months 2+)
1. Complete Phase 4 (Dashboard)
2. Performance optimization
3. Testing with biological validation
4. User documentation and publication

---

**Document Author**: Bioinformatics Platform Review  
**Last Updated**: 2026-04-14  
**Status**: Ready for implementation planning

