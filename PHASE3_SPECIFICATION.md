# Phase 3: Cross-Scale Parameter Propagation & SIS-as-Network-Biomarker-Discovery

**Status**: Implementation Roadmap  
**Date**: 2026-04-14  
**Objective**: Implement SIS dynamics on multi-scale biological networks for biomarker discovery  
**Timeline**: 2-3 weeks  
**Effort**: 40-60 hours

---

## Executive Summary

Phase 3 bridges tissue-level insights (pathways, modules, co-expression networks) to population-scale disease dynamics by:

1. **Extracting SIS parameters** from molecular-scale networks
2. **Running SIS epidemic model** on disease-relevant network structures
3. **Identifying persistent biomarkers** through network propagation
4. **Validating** predictions against clinical outcomes
5. **Integrating** into Gradio app (Tab 7: Multi-Scale Analysis)

This creates a novel application where epidemiological dynamics identify disease biomarkers through network propagation mechanisms.

---

## Part 1: Conceptual Framework

### The SIS-as-Biomarker-Discovery Model

**Biological Interpretation**:
- **Nodes** = Genes in a disease-relevant network
- **Edges** = Regulatory or co-expression relationships (from Phase 2)
- **"Infection"** = Dysregulation (high expression or abnormal activity)
- **Transmission** = Propagation of dysregulation through network
- **Recovery** = Normalization of gene expression

**Key Insight**: Genes that remain "infected" (dysregulated) despite recovery attempts are **robust biomarkers** - they have network-level importance that makes them central to disease progression.

### Multi-Scale Integration

```
Phase 1: Gene/Pathway Level
  ↓ [Pathway activity scores + hub genes]
  ↓
Phase 2: Tissue/Module Level  
  ↓ [Disease modules + co-expression networks]
  ↓
Phase 3: Population/Dynamics Level
  ↓ [SIS parameters extracted from modules]
  ↓
Population Dynamics Simulation
  ↓ [Network propagation identifies biomarkers]
  ↓
Clinical Validation
  ↓ [Biomarkers correlate with outcomes]
```

---

## Part 2: Technical Architecture

### 2.1 Core Components

#### Component 1: Parameter Extraction Module
**File**: `parameter_extraction.py`

**Purpose**: Extract SIS parameters (β, γ, initial infection) from molecular networks

**Key Methods**:

```python
class ParameterExtractor:
    """
    Extract SIS model parameters from disease modules and networks
    """
    
    def __init__(self, disease_modules, expression_data, clinical_data):
        """
        Initialize with Phase 2 outputs
        
        Parameters:
        -----------
        disease_modules : dict
            Module name → list of genes
        expression_data : pd.DataFrame
            Gene expression (genes × samples)
        clinical_data : pd.DataFrame
            Clinical traits (samples × traits)
        """
    
    def extract_transmission_rate(self, module: str) -> float:
        """
        Estimate β (transmission rate) from module structure
        
        Theory: Highly connected hubs transmit dysregulation more effectively
        
        Method:
        1. Compute avg degree of genes in module
        2. Scale by expression variance (high variance = more active)
        3. β = f(connectivity, expression_variability)
        
        Returns: β ∈ [0, 1]
        """
    
    def extract_recovery_rate(self, module: str) -> float:
        """
        Estimate γ (recovery rate) from module robustness
        
        Theory: Modules with high pathway redundancy recover faster
        
        Method:
        1. Count alternative regulatory paths
        2. Measure feedback loop density
        3. γ = 1 / (1 + pathway_redundancy)
        
        Returns: γ ∈ [0, 1]
        """
    
    def extract_initial_infection(self, module: str, disease: str) -> np.ndarray:
        """
        Estimate initial "infection" from differential expression
        
        Method:
        1. Compare healthy vs diseased samples
        2. Genes with high |log2FC| and low p-value = high initial infection
        3. Scale to [0, 1]
        
        Returns: I ∈ [0, 1]^n where n = module size
        """
    
    def extract_network_structure(self, module: str) -> np.ndarray:
        """
        Build adjacency matrix from co-expression or PPI data
        
        Method:
        1. Get genes in module
        2. Extract edges from Phase 2 co-expression / PPI networks
        3. Weight edges by interaction strength
        
        Returns: Adjacency matrix (n × n)
        """
```

**Outputs**:
- Per-module parameters: (β, γ, I₀, A)
- Per-disease parameters (aggregate across modules)
- Parameter distributions for sensitivity analysis

---

#### Component 2: Network Propagation Engine
**File**: `sis_network_propagation.py`

**Purpose**: Run SIS dynamics and track biomarker persistence

**Key Methods**:

```python
class SISNetworkPropagation:
    """
    Run SIS dynamics on disease networks for biomarker discovery
    """
    
    def __init__(self, adjacency_matrix, parameters):
        """
        Initialize SIS model on network
        
        Parameters:
        -----------
        adjacency_matrix : np.ndarray
            Network structure (n × n)
        parameters : dict
            {'beta': float, 'gamma': float, 'initial_infection': np.ndarray}
        """
    
    def run_dynamics(self, n_steps: int = 1000, n_runs: int = 100) -> dict:
        """
        Run SIS dynamics for multiple timesteps and runs
        
        Method:
        1. Initialize with I₀
        2. At each step: For each infected node i
           - With prob β*I_neighbors: node j becomes infected
           - With prob γ: node i recovers
        3. Track persistence of infections
        
        Parameters:
        -----------
        n_steps : int
            Number of simulation steps
        n_runs : int
            Number of stochastic runs
        
        Returns:
        --------
        dict : {
            'infection_history': np.ndarray (n_runs × n_steps × n_nodes),
            'persistence_scores': np.ndarray (n_nodes,),  # % of steps infected
            'biomarkers': list,  # genes with high persistence
            'network_metrics': dict  # centrality, clustering, etc.
        }
        """
    
    def compute_persistence_score(self, node: int, infection_history: np.ndarray) -> float:
        """
        Compute persistence score for a single gene
        
        Method:
        1. Sum timesteps where node is infected across all runs
        2. Normalize by (n_runs × n_steps)
        3. Score = persistence + betweenness centrality weighted
        
        Returns: score ∈ [0, 1]
        """
    
    def identify_biomarkers(self, persistence_threshold: float = 0.5) -> pd.DataFrame:
        """
        Identify genes as biomarkers based on persistence
        
        Logic: Genes that remain "infected" despite recovery mechanisms
               are robust biomarkers
        
        Returns: DataFrame with genes, persistence scores, network metrics
        """
```

**Key Algorithm**:

```
SIS Dynamics on Network:
  For each run r = 1 to n_runs:
    I = I₀  # Initial infection vector
    For each timestep t = 1 to n_steps:
      I_new = I.copy()
      For each infected node i:
        # Transmission: can infect neighbors
        For each neighbor j of i:
          if rand() < β * w_ij * I[i]:
            I_new[j] = 1
        # Recovery: can recover
        if rand() < γ:
          I_new[i] = 0
      I = I_new
      store(I)  # Record state
  
  # Compute persistence
  persistence[i] = (# timesteps i was infected) / (n_runs × n_steps)
```

---

#### Component 3: Validation & Clinical Integration
**File**: `biomarker_validation.py`

**Purpose**: Validate predicted biomarkers against clinical outcomes

**Key Methods**:

```python
class BiomarkerValidator:
    """
    Validate SIS-predicted biomarkers against clinical data
    """
    
    def __init__(self, predicted_biomarkers, expression_data, clinical_data):
        """
        Initialize validator
        
        Parameters:
        -----------
        predicted_biomarkers : pd.DataFrame
            Genes with persistence scores from Phase 3
        expression_data : pd.DataFrame
            Gene expression (genes × samples)
        clinical_data : pd.DataFrame
            Clinical outcomes (samples × traits)
        """
    
    def validate_correlation_with_outcome(self, outcome_variable: str) -> pd.DataFrame:
        """
        Test if predicted biomarkers correlate with clinical outcome
        
        Method:
        1. For each predicted biomarker gene
        2. Compute correlation with outcome (stage, survival, etc.)
        3. Test significance (p-value)
        4. Compare with random genes (null distribution)
        
        Returns:
        --------
        pd.DataFrame : Validation metrics per biomarker
        """
    
    def compute_biomarker_signature_score(self, genes: List[str]) -> np.ndarray:
        """
        Compute multi-gene biomarker signature score
        
        Method:
        1. Average expression of predicted biomarkers
        2. Scale to [0, 1]
        3. Test association with clinical outcomes
        
        Returns: Signature score per sample
        """
    
    def compare_with_literature(self, genes: List[str]) -> pd.DataFrame:
        """
        Compare predicted biomarkers with known disease signatures
        
        Method:
        1. Query biomarker databases
        2. Check overlap with literature
        3. Report literature validation rate
        
        Returns: Literature validation metrics
        """
```

---

#### Component 4: Gradio Integration
**File**: Update `gradio_phase3_integration.py`

**UI Structure**:

```python
def create_sis_biomarker_tab():
    """
    Creates Phase 3 tab for SIS-as-Biomarker-Discovery
    
    Interface:
    -----------
    [Disease Selection] [Module Selection] [Network Type]
         ↓                    ↓                    ↓
    [Disease Module]  [Co-expression]  [PPI Network]
         ↓
    [Parameter Extraction]
         ↓
    [Set SIS Parameters: β, γ, I₀]
         ↓
    [Run Propagation] → Progress: 0-100%
         ↓
    [Results Tabs]:
    ├─ Infection Dynamics (heatmap over time)
    ├─ Persistence Scores (biomarker ranking)
    ├─ Network Visualization (highlight persistent genes)
    ├─ Clinical Validation (biomarker-outcome correlation)
    └─ Comparison with Literature
    """
```

---

### 2.2 Data Flow

```
Phase 1 Output
  ├─ Pathway activity scores
  ├─ Hub genes per pathway
  └─ Differential pathway analysis
    ↓
Phase 2 Output
  ├─ Disease modules
  ├─ Co-expression network (WGCNA)
  ├─ miRNA regulatory network
  └─ Module separation metrics
    ↓
Phase 3 Input Processing
  ├─ Extract network structure from Phase 2 modules
  ├─ Extract initial condition from Phase 1 hub genes
  └─ Compute parameters from network metrics
    ↓
SIS Dynamics
  ├─ Run propagation 100 times
  ├─ Track infection persistence
  └─ Compute biomarker scores
    ↓
Phase 3 Output
  ├─ Biomarker rankings
  ├─ Persistence profiles
  ├─ Clinical validation
  └─ Integrated biomarker signatures
    ↓
Gradio Tab 7: Visualization & Interaction
```

---

## Part 3: Implementation Details

### 3.1 File Organization

```
webapp/
├── parameter_extraction.py     (NEW, 200-250 lines)
├── sis_network_propagation.py  (NEW, 250-300 lines)
├── biomarker_validation.py     (NEW, 150-200 lines)
├── gradio_phase3_integration.py (NEW, 300-350 lines)
├── test_phase3_*.py            (NEW, 3 test files, 200-250 lines each)
├── app_full.py                 (MODIFY - Add Phase 3 tab)
└── PHASE3_*.md                 (NEW, documentation)
```

### 3.2 Dependencies

**Required Libraries**:
```python
# Network analysis (already available)
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# Visualization (already available)
import plotly.graph_objects as go
import plotly.express as px

# New dependencies
import seaborn as sns  # Heatmaps for dynamics
import matplotlib.pyplot as plt  # Additional visualization
```

**No additional pip installs needed** - leveraging existing libraries!

---

## Part 4: Algorithm Details

### 4.1 Parameter Extraction from Disease Modules

#### β (Transmission Rate)

**Principle**: Higher connectivity and expression variability → higher transmission

**Formula**:
```
β = min(1, α · connectivity · expression_variability)

where:
  connectivity = avg_degree / max_degree
  expression_variability = std(expression) / mean(expression)
  α ∈ [0.3, 0.7] (calibration parameter)
```

**Interpretation**: 
- β ≈ 0.1-0.3 for isolated modules (low transmission)
- β ≈ 0.5-0.8 for hub-rich modules (high transmission)

---

#### γ (Recovery Rate)

**Principle**: Pathway redundancy limits recovery (high redundancy = slow recovery)

**Formula**:
```
γ = 1 / (1 + redundancy_factor)

where:
  redundancy_factor = n_alternative_paths / n_genes
  high redundancy (e.g., 5 alternative paths per gene) → low γ ≈ 0.1
  low redundancy (e.g., 1 alternative path per gene) → high γ ≈ 0.5
```

**Interpretation**:
- γ ≈ 0.1-0.2 for robust modules with high redundancy
- γ ≈ 0.5-0.8 for fragile modules with low redundancy

---

#### I₀ (Initial Infection)

**Principle**: Highly dysregulated genes start "infected"

**Formula**:
```
For each gene g in disease_samples vs control_samples:
  log2FC[g] = mean(expr_disease[g]) - mean(expr_control[g])
  p-value[g] = t-test(expr_disease[g], expr_control[g])
  
  I₀[g] = {
    max(0, min(1, |log2FC[g]| / 2))   if p-value < 0.05
    0                                  otherwise
  }
```

**Interpretation**:
- I₀[g] ≈ 0 for stable genes
- I₀[g] ≈ 0.8-1.0 for highly dysregulated hub genes

---

### 4.2 SIS Dynamics Algorithm (Pseudocode)

```python
def sis_dynamics_stochastic(A, beta, gamma, I0, n_steps=1000, n_runs=100):
    """
    Stochastic SIS dynamics on network
    
    Parameters:
    -----------
    A : np.ndarray
        Adjacency matrix (weighted)
    beta : float
        Transmission rate
    gamma : float
        Recovery rate
    I0 : np.ndarray
        Initial infection vector
    n_steps : int
        Number of timesteps per run
    n_runs : int
        Number of stochastic realizations
    """
    
    n = A.shape[0]  # Number of nodes
    infection_history = np.zeros((n_runs, n_steps, n))
    
    for run in range(n_runs):
        I = I0.copy()  # Current infection state
        
        for t in range(n_steps):
            infection_history[run, t, :] = I
            
            # Transmission step: infected nodes infect neighbors
            I_new = I.copy()
            for i in range(n):
                if I[i] > 0:  # i is infected
                    for j in range(n):
                        if A[i, j] > 0 and I[j] == 0:  # i-j edge exists, j is susceptible
                            # Transmission probability
                            p_transmit = beta * A[i, j] * I[i]
                            if np.random.rand() < p_transmit:
                                I_new[j] = 1
            
            # Recovery step: infected nodes recover
            for i in range(n):
                if I[i] > 0:  # i is infected
                    if np.random.rand() < gamma:
                        I_new[i] = 0  # Recovery
            
            I = I_new
    
    return infection_history
```

---

### 4.3 Biomarker Identification

**Persistence Score**:
```
persistence[i] = (# timesteps i is infected across all runs) / (n_runs × n_steps)
                = (sum over all runs and timesteps of I[run, t, i]) / (n_runs × n_steps)
```

**Weighted Biomarker Score** (combining persistence + centrality):
```
biomarker_score[i] = 0.6 × persistence[i] + 0.4 × centrality[i]

where:
  centrality[i] = normalized degree or betweenness centrality
```

**Biomarker Selection**:
```
Biomarkers = genes with biomarker_score > percentile_75
           ∪ genes with high persistence that are hubs
           ∪ genes with early infection times
```

---

## Part 5: Expected Results

### 5.1 Per-Disease Module

For a disease module with 20-50 genes:

- **Top biomarkers**: 3-5 genes with persistence > 0.7
- **Runner-up biomarkers**: 5-10 genes with persistence > 0.5
- **Robust signatures**: Multi-gene signatures with clinical correlation

### 5.2 Validation Metrics

**Literature Validation**:
- % of predicted biomarkers in disease literature: Target > 60%
- Comparison with known drug targets: Overlap > 40%

**Clinical Validation** (on TCGA-COAD):
- Biomarker expression correlation with disease stage: |r| > 0.3, p < 0.05
- Survival analysis (if available): HR > 1.5, p < 0.05

---

## Part 6: Implementation Schedule

### Week 1: Core Infrastructure
- [ ] Implement `parameter_extraction.py` (8 hours)
- [ ] Implement `sis_network_propagation.py` (10 hours)
- [ ] Unit tests: parameter extraction (6 tests, 4 hours)
- [ ] Unit tests: SIS dynamics (5 tests, 4 hours)

### Week 2: Validation & Integration
- [ ] Implement `biomarker_validation.py` (6 hours)
- [ ] Implement `gradio_phase3_integration.py` (10 hours)
- [ ] Integration with app_full.py (4 hours)
- [ ] Unit tests: validation + integration (5 tests, 4 hours)

### Week 3: Polish & Documentation
- [ ] Sensitivity analysis (parameter sweeps) (6 hours)
- [ ] Create comprehensive test suite (10 tests, 8 hours)
- [ ] Documentation & examples (8 hours)
- [ ] Performance optimization (4 hours)

**Total**: 40-50 hours of implementation

---

## Part 7: Testing Strategy

### Unit Tests (20-25 tests)

**Test Classes**:
1. `TestParameterExtraction` (6 tests)
   - β extraction from different module types
   - γ extraction with varying redundancy
   - I₀ extraction from expression data
   - Edge cases (empty modules, single genes)

2. `TestSISNetworkPropagation` (6 tests)
   - Basic SIS dynamics on simple networks
   - Parameter sweeps (β, γ variations)
   - Persistence score computation
   - Biomarker identification consistency

3. `TestBiomarkerValidation` (5 tests)
   - Correlation computation
   - Clinical outcome association
   - Literature validation
   - Signature score aggregation

4. `TestPhase3Integration` (4 tests)
   - Gradio tab creation
   - Data flow integration
   - Phase 1-2-3 data compatibility
   - End-to-end workflow

### Integration Tests (5 tests)

1. Full pipeline: Phase 2 modules → Phase 3 biomarkers
2. Multiple diseases: Biomarker comparison across diseases
3. Sensitivity: Parameter sweep effects on results
4. Clinical validation: Biomarker-outcome correlation on TCGA
5. Cross-scale validation: Biomarkers match Phase 1 hub genes

### Validation Tests (3 tests)

1. **Known Gene Test**: Compare with TCGA-COAD known biomarkers
2. **Literature Overlap**: Biomarkers match disease literature
3. **Network Structure**: Persistent genes are high-degree hubs

---

## Part 8: Success Criteria

Phase 3 is complete when:

1. **Functionality** ✅
   - [ ] Extract SIS parameters for 10+ disease modules
   - [ ] Run SIS dynamics: 100+ runs × 1000 steps in < 5 minutes
   - [ ] Identify 20-50 biomarkers per disease
   - [ ] Clinical validation: > 50% biomarkers correlate with outcome

2. **Validation** ✅
   - [ ] Predicted biomarkers overlap with TCGA literature: > 60%
   - [ ] Biomarker-outcome correlation: p < 0.05 for > 70% of markers
   - [ ] Network structure validation: persistent genes are hubs (r > 0.7)

3. **Performance** ✅
   - [ ] Parameter extraction: < 1 sec per module
   - [ ] SIS propagation: < 5 min for 100 runs
   - [ ] Full pipeline: < 10 min for disease module
   - [ ] Gradio interface: < 3 sec response time

4. **Integration** ✅
   - [ ] Seamless Phase 1-2-3 data flow
   - [ ] Tab 7 displays all results interactively
   - [ ] Export biomarker predictions (CSV, JSON)
   - [ ] No breaking changes to Phases 1-2

5. **Documentation** ✅
   - [ ] 100% code documentation
   - [ ] 20-25 unit tests passing
   - [ ] Integration tests passing
   - [ ] User guide with examples
   - [ ] Biological interpretation guide

---

## Part 9: Known Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Parameter extraction may be suboptimal | Biomarkers unreliable | Run sensitivity analysis, calibrate against known biomarkers |
| SIS stochasticity adds variance | Results inconsistent | Run 100+ replicates, report confidence intervals |
| Biomarker validation limited by TCGA sample size | Low statistical power | Use permutation testing, cross-validation |
| Computational cost high for all diseases | Runtime > 1 hour | Cache results, parallelize across modules |
| Literature overlap limited by database coverage | Validation incomplete | Use multiple databases, manual verification |

---

## Part 10: Phase 3 → Future Phases

### Post-Phase 3 Enhancements

1. **Advanced Parameter Learning** (Phase 3.5)
   - Use machine learning to optimize β, γ from clinical outcomes
   - Personalized parameters per patient subgroup

2. **Multi-Network Integration** (Phase 4)
   - Run SIS on multiple networks simultaneously (PPI + co-expression)
   - Aggregate biomarker predictions

3. **Temporal Dynamics** (Phase 4.5)
   - Time-varying parameters reflecting disease progression
   - Disease stage-specific SIS models

4. **Drug Target Discovery** (Phase 5)
   - Map biomarkers to approved drugs
   - Predict drug efficacy based on network structure
   - Suggest combination therapies

---

## Success Metrics on Completion

### Biomarker Discovery Performance

- **Per-disease modules**: 20-30 robust biomarkers identified
- **Literature match**: 60-75% biomarkers in disease literature  
- **Clinical correlation**: > 70% biomarkers correlate with stage/outcome
- **Novel predictions**: 20-40% of biomarkers are novel candidates

### Cross-Scale Integration

- **Phase 1 Match**: Biomarkers overlap with Hub genes (r > 0.6)
- **Phase 2 Match**: Biomarkers located in disease module cores
- **Network Robustness**: Persistent genes have avg degree > 2× random

### User Experience

- **Interactive Exploration**: Users can adjust β, γ, visualize dynamics
- **Clinical Actionability**: Biomarkers translate to treatment targets
- **Reproducibility**: Results consistent across multiple runs

---

## Deliverables

### Code (4 new modules)
- [ ] `parameter_extraction.py` - ~250 lines
- [ ] `sis_network_propagation.py` - ~300 lines
- [ ] `biomarker_validation.py` - ~200 lines
- [ ] `gradio_phase3_integration.py` - ~350 lines
- [ ] Update `app_full.py` - Add Phase 3 tab (Tab 7)

### Documentation
- [ ] `PHASE3_SPECIFICATION.md` (this document)
- [ ] `PHASE3_IMPLEMENTATION_GUIDE.md` - Technical details
- [ ] `PHASE3_USER_GUIDE.md` - How to use Phase 3
- [ ] `PHASE3_TEST_REPORT.md` - Test results

### Tests
- [ ] 20-25 unit tests (100% pass rate)
- [ ] 5 integration tests
- [ ] 3 validation tests
- [ ] Comprehensive test coverage

### Results & Visualizations
- [ ] Biomarker rankings per disease
- [ ] Persistence dynamics heatmaps
- [ ] Network visualization with biomarker highlighting
- [ ] Clinical validation plots
- [ ] Cross-scale comparison matrices

---

**Status**: Ready for Implementation  
**Next Step**: Begin Week of 2026-05-01  
**Estimated Completion**: 2026-05-14

---

## Document History

- **2026-04-14**: Initial specification created
- **Status**: Pending implementation
- **Author**: Bioinformatics Platform Development
