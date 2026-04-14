# Comprehensive Project Data & Functionality Inventory

**Project**: Cross-Scale Gene & Social Network Simulation Platform  
**Framework**: Gradio Web Application  
**Analysis Date**: 2025-12-20

---

## 📊 PART 1: DATA SOURCES & STRUCTURE

### A. DATA DIRECTORY FILES

#### 1. **gene_disease.tsv** (8,498 rows)
- **Biological Entity**: Gene-Disease-Pathway-Drug relationships
- **Format**: Tab-separated values
- **Scope**: All 2,503 diseases in knowledge base
- **Key Columns**:
  - `disease_id` - KEGG disease ID
  - `disease_name` - Human-readable disease name
  - `disease_category` - Disease classification
  - `disease_pathway` - Associated KEGG pathways
  - `disease_gene` - KEGG gene IDs and names for that disease
  - `disease_drug` - Treatment drugs (DR:codes)
  - `gene_ENSG` - Ensembl gene ID
  - `gene_pathway` - Pathways the gene participates in
  - `gene_disease` - Diseases the gene is involved in
  - `gene_symbol` - Gene symbol (may include comma-separated aliases like "PTEN, PTEN1, PTENbeta")

**Biological Scale**: **MOLECULAR SCALE** - Gene/Pathway level
**Entity Types Covered**: 
- Genes (4,322 unique)
- Diseases (2,503 unique)
- Pathways (347 unique)
- Drugs (implicit, in disease_drug field)

---

#### 2. **pathway(最终版).tsv** (360 rows)
- **Biological Entity**: Metabolic and signaling pathways
- **Format**: Tab-separated values
- **Key Columns**:
  - `Pathway_ID` - KEGG pathway ID (e.g., "hsa01100")
  - `Pathway_Name` - Full pathway name (e.g., "Glycolysis / Gluconeogenesis - Homo sapiens")
  - `Class` - Pathway classification (mostly NA)
  - `Gene` - Comma/semicolon-separated gene list (up to ~30 genes per pathway)
  - `Drug` - Associated drug IDs (semicolon-separated)
  - `Related_Pathway` - Related pathway names (semicolon-separated)

**Biological Scale**: **MOLECULAR SCALE** - Pathway/Network level
**Entity Types Covered**:
- Pathways (347 total, KEGG Homo sapiens)
- Gene-pathway relationships
- Pathway-pathway relationships

**Example Pathways**:
- Glycolysis / Gluconeogenesis
- Citrate cycle (TCA cycle)
- Pentose phosphate pathway
- Oxidative phosphorylation
- Signal transduction pathways

---

#### 3. **pathway(基因名映射版).tsv** (Same as above, alternate version)
- Gene names already mapped
- Used as primary source in data_loader.py

---

#### 4. **Related Pathway.txt** (1,585 pathway pairs)
- **Biological Entity**: Pathway-pathway network relationships
- **Format**: Two-column tab-separated
- **Columns**:
  - `Pathway` - Source pathway ID
  - `Related Pathway` - Target pathway ID
- **Network Type**: Undirected pathway relationships
- **Example**: hsa00010 connects to hsa00020, hsa00030, hsa00500 (metabolic cascade)

**Biological Scale**: **MOLECULAR SCALE** - Metanetwork level

---

#### 5. **merged_output.tsv** (129,521 rows)
- **Biological Entity**: Gene-Gene relationships (network edges)
- **Format**: Tab-separated values
- **Columns**:
  - `source` - First gene/disease name
  - `target` - Second gene/disease name  
  - `weight` - Edge weight (interaction strength)
  - `disease` - Associated disease context
- **Purpose**: Pre-computed gene interaction network for visualization

**Biological Scale**: **MOLECULAR SCALE** - Gene network interactions

---

#### 6. **preprocessed_summary.txt**
- Summary statistics file (auto-generated)
- Contains: Disease count (2503), Pathway count (347), Gene count (4322)
- Lists top 20 diseases and pathways

---

### B. TCGA-COAD DIRECTORY FILES

#### 1. **clinical.tsv** (462 rows)
- **Biological Entity**: Patient clinical metadata
- **Format**: Tab-separated values
- **Subject Count**: 461 TCGA-COAD patients
- **Key Columns**:
  - `case_submitter_id` - Patient ID (e.g., "TCGA-F4-6854")
  - `age_at_index` - Age at diagnosis (numeric)
  - `gender` - Biological sex
  - `race` - Ethnic/racial category
  - `ajcc_pathologic_stage` - Tumor stage (e.g., "Stage IIA", "Stage IIIB")

**Biological Scale**: **CLINICAL/POPULATION SCALE** - Patient phenotypes
**Entity Types**: Patients with tumor type: Colon Adenocarcinoma (COAD)

---

#### 2. **filtered_hiseq_data.csv** (14,521 rows × ~255 columns)
- **Biological Entity**: Gene expression measurements
- **Format**: CSV matrix (genes × samples)
- **Data Type**: RNA-seq from Illumina HiSeq
- **Rows**: Genes (14,520 genes + 1 header)
- **Columns**: Patient samples (255 TCGA-COAD samples)
- **Structure**: 
  - First column: "gene" (header) followed by sample barcodes
  - Row 2+: Sample TCGA IDs (e.g., "TCGA-D5-6920-01")
  - Values: Normalized gene expression counts

**Biological Scale**: **CELLULAR/MOLECULAR SCALE** - Gene expression level
**Entity Types**: 
- Genes (14,520 measurable transcripts)
- Samples (255 patient tissue samples)
- Expression values (continuous)

---

#### 3. **filtered_methylation_data.csv** (14,521 rows × ~255 columns)
- **Biological Entity**: DNA methylation measurements
- **Format**: CSV matrix (CpG sites × samples)
- **Data Type**: Illumina 450K methylation arrays
- **Rows**: CpG/methylation sites (14,520 sites)
- **Columns**: Patient samples (255 samples, matched to hiseq)
- **Values**: Beta values (0-1, methylation ratio)

**Biological Scale**: **CELLULAR/MOLECULAR SCALE** - Epigenetic regulation
**Entity Types**:
- Methylation sites (14,520 CpG positions)
- Samples (255 matched to expression data)
- Methylation beta values

---

#### 4. **filtered_miRNA_with_names.csv** (620 rows × ~255 columns)
- **Biological Entity**: microRNA expression
- **Format**: CSV matrix (miRNAs × samples)
- **Data Type**: Small RNA sequencing
- **Rows**: miRNAs (619 mature miRNAs)
- **Columns**: Patient samples (255 samples)
- **Values**: Normalized miRNA expression counts

**Biological Scale**: **MOLECULAR SCALE** - Post-transcriptional regulation
**Entity Types**:
- miRNAs (619 mature microRNAs)
- Samples (255 patient tissues)
- Expression values

---

#### 5. **filtered_clinical.csv** (256 rows)
- **Biological Entity**: Filtered clinical metadata
- **Format**: CSV
- **Subject Count**: 255 samples (matched to omics data)
- **Key Columns**:
  - `case_submitter_id`
  - `gender`
  - `race`
  - `ajcc_pathologic_stage`
- **Purpose**: Clinical metadata for samples with complete omics data

**Biological Scale**: **CLINICAL SCALE**

---

### DATA SUMMARY TABLE

| Data File | Rows | Cols | Entity | Biological Scale | Format |
|-----------|------|------|--------|------------------|--------|
| gene_disease.tsv | 8,498 | 10 | Gene-Disease-Drug | Molecular | TSV |
| pathway(最终版).tsv | 360 | 6 | Pathway-Gene | Molecular | TSV |
| Related Pathway.txt | 1,585 | 2 | Pathway-Pathway | Molecular | TSV |
| merged_output.tsv | 129,521 | 4 | Gene-Gene edges | Molecular | TSV |
| clinical.tsv | 462 | 5 | Patient metadata | Clinical/Population | TSV |
| filtered_hiseq_data.csv | 14,521 | 255 | Gene expression | Molecular | CSV |
| filtered_methylation_data.csv | 14,521 | 255 | DNA methylation | Molecular | CSV |
| filtered_miRNA_with_names.csv | 620 | 255 | miRNA expression | Molecular | CSV |
| filtered_clinical.csv | 256 | 4 | Filtered patient data | Clinical | CSV |

---

## 🎯 PART 2: APPLICATION TABS & FUNCTIONALITY

### **Tab 0: 🔬 Gene Network Visualization**
**Biological Scale**: **MOLECULAR SCALE**

**Functionality**:
1. Select disease → visualize gene interaction/regulation networks
2. Gene-to-pathway query: Input gene name → get associated pathways
3. Pathway network overlay: Show genes in disease context vs. pathway context
4. Drug listing: Show treatment options for selected disease

**Visualizations**:
- Force-directed gene network graph (Plotly)
- Node colors indicate: Hub genes, regulators, targets
- Edge types: Interaction (undirected) vs. Regulation (directed)
- Pathway network with color-coding (pink = disease-pathway overlap, blue = pathway-only)

**Data Sources**: gene_disease.tsv, pathway files, merged_output.tsv
**Network Operations**: 
- Network density calculation
- Clustering coefficient
- Degree centrality (Hub gene detection)

---

### **Tab 1: 📊 Network Model Calculation (IS Coefficient)**
**Biological Scale**: **MOLECULAR SCALE**

**Functionality**:
1. Select disease
2. Multi-select pathways (auto-filters to overlapping pathways)
3. Calculate Influence Score (IS) = f(coverage, pathway-disease overlap, connection density)
4. Rank pathways by IS coefficient

**Computation**:
- Coverage rate: % of disease genes in pathway
- Pathway overlap: % of pathway genes in disease
- Connection density: Interaction density among overlapping genes
- IS = Weighted combination of above metrics

**Outputs**:
- Bar chart: IS scores ranked by pathway
- Data table: Pathway name, IS coefficient, statistical metrics
- Status updates: Computation progress

**Data Sources**: Preprocessed disease/pathway objects from pickle

---

### **Tab 2: 📈 Data Statistics**
**Biological Scale**: **MOLECULAR SCALE** & **POPULATION SCALE**

**Functionality**:
1. Database type selector: "Gene Network" vs. "Social Network"
2. Display aggregate statistics

**Gene Network Statistics**:
- Total connections: 404,830,036
- Interaction connections: Gene-gene undirected edges
- Regulation connections: Gene-gene directed edges
- Network density, clustering coefficient

**Social Network Statistics**:
- Total social connections: 75,134,767
- Community structure statistics
- Node degree distribution

**Visualizations**: HTML cards with styled statistics boxes

---

### **Tab 3: 🌐 Social Network Simulation (SIS Model)**
**Biological Scale**: **POPULATION SCALE** (Analogous to disease spread)

**Functionality**:
1. Build network with community structure (LFR benchmark networks)
2. Run SIS epidemic simulation
3. Visualize network structure and infection dynamics

**Parameters**:
- **Network**: N (50-200 nodes), c (2-10 communities), k (avg degree), Z_in (community cohesion)
- **Transmission**: β (infection rate, 0.01-0.2), γ (recovery rate, 0.05-0.5)
- **Dynamics**: Initial infection ratio, simulation steps (50-200)

**Outputs**:
- Subtab "Network Structure": Community network visualization with node colors
- Subtab "Infection Dynamics": Time-series line plot of infection density
- Subtab "Infection Snapshot": State at specific timestep (susceptible/infected coloring)

**Data Sources**: Simulator generates synthetic networks and runs dynamics

**Models Implemented**:
- SIS: P(S→I) = 1-(1-β)^k, P(I→S) = γ
- Community detection via LFR model

---

### **Tab 4: 🧬 Gene Network Simulation (TCGA-COAD Analysis)**
**Biological Scale**: **MOLECULAR/CELLULAR SCALE**

**Functionality**:
1. Load TCGA-COAD real patient data
2. Select analysis dimension: Age, Sex, Disease Stage
3. Build gene/miRNA networks using MRNetB algorithm
4. Stratified analysis: Separate networks for age groups / sexes / stages

**Analysis Types**:
- **Age Grouping**: Young (0-50), Middle (50-70), Old (70+)
- **Sex Grouping**: Male vs. Female networks
- **Disease Stage**: Stage IA, IB, IIA, IIB, IIIA, IIIB, IV (+ optional normal samples)

**Data Types**:
- Gene expression (hiseq)
- miRNA expression
- Methylation patterns (optional)

**Algorithm**:
- **MRNetB** (Mutual Information-based Network Reconstruction):
  - Weight(i,j) = MI(i,j) - max(min(MI(i,:), MI(j,:)))
  - Filters out indirect associations
- Feature selection: Variance/mean-based or random sampling
- Max features: Adjustable (50-300) for performance tuning

**Outputs**:
- Subtab "Network Statistics": Density, clustering coefficient, degree distribution
- Subtab "Network Visualization": Interactive node-link diagram
- Subtab "Result Data": Edge list, centrality measures, module detection

**Data Sources**: TCGA-COAD filtered gene expression, methylation, miRNA, clinical metadata

---

### **Tab 5: 📚 Model Library**
**Biological Scale**: **MULTI-SCALE** (Molecular → Cellular → Population)

**Functionality**:
1. Display catalog of 6 available models
2. Model cards with: Name, scale, description, algorithm, input/output specs
3. Summary table: Model name, scale, algorithm, I/O
4. Scale distribution pie/bar chart

**6 Models Defined** (from model_library.py):

| Model | Scale | Algorithm | Input | Output |
|-------|-------|-----------|-------|--------|
| **Gene Interaction Network** | 🧬 Molecular | Graph theory + degree centrality | Gene-disease DB | Network topology, Hub genes |
| **Gene Regulatory Network** | 🧬 Molecular | Directed graph + regulatory cascade | Gene-disease DB | Regulatory networks, TF rankings |
| **IS Coefficient** | 🧬 Molecular | Influence score computation | Gene-pathway mapping | IS rankings, pathway impact |
| **MRNetB** | 🔬 Cellular/Tissue | Mutual information network inference | TCGA gene/miRNA expression | Gene association networks, edge weights |
| **Network Flow Entropy (MFE)** | 🔬 Cellular/Tissue | Information-theoretic entropy on networks | Gene expression + network topology | Node entropy values, information flow |
| **SIS Epidemic Dynamics** | 👥 Population | SIS transmission dynamics on networks | Network params + β/γ rates | Infection curves, community effects |

**Visualizations**: 
- HTML cards (2-column grid layout)
- Interactive parameter details (expandable sections)

---

### **Tab 6: 🔗 Multi-Scale Linkage Analysis**
**Biological Scale**: **MULTI-SCALE** (Molecular → Cellular → Population)

**Functionality**: Cross-scale parameter cascade
- **Molecular Layer** (Tab 0): Select disease → extract Hub genes + network density
- **Cellular Layer** (Tab 4): Hub genes as MRNetB seed nodes → TCGA networks → extract clustering coefficient
- **Population Layer** (Tab 3): Use molecular density + cellular clustering as initial conditions for SIS model

**Subtabs**:
1. **"Cascade Analysis"**: Three-layer parameter passing
   - Layer 1→2: Hub genes from disease network seed MRNetB selection
   - Layer 2→3: Network density + clustering coefficient inform β and c parameters
   - Outputs: 3-layer network comparison, radar chart (multi-dimensional metrics), summary table

2. **"Multi-Disease Comparison"** (optional subtab): Compare IS scores across diseases
   
3. **"Gene Tracking"** (optional subtab): Trace specific genes through multi-scale model

**Data Sources**: All above (gene_disease, TCGA-COAD, social network simulator)

---

## 🔍 PART 3: BIOLOGICAL SCALES COVERED

### **Scale 1: MOLECULAR SCALE** ✓ Comprehensive
- **Gene Networks**: Interaction and regulatory networks (4,322 genes)
- **Pathways**: 347 metabolic and signaling pathways  
- **Diseases**: 2,503 diseases linked to genes/pathways
- **Data**: Gene-disease-drug associations, pathway memberships
- **Models**: Gene Interaction, Gene Regulatory, IS Coefficient, MRNetB
- **Tabs**: 0 (Gene Network Viz), 1 (IS Coefficient), 5 (Model Library)

### **Scale 2: CELLULAR/TISSUE SCALE** ✓ Partial
- **TCGA-COAD Data**: 
  - Gene expression (14,520 genes × 255 samples)
  - miRNA expression (619 miRNAs × 255 samples)
  - DNA methylation (14,520 sites × 255 samples)
- **Patient Stratification**: Age, sex, disease stage
- **Models**: MRNetB network inference, MFE entropy
- **Tabs**: 4 (TCGA-COAD Gene Network Simulation)

### **Scale 3: POPULATION SCALE** ✓ Implemented
- **Social Networks**: Community-structured synthetic networks
- **Transmission Models**: SIS epidemic dynamics
- **Parameters**: Node count, community structure, infection rates
- **Models**: SIS Epidemic Dynamics
- **Tabs**: 3 (Social Network Simulation)

### **Scale 4: MULTI-SCALE INTEGRATION** ✓ Emerging
- **Cross-Scale Engine**: Automated parameter cascade
- **Tab 6**: Multi-scale linkage (molecular → cellular → population)

---

## 📈 PART 4: DATA STATISTICS SUMMARY

| Category | Count | Notes |
|----------|-------|-------|
| **Diseases** | 2,503 | KEGG disease database |
| **Genes** | 4,322 | Unique gene symbols (after alias handling) |
| **Pathways** | 347 | KEGG Homo sapiens pathways |
| **Pathway-Pathway Edges** | 1,585 | Related pathway network |
| **Gene-Disease Records** | 8,498 | Multiple genes per disease, multiple pathways per gene |
| **Gene-Gene Edges** | 129,521 | Pre-computed interaction network |
| **TCGA-COAD Patients** | 461 | Total patient records |
| **TCGA-COAD Samples with Complete Omics** | 255 | Matched expression + methylation + miRNA + clinical |
| **Genes in TCGA Expression Data** | 14,520 | RNA-seq from HiSeq |
| **miRNAs in TCGA** | 619 | Small RNA-seq |
| **Methylation Sites in TCGA** | 14,520 | CpG sites (450K array) |

---

## 🛠️ PART 5: DATA PROCESSING PIPELINE

### Load Order:
1. **Real Data Loader** (data_loader.py):
   - Load gene_disease.tsv
   - Load pathway(基因名映射版).tsv
   - Load Related Pathway.txt
   - Build disease/pathway/gene indices

2. **Preprocessing** (preprocess_data.py):
   - Parse gene_disease.tsv → Disease objects (with genes, drugs, interactions, regulations)
   - Parse pathway TSV → Pathway objects (with genes, interactions)
   - Generate synthetic network relationships (gene-gene for visualization)
   - Save as pickle (data/preprocessed_data.pkl) for fast loading

3. **TCGA-COAD Simulator** (tcga_coad_simulator.py):
   - Load clinical.tsv + filtered_hiseq_data.csv + filtered_methylation_data.csv + filtered_miRNA_with_names.csv
   - Stratify by age/sex/stage
   - Apply MRNetB algorithm on expression data
   - Generate networks per group

4. **Social Network Simulator** (social_network_sim.py):
   - Generate LFR benchmark networks with community structure
   - Simulate SIS epidemic dynamics

---

## 💡 PART 6: KEY INSIGHTS

**Multi-Scale Integration**:
- Molecular layer (genes/pathways/diseases) is comprehensive (2,503 diseases, 4,322 genes)
- Cellular layer (TCGA-COAD) provides real patient omics data (255 matched samples)
- Population layer (SIS model) simulates transmission at network level
- Cross-scale engine attempts to link: Hub genes → cellular networks → population dynamics

**Data Coverage**:
- Strong in knowledge base (KEGG genes, diseases, pathways)
- Real-world transcriptomics data available (RNA-seq, miRNA, methylation)
- Patient stratification possible (age, sex, stage)
- Synthetic social networks for population-level modeling

**Visualization & Models**:
- 6 computational models spanning 3 scales
- Interactive Gradio web interface
- Network visualization, statistics, time-series dynamics
- Stratified analysis capabilities

