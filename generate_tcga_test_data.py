"""Generate test TCGA-COAD data for WGCNA analysis"""
import pandas as pd
import numpy as np
import os

np.random.seed(42)

# Create data directory
os.makedirs('data/TCGA-COAD', exist_ok=True)

# Generate gene expression data (14,520 genes × 255 samples)
print("Generating gene expression data...")
n_genes = 14520
n_samples = 255

# Create realistic gene expression with some correlation structure
gene_expr = np.random.negative_binomial(5, 0.3, size=(n_genes, n_samples))
gene_expr = gene_expr.astype(np.float32)
# Log transform (typical for RNA-seq)
gene_expr = np.log2(gene_expr + 1)

# Create gene names
gene_names = [f"ENSG{i:08d}" for i in range(n_genes)]

# Create sample names
sample_names = [f"TCGA-{i:04d}" for i in range(n_samples)]

gene_expr_df = pd.DataFrame(gene_expr, index=gene_names, columns=sample_names)
gene_expr_df.to_csv('data/TCGA-COAD/filtered_hiseq_data.csv')
print(f"✓ Saved gene expression: {gene_expr_df.shape}")

# Generate clinical data (255 samples × 5 clinical variables)
print("Generating clinical data...")
clinical_data = pd.DataFrame({
    'Age': np.random.randint(30, 85, n_samples),
    'Gender': np.random.choice(['M', 'F'], n_samples),
    'Stage': np.random.choice(['I', 'II', 'III', 'IV'], n_samples),
    'Grade': np.random.choice(['G1', 'G2', 'G3'], n_samples),
    'Status': np.random.choice(['Alive', 'Dead'], n_samples)
}, index=sample_names)

# Create age groups
clinical_data['Age_Group'] = pd.cut(clinical_data['Age'], bins=[0, 50, 65, 100], labels=['Young', 'Middle', 'Old'])

clinical_data.to_csv('data/TCGA-COAD/filtered_clinical.csv')
print(f"✓ Saved clinical data: {clinical_data.shape}")

# Generate miRNA expression data (619 miRNAs × 255 samples)
print("Generating miRNA expression data...")
n_mirnas = 619
mirna_expr = np.random.negative_binomial(3, 0.4, size=(n_mirnas, n_samples))
mirna_expr = mirna_expr.astype(np.float32)
mirna_expr = np.log2(mirna_expr + 1)

# Create miRNA names
mirna_names = [f"hsa-miR-{i}" for i in range(n_mirnas)]

mirna_expr_df = pd.DataFrame(mirna_expr, index=mirna_names, columns=sample_names)
mirna_expr_df.to_csv('data/TCGA-COAD/filtered_miRNA_with_names.csv')
print(f"✓ Saved miRNA expression: {mirna_expr_df.shape}")

print("\n✅ All test data generated successfully!")
print(f"   Expression data: {gene_expr_df.shape}")
print(f"   Clinical data: {clinical_data.shape}")
print(f"   miRNA data: {mirna_expr_df.shape}")
