# scripts/04_celltype_annotation.py

"""
Step 4: Cell Type Annotation using CellTypist and Manual Curation
"""

import scanpy as sc
import pandas as pd
import numpy as np
import random
import celltypist

# Set seeds
random.seed(42)
np.random.seed(42)

# Load clustered AnnData
adata = sc.read("../results/anndata/clustered.h5ad")

# Step 1: Run CellTypist with immune + epithelial model
model = celltypist.models.download_model('Immune_All_Low')
predictions = celltypist.annotate(adata, model=model, majority_voting=True)

# Add prediction to AnnData
adata.obs['celltypist_labels'] = predictions.predicted_labels

# Step 2: Assign dominant cell type per cluster (Leiden)
dominant = (
    adata.obs.groupby('leiden')['celltypist_labels']
    .agg(lambda x: x.value_counts().idxmax())
    .rename("dominant_cell_type")
)
adata.obs = adata.obs.join(dominant, on='leiden')

# Step 3: Manual curation (placeholder — replace with curated mapping)
cluster_to_type = {
    '0': 'Alveolar macrophages',
    '1': 'Epithelial cells',
    '2': 'NK cells',
    # Add all cluster IDs as needed
}
adata.obs['final_celltype'] = adata.obs['leiden'].map(cluster_to_type)

# Save final annotated object
adata.write("../results/anndata/clustered_final_celltypes.h5ad")
