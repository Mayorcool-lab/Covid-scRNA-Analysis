# scripts/02_preprocessing_and_clustering.py

"""
Step 2: Preprocessing, Dimensionality Reduction, Clustering, and UMAP
"""

import scanpy as sc
import numpy as np
import pandas as pd
import random

# Set seeds for reproducibility
random.seed(42)
np.random.seed(42)

# Load quality-controlled data
adata = sc.read("../results/anndata/qc_nodoublets.h5ad")

# Step 1: Normalize & log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Step 2: Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable].copy()

# Step 3: Scale and regress out unwanted sources of variation
sc.pp.scale(adata, max_value=10)

# Step 4: PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, show=False)

# Step 5: Neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# Step 6: Clustering
sc.tl.leiden(adata, resolution=0.5)

# Save the object for next steps
# adata.write("../results/anndata/clustered.h5ad")
