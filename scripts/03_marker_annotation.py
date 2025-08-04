# scripts/03_marker_annotation.py

"""
Step 3: Marker Gene Identification for Each Cluster
"""

import scanpy as sc
import pandas as pd
import numpy as np
import random

# Set seeds for reproducibility
random.seed(42)
np.random.seed(42)

# Load clustered AnnData
adata = sc.read("../results/anndata/clustered.h5ad")

# Step 1: Identify marker genes
sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    method='wilcoxon',
    use_raw=True,
    pts=True
)

# Step 2: Export top 5 marker genes per cluster
markers = sc.get.rank_genes_groups_df(adata, group=None)
top5 = markers.groupby('group').head(5)
top5.to_csv("../results/markers/top5_marker_genes.csv", index=False)

# Optional: Plot (disabled in script mode)
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, show=True)

# Save AnnData with markers (optional)
# adata.write("../results/anndata/with_markers.h5ad")
