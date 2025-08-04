# scripts/05_differential_expression.py

"""
Step 5: Differential Expression Analysis (COVID vs Control)
for selected cell types: Alveolar macrophages and NK cells
"""

import scanpy as sc
import pandas as pd
import numpy as np
import random

# Set seeds for reproducibility
random.seed(42)
np.random.seed(42)

# Load annotated dataset
adata = sc.read("../results/anndata/clustered_final_celltypes.h5ad")

# Define cell types of interest
celltypes_of_interest = ['Alveolar macrophages', 'NK cells']

# Loop through each cell type
for celltype in celltypes_of_interest:
    print(f"Running DE analysis for: {celltype}")
    
    # Subset the data
    adata_subset = adata[adata.obs['final_celltype'] == celltype].copy()
    
    # Ensure 'condition' is categorical
    adata_subset.obs['condition'] = adata_subset.obs['condition'].astype('category')
    
    # Run differential expression
    sc.tl.rank_genes_groups(
        adata_subset,
        groupby='condition',
        method='wilcoxon',
        use_raw=True,
        pts=True
    )
    
    # Save full DEG result
    df = sc.get.rank_genes_groups_df(adata_subset, group='COVID')
    outpath = f"../results/deg/DEG_{celltype.replace(' ', '_')}_COVID_vs_Control.csv"
    df.to_csv(outpath, index=False)
    print(f"Saved: {outpath}")

    # Filter for significantly upregulated genes in COVID
    df_sig = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > 0)].sort_values('pvals_adj')
    df_sig.to_csv(
        f"../results/deg/Significant_DEG_{celltype.replace(' ', '_')}_COVID.csv",
        index=False
    )
