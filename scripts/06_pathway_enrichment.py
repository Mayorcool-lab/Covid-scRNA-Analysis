# scripts/06_pathway_enrichment.py

"""
Step 6: Pathway Enrichment Analysis using gseapy (Enrichr)
Analyzes upregulated DEGs in COVID lungs for:
- Alveolar Macrophages
- NK Cells
Databases: GO BP, KEGG, Reactome
"""

import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import random

# Set seed
random.seed(42)
np.random.seed(42)

# Load DEG files
sig_alveolar = pd.read_csv("../results/deg/Significant_DEG_Alveolar_macrophages_COVID.csv")
sig_nk = pd.read_csv("../results/deg/Significant_DEG_NK_cells_COVID.csv")

# Define gene sets
gene_sets = [
    "GO_Biological_Process_2021",
    "KEGG_2021_Human",
    "Reactome_2022"
]

# Dictionary of cell type → DEG dataframe
deg_dict = {
    "Alveolar_Macrophages": sig_alveolar,
    "NK_Cells": sig_nk
}

# Loop over both cell types
for celltype, df in deg_dict.items():
    gene_list = df["names"].tolist()
    
    for gs in gene_sets:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gs,
            organism='Human',
            outdir=None,
            cutoff=0.05
        )
        result_df = enr.results.sort_values("Adjusted P-value").head(10)
        
        # Plot barplot
        plt.figure(figsize=(8, 5))
        sns.barplot(
            x=-np.log10(result_df["Adjusted P-value"]),
            y=result_df["Term"],
            color='teal'
        )
        plt.xlabel("-log₁₀(Adjusted P-value)")
        plt.title(f"{celltype.replace('_', ' ')} – {gs.split('_')[0]}")
        plt.tight_layout()
        fname = f"../results/figures/{celltype.lower()}_{gs.split('_')[0].lower()}.png"
        plt.savefig(fname, dpi=300)
        plt.close()
        print(f"Saved: {fname}")
