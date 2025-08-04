# scripts/merge_counts.py

import scanpy as sc
import pandas as pd
import os
import argparse

# -------------------
# Argument Parser
# -------------------
parser = argparse.ArgumentParser(description="Merge raw count matrices into one AnnData object")
parser.add_argument("--input_dir", type=str, required=True, help="Directory with raw CSV count files")
parser.add_argument("--metadata_file", type=str, required=True, help="Path to metadata CSV")
parser.add_argument("--output", type=str, required=True, help="Path to save merged AnnData (.h5ad)")
args = parser.parse_args()

# -------------------
# Load Metadata
# -------------------
metadata = pd.read_csv(args.metadata_file)
print(f"Loaded metadata for {len(metadata)} samples")

# -------------------
# Load and Merge Counts
# -------------------
adatas = []
for _, row in metadata.iterrows():
    sample_id = row['sample']
    condition = row['condition']
    filepath = os.path.join(args.input_dir, f"{sample_id}.csv")

    print(f"Loading: {filepath}")
    counts = sc.read_csv(filepath).T
    counts.obs['sample'] = sample_id
    counts.obs['condition'] = condition
    adatas.append(counts)

adata_merged = adatas[0].concatenate(adatas[1:], batch_key="sample_id", batch_categories=metadata['sample'].tolist())
print(f"Merged AnnData shape: {adata_merged.shape}")

# -------------------
# Save Output
# -------------------
adata_merged.write(args.output)
print(f"Saved merged data to {args.output}")
