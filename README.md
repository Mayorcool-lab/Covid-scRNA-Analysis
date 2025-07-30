 # 🦠 Single-Cell Transcriptomic Analysis of COVID-19 Infected Human Lungs | Scanpy + Snakemake

This project presents a **reproducible single-cell RNA-seq (scRNA-seq) analysis** pipeline for **post-mortem lung samples** from COVID-19 patients using **Scanpy**, **Snakemake**, **Scrublet**, and Python data science tools. The pipeline is modular, version-controlled, and optimized for reproducibility and scalability.

---

## 🎯 Objective

To analyze cellular heterogeneity in COVID-19-infected lungs, identify **highly variable genes (HVGs)**, perform **clustering**, **UMAP visualization**, and prepare the data for **downstream differential expression** and **cell–cell interaction analysis**. The project aims to pinpoint candidate genes and pathways involved in **inflammation, impaired regeneration**, and **fibrosis**.

---

## 📚 Dataset

- **Source**: [COVID-19 lung single-nucleus RNA-seq study (Nature, 2021)](https://doi.org/10.1038/s41586-021-03504-w)
- **Samples**: Autopsy-derived lung tissues from COVID-19 fatalities and healthy controls
- **Format**: Preprocessed gene-cell count matrices (`.csv`) per sample

---

## 🛠️ Tools & Packages

- `Scanpy`, `anndata`, `matplotlib`, `numpy`, `pandas`
- `Scrublet` (doublet detection)
- `Snakemake` (pipeline automation)
- `Git LFS` (for tracking large `.h5ad` files)
- Reproducibility via `conda` (`envs/scanpy_env.yml`)

---

## 🔄 Analysis Workflow

### 📁 Phase 1: Data Import & Merging
- Loaded raw count matrices per sample
- Merged into single `AnnData` object
- Tracked sample and condition metadata

### 🧹 Phase 2: Quality Control
- Filtered out low-quality cells and genes
- Removed cells with high mitochondrial/ribosomal content
- Detected and excluded doublets with Scrublet

### 🧪 Phase 3: Normalization & HVG Selection
- Total-count normalization and log1p transform
- Identified highly variable genes (HVGs)
- Saved HVG plot and filtered object

### 📊 Phase 4: Dimensionality Reduction & Clustering
- Scaled data and ran PCA
- Computed neighborhood graph
- Performed Leiden clustering (resolution = 0.5)
- Visualized results with UMAP
- Plotted clusters, QC metrics, and batch effects

📌 *Status: downstream steps such as DE analysis, marker identification, and cell-type annotation are in progress.*

---

## 📁 Project Structure

