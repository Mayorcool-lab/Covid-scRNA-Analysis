# COVID-19 Lung scRNA-seq Analysis (Scanpy Pipeline)

![scRNA-seq Workflow](https://img.shields.io/badge/Analysis-scanpy-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)

## 🔍 Project Overview
This repository presents a complete Python-based single-cell RNA-seq (scRNA-seq) pipeline analyzing COVID-19 lung tissue data from the study:
**"A molecular single-cell lung atlas of lethal COVID-19"**

- Annotate immune and epithelial populations in lethal COVID-19 lungs.
- Discover COVID-driven changes in cell proportions and gene expression.
- Identify cell-type-specific biomarkers and enriched immune pathways.

### Key Features
- ✅ **Modernized workflow**: Replaces original Seurat/R implementation with Scanpy/Python
- ✅ **End-to-end analysis**: From raw data merging to pathway enrichment
- ✅ **Reproducible**: Showcase best practices in scRNA-seq analysis using Scanpy, Git, and Snakemake.


## 🧪 Biological Objectives
1. Identify cell populations enriched in COVID-19 lungs
2. Discover cell-type-specific differentially expressed genes (DEGs)
3. Characterize altered biological pathways in:
   - Alveolar Macrophages
   - NK Cells

## 🛠️ Complete Analysis Pipeline

The analysis pipeline is structured as a series of Jupyter notebooks, which are also integrated into an automated **Snakemake** workflow for full reproducibility.

### 📁 0. Data Merging (`merge_counts.ipynb`)
**Goal**: Combine raw count matrices into unified AnnData object
**Key Steps**:
- Load multiple sample-specific CSV files
- Concatenate with sample/condition metadata
- Output: `merged_counts.h5ad`

### 📁 1. Quality Control (`01_quality_control.ipynb`)
**Goal**: Filter low-quality cells and genes
**Methods**:
- Remove cells with:
  - High mitochondrial content (>10%)
  - Low gene counts (<200 genes/cell)
- Doublet detection (Scrublet)
**Output**: `qc_nodoublets.h5ad`

### 📁 2. Preprocessing & Clustering (`02_preprocessing_and_clustering.ipynb`)
**Goal**: Normalize and cluster cells
**Workflow**:
1. Normalization: `sc.pp.normalize_total()` + `sc.pp.log1p()`
2. Feature selection: HVGs (`sc.pp.highly_variable_genes`)
3. Dimensionality reduction: PCA → UMAP
4. Clustering: Leiden algorithm
**Output**: Cluster-labeled UMAP plots

### 📁 3. Marker Identification (`03_marker_annotation.ipynb`)
**Goal**: Find cluster biomarkers
**Methods**:
- Wilcoxon rank-sum test (`sc.tl.rank_genes_groups`)
- Top 5 markers per cluster exported to CSV

### 📁 4. Cell Type Annotation (`04_celltype_annotation.ipynb`)
**Goal**: Assign biological identities
**Approach**:
1. Automated annotation: CellTypist (immune/epithelial models)
2. Manual validation: Canonical marker expression
3. Finalized 10 major lung cell types
**Output**: Annotated UMAPs + cell type proportion plots

### 📁 5. Differential Expression (`05_differential_expression.ipynb`)
**Goal**: Identify COVID-associated DEGs
**Focus Populations**:
- Alveolar Macrophages
- NK Cells
**Thresholds**: adj-p < 0.05, logFC > 0
**Visualizations**: Dot plots, heatmaps

### 📁 6. Pathway Enrichment (`06_pathway_enrichment.ipynb`)
**Goal**: Uncover altered biological pathways
**Databases**:
- GO Biological Processes
- KEGG
- Reactome
**Key Findings**:
- Macrophages: Inflammatory cytokine signaling
- NK Cells: Cytotoxic effector functions


## 📁 Repo Structure (Planned)

```bash
Covid-scRNA-Analysis/
├── notebook/           # ✅ All updated notebooks (merge_counts → 06)
├── scripts/            # 🔜 Python modules for each phase
├── Snakefile           # 🔜 Snakemake pipeline for automation
├── envs/               # ✅ Conda environment definitions
├── .gitignore          # ✅ Ignore .csv/.png/.h5ad outputs
├── README.md           # ✅ You’re reading it!
