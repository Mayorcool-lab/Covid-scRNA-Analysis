# 🧬 Single-Cell RNA-Seq Analysis of COVID-19 Lung Cells | Seurat + SingleR

This project delivers an end-to-end single-cell transcriptomic analysis comparing **COVID-19 patients** and **healthy controls** using the **Seurat** framework in R. It includes quality control, batch correction, clustering, differential expression, and cell-type annotation with **Azimuth** (lung reference).

---

## 🎯 Objective

To explore cell-type-specific transcriptional changes in lung tissue during COVID-19 infection and identify key differentially expressed genes (DEGs) across conditions.

---

## 📁 Dataset

- **Samples**: 10 single-cell RNA-seq samples  
  - `5 COVID-19 patient lungs`  
  - `5 healthy control lungs`  
- **Source**: GEO (e.g., GSM5226574, GSM5226581, etc.)
- **Input format**: Raw gene expression count matrices (CSV)

---

## 🛠️ Tools & Packages Used

- `Seurat` for scRNA-seq analysis  
- `SingleR`, `celldex` for automated cell annotation  
- `pheatmap`, `ggplot2`, `ggrepel` for visualization  
- `Azimuth (HumanPrimaryCellAtlasData)` for lung reference  
- `DAVID` (external) for enrichment post-DEG analysis

---

## 🧪 Analysis Pipeline

### 1️⃣ Data Ingestion & Merging
- Loaded 10 raw count matrices
- Created Seurat objects with metadata labels (`Condition = covid/normal`)
- Merged all samples into one object

### 2️⃣ Quality Control
- Filtered cells with:
  - `nFeature_RNA > 200`
  - `nFeature_RNA < 3000`
  - `percent.mt < 10%`
- Visualized metrics using violin plots

### 3️⃣ Normalization & Variable Gene Detection
- Performed log normalization
- Identified highly variable genes
- Highlighted top 10 variable genes

### 4️⃣ Dimensionality Reduction
- Performed PCA
- Determined optimal PCs via JackStraw + ElbowPlot

### 5️⃣ Clustering
- Ran UMAP for nonlinear projection
- Tested resolutions from `0.2` to `2.0`
- Final clustering at `resolution = 0.2`

### 6️⃣ Batch Correction & Integration
- Split by condition
- Performed data integration using Canonical Correlation Analysis (CCA)
- Re-ran PCA and UMAP on integrated data

### 7️⃣ Cell Annotation
- Annotated cells using `SingleR` and `HumanPrimaryCellAtlasData` lung reference
- Added cell-type metadata
- Visualized clusters by predicted cell type on UMAP

### 8️⃣ Differential Expression Analysis
- Identified cluster-specific markers using `FindAllMarkers()`
- Compared gene expression between:
  - COVID vs Normal
  - Clusters across each condition
  - SFTPB-high vs SFTPB-low cells

### 9️⃣ Marker Export
- Saved:
  - Cluster-specific markers
  - DEGs by condition
  - SFTPB-specific markers
- Exported to CSVs for further use (e.g., DAVID enrichment)

---

## 📊 Outputs

| File | Description |
|------|-------------|
| `top_20_upregulated_genes.csv` | DEGs highly expressed in COVID samples |
| `top_20_downregulated_genes.csv` | Genes suppressed in COVID lungs |
| `conditions_marker_files/*.csv` | Per-cluster markers under each condition |
| `covid.Rdata`, `covid.Rhistory` | R workspace and history |

---

## 📈 Visualizations

- PCA & UMAP plots colored by:
  - Condition
  - Cell type (SingleR)
  - Clusters
- Heatmaps for top DEGs
- Violin plots of QC metrics
- Pie charts for SFTPB-expressing cells

---

## 🧠 Skills Demonstrated

- Advanced use of Seurat and scRNA-seq workflows  
- Batch effect correction and sample integration  
- Multi-level DEG analysis (cluster, condition, marker-specific)  
- Automated cell-type annotation with Azimuth  
- Pipeline modularization and scalable design  
- Data export and downstream compatibility (e.g., DAVID)

---

## 💡 Key Insights

- Identified **distinct transcriptional signatures** in COVID-infected lung cells  
- Found **SFTPB-overexpressing cells** enriched in COVID samples  
- Annotated a diverse population of immune and epithelial cells  
- Built a scalable workflow for large multi-condition scRNA-seq studies

---
