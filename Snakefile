rule all:
    input:
        "results/anndata/clustered_final_celltypes.h5ad",
        "results/deg/DEG_Alveolar_macrophages_COVID_vs_Control.csv",
        "results/deg/DEG_NK_cells_COVID_vs_Control.csv",
        "results/deg/Significant_DEG_Alveolar_macrophages_COVID.csv",
        "results/deg/Significant_DEG_NK_cells_COVID.csv",
        "figures/dotplot_COVID_DEG_dotplot.png",
        "figures/heatmap_COVID_DEG_heatmap.png",
        "figures/alveolar_macrophages_go.png",
        "figures/nk_cells_go.png",
        "figures/alveolar_macrophages_kegg.png",
        "figures/nk_cells_kegg.png",
        "figures/alveolar_macrophages_reactome.png",
        "figures/nk_cells_reactome.png"

rule merge_counts:
    output:
        "results/anndata/merged_counts.h5ad"
    script:
        "scripts/merge_counts.py"

rule quality_control:
    input:
        "results/anndata/merged_counts.h5ad"
    output:
        "results/anndata/qc_nodoublets.h5ad"
    script:
        "scripts/01_quality_control.py"

rule preprocessing_and_clustering:
    input:
        "results/anndata/qc_nodoublets.h5ad"
    output:
        "results/anndata/clustered.h5ad"
    script:
        "scripts/02_preprocessing_and_clustering.py"

rule marker_annotation:
    input:
        "results/anndata/clustered.h5ad"
    output:
        "results/markers/top5_marker_genes_per_cluster.csv"
    script:
        "scripts/03_marker_annotation.py"

rule celltype_annotation:
    input:
        "results/anndata/clustered.h5ad",
        "results/markers/top5_marker_genes_per_cluster.csv"
    output:
        "results/anndata/clustered_final_celltypes.h5ad"
    script:
        "scripts/04_celltype_annotation.py"

rule differential_expression:
    input:
        "results/anndata/clustered_final_celltypes.h5ad"
    output:
        "results/deg/DEG_Alveolar_macrophages_COVID_vs_Control.csv",
        "results/deg/DEG_NK_cells_COVID_vs_Control.csv",
        "results/deg/Significant_DEG_Alveolar_macrophages_COVID.csv",
        "results/deg/Significant_DEG_NK_cells_COVID.csv",
        "figures/dotplot_COVID_DEG_dotplot.png",
        "figures/heatmap_COVID_DEG_heatmap.png"
    script:
        "scripts/05_differential_expression.py"

rule pathway_enrichment:
    input:
        "results/deg/Significant_DEG_Alveolar_macrophages_COVID.csv",
        "results/deg/Significant_DEG_NK_cells_COVID.csv"
    output:
        "figures/alveolar_macrophages_go.png",
        "figures/nk_cells_go.png",
        "figures/alveolar_macrophages_kegg.png",
        "figures/nk_cells_kegg.png",
        "figures/alveolar_macrophages_reactome.png",
        "figures/nk_cells_reactome.png"
    script:
        "scripts/06_pathway_enrichment.py"
 
