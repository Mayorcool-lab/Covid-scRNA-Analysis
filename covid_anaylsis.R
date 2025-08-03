#### Single-cell RNA analysis of COVID and save workplace####

save.image("covid.Rdata")
savehistory("covid.Rhistory")
loadhistory("covid.Rdata")
load("covid.Rdata")

####Load packages####
install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("celldex")

BiocManager::install("celldex", update = TRUE)
install.packages("pheatmap")
install.packages("viridis")

library(celldex)
library(SingleR)
library(celldex)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(viridis)
library(ggplot2)
library(gridExtra)
library(ggrepel)

####Load data ####
# Define a vector with the paths to your CSV files
files <- c("GSM5226574_C51ctr_raw_counts.csv", "GSM5226576_C53ctr_raw_counts.csv", "GSM5226578_C55ctr_raw_counts.csv", "GSM5226575_C52ctr_raw_counts.csv", "GSM5226577_C54ctr_raw_counts.csv",
           "GSM5226582_L03cov_raw_counts.csv", "GSM5226584_L04covaddon_raw_counts.csv", "GSM5226581_L01cov_raw_counts.csv", "GSM5226583_L04cov_raw_counts.csv", "GSM5226585_L05cov_raw_counts.csv")

# Define a vector with labels for your CSV files
labels <- c(rep("covid", 5), rep("normal", 5))

seurat_list <- list()

for (i in seq_along(files)) {
  expression_data <- read.csv(files[i], row.names = 1)
  seurat_object <- CreateSeuratObject(counts = expression_data, project = paste0("Project", i), min.cells = 3, min.features = 200)
  
  # Add a metadata column
  seurat_object[["Condition"]] <- labels[i]
  
  seurat_list[[i]] <- seurat_object
}

  # Merge Seurat objects
merged_seurat_object <- merge(seurat_list[[1]], y = seurat_list[-1])

# Now we check the metadata
tail(merged_seurat_object@meta.data, 10) # we get the barcodes of the cells etc.

# Check if the data contain both normal and covid.
unique(merged_seurat_object@meta.data$Condition)

#### Pre-processing####
## We check the mitochondrial cell and check the number of unique genes i.e filter out low quality cell
#Mitochondrial percentage calculation

merged_seurat_object$percent_mt = PercentageFeatureSet(merged_seurat_object, pattern = "^MT-")
head(merged_seurat_object@meta.data, 5) # Reprint the metadata

#Visualization of with violin plot
VlnPlot(merged_seurat_object, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)

# Now filter assigend to a new varibale called filtered_data
filtered_data <- subset(merged_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt < 10)

# Check cells left with and recheck with violin plot
filtered_data 
VlnPlot(filtered_data, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)

####Normalization of the data on the basis of cells.Log Normalization ####

filtered_data = NormalizeData(filtered_data)
#### Variable feature selection.Choose highly variable genes####

filtered_data= FindVariableFeatures(filtered_data)
#Check top10 genes

top10 <- head(VariableFeatures(filtered_data), 10)
top10

#Visualize this top genes

plot1 <- plot(VariableFeaturePlot(filtered_data))

plot2 <- LabelPoints(plot = plot1, points = top10, repel = F)#Label the plot
plot1 + plot2

####Scale the Data scaling####

filtered_data <- ScaleData(filtered_data)
#Now assess the scaled data

head(filtered_data[["RNA"]]@scale.data, 5)

#### Linear dimension reduction PCA ####

filtered_data = RunPCA(filtered_data , features = VariableFeatures(filtered_data)) # we only consider the variable features

#Visualization of PCA with VizDim,  and Dimplot 

VizDimLoadings(filtered_data, dims = 1:2, reduction = "pca")

DimPlot(filtered_data, reduction = "pca")

DimHeatmap(filtered_data, dims = 1, cells = 500, balanced = T)

####Determining the dimentionality of the dataset####

filtered_data <- JackStraw(filtered_data, num.replicate = 100) #num.replicate can less or more
filtered_data <- ScoreJackStraw(filtered_data, dims = 1:20) #Get the score for the first 20 to check if they are significant

#Then visualize the significance
#Option 1
JackStrawPlot(filtered_data, dims = 1:20) 
#Option 2 using elbowplot

ElbowPlot(filtered_data)

#### Clustering #### 
#KNN graph using the Find Neighbor's function of the Seurat

filtered_data <- FindNeighbors(filtered_data, dims = 1:20)

#Find the best resolution
resolutions <- c(0.2, 0.5, 1.0, 1.5, 2.0)

for (res in resolutions) {
  filtered_data <- FindClusters(filtered_data, resolution = res)
  
  # Number of clusters
  num_clusters <- length(levels(Idents(filtered_data)))
  
  print(paste("Resolution:", res, "Number of clusters:", num_clusters))
  
  # Plot
  DimPlot(filtered_data, group.by = "ident") + ggtitle(paste("Resolution:", res))
}

# 0.2 = 14 clus,0.5 = 19 clus,1.0 = 27 clus,1.5 = 33 clus,2.0 = 39 clus

#Find the clusters based on the best resolution

filtered_data <- FindClusters(filtered_data, resolution = 0.2)

#Check the Idents or clusters

head(Idents(filtered_data),5)

####Non-linear dimension reduction####
filtered_data <- RunUMAP(filtered_data, dims = 1:20)

#Use dimplot to visualize and color to the conditions
umap_plot <- DimPlot(filtered_data, reduction = "umap", group.by = "Condition")
umap_plot

#### Perform integration to correct for batch effects####
obj.list <- SplitObject(filtered_data, split.by = "Condition")
for (i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
}

#Now we select the integration features and find integration anchors

features <-  SelectIntegrationFeatures(obj.list)

#Canonical Correlation Analysis (CCA) statistical method to understand the relationships between the variables

anchors <- FindIntegrationAnchors(object.list = obj.list,
                       anchor.features = features) 
  
#Now integrate the data

seurat.integrated <- IntegrateData(anchorset = anchors)  

#### Now perform scaling, PCA and UMAP on the integrated data####

seurat.integrated <- ScaleData(seurat.integrated)
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:20)

# Now visualize the data

DimPlot(seurat.integrated, group.by = "ident") + ggtitle("UMAP plot colored by cluster")
DimPlot(seurat.integrated, group.by = "Condition") + ggtitle("UMAP plot colored by condition")

# Check before and after integration

grid.arrange(umap_plot, umap_plot_integrated, ncol = 2, nrow = 2)

####Identify marker genes for each cluster####
cluster.markers <- FindAllMarkers(seurat.integrated) # i.e the 14 clusters

levels(seurat.integrated)

# Perform differential expression analysis between conditions within each cluster
# Fetch cluster identities
clusters <- levels(seurat.integrated)

#Option 1: FindMarkers() function is for between each cluster and all other cells,####
#irrespective of the condition ("covid" or "control")
# Initialize an empty list to store the markers for each cluster
marker.list <- list()

# Loop through each cluster
for (cluster in clusters) {
  
  # Find markers for current cluster against all other cells
  markers <- FindMarkers(seurat.integrated, ident.1 = cluster)
  
  # Store the markers in the list
  marker.list[[cluster]] <- markers
}

head(marker.list[["0"]]) # for cluster 1 or zero

# Loop through each cluster in the marker list for all cluster
for (cluster in names(marker.list)) {
  # Print the first six lines of the marker genes for the current cluster
  print(head(marker.list[[cluster]]))
  
  # Save the marker genes for the current cluster to a CSV file
  write.csv(marker.list[[cluster]], paste0(cluster, "new_clusters_markers.csv"))
}


#Option 2:Find markers for each cluster within each condition####
conditions <- unique(seurat.integrated$Condition)
clusters <- levels(seurat.integrated)

# Loop through each condition
for (cond in conditions) {
  # Subset Seurat object by condition
  seurat.sub <- subset(seurat.integrated, Condition == cond)
  
  # Loop through each cluster
  for (cluster in clusters) {
    # Find markers for current cluster against all other cells
    markers <- FindMarkers(seurat.sub, ident.1 = cluster)
    
    # Store the markers in the list, naming the list element with the cluster and condition
    marker.list[[paste0("Condition_", cond, "_Cluster_", cluster)]] <- markers
  }
}

####save all the markers for each condition

# Create a directory to store the marker files
dir.create("conditions_marker_files")

# Loop through each element of the marker list
for (name in names(marker.list)) {
  # Write the markers for the current element to a CSV file
  write.csv(marker.list[[name]], file = paste0("conditions_marker_files/", name, "condition_markers.csv"))
}




####to get the differentially expressed genes between conditions####
# Set the active idents to "Condition"
Idents(seurat.integrated) <- "Condition"

# Now perform differential expression analysis
de.genes <- FindMarkers(seurat.integrated, ident.1 = "covid", ident.2 = "normal", min.pct = 0.25)

# Sort by avg_log2FC to get the most upregulated genes
upregulated_genes <- de.genes[order(de.genes$avg_log2FC, decreasing = TRUE), ]
top_20_upregulated <- head(upregulated_genes, n = 20)

# Sort by avg_log2FC to get the most downregulated genes
downregulated_genes <- de.genes[order(de.genes$avg_log2FC), ]
top_20_downregulated <- head(downregulated_genes, n = 20)

##### Save both up and down regulated in both conditions####
# Check if the directory "results" exists, if not, create it
if (!dir.exists("results")) {
  dir.create("results")
}

write.csv(top_20_upregulated, file = "results/top_20_upregulated_genes.csv")
write.csv(top_20_downregulated, file = "results/top_20_downregulated_genes.csv")
# Then export for functional enrichment in DAVID


#### Cell annotation using Azimuth (Human - Lung reference dataset ####
# Extract normalized count data
counts_2 <- seurat.integrated@assays$RNA@counts
#counts <- seurat.integrated@assays$RNA@counts

# Load the reference dataset
ref <- HumanPrimaryCellAtlasData()
# Run SingleR
predictions_2 <- SingleR(test = counts_2, ref = ref, labels = ref$label.main)
#predictions <- SingleR(test = counts, ref = ref, labels = ref$label.main)

# Add SingleR predictions to Seurat object metadata
seurat.integrated$SingleR <- predictions_2$labels
#seurat.integrated$SingleR <- predictions$labels

# Visualize the cell clusters and their annotations
DimPlot(seurat.integrated, group.by = "SingleR", label = T)

# Get the UMAP plot data
umap.data <- as.data.frame(seurat.integrated@reductions$umap@cell.embeddings)
umap.data$Cluster <- seurat.integrated$SingleR

# Generate the plot
ggplot(umap.data, aes(UMAP_1, UMAP_2, color = Cluster)) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = Cluster), size = 3) +
  theme_classic() +
  theme(legend.position = "bottom")


#### Find the cells where SFTPB is highly expressed####
# Calculate the 90th percentile of SFTPB expression
percentile90 <- quantile(seurat.integrated[['RNA']]@data["SFTPB", ], 0.90)

# Use this as the threshold for high expression
cells <- WhichCells(seurat.integrated, expression = SFTPB > percentile90)

# Proceed as before
cell.types <- seurat.integrated$SingleR[cells]
table(cell.types)

 
pie(table(cell.types), labels = names(table(cell.types)))

# Create a new metadata column to mark the cells
seurat.integrated$HighSFTPB <- "No"
seurat.integrated$HighSFTPB[cells] <- "Yes"

# Plot on a UMAP plot
DimPlot(seurat.integrated, group.by = "HighSFTPB")

de.genes <- FindMarkers(seurat.integrated, ident.1 = "Yes", ident.2 = "No", min.pct = 0.25, group.by = "HighSFTPB")



####Heat map####
# Perform differential expression analysis
de.genes <- FindMarkers(seurat.integrated, ident.1 = "covid", ident.2 = "normal", min.pct = 0.25)

 
















