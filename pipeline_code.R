library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)


library(data.table)
library(Seurat)

# Set working directory to where CSVs are extracted
setwd("D:/cell-type-annotation-from-scRNAseq-brain-dataset")
library(data.table)


files <- list.files(pattern = "*.csv.gz")

# Initialize list to hold all cell columns
all_expr <- list()
gene_names <- NULL

for (f in files) {
  df <- fread(f)
  sample_id <- sub(".csv.gz", "", f)
  
  # Store gene names only once
  if (is.null(gene_names)) {
    gene_names <- df[[1]]
  }
  
  # Extract expression vector
  expr_vec <- df[[2]]
  all_expr[[sample_id]] <- expr_vec
}

# Combine all columns into a matrix
expr_matrix <- do.call(cbind, all_expr)
rownames(expr_matrix) <- gene_names

# Replace NA with 0
expr_matrix[is.na(expr_matrix)] <- 0
storage.mode(expr_matrix) <- "numeric"

# Create Seurat object
library(Seurat)
hippo <- CreateSeuratObject(counts = expr_matrix, project = "Hippocampus")
hippo <- NormalizeData(hippo)
hippo <- FindVariableFeatures(hippo)
hippo <- ScaleData(hippo)
hippo <- RunPCA(hippo)
hippo <- FindNeighbors(hippo, dims = 1:10)
hippo <- FindClusters(hippo, resolution = 0.5)
hippo <- RunUMAP(hippo, dims = 1:10)
Idents(hippo) <- "seurat_clusters"
table(Idents(hippo))  # Should show cluster counts (e.g., 0, 1, 2, ...)
FeaturePlot(hippo, features = c("SNAP25", "GFAP", "CX3CR1", "MBP", "PDGFRA", "CLDN5"))

library(plyr)

hippo$celltype <- plyr::mapvalues(
  hippo$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6"),
  to   = c("Oligodendrocyte", "Astrocyte", "Neuron", "OPC", "Microglia", "Unknown", "Endothelial")
)
DimPlot(hippo, reduction = "umap", group.by = "celltype", label = TRUE) +
  ggtitle("UMAP with Annotated Cell Types")
