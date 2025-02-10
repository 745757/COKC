# Load required libraries
library(Seurat)  # For single-cell RNA-seq analysis
library(DoubletFinder)  # For doublet detection
library(Harmony)  # For batch effect correction
library(ggplot2)  # For plotting
library(ggsignif)  # For statistical tests in plots
library(ggdist)  # For raincloud plots
library(clustree)
# Preprocessing of scRNA-seq Data
sample_dirs <- list.dirs(data_dir, recursive = FALSE)
seurat_list <- lapply(sample_dirs, function(sample_dir) {
  sample_name <- basename(sample_dir) 
  sample_data <- Read10X(data.dir = sample_dir)
  seurat_obj <- CreateSeuratObject(counts = sample_data,min.cells = 3, min.features = 200)
  seurat_obj$patient <- sample_name
  return(seurat_obj)
})
OKC <- Reduce(function(x, y) merge(x, y), seurat_list)
# Add mitochondrial and ribosomal gene percentage to metadata
OKC[["percent.mt"]] <- PercentageFeatureSet(OKC, pattern = "^MT-")  
OKC <- subset(OKC, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
# Identify doublets using DoubletFinder (adjust the parameters as needed)
OKC <- DoubletFinder::doubletFinder_v3(OKC, pN = 0.25, pK = 0.05, nExp = 500, reuse.pANN = FALSE)
OKC <- OKC[,OKC$doublet_info%in%c("Singlet")]
# Normalize the data using LogNormalize method
OKC <- NormalizeData(object = OKC, normalization.method = "LogNormalize", scale.factor = 10000)
# Identify the 2000 most variable genes
OKC <- FindVariableFeatures(OKC, selection.method = "vst", nfeatures = 2000)
# Perform cell cycle scoring using the predefined S and G2M genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
OKC <- CellCycleScoring(OKC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regress out unwanted sources of variation (cell cycle scores, RNA count, mitochondrial and ribosomal genes)
OKC <- ScaleData(OKC, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"), features = rownames(OKC))
# Integration of datasets from multiple samples using Harmony to correct batch effects
OKC <- RunHarmony(OKC, group.by.vars = "orig.ident")
# Run PCA for dimensionality reduction
OKC <- RunPCA(OKC, features = VariableFeatures(object = OKC), verbose = FALSE)
# Visualize elbow plot to decide the number of principal components to use
ElbowPlot(OKC, 50)
# Run UMAP for visualization in two dimensions (using 1:25 principal components)
OKC <- RunUMAP(OKC, reduction = "harmony", dims = 1:25, min.dist = 1)
# Perform clustering based on the first 25 principal components and find clusters
OKC <- FindNeighbors(OKC, reduction = "harmony", dims = 1:25) %>% FindClusters(resolution = seq(0.1, 1.6, 0.1))
# Visualize the clustering tree
clustree(OKC@meta.data, prefix = "RNA_snn_res.")
# Final clustering with resolution = 0.4
OKC <- FindClusters(OKC, resolution = 0.4)
# Visualize the UMAP with custom cell type colors and labels
DimPlot(OKC, label = FALSE, cols = cell_type_cols, label.size = 6) + 
  labs(x = "UMAP1", y = "UMAP2") +   
  theme(axis.text.y = element_blank(),   
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank()) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))
