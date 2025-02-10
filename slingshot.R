# Load necessary libraries
library(Seurat)
library(SingleCellExperiment)
library(Slingshot)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# Step 1: Convert Seurat Object to SingleCellExperiment and Perform Slingshot Analysis
sce <- as.SingleCellExperiment(OKCepi, assay = "RNA")
sce_slingshot <- slingshot(sce, 
                           reducedDim = "UMAP", 
                           clusterLabels = sce$min.celltype, 
                           start.clus = "Basal-C1-ALDH3A1",
                           approx_points = 150)

# Step 2: Add Pseudotime to Seurat Object
pseudotime <- as.data.frame(slingPseudotime(sce_slingshot))
colnames(pseudotime) <- paste0("Lineage", 1:ncol(pseudotime))
OKCepi <- AddMetaData(OKCepi, metadata = pseudotime)

# Step 3: Pseudotime Density Plot Function
pseudotime_density <- function(seurat_obj, Lineage, cluster_label, colors) {
  df <- na.omit(data.frame(cluster = seurat_obj[[cluster_label]], pseudotime = seurat_obj[[Lineage]]))
  ggplot(df, aes(x = pseudotime, fill = cluster)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = colors)
}

# Generate density plots
library(dittoSeq)
p1 <- pseudotime_density(OKCepi, "Lineage1", "Tissue", dittoColors()) + theme(axis.title.x = element_blank())
p2 <- pseudotime_density(OKCepi, "Lineage1", "min.celltype", c("#EAAA60", "#E68B81", "#B7B2D0", "#7DA6C6", "#84C3B7"))

# Combine density plots
library(patchwork)
p1 + p2 + plot_layout(ncol = 1, heights = c(1, 1), guides = 'collect')

# Step 4: Fit GAM Model for Pseudotime Gene Expression Analysis
slingsce <- SlingshotDataSet(sce_slingshot)
pseudotime_vals <- slingPseudotime(slingsce, na = FALSE)
cellWeights_vals <- slingCurveWeights(slingsce)
counts_vals <- assays(sce_slingshot)$counts

sce_tradeSeq <- fitGAM(counts = counts_vals, 
                       pseudotime = pseudotime_vals, 
                       cellWeights = cellWeights_vals, 
                       nknots = 5, verbose = TRUE)

# Step 5: Differential Expression Analysis Along Pseudotime
assoc_res <- associationTest(sce_tradeSeq, lineages = TRUE, l2fc = log2(2))
gene_list <- rownames(assoc_res[assoc_res$pvalue < 0.05, ])  # Extract significant genes

# Step 6: Generate Pseudotime Heatmap Matrix
slingshot_for_plotMatrix <- function(seurat_obj, n_bins = 50, min_exp = 0.2) {
  seurat_meta <- seurat_obj@meta.data[order(seurat_obj$sling_pseudotime), ]
  expr_mat <- as.matrix(seurat_obj@assays$RNA@data[rownames(seurat_obj) %in% gene_list, ])
  expr_mat <- expr_mat[, order(match(colnames(expr_mat), rownames(seurat_meta)))]
  
  # Bin pseudotime and calculate average expression
  max_pseudotime <- max(seurat_meta$sling_pseudotime, na.rm = TRUE)
  pseudotime_bin_size <- max_pseudotime / n_bins
  binned_expr <- matrix(nrow = nrow(expr_mat), ncol = n_bins, dimnames = list(rownames(expr_mat), 1:n_bins))
  
  for (i in 1:n_bins) {
    bin_cells <- rownames(seurat_meta)[seurat_meta$sling_pseudotime > (i - 1) * pseudotime_bin_size & seurat_meta$sling_pseudotime <= i * pseudotime_bin_size]
    if (length(bin_cells) > 10) {
      binned_expr[, i] <- rowMeans(expr_mat[, colnames(expr_mat) %in% bin_cells, drop = FALSE], na.rm = TRUE)
    }
  }
  
  # Scale expression data and remove low-expressed genes
  scaled_expr <- t(scale(t(binned_expr)))
  scaled_expr[apply(abs(scaled_expr), 1, max, na.rm = TRUE) > min_exp, ]
}

# Generate heatmap matrix
heatmap_matrix <- slingshot_for_plotMatrix(OKCepi)

# Step 7: Pseudotime Heatmap Visualization
pheatmap(heatmap_matrix, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("#61AACF", "#EAEFF6", "#DA9599"))(250))

# Step 8: Functional Enrichment Analysis (GO & KEGG)
module_gene <- data.frame(Module = cutree(hclust(dist(heatmap_matrix)), k = 4), gene = rownames(heatmap_matrix))

# GO Enrichment
Module_GO <- do.call(rbind, lapply(unique(module_gene$Module), function(i) {
  genes <- filter(module_gene, Module == i)$gene
  df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  go_res <- enrichGO(gene = unique(df$ENTREZID), OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  go_res@result$cluster <- i
  go_res@result
}))

# KEGG Enrichment
Module_KEGG <- do.call(rbind, lapply(unique(module_gene$Module), function(i) {
  genes <- filter(module_gene, Module == i)$gene
  df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  kegg_res <- enrichKEGG(gene = unique(df$ENTREZID), organism = "hsa",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  kegg_res@result$cluster <- i
  kegg_res@result
}))

# Save enrichment results
write.csv(Module_GO, file = "Module_GO.csv")
write.csv(Module_KEGG, file = "Module_KEGG.csv")

# Step 9: Annotate and Visualize Genes in Heatmap
selected_genes <- c("MYC", "Jun", "KRT14", "TCF12", "AQP3", "HS6ST1", "RCOR1", "LIN7C", 
                    "CENPP", "SOD1", "FKBP1A", "CDCA7L", "TRIM29", "CENPN", "PTGES",
                    "HSPD1", "DTYMK", "SFR1", "H2AFZ", "UBXN4", "PFDN2", "YPLAL1", 
                    "ITGA2", "RACK1", "EEF1A1")

source('/home/guile/角化囊肿癌变/add.flag.R')
add.flag(pheatmap(heatmap_matrix), kept.labels = selected_genes, repel.degree = 0.2)
