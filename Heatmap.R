# Load required libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggsci)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

# Get mean gene expression for each cell group (Tissue)
mean_gene_exp <- AverageExpression(Mac,
                                   features = gene,
                                   group.by = 'Tissue',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

# Z-score normalization (mean centering)
htdf <- t(scale(t(mean_gene_exp), scale = FALSE, center = TRUE))

# Define color function for heatmap
col_fun = colorRamp2(c(-2, 0, 2), c( "#BCDCE4", "white", "#766DA7"))

# Define the genes of interest and their associated groups (M1 vs M2 signature)
gene <- c("IL1B", "CCL2", "PTGS2", "CD163", "MRC1", "CHI3L1")
gene_groups <- c("M1 signature", "M1 signature", "M1 signature", "M2 signature", "M2 signature", "M2 signature")

# Row annotation based on gene signature
row_ha <- rowAnnotation(
  Signature = factor(gene_groups, levels = c("M1 signature", "M2 signature")),
  col = list(Signature = c("M1 signature" = "#66C2A5", "M2 signature" = "#BBD2A0"))  # You can modify colors
)

# Define cluster colors for OKC and COKC groups
cluster_colors <- c("OKC" = "#C2ADCE", "COKC" = "#ECA8A9")

# Column annotation for clusters (OKC vs COKC)
column_ha <- HeatmapAnnotation(
  Cluster = factor(colnames(htdf), levels = names(cluster_colors)),  # Ensure correct matching of cluster colors
  col = list(Cluster = cluster_colors)
)

# Generate the heatmap
Heatmap(htdf,
        name = "Average Expression",
        cluster_columns = FALSE,  # Do not cluster columns
        cluster_rows = FALSE,  # Do not cluster rows
        row_names_gp = gpar(fontface = 'italic', fontsize = 10),  # Italicize and set font size for row names
        row_names_side = "right",  # Place gene names on the right side of the heatmap
        border = TRUE,
        rect_gp = gpar(col = "white", lwd = 1),  # Add borders between cells
        column_names_side = "top",  # Place column names at the top
        column_names_rot = 45,  # Rotate column names for readability
        top_annotation = column_ha,  # Add top annotation for cluster
        left_annotation = row_ha,  # Add left annotation for gene signature
        col = col_fun)  # Apply color scale defined earlier


###Plotting a Single Group Circular Heatmap

# If the matrix data is grouped, use the 'split' parameter to specify the grouping variable
gene_groups <- Gene$group  # Assigning group information to 'gene_groups'
names(gene_groups) <- rownames(cir1)  # Assigning the row names of 'cir1' as names for gene_groups

# Convert the gene groups into factors required for splitting the heatmap
split_factors <- factor(gene_groups)

# Define color scale for the heatmap
mycol1 = colorRamp2(c(-2, 0, 2), c("#F28E2B", "white", "#E15759"))

# Adjust circular plot parameters
circos.par(gap.after = c(2, 2, 2, 2, 50))  # Adjusts the space between the sectors in the circular plot
# Increasing the gap size allows for a wider separation between sectors (for visual clarity)

# Create the circular heatmap
circos.heatmap(cir1, 
               col = mycol1,  # Apply the color scale
               dend.side = "inside",  # Position the dendrogram (row clustering) inside the circular plot
               rownames.side = "outside",  # Position the row names outside the circle
               track.height = 0.22,  # Height of each track (thickness of the circle)
               rownames.col = "black",  # Color for row names
               bg.border = "black",  # Background border color
               split = split_factors,  # Split the heatmap by row annotation
               show.sector.labels = TRUE,  # Show labels for each sector (group)
               rownames.cex = 0.9,  # Adjust font size of row names
               rownames.font = 1,  # Font style for row names (1 = normal)
               cluster = TRUE,  # Perform hierarchical clustering for rows (TRUE = clustering, FALSE = no clustering)
               dend.track.height = 0.18,  # Adjust the height of the dendrogram track
               dend.callback = function(dend, m, si) {  # Callback for dendrogram customization
                 # Modify the dendrogram color branches based on cluster
                 color_branches(dend, k = 5, col = 1:5)  # Color the branches in 5 different colors
               }
)

# Legend and column name settings
lg = Legend(title = "Legend", col_fun = mycol1, direction = c("horizontal"))  # Define legend with horizontal orientation
grid.draw(lg)  # Draw the legend

# Additional text labels for columns (at the end of the plot)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 3) {  # The last sector (for the final group)
    cn = colnames(cir1)  # Get the column names of the matrix
    n = length(cn)  # Number of columns
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),  # X coordinates for text
                (1:n) * 0.2 - 1.3,  # Y coordinates for text, adjust based on desired spacing
                cn, cex = 0.8, adj = c(0, 1), facing = "inside")  # Adjust text position and style
  }
}, bg.border = NA)

# Clear the circular plot after completion
circos.clear()


