# Load required libraries
library(CellChat)
library(ggplot2)

# Prepare data for CellChat
data.input <- cellchat[["RNA"]]@data  # Extract the RNA expression data from the CellChat object
meta <- cellchat@meta.data  # Extract the metadata (cell annotations)
cellchatC <- createCellChat(object = data.input, meta = meta, group.by = "celltype")  # Create a CellChat object based on the data and metadata

# Load the human CellChat database and display available categories
CellChatDB <- CellChatDB.human  # Human CellChat database
showDatabaseCategory(CellChatDB)  # Display the categories in the database

# Subset the CellChat database to focus on secreted signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchatC@DB <- CellChatDB.use  # Assign the subset database to the CellChat object

# Subset data by selected features (set to NULL to keep all features)
cellchatC <- subsetData(cellchatC, features = NULL)

# Identify over-expressed genes and interactions
cellchatC <- identifyOverExpressedGenes(cellchatC)  # Identify genes that are over-expressed
cellchatC <- identifyOverExpressedInteractions(cellchatC)  # Identify interactions between genes that are over-expressed

# Project data onto the protein-protein interaction (PPI) network (mouse in this case)
cellchatC <- projectData(cellchatC, PPI.mouse)

# Compute the communication probability matrix
cellchatC <- computeCommunProb(cellchatC, raw.use = TRUE)  # Compute the raw communication probability

# Filter communication pathways with a minimum of 10 cells involved
cellchatC <- filterCommunication(cellchatC, min.cells = 10)

# Compute communication probability for individual pathways
cellchatC <- computeCommunProbPathway(cellchatC)

# Aggregate the network for visualization
cellchatC <- aggregateNet(cellchatC)

# Calculate the size of each group (number of cells per cluster)
groupSize <- as.numeric(table(cellchatC@idents))  # Create a numeric vector with group sizes

# Visualize the network as a circle plot, scaling edge weights by group size
netVisual_circle(cellchatC@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")

# Compute the centrality of the network (importance of each node)
cellchatC <- netAnalysis_computeCentrality(cellchatC, slot.name = "netP")

# 1. Generate a basic bubble plot of communication between specified sources and targets
bubble_plot <- netVisual_bubble(cellchatC, 
                                sources.use = c("Macrophage-APOE", "Macrophage-C1QA", "Macrophage-CD209", "Macrophage-IL1B"), 
                                targets.use = c("Endothelial"), 
                                remove.isolate = TRUE, angle.x = 45)  # Create a bubble plot for the interaction between macrophage subtypes and endothelial cells
