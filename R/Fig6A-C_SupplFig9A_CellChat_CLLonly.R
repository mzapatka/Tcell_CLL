library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

# Load CLL data set
ptm = Sys.time()
CLL_seurat_harmony_umap_obj <- readRDS("CLL_seurat_harmony_umap_obj.rds")

new.cluster.ids <- c("CD4 TN", "CD4 TH KLF2", "CD8 TPEX", "CD8 TEM", "CD4 TH CD69", "CD4 TREG",
                     "CD4 TFH", "CD8 TN", "CD8 TEX", "CLL", "CD4 TEM GZMK", "TPR", "DC")
names(new.cluster.ids) <- levels(CLL_seurat_harmony_umap_obj)
CLL_seurat_harmony_umap_obj <- RenameIdents(CLL_seurat_harmony_umap_obj, new.cluster.ids)
CLL_seurat_harmony_umap_obj@meta.data$seurat_clusters <- CLL_seurat_harmony_umap_obj@meta.data$new.cluster.ids

# Create a CellChat object
cellchat <- createCellChat(object = CLL_seurat_harmony_umap_obj, group.by = "ident")

# Set the ligand-receptor interactions base
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Compute the communication probability and infer cellular communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

###### Figure 6A: Heatmap depicting the number of interactions between cell subsets in CLL LNs ######
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = NULL, color.heatmap = "YlGnBu")

##### Figure 6B: Scatter plot showing the dominant sender and receiver cell subsets #####
netAnalysis_signalingRole_scatter(cellchat)

##### Figure 6C: Heatmap depicting the list of significant ligand-receptor pairs between CLL cells and all the other cell subsets #####
# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# show all the significant interactions (L-R pairs) from CLL cells (defined by 'sources.use') to all T cell clusters (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 10, targets.use = c(1:13), remove.isolate = FALSE)

##### Suppl. Figure 9A: 
# show all the significant interactions (L-R pairs) from CLL cells (defined by 'sources.use') to all T cell clusters (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:13), remove.isolate = FALSE)


