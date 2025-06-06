
library(Seurat)
library(ggplot2)
library(viridis)

##### Suppl. Figure 9E: Violin plots showing the expression of LGALS9 and HAVCR2 in CLL LN and rLN (from Aoki et al.) #####

## CLL LNs
CLL_seurat_harmony_umap_obj <- readRDS("CLL_seurat_harmony_umap_obj.rds")

# Violin plot
VlnPlot(CLL_seurat_harmony_umap_obj, features = c("HAVCR2", "LGALS9"))

## rLN 
rLN_seurat_harmony_umap_obj <- readRDS("rLN_seurat_harmony_umap_obj.rds")

# Violin plot
VlnPlot(rLN_seurat_harmony_umap_obj, features = c("HAVCR2", "LGALS9"))
