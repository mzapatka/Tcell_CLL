library(spatstat.utils)
library(Seurat)
library(SeuratObject)
library("devtools")
library(withr)
library(ggplot2)
library(viridis)
library("harmony")
library(qs)


#Sample merging and integration using harmony

filenames <- list.files(pattern="*.qs", full.names=TRUE)
ldf <- lapply(filenames, qread)
qs_names <- sub('.qs$', '', basename(filenames))


for (i in seq_along(ldf)) {
  DefaultAssay(ldf[[i]]) <- "RNA" 
  ldf[[i]]$sample_name <- qs_names[i]
  ldf[[i]][["orig.ident"]] <- qs_names[i]
}

scDataMerged <- merge(x = ldf[[1]], 
                      y = ldf[2:length(ldf)], 
                      add.cell.ids = (qs_names), 
                      merge.data = TRUE)
DefaultAssay(scDataMerged) = "RNA"
scDataMerged=NormalizeData(scDataMerged)
scDataMerged=ScaleData(scDataMerged)
scDataMerged <- FindVariableFeatures(scDataMerged, selection.method = "vst", nfeatures = 2000)
scDataMerged <- RunPCA(scDataMerged, verbose = TRUE, features = VariableFeatures(object = scDataMerged))
scDataMerged_harmony <- RunHarmony(scDataMerged,c("orig.ident"),plot_convergence = TRUE,theta=c(3))
ClustersTarget=16
scDataMerged_harmony <- RunUMAP(scDataMerged_harmony, reduction = "harmony", dims = 1:ClustersTarget)
scDataMerged_harmony <- FindNeighbors(scDataMerged_harmony,reduction = "harmony", dims = 1:ClustersTarget)
scDataMerged_harmony <- FindClusters(scDataMerged_harmony, resolution = (0.8))

pdf(file = "CLL_LN_UMAP_16CLres0.8_27112024.pdf",   # The directory you want to save the file in
    width = 10, height = 6) 
DimPlot(scDataMerged_harmony, reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE, pt.size = 0.5, raster = FALSE)
dev.off()

# Set identity classes to res.0.8
Idents(scDataMerged_harmony) <- "RNA_snn_res.0.8"

qsave(scDataMerged_harmony, "scDataMerged_harmony_23102024.qs")

# Find markers for every cluster compared to all remaining cells, report only the positive ones
scDataMerged_harmony.markers <- FindAllMarkers(scDataMerged_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

write.csv(scDataMerged_harmony.markers, "CLL_DEG_23102024.csv")


### Dot plot for marker genes per cluster (Figure 4B)
scDataMerged_harmony@active.ident <- factor(scDataMerged_harmony@active.ident, 
                                            levels=c("2",  "11", "1",   "4", "6",  "3",  "5", "9",
                                                     "0",   "14", "15", "8", "7", "10",  "12", "13" ))

features <- c("CD3E", "CD4", "CD8A", 	 "LEF1",	"CCR7",	"SELL", "IL7R","LTB", "CD40LG","CD69",
              "GNLY", "FCGR3A", "NKG7", "CCL5",	"PRF1", "GZMA", "GZMK",	"EOMES", 		
              "HAVCR2", "LAG3","ENTPD1", 	"PDCD1",	"CTLA4",	"TIGIT","TOX","MKI67", "TOP2A",
              "CXCR5", 	"CD200",	 	"ICOS",	"FOXP3",	"IL2RA"	, "IKZF2",  "SLC4A10", "TRAV1-2",
              "CD19",	"CD79A")
pdf(file = "CLL_LN_UMAP_DotPlot_271120224.pdf",   # The directory you want to save the file in
    width = 10, height = 5) 
DotPlot(scDataMerged_harmony, features = features) + RotatedAxis() + aes(colour = z1) +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red')
dev.off()


pdf(file = "CLL_LN_UMAP_16CLres0.8_orig.ident.pdf",   # The directory you want to save the file in
    width = 10, height = 6) 
DimPlot(scDataMerged_harmony, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5, raster = FALSE)
dev.off()

pdf(file = "CLL_LN_UMAP_PBMCpredicted.celltype.l2.pdf",   # The directory you want to save the file in
    width = 12, height = 6) 
DimPlot(scDataMerged_harmony, reduction = "umap", group.by = "PBMCpredicted.celltype.l2", label = TRUE, pt.size = 0.5, raster = FALSE)
dev.off()
