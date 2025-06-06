library(HDCytoData)
library('tidyverse')
library('readxl')
library('flowCore')
library('Cairo')
library("CATALYST")
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(dichromat)
library(RColorBrewer)
library(scales)
library(ggplot2)

####  load object for T cells for all the plots below
daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/TcellsDaf2.rds")

###converted this sce to seurat object ### for better visualization of the UMAP
counts <- assay(daf2, "exprs")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(daf2) <- log2(t(t(counts)/size.factors) + 1)
assayNames(daf2)

manno.seurat <- as.Seurat(daf2, counts = "exprs", data = "logcounts")


##############################################   GENE EXPRESSION ON UMAP      ###############################################################################


######   this needs to be visualized as one image after rasterizing????  ########  (remaining work)  ###################################################

pdf("cytofALLMarkers_TcellsDaf2.pdf",useDingbats = FALSE, width = 13)

plotDR(daf2, "UMAP", color_by = "TBET", a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "FOXP3", a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "TCF1",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD45RO",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD44",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CCR7",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD73",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CTLA4",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "TIGIT",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD27",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "HELIOS",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "x2B4",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "KLRG1",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD19",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "EOMES",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "GZMK",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD45RA",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD38",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD4",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD8a",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "ICOS",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "TOX",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD3",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD7",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "IL7Ra",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "FAS",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "TIM3",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD25",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD45",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CXCR5",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "Ki67",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "x4_1BB",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "PD1",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "NCAM",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
plotDR(daf2, "UMAP", color_by = "CD39",a_pal= c('navyblue','darkgreen','yellow','firebrick4')) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


