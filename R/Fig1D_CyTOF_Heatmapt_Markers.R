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

manno.seurat <- as.Seurat(daf2, counts = "exprs", data = "logcounts")

manno.seurat <- readRDS('CyTOFfcsfileslivecells/CYTOFanalysis_REFER/mannoSeurat.rds')

manno.seurat@meta.data$seurat_annotations='CD8 TEF2 NCAM+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('1','2','3','31','39',
                                                                                   '11','12','13','30','40',
                                                                                   '42','14','18','21','29')]='CD4 TCM1'

manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('4')]='CD4 TCM2 (low)'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('5')]='CD4 TEM NCAM+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('6','15','16')]='CD4 TEM TBET+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('7','8','9','10',
                                                                                   '19','20')]='CD4 TN'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('17','45','54','55',
                                                                                   '64','74')]='CD4 TR1'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('22','23','24','25','33','35',
                                                                                   '44')]='CD4 TCM3'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('26','27','36','37','38','46',
                                                                                   '47')]='CD4 TEM CD39+ CD38+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('28')]='CD4 TEM CD38+ CTLA4+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('32','71')]='CD4 Treg1 Helios-'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('34','83','92',
                                                                                   '93','94')]='CD4 rTregs2'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('41','51','53',
                                                                                   '61')]='CD4 TFH PD-1hi'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('43','52')]='CD4 TEM PD-1+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('48')]='Mixed ICOS+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('49','84')]='CD8 TEX1 CD39+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('50','59','60','70',
                                                                                   '80')]='CD8 TN'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('56','73')]='DN EF1'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('57')]='DN TEM1 KLRG1+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('58','68')]='CD8 TCM'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('62','72')]='CD4 Treg3 CD39+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('63')]='Mixed TEX Proliferation'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('65','100')]='DN TEF2 HELIOS+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('66')]='DN TEM2 HELIOS+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('67','76','77','87')]='CD8 TEM GZMK+'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('69','86','96')]='CD8 TEM KLRG1 TBET'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('75','85','95')]='CD8 TEX2'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('78','88','89','97','98')]='CD8 TEM TBET'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('79','90')]='CD8 TEF1'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('81','82','91')]='CD4 aTreg4'
manno.seurat@meta.data$seurat_annotations[manno.seurat@meta.data$cluster_id %in% c('99')]='CD8 TEF2 NCAM+'



dot_genes = rev(c("CD4","CD8a","CD45RA","CCR7","TCF1","IL7Ra","CD45RO","CD25","TBET","KLRG1","NCAM","ICOS","GZMK","TOX","TIGIT","PD1","CD38","CD39","CXCR5","CTLA4","HELIOS","FOXP3","Ki67"))
identsuse =  c(
  "CD4 TN",
  "CD8 TN", 
  "CD4 TCM1", 
  "CD4 TCM2 (low)",
  "CD4 TCM3",  
  "CD8 TCM", 
  "CD4 TEM TBET+", 
  "CD4 TEM NCAM+", 
  "Mixed ICOS+", 
  "CD8 TEM GZMK+", 
  "CD8 TEX2", 
  "CD8 TEX1 CD39+",
  "CD8 TEM KLRG1 TBET",
  'DN TEF2 HELIOS+', 
  "DN TEM1 KLRG1+", 
  "CD8 TEM TBET", 
  "CD8 TEF2 NCAM+",
  "CD8 TEF1", 
  "CD4 TEM CD39+ CD38+", 
  "CD4 TR1", 
  "CD4 TEM PD-1+", 
  "CD4 TFH PD-1hi", 
  "CD4 TEM CD38+ CTLA4+", 
  "CD4 Treg1 Helios-",
  "CD4 rTregs2", 
  "CD4 Treg3 CD39+",  
  "CD4 aTreg4", 
  'DN TEM2 HELIOS+',
  "DN EF1", 
  "Mixed TEX Proliferation")
manno.seurat$seurat_annotations = factor(manno.seurat$seurat_annotations, levels=identsuse)

#gene gates decision
pdf(paste(sep="","dotplot_Ridge_GateSelector.pdf"),width=13,height=20)
RidgePlot(manno.seurat, features = dot_genes, ncol = 2)
dev.off()

DefaultAssay(manno.seurat) = "RNA"
manno.seurat = NormalizeData(manno.seurat,normalization.method = 'CLR')

manno.seurat@assays['RNAGated'] = manno.seurat@assays['RNA']
DefaultAssay(manno.seurat) = "RNAGated"

#Gating for CyTof intensity data
manno.seurat@assays$RNAGated@data[which(manno.seurat@assays$RNAGated@data < 0.5)] = 0

#col.min
g <- DotPlot(object = manno.seurat, features = dot_genes, cols = "RdBu",
             group.by = 'seurat_annotations',cluster.idents =FALSE, scale=TRUE) + 
  coord_flip() +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

g

ggplot2::ggsave(filename = paste(sep="","dotplot_ctyof_normalized_CLR_gated.pdf"),width=12,height=6.7, plot=g)

