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
library(circlize)

counts <- assay(daf2, "exprs")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(daf2) <- log2(t(t(counts)/size.factors) + 1)
assayNames(daf2)


manno.seurat <- as.Seurat(daf2, counts = "exprs", data = "logcounts")


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


##below is the cluster annotation in the order 1-30. add this as metadata in seurat object
merging_table2 <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/Figures/ClusterAnnotation.txt", header = T, sep = "\t")

head(data.frame(merging_table2)) 
merging_table2$new_cluster <- factor(merging_table2$new_cluster, 
                                     levels = c("CD4 TCM1", "CD4 TCM2 (low)", "CD4 TEM NCAM+", "CD4 TEM TBET+","CD4 TN",
                                                "CD4 TR1","CD4 TCM3","CD4 TEM CD39+ CD38+", "CD4 TEM CD38+ CTLA4+", "CD4 Treg1 Helios-",
                                                "CD4 rTregs2", "CD4 TFH PD-1hi","CD4 TEM PD-1+","Mixed ICOS+","CD8 TEX1 CD39+",
                                                "CD8 TN","DN EF1", "DN TEM1 KLRG1+", "CD8 TCM","CD4 Treg3 CD39+", 
                                                "Mixed TEX Proliferation",'DN TEF2 HELIOS+','DN TEM2 HELIOS+',"CD8 TEM GZMK+","CD8 TEM KLRG1 TBET",
                                                "CD8 TEX2","CD8 TEM TBET","CD8 TEF1","CD4 aTreg4","CD8 TEF2 NCAM+"))


####  UMAP plot by cluster
g<-DimPlot(manno.seurat, group.by='seurat_annotations', 
        cols=c("CD4 TCM1"="#0570b0","CD4 TCM2 (low)"="#74a9cf",
               "CD4 TEM NCAM+"="#FF9966","CD4 TEM TBET+"="#FFFF33","CD4 TN"="#0099FF",
               "CD4 TR1"="#CC0000","CD4 TCM3"="#a6bddb","CD4 TEM CD39+ CD38+"="#FF6633",
               "CD4 TEM CD38+ CTLA4+"= "#CC6633","CD4 Treg1 Helios-"= "#996666",
               "CD4 rTregs2"="#CC9999","CD4 TFH PD-1hi"= "#339900","CD4 TEM PD-1+"= "#CC6666",
               "Mixed ICOS"="#666699","CD8 TEX1 CD39+"="#9933FF",
               "CD8 TN"="#0066CC","DN EF1"="#666600","DN TEM1 KLRG1+"="#FF99CC",
               "CD8 TCM"= "#339999","CD4 Treg3 CD39+"="#999999",
               "Mixed TEX Proliferation"="#660000",'DN TEF2 HELIOS+'="#336600",
               'DN TEM2 HELIOS+'="#663366","CD8 TEM GZMK+"="#CC99FF","CD8 TEM KLRG1 TBET"="#CC6699",
               "CD8 TEX2"="#6600CC", "CD8 TEM TBET"="#CC0066","CD8 TEF1"="#99CC99",
               "CD4 aTreg4"="#333333","CD8 TEF2 NCAM+"="#66CC00"))  # + theme(legend.position = 'none')   ##PLOT WITH OR WITHOUT LEGEND

ggplot2::ggsave(filename = "umap_test.pdf",width=13,height=8, plot=g)
