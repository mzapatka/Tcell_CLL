#library(HDCytoData)
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


####  load object for T cells for all the plots below
daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/TcellsDaf2.rds")

############################################### UMAP to show sample distribution   ###############################################################

###converted this sce to seurat object ### for better visualization of the UMAP

counts <- assay(daf2, "exprs")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(daf2) <- log2(t(t(counts)/size.factors) + 1)
assayNames(daf2)

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




manno.seurat@meta.data$sample_id = factor(manno.seurat@meta.data$sample_id,
                                          levels=c('RLN1','RLN2','RLN3','RLN4','RLN5','RLN6','RLN7','RLN8','RLN9','RLN10','RLN11',
                                                   'RLN12','RLN13',"BC1LNA","BC2LNA", "BC3LN","BC4LN", "BC5LN1","BC5LN2", 
                                                   "BC8LN","BC9LN", "BC12LN","BC13LN", "BC14LN","BC15LN1", "BC15LN2",
                                                   "HD1LN","HD2LN","HD4LN","HD5LN","HD6LN", "HD8LN","HD9LN","HD10LN", 
                                                   "HD11LN",'BC1PB','BC05PB','BC8PB','BC09PB','BC12PB', 'BC14PB',
                                                   'BC15PB', 'BC3BM','BC09BM','BC11BM'))

##umap plot by sample
DimPlot(manno.seurat, group.by='sample_id',
        cols = c(
          "HD8LN"='#333399',"HD9LN"='#6699CC',   "BC12LN"='#3399FF',  "BC5LN2"='#6666FF', 
          'RLN12'='#FF9966','RLN13'='#FF6633', 'BC3BM'='#99FF66','BC09BM'='#66FF33','BC11BM'='#00FF66',
          'BC12PB'='#CC99FF','BC05PB'='#9999FF','BC09PB'='#9966CC','BC14PB'='#663399',
          'BC15PB'='#9933FF','BC1PB'='#6600CC','BC8PB'='#9900CC',
          "BC13LN" ='#3333FF',
          "BC14LN"='#6666CC',  "BC15LN1"='#99CCFF', "BC15LN2"='#0066CC', "HD10LN"='#3333CC',  "HD11LN" ='#0000CC',
          'RLN1'='#FF9999','RLN2'='#FF6666','RLN3'='#FF3333',
          "BC1LNA"='#0099CC',"BC2LNA"='#33CCFF',"BC3LN"='#0099FF',"BC4LN"='#0066FF',
          'RLN4'='#FF0000','RLN5'='#CC0000','RLN6'='#CC3333','RLN7'='#CC6666',
          "BC5LN1"='#0033FF',"BC8LN"='#0033CC',
          "BC9LN"='#3366FF',"HD1LN"='#6699FF',"HD2LN"='#3399CC',"HD4LN"='#0066CC',
          'RLN8'='#FFCCCC','RLN9'='#FF3366','RLN10'='#FF0033','RLN11'='#996666',
          "HD5LN"='#006699',"HD6LN"='#003399'))    # + theme(legend.position = 'none')   ##PLOT WITH OR WITHOUT LEGEND
