library(HDCytoData)
library('tidyverse')
library('readxl')
library('flowCore')
library('Cairo')
library("CATALYST")
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(circlize)
library("RColorBrewer")

####  load object for T cells for all the plots below
daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/TcellsDaf2.rds")

##############   Abundance heatmap/barplot with pathological features like age, tumor load etc  ######################################

library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")



freqs2 <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/Figures/DafN_Tcells_Table_perSampleClusterCells", header = T, sep = "\t", row.names = 1)

freqsN <- freqs2[c("CD4.TN", "CD8.TN", "CD4.TCM1", "CD4.TCM2..low.", "CD4.TCM3", 
                   "CD8.TCM", "CD4.TEM.TBET.", "CD4.TEM.NCAM.", "CD4.TEM.CD38..CTLA4.", "CD4.TEM.CD39..CD38.", 
                   "CD4.TEM.PD.1.",  "CD4.TR1", "Mixed.TEX.Proliferation", "CD8.TEM.KLRG1.TBET", "CD8.TEM.TBET", 
                   "CD8.TEM.GZMK.", "CD8.TEX1.CD39.", "CD8.TEX2", "DN.TEM1.KLRG1.", 'DN.TEM2.HELIOS.',
                   "Mixed.ICOS.", "CD8.TEF1", "CD8.TEF2.NCAM.", "DN.EF1", 'DN.TEF2.HELIOS.',
                   "CD4.TFH.PD.1hi", "CD4.Treg1.Helios.",  "CD4.rTregs2",  "CD4.Treg3.CD39.", "CD4.aTreg4")]
my_color_palette <- c("#0099FF","#0066CC","#0570b0","#74a9cf","#a6bddb",
                      "#339999","#FFFF33","#FF9966", "#CC6633", "#FF6633",
                      "#CC6666", "#CC0000", "#660000", "#CC6699","#CC0066",
                      "#CC99FF","#9933FF", "#6600CC","#FF99CC","#663366",
                      "#666699","#99CC99","#66CC00","#666600","#336600",
                      "#339900", "#996666","#CC9999","#999999","#333333")
freqs22 <- freqsN/apply(freqsN,1,sum)
freqs22$Total <- NULL

hc <- hclust(dist(freqs22),method="ward.D")
meta <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/Figures/files", header = T, sep = '\t')

haB <- rowAnnotation(Condition = meta$type,
                     TumorSite = meta$TumorSite,
                     IGHVstatus = meta$IGHVstatus,
                     Age = meta$age,
                     PercentTL = meta$TumorLoad,
                     col=list(Condition = c("tumor"="indianred2","control"="olivedrab3"),
                              TumorSite = c("LymphNode"="orange","PeripheralBlood"="grey","BoneMarrow"="pink"),
                              IGHVstatus = c("mutated"="coral1","unmutated"="mediumpurple1","NA"='grey'),
                              PercentTL = colorRamp2(c(30, 60, 90), c("lightcyan2", "lightblue1", "lightblue3")),
                              Age = colorRamp2(c(0,50 ,100), c("white", "pink2","pink4"))),
                     width = unit(4,"cm"),
                     show_legend =TRUE,
                     rn = anno_text(rownames(freqs2),gp = gpar(fontsize = 7, fontface="bold")),
                     annotation_name_gp = gpar(fontsize = 10, fontface="bold"))

ha_list2 = rowAnnotation(ClusterAbundance = anno_barplot(freqs22, bar_width = 1.0,
                                                         width = unit(6, "cm"),
                                                         height = unit(10, "cm"),gp=gpar(fill=my_color_palette),
                                                         which = "row",show_legend = TRUE,
                                                         rn = anno_text(colnames(freqs22)),
                                                         annotation_name_gp = gpar(fontsize = 10, fontface="bold")))
lgd_boxplot = Legend(labels = c("CD4 TN", "CD8 TN", "CD4 TCM1", "CD4 TCM2 (low)", "CD4 TCM3",
                                "CD8 TCM", "CD4 TEM TBET+", "CD4 TEM NCAM+", "CD4 TEM CD38+ CTLA4+", "CD4 TEM CD39+ CD38+",
                                "CD4 TEM PD-1+",  "CD4 TR1", "Mixed TEX Proliferation", "CD8 TEM KLRG1 TBET", "CD8 TEM TBET", 
                                "CD8 TEM GZMK+", "CD8 TEX1 CD39+", "CD8 TEX2", "DN TEM1 KLRG1+", 'DN TEM2 HELIOS+',
                                "Mixed ICOS+", "CD8 TEF1", "CD8 TEF2 NCAM+", "DN EF1",  'DN TEF2 HELIOS+', "CD4 TFH PD-1hi",
                                "CD4 Treg1 Helios-",  "CD4 rTregs2",  "CD4 Treg3 CD39+", "CD4 aTreg4"), title = "Clusters",
                     legend_gp = gpar(fill = c("#0099FF","#0066CC","#0570b0","#74a9cf","#a6bddb",
                                               "#339999","#FFFF33","#FF9966", "#CC6633", "#FF6633",
                                               "#CC6666", "#CC0000", "#660000", "#CC6699","#CC0066",
                                               "#CC99FF","#9933FF", "#6600CC","#FF99CC","#663366",
                                               "#666699","#99CC99","#66CC00","#666600","#336600",
                                               "#339900", "#996666","#CC9999","#999999","#333333")))

pdf("CyToF/cell_types_v2_CyToF_HM.pdf",width = 7,height=7)
ss <- Heatmap(matrix(nc = 0, nr = 45), cluster_rows = hc, 
              left_annotation = ha_list2, show_row_dend = T, 
              right_annotation = haB,show_heatmap_legend = TRUE) 

draw(ss,heatmap_legend_list = list(lgd_boxplot))
dev.off()