#  library(HDCytoData)
library('tidyverse')
library('readxl')
library('flowCore')
library('Cairo')
library("CATALYST")
library(corrplot)
library(ComplexHeatmap)
library(affycoretools)

daf <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/LiveCellsCytof.rds")


daf2 <- filterSCE(daf, cluster_id %in% c(1,2,3,4,5,6,7,8,10,
                                          13,14,15,22,
                                          20,21,23,26,27,28,30), k = 'meta30')


#########... median  expression _marker_cluster_meta30: T cells ######FIGURE 1B
#####new annotation is added here
merging_table2 <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/ClusterAnnotation.csv", header = T, sep = ",")

head(data.frame(merging_table2)) 

setwd("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/")
merging_table2$new_cluster <- factor(merging_table2$new_cluster, levels = c("CD4 Tnaive", "CD8 Tnaive", "CD4 Tcm", "CD4 Tcm Low","CD8 Tcm",
                                                                            "CD4 Tem","CD4 Tem TBET+","CD4 Tem NCAM+", "CD8 Tem KLRG1 TBET", "CD4 Tem CD38+","CD4 Tem CD39+ CD38+", "CD4 Tem PD-1+","CD4 Tex (TR1)","CD4 Tex Prolif",
                                                                            "CD8 Tem GZM+","CD8 Tex","CD8 Tex Prolif", "DN Tem CD39+", "DN ICOS+","DN Tem HELIOS+", "DP ICOS+",
                                                                            'CD4 Tef PD-1hi CXCR5','CD8 Tef',"DN Tef","DN Tef/HELIOS+",
                                                                            "CD4 rTregs","CD4 aTregs","CD4 Treg HELIOS-",
                                                                            "CD4 TFH","CD4 CD45RO+ CD45RA+"))
dafN <- mergeClusters(daf2, k = "meta30",
                      table = merging_table2, id = "merging4")
colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(200)
hm5 <- read.table("med_exprsMarkersPerClusterDaf2", header = T, sep = "\t", row.names = 1)
normalized = (hm5-min(hm5))/(max(hm5)-min(hm5))
normalized <- as.matrix(normalized)
Heatmap(normalized,
        cluster_rows = TRUE,column_names_gp = gpar(fontsize = 9),
        cluster_columns = TRUE, row_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = TRUE,col = colors,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson"
)


daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsDaf2.rds")
                      
freqs2 <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/Daf2_Tcells_Table_perSampleClusterCells", header = T, sep = "\t", row.names = 1)
freqs22 <- freqs2/apply(freqs2,1,sum)

#remove samples
freqs22 = freqs22[which(rownames(freqs22) != "HD3LN"),]

#pch annotation list
pchlist = rownames(freqs22)
pchlist[which(grepl("^RLN",pchlist))] = 17
pchlist[which(grepl("PB$",pchlist))] = 16
pchlist[which(grepl("BM$",pchlist))] = 15
pchlist[which(grepl("LN$",pchlist))] = 18
pchlist[which(grepl("LN.$",pchlist))] = 18
pchlist = as.numeric(pchlist)

pdf("cell_types_CyToF_PCA.pdf",width = 5,height=5)

affycoretools::plotPCA(t(freqs22),addtext=rownames(freqs22),pch=pchlist,col=pchlist)
title(xlab="                             (22.9%)",ylab="                             (12.4%)")
#,pch=20,labels=rownames(freqs22))
dev.off()
pdf("cell_types_CyToF_PCA_nolabel.pdf",width = 5,height=5)

affycoretools::plotPCA(t(freqs22),pch=pchlist,col=pchlist,cex=1.5)
legend(-0.4, 0.3, legend=c("RLN", "PB","BM","LN"),
       col=c("black", "gray", "gold", "#e57588"),
       pch=c(17,16,15,18) )
title(xlab="                             (18.32%)",ylab="                             (31.23%)")
#,pch=20,labels=rownames(freqs22))
dev.off()

