
setwd("/CYTOFanalysis_REFER")

library(HDCytoData)
library('tidyverse')
library('readxl')
library('flowCore')
library('Cairo')
library("CATALYST")
library(corrplot)
library(ComplexHeatmap)



####### Read in the files
#######################################
mdH <- read.table("files", header = T, sep = ",")
mdH <- as.data.frame(mdH)
mdH$file_name <- as.character(mdH$file_name)
mdH$sample_id <- as.character(mdH$sample_id)
mdH$patient_id <- as.character(mdH$patient_id)
mdH$TumorSite <- as.factor(mdH$TumorSite)
mdH$condition <- as.character(mdH$condition)
mdH$treatment <- as.character(mdH$treatment)
mdH$IGHVstatus <- as.character(mdH$IGHVstatus)
mdH$type <- as.character(mdH$type)

head(data.frame(mdH))  

##load .fcs files 
##this steps combines the intensitities across all samples into one combined exprs matrix.
##this by default transforms marker intensities and removes cells with extreme positive values
fcs_raw <- read.flowSet(mdH$file_name, transformation = FALSE,truncate_max_range = FALSE,
                        alter.names=FALSE, as.is = T, invert.pattern = FALSE,ignore.text.offset = TRUE)
                        


fcs_raw


####now loading isotope/marker information
panel <- read.csv("panelNEW.csv", header = TRUE, sep = ",")
head(data.frame(panel)) 

panel <- as_tibble(panel)
panel$fcs_colname <- as.character(panel$fcs_colname)
panel$antigen <- as.character(panel$antigen)
panel$marker_class <- as.character(panel$marker_class) #state or type 



#####this checks if column name in .fcs file matches the row names in panel file
###both these files should have same number of markers and the code line should return TRUE
all(panel$fcs_colname %in% colnames(fcs_raw))


####reordering levels as factor
mdH$sample_id <- factor(mdH$sample_id, levels = mdH$sample_id[order(mdH$condition)])    


##preparing the dataObject that will contain exprs matrix, metaData, transformed values and other parameters
###used during plotting and downstream analysis
daf <- CATALYST::prepData(fcs_raw, panel, mdH, features = panel$fcs_colname)

##provinding my own colors
condition_colors = c("tumorAccLN"="red","tumorLN"="cyan","tumorPB"="cornflowerblue","tumorBM"="darkcyan","tumorHodgkinLN"="yellow", "controlLN" = "violet")

setwd("CyTOFfcsfileslivecells/CYTOFanalysis_REFER")
pdf("Counts.pdf", useDingbats = FALSE, width = 10)
plotCounts(daf, color_by = "condition") + scale_color_manual(values = condition_colors)
dev.off()

write.table(daf@metadata$experiment_info, file="sampleCounts")

pdf("SampleMDS.pdf", useDingbats = FALSE, width = 10)
plotMDS(daf, color_by = "condition") + scale_color_manual(values = condition_colors)
dev.off()

ss <- plotMDS(daf, color_by = "condition")  ###daf only for creating this MDS with just points
ss$layers[[1]] <- NULL
pdf("plotMDS_points.pdf", useDingbats = FALSE, width = 7, height = 5)
ss + geom_point(size = 5) + theme(panel.background = element_blank())

dev.off()
#


##stateMarkersClustering
set.seed(1342)
daf <- cluster(daf, features = state_markers(daf), 
               xdim = 10, ydim = 10, maxK = 30, seed = 1234)   



daf <- runDR(daf, dr = "UMAP", features = "state", cells=4000) 

pdf("clusterAllCellsUMAP.pdf",useDingbats = FALSE, width = 13)
plotDR(daf, "UMAP", color_by = "meta30")
dev.off()



pdf("cluster4kCellsmarkers.pdf",useDingbats = FALSE, width = 13)
plotDR(daf, "UMAP", color_by = "meta30")
plotDR(daf, "UMAP", color_by = "CD3")
plotDR(daf, "UMAP", color_by = "CD4")
plotDR(daf, "UMAP", color_by = "CD8a")
plotDR(daf, "UMAP", color_by = "CD19")
plotDR(daf, "UMAP", color_by = "NCAM")
plotDR(daf, "UMAP", color_by = "CD45RA")
dev.off()

pdf("cytofALLMarkers_TcellsDaf2.pdf",useDingbats = FALSE, width = 13)
plotDR(daf2, "UMAP", color_by = "TBET")
plotDR(daf2, "UMAP", color_by = "FOXP3")
plotDR(daf2, "UMAP", color_by = "TCF1")
plotDR(daf2, "UMAP", color_by = "CD45RO")
plotDR(daf2, "UMAP", color_by = "CD44")
plotDR(daf2, "UMAP", color_by = "CCR7")
plotDR(daf2, "UMAP", color_by = "CD73")
plotDR(daf2, "UMAP", color_by = "CTLA4")
plotDR(daf2, "UMAP", color_by = "TIGIT")
plotDR(daf2, "UMAP", color_by = "CD27")
plotDR(daf2, "UMAP", color_by = "HELIOS")
plotDR(daf2, "UMAP", color_by = "x2B4")
plotDR(daf2, "UMAP", color_by = "KLRG1")
plotDR(daf2, "UMAP", color_by = "CD19")
plotDR(daf2, "UMAP", color_by = "EOMES")
plotDR(daf2, "UMAP", color_by = "GZMK")
plotDR(daf2, "UMAP", color_by = "CD45RA")
plotDR(daf2, "UMAP", color_by = "CD38")
plotDR(daf2, "UMAP", color_by = "CD4")
plotDR(daf2, "UMAP", color_by = "CD8a")
plotDR(daf2, "UMAP", color_by = "ICOS")
plotDR(daf2, "UMAP", color_by = "TOX")
plotDR(daf2, "UMAP", color_by = "CD3")
plotDR(daf2, "UMAP", color_by = "CD7")
plotDR(daf2, "UMAP", color_by = "IL7Ra")
plotDR(daf2, "UMAP", color_by = "FAS")
plotDR(daf2, "UMAP", color_by = "TIM3")
plotDR(daf2, "UMAP", color_by = "CD25")
plotDR(daf2, "UMAP", color_by = "CD45")
plotDR(daf2, "UMAP", color_by = "CXCR5")
plotDR(daf2, "UMAP", color_by = "Ki67")
plotDR(daf2, "UMAP", color_by = "x4_1BB")
plotDR(daf2, "UMAP", color_by = "PD1")
plotDR(daf2, "UMAP", color_by = "NCAM")
dev.off()

saveRDS(daf, file = "CyTOFfcsfileslivecells/LiveCellsCytof.rds")

daf <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/LiveCellsCytof.rds")

##expression per cluster    
pdf("clusterExprsHeatmap.pdf", useDingbats = FALSE, width = 10)
s <- plotClusterHeatmap(daf, hm2 = 'state_markers', k = "meta30", m = NULL,               
                   cluster_anno = TRUE, draw_freqs = TRUE)
dev.off()


#condition_colors <- c("controlLN" = "violet","tumorAccLN"="red","tumorLN"="aquamarine2","tumorPB"="cornflowerblue","tumorBM"="darkcyan","tumorHodgkinLN"="yellow")


## Facet per sample ## correct for pointSize
pdf("conditionSampleFacets.pdf", useDingbats = FALSE)
plotDR(daf1, dr="UMAP", color_by = "condition") + facet_wrap("condition")
plotDR(daf1, dr="UMAP", color_by = "meta30") +  facet_wrap("sample_id")
plotDR(daf1, dr="UMAP", color_by = "meta30") +  facet_wrap("condition")
dev.off()


### sample wise and marker heatmap ##ONe TIme###per sample_id ####meta30
med_exprs <- data.frame(t(daf1@assays@data$exprs), sample_id=sample_ids(daf1)) %>%
  group_by_(~sample_id) %>% summarize_all(funs(median))
med_exprs <- data.frame(med_exprs, row.names=1)
med_exprsScaled <- scale(med_exprs)
write.table(med_exprs, file = "med_exprsMarkers")
write.table(med_exprsScaled, file = "med_exprsMarkersScaled")

### sample wise and marker heatmap ##ONe TIme###per cluster ####meta30
med_exprs1 <- data.frame(t(daf2@assays@data$exprs), cluster_id=cluster_ids(daf2, k='meta30')) %>%
  group_by_(~cluster_id) %>% summarize_all(funs(median))
med_exprs1 <- data.frame(med_exprs1, row.names=1)
med_exprs1Scaled <- scale(med_exprs1)
write.table(med_exprs1, file = "med_exprsMarkersPerClusterDaf2")
write.table(med_exprs1Scaled, file = "med_exprsMarkersPerClusterScaled")


######table of cluster ids and sample ids -> number of cells in each cluster per sample used as input for abundance analysis
hh <- table(sample_ids(daf), cluster_ids(daf, k='meta30'))


#########1... median scaled expression _marker_sample : ALl cells
hm <- read.table("med_exprsMarkersPerCluster", header = T, sep = "\t", row.names = 1)

hm <- as.matrix(hm)


library(circlize)
library("RColorBrewer")


colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(200)

pdf("med_exprsMarkersPerCLuster.pdf", useDingbats = FALSE, width = 13)
Heatmap(hm, 
        cluster_rows = TRUE,column_names_gp = gpar(fontsize = 9),
        cluster_columns = TRUE, row_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = TRUE,col = colors,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson"
)
dev.off()


pdf('abundancesPerSamples.pdf', useDingbats = FALSE, width = 10)
plotAbundances(daf1, k = "meta30", by = "sample_id")
dev.off()

pdf('abundancesPerSamplesBoxPlots.pdf', useDingbats = FALSE, width = 10)
plotAbundances(daf1, k = "meta30", by = "cluster_id", shape = "patient_id")
dev.off()

#############subset T cells#########
setwd("/CyTOFfcsfileslivecells/CYTOFanalysis_REFER")

####did it first time

daf2 <- filterSCE(daf, cluster_id %in% c(1,2,3,4,5,6,7,8,10,
                                          13,14,15,22,
                                          20,21,23,26,27,28,30), k = 'meta30')


###cluster anew only Tcells 

set.seed(1347)




daf2 <- cluster(daf2, features = state_markers(daf2[c(2:19,23:38,40:43,45,46)]), 
                xdim = 10, ydim = 10, maxK = 30, seed = 1219) 

condition_colors = c("tumorAccLN"="red","tumorLN"="cyan","tumorPB"="cornflowerblue","tumorBM"="darkcyan","tumorHodgkinLN"="yellow", "controlLN" = "violet")


pdf("Counts_TcellsDaf2.pdf", useDingbats = FALSE, width = 10)
plotCounts(daf2, color_by = "condition") + scale_color_manual(values = condition_colors)
dev.off()

write.table(daf2@metadata$experiment_info, file="sampleCountsTcellsDaf2")

pdf("SampleMDS_TcellsDaf2.pdf", useDingbats = FALSE, width = 10)
plotMDS(daf2, color_by = "condition") + scale_color_manual(values = condition_colors)
dev.off()


daf2 <- runDR(daf2, dr = "UMAP", features = "state", cells=2500) 

pdf("clusterCellsUMAP_Tcells.pdf_Daf2.pdf",useDingbats = FALSE, width = 13)
plotDR(daf2, "UMAP", color_by = "meta30")  #+ scale_color_manual(values = my_color_palette)
dev.off()

pdf("cluster4kCellsmarkers_TcellsDaf2.pdf",useDingbats = FALSE, width = 13)
plotDR(daf2, "UMAP", color_by = "meta30")
plotDR(daf2, "UMAP", color_by = "CD3")
plotDR(daf2, "UMAP", color_by = "CD4")
plotDR(daf2, "UMAP", color_by = "CD8a")
plotDR(daf2, "UMAP", color_by = "CD19")
plotDR(daf2, "UMAP", color_by = "NCAM")
plotDR(daf2, "UMAP", color_by = "CD45RA")
dev.off()

hh <- table(sample_ids(daf2), cluster_ids(daf2, k='meta30'))
write.table(hh, file = "Daf2_Tcells_Table_perSampleClusterCells")

saveRDS(daf2, file = "CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsDaf2.rds")

daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsDaf2.rds")


pdf("clusterExprsHeatmap_TcellsDaf2check.pdf", useDingbats = FALSE, width = 10)
plotClusterHeatmap(daf2, hm2 = 'state_markers', k = "meta30", m = NULL,               
                        cluster_anno = TRUE, draw_freqs = TRUE)
dev.off()


pdf("conditionSampleFacets_TcellsDaf2.pdf", useDingbats = FALSE)
plotDR(daf2, dr="UMAP", color_by = "condition") + facet_wrap("condition")
plotDR(daf2, dr="UMAP", color_by = "meta30") +  facet_wrap("sample_id")
plotDR(daf2, dr="UMAP", color_by = "meta30") +  facet_wrap("condition")
dev.off()

### sample wise and marker heatmap ##ONe TIme###per sample_id ####meta30
med_exprs_Tcells <- data.frame(t(daf2@assays@data$exprs), sample_id=sample_ids(daf2)) %>%
  group_by_(~sample_id) %>% summarize_all(funs(median))
med_exprs_Tcells <- data.frame(med_exprs_Tcells, row.names=1)
med_exprs_TcellsScaled <- scale(med_exprs_Tcells)
write.table(med_exprs_Tcells, file = "med_exprsMarkers_TcellsDaf2")
write.table(med_exprs_TcellsScaled, file = "med_exprsMarkersScaled_TcellsDaf2")

### sample wise and marker heatmap ##ONe TIme###per cluster ####meta30
med_exprs1 <- data.frame(t(daf2@assays@data$exprs), cluster_id=cluster_ids(daf2, k='meta30')) %>%
  group_by_(~cluster_id) %>% summarize_all(funs(median))
med_exprs1 <- data.frame(med_exprs1, row.names=1)
med_exprs1Scaled <- scale(med_exprs1)
write.table(med_exprs1, file = "med_exprsMarkersPerClusterDaf2")
write.table(med_exprs1Scaled, file = "med_exprsMarkersPerClusterScaled")


pdf('abundancesPerSamples_TcellsDaf2.pdf', useDingbats = FALSE, width = 10)
plotAbundances(daf2, k = "meta30", by = "sample_id")
dev.off()

pdf('abundancesPerSamples_TcellsDaf2.pdf', useDingbats = FALSE, width = 10) ###with clusterAnno
plotAbundances(daf2, k = "meta30", by = "sample_id")
dev.off()

pdf('abundancesPerSamplesBoxPlots_TcellsDaf2.pdf', useDingbats = FALSE, width = 10)
plotAbundances(daf2, k = "meta30", by = "cluster_id", shape = "patient_id")
dev.off()

dd <- plotAbundances(daf2, k = "meta30", by = "sample_id")
write.table(dd$data, file = "cluster_freqs_TcellsDaf2")
##############heatmaps


library(circlize)
library("RColorBrewer")


colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(200)


#########.. median scaled expression per cluster

hm <- read.table("med_exprsMarkersPerClusterDaf2", header = T, sep = "\t", row.names = 1)

normalized = (hm-min(hm))/(max(hm)-min(hm))
normalized <- as.matrix(normalized)

pdf("med_exprsMarkersPerClusterDaf2.pdf", useDingbats = FALSE, width = 13)
Heatmap(hm, 
        cluster_rows = TRUE,column_names_gp = gpar(fontsize = 9),
        cluster_columns = TRUE, row_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = TRUE,col = colors,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson"
)
dev.off()

############MDS FIGURE 1A
daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsDaf2.rds")

condition_colors = c("tumorAccLN"="red","tumorLN"="cyan","tumorPB"="cornflowerblue","tumorBM"="darkcyan","tumorHodgkinLN"="yellow", "controlLN" = "violet")

setwd("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/Figures/")

pdf("SampleMDS_TcellsDaf2_Figure1A.pdf", useDingbats = FALSE, width = 10)
plotMDS(daf2, color_by = "condition") + scale_color_manual(values = condition_colors)
dev.off()

gg <- plotMDS(daf2, color_by = "condition") + scale_color_manual(values = condition_colors)
saveRDS(gg, file = 'figure1A_SampleMDS_TcellsDaf2.rds')

#########5... median  expression _marker_cluster_meta30: T cells ######FIGURE 1B
#####new annotation is added here
merging_table2 <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/ClusterAnnotation.csv", header = T, sep = ",")

head(data.frame(merging_table2)) 
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

pdf("med_exprsMarkersPerClusterDaf2.pdf", useDingbats = FALSE, width = 15, height = 9)
Heatmap(normalized, 
        cluster_rows = TRUE,column_names_gp = gpar(fontsize = 9),
        cluster_columns = TRUE, row_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = TRUE,col = colors,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson"
)
dev.off()

ll <- Heatmap(normalized, 
        cluster_rows = TRUE,column_names_gp = gpar(fontsize = 9),
        cluster_columns = TRUE, row_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = TRUE,col = colors,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson")

saveRDS(ll, file = 'med_exprsMarkersPerClusterDaf2Figure1B.rds')

saveRDS(dafN, file = 'TcellsDafN_withAnnotations.rds')
##############################cluster annotations 101####### FIGURE 1C



daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsDaf2.rds")


my_color_palette <- c("peru","#A6761D","plum1","orchid2","mediumorchid3",
                      "lightcyan2","lightblue2","lightskyblue", "skyblue3", "steelblue1", "dodgerblue1", "slateblue2", "royalblue2", "blue2",
                      "darkseagreen2","palegreen3", "aquamarine2","seagreen2","lawngreen", "limegreen", "forestgreen",
                      "tan1","sienna2","darkorange","tomato",
                      "lightpink1","palevioletred2", "hotpink1",
                      "yellow1","burlywood3")

pdf("UMAPannottaion_Tcells_2.pdf", useDingbats = FALSE)
plotDR(dafN, "UMAP", color_by = "merging4") + scale_color_manual(values = my_color_palette)
dev.off()

plot2 <- plotDR(dafN, "UMAP", color_by = "merging4") + scale_color_manual(values = my_color_palette)

pdf("UMAPannottaion_TcellsWithoutLegend.pdf", useDingbats = FALSE)
plotDR(dafN, "UMAP", color_by = "merging4") + scale_color_manual(values = my_color_palette) +
  theme(legend.position = "none")
dev.off()

plot3 <- plotDR(dafN, "UMAP", color_by = "merging4") + scale_color_manual(values = my_color_palette) +
  theme(legend.position = "none") 

saveRDS(plot2, file = 'figure1C_umapAnnotations.rds')
saveRDS(plot3, file = 'figure1C_umapAnnotationsWithoutLegend.rds') ###  more space for the clusters


#theme(legend.position = "none")

##############hm6 abundances and clinical parameters ########FIGURE 1D
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

freqs2 <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/Daf2_Tcells_Table_perSampleClusterCells", header = T, sep = "\t", row.names = 1)
freqs22 <- freqs2/apply(freqs2,1,sum)

hc <- hclust(dist(freqs22),method="ward.D")

meta <- read.table("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/files", header = T, sep = '\t')

haB <- rowAnnotation(condition = meta$condition,
                     TumorSite = meta$TumorSite,
                     IGHVstatus = meta$IGHVstatus,
                     gender = meta$gender,
                     col=list(condition = c("tumorLN"="yellow",
                                            "tumorPB"="green",
                                            "tumorBM"="darkcyan",
                                            "controlLN" = "violet"),
                              TumorSite = c("LymphNode"="orange","PeripheralBlood"="grey","BoneMarrow"="pink"),
                              IGHVstatus=c("mutated"="coral1","unmutated"="darkolivegreen2","NA"="grey"),
                              gender=c("male"="lightblue4","female"="hotpink4")),
                     width = unit(4,"cm"),
                     show_legend =TRUE,
                     rn = anno_text(rownames(freqs2),gp = gpar(fontsize = 7, fontface="bold")),
                     annotation_name_gp = gpar(fontsize = 10, fontface="bold"))

####k_pal = CATALYST:::.cluster_cols (use this identify default colors in the daf2 object)

my_color_palette <- c("peru","#A6761D","plum1","orchid2","mediumorchid3",
                      "lightcyan2","lightblue2","lightskyblue", "skyblue3", "steelblue1", "dodgerblue1", "slateblue2", "royalblue2", "blue2",
                      "darkseagreen2","palegreen3", "aquamarine2","seagreen2","lawngreen", "limegreen", "forestgreen",
                      "tan1","sienna2","darkorange","tomato",
                      "lightpink1","palevioletred2", "hotpink1",
                      "yellow1","burlywood3")

ha_list2 = rowAnnotation(ClusterAbundance = anno_barplot(freqs22, bar_width = 1.0,
                                                         width = unit(6, "cm"),
                                                         height = unit(10, "cm"),gp=gpar(fill=my_color_palette),
                                                         which = "row"),
                         show_legend = TRUE, annotation_name_gp = gpar(fontsize = 10, fontface="bold"))





pdf("abundanceHM3_WithAnno_FIGURE1D.pdf", useDingbats = FALSE)
 Heatmap(matrix(nc = 0, nr = 46), cluster_rows = hc, 
        left_annotation = ha_list2, show_row_dend = T, 
        right_annotation = haB,show_heatmap_legend = TRUE) 
dev.off()

kk <- Heatmap(matrix(nc = 0, nr = 46), cluster_rows = hc$order, 
              left_annotation = ha_list2, show_row_dend = T, 
              right_annotation = haB,show_heatmap_legend = TRUE) 
saveRDS(kk, file = 'figure1D_abundanceHM3_WithAnno.rds')






