scDataMerged_harmony <- qread("scDataMerged_harmony_23102024.qs")


#List of data frames
BC0_LN <- read.csv("BC0_LN.csv")
BC1_LN <- read.csv("BC1_LN.csv")
BC2_LN <- read.csv("BC2_LN.csv")
BC3_LN <- read.csv("BC3_LN.csv")
BC9_LN <- read.csv("BC9_LN.csv")
BC0_PB <- read.csv("BC0_PB.csv")
BC1_PB <- read.csv("BC1_PB.csv")
BC2_PB <- read.csv("BC2_PB.csv")
BC3_PB <- read.csv("BC3_PB.csv")
BC9_PB <- read.csv("BC9_PB.csv")


contig_list <- list(BC0_LN, BC1_LN, BC2_LN, BC3_LN, BC9_LN, BC0_PB, 
                    BC1_PB, BC2_PB, BC3_PB, BC9_PB)
contig_list <- loadContigs(input = contig_list, 
                           format = "WAT3R")

head(contig_list[[1]])

combined <- combineTCR(contig_list, 
                       samples = c("BC0_LN", "BC1_LN", "BC2_LN", "BC3_LN", 
                                   "BC9_LN","BC0_PB", "BC1_PB", "BC2_PB", "BC3_PB", 
                                   "BC9_PB"))

## Rename cell barcodes from Seurat Object
barcodes_scDataMerged <- colnames(scDataMerged_harmony)

barcodes_scDataMerged_modified <- gsub("_BC0_lymph-node-tumor-00", "", barcodes_scDataMerged)
barcodes_scDataMerged_modified <- gsub("_BC1_lymph-node-tumor-01", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_BC2_lymph-node-tumor-02", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_BC3_lymph-node-tumor-01", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_BC9_lymph-node-tumor09", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_15-8601_pbmc-bc3", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_17-9411_pbmc-bc2", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_10-7433_pbmc-bc0", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_09-316_pbmc-bc1", "", barcodes_scDataMerged_modified)
barcodes_scDataMerged_modified <- gsub("_02-97_pbmc-bc9", "", barcodes_scDataMerged_modified)

new_barcodes <- setNames(barcodes_scDataMerged_modified, barcodes_scDataMerged)
scDataMerged_harmony <- RenameCells(scDataMerged_harmony, new.names = new_barcodes)


qsave(scDataMerged_harmony, "scDataMerged_harmony2_23102024.qs")


##### Figure 5A: Clone size frequency per sample in PB and LN for CD8 and CD4 T cells ####
library(RColorBrewer)
library(viridis)
library(reshape2) 
library(data.table) 
library(dplyr)
######classify clusters
CD4clusters = c(1, 2, 3, 4,7, 8)
CD8clusters = c(0, 5, 6, 10, 11, 14, 15)
# Use original clustering column (replace "subcluster" with the main cluster column)
clonerecurrCD4 = seurat@meta.data[which(seurat@meta.data$seurat_clusters %in% CD4clusters), c("orig.ident", "cloneSize")]
clonerecurrCD8 = seurat@meta.data[which(seurat@meta.data$seurat_clusters %in% CD8clusters), c("orig.ident", "cloneSize")]
clonerecurrCD4$marker = "CD4"
clonerecurrCD8$marker = "CD8"
clonerecurrCD4_CD8 = rbind(clonerecurrCD4,clonerecurrCD8)
clonerecurrCD4_CD8 = clonerecurrCD4_CD8[which(!is.na(clonerecurrCD4_CD8$cloneSize)),]
clonerecurrCD4_CD8$cloneSize = paste(clonerecurrCD4_CD8$marker,clonerecurrCD4_CD8$cloneSize )
clonerecurrCD4_CD8sum =  melt(table(clonerecurrCD4_CD8))
#clonerecurrCD4_CD8sum <- melt(table(clonerecurrCD4_CD8), na.rm = TRUE)
setDT(clonerecurrCD4_CD8sum)[, frac := value / sum(value), by=orig.ident]
clonerecurrCD4_CD8sum = clonerecurrCD4_CD8sum[which(clonerecurrCD4_CD8sum$frac > 0),]
clonerecurrCD4_CD8sum$cloneSize = factor(clonerecurrCD4_CD8sum$cloneSize, levels=c("CD4 Hyperexpanded (100 < X <= 1250)",
                                                                                   "CD4 Large (20 < X <= 100)",
                                                                                   "CD4 Medium (5 < X <= 20)",
                                                                                   "CD4 Small (1 < X <= 5)",
                                                                                   "CD4 Single (0 < X <= 1)",
                                                                                   "CD8 Hyperexpanded (100 < X <= 1250)",
                                                                                   "CD8 Large (20 < X <= 100)",
                                                                                   "CD8 Medium (5 < X <= 20)",
                                                                                   "CD8 Small (1 < X <= 5)",
                                                                                   "CD8 Single (0 < X <= 1)"))

# Set cohort name
scCohort <- "your_sc_cohort_name"

clonebarplotf <- ggplot(data=clonerecurrCD4_CD8sum, aes(x=orig.ident, y=frac, fill=cloneSize)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(10, "Set1"))(12))+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  ylab("Fraction of cells")+ 
  xlab("Sample ID") 

pdf(file = "CLL_LN_UMAP_clonalHomeostasis.pdf", width = 12, height = 10)
clonebarplotf
dev.off()

ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/TCR/",scCohort,"/TCR_cell_CD4_CD8_fraction_barplot.pdf"),width=7,height=6, plot=clonebarplotf)

#CD4 only
clonerecurrCD4sum =  melt(table(clonerecurrCD4))
setDT(clonerecurrCD4sum)[, frac := value / sum(value), by=orig.ident]
clonerecurrCD4sum$cloneSize = paste("CD4",clonerecurrCD4sum$cloneSize)
clonerecurrCD4sum$cloneSize = factor(clonerecurrCD4sum$cloneSize, levels=c("CD4 Hyperexpanded (100 < X <= 1250)","CD4 Large (20 < X <= 100)","CD4 Medium (5 < X <= 20)","CD4 Small (1 < X <= 5)","CD4 Single (0 < X <= 1)"))
clonerecurrCD4sum$pct = paste(round(clonerecurrCD4sum$frac * 100),"%",sep="")
clonerecurrCD4sum$pct[which(clonerecurrCD4sum$frac < 0.01)] = NA

clonebarplotf <- ggplot(data=clonerecurrCD4sum, aes(x=orig.ident, y=frac, fill=cloneSize)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(10, "Set1"))(12))+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  ylab("Fraction of cells")+ 
  xlab("Sample ID") +
  geom_text(
    aes(label = pct), 
    position = position_stack(vjust = 0.5)
  ) 


pdf(file = "CLL_LN_UMAP_clonalHomeostasis_CD4.pdf", width = 12, height = 10)
clonebarplotf
dev.off()


ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/TCR/",scCohort,"/TCR_cell_CD4_fraction_barplot.pdf"),width=6,height=6, plot=clonebarplotf)


clonerecurrCD4sum =  melt(table(clonerecurrCD4))
setDT(clonerecurrCD4sum)[, frac := value / sum(value), by=orig.ident]
clonerecurrCD4sum$cloneSize = paste("CD4",clonerecurrCD4sum$cloneSize)
clonerecurrCD4sum$cloneSize = factor(clonerecurrCD4sum$cloneSize, levels=c("CD4 Hyperexpanded (100 < X <= 500)","CD4 Large (20 < X <= 100)","CD4 Medium (5 < X <= 20)","CD4 Small (1 < X <= 5)","CD4 Single (0 < X <= 1)"))

#CD8
clonerecurrCD8sum =  melt(table(clonerecurrCD8))
setDT(clonerecurrCD8sum)[, frac := value / sum(value), by=orig.ident]
clonerecurrCD8sum$cloneSize = paste("CD8",clonerecurrCD8sum$cloneSize)
clonerecurrCD8sum$cloneSize = factor(clonerecurrCD8sum$cloneSize, levels=c("CD8 Hyperexpanded (100 < X <= 500)","CD8 Large (20 < X <= 100)","CD8 Medium (5 < X <= 20)","CD8 Small (1 < X <= 5)","CD8 Single (0 < X <= 1)"))
clonerecurrCD8sum$pct = paste(round(clonerecurrCD8sum$frac * 100),"%",sep="")
clonerecurrCD8sum$pct[which(clonerecurrCD8sum$frac == 0)] = NA

clonebarplotf <- ggplot(data=clonerecurrCD8sum, aes(x=orig.ident, y=frac, fill=cloneSize)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(10, "Set1"))(12)[6:12])+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  ylab("Fraction of cells")+ 
  xlab("Sample ID") +
  geom_text(
    aes(label = pct), 
    position = position_stack(vjust = 0.5)
  )

pdf(file = "CLL_LN_UMAP_clonalHomeostasis_CD8.pdf", width = 12, height = 10)
clonebarplotf
dev.off()

ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/TCR/",scCohort,"/TCR_cell_CD8_fraction_barplot.pdf"),width=6,height=6, plot=clonebarplotf)


##### Figure 5B: UMAP colored accoridng to Clone Size ######

colorblind_vector <- hcl.colors(n=9, palette = "inferno", fixup = TRUE)

scDataMerged_harmony <- combineExpression(combined, scDataMerged_harmony, cloneCall="strict", 
                                          group.by = "sample", proportion = FALSE, filterNA = FALSE,
                                          cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=1250))


pdf(file = "CLL_LN_UMAP_cloneSize.pdf",
    width = 10, height = 6) 
Seurat::DimPlot(scDataMerged_harmony, group.by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,5,6,7,9)]))
dev.off()


##### Figure 5C: Shannon Index per cluster in PB and LN #####
# Calculate Shannon diversity for PB
plot <- clonalDiversity(PB, cloneCall = "strict", group.by = "RNA_snn_res.0.8")
ShannonData <- plot$data[which(plot$data$variable == "shannon"),]
ShannonData$norm <- ShannonData$value / max(ShannonData$value)
ShannonData <- ShannonData[order(ShannonData$value),]
ShannonData$samples <- ShannonData$RNA_snn_res.0.8
clusternames <- data.frame(row.names=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                           cluster=c("CD8 TPEX","CD4 TCM1","CD4 TN","CD4 TEM","CD4 TCM2","CD8 TEM",
                                     "CD8 TEM2", "CD4 TREG", "T FH","NK","MAIT","CD8 TN","CLL","CLL2","CD8 TEX", "TPR"))
ShannonData$cluster <- clusternames[ShannonData$samples, ]
ShannonData$cluster <- factor(ShannonData$cluster, levels=rev(clusternames$cluster))

# Filter out CLL clusters
ShannonData <- ShannonData[!ShannonData$cluster %in% c("CLL", "CLL2"), ]
ShannonData$cluster <- factor(ShannonData$cluster, levels=rev(unique(ShannonData$cluster)))

ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=length(unique(ShannonData$cluster)))
ShannonData$color <- color_list[as.numeric(ShannonData$cluster)]

Shannonmin <- 0.5 

# Plot
shaplot_norm <- ggplot(ShannonData, aes(x=norm, y=cluster)) +
  geom_segment(aes(x=Shannonmin, xend=norm, y=cluster, yend=cluster), color="grey") + 
  geom_point(aes(x=norm, y=cluster, colour=cluster), size=4) +
  theme_classic() +
  scale_color_manual(breaks=as.character(ShannonData$cluster), values=ShannonData$color) +
  xlim(Shannonmin, 1)

# Save normalized plot
pdf(file = "CLL_LN_ShannonIndex_normalized_PBplot.pdf", width = 12, height = 10)
print(shaplot_norm)
dev.off()

# Calculate Shannon diversity for LN
plot <- clonalDiversity(LN, cloneCall = "strict", group.by = "RNA_snn_res.0.8")
ShannonData <- plot$data[which(plot$data$variable == "shannon"),]
ShannonData$norm <- ShannonData$value / max(ShannonData$value)
ShannonData <- ShannonData[order(ShannonData$value),]
ShannonData$samples <- ShannonData$RNA_snn_res.0.8
clusternames <- data.frame(row.names=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                           cluster=c("CD8 TPEX","CD4 TCM1","CD4 TN","CD4 TEM","CD4 TCM2","CD8 TEM",
                                     "CD8 TEM2", "CD4 TREG", "T FH","NK","MAIT","CD8 TN","CLL","CLL2","CD8 TEX", "TPR"))
ShannonData$cluster <- clusternames[ShannonData$samples, ]
ShannonData$cluster <- factor(ShannonData$cluster, levels=rev(clusternames$cluster))

# Filter out CLL clusters
ShannonData <- ShannonData[!ShannonData$cluster %in% c("CLL", "CLL2"), ]
ShannonData$cluster <- factor(ShannonData$cluster, levels=rev(unique(ShannonData$cluster)))

ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=length(unique(ShannonData$cluster)))
ShannonData$color <- color_list[as.numeric(ShannonData$cluster)]

Shannonmin <- 0.5 

# Plot
shaplot_norm <- ggplot(ShannonData, aes(x=norm, y=cluster)) +
  geom_segment(aes(x=Shannonmin, xend=norm, y=cluster, yend=cluster), color="grey") + 
  geom_point(aes(x=norm, y=cluster, colour=cluster), size=4) +
  theme_classic() +
  scale_color_manual(breaks=as.character(ShannonData$cluster), values=ShannonData$color) +
  xlim(Shannonmin, 1)

# Save normalized plot
pdf(file = "CLL_LN_ShannonIndex_normalized_LNplot.pdf", width = 12, height = 10)
print(shaplot_norm)
dev.off()


###### Figure 5D: Alluvial plot displaying the top 10 most frequent clones for LN and PB for each patient #####
# Example for BC0
pdf(file = "CLL_LN_UMAP_clonalCompareBC0.pdf", width = 12, height = 10)
clonalCompare(scDataMerged_harmony, 
              top.clones = 10, 
              samples = c("BC0_PB", "BC0_LN"), 
              cloneCall="strict", group.by = "orig.ident",
              graph = "alluvial", relabel.clones = TRUE)
dev.off()



##### Addition of predictTCR information into Seurat Object
predictTCR <- read.csv("/slgpfs/scratch/cli79/cli79431/Projects/Nat Comm/CLL predicTCR.csv")

# Apply multiple gsub replacements to standardize the barcodes in predictTCR$X
predictTCR$X <- gsub("_BC0_lymph-node-tumor-00", "", predictTCR$X)
predictTCR$X <- gsub("_BC1_lymph-node-tumor-01", "", predictTCR$X)
predictTCR$X <- gsub("_BC2_lymph-node-tumor-02", "", predictTCR$X)
predictTCR$X <- gsub("_BC3_lymph-node-tumor-01", "", predictTCR$X)
predictTCR$X <- gsub("_BC9_lymph-node-tumor09", "", predictTCR$X)
predictTCR$X <- gsub("_15-8601_pbmc-bc3", "", predictTCR$X)
predictTCR$X <- gsub("_17-9411_pbmc-bc2", "", predictTCR$X)
predictTCR$X <- gsub("_10-7433_pbmc-bc0", "", predictTCR$X)
predictTCR$X <- gsub("_09-316_pbmc-bc1", "", predictTCR$X)
predictTCR$X <- gsub("_02-97_pbmc-bc9", "", predictTCR$X)

# Set the modified barcodes as row names for predictTCR
rownames(predictTCR) <- predictTCR$X
predictTCR$X <- NULL  # Remove the barcode column after setting row names

# Find matching barcodes between predictTCR and scDataMerged_harmony
matching_barcodes <- intersect(rownames(predictTCR), colnames(scDataMerged_harmony))

# Merge the metadata into the Seurat object
scDataMerged_harmony <- AddMetaData(object = scDataMerged_harmony, metadata = predictTCR[matching_barcodes, ])


##### Figure 5E: Proportion of reactive, non-reactive and unknown/ NA T cell clonotypes out of total T cells in LN and PB #####
library(dplyr)
library(ggplot2)

# Create a table of counts for each Reactivity category within each seurat_cluster and Tissue
reactivity_summary <- scDataMerged_harmony@meta.data %>%
  group_by(Tissue, seurat_clusters, Reactivity) %>%
  summarise(Frequency = n()) %>%
  ungroup()

# Filter out CLL clusters 12 and 13
reactivity_summary_filtered_no_clusters <- reactivity_summary_filtered %>%
  filter(!seurat_clusters %in% c(12, 13))


pdf(file = "CLL_LN_UMAP_TumorReactive_ClusterFreq_noNA_no12_13.pdf", 
    width = 10, height = 6)
ggplot(reactivity_summary_filtered_no_clusters, aes(x = seurat_clusters, y = Frequency, fill = as.factor(Reactivity))) +
  geom_bar(stat = "identity", position = "fill") +  
  facet_wrap(~ Tissue) +  
  scale_fill_manual(values = c("Unknown" = "#4d4d4d", 
                               "Non_Reactive" = "#D8B365", 
                               "Reactive" = "#5AB4AC")) +  
  labs(x = "Seurat Clusters", y = "Proportion", fill = "Reactivity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

##### Figure 5F: Scatter plot showing LN and PB clone sizes from all 5 CLL patients #####
clonal_scatter_output <- clonalScatter(scDataMerged_harmony, group.by = "Tissue",
                                       cloneCall = "strict", 
                                       x.axis = "PB", 
                                       y.axis = "LN",
                                       dot.size = "total", exportTable = "TRUE",
                                       graph = "proportion")

clonal_scatter_output$PB.fraction.plot <- ifelse(clonal_scatter_output$PB.fraction == 0, 1e-5, clonal_scatter_output$PB.fraction)
clonal_scatter_output$LN.fraction.plot <- ifelse(clonal_scatter_output$LN.fraction == 0, 1e-5, clonal_scatter_output$LN.fraction)

fraction_limits <- c(1e-5, max(c(clonal_scatter_output$PB.fraction.plot, clonal_scatter_output$LN.fraction.plot), na.rm = TRUE) * 1.1)

reactivity_colors <- c("NA" = "#bababa", "Unknown" = "#4d4d4d", "Non_Reactive" = "#d8b365", "Reactive" = "#5ab4ac")

# Make sure Reactivity is mapped correctly to clonal_scatter_output
clonal_scatter_output$Reactivity <- scDataMerged_harmony$Reactivity[match(clonal_scatter_output$Var1, scDataMerged_harmony$CTstrict)]

# Plot
pdf(file = "CLL_LN_UMAP_clonalScatter_log_27112024.pdf", width = 12, height = 10)
ggplot_clonal_scatter <- ggplot(clonal_scatter_output, aes(x = PB.fraction.plot, y = LN.fraction.plot, size = PB + LN, color = Reactivity)) +
  geom_point(alpha = 1/2) +  # Make points slightly transparent
  scale_size_continuous(name = "Total PB + LN", range = c(5, 30)) +  # Adjust point size for better visibility
  labs(title = "PB vs LN", x = "PB Fraction", y = "LN Fraction") +
  theme_light() +
  theme(
    axis.title = element_text(size = 20),  
    axis.text = element_text(size = 20),   
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 20),  
    plot.title = element_text(size = 24),  
    strip.text = element_text(size = 20)   
  ) +
  scale_x_log10(limits = fraction_limits, expand = expansion(mult = 0.1)) +  
  scale_y_log10(limits = fraction_limits, expand = expansion(mult = 0.1)) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +  
  scale_color_manual(values = reactivity_colors) +  
  guides(color = guide_legend(override.aes = list(size = 6)))  
dev.off()



