library(scRepertoire)

## Using the Seurat object that includes in meta.data the predictTCR information
scDataMerged_harmony <- qread("scDataMerged_harmony_23102024.qs")

######### Suppl. Figure 7A: % expanded clones per sample LN vs PB ####
library(ggplot2)
library(dplyr)

metadata <- scDataMerged_harmony@meta.data

# Define the clone categories for expanded clones
expanded_clone_categories <- c("Small (1 < X <= 5)", "Medium (5 < X <= 20)", 
                               "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 1250)")
metadata$expandedClone <- metadata$cloneSize %in% expanded_clone_categories

# Group by sample and tissue type, calculate the percentage of expanded clones
clone_expansion_percent <- metadata %>%
  group_by(orig.ident, Tissue) %>%
  summarise(expandedPercent = mean(expandedClone) * 100) %>%
  ungroup()

# Plot
pdf(file = "CLL_LN_CloneFreq_PBvsLN.pdf", width = 6, height = 6)
ggplot(clone_expansion_percent, aes(x = Tissue, y = expandedPercent, color = Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, aes(color = Tissue), size = 2, alpha = 0.6) +
  scale_color_manual(values = c("PB" = "#560BAD", "LN" = "#00A6FB")) +
  labs(x = "Tissue", y = "Percentage of Expanded Clones (%)", title = "Percentage of Expanded Clones by Sample") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
dev.off()


##### Suppl. Figure 7B: Clone Size per cluster in PB and LN #####
PB <- subset(scDataMerged_harmony, subset = Tissue == "PB")
LN <- subset(scDataMerged_harmony, subset = Tissue == "LN")

pdf(file = "CLL_LN_ClonalOccupy_perClusterLN_27112024.pdf",
    width = 10, height = 6) 
clonalOccupy(LN,
             x.axis = "RNA_snn_res.0.8", proportion = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(file = "CLL_LN_ClonalOccupy_perClusterPB_27112024.pdf",
    width = 10, height = 6) 
clonalOccupy(PB,
             x.axis = "RNA_snn_res.0.8", proportion = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


##### Suppl. Figure 7C: Cluster distribution of shared expanded T-cell clones #####
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggalluvial)
library(tidyverse)

# Define colors for each cluster
cluster_colors <- c(
  "0" = "#F8766D", "1" = "#E68613", "2" = "#CD9600", "3" = "#ABA300", 
  "4" = "#7CAE00", "5" = "#0BB702", "6" = "#00BE67", "7" = "#00C19A", 
  "8" = "#00BFC4", "9" = "#00B8E7", "10" = "#00A9FF", "11" = "#8494FF", 
  "14" = "#FF61CC", "15" = "#FF68A1"
)


# Filter CLL clusters 12 and 13 
scDataFiltered <- subset(scDataMerged_harmony, 
                         seurat_clusters %in% setdiff(0:15, c(12, 13)))

# Identify expanded clones across PB and LN
expanded_clones <- scDataFiltered@meta.data %>%
  filter(cloneSize %in% c("Small (1 < X <= 5)", "Medium (5 < X <= 20)", 
                          "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 1250)")) %>%
  group_by(CTstrict) %>%
  filter(n_distinct(Tissue) > 1) %>%
  ungroup()
tissue_cluster_freq <- top_clone_data %>%
  group_by(Tissue, seurat_clusters) %>%
  summarize(Frequency = n(), .groups = "drop") %>%
  group_by(Tissue) %>%
  mutate(Proportion = 100 * Frequency / sum(Frequency)) 

alluvial_data <- tissue_cluster_freq %>%
  mutate(Var1 = factor(seurat_clusters), Var2 = Tissue)  

# Plot the alluvial diagram
g <- ggplot(alluvial_data, aes(y = Proportion, x = Var2, fill = Var1)) +
  geom_flow(aes(alluvium = Var1), alpha = .5, color = "black", curve_type = "linear", width = .8) +
  geom_col(width = .8, color = "black") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15))
ggsave(filename = "CLL_LN_SharedClones_Alluvial_PBvsLN.pdf", 
       width = 6, height = 6, plot = g)


##### Suppl. Figure 7F: UMAP plot colored according to tumor reactivity as predicted by predictTCR #####
custom_colors <- c("Unknown" = "#4d4d4d", "Non_Reactive" = "#D8B365", "Reactive" = "#5AB4AC", "NA" = "#bababa")

pdf(file = "CLL_LN_UMAP_TumorReactive.pdf", width = 10, height = 6)
DimPlot(scDataMerged_harmony, group.by = "Reactivity", raster = FALSE, label = TRUE) +
  scale_color_manual(values = custom_colors) + 
  theme(plot.title = element_blank())

dev.off()

##### Suppl. Figure 7G: UMAP with tumor-reactive cells colored according to the clone size #####
library(Seurat)
library(ggplot2)

# Set the color palette for cloneSize
clone_size_colors <- c("Single (0 < X <= 1)" = "#0D0887FF",                         
                       "Small (1 < X <= 5)" = "#7301A8FF",                         
                       "Medium (5 < X <= 20)" = "#BD3786FF",                         
                       "Large (20 < X <= 100)" = "#FA9E3BFF",                         
                       "Hyperexpanded (100 < X <= 1250)" = "#F0F921FF",                        
                       "NA" = "grey", 
                       "Not Reactive" = "grey") 


umap_data <- Embeddings(scDataMerged_harmony, reduction = "umap")
plot_data <- data.frame(UMAP_1 = umap_data[, 1], 
                        UMAP_2 = umap_data[, 2],
                        Reactivity = scDataMerged_harmony$Reactivity,
                        cloneSize = scDataMerged_harmony$cloneSize)

# Create a new column in metadata for the color group
plot_data$color_group <- ifelse(plot_data$Reactivity == "Reactive", 
                                as.character(plot_data$cloneSize), 
                                "Not Reactive")

# Update NAs in cloneSize to "NA" for consistent coloring
plot_data$color_group[is.na(plot_data$cloneSize)] <- "NA"


plot_data$color_group <- factor(plot_data$color_group, 
                                levels = c("Not Reactive", "NA", 
                                           "Single (0 < X <= 1)", 
                                           "Small (1 < X <= 5)", 
                                           "Medium (5 < X <= 20)", 
                                           "Large (20 < X <= 100)", 
                                           "Hyperexpanded (100 < X <= 1250)"))

pdf(file = "CLL_LN_UMAP_Tumor-Reactive_ClonalHomeostasis_customized.pdf", 
    width = 10, height = 6)
p_custom <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = color_group)) +
  geom_point(data = subset(plot_data, Reactivity != "Reactive"), size = 0.1) + 
  geom_point(data = subset(plot_data, Reactivity == "Reactive"), size = 0.1) + 
  scale_color_manual(values = clone_size_colors) 
dev.off()


##### Suppl. Figure 7H: Bar plot displaying the CLL-reactive T cell clone size percentage

custom_colors <- c(
  "Single (0 < X <= 1)" = "#0D0887FF",
  "Small (1 < X <= 5)" = "#7301A8FF",
  "Medium (5 < X <= 20)" = "#BD3786FF",
  "Large (20 < X <= 100)" = "#FA9E3BFF",
  "Hyperexpanded (100 < X <= 1250)" = "#F0F921FF"
)

pdf(file = "CLL_LN_ClonalOccupy_perCluster_CLLReactive_29112024.pdf", 
    width = 10, height = 6)
clonalOccupy(Reactive, 
             x.axis = "RNA_snn_res.0.8", 
             proportion = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_colors)
dev.off()


