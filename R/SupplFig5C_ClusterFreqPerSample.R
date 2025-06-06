###### Cluster Freq per sample (as John did)
library(ggalluvial)
library(tidyverse)
library(ggplot2)

set.seed(0)

# Extract cluster and sample information from scDataMerged_harmony
clustersstats <- scDataMerged_harmony$seurat_clusters
clustersstatsTable <- data.frame(table(scDataMerged_harmony$seurat_clusters, scDataMerged_harmony$orig.ident))
clustersstatsTableTotalPerSample <- data.frame(table(scDataMerged_harmony$orig.ident))
rownames(clustersstatsTableTotalPerSample) <- clustersstatsTableTotalPerSample$Var1
clustersstatsTable$prop <- clustersstatsTable$Freq / clustersstatsTableTotalPerSample[clustersstatsTable$Var2,]$Freq

# Generate colors for each of the 16 clusters (0 to 15)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colorclust <- gg_color_hue(16)  # Set color for 16 clusters

# No renaming of samples needed, so directly set factors
clustersstatsTable$Var2 <- as.factor(clustersstatsTable$Var2)

# Optional: Add a dummy reference sample if you need a reference across clusters (adjust as necessary)
DummyRef <- clustersstatsTable[which(clustersstatsTable$Var2 == levels(clustersstatsTable$Var2)[1]), ]
DummyRef$Var2 <- "Ref"
DummyRef$prop <- rep(1 / 32, 16) + ((15 - as.numeric(as.character(DummyRef$Var1))) / 240)
DummyRef$name <- DummyRef$Var1
clustersstatsTable$name <- ""
clustersstatsTable <- rbind(clustersstatsTable, DummyRef)

# Plot without renaming samples
g <- ggplot(clustersstatsTable, aes(y = prop, x = Var2, fill = (Var1))) +
  geom_flow(aes(alluvium = (Var1)), alpha = .5, color = "black", curve_type = "linear", width = .8) +
  geom_col(width = .8, color = "black") +
  scale_fill_manual(values = colorclust) +
  theme_bw() +
  geom_text(aes(label = name, y = prop), position = position_stack(0.5)) +
  xlab("") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15))

# Save plot with an appropriate file path
ggsave(filename = paste0("barplot_harmony_cluster_proportion.pdf"), 
       width = 5, height = 5, plot = g)
