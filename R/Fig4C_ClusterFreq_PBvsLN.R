library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(qs)
library("paletteer")
library("RColorBrewer")
library(tidyverse)
library(circlize)
library(scCustomize)


###### T cell subset frequency out of total T cells
# Load data
freqs1 <- read.csv("my_data 2.csv", header = TRUE, row.names = 1)

# Exclude columns "X12" and "X13" (CLL clusters) from the frequency calculations
freqs1_filtered <- freqs1[, !colnames(freqs1) %in% c("X12", "X13")]
freqs11 <- freqs1_filtered / rowSums(freqs1_filtered) 
freqs11$Samples <- rownames(freqs11)  

# Load metadata
meta <- read.csv("/meta.data.csv", header = TRUE, row.names = 1)
meta$Samples <- rownames(meta)

freqs11 <- freqs11[, c("Samples", colnames(freqs11)[-ncol(freqs11)])]
meta <- meta[, c("Samples", colnames(meta)[-ncol(meta)])]
merged_data <- merge(meta, freqs11, by = "Samples")

# Define tissue order and colors
tissue_order <- c("PB", "LN")
merged_data$Tissue <- factor(merged_data$Tissue, levels = tissue_order)
tissue_colors <- c(PB = "#560bad", LN = "#00a6fb")

tcell_subset_names <- setdiff(names(merged_data), c("Samples", "Tissue", "Patient", "X12", "X13"))

t_test_results <- list()

# Perform a t-test without multiple testing correction
for (subset_name in tcell_subset_names) {
  formula_str <- as.formula(paste(subset_name, "~ Tissue"))
  test_result <- t.test(
    formula_str,
    data = merged_data
  )
  t_test_results[[subset_name]] <- test_result
}

####### Boxplots #######
# Convert `merged_data` to long format for grouped boxplots
long_data <- melt(merged_data, id.vars = c("Samples", "Tissue"), 
                  measure.vars = tcell_subset_names,
                  variable.name = "Subset", value.name = "Frequency")

# Create boxplot with T cell subsets on the x-axis, grouped by tissue type
p <- ggplot(long_data, aes(x = Subset, y = Frequency, fill = Tissue)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             color = "black", size = 1) +
  scale_fill_manual(values = tissue_colors, name = "Tissue") +
  labs(x = "T Cell Subset", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "top")
print(p)
