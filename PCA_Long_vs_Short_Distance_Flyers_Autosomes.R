#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

# ---------------------------
# CONFIG PATHS
# ---------------------------
pca_file <- "pigeons_pca.eigenvec"  # Output from PLINK PCA
metadata_file <- "metadata.txt"      # Must include IID and FlightType
output_file <- "pca_long_vs_short_all_chr.png"

# ---------------------------
# LOAD DATA
# ---------------------------
pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# Ensure IID columns are character for merging
pca$IID <- as.character(pca$IID)
metadata$IID <- as.character(metadata$IID)

# FIX 1: Reconstruct full ID (e.g., NPO_146)
pca$IID <- paste0("NPO_", pca$IID)

# ---------------------------
# MERGE PCA WITH METADATA
# ---------------------------
pca <- left_join(pca, metadata, by = "IID")

# Check for missing FlightType
if (any(is.na(pca$FlightType))) warning("Some PCA samples did not match metadata.")

# ---------------------------
# FILTER TO LONG VS SHORT
# ---------------------------
pca <- pca %>% filter(FlightType %in% c("Long", "Short"))

cat("Number of samples per flight type:\n")
print(table(pca$FlightType))

# ---------------------------
# PLOT PCA (PC1 vs PC2)
# ---------------------------
png(output_file, width = 900, height = 700)

ggplot(pca, aes(x = PC1, y = PC2, color = FlightType)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +      # 95% confidence ellipse
  theme_minimal() +
  labs(
    title = "PCA of Columba livia domestica (Long vs Short Distance Flyers, Autosomes)",
    x = "PC1 (16.4%)",
    y = "PC2 (12.9%)"
  ) +
  theme(
    # Add solid axis lines
    axis.line = element_line(color = "black", linewidth = 0.8),

    # Increase font sizes
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 18, face = "bold"),

    legend.position = "right"
  ) +
  scale_color_brewer(palette = "Set1")

dev.off()

cat("PCA plot saved to", output_file, "\n")

