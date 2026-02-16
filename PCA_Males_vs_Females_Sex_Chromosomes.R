#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

# ---------------------------
# CONFIG PATHS
# ---------------------------
pca_file <- "pigeons_sexchr_pca.eigenvec"   # PCA output from sex chromosomes
metadata_file <- "metadata.txt"              # Must include IID and Sex
output_file <- "pca_sexchr_male_vs_female.png"

# ---------------------------
# LOAD DATA
# ---------------------------
pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# Ensure IID columns are character for merging
pca$IID <- as.character(pca$IID)
metadata$IID <- as.character(metadata$IID)

# ---------------------------
# FIX: Reconstruct full ID to match metadata
# ---------------------------
pca$IID <- paste0("NPO_", pca$IID)

# ---------------------------
# MERGE PCA WITH METADATA
# ---------------------------
pca <- left_join(pca, metadata, by = "IID")

# ---------------------------
# FILTER TO MALES AND FEMALES
# ---------------------------
pca <- pca %>% filter(Sex %in% c("Male", "Female"))

# ---------------------------
# CHECK SAMPLE COUNTS
# ---------------------------
cat("Number of samples per sex:\n")
print(table(pca$Sex))

# ---------------------------
# PLOT PCA (PC1 vs PC2)
# ---------------------------
png(output_file, width = 900, height = 700)

ggplot(pca, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +       # 95% confidence ellipse
  theme_minimal() +
  labs(
    title = "PCA of Columba livia domestica (Males vs Females, Sex Chromosomes)",
    x = "PC1 (22.6%)",
    y = "PC2 (12.8%)"
  ) +
  theme(
    # Solid axes
    axis.line = element_line(color = "black", linewidth = 0.8),

    # Font size improvements
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
