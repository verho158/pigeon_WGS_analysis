#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

pca_file <- "pigeons_sexchr_pca.eigenvec"
metadata_file <- "metadata.txt"
output_file <- "pca_sexchr_male_vs_female.png"

pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

pca$IID <- as.character(pca$IID)
metadata$IID <- as.character(metadata$IID)
pca <- left_join(pca, metadata, by = "IID")

pca <- pca %>% filter(Sex %in% c("Male", "Female"))

cat("Number of samples per sex:\n")
print(table(pca$Sex))

png(output_file, width = 900, height = 700)
ggplot(pca, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(title = "Sex Chromosomes PCA: Male vs Female",
       x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Set1")
dev.off()

cat("PCA plot saved to", output_file, "\n")
