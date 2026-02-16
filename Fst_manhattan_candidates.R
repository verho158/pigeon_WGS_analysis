# fst_manhattan_candidates.R
# Generate a Manhattan plot of Weir & Cockerham FST values with autosome labels
# Also extract candidate SNPs based on top FST values

library(data.table)
library(ggplot2)

# ---- Step 1: Load FST data ----
fst_file <- "long_vs_short.weir.fst"
fst <- fread(fst_file)

# Remove missing or negative FST values
fst <- fst[!is.na(WEIR_AND_COCKERHAM_FST) & WEIR_AND_COCKERHAM_FST >= 0]

# ---- Step 2: Define chromosome order and sizes ----
chrom_order <- c(
  "NC_088602.1","NC_088603.1","NC_088604.1","NC_088605.1","NC_088606.1",
  "NC_088607.1","NC_088608.1","NC_088609.1","NC_088610.1","NC_088611.1",
  "NC_088612.1","NC_088613.1","NC_088614.1","NC_088615.1","NC_088616.1",
  "NC_088617.1","NC_088618.1","NC_088619.1","NC_088620.1","NC_088621.1",
  "NC_088622.1","NC_088623.1","NC_088624.1","NC_088625.1","NC_088626.1",
  "NC_088627.1","NC_088628.1","NC_088629.1","NC_088630.1","NC_088631.1",
  "NC_088632.1","NC_088633.1","NC_088634.1","NC_088635.1","NC_088636.1",
  "NC_088637.1","NC_088638.1","NC_088639.1","NC_088640.1","NC_088641.1","NC_088642.1"
)

chrom_sizes <- c(
  212386202, 163726572, 122092291, 78855516, 68980706, 41865730, 41451919,
  34522999, 29946398, 23263835, 22349698, 21496342, 21485835, 20651115,
  17800093, 15524004, 14482990, 14042120, 11817371, 11100612, 9002469,
  7344926, 6871264, 6290118, 6129807, 5959470, 4741841, 3632745, 3051461,
  1359805, 1272592, 689683, 564255, 538117, 489914, 459057, 302372, 276054,
  204328, 238678768, 84824678
)

# ---- Step 2b: Keep only autosomes ----
autosomes <- chrom_order[1:39]       # only true autosomes (exclude W/Z)
fst <- fst[CHROM %in% autosomes]

# ---- Step 3: Compute cumulative positions ----
fst$CHROM <- factor(fst$CHROM, levels = autosomes)
cumsum_sizes <- cumsum(c(0, chrom_sizes[1:39][-length(autosomes)]))
names(cumsum_sizes) <- autosomes
fst[, BPcum := POS + cumsum_sizes[as.character(CHROM)]]

axisdf <- fst[, .(center = mean(BPcum)), by = CHROM]
axisdf$label <- paste0("Chr", 1:39)  # relabel for plotting

# ---- Step 4: Select candidate SNPs (top 1%) ----
fst_threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, 0.99)
candidates <- fst[WEIR_AND_COCKERHAM_FST >= fst_threshold]

# Save candidate SNPs to CSV
fwrite(candidates, "candidate_SNPs_top1pct.csv")

cat("Number of candidate SNPs (top 1%):", nrow(candidates), "\n")
cat("Threshold FST for top 1%:", round(fst_threshold, 3), "\n")

# ---- Step 5: Generate Manhattan plot ----
manhattan_plot <- ggplot(fst, aes(x = BPcum, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(aes(color = CHROM), alpha = 0.6, size = 0.8) +
  geom_point(data = candidates, color = "red", size = 1) + # highlight top SNPs
  scale_x_continuous(label = axisdf$label, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = rep(c("grey30","grey70"), length.out = length(autosomes))) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # --- Improved aesthetics ---
    axis.line = element_line(color = "black", linewidth = 0.8),  # solid axes
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  xlab("Chromosome") +
  ylab("Weir & Cockerham FST") +
  ggtitle("FST Manhattan Plot (Autosomes Only) with Candidate SNPs Highlighted")

# Save plot
ggsave("FST_Manhattan_candidates_autosomes.png", plot = manhattan_plot, width = 12, height = 6)
