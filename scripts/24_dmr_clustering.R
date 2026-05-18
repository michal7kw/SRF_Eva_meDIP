#!/usr/bin/env Rscript

# ============================================================================
# DMR Clustering Analysis - R Statistics and Visualization
# ============================================================================
#
# Purpose: Analyze DMR clusters and create visualizations
#
# Input:
#   - cluster_stats_*.tsv files (from bash script)
#   - mega_dmrs.tsv (from bash script)
#
# Output:
#   - Cluster size distributions
#   - Comparison across distance thresholds
#   - Mega-DMR characterization
#   - Genome-wide ideogram (if possible)
#
# ============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(gridExtra)
})

cat("============================================\n")
cat("DMR Clustering Analysis - R Script\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP")

OUTDIR <- "results/24_dmr_clustering"
PLOTDIR <- file.path(OUTDIR, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load data
# ============================================================================

cat("=== Loading data ===\n")

# Load cluster statistics at different thresholds
distances <- c("500bp", "1kb", "2kb", "5kb")
cluster_data <- list()

for (dist in distances) {
    file_path <- file.path(OUTDIR, paste0("cluster_stats_", dist, ".tsv"))
    if (file.exists(file_path)) {
        cluster_data[[dist]] <- read.delim(file_path, stringsAsFactors = FALSE)
        cluster_data[[dist]]$distance_threshold <- dist
        cat("  Loaded", nrow(cluster_data[[dist]]), "clusters at", dist, "\n")
    }
}

# Combine all cluster data
all_clusters <- do.call(rbind, cluster_data)
all_clusters$distance_threshold <- factor(all_clusters$distance_threshold,
                                          levels = distances)

# Load mega-DMRs
mega_dmrs <- read.delim(file.path(OUTDIR, "mega_dmrs.tsv"),
                        stringsAsFactors = FALSE)
cat("  Loaded", nrow(mega_dmrs), "mega-DMRs\n\n")

# ============================================================================
# Summary statistics
# ============================================================================

cat("=== Calculating summary statistics ===\n")

# Summary by distance threshold
summary_by_dist <- all_clusters %>%
    group_by(distance_threshold) %>%
    summarise(
        n_clusters = n(),
        total_dmrs = sum(n_dmrs),
        mean_cluster_size = mean(n_dmrs),
        median_cluster_size = median(n_dmrs),
        max_cluster_size = max(n_dmrs),
        mean_span_bp = mean(span_bp),
        median_span_bp = median(span_bp),
        max_span_bp = max(span_bp),
        n_consistent = sum(direction_consistent == "yes"),
        pct_consistent = mean(direction_consistent == "yes") * 100,
        .groups = "drop"
    )

print(summary_by_dist)

write.csv(summary_by_dist, file.path(OUTDIR, "clustering_summary.csv"),
          row.names = FALSE)
cat("\n  Saved: clustering_summary.csv\n\n")

# ============================================================================
# Visualization 1: Cluster count by distance threshold
# ============================================================================

cat("=== Generating visualizations ===\n")

p1 <- ggplot(summary_by_dist, aes(x = distance_threshold, y = n_clusters)) +
    geom_bar(stat = "identity", fill = "#3182BD", color = "black", linewidth = 0.5) +
    geom_text(aes(label = n_clusters), vjust = -0.3, size = 4) +
    labs(title = "Number of DMR Clusters by Distance Threshold",
         subtitle = "Clusters formed by merging adjacent DMRs",
         x = "Merge Distance", y = "Number of Clusters") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(PLOTDIR, "cluster_count_by_distance.png"), p1,
       width = 8, height = 6, dpi = 300)
cat("  Saved: cluster_count_by_distance.png\n")

# ============================================================================
# Visualization 2: Cluster size distribution
# ============================================================================

# Use 1kb as default
clusters_1kb <- all_clusters %>% filter(distance_threshold == "1kb")

p2 <- ggplot(clusters_1kb, aes(x = n_dmrs)) +
    geom_histogram(binwidth = 1, fill = "#3182BD", color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 5, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = 5.5, y = Inf, label = "Mega-DMR threshold",
             vjust = 2, hjust = 0, color = "red", size = 3.5) +
    scale_x_continuous(breaks = seq(0, max(clusters_1kb$n_dmrs, na.rm = TRUE), by = 5)) +
    labs(title = "Distribution of DMRs per Cluster (1kb merge distance)",
         x = "Number of DMRs in Cluster", y = "Frequency") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(PLOTDIR, "cluster_size_distribution_1kb.png"), p2,
       width = 10, height = 6, dpi = 300)
cat("  Saved: cluster_size_distribution_1kb.png\n")

# ============================================================================
# Visualization 3: Cluster span distribution
# ============================================================================

p3 <- ggplot(clusters_1kb, aes(x = span_bp / 1000)) +
    geom_histogram(bins = 50, fill = "#3182BD", color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = 11, y = Inf, label = "Mega-DMR threshold",
             vjust = 2, hjust = 0, color = "red", size = 3.5) +
    scale_x_log10(labels = comma) +
    labs(title = "Distribution of Cluster Span (1kb merge distance)",
         x = "Cluster Span (kb, log scale)", y = "Frequency") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(PLOTDIR, "cluster_span_distribution_1kb.png"), p3,
       width = 10, height = 6, dpi = 300)
cat("  Saved: cluster_span_distribution_1kb.png\n")

# ============================================================================
# Visualization 4: Comparison across distance thresholds
# ============================================================================

# Median cluster size
p4a <- ggplot(summary_by_dist, aes(x = distance_threshold, y = median_cluster_size)) +
    geom_bar(stat = "identity", fill = "#9ECAE1", color = "black", linewidth = 0.5) +
    geom_text(aes(label = round(median_cluster_size, 1)), vjust = -0.3, size = 4) +
    labs(title = "Median Cluster Size",
         x = "Merge Distance", y = "Median DMRs per Cluster") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# % direction-consistent
p4b <- ggplot(summary_by_dist, aes(x = distance_threshold, y = pct_consistent)) +
    geom_bar(stat = "identity", fill = "#A1D99B", color = "black", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(pct_consistent, 1), "%")), vjust = -0.3, size = 4) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray40") +
    labs(title = "Direction Consistency",
         x = "Merge Distance", y = "% Clusters with Consistent Direction") +
    ylim(0, 100) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p4_combined <- grid.arrange(p4a, p4b, ncol = 2,
                            top = grid::textGrob("Cluster Properties by Merge Distance",
                                                 gp = grid::gpar(fontsize = 14, fontface = "bold")))

ggsave(file.path(PLOTDIR, "cluster_properties_by_distance.png"), p4_combined,
       width = 12, height = 5, dpi = 300)
cat("  Saved: cluster_properties_by_distance.png\n")

# ============================================================================
# Visualization 5: Mean logFC vs cluster size
# ============================================================================

p5 <- ggplot(clusters_1kb, aes(x = n_dmrs, y = mean_logFC)) +
    geom_point(aes(color = direction_consistent), alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_smooth(method = "loess", color = "red", se = FALSE, linewidth = 1) +
    scale_color_manual(values = c("yes" = "#2ECC40", "no" = "#FF4136"),
                       labels = c("yes" = "Consistent", "no" = "Mixed"),
                       name = "Direction") +
    labs(title = "Mean Methylation Change vs Cluster Size",
         subtitle = "1kb merge distance",
         x = "Number of DMRs in Cluster", y = "Mean log2(Fold Change)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right")

ggsave(file.path(PLOTDIR, "logFC_vs_cluster_size.png"), p5,
       width = 10, height = 7, dpi = 300)
cat("  Saved: logFC_vs_cluster_size.png\n")

# ============================================================================
# Visualization 6: Mega-DMR analysis
# ============================================================================

if (nrow(mega_dmrs) > 0) {
    # Mega-DMR characteristics
    mega_dmrs$qualifies_by <- factor(mega_dmrs$qualifies_by,
                                     levels = c("n_dmrs", "span", "n_dmrs;span"))

    p6a <- ggplot(mega_dmrs, aes(x = n_dmrs, y = span_bp / 1000, color = mean_logFC)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 10, linetype = "dashed", color = "gray60") +
        geom_vline(xintercept = 5, linetype = "dashed", color = "gray60") +
        scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                              midpoint = 0, name = "Mean logFC") +
        labs(title = "Mega-DMR Characteristics",
             subtitle = paste0("n = ", nrow(mega_dmrs), " mega-DMRs"),
             x = "Number of DMRs", y = "Span (kb)") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))

    ggsave(file.path(PLOTDIR, "mega_dmr_characteristics.png"), p6a,
           width = 9, height = 7, dpi = 300)
    cat("  Saved: mega_dmr_characteristics.png\n")

    # Mega-DMR logFC distribution
    p6b <- ggplot(mega_dmrs, aes(x = mean_logFC)) +
        geom_histogram(bins = 30, fill = "#3182BD", color = "black", linewidth = 0.3) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        labs(title = "Methylation Change in Mega-DMRs",
             x = "Mean log2(Fold Change)", y = "Frequency") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(file.path(PLOTDIR, "mega_dmr_logFC_distribution.png"), p6b,
           width = 8, height = 6, dpi = 300)
    cat("  Saved: mega_dmr_logFC_distribution.png\n")
} else {
    cat("  No mega-DMRs found, skipping mega-DMR visualizations\n")
}

# ============================================================================
# Visualization 7: Chromosome distribution of clusters
# ============================================================================

# Count clusters per chromosome
chr_dist <- clusters_1kb %>%
    mutate(chr = gsub("chr", "", chr)) %>%
    group_by(chr) %>%
    summarise(
        n_clusters = n(),
        n_mega = sum(n_dmrs >= 5 | span_bp >= 10000),
        total_dmrs = sum(n_dmrs),
        .groups = "drop"
    )

# Order chromosomes properly
chr_order <- c(as.character(1:22), "X", "Y", "M")
chr_dist$chr <- factor(chr_dist$chr, levels = chr_order)
chr_dist <- chr_dist %>% filter(!is.na(chr))

p7 <- ggplot(chr_dist, aes(x = chr, y = n_clusters)) +
    geom_bar(stat = "identity", fill = "#3182BD", color = "black", linewidth = 0.3) +
    geom_bar(aes(y = n_mega), stat = "identity", fill = "#D73027",
             color = "black", linewidth = 0.3) +
    labs(title = "Cluster Distribution by Chromosome",
         subtitle = "Blue = all clusters, Red = mega-DMRs",
         x = "Chromosome", y = "Number of Clusters") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(PLOTDIR, "cluster_distribution_by_chr.png"), p7,
       width = 12, height = 6, dpi = 300)
cat("  Saved: cluster_distribution_by_chr.png\n")

# ============================================================================
# Visualization 8: Direction consistency by cluster size
# ============================================================================

clusters_1kb$size_category <- cut(clusters_1kb$n_dmrs,
                                  breaks = c(0, 1, 2, 5, 10, Inf),
                                  labels = c("1", "2", "3-5", "6-10", ">10"),
                                  include.lowest = TRUE)

consistency_by_size <- clusters_1kb %>%
    group_by(size_category) %>%
    summarise(
        n_clusters = n(),
        pct_consistent = mean(direction_consistent == "yes") * 100,
        .groups = "drop"
    )

p8 <- ggplot(consistency_by_size, aes(x = size_category, y = pct_consistent)) +
    geom_bar(stat = "identity", fill = "#A1D99B", color = "black", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(pct_consistent, 1), "%\n(n=", n_clusters, ")")),
              vjust = -0.1, size = 3.5) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray40") +
    labs(title = "Direction Consistency by Cluster Size",
         subtitle = "Dashed line = 80% threshold for consistent direction",
         x = "DMRs per Cluster", y = "% Consistent Direction") +
    ylim(0, 105) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(PLOTDIR, "consistency_by_cluster_size.png"), p8,
       width = 9, height = 6, dpi = 300)
cat("  Saved: consistency_by_cluster_size.png\n")

# ============================================================================
# Create summary figure
# ============================================================================

p_summary <- grid.arrange(
    p1 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    p2 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    p5 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    p8 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    ncol = 2, nrow = 2,
    top = grid::textGrob("DMR Clustering Analysis Summary",
                         gp = grid::gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(PLOTDIR, "clustering_summary_figure.png"), p_summary,
       width = 14, height = 12, dpi = 300)
cat("  Saved: clustering_summary_figure.png\n")

# ============================================================================
# Key findings
# ============================================================================

cat("\n============================================\n")
cat("Analysis Complete\n")
cat("============================================\n\n")

cat("Key Findings:\n")
cat("-------------\n")

# Optimal distance threshold (based on direction consistency)
best_dist <- summary_by_dist %>%
    filter(pct_consistent >= 80) %>%
    arrange(desc(median_cluster_size)) %>%
    slice(1)

if (nrow(best_dist) > 0) {
    cat("  Recommended merge distance:", as.character(best_dist$distance_threshold), "\n")
    cat("    Direction consistency:", round(best_dist$pct_consistent, 1), "%\n")
    cat("    Median cluster size:", round(best_dist$median_cluster_size, 1), "DMRs\n")
}

cat("\n  At 1kb merge distance:\n")
cat("    Total clusters:", nrow(clusters_1kb), "\n")
cat("    Singleton clusters (1 DMR):", sum(clusters_1kb$n_dmrs == 1), "\n")
cat("    Multi-DMR clusters:", sum(clusters_1kb$n_dmrs > 1), "\n")
cat("    Mega-DMRs:", nrow(mega_dmrs), "\n")

if (nrow(mega_dmrs) > 0) {
    cat("\n  Mega-DMR summary:\n")
    cat("    Mean size:", round(mean(mega_dmrs$n_dmrs), 1), "DMRs\n")
    cat("    Mean span:", round(mean(mega_dmrs$span_bp) / 1000, 1), "kb\n")
    cat("    Mean logFC:", round(mean(mega_dmrs$mean_logFC), 3), "\n")
    cat("    % hypermethylated:", round(mean(mega_dmrs$pct_hyper), 1), "%\n")
}

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
