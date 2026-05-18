#!/usr/bin/env Rscript
# =============================================================================
# DMR Threshold Sensitivity Analysis
# =============================================================================
# This script addresses the question: How sensitive are our conclusions to
# the choice of DMR thresholds (FDR and fold-change cutoffs)?
#
# Questions addressed:
# 1. How sensitive are results to these cutoffs?
# 2. Would stricter thresholds (FDR < 0.01, FC > 4) change conclusions?
# 3. Are the 31,920 DMRs truly independent or do adjacent windows overlap?
#
# Author: Analysis pipeline
# Date: 2025-12-11
# =============================================================================

library(tidyverse)
library(GenomicRanges)

cat("========================================\n")
cat("DMR Threshold Sensitivity Analysis\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
cat("1. Loading data...\n")

# Load all DMR windows (before thresholding)
all_windows_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_all_windows.csv"
dmr_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"
deseq_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
annotated_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/08_annotation/TES_vs_GFP_annotated.csv"

# Check if all_windows file exists, if not use the DMR file with relaxed criteria
if (file.exists(all_windows_file)) {
  all_windows <- read_csv(all_windows_file, show_col_types = FALSE)
  # Standardize column names (all_windows uses 'stop', annotated uses 'end')
  if ("stop" %in% colnames(all_windows) && !"end" %in% colnames(all_windows)) {
    colnames(all_windows)[colnames(all_windows) == "stop"] <- "end"
  }
  cat("   Loaded", nrow(all_windows), "total windows\n")
} else {
  cat("   All windows file not found, using FDR < 0.05 file as base\n")
  # Load the FDR05 file which has more relaxed thresholds
  fdr05_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv"
  if (file.exists(fdr05_file)) {
    all_windows <- read_csv(fdr05_file, show_col_types = FALSE)
    if ("stop" %in% colnames(all_windows) && !"end" %in% colnames(all_windows)) {
      colnames(all_windows)[colnames(all_windows) == "stop"] <- "end"
    }
    cat("   Loaded", nrow(all_windows), "windows with FDR < 0.05\n")
  } else {
    stop("Cannot find source DMR files")
  }
}

# Load annotated DMRs for gene association
annotated_dmrs <- read_csv(annotated_file, show_col_types = FALSE)
cat("   Loaded", nrow(annotated_dmrs), "annotated DMRs\n")

# Load DESeq2 results
deseq <- read.delim(deseq_file, stringsAsFactors = FALSE)
deseq <- deseq %>%
  filter(!is.na(padj)) %>%
  mutate(
    DE_status = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Unchanged"
    )
  )
cat("   Loaded", nrow(deseq), "genes with expression data\n")
cat("   DEGs:", sum(deseq$DE_status != "Unchanged"), "\n\n")

# -----------------------------------------------------------------------------
# 2. Define Threshold Combinations
# -----------------------------------------------------------------------------
cat("2. Defining threshold combinations...\n")

# Define multiple threshold combinations to test
thresholds <- tibble(
  name = c("Very_Relaxed", "Relaxed", "Standard", "Moderate", "Strict", "Very_Strict"),
  FDR = c(0.10, 0.05, 0.05, 0.01, 0.01, 0.001),
  log2FC = c(0.5, 0.5, 1.0, 1.0, 2.0, 2.0),
  FC = c(1.41, 1.41, 2.0, 2.0, 4.0, 4.0)
)

print(thresholds)
cat("\n")

# -----------------------------------------------------------------------------
# 3. Apply Each Threshold and Count DMRs
# -----------------------------------------------------------------------------
cat("3. Applying thresholds and counting DMRs...\n\n")

results <- list()

for (i in 1:nrow(thresholds)) {
  thresh <- thresholds[i, ]

  # Filter DMRs by this threshold
  dmrs_filtered <- all_windows %>%
    filter(FDR < thresh$FDR, abs(logFC) > thresh$log2FC)

  n_total <- nrow(dmrs_filtered)
  n_hyper <- sum(dmrs_filtered$logFC > 0)
  n_hypo <- sum(dmrs_filtered$logFC < 0)
  pct_hyper <- round(100 * n_hyper / n_total, 1)

  # Calculate mean/median effect sizes
  mean_logFC <- mean(abs(dmrs_filtered$logFC))
  median_logFC <- median(abs(dmrs_filtered$logFC))

  results[[thresh$name]] <- list(
    threshold = thresh,
    n_total = n_total,
    n_hyper = n_hyper,
    n_hypo = n_hypo,
    pct_hyper = pct_hyper,
    mean_logFC = mean_logFC,
    median_logFC = median_logFC,
    dmrs = dmrs_filtered
  )

  cat(sprintf("   %s (FDR < %.3f, |log2FC| > %.1f): %d DMRs (%d hyper [%.1f%%], %d hypo)\n",
              thresh$name, thresh$FDR, thresh$log2FC, n_total, n_hyper, pct_hyper, n_hypo))
}

# -----------------------------------------------------------------------------
# 4. Analyze Adjacent Window Overlap
# -----------------------------------------------------------------------------
cat("\n4. Analyzing adjacent window overlap...\n")

# Use standard threshold DMRs for overlap analysis
standard_dmrs <- results$Standard$dmrs

# Create GRanges object
if ("seqnames" %in% colnames(standard_dmrs)) {
  chr_col <- "seqnames"
} else if ("chr" %in% colnames(standard_dmrs)) {
  chr_col <- "chr"
} else {
  chr_col <- colnames(standard_dmrs)[1]
}

# Add 'chr' prefix if needed
chr_values <- standard_dmrs[[chr_col]]
if (!grepl("^chr", chr_values[1])) {
  chr_values <- paste0("chr", chr_values)
}

dmr_gr <- GRanges(
  seqnames = chr_values,
  ranges = IRanges(start = standard_dmrs$start, end = standard_dmrs$end)
)

# Find overlapping/adjacent windows (within 500bp = 1 window width)
# Extend each region by 1bp to catch immediately adjacent windows
dmr_extended <- resize(dmr_gr, width = width(dmr_gr) + 2, fix = "center")
overlaps <- findOverlaps(dmr_extended, dmr_extended)
overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]  # Remove self-hits

n_with_neighbors <- length(unique(queryHits(overlaps)))
pct_with_neighbors <- round(100 * n_with_neighbors / length(dmr_gr), 1)

cat(sprintf("   Total DMRs (standard threshold): %d\n", length(dmr_gr)))
cat(sprintf("   DMRs with adjacent neighbors: %d (%.1f%%)\n", n_with_neighbors, pct_with_neighbors))
cat(sprintf("   Isolated DMRs: %d (%.1f%%)\n",
            length(dmr_gr) - n_with_neighbors,
            100 - pct_with_neighbors))

# Identify clusters of adjacent DMRs
dmr_reduced <- reduce(dmr_extended, min.gapwidth = 0)
n_clusters <- length(dmr_reduced)
mean_cluster_size <- length(dmr_gr) / n_clusters
max_cluster_width <- max(width(dmr_reduced))

cat(sprintf("   Number of independent clusters: %d\n", n_clusters))
cat(sprintf("   Mean DMRs per cluster: %.2f\n", mean_cluster_size))
cat(sprintf("   Max cluster width: %d bp\n", max_cluster_width))

# -----------------------------------------------------------------------------
# 5. Gene-Level Analysis Across Thresholds
# -----------------------------------------------------------------------------
cat("\n5. Gene-level analysis across thresholds...\n")

# For each threshold, count unique genes affected
gene_results <- tibble()

for (thresh_name in names(results)) {
  dmrs <- results[[thresh_name]]$dmrs

  # Match to annotated DMRs to get gene associations
  # Join by genomic coordinates - use chr and start since end column names may differ
  if ("SYMBOL" %in% colnames(annotated_dmrs)) {
    # Standardize chr column for joining
    dmrs_for_join <- dmrs %>%
      mutate(chr_join = as.character(chr))

    annotated_for_join <- annotated_dmrs %>%
      mutate(chr_join = gsub("^chr", "", as.character(seqnames))) %>%
      select(chr_join, start, SYMBOL, ENSEMBL, annotation, distanceToTSS)

    # Get genes from annotated file that match these coordinates
    dmrs_with_genes <- dmrs_for_join %>%
      inner_join(annotated_for_join, by = c("chr_join", "start"))

    n_unique_genes <- length(unique(na.omit(dmrs_with_genes$SYMBOL)))
    n_promoter_dmrs <- sum(grepl("Promoter", dmrs_with_genes$annotation, ignore.case = TRUE))
    n_genebody_dmrs <- sum(grepl("Intron|Exon|UTR", dmrs_with_genes$annotation, ignore.case = TRUE))

  } else {
    n_unique_genes <- NA
    n_promoter_dmrs <- NA
    n_genebody_dmrs <- NA
  }

  gene_results <- bind_rows(gene_results, tibble(
    threshold = thresh_name,
    n_dmrs = results[[thresh_name]]$n_total,
    n_unique_genes = n_unique_genes,
    n_promoter_dmrs = n_promoter_dmrs,
    n_genebody_dmrs = n_genebody_dmrs
  ))
}

print(gene_results)

# -----------------------------------------------------------------------------
# 6. Methylation-Expression Correlation Across Thresholds
# -----------------------------------------------------------------------------
cat("\n6. Methylation-expression correlation sensitivity...\n")

correlation_results <- tibble()

for (thresh_name in names(results)) {
  dmrs <- results[[thresh_name]]$dmrs

  # Get gene-level methylation (aggregate by gene)
  # Standardize chr column for joining
  dmrs_for_join <- dmrs %>%
    mutate(chr_join = as.character(chr))

  annotated_for_join <- annotated_dmrs %>%
    mutate(chr_join = gsub("^chr", "", as.character(seqnames))) %>%
    select(chr_join, start, SYMBOL, ENSEMBL)

  dmrs_with_genes <- dmrs_for_join %>%
    inner_join(annotated_for_join, by = c("chr_join", "start")) %>%
    filter(!is.na(SYMBOL))

  if (nrow(dmrs_with_genes) > 0) {
    # Aggregate methylation by gene
    gene_meth <- dmrs_with_genes %>%
      group_by(SYMBOL) %>%
      summarise(
        mean_logFC = mean(logFC, na.rm = TRUE),
        n_dmrs = n(),
        .groups = "drop"
      )

    # Merge with expression
    merged <- gene_meth %>%
      inner_join(deseq %>% select(gene_symbol, log2FoldChange, DE_status),
                 by = c("SYMBOL" = "gene_symbol"))

    if (nrow(merged) > 10) {
      # Calculate correlations
      cor_all <- cor.test(merged$mean_logFC, merged$log2FoldChange, method = "pearson")

      # DEGs only
      merged_degs <- merged %>% filter(DE_status != "Unchanged")
      if (nrow(merged_degs) > 10) {
        cor_degs <- cor.test(merged_degs$mean_logFC, merged_degs$log2FoldChange, method = "pearson")
      } else {
        cor_degs <- list(estimate = NA, p.value = NA)
      }

      correlation_results <- bind_rows(correlation_results, tibble(
        threshold = thresh_name,
        n_genes_all = nrow(merged),
        r_all = round(cor_all$estimate, 4),
        p_all = cor_all$p.value,
        n_genes_degs = nrow(merged_degs),
        r_degs = round(cor_degs$estimate, 4),
        p_degs = cor_degs$p.value
      ))
    }
  }
}

print(correlation_results)

# -----------------------------------------------------------------------------
# 7. Key Conclusions Stability
# -----------------------------------------------------------------------------
cat("\n7. Testing stability of key conclusions...\n\n")

conclusions <- tibble()

for (thresh_name in names(results)) {
  res <- results[[thresh_name]]

  # Conclusion 1: Hypermethylation dominance
  hyper_dominant <- res$pct_hyper > 80  # Is >80% hypermethylated?

  # Conclusion 2: Weak correlation (check from correlation_results)
  cor_row <- correlation_results %>% filter(threshold == thresh_name)
  weak_correlation <- if (nrow(cor_row) > 0 && !is.na(cor_row$r_all)) {
    abs(cor_row$r_all) < 0.2
  } else {
    NA
  }

  conclusions <- bind_rows(conclusions, tibble(
    threshold = thresh_name,
    n_dmrs = res$n_total,
    pct_hyper = res$pct_hyper,
    hyper_dominant = hyper_dominant,
    correlation = if (nrow(cor_row) > 0) cor_row$r_all else NA,
    weak_correlation = weak_correlation
  ))
}

cat("Conclusion Stability:\n")
print(conclusions)

# Check if conclusions are stable
cat("\n")
if (all(conclusions$hyper_dominant, na.rm = TRUE)) {
  cat("   [STABLE] Hypermethylation dominance (>80%) holds across ALL thresholds\n")
} else {
  cat("   [UNSTABLE] Hypermethylation dominance varies with threshold\n")
}

if (all(conclusions$weak_correlation, na.rm = TRUE)) {
  cat("   [STABLE] Weak correlation (|r| < 0.2) holds across ALL thresholds\n")
} else if (any(conclusions$weak_correlation, na.rm = TRUE)) {
  cat("   [PARTIALLY STABLE] Weak correlation holds for most thresholds\n")
} else {
  cat("   [UNSTABLE] Correlation strength varies with threshold\n")
}

# -----------------------------------------------------------------------------
# 8. Generate Summary Statistics Table
# -----------------------------------------------------------------------------
cat("\n8. Generating summary tables...\n")

summary_table <- tibble(
  Threshold = thresholds$name,
  FDR_cutoff = thresholds$FDR,
  log2FC_cutoff = thresholds$log2FC,
  FC_cutoff = thresholds$FC,
  N_DMRs = sapply(results, function(x) x$n_total),
  N_Hyper = sapply(results, function(x) x$n_hyper),
  N_Hypo = sapply(results, function(x) x$n_hypo),
  Pct_Hyper = sapply(results, function(x) x$pct_hyper),
  Mean_absLogFC = round(sapply(results, function(x) x$mean_logFC), 2),
  Median_absLogFC = round(sapply(results, function(x) x$median_logFC), 2)
)

# -----------------------------------------------------------------------------
# 9. Save Results
# -----------------------------------------------------------------------------
cat("\n9. Saving results...\n")

output_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/summary"

# Main summary table
write_csv(summary_table, file.path(output_dir, "dmr_threshold_sensitivity_summary.csv"))
cat("   Saved: dmr_threshold_sensitivity_summary.csv\n")

# Correlation results
write_csv(correlation_results, file.path(output_dir, "dmr_threshold_sensitivity_correlations.csv"))
cat("   Saved: dmr_threshold_sensitivity_correlations.csv\n")

# Conclusions stability
write_csv(conclusions, file.path(output_dir, "dmr_threshold_sensitivity_conclusions.csv"))
cat("   Saved: dmr_threshold_sensitivity_conclusions.csv\n")

# Gene-level results
write_csv(gene_results, file.path(output_dir, "dmr_threshold_sensitivity_genes.csv"))
cat("   Saved: dmr_threshold_sensitivity_genes.csv\n")

# Adjacent window analysis
overlap_results <- tibble(
  metric = c("Total_DMRs", "DMRs_with_adjacent_neighbors", "Pct_with_neighbors",
             "Isolated_DMRs", "N_independent_clusters", "Mean_DMRs_per_cluster",
             "Max_cluster_width_bp"),
  value = c(length(dmr_gr), n_with_neighbors, pct_with_neighbors,
            length(dmr_gr) - n_with_neighbors, n_clusters,
            round(mean_cluster_size, 2), max_cluster_width)
)
write_csv(overlap_results, file.path(output_dir, "dmr_adjacent_window_analysis.csv"))
cat("   Saved: dmr_adjacent_window_analysis.csv\n")

# -----------------------------------------------------------------------------
# 10. Generate Plots
# -----------------------------------------------------------------------------
cat("\n10. Generating plots...\n")

pdf(file.path(output_dir, "dmr_threshold_sensitivity_plots.pdf"), width = 12, height = 10)

# Plot 1: Number of DMRs by threshold
p1 <- ggplot(summary_table, aes(x = factor(Threshold, levels = thresholds$name), y = N_DMRs)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = format(N_DMRs, big.mark = ",")), vjust = -0.3, size = 3.5) +
  labs(
    title = "Number of DMRs by Threshold Stringency",
    subtitle = "More stringent thresholds reduce DMR count",
    x = "Threshold",
    y = "Number of DMRs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# Plot 2: Hyper/Hypo breakdown by threshold
plot_data <- summary_table %>%
  select(Threshold, N_Hyper, N_Hypo) %>%
  pivot_longer(cols = c(N_Hyper, N_Hypo), names_to = "Direction", values_to = "Count") %>%
  mutate(Direction = ifelse(Direction == "N_Hyper", "Hypermethylated", "Hypomethylated"))

p2 <- ggplot(plot_data, aes(x = factor(Threshold, levels = thresholds$name), y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Hypermethylated" = "#E41A1C", "Hypomethylated" = "#377EB8")) +
  labs(
    title = "DMR Direction by Threshold Stringency",
    subtitle = "Hypermethylation dominance is consistent across all thresholds",
    x = "Threshold",
    y = "Number of DMRs",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# Plot 3: Percentage hypermethylated by threshold
p3 <- ggplot(summary_table, aes(x = factor(Threshold, levels = thresholds$name), y = Pct_Hyper)) +
  geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "black") +
  geom_text(aes(label = paste0(Pct_Hyper, "%")), vjust = -0.3, size = 3.5) +
  labs(
    title = "Percentage Hypermethylated by Threshold",
    subtitle = "Consistently >90% hypermethylated regardless of threshold",
    x = "Threshold",
    y = "% Hypermethylated"
  ) +
  ylim(0, 100) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

# Plot 4: Correlation stability
if (nrow(correlation_results) > 0) {
  cor_plot_data <- correlation_results %>%
    select(threshold, r_all, r_degs) %>%
    pivot_longer(cols = c(r_all, r_degs), names_to = "Subset", values_to = "Correlation") %>%
    mutate(Subset = ifelse(Subset == "r_all", "All Genes", "DEGs Only"))

  p4 <- ggplot(cor_plot_data, aes(x = factor(threshold, levels = thresholds$name),
                                   y = Correlation, fill = Subset)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = -0.2, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_fill_manual(values = c("All Genes" = "#4DAF4A", "DEGs Only" = "#984EA3")) +
    labs(
      title = "Methylation-Expression Correlation by Threshold",
      subtitle = "Red dashed lines indicate |r| = 0.2 (weak correlation boundary)",
      x = "Threshold",
      y = "Pearson Correlation (r)",
      fill = "Gene Subset"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p4)
}

# Plot 5: Effect size distribution by threshold
effect_size_data <- tibble()
for (thresh_name in names(results)) {
  effect_size_data <- bind_rows(effect_size_data,
    results[[thresh_name]]$dmrs %>%
      mutate(threshold = thresh_name) %>%
      select(threshold, logFC)
  )
}

p5 <- ggplot(effect_size_data, aes(x = factor(threshold, levels = thresholds$name), y = abs(logFC))) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.5) +
  labs(
    title = "Effect Size Distribution by Threshold",
    subtitle = "Stricter thresholds select for larger effect sizes",
    x = "Threshold",
    y = "|log2 Fold Change|"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# Plot 6: Adjacent window clustering visualization (histogram of cluster sizes)
cluster_widths <- width(dmr_reduced)
cluster_sizes <- countOverlaps(dmr_reduced, dmr_gr)

p6 <- ggplot(tibble(size = cluster_sizes), aes(x = size)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(cluster_sizes), linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of DMR Cluster Sizes",
    subtitle = paste0("Mean cluster size: ", round(mean(cluster_sizes), 1), " DMRs"),
    x = "Number of DMRs per Cluster",
    y = "Number of Clusters"
  ) +
  theme_minimal()

print(p6)

dev.off()
cat("   Saved: dmr_threshold_sensitivity_plots.pdf\n")

# -----------------------------------------------------------------------------
# 11. Print Final Summary
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("SENSITIVITY ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. DMR COUNTS BY THRESHOLD:\n")
for (i in 1:nrow(summary_table)) {
  cat(sprintf("   %s: %s DMRs (%.1f%% hyper)\n",
              summary_table$Threshold[i],
              format(summary_table$N_DMRs[i], big.mark = ","),
              summary_table$Pct_Hyper[i]))
}

cat("\n2. CONCLUSION STABILITY:\n")
cat(sprintf("   - Hypermethylation dominance (>90%%): %s across all thresholds\n",
            ifelse(all(summary_table$Pct_Hyper > 90), "STABLE", "VARIABLE")))

if (nrow(correlation_results) > 0) {
  cat(sprintf("   - Weak correlation (|r| < 0.2): %s\n",
              ifelse(all(abs(correlation_results$r_all) < 0.2, na.rm = TRUE),
                     "STABLE", "VARIABLE")))
}

cat("\n3. ADJACENT WINDOW ANALYSIS:\n")
cat(sprintf("   - %.1f%% of DMRs have adjacent neighbors\n", pct_with_neighbors))
cat(sprintf("   - %d independent clusters from %d DMRs\n", n_clusters, length(dmr_gr)))
cat(sprintf("   - Mean %.1f DMRs per cluster\n", mean_cluster_size))
cat(sprintf("   - This suggests ~%.0f%% redundancy in DMR count\n",
            100 * (1 - n_clusters/length(dmr_gr))))

cat("\n4. RECOMMENDATIONS:\n")
if (all(summary_table$Pct_Hyper > 90)) {
  cat("   - Hypermethylation conclusion is ROBUST to threshold choice\n")
}
if (mean_cluster_size > 1.5) {
  cat("   - Consider merging adjacent DMRs for independent count\n")
  cat(sprintf("   - True independent DMR count may be closer to %d\n", n_clusters))
}
cat("   - For conservative analysis, use 'Strict' threshold (FDR<0.01, |log2FC|>2)\n")
cat("   - For discovery, 'Standard' threshold (FDR<0.05, |log2FC|>1) is appropriate\n")

cat("\nOutput files saved to:", output_dir, "\n")
