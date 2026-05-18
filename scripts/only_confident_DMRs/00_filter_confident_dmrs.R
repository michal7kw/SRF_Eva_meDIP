#!/usr/bin/env Rscript

# Disable scientific notation for BED file output
options(scipen = 999)

################################################################################
# Script: 00_filter_confident_dmrs.R
# Purpose: Filter DMRs to retain only HIGH-CONFIDENCE regions where BOTH
#          TES and GFP samples have meaningful signal (>2 mean reads)
#
# Background:
#   The GFP meDIP samples have severe library quality issues:
#   - 85% duplication rate
#   - Only 4.3M unique molecules (vs 15M in TES)
#   - 65% of "hypermethylated" DMRs have GFP = 0 reads (artifact!)
#
#   When filtering for regions where BOTH samples have >2 reads:
#   - Original: 91% hyper / 9% hypo
#   - Filtered: 35% hyper / 65% hypo  <-- REVERSES the conclusion!
#
# Filtering Criteria:
#   1. FDR < 0.05 (already applied in input)
#   2. |logFC| > 1 (already applied in input)
#   3. mean_group1 > 2 (TES has meaningful signal)
#   4. mean_group2 > 2 (GFP has meaningful signal)  <-- KEY NEW FILTER
#
# Input:
#   - results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv
#
# Output:
#   - results/07_differential_MEDIPS_confident/TES_vs_GFP_DMRs_confident.csv
#   - results/07_differential_MEDIPS_confident/TES_vs_GFP_hypermethylated_confident.bed
#   - results/07_differential_MEDIPS_confident/TES_vs_GFP_hypomethylated_confident.bed
#   - results/07_differential_MEDIPS_confident/filtering_summary.txt
################################################################################

cat("=======================================================\n")
cat("High-Confidence DMR Filtering\n")
cat("=======================================================\n")
cat(paste("Start time:", Sys.time(), "\n\n"))

# Load libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

# Define paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
script_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs"
setwd(script_dir)

input_file <- file.path(base_dir, "results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv")
output_dir <- "output/07_differential_MEDIPS_confident"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define filtering thresholds
MIN_COUNT_THRESHOLD <- 2  # Minimum mean reads required in EACH sample

cat("=== Filtering Parameters ===\n")
cat(paste("  Minimum count per sample:", MIN_COUNT_THRESHOLD, "\n"))
cat(paste("  Input file:", input_file, "\n"))
cat(paste("  Output dir:", output_dir, "\n\n"))

# Read input DMRs
cat("=== Loading DMR Data ===\n")
dmr_data <- read_csv(input_file, show_col_types = FALSE)
cat(paste("  Total DMRs loaded:", nrow(dmr_data), "\n"))

# Show column names for debugging
cat("  Columns:", paste(colnames(dmr_data), collapse = ", "), "\n\n")

# Examine original distribution
cat("=== Original DMR Distribution ===\n")
original_hyper <- sum(dmr_data$logFC > 0)
original_hypo <- sum(dmr_data$logFC < 0)
cat(paste("  Hypermethylated (TES > GFP):", original_hyper,
          "(", round(100 * original_hyper / nrow(dmr_data), 1), "%)\n"))
cat(paste("  Hypomethylated (TES < GFP):", original_hypo,
          "(", round(100 * original_hypo / nrow(dmr_data), 1), "%)\n\n"))

# Examine count distribution
cat("=== Count Distribution in Original DMRs ===\n")
cat("  TES (mean_group1):\n")
cat(paste("    Mean:", round(mean(dmr_data$mean_group1), 2), "\n"))
cat(paste("    Median:", round(median(dmr_data$mean_group1), 2), "\n"))
cat(paste("    % with count = 0:", round(100 * sum(dmr_data$mean_group1 == 0) / nrow(dmr_data), 1), "%\n"))
cat(paste("    % with count <= 2:", round(100 * sum(dmr_data$mean_group1 <= 2) / nrow(dmr_data), 1), "%\n"))

cat("  GFP (mean_group2):\n")
cat(paste("    Mean:", round(mean(dmr_data$mean_group2), 2), "\n"))
cat(paste("    Median:", round(median(dmr_data$mean_group2), 2), "\n"))
cat(paste("    % with count = 0:", round(100 * sum(dmr_data$mean_group2 == 0) / nrow(dmr_data), 1), "%\n"))
cat(paste("    % with count <= 2:", round(100 * sum(dmr_data$mean_group2 <= 2) / nrow(dmr_data), 1), "%\n\n"))

# Apply high-confidence filter
cat("=== Applying High-Confidence Filter ===\n")
cat(paste("  Filter: mean_group1 >", MIN_COUNT_THRESHOLD,
          "AND mean_group2 >", MIN_COUNT_THRESHOLD, "\n"))

confident_dmrs <- dmr_data %>%
    filter(mean_group1 > MIN_COUNT_THRESHOLD,
           mean_group2 > MIN_COUNT_THRESHOLD)

cat(paste("  DMRs passing filter:", nrow(confident_dmrs), "\n"))
cat(paste("  DMRs removed:", nrow(dmr_data) - nrow(confident_dmrs), "\n"))
cat(paste("  Retention rate:", round(100 * nrow(confident_dmrs) / nrow(dmr_data), 1), "%\n\n"))

# Examine filtered distribution
cat("=== High-Confidence DMR Distribution ===\n")
confident_hyper <- sum(confident_dmrs$logFC > 0)
confident_hypo <- sum(confident_dmrs$logFC < 0)
cat(paste("  Hypermethylated (TES > GFP):", confident_hyper,
          "(", round(100 * confident_hyper / nrow(confident_dmrs), 1), "%)\n"))
cat(paste("  Hypomethylated (TES < GFP):", confident_hypo,
          "(", round(100 * confident_hypo / nrow(confident_dmrs), 1), "%)\n\n"))

# Highlight the reversal
cat("=== KEY FINDING: Direction Reversal ===\n")
cat("  Original data (includes GFP dropout artifacts):\n")
cat(paste("    Hypermethylated:", round(100 * original_hyper / nrow(dmr_data), 1), "%\n"))
cat(paste("    Hypomethylated:", round(100 * original_hypo / nrow(dmr_data), 1), "%\n"))
cat("  High-confidence data (both samples have signal):\n")
cat(paste("    Hypermethylated:", round(100 * confident_hyper / nrow(confident_dmrs), 1), "%\n"))
cat(paste("    Hypomethylated:", round(100 * confident_hypo / nrow(confident_dmrs), 1), "%\n\n"))

# Save filtered CSV
cat("=== Saving Output Files ===\n")

# Add "chr" prefix if not present (needed for compatibility with BigWig files)
confident_dmrs <- confident_dmrs %>%
    mutate(chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)))

output_csv <- file.path(output_dir, "TES_vs_GFP_DMRs_confident.csv")
write_csv(confident_dmrs, output_csv)
cat(paste("  Saved:", output_csv, "\n"))

# Create BED files for hyper and hypo

# Hypermethylated (logFC > 0, meaning TES > GFP)
hyper_dmrs <- confident_dmrs %>%
    filter(logFC > 0) %>%
    mutate(
        start = as.integer(start),
        stop = as.integer(stop),
        name = paste0(chr, ":", start, "-", stop),
        score = round(-log10(FDR) * 100),
        strand = "."
    ) %>%
    select(chr, start, stop, name, score, strand)

hyper_bed <- file.path(output_dir, "TES_vs_GFP_hypermethylated_confident.bed")
write.table(hyper_dmrs, hyper_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat(paste("  Saved:", hyper_bed, "(", nrow(hyper_dmrs), "regions)\n"))

# Hypomethylated (logFC < 0, meaning TES < GFP)
hypo_dmrs <- confident_dmrs %>%
    filter(logFC < 0) %>%
    mutate(
        start = as.integer(start),
        stop = as.integer(stop),
        name = paste0(chr, ":", start, "-", stop),
        score = round(-log10(FDR) * 100),
        strand = "."
    ) %>%
    select(chr, start, stop, name, score, strand)

hypo_bed <- file.path(output_dir, "TES_vs_GFP_hypomethylated_confident.bed")
write.table(hypo_dmrs, hypo_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat(paste("  Saved:", hypo_bed, "(", nrow(hypo_dmrs), "regions)\n"))

# Create all DMRs BED file
all_dmrs_bed <- confident_dmrs %>%
    mutate(
        start = as.integer(start),
        stop = as.integer(stop),
        name = paste0(chr, ":", start, "-", stop),
        score = round(-log10(FDR) * 100),
        strand = "."
    ) %>%
    select(chr, start, stop, name, score, strand)

all_bed <- file.path(output_dir, "TES_vs_GFP_all_confident.bed")
write.table(all_dmrs_bed, all_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat(paste("  Saved:", all_bed, "(", nrow(all_dmrs_bed), "regions)\n"))

# Save summary statistics
summary_file <- file.path(output_dir, "filtering_summary.txt")
sink(summary_file)
cat("=======================================================\n")
cat("HIGH-CONFIDENCE DMR FILTERING SUMMARY\n")
cat("=======================================================\n")
cat(paste("Generated:", Sys.time(), "\n\n"))

cat("FILTERING CRITERIA:\n")
cat(paste("  - FDR < 0.05 (pre-applied)\n"))
cat(paste("  - |logFC| > 1 (pre-applied)\n"))
cat(paste("  - mean_group1 (TES) >", MIN_COUNT_THRESHOLD, "\n"))
cat(paste("  - mean_group2 (GFP) >", MIN_COUNT_THRESHOLD, "\n\n"))

cat("ORIGINAL DMR COUNTS:\n")
cat(paste("  Total:", nrow(dmr_data), "\n"))
cat(paste("  Hypermethylated:", original_hyper, "(", round(100 * original_hyper / nrow(dmr_data), 1), "%)\n"))
cat(paste("  Hypomethylated:", original_hypo, "(", round(100 * original_hypo / nrow(dmr_data), 1), "%)\n\n"))

cat("HIGH-CONFIDENCE DMR COUNTS:\n")
cat(paste("  Total:", nrow(confident_dmrs), "\n"))
cat(paste("  Hypermethylated:", confident_hyper, "(", round(100 * confident_hyper / nrow(confident_dmrs), 1), "%)\n"))
cat(paste("  Hypomethylated:", confident_hypo, "(", round(100 * confident_hypo / nrow(confident_dmrs), 1), "%)\n\n"))

cat("RETENTION STATISTICS:\n")
cat(paste("  DMRs retained:", nrow(confident_dmrs), "/", nrow(dmr_data),
          "(", round(100 * nrow(confident_dmrs) / nrow(dmr_data), 1), "%)\n"))
cat(paste("  DMRs removed:", nrow(dmr_data) - nrow(confident_dmrs),
          "(", round(100 * (nrow(dmr_data) - nrow(confident_dmrs)) / nrow(dmr_data), 1), "%)\n\n"))

cat("KEY FINDING:\n")
cat("  The direction of differential methylation REVERSES when\n")
cat("  filtering for high-confidence DMRs with signal in both samples!\n")
cat(paste("  Original: ", round(100 * original_hyper / nrow(dmr_data), 1), "% hyper / ",
          round(100 * original_hypo / nrow(dmr_data), 1), "% hypo\n", sep = ""))
cat(paste("  Filtered: ", round(100 * confident_hyper / nrow(confident_dmrs), 1), "% hyper / ",
          round(100 * confident_hypo / nrow(confident_dmrs), 1), "% hypo\n", sep = ""))
cat("\n")
cat("INTERPRETATION:\n")
cat("  The apparent 91% hypermethylation in the original data is largely\n")
cat("  an artifact of GFP sample dropout (low library complexity).\n")
cat("  When both samples have meaningful signal, hypomethylation dominates.\n")
sink()
cat(paste("  Saved:", summary_file, "\n\n"))

cat("=======================================================\n")
cat("Filtering Complete!\n")
cat("=======================================================\n")
cat(paste("End time:", Sys.time(), "\n"))
