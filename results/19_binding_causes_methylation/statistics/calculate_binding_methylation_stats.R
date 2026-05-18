#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

cat("Calculating methylation statistics at binding sites...\n\n")

# ============================================================================
# Load data
# ============================================================================

# Load DMRs
dmr_file <- "results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv"
if (file.exists(dmr_file)) {
    dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)
    dmrs <- dmrs %>%
        mutate(
            dmr_id = paste(chr, start, end, sep = "_"),
            direction = ifelse(logFC > 0, "hyper", "hypo")
        )
    cat(sprintf("Loaded %d DMRs\n", nrow(dmrs)))
} else {
    stop("DMR file not found")
}

# Function to load BED file
load_bed <- function(filepath) {
    if (!file.exists(filepath)) return(NULL)
    df <- read.table(filepath, stringsAsFactors = FALSE)
    colnames(df)[1:3] <- c("chr", "start", "end")
    df$peak_id <- paste(df$chr, df$start, df$end, sep = "_")
    return(df)
}

# Load binding site categories
beds_dir <- file.path(outdir, "beds")
tes_all <- load_bed(file.path(beds_dir, "TES_all_peaks.bed"))
tead1_all <- load_bed(file.path(beds_dir, "TEAD1_all_peaks.bed"))
tes_unique <- load_bed(file.path(beds_dir, "TES_unique_peaks.bed"))
tead1_unique <- load_bed(file.path(beds_dir, "TEAD1_unique_peaks.bed"))
shared <- load_bed(file.path(beds_dir, "shared_peaks.bed"))
tes_with_dmr <- load_bed(file.path(beds_dir, "TES_peaks_with_DMR.bed"))
tes_no_dmr <- load_bed(file.path(beds_dir, "TES_peaks_no_DMR.bed"))

# ============================================================================
# Calculate statistics
# ============================================================================

stats_list <- list()

add_stat <- function(category, metric, value) {
    stats_list[[length(stats_list) + 1]] <<- data.frame(
        Category = category,
        Metric = metric,
        Value = as.character(value),
        stringsAsFactors = FALSE
    )
}

# Overview
add_stat("Overview", "Total_TES_peaks", nrow(tes_all))
add_stat("Overview", "Total_TEAD1_peaks", nrow(tead1_all))
add_stat("Overview", "Total_DMRs", nrow(dmrs))
add_stat("Overview", "Hypermethylated_DMRs", sum(dmrs$direction == "hyper"))
add_stat("Overview", "Hypomethylated_DMRs", sum(dmrs$direction == "hypo"))

# Binding site categories
add_stat("Binding_Categories", "TES_unique_peaks", nrow(tes_unique))
add_stat("Binding_Categories", "TEAD1_unique_peaks", nrow(tead1_unique))
add_stat("Binding_Categories", "Shared_peaks", nrow(shared))

# DMR overlap
if (!is.null(tes_with_dmr) && !is.null(tes_no_dmr)) {
    n_tes_dmr <- nrow(tes_with_dmr)
    n_tes_no_dmr <- nrow(tes_no_dmr)
    pct_tes_dmr <- round(100 * n_tes_dmr / nrow(tes_all), 1)

    add_stat("DMR_Overlap", "TES_peaks_with_DMR", n_tes_dmr)
    add_stat("DMR_Overlap", "TES_peaks_without_DMR", n_tes_no_dmr)
    add_stat("DMR_Overlap", "Pct_TES_peaks_with_DMR", paste0(pct_tes_dmr, "%"))

    # TEAD1 DMR overlap (if available)
    tead1_with_dmr <- load_bed(file.path(beds_dir, "TEAD1_peaks_with_DMR.bed"))
    if (!is.null(tead1_with_dmr)) {
        n_tead1_dmr <- nrow(tead1_with_dmr)
        pct_tead1_dmr <- round(100 * n_tead1_dmr / nrow(tead1_all), 1)
        add_stat("DMR_Overlap", "TEAD1_peaks_with_DMR", n_tead1_dmr)
        add_stat("DMR_Overlap", "Pct_TEAD1_peaks_with_DMR", paste0(pct_tead1_dmr, "%"))
    }
}

# Fisher's exact test: Is TES binding enriched at DMRs compared to TEAD1?
# This tests whether TES specifically recruits methylation machinery

# We need genome-wide background for proper enrichment
# Simplified: compare TES vs TEAD1 DMR overlap rates
if (!is.null(tes_with_dmr) && !is.null(tead1_all)) {
    tead1_with_dmr <- load_bed(file.path(beds_dir, "TEAD1_peaks_with_DMR.bed"))
    if (!is.null(tead1_with_dmr)) {
        # Contingency table
        tes_dmr <- nrow(tes_with_dmr)
        tes_no <- nrow(tes_all) - tes_dmr
        tead1_dmr <- nrow(tead1_with_dmr)
        tead1_no <- nrow(tead1_all) - tead1_dmr

        contingency <- matrix(c(tes_dmr, tes_no, tead1_dmr, tead1_no), nrow = 2,
                              dimnames = list(DMR = c("Yes", "No"), TF = c("TES", "TEAD1")))

        fisher_result <- fisher.test(contingency)

        add_stat("Enrichment_Test", "Fisher_OR_TES_vs_TEAD1_at_DMRs", round(fisher_result$estimate, 3))
        add_stat("Enrichment_Test", "Fisher_pvalue", format(fisher_result$p.value, scientific = TRUE, digits = 3))
        add_stat("Enrichment_Test", "Fisher_95CI_low", round(fisher_result$conf.int[1], 3))
        add_stat("Enrichment_Test", "Fisher_95CI_high", round(fisher_result$conf.int[2], 3))

        cat("\n=== Fisher's Exact Test: TES vs TEAD1 DMR Overlap ===\n")
        cat(sprintf("  TES peaks with DMR: %d / %d (%.1f%%)\n", tes_dmr, nrow(tes_all), 100*tes_dmr/nrow(tes_all)))
        cat(sprintf("  TEAD1 peaks with DMR: %d / %d (%.1f%%)\n", tead1_dmr, nrow(tead1_all), 100*tead1_dmr/nrow(tead1_all)))
        cat(sprintf("  Odds Ratio: %.3f (95%% CI: %.3f - %.3f)\n",
                    fisher_result$estimate, fisher_result$conf.int[1], fisher_result$conf.int[2]))
        cat(sprintf("  p-value: %s\n", format(fisher_result$p.value, scientific = TRUE, digits = 3)))

        if (fisher_result$estimate > 1 && fisher_result$p.value < 0.05) {
            cat("  Interpretation: TES binding is ENRICHED at DMRs compared to TEAD1\n")
        } else if (fisher_result$estimate < 1 && fisher_result$p.value < 0.05) {
            cat("  Interpretation: TES binding is DEPLETED at DMRs compared to TEAD1\n")
        } else {
            cat("  Interpretation: No significant difference in DMR overlap\n")
        }
    }
}

# ============================================================================
# Write output
# ============================================================================

all_stats <- do.call(rbind, stats_list)
output_file <- file.path(outdir, "statistics", "binding_methylation_statistics.csv")
write.csv(all_stats, output_file, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", output_file))

# Summary
cat("\n=== SUMMARY ===\n")
cat(sprintf("TES binding sites: %d\n", nrow(tes_all)))
cat(sprintf("  - TES-unique: %d\n", nrow(tes_unique)))
cat(sprintf("  - Shared with TEAD1: %d\n", nrow(shared)))
if (!is.null(tes_with_dmr)) {
    cat(sprintf("  - Overlapping DMRs: %d (%.1f%%)\n",
                nrow(tes_with_dmr), 100*nrow(tes_with_dmr)/nrow(tes_all)))
}
cat(sprintf("\nTEAD1 binding sites: %d\n", nrow(tead1_all)))
cat(sprintf("  - TEAD1-unique: %d\n", nrow(tead1_unique)))

