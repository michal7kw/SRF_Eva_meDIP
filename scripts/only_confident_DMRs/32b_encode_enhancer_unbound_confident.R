#!/usr/bin/env Rscript
#
# ENCODE ENHANCER METHYLATION - UNBOUND ENHANCERS VERSION (32b)
# =============================================================================
#
# Purpose: Prepare ENCODE enhancer regions stratified by:
#   1. TES/TEAD1 binding status
#   2. High-confidence DMR overlap (where BOTH samples have >2 reads)
#
# This version is identical to 32 but outputs to a different directory
# for the unbound-focused analysis (which has more regions with DMR overlap).
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

# Disable scientific notation for BED output
options(scipen = 999)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")

cat("==========================================================\n")
cat("ENCODE ENHANCER PREPARATION - UNBOUND VERSION (32b)\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# ENCODE enhancers (from standard download step in original script 32)
# Check if exists, if not we'll need the original script to run first
ENCODE_ENHANCERS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1/output/32_encode_enhancer/ENCODE_distal_enhancers.bed"

# If not available there, use alternative path
if (!file.exists(ENCODE_ENHANCERS)) {
    ENCODE_ENHANCERS <- "output/32_encode_enhancer_confident/ENCODE_distal_enhancers.bed"
}

# Peak files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# HIGH-CONFIDENCE DMRs
DMR_FILE <- "output/07_differential_MEDIPS_confident/TES_vs_GFP_DMRs_confident.csv"

# Output directory - 32b version
OUTPUT_DIR <- "output/32b_encode_enhancer_unbound_confident"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD ENCODE ENHANCERS
# =============================================================================

cat("=== PHASE 1: Loading ENCODE Enhancers ===\n")

if (!file.exists(ENCODE_ENHANCERS)) {
    cat("  ERROR: ENCODE enhancers file not found at:\n")
    cat(sprintf("    %s\n", ENCODE_ENHANCERS))
    cat("  Please run the original 32_encode_enhancer_methylation.sh first\n")
    cat("  to download and prepare ENCODE cCRE annotations.\n")
    quit(status = 1)
}

encode_raw <- read.table(ENCODE_ENHANCERS, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
cat(sprintf("  Raw ENCODE enhancers loaded: %d\n", nrow(encode_raw)))

encode_gr <- GRanges(
    seqnames = encode_raw$V1,
    ranges = IRanges(start = encode_raw$V2 + 1, end = encode_raw$V3)
)
encode_gr$accession <- encode_raw$V4

cat(sprintf("  ENCODE dELS enhancers: %d\n", length(encode_gr)))

# =============================================================================
# PHASE 2: LOAD TES AND TEAD1 PEAKS
# =============================================================================

cat("\n=== PHASE 2: Loading TES/TEAD1 Peaks ===\n")

load_peaks <- function(peak_file) {
    peaks <- read.table(peak_file,
        header = FALSE, stringsAsFactors = FALSE,
        col.names = c("chr", "start", "end", "name", "score",
                      "strand", "signalValue", "pValue", "qValue", "peak"))
    if (!grepl("^chr", peaks$chr[1])) peaks$chr <- paste0("chr", peaks$chr)
    GRanges(seqnames = peaks$chr, ranges = IRanges(peaks$start + 1, peaks$end))
}

tes_gr <- load_peaks(TES_PEAKS)
tead1_gr <- load_peaks(TEAD1_PEAKS)

cat(sprintf("  TES peaks: %d\n", length(tes_gr)))
cat(sprintf("  TEAD1 peaks: %d\n", length(tead1_gr)))

# =============================================================================
# PHASE 3: LOAD HIGH-CONFIDENCE DMRs
# =============================================================================

cat("\n=== PHASE 3: Loading High-Confidence DMRs ===\n")

dmr_data <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
cat(sprintf("  High-confidence DMRs loaded: %d\n", nrow(dmr_data)))
cat(sprintf("    (These have >2 reads in BOTH TES and GFP samples)\n"))

# Add "chr" prefix to DMR chromosomes if not present
dmr_chr <- as.character(dmr_data$chr)
if (!grepl("^chr", dmr_chr[1])) {
    dmr_chr <- paste0("chr", dmr_chr)
}

dmr_gr <- GRanges(
    seqnames = dmr_chr,
    ranges = IRanges(start = dmr_data$start, end = dmr_data$stop),
    logFC = dmr_data$logFC
)

hyper_dmrs <- dmr_gr[dmr_gr$logFC > 0]
hypo_dmrs <- dmr_gr[dmr_gr$logFC < 0]
cat(sprintf("  Hypermethylated: %d (%.1f%%)\n", length(hyper_dmrs), 100*length(hyper_dmrs)/length(dmr_gr)))
cat(sprintf("  Hypomethylated: %d (%.1f%%)\n", length(hypo_dmrs), 100*length(hypo_dmrs)/length(dmr_gr)))

# =============================================================================
# PHASE 4: OVERLAP ENHANCERS WITH PEAKS AND DMRs
# =============================================================================

cat("\n=== PHASE 4: Overlapping Enhancers with Peaks and DMRs ===\n")

# Peak overlaps
tes_overlaps <- findOverlaps(encode_gr, tes_gr)
tead1_overlaps <- findOverlaps(encode_gr, tead1_gr)
tes_bound_idx <- unique(queryHits(tes_overlaps))
tead1_bound_idx <- unique(queryHits(tead1_overlaps))

cat(sprintf("  Enhancers overlapping TES peaks: %d (%.1f%%)\n",
    length(tes_bound_idx), 100*length(tes_bound_idx)/length(encode_gr)))
cat(sprintf("  Enhancers overlapping TEAD1 peaks: %d (%.1f%%)\n",
    length(tead1_bound_idx), 100*length(tead1_bound_idx)/length(encode_gr)))

# DMR overlaps
dmr_overlaps <- findOverlaps(encode_gr, dmr_gr)
dmr_bound_idx <- unique(queryHits(dmr_overlaps))

cat(sprintf("  Enhancers overlapping high-conf DMRs: %d (%.1f%%)\n",
    length(dmr_bound_idx), 100*length(dmr_bound_idx)/length(encode_gr)))

# Create binding categories
both_idx <- intersect(tes_bound_idx, tead1_bound_idx)
tes_only_idx <- setdiff(tes_bound_idx, tead1_bound_idx)
tead1_only_idx <- setdiff(tead1_bound_idx, tes_bound_idx)
unbound_idx <- setdiff(1:length(encode_gr), union(tes_bound_idx, tead1_bound_idx))

# Cross with DMR status
tes_bound_dmr <- intersect(tes_bound_idx, dmr_bound_idx)
tes_bound_nodmr <- setdiff(tes_bound_idx, dmr_bound_idx)
unbound_dmr <- intersect(unbound_idx, dmr_bound_idx)
unbound_nodmr <- setdiff(unbound_idx, dmr_bound_idx)

cat(sprintf("\n  Category breakdown:\n"))
cat(sprintf("    TES-bound WITH DMR: %d (%.1f%%)\n",
    length(tes_bound_dmr), 100*length(tes_bound_dmr)/length(encode_gr)))
cat(sprintf("    TES-bound NO DMR: %d (%.1f%%)\n",
    length(tes_bound_nodmr), 100*length(tes_bound_nodmr)/length(encode_gr)))
cat(sprintf("    Unbound WITH DMR: %d (%.1f%%)\n",
    length(unbound_dmr), 100*length(unbound_dmr)/length(encode_gr)))
cat(sprintf("    Unbound NO DMR: %d (%.1f%%)\n",
    length(unbound_nodmr), 100*length(unbound_nodmr)/length(encode_gr)))

# =============================================================================
# PHASE 5: EXPORT BED FILES
# =============================================================================

cat("\n=== PHASE 5: Exporting BED Files ===\n")

write_bed <- function(gr, filename) {
    if (length(gr) == 0) {
        cat(sprintf("  WARNING: No regions for %s (skipping)\n", basename(filename)))
        return(invisible(NULL))
    }

    bed_df <- data.frame(
        chr = as.character(seqnames(gr)),
        start = as.integer(start(gr) - 1),
        end = as.integer(end(gr)),
        name = if (!is.null(gr$accession)) gr$accession else paste0("enh_", seq_len(length(gr))),
        score = rep(0, length(gr)),
        strand = rep(".", length(gr))
    )
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]
    write.table(bed_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s (%d regions)\n", basename(filename), nrow(bed_df)))
}

# DMR-stratified exports
write_bed(encode_gr[tes_bound_dmr], file.path(OUTPUT_DIR, "TES_bound_enhancers_with_DMR.bed"))
write_bed(encode_gr[tes_bound_nodmr], file.path(OUTPUT_DIR, "TES_bound_enhancers_no_DMR.bed"))
write_bed(encode_gr[tes_bound_idx], file.path(OUTPUT_DIR, "TES_bound_enhancers_all.bed"))
write_bed(encode_gr[tead1_only_idx], file.path(OUTPUT_DIR, "TEAD1_only_enhancers.bed"))
write_bed(encode_gr[unbound_dmr], file.path(OUTPUT_DIR, "Unbound_enhancers_with_DMR.bed"))
write_bed(encode_gr[unbound_nodmr], file.path(OUTPUT_DIR, "Unbound_enhancers_no_DMR.bed"))

# Subsample unbound NO DMR for comparison (match Unbound WITH DMR count)
set.seed(42)
n_sample <- length(unbound_dmr)
if (length(unbound_nodmr) > n_sample) {
    unbound_nodmr_subsample <- sample(unbound_nodmr, n_sample)
} else {
    unbound_nodmr_subsample <- unbound_nodmr
}
write_bed(encode_gr[unbound_nodmr_subsample], file.path(OUTPUT_DIR, "Unbound_enhancers_no_DMR_subsampled.bed"))

# =============================================================================
# PHASE 6: SUMMARY STATISTICS
# =============================================================================

cat("\n=== PHASE 6: Summary Statistics ===\n")

summary_df <- data.frame(
    Category = c("TES_bound_with_DMR", "TES_bound_no_DMR", "TES_bound_all",
                 "TEAD1_only", "Unbound_with_DMR", "Unbound_no_DMR", "Total"),
    N_enhancers = c(
        length(tes_bound_dmr), length(tes_bound_nodmr), length(tes_bound_idx),
        length(tead1_only_idx), length(unbound_dmr), length(unbound_nodmr),
        length(encode_gr)
    ),
    Percentage = c(
        100*length(tes_bound_dmr)/length(encode_gr),
        100*length(tes_bound_nodmr)/length(encode_gr),
        100*length(tes_bound_idx)/length(encode_gr),
        100*length(tead1_only_idx)/length(encode_gr),
        100*length(unbound_dmr)/length(encode_gr),
        100*length(unbound_nodmr)/length(encode_gr),
        100
    ),
    Note = c(
        "MOST RELIABLE for methylation",
        "May have GFP dropout issues",
        "All TES-bound",
        "TEAD1 but not TES",
        "RELIABLE for methylation - PRIMARY FOCUS",
        "May have GFP dropout issues",
        "Total ENCODE dELS enhancers"
    )
)

print(summary_df)

write.csv(summary_df, file.path(OUTPUT_DIR, "enhancer_binding_dmr_summary.csv"), row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("PREPARATION COMPLETE - UNBOUND VERSION (32b)\n")
cat("==========================================================\n")

cat("\nKey insight:\n")
cat("  This version focuses on UNBOUND enhancers with DMR overlap\n")
cat("  which provides more regions (454) for reliable analysis.\n\n")
cat("  TES-bound enhancers with DMR overlap are very few (~3),\n")
cat("  so use the Unbound comparison for robust statistics.\n")
