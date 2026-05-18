#!/usr/bin/env Rscript
#
# ENCODE ENHANCERS OF DEGs DOWN - HIGH-CONFIDENCE VERSION
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE enhancers specifically associated
#          with downregulated DEGs, stratified by:
#          1. TES/TEAD1 binding status
#          2. High-confidence DMR overlap (where BOTH samples have >2 reads)
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
})

options(scipen = 999)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")

cat("==========================================================\n")
cat("ENCODE ENHANCERS OF DEGs DOWN - HIGH-CONFIDENCE VERSION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# ENCODE enhancers (check multiple locations)
ENCODE_SOURCE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1/output/32_encode_enhancer/ENCODE_distal_enhancers.bed"
if (!file.exists(ENCODE_SOURCE)) {
    ENCODE_SOURCE <- "output/32_encode_enhancer_confident/ENCODE_distal_enhancers.bed"
}

# Peak files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# HIGH-CONFIDENCE DMRs
DMR_FILE <- "output/07_differential_MEDIPS_confident/TES_vs_GFP_DMRs_confident.csv"

# Output directory
OUTPUT_DIR <- "output/33_encode_enhancer_degs_down_confident"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD ENCODE ENHANCERS
# =============================================================================

cat("=== PHASE 1: Loading ENCODE Enhancers ===\n")

if (!file.exists(ENCODE_SOURCE)) {
    cat("  ERROR: ENCODE enhancers file not found.\n")
    cat("  Please run 32_encode_enhancer_methylation_confident.sh first.\n")
    quit(status = 1)
}

encode_raw <- read.table(ENCODE_SOURCE, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
cat(sprintf("  Raw ENCODE enhancers loaded: %d\n", nrow(encode_raw)))

encode_gr <- GRanges(
    seqnames = encode_raw$V1,
    ranges = IRanges(start = encode_raw$V2 + 1, end = encode_raw$V3)
)
encode_gr$accession <- encode_raw$V4

cat(sprintf("  ENCODE dELS enhancers: %d\n", length(encode_gr)))

# =============================================================================
# PHASE 2: ANNOTATE ENHANCERS TO NEAREST GENES
# =============================================================================

cat("\n=== PHASE 2: Annotating Enhancers to Nearest Genes ===\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cat("  Running ChIPseeker annotation...\n")
peak_anno <- annotatePeak(encode_gr, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
peak_anno_df <- as.data.frame(peak_anno)

cat(sprintf("  Annotated enhancers: %d\n", nrow(peak_anno_df)))

# =============================================================================
# PHASE 3: LOAD HIGH-CONFIDENCE DMRs
# =============================================================================

cat("\n=== PHASE 3: Loading High-Confidence DMRs ===\n")

dmr_data <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
cat(sprintf("  High-confidence DMRs loaded: %d\n", nrow(dmr_data)))

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

# =============================================================================
# PHASE 4: LOAD DESeq2 RESULTS
# =============================================================================

cat("\n=== PHASE 4: Loading DESeq2 Results ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total genes in DESeq2: %d\n", nrow(deseq2)))

degs_down_genes <- deseq2 %>%
    filter(padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_symbol)
cat(sprintf("  DEGs DOWN genes: %d\n", length(degs_down_genes)))

control_genes <- deseq2 %>%
    filter(!is.na(padj), padj >= 0.05) %>%
    pull(gene_symbol)
cat(sprintf("  Control genes (non-DE): %d\n", length(control_genes)))

# =============================================================================
# PHASE 5: FILTER ENHANCERS FOR DEGs DOWN
# =============================================================================

cat("\n=== PHASE 5: Filtering Enhancers for DEGs DOWN ===\n")

enhancers_degs_down <- peak_anno_df %>%
    filter(SYMBOL %in% degs_down_genes)
cat(sprintf("  Enhancers linked to DEGs DOWN: %d\n", nrow(enhancers_degs_down)))

enhancers_control_pool <- peak_anno_df %>%
    filter(SYMBOL %in% control_genes)
cat(sprintf("  Enhancers linked to Control genes: %d\n", nrow(enhancers_control_pool)))

# =============================================================================
# PHASE 6: LOAD TES/TEAD1 PEAKS AND OVERLAP WITH ENHANCERS
# =============================================================================

cat("\n=== PHASE 6: Overlapping with TES/TEAD1 Peaks and DMRs ===\n")

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

# Convert enhancers to GRanges
enhancers_degs_down_gr <- GRanges(
    seqnames = enhancers_degs_down$seqnames,
    ranges = IRanges(start = enhancers_degs_down$start, end = enhancers_degs_down$end)
)
enhancers_degs_down_gr$SYMBOL <- enhancers_degs_down$SYMBOL
enhancers_degs_down_gr$accession <- enhancers_degs_down$accession

# Find overlaps with peaks
tes_overlaps <- findOverlaps(enhancers_degs_down_gr, tes_gr)
tead1_overlaps <- findOverlaps(enhancers_degs_down_gr, tead1_gr)
tes_bound_idx <- unique(queryHits(tes_overlaps))
tead1_bound_idx <- unique(queryHits(tead1_overlaps))

# Find overlaps with DMRs
dmr_overlaps <- findOverlaps(enhancers_degs_down_gr, dmr_gr)
dmr_idx <- unique(queryHits(dmr_overlaps))

# Create categories
both_idx <- intersect(tes_bound_idx, tead1_bound_idx)
tes_only_idx <- setdiff(tes_bound_idx, tead1_bound_idx)
tead1_only_idx <- setdiff(tead1_bound_idx, tes_bound_idx)
unbound_idx <- setdiff(1:length(enhancers_degs_down_gr), union(tes_bound_idx, tead1_bound_idx))

# Cross binding with DMR status
tes_bound_dmr_idx <- intersect(tes_bound_idx, dmr_idx)
tes_bound_nodmr_idx <- setdiff(tes_bound_idx, dmr_idx)
unbound_dmr_idx <- intersect(unbound_idx, dmr_idx)
unbound_nodmr_idx <- setdiff(unbound_idx, dmr_idx)

cat(sprintf("\n  DEGs DOWN Enhancers by binding status:\n"))
cat(sprintf("    TES-bound: %d\n", length(tes_bound_idx)))
cat(sprintf("      - with DMR: %d (MOST RELIABLE)\n", length(tes_bound_dmr_idx)))
cat(sprintf("      - no DMR: %d\n", length(tes_bound_nodmr_idx)))
cat(sprintf("    Unbound: %d\n", length(unbound_idx)))
cat(sprintf("      - with DMR: %d (RELIABLE)\n", length(unbound_dmr_idx)))
cat(sprintf("      - no DMR: %d\n", length(unbound_nodmr_idx)))

# =============================================================================
# PHASE 7: CREATE CONTROL SET (MATCHED SIZE)
# =============================================================================

cat("\n=== PHASE 7: Creating Control Set ===\n")

n_target <- length(tes_bound_dmr_idx)
set.seed(42)

if (nrow(enhancers_control_pool) > n_target && n_target > 0) {
    control_sample_idx <- sample(1:nrow(enhancers_control_pool), n_target)
    enhancers_control <- enhancers_control_pool[control_sample_idx, ]
    cat(sprintf("  Subsampled Control enhancers: %d\n", nrow(enhancers_control)))
} else {
    enhancers_control <- enhancers_control_pool
    cat(sprintf("  Using all Control enhancers: %d\n", nrow(enhancers_control)))
}

# =============================================================================
# PHASE 8: EXPORT BED FILES
# =============================================================================

cat("\n=== PHASE 8: Exporting BED Files ===\n")

write_bed <- function(gr, filename) {
    if (length(gr) == 0) {
        cat(sprintf("  WARNING: No regions for %s\n", basename(filename)))
        return()
    }

    bed_df <- data.frame(
        chr = as.character(seqnames(gr)),
        start = as.integer(start(gr) - 1),
        end = as.integer(end(gr)),
        name = if (!is.null(gr$SYMBOL)) gr$SYMBOL else if (!is.null(gr$accession)) gr$accession else paste0("enh_", 1:length(gr)),
        score = 0,
        strand = "."
    )
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]
    write.table(bed_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s (%d regions)\n", basename(filename), nrow(bed_df)))
}

# DMR-stratified DEGs DOWN enhancers
write_bed(enhancers_degs_down_gr[tes_bound_dmr_idx], file.path(OUTPUT_DIR, "TES_bound_enhancers_DEGs_DOWN_with_DMR.bed"))
write_bed(enhancers_degs_down_gr[tes_bound_nodmr_idx], file.path(OUTPUT_DIR, "TES_bound_enhancers_DEGs_DOWN_no_DMR.bed"))
write_bed(enhancers_degs_down_gr[tes_bound_idx], file.path(OUTPUT_DIR, "TES_bound_enhancers_DEGs_DOWN_all.bed"))
write_bed(enhancers_degs_down_gr[tead1_only_idx], file.path(OUTPUT_DIR, "TEAD1_only_enhancers_DEGs_DOWN.bed"))
write_bed(enhancers_degs_down_gr[unbound_dmr_idx], file.path(OUTPUT_DIR, "Unbound_enhancers_DEGs_DOWN_with_DMR.bed"))
write_bed(enhancers_degs_down_gr[unbound_nodmr_idx], file.path(OUTPUT_DIR, "Unbound_enhancers_DEGs_DOWN_no_DMR.bed"))

# Control enhancers
control_gr <- GRanges(
    seqnames = enhancers_control$seqnames,
    ranges = IRanges(start = enhancers_control$start, end = enhancers_control$end)
)
control_gr$SYMBOL <- enhancers_control$SYMBOL
write_bed(control_gr, file.path(OUTPUT_DIR, "Control_enhancers.bed"))

# =============================================================================
# PHASE 9: SUMMARY STATISTICS
# =============================================================================

cat("\n=== PHASE 9: Summary Statistics ===\n")

summary_df <- data.frame(
    Category = c("TES_bound_DEGs_DOWN_with_DMR", "TES_bound_DEGs_DOWN_no_DMR",
                 "TES_bound_DEGs_DOWN_all", "TEAD1_only_DEGs_DOWN",
                 "Unbound_DEGs_DOWN_with_DMR", "Unbound_DEGs_DOWN_no_DMR", "Control"),
    N_enhancers = c(length(tes_bound_dmr_idx), length(tes_bound_nodmr_idx),
                    length(tes_bound_idx), length(tead1_only_idx),
                    length(unbound_dmr_idx), length(unbound_nodmr_idx),
                    nrow(enhancers_control)),
    Reliability = c("HIGH", "LOW (GFP dropout possible)", "MIXED",
                    "MEDIUM", "HIGH", "LOW (GFP dropout possible)", "MIXED")
)

print(summary_df)
write.csv(summary_df, file.path(OUTPUT_DIR, "enhancer_summary_dmr.csv"), row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("PREPARATION COMPLETE (HIGH-CONFIDENCE VERSION)\n")
cat("==========================================================\n")

cat("\nKey insight:\n")
cat("  Enhancers with DMR overlap have RELIABLE methylation data.\n")
cat("  Regions WITHOUT DMR overlap may show false hypermethylation\n")
cat("  due to GFP library dropout (65% of original 'hyper' DMRs had GFP=0).\n")
