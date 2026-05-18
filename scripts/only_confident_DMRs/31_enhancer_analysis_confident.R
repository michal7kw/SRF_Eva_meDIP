#!/usr/bin/env Rscript
#
# ENHANCER ANALYSIS: DEGs DOWN vs CONTROL - HIGH-CONFIDENCE VERSION
# =============================================================================
#
# Purpose: Define "Enhancer" regions as TES/TEAD1 peaks that are NOT in promoters.
#          Stratify by expression status AND high-confidence DMR overlap.
#
# Key insight: High-confidence DMRs show 65% hypomethylation (vs 91% hyper in
# original analysis). Only regions with high-conf DMRs have reliable methylation.
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

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")

cat("==========================================================\n")
cat("ENHANCER ANALYSIS - HIGH-CONFIDENCE VERSION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# HIGH-CONFIDENCE DMRs
DMR_FILE <- "output/07_differential_MEDIPS_confident/TES_vs_GFP_DMRs_confident.csv"

# Output directory
OUTPUT_BASE <- "output/31_enhancer_analysis_confident"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD AND ANNOTATE PEAKS
# =============================================================================

cat("=== PHASE 1: Loading and Annotating Peaks ===\n")

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

all_peaks <- unique(c(tes_gr, tead1_gr))
cat(sprintf("  Total unique TES/TEAD1 peaks: %d\n", length(all_peaks)))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cat("  Annotating peaks to nearest genes...\n")
peak_anno <- annotatePeak(all_peaks, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
peak_anno_df <- as.data.frame(peak_anno)

# Filter for Enhancers (Exclude Promoters)
cat("  Annotation categories found:\n")
print(table(gsub(" \\(.*", "", peak_anno_df$annotation)))

enhancer_idx <- !grepl("^Promoter", peak_anno_df$annotation)
enhancer_peaks <- peak_anno_df[enhancer_idx, ]

cat(sprintf("\n  Defined 'Enhancers' (Non-Promoter Peaks): %d\n", nrow(enhancer_peaks)))

# =============================================================================
# PHASE 2: LOAD HIGH-CONFIDENCE DMRs
# =============================================================================

cat("\n=== PHASE 2: Loading High-Confidence DMRs ===\n")

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

# Find enhancers overlapping DMRs
enhancer_gr <- GRanges(
    seqnames = enhancer_peaks$seqnames,
    ranges = IRanges(start = enhancer_peaks$start, end = enhancer_peaks$end)
)

dmr_overlaps <- findOverlaps(enhancer_gr, dmr_gr)
enhancers_with_dmr_idx <- unique(queryHits(dmr_overlaps))
enhancers_no_dmr_idx <- setdiff(1:length(enhancer_gr), enhancers_with_dmr_idx)

cat(sprintf("  Enhancers overlapping high-conf DMRs: %d (%.1f%%)\n",
            length(enhancers_with_dmr_idx), 100*length(enhancers_with_dmr_idx)/length(enhancer_gr)))
cat(sprintf("  Enhancers without high-conf DMR: %d\n", length(enhancers_no_dmr_idx)))

# =============================================================================
# PHASE 3: INTEGRATE RNA-SEQ (DESEQ2)
# =============================================================================

cat("\n=== PHASE 3: Integrating RNA-seq Data ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total genes in DESeq2: %d\n", nrow(deseq2)))

degs_down_genes <- deseq2 %>%
    filter(padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_symbol)

control_genes <- deseq2 %>%
    filter(!is.na(padj), padj >= 0.05) %>%
    pull(gene_symbol)

cat(sprintf("  DEGs DOWN genes: %d\n", length(degs_down_genes)))
cat(sprintf("  Control genes (Expressed, Non-DE): %d\n", length(control_genes)))

# =============================================================================
# PHASE 4: ASSIGN ENHANCERS TO GROUPS (WITH DMR STRATIFICATION)
# =============================================================================

cat("\n=== PHASE 4: Assigning Enhancers to Groups ===\n")

# Enhancers of DEGs DOWN - with and without DMR
enhancers_degs_down <- enhancer_peaks %>%
    filter(SYMBOL %in% degs_down_genes)
cat(sprintf("  Enhancers assigned to DEGs DOWN: %d\n", nrow(enhancers_degs_down)))

# Create GRanges for DEGs DOWN enhancers
enhancers_degs_down_gr <- GRanges(
    seqnames = enhancers_degs_down$seqnames,
    ranges = IRanges(start = enhancers_degs_down$start, end = enhancers_degs_down$end)
)

# Find which DEGs DOWN enhancers have DMR overlap
dmr_overlaps_down <- findOverlaps(enhancers_degs_down_gr, dmr_gr)
degs_down_with_dmr_idx <- unique(queryHits(dmr_overlaps_down))
degs_down_no_dmr_idx <- setdiff(1:nrow(enhancers_degs_down), degs_down_with_dmr_idx)

enhancers_degs_down_dmr <- enhancers_degs_down[degs_down_with_dmr_idx, ]
enhancers_degs_down_nodmr <- enhancers_degs_down[degs_down_no_dmr_idx, ]

cat(sprintf("    DEGs DOWN enhancers WITH DMR: %d\n", nrow(enhancers_degs_down_dmr)))
cat(sprintf("    DEGs DOWN enhancers NO DMR: %d\n", nrow(enhancers_degs_down_nodmr)))

# Enhancers of Control Genes
enhancers_control_pool <- enhancer_peaks %>%
    filter(SYMBOL %in% control_genes)
cat(sprintf("  Enhancers assigned to Control Genes (Pool): %d\n", nrow(enhancers_control_pool)))

# Match sample size to DEGs DOWN with DMR
n_target <- nrow(enhancers_degs_down_dmr)
set.seed(42)

if (nrow(enhancers_control_pool) > n_target && n_target > 0) {
    enhancers_control <- enhancers_control_pool[sample(nrow(enhancers_control_pool), n_target), ]
    cat(sprintf("  Subsampled Control Enhancers to match N: %d\n", nrow(enhancers_control)))
} else {
    enhancers_control <- enhancers_control_pool
    cat("  Using all Control enhancers\n")
}

# =============================================================================
# PHASE 5: EXPORT BED FILES
# =============================================================================

cat("\n=== PHASE 5: Exporting BED Files ===\n")

write_bed <- function(df, filename) {
    if (nrow(df) == 0) {
        cat(sprintf("  WARNING: No regions for %s\n", basename(filename)))
        return()
    }

    bed_df <- data.frame(
        chr = df$seqnames,
        start = df$start - 1,
        end = df$end,
        name = df$SYMBOL,
        score = 0,
        strand = df$strand
    )
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]
    write.table(bed_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s (%d regions)\n", basename(filename), nrow(bed_df)))
}

# Export DMR-stratified DEGs DOWN enhancers
write_bed(enhancers_degs_down_dmr, file.path(OUTPUT_BASE, "Enhancers_DEGs_DOWN_with_DMR.bed"))
write_bed(enhancers_degs_down_nodmr, file.path(OUTPUT_BASE, "Enhancers_DEGs_DOWN_no_DMR.bed"))
write_bed(enhancers_degs_down, file.path(OUTPUT_BASE, "Enhancers_DEGs_DOWN_all.bed"))
write_bed(enhancers_control, file.path(OUTPUT_BASE, "Enhancers_Random_Control.bed"))

# Save summary stats
summary_df <- data.frame(
    Group = c("Enhancers_DEGs_DOWN_with_DMR", "Enhancers_DEGs_DOWN_no_DMR",
              "Enhancers_DEGs_DOWN_all", "Enhancers_Random_Control"),
    N_Peaks = c(nrow(enhancers_degs_down_dmr), nrow(enhancers_degs_down_nodmr),
                nrow(enhancers_degs_down), nrow(enhancers_control)),
    Description = c(
        "Non-promoter peaks linked to DEGs DOWN, WITH high-conf DMR (MOST RELIABLE)",
        "Non-promoter peaks linked to DEGs DOWN, NO high-conf DMR",
        "All non-promoter peaks linked to DEGs DOWN",
        "Non-promoter peaks linked to Non-DE genes (matched control)"
    )
)
write.csv(summary_df, file.path(OUTPUT_BASE, "enhancer_counts.csv"), row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("ENHANCER PREPARATION COMPLETE (HIGH-CONFIDENCE VERSION)\n")
cat("==========================================================\n")

cat("\nKey insight:\n")
cat("  Enhancers WITH high-confidence DMR overlap have reliable methylation\n")
cat("  signals because BOTH TES and GFP samples have >2 reads.\n")
cat("  The GFP library has quality issues causing 65% of 'hyper' DMRs\n")
cat("  to have GFP=0 reads (artifact).\n")
