#!/usr/bin/env Rscript
#
# PEAKED GENES METAGENE PROFILES - HIGH-CONFIDENCE VERSION
# =============================================================================
#
# Creates BED files with gene body coordinates for deepTools scale-regions
# Includes genes with Cut&Tag peaks AND integrates high-confidence DMR info.
#
# Version 1: All genes with peaks (stratified by DMR status)
# Version 2: DEGs with peaks (UP/DOWN stratified by DMR status)
#
# Key insight: High-confidence DMRs show 65% hypomethylation vs 91% hyper
# in original analysis. This reversal is due to GFP library quality issues.
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(dplyr)
    library(rtracklayer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")

cat("==========================================================\n")
cat("PEAKED GENES METAGENE - HIGH-CONFIDENCE VERSION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak annotation files
TES_PEAKS_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv"
TEAD1_PEAKS_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# HIGH-CONFIDENCE DMRs
DMR_FILE <- "output/07_differential_MEDIPS_confident/TES_vs_GFP_DMRs_confident.csv"

# Output directory
OUTPUT_BASE <- "output/30_peaked_genes_metagene_confident"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

GENE_BODY_EXTENSION <- 5000  # bp extension for DMR overlap

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS FROM GTF
# =============================================================================

cat("=== PHASE 1: Loading Gene Annotations ===\n")

gtf <- rtracklayer::import(GTF_FILE)
genes_gtf <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes_gtf)

genes_df <- genes_df %>%
    filter(gene_type == "protein_coding") %>%
    filter(seqnames %in% paste0("chr", c(1:22, "X", "Y")))

cat(sprintf("  Protein-coding genes on standard chromosomes: %d\n", nrow(genes_df)))

genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

genes_lookup <- genes_df %>%
    select(gene_name, seqnames, start, end, strand, gene_id_clean) %>%
    distinct(gene_name, .keep_all = TRUE)

cat(sprintf("  Unique gene symbols with coordinates: %d\n", nrow(genes_lookup)))

# Create GRanges for genes
genes_gr <- GRanges(
    seqnames = genes_lookup$seqnames,
    ranges = IRanges(start = genes_lookup$start, end = genes_lookup$end),
    strand = genes_lookup$strand,
    gene_name = genes_lookup$gene_name,
    gene_id = genes_lookup$gene_id_clean
)

# Create extended gene regions for DMR overlap
extended_genes_gr <- GRanges(
    seqnames = seqnames(genes_gr),
    ranges = IRanges(
        start = pmax(1, start(genes_gr) - GENE_BODY_EXTENSION),
        end = end(genes_gr) + GENE_BODY_EXTENSION
    ),
    strand = strand(genes_gr),
    gene_name = genes_gr$gene_name,
    gene_id = genes_gr$gene_id
)

cat("\n")

# =============================================================================
# PHASE 2: LOAD PEAK ANNOTATIONS AND EXTRACT PEAKED GENES
# =============================================================================

cat("=== PHASE 2: Loading Peak Annotations ===\n")

tes_peaks <- read.csv(TES_PEAKS_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total TES peaks: %d\n", nrow(tes_peaks)))

tead1_peaks <- read.csv(TEAD1_PEAKS_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total TEAD1 peaks: %d\n", nrow(tead1_peaks)))

# Extract unique gene IDs from peak annotations (Entrez IDs)
tes_entrez_ids <- unique(tes_peaks$geneId[!is.na(tes_peaks$geneId)])
tead1_entrez_ids <- unique(tead1_peaks$geneId[!is.na(tead1_peaks$geneId)])

# Convert Entrez IDs to gene symbols
tes_symbols <- mapIds(org.Hs.eg.db,
                      keys = as.character(tes_entrez_ids),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")
tes_symbols <- unique(tes_symbols[!is.na(tes_symbols)])

tead1_symbols <- mapIds(org.Hs.eg.db,
                        keys = as.character(tead1_entrez_ids),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
tead1_symbols <- unique(tead1_symbols[!is.na(tead1_symbols)])

cat(sprintf("  Gene symbols with TES peaks: %d\n", length(tes_symbols)))
cat(sprintf("  Gene symbols with TEAD1 peaks: %d\n", length(tead1_symbols)))

cat("\n")

# =============================================================================
# PHASE 3: LOAD HIGH-CONFIDENCE DMRs
# =============================================================================

cat("=== PHASE 3: Loading High-Confidence DMRs ===\n")

dmr_data <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
cat(sprintf("  High-confidence DMRs loaded: %d\n", nrow(dmr_data)))
cat(sprintf("    (These DMRs have >2 reads in BOTH TES and GFP samples)\n"))

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

# Separate hyper and hypo
hyper_dmrs <- dmr_gr[dmr_gr$logFC > 0]
hypo_dmrs <- dmr_gr[dmr_gr$logFC < 0]
cat(sprintf("  Hypermethylated (TES > GFP): %d (%.1f%%)\n",
            length(hyper_dmrs), 100*length(hyper_dmrs)/length(dmr_gr)))
cat(sprintf("  Hypomethylated (TES < GFP): %d (%.1f%%)\n",
            length(hypo_dmrs), 100*length(hypo_dmrs)/length(dmr_gr)))

# Find genes with DMR overlap
dmr_overlaps <- findOverlaps(extended_genes_gr, dmr_gr)
dmr_gene_names <- unique(genes_gr$gene_name[queryHits(dmr_overlaps)])
cat(sprintf("  Genes with high-confidence DMR overlap: %d\n", length(dmr_gene_names)))

cat("\n")

# =============================================================================
# PHASE 4: LOAD DESEQ2 FOR DEG FILTERING
# =============================================================================

cat("=== PHASE 4: Loading DESeq2 Results ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

degs_up <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange > 0)
cat(sprintf("  DEGs UP: %d\n", nrow(degs_up)))

degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0)
cat(sprintf("  DEGs DOWN: %d\n", nrow(degs_down)))

degs_up_symbols <- unique(degs_up$gene_symbol[!is.na(degs_up$gene_symbol)])
degs_down_symbols <- unique(degs_down$gene_symbol[!is.na(degs_down$gene_symbol)])

cat("\n")

# =============================================================================
# PHASE 5: CREATE GENE SETS (STRATIFIED BY DMR STATUS)
# =============================================================================

cat("=== PHASE 5: Creating Gene Sets (Stratified by DMR) ===\n")

# VERSION 1: All genes with peaks, stratified by DMR
tes_peaked_genes <- intersect(tes_symbols, genes_lookup$gene_name)
tead1_peaked_genes <- intersect(tead1_symbols, genes_lookup$gene_name)

# Stratify by DMR presence
tes_peaked_with_dmr <- intersect(tes_peaked_genes, dmr_gene_names)
tes_peaked_no_dmr <- setdiff(tes_peaked_genes, dmr_gene_names)
tead1_peaked_with_dmr <- intersect(tead1_peaked_genes, dmr_gene_names)
tead1_peaked_no_dmr <- setdiff(tead1_peaked_genes, dmr_gene_names)

cat(sprintf("  VERSION 1 (All peaked genes):\n"))
cat(sprintf("    TES-peaked with DMR: %d\n", length(tes_peaked_with_dmr)))
cat(sprintf("    TES-peaked no DMR: %d\n", length(tes_peaked_no_dmr)))
cat(sprintf("    TEAD1-peaked with DMR: %d\n", length(tead1_peaked_with_dmr)))
cat(sprintf("    TEAD1-peaked no DMR: %d\n", length(tead1_peaked_no_dmr)))

# VERSION 2: DEGs with peaks, stratified by DMR
tes_degs_down_peaked <- intersect(intersect(tes_symbols, degs_down_symbols), genes_lookup$gene_name)
tes_degs_down_dmr <- intersect(tes_degs_down_peaked, dmr_gene_names)
tes_degs_down_no_dmr <- setdiff(tes_degs_down_peaked, dmr_gene_names)

tead1_degs_down_peaked <- intersect(intersect(tead1_symbols, degs_down_symbols), genes_lookup$gene_name)
tead1_degs_down_dmr <- intersect(tead1_degs_down_peaked, dmr_gene_names)
tead1_degs_down_no_dmr <- setdiff(tead1_degs_down_peaked, dmr_gene_names)

cat(sprintf("\n  VERSION 2 (DEGs DOWN with peaks):\n"))
cat(sprintf("    TES DEGs DOWN with DMR: %d\n", length(tes_degs_down_dmr)))
cat(sprintf("    TES DEGs DOWN no DMR: %d\n", length(tes_degs_down_no_dmr)))
cat(sprintf("    TEAD1 DEGs DOWN with DMR: %d\n", length(tead1_degs_down_dmr)))
cat(sprintf("    TEAD1 DEGs DOWN no DMR: %d\n", length(tead1_degs_down_no_dmr)))

cat("\n")

# =============================================================================
# PHASE 6: CREATE BED FILES FOR DEEPTOOLS
# =============================================================================

cat("=== PHASE 6: Creating BED Files ===\n")

create_gene_body_bed <- function(gene_symbols, gene_lookup, output_file) {
    gene_data <- gene_lookup %>%
        filter(gene_name %in% gene_symbols)

    if (nrow(gene_data) == 0) {
        cat(sprintf("  WARNING: No genes found for %s\n", basename(output_file)))
        return(0)
    }

    bed <- data.frame(
        chr = gene_data$seqnames,
        start = gene_data$start - 1,
        end = gene_data$end,
        name = gene_data$gene_name,
        score = 0,
        strand = gene_data$strand,
        stringsAsFactors = FALSE
    )

    bed <- bed[!is.na(bed$start) & !is.na(bed$end), ]
    bed <- bed[bed$start >= 0, ]
    bed <- bed[order(bed$chr, bed$start), ]
    bed <- bed[!duplicated(paste(bed$chr, bed$start, bed$end)), ]

    write.table(bed, output_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    cat(sprintf("  Created: %s (%d genes)\n", basename(output_file), nrow(bed)))
    return(nrow(bed))
}

# VERSION 1: All peaked genes (stratified by DMR)
cat("\n  VERSION 1 - All peaked genes:\n")
n_tes_dmr <- create_gene_body_bed(tes_peaked_with_dmr, genes_lookup,
                                   file.path(OUTPUT_BASE, "TES_peaked_with_DMR.bed"))
n_tes_nodmr <- create_gene_body_bed(tes_peaked_no_dmr, genes_lookup,
                                     file.path(OUTPUT_BASE, "TES_peaked_no_DMR.bed"))
n_tead1_dmr <- create_gene_body_bed(tead1_peaked_with_dmr, genes_lookup,
                                     file.path(OUTPUT_BASE, "TEAD1_peaked_with_DMR.bed"))
n_tead1_nodmr <- create_gene_body_bed(tead1_peaked_no_dmr, genes_lookup,
                                       file.path(OUTPUT_BASE, "TEAD1_peaked_no_DMR.bed"))

# VERSION 2: DEGs DOWN with peaks (stratified by DMR)
cat("\n  VERSION 2 - DEGs DOWN with peaks:\n")
n_tes_down_dmr <- create_gene_body_bed(tes_degs_down_dmr, genes_lookup,
                                        file.path(OUTPUT_BASE, "TES_DEGs_DOWN_with_DMR.bed"))
n_tes_down_nodmr <- create_gene_body_bed(tes_degs_down_no_dmr, genes_lookup,
                                          file.path(OUTPUT_BASE, "TES_DEGs_DOWN_no_DMR.bed"))
n_tead1_down_dmr <- create_gene_body_bed(tead1_degs_down_dmr, genes_lookup,
                                          file.path(OUTPUT_BASE, "TEAD1_DEGs_DOWN_with_DMR.bed"))
n_tead1_down_nodmr <- create_gene_body_bed(tead1_degs_down_no_dmr, genes_lookup,
                                            file.path(OUTPUT_BASE, "TEAD1_DEGs_DOWN_no_DMR.bed"))

# Also create original-style files for comparison
cat("\n  Original-style files:\n")
n_tes_all <- create_gene_body_bed(tes_peaked_genes, genes_lookup,
                                   file.path(OUTPUT_BASE, "TES_peaked_genes_all.bed"))
n_tead1_all <- create_gene_body_bed(tead1_peaked_genes, genes_lookup,
                                     file.path(OUTPUT_BASE, "TEAD1_peaked_genes_all.bed"))
n_tes_down_all <- create_gene_body_bed(tes_degs_down_peaked, genes_lookup,
                                        file.path(OUTPUT_BASE, "TES_DEGs_DOWN_all.bed"))
n_tead1_down_all <- create_gene_body_bed(tead1_degs_down_peaked, genes_lookup,
                                          file.path(OUTPUT_BASE, "TEAD1_DEGs_DOWN_all.bed"))

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("==========================================================\n")
cat("BED FILE PREPARATION COMPLETE (HIGH-CONFIDENCE VERSION)\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))

cat("VERSION 1 - All peaked genes (stratified by DMR):\n")
cat(sprintf("  - TES_peaked_with_DMR.bed (%d genes) - MOST RELIABLE\n", n_tes_dmr))
cat(sprintf("  - TES_peaked_no_DMR.bed (%d genes)\n", n_tes_nodmr))
cat(sprintf("  - TEAD1_peaked_with_DMR.bed (%d genes) - MOST RELIABLE\n", n_tead1_dmr))
cat(sprintf("  - TEAD1_peaked_no_DMR.bed (%d genes)\n", n_tead1_nodmr))

cat("\nVERSION 2 - DEGs DOWN with peaks (stratified by DMR):\n")
cat(sprintf("  - TES_DEGs_DOWN_with_DMR.bed (%d genes) - MOST RELIABLE\n", n_tes_down_dmr))
cat(sprintf("  - TES_DEGs_DOWN_no_DMR.bed (%d genes)\n", n_tes_down_nodmr))
cat(sprintf("  - TEAD1_DEGs_DOWN_with_DMR.bed (%d genes) - MOST RELIABLE\n", n_tead1_down_dmr))
cat(sprintf("  - TEAD1_DEGs_DOWN_no_DMR.bed (%d genes)\n", n_tead1_down_nodmr))

cat("\nKey insight: Genes with high-confidence DMRs have reliable methylation\n")
cat("data because BOTH TES and GFP samples have >2 reads at these regions.\n")

cat("\nNext: Run deepTools to generate metagene profiles\n")
