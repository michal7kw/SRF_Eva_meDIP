#!/usr/bin/env Rscript
#
# DEGs DOWN STRATIFIED BY TES/TEAD1 BINDING - HIGH-CONFIDENCE DMR VERSION
# =============================================================================
#
# Creates gene sets based on:
#   1. DEG status (downregulated in TES vs GFP)
#   2. TES/TEAD1 binding at promoter (±2kb from TSS)
#   3. Presence of HIGH-CONFIDENCE DMRs (where both samples have >2 reads)
#
# Gene Groups:
#   1. DEGs DOWN WITH binding + WITH high-conf DMR
#   2. DEGs DOWN WITH binding - NO high-conf DMR
#   3. DEGs DOWN WITHOUT binding + WITH high-conf DMR
#   4. DEGs DOWN WITHOUT binding - NO high-conf DMR
#   5. Random control (no binding, no DMR)
#
# This version addresses GFP library quality issues by focusing on regions
# where BOTH samples have meaningful signal.
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")

cat("==========================================================\n")
cat("DEGs DOWN STRATIFIED BY BINDING - HIGH-CONFIDENCE VERSION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak files (narrow peaks from merged replicates)
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# HIGH-CONFIDENCE DMRs (filtered for both samples >2 reads)
DMR_FILE <- "output/07_differential_MEDIPS_confident/TES_vs_GFP_DMRs_confident.csv"

# Output directory
OUTPUT_BASE <- "output/29_degs_down_binding_methylation_confident"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Parameters
PROMOTER_UPSTREAM <- 2000   # bp upstream of TSS
PROMOTER_DOWNSTREAM <- 2000 # bp downstream of TSS
GENE_BODY_EXTENSION <- 5000 # bp extension beyond gene body for DMR overlap

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS
# =============================================================================

cat("=== PHASE 1: Loading Gene Annotations ===\n")

gtf <- rtracklayer::import(GTF_FILE)
genes_gtf <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes_gtf)

# Filter for protein-coding genes on standard chromosomes
genes_df <- genes_df %>%
    filter(gene_type == "protein_coding") %>%
    filter(seqnames %in% paste0("chr", c(1:22, "X", "Y")))

cat(sprintf("  Protein-coding genes on standard chromosomes: %d\n", nrow(genes_df)))

# Clean gene IDs (remove version)
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

# Create GRanges for genes
genes_gr <- GRanges(
    seqnames = genes_df$seqnames,
    ranges = IRanges(start = genes_df$start, end = genes_df$end),
    strand = genes_df$strand,
    gene_name = genes_df$gene_name,
    gene_id = genes_df$gene_id_clean
)

# Create promoter regions (±2kb from TSS)
promoters_gr <- promoters(genes_gr, upstream = PROMOTER_UPSTREAM, downstream = PROMOTER_DOWNSTREAM)

# Create extended gene regions (for DMR overlap)
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

cat(sprintf("  Created promoter regions: +/-%d bp from TSS\n", PROMOTER_UPSTREAM))
cat(sprintf("  Created extended gene regions: +/-%d bp beyond gene body\n", GENE_BODY_EXTENSION))
cat("\n")

# =============================================================================
# PHASE 2: LOAD TES/TEAD1 PEAKS
# =============================================================================

cat("=== PHASE 2: Loading TES/TEAD1 Binding Peaks ===\n")

# Function to load narrowPeak files
load_peaks <- function(peak_file, name) {
    peaks <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("chr", "start", "end", "name", "score",
                                      "strand", "signalValue", "pValue", "qValue", "peak"))
    cat(sprintf("  %s peaks loaded: %d\n", name, nrow(peaks)))

    # Add "chr" prefix if not present
    chr_names <- peaks$chr
    if (!grepl("^chr", chr_names[1])) {
        chr_names <- paste0("chr", chr_names)
    }

    GRanges(
        seqnames = chr_names,
        ranges = IRanges(start = peaks$start + 1, end = peaks$end),  # Convert to 1-based
        score = peaks$score,
        signalValue = peaks$signalValue
    )
}

tes_peaks_gr <- load_peaks(TES_PEAKS, "TES")
tead1_peaks_gr <- load_peaks(TEAD1_PEAKS, "TEAD1")

# Combine all TES/TEAD1 peaks
all_binding_peaks <- c(tes_peaks_gr, tead1_peaks_gr)
cat(sprintf("  Combined TES+TEAD1 peaks: %d\n", length(all_binding_peaks)))
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

# Create GRanges
dmr_gr <- GRanges(
    seqnames = dmr_chr,
    ranges = IRanges(start = dmr_data$start, end = dmr_data$stop),
    logFC = dmr_data$logFC,
    FDR = dmr_data$FDR
)

# Separate hyper and hypo
hyper_dmrs <- dmr_gr[dmr_gr$logFC > 0]
hypo_dmrs <- dmr_gr[dmr_gr$logFC < 0]
cat(sprintf("  Hypermethylated (TES > GFP): %d (%.1f%%)\n",
            length(hyper_dmrs), 100*length(hyper_dmrs)/length(dmr_gr)))
cat(sprintf("  Hypomethylated (TES < GFP): %d (%.1f%%)\n",
            length(hypo_dmrs), 100*length(hypo_dmrs)/length(dmr_gr)))
cat("\n")

# =============================================================================
# PHASE 4: IDENTIFY GENES WITH/WITHOUT BINDING AND DMRs
# =============================================================================

cat("=== PHASE 4: Classifying Genes by Binding and DMR Status ===\n")

# Find genes with binding at promoter
promoter_overlaps <- findOverlaps(promoters_gr, all_binding_peaks)
genes_with_binding_idx <- unique(queryHits(promoter_overlaps))
genes_without_binding_idx <- setdiff(1:length(genes_gr), genes_with_binding_idx)

cat(sprintf("  Genes with TES/TEAD1 binding at promoter: %d\n", length(genes_with_binding_idx)))
cat(sprintf("  Genes without TES/TEAD1 binding at promoter: %d\n", length(genes_without_binding_idx)))

# Find genes with high-confidence DMR overlap (extended gene body)
dmr_overlaps <- findOverlaps(extended_genes_gr, dmr_gr)
genes_with_dmr_idx <- unique(queryHits(dmr_overlaps))
genes_without_dmr_idx <- setdiff(1:length(genes_gr), genes_with_dmr_idx)

cat(sprintf("  Genes with high-confidence DMR overlap: %d\n", length(genes_with_dmr_idx)))
cat(sprintf("  Genes without high-confidence DMR overlap: %d\n", length(genes_without_dmr_idx)))

# Get gene IDs for each category
bound_gene_ids <- genes_gr$gene_id[genes_with_binding_idx]
unbound_gene_ids <- genes_gr$gene_id[genes_without_binding_idx]
dmr_gene_ids <- genes_gr$gene_id[genes_with_dmr_idx]
no_dmr_gene_ids <- genes_gr$gene_id[genes_without_dmr_idx]

cat("\n")

# =============================================================================
# PHASE 5: LOAD DESEQ2 AND DEFINE GENE SETS
# =============================================================================

cat("=== PHASE 5: Loading DESeq2 Results and Defining Gene Sets ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

cat(sprintf("  Total genes in DESeq2: %d\n", nrow(deseq2)))

# All expressed genes (non-NA padj)
expressed_genes <- deseq2 %>%
    filter(!is.na(padj)) %>%
    pull(gene_id_clean)
cat(sprintf("  Expressed genes (non-NA padj): %d\n", length(expressed_genes)))

# DEGs DOWN: significantly downregulated (padj < 0.05 & log2FC < 0)
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_id_clean)
cat(sprintf("  DEGs DOWN (padj<0.05, log2FC<0): %d\n", length(degs_down)))

# =============================================================================
# PHASE 6: CREATE 5 GENE SETS
# =============================================================================

cat("\n=== PHASE 6: Creating 5 Gene Sets ===\n")

# GROUP 1: DEGs DOWN WITH binding AND WITH high-confidence DMR
degs_down_bound_dmr <- intersect(degs_down, intersect(bound_gene_ids, dmr_gene_ids))
cat(sprintf("  GROUP 1 - DEGs DOWN + binding + DMR: %d genes\n", length(degs_down_bound_dmr)))

# GROUP 2: DEGs DOWN WITH binding BUT WITHOUT high-confidence DMR
degs_down_bound_no_dmr <- intersect(degs_down, intersect(bound_gene_ids, no_dmr_gene_ids))
cat(sprintf("  GROUP 2 - DEGs DOWN + binding - no DMR: %d genes\n", length(degs_down_bound_no_dmr)))

# GROUP 3: DEGs DOWN WITHOUT binding BUT WITH high-confidence DMR
degs_down_no_bound_dmr <- intersect(degs_down, intersect(unbound_gene_ids, dmr_gene_ids))
cat(sprintf("  GROUP 3 - DEGs DOWN - no binding + DMR: %d genes\n", length(degs_down_no_bound_dmr)))

# GROUP 4: DEGs DOWN WITHOUT binding AND WITHOUT high-confidence DMR
degs_down_no_bound_no_dmr <- intersect(degs_down, intersect(unbound_gene_ids, no_dmr_gene_ids))
cat(sprintf("  GROUP 4 - DEGs DOWN - no binding - no DMR: %d genes\n", length(degs_down_no_bound_no_dmr)))

# GROUP 5: Random expressed genes (same N as group 1) WITHOUT binding and WITHOUT DMR
all_de_genes <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05) %>%
    pull(gene_id_clean)

non_de_expressed_unbound_no_dmr <- setdiff(
    intersect(intersect(expressed_genes, unbound_gene_ids), no_dmr_gene_ids),
    all_de_genes
)
cat(sprintf("  Non-DE expressed genes without binding and without DMR (pool): %d\n",
            length(non_de_expressed_unbound_no_dmr)))

# Sample same number as group 1
set.seed(42)  # For reproducibility
n_to_sample <- length(degs_down_bound_dmr)
if (n_to_sample == 0) {
    n_to_sample <- min(100, length(degs_down))  # fallback
}
if (length(non_de_expressed_unbound_no_dmr) >= n_to_sample) {
    random_control <- sample(non_de_expressed_unbound_no_dmr, n_to_sample)
} else {
    cat("  WARNING: Not enough non-DE unbound genes. Using all available.\n")
    random_control <- non_de_expressed_unbound_no_dmr
}
cat(sprintf("  GROUP 5 - Random control (matched N, no binding, no DMR): %d genes\n", length(random_control)))

cat("\n")

# =============================================================================
# PHASE 7: CREATE BED FILES
# =============================================================================

cat("=== PHASE 7: Creating BED Files ===\n")

# Function to create gene body BED file
create_gene_body_bed <- function(gene_ids_to_include, genes_df, output_file, description) {
    selected <- genes_df[genes_df$gene_id_clean %in% gene_ids_to_include, ]

    if (nrow(selected) == 0) {
        cat(sprintf("  WARNING: No genes matched for %s\n", description))
        return(0)
    }

    bed <- data.frame(
        chr = selected$seqnames,
        start = selected$start - 1,  # BED is 0-based
        end = selected$end,
        name = selected$gene_name,
        score = 0,
        strand = selected$strand,
        stringsAsFactors = FALSE
    )

    # Filter valid entries
    bed <- bed[!is.na(bed$start) & !is.na(bed$end), ]
    bed <- bed[bed$start >= 0, ]

    # Sort by chromosome and position
    bed <- bed[order(bed$chr, bed$start), ]

    # Remove duplicates
    bed <- bed[!duplicated(paste(bed$chr, bed$start, bed$end)), ]

    write.table(bed, output_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    cat(sprintf("  Created: %s (%d genes) - %s\n",
                basename(output_file), nrow(bed), description))
    return(nrow(bed))
}

# Create BED files for each group
n_group1 <- create_gene_body_bed(
    degs_down_bound_dmr, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_binding_DMR.bed"),
    "DEGs DOWN + binding + high-conf DMR"
)

n_group2 <- create_gene_body_bed(
    degs_down_bound_no_dmr, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_binding_noDMR.bed"),
    "DEGs DOWN + binding - no DMR"
)

n_group3 <- create_gene_body_bed(
    degs_down_no_bound_dmr, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_noBinding_DMR.bed"),
    "DEGs DOWN - no binding + DMR"
)

n_group4 <- create_gene_body_bed(
    degs_down_no_bound_no_dmr, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_noBinding_noDMR.bed"),
    "DEGs DOWN - no binding - no DMR"
)

n_group5 <- create_gene_body_bed(
    random_control, genes_df,
    file.path(OUTPUT_BASE, "random_control.bed"),
    "Random control (no binding, no DMR)"
)

# Also create simplified version for original-style analysis (3 groups)
# Group A: DEGs DOWN with binding (any DMR status)
degs_down_with_binding <- intersect(degs_down, bound_gene_ids)
n_groupA <- create_gene_body_bed(
    degs_down_with_binding, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_with_binding.bed"),
    "DEGs DOWN with binding (all)"
)

# Group B: DEGs DOWN without binding (any DMR status)
degs_down_without_binding <- intersect(degs_down, unbound_gene_ids)
n_groupB <- create_gene_body_bed(
    degs_down_without_binding, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_without_binding.bed"),
    "DEGs DOWN without binding (all)"
)

# Save gene lists as CSV for reference
gene_list_df <- data.frame(
    gene_id = c(
        degs_down_bound_dmr, degs_down_bound_no_dmr,
        degs_down_no_bound_dmr, degs_down_no_bound_no_dmr,
        random_control
    ),
    group = c(
        rep("DEGs_DOWN_binding_DMR", length(degs_down_bound_dmr)),
        rep("DEGs_DOWN_binding_noDMR", length(degs_down_bound_no_dmr)),
        rep("DEGs_DOWN_noBinding_DMR", length(degs_down_no_bound_dmr)),
        rep("DEGs_DOWN_noBinding_noDMR", length(degs_down_no_bound_no_dmr)),
        rep("random_control", length(random_control))
    )
)
write.csv(gene_list_df, file.path(OUTPUT_BASE, "gene_lists.csv"), row.names = FALSE)

# Save summary counts
counts_df <- data.frame(
    group = c("DEGs_DOWN_binding_DMR", "DEGs_DOWN_binding_noDMR",
              "DEGs_DOWN_noBinding_DMR", "DEGs_DOWN_noBinding_noDMR",
              "random_control"),
    n_genes = c(n_group1, n_group2, n_group3, n_group4, n_group5),
    description = c(
        "DEGs DOWN with binding AND high-confidence DMR",
        "DEGs DOWN with binding but NO high-confidence DMR",
        "DEGs DOWN without binding but WITH high-confidence DMR",
        "DEGs DOWN without binding AND NO high-confidence DMR",
        "Random expressed non-DE genes (matched N)"
    )
)
write.csv(counts_df, file.path(OUTPUT_BASE, "gene_counts.csv"), row.names = FALSE)

cat("\n")

# =============================================================================
# PHASE 8: ADDITIONAL STATISTICS
# =============================================================================

cat("=== PHASE 8: Additional Statistics ===\n")

# Get log2FC values for DEGs DOWN stratified by binding+DMR
groups_of_interest <- c(degs_down_bound_dmr, degs_down_bound_no_dmr,
                        degs_down_no_bound_dmr, degs_down_no_bound_no_dmr)

deseq2_subset <- deseq2[deseq2$gene_id_clean %in% groups_of_interest, ]
deseq2_subset$group <- case_when(
    deseq2_subset$gene_id_clean %in% degs_down_bound_dmr ~ "binding_DMR",
    deseq2_subset$gene_id_clean %in% degs_down_bound_no_dmr ~ "binding_noDMR",
    deseq2_subset$gene_id_clean %in% degs_down_no_bound_dmr ~ "noBinding_DMR",
    deseq2_subset$gene_id_clean %in% degs_down_no_bound_no_dmr ~ "noBinding_noDMR"
)

# Summary statistics
group_stats <- deseq2_subset %>%
    group_by(group) %>%
    summarise(
        n = n(),
        mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
        median_log2FC = median(log2FoldChange, na.rm = TRUE),
        sd_log2FC = sd(log2FoldChange, na.rm = TRUE),
        .groups = "drop"
    )

cat("\n  Log2FC statistics for DEGs DOWN by Binding+DMR status:\n")
print(as.data.frame(group_stats), row.names = FALSE)

# Save statistics
write.csv(group_stats, file.path(OUTPUT_BASE, "log2FC_statistics.csv"), row.names = FALSE)

cat("\n")

# =============================================================================
# PHASE 9: HIGH-CONFIDENCE DMR SUMMARY AT GENES
# =============================================================================

cat("=== PHASE 9: DMR Direction Summary at DEGs DOWN ===\n")

# Function to count DMR directions for a gene set
count_dmr_directions <- function(gene_ids, genes_gr, dmr_gr, ext_genes_gr, label) {
    idx <- which(genes_gr$gene_id %in% gene_ids)
    if (length(idx) == 0) return(NULL)

    ext_subset <- ext_genes_gr[idx]
    overlaps <- findOverlaps(ext_subset, dmr_gr)

    if (length(overlaps) == 0) {
        return(data.frame(
            group = label,
            n_genes = length(gene_ids),
            n_genes_with_dmr = 0,
            n_hyper_dmrs = 0,
            n_hypo_dmrs = 0,
            pct_hyper = NA,
            pct_hypo = NA
        ))
    }

    dmr_subset <- dmr_gr[subjectHits(overlaps)]
    n_hyper <- sum(dmr_subset$logFC > 0)
    n_hypo <- sum(dmr_subset$logFC < 0)

    data.frame(
        group = label,
        n_genes = length(gene_ids),
        n_genes_with_dmr = length(unique(queryHits(overlaps))),
        n_hyper_dmrs = n_hyper,
        n_hypo_dmrs = n_hypo,
        pct_hyper = 100 * n_hyper / (n_hyper + n_hypo),
        pct_hypo = 100 * n_hypo / (n_hyper + n_hypo)
    )
}

dmr_direction_summary <- rbind(
    count_dmr_directions(degs_down_with_binding, genes_gr, dmr_gr, extended_genes_gr,
                         "DEGs_DOWN_with_binding"),
    count_dmr_directions(degs_down_without_binding, genes_gr, dmr_gr, extended_genes_gr,
                         "DEGs_DOWN_without_binding")
)

cat("\n  High-confidence DMR directions at DEGs DOWN:\n")
print(dmr_direction_summary, row.names = FALSE)

write.csv(dmr_direction_summary, file.path(OUTPUT_BASE, "dmr_direction_at_degs.csv"),
          row.names = FALSE)

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("==========================================================\n")
cat("GENE SET CREATION COMPLETE (HIGH-CONFIDENCE VERSION)\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))

cat("Created files (5-group stratification):\n")
cat(sprintf("  1. DEGs_DOWN_binding_DMR.bed (%d genes)\n", n_group1))
cat(sprintf("     -> DEGs DOWN with TES/TEAD1 binding AND high-conf DMR\n"))
cat(sprintf("     -> MOST RELIABLE for methylation analysis\n\n"))

cat(sprintf("  2. DEGs_DOWN_binding_noDMR.bed (%d genes)\n", n_group2))
cat(sprintf("     -> DEGs DOWN with binding but NO high-conf DMR detected\n"))
cat(sprintf("     -> May have unreliable GFP signal in original analysis\n\n"))

cat(sprintf("  3. DEGs_DOWN_noBinding_DMR.bed (%d genes)\n", n_group3))
cat(sprintf("     -> DEGs DOWN without binding but WITH high-conf DMR\n"))
cat(sprintf("     -> Indirect methylation targets?\n\n"))

cat(sprintf("  4. DEGs_DOWN_noBinding_noDMR.bed (%d genes)\n", n_group4))
cat(sprintf("     -> DEGs DOWN without binding AND without DMR\n"))
cat(sprintf("     -> Purely expression-based regulation\n\n"))

cat(sprintf("  5. random_control.bed (%d genes)\n", n_group5))
cat(sprintf("     -> Random control (no binding, no DMR)\n\n"))

cat("Also created (for original-style 3-group analysis):\n")
cat(sprintf("  - DEGs_DOWN_with_binding.bed (%d genes)\n", n_groupA))
cat(sprintf("  - DEGs_DOWN_without_binding.bed (%d genes)\n", n_groupB))

cat("\nKey insight from high-confidence filtering:\n")
cat("  Original DMRs: 91% hypermethylated (artifact of GFP dropout)\n")
cat("  High-conf DMRs: 35% hyper / 65% hypo (real biology!)\n")

cat("\nNext step: Run the shell script to compute metagene profiles\n")
cat("  sbatch 29_degs_down_binding_methylation_confident.sh\n")
