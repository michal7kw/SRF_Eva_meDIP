#!/usr/bin/env Rscript

# ============================================================================
# 23_prepare_binding_stratified_beds.R
# ============================================================================
# Purpose: Create gene sets for binding-stratified MeDIP metaprofile analysis
#
# Gene Sets Created:
#   1. Direct targets: DEGs down + TES/TEAD binding (tes_bound OR tead1_bound)
#   2. Indirect targets: DEGs down + NO binding (Neither_bound)
#   3. Random control: Unchanged + NO binding (sampled to match direct size)
#
# Output: BED files with promoter coordinates for deepTools visualization
# ============================================================================

cat("============================================\n")
cat("Binding-Stratified Gene Set Preparation\n")
cat("Started:", format(Sys.time()), "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicFeatures)
    library(rtracklayer)
})

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"

# Input: Gene classification from integrative analysis
CLASSIFICATION_FILE <- file.path(
    BASE_DIR,
    "SRF_Eva_integrated_analysis/results/20_methylation_binding_expression/gene_classification_summary.csv"
)

# Output directory
OUT_DIR <- file.path(BASE_DIR, "meDIP/results/23_binding_stratified_metaprofiles")
BED_DIR <- file.path(OUT_DIR, "beds")

# GTF for gene coordinates
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Random seed for reproducibility
set.seed(42)

# Create output directories
dir.create(BED_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load Data
# ============================================================================

cat("=== Loading gene classification data ===\n\n")

if (!file.exists(CLASSIFICATION_FILE)) {
    stop("Classification file not found: ", CLASSIFICATION_FILE)
}

classification <- read.csv(CLASSIFICATION_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Loaded classification for %d genes\n", nrow(classification)))

# Summary of current classification
cat("\nExpression status summary:\n")
print(table(classification$expression_status))

cat("\nBinding status summary:\n")
print(table(classification$binding_status))

# ============================================================================
# Define Gene Sets
# ============================================================================

cat("\n=== Defining gene sets ===\n\n")

# 1. DIRECT TARGETS: Downregulated AND (TES-bound OR TEAD1-bound)
# This includes: TES_only_bound, TEAD1_only_bound, TES_TEAD1_bound
direct_targets <- classification %>%
    filter(expression_status == "downregulated") %>%
    filter(tes_bound | tead1_bound)

cat(sprintf("Direct targets (DEGs down + TES/TEAD binding): %d genes\n", nrow(direct_targets)))
cat(sprintf("  - TES-only bound: %d\n", sum(direct_targets$binding_status == "TES_only_bound")))
cat(sprintf("  - TEAD1-only bound: %d\n", sum(direct_targets$binding_status == "TEAD1_only_bound")))
cat(sprintf("  - Both bound: %d\n", sum(direct_targets$binding_status == "TES_TEAD1_bound")))

# 2. INDIRECT TARGETS: Downregulated AND Neither_bound
indirect_targets <- classification %>%
    filter(expression_status == "downregulated") %>%
    filter(binding_status == "Neither_bound")

cat(sprintf("\nIndirect targets (DEGs down + NO binding): %d genes\n", nrow(indirect_targets)))

# 3. RANDOM CONTROL: Unchanged AND Neither_bound, sampled to match direct size
random_pool <- classification %>%
    filter(expression_status == "unchanged") %>%
    filter(binding_status == "Neither_bound")

cat(sprintf("\nRandom pool (Unchanged + NO binding): %d genes available\n", nrow(random_pool)))

# Sample to match direct targets size
target_size <- nrow(direct_targets)
if (nrow(random_pool) < target_size) {
    warning("Random pool smaller than direct targets, using all available genes")
    random_control <- random_pool
} else {
    random_control <- random_pool %>%
        sample_n(size = target_size, replace = FALSE)
}

cat(sprintf("Random control (sampled): %d genes\n", nrow(random_control)))

# ============================================================================
# Load Gene Coordinates from GTF
# ============================================================================

cat("\n=== Loading gene coordinates ===\n\n")

if (!file.exists(GTF_FILE)) {
    stop("GTF file not found: ", GTF_FILE)
}

cat("  Loading GTF annotation...\n")
txdb <- makeTxDbFromGFF(GTF_FILE, format = "gtf")
genes_gr <- genes(txdb)

cat(sprintf("  Loaded %d gene annotations\n", length(genes_gr)))

# Clean Ensembl IDs (remove version numbers)
genes_gr$clean_id <- gsub("\\..*", "", names(genes_gr))

# ============================================================================
# Create BED Files
# ============================================================================

cat("\n=== Creating BED files ===\n\n")

# Function to create BED file from gene set
create_bed_file <- function(gene_data, output_name, genes_ref = genes_gr,
                             promoter_size = 2000) {

    # Get Ensembl IDs (clean, no version)
    gene_ids <- unique(gene_data$ensembl_id)

    # Match to reference
    matched_idx <- which(genes_ref$clean_id %in% gene_ids)

    if (length(matched_idx) == 0) {
        warning(sprintf("No genes matched for %s", output_name))
        return(NULL)
    }

    matched_genes <- genes_ref[matched_idx]

    # Get coordinates
    chrs <- as.character(seqnames(matched_genes))
    starts <- start(matched_genes)
    ends <- end(matched_genes)
    strands <- as.character(strand(matched_genes))

    # Use gene symbols where available, otherwise use Ensembl ID
    names_col <- ifelse(
        matched_genes$clean_id %in% gene_data$ensembl_id,
        gene_data$gene_symbol[match(matched_genes$clean_id, gene_data$ensembl_id)],
        matched_genes$clean_id
    )
    # Fall back to Ensembl ID if symbol is NA
    names_col <- ifelse(is.na(names_col), matched_genes$clean_id, names_col)

    # Convert UCSC chr names (chr1) to Ensembl (1) for compatibility with BigWig
    chrs <- gsub("^chr", "", chrs)

    # Calculate TSS based on strand
    tss <- ifelse(strands == "+", starts, ends)

    # Create promoter BED (TSS +/- promoter_size)
    promoter_bed <- data.frame(
        chr = chrs,
        start = pmax(0, tss - promoter_size),
        end = tss + promoter_size,
        name = names_col,
        score = 0,
        strand = strands,
        stringsAsFactors = FALSE
    )

    # Sort by chromosome and position
    promoter_bed <- promoter_bed %>%
        arrange(chr, start)

    # Write BED file
    bed_file <- file.path(BED_DIR, sprintf("%s_promoter.bed", output_name))
    write.table(
        promoter_bed,
        bed_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    cat(sprintf("  %s: %d genes -> %s\n", output_name, nrow(promoter_bed), bed_file))

    return(list(file = bed_file, n_genes = nrow(promoter_bed)))
}

# Create BED files for each gene set
cat("Creating promoter BED files:\n")

direct_result <- create_bed_file(
    direct_targets,
    "direct_bound_downregulated"
)

indirect_result <- create_bed_file(
    indirect_targets,
    "indirect_unbound_downregulated"
)

random_result <- create_bed_file(
    random_control,
    "random_unbound_unchanged"
)

# ============================================================================
# Save Gene Lists (for reference)
# ============================================================================

cat("\n=== Saving gene lists ===\n\n")

# Save CSV files with full gene info
write.csv(
    direct_targets,
    file.path(OUT_DIR, "direct_bound_downregulated_genes.csv"),
    row.names = FALSE
)

write.csv(
    indirect_targets,
    file.path(OUT_DIR, "indirect_unbound_downregulated_genes.csv"),
    row.names = FALSE
)

write.csv(
    random_control,
    file.path(OUT_DIR, "random_unbound_unchanged_genes.csv"),
    row.names = FALSE
)

cat("  Saved gene lists to CSV\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n============================================\n")
cat("Gene Set Preparation Complete\n")
cat("============================================\n\n")

summary_df <- data.frame(
    Category = c(
        "Direct (DEGs down + bound)",
        "Indirect (DEGs down + unbound)",
        "Random control (unchanged + unbound)"
    ),
    N_genes = c(
        nrow(direct_targets),
        nrow(indirect_targets),
        nrow(random_control)
    ),
    Description = c(
        "TES and/or TEAD1 binding at promoter",
        "No TF binding at promoter",
        "Expression unchanged, no binding (control)"
    )
)

print(summary_df, row.names = FALSE)

cat("\n\nOutput files:\n")
cat(sprintf("  BED files: %s/\n", BED_DIR))
cat(sprintf("  Gene lists: %s/\n", OUT_DIR))

cat("\nNext steps:\n")
cat("  1. Run 23_compute_binding_matrices.sh to compute deepTools matrices\n")
cat("  2. Run 23_plot_binding_metaprofiles.sh to generate metaprofiles\n")

cat("\nFinished:", format(Sys.time()), "\n")
