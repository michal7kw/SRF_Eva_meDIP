#!/usr/bin/env Rscript

# Quick diagnostic to check gene ID format mismatch

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(rtracklayer)
})

base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
annotation_file <- file.path(base_dir, "../COMMONS/annotation/gencode.v44.annotation.gtf")
rnaseq_dir <- file.path(base_dir, "SRF_Eva_RNA/results/05_deseq2")

cat("Checking gene ID formats...\n\n")

# 1. Check RNA-seq gene IDs
cat("=== RNA-seq gene IDs (first 10) ===\n")
counts <- read.delim(file.path(rnaseq_dir, "normalized_counts.txt"), stringsAsFactors = FALSE)
cat(head(counts$gene_id, 10), sep = "\n")
cat("\n")

# 2. Check what TxDb produces
cat("=== TxDb gene IDs (first 10) ===\n")
txdb <- makeTxDbFromGFF(annotation_file, format = "gtf")
genes_gr <- genes(txdb)
cat(head(names(genes_gr), 10), sep = "\n")
cat("\n")

# 3. Check if stripping version numbers helps
cat("=== RNA-seq IDs without version (first 10) ===\n")
rna_ids_no_version <- gsub("\\..*", "", counts$gene_id)
cat(head(rna_ids_no_version, 10), sep = "\n")
cat("\n")

# 4. Test matching
cat("=== Matching Test ===\n")
cat(paste("Total RNA-seq genes:", length(counts$gene_id), "\n"))
cat(paste("Total TxDb genes:", length(genes_gr), "\n"))

# Direct match
match_direct <- sum(counts$gene_id %in% names(genes_gr))
cat(paste("Direct matches (with version):", match_direct, "\n"))

# Match without version
match_no_version <- sum(rna_ids_no_version %in% names(genes_gr))
cat(paste("Matches (RNA without version):", match_no_version, "\n"))

# Try stripping TxDb version
txdb_ids_no_version <- gsub("\\..*", "", names(genes_gr))
match_both_no_version <- sum(rna_ids_no_version %in% txdb_ids_no_version)
cat(paste("Matches (both without version):", match_both_no_version, "\n"))

cat("\n=== Sample Matches ===\n")
matching_indices <- which(rna_ids_no_version %in% txdb_ids_no_version)[1:5]
for (i in matching_indices) {
    rna_id <- counts$gene_id[i]
    rna_id_clean <- rna_ids_no_version[i]
    txdb_match <- names(genes_gr)[match(rna_id_clean, txdb_ids_no_version)]
    cat(paste("RNA:", rna_id, "->", rna_id_clean, "-> TxDb:", txdb_match, "\n"))
}
