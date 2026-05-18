#!/usr/bin/env Rscript

################################################################################
# Script: 12_define_gene_sets.R
# Purpose: Define gene sets based on RNA-seq expression data for meDIP analysis
#
# Description:
#   Creates gene lists stratified by expression level and differential expression
#   for use in meDIP heatmap and metaprofile analysis. Generates BED files with
#   promoter coordinates (TSS ±2kb) for each gene set.
#
# Gene Sets Created:
#   1. All genes - All genes with valid expression data
#   2. Highly expressed - Top 25% by mean normalized counts
#   3. Lowly expressed - Bottom 25% by mean normalized counts
#   4. Upregulated - DEGs with log2FC > 1, FDR < 0.05 (TES vs GFP)
#   5. Downregulated - DEGs with log2FC < -1, FDR < 0.05 (TES vs GFP)
#   6. Unchanged - Not significantly differentially expressed
#
# Input:
#   - RNA-seq normalized counts
#   - DESeq2 differential expression results
#   - GTF annotation for promoter coordinates
#
# Output:
#   - BED files for each gene set (promoter TSS points)
#   - Gene lists (CSV format)
#   - Summary statistics
#
# Runtime: ~5-10 minutes
################################################################################

cat("====================================================\n")
cat("Gene Set Definition for meDIP Heatmap Analysis\n")
cat("====================================================\n")
cat(paste("Start:", Sys.time(), "\n\n"))

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
})

# Define paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
rnaseq_dir <- file.path(base_dir, "SRF_Eva_RNA/results/05_deseq2")
annotation_file <- file.path(base_dir, "../COMMONS/annotation/gencode.v44.annotation.gtf")
out_dir <- file.path(base_dir, "meDIP/results/12_gene_sets")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

################################################################################
# 1. Load RNA-seq Data
################################################################################

cat("========================================\n")
cat("1. Loading RNA-seq Data\n")
cat("========================================\n")

# Load normalized counts
counts_file <- file.path(rnaseq_dir, "normalized_counts.txt")
if (!file.exists(counts_file)) {
    stop("ERROR: Normalized counts file not found at: ", counts_file)
}
counts <- read.delim(counts_file, stringsAsFactors = FALSE)

cat(paste("  Loaded normalized counts for", nrow(counts), "genes\n"))
cat(paste("  Samples:", paste(colnames(counts)[-1], collapse = ", "), "\n\n"))

# Load differential expression results
deseq_file <- file.path(rnaseq_dir, "deseq2_results_TES_vs_GFP.txt")
if (!file.exists(deseq_file)) {
    stop("ERROR: DESeq2 results file not found at: ", deseq_file)
}
degs <- read.delim(deseq_file, stringsAsFactors = FALSE)

cat(paste("  Loaded DESeq2 results for", nrow(degs), "genes\n"))
cat(paste(
    "  Significant DEGs (FDR<0.05, |FC|>2):",
    sum(!is.na(degs$padj) & degs$padj < 0.05 & abs(degs$log2FoldChange) > 1), "\n\n"
))

################################################################################
# 2. Calculate Mean Expression Across All Samples
################################################################################

cat("========================================\n")
cat("2. Calculating Expression Statistics\n")
cat("========================================\n")

# Calculate mean expression across all samples
count_cols <- grep("^ASE01_", colnames(counts), value = TRUE)
counts$mean_expression <- rowMeans(counts[, count_cols], na.rm = TRUE)

# Calculate mean for GFP and TES separately
gfp_cols <- grep("GFP", count_cols, value = TRUE)
tes_cols <- grep("TES", count_cols, value = TRUE)
counts$mean_GFP <- rowMeans(counts[, gfp_cols], na.rm = TRUE)
counts$mean_TES <- rowMeans(counts[, tes_cols], na.rm = TRUE)

cat(paste(
    "  Mean expression range:",
    round(min(counts$mean_expression), 2), "-",
    round(max(counts$mean_expression), 2), "\n"
))

# Remove genes with zero or near-zero expression
counts_filtered <- counts[counts$mean_expression > 1, ]
cat(paste("  Genes with expression > 1:", nrow(counts_filtered), "\n\n"))

################################################################################
# 3. Define Gene Sets by Expression Level
################################################################################

cat("========================================\n")
cat("3. Defining Gene Sets by Expression\n")
cat("========================================\n")

# Calculate quartiles
expr_quantiles <- quantile(counts_filtered$mean_expression, probs = c(0.25, 0.75))
cat(paste("  25th percentile:", round(expr_quantiles[1], 2), "\n"))
cat(paste("  75th percentile:", round(expr_quantiles[2], 2), "\n\n"))

# Define expression-based gene sets
highly_expressed <- counts_filtered[counts_filtered$mean_expression >= expr_quantiles[2], ]
lowly_expressed <- counts_filtered[counts_filtered$mean_expression <= expr_quantiles[1], ]

cat(paste("  Highly expressed genes (top 25%):", nrow(highly_expressed), "\n"))
cat(paste("  Lowly expressed genes (bottom 25%):", nrow(lowly_expressed), "\n\n"))

################################################################################
# 4. Define Gene Sets by Differential Expression
################################################################################

cat("========================================\n")
cat("4. Defining Gene Sets by Regulation\n")
cat("========================================\n")

# Merge DESeq2 results with expression data
gene_data <- merge(counts_filtered, degs, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)

# Define DEG thresholds
fdr_threshold <- 0.05
fc_threshold <- 1 # log2 scale (2-fold)

# Identify DEGs
upregulated <- gene_data[
    !is.na(gene_data$padj) &
        gene_data$padj < fdr_threshold &
        gene_data$log2FoldChange > fc_threshold,
]

downregulated <- gene_data[
    !is.na(gene_data$padj) &
        gene_data$padj < fdr_threshold &
        gene_data$log2FoldChange < -fc_threshold,
]

unchanged <- gene_data[
    is.na(gene_data$padj) |
        gene_data$padj >= fdr_threshold |
        abs(gene_data$log2FoldChange) <= fc_threshold,
]

cat(paste("  Upregulated genes (log2FC > 1, FDR < 0.05):", nrow(upregulated), "\n"))
cat(paste("  Downregulated genes (log2FC < -1, FDR < 0.05):", nrow(downregulated), "\n"))
cat(paste("  Unchanged genes:", nrow(unchanged), "\n\n"))

################################################################################
# 5. Load GTF Annotation and Extract Promoter Coordinates
################################################################################

cat("========================================\n")
cat("5. Extracting Promoter Coordinates\n")
cat("========================================\n")

if (!file.exists(annotation_file)) {
    stop("ERROR: GTF annotation file not found at: ", annotation_file)
}

cat("  Loading GTF annotation...\n")
txdb <- makeTxDbFromGFF(annotation_file, format = "gtf")

# Extract gene coordinates
genes_gr <- genes(txdb)
cat(paste("  Loaded", length(genes_gr), "gene annotations\n"))

# Define promoter regions (TSS points)
# We define promoters as just the 1bp TSS so that computeMatrix centers exactly on it
promoters_gr <- promoters(genes_gr, upstream = 0, downstream = 1)

cat(paste("  Defined promoters as TSS points (1bp)\n\n"))

################################################################################
# 6. Create BED Files for Each Gene Set
################################################################################

cat("========================================\n")
cat("6. Generating BED Files\n")
cat("========================================\n")

# Function to create BED file for a gene set
create_bed_file <- function(gene_set, set_name, promoters) {
    # Use gene IDs with version numbers intact (matches TxDb format)
    gene_ids <- gene_set$gene_id

    # Match with promoter coordinates
    promoters_subset <- promoters[names(promoters) %in% gene_ids]

    if (length(promoters_subset) == 0) {
        cat(paste("  WARNING: No promoters found for", set_name, "\n"))
        return(NULL)
    }

    # Convert to data frame for BED format (0-based start, 1-based end)
    bed_df <- data.frame(
        chr = as.character(seqnames(promoters_subset)),
        start = start(promoters_subset) - 1,
        end = end(promoters_subset),
        name = names(promoters_subset),
        score = 0,
        strand = as.character(strand(promoters_subset))
    )

    # Sort by chromosome and position
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]

    # Write BED file
    bed_file <- paste0(set_name, "_promoters.bed")
    write.table(bed_df, bed_file,
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    cat(paste("  Saved:", nrow(bed_df), "promoters ->", bed_file, "\n"))

    # Also save gene list with expression data
    gene_list_file <- paste0(set_name, "_genes.csv")
    write.csv(gene_set, gene_list_file, row.names = FALSE)

    return(bed_df)
}

# Create BED files for all gene sets
all_genes_bed <- create_bed_file(counts_filtered, "all_genes", promoters_gr)
highly_expressed_bed <- create_bed_file(highly_expressed, "highly_expressed", promoters_gr)
lowly_expressed_bed <- create_bed_file(lowly_expressed, "lowly_expressed", promoters_gr)
upregulated_bed <- create_bed_file(upregulated, "upregulated", promoters_gr)
downregulated_bed <- create_bed_file(downregulated, "downregulated", promoters_gr)
unchanged_bed <- create_bed_file(unchanged, "unchanged", promoters_gr)

cat("\n")

################################################################################
# 7. Summary Visualization
################################################################################

cat("========================================\n")
cat("7. Generating Summary Plots\n")
cat("========================================\n")

# Expression distribution plot
pdf("expression_distribution.pdf", width = 10, height = 6)

p1 <- ggplot(counts_filtered, aes(x = log10(mean_expression + 1))) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(
        xintercept = log10(expr_quantiles[1] + 1),
        linetype = "dashed", color = "red", size = 1
    ) +
    geom_vline(
        xintercept = log10(expr_quantiles[2] + 1),
        linetype = "dashed", color = "red", size = 1
    ) +
    annotate("text",
        x = log10(expr_quantiles[1] + 1), y = Inf,
        label = "25th percentile", vjust = 2, hjust = 1.1, color = "red"
    ) +
    annotate("text",
        x = log10(expr_quantiles[2] + 1), y = Inf,
        label = "75th percentile", vjust = 2, hjust = -0.1, color = "red"
    ) +
    labs(
        title = "Expression Distribution Across All Genes",
        x = "Log10(Mean Normalized Counts + 1)",
        y = "Number of Genes"
    ) +
    theme_bw(base_size = 12)

print(p1)

# MA plot with gene categories
if (nrow(gene_data[!is.na(gene_data$log2FoldChange), ]) > 0) {
    gene_data$category <- "Unchanged"
    gene_data$category[!is.na(gene_data$padj) &
        gene_data$padj < fdr_threshold &
        gene_data$log2FoldChange > fc_threshold] <- "Upregulated"
    gene_data$category[!is.na(gene_data$padj) &
        gene_data$padj < fdr_threshold &
        gene_data$log2FoldChange < -fc_threshold] <- "Downregulated"

    p2 <- ggplot(
        gene_data[!is.na(gene_data$log2FoldChange), ],
        aes(x = log10(baseMean + 1), y = log2FoldChange, color = category)
    ) +
        geom_point(alpha = 0.5, size = 0.8) +
        scale_color_manual(values = c(
            "Upregulated" = "red",
            "Downregulated" = "blue",
            "Unchanged" = "grey60"
        )) +
        geom_hline(
            yintercept = c(-fc_threshold, fc_threshold),
            linetype = "dashed", color = "black"
        ) +
        labs(
            title = "MA Plot: Gene Categories",
            x = "Log10(Mean Expression + 1)",
            y = "Log2 Fold Change (TES vs GFP)",
            color = "Category"
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "top")

    print(p2)
}

dev.off()
cat("  Saved: expression_distribution.pdf\n\n")

################################################################################
# 8. Save Summary Statistics
################################################################################

cat("========================================\n")
cat("8. Saving Summary Statistics\n")
cat("========================================\n")

summary_stats <- data.frame(
    Category = c(
        "All genes", "Highly expressed (top 25%)",
        "Lowly expressed (bottom 25%)", "Upregulated DEGs",
        "Downregulated DEGs", "Unchanged genes"
    ),
    Number_of_genes = c(
        nrow(counts_filtered), nrow(highly_expressed),
        nrow(lowly_expressed), nrow(upregulated),
        nrow(downregulated), nrow(unchanged)
    ),
    Mean_expression = c(
        mean(counts_filtered$mean_expression),
        mean(highly_expressed$mean_expression),
        mean(lowly_expressed$mean_expression),
        mean(upregulated$mean_expression, na.rm = TRUE),
        mean(downregulated$mean_expression, na.rm = TRUE),
        mean(unchanged$mean_expression, na.rm = TRUE)
    )
)

write.csv(summary_stats, "gene_set_summary.csv", row.names = FALSE)
cat("  Saved: gene_set_summary.csv\n\n")

# Print summary table
cat("\nGene Set Summary:\n")
cat("==================\n")
print(summary_stats, row.names = FALSE)
cat("\n")

################################################################################
# Final Output
################################################################################

cat("====================================================\n")
cat("Gene Set Definition Complete!\n")
cat("====================================================\n")
cat(paste("End:", Sys.time(), "\n\n"))

cat("Output files:\n")
cat("  BED files (promoter coordinates TSS point):\n")
cat("    - all_genes_promoters.bed\n")
cat("    - highly_expressed_promoters.bed\n")
cat("    - lowly_expressed_promoters.bed\n")
cat("    - upregulated_promoters.bed\n")
cat("    - downregulated_promoters.bed\n")
cat("    - unchanged_promoters.bed\n\n")
cat("  Gene lists (CSV with expression data):\n")
cat("    - all_genes_genes.csv\n")
cat("    - highly_expressed_genes.csv\n")
cat("    - lowly_expressed_genes.csv\n")
cat("    - upregulated_genes.csv\n")
cat("    - downregulated_genes.csv\n")
cat("    - unchanged_genes.csv\n\n")
cat("  Visualizations:\n")
cat("    - expression_distribution.pdf\n")
cat("    - gene_set_summary.csv\n\n")

cat("Next steps:\n")
cat("  1. Run 13_compute_matrices.sh to calculate meDIP signal at promoters\n")
cat("  2. Run 14_plot_heatmaps.sh to generate heatmaps\n")
cat("  3. Run 15_plot_metaprofiles.sh to generate metaprofiles\n\n")
