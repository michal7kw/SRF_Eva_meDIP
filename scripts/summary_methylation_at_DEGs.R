#!/usr/bin/env Rscript
################################################################################
# Summary Statistics: DNA Methylation at Differentially Expressed Genes
#
# This script extracts methylation statistics specifically for DEGs, including:
# 1. Promoter-specific methylation at DEGs
# 2. Gene body methylation at DEGs
# 3. Methylation-expression correlations (properly stratified)
# 4. Direction analysis (hyper/hypo vs up/down)
#
# Author: Generated for SRF_Eva project
# Date: 2024-12
#
# Usage: Rscript summary_methylation_at_DEGs.R
# Output: results/summary/methylation_at_DEGs_summary.csv
#         results/summary/methylation_at_DEGs_top_genes.csv
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

cat("========================================\n")
cat("Methylation at DEGs Summary Statistics\n")
cat("========================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
MEDIP_DIR <- file.path(BASE_DIR, "meDIP")
RNA_DIR <- file.path(BASE_DIR, "SRF_Eva_RNA")
OUTPUT_DIR <- file.path(MEDIP_DIR, "results/summary")

# Thresholds
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Load Data
################################################################################

cat("Loading data files...\n")

# 1. Load DESeq2 results
deseq_file <- file.path(RNA_DIR, "results/05_deseq2/deseq2_results_TES_vs_GFP.txt")
if (!file.exists(deseq_file)) {
    stop("DESeq2 results not found: ", deseq_file)
}
deseq <- read.delim(deseq_file, stringsAsFactors = FALSE)
deseq$gene_id_clean <- gsub("\\..*", "", deseq$gene_id)
cat(sprintf("  Loaded %d genes from DESeq2\n", nrow(deseq)))

# 2. Load annotated DMRs
dmr_file <- file.path(MEDIP_DIR, "results/08_annotation/TES_vs_GFP_annotated.csv")
if (!file.exists(dmr_file)) {
    stop("Annotated DMRs not found: ", dmr_file)
}
dmr_annot <- read.csv(dmr_file, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d annotated DMRs\n", nrow(dmr_annot)))

# 3. Load integrated gene-level data (TSS ±5kb)
integrated_file <- file.path(MEDIP_DIR, "results/16_advanced_visualization/integrated_medip_rnaseq_data.csv")
if (!file.exists(integrated_file)) {
    stop("Integrated data not found: ", integrated_file)
}
integrated <- read.csv(integrated_file, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d genes from integrated data\n", nrow(integrated)))

################################################################################
# Define Gene Categories
################################################################################

cat("\nDefining gene categories...\n")

# Define DEGs in DESeq2 data
deseq <- deseq %>%
    mutate(
        is_DEG = !is.na(padj) & padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD,
        expression_status = case_when(
            is_DEG & log2FoldChange > 0 ~ "upregulated",
            is_DEG & log2FoldChange < 0 ~ "downregulated",
            TRUE ~ "unchanged"
        )
    )

# Define DEGs in integrated data
integrated <- integrated %>%
    mutate(
        is_DEG = !is.na(padj) & padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD,
        expression_status = case_when(
            is_DEG & log2FoldChange > 0 ~ "upregulated",
            is_DEG & log2FoldChange < 0 ~ "downregulated",
            TRUE ~ "unchanged"
        )
    )

n_up <- sum(deseq$expression_status == "upregulated", na.rm = TRUE)
n_down <- sum(deseq$expression_status == "downregulated", na.rm = TRUE)
n_unchanged <- sum(deseq$expression_status == "unchanged", na.rm = TRUE)

cat(sprintf("  Upregulated DEGs: %d\n", n_up))
cat(sprintf("  Downregulated DEGs: %d\n", n_down))
cat(sprintf("  Unchanged genes: %d\n", n_unchanged))

################################################################################
# Initialize Statistics Collection
################################################################################

stats_list <- list()
add_stat <- function(category, metric, value) {
    stats_list[[length(stats_list) + 1]] <<- data.frame(
        Category = category,
        Metric = metric,
        Value = as.character(value),
        stringsAsFactors = FALSE
    )
}

################################################################################
# Section 1: Overview Statistics
################################################################################

cat("\n--- Section 1: Overview ---\n")

add_stat("Overview", "Total_genes_analyzed", nrow(deseq))
add_stat("Overview", "Total_DEGs", n_up + n_down)
add_stat("Overview", "Upregulated_DEGs", n_up)
add_stat("Overview", "Downregulated_DEGs", n_down)
add_stat("Overview", "Unchanged_genes", n_unchanged)
add_stat("Overview", "Total_DMRs", nrow(dmr_annot))

################################################################################
# Section 2: Promoter DMRs at DEGs
################################################################################

cat("\n--- Section 2: Promoter DMRs at DEGs ---\n")

# Extract promoter DMRs
promoter_dmrs <- dmr_annot %>%
    filter(grepl("Promoter", annotation, ignore.case = TRUE))

cat(sprintf("  Total promoter DMRs: %d\n", nrow(promoter_dmrs)))

# Aggregate by gene
promoter_by_gene <- promoter_dmrs %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarize(
        n_promoter_dmrs = n(),
        mean_logFC = mean(logFC, na.rm = TRUE),
        max_logFC = max(logFC, na.rm = TRUE),
        min_logFC = min(logFC, na.rm = TRUE),
        mean_FDR = mean(FDR, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        methylation_direction = case_when(
            mean_logFC > 0 ~ "hypermethylated",
            mean_logFC < 0 ~ "hypomethylated",
            TRUE ~ "no_change"
        )
    )

# Merge with expression data
promoter_expr <- promoter_by_gene %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )

cat(sprintf("  Genes with promoter DMRs matched to expression: %d\n", nrow(promoter_expr)))

# Statistics
n_genes_promoter_dmr <- nrow(promoter_by_gene)
n_degs_promoter_dmr <- sum(promoter_expr$is_DEG, na.rm = TRUE)
n_up_promoter_dmr <- sum(promoter_expr$expression_status == "upregulated", na.rm = TRUE)
n_down_promoter_dmr <- sum(promoter_expr$expression_status == "downregulated", na.rm = TRUE)

add_stat("Promoter_DMRs", "Total_promoter_DMRs", nrow(promoter_dmrs))
add_stat("Promoter_DMRs", "Unique_genes_with_promoter_DMRs", n_genes_promoter_dmr)
add_stat("Promoter_DMRs", "Genes_matched_to_expression", nrow(promoter_expr))
add_stat("Promoter_DMRs", "DEGs_with_promoter_DMRs", n_degs_promoter_dmr)
add_stat("Promoter_DMRs", "Upregulated_DEGs_with_promoter_DMRs", n_up_promoter_dmr)
add_stat("Promoter_DMRs", "Downregulated_DEGs_with_promoter_DMRs", n_down_promoter_dmr)
add_stat("Promoter_DMRs", "Pct_DEGs_with_promoter_DMRs",
         sprintf("%.1f%%", 100 * n_degs_promoter_dmr / (n_up + n_down)))

################################################################################
# Section 3: Direction Analysis (Promoter Methylation vs Expression)
################################################################################

cat("\n--- Section 3: Direction Analysis ---\n")

# DEGs with promoter DMRs
deg_promoter <- promoter_expr %>% filter(is_DEG)

# Cross-tabulation
hyper_down <- sum(deg_promoter$methylation_direction == "hypermethylated" &
                  deg_promoter$expression_status == "downregulated", na.rm = TRUE)
hyper_up <- sum(deg_promoter$methylation_direction == "hypermethylated" &
                deg_promoter$expression_status == "upregulated", na.rm = TRUE)
hypo_down <- sum(deg_promoter$methylation_direction == "hypomethylated" &
                 deg_promoter$expression_status == "downregulated", na.rm = TRUE)
hypo_up <- sum(deg_promoter$methylation_direction == "hypomethylated" &
               deg_promoter$expression_status == "upregulated", na.rm = TRUE)

cat(sprintf("  Hypermethylated + Downregulated: %d\n", hyper_down))
cat(sprintf("  Hypermethylated + Upregulated: %d\n", hyper_up))
cat(sprintf("  Hypomethylated + Downregulated: %d\n", hypo_down))
cat(sprintf("  Hypomethylated + Upregulated: %d\n", hypo_up))

add_stat("Direction_Analysis", "Hypermethylated_Downregulated", hyper_down)
add_stat("Direction_Analysis", "Hypermethylated_Upregulated", hyper_up)
add_stat("Direction_Analysis", "Hypomethylated_Downregulated", hypo_down)
add_stat("Direction_Analysis", "Hypomethylated_Upregulated", hypo_up)

# Percentages among hypermethylated DEGs
total_hyper_deg <- hyper_down + hyper_up
if (total_hyper_deg > 0) {
    add_stat("Direction_Analysis", "Pct_Hypermethylated_that_are_Down",
             sprintf("%.1f%%", 100 * hyper_down / total_hyper_deg))
    add_stat("Direction_Analysis", "Pct_Hypermethylated_that_are_Up",
             sprintf("%.1f%%", 100 * hyper_up / total_hyper_deg))
}

# Classical expectation test (Fisher's exact)
# Classical: hypermethylation -> downregulation
# H0: no association between methylation direction and expression direction
contingency <- matrix(c(hyper_down, hyper_up, hypo_down, hypo_up), nrow = 2,
                      dimnames = list(
                          Methylation = c("Hyper", "Hypo"),
                          Expression = c("Down", "Up")
                      ))

if (all(contingency > 0)) {
    fisher_result <- fisher.test(contingency)
    add_stat("Direction_Analysis", "Fisher_OR_HyperDown_vs_HyperUp",
             round(fisher_result$estimate, 3))
    add_stat("Direction_Analysis", "Fisher_pvalue",
             format(fisher_result$p.value, scientific = TRUE, digits = 3))
} else {
    # Use chi-squared with correction for small counts
    add_stat("Direction_Analysis", "Fisher_OR_HyperDown_vs_HyperUp", "N/A (zero cells)")
    add_stat("Direction_Analysis", "Fisher_pvalue", "N/A")
}

################################################################################
# Section 4: Methylation-Expression Correlations
################################################################################

cat("\n--- Section 4: Methylation-Expression Correlations ---\n")

# Function to calculate and store correlations
calc_correlations <- function(data, meth_col, expr_col, prefix) {
    valid <- !is.na(data[[meth_col]]) & !is.na(data[[expr_col]])
    n <- sum(valid)

    if (n < 10) {
        add_stat(prefix, "N_genes", n)
        add_stat(prefix, "Pearson_r", "N/A (n<10)")
        add_stat(prefix, "Spearman_rho", "N/A (n<10)")
        return(NULL)
    }

    x <- data[[meth_col]][valid]
    y <- data[[expr_col]][valid]

    pearson <- cor.test(x, y, method = "pearson")
    spearman <- suppressWarnings(cor.test(x, y, method = "spearman"))

    add_stat(prefix, "N_genes", n)
    add_stat(prefix, "Pearson_r", round(pearson$estimate, 4))
    add_stat(prefix, "Pearson_pvalue", format(pearson$p.value, scientific = TRUE, digits = 3))
    add_stat(prefix, "Spearman_rho", round(spearman$estimate, 4))
    add_stat(prefix, "Spearman_pvalue", format(spearman$p.value, scientific = TRUE, digits = 3))

    cat(sprintf("  %s (n=%d): r=%.4f, rho=%.4f\n", prefix, n,
                pearson$estimate, spearman$estimate))

    return(list(pearson = pearson, spearman = spearman, n = n))
}

# 4a. TSS ±5kb methylation (from integrated data)
cat("\n  TSS ±5kb Region:\n")

# All genes
calc_correlations(integrated, "meDIP_delta", "log2FoldChange",
                  "Correlation_TSS5kb_AllGenes")

# DEGs only
calc_correlations(integrated %>% filter(is_DEG), "meDIP_delta", "log2FoldChange",
                  "Correlation_TSS5kb_DEGs")

# Upregulated DEGs
calc_correlations(integrated %>% filter(expression_status == "upregulated"),
                  "meDIP_delta", "log2FoldChange",
                  "Correlation_TSS5kb_Upregulated")

# Downregulated DEGs
calc_correlations(integrated %>% filter(expression_status == "downregulated"),
                  "meDIP_delta", "log2FoldChange",
                  "Correlation_TSS5kb_Downregulated")

# 4b. Promoter-specific methylation
cat("\n  Promoter DMRs:\n")

# All genes with promoter DMRs
calc_correlations(promoter_expr, "mean_logFC", "log2FoldChange",
                  "Correlation_Promoter_AllGenes")

# DEGs with promoter DMRs
calc_correlations(promoter_expr %>% filter(is_DEG), "mean_logFC", "log2FoldChange",
                  "Correlation_Promoter_DEGs")

# Upregulated DEGs with promoter DMRs
calc_correlations(promoter_expr %>% filter(expression_status == "upregulated"),
                  "mean_logFC", "log2FoldChange",
                  "Correlation_Promoter_Upregulated")

# Downregulated DEGs with promoter DMRs
calc_correlations(promoter_expr %>% filter(expression_status == "downregulated"),
                  "mean_logFC", "log2FoldChange",
                  "Correlation_Promoter_Downregulated")

################################################################################
# Section 5: Methylation Signal by Expression Category
################################################################################

cat("\n--- Section 5: Methylation by Expression Category ---\n")

# TSS ±5kb region
meth_by_expr <- integrated %>%
    group_by(expression_status) %>%
    summarize(
        n = n(),
        mean_meDIP_delta = mean(meDIP_delta, na.rm = TRUE),
        median_meDIP_delta = median(meDIP_delta, na.rm = TRUE),
        sd_meDIP_delta = sd(meDIP_delta, na.rm = TRUE),
        .groups = "drop"
    )

for (i in 1:nrow(meth_by_expr)) {
    status <- meth_by_expr$expression_status[i]
    add_stat("Methylation_by_Expression", paste0("N_", status), meth_by_expr$n[i])
    add_stat("Methylation_by_Expression", paste0("Mean_meDIP_delta_", status),
             round(meth_by_expr$mean_meDIP_delta[i], 4))
    add_stat("Methylation_by_Expression", paste0("Median_meDIP_delta_", status),
             round(meth_by_expr$median_meDIP_delta[i], 4))
}

# Wilcoxon test: upregulated vs downregulated
up_meth <- integrated$meDIP_delta[integrated$expression_status == "upregulated"]
down_meth <- integrated$meDIP_delta[integrated$expression_status == "downregulated"]
wilcox_result <- wilcox.test(up_meth, down_meth)
add_stat("Methylation_by_Expression", "Wilcoxon_Up_vs_Down_pvalue",
         format(wilcox_result$p.value, scientific = TRUE, digits = 3))

cat(sprintf("  Upregulated: mean delta = %.4f (n=%d)\n",
            mean(up_meth, na.rm=TRUE), length(up_meth)))
cat(sprintf("  Downregulated: mean delta = %.4f (n=%d)\n",
            mean(down_meth, na.rm=TRUE), length(down_meth)))

################################################################################
# Section 6: Gene Body DMRs at DEGs
################################################################################

cat("\n--- Section 6: Gene Body DMRs at DEGs ---\n")

# Extract gene body DMRs (Intron + Exon)
genebody_dmrs <- dmr_annot %>%
    filter(grepl("Intron|Exon", annotation, ignore.case = TRUE))

cat(sprintf("  Total gene body DMRs: %d\n", nrow(genebody_dmrs)))

# Aggregate by gene
genebody_by_gene <- genebody_dmrs %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarize(
        n_genebody_dmrs = n(),
        mean_logFC = mean(logFC, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        methylation_direction = ifelse(mean_logFC > 0, "hypermethylated", "hypomethylated")
    )

# Merge with expression
genebody_expr <- genebody_by_gene %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )

n_degs_genebody <- sum(genebody_expr$is_DEG, na.rm = TRUE)

add_stat("GeneBody_DMRs", "Total_genebody_DMRs", nrow(genebody_dmrs))
add_stat("GeneBody_DMRs", "Unique_genes_with_genebody_DMRs", nrow(genebody_by_gene))
add_stat("GeneBody_DMRs", "DEGs_with_genebody_DMRs", n_degs_genebody)

# Correlation for gene body
cat("\n  Gene Body Methylation Correlations:\n")
calc_correlations(genebody_expr %>% filter(is_DEG), "mean_logFC", "log2FoldChange",
                  "Correlation_GeneBody_DEGs")

################################################################################
# Section 7: DMR Enrichment at DEGs vs Non-DEGs
################################################################################

cat("\n--- Section 7: DMR Enrichment Analysis ---\n")

# Genes with ANY DMR
genes_with_dmr <- unique(c(
    dmr_annot$SYMBOL[!is.na(dmr_annot$SYMBOL) & dmr_annot$SYMBOL != ""]
))

# Build contingency table
degs_with_dmr <- sum(deseq$gene_symbol %in% genes_with_dmr & deseq$is_DEG, na.rm = TRUE)
degs_without_dmr <- sum(!deseq$gene_symbol %in% genes_with_dmr & deseq$is_DEG, na.rm = TRUE)
nondegs_with_dmr <- sum(deseq$gene_symbol %in% genes_with_dmr & !deseq$is_DEG, na.rm = TRUE)
nondegs_without_dmr <- sum(!deseq$gene_symbol %in% genes_with_dmr & !deseq$is_DEG, na.rm = TRUE)

add_stat("DMR_Enrichment", "DEGs_with_any_DMR", degs_with_dmr)
add_stat("DMR_Enrichment", "DEGs_without_DMR", degs_without_dmr)
add_stat("DMR_Enrichment", "NonDEGs_with_any_DMR", nondegs_with_dmr)
add_stat("DMR_Enrichment", "NonDEGs_without_DMR", nondegs_without_dmr)

pct_degs_dmr <- 100 * degs_with_dmr / (degs_with_dmr + degs_without_dmr)
pct_nondegs_dmr <- 100 * nondegs_with_dmr / (nondegs_with_dmr + nondegs_without_dmr)

add_stat("DMR_Enrichment", "Pct_DEGs_with_DMR", sprintf("%.1f%%", pct_degs_dmr))
add_stat("DMR_Enrichment", "Pct_NonDEGs_with_DMR", sprintf("%.1f%%", pct_nondegs_dmr))

# Fisher's exact test
fisher_dmr <- fisher.test(matrix(c(degs_with_dmr, degs_without_dmr,
                                    nondegs_with_dmr, nondegs_without_dmr), nrow = 2))
add_stat("DMR_Enrichment", "Fisher_OR", round(fisher_dmr$estimate, 3))
add_stat("DMR_Enrichment", "Fisher_pvalue", format(fisher_dmr$p.value, scientific = TRUE, digits = 3))
add_stat("DMR_Enrichment", "Fisher_95CI_low", round(fisher_dmr$conf.int[1], 3))
add_stat("DMR_Enrichment", "Fisher_95CI_high", round(fisher_dmr$conf.int[2], 3))

cat(sprintf("  DEGs with DMR: %.1f%% (%d/%d)\n", pct_degs_dmr,
            degs_with_dmr, degs_with_dmr + degs_without_dmr))
cat(sprintf("  Non-DEGs with DMR: %.1f%% (%d/%d)\n", pct_nondegs_dmr,
            nondegs_with_dmr, nondegs_with_dmr + nondegs_without_dmr))
cat(sprintf("  Fisher's OR: %.3f (p = %s)\n", fisher_dmr$estimate,
            format(fisher_dmr$p.value, scientific = TRUE, digits = 3)))

################################################################################
# Section 8: Promoter DMR Enrichment at DEGs
################################################################################

cat("\n--- Section 8: Promoter DMR Enrichment ---\n")

genes_with_promoter_dmr <- unique(promoter_by_gene$SYMBOL)

degs_prom_dmr <- sum(deseq$gene_symbol %in% genes_with_promoter_dmr & deseq$is_DEG, na.rm = TRUE)
degs_no_prom_dmr <- sum(!deseq$gene_symbol %in% genes_with_promoter_dmr & deseq$is_DEG, na.rm = TRUE)
nondegs_prom_dmr <- sum(deseq$gene_symbol %in% genes_with_promoter_dmr & !deseq$is_DEG, na.rm = TRUE)
nondegs_no_prom_dmr <- sum(!deseq$gene_symbol %in% genes_with_promoter_dmr & !deseq$is_DEG, na.rm = TRUE)

add_stat("Promoter_DMR_Enrichment", "DEGs_with_promoter_DMR", degs_prom_dmr)
add_stat("Promoter_DMR_Enrichment", "DEGs_without_promoter_DMR", degs_no_prom_dmr)
add_stat("Promoter_DMR_Enrichment", "NonDEGs_with_promoter_DMR", nondegs_prom_dmr)
add_stat("Promoter_DMR_Enrichment", "NonDEGs_without_promoter_DMR", nondegs_no_prom_dmr)

pct_degs_prom <- 100 * degs_prom_dmr / (degs_prom_dmr + degs_no_prom_dmr)
pct_nondegs_prom <- 100 * nondegs_prom_dmr / (nondegs_prom_dmr + nondegs_no_prom_dmr)

add_stat("Promoter_DMR_Enrichment", "Pct_DEGs_with_promoter_DMR", sprintf("%.1f%%", pct_degs_prom))
add_stat("Promoter_DMR_Enrichment", "Pct_NonDEGs_with_promoter_DMR", sprintf("%.1f%%", pct_nondegs_prom))

fisher_prom <- fisher.test(matrix(c(degs_prom_dmr, degs_no_prom_dmr,
                                     nondegs_prom_dmr, nondegs_no_prom_dmr), nrow = 2))
add_stat("Promoter_DMR_Enrichment", "Fisher_OR", round(fisher_prom$estimate, 3))
add_stat("Promoter_DMR_Enrichment", "Fisher_pvalue", format(fisher_prom$p.value, scientific = TRUE, digits = 3))

cat(sprintf("  DEGs with promoter DMR: %.1f%%\n", pct_degs_prom))
cat(sprintf("  Non-DEGs with promoter DMR: %.1f%%\n", pct_nondegs_prom))
cat(sprintf("  Fisher's OR: %.3f (p = %s)\n", fisher_prom$estimate,
            format(fisher_prom$p.value, scientific = TRUE, digits = 3)))

################################################################################
# Combine and Write Statistics
################################################################################

cat("\n--- Writing Output Files ---\n")

# Combine all statistics
all_stats <- do.call(rbind, stats_list)

# Write long format
output_file <- file.path(OUTPUT_DIR, "methylation_at_DEGs_summary.csv")
write.csv(all_stats, output_file, row.names = FALSE)
cat(sprintf("  Written: %s\n", output_file))

# Write wide format
stats_wide <- all_stats %>%
    mutate(Full_Metric = paste(Category, Metric, sep = "_")) %>%
    select(Full_Metric, Value) %>%
    pivot_wider(names_from = Full_Metric, values_from = Value)

output_wide <- file.path(OUTPUT_DIR, "methylation_at_DEGs_summary_wide.csv")
write.csv(stats_wide, output_wide, row.names = FALSE)
cat(sprintf("  Written: %s\n", output_wide))

################################################################################
# Generate Top Genes Tables
################################################################################

cat("\n--- Generating Top Genes Tables ---\n")

# Top DEGs with strongest promoter hypermethylation
top_hyper_down <- deg_promoter %>%
    filter(methylation_direction == "hypermethylated" & expression_status == "downregulated") %>%
    arrange(desc(mean_logFC)) %>%
    head(50) %>%
    mutate(category = "Hypermethylated_Downregulated")

top_hyper_up <- deg_promoter %>%
    filter(methylation_direction == "hypermethylated" & expression_status == "upregulated") %>%
    arrange(desc(mean_logFC)) %>%
    head(50) %>%
    mutate(category = "Hypermethylated_Upregulated")

top_hypo_down <- deg_promoter %>%
    filter(methylation_direction == "hypomethylated" & expression_status == "downregulated") %>%
    arrange(mean_logFC) %>%
    head(50) %>%
    mutate(category = "Hypomethylated_Downregulated")

top_hypo_up <- deg_promoter %>%
    filter(methylation_direction == "hypomethylated" & expression_status == "upregulated") %>%
    arrange(mean_logFC) %>%
    head(50) %>%
    mutate(category = "Hypomethylated_Upregulated")

# Combine
top_genes <- bind_rows(top_hyper_down, top_hyper_up, top_hypo_down, top_hypo_up) %>%
    select(category, SYMBOL, n_promoter_dmrs, mean_logFC, log2FoldChange, padj)

output_top <- file.path(OUTPUT_DIR, "methylation_at_DEGs_top_genes.csv")
write.csv(top_genes, output_top, row.names = FALSE)
cat(sprintf("  Written: %s (%d genes)\n", output_top, nrow(top_genes)))

################################################################################
# Summary
################################################################################

cat("\n========================================\n")
cat("SUMMARY COMPLETE\n")
cat("========================================\n\n")

cat("Key Findings:\n")
cat(sprintf("  - Total DEGs: %d (Up: %d, Down: %d)\n", n_up + n_down, n_up, n_down))
cat(sprintf("  - DEGs with promoter DMRs: %d (%.1f%%)\n",
            n_degs_promoter_dmr, 100 * n_degs_promoter_dmr / (n_up + n_down)))
cat(sprintf("  - Hypermethylated promoter + Down: %d\n", hyper_down))
cat(sprintf("  - Hypermethylated promoter + Up: %d\n", hyper_up))
cat(sprintf("  - Promoter DMR-expression correlation (DEGs): shown above\n"))
cat(sprintf("\nOutput files:\n"))
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", output_wide))
cat(sprintf("  - %s\n", output_top))

cat("\n========================================\n")
