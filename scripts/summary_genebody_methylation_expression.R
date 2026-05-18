#!/usr/bin/env Rscript
################################################################################
# Gene Body Methylation-Expression Analysis
#
# This script analyzes the relationship between gene body methylation and
# gene expression, including:
# 1. Gene body correlation split by upregulated vs downregulated DEGs
# 2. Direction analysis (hypermethylated/hypomethylated gene body + up/down)
# 3. Comparison with promoter methylation effects
#
# Author: Generated for SRF_Eva project
# Date: 2024-12
#
# Usage: Rscript summary_genebody_methylation_expression.R
# Output: results/summary/genebody_methylation_expression_summary.csv
#         results/summary/genebody_methylation_expression_plots.pdf
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(gridExtra)
})

cat("========================================\n")
cat("Gene Body Methylation-Expression Analysis\n")
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

cat("1. Loading data files...\n")

# Load DESeq2 results
deseq_file <- file.path(RNA_DIR, "results/05_deseq2/deseq2_results_TES_vs_GFP.txt")
if (!file.exists(deseq_file)) {
    stop("DESeq2 results not found: ", deseq_file)
}
deseq <- read.delim(deseq_file, stringsAsFactors = FALSE)
deseq$gene_id_clean <- gsub("\\..*", "", deseq$gene_id)
cat(sprintf("   Loaded %d genes from DESeq2\n", nrow(deseq)))

# Load annotated DMRs
dmr_file <- file.path(MEDIP_DIR, "results/08_annotation/TES_vs_GFP_annotated.csv")
if (!file.exists(dmr_file)) {
    stop("Annotated DMRs not found: ", dmr_file)
}
dmr_annot <- read.csv(dmr_file, stringsAsFactors = FALSE)
cat(sprintf("   Loaded %d annotated DMRs\n", nrow(dmr_annot)))

################################################################################
# Define Gene Categories
################################################################################

cat("\n2. Defining gene categories...\n")

deseq <- deseq %>%
    mutate(
        is_DEG = !is.na(padj) & padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD,
        expression_status = case_when(
            is_DEG & log2FoldChange > 0 ~ "Upregulated",
            is_DEG & log2FoldChange < 0 ~ "Downregulated",
            TRUE ~ "Unchanged"
        )
    )

n_up <- sum(deseq$expression_status == "Upregulated", na.rm = TRUE)
n_down <- sum(deseq$expression_status == "Downregulated", na.rm = TRUE)
cat(sprintf("   Upregulated DEGs: %d\n", n_up))
cat(sprintf("   Downregulated DEGs: %d\n", n_down))

################################################################################
# Extract Gene Body and Promoter DMRs
################################################################################

cat("\n3. Extracting genomic region DMRs...\n")

# Gene body DMRs (Intron + Exon)
genebody_dmrs <- dmr_annot %>%
    filter(grepl("Intron|Exon", annotation, ignore.case = TRUE))
cat(sprintf("   Gene body DMRs: %d\n", nrow(genebody_dmrs)))

# Promoter DMRs
promoter_dmrs <- dmr_annot %>%
    filter(grepl("Promoter", annotation, ignore.case = TRUE))
cat(sprintf("   Promoter DMRs: %d\n", nrow(promoter_dmrs)))

################################################################################
# Aggregate Methylation by Gene
################################################################################

cat("\n4. Aggregating methylation by gene...\n")

# Gene body aggregation
genebody_by_gene <- genebody_dmrs %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarize(
        n_genebody_dmrs = n(),
        mean_logFC_genebody = mean(logFC, na.rm = TRUE),
        median_logFC_genebody = median(logFC, na.rm = TRUE),
        max_logFC_genebody = max(logFC, na.rm = TRUE),
        min_logFC_genebody = min(logFC, na.rm = TRUE),
        sum_logFC_genebody = sum(logFC, na.rm = TRUE),
        mean_FDR_genebody = mean(FDR, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        genebody_meth_direction = case_when(
            mean_logFC_genebody > 0.5 ~ "Hypermethylated",
            mean_logFC_genebody < -0.5 ~ "Hypomethylated",
            TRUE ~ "No_change"
        )
    )
cat(sprintf("   Genes with gene body DMRs: %d\n", nrow(genebody_by_gene)))

# Promoter aggregation
promoter_by_gene <- promoter_dmrs %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarize(
        n_promoter_dmrs = n(),
        mean_logFC_promoter = mean(logFC, na.rm = TRUE),
        median_logFC_promoter = median(logFC, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        promoter_meth_direction = case_when(
            mean_logFC_promoter > 0.5 ~ "Hypermethylated",
            mean_logFC_promoter < -0.5 ~ "Hypomethylated",
            TRUE ~ "No_change"
        )
    )
cat(sprintf("   Genes with promoter DMRs: %d\n", nrow(promoter_by_gene)))

################################################################################
# Merge with Expression Data
################################################################################

cat("\n5. Merging methylation with expression data...\n")

# Gene body + expression
genebody_expr <- genebody_by_gene %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, baseMean,
                         expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )
cat(sprintf("   Gene body DMRs matched to expression: %d genes\n", nrow(genebody_expr)))

# Promoter + expression
promoter_expr <- promoter_by_gene %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, baseMean,
                         expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )
cat(sprintf("   Promoter DMRs matched to expression: %d genes\n", nrow(promoter_expr)))

# Combined (genes with both promoter and gene body DMRs)
combined_expr <- genebody_by_gene %>%
    inner_join(promoter_by_gene, by = "SYMBOL") %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, baseMean,
                         expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )
cat(sprintf("   Genes with both promoter AND gene body DMRs: %d\n", nrow(combined_expr)))

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
# Section 1: Overview
################################################################################

cat("\n6. Calculating statistics...\n")
cat("   --- Overview ---\n")

add_stat("Overview", "Total_DEGs", n_up + n_down)
add_stat("Overview", "Upregulated_DEGs", n_up)
add_stat("Overview", "Downregulated_DEGs", n_down)
add_stat("Overview", "Total_genebody_DMRs", nrow(genebody_dmrs))
add_stat("Overview", "Total_promoter_DMRs", nrow(promoter_dmrs))
add_stat("Overview", "Genes_with_genebody_DMRs", nrow(genebody_by_gene))
add_stat("Overview", "Genes_with_promoter_DMRs", nrow(promoter_by_gene))

################################################################################
# Section 2: Gene Body Correlation by Expression Direction
################################################################################

cat("   --- Gene Body Correlations ---\n")

calc_and_store_correlation <- function(data, meth_col, expr_col, prefix) {
    valid <- !is.na(data[[meth_col]]) & !is.na(data[[expr_col]])
    n <- sum(valid)

    if (n < 10) {
        add_stat(prefix, "N_genes", n)
        add_stat(prefix, "Pearson_r", "N/A")
        add_stat(prefix, "Spearman_rho", "N/A")
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

    cat(sprintf("      %s: n=%d, r=%.4f (p=%s), rho=%.4f (p=%s)\n",
                prefix, n, pearson$estimate,
                format(pearson$p.value, scientific = TRUE, digits = 2),
                spearman$estimate,
                format(spearman$p.value, scientific = TRUE, digits = 2)))

    return(list(pearson = pearson, spearman = spearman, n = n))
}

# Gene body correlations
calc_and_store_correlation(genebody_expr, "mean_logFC_genebody", "log2FoldChange",
                           "GeneBody_AllGenes")

calc_and_store_correlation(genebody_expr %>% filter(is_DEG),
                           "mean_logFC_genebody", "log2FoldChange",
                           "GeneBody_DEGs")

calc_and_store_correlation(genebody_expr %>% filter(expression_status == "Upregulated"),
                           "mean_logFC_genebody", "log2FoldChange",
                           "GeneBody_Upregulated")

calc_and_store_correlation(genebody_expr %>% filter(expression_status == "Downregulated"),
                           "mean_logFC_genebody", "log2FoldChange",
                           "GeneBody_Downregulated")

# Promoter correlations for comparison
cat("   --- Promoter Correlations (for comparison) ---\n")

calc_and_store_correlation(promoter_expr, "mean_logFC_promoter", "log2FoldChange",
                           "Promoter_AllGenes")

calc_and_store_correlation(promoter_expr %>% filter(is_DEG),
                           "mean_logFC_promoter", "log2FoldChange",
                           "Promoter_DEGs")

calc_and_store_correlation(promoter_expr %>% filter(expression_status == "Upregulated"),
                           "mean_logFC_promoter", "log2FoldChange",
                           "Promoter_Upregulated")

calc_and_store_correlation(promoter_expr %>% filter(expression_status == "Downregulated"),
                           "mean_logFC_promoter", "log2FoldChange",
                           "Promoter_Downregulated")

################################################################################
# Section 3: Direction Analysis - Gene Body
################################################################################

cat("   --- Gene Body Direction Analysis ---\n")

# DEGs with gene body DMRs
deg_genebody <- genebody_expr %>% filter(is_DEG)

# Cross-tabulation: gene body methylation direction vs expression direction
gb_hyper_down <- sum(deg_genebody$genebody_meth_direction == "Hypermethylated" &
                     deg_genebody$expression_status == "Downregulated", na.rm = TRUE)
gb_hyper_up <- sum(deg_genebody$genebody_meth_direction == "Hypermethylated" &
                   deg_genebody$expression_status == "Upregulated", na.rm = TRUE)
gb_hypo_down <- sum(deg_genebody$genebody_meth_direction == "Hypomethylated" &
                    deg_genebody$expression_status == "Downregulated", na.rm = TRUE)
gb_hypo_up <- sum(deg_genebody$genebody_meth_direction == "Hypomethylated" &
                  deg_genebody$expression_status == "Upregulated", na.rm = TRUE)
gb_nochange_down <- sum(deg_genebody$genebody_meth_direction == "No_change" &
                        deg_genebody$expression_status == "Downregulated", na.rm = TRUE)
gb_nochange_up <- sum(deg_genebody$genebody_meth_direction == "No_change" &
                      deg_genebody$expression_status == "Upregulated", na.rm = TRUE)

add_stat("GeneBody_Direction", "Hypermethylated_Downregulated", gb_hyper_down)
add_stat("GeneBody_Direction", "Hypermethylated_Upregulated", gb_hyper_up)
add_stat("GeneBody_Direction", "Hypomethylated_Downregulated", gb_hypo_down)
add_stat("GeneBody_Direction", "Hypomethylated_Upregulated", gb_hypo_up)
add_stat("GeneBody_Direction", "NoChange_Downregulated", gb_nochange_down)
add_stat("GeneBody_Direction", "NoChange_Upregulated", gb_nochange_up)

cat(sprintf("      Hypermethylated + Down: %d\n", gb_hyper_down))
cat(sprintf("      Hypermethylated + Up: %d\n", gb_hyper_up))
cat(sprintf("      Hypomethylated + Down: %d\n", gb_hypo_down))
cat(sprintf("      Hypomethylated + Up: %d\n", gb_hypo_up))

# Percentages
total_gb_hyper <- gb_hyper_down + gb_hyper_up
total_gb_hypo <- gb_hypo_down + gb_hypo_up

if (total_gb_hyper > 0) {
    add_stat("GeneBody_Direction", "Pct_Hypermethylated_Down",
             sprintf("%.1f%%", 100 * gb_hyper_down / total_gb_hyper))
    add_stat("GeneBody_Direction", "Pct_Hypermethylated_Up",
             sprintf("%.1f%%", 100 * gb_hyper_up / total_gb_hyper))
}

if (total_gb_hypo > 0) {
    add_stat("GeneBody_Direction", "Pct_Hypomethylated_Down",
             sprintf("%.1f%%", 100 * gb_hypo_down / total_gb_hypo))
    add_stat("GeneBody_Direction", "Pct_Hypomethylated_Up",
             sprintf("%.1f%%", 100 * gb_hypo_up / total_gb_hypo))
}

# Fisher's exact test for association
if (gb_hyper_down > 0 && gb_hyper_up > 0 && gb_hypo_down > 0 && gb_hypo_up > 0) {
    contingency_gb <- matrix(c(gb_hyper_down, gb_hyper_up, gb_hypo_down, gb_hypo_up),
                             nrow = 2,
                             dimnames = list(
                                 GeneBody_Meth = c("Hyper", "Hypo"),
                                 Expression = c("Down", "Up")
                             ))
    fisher_gb <- fisher.test(contingency_gb)
    add_stat("GeneBody_Direction", "Fisher_OR", round(fisher_gb$estimate, 3))
    add_stat("GeneBody_Direction", "Fisher_pvalue",
             format(fisher_gb$p.value, scientific = TRUE, digits = 3))
    cat(sprintf("      Fisher's OR: %.3f (p = %s)\n",
                fisher_gb$estimate, format(fisher_gb$p.value, scientific = TRUE, digits = 2)))
} else {
    add_stat("GeneBody_Direction", "Fisher_OR", "N/A (zero cells)")
    add_stat("GeneBody_Direction", "Fisher_pvalue", "N/A")
}

################################################################################
# Section 4: Direction Analysis - Promoter (for comparison)
################################################################################

cat("   --- Promoter Direction Analysis (for comparison) ---\n")

deg_promoter <- promoter_expr %>% filter(is_DEG)

pr_hyper_down <- sum(deg_promoter$promoter_meth_direction == "Hypermethylated" &
                     deg_promoter$expression_status == "Downregulated", na.rm = TRUE)
pr_hyper_up <- sum(deg_promoter$promoter_meth_direction == "Hypermethylated" &
                   deg_promoter$expression_status == "Upregulated", na.rm = TRUE)
pr_hypo_down <- sum(deg_promoter$promoter_meth_direction == "Hypomethylated" &
                    deg_promoter$expression_status == "Downregulated", na.rm = TRUE)
pr_hypo_up <- sum(deg_promoter$promoter_meth_direction == "Hypomethylated" &
                  deg_promoter$expression_status == "Upregulated", na.rm = TRUE)

add_stat("Promoter_Direction", "Hypermethylated_Downregulated", pr_hyper_down)
add_stat("Promoter_Direction", "Hypermethylated_Upregulated", pr_hyper_up)
add_stat("Promoter_Direction", "Hypomethylated_Downregulated", pr_hypo_down)
add_stat("Promoter_Direction", "Hypomethylated_Upregulated", pr_hypo_up)

cat(sprintf("      Hypermethylated + Down: %d\n", pr_hyper_down))
cat(sprintf("      Hypermethylated + Up: %d\n", pr_hyper_up))

total_pr_hyper <- pr_hyper_down + pr_hyper_up
if (total_pr_hyper > 0) {
    add_stat("Promoter_Direction", "Pct_Hypermethylated_Down",
             sprintf("%.1f%%", 100 * pr_hyper_down / total_pr_hyper))
    add_stat("Promoter_Direction", "Pct_Hypermethylated_Up",
             sprintf("%.1f%%", 100 * pr_hyper_up / total_pr_hyper))
}

################################################################################
# Section 5: Mean Methylation by Expression Category
################################################################################

cat("   --- Mean Methylation by Expression Category ---\n")

# Gene body methylation by expression status
gb_meth_by_expr <- genebody_expr %>%
    group_by(expression_status) %>%
    summarize(
        n = n(),
        mean_logFC = mean(mean_logFC_genebody, na.rm = TRUE),
        median_logFC = median(mean_logFC_genebody, na.rm = TRUE),
        sd_logFC = sd(mean_logFC_genebody, na.rm = TRUE),
        .groups = "drop"
    )

for (i in 1:nrow(gb_meth_by_expr)) {
    status <- gb_meth_by_expr$expression_status[i]
    add_stat("GeneBody_MethByExpr", paste0("N_", status), gb_meth_by_expr$n[i])
    add_stat("GeneBody_MethByExpr", paste0("Mean_logFC_", status),
             round(gb_meth_by_expr$mean_logFC[i], 4))
    add_stat("GeneBody_MethByExpr", paste0("Median_logFC_", status),
             round(gb_meth_by_expr$median_logFC[i], 4))
}

# Wilcoxon test
up_gb_meth <- genebody_expr$mean_logFC_genebody[genebody_expr$expression_status == "Upregulated"]
down_gb_meth <- genebody_expr$mean_logFC_genebody[genebody_expr$expression_status == "Downregulated"]

if (length(up_gb_meth) > 5 && length(down_gb_meth) > 5) {
    wilcox_gb <- wilcox.test(up_gb_meth, down_gb_meth)
    add_stat("GeneBody_MethByExpr", "Wilcoxon_Up_vs_Down_pvalue",
             format(wilcox_gb$p.value, scientific = TRUE, digits = 3))
    cat(sprintf("      Wilcoxon Up vs Down: p = %s\n",
                format(wilcox_gb$p.value, scientific = TRUE, digits = 2)))
}

# Promoter methylation by expression status (for comparison)
pr_meth_by_expr <- promoter_expr %>%
    group_by(expression_status) %>%
    summarize(
        n = n(),
        mean_logFC = mean(mean_logFC_promoter, na.rm = TRUE),
        median_logFC = median(mean_logFC_promoter, na.rm = TRUE),
        .groups = "drop"
    )

for (i in 1:nrow(pr_meth_by_expr)) {
    status <- pr_meth_by_expr$expression_status[i]
    add_stat("Promoter_MethByExpr", paste0("N_", status), pr_meth_by_expr$n[i])
    add_stat("Promoter_MethByExpr", paste0("Mean_logFC_", status),
             round(pr_meth_by_expr$mean_logFC[i], 4))
}

################################################################################
# Section 6: Comparison Summary
################################################################################

cat("   --- Comparison Summary ---\n")

# DEGs with gene body DMRs
n_degs_gb <- sum(genebody_expr$is_DEG)
n_up_gb <- sum(genebody_expr$expression_status == "Upregulated")
n_down_gb <- sum(genebody_expr$expression_status == "Downregulated")

add_stat("Comparison", "DEGs_with_genebody_DMRs", n_degs_gb)
add_stat("Comparison", "Upregulated_with_genebody_DMRs", n_up_gb)
add_stat("Comparison", "Downregulated_with_genebody_DMRs", n_down_gb)
add_stat("Comparison", "Pct_DEGs_with_genebody_DMRs",
         sprintf("%.1f%%", 100 * n_degs_gb / (n_up + n_down)))

# DEGs with promoter DMRs
n_degs_pr <- sum(promoter_expr$is_DEG)
add_stat("Comparison", "DEGs_with_promoter_DMRs", n_degs_pr)
add_stat("Comparison", "Pct_DEGs_with_promoter_DMRs",
         sprintf("%.1f%%", 100 * n_degs_pr / (n_up + n_down)))

################################################################################
# Save Statistics
################################################################################

cat("\n7. Saving results...\n")

all_stats <- do.call(rbind, stats_list)

# Long format
output_file <- file.path(OUTPUT_DIR, "genebody_methylation_expression_summary.csv")
write.csv(all_stats, output_file, row.names = FALSE)
cat(sprintf("   Saved: %s\n", output_file))

# Wide format
stats_wide <- all_stats %>%
    mutate(Full_Metric = paste(Category, Metric, sep = "_")) %>%
    select(Full_Metric, Value) %>%
    pivot_wider(names_from = Full_Metric, values_from = Value)

output_wide <- file.path(OUTPUT_DIR, "genebody_methylation_expression_wide.csv")
write.csv(stats_wide, output_wide, row.names = FALSE)
cat(sprintf("   Saved: %s\n", output_wide))

################################################################################
# Generate Plots
################################################################################

cat("\n8. Generating plots...\n")

pdf_file <- file.path(OUTPUT_DIR, "genebody_methylation_expression_plots.pdf")
pdf(pdf_file, width = 12, height = 10)

# Plot 1: Scatter plot - Gene body methylation vs expression (DEGs)
p1 <- ggplot(genebody_expr %>% filter(is_DEG),
             aes(x = mean_logFC_genebody, y = log2FoldChange, color = expression_status)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    labs(title = "Gene Body Methylation vs Expression (DEGs)",
         subtitle = sprintf("n = %d genes", sum(genebody_expr$is_DEG)),
         x = "Gene Body Methylation (mean log2FC)",
         y = "Expression (log2FC TES/GFP)",
         color = "Expression\nStatus") +
    theme_bw() +
    theme(legend.position = "right")

# Plot 2: Scatter plot - Promoter methylation vs expression (DEGs)
p2 <- ggplot(promoter_expr %>% filter(is_DEG),
             aes(x = mean_logFC_promoter, y = log2FoldChange, color = expression_status)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    labs(title = "Promoter Methylation vs Expression (DEGs)",
         subtitle = sprintf("n = %d genes", sum(promoter_expr$is_DEG)),
         x = "Promoter Methylation (mean log2FC)",
         y = "Expression (log2FC TES/GFP)",
         color = "Expression\nStatus") +
    theme_bw() +
    theme(legend.position = "right")

grid.arrange(p1, p2, ncol = 2)

# Plot 3: Direction analysis bar plot - Gene body
direction_data_gb <- data.frame(
    Methylation = c("Hypermethylated", "Hypermethylated", "Hypomethylated", "Hypomethylated"),
    Expression = c("Downregulated", "Upregulated", "Downregulated", "Upregulated"),
    Count = c(gb_hyper_down, gb_hyper_up, gb_hypo_down, gb_hypo_up)
)

p3 <- ggplot(direction_data_gb, aes(x = Methylation, y = Count, fill = Expression)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3) +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    labs(title = "Gene Body Methylation vs Expression Direction",
         subtitle = "DEGs with gene body DMRs",
         x = "Gene Body Methylation Direction",
         y = "Number of Genes") +
    theme_bw() +
    theme(legend.position = "right")

# Plot 4: Direction analysis bar plot - Promoter
direction_data_pr <- data.frame(
    Methylation = c("Hypermethylated", "Hypermethylated", "Hypomethylated", "Hypomethylated"),
    Expression = c("Downregulated", "Upregulated", "Downregulated", "Upregulated"),
    Count = c(pr_hyper_down, pr_hyper_up, pr_hypo_down, pr_hypo_up)
)

p4 <- ggplot(direction_data_pr, aes(x = Methylation, y = Count, fill = Expression)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3) +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    labs(title = "Promoter Methylation vs Expression Direction",
         subtitle = "DEGs with promoter DMRs",
         x = "Promoter Methylation Direction",
         y = "Number of Genes") +
    theme_bw() +
    theme(legend.position = "right")

grid.arrange(p3, p4, ncol = 2)

# Plot 5: Boxplot - Gene body methylation by expression category
p5 <- ggplot(genebody_expr %>% filter(expression_status != "Unchanged"),
             aes(x = expression_status, y = mean_logFC_genebody, fill = expression_status)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    labs(title = "Gene Body Methylation by Expression Status",
         x = "Expression Status",
         y = "Gene Body Methylation (mean log2FC)") +
    theme_bw() +
    theme(legend.position = "none")

# Plot 6: Boxplot - Promoter methylation by expression category
p6 <- ggplot(promoter_expr %>% filter(expression_status != "Unchanged"),
             aes(x = expression_status, y = mean_logFC_promoter, fill = expression_status)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    labs(title = "Promoter Methylation by Expression Status",
         x = "Expression Status",
         y = "Promoter Methylation (mean log2FC)") +
    theme_bw() +
    theme(legend.position = "none")

grid.arrange(p5, p6, ncol = 2)

# Plot 7: Correlation comparison bar plot
cor_data <- data.frame(
    Region = c("Gene Body", "Gene Body", "Gene Body", "Gene Body",
               "Promoter", "Promoter", "Promoter", "Promoter"),
    Subset = c("All Genes", "DEGs", "Upregulated", "Downregulated",
               "All Genes", "DEGs", "Upregulated", "Downregulated"),
    stringsAsFactors = FALSE
)

# Extract correlation values from stats
get_cor <- function(prefix) {
    idx <- which(all_stats$Category == prefix & all_stats$Metric == "Spearman_rho")
    if (length(idx) > 0) {
        val <- all_stats$Value[idx]
        if (val != "N/A") return(as.numeric(val))
    }
    return(NA)
}

cor_data$Correlation <- c(
    get_cor("GeneBody_AllGenes"),
    get_cor("GeneBody_DEGs"),
    get_cor("GeneBody_Upregulated"),
    get_cor("GeneBody_Downregulated"),
    get_cor("Promoter_AllGenes"),
    get_cor("Promoter_DEGs"),
    get_cor("Promoter_Upregulated"),
    get_cor("Promoter_Downregulated")
)

cor_data$Subset <- factor(cor_data$Subset,
                          levels = c("All Genes", "DEGs", "Upregulated", "Downregulated"))

p7 <- ggplot(cor_data %>% filter(!is.na(Correlation)),
             aes(x = Subset, y = Correlation, fill = Region)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Gene Body" = "darkgreen", "Promoter" = "purple")) +
    labs(title = "Methylation-Expression Correlation Comparison",
         subtitle = "Spearman correlation (rho)",
         x = "Gene Subset",
         y = "Correlation (rho)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p7)

dev.off()
cat(sprintf("   Saved: %s\n", pdf_file))

################################################################################
# Top Genes Tables
################################################################################

cat("\n9. Generating top genes tables...\n")

# Top genes with hypermethylated gene body + downregulated
top_gb_hyper_down <- deg_genebody %>%
    filter(genebody_meth_direction == "Hypermethylated" &
           expression_status == "Downregulated") %>%
    arrange(desc(mean_logFC_genebody)) %>%
    head(50) %>%
    select(SYMBOL, n_genebody_dmrs, mean_logFC_genebody, log2FoldChange, padj) %>%
    mutate(Category = "GeneBody_Hyper_Down")

# Top genes with hypermethylated gene body + upregulated
top_gb_hyper_up <- deg_genebody %>%
    filter(genebody_meth_direction == "Hypermethylated" &
           expression_status == "Upregulated") %>%
    arrange(desc(mean_logFC_genebody)) %>%
    head(50) %>%
    select(SYMBOL, n_genebody_dmrs, mean_logFC_genebody, log2FoldChange, padj) %>%
    mutate(Category = "GeneBody_Hyper_Up")

# Top genes with hypomethylated gene body
top_gb_hypo <- deg_genebody %>%
    filter(genebody_meth_direction == "Hypomethylated") %>%
    arrange(mean_logFC_genebody) %>%
    head(50) %>%
    select(SYMBOL, n_genebody_dmrs, mean_logFC_genebody, log2FoldChange, padj) %>%
    mutate(Category = "GeneBody_Hypo")

top_genes <- bind_rows(top_gb_hyper_down, top_gb_hyper_up, top_gb_hypo)

top_genes_file <- file.path(OUTPUT_DIR, "genebody_methylation_top_genes.csv")
write.csv(top_genes, top_genes_file, row.names = FALSE)
cat(sprintf("   Saved: %s (%d genes)\n", top_genes_file, nrow(top_genes)))

################################################################################
# Summary
################################################################################

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Key Findings:\n")
cat(sprintf("  Gene Body Correlations:\n"))
cat(sprintf("    - All genes: rho = %.4f\n", get_cor("GeneBody_AllGenes")))
cat(sprintf("    - DEGs: rho = %.4f\n", get_cor("GeneBody_DEGs")))
cat(sprintf("    - Upregulated: rho = %.4f\n", get_cor("GeneBody_Upregulated")))
cat(sprintf("    - Downregulated: rho = %.4f\n", get_cor("GeneBody_Downregulated")))
cat(sprintf("\n  Gene Body Direction (DEGs):\n"))
cat(sprintf("    - Hypermethylated + Down: %d\n", gb_hyper_down))
cat(sprintf("    - Hypermethylated + Up: %d\n", gb_hyper_up))
cat(sprintf("    - Hypomethylated + Down: %d\n", gb_hypo_down))
cat(sprintf("    - Hypomethylated + Up: %d\n", gb_hypo_up))

cat("\nOutput files:\n")
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", output_wide))
cat(sprintf("  - %s\n", pdf_file))
cat(sprintf("  - %s\n", top_genes_file))

cat("\n========================================\n")
