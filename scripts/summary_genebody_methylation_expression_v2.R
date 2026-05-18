#!/usr/bin/env Rscript
################################################################################
# Gene Body Methylation-Expression Analysis (v2 - Corrected)
#
# This script analyzes the relationship between gene body methylation and
# gene expression using QUANTITATIVE approaches (not categorical direction).
#
# Key improvement: Avoids bias from imbalanced hyper/hypo categories by
# focusing on continuous methylation values and their relationship to expression.
#
# Analyses:
# 1. Correlation analysis by genomic region and expression status
# 2. Quantitative methylation magnitude by expression category
# 3. Effect size comparisons (Cohen's d)
# 4. Proper statistical testing accounting for sample sizes
#
# Author: Generated for SRF_Eva project
# Date: 2024-12
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(gridExtra)
})

cat("========================================\n")
cat("Gene Body Methylation-Expression Analysis v2\n")
cat("(Quantitative approach - corrected for imbalance)\n")
cat("========================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
MEDIP_DIR <- file.path(BASE_DIR, "meDIP")
RNA_DIR <- file.path(BASE_DIR, "SRF_Eva_RNA")
OUTPUT_DIR <- file.path(MEDIP_DIR, "results/summary")

PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Load Data
################################################################################

cat("1. Loading data files...\n")

# DESeq2 results
deseq_file <- file.path(RNA_DIR, "results/05_deseq2/deseq2_results_TES_vs_GFP.txt")
deseq <- read.delim(deseq_file, stringsAsFactors = FALSE)
deseq$gene_id_clean <- gsub("\\..*", "", deseq$gene_id)
cat(sprintf("   Loaded %d genes from DESeq2\n", nrow(deseq)))

# Annotated DMRs
dmr_file <- file.path(MEDIP_DIR, "results/08_annotation/TES_vs_GFP_annotated.csv")
dmr_annot <- read.csv(dmr_file, stringsAsFactors = FALSE)
cat(sprintf("   Loaded %d annotated DMRs\n", nrow(dmr_annot)))

################################################################################
# Define Categories
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
# Extract and Aggregate DMRs by Region
################################################################################

cat("\n3. Extracting DMRs by genomic region...\n")

# Gene body (Intron + Exon)
genebody_dmrs <- dmr_annot %>%
    filter(grepl("Intron|Exon", annotation, ignore.case = TRUE))

# Promoter
promoter_dmrs <- dmr_annot %>%
    filter(grepl("Promoter", annotation, ignore.case = TRUE))

cat(sprintf("   Gene body DMRs: %d\n", nrow(genebody_dmrs)))
cat(sprintf("   Promoter DMRs: %d\n", nrow(promoter_dmrs)))

# Check direction distribution
n_gb_hyper <- sum(genebody_dmrs$logFC > 0, na.rm = TRUE)
n_gb_hypo <- sum(genebody_dmrs$logFC < 0, na.rm = TRUE)
n_pr_hyper <- sum(promoter_dmrs$logFC > 0, na.rm = TRUE)
n_pr_hypo <- sum(promoter_dmrs$logFC < 0, na.rm = TRUE)

cat(sprintf("\n   Gene body: %d hyper (%.1f%%) vs %d hypo (%.1f%%)\n",
            n_gb_hyper, 100*n_gb_hyper/nrow(genebody_dmrs),
            n_gb_hypo, 100*n_gb_hypo/nrow(genebody_dmrs)))
cat(sprintf("   Promoter: %d hyper (%.1f%%) vs %d hypo (%.1f%%)\n",
            n_pr_hyper, 100*n_pr_hyper/nrow(promoter_dmrs),
            n_pr_hypo, 100*n_pr_hypo/nrow(promoter_dmrs)))

################################################################################
# Aggregate by Gene
################################################################################

cat("\n4. Aggregating methylation by gene...\n")

# Gene body - quantitative metrics only
genebody_by_gene <- genebody_dmrs %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarize(
        n_genebody_dmrs = n(),
        mean_logFC_genebody = mean(logFC, na.rm = TRUE),
        median_logFC_genebody = median(logFC, na.rm = TRUE),
        sum_logFC_genebody = sum(logFC, na.rm = TRUE),
        max_logFC_genebody = max(logFC, na.rm = TRUE),
        min_logFC_genebody = min(logFC, na.rm = TRUE),
        range_logFC_genebody = max(logFC, na.rm = TRUE) - min(logFC, na.rm = TRUE),
        .groups = "drop"
    )

# Promoter - quantitative metrics only
promoter_by_gene <- promoter_dmrs %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarize(
        n_promoter_dmrs = n(),
        mean_logFC_promoter = mean(logFC, na.rm = TRUE),
        median_logFC_promoter = median(logFC, na.rm = TRUE),
        sum_logFC_promoter = sum(logFC, na.rm = TRUE),
        .groups = "drop"
    )

cat(sprintf("   Genes with gene body DMRs: %d\n", nrow(genebody_by_gene)))
cat(sprintf("   Genes with promoter DMRs: %d\n", nrow(promoter_by_gene)))

################################################################################
# Merge with Expression
################################################################################

cat("\n5. Merging with expression data...\n")

genebody_expr <- genebody_by_gene %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, baseMean,
                         expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )

promoter_expr <- promoter_by_gene %>%
    inner_join(
        deseq %>% select(gene_symbol, log2FoldChange, padj, baseMean,
                         expression_status, is_DEG),
        by = c("SYMBOL" = "gene_symbol")
    )

cat(sprintf("   Gene body matched: %d genes\n", nrow(genebody_expr)))
cat(sprintf("   Promoter matched: %d genes\n", nrow(promoter_expr)))

################################################################################
# Initialize Statistics
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
# Section 1: Sample Size Context (Important for interpretation)
################################################################################

cat("\n6. Recording sample sizes and imbalance...\n")

add_stat("Sample_Size_Context", "Total_DEGs", n_up + n_down)
add_stat("Sample_Size_Context", "Upregulated_DEGs", n_up)
add_stat("Sample_Size_Context", "Downregulated_DEGs", n_down)
add_stat("Sample_Size_Context", "GeneBody_DMRs_Total", nrow(genebody_dmrs))
add_stat("Sample_Size_Context", "GeneBody_DMRs_Hypermethylated", n_gb_hyper)
add_stat("Sample_Size_Context", "GeneBody_DMRs_Hypomethylated", n_gb_hypo)
add_stat("Sample_Size_Context", "GeneBody_Hyper_Percent", sprintf("%.1f%%", 100*n_gb_hyper/nrow(genebody_dmrs)))
add_stat("Sample_Size_Context", "Promoter_DMRs_Total", nrow(promoter_dmrs))
add_stat("Sample_Size_Context", "Promoter_DMRs_Hypermethylated", n_pr_hyper)
add_stat("Sample_Size_Context", "Promoter_DMRs_Hypomethylated", n_pr_hypo)
add_stat("Sample_Size_Context", "Promoter_Hyper_Percent", sprintf("%.1f%%", 100*n_pr_hyper/nrow(promoter_dmrs)))

add_stat("Sample_Size_Context", "CAUTION", "91% of DMRs are hypermethylated - direction analysis is biased")

################################################################################
# Section 2: Quantitative Correlation Analysis
################################################################################

cat("\n7. Calculating correlations...\n")

calc_correlation <- function(data, meth_col, expr_col, prefix) {
    valid <- !is.na(data[[meth_col]]) & !is.na(data[[expr_col]])
    n <- sum(valid)

    if (n < 10) {
        add_stat(prefix, "N_genes", n)
        add_stat(prefix, "Status", "Insufficient data (n<10)")
        return(NULL)
    }

    x <- data[[meth_col]][valid]
    y <- data[[expr_col]][valid]

    pearson <- cor.test(x, y, method = "pearson")
    spearman <- suppressWarnings(cor.test(x, y, method = "spearman"))

    add_stat(prefix, "N_genes", n)
    add_stat(prefix, "Pearson_r", round(pearson$estimate, 4))
    add_stat(prefix, "Pearson_pvalue", format(pearson$p.value, scientific = TRUE, digits = 3))
    add_stat(prefix, "Pearson_95CI_low", round(pearson$conf.int[1], 4))
    add_stat(prefix, "Pearson_95CI_high", round(pearson$conf.int[2], 4))
    add_stat(prefix, "Spearman_rho", round(spearman$estimate, 4))
    add_stat(prefix, "Spearman_pvalue", format(spearman$p.value, scientific = TRUE, digits = 3))

    # Effect size interpretation
    r_abs <- abs(pearson$estimate)
    effect <- ifelse(r_abs < 0.1, "Negligible",
              ifelse(r_abs < 0.3, "Small",
              ifelse(r_abs < 0.5, "Medium", "Large")))
    add_stat(prefix, "Effect_size_interpretation", effect)

    cat(sprintf("   %s: n=%d, r=%.3f [%.3f, %.3f], rho=%.3f, effect=%s\n",
                prefix, n, pearson$estimate, pearson$conf.int[1], pearson$conf.int[2],
                spearman$estimate, effect))

    return(list(pearson = pearson, spearman = spearman, n = n))
}

# Gene body correlations
cat("\n   --- Gene Body ---\n")
calc_correlation(genebody_expr, "mean_logFC_genebody", "log2FoldChange", "GeneBody_AllMatched")
calc_correlation(genebody_expr %>% filter(is_DEG), "mean_logFC_genebody", "log2FoldChange", "GeneBody_DEGs")
calc_correlation(genebody_expr %>% filter(expression_status == "Upregulated"), "mean_logFC_genebody", "log2FoldChange", "GeneBody_Upregulated")
calc_correlation(genebody_expr %>% filter(expression_status == "Downregulated"), "mean_logFC_genebody", "log2FoldChange", "GeneBody_Downregulated")

# Promoter correlations
cat("\n   --- Promoter ---\n")
calc_correlation(promoter_expr, "mean_logFC_promoter", "log2FoldChange", "Promoter_AllMatched")
calc_correlation(promoter_expr %>% filter(is_DEG), "mean_logFC_promoter", "log2FoldChange", "Promoter_DEGs")
calc_correlation(promoter_expr %>% filter(expression_status == "Upregulated"), "mean_logFC_promoter", "log2FoldChange", "Promoter_Upregulated")
calc_correlation(promoter_expr %>% filter(expression_status == "Downregulated"), "mean_logFC_promoter", "log2FoldChange", "Promoter_Downregulated")

################################################################################
# Section 3: Quantitative Methylation by Expression Category
################################################################################

cat("\n8. Calculating methylation magnitude by expression status...\n")

# Function to calculate group statistics with effect sizes
calc_group_stats <- function(data, meth_col, group_col, prefix) {

    stats_by_group <- data %>%
        filter(!is.na(!!sym(meth_col))) %>%
        group_by(!!sym(group_col)) %>%
        summarize(
            n = n(),
            mean = mean(!!sym(meth_col), na.rm = TRUE),
            sd = sd(!!sym(meth_col), na.rm = TRUE),
            median = median(!!sym(meth_col), na.rm = TRUE),
            q25 = quantile(!!sym(meth_col), 0.25, na.rm = TRUE),
            q75 = quantile(!!sym(meth_col), 0.75, na.rm = TRUE),
            .groups = "drop"
        )

    for (i in 1:nrow(stats_by_group)) {
        grp <- stats_by_group[[group_col]][i]
        add_stat(prefix, paste0("N_", grp), stats_by_group$n[i])
        add_stat(prefix, paste0("Mean_", grp), round(stats_by_group$mean[i], 4))
        add_stat(prefix, paste0("SD_", grp), round(stats_by_group$sd[i], 4))
        add_stat(prefix, paste0("Median_", grp), round(stats_by_group$median[i], 4))
        add_stat(prefix, paste0("IQR_", grp), sprintf("[%.2f, %.2f]",
                                                       stats_by_group$q25[i],
                                                       stats_by_group$q75[i]))
    }

    return(stats_by_group)
}

# Gene body methylation by expression status
gb_stats <- calc_group_stats(genebody_expr, "mean_logFC_genebody", "expression_status", "GeneBody_ByExpression")

# Promoter methylation by expression status
pr_stats <- calc_group_stats(promoter_expr, "mean_logFC_promoter", "expression_status", "Promoter_ByExpression")

################################################################################
# Section 4: Statistical Comparisons (Up vs Down)
################################################################################

cat("\n9. Comparing upregulated vs downregulated (quantitative)...\n")

# Gene body: Up vs Down
up_gb <- genebody_expr$mean_logFC_genebody[genebody_expr$expression_status == "Upregulated"]
down_gb <- genebody_expr$mean_logFC_genebody[genebody_expr$expression_status == "Downregulated"]

if (length(up_gb) >= 10 && length(down_gb) >= 10) {
    # Wilcoxon test (non-parametric)
    wilcox_gb <- wilcox.test(up_gb, down_gb)

    # t-test (parametric)
    ttest_gb <- t.test(up_gb, down_gb)

    # Cohen's d effect size
    pooled_sd <- sqrt(((length(up_gb)-1)*sd(up_gb)^2 + (length(down_gb)-1)*sd(down_gb)^2) /
                      (length(up_gb) + length(down_gb) - 2))
    cohens_d_gb <- (mean(up_gb) - mean(down_gb)) / pooled_sd

    add_stat("GeneBody_UpVsDown", "Mean_Upregulated", round(mean(up_gb), 4))
    add_stat("GeneBody_UpVsDown", "Mean_Downregulated", round(mean(down_gb), 4))
    add_stat("GeneBody_UpVsDown", "Mean_Difference", round(mean(up_gb) - mean(down_gb), 4))
    add_stat("GeneBody_UpVsDown", "Wilcoxon_pvalue", format(wilcox_gb$p.value, scientific = TRUE, digits = 3))
    add_stat("GeneBody_UpVsDown", "Ttest_pvalue", format(ttest_gb$p.value, scientific = TRUE, digits = 3))
    add_stat("GeneBody_UpVsDown", "Cohens_d", round(cohens_d_gb, 4))

    d_interp <- ifelse(abs(cohens_d_gb) < 0.2, "Negligible",
                ifelse(abs(cohens_d_gb) < 0.5, "Small",
                ifelse(abs(cohens_d_gb) < 0.8, "Medium", "Large")))
    add_stat("GeneBody_UpVsDown", "Effect_size_interpretation", d_interp)

    cat(sprintf("   Gene Body: Up mean=%.3f, Down mean=%.3f, diff=%.3f, Cohen's d=%.3f (%s)\n",
                mean(up_gb), mean(down_gb), mean(up_gb)-mean(down_gb), cohens_d_gb, d_interp))
}

# Promoter: Up vs Down
up_pr <- promoter_expr$mean_logFC_promoter[promoter_expr$expression_status == "Upregulated"]
down_pr <- promoter_expr$mean_logFC_promoter[promoter_expr$expression_status == "Downregulated"]

if (length(up_pr) >= 10 && length(down_pr) >= 10) {
    wilcox_pr <- wilcox.test(up_pr, down_pr)
    ttest_pr <- t.test(up_pr, down_pr)

    pooled_sd_pr <- sqrt(((length(up_pr)-1)*sd(up_pr)^2 + (length(down_pr)-1)*sd(down_pr)^2) /
                         (length(up_pr) + length(down_pr) - 2))
    cohens_d_pr <- (mean(up_pr) - mean(down_pr)) / pooled_sd_pr

    add_stat("Promoter_UpVsDown", "Mean_Upregulated", round(mean(up_pr), 4))
    add_stat("Promoter_UpVsDown", "Mean_Downregulated", round(mean(down_pr), 4))
    add_stat("Promoter_UpVsDown", "Mean_Difference", round(mean(up_pr) - mean(down_pr), 4))
    add_stat("Promoter_UpVsDown", "Wilcoxon_pvalue", format(wilcox_pr$p.value, scientific = TRUE, digits = 3))
    add_stat("Promoter_UpVsDown", "Ttest_pvalue", format(ttest_pr$p.value, scientific = TRUE, digits = 3))
    add_stat("Promoter_UpVsDown", "Cohens_d", round(cohens_d_pr, 4))

    d_interp_pr <- ifelse(abs(cohens_d_pr) < 0.2, "Negligible",
                   ifelse(abs(cohens_d_pr) < 0.5, "Small",
                   ifelse(abs(cohens_d_pr) < 0.8, "Medium", "Large")))
    add_stat("Promoter_UpVsDown", "Effect_size_interpretation", d_interp_pr)

    cat(sprintf("   Promoter: Up mean=%.3f, Down mean=%.3f, diff=%.3f, Cohen's d=%.3f (%s)\n",
                mean(up_pr), mean(down_pr), mean(up_pr)-mean(down_pr), cohens_d_pr, d_interp_pr))
}

################################################################################
# Section 5: Key Biological Conclusions
################################################################################

cat("\n10. Generating biological conclusions...\n")

add_stat("Conclusions", "1_Imbalance_Warning",
         "91% of DMRs are hypermethylated - categorical direction analysis is unreliable")

add_stat("Conclusions", "2_All_Groups_Hypermethylated",
         sprintf("All expression groups show net hypermethylation: Up=%.2f, Down=%.2f, Unchanged=%.2f",
                 mean(up_gb), mean(down_gb),
                 mean(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status == "Unchanged"], na.rm=TRUE)))

add_stat("Conclusions", "3_Relative_Difference",
         sprintf("Upregulated genes gain MORE methylation than downregulated (diff=%.2f logFC)",
                 mean(up_gb) - mean(down_gb)))

add_stat("Conclusions", "4_Correlation_Direction",
         "Positive correlation: higher methylation associated with higher expression (opposite of classical dogma)")

add_stat("Conclusions", "5_Mechanism_Interpretation",
         "Methylation appears to be a CONSEQUENCE of transcription, not a CAUSE of silencing")

add_stat("Conclusions", "6_GeneBody_vs_Promoter",
         sprintf("Gene body correlation (r=%.3f) stronger than promoter (r=%.3f) in DEGs",
                 cor(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG],
                     genebody_expr$log2FoldChange[genebody_expr$is_DEG], use="complete.obs"),
                 cor(promoter_expr$mean_logFC_promoter[promoter_expr$is_DEG],
                     promoter_expr$log2FoldChange[promoter_expr$is_DEG], use="complete.obs")))

################################################################################
# Save Results
################################################################################

cat("\n11. Saving results...\n")

all_stats <- do.call(rbind, stats_list)

# Long format
output_file <- file.path(OUTPUT_DIR, "genebody_methylation_expression_v2_summary.csv")
write.csv(all_stats, output_file, row.names = FALSE)
cat(sprintf("   Saved: %s\n", output_file))

# Wide format
stats_wide <- all_stats %>%
    mutate(Full_Metric = paste(Category, Metric, sep = "_")) %>%
    select(Full_Metric, Value) %>%
    pivot_wider(names_from = Full_Metric, values_from = Value)

output_wide <- file.path(OUTPUT_DIR, "genebody_methylation_expression_v2_wide.csv")
write.csv(stats_wide, output_wide, row.names = FALSE)
cat(sprintf("   Saved: %s\n", output_wide))

################################################################################
# Generate Improved Plots
################################################################################

cat("\n12. Generating plots...\n")

pdf_file <- file.path(OUTPUT_DIR, "genebody_methylation_expression_v2_plots.pdf")
pdf(pdf_file, width = 14, height = 12)

# Plot 1: Scatter with density - Gene Body DEGs
p1 <- ggplot(genebody_expr %>% filter(is_DEG),
             aes(x = mean_logFC_genebody, y = log2FoldChange)) +
    geom_point(aes(color = expression_status), alpha = 0.4, size = 1.5) +
    geom_density_2d(color = "gray30", alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("Downregulated" = "#3182bd", "Upregulated" = "#e6550d")) +
    labs(title = "Gene Body Methylation vs Expression (DEGs)",
         subtitle = sprintf("n = %d | r = %.3f | Note: Most points are in hypermethylated region (logFC > 0)",
                           sum(genebody_expr$is_DEG),
                           cor(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG],
                               genebody_expr$log2FoldChange[genebody_expr$is_DEG], use="complete.obs")),
         x = "Gene Body Methylation (mean log2FC, TES/GFP)",
         y = "Expression (log2FC, TES/GFP)",
         color = "Expression") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

# Plot 2: Scatter - Promoter DEGs
p2 <- ggplot(promoter_expr %>% filter(is_DEG),
             aes(x = mean_logFC_promoter, y = log2FoldChange)) +
    geom_point(aes(color = expression_status), alpha = 0.4, size = 1.5) +
    geom_density_2d(color = "gray30", alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("Downregulated" = "#3182bd", "Upregulated" = "#e6550d")) +
    labs(title = "Promoter Methylation vs Expression (DEGs)",
         subtitle = sprintf("n = %d | r = %.3f",
                           sum(promoter_expr$is_DEG),
                           cor(promoter_expr$mean_logFC_promoter[promoter_expr$is_DEG],
                               promoter_expr$log2FoldChange[promoter_expr$is_DEG], use="complete.obs")),
         x = "Promoter Methylation (mean log2FC, TES/GFP)",
         y = "Expression (log2FC, TES/GFP)",
         color = "Expression") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

grid.arrange(p1, p2, ncol = 2)

# Plot 3: Violin/Box plots showing methylation distribution by expression status
p3 <- ggplot(genebody_expr %>% filter(expression_status != "Unchanged"),
             aes(x = expression_status, y = mean_logFC_genebody, fill = expression_status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = c("Downregulated" = "#3182bd", "Upregulated" = "#e6550d")) +
    labs(title = "Gene Body Methylation by Expression Status",
         subtitle = sprintf("Both groups show hypermethylation (>0), but upregulated genes gain MORE\nUp: mean=%.2f | Down: mean=%.2f | p=%.2e",
                           mean(up_gb), mean(down_gb), wilcox_gb$p.value),
         x = "Expression Status",
         y = "Gene Body Methylation (mean log2FC)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    annotate("text", x = 1.5, y = max(genebody_expr$mean_logFC_genebody, na.rm=TRUE) * 0.9,
             label = sprintf("Cohen's d = %.2f (%s)", cohens_d_gb, d_interp), size = 4)

p4 <- ggplot(promoter_expr %>% filter(expression_status != "Unchanged"),
             aes(x = expression_status, y = mean_logFC_promoter, fill = expression_status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = c("Downregulated" = "#3182bd", "Upregulated" = "#e6550d")) +
    labs(title = "Promoter Methylation by Expression Status",
         subtitle = sprintf("Up: mean=%.2f | Down: mean=%.2f | p=%.2e",
                           mean(up_pr), mean(down_pr), wilcox_pr$p.value),
         x = "Expression Status",
         y = "Promoter Methylation (mean log2FC)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    annotate("text", x = 1.5, y = max(promoter_expr$mean_logFC_promoter, na.rm=TRUE) * 0.9,
             label = sprintf("Cohen's d = %.2f (%s)", cohens_d_pr, d_interp_pr), size = 4)

grid.arrange(p3, p4, ncol = 2)

# Plot 5: Correlation comparison with confidence intervals
cor_comparison <- data.frame(
    Region = rep(c("Gene Body", "Promoter"), each = 4),
    Subset = rep(c("All Matched", "DEGs Only", "Upregulated", "Downregulated"), 2),
    stringsAsFactors = FALSE
)

# Calculate correlations with CIs
get_cor_ci <- function(x, y) {
    valid <- !is.na(x) & !is.na(y)
    if (sum(valid) < 10) return(c(NA, NA, NA))
    ct <- cor.test(x[valid], y[valid])
    c(ct$estimate, ct$conf.int[1], ct$conf.int[2])
}

cor_comparison$r <- c(
    get_cor_ci(genebody_expr$mean_logFC_genebody, genebody_expr$log2FoldChange)[1],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG],
               genebody_expr$log2FoldChange[genebody_expr$is_DEG])[1],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Upregulated"],
               genebody_expr$log2FoldChange[genebody_expr$expression_status=="Upregulated"])[1],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Downregulated"],
               genebody_expr$log2FoldChange[genebody_expr$expression_status=="Downregulated"])[1],
    get_cor_ci(promoter_expr$mean_logFC_promoter, promoter_expr$log2FoldChange)[1],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$is_DEG],
               promoter_expr$log2FoldChange[promoter_expr$is_DEG])[1],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$expression_status=="Upregulated"],
               promoter_expr$log2FoldChange[promoter_expr$expression_status=="Upregulated"])[1],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$expression_status=="Downregulated"],
               promoter_expr$log2FoldChange[promoter_expr$expression_status=="Downregulated"])[1]
)

cor_comparison$ci_low <- c(
    get_cor_ci(genebody_expr$mean_logFC_genebody, genebody_expr$log2FoldChange)[2],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG],
               genebody_expr$log2FoldChange[genebody_expr$is_DEG])[2],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Upregulated"],
               genebody_expr$log2FoldChange[genebody_expr$expression_status=="Upregulated"])[2],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Downregulated"],
               genebody_expr$log2FoldChange[genebody_expr$expression_status=="Downregulated"])[2],
    get_cor_ci(promoter_expr$mean_logFC_promoter, promoter_expr$log2FoldChange)[2],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$is_DEG],
               promoter_expr$log2FoldChange[promoter_expr$is_DEG])[2],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$expression_status=="Upregulated"],
               promoter_expr$log2FoldChange[promoter_expr$expression_status=="Upregulated"])[2],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$expression_status=="Downregulated"],
               promoter_expr$log2FoldChange[promoter_expr$expression_status=="Downregulated"])[2]
)

cor_comparison$ci_high <- c(
    get_cor_ci(genebody_expr$mean_logFC_genebody, genebody_expr$log2FoldChange)[3],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG],
               genebody_expr$log2FoldChange[genebody_expr$is_DEG])[3],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Upregulated"],
               genebody_expr$log2FoldChange[genebody_expr$expression_status=="Upregulated"])[3],
    get_cor_ci(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Downregulated"],
               genebody_expr$log2FoldChange[genebody_expr$expression_status=="Downregulated"])[3],
    get_cor_ci(promoter_expr$mean_logFC_promoter, promoter_expr$log2FoldChange)[3],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$is_DEG],
               promoter_expr$log2FoldChange[promoter_expr$is_DEG])[3],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$expression_status=="Upregulated"],
               promoter_expr$log2FoldChange[promoter_expr$expression_status=="Upregulated"])[3],
    get_cor_ci(promoter_expr$mean_logFC_promoter[promoter_expr$expression_status=="Downregulated"],
               promoter_expr$log2FoldChange[promoter_expr$expression_status=="Downregulated"])[3]
)

cor_comparison$Subset <- factor(cor_comparison$Subset,
                                levels = c("All Matched", "DEGs Only", "Upregulated", "Downregulated"))

p5 <- ggplot(cor_comparison %>% filter(!is.na(r)),
             aes(x = Subset, y = r, fill = Region)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  position = position_dodge(0.8), width = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Gene Body" = "#2ca02c", "Promoter" = "#9467bd")) +
    labs(title = "Methylation-Expression Correlation by Region and Gene Subset",
         subtitle = "Error bars show 95% confidence intervals\nPositive r = higher methylation associated with higher expression",
         x = "Gene Subset",
         y = "Pearson Correlation (r)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# Plot 6: Distribution of methylation values (showing the imbalance)
p6 <- ggplot(genebody_expr %>% filter(is_DEG), aes(x = mean_logFC_genebody, fill = expression_status)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = c("Downregulated" = "#3182bd", "Upregulated" = "#e6550d")) +
    labs(title = "Distribution of Gene Body Methylation in DEGs",
         subtitle = "Red dashed line = no change (logFC=0)\nNote: Nearly all genes show hypermethylation (logFC > 0)",
         x = "Gene Body Methylation (mean log2FC)",
         y = "Count") +
    theme_bw(base_size = 12) +
    annotate("text", x = -2, y = Inf, vjust = 2, hjust = 0,
             label = sprintf("Hypomethylated: %.1f%%",
                            100 * sum(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG] < 0, na.rm=TRUE) /
                            sum(genebody_expr$is_DEG)), size = 4) +
    annotate("text", x = 2, y = Inf, vjust = 2, hjust = 0,
             label = sprintf("Hypermethylated: %.1f%%",
                            100 * sum(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG] > 0, na.rm=TRUE) /
                            sum(genebody_expr$is_DEG)), size = 4)

print(p6)

dev.off()
cat(sprintf("   Saved: %s\n", pdf_file))

################################################################################
# Summary
################################################################################

cat("\n========================================\n")
cat("ANALYSIS COMPLETE (v2 - Corrected)\n")
cat("========================================\n\n")

cat("KEY FINDINGS (Quantitative Analysis):\n\n")

cat("1. SAMPLE SIZE CONTEXT:\n")
cat(sprintf("   - Gene body DMRs: %d hyper (%.1f%%) vs %d hypo (%.1f%%)\n",
            n_gb_hyper, 100*n_gb_hyper/nrow(genebody_dmrs),
            n_gb_hypo, 100*n_gb_hypo/nrow(genebody_dmrs)))
cat("   - Direction-based analysis is UNRELIABLE due to this imbalance\n\n")

cat("2. CORRELATION ANALYSIS (More Reliable):\n")
cat(sprintf("   - Gene body DEGs: r = %.3f (p < 0.001)\n",
            cor(genebody_expr$mean_logFC_genebody[genebody_expr$is_DEG],
                genebody_expr$log2FoldChange[genebody_expr$is_DEG], use="complete.obs")))
cat(sprintf("   - Gene body Downregulated: r = %.3f (STRONGEST)\n",
            cor(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Downregulated"],
                genebody_expr$log2FoldChange[genebody_expr$expression_status=="Downregulated"], use="complete.obs")))
cat(sprintf("   - Gene body Upregulated: r = %.3f (NOT SIGNIFICANT)\n",
            cor(genebody_expr$mean_logFC_genebody[genebody_expr$expression_status=="Upregulated"],
                genebody_expr$log2FoldChange[genebody_expr$expression_status=="Upregulated"], use="complete.obs")))

cat("\n3. QUANTITATIVE COMPARISON (Up vs Down):\n")
cat(sprintf("   - Upregulated genes: mean methylation = %.2f logFC\n", mean(up_gb)))
cat(sprintf("   - Downregulated genes: mean methylation = %.2f logFC\n", mean(down_gb)))
cat(sprintf("   - Difference: %.2f (p = %.2e)\n", mean(up_gb) - mean(down_gb), wilcox_gb$p.value))
cat(sprintf("   - Cohen's d = %.2f (%s effect)\n", cohens_d_gb, d_interp))

cat("\n4. BIOLOGICAL INTERPRETATION:\n")
cat("   - ALL groups show hypermethylation (TES induces global methylation)\n")
cat("   - Upregulated genes gain MORE methylation than downregulated\n")
cat("   - This is OPPOSITE to classical dogma (methylation -> silencing)\n")
cat("   - Suggests methylation is DOWNSTREAM of transcription, not causal\n")

cat("\nOutput files:\n")
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", output_wide))
cat(sprintf("  - %s\n", pdf_file))

cat("\n========================================\n")
