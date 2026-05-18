#!/usr/bin/env Rscript

# ============================================================================
# CpG Island Context Analysis - R Statistics and Visualization
# ============================================================================
#
# Purpose: Analyze DMR distribution across CpG island geography and
#          calculate enrichment statistics with proper controls
#
# Input:
#   - dmr_cpg_annotation.tsv (from bash script)
#   - genome_background.tsv (from bash script)
#   - Original DMR file (for additional statistics)
#
# Output:
#   - cpg_context_enrichment.csv (Fisher's exact test results)
#   - cpg_context_summary.csv (summary statistics)
#   - Visualization plots
#
# ============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
})

cat("============================================\n")
cat("CpG Context Analysis - R Script\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP")

OUTDIR <- "results/22_cpg_context"
PLOTDIR <- file.path(OUTDIR, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load data
# ============================================================================

cat("=== Loading data ===\n")

# Load annotated DMRs
dmr_annot <- read.delim(file.path(OUTDIR, "dmr_cpg_annotation.tsv"),
                        stringsAsFactors = FALSE)
cat("  Loaded", nrow(dmr_annot), "annotated DMRs\n")

# Load genome background
background <- read.delim(file.path(OUTDIR, "genome_background.tsv"),
                         stringsAsFactors = FALSE)
cat("  Loaded genome background\n")

# Load original DMR file for full logFC data
dmr_full <- read.csv("results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv",
                     stringsAsFactors = FALSE)
cat("  Loaded", nrow(dmr_full), "original DMRs\n\n")

# ============================================================================
# Summary statistics by CpG context
# ============================================================================

cat("=== Calculating summary statistics ===\n")

# Define context order (from most to least CpG-dense)
context_order <- c("Island", "Shore", "Shelf", "OpenSea")
dmr_annot$cpg_context <- factor(dmr_annot$cpg_context, levels = context_order)

# Count DMRs per context
dmr_counts <- dmr_annot %>%
    group_by(cpg_context) %>%
    summarise(
        n_dmrs = n(),
        mean_logFC = mean(logFC, na.rm = TRUE),
        median_logFC = median(logFC, na.rm = TRUE),
        sd_logFC = sd(logFC, na.rm = TRUE),
        n_hyper = sum(logFC > 0, na.rm = TRUE),
        n_hypo = sum(logFC < 0, na.rm = TRUE),
        pct_hyper = n_hyper / n() * 100,
        .groups = "drop"
    ) %>%
    mutate(pct_total = n_dmrs / sum(n_dmrs) * 100)

# Add genome background percentages
background_pct <- background %>%
    select(context, pct_genome) %>%
    rename(cpg_context = context, genome_pct = pct_genome)

dmr_counts <- dmr_counts %>%
    left_join(background_pct, by = "cpg_context")

print(dmr_counts)

# Save summary
write.csv(dmr_counts, file.path(OUTDIR, "cpg_context_summary.csv"),
          row.names = FALSE)
cat("\n  Saved: cpg_context_summary.csv\n\n")

# ============================================================================
# Enrichment analysis (Fisher's exact test)
# ============================================================================

cat("=== Enrichment analysis ===\n")

# For enrichment, we need to compare:
# - Proportion of DMRs in each context vs
# - Proportion of genome in each context

# Total genome size and DMR count
total_genome_bp <- sum(background$total_bp)
total_dmrs <- nrow(dmr_annot)

# Calculate enrichment for each context
enrichment_results <- data.frame()

for (ctx in context_order) {
    ctx_dmrs <- sum(dmr_annot$cpg_context == ctx, na.rm = TRUE)
    ctx_genome_bp <- background$total_bp[background$context == ctx]
    ctx_genome_pct <- ctx_genome_bp / total_genome_bp

    # Expected DMRs based on genome coverage
    expected_dmrs <- total_dmrs * ctx_genome_pct

    # Fold enrichment
    fold_enrichment <- ctx_dmrs / expected_dmrs

    # Fisher's exact test
    # Contingency table:
    #                    In Context    Not In Context
    # DMR                ctx_dmrs      total_dmrs - ctx_dmrs
    # Expected (genome)  expected      total - expected

    # Use a chi-squared approximation for large numbers
    observed <- ctx_dmrs
    expected <- expected_dmrs

    # Binomial test (more appropriate for this setting)
    binom_test <- binom.test(ctx_dmrs, total_dmrs, p = ctx_genome_pct)

    enrichment_results <- rbind(enrichment_results, data.frame(
        cpg_context = ctx,
        n_dmrs = ctx_dmrs,
        pct_dmrs = ctx_dmrs / total_dmrs * 100,
        genome_bp = ctx_genome_bp,
        pct_genome = ctx_genome_pct * 100,
        expected_dmrs = round(expected_dmrs, 1),
        fold_enrichment = round(fold_enrichment, 3),
        log2_enrichment = round(log2(fold_enrichment), 3),
        binom_pvalue = binom_test$p.value,
        binom_95CI_low = binom_test$conf.int[1] * 100,
        binom_95CI_high = binom_test$conf.int[2] * 100
    ))
}

# Adjust p-values for multiple testing
enrichment_results$padj <- p.adjust(enrichment_results$binom_pvalue, method = "BH")

# Add significance stars
enrichment_results$significance <- case_when(
    enrichment_results$padj < 0.001 ~ "***",
    enrichment_results$padj < 0.01 ~ "**",
    enrichment_results$padj < 0.05 ~ "*",
    TRUE ~ "ns"
)

print(enrichment_results)

write.csv(enrichment_results, file.path(OUTDIR, "cpg_context_enrichment.csv"),
          row.names = FALSE)
cat("\n  Saved: cpg_context_enrichment.csv\n\n")

# ============================================================================
# Statistical tests for logFC differences between contexts
# ============================================================================

cat("=== Testing logFC differences between contexts ===\n")

# Kruskal-Wallis test (non-parametric ANOVA)
kw_test <- kruskal.test(logFC ~ cpg_context, data = dmr_annot)
cat("  Kruskal-Wallis test:\n")
cat("    Chi-squared =", round(kw_test$statistic, 2), "\n")
cat("    df =", kw_test$parameter, "\n")
cat("    p-value =", format(kw_test$p.value, digits = 3, scientific = TRUE), "\n\n")

# Pairwise Wilcoxon tests
pairwise_tests <- pairwise.wilcox.test(dmr_annot$logFC, dmr_annot$cpg_context,
                                        p.adjust.method = "BH")
cat("  Pairwise Wilcoxon tests (BH-adjusted p-values):\n")
print(pairwise_tests$p.value)

# Save statistical test results
sink(file.path(OUTDIR, "statistical_tests.txt"))
cat("CpG Context Analysis - Statistical Tests\n")
cat("=========================================\n\n")
cat("Kruskal-Wallis Test (logFC ~ CpG Context):\n")
print(kw_test)
cat("\nPairwise Wilcoxon Tests (BH-adjusted):\n")
print(pairwise_tests)
sink()
cat("\n  Saved: statistical_tests.txt\n\n")

# ============================================================================
# Visualization 1: DMR distribution by CpG context
# ============================================================================

cat("=== Generating visualizations ===\n")

# Bar plot of DMR counts
p1 <- ggplot(dmr_counts, aes(x = cpg_context, y = n_dmrs, fill = cpg_context)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
    geom_text(aes(label = paste0(n_dmrs, "\n(", round(pct_total, 1), "%)")),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Island" = "#1B9E77", "Shore" = "#D95F02",
                                 "Shelf" = "#7570B3", "OpenSea" = "#E7298A")) +
    labs(title = "DMR Distribution by CpG Island Context",
         subtitle = paste0("Total DMRs: ", total_dmrs),
         x = "CpG Context", y = "Number of DMRs") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(PLOTDIR, "dmr_count_by_cpg_context.png"), p1,
       width = 8, height = 6, dpi = 300)
cat("  Saved: dmr_count_by_cpg_context.png\n")

# ============================================================================
# Visualization 2: Enrichment plot
# ============================================================================

# Enrichment bar plot with error bars
enrichment_results$cpg_context <- factor(enrichment_results$cpg_context,
                                         levels = context_order)

p2 <- ggplot(enrichment_results, aes(x = cpg_context, y = log2_enrichment,
                                      fill = cpg_context)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_text(aes(label = significance,
                  y = ifelse(log2_enrichment >= 0, log2_enrichment + 0.1, log2_enrichment - 0.15)),
              size = 5) +
    scale_fill_manual(values = c("Island" = "#1B9E77", "Shore" = "#D95F02",
                                 "Shelf" = "#7570B3", "OpenSea" = "#E7298A")) +
    labs(title = "DMR Enrichment by CpG Context",
         subtitle = "Compared to genomic background",
         x = "CpG Context", y = "log2(Fold Enrichment)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(PLOTDIR, "dmr_enrichment_by_cpg_context.png"), p2,
       width = 8, height = 6, dpi = 300)
cat("  Saved: dmr_enrichment_by_cpg_context.png\n")

# ============================================================================
# Visualization 3: logFC distribution by context (boxplot)
# ============================================================================

p3 <- ggplot(dmr_annot, aes(x = cpg_context, y = logFC, fill = cpg_context)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("Island" = "#1B9E77", "Shore" = "#D95F02",
                                 "Shelf" = "#7570B3", "OpenSea" = "#E7298A")) +
    labs(title = "Methylation Change (logFC) by CpG Context",
         subtitle = paste0("Kruskal-Wallis p = ",
                          format(kw_test$p.value, digits = 3, scientific = TRUE)),
         x = "CpG Context", y = "log2(Fold Change)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(PLOTDIR, "logFC_by_cpg_context_boxplot.png"), p3,
       width = 8, height = 6, dpi = 300)
cat("  Saved: logFC_by_cpg_context_boxplot.png\n")

# ============================================================================
# Visualization 4: Hyper vs Hypomethylation by context
# ============================================================================

# Prepare data for stacked bar
hyper_hypo_data <- dmr_counts %>%
    select(cpg_context, n_hyper, n_hypo) %>%
    pivot_longer(cols = c(n_hyper, n_hypo),
                 names_to = "direction", values_to = "count") %>%
    mutate(direction = ifelse(direction == "n_hyper", "Hypermethylated", "Hypomethylated"))

p4 <- ggplot(hyper_hypo_data, aes(x = cpg_context, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Hypermethylated" = "#D73027",
                                 "Hypomethylated" = "#4575B4"),
                      name = "Direction") +
    labs(title = "DMR Direction by CpG Context",
         x = "CpG Context", y = "Number of DMRs") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")

ggsave(file.path(PLOTDIR, "dmr_direction_by_cpg_context.png"), p4,
       width = 9, height = 6, dpi = 300)
cat("  Saved: dmr_direction_by_cpg_context.png\n")

# ============================================================================
# Visualization 5: Comparison of observed vs expected DMRs
# ============================================================================

comparison_data <- enrichment_results %>%
    select(cpg_context, n_dmrs, expected_dmrs) %>%
    pivot_longer(cols = c(n_dmrs, expected_dmrs),
                 names_to = "type", values_to = "count") %>%
    mutate(type = ifelse(type == "n_dmrs", "Observed", "Expected"))

p5 <- ggplot(comparison_data, aes(x = cpg_context, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Observed" = "#66C2A5", "Expected" = "#FC8D62"),
                      name = "") +
    labs(title = "Observed vs Expected DMRs by CpG Context",
         subtitle = "Expected based on genomic coverage",
         x = "CpG Context", y = "Number of DMRs") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "top")

ggsave(file.path(PLOTDIR, "observed_vs_expected_dmrs.png"), p5,
       width = 9, height = 6, dpi = 300)
cat("  Saved: observed_vs_expected_dmrs.png\n")

# ============================================================================
# Visualization 6: Combined summary figure
# ============================================================================

# Create a multi-panel figure
library(gridExtra)

p_combined <- grid.arrange(
    p1 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    p2 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    p3 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    p4 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    ncol = 2, nrow = 2,
    top = grid::textGrob("CpG Island Context Analysis Summary",
                         gp = grid::gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(PLOTDIR, "cpg_context_summary_figure.png"), p_combined,
       width = 14, height = 12, dpi = 300)
cat("  Saved: cpg_context_summary_figure.png\n")

# ============================================================================
# Summary output
# ============================================================================

cat("\n============================================\n")
cat("Analysis Complete\n")
cat("============================================\n\n")

cat("Key Findings:\n")
cat("-------------\n")

# Most enriched context
max_enrich <- enrichment_results[which.max(enrichment_results$fold_enrichment), ]
cat("  Most enriched context:", as.character(max_enrich$cpg_context), "\n")
cat("    Fold enrichment:", max_enrich$fold_enrichment, "\n")
cat("    Adjusted p-value:", format(max_enrich$padj, digits = 3, scientific = TRUE), "\n\n")

# Context with highest hypermethylation
max_hyper <- dmr_counts[which.max(dmr_counts$pct_hyper), ]
cat("  Highest % hypermethylated:", as.character(max_hyper$cpg_context), "\n")
cat("    Percentage:", round(max_hyper$pct_hyper, 1), "%\n\n")

# Overall pattern
cat("  Overall logFC by context:\n")
for (i in 1:nrow(dmr_counts)) {
    cat("    ", as.character(dmr_counts$cpg_context[i]), ": median logFC =",
        round(dmr_counts$median_logFC[i], 3), "\n")
}

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
