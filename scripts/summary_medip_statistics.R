#!/usr/bin/env Rscript
################################################################################
# MeDIP Analysis Summary Statistics
# Extracts quantitative results from existing meDIP analysis outputs
#
# Author: Generated for SRF_Eva project
# Date: 2024-12
#
# Usage: Rscript summary_medip_statistics.R
# Output: results/summary/medip_summary_statistics.csv
#         results/summary/medip_top_dmrs_by_category.csv
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

cat("========================================\n")
cat("MeDIP Summary Statistics Generator\n")
cat("========================================\n\n")

# Define paths
# BASE_DIR: The root directory for the MeDIP analysis
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
RESULTS_DIR <- file.path(BASE_DIR, "results")
# OUTPUT_DIR: Where the summary CSV files will be saved
OUTPUT_DIR <- file.path(RESULTS_DIR, "summary")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1. Load DMR Data
################################################################################

cat("1. Loading DMR data...\n")

# All significant DMRs (Differentially Methylated Regions)
# We expect a CSV file containing DMRs with FDR < 0.05
dmr_file <- file.path(RESULTS_DIR, "07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv")
if (!file.exists(dmr_file)) {
    stop("DMR file not found: ", dmr_file)
}
dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)
cat(sprintf("   Loaded %d DMRs from FDR<0.05 file\n", nrow(dmrs)))

# Annotated DMRs
annotated_file <- file.path(RESULTS_DIR, "08_annotation/TES_vs_GFP_annotated.csv")
if (file.exists(annotated_file)) {
    annotated_dmrs <- read.csv(annotated_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded %d annotated DMRs\n", nrow(annotated_dmrs)))
} else {
    annotated_dmrs <- NULL
    cat("   WARNING: Annotated DMR file not found\n")
}

# Promoter DMRs
promoter_file <- file.path(RESULTS_DIR, "08_annotation/TES_vs_GFP_promoter_DMRs.csv")
if (file.exists(promoter_file)) {
    promoter_dmrs <- read.csv(promoter_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded %d promoter DMRs\n", nrow(promoter_dmrs)))
} else {
    promoter_dmrs <- NULL
    cat("   WARNING: Promoter DMR file not found\n")
}

################################################################################
# 2. Load Integration Data
################################################################################

cat("\n2. Loading integration data...\n")

# Methylation-expression integration
integration_file <- file.path(RESULTS_DIR, "16_advanced_visualization/integrated_medip_rnaseq_data.csv")
if (file.exists(integration_file)) {
    integration_data <- read.csv(integration_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded %d genes from integration file\n", nrow(integration_data)))
} else {
    integration_data <- NULL
    cat("   WARNING: Integration file not found\n")
}

# Gene set summary
geneset_file <- file.path(RESULTS_DIR, "12_gene_sets/gene_set_summary.csv")
if (file.exists(geneset_file)) {
    geneset_summary <- read.csv(geneset_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded gene set summary\n"))
} else {
    geneset_summary <- NULL
    cat("   WARNING: Gene set summary not found\n")
}

################################################################################
# 3. Calculate Summary Statistics
################################################################################

cat("\n3. Calculating summary statistics...\n")

# Initialize an empty data frame to store key-value pairs of statistics
summary_stats <- data.frame(
    Category = character(),
    Metric = character(),
    Value = character(),
    stringsAsFactors = FALSE
)

# Helper function to easily add a new statistic row to the dataframe
# category: Grouping for the stat (e.g., "DMR_Overview")
# metric: Name of the stat (e.g., "Total_DMRs")
# value: The calculated value
add_stat <- function(category, metric, value) {
    summary_stats <<- rbind(summary_stats, data.frame(
        Category = category,
        Metric = metric,
        Value = as.character(value),
        stringsAsFactors = FALSE
    ))
}

# --- DMR Statistics ---
cat("   Calculating DMR statistics...\n")

add_stat("DMR_Overview", "Total_DMRs_FDR05", nrow(dmrs))
add_stat("DMR_Overview", "Total_windows_analyzed", "6,176,550")

# Directionality based on logFC (Log Fold Change)
# logFC > 0 implies higher methylation in TES compared to GFP (Hypermethylated)
# logFC < 0 implies lower methylation in TES compared to GFP (Hypomethylated)
hyper_count <- sum(dmrs$logFC > 0, na.rm = TRUE)
hypo_count <- sum(dmrs$logFC < 0, na.rm = TRUE)
add_stat("DMR_Direction", "Hypermethylated_DMRs", hyper_count)
add_stat("DMR_Direction", "Hypomethylated_DMRs", hypo_count)
add_stat("DMR_Direction", "Percent_Hypermethylated", sprintf("%.1f%%", 100 * hyper_count / nrow(dmrs)))
add_stat("DMR_Direction", "Percent_Hypomethylated", sprintf("%.1f%%", 100 * hypo_count / nrow(dmrs)))
add_stat("DMR_Direction", "Hyper_to_Hypo_Ratio", sprintf("%.1f:1", hyper_count / hypo_count))

# Fold change statistics
add_stat("DMR_Effect_Size", "Mean_logFC", sprintf("%.3f", mean(dmrs$logFC, na.rm = TRUE)))
add_stat("DMR_Effect_Size", "Median_logFC", sprintf("%.3f", median(dmrs$logFC, na.rm = TRUE)))
add_stat("DMR_Effect_Size", "Max_logFC_hyper", sprintf("%.3f", max(dmrs$logFC, na.rm = TRUE)))
add_stat("DMR_Effect_Size", "Min_logFC_hypo", sprintf("%.3f", min(dmrs$logFC, na.rm = TRUE)))
add_stat("DMR_Effect_Size", "Max_fold_change", sprintf("%.1f", max(dmrs$fold_change, na.rm = TRUE)))

# With FC > 2 threshold
fc2_count <- sum(abs(dmrs$logFC) > 1, na.rm = TRUE)
add_stat("DMR_Stringent", "DMRs_with_FC_gt_2", fc2_count)
add_stat("DMR_Stringent", "Percent_FC_gt_2", sprintf("%.1f%%", 100 * fc2_count / nrow(dmrs)))

# Statistical significance breakdown
# Counting how many DMRs pass stricter significance thresholds
add_stat("DMR_Significance", "DMRs_FDR_lt_0.01", sum(dmrs$FDR < 0.01, na.rm = TRUE))
add_stat("DMR_Significance", "DMRs_FDR_lt_0.001", sum(dmrs$FDR < 0.001, na.rm = TRUE))
add_stat("DMR_Significance", "DMRs_Pvalue_lt_1e-4", sum(dmrs$PValue < 1e-4, na.rm = TRUE))

# CpG content
add_stat("DMR_CpG_Content", "Mean_CpG_per_window", sprintf("%.1f", mean(dmrs$CpG_count, na.rm = TRUE)))
add_stat("DMR_CpG_Content", "Median_CpG_per_window", sprintf("%.1f", median(dmrs$CpG_count, na.rm = TRUE)))

# --- Genomic Distribution (from annotated DMRs) ---
if (!is.null(annotated_dmrs)) {
    cat("   Calculating genomic distribution...\n")

    # Parse annotation categories
    # Simplifies complex annotation strings from HOMER/ChIPseeker into broad categories
    # e.g., "Intron (NM_001...)" becomes "Intron"
    annotated_dmrs$annotation_simple <- sapply(annotated_dmrs$annotation, function(x) {
        if (grepl("Promoter", x)) {
            return("Promoter")
        }
        if (grepl("Exon", x)) {
            return("Exon")
        }
        if (grepl("Intron", x)) {
            return("Intron")
        }
        if (grepl("3' UTR", x)) {
            return("3_UTR")
        }
        if (grepl("5' UTR", x)) {
            return("5_UTR")
        }
        if (grepl("Downstream", x)) {
            return("Downstream")
        }
        if (grepl("Distal Intergenic", x)) {
            return("Intergenic")
        }
        return("Other")
    })

    anno_counts <- table(annotated_dmrs$annotation_simple)
    for (anno in names(anno_counts)) {
        add_stat("Genomic_Distribution", paste0("DMRs_", anno), anno_counts[anno])
        add_stat(
            "Genomic_Distribution", paste0("Percent_", anno),
            sprintf("%.1f%%", 100 * anno_counts[anno] / nrow(annotated_dmrs))
        )
    }

    # Unique genes affected
    unique_genes <- length(unique(na.omit(annotated_dmrs$SYMBOL)))
    add_stat("Gene_Coverage", "Unique_genes_with_DMRs", unique_genes)
}

# --- Promoter DMRs ---
if (!is.null(promoter_dmrs)) {
    cat("   Calculating promoter DMR statistics...\n")

    add_stat("Promoter_DMRs", "Total_promoter_DMRs", nrow(promoter_dmrs))
    add_stat(
        "Promoter_DMRs", "Percent_of_all_DMRs",
        sprintf("%.1f%%", 100 * nrow(promoter_dmrs) / nrow(dmrs))
    )

    # Direction at promoters
    promo_hyper <- sum(promoter_dmrs$logFC > 0, na.rm = TRUE)
    promo_hypo <- sum(promoter_dmrs$logFC < 0, na.rm = TRUE)
    add_stat("Promoter_DMRs", "Promoter_hypermethylated", promo_hyper)
    add_stat("Promoter_DMRs", "Promoter_hypomethylated", promo_hypo)

    # Unique genes with promoter DMRs
    promo_genes <- length(unique(na.omit(promoter_dmrs$SYMBOL)))
    add_stat("Promoter_DMRs", "Unique_genes_promoter_DMRs", promo_genes)
}

# --- Integration Statistics ---
if (!is.null(integration_data)) {
    cat("   Calculating integration statistics...\n")

    add_stat("Integration", "Total_genes_analyzed", nrow(integration_data))

    # DEG counts
    if ("padj" %in% colnames(integration_data) && "log2FoldChange" %in% colnames(integration_data)) {
        up_degs <- sum(integration_data$padj < 0.05 & integration_data$log2FoldChange > 1, na.rm = TRUE)
        down_degs <- sum(integration_data$padj < 0.05 & integration_data$log2FoldChange < -1, na.rm = TRUE)
        add_stat("Integration", "Upregulated_DEGs", up_degs)
        add_stat("Integration", "Downregulated_DEGs", down_degs)
        add_stat("Integration", "Total_DEGs", up_degs + down_degs)
    }

    # Methylation signal
    if ("meDIP_GFP_mean" %in% colnames(integration_data)) {
        add_stat(
            "Integration", "Mean_meDIP_signal_GFP",
            sprintf("%.2f", mean(integration_data$meDIP_GFP_mean, na.rm = TRUE))
        )
        add_stat(
            "Integration", "Mean_meDIP_signal_TES",
            sprintf("%.2f", mean(integration_data$meDIP_TES_mean, na.rm = TRUE))
        )
    }

    # Correlation analysis
    # We check if migration/expression correlation exists (e.g. does low methylation = high expression?)
    if ("meDIP_log2FC" %in% colnames(integration_data) && "log2FoldChange" %in% colnames(integration_data)) {
        valid_idx <- !is.na(integration_data$meDIP_log2FC) & !is.na(integration_data$log2FoldChange)
        # Require a minimum number of points to calculate correlation
        if (sum(valid_idx) > 10) {
            pearson_cor <- cor(integration_data$meDIP_log2FC[valid_idx],
                integration_data$log2FoldChange[valid_idx],
                method = "pearson"
            )
            spearman_cor <- cor(integration_data$meDIP_log2FC[valid_idx],
                integration_data$log2FoldChange[valid_idx],
                method = "spearman"
            )
            pearson_test <- cor.test(integration_data$meDIP_log2FC[valid_idx],
                integration_data$log2FoldChange[valid_idx],
                method = "pearson"
            )
            spearman_test <- cor.test(integration_data$meDIP_log2FC[valid_idx],
                integration_data$log2FoldChange[valid_idx],
                method = "spearman"
            )

            add_stat("Correlation", "Pearson_r", sprintf("%.4f", pearson_cor))
            add_stat("Correlation", "Pearson_pvalue", sprintf("%.2e", pearson_test$p.value))
            add_stat("Correlation", "Spearman_rho", sprintf("%.4f", spearman_cor))
            add_stat("Correlation", "Spearman_pvalue", sprintf("%.2e", spearman_test$p.value))
        }
    }
}

# --- Gene Set Summary ---
if (!is.null(geneset_summary)) {
    cat("   Adding gene set statistics...\n")
    for (i in 1:nrow(geneset_summary)) {
        add_stat(
            "Gene_Sets",
            paste0("N_", gsub(" ", "_", geneset_summary$Category[i])),
            geneset_summary$Number_of_genes[i]
        )
    }
}

################################################################################
# 4. Generate Top DMRs by Category
################################################################################

cat("\n4. Generating top DMRs by category...\n")

# List to collect different categories of top DMRs
top_dmrs_list <- list()

# Top hypermethylated DMRs
if (!is.null(annotated_dmrs)) {
    top_hyper <- annotated_dmrs %>%
        filter(logFC > 0) %>%
        arrange(desc(logFC)) %>%
        head(50) %>%
        select(
            chr = seqnames, start, end, logFC, FDR, fold_change,
            annotation, SYMBOL, GENENAME, direction
        ) %>%
        mutate(Category = "Top_Hypermethylated")

    top_hypo <- annotated_dmrs %>%
        filter(logFC < 0) %>%
        arrange(logFC) %>%
        head(50) %>%
        select(
            chr = seqnames, start, end, logFC, FDR, fold_change,
            annotation, SYMBOL, GENENAME, direction
        ) %>%
        mutate(Category = "Top_Hypomethylated")

    # Top promoter DMRs
    top_promoter <- annotated_dmrs %>%
        filter(grepl("Promoter", annotation)) %>%
        arrange(desc(abs(logFC))) %>%
        head(50) %>%
        select(
            chr = seqnames, start, end, logFC, FDR, fold_change,
            annotation, SYMBOL, GENENAME, direction
        ) %>%
        mutate(Category = "Top_Promoter_DMRs")

    # Most significant DMRs
    top_significant <- annotated_dmrs %>%
        arrange(FDR) %>%
        head(50) %>%
        select(
            chr = seqnames, start, end, logFC, FDR, fold_change,
            annotation, SYMBOL, GENENAME, direction
        ) %>%
        mutate(Category = "Most_Significant")

    top_dmrs <- bind_rows(top_hyper, top_hypo, top_promoter, top_significant)

    cat(sprintf("   Generated top DMR lists (%d entries)\n", nrow(top_dmrs)))
}

################################################################################
# 5. Save Results
################################################################################

cat("\n5. Saving results...\n")

# Save summary statistics
summary_file <- file.path(OUTPUT_DIR, "medip_summary_statistics.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(sprintf("   Saved: %s\n", summary_file))

# Save top DMRs
if (exists("top_dmrs")) {
    top_dmrs_file <- file.path(OUTPUT_DIR, "medip_top_dmrs_by_category.csv")
    write.csv(top_dmrs, top_dmrs_file, row.names = FALSE)
    cat(sprintf("   Saved: %s\n", top_dmrs_file))
}

# Create a wide-format summary for quick viewing
# Flattens the Category+Metric into a single column name for easier reading in Excel
summary_wide <- summary_stats %>%
    unite(Full_Metric, Category, Metric, sep = "__") %>%
    pivot_wider(names_from = Full_Metric, values_from = Value)

summary_wide_file <- file.path(OUTPUT_DIR, "medip_summary_wide.csv")
write.csv(summary_wide, summary_wide_file, row.names = FALSE)
cat(sprintf("   Saved: %s\n", summary_wide_file))

################################################################################
# 6. Print Summary to Console
################################################################################

cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n\n")

# Print key statistics
key_stats <- summary_stats %>%
    filter(Category %in% c("DMR_Overview", "DMR_Direction", "Promoter_DMRs", "Correlation"))

for (cat in unique(key_stats$Category)) {
    cat(sprintf("\n%s:\n", cat))
    cat_stats <- key_stats %>% filter(Category == cat)
    for (i in 1:nrow(cat_stats)) {
        cat(sprintf("  %s: %s\n", cat_stats$Metric[i], cat_stats$Value[i]))
    }
}

cat("\n========================================\n")
cat("Summary generation complete!\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("========================================\n")
