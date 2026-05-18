#!/usr/bin/env Rscript
################################################################################
# Generate Markdown Summary Report
# Creates a formatted markdown file with key statistics from summary CSVs
#
# Author: Generated for SRF_Eva project
# Date: 2024-12
#
# Usage: Rscript generate_markdown_summary.R
# Output: results_summary/KEY_STATISTICS_SUMMARY.md
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

cat("========================================\n")
cat("Generating Markdown Summary Report\n")
cat("========================================\n\n")

# Define paths
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
MEDIP_SUMMARY <- file.path(BASE_DIR, "meDIP/results/summary/medip_summary_statistics.csv")
INTEGRATION_SUMMARY <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/results/summary/integration_summary_statistics.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "results_summary")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "KEY_STATISTICS_SUMMARY.md")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Load Summary Data
################################################################################

cat("Loading summary files...\n")

if (!file.exists(MEDIP_SUMMARY)) {
    stop("meDIP summary not found. Run summary_medip_statistics.R first.")
}
if (!file.exists(INTEGRATION_SUMMARY)) {
    stop("Integration summary not found. Run summary_integration_statistics.R first.")
}

medip_stats <- read.csv(MEDIP_SUMMARY, stringsAsFactors = FALSE)
integration_stats <- read.csv(INTEGRATION_SUMMARY, stringsAsFactors = FALSE)

cat(sprintf("  Loaded %d meDIP statistics\n", nrow(medip_stats)))
cat(sprintf("  Loaded %d integration statistics\n", nrow(integration_stats)))

################################################################################
# Helper Function to Extract Values
################################################################################

get_stat <- function(df, category, metric) {
    val <- df %>%
        filter(Category == category, Metric == metric) %>%
        pull(Value)
    if (length(val) == 0) return("N/A")
    return(val[1])
}

# Shorthand functions
m <- function(cat, met) get_stat(medip_stats, cat, met)
i <- function(cat, met) get_stat(integration_stats, cat, met)

################################################################################
# Generate Markdown Content
################################################################################

cat("Generating markdown content...\n")

md_content <- paste0(
'# Multi-Omics Analysis: Key Statistics Summary

**Project:** TES vs GFP Comparison in SNB19 Glioblastoma Cells
**Generated:** ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '
**Data Sources:** Cut&Tag (Binding), meDIP-seq (Methylation), RNA-seq (Expression)

---

## 1. Differential Methylation Analysis (meDIP-seq)

### Overview
| Metric | Value |
|--------|-------|
| Total windows analyzed | ', m("DMR_Overview", "Total_windows_analyzed"), ' |
| Significant DMRs (FDR < 0.05) | ', m("DMR_Overview", "Total_DMRs_FDR05"), ' |
| Hypermethylated DMRs | ', m("DMR_Direction", "Hypermethylated_DMRs"), ' (', m("DMR_Direction", "Percent_Hypermethylated"), ') |
| Hypomethylated DMRs | ', m("DMR_Direction", "Hypomethylated_DMRs"), ' (', m("DMR_Direction", "Percent_Hypomethylated"), ') |
| Hyper:Hypo ratio | ', m("DMR_Direction", "Hyper_to_Hypo_Ratio"), ' |

### Effect Sizes
| Metric | Value |
|--------|-------|
| Mean log2FC | ', m("DMR_Effect_Size", "Mean_logFC"), ' |
| Median log2FC | ', m("DMR_Effect_Size", "Median_logFC"), ' |
| Max log2FC (hypermethylated) | ', m("DMR_Effect_Size", "Max_logFC_hyper"), ' |
| Min log2FC (hypomethylated) | ', m("DMR_Effect_Size", "Min_logFC_hypo"), ' |
| Max fold change | ', m("DMR_Effect_Size", "Max_fold_change"), 'x |

### Genomic Distribution
| Location | Count | Percentage |
|----------|-------|------------|
| Promoter | ', m("Genomic_Distribution", "DMRs_Promoter"), ' | ', m("Genomic_Distribution", "Percent_Promoter"), ' |
| Intron | ', m("Genomic_Distribution", "DMRs_Intron"), ' | ', m("Genomic_Distribution", "Percent_Intron"), ' |
| Intergenic | ', m("Genomic_Distribution", "DMRs_Intergenic"), ' | ', m("Genomic_Distribution", "Percent_Intergenic"), ' |
| Exon | ', m("Genomic_Distribution", "DMRs_Exon"), ' | ', m("Genomic_Distribution", "Percent_Exon"), ' |
| 3\' UTR | ', m("Genomic_Distribution", "DMRs_3_UTR"), ' | ', m("Genomic_Distribution", "Percent_3_UTR"), ' |
| 5\' UTR | ', m("Genomic_Distribution", "DMRs_5_UTR"), ' | ', m("Genomic_Distribution", "Percent_5_UTR"), ' |

### Promoter DMRs
| Metric | Value |
|--------|-------|
| Total promoter DMRs | ', m("Promoter_DMRs", "Total_promoter_DMRs"), ' |
| Promoter hypermethylated | ', m("Promoter_DMRs", "Promoter_hypermethylated"), ' |
| Promoter hypomethylated | ', m("Promoter_DMRs", "Promoter_hypomethylated"), ' |
| Unique genes with promoter DMRs | ', m("Promoter_DMRs", "Unique_genes_promoter_DMRs"), ' |

### CpG Content
| Metric | Value |
|--------|-------|
| Mean CpGs per window | ', m("DMR_CpG_Content", "Mean_CpG_per_window"), ' |
| Median CpGs per window | ', m("DMR_CpG_Content", "Median_CpG_per_window"), ' |

---

## 2. Differential Expression Analysis (RNA-seq)

### Overview
| Metric | Value |
|--------|-------|
| Total genes analyzed | ', i("Overview", "Total_genes_analyzed"), ' |
| Upregulated DEGs | ', i("Expression_Status", "N_upregulated"), ' (', i("Expression_Status", "Pct_upregulated"), ') |
| Downregulated DEGs | ', i("Expression_Status", "N_downregulated"), ' (', i("Expression_Status", "Pct_downregulated"), ') |
| Total DEGs | ', i("DEG_Summary", "Total_DEGs"), ' |
| Unchanged genes | ', i("Expression_Status", "N_unchanged"), ' (', i("Expression_Status", "Pct_unchanged"), ') |
| Up:Down ratio | ', i("Direction_Summary", "Up_to_Down_ratio"), ' |

---

## 3. Transcription Factor Binding (Cut&Tag)

### Binding Status (All Genes)
| Category | Count | Percentage |
|----------|-------|------------|
| Neither bound | ', i("Binding_Status", "N_Neither_bound"), ' | ', i("Binding_Status", "Pct_Neither_bound"), ' |
| TEAD1 only bound | ', i("Binding_Status", "N_TEAD1_only_bound"), ' | ', i("Binding_Status", "Pct_TEAD1_only_bound"), ' |
| TES only bound | ', i("Binding_Status", "N_TES_only_bound"), ' | ', i("Binding_Status", "Pct_TES_only_bound"), ' |
| TES + TEAD1 bound | ', i("Binding_Status", "N_TES_TEAD1_bound"), ' | ', i("Binding_Status", "Pct_TES_TEAD1_bound"), ' |

### Direct Targets (DEGs with Binding)
| Metric | Value |
|--------|-------|
| DEGs with TES binding | ', i("Direct_Targets", "DEGs_with_TES_binding"), ' (', i("Direct_Targets", "Pct_DEGs_TES_bound"), ') |
| DEGs with TEAD1 binding | ', i("Direct_Targets", "DEGs_with_TEAD1_binding"), ' |
| DEGs with both TES + TEAD1 | ', i("Direct_Targets", "DEGs_with_TES_and_TEAD1_binding"), ' |
| DEGs with any binding | ', i("Direct_Targets", "DEGs_with_any_binding"), ' (', i("Direct_Targets", "Pct_DEGs_any_bound"), ') |
| Indirect DEGs (no binding) | ', i("Direct_Targets", "DEGs_without_binding_indirect"), ' (', i("Direct_Targets", "Pct_DEGs_indirect"), ') |

### Direction by Binding Status
| Category | Upregulated | Downregulated |
|----------|-------------|---------------|
| TES-bound | ', i("Direction_by_Binding", "TES_bound_upregulated"), ' | ', i("Direction_by_Binding", "TES_bound_downregulated"), ' |
| TEAD1-bound | ', i("Direction_by_Binding", "TEAD1_bound_upregulated"), ' | ', i("Direction_by_Binding", "TEAD1_bound_downregulated"), ' |

---

## 4. Combined Classification (Binding + Expression)

| Category | Count |
|----------|-------|
| TES+TEAD1 bound + downregulated | ', i("Combined_Classification", "TES_TEAD1_bound_downregulated"), ' |
| TES+TEAD1 bound + upregulated | ', i("Combined_Classification", "TES_TEAD1_bound_upregulated"), ' |
| TES+TEAD1 bound + unchanged | ', i("Combined_Classification", "TES_TEAD1_bound_unchanged"), ' |
| TES only bound + downregulated | ', i("Combined_Classification", "TES_only_bound_downregulated"), ' |
| TES only bound + upregulated | ', i("Combined_Classification", "TES_only_bound_upregulated"), ' |
| TEAD1 only bound + downregulated | ', i("Combined_Classification", "TEAD1_only_bound_downregulated"), ' |
| TEAD1 only bound + upregulated | ', i("Combined_Classification", "TEAD1_only_bound_upregulated"), ' |
| Neither bound + downregulated | ', i("Combined_Classification", "Neither_bound_downregulated"), ' |
| Neither bound + upregulated | ', i("Combined_Classification", "Neither_bound_upregulated"), ' |

---

## 5. Methylation-Expression Correlation

### Global Correlation
| Metric | Value |
|--------|-------|
| Pearson r | ', i("Meth_Expr_Correlation", "Pearson_r"), ' |
| Pearson p-value | ', i("Meth_Expr_Correlation", "Pearson_pvalue"), ' |
| Spearman rho | ', i("Meth_Expr_Correlation", "Spearman_rho"), ' |
| Spearman p-value | ', i("Meth_Expr_Correlation", "Spearman_pvalue"), ' |

### Methylation Change by Expression Category
| Category | Mean meDIP delta (TES - GFP) |
|----------|------------------------------|
| Upregulated DEGs | ', i("Meth_by_Expression", "Mean_meDIP_delta_upregulated"), ' |
| Downregulated DEGs | ', i("Meth_by_Expression", "Mean_meDIP_delta_downregulated"), ' |
| Unchanged genes | ', i("Meth_by_Expression", "Mean_meDIP_delta_unchanged"), ' |
| Wilcoxon Up vs Down p-value | ', i("Meth_by_Expression", "Wilcoxon_Up_vs_Down_pvalue"), ' |

---

## 6. DMR-DEG Enrichment Analysis

### Contingency Table
|  | With DMRs | Without DMRs |
|--|-----------|--------------|
| **DEGs** | ', i("DMR_DEG_Enrichment", "DEGs_with_DMRs"), ' | ', i("DMR_DEG_Enrichment", "DEGs_without_DMRs"), ' |
| **Non-DEGs** | ', i("DMR_DEG_Enrichment", "NonDEGs_with_DMRs"), ' | ', i("DMR_DEG_Enrichment", "NonDEGs_without_DMRs"), ' |

### Enrichment Statistics
| Metric | Value |
|--------|-------|
| % DEGs with DMRs | ', i("DMR_DEG_Enrichment", "Pct_DEGs_with_DMRs"), ' |
| % Non-DEGs with DMRs | ', i("DMR_DEG_Enrichment", "Pct_NonDEGs_with_DMRs"), ' |
| Fisher\'s Odds Ratio | ', i("DMR_DEG_Enrichment", "Fisher_odds_ratio"), ' |
| Fisher\'s p-value | ', i("DMR_DEG_Enrichment", "Fisher_pvalue"), ' |
| 95% CI (lower) | ', i("DMR_DEG_Enrichment", "Fisher_95CI_low"), ' |
| 95% CI (upper) | ', i("DMR_DEG_Enrichment", "Fisher_95CI_high"), ' |

---

## 7. Gene Set Summary

| Category | Number of Genes |
|----------|-----------------|
| All genes | ', m("Gene_Sets", "N_All_genes"), ' |
| Highly expressed (top 25%) | ', m("Gene_Sets", "N_Highly_expressed_(top_25%)"), ' |
| Lowly expressed (bottom 25%) | ', m("Gene_Sets", "N_Lowly_expressed_(bottom_25%)"), ' |
| Upregulated DEGs | ', m("Gene_Sets", "N_Upregulated_DEGs"), ' |
| Downregulated DEGs | ', m("Gene_Sets", "N_Downregulated_DEGs"), ' |
| Unchanged genes | ', m("Gene_Sets", "N_Unchanged_genes"), ' |

---

## 8. Key Biological Conclusions

### TES Promotes DNA Hypermethylation
- **91% of DMRs are hypermethylated** in TES-expressing cells
- Hyper:Hypo ratio of **', m("DMR_Direction", "Hyper_to_Hypo_Ratio"), '** indicates strong methyltransferase recruitment
- Consistent with TES construct containing DNMT3A/3L methyltransferase domains

### Most DEGs Are Indirect Targets
- Only **', i("Direct_Targets", "Pct_DEGs_any_bound"), '** of DEGs have direct TF binding
- **', i("Direct_Targets", "Pct_DEGs_indirect"), '** of DEGs are indirect (downstream cascade effects)
- TES binding preferentially associated with **downregulation** (', i("Direction_by_Binding", "TES_bound_downregulated"), ' down vs ', i("Direction_by_Binding", "TES_bound_upregulated"), ' up)

### Weak But Significant Methylation-Expression Coupling
- Pearson correlation r = ', i("Meth_Expr_Correlation", "Pearson_r"), ' (p = ', i("Meth_Expr_Correlation", "Pearson_pvalue"), ')
- Downregulated genes show **negative methylation change** (mean = ', i("Meth_by_Expression", "Mean_meDIP_delta_downregulated"), ')
- Suggests methylation is one of multiple regulatory mechanisms

### Modest DMR Enrichment at DEGs
- ', i("DMR_DEG_Enrichment", "Pct_DEGs_with_DMRs"), ' of DEGs have DMRs vs ', i("DMR_DEG_Enrichment", "Pct_NonDEGs_with_DMRs"), ' of non-DEGs
- Fisher\'s OR = ', i("DMR_DEG_Enrichment", "Fisher_odds_ratio"), ' (p = ', i("DMR_DEG_Enrichment", "Fisher_pvalue"), ')
- Significant but modest enrichment supports indirect methylation hypothesis

---

## Source Files

### meDIP Summary
- `meDIP/results/summary/medip_summary_statistics.csv`
- `meDIP/results/summary/medip_top_dmrs_by_category.csv`

### Integration Summary
- `SRF_Eva_integrated_analysis/results/summary/integration_summary_statistics.csv`
- `SRF_Eva_integrated_analysis/results/summary/integration_top_genes_by_category.csv`

---

*Report generated automatically from analysis outputs.*
')

################################################################################
# Write Output
################################################################################

cat("Writing markdown file...\n")

writeLines(md_content, OUTPUT_FILE)

cat(sprintf("\nSUCCESS: Markdown summary written to:\n  %s\n", OUTPUT_FILE))

# Also copy to both result directories for convenience
file.copy(OUTPUT_FILE, file.path(BASE_DIR, "meDIP/results/summary/KEY_STATISTICS_SUMMARY.md"), overwrite = TRUE)
file.copy(OUTPUT_FILE, file.path(BASE_DIR, "SRF_Eva_integrated_analysis/results/summary/KEY_STATISTICS_SUMMARY.md"), overwrite = TRUE)

cat("Copies also saved to meDIP and integration summary directories.\n")

cat("\n========================================\n")
cat("Markdown generation complete!\n")
cat("========================================\n")
