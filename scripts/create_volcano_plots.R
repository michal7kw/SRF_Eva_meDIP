#!/usr/bin/env Rscript
################################################################################
# Volcano Plots for meDIP Differential Methylation Analysis
#
# Creates publication-quality volcano plots from annotated DMR data using
# EnhancedVolcano package
#
# Input: results/08_annotation/TES_vs_GFP_annotated.csv
# Output: results/volcano_plots/*.pdf
################################################################################

suppressPackageStartupMessages({
    library(EnhancedVolcano)
    library(dplyr)
    library(ggplot2)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
INPUT_FILE <- file.path(BASE_DIR, "results/08_annotation/TES_vs_GFP_annotated.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "results/volcano_plots")

# Thresholds
FDR_CUTOFF <- 0.05
FC_CUTOFF <- 2  # Linear fold change (corresponds to logFC = 1)
LOGFC_CUTOFF <- log2(FC_CUTOFF)
TOP_GENES_TO_LABEL <- 20

cat("========================================\n")
cat("meDIP Volcano Plot Generation\n")
cat("========================================\n\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", OUTPUT_DIR, "\n\n")

# Read data
cat("Reading annotated DMR data...\n")
dmr_data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
cat("  Total DMRs loaded:", nrow(dmr_data), "\n")

# Data cleaning
dmr_data <- dmr_data %>%
    mutate(
        SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", NA_character_, SYMBOL),
        # Create simplified annotation categories
        annotation_simple = case_when(
            grepl("^Promoter", annotation) ~ "Promoter",
            grepl("^Exon", annotation) ~ "Exon",
            grepl("^Intron", annotation) ~ "Intron",
            grepl("^3' UTR", annotation) ~ "3' UTR",
            grepl("^5' UTR", annotation) ~ "5' UTR",
            grepl("^Downstream", annotation) ~ "Downstream",
            grepl("^Distal Intergenic", annotation) ~ "Intergenic",
            TRUE ~ "Other"
        )
    )

# Summary statistics
cat("\n--- DMR Summary ---\n")
cat("By significance:\n")
sig_dmrs <- dmr_data %>% filter(FDR < FDR_CUTOFF)
strong_dmrs <- dmr_data %>% filter(FDR < FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)
cat("  FDR < 0.05:", nrow(sig_dmrs), "\n")
cat("  FDR < 0.05 & |FC| >= 2:", nrow(strong_dmrs), "\n")
cat("    Hypermethylated:", sum(strong_dmrs$logFC > 0), "\n")
cat("    Hypomethylated:", sum(strong_dmrs$logFC < 0), "\n")

cat("\nBy annotation:\n")
print(table(dmr_data$annotation_simple))

# Save summary statistics
summary_file <- file.path(OUTPUT_DIR, "volcano_summary.txt")
sink(summary_file)
cat("meDIP Volcano Plot Summary\n")
cat("Generated:", format(Sys.time()), "\n")
cat("========================================\n\n")
cat("Input file:", INPUT_FILE, "\n")
cat("Total DMRs:", nrow(dmr_data), "\n\n")
cat("Thresholds:\n")
cat("  FDR cutoff:", FDR_CUTOFF, "\n")
cat("  Fold change cutoff:", FC_CUTOFF, "(logFC =", LOGFC_CUTOFF, ")\n\n")
cat("DMR counts by significance:\n")
cat("  All DMRs:", nrow(dmr_data), "\n")
cat("  FDR < 0.05:", nrow(sig_dmrs), "\n")
cat("  FDR < 0.05 & |FC| >= 2:", nrow(strong_dmrs), "\n")
cat("    Hypermethylated:", sum(strong_dmrs$logFC > 0), "\n")
cat("    Hypomethylated:", sum(strong_dmrs$logFC < 0), "\n\n")
cat("DMR counts by annotation:\n")
print(table(dmr_data$annotation_simple))
sink()
cat("\nSummary saved to:", summary_file, "\n")

################################################################################
# Plot 1: Main Volcano Plot (All DMRs)
################################################################################

cat("\n--- Creating Main Volcano Plot ---\n")

# Select top genes to label (by p-value, with valid symbols)
top_genes <- dmr_data %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(PValue) %>%
    head(TOP_GENES_TO_LABEL) %>%
    pull(SYMBOL)

cat("Top", length(top_genes), "genes to label:", paste(head(top_genes, 10), collapse = ", "), "...\n")

# Create main volcano plot
p_main <- EnhancedVolcano(
    dmr_data,
    lab = dmr_data$SYMBOL,
    selectLab = top_genes,
    x = 'logFC',
    y = 'PValue',
    title = 'TES vs GFP: Differential DNA Methylation',
    subtitle = 'meDIP-seq Analysis (All DMRs)',
    caption = paste0('FDR < ', FDR_CUTOFF, ' & |FC| >= ', FC_CUTOFF),
    pCutoff = FDR_CUTOFF,
    FCcutoff = LOGFC_CUTOFF,
    pointSize = 1.5,
    labSize = 3.5,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.6,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3,
    xlim = c(min(dmr_data$logFC, na.rm = TRUE) - 0.5,
             max(dmr_data$logFC, na.rm = TRUE) + 0.5),
    ylim = c(0, max(-log10(dmr_data$PValue), na.rm = TRUE) + 1)
)

# Save main volcano
pdf(file.path(OUTPUT_DIR, "volcano_main.pdf"), width = 10, height = 8)
print(p_main)
dev.off()
cat("Saved: volcano_main.pdf\n")

################################################################################
# Plot 2: Promoter-only Volcano Plot
################################################################################

cat("\n--- Creating Promoter Volcano Plot ---\n")

promoter_dmrs <- dmr_data %>%
    filter(annotation_simple == "Promoter")

cat("Promoter DMRs:", nrow(promoter_dmrs), "\n")

# Top genes for promoter plot
top_promoter_genes <- promoter_dmrs %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(PValue) %>%
    head(TOP_GENES_TO_LABEL) %>%
    pull(SYMBOL)

cat("Top promoter genes:", paste(head(top_promoter_genes, 10), collapse = ", "), "...\n")

p_promoter <- EnhancedVolcano(
    promoter_dmrs,
    lab = promoter_dmrs$SYMBOL,
    selectLab = top_promoter_genes,
    x = 'logFC',
    y = 'PValue',
    title = 'TES vs GFP: Promoter Methylation Changes',
    subtitle = 'meDIP-seq Analysis (Promoter DMRs only)',
    caption = paste0('FDR < ', FDR_CUTOFF, ' & |FC| >= ', FC_CUTOFF),
    pCutoff = FDR_CUTOFF,
    FCcutoff = LOGFC_CUTOFF,
    pointSize = 2.0,
    labSize = 4.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.7,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3
)

# Save promoter volcano
pdf(file.path(OUTPUT_DIR, "volcano_promoters.pdf"), width = 10, height = 8)
print(p_promoter)
dev.off()
cat("Saved: volcano_promoters.pdf\n")

################################################################################
# Plot 3: Volcano by Methylation Direction
################################################################################

cat("\n--- Creating Direction-colored Volcano ---\n")

# Custom color scheme for hyper/hypo methylation
dmr_data <- dmr_data %>%
    mutate(
        sig_direction = case_when(
            FDR < FDR_CUTOFF & logFC >= LOGFC_CUTOFF ~ "Hypermethylated",
            FDR < FDR_CUTOFF & logFC <= -LOGFC_CUTOFF ~ "Hypomethylated",
            FDR < FDR_CUTOFF ~ "Significant (|FC| < 2)",
            TRUE ~ "Not significant"
        )
    )

p_direction <- ggplot(dmr_data, aes(x = logFC, y = -log10(PValue), color = sig_direction)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(
        values = c(
            "Hypermethylated" = "#D73027",
            "Hypomethylated" = "#4575B4",
            "Significant (|FC| < 2)" = "#FEE090",
            "Not significant" = "grey70"
        ),
        name = "DMR Status"
    ) +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40") +
    labs(
        title = "TES vs GFP: Differential DNA Methylation",
        subtitle = "Colored by methylation direction",
        x = "log2 Fold Change (TES / GFP)",
        y = "-log10(P-value)",
        caption = paste0("Dashed lines: FDR = ", FDR_CUTOFF, ", |FC| = ", FC_CUTOFF)
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
    )

# Save direction volcano
pdf(file.path(OUTPUT_DIR, "volcano_direction.pdf"), width = 10, height = 8)
print(p_direction)
dev.off()
cat("Saved: volcano_direction.pdf\n")

################################################################################
# Summary
################################################################################

cat("\n========================================\n")
cat("Volcano Plot Generation Complete\n")
cat("========================================\n")
cat("\nOutput files in:", OUTPUT_DIR, "\n")
cat("  - volcano_main.pdf\n")
cat("  - volcano_promoters.pdf\n")
cat("  - volcano_direction.pdf\n")
cat("  - volcano_summary.txt\n")
cat("\nDone!\n")
