#!/usr/bin/env Rscript

################################################################################
# Script: 07_differential_methylation_MEDIPS_with_INPUT.R
# Purpose: Quantitative differential methylation analysis using MEDIPS
#          WITH INPUT NORMALIZATION for proper meDIP-seq analysis
#
# Description:
#   This version properly incorporates INPUT controls for normalization,
#   which corrects for:
#   - Genomic copy number variations
#   - Chromatin accessibility differences
#   - Sequencing depth variations
#
# Key Difference from Original Script:
#   - ISet1 and ISet2 are properly populated with INPUT samples
#   - INPUT samples are matched to their corresponding IP samples by group
#   - More rigorous normalization for meDIP-seq data
#
# Input:
#   - BAM files: ../results/04_filtered/*_filtered_dedup.bam (paired-end)
#   - Design matrix: ../config/design.txt
#   - Reference genome: BSgenome.Hsapiens.NCBI.GRCh38
#
# Output:
#   - Results saved to: results/07_differential_MEDIPS_INPUT_normalized/
#   - Does NOT overwrite original results in 07_differential_MEDIPS/
#
# Statistical Model:
#   - Window-based read counts (500bp windows, genome-wide)
#   - CpG density normalization (accounts for meDIP enrichment bias)
#   - INPUT normalization (accounts for copy number and accessibility)
#   - edgeR negative binomial GLM
#   - Significance: FDR < 0.05, |logFC| > 1
#
# Runtime: ~2-4 hours (genome-wide analysis)
# Memory: 32-64 GB
################################################################################

cat("=======================================================\n")
cat("meDIP-seq: MEDIPS Analysis WITH INPUT Normalization\n")
cat("=======================================================\n")
cat(paste("Start time:", Sys.time(), "\n\n"))

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
    library(MEDIPS)
    library(BSgenome.Hsapiens.NCBI.GRCh38)
    library(edgeR)
    library(GenomicRanges)
    library(rtracklayer)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
})

cat("Loaded R packages:\n")
cat(paste("  MEDIPS version:", packageVersion("MEDIPS"), "\n"))
cat(paste("  BSgenome.Hsapiens.NCBI.GRCh38 version:", packageVersion("BSgenome.Hsapiens.NCBI.GRCh38"), "\n"))
cat(paste("  edgeR version:", packageVersion("edgeR"), "\n\n"))

# Define paths - NEW OUTPUT DIRECTORY
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
bam_dir <- file.path(base_dir, "results/04_filtered")
out_dir <- file.path(base_dir, "results/07_differential_MEDIPS_INPUT_normalized") # NEW DIRECTORY
design_file <- file.path(base_dir, "config/design.txt")

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

cat("=======================================================\n")
cat("OUTPUT DIRECTORY: ", out_dir, "\n")
cat("(Original results preserved in 07_differential_MEDIPS/)\n")
cat("=======================================================\n\n")

# Read design matrix
cat("Reading experimental design...\n")
design <- read.delim(design_file, stringsAsFactors = FALSE, comment.char = "#")
cat(paste("Loaded", nrow(design), "samples\n"))
print(design)
cat("\n")

# MEDIPS parameters
cat("MEDIPS Analysis Parameters:\n")
window_size <- 500
bsg <- "BSgenome.Hsapiens.NCBI.GRCh38"
uniq <- 1e-3
extend <- 250
shift <- 0
chr_select <- c(1:22, "X", "Y")

cat(paste("  Window size:", window_size, "bp\n"))
cat(paste("  Genome:", bsg, "\n"))
cat(paste("  Extend reads to:", extend, "bp\n"))
cat(paste("  Chromosomes:", paste(chr_select, collapse = ", "), "\n"))
cat("  INPUT normalization: ENABLED\n\n")

################################################################################
# Step 1: Create MEDIPS sets for ALL samples (IP + INPUT)
################################################################################

cat("========================================\n")
cat("Step 1: Creating MEDIPS sets\n")
cat("========================================\n\n")

# TEACHING NOTE:
# A "MEDIPS set" is an object that stores the read counts for every window in the genome.
# It's the large data structure that holds our raw data before any fancy statistics.
medips_sets <- list()

for (i in 1:nrow(design)) {
    sample_name <- design$sample_name[i]
    sample_type <- design$sample_type[i]
    group <- design$group[i]

    bam_file <- file.path(bam_dir, paste0(sample_name, "_filtered_dedup.bam"))

    cat(paste("Processing:", sample_name, "(", group, "-", sample_type, ")\n"))

    if (!file.exists(bam_file)) {
        cat(paste("  ERROR: BAM file not found:", bam_file, "\n"))
        next
    }

    # Create MEDIPS set
    # This function scans the BAM file and counts reads in every 500bp window.
    # - extend: Extends reads to 250bp (av. fragment size) to cover the CpG site.
    # - uniq: Keeps only 1 read per position if duplicates weren't already removed.
    # - shift: 0 (we center reads by extending them, not shifting them).
    cat("  Creating MEDIPS set...\n")
    medips_set <- tryCatch(
        {
            MEDIPS.createSet(
                file = bam_file,
                BSgenome = bsg,
                extend = extend,
                shift = shift,
                uniq = uniq,
                window_size = window_size,
                chr.select = chr_select,
                paired = TRUE
            )
        },
        error = function(e) {
            cat(paste("  ERROR creating MEDIPS set:", e$message, "\n"))
            return(NULL)
        }
    )

    if (!is.null(medips_set)) {
        medips_sets[[sample_name]] <- medips_set
        cat(paste("  Successfully created MEDIPS set\n"))
        cat(paste("  Total windows:", length(medips_set@genome_count), "\n"))
        cat(paste("  Total reads:", sum(medips_set@genome_count), "\n\n"))
    }
}

cat(paste("Created", length(medips_sets), "MEDIPS sets\n\n"))

# Check if any MEDIPS sets were successfully created
if (length(medips_sets) == 0) {
    cat("ERROR: No MEDIPS sets were successfully created!\n")
    quit(status = 1)
}

# Save MEDIPS sets
cat("Saving MEDIPS sets to RData...\n")
save(medips_sets, file = "MEDIPS_sets.RData")
cat("Saved: MEDIPS_sets.RData\n\n")

################################################################################
# Step 2: Quality Control
################################################################################

cat("========================================\n")
cat("Step 2: Quality Control\n")
cat("========================================\n\n")

# CpG enrichment for IP samples
# This checks if the IP actually worked.
# It compares the density of CpGs in the pulled-down regions vs the whole genome.
# A high score (>1.5) means we successfully enriched for methylated (CpG-rich) DNA.
cat("Generating CpG enrichment analysis...\n")
enrichment_scores <- list()
for (sample_name in names(medips_sets)) {
    if (grepl("IP", sample_name)) {
        cat(paste("  CpG enrichment for:", sample_name, "\n"))
        bam_file <- file.path(bam_dir, paste0(sample_name, "_filtered_dedup.bam"))
        if (file.exists(bam_file)) {
            tryCatch(
                {
                    er <- MEDIPS.CpGenrich(
                        file = bam_file,
                        BSgenome = bsg,
                        extend = 0,
                        shift = 0,
                        uniq = uniq,
                        chr.select = chr_select,
                        paired = TRUE
                    )
                    score <- er$enrichment.score.GoGe
                    enrichment_scores[[sample_name]] <- score
                    cat(paste("    Enrichment score:", round(score, 2), "\n"))
                },
                error = function(e) {
                    cat(paste("    WARNING: Could not calculate enrichment:", e$message, "\n"))
                }
            )
        }
    }
}

if (length(enrichment_scores) > 0) {
    enrichment_df <- data.frame(
        sample = names(enrichment_scores),
        enrichment_score = unlist(enrichment_scores)
    )
    write.csv(enrichment_df, "QC_CpG_enrichment_scores.csv", row.names = FALSE)
    cat("\nCpG enrichment scores saved: QC_CpG_enrichment_scores.csv\n")
    cat("Mean enrichment score:", round(mean(unlist(enrichment_scores)), 2), "\n")
    cat("(Score > 1.5 indicates good IP quality)\n\n")
}

################################################################################
# Step 3: Calculate genome-wide CpG density
################################################################################

cat("========================================\n")
cat("Step 3: Calculate CpG density\n")
cat("========================================\n\n")

# TEACHING NOTE:
# The "Coupling Vector" is the key to MEDIPS normalization.
# It calculates the local density of CpG sites for every single window in the genome.
# Why? Because meDIP is biased: it pulls down CpG-rich regions more easily.
# We need this vector to mathmatically "correct" the read counts later,
# so a high count means "high methylation", not just "high CpG density".
cat("Calculating CpG density (coupling vector) for normalization...\n")
template_set <- medips_sets[[1]]

CpG <- MEDIPS.couplingVector(
    pattern = "CG",
    refObj = template_set
)

cat(paste("CpG coupling vector calculated for", length(CpG@genome_CF), "windows\n\n"))

################################################################################
# Step 4: Sample Organization and IP-INPUT Pairing
################################################################################

cat("========================================\n")
cat("Step 4: Sample Organization & INPUT Pairing\n")
cat("========================================\n\n")

# Separate IP and INPUT samples
ip_samples <- design[design$sample_type == "IP", "sample_name"]
input_samples <- design[design$sample_type == "INPUT", "sample_name"]

cat(paste("IP samples in design:", length(ip_samples), "\n"))
cat(paste("INPUT samples in design:", length(input_samples), "\n"))

# Filter for successfully created MEDIPS sets
ip_samples <- ip_samples[ip_samples %in% names(medips_sets)]
input_samples <- input_samples[input_samples %in% names(medips_sets)]

cat(paste("IP samples with data:", length(ip_samples), "\n"))
cat(paste("INPUT samples with data:", length(input_samples), "\n\n"))

# Use COMMON INPUT for all samples
# TESmut-1-INPUT is used as the common background control for all IP samples
common_input <- "TESmut-1-INPUT"

cat("========================================\n")
cat("USING COMMON INPUT FOR ALL SAMPLES\n")
cat("========================================\n")
cat(paste("Common INPUT sample:", common_input, "\n"))
cat("This INPUT will be used for both TES and GFP groups\n")
cat("Rationale: TES-1-INPUT and GFP-1-INPUT were not processed\n\n")

# Create IP-INPUT pairing map (all samples use common INPUT)
cat("Creating IP-INPUT sample pairs:\n")
sample_pairs <- list()
for (ip in ip_samples) {
    if (common_input %in% names(medips_sets)) {
        sample_pairs[[ip]] <- common_input
        cat(paste("  ", ip, " <-> ", common_input, "\n"))
    } else {
        cat(paste("  WARNING: Common INPUT", common_input, "not available for", ip, "\n"))
        sample_pairs[[ip]] <- NA
    }
}
cat("\n")

# Verify we have the common INPUT sample
if (!(common_input %in% names(medips_sets))) {
    cat(paste("ERROR: Common INPUT sample", common_input, "not available!\n"))
    cat("Cannot proceed with INPUT-normalized analysis.\n")
    quit(status = 1)
}

################################################################################
# Step 5: Differential Methylation Analysis WITH INPUT
################################################################################

cat("========================================\n")
cat("Step 5: Differential Methylation (INPUT Normalized)\n")
cat("========================================\n\n")

# Define contrasts - TES vs GFP only (TESmut excluded)
comparisons <- list(
    list(
        name = "TES_vs_GFP",
        group1 = "TES", group1_samples = grep("^TES-[0-9]+-IP$", ip_samples, value = TRUE),
        group2 = "GFP", group2_samples = grep("^GFP", ip_samples, value = TRUE)
    )
)

for (comp in comparisons) {
    cat(paste("\n==========================================\n"))
    cat(paste("Contrast:", comp$name, "\n"))
    cat(paste("==========================================\n\n"))

    cat(paste("Group 1 (", comp$group1, "):", paste(comp$group1_samples, collapse = ", "), "\n"))
    cat(paste("Group 2 (", comp$group2, "):", paste(comp$group2_samples, collapse = ", "), "\n\n"))

    # Get MEDIPS sets for IP samples
    group1_sets <- medips_sets[comp$group1_samples]
    group2_sets <- medips_sets[comp$group2_samples]

    if (length(group1_sets) == 0 || length(group2_sets) == 0) {
        cat("ERROR: Missing IP samples for this comparison\n")
        next
    }

    # Use COMMON INPUT for both groups
    # Both TES and GFP use TESmut-1-INPUT as the background control
    cat("Using common INPUT for both groups:\n")
    cat(paste("  Common INPUT:", common_input, "\n\n"))

    # Get INPUT MEDIPS sets (same for both groups)
    common_input_set <- medips_sets[[common_input]]

    if (is.null(common_input_set)) {
        cat(paste("ERROR: Common INPUT", common_input, "not found in MEDIPS sets!\n"))
        next
    }

    # Both groups use the same common INPUT
    group1_input_sets <- list(common_input_set)
    group2_input_sets <- list(common_input_set)
    names(group1_input_sets) <- common_input
    names(group2_input_sets) <- common_input

    cat(paste("  Group 1 (", comp$group1, ") INPUT:", common_input, "\n"))
    cat(paste("  Group 2 (", comp$group2, ") INPUT:", common_input, "\n\n"))

    # Perform differential methylation analysis WITH INPUT NORMALIZATION
    cat("Running MEDIPS differential methylation analysis WITH INPUT normalization...\n")
    cat("This may take 30-60 minutes...\n\n")

    # TEACHING NOTE:
    # This is the core function. It calls 'edgeR' to do the statistics.
    # - MSet: Methylation Sets (IP samples)
    # - ISet: Input Sets (Control samples). Used to correct for local background.
    # - CSet: CpG Coupling Vector. Used to correct for CpG density bias.
    # - diff.method = "edgeR": Uses Negative Binomial statistics (best for this data).
    # - MeDIP = TRUE: Tells it to apply the specific meDIP biases correction.
    diff_meth <- tryCatch(
        {
            MEDIPS.meth(
                MSet1 = group1_sets,
                MSet2 = group2_sets,
                CSet = CpG,
                ISet1 = group1_input_sets, # INPUT for group 1
                ISet2 = group2_input_sets, # INPUT for group 2
                p.adj = "BH", # Benjamini-Hochberg FDR correction
                diff.method = "edgeR",
                MeDIP = TRUE,
                CNV = FALSE,
                minRowSum = 10, # Ignore windows with <10 reads total
                diffnorm = "tmm" # TMM Normalization (Scales samples so libraries are comparable)
            )
        },
        error = function(e) {
            cat(paste("ERROR in differential analysis:", e$message, "\n"))
            return(NULL)
        }
    )

    if (is.null(diff_meth)) {
        cat("Skipping this comparison due to error\n")
        next
    }

    cat("Differential analysis completed WITH INPUT normalization\n\n")

    # Extract results
    # The result is a table with info for every window.
    # logFC: Log2 Fold Change (positive = higher in Group 1)
    # logCPM: Log Counts Per Million (how much signal is there?)
    # FDR: False Discovery Rate (adjusted p-value)
    dmr_results <- data.frame(
        chr = diff_meth[, "chr"],
        start = diff_meth[, "start"],
        stop = diff_meth[, "stop"],
        CpG_count = diff_meth[, "CF"],
        mean_group1 = diff_meth[, "MSets1.counts.mean"],
        mean_group2 = diff_meth[, "MSets2.counts.mean"],
        logFC = diff_meth[, "edgeR.logFC"],
        logCPM = diff_meth[, "edgeR.logCPM"],
        PValue = diff_meth[, "edgeR.p.value"],
        FDR = diff_meth[, "edgeR.adj.p.value"]
    )

    # Add fold change
    dmr_results$fold_change <- 2^dmr_results$logFC

    # Save all results
    out_file_all <- paste0(comp$name, "_all_windows.csv")
    write.csv(dmr_results, out_file_all, row.names = FALSE)
    cat(paste("All windows saved:", out_file_all, "\n"))

    # Filter for significant DMRs (FDR < 0.05)
    sig_dmr <- dmr_results[dmr_results$FDR < 0.05 & !is.na(dmr_results$FDR), ]
    out_file_sig <- paste0(comp$name, "_DMRs_FDR05.csv")
    write.csv(sig_dmr, out_file_sig, row.names = FALSE)
    cat(paste("Significant DMRs (FDR<0.05):", nrow(sig_dmr), "->", out_file_sig, "\n"))

    # Filter for strong effect (FDR < 0.05, |logFC| > 1)
    strong_dmr <- sig_dmr[abs(sig_dmr$logFC) > 1, ]
    out_file_strong <- paste0(comp$name, "_DMRs_FDR05_FC2.csv")
    write.csv(strong_dmr, out_file_strong, row.names = FALSE)
    cat(paste("Strong DMRs (FDR<0.05, FC>2):", nrow(strong_dmr), "->", out_file_strong, "\n"))

    # Export BED files for IGV
    if (nrow(strong_dmr) > 0) {
        # Hypermethylated (increased in group1)
        hyper_dmr <- strong_dmr[strong_dmr$logFC > 1, ]
        if (nrow(hyper_dmr) > 0) {
            hyper_bed <- data.frame(
                chrom = hyper_dmr$chr,
                start = hyper_dmr$start,
                end = hyper_dmr$stop,
                name = paste0("hyper_", 1:nrow(hyper_dmr)),
                score = pmin(1000, -log10(hyper_dmr$FDR) * 100),
                strand = "."
            )
            write.table(hyper_bed, paste0(comp$name, "_hypermethylated.bed"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
            )
            cat(paste(
                "Hypermethylated regions:", nrow(hyper_dmr), "->",
                paste0(comp$name, "_hypermethylated.bed\n")
            ))
        }

        # Hypomethylated (decreased in group1)
        hypo_dmr <- strong_dmr[strong_dmr$logFC < -1, ]
        if (nrow(hypo_dmr) > 0) {
            hypo_bed <- data.frame(
                chrom = hypo_dmr$chr,
                start = hypo_dmr$start,
                end = hypo_dmr$stop,
                name = paste0("hypo_", 1:nrow(hypo_dmr)),
                score = pmin(1000, -log10(hypo_dmr$FDR) * 100),
                strand = "."
            )
            write.table(hypo_bed, paste0(comp$name, "_hypomethylated.bed"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
            )
            cat(paste(
                "Hypomethylated regions:", nrow(hypo_dmr), "->",
                paste0(comp$name, "_hypomethylated.bed\n")
            ))
        }
    }

    # Generate plots
    cat("\nGenerating plots...\n")

    # MA plot
    # TEACHING NOTE:
    # MA Plot checks for bias.
    # X-axis: Average strength of signal (Mean Abundance).
    # Y-axis: Log Fold Change (M).
    # You want the cloud to be centered on Y=0. If it curves up or down, normalization failed.
    # Red dots are significant changes.
    pdf(paste0(comp$name, "_MA_plot.pdf"), width = 8, height = 6)
    plot(dmr_results$logCPM, dmr_results$logFC,
        pch = 16, cex = 0.3, col = ifelse(dmr_results$FDR < 0.05, "red", "gray50"),
        xlab = "Average log2 CPM", ylab = paste0("log2 FC (", comp$group1, " vs ", comp$group2, ")"),
        main = paste("MA Plot:", comp$name, "(INPUT Normalized)")
    )
    abline(h = c(-1, 0, 1), col = c("blue", "black", "blue"), lty = c(2, 1, 2))
    legend("topright",
        legend = c("FDR < 0.05", "Not significant"),
        col = c("red", "gray50"), pch = 16
    )
    dev.off()
    cat(paste("  MA plot:", paste0(comp$name, "_MA_plot.pdf\n")))

    # Volcano plot
    # TEACHING NOTE:
    # Visualizes significance vs effect size.
    # Top-Right: Significantly Hyper-methylated.
    # Top-Left: Significantly Hypo-methylated.
    pdf(paste0(comp$name, "_volcano_plot.pdf"), width = 8, height = 6)
    plot(dmr_results$logFC, -log10(dmr_results$PValue),
        pch = 16, cex = 0.3, col = ifelse(dmr_results$FDR < 0.05 & abs(dmr_results$logFC) > 1, "red", "gray50"),
        xlab = paste0("log2 FC (", comp$group1, " vs ", comp$group2, ")"),
        ylab = "-log10(P-value)",
        main = paste("Volcano Plot:", comp$name, "(INPUT Normalized)")
    )
    abline(v = c(-1, 1), col = "blue", lty = 2)
    abline(h = -log10(0.05), col = "blue", lty = 2)
    legend("topright",
        legend = c("FDR<0.05 & |FC|>2", "Not significant"),
        col = c("red", "gray50"), pch = 16
    )
    dev.off()
    cat(paste("  Volcano plot:", paste0(comp$name, "_volcano_plot.pdf\n")))

    # Summary statistics
    cat("\n--- Summary Statistics ---\n")
    cat(paste("Total windows analyzed:", nrow(dmr_results), "\n"))
    cat(paste("Significant DMRs (FDR<0.05):", nrow(sig_dmr), "\n"))
    cat(paste(
        "  Hypermethylated (", comp$group1, " > ", comp$group2, "):",
        sum(sig_dmr$logFC > 0), "\n"
    ))
    cat(paste(
        "  Hypomethylated (", comp$group1, " < ", comp$group2, "):",
        sum(sig_dmr$logFC < 0), "\n"
    ))
    cat(paste("Strong DMRs (FDR<0.05, |FC|>2):", nrow(strong_dmr), "\n\n"))
}

################################################################################
# Step 6: Sample correlation and PCA
################################################################################

cat("========================================\n")
cat("Step 6: Sample QC and Clustering\n")
cat("========================================\n\n")

# Create count matrix for all IP samples
cat("Creating count matrix for correlation analysis...\n")
count_matrix <- sapply(ip_samples, function(s) {
    if (s %in% names(medips_sets)) {
        medips_sets[[s]]@genome_count
    } else {
        rep(NA, length(template_set@genome_count))
    }
})

# Remove windows with all zeros
count_matrix <- count_matrix[rowSums(count_matrix, na.rm = TRUE) > 0, ]
cat(paste("Count matrix:", nrow(count_matrix), "windows x", ncol(count_matrix), "samples\n\n"))

# Correlation heatmap
# TEACHING NOTE:
# Shows how similar the samples are to each other.
# We expect replicates (TES-1, TES-2) to cluster together.
# If they don't, we have a "batch effect" or a bad sample.
cat("Generating correlation heatmap...\n")
cor_matrix <- cor(count_matrix, method = "spearman")
pdf("sample_correlation_heatmap.pdf", width = 8, height = 8)
pheatmap(cor_matrix,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    main = "Sample Correlation (Spearman) - INPUT Normalized Analysis",
    display_numbers = TRUE,
    number_format = "%.2f"
)
dev.off()
cat("Correlation heatmap saved: sample_correlation_heatmap.pdf\n\n")

# PCA
cat("Performing PCA...\n")
norm_counts <- log2(t(t(count_matrix) / colSums(count_matrix) * 1e6) + 1)

pca_result <- prcomp(t(norm_counts), scale. = TRUE)
pca_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample = rownames(pca_result$x),
    group = sapply(rownames(pca_result$x), function(s) {
        if (grepl("^TES-[0-9]+-IP$", s)) {
            return("TES")
        }
        if (grepl("^GFP", s)) {
            return("GFP")
        }
        return("Other")
    })
)

pdf("PCA_plot.pdf", width = 8, height = 6)
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    labs(
        title = "PCA of meDIP-seq Samples (INPUT Normalized Analysis)",
        x = paste0("PC1 (", pca_var[1], "% variance)"),
        y = paste0("PC2 (", pca_var[2], "% variance)")
    ) +
    theme_bw() +
    theme(legend.position = "right")
dev.off()
cat("PCA plot saved: PCA_plot.pdf\n\n")

################################################################################
# Final Summary
################################################################################

cat("========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n\n")

cat(paste("End time:", Sys.time(), "\n\n"))

cat("OUTPUT DIRECTORY:", out_dir, "\n")
cat("(Original results preserved in 07_differential_MEDIPS/)\n\n")

cat("Key differences from original analysis:\n")
cat("  - INPUT normalization: ENABLED\n")
cat("  - COMMON INPUT used: TESmut-1-INPUT\n")
cat("  - TES samples paired with: TESmut-1-INPUT (TES-1-INPUT not available)\n")
cat("  - GFP samples paired with: TESmut-1-INPUT (GFP-1-INPUT not processed)\n")
cat("  - This corrects for copy number and accessibility bias\n")
cat("  - Using common INPUT is valid as all samples are from the same cell line\n\n")

cat("Output files:\n")
cat("1. MEDIPS sets:\n")
cat("   - MEDIPS_sets.RData: Saved MEDIPS objects\n\n")

cat("2. Quality control:\n")
cat("   - QC_CpG_enrichment_scores.csv: IP enrichment validation\n")
cat("   - sample_correlation_heatmap.pdf: Sample similarity\n")
cat("   - PCA_plot.pdf: Sample clustering\n\n")

cat("3. Differential methylation results:\n")
cat("   - *_all_windows.csv: All genomic windows with statistics\n")
cat("   - *_DMRs_FDR05.csv: Significant DMRs (FDR < 0.05)\n")
cat("   - *_DMRs_FDR05_FC2.csv: Strong DMRs (FDR < 0.05, FC > 2)\n")
cat("   - *_hypermethylated.bed: BED file for IGV\n")
cat("   - *_hypomethylated.bed: BED file for IGV\n")
cat("   - *_MA_plot.pdf: MA plot\n")
cat("   - *_volcano_plot.pdf: Volcano plot\n\n")

cat("Comparison with original analysis:\n")
cat("  To compare INPUT-normalized vs non-normalized results:\n")
cat("  1. Compare DMR counts between analyses\n")
cat("  2. Check overlap of significant DMRs\n")
cat("  3. INPUT normalization typically reduces false positives\n\n")

cat("=======================================================\n")
cat("MEDIPS analysis WITH INPUT completed successfully!\n")
cat("=======================================================\n")
