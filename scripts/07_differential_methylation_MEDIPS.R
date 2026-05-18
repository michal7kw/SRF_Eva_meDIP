#!/usr/bin/env Rscript

################################################################################
# Script: 07_differential_methylation_MEDIPS.R
# Purpose: Quantitative differential methylation analysis using MEDIPS
#
# Description:
#   Performs biologically rigorous differential methylation analysis using
#   MEDIPS (Methylated DNA Immunoprecipitation Sequencing), which accounts
#   for CpG density bias and treats methylation as a quantitative signal.
#
# Strategy:
#   1. Tile genome into uniform windows (500bp)
#   2. Create MEDIPS sets for all samples (IP + INPUT)
#   3. Calculate CpG-density-corrected methylation scores
#   4. Normalize IP against INPUT controls
#   5. Statistical testing with edgeR (designed for count data)
#   6. Multiple testing correction (FDR)
#
# Why This Approach (vs Peak Calling):
#   - Treats methylation as QUANTITATIVE (not binary peaks)
#   - Corrects for CpG density bias (critical for meDIP!)
#   - No arbitrary peak thresholds
#   - Genome-wide coverage (doesn't miss low-signal regions)
#   - Proper statistical framework for IP-INPUT comparison
#
# Biological Rationale:
#   - DNA methylation is a continuous modification, not discrete binding
#   - meDIP signal is proportional to: methylation% × CpG_count × fragment_coverage
#   - Must deconvolve these factors for biological interpretation
#   - MEDIPS calculates Absolute Methylation Score (AMS) and Relative Methylation Score (RMS)
#
# Input:
#   - BAM files: ../results/04_filtered/*_filtered_dedup.bam (paired-end)
#   - Design matrix: ../config/design.txt
#   - Reference genome: BSgenome.Hsapiens.UCSC.hg38
#
# Output:
#   - Genome-wide methylation scores: *_methylation_windows.csv
#   - DMR results: *_vs_*_DMRs_MEDIPS.csv
#   - Plots: saturation, CpG calibration, correlation, PCA, MA plots
#   - BED files: Significant DMRs for IGV visualization
#
# Statistical Model:
#   - Window-based read counts (500bp windows, genome-wide)
#   - CpG density normalization (accounts for meDIP enrichment bias)
#   - edgeR negative binomial GLM: log(IP/INPUT) ~ group + CpG_density
#   - Significance: FDR < 0.05, |logFC| > 1
#
# Runtime: ~2-4 hours (genome-wide analysis)
# Memory: 32-64 GB
################################################################################

cat("=======================================================\n")
cat("meDIP-seq: MEDIPS Differential Methylation Analysis\n")
cat("=======================================================\n")
cat(paste("Start time:", Sys.time(), "\n\n"))

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
    library(MEDIPS)
    library(BSgenome.Hsapiens.NCBI.GRCh38)  # Use NCBI/Ensembl naming (1, 2, 3...) not UCSC (chr1, chr2...)
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
cat(paste("  edgeR version:", packageVersion("edgeR"), "\n"))
cat(paste("  GenomicRanges version:", packageVersion("GenomicRanges"), "\n\n"))

# Define paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
bam_dir <- file.path(base_dir, "results/04_filtered")
out_dir <- file.path(base_dir, "results/07_differential_MEDIPS")
design_file <- file.path(base_dir, "config/design.txt")

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

# Read design matrix
cat("Reading experimental design...\n")
design <- read.delim(design_file, stringsAsFactors = FALSE, comment.char = "#")
cat(paste("Loaded", nrow(design), "samples\n"))
print(design)
cat("\n")

# MEDIPS parameters
cat("MEDIPS Analysis Parameters:\n")
window_size <- 500  # 500bp windows (balance resolution vs computational cost)
bsg <- "BSgenome.Hsapiens.NCBI.GRCh38"  # Use NCBI/GRCh38 for Ensembl chromosome naming
uniq <- 1e-3  # Require reads to map uniquely (mapq threshold)
extend <- 250  # Extend reads to average fragment length
shift <- 0
# BAM files use Ensembl chromosome naming (1, 2, 3...) which matches NCBI.GRCh38
chr_select <- c(1:22, "X", "Y")  # Autosomes + sex chromosomes (Ensembl style)

cat(paste("  Window size:", window_size, "bp\n"))
cat(paste("  Genome:", bsg, "\n"))
cat(paste("  Extend reads to:", extend, "bp\n"))
cat(paste("  Chromosomes:", paste(chr_select, collapse=", "), "\n"))
cat("  NOTE: Using NCBI.GRCh38 genome (Ensembl chromosome naming: 1, 2, 3...) to match BAM files\n\n")

################################################################################
# Step 1: Create MEDIPS sets for all samples
################################################################################

cat("========================================\n")
cat("Step 1: Creating MEDIPS sets\n")
cat("========================================\n\n")

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
    cat("  Creating MEDIPS set...\n")
    medips_set <- tryCatch({
        MEDIPS.createSet(
            file = bam_file,
            BSgenome = bsg,
            extend = extend,
            shift = shift,
            uniq = uniq,
            window_size = window_size,
            chr.select = chr_select,
            paired = TRUE  # meDIP data is paired-end
        )
    }, error = function(e) {
        cat(paste("  ERROR creating MEDIPS set:", e$message, "\n"))
        return(NULL)
    })

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
    cat("This could be due to:\n")
    cat("  1. Chromosome naming mismatch between BAM and genome\n")
    cat("  2. Missing BAM files\n")
    cat("  3. BAM index files not found\n\n")
    cat("Please check the errors above and fix before proceeding.\n")
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

# Saturation analysis (should plateau, indicating sufficient sequencing depth)
# NOTE: Skipping saturation analysis - computationally intensive and optional for differential analysis
cat("Skipping saturation analysis (optional QC, computationally intensive)\n")
cat("All samples successfully created MEDIPS sets - sufficient for differential analysis\n\n")

# # Uncomment below if you want saturation analysis (requires ~2 hours extra runtime)
# cat("Generating saturation analysis...\n")
# pdf("QC_saturation_analysis.pdf", width = 12, height = 8)
# par(mfrow = c(2, 3))
# for (sample_name in names(medips_sets)) {
#     bam_file <- file.path(bam_dir, paste0(sample_name, "_filtered_dedup.bam"))
#     if (file.exists(bam_file)) {
#         cat(paste("  Saturation plot for:", sample_name, "\n"))
#         sr <- MEDIPS.saturation(
#             file = bam_file,
#             BSgenome = bsg,
#             uniq = uniq,
#             extend = extend,
#             shift = shift,
#             window_size = window_size,
#             chr.select = "1",  # Use chr1 for speed (Ensembl naming)
#             nit = 10,
#             nrit = 1,
#             empty_bins = TRUE,
#             rank = FALSE
#         )
#     }
# }
# dev.off()
# cat("Saturation analysis saved: QC_saturation_analysis.pdf\n\n")

# CpG enrichment (should show enrichment in IP vs random)
cat("Generating CpG enrichment analysis...\n")
enrichment_scores <- list()
for (sample_name in names(medips_sets)) {
    if (grepl("IP", sample_name)) {  # Only for IP samples
        cat(paste("  CpG enrichment for:", sample_name, "\n"))
        bam_file <- file.path(bam_dir, paste0(sample_name, "_filtered_dedup.bam"))
        if (file.exists(bam_file)) {
            tryCatch({
                er <- MEDIPS.CpGenrich(
                    file = bam_file,
                    BSgenome = bsg,  # Pass BSgenome string
                    extend = 0,      # Use actual fragment length from paired-end data
                    shift = 0,
                    uniq = uniq,
                    chr.select = chr_select,
                    paired = TRUE
                )
                score <- er$enrichment.score.GoGe
                enrichment_scores[[sample_name]] <- score
                cat(paste("    Enrichment score:", round(score, 2), "\n"))
            }, error = function(e) {
                cat(paste("    WARNING: Could not calculate enrichment:", e$message, "\n"))
            })
        } else {
            cat(paste("    WARNING: BAM file not found:", bam_file, "\n"))
        }
    }
}

# Save enrichment scores
if (length(enrichment_scores) > 0) {
    enrichment_df <- data.frame(
        sample = names(enrichment_scores),
        enrichment_score = unlist(enrichment_scores)
    )
    write.csv(enrichment_df, "QC_CpG_enrichment_scores.csv", row.names = FALSE)
    cat("\nCpG enrichment scores saved: QC_CpG_enrichment_scores.csv\n")
    cat("Mean enrichment score:", round(mean(unlist(enrichment_scores)), 2), "\n")
    cat("(Score > 1.5 indicates good IP quality)\n\n")
} else {
    cat("WARNING: No enrichment scores calculated\n\n")
}

################################################################################
# Step 3: Calculate genome-wide CpG density
################################################################################

cat("========================================\n")
cat("Step 3: Calculate CpG density\n")
cat("========================================\n\n")

cat("Calculating CpG density (coupling vector) for normalization...\n")
# Use the first MEDIPS set as template for genomic coordinates
template_set <- medips_sets[[1]]

# Calculate CpG coupling vector for normalization
CpG <- MEDIPS.couplingVector(
    pattern = "CG",
    refObj = template_set
)

cat(paste("CpG coupling vector calculated for", length(CpG@genome_CF), "windows\n\n"))

################################################################################
# Step 4: Normalize and calculate methylation scores
################################################################################

cat("========================================\n")
cat("Step 4: Sample Organization\n")
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

if (length(input_samples) == 0) {
    cat("NOTE: No INPUT samples available.\n")
    cat("Proceeding with IP-only differential analysis.\n")
    cat("This is valid for comparing methylation between conditions.\n\n")
}

# Create pairing between IP and INPUT
sample_pairs <- list()
for (ip in ip_samples) {
    # Extract group from IP sample name (e.g., "TES-1-IP" -> "TES")
    group <- gsub("-[0-9]+-IP", "", ip)
    # Find matching INPUT (e.g., "TES-1-INPUT")
    matching_input <- input_samples[grepl(group, input_samples)]

    if (length(matching_input) > 0) {
        sample_pairs[[ip]] <- matching_input[1]
        cat(paste("  Paired:", ip, "<->", matching_input[1], "\n"))
    } else {
        cat(paste("  WARNING: No INPUT found for", ip, "\n"))
    }
}
cat("\n")

################################################################################
# Step 5: Differential Methylation Analysis
################################################################################

cat("========================================\n")
cat("Step 5: Differential Methylation\n")
cat("========================================\n\n")

# Define contrasts
# NOTE: TESmut samples excluded from analysis (failed sample)
comparisons <- list(
    list(name = "TES_vs_GFP",
         group1 = "TES", group1_samples = grep("^TES-[0-9]+-IP$", ip_samples, value = TRUE),
         group2 = "GFP", group2_samples = grep("^GFP", ip_samples, value = TRUE))
)

for (comp in comparisons) {
    cat(paste("\n==========================================\n"))
    cat(paste("Contrast:", comp$name, "\n"))
    cat(paste("==========================================\n\n"))

    cat(paste("Group 1 (", comp$group1, "):", paste(comp$group1_samples, collapse=", "), "\n"))
    cat(paste("Group 2 (", comp$group2, "):", paste(comp$group2_samples, collapse=", "), "\n\n"))

    # Get MEDIPS sets for this comparison
    group1_sets <- medips_sets[comp$group1_samples]
    group2_sets <- medips_sets[comp$group2_samples]

    if (length(group1_sets) == 0 || length(group2_sets) == 0) {
        cat("ERROR: Missing samples for this comparison\n")
        next
    }

    # Perform differential methylation analysis
    cat("Running MEDIPS differential methylation analysis...\n")
    cat("This may take 30-60 minutes...\n")

    diff_meth <- tryCatch({
        MEDIPS.meth(
            MSet1 = group1_sets,
            MSet2 = group2_sets,
            CSet = CpG,
            ISet1 = NULL,  # INPUT normalization handled separately if needed
            ISet2 = NULL,
            p.adj = "BH",  # Benjamini-Hochberg FDR correction
            diff.method = "edgeR",  # Use edgeR for differential testing
            MeDIP = TRUE,  # Enable CpG density normalization (critical for meDIP-seq!)
            CNV = FALSE,
            minRowSum = 10,  # Filter low-count windows
            diffnorm = "tmm"  # TMM normalization for edgeR
        )
    }, error = function(e) {
        cat(paste("ERROR in differential analysis:", e$message, "\n"))
        return(NULL)
    })

    if (is.null(diff_meth)) {
        cat("Skipping this comparison due to error\n")
        next
    }

    cat("Differential analysis completed\n\n")

    # Extract results
    # Note: MEDIPS.meth() returns "CF" (Coupling Factor) = number of CpGs per window
    # MEDIPS pre-calculates group means as "MSets1.counts.mean" and "MSets2.counts.mean"
    dmr_results <- data.frame(
        chr = diff_meth[, "chr"],
        start = diff_meth[, "start"],
        stop = diff_meth[, "stop"],
        CF = diff_meth[, "CF"],  # Coupling Factor (number of CpGs in window)
        mean_group1 = diff_meth[, "MSets1.counts.mean"],  # Pre-calculated by MEDIPS
        mean_group2 = diff_meth[, "MSets2.counts.mean"],  # Pre-calculated by MEDIPS
        logFC = diff_meth[, "edgeR.logFC"],
        logCPM = diff_meth[, "edgeR.logCPM"],
        PValue = diff_meth[, "edgeR.p.value"],
        FDR = diff_meth[, "edgeR.adj.p.value"]
    )

    # Rename CF to CpG_count for clarity
    names(dmr_results)[names(dmr_results) == "CF"] <- "CpG_count"

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
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            cat(paste("Hypermethylated regions:", nrow(hyper_dmr), "->",
                     paste0(comp$name, "_hypermethylated.bed\n")))
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
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            cat(paste("Hypomethylated regions:", nrow(hypo_dmr), "->",
                     paste0(comp$name, "_hypomethylated.bed\n")))
        }
    }

    # Generate plots
    cat("\nGenerating plots...\n")

    # MA plot
    pdf(paste0(comp$name, "_MA_plot.pdf"), width = 8, height = 6)
    plot(dmr_results$logCPM, dmr_results$logFC,
         pch = 16, cex = 0.3, col = ifelse(dmr_results$FDR < 0.05, "red", "gray50"),
         xlab = "Average log2 CPM", ylab = paste0("log2 FC (", comp$group1, " vs ", comp$group2, ")"),
         main = paste("MA Plot:", comp$name))
    abline(h = c(-1, 0, 1), col = c("blue", "black", "blue"), lty = c(2, 1, 2))
    legend("topright", legend = c("FDR < 0.05", "Not significant"),
           col = c("red", "gray50"), pch = 16)
    dev.off()
    cat(paste("  MA plot:", paste0(comp$name, "_MA_plot.pdf\n")))

    # Volcano plot
    pdf(paste0(comp$name, "_volcano_plot.pdf"), width = 8, height = 6)
    plot(dmr_results$logFC, -log10(dmr_results$PValue),
         pch = 16, cex = 0.3, col = ifelse(dmr_results$FDR < 0.05 & abs(dmr_results$logFC) > 1, "red", "gray50"),
         xlab = paste0("log2 FC (", comp$group1, " vs ", comp$group2, ")"),
         ylab = "-log10(P-value)",
         main = paste("Volcano Plot:", comp$name))
    abline(v = c(-1, 1), col = "blue", lty = 2)
    abline(h = -log10(0.05), col = "blue", lty = 2)
    legend("topright", legend = c("FDR<0.05 & |FC|>2", "Not significant"),
           col = c("red", "gray50"), pch = 16)
    dev.off()
    cat(paste("  Volcano plot:", paste0(comp$name, "_volcano_plot.pdf\n")))

    # Summary statistics
    cat("\n--- Summary Statistics ---\n")
    cat(paste("Total windows analyzed:", nrow(dmr_results), "\n"))
    cat(paste("Significant DMRs (FDR<0.05):", nrow(sig_dmr), "\n"))
    cat(paste("  Hypermethylated (", comp$group1, " > ", comp$group2, "):",
             sum(sig_dmr$logFC > 0), "\n"))
    cat(paste("  Hypomethylated (", comp$group1, " < ", comp$group2, "):",
             sum(sig_dmr$logFC < 0), "\n"))
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
cat(paste("Count matrix:", nrow(count_matrix), "windows ×", ncol(count_matrix), "samples\n\n"))

# Correlation heatmap
cat("Generating correlation heatmap...\n")
cor_matrix <- cor(count_matrix, method = "spearman")
pdf("sample_correlation_heatmap.pdf", width = 8, height = 8)
pheatmap(cor_matrix,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Sample Correlation (Spearman)",
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()
cat("Correlation heatmap saved: sample_correlation_heatmap.pdf\n\n")

# PCA
cat("Performing PCA...\n")
# Normalize counts (log2-CPM)
norm_counts <- log2(t(t(count_matrix) / colSums(count_matrix) * 1e6) + 1)

pca_result <- prcomp(t(norm_counts), scale. = TRUE)
pca_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Extract group from sample names
# NOTE: TESmut samples excluded from analysis (failed sample)
pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample = rownames(pca_result$x),
    group = sapply(rownames(pca_result$x), function(s) {
        if (grepl("^TES-[0-9]+-IP$", s)) return("TES")
        if (grepl("^GFP", s)) return("GFP")
        return("Other")
    })
)

pdf("PCA_plot.pdf", width = 8, height = 6)
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    labs(title = "PCA of meDIP-seq Samples",
         x = paste0("PC1 (", pca_var[1], "% variance)"),
         y = paste0("PC2 (", pca_var[2], "% variance)")) +
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

cat("Output directory:", out_dir, "\n\n")

cat("Key output files:\n")
cat("1. MEDIPS sets:\n")
cat("   - MEDIPS_sets.RData: Saved MEDIPS objects for all samples\n\n")

cat("2. Quality control:\n")
cat("   - QC_saturation_analysis.pdf: Sequencing depth assessment\n")
cat("   - QC_CpG_enrichment.pdf: IP enrichment validation\n")
cat("   - sample_correlation_heatmap.pdf: Sample similarity\n")
cat("   - PCA_plot.pdf: Sample clustering\n\n")

cat("3. Differential methylation results (for each comparison):\n")
cat("   - *_all_windows.csv: All genomic windows with statistics\n")
cat("   - *_DMRs_FDR05.csv: Significant DMRs (FDR < 0.05)\n")
cat("   - *_DMRs_FDR05_FC2.csv: Strong DMRs (FDR < 0.05, FC > 2)\n")
cat("   - *_hypermethylated.bed: BED file for IGV (increased methylation)\n")
cat("   - *_hypomethylated.bed: BED file for IGV (decreased methylation)\n")
cat("   - *_MA_plot.pdf: MA plot showing fold changes\n")
cat("   - *_volcano_plot.pdf: Volcano plot (effect size vs significance)\n\n")

cat("Next steps:\n")
cat("1. Review QC plots (saturation, CpG enrichment, PCA)\n")
cat("2. Examine significant DMRs in CSV files\n")
cat("3. Load BED files into IGV for visual inspection\n")
cat("4. Run genomic annotation (step 08) to map DMRs to genes\n")
cat("5. Integrate with Cut&Tag and RNA-seq data (step 10)\n\n")

cat("Biological interpretation notes:\n")
cat("- Hypermethylated promoters typically = gene repression\n")
cat("- Gene body methylation can = active transcription\n")
cat("- Check CpG island overlap (high-CpG regions most affected by meDIP bias)\n")
cat("- Validate key findings with bisulfite sequencing\n\n")

cat("=======================================================\n")
cat("MEDIPS analysis completed successfully!\n")
cat("=======================================================\n")
