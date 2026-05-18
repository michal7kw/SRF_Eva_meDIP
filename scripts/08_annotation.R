#!/usr/bin/env Rscript

################################################################################
# Script: 08_annotation.R
# Purpose: Annotate DMRs from MEDIPS with genomic features and CpG islands
#
# Description:
#   Annotates differentially methylated regions (DMRs) from MEDIPS analysis
#   with genomic context (promoter, exon, intron, intergenic) and overlaps
#   with CpG islands. Generates plots showing genomic distribution.
#
# Why Annotation Is Essential for meDIP-seq:
#   - Promoter methylation: Strongest gene expression effects (gene silencing)
#   - Gene body methylation: Often associated with active transcription
#   - CpG islands: Regulatory significance, typically at promoters
#   - Enhancers/intergenic: Distal regulatory elements
#   - Maps DMRs to genes: Required for RNA-seq integration
#
# Input:
#   - MEDIPS DMR files: ../results/07_differential_MEDIPS/*_DMRs_FDR05.csv
#   - Format: chr, start, stop, CpG_count, logFC, FDR, etc.
#   - TxDb: TxDb.Hsapiens.UCSC.hg38.knownGene (for genomic features)
#
# Output:
#   - Annotated DMRs: *_annotated.csv (with gene names, distances, features)
#   - Distribution plots: *_genomic_distribution.pdf
#   - Summary statistics: annotation_summary.txt
#   - Promoter DMRs: *_promoter_DMRs.csv (for motif analysis and integration)
#
# Tools: ChIPseeker for genomic annotation
#
# Runtime: ~10-30 minutes
################################################################################

cat("=================================================\n")
cat("meDIP-seq: MEDIPS DMR Annotation\n")
cat("=================================================\n")
cat(paste("Start:", Sys.time(), "\n\n"))

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(ggplot2)
    library(dplyr)
})

# Paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
in_dir <- file.path(base_dir, "results/07_differential_MEDIPS")  # FIXED: Use MEDIPS output directory
out_dir <- file.path(base_dir, "results/08_annotation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

# Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cat("Loaded hg38 transcript database (UCSC naming: chr1, chr2...)\n\n")

################################################################################
# Annotate DMRs for each contrast
################################################################################

# Find all MEDIPS DMR files
dmr_files <- list.files(in_dir, pattern = "*_DMRs_FDR05_FC2.csv$", full.names = TRUE)
cat(paste("Found", length(dmr_files), "MEDIPS DMR files to annotate\n\n"))

if (length(dmr_files) == 0) {
    cat("ERROR: No DMR files found in", in_dir, "\n")
    cat("Expected files matching pattern: *_DMRs_FDR05_FC2.csv\n")
    cat("Please check that MEDIPS analysis (step 07) completed successfully.\n")
    quit(status = 1)
}

annotation_summary <- list()

for (dmr_file in dmr_files) {
    contrast_name <- gsub("_DMRs_FDR05_FC2.csv$", "", basename(dmr_file))
    cat(paste("==========================================\n"))
    cat(paste("Processing:", contrast_name, "\n"))
    cat(paste("==========================================\n"))

    dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)
    cat(paste("Loaded", nrow(dmrs), "DMRs\n"))

    if (nrow(dmrs) == 0) {
        cat("  No DMRs found, skipping\n\n")
        next
    }

    # IMPORTANT: MEDIPS uses Ensembl chromosome naming (1, 2, 3...)
    # ChIPseeker/UCSC uses "chr1", "chr2"... format
    # Need to add "chr" prefix for compatibility
    cat("Converting chromosome names: Ensembl (1,2,3...) -> UCSC (chr1,chr2...)\n")
    dmrs$seqnames <- paste0("chr", dmrs$chr)

    # Handle sex chromosomes
    dmrs$seqnames <- gsub("chrX", "chrX", dmrs$seqnames)
    dmrs$seqnames <- gsub("chrY", "chrY", dmrs$seqnames)

    # Convert to GRanges
    dmrs_gr <- makeGRangesFromDataFrame(dmrs,
        seqnames.field = "seqnames",
        start.field = "start",
        end.field = "stop",  # MEDIPS uses "stop" not "end"
        keep.extra.columns = TRUE
    )

    cat(paste("Created GRanges object with", length(dmrs_gr), "regions\n"))

    # Annotate with ChIPseeker
    cat("Annotating with genomic features...\n")
    dmr_anno <- annotatePeak(dmrs_gr,
                             tssRegion = c(-2000, 500),  # Promoter = 2kb upstream to 500bp downstream
                             TxDb = txdb,
                             annoDb = "org.Hs.eg.db")

    # Convert to data frame
    dmr_anno_df <- as.data.frame(dmr_anno)

    # Add original MEDIPS columns
    dmr_anno_df <- dmr_anno_df %>%
        mutate(
            contrast = contrast_name,
            direction = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated")
        )

    # Save annotated DMRs
    out_name <- paste0(contrast_name, "_annotated.csv")
    write.csv(dmr_anno_df, out_name, row.names = FALSE)
    cat(paste("  Saved annotated DMRs:", out_name, "\n"))

    # Extract promoter DMRs (most functionally relevant)
    promoter_dmrs <- dmr_anno_df[grepl("Promoter", dmr_anno_df$annotation), ]
    promoter_file <- paste0(contrast_name, "_promoter_DMRs.csv")
    write.csv(promoter_dmrs, promoter_file, row.names = FALSE)
    cat(paste("  Saved promoter DMRs:", nrow(promoter_dmrs), "regions ->", promoter_file, "\n"))

    # Extract gene body DMRs
    genebody_dmrs <- dmr_anno_df[grepl("Exon|Intron", dmr_anno_df$annotation), ]
    cat(paste("  Gene body DMRs:", nrow(genebody_dmrs), "\n"))

    # Extract intergenic DMRs
    intergenic_dmrs <- dmr_anno_df[grepl("Intergenic", dmr_anno_df$annotation), ]
    cat(paste("  Intergenic DMRs:", nrow(intergenic_dmrs), "\n\n"))

    # Summary statistics
    annotation_summary[[contrast_name]] <- data.frame(
        contrast = contrast_name,
        total_dmrs = nrow(dmr_anno_df),
        promoter = nrow(promoter_dmrs),
        genebody = nrow(genebody_dmrs),
        intergenic = nrow(intergenic_dmrs),
        hypermethylated = sum(dmr_anno_df$logFC > 0),
        hypomethylated = sum(dmr_anno_df$logFC < 0)
    )

    # Generate plots
    cat("Generating annotation plots...\n")
    pdf_name <- paste0(contrast_name, "_genomic_distribution.pdf")
    pdf(pdf_name, width = 12, height = 10)

    # Pie chart
    print(plotAnnoPie(dmr_anno, main = paste(contrast_name, "- Genomic Distribution")))

    # Bar chart
    print(plotAnnoBar(dmr_anno, main = paste(contrast_name, "- Feature Distribution")))

    # Distance to TSS
    print(plotDistToTSS(dmr_anno,
                       title = paste(contrast_name, "- Distance to TSS"),
                       xlab = "Distance to TSS (bp)"))

    dev.off()
    cat(paste("  Plots saved:", pdf_name, "\n\n"))
}

################################################################################
# Generate combined summary
################################################################################

cat("==========================================\n")
cat("Generating summary report\n")
cat("==========================================\n\n")

if (length(annotation_summary) > 0) {
    summary_df <- do.call(rbind, annotation_summary)

    # Calculate percentages
    summary_df <- summary_df %>%
        mutate(
            promoter_pct = round(100 * promoter / total_dmrs, 1),
            genebody_pct = round(100 * genebody / total_dmrs, 1),
            intergenic_pct = round(100 * intergenic / total_dmrs, 1),
            hyper_pct = round(100 * hypermethylated / total_dmrs, 1),
            hypo_pct = round(100 * hypomethylated / total_dmrs, 1)
        )

    # Save summary
    write.csv(summary_df, "annotation_summary.csv", row.names = FALSE)
    cat("Summary saved: annotation_summary.csv\n\n")

    # Print summary
    cat("==========================================\n")
    cat("Annotation Summary\n")
    cat("==========================================\n\n")
    print(summary_df)
    cat("\n")

    # Interpretation notes
    cat("==========================================\n")
    cat("Interpretation Guide\n")
    cat("==========================================\n\n")
    cat("Promoter DMRs (±2kb from TSS):\n")
    cat("  - Strongest functional impact on gene expression\n")
    cat("  - Hypermethylation typically = gene silencing\n")
    cat("  - Priority for validation and integration with RNA-seq\n\n")

    cat("Gene body DMRs (exons/introns):\n")
    cat("  - Often associated with active transcription (opposite of promoters)\n")
    cat("  - May regulate alternative splicing\n")
    cat("  - Less predictable effects on expression\n\n")

    cat("Intergenic DMRs:\n")
    cat("  - Potential enhancers or silencers\n")
    cat("  - May regulate distant genes\n")
    cat("  - Require additional evidence for functional assignment\n\n")
}

cat("==========================================\n")
cat("Annotation complete!\n")
cat("==========================================\n")
cat(paste("End:", Sys.time(), "\n\n"))

cat("Output files in:", out_dir, "\n\n")

cat("Key outputs:\n")
cat("1. *_annotated.csv: Full annotation for all DMRs\n")
cat("2. *_promoter_DMRs.csv: Promoter-associated DMRs (for motif analysis)\n")
cat("3. *_genomic_distribution.pdf: Distribution plots\n")
cat("4. annotation_summary.csv: Summary statistics\n\n")

cat("Next steps:\n")
cat("1. Review promoter DMR counts (expect 20-40% of total DMRs)\n")
cat("2. Check if promoter methylation enrichment makes biological sense\n")
cat("3. OPTIONAL: Run motif analysis (step 09) on promoter DMRs\n")
cat("4. ESSENTIAL: Run integration analysis (step 10) to combine with Cut&Tag + RNA-seq\n\n")

cat("Expected patterns for TES project:\n")
cat("- IF TES recruits methylation to TEAD1 sites → expect promoter enrichment\n")
cat("- Compare TES vs GFP promoter DMRs with TES Cut&Tag binding sites\n")
cat("- Genes with promoter hypermethylation should show decreased expression (RNA-seq)\n\n")
