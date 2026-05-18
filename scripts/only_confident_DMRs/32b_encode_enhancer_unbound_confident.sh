#!/bin/bash
#SBATCH --job-name=32b_encode_unbound
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/32b_encode_enhancer_unbound_confident.out
#SBATCH --error=logs/32b_encode_enhancer_unbound_confident.err

# =============================================================================
# ENCODE ENHANCER METHYLATION - UNBOUND ENHANCERS VERSION
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE-defined UNBOUND enhancers
#          stratified by high-confidence DMR overlap.
#
# Strategy: Since TES-bound enhancers with DMR overlap are too few (~3),
#           this version focuses on UNBOUND enhancers with DMR (454 regions)
#           vs matched UNBOUND enhancers without DMR.
#
# =============================================================================

echo "=========================================="
echo "ENCODE ENHANCER METHYLATION"
echo "(UNBOUND ENHANCERS - HIGH-CONFIDENCE)"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

OUTDIR="output/32b_encode_enhancer_unbound_confident"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: CHECK/DOWNLOAD ENCODE cCRE ANNOTATIONS
# =============================================================================

echo "=== STEP 1: Checking ENCODE cCRE Annotations ==="
echo ""

ENCODE_SOURCE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1/output/32_encode_enhancer/ENCODE_distal_enhancers.bed"

if [ ! -f "${ENCODE_SOURCE}" ]; then
    echo "ENCODE enhancers not found. Downloading..."

    CCRE_FILE="${OUTDIR}/GRCh38-cCREs.bed"

    if [ ! -f "${CCRE_FILE}" ]; then
        # Download BigBed from UCSC
        wget -q -O "${OUTDIR}/encodeCcreCombined.bb" \
            "http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"

        # Download bigBedToBed if not present
        if [ ! -f "${OUTDIR}/bigBedToBed" ]; then
            wget -q -O "${OUTDIR}/bigBedToBed" \
                "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
            chmod +x "${OUTDIR}/bigBedToBed"
        fi

        # Convert BigBed to BED
        "${OUTDIR}/bigBedToBed" "${OUTDIR}/encodeCcreCombined.bb" "${CCRE_FILE}"
    fi

    # Extract distal enhancers (dELS)
    awk -F'\t' '$11 == "dELS"' "${CCRE_FILE}" > "${OUTDIR}/ENCODE_distal_enhancers.bed"
    ENCODE_SOURCE="${OUTDIR}/ENCODE_distal_enhancers.bed"

    echo "Downloaded and filtered ENCODE dELS enhancers"
else
    echo "Using existing ENCODE enhancers: ${ENCODE_SOURCE}"
    # Copy to our output dir for consistency
    cp "${ENCODE_SOURCE}" "${OUTDIR}/ENCODE_distal_enhancers.bed"
fi

N_ENCODE=$(wc -l < "${OUTDIR}/ENCODE_distal_enhancers.bed")
echo "Total ENCODE dELS enhancers: ${N_ENCODE}"
echo ""

# =============================================================================
# STEP 2: PREPARE ENHANCER SETS BY BINDING + DMR STATUS
# =============================================================================

echo "=== STEP 2: Preparing Enhancer Sets ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 32b_encode_enhancer_unbound_confident.R

# Check if essential BED files were created
# Note: DMR-containing files might not exist if no overlaps found
if [ ! -f "${OUTDIR}/TES_bound_enhancers_no_DMR.bed" ]; then
    echo "ERROR: R script failed to create essential BED files"
    exit 1
fi

# Get counts (with defaults for missing files)
N_TES_DMR=$(wc -l < "${OUTDIR}/TES_bound_enhancers_with_DMR.bed" 2>/dev/null || echo 0)
N_TES_NODMR=$(wc -l < "${OUTDIR}/TES_bound_enhancers_no_DMR.bed" 2>/dev/null || echo 0)
N_UNBOUND_DMR=$(wc -l < "${OUTDIR}/Unbound_enhancers_with_DMR.bed" 2>/dev/null || echo 0)

echo ""
echo "Enhancer counts by binding + DMR status:"
echo "  TES-bound WITH DMR:     ${N_TES_DMR} (MOST RELIABLE)"
echo "  TES-bound NO DMR:       ${N_TES_NODMR}"
echo "  Unbound WITH DMR:       ${N_UNBOUND_DMR} (RELIABLE)"
echo ""

# =============================================================================
# STEP 3: COMPUTE METHYLATION MATRIX
# =============================================================================

echo "=== STEP 3: Computing Methylation Matrix ==="
echo ""

conda activate tg

# BigWig files
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_average.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_average.bw"
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# Verify inputs
for f in "$TES_METH" "$GFP_METH" "$TES_BIND" "$TEAD1_BIND"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

MIN_REGIONS=50

# =============================================================================
# ANALYSIS STRATEGY:
# Since TES-bound + DMR has very few regions (~3), we focus on:
# 1. Unbound WITH DMR (454 regions) - reliable methylation data
# 2. Unbound WITHOUT DMR (subsampled) - for comparison
# 3. TES-bound ALL (for binding signal reference)
# =============================================================================

if [ "$N_UNBOUND_DMR" -ge $MIN_REGIONS ]; then
    echo "Using Unbound enhancers as primary comparison (N=${N_UNBOUND_DMR})"
    echo ""

    # Create subsampled control from Unbound NO DMR (match N to Unbound WITH DMR)
    shuf -n ${N_UNBOUND_DMR} ${OUTDIR}/Unbound_enhancers_no_DMR.bed > ${OUTDIR}/Unbound_enhancers_no_DMR_matched.bed

    # Compute matrix: Unbound WITH DMR vs Unbound WITHOUT DMR vs TES-bound ALL
    computeMatrix reference-point \
        --referencePoint center \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${OUTDIR}/Unbound_enhancers_with_DMR.bed \
           ${OUTDIR}/Unbound_enhancers_no_DMR_matched.bed \
           ${OUTDIR}/TES_bound_enhancers_all.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/encode_enhancer_dmr_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    echo "Matrix created: ${OUTDIR}/encode_enhancer_dmr_matrix.gz"
    echo ""

    # =============================================================================
    # STEP 4: GENERATE PROFILE PLOTS
    # =============================================================================

    echo "=== STEP 4: Generating Profile Plots ==="
    echo ""

    N_TES_ALL=$(wc -l < "${OUTDIR}/TES_bound_enhancers_all.bed")

    # Main comparison plot
    plotProfile -m ${OUTDIR}/encode_enhancer_dmr_matrix.gz \
        -out ${OUTDIR}/MAIN_ENCODE_Enhancer_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --refPointLabel "Enhancer Center" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "Unbound + DMR (n=${N_UNBOUND_DMR}, RELIABLE)" \
                       "Unbound - no DMR (n=${N_UNBOUND_DMR}, matched)" \
                       "TES-bound ALL (n=${N_TES_ALL})" \
        --plotTitle "ENCODE Enhancers: High-Confidence DMR Comparison" \
        --plotHeight 14 \
        --plotWidth 18 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Created: MAIN_ENCODE_Enhancer_DMR_Stratified.png"

    # Methylation only plot - FOCUS ON RELIABLE DATA
    plotProfile -m ${OUTDIR}/encode_enhancer_dmr_matrix.gz \
        -out ${OUTDIR}/METHYLATION_ENCODE_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --samplesLabel "TES meth" "GFP meth" \
        --regionsLabel "Unbound + DMR (RELIABLE)" "Unbound - no DMR" "TES-bound ALL" \
        --plotTitle "Methylation at ENCODE Enhancers (DMR-stratified)" \
        --yMin 0 \
        --dpi 300

    echo "  Created: METHYLATION_ENCODE_DMR_Stratified.png"

    # Heatmap
    plotHeatmap -m ${OUTDIR}/encode_enhancer_dmr_matrix.gz \
        -out ${OUTDIR}/ENCODE_Enhancer_DMR_Heatmap.png \
        --colorMap RdBu_r \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "Unbound+DMR" "Unbound-noDMR" "TES-bound" \
        --sortUsing mean \
        --sortUsingSamples 1 \
        --zMin 0 \
        --dpi 200

    echo "  Created: ENCODE_Enhancer_DMR_Heatmap.png"

    # Additional: Methylation difference at DMR-containing vs non-DMR enhancers
    echo ""
    echo "=== Creating Unbound DMR-only comparison ==="

    computeMatrix reference-point \
        --referencePoint center \
        -S $TES_METH $GFP_METH \
        -R ${OUTDIR}/Unbound_enhancers_with_DMR.bed \
           ${OUTDIR}/Unbound_enhancers_no_DMR_matched.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/unbound_dmr_comparison_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/unbound_dmr_comparison_matrix.gz \
        -out ${OUTDIR}/UNBOUND_DMR_vs_noDMR_profile.png \
        --perGroup \
        --colors "#E31A1C" "#377EB8" \
        --samplesLabel "TES meth" "GFP meth" \
        --regionsLabel "WITH DMR (RELIABLE, n=${N_UNBOUND_DMR})" "NO DMR (matched, n=${N_UNBOUND_DMR})" \
        --plotTitle "Unbound Enhancers: DMR vs No-DMR Comparison" \
        --plotHeight 10 \
        --plotWidth 14 \
        --yMin 0 \
        --dpi 300

    echo "  Created: UNBOUND_DMR_vs_noDMR_profile.png"

    plotHeatmap -m ${OUTDIR}/unbound_dmr_comparison_matrix.gz \
        -out ${OUTDIR}/UNBOUND_DMR_vs_noDMR_heatmap.png \
        --colorMap RdBu_r \
        --samplesLabel "TES meth" "GFP meth" \
        --regionsLabel "WITH DMR" "NO DMR" \
        --sortUsing mean \
        --sortUsingSamples 1 \
        --zMin 0 \
        --heatmapHeight 15 \
        --dpi 200

    echo "  Created: UNBOUND_DMR_vs_noDMR_heatmap.png"
else
    echo "Skipping analysis - insufficient Unbound regions with DMR overlap"
    echo "  Unbound with DMR: ${N_UNBOUND_DMR}"
fi

# =============================================================================
# STEP 5: STATISTICAL QUANTIFICATION
# =============================================================================

echo ""
echo "=== STEP 5: Statistical Quantification ==="
echo ""

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")
OUTPUT_DIR <- "output/32b_encode_enhancer_unbound_confident"

matrix_file <- file.path(OUTPUT_DIR, "encode_enhancer_dmr_matrix.gz")

if (!file.exists(matrix_file)) {
    cat("Matrix file not found - skipping quantification\n")
    quit(status = 0)
}

con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)
header_json <- fromJSON(gsub("^@", "", header_line))

cat("Matrix info:\n")
cat(sprintf("  Samples: %s\n", paste(header_json$sample_labels, collapse=", ")))
cat(sprintf("  Group boundaries: %s\n", paste(header_json$group_boundaries, collapse=", ")))

mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)

# Groups (NEW: Unbound+DMR, Unbound-noDMR, TES-bound ALL)
group_bounds <- header_json$group_boundaries
n_g1 <- group_bounds[2]
n_g2 <- group_bounds[3] - group_bounds[2]
n_g3 <- nrow(mat) - group_bounds[3]

g1_rows <- 1:n_g1  # Unbound + DMR (RELIABLE)
g2_rows <- (n_g1 + 1):(n_g1 + n_g2)  # Unbound - no DMR (matched control)
g3_rows <- (n_g1 + n_g2 + 1):nrow(mat)  # TES-bound ALL

# Center bins
n_bins_per_sample <- as.integer(header_json$upstream[1] + header_json$downstream[1]) / header_json$`bin size`[1]
center_bins <- 90:110

tes_meth_cols <- 6 + center_bins
gfp_meth_cols <- 6 + n_bins_per_sample + center_bins

calc_diff <- function(rows) {
    tes <- rowMeans(as.matrix(mat[rows, ..tes_meth_cols]), na.rm=TRUE)
    gfp <- rowMeans(as.matrix(mat[rows, ..gfp_meth_cols]), na.rm=TRUE)
    return(tes - gfp)
}

g1_diff <- calc_diff(g1_rows)
g2_diff <- calc_diff(g2_rows)
g3_diff <- calc_diff(g3_rows)

cat("\n========================================\n")
cat("METHYLATION at ENCODE ENHANCERS (HIGH-CONFIDENCE)\n")
cat("(Center ±500bp)\n")
cat("========================================\n\n")

cat(sprintf("Unbound WITH DMR (n=%d) - RELIABLE DATA:\n", n_g1))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g1_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g1_diff < 0, na.rm=TRUE)))

cat(sprintf("Unbound NO DMR (n=%d) - MATCHED CONTROL:\n", n_g2))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g2_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g2_diff < 0, na.rm=TRUE)))

cat(sprintf("TES-bound ALL (n=%d) - BINDING REFERENCE:\n", n_g3))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g3_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g3_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g3_diff < 0, na.rm=TRUE)))

# Test: Unbound WITH DMR vs Unbound NO DMR (most reliable comparison)
res_1v2 <- wilcox.test(g1_diff, g2_diff)
cat("========================================\n")
cat("STATISTICAL TEST (RELIABLE COMPARISON)\n")
cat("========================================\n\n")
cat("Unbound + DMR vs Unbound - no DMR:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v2$p.value))
cat(sprintf("  Effect: WITH DMR mean=%.4f, NO DMR mean=%.4f\n",
    mean(g1_diff, na.rm=TRUE), mean(g2_diff, na.rm=TRUE)))
cat("\n")

cat("KEY INSIGHT:\n")
cat("  Enhancers WITH high-confidence DMR overlap show the TRUE methylation\n")
cat("  difference. Enhancers WITHOUT DMR overlap may have artifactual\n")
cat("  hypermethylation due to GFP library dropout.\n\n")

# Save results
results_df <- data.frame(
    Group = c("Unbound_with_DMR", "Unbound_no_DMR", "TES_bound_all"),
    N = c(n_g1, n_g2, n_g3),
    Mean_diff = c(mean(g1_diff, na.rm=TRUE), mean(g2_diff, na.rm=TRUE), mean(g3_diff, na.rm=TRUE)),
    Pct_hypermethylated = c(100*mean(g1_diff > 0, na.rm=TRUE), 100*mean(g2_diff > 0, na.rm=TRUE), 100*mean(g3_diff > 0, na.rm=TRUE)),
    Pct_hypomethylated = c(100*mean(g1_diff < 0, na.rm=TRUE), 100*mean(g2_diff < 0, na.rm=TRUE), 100*mean(g3_diff < 0, na.rm=TRUE)),
    Reliability = c("HIGH (DMR overlap)", "LOW (may have GFP dropout)", "MIXED")
)

write.csv(results_df, file.path(OUTPUT_DIR, "methylation_statistics_dmr.csv"), row.names=FALSE)
cat("Results saved to: methylation_statistics_dmr.csv\n")

RSCRIPT_QUANTIFY

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE (HIGH-CONFIDENCE)"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Key insight:"
echo "  Enhancers with DMR overlap have RELIABLE methylation data."
echo "  Comparing TES-bound vs Unbound WITHIN DMR-overlapping enhancers"
echo "  provides the most accurate assessment of binding-methylation link."
echo ""
