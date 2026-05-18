#!/bin/bash
#SBATCH --job-name=33_degs_enh_conf
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/33_encode_enhancer_degs_down_confident.out
#SBATCH --error=logs/33_encode_enhancer_degs_down_confident.err

# =============================================================================
# ENCODE ENHANCERS OF DEGs DOWN - TES-BOUND FOCUS (HIGH-CONFIDENCE VERSION)
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE enhancers associated with
#          downregulated DEGs, stratified by:
#   1. TES/TEAD1 binding status
#   2. High-confidence DMR overlap (where BOTH samples have >2 reads)
#
# Strategy: Compare TES-bound enhancers WITH DMR vs TES-bound WITHOUT DMR
#           vs Unbound enhancers WITH DMR (subsampled control)
#
# Note: This version may skip analysis if TES-bound + DMR regions are too few.
#       See 33b version for Unbound-focused analysis with more regions.
#
# =============================================================================

echo "=========================================="
echo "ENCODE ENHANCERS OF DEGs DOWN"
echo "(TES-BOUND FOCUS - HIGH-CONFIDENCE)"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

OUTDIR="output/33_encode_enhancer_degs_down_confident"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: PREPARE BED FILES (R SCRIPT)
# =============================================================================

echo "=== STEP 1: Preparing Enhancer BED Files ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run R script (same as 33b - produces all categories)
Rscript 33_encode_enhancer_degs_down_confident.R

# Check if BED files were created
BED_TES_DMR="${OUTDIR}/TES_bound_enhancers_DEGs_DOWN_with_DMR.bed"
BED_TES_NODMR="${OUTDIR}/TES_bound_enhancers_DEGs_DOWN_no_DMR.bed"
BED_UNBOUND_DMR="${OUTDIR}/Unbound_enhancers_DEGs_DOWN_with_DMR.bed"
BED_CONTROL="${OUTDIR}/Control_enhancers.bed"

# Note: DMR files might not exist if no enhancers overlap DMRs - this is OK
if [[ ! -f "$BED_TES_NODMR" && ! -f "$BED_CONTROL" ]]; then
    echo "ERROR: Essential BED files not created. Check R script output."
    exit 1
fi

echo ""
echo "BED files created successfully."
echo ""

# Get counts
N_TES_DMR=$(wc -l < ${BED_TES_DMR} 2>/dev/null || echo 0)
N_TES_NODMR=$(wc -l < ${BED_TES_NODMR} 2>/dev/null || echo 0)
N_UNBOUND_DMR=$(wc -l < ${BED_UNBOUND_DMR} 2>/dev/null || echo 0)
N_CONTROL=$(wc -l < ${BED_CONTROL} 2>/dev/null || echo 0)

echo "Enhancer Counts (DEGs DOWN):"
echo "  TES-bound + DMR:      ${N_TES_DMR} (MOST RELIABLE)"
echo "  TES-bound - no DMR:   ${N_TES_NODMR}"
echo "  Unbound + DMR:        ${N_UNBOUND_DMR} (RELIABLE)"
echo "  Control enhancers:    ${N_CONTROL}"
echo ""

# =============================================================================
# STEP 2: SETUP DEEPTOOLS
# =============================================================================

echo "=== STEP 2: Setting up deepTools ==="

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

# =============================================================================
# STEP 3: COMPUTE METHYLATION MATRIX
# =============================================================================

echo "=== STEP 3: Computing Methylation Matrix ==="
echo ""

MIN_REGIONS=50

# =============================================================================
# ANALYSIS STRATEGY (TES-BOUND FOCUS):
# Compare TES-bound + DMR vs TES-bound - no DMR vs Unbound + DMR (subsampled)
# This tests whether binding affects methylation at enhancers with reliable data
# =============================================================================

if [ "$N_TES_DMR" -ge $MIN_REGIONS ] && [ "$N_UNBOUND_DMR" -ge $MIN_REGIONS ]; then
    echo "Sufficient regions for TES-bound DMR comparison"
    echo ""

    # Create subsampled control from Unbound WITH DMR (match N to TES-bound WITH DMR)
    if [ "$N_UNBOUND_DMR" -gt "$N_TES_DMR" ]; then
        shuf -n ${N_TES_DMR} ${BED_UNBOUND_DMR} > ${OUTDIR}/Unbound_enhancers_DMR_subsampled.bed
    else
        cp ${BED_UNBOUND_DMR} ${OUTDIR}/Unbound_enhancers_DMR_subsampled.bed
    fi

    # Compute matrix: TES-bound WITH DMR vs TES-bound NO DMR vs Unbound WITH DMR
    computeMatrix reference-point \
        --referencePoint center \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${BED_TES_DMR} \
           ${BED_TES_NODMR} \
           ${OUTDIR}/Unbound_enhancers_DMR_subsampled.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    echo "Matrix created: ${OUTDIR}/encode_degs_down_dmr_matrix.gz"
    echo ""

    # =============================================================================
    # STEP 4: GENERATE PROFILE PLOTS
    # =============================================================================

    echo "=== STEP 4: Generating Profile Plots ==="
    echo ""

    N_UNBOUND_SUB=$(wc -l < "${OUTDIR}/Unbound_enhancers_DMR_subsampled.bed")

    # Main comparison plot
    plotProfile -m ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -out ${OUTDIR}/MAIN_DEGs_DOWN_Enhancer_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --refPointLabel "Enhancer Center" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "TES-bound + DMR (n=${N_TES_DMR})" \
                       "TES-bound - no DMR (n=${N_TES_NODMR})" \
                       "Unbound + DMR (n=${N_UNBOUND_SUB})" \
        --plotTitle "DEGs DOWN Enhancers: Methylation (DMR-stratified)" \
        --plotHeight 14 \
        --plotWidth 18 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Created: MAIN_DEGs_DOWN_Enhancer_DMR_Stratified.png"

    # Methylation only plot
    plotProfile -m ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -out ${OUTDIR}/METHYLATION_DEGs_DOWN_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --samplesLabel "TES meth" "GFP meth" \
        --regionsLabel "TES-bound + DMR" "TES-bound - no DMR" "Unbound + DMR" \
        --plotTitle "Methylation at DEGs DOWN Enhancers (DMR-stratified)" \
        --yMin 0 \
        --dpi 300

    echo "  Created: METHYLATION_DEGs_DOWN_DMR_Stratified.png"

    # Heatmap
    plotHeatmap -m ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -out ${OUTDIR}/DEGs_DOWN_Enhancer_DMR_Heatmap.png \
        --colorMap RdBu_r \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "TES+DMR" "TES-noDMR" "Unbound+DMR" \
        --sortUsing mean \
        --sortUsingSamples 1 \
        --zMin 0 \
        --heatmapHeight 15 \
        --dpi 200

    echo "  Created: DEGs_DOWN_Enhancer_DMR_Heatmap.png"

else
    echo "Skipping analysis - insufficient regions with DMR overlap"
    echo "  TES-bound with DMR: ${N_TES_DMR} (need >= ${MIN_REGIONS})"
    echo "  Unbound with DMR: ${N_UNBOUND_DMR} (need >= ${MIN_REGIONS})"
    echo ""
    echo "NOTE: Use 33b_encode_enhancer_degs_unbound_confident.sh for Unbound-focused analysis"
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
OUTPUT_DIR <- "output/33_encode_enhancer_degs_down_confident"

matrix_file <- file.path(OUTPUT_DIR, "encode_degs_down_dmr_matrix.gz")

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

# Groups: TES-bound+DMR, TES-bound-noDMR, Unbound+DMR
group_bounds <- header_json$group_boundaries
n_g1 <- group_bounds[2]
n_g2 <- group_bounds[3] - group_bounds[2]
n_g3 <- nrow(mat) - group_bounds[3]

g1_rows <- 1:n_g1
g2_rows <- (n_g1 + 1):(n_g1 + n_g2)
g3_rows <- (n_g1 + n_g2 + 1):nrow(mat)

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
cat("METHYLATION at DEGs DOWN ENHANCERS (HIGH-CONFIDENCE)\n")
cat("(Center ±500bp)\n")
cat("========================================\n\n")

cat(sprintf("TES-bound WITH DMR (n=%d) - MOST RELIABLE:\n", n_g1))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g1_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g1_diff < 0, na.rm=TRUE)))

cat(sprintf("TES-bound NO DMR (n=%d) - may have GFP dropout:\n", n_g2))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g2_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g2_diff < 0, na.rm=TRUE)))

cat(sprintf("Unbound WITH DMR (n=%d) - RELIABLE:\n", n_g3))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g3_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g3_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g3_diff < 0, na.rm=TRUE)))

res_1v3 <- wilcox.test(g1_diff, g3_diff)
cat("========================================\n")
cat("STATISTICAL TEST (RELIABLE COMPARISON)\n")
cat("========================================\n\n")
cat("TES-bound + DMR vs Unbound + DMR:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v3$p.value))
if (mean(g1_diff, na.rm=TRUE) > mean(g3_diff, na.rm=TRUE)) {
    cat("  Direction: TES binding causes MORE methylation\n\n")
} else {
    cat("  Direction: TES binding causes LESS methylation\n\n")
}

results_df <- data.frame(
    Group = c("TES_bound_with_DMR", "TES_bound_no_DMR", "Unbound_with_DMR"),
    N = c(n_g1, n_g2, n_g3),
    Mean_diff = c(mean(g1_diff, na.rm=TRUE), mean(g2_diff, na.rm=TRUE), mean(g3_diff, na.rm=TRUE)),
    Pct_hypermethylated = c(100*mean(g1_diff > 0, na.rm=TRUE), 100*mean(g2_diff > 0, na.rm=TRUE), 100*mean(g3_diff > 0, na.rm=TRUE)),
    Pct_hypomethylated = c(100*mean(g1_diff < 0, na.rm=TRUE), 100*mean(g2_diff < 0, na.rm=TRUE), 100*mean(g3_diff < 0, na.rm=TRUE)),
    Reliability = c("HIGH", "LOW (GFP dropout possible)", "HIGH")
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

