#!/bin/bash
#SBATCH --job-name=31_enhancer_confident
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/31_enhancer_analysis_confident.out
#SBATCH --error=logs/31_enhancer_analysis_confident.err

# =============================================================================
# ENHANCER LEVEL ANALYSIS - HIGH-CONFIDENCE VERSION
# =============================================================================
#
# Purpose: Analyze Binding and Methylation profiles at Enhancer regions.
#          Stratified by target gene expression AND high-confidence DMR status.
#
# Key insight: Only enhancers overlapping high-conf DMRs have reliable
# methylation data (both TES and GFP have >2 reads).
#
# =============================================================================

echo "=========================================="
echo "ENHANCER LEVEL ANALYSIS"
echo "(HIGH-CONFIDENCE VERSION)"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

OUTDIR="output/31_enhancer_analysis_confident"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: CREATE BED FILES (R SCRIPT)
# =============================================================================

echo "=== STEP 1: Creating Enhancer BED files ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 31_enhancer_analysis_confident.R

# Check if BED files were created
BED_DMR="${OUTDIR}/Enhancers_DEGs_DOWN_with_DMR.bed"
BED_NODMR="${OUTDIR}/Enhancers_DEGs_DOWN_no_DMR.bed"
BED_CONTROL="${OUTDIR}/Enhancers_Random_Control.bed"

# Note: DMR file might not exist if no enhancers overlap DMRs - this is OK
if [[ ! -f "$BED_NODMR" && ! -f "$BED_CONTROL" ]]; then
    echo "ERROR: Essential BED files not created. Check R script output."
    exit 1
fi

echo ""
echo "BED files created successfully."
echo ""

# Get counts
N_DMR=$(wc -l < ${BED_DMR} 2>/dev/null || echo 0)
N_NODMR=$(wc -l < ${BED_NODMR} 2>/dev/null || echo 0)
N_CONTROL=$(wc -l < ${BED_CONTROL})

echo "Enhancer Counts:"
echo "  DEGs DOWN with DMR:    ${N_DMR} (MOST RELIABLE)"
echo "  DEGs DOWN no DMR:      ${N_NODMR}"
echo "  Control Enhancers:     ${N_CONTROL}"
echo ""

# =============================================================================
# STEP 2: SETUP DEEPTOOLS
# =============================================================================

echo "=== STEP 2: Setting up deepTools ==="

conda activate tg

# BigWig files
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_average.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_average.bw"

# Verify inputs
for f in "$TES_BIND" "$TEAD1_BIND" "$TES_METH" "$GFP_METH"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

# =============================================================================
# STEP 3: COMPUTE MATRIX (DMR-STRATIFIED)
# =============================================================================

echo "=== STEP 3: Computing Matrix ==="

MIN_REGIONS=10

if [ "$N_DMR" -ge $MIN_REGIONS ] && [ "$N_NODMR" -ge $MIN_REGIONS ]; then
    computeMatrix reference-point \
        --referencePoint center \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${BED_DMR} ${BED_NODMR} ${BED_CONTROL} \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/enhancer_dmr_stratified_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    # =============================================================================
    # STEP 4: PLOT PROFILES
    # =============================================================================

    echo "=== STEP 4: Plotting Profiles ==="

    # Main Comparison (All Signals)
    plotProfile -m ${OUTDIR}/enhancer_dmr_stratified_matrix.gz \
        -out ${OUTDIR}/MAIN_Enhancer_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --refPointLabel "Peak Center" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "DEGs DOWN + DMR (n=${N_DMR})" \
                       "DEGs DOWN - no DMR (n=${N_NODMR})" \
                       "Control (n=${N_CONTROL})" \
        --plotTitle "Enhancer Profile: DEGs DOWN (DMR-stratified)" \
        --plotHeight 14 \
        --plotWidth 18 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Created: MAIN_Enhancer_DMR_Stratified.png"

    # Methylation Only
    plotProfile -m ${OUTDIR}/enhancer_dmr_stratified_matrix.gz \
        -out ${OUTDIR}/METHYLATION_Enhancer_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --samplesLabel "TES meth" "GFP meth" \
        --regionsLabel "DEGs DOWN + DMR" "DEGs DOWN - no DMR" "Control" \
        --plotTitle "Methylation at Enhancers (DMR-stratified)" \
        --yMin 0 \
        --dpi 300

    echo "  Created: METHYLATION_Enhancer_DMR_Stratified.png"
else
    echo "  Skipping DMR-stratified analysis (insufficient regions)"
    echo "  DEGs DOWN with DMR: ${N_DMR}, no DMR: ${N_NODMR}"
fi

# =============================================================================
# STEP 5: QUANTIFY METHYLATION DIFFERENCES
# =============================================================================

echo "=== STEP 5: Quantifying Methylation via R ==="

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs")
OUTPUT_DIR <- "output/31_enhancer_analysis_confident"

matrix_file <- file.path(OUTPUT_DIR, "enhancer_dmr_stratified_matrix.gz")

if (!file.exists(matrix_file)) {
    cat("Matrix file not found - skipping quantification\n")
    quit(status = 0)
}

con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)
header_json <- fromJSON(gsub("^@", "", header_line))

mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)

# Groups
group_bounds <- header_json$group_boundaries
n_g1 <- group_bounds[2]
n_g2 <- group_bounds[3] - group_bounds[2]
n_g3 <- nrow(mat) - group_bounds[3]

g1_rows <- 1:n_g1
g2_rows <- (n_g1 + 1):(n_g1 + n_g2)
g3_rows <- (n_g1 + n_g2 + 1):nrow(mat)

# Center bins
n_bins_per_sample <- (header_json$upstream + header_json$downstream) / header_json$`bin size`
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
cat("METHYLATION DIFFERENCE at ENHANCERS (HIGH-CONFIDENCE)\n")
cat("========================================\n\n")

cat(sprintf("DEGs DOWN + DMR (n=%d):\n", n_g1))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g1_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g1_diff < 0, na.rm=TRUE)))

cat(sprintf("DEGs DOWN - no DMR (n=%d):\n", n_g2))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g2_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g2_diff < 0, na.rm=TRUE)))

cat(sprintf("Control (n=%d):\n", n_g3))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g3_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g3_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g3_diff < 0, na.rm=TRUE)))

# Test
res <- wilcox.test(g1_diff, g3_diff)
cat(sprintf("Wilcoxon test (DEGs DOWN + DMR vs Control): p=%.4e\n", res$p.value))

RSCRIPT_QUANTIFY

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE (HIGH-CONFIDENCE)"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Key insight:"
echo "  Enhancers WITH high-confidence DMRs have RELIABLE methylation"
echo "  (both TES and GFP have >2 reads at overlapping regions)"
echo ""
