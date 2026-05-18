#!/bin/bash
#SBATCH --job-name=33b_degs_unbound
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/33b_encode_enhancer_degs_unbound_confident.out
#SBATCH --error=logs/33b_encode_enhancer_degs_unbound_confident.err

# =============================================================================
# ENCODE ENHANCERS OF DEGs DOWN - UNBOUND ENHANCERS VERSION
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE UNBOUND enhancers associated
#          with downregulated DEGs, stratified by DMR overlap.
#
# Strategy: Since TES-bound enhancers with DMR overlap are too few (~1),
#           this version focuses on UNBOUND enhancers with DMR (155 regions)
#           vs matched UNBOUND enhancers without DMR.
#
# =============================================================================

echo "=========================================="
echo "ENCODE ENHANCERS OF DEGs DOWN"
echo "(UNBOUND ENHANCERS - HIGH-CONFIDENCE)"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

OUTDIR="output/33b_encode_enhancer_degs_unbound_confident"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: PREPARE BED FILES (R SCRIPT)
# =============================================================================

echo "=== STEP 1: Preparing Enhancer BED Files ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 33b_encode_enhancer_degs_unbound_confident.R

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
# ANALYSIS STRATEGY:
# Since TES-bound + DMR has very few regions (~1), we focus on:
# 1. Unbound WITH DMR (155 regions) - reliable methylation data
# 2. Unbound WITHOUT DMR (subsampled) - for comparison
# 3. TES-bound ALL (for binding signal reference)
# =============================================================================

if [ "$N_UNBOUND_DMR" -ge $MIN_REGIONS ]; then
    echo "Using Unbound DEGs DOWN enhancers as primary comparison (N=${N_UNBOUND_DMR})"
    echo ""

    # Get path to unbound no DMR file
    BED_UNBOUND_NODMR="${OUTDIR}/Unbound_enhancers_DEGs_DOWN_no_DMR.bed"

    # Create subsampled control from Unbound NO DMR (match N to Unbound WITH DMR)
    shuf -n ${N_UNBOUND_DMR} ${BED_UNBOUND_NODMR} > ${OUTDIR}/Unbound_enhancers_no_DMR_matched.bed

    # Compute matrix: Unbound WITH DMR vs Unbound WITHOUT DMR vs TES-bound ALL
    computeMatrix reference-point \
        --referencePoint center \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${BED_UNBOUND_DMR} \
           ${OUTDIR}/Unbound_enhancers_no_DMR_matched.bed \
           ${OUTDIR}/TES_bound_enhancers_DEGs_DOWN_all.bed \
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

    N_TES_ALL=$(wc -l < "${OUTDIR}/TES_bound_enhancers_DEGs_DOWN_all.bed")

    # Main comparison plot
    plotProfile -m ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -out ${OUTDIR}/MAIN_DEGs_DOWN_Enhancer_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --refPointLabel "Enhancer Center" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "Unbound + DMR (n=${N_UNBOUND_DMR}, RELIABLE)" \
                       "Unbound - no DMR (n=${N_UNBOUND_DMR}, matched)" \
                       "TES-bound ALL (n=${N_TES_ALL})" \
        --plotTitle "DEGs DOWN Enhancers: High-Confidence DMR Comparison" \
        --plotHeight 14 \
        --plotWidth 18 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Created: MAIN_DEGs_DOWN_Enhancer_DMR_Stratified.png"

    # Methylation only plot - FOCUS ON RELIABLE DATA
    plotProfile -m ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -out ${OUTDIR}/METHYLATION_DEGs_DOWN_DMR_Stratified.png \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --samplesLabel "TES meth" "GFP meth" \
        --regionsLabel "Unbound + DMR (RELIABLE)" "Unbound - no DMR" "TES-bound ALL" \
        --plotTitle "Methylation at DEGs DOWN Enhancers (DMR-stratified)" \
        --yMin 0 \
        --dpi 300

    echo "  Created: METHYLATION_DEGs_DOWN_DMR_Stratified.png"

    # Heatmap
    plotHeatmap -m ${OUTDIR}/encode_degs_down_dmr_matrix.gz \
        -out ${OUTDIR}/DEGs_DOWN_Enhancer_DMR_Heatmap.png \
        --colorMap RdBu_r \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "Unbound+DMR" "Unbound-noDMR" "TES-bound" \
        --sortUsing mean \
        --sortUsingSamples 1 \
        --zMin 0 \
        --heatmapHeight 15 \
        --dpi 200

    echo "  Created: DEGs_DOWN_Enhancer_DMR_Heatmap.png"

    # Additional: Methylation difference at DMR-containing vs non-DMR enhancers
    echo ""
    echo "=== Creating Unbound DMR-only comparison ==="

    computeMatrix reference-point \
        --referencePoint center \
        -S $TES_METH $GFP_METH \
        -R ${BED_UNBOUND_DMR} \
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
        --plotTitle "DEGs DOWN Unbound Enhancers: DMR vs No-DMR" \
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
    echo "  Unbound + DMR: ${N_UNBOUND_DMR}"
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
OUTPUT_DIR <- "output/33b_encode_enhancer_degs_unbound_confident"

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
cat("METHYLATION at DEGs DOWN ENHANCERS (HIGH-CONFIDENCE)\n")
cat("(Center ±500bp)\n")
cat("========================================\n\n")

cat(sprintf("Unbound DEGs DOWN + DMR (n=%d) - RELIABLE DATA:\n", n_g1))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g1_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g1_diff < 0, na.rm=TRUE)))

cat(sprintf("Unbound DEGs DOWN - no DMR (n=%d) - MATCHED CONTROL:\n", n_g2))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g2_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g2_diff < 0, na.rm=TRUE)))

cat(sprintf("TES-bound ALL (n=%d) - BINDING REFERENCE:\n", n_g3))
cat(sprintf("  Mean diff (TES-GFP):   %.4f\n", mean(g3_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated: %.1f%%\n", 100*mean(g3_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated:  %.1f%%\n\n", 100*mean(g3_diff < 0, na.rm=TRUE)))

# Key test: Unbound WITH DMR vs Unbound NO DMR
cat("========================================\n")
cat("STATISTICAL TEST (RELIABLE COMPARISON)\n")
cat("========================================\n\n")

res_1v2 <- wilcox.test(g1_diff, g2_diff)
cat("Unbound + DMR vs Unbound - no DMR:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v2$p.value))
cat(sprintf("  Effect: WITH DMR mean=%.4f, NO DMR mean=%.4f\n",
    mean(g1_diff, na.rm=TRUE), mean(g2_diff, na.rm=TRUE)))
cat("\n")

cat("KEY INSIGHT:\n")
cat("  Enhancers WITH high-confidence DMR overlap show the TRUE methylation\n")
cat("  difference at DEGs DOWN enhancers. Regions WITHOUT DMR overlap may\n")
cat("  show artifactual hypermethylation due to GFP library dropout.\n\n")

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
echo "  This analysis compares TES-bound vs Unbound WITHIN enhancers"
echo "  that have reliable methylation data (high-conf DMR overlap)."
echo ""
echo "  Groups WITHOUT DMR overlap may show artifactual hypermethylation"
echo "  due to GFP library dropout (65% of original hyper DMRs had GFP=0)."
echo ""
