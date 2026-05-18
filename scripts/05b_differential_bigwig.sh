#!/bin/bash
#SBATCH --job-name=05b_differential_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/05b_differential_bigwig.out
#SBATCH --error=logs/05b_differential_bigwig.err

################################################################################
# Script: 05b_differential_bigwig.sh
# Purpose: Generate differential methylation BigWig tracks
#
# Description:
#   Creates BigWig files showing the difference in methylation signal between
#   conditions (TES vs GFP). These tracks visualize where methylation is
#   increased (hypermethylated) or decreased (hypomethylated) in TES samples.
#
# Output Types:
#   1. Log2 Ratio Tracks:
#      - Positive values = hypermethylated in TES (increased methylation)
#      - Negative values = hypomethylated in TES (decreased methylation)
#      - Zero = no change
#      - Best for: Fold-change visualization, identifying strong differences
#
#   2. Subtraction Tracks:
#      - Positive values = higher signal in TES
#      - Negative values = higher signal in GFP
#      - Best for: Absolute signal differences, when log2 would be extreme
#
# Strategy:
#   - Average replicates within each condition using bigwigAverage
#   - Then compare averaged tracks using bigwigCompare
#   - This reduces noise and provides robust differential signal
#
# Input:
#   - RPKM-normalized BigWig files from 05_bigwig.sh
#
# Output:
#   - TES_vs_GFP_log2ratio.bw   (log2 of TES/GFP)
#   - TES_vs_GFP_subtract.bw    (TES - GFP)
#   - TES_averaged.bw           (mean of TES replicates)
#   - GFP_averaged.bw           (mean of GFP replicates)
#
# Visualization:
#   In IGV, load differential tracks with diverging color scale:
#   - Blue-White-Red: Blue=hypomethylated, White=no change, Red=hypermethylated
#
# Runtime: ~30-60 minutes
# Memory: ~16-32 GB
################################################################################

set -euo pipefail

echo "========================================="
echo "meDIP-seq Pipeline - Step 05b: Differential Methylation Tracks"
echo "========================================="
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig-generation

echo "deepTools version:"
bigwigCompare --version 2>&1 | head -1 || echo "bigwigCompare loaded"
echo ""

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
INPUT_DIR="${BASE_DIR}/results/05_bigwig"
OUTPUT_DIR="${BASE_DIR}/results/05b_differential_bigwig"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/individual_comparisons"

# Parameters
THREADS=${SLURM_CPUS_PER_TASK:-16}
BIN_SIZE=50  # Slightly larger bins for smoother differential tracks
PSEUDOCOUNT=1  # Add pseudocount to avoid log2(0) issues

echo "Parameters:"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Pseudocount: ${PSEUDOCOUNT}"
echo "  Threads: ${THREADS}"
echo ""

################################################################################
# Define sample groups
################################################################################

# TES group (treatment)
TES_SAMPLES=("TES-1-IP" "TES-2-IP")

# GFP group (control)
GFP_SAMPLES=("GFP-1-IP" "GFP-2-IP")

echo "Sample groups:"
echo "  TES (treatment): ${TES_SAMPLES[*]}"
echo "  GFP (control): ${GFP_SAMPLES[*]}"
echo ""

################################################################################
# Step 1: Verify all input BigWig files exist
################################################################################

echo "========================================="
echo "Step 1: Checking input files"
echo "========================================="

ALL_SAMPLES=("${TES_SAMPLES[@]}" "${GFP_SAMPLES[@]}")
MISSING=0

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    BW_FILE="${INPUT_DIR}/${SAMPLE}_RPKM.bw"
    if [ -f "${BW_FILE}" ]; then
        echo "  [OK] ${SAMPLE}_RPKM.bw"
    else
        echo "  [MISSING] ${BW_FILE}"
        MISSING=$((MISSING + 1))
    fi
done

if [ ${MISSING} -gt 0 ]; then
    echo ""
    echo "ERROR: ${MISSING} input file(s) missing!"
    echo "Please run 05_bigwig.sh first to generate all BigWig files."
    exit 1
fi

echo ""
echo "All input files present."
echo ""

################################################################################
# Step 2: Create averaged BigWig tracks per condition
################################################################################

echo "========================================="
echo "Step 2: Averaging replicates within conditions"
echo "========================================="

# Average TES replicates
TES_BW_FILES=""
for SAMPLE in "${TES_SAMPLES[@]}"; do
    TES_BW_FILES="${TES_BW_FILES} ${INPUT_DIR}/${SAMPLE}_RPKM.bw"
done

TES_AVERAGED="${OUTPUT_DIR}/TES_averaged.bw"
echo "Creating TES averaged track..."
echo "  Input files: ${TES_BW_FILES}"

bigwigAverage \
    --bigwigs ${TES_BW_FILES} \
    --outFileName "${TES_AVERAGED}" \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS} \
    --verbose

echo "  Output: ${TES_AVERAGED}"
echo ""

# Average GFP replicates
GFP_BW_FILES=""
for SAMPLE in "${GFP_SAMPLES[@]}"; do
    GFP_BW_FILES="${GFP_BW_FILES} ${INPUT_DIR}/${SAMPLE}_RPKM.bw"
done

GFP_AVERAGED="${OUTPUT_DIR}/GFP_averaged.bw"
echo "Creating GFP averaged track..."
echo "  Input files: ${GFP_BW_FILES}"

bigwigAverage \
    --bigwigs ${GFP_BW_FILES} \
    --outFileName "${GFP_AVERAGED}" \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS} \
    --verbose

echo "  Output: ${GFP_AVERAGED}"
echo ""

################################################################################
# Step 3: Generate differential tracks (TES vs GFP)
################################################################################

echo "========================================="
echo "Step 3: Generating differential methylation tracks"
echo "========================================="

# Log2 ratio track (TES / GFP)
LOG2_OUTPUT="${OUTPUT_DIR}/TES_vs_GFP_log2ratio.bw"
echo "Creating log2 ratio track (TES/GFP)..."
echo "  Interpretation: Positive = hypermethylated in TES"

bigwigCompare \
    --bigwig1 "${TES_AVERAGED}" \
    --bigwig2 "${GFP_AVERAGED}" \
    --outFileName "${LOG2_OUTPUT}" \
    --operation log2 \
    --pseudocount ${PSEUDOCOUNT} \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS}

echo "  Output: ${LOG2_OUTPUT}"
echo ""

# Subtraction track (TES - GFP)
SUBTRACT_OUTPUT="${OUTPUT_DIR}/TES_vs_GFP_subtract.bw"
echo "Creating subtraction track (TES - GFP)..."
echo "  Interpretation: Positive = higher signal in TES"

bigwigCompare \
    --bigwig1 "${TES_AVERAGED}" \
    --bigwig2 "${GFP_AVERAGED}" \
    --outFileName "${SUBTRACT_OUTPUT}" \
    --operation subtract \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS}

echo "  Output: ${SUBTRACT_OUTPUT}"
echo ""

################################################################################
# Step 4: Generate individual replicate comparisons (optional, for QC)
################################################################################

echo "========================================="
echo "Step 4: Individual replicate comparisons (for QC)"
echo "========================================="

# Compare each TES replicate to GFP average
for TES_REP in "${TES_SAMPLES[@]}"; do
    TES_BW="${INPUT_DIR}/${TES_REP}_RPKM.bw"
    OUT_FILE="${OUTPUT_DIR}/individual_comparisons/${TES_REP}_vs_GFP_log2ratio.bw"

    echo "  ${TES_REP} vs GFP..."

    bigwigCompare \
        --bigwig1 "${TES_BW}" \
        --bigwig2 "${GFP_AVERAGED}" \
        --outFileName "${OUT_FILE}" \
        --operation log2 \
        --pseudocount ${PSEUDOCOUNT} \
        --binSize ${BIN_SIZE} \
        --numberOfProcessors ${THREADS} 2>/dev/null

    echo "    Output: ${OUT_FILE}"
done

echo ""

################################################################################
# Step 5: Generate summary statistics
################################################################################

echo "========================================="
echo "Step 5: Track statistics"
echo "========================================="

STATS_FILE="${OUTPUT_DIR}/differential_track_stats.txt"

cat > "${STATS_FILE}" << EOF
Differential Methylation Track Statistics
==========================================
Generated: $(date)
Script: 05b_differential_bigwig.sh

Sample Groups:
  TES (treatment): ${TES_SAMPLES[*]}
  GFP (control): ${GFP_SAMPLES[*]}

Parameters:
  Bin size: ${BIN_SIZE} bp
  Pseudocount (for log2): ${PSEUDOCOUNT}
  Tool: deepTools bigwigCompare / bigwigAverage

Output Files:
  1. TES_averaged.bw - Mean of TES replicates
  2. GFP_averaged.bw - Mean of GFP replicates
  3. TES_vs_GFP_log2ratio.bw - log2(TES/GFP) differential track
  4. TES_vs_GFP_subtract.bw - (TES - GFP) differential track

Interpretation Guide:

  Log2 Ratio Track:
    +2 = 4-fold hypermethylated in TES
    +1 = 2-fold hypermethylated in TES
     0 = No change
    -1 = 2-fold hypomethylated in TES
    -2 = 4-fold hypomethylated in TES

  Subtraction Track:
    Positive values = Higher RPKM signal in TES
    Negative values = Higher RPKM signal in GFP
    Units: RPKM difference

File Sizes:
EOF

# Add file sizes
for FILE in "${OUTPUT_DIR}"/*.bw; do
    if [ -f "${FILE}" ]; then
        SIZE=$(ls -lh "${FILE}" | awk '{print $5}')
        NAME=$(basename "${FILE}")
        echo "  ${NAME}: ${SIZE}" >> "${STATS_FILE}"
    fi
done

cat "${STATS_FILE}"
echo ""

################################################################################
# Summary
################################################################################

echo "========================================="
echo "Differential Methylation Tracks Summary"
echo "========================================="

echo ""
echo "Generated tracks:"
ls -lh "${OUTPUT_DIR}"/*.bw
echo ""

echo "Individual replicate comparisons (QC):"
ls -lh "${OUTPUT_DIR}/individual_comparisons/"*.bw 2>/dev/null || echo "  None generated"
echo ""

echo "========================================="
echo "Visualization Instructions"
echo "========================================="
echo ""
echo "Loading in IGV:"
echo "  1. File > Load from File > Select differential .bw files"
echo "  2. Right-click track > Set Data Range"
echo "     - Log2 ratio: -2 to +2 (captures most changes)"
echo "     - Subtraction: -50 to +50 (adjust based on signal)"
echo ""
echo "Recommended Color Scheme (in IGV):"
echo "  - Right-click track > Change Track Color"
echo "  - Use diverging scale: Blue (neg) - White (0) - Red (pos)"
echo "  - Or load as 'Heatmap' track type"
echo ""
echo "What to look for:"
echo "  1. Regions with consistent positive signal = TES hypermethylation"
echo "  2. Regions with consistent negative signal = TES hypomethylation"
echo "  3. Compare with DMRs from 07_differential_methylation_MEDIPS.R"
echo "  4. Navigate to DEG promoters to see methylation changes"
echo ""
echo "Quality Control:"
echo "  - Compare individual replicate tracks in individual_comparisons/"
echo "  - Replicates should show similar patterns"
echo "  - Check consistency at known DMRs"
echo ""

echo "========================================="
echo "Next steps:"
echo "  1. Load differential tracks in IGV alongside:"
echo "     - Individual sample BigWigs"
echo "     - DMR BED files from step 07"
echo "     - Cut&Tag binding peaks"
echo "  2. Validate DMR calls visually"
echo "  3. Navigate to candidate genes"
echo "========================================="

echo ""
echo "End time: $(date)"
echo "========================================="

################################################################################
# Understanding Differential Tracks:
#
# Why Average Replicates First?
#   - Reduces biological and technical noise
#   - Provides more stable denominator for log2 ratio
#   - Single differential track is easier to interpret than multiple
#
# Log2 Ratio vs Subtraction:
#   - Log2: Better for fold-change interpretation, normalized scale
#     * Symmetric around 0 (2-fold up = +1, 2-fold down = -1)
#     * Issues: Can be extreme when denominator is low
#   - Subtract: Better for absolute differences
#     * Not affected by low-signal denominator
#     * Harder to compare across different baseline levels
#
# Pseudocount Purpose:
#   - Prevents log2(x/0) = infinity
#   - Value of 1 is standard (adds 1 RPKM to both tracks)
#   - Larger pseudocount = more conservative (smaller log2 ratios)
#
# Bin Size Choice:
#   - 50bp: Smoother tracks, faster computation
#   - Smaller bins (10bp): More detail but noisier
#   - Match to DMR calling window size for consistency
#
# Visualization Tips:
#   - Always view differential tracks alongside original signals
#   - Artifacts: High differential at very low signal regions
#   - Use consistent scale across all differential tracks
################################################################################
