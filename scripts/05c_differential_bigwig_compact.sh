#!/bin/bash
#SBATCH --job-name=05c_diff_bw_compact
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/05c_differential_bigwig_compact.out
#SBATCH --error=logs/05c_differential_bigwig_compact.err

################################################################################
# Script: 05c_differential_bigwig_compact.sh
# Purpose: Generate compact differential methylation BigWig tracks
#
# Description:
#   Same as 05b but with optimized settings for smaller file sizes:
#   - 200bp bin size (vs 50bp) - ~4x smaller files
#   - skipZeroOverZero flag - skips regions with no signal in both tracks
#
# Output:
#   results/05c_differential_bigwig_compact/
#     - TES_vs_GFP_log2ratio.bw   (~50-60MB vs ~240MB)
#     - TES_vs_GFP_subtract.bw    (~50-60MB vs ~240MB)
#     - TES_averaged.bw           (~25-40MB vs ~175MB)
#     - GFP_averaged.bw           (~25-40MB vs ~110MB)
#
# Resolution: 200bp bins - sufficient for DMR-level visualization
################################################################################

set -euo pipefail

echo "========================================="
echo "meDIP-seq - Step 05c: Compact Differential Tracks"
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
OUTPUT_DIR="${BASE_DIR}/results/05c_differential_bigwig_compact"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/individual_comparisons"

# Parameters - OPTIMIZED FOR SMALLER FILES
THREADS=${SLURM_CPUS_PER_TASK:-16}
BIN_SIZE=200        # Larger bins = smaller files (was 50)
PSEUDOCOUNT=1

echo "Parameters (compact mode):"
echo "  Bin size: ${BIN_SIZE} bp (4x larger than standard)"
echo "  Pseudocount: ${PSEUDOCOUNT}"
echo "  Skip zero regions: Yes"
echo "  Threads: ${THREADS}"
echo ""

################################################################################
# Define sample groups
################################################################################

TES_SAMPLES=("TES-1-IP" "TES-2-IP")
GFP_SAMPLES=("GFP-1-IP" "GFP-2-IP")

echo "Sample groups:"
echo "  TES (treatment): ${TES_SAMPLES[*]}"
echo "  GFP (control): ${GFP_SAMPLES[*]}"
echo ""

################################################################################
# Step 1: Verify input files
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
    echo "ERROR: ${MISSING} input file(s) missing!"
    exit 1
fi

echo ""

################################################################################
# Step 2: Create averaged BigWig tracks per condition
################################################################################

echo "========================================="
echo "Step 2: Averaging replicates (compact)"
echo "========================================="

# Average TES replicates
TES_BW_FILES=""
for SAMPLE in "${TES_SAMPLES[@]}"; do
    TES_BW_FILES="${TES_BW_FILES} ${INPUT_DIR}/${SAMPLE}_RPKM.bw"
done

TES_AVERAGED="${OUTPUT_DIR}/TES_averaged.bw"
echo "Creating TES averaged track..."

bigwigAverage \
    --bigwigs ${TES_BW_FILES} \
    --outFileName "${TES_AVERAGED}" \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS}

echo "  Output: ${TES_AVERAGED}"
echo "  Size: $(ls -lh "${TES_AVERAGED}" | awk '{print $5}')"
echo ""

# Average GFP replicates
GFP_BW_FILES=""
for SAMPLE in "${GFP_SAMPLES[@]}"; do
    GFP_BW_FILES="${GFP_BW_FILES} ${INPUT_DIR}/${SAMPLE}_RPKM.bw"
done

GFP_AVERAGED="${OUTPUT_DIR}/GFP_averaged.bw"
echo "Creating GFP averaged track..."

bigwigAverage \
    --bigwigs ${GFP_BW_FILES} \
    --outFileName "${GFP_AVERAGED}" \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS}

echo "  Output: ${GFP_AVERAGED}"
echo "  Size: $(ls -lh "${GFP_AVERAGED}" | awk '{print $5}')"
echo ""

################################################################################
# Step 3: Generate differential tracks
################################################################################

echo "========================================="
echo "Step 3: Generating differential tracks (compact)"
echo "========================================="

# Log2 ratio track
LOG2_OUTPUT="${OUTPUT_DIR}/TES_vs_GFP_log2ratio.bw"
echo "Creating log2 ratio track..."

bigwigCompare \
    --bigwig1 "${TES_AVERAGED}" \
    --bigwig2 "${GFP_AVERAGED}" \
    --outFileName "${LOG2_OUTPUT}" \
    --operation log2 \
    --pseudocount ${PSEUDOCOUNT} \
    --skipZeroOverZero \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS}

echo "  Output: ${LOG2_OUTPUT}"
echo "  Size: $(ls -lh "${LOG2_OUTPUT}" | awk '{print $5}')"
echo ""

# Subtraction track
SUBTRACT_OUTPUT="${OUTPUT_DIR}/TES_vs_GFP_subtract.bw"
echo "Creating subtraction track..."

bigwigCompare \
    --bigwig1 "${TES_AVERAGED}" \
    --bigwig2 "${GFP_AVERAGED}" \
    --outFileName "${SUBTRACT_OUTPUT}" \
    --operation subtract \
    --skipZeroOverZero \
    --binSize ${BIN_SIZE} \
    --numberOfProcessors ${THREADS}

echo "  Output: ${SUBTRACT_OUTPUT}"
echo "  Size: $(ls -lh "${SUBTRACT_OUTPUT}" | awk '{print $5}')"
echo ""

################################################################################
# Step 4: Individual replicate comparisons
################################################################################

echo "========================================="
echo "Step 4: Individual replicate comparisons"
echo "========================================="

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
        --skipZeroOverZero \
        --binSize ${BIN_SIZE} \
        --numberOfProcessors ${THREADS} 2>/dev/null

    echo "    Size: $(ls -lh "${OUT_FILE}" | awk '{print $5}')"
done

echo ""

################################################################################
# Summary
################################################################################

echo "========================================="
echo "Summary - Compact Differential Tracks"
echo "========================================="

echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/*.bw

echo ""
echo "Individual comparisons:"
ls -lh "${OUTPUT_DIR}/individual_comparisons/"*.bw 2>/dev/null || echo "  None"

echo ""
echo "Size comparison (compact vs full resolution):"
echo "  Standard (50bp):  ~750MB total"
echo "  Compact (200bp):  ~$(du -sh "${OUTPUT_DIR}" | cut -f1) total"

echo ""
echo "Note: 200bp resolution is sufficient for:"
echo "  - DMR visualization (500bp windows)"
echo "  - Regional methylation patterns"
echo "  - Genome browser viewing"
echo ""

echo "End time: $(date)"
echo "========================================="
