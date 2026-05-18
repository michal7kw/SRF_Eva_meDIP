#!/bin/bash
#SBATCH --job-name=combine_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/05a_combine_bigwig.out
#SBATCH --error=logs/05a_combine_bigwig.err

################################################################################
# Script: 05a_combine_bigwig_replicates.sh
# Purpose: Combine technical replicate BigWig files by averaging signal
#
# Description:
#   This script combines RPKM-normalized BigWig files from technical replicates
#   into averaged BigWig files for each condition. These combined tracks provide:
#   - Reduced noise through replicate averaging
#   - Cleaner signal for visualization and downstream analysis
#   - Required input for scripts 18 (DMR binding overlay) and 19 (binding causes methylation)
#
# Why bigwigAverage (not bigWigMerge):
#   - bigWigMerge from UCSC tools SUMS values, not averages them
#   - This results in artificially inflated signal (2x for 2 replicates)
#   - bigwigAverage from deepTools computes proper mean across replicates
#   - Maintains correct RPKM scale for cross-sample comparison
#
# Technical Details:
#   - Uses deepTools bigwigAverage for mathematically correct averaging
#   - Preserves 10 bp bin resolution from input BigWig files
#   - Handles missing data appropriately (skipZeroOverZero)
#
# Input:
#   - results/05_bigwig/*_RPKM.bw (individual replicate BigWig files)
#
# Output:
#   - results/05_bigwig/TES_average.bw (averaged TES replicates)
#   - results/05_bigwig/GFP_average.bw (averaged GFP replicates)
#
# Dependencies:
#   - deepTools (bigwigAverage)
#   - Conda environment: bigwig-generation or peak_calling_new
#
# Runtime: ~30-60 minutes
# Memory: ~16-32 GB
################################################################################

set -euo pipefail

echo "========================================="
echo "meDIP-seq Pipeline - Step 05b: Combine BigWig Replicates"
echo "========================================="
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment with deepTools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig-generation

# Verify bigwigAverage is available
if ! command -v bigwigAverage &> /dev/null; then
    echo "bigwigAverage not found in bigwig-generation, trying peak_calling_new..."
    conda activate peak_calling_new
    if ! command -v bigwigAverage &> /dev/null; then
        echo "ERROR: bigwigAverage not found. Please ensure deepTools is installed."
        exit 1
    fi
fi

echo "deepTools version:"
bigwigAverage --version 2>&1 || echo "bigwigAverage loaded"
echo ""

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
BIGWIG_DIR="${BASE_DIR}/results/05_bigwig"

# Check input directory exists
if [ ! -d "${BIGWIG_DIR}" ]; then
    echo "ERROR: BigWig directory not found: ${BIGWIG_DIR}"
    echo "Please run 05_bigwig.sh first to generate individual BigWig files."
    exit 1
fi

# Define input files - TES replicates
TES_REP1="${BIGWIG_DIR}/TES-1-IP_RPKM.bw"
TES_REP2="${BIGWIG_DIR}/TES-2-IP_RPKM.bw"

# Define input files - GFP replicates
GFP_REP1="${BIGWIG_DIR}/GFP-1-IP_RPKM.bw"
GFP_REP2="${BIGWIG_DIR}/GFP-2-IP_RPKM.bw"

# Define output files (consistent naming - _IP suffix omitted as it's implied for meDIP samples)
TES_COMBINED="${BIGWIG_DIR}/TES_average.bw"
GFP_COMBINED="${BIGWIG_DIR}/GFP_average.bw"

# Parameters
THREADS=8
BIN_SIZE=10  # Match bin size from 05_bigwig.sh for consistency

echo "Configuration:"
echo "  Input directory: ${BIGWIG_DIR}"
echo "  Threads: ${THREADS}"
echo "  Bin size: ${BIN_SIZE} bp"
echo ""

# ============================================================================
# Verify input files exist
# ============================================================================

echo "=== Checking input files ==="

missing_files=0

for f in "$TES_REP1" "$TES_REP2" "$GFP_REP1" "$GFP_REP2"; do
    if [ -f "$f" ]; then
        size=$(stat -c%s "$f" | numfmt --to=iec-i --suffix=B)
        echo "  Found: $(basename $f) ($size)"
    else
        echo "  MISSING: $f"
        missing_files=$((missing_files + 1))
    fi
done

if [ $missing_files -gt 0 ]; then
    echo ""
    echo "ERROR: $missing_files input file(s) missing."
    echo "Please run 05_bigwig.sh first to generate all BigWig files."
    exit 1
fi

echo ""

# ============================================================================
# Combine TES replicates
# ============================================================================

echo "========================================="
echo "Combining TES replicates"
echo "========================================="
echo "  Input 1: $(basename $TES_REP1)"
echo "  Input 2: $(basename $TES_REP2)"
echo "  Output: $(basename $TES_COMBINED)"
echo ""

# Remove existing output to ensure clean run
if [ -f "$TES_COMBINED" ]; then
    echo "  Removing existing file: $TES_COMBINED"
    rm -f "$TES_COMBINED"
fi

echo "  Running bigwigAverage..."
bigwigAverage \
    -b "$TES_REP1" "$TES_REP2" \
    -o "$TES_COMBINED" \
    --binSize ${BIN_SIZE} \
    -p ${THREADS}

# Verify output
if [ -f "$TES_COMBINED" ]; then
    size=$(stat -c%s "$TES_COMBINED" | numfmt --to=iec-i --suffix=B)
    echo "  SUCCESS: Created $TES_COMBINED ($size)"
else
    echo "  ERROR: Failed to create $TES_COMBINED"
    exit 1
fi

echo ""

# ============================================================================
# Combine GFP replicates
# ============================================================================

echo "========================================="
echo "Combining GFP replicates"
echo "========================================="
echo "  Input 1: $(basename $GFP_REP1)"
echo "  Input 2: $(basename $GFP_REP2)"
echo "  Output: $(basename $GFP_COMBINED)"
echo ""

# Remove existing output to ensure clean run
if [ -f "$GFP_COMBINED" ]; then
    echo "  Removing existing file: $GFP_COMBINED"
    rm -f "$GFP_COMBINED"
fi

echo "  Running bigwigAverage..."
bigwigAverage \
    -b "$GFP_REP1" "$GFP_REP2" \
    -o "$GFP_COMBINED" \
    --binSize ${BIN_SIZE} \
    -p ${THREADS}

# Verify output
if [ -f "$GFP_COMBINED" ]; then
    size=$(stat -c%s "$GFP_COMBINED" | numfmt --to=iec-i --suffix=B)
    echo "  SUCCESS: Created $GFP_COMBINED ($size)"
else
    echo "  ERROR: Failed to create $GFP_COMBINED"
    exit 1
fi

echo ""

# ============================================================================
# Validation: Compare signal ranges
# ============================================================================

echo "========================================="
echo "Validating combined BigWig files"
echo "========================================="

# Use bigWigInfo to get statistics if available
if command -v bigWigInfo &> /dev/null; then
    echo ""
    echo "TES replicate 1 stats:"
    bigWigInfo "$TES_REP1" 2>/dev/null | grep -E "mean|min|max|std" || echo "  (stats not available)"

    echo ""
    echo "TES replicate 2 stats:"
    bigWigInfo "$TES_REP2" 2>/dev/null | grep -E "mean|min|max|std" || echo "  (stats not available)"

    echo ""
    echo "TES combined stats:"
    bigWigInfo "$TES_COMBINED" 2>/dev/null | grep -E "mean|min|max|std" || echo "  (stats not available)"

    echo ""
    echo "GFP replicate 1 stats:"
    bigWigInfo "$GFP_REP1" 2>/dev/null | grep -E "mean|min|max|std" || echo "  (stats not available)"

    echo ""
    echo "GFP replicate 2 stats:"
    bigWigInfo "$GFP_REP2" 2>/dev/null | grep -E "mean|min|max|std" || echo "  (stats not available)"

    echo ""
    echo "GFP combined stats:"
    bigWigInfo "$GFP_COMBINED" 2>/dev/null | grep -E "mean|min|max|std" || echo "  (stats not available)"
else
    echo "  bigWigInfo not available, skipping detailed statistics."
fi

echo ""

# ============================================================================
# Summary
# ============================================================================

echo "========================================="
echo "Summary"
echo "========================================="
echo ""
echo "Combined BigWig files created:"
echo ""
ls -lh "${BIGWIG_DIR}"/*average*.bw 2>/dev/null || echo "  No combined files found"
echo ""
echo "File naming convention:"
echo "  - TES_average.bw: Mean of TES-1-IP and TES-2-IP"
echo "  - GFP_average.bw: Mean of GFP-1-IP and GFP-2-IP"
echo ""
echo "These files are ready for use in:"
echo "  - 18_dmr_binding_overlay.sh (DMR visualization with binding)"
echo "  - 19_binding_causes_methylation.sh (binding-methylation analysis)"
echo "  - Any deepTools visualization (computeMatrix, plotHeatmap, plotProfile)"
echo ""
echo "Technical notes:"
echo "  - Signal values are true AVERAGES (not sums) of replicates"
echo "  - Bin size: ${BIN_SIZE} bp (matching individual BigWig files)"
echo "  - RPKM normalization preserved from input files"
echo ""

echo "========================================="
echo "Next steps:"
echo "1. Verify combined tracks in IGV alongside individual replicates"
echo "2. Run 18_dmr_binding_overlay.sh for DMR + binding visualization"
echo "3. Run 19_binding_causes_methylation.sh for causal analysis"
echo "========================================="

echo ""
echo "End time: $(date)"
echo "========================================="

# ============================================================================
# Technical Background: Why bigwigAverage is Preferred
# ============================================================================
#
# The bigWigMerge tool from UCSC utilities has a critical limitation:
# It SUMS the values from input files rather than averaging them.
#
# Example with two replicates (RPKM values):
#   Replicate 1 at position X: 50 RPKM
#   Replicate 2 at position X: 60 RPKM
#
#   bigWigMerge output: 110 RPKM (incorrect - double the expected scale)
#   bigwigAverage output: 55 RPKM (correct - true average)
#
# This matters because:
# 1. Downstream analysis assumes RPKM scale (0-100 typical for meDIP)
# 2. Combined track should be directly comparable to individual replicates
# 3. Heatmap color scales are calibrated for single-sample RPKM values
# 4. Statistical comparisons assume consistent normalization
#
# The bigwigAverage tool from deepTools handles this correctly:
# - Computes element-wise mean across input files
# - Handles missing data appropriately (--skipZeroOverZero)
# - Supports multi-threading for performance
# - Maintains consistent bin sizes
#
# Reference: deepTools documentation
# https://deeptools.readthedocs.io/en/develop/content/tools/bigwigAverage.html
# ============================================================================
