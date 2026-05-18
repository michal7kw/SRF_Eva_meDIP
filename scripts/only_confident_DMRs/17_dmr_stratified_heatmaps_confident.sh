#!/bin/bash
#SBATCH --job-name=dmr_heatmaps_confident
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/17_dmr_stratified_heatmaps_confident.out
#SBATCH --error=logs/17_dmr_stratified_heatmaps_confident.err

# ============================================================================
# DMR Stratified Heatmaps - HIGH-CONFIDENCE VERSION
# ============================================================================
# Purpose: Create heatmaps for HIGH-CONFIDENCE DMR categories where BOTH
#          TES and GFP samples have meaningful signal (>2 reads)
#
# Categories:
#   1. All high-confidence DMRs
#   2. Hypermethylated only (high-confidence)
#   3. Hypomethylated only (high-confidence)
#
# NOTE: This version uses filtered DMRs that address the GFP library
#       quality issues (85% duplication, dropout artifacts).
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

# Define base meDIP directory for BigWig files
MEDIP_BASE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"

# Create output directories
OUTDIR="output/17_dmr_heatmaps_confident"
mkdir -p "$OUTDIR"/{matrices,heatmaps,beds}
mkdir -p logs

echo "=============================================="
echo "DMR Stratified Heatmaps - HIGH-CONFIDENCE"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment with deepTools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# Verify deepTools is available
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found"
    exit 1
fi

# Define input files - USE CONFIDENT DMR FILES
DMR_DIR="output/07_differential_MEDIPS_confident"
DMR_CSV="$DMR_DIR/TES_vs_GFP_DMRs_confident.csv"
HYPER_BED="$DMR_DIR/TES_vs_GFP_hypermethylated_confident.bed"
HYPO_BED="$DMR_DIR/TES_vs_GFP_hypomethylated_confident.bed"
ALL_BED="$DMR_DIR/TES_vs_GFP_all_confident.bed"

# Verify input files exist
echo ""
echo "=== Checking input files ==="
for f in "$DMR_CSV" "$HYPER_BED" "$HYPO_BED" "$ALL_BED"; do
    if [[ -f "$f" ]]; then
        if [[ "$f" == *.bed ]]; then
            count=$(wc -l < "$f")
            echo "  Found: $f ($count regions)"
        else
            echo "  Found: $f"
        fi
    else
        echo "  ERROR: Not found: $f"
        echo "  Please run 00_filter_confident_dmrs.sh first!"
        exit 1
    fi
done

# BigWig files - individual replicates (use absolute paths)
GFP_BW1="$MEDIP_BASE/results/05_bigwig/GFP-1-IP_RPKM.bw"
GFP_BW2="$MEDIP_BASE/results/05_bigwig/GFP-2-IP_RPKM.bw"
TES_BW1="$MEDIP_BASE/results/05_bigwig/TES-1-IP_RPKM.bw"
TES_BW2="$MEDIP_BASE/results/05_bigwig/TES-2-IP_RPKM.bw"

# Combined BigWig files
GFP_COMBINED="$MEDIP_BASE/results/05_bigwig/GFP_average.bw"
TES_COMBINED="$MEDIP_BASE/results/05_bigwig/TES_average.bw"

# Verify BigWig files
echo ""
echo "=== Checking BigWig files ==="
for f in "$GFP_BW1" "$GFP_BW2" "$TES_BW1" "$TES_BW2"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  WARNING: Not found: $f"
    fi
done

if [[ -f "$GFP_COMBINED" ]]; then
    echo "  Found: $GFP_COMBINED"
else
    echo "  ERROR: GFP average BigWig not found: $GFP_COMBINED"
    exit 1
fi
if [[ -f "$TES_COMBINED" ]]; then
    echo "  Found: $TES_COMBINED"
else
    echo "  ERROR: TES average BigWig not found: $TES_COMBINED"
    exit 1
fi

# ============================================================================
# Copy BED files to output directory for consistency with downstream scripts
# ============================================================================
echo ""
echo "=== Setting up BED files ==="

cp "$ALL_BED" "$OUTDIR/beds/all_dmrs.bed"
cp "$HYPER_BED" "$OUTDIR/beds/hypermethylated_dmrs.bed"
cp "$HYPO_BED" "$OUTDIR/beds/hypomethylated_dmrs.bed"

ALL_DMR_COUNT=$(wc -l < "$OUTDIR/beds/all_dmrs.bed")
HYPER_COUNT=$(wc -l < "$OUTDIR/beds/hypermethylated_dmrs.bed")
HYPO_COUNT=$(wc -l < "$OUTDIR/beds/hypomethylated_dmrs.bed")

echo "  All high-confidence DMRs: $ALL_DMR_COUNT"
echo "  Hypermethylated: $HYPER_COUNT"
echo "  Hypomethylated: $HYPO_COUNT"

# ============================================================================
# Function to create heatmap for a DMR category
# ============================================================================

create_dmr_heatmap() {
    local category=$1
    local bed_file=$2
    local title=$3

    local region_count=$(wc -l < "$bed_file")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $category - too few regions ($region_count)"
        return
    fi

    echo ""
    echo "Processing $category ($region_count regions)..."

    # Compute matrix - center on DMR, extend 5kb each side
    echo "  Computing matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$GFP_BW1" "$GFP_BW2" "$TES_BW1" "$TES_BW2" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${category}_matrix.gz" \
        --outFileNameMatrix "$OUTDIR/matrices/${category}_matrix.tab" \
        2>/dev/null

    # Generate heatmap - sorted by mean signal
    echo "  Creating heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${category}_matrix.gz" \
        -o "$OUTDIR/heatmaps/${category}_heatmap.png" \
        --colorMap RdYlBu_r \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 5 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "$title" \
        --legendLocation best \
        --dpi 300 \
        --outFileSortedRegions "$OUTDIR/beds/${category}_sorted.bed" \
        2>/dev/null

    # Generate profile plot
    echo "  Creating profile..."
    plotProfile \
        -m "$OUTDIR/matrices/${category}_matrix.gz" \
        -o "$OUTDIR/heatmaps/${category}_profile.png" \
        --plotTitle "$title" \
        --perGroup \
        --colors "#2166AC" "#4393C3" "#B2182B" "#D6604D" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null

    echo "  Done with $category"
}

# Function to create heatmap using combined replicates
create_dmr_heatmap_combined() {
    local category=$1
    local bed_file=$2
    local title=$3

    local region_count=$(wc -l < "$bed_file")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $category (combined) - too few regions ($region_count)"
        return
    fi

    echo ""
    echo "Processing $category - Combined Replicates ($region_count regions)..."

    # Compute matrix
    echo "  Computing matrix (combined)..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$GFP_COMBINED" "$TES_COMBINED" \
        --samplesLabel "GFP (combined)" "TES (combined)" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${category}_combined_matrix.gz" \
        --outFileNameMatrix "$OUTDIR/matrices/${category}_combined_matrix.tab" \
        2>/dev/null

    # Generate heatmap
    echo "  Creating heatmap (combined)..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${category}_combined_matrix.gz" \
        -o "$OUTDIR/heatmaps/${category}_combined_heatmap.png" \
        --colorMap RdYlBu_r \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "$title (Combined Replicates)" \
        --legendLocation best \
        --dpi 300 \
        --outFileSortedRegions "$OUTDIR/beds/${category}_combined_sorted.bed" \
        2>/dev/null

    # Generate profile plot
    echo "  Creating profile (combined)..."
    plotProfile \
        -m "$OUTDIR/matrices/${category}_combined_matrix.gz" \
        -o "$OUTDIR/heatmaps/${category}_combined_profile.png" \
        --plotTitle "$title (Combined Replicates)" \
        --perGroup \
        --colors "#2166AC" "#B2182B" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null

    echo "  Done with $category (combined)"
}

# ============================================================================
# Generate heatmaps for each DMR category
# ============================================================================

echo ""
echo "=== Generating heatmaps for HIGH-CONFIDENCE DMR categories ==="

# Individual Replicates
echo ""
echo "--- Individual Replicate Heatmaps ---"
create_dmr_heatmap "all_dmrs" "$OUTDIR/beds/all_dmrs.bed" "High-Confidence DMRs (Both samples >2 reads)"
create_dmr_heatmap "hypermethylated" "$OUTDIR/beds/hypermethylated_dmrs.bed" "High-Conf Hypermethylated (TES > GFP)"
create_dmr_heatmap "hypomethylated" "$OUTDIR/beds/hypomethylated_dmrs.bed" "High-Conf Hypomethylated (TES < GFP)"

# Combined Replicates
echo ""
echo "--- Combined Replicate Heatmaps ---"
create_dmr_heatmap_combined "all_dmrs" "$OUTDIR/beds/all_dmrs.bed" "High-Confidence DMRs"
create_dmr_heatmap_combined "hypermethylated" "$OUTDIR/beds/hypermethylated_dmrs.bed" "High-Conf Hypermethylated"
create_dmr_heatmap_combined "hypomethylated" "$OUTDIR/beds/hypomethylated_dmrs.bed" "High-Conf Hypomethylated"

# ============================================================================
# Create side-by-side comparison
# ============================================================================

echo ""
echo "=== Creating Hyper vs Hypo comparison ==="

if [[ $HYPER_COUNT -ge 10 ]] && [[ $HYPO_COUNT -ge 10 ]]; then
    echo "Creating comparison plots..."

    # Combined replicates comparison
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$OUTDIR/beds/hypermethylated_dmrs.bed" "$OUTDIR/beds/hypomethylated_dmrs.bed" \
        -S "$GFP_COMBINED" "$TES_COMBINED" \
        --samplesLabel "GFP (combined)" "TES (combined)" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/hyper_vs_hypo_combined_matrix.gz" \
        2>/dev/null

    plotProfile \
        -m "$OUTDIR/matrices/hyper_vs_hypo_combined_matrix.gz" \
        -o "$OUTDIR/heatmaps/hyper_vs_hypo_combined_comparison.png" \
        --plotTitle "High-Confidence: Hyper vs Hypo DMRs" \
        --perGroup \
        --colors "#2166AC" "#B2182B" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null

    plotHeatmap \
        -m "$OUTDIR/matrices/hyper_vs_hypo_combined_matrix.gz" \
        -o "$OUTDIR/heatmaps/hyper_vs_hypo_combined_heatmap.png" \
        --colorMap RdYlBu_r \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "High-Confidence: Hyper vs Hypo DMRs" \
        --legendLocation best \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --dpi 300 \
        2>/dev/null

    echo "  Done with comparison"
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="
echo ""
echo "HIGH-CONFIDENCE DMR Categories:"
echo "  All DMRs: $ALL_DMR_COUNT regions"
echo "  Hypermethylated: $HYPER_COUNT regions ($(echo "scale=1; 100*$HYPER_COUNT/$ALL_DMR_COUNT" | bc)%)"
echo "  Hypomethylated: $HYPO_COUNT regions ($(echo "scale=1; 100*$HYPO_COUNT/$ALL_DMR_COUNT" | bc)%)"
echo ""
echo "NOTE: These are DMRs where BOTH TES and GFP have >2 mean reads,"
echo "      filtering out artifactual calls from GFP sample dropout."
echo ""
echo "Output files in: $OUTDIR/"
echo "  - beds/: BED files for each category"
echo "  - matrices/: deepTools matrices"
echo "  - heatmaps/: PNG heatmaps and profiles (300 DPI)"
echo ""
echo "Finished: $(date)"
