#!/bin/bash
#SBATCH --job-name=medip_union_binding
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --output=logs/20_medip_at_union_binding.out
#SBATCH --error=logs/20_medip_at_union_binding.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# meDIP Signal at Union of TES + TEAD1 Binding Sites
# ============================================================================
#
# Purpose: Visualize ABSOLUTE meDIP signal at ANY TF binding site
#          (union of TES and TEAD1 peaks)
#
# This shows: What is the methylation level at regions bound by either TF?
#
# Output:
#   - meDIP signal profiles at union binding sites
#   - Heatmaps showing GFP and TES meDIP signal separately
#   - Summary statistics
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/20_medip_union_binding"
mkdir -p "$OUTDIR"/{beds,matrices,heatmaps,profiles}
mkdir -p logs

echo "=============================================="
echo "meDIP Signal at Union of TES + TEAD1 Binding"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ============================================================================
# Define input files
# ============================================================================

# Cut&Tag peak files
CUTNTAG_PEAKS="../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/consensus_peaks"
TES_PEAKS="$CUTNTAG_PEAKS/TES_consensus_peaks.bed"
TEAD1_PEAKS="$CUTNTAG_PEAKS/TEAD1_consensus_peaks.bed"

# Alternative locations
CUTNTAG_PEAKS_ALT="../SRF_Eva_CUTandTAG/results/05_peaks_narrow"

# meDIP BigWig files (ABSOLUTE signal)
MEDIP_DIR="results/05_bigwig"
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
GFP_MEDIP_BW2="$MEDIP_DIR/GFP-2-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"
TES_MEDIP_BW2="$MEDIP_DIR/TES-2-IP_RPKM.bw"

# ============================================================================
# Locate and prepare input files
# ============================================================================

echo ""
echo "=== Checking input files ==="

# Find TES peaks
if [[ -f "$TES_PEAKS" ]]; then
    echo "  Found TES consensus peaks: $TES_PEAKS"
    TES_PEAKS_USE="$TES_PEAKS"
elif [[ -f "$CUTNTAG_PEAKS_ALT/TES-1_peaks.narrowPeak" ]]; then
    echo "  Merging individual TES peaks..."
    TES_PEAKS_USE="$OUTDIR/beds/TES_merged_peaks.bed"
    cat "$CUTNTAG_PEAKS_ALT"/TES-*_peaks.narrowPeak | \
        cut -f1-3 | sort -k1,1 -k2,2n | \
        bedtools merge -i - > "$TES_PEAKS_USE"
else
    echo "  ERROR: No TES peaks found"
    exit 1
fi

# Find TEAD1 peaks
if [[ -f "$TEAD1_PEAKS" ]]; then
    echo "  Found TEAD1 consensus peaks: $TEAD1_PEAKS"
    TEAD1_PEAKS_USE="$TEAD1_PEAKS"
elif [[ -f "$CUTNTAG_PEAKS_ALT/TEAD1-1_peaks.narrowPeak" ]]; then
    echo "  Merging individual TEAD1 peaks..."
    TEAD1_PEAKS_USE="$OUTDIR/beds/TEAD1_merged_peaks.bed"
    cat "$CUTNTAG_PEAKS_ALT"/TEAD1-*_peaks.narrowPeak | \
        cut -f1-3 | sort -k1,1 -k2,2n | \
        bedtools merge -i - > "$TEAD1_PEAKS_USE"
else
    echo "  ERROR: No TEAD1 peaks found"
    exit 1
fi

# Check meDIP files
for f in "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  ERROR: Not found: $f"
        exit 1
    fi
done

# ============================================================================
# Create union of binding sites
# ============================================================================

echo ""
echo "=== Creating union of TES + TEAD1 binding sites ==="

# Standardize peak files (first 3 columns)
cut -f1-3 "$TES_PEAKS_USE" | sort -k1,1 -k2,2n > "$OUTDIR/beds/TES_peaks.bed"
cut -f1-3 "$TEAD1_PEAKS_USE" | sort -k1,1 -k2,2n > "$OUTDIR/beds/TEAD1_peaks.bed"

TES_COUNT=$(wc -l < "$OUTDIR/beds/TES_peaks.bed")
TEAD1_COUNT=$(wc -l < "$OUTDIR/beds/TEAD1_peaks.bed")
echo "  TES peaks: $TES_COUNT"
echo "  TEAD1 peaks: $TEAD1_COUNT"

# Create UNION (merge overlapping peaks from both TFs)
cat "$OUTDIR/beds/TES_peaks.bed" "$OUTDIR/beds/TEAD1_peaks.bed" | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > "$OUTDIR/beds/union_TES_TEAD1_peaks.bed"

UNION_COUNT=$(wc -l < "$OUTDIR/beds/union_TES_TEAD1_peaks.bed")
echo "  Union peaks (merged): $UNION_COUNT"

# Also create simple concatenation (non-merged) for comparison
cat "$OUTDIR/beds/TES_peaks.bed" "$OUTDIR/beds/TEAD1_peaks.bed" | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/all_binding_sites_unmerged.bed"

UNMERGED_COUNT=$(wc -l < "$OUTDIR/beds/all_binding_sites_unmerged.bed")
echo "  All binding sites (unmerged): $UNMERGED_COUNT"

# ============================================================================
# Compute meDIP signal matrices at union binding sites
# ============================================================================

echo ""
echo "=== Computing meDIP signal at union binding sites ==="

# Matrix with all 4 meDIP samples (absolute signal)
echo "  Computing matrix for all meDIP replicates..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/union_TES_TEAD1_peaks.bed" \
    -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
    --samplesLabel "GFP_rep1" "GFP_rep2" "TES_rep1" "TES_rep2" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/union_binding_medip_all_reps.gz" \
    2>/dev/null || { echo "ERROR: computeMatrix failed"; exit 1; }

echo "  Matrix computed successfully"

# ============================================================================
# Generate heatmaps
# ============================================================================

echo ""
echo "=== Generating heatmaps ==="

# Heatmap showing all 4 replicates
echo "  Creating heatmap (all replicates)..."
plotHeatmap \
    -m "$OUTDIR/matrices/union_binding_medip_all_reps.gz" \
    -o "$OUTDIR/heatmaps/union_binding_medip_all_reps_heatmap.png" \
    --colorMap YlGn YlGn Purples Purples \
    --sortRegions descend \
    --sortUsing mean \
    --heatmapHeight 15 \
    --heatmapWidth 2.5 \
    --xAxisLabel "" \
    --refPointLabel "Peak" \
    --plotTitle "meDIP Signal at Union TES+TEAD1 Binding Sites (n=$UNION_COUNT)" \
    --legendLocation upper-right \
    --dpi 300 \
    2>/dev/null || echo "WARNING: heatmap failed"

# Heatmap sorted by TES meDIP signal (to show where methylation is highest)
echo "  Creating heatmap (sorted by TES signal)..."
plotHeatmap \
    -m "$OUTDIR/matrices/union_binding_medip_all_reps.gz" \
    -o "$OUTDIR/heatmaps/union_binding_medip_sorted_by_TES.png" \
    --colorMap YlGn YlGn Purples Purples \
    --sortRegions descend \
    --sortUsing mean \
    --sortUsingSamples 3 4 \
    --heatmapHeight 15 \
    --heatmapWidth 2.5 \
    --xAxisLabel "" \
    --refPointLabel "Peak" \
    --plotTitle "meDIP at Union Sites (Sorted by TES meDIP)" \
    --legendLocation upper-right \
    --dpi 300 \
    --outFileSortedRegions "$OUTDIR/beds/union_sorted_by_TES_medip.bed" \
    2>/dev/null || echo "WARNING: heatmap failed"

# ============================================================================
# Generate profile plots
# ============================================================================

echo ""
echo "=== Generating profile plots ==="

# Profile with all 4 replicates overlaid
echo "  Creating profile (all replicates)..."
plotProfile \
    -m "$OUTDIR/matrices/union_binding_medip_all_reps.gz" \
    -o "$OUTDIR/profiles/union_binding_medip_all_reps_profile.png" \
    --plotTitle "meDIP Signal at Union TES+TEAD1 Binding Sites" \
    --perGroup \
    --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: profile failed"

# Profile per sample (separate Y-axes)
echo "  Creating profile (per sample panels)..."
plotProfile \
    -m "$OUTDIR/matrices/union_binding_medip_all_reps.gz" \
    -o "$OUTDIR/profiles/union_binding_medip_per_sample_profile.png" \
    --plotTitle "meDIP at Union Binding Sites (Per Sample)" \
    --perSample \
    --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 8 \
    --plotWidth 14 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: per-sample profile failed"

# ============================================================================
# Create averaged GFP vs TES comparison
# ============================================================================

echo ""
echo "=== Creating GFP vs TES averaged comparison ==="

# Create matrix with just 2 samples (GFP average, TES average would require merging)
# For now, use representative replicates
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/union_TES_TEAD1_peaks.bed" \
    -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
    --samplesLabel "GFP_meDIP" "TES_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/union_binding_GFP_vs_TES.gz" \
    2>/dev/null

# Clean comparison heatmap
echo "  Creating GFP vs TES comparison heatmap..."
plotHeatmap \
    -m "$OUTDIR/matrices/union_binding_GFP_vs_TES.gz" \
    -o "$OUTDIR/heatmaps/union_binding_GFP_vs_TES_heatmap.png" \
    --colorMap Greens Purples \
    --sortRegions descend \
    --sortUsing mean \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --xAxisLabel "" \
    --refPointLabel "Peak" \
    --plotTitle "meDIP: GFP vs TES at Union Binding Sites" \
    --legendLocation upper-right \
    --dpi 300 \
    2>/dev/null || echo "WARNING: comparison heatmap failed"

# Clean comparison profile
echo "  Creating GFP vs TES comparison profile..."
plotProfile \
    -m "$OUTDIR/matrices/union_binding_GFP_vs_TES.gz" \
    -o "$OUTDIR/profiles/union_binding_GFP_vs_TES_profile.png" \
    --plotTitle "meDIP Signal: GFP vs TES at Union Binding Sites" \
    --perGroup \
    --colors "#1B7837" "#762A83" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: comparison profile failed"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Summary:"
echo "  TES binding sites: $TES_COUNT"
echo "  TEAD1 binding sites: $TEAD1_COUNT"
echo "  Union (merged): $UNION_COUNT"
echo ""
echo "Key outputs:"
echo "  beds/union_TES_TEAD1_peaks.bed - Union of all binding sites"
echo ""
echo "Heatmaps:"
ls -1 "$OUTDIR/heatmaps/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Profiles:"
ls -1 "$OUTDIR/profiles/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Interpretation:"
echo "  - Compare GFP vs TES meDIP signal at binding sites"
echo "  - Higher TES signal indicates methylation gain upon TES expression"
echo "  - Profile shape shows methylation distribution around binding sites"
echo ""
echo "Finished: $(date)"
