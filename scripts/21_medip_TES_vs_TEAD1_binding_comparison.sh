#!/bin/bash
#SBATCH --job-name=medip_TES_vs_TEAD1
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --output=logs/21_medip_TES_vs_TEAD1_comparison.out
#SBATCH --error=logs/21_medip_TES_vs_TEAD1_comparison.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# Profile Overlay: meDIP Signal at TES vs TEAD1 Binding Sites
# ============================================================================
#
# Purpose: Compare ABSOLUTE meDIP signal between TES and TEAD1 binding sites
#          using profile overlays
#
# Key question: Is methylation different at TES-bound vs TEAD1-bound regions?
#
# Comparisons:
#   1. All TES sites vs All TEAD1 sites (same plot, different region groups)
#   2. TES-unique vs TEAD1-unique vs Shared sites
#   3. GFP meDIP and TES meDIP at each binding category
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/21_medip_TES_vs_TEAD1_comparison"
mkdir -p "$OUTDIR"/{beds,matrices,profiles,heatmaps}
mkdir -p logs

echo "=============================================="
echo "meDIP Profile Comparison: TES vs TEAD1 Binding"
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
# Create binding site categories
# ============================================================================

echo ""
echo "=== Creating binding site categories ==="

# Standardize peak files
cut -f1-3 "$TES_PEAKS_USE" | sort -k1,1 -k2,2n > "$OUTDIR/beds/TES_all.bed"
cut -f1-3 "$TEAD1_PEAKS_USE" | sort -k1,1 -k2,2n > "$OUTDIR/beds/TEAD1_all.bed"

TES_COUNT=$(wc -l < "$OUTDIR/beds/TES_all.bed")
TEAD1_COUNT=$(wc -l < "$OUTDIR/beds/TEAD1_all.bed")
echo "  TES peaks: $TES_COUNT"
echo "  TEAD1 peaks: $TEAD1_COUNT"

# TES-unique (TES but not TEAD1)
bedtools intersect -a "$OUTDIR/beds/TES_all.bed" \
    -b "$OUTDIR/beds/TEAD1_all.bed" -v > "$OUTDIR/beds/TES_unique.bed"
TES_UNIQUE=$(wc -l < "$OUTDIR/beds/TES_unique.bed")
echo "  TES-unique: $TES_UNIQUE"

# TEAD1-unique (TEAD1 but not TES)
bedtools intersect -a "$OUTDIR/beds/TEAD1_all.bed" \
    -b "$OUTDIR/beds/TES_all.bed" -v > "$OUTDIR/beds/TEAD1_unique.bed"
TEAD1_UNIQUE=$(wc -l < "$OUTDIR/beds/TEAD1_unique.bed")
echo "  TEAD1-unique: $TEAD1_UNIQUE"

# Shared (TES overlapping TEAD1)
bedtools intersect -a "$OUTDIR/beds/TES_all.bed" \
    -b "$OUTDIR/beds/TEAD1_all.bed" -u > "$OUTDIR/beds/shared.bed"
SHARED=$(wc -l < "$OUTDIR/beds/shared.bed")
echo "  Shared: $SHARED"

# ============================================================================
# COMPARISON 1: All TES sites vs All TEAD1 sites (Profile Overlay)
# ============================================================================

echo ""
echo "=== Comparison 1: All TES vs All TEAD1 sites ==="

# GFP meDIP at TES vs TEAD1 sites
echo "  Computing matrix: GFP meDIP at TES vs TEAD1..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_all.bed" "$OUTDIR/beds/TEAD1_all.bed" \
    -S "$GFP_MEDIP_BW1" \
    --samplesLabel "GFP_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/TES_vs_TEAD1_GFP_medip.gz" \
    2>/dev/null

# Profile overlay: GFP meDIP
echo "  Creating profile: GFP meDIP at TES vs TEAD1..."
plotProfile \
    -m "$OUTDIR/matrices/TES_vs_TEAD1_GFP_medip.gz" \
    -o "$OUTDIR/profiles/GFP_medip_TES_vs_TEAD1_overlay.png" \
    --plotTitle "GFP meDIP: TES Binding Sites vs TEAD1 Binding Sites" \
    --perGroup \
    --colors "#2166AC" "#B2182B" \
    --regionsLabel "TES sites (n=$TES_COUNT)" "TEAD1 sites (n=$TEAD1_COUNT)" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: GFP profile failed"

# TES meDIP at TES vs TEAD1 sites
echo "  Computing matrix: TES meDIP at TES vs TEAD1..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_all.bed" "$OUTDIR/beds/TEAD1_all.bed" \
    -S "$TES_MEDIP_BW1" \
    --samplesLabel "TES_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/TES_vs_TEAD1_TES_medip.gz" \
    2>/dev/null

# Profile overlay: TES meDIP
echo "  Creating profile: TES meDIP at TES vs TEAD1..."
plotProfile \
    -m "$OUTDIR/matrices/TES_vs_TEAD1_TES_medip.gz" \
    -o "$OUTDIR/profiles/TES_medip_TES_vs_TEAD1_overlay.png" \
    --plotTitle "TES meDIP: TES Binding Sites vs TEAD1 Binding Sites" \
    --perGroup \
    --colors "#2166AC" "#B2182B" \
    --regionsLabel "TES sites (n=$TES_COUNT)" "TEAD1 sites (n=$TEAD1_COUNT)" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: TES profile failed"

# Combined: Both GFP and TES meDIP at TES vs TEAD1 sites
echo "  Computing combined matrix..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_all.bed" "$OUTDIR/beds/TEAD1_all.bed" \
    -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
    --samplesLabel "GFP_meDIP" "TES_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/TES_vs_TEAD1_both_medip.gz" \
    2>/dev/null

# Heatmap showing both conditions at both region types
echo "  Creating heatmap: GFP & TES meDIP at TES vs TEAD1..."
plotHeatmap \
    -m "$OUTDIR/matrices/TES_vs_TEAD1_both_medip.gz" \
    -o "$OUTDIR/heatmaps/TES_vs_TEAD1_both_medip_heatmap.png" \
    --colorMap Greens Purples \
    --sortRegions descend \
    --sortUsing mean \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --xAxisLabel "" \
    --refPointLabel "Peak" \
    --regionsLabel "TES sites" "TEAD1 sites" \
    --plotTitle "meDIP at TES vs TEAD1 Binding Sites" \
    --legendLocation upper-right \
    --dpi 300 \
    2>/dev/null || echo "WARNING: combined heatmap failed"

# Profile with per-sample panels
echo "  Creating per-sample profile..."
plotProfile \
    -m "$OUTDIR/matrices/TES_vs_TEAD1_both_medip.gz" \
    -o "$OUTDIR/profiles/TES_vs_TEAD1_both_medip_perSample.png" \
    --plotTitle "meDIP at TES vs TEAD1 Binding Sites" \
    --perSample \
    --colors "#2166AC" "#B2182B" \
    --regionsLabel "TES sites" "TEAD1 sites" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 12 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: per-sample profile failed"

# ============================================================================
# COMPARISON 2: TES-unique vs TEAD1-unique vs Shared (3-way)
# ============================================================================

echo ""
echo "=== Comparison 2: TES-unique vs Shared vs TEAD1-unique ==="

# GFP meDIP at 3 categories
echo "  Computing matrix: GFP meDIP at 3 categories..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_unique.bed" \
       "$OUTDIR/beds/shared.bed" \
       "$OUTDIR/beds/TEAD1_unique.bed" \
    -S "$GFP_MEDIP_BW1" \
    --samplesLabel "GFP_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/3way_GFP_medip.gz" \
    2>/dev/null

# Profile: GFP meDIP at 3 categories
echo "  Creating 3-way profile: GFP meDIP..."
plotProfile \
    -m "$OUTDIR/matrices/3way_GFP_medip.gz" \
    -o "$OUTDIR/profiles/GFP_medip_3way_overlay.png" \
    --plotTitle "GFP meDIP: TES-unique vs Shared vs TEAD1-unique" \
    --perGroup \
    --colors "#2166AC" "#7570B3" "#B2182B" \
    --regionsLabel "TES-unique (n=$TES_UNIQUE)" "Shared (n=$SHARED)" "TEAD1-unique (n=$TEAD1_UNIQUE)" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: 3-way GFP profile failed"

# TES meDIP at 3 categories
echo "  Computing matrix: TES meDIP at 3 categories..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_unique.bed" \
       "$OUTDIR/beds/shared.bed" \
       "$OUTDIR/beds/TEAD1_unique.bed" \
    -S "$TES_MEDIP_BW1" \
    --samplesLabel "TES_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/3way_TES_medip.gz" \
    2>/dev/null

# Profile: TES meDIP at 3 categories
echo "  Creating 3-way profile: TES meDIP..."
plotProfile \
    -m "$OUTDIR/matrices/3way_TES_medip.gz" \
    -o "$OUTDIR/profiles/TES_medip_3way_overlay.png" \
    --plotTitle "TES meDIP: TES-unique vs Shared vs TEAD1-unique" \
    --perGroup \
    --colors "#2166AC" "#7570B3" "#B2182B" \
    --regionsLabel "TES-unique (n=$TES_UNIQUE)" "Shared (n=$SHARED)" "TEAD1-unique (n=$TEAD1_UNIQUE)" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: 3-way TES profile failed"

# Combined: Both meDIP conditions at 3 categories
echo "  Computing combined 3-way matrix..."
computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_unique.bed" \
       "$OUTDIR/beds/shared.bed" \
       "$OUTDIR/beds/TEAD1_unique.bed" \
    -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
    --samplesLabel "GFP_meDIP" "TES_meDIP" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/3way_both_medip.gz" \
    2>/dev/null

# Heatmap: 3-way comparison
echo "  Creating 3-way heatmap..."
plotHeatmap \
    -m "$OUTDIR/matrices/3way_both_medip.gz" \
    -o "$OUTDIR/heatmaps/3way_both_medip_heatmap.png" \
    --colorMap Greens Purples \
    --sortRegions descend \
    --sortUsing mean \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --xAxisLabel "" \
    --refPointLabel "Peak" \
    --regionsLabel "TES-unique" "Shared" "TEAD1-unique" \
    --plotTitle "meDIP at Binding Site Categories" \
    --legendLocation upper-right \
    --dpi 300 \
    2>/dev/null || echo "WARNING: 3-way heatmap failed"

# Per-sample profile
plotProfile \
    -m "$OUTDIR/matrices/3way_both_medip.gz" \
    -o "$OUTDIR/profiles/3way_both_medip_perSample.png" \
    --plotTitle "GFP vs TES meDIP at Binding Categories" \
    --perSample \
    --colors "#2166AC" "#7570B3" "#B2182B" \
    --regionsLabel "TES-unique" "Shared" "TEAD1-unique" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 12 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: per-sample 3-way profile failed"

# ============================================================================
# COMPARISON 3: All replicates for comprehensive view
# ============================================================================

echo ""
echo "=== Comparison 3: All replicates at TES vs TEAD1 sites ==="

computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$OUTDIR/beds/TES_all.bed" "$OUTDIR/beds/TEAD1_all.bed" \
    -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
    --samplesLabel "GFP_rep1" "GFP_rep2" "TES_rep1" "TES_rep2" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/TES_vs_TEAD1_all_reps.gz" \
    2>/dev/null

# Profile showing all replicates
echo "  Creating all-replicates profile..."
plotProfile \
    -m "$OUTDIR/matrices/TES_vs_TEAD1_all_reps.gz" \
    -o "$OUTDIR/profiles/TES_vs_TEAD1_all_reps_overlay.png" \
    --plotTitle "All meDIP Replicates: TES vs TEAD1 Sites" \
    --perGroup \
    --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
    --regionsLabel "TES sites" "TEAD1 sites" \
    --legendLocation upper-right \
    --refPointLabel "Peak" \
    --plotHeight 6 \
    --plotWidth 10 \
    --dpi 300 \
    2>/dev/null || echo "WARNING: all-reps profile failed"

# Heatmap with all replicates
echo "  Creating all-replicates heatmap..."
plotHeatmap \
    -m "$OUTDIR/matrices/TES_vs_TEAD1_all_reps.gz" \
    -o "$OUTDIR/heatmaps/TES_vs_TEAD1_all_reps_heatmap.png" \
    --colorMap YlGn YlGn Purples Purples \
    --sortRegions descend \
    --sortUsing mean \
    --heatmapHeight 15 \
    --heatmapWidth 2.5 \
    --xAxisLabel "" \
    --refPointLabel "Peak" \
    --regionsLabel "TES sites" "TEAD1 sites" \
    --plotTitle "All meDIP Replicates at TES vs TEAD1 Sites" \
    --legendLocation upper-right \
    --dpi 300 \
    2>/dev/null || echo "WARNING: all-reps heatmap failed"

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
echo "Binding site categories:"
echo "  TES all: $TES_COUNT"
echo "  TEAD1 all: $TEAD1_COUNT"
echo "  TES-unique: $TES_UNIQUE"
echo "  TEAD1-unique: $TEAD1_UNIQUE"
echo "  Shared: $SHARED"
echo ""
echo "Profile overlays generated:"
echo ""
echo "  COMPARISON 1: TES sites vs TEAD1 sites"
ls -1 "$OUTDIR/profiles/"*TES_vs_TEAD1*.png 2>/dev/null | while read f; do echo "    $(basename $f)"; done
echo ""
echo "  COMPARISON 2: 3-way (TES-unique vs Shared vs TEAD1-unique)"
ls -1 "$OUTDIR/profiles/"*3way*.png 2>/dev/null | while read f; do echo "    $(basename $f)"; done
echo ""
echo "Heatmaps:"
ls -1 "$OUTDIR/heatmaps/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "=============================================="
echo "INTERPRETATION GUIDE"
echo "=============================================="
echo ""
echo "Key profiles to examine:"
echo ""
echo "1. GFP_medip_TES_vs_TEAD1_overlay.png"
echo "   - Shows baseline methylation (before TES expression)"
echo "   - If lines overlap: similar baseline methylation at both site types"
echo ""
echo "2. TES_medip_TES_vs_TEAD1_overlay.png"
echo "   - Shows methylation after TES expression"
echo "   - If TES sites > TEAD1 sites: TES binding associated with more methylation"
echo ""
echo "3. 3way profiles (GFP and TES)"
echo "   - Compare methylation across binding categories"
echo "   - TES-unique vs Shared vs TEAD1-unique"
echo ""
echo "Expected results if TES recruits methylation machinery:"
echo "   - TES meDIP higher at TES-unique and Shared sites"
echo "   - Difference between GFP and TES conditions at TES sites"
echo ""
echo "Finished: $(date)"
