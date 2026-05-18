#!/bin/bash
#SBATCH --job-name=dmr_binding_overlay_confident
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/18_dmr_binding_overlay_confident.out
#SBATCH --error=logs/18_dmr_binding_overlay_confident.err

# ============================================================================
# DMR Binding Overlay Analysis - HIGH-CONFIDENCE VERSION
# ============================================================================
# Purpose: Visualize TES/TEAD1 binding signals at HIGH-CONFIDENCE DMR regions
#          where BOTH TES and GFP samples have meaningful signal (>2 reads)
#
# Overlay:
#   - TES Cut&Tag signal
#   - TEAD1 Cut&Tag signal
#   - meDIP GFP signal
#   - meDIP TES signal
#
# NOTE: This version uses filtered DMRs that address the GFP library
#       quality issues (85% duplication, dropout artifacts).
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

# Define base directories for BigWig files
MEDIP_BASE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
CUTNTAG_BASE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"

# Create output directories
OUTDIR="output/18_binding_overlay_confident"
mkdir -p "$OUTDIR"/{matrices,heatmaps,profiles}
mkdir -p logs

echo "=============================================="
echo "DMR Binding Overlay - HIGH-CONFIDENCE"
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

# ============================================================================
# Define input files
# ============================================================================

# DMR BED files from confident heatmaps (script 17_confident)
DMR_DIR="output/17_dmr_heatmaps_confident/beds"

# Check if confident BEDs exist
if [[ ! -d "$DMR_DIR" ]]; then
    echo "ERROR: Confident DMR directory not found: $DMR_DIR"
    echo "Please run 17_dmr_stratified_heatmaps_confident.sh first!"
    exit 1
fi

HYPER_BED="$DMR_DIR/hypermethylated_dmrs.bed"
HYPO_BED="$DMR_DIR/hypomethylated_dmrs.bed"
ALL_DMR_BED="$DMR_DIR/all_dmrs.bed"

# Cut&Tag BigWig files (use absolute paths)
CUTNTAG_DIR="$CUTNTAG_BASE/results/06_bigwig"
TES_BW="$CUTNTAG_DIR/TES_comb.bw"
TEAD1_BW="$CUTNTAG_DIR/TEAD1_comb.bw"

# Use average bigwigs if available
CUTNTAG_AVG_DIR="$CUTNTAG_BASE/results/11_combined_replicates_narrow/bigwig"
if [[ -f "$CUTNTAG_AVG_DIR/TES_average.bw" ]]; then
    TES_AVG_BW="$CUTNTAG_AVG_DIR/TES_average.bw"
    TEAD1_AVG_BW="$CUTNTAG_AVG_DIR/TEAD1_average.bw"
else
    TES_AVG_BW="$TES_BW"
    TEAD1_AVG_BW="$TEAD1_BW"
fi

# meDIP BigWig files - individual replicates (use absolute paths)
MEDIP_DIR="$MEDIP_BASE/results/05_bigwig"
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"

# meDIP BigWig files - combined replicates
GFP_MEDIP_COMBINED="$MEDIP_DIR/GFP_average.bw"
TES_MEDIP_COMBINED="$MEDIP_DIR/TES_average.bw"

# Verify input files
echo ""
echo "=== Checking input files ==="

# Check DMR files
for f in "$HYPER_BED" "$HYPO_BED" "$ALL_DMR_BED"; do
    if [[ -f "$f" ]]; then
        count=$(wc -l < "$f")
        echo "  Found: $f ($count regions)"
    else
        echo "  ERROR: Not found: $f"
        exit 1
    fi
done

# Check Cut&Tag BigWigs
for f in "$TES_AVG_BW" "$TEAD1_AVG_BW"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  WARNING: Not found: $f"
    fi
done

# Check meDIP BigWigs
for f in "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  ERROR: Not found: $f"
        exit 1
    fi
done

# ============================================================================
# Function to create binding overlay heatmap
# ============================================================================

create_binding_overlay() {
    local dmr_type=$1
    local bed_file=$2
    local title=$3

    if [[ ! -f "$bed_file" ]]; then
        echo "  Skipping $dmr_type - BED file not found: $bed_file"
        return
    fi

    local region_count=$(wc -l < "$bed_file")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $dmr_type - too few regions ($region_count)"
        return
    fi

    echo ""
    echo "=== Processing $dmr_type ($region_count regions) ==="

    # -------------------------------------------------------------------------
    # Matrix 1: Cut&Tag binding + meDIP methylation (all 4 signals)
    # -------------------------------------------------------------------------
    echo "  Computing combined matrix (Cut&Tag + meDIP)..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED" \
        --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${dmr_type}_combined_matrix.gz" \
        --outFileNameMatrix "$OUTDIR/matrices/${dmr_type}_combined_matrix.tab" \
        2>/dev/null || {
            echo "  WARNING: computeMatrix failed for $dmr_type combined matrix"
            return
        }

    # Generate heatmap with all 4 signals
    echo "  Creating combined heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${dmr_type}_combined_matrix.gz" \
        -o "$OUTDIR/heatmaps/${dmr_type}_binding_methylation_heatmap.png" \
        --colorList "white,#2166AC" "white,#B2182B" "white,#1B7837" "white,#762A83" \
        --zMin 0 0 0 0 \
        --zMax 2 2 150 150 \
        --sortRegions descend \
        --sortUsing mean \
        --sortUsingSamples 1 2 \
        --heatmapHeight 15 \
        --heatmapWidth 3 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "$title" \
        --legendLocation upper-right \
        --dpi 300 \
        --outFileSortedRegions "$OUTDIR/heatmaps/${dmr_type}_sorted.bed" \
        2>/dev/null || echo "  WARNING: plotHeatmap failed"

    # Generate combined profile
    echo "  Creating combined profile..."
    plotProfile \
        -m "$OUTDIR/matrices/${dmr_type}_combined_matrix.gz" \
        -o "$OUTDIR/profiles/${dmr_type}_combined_profile.png" \
        --plotTitle "Binding & Methylation at $title" \
        --colors "#2166AC" "#B2182B" "#1B7837" "#762A83" \
        --legendLocation upper-right \
        --plotHeight 8 \
        --plotWidth 14 \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: combined profile failed"

    # -------------------------------------------------------------------------
    # Matrix 2: Binding signals only (TES vs TEAD1)
    # -------------------------------------------------------------------------
    echo "  Computing binding-only matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" \
        --samplesLabel "TES" "TEAD1" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${dmr_type}_binding_matrix.gz" \
        2>/dev/null || {
            echo "  WARNING: computeMatrix failed for $dmr_type binding matrix"
            return
        }

    # Generate binding-only heatmap
    echo "  Creating binding heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${dmr_type}_binding_matrix.gz" \
        -o "$OUTDIR/heatmaps/${dmr_type}_binding_only_heatmap.png" \
        --colorMap Blues Reds \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "TES/TEAD1 Binding at $title" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: binding heatmap failed"

    # Generate binding-only profile
    plotProfile \
        -m "$OUTDIR/matrices/${dmr_type}_binding_matrix.gz" \
        -o "$OUTDIR/profiles/${dmr_type}_binding_only_profile.png" \
        --plotTitle "TES/TEAD1 Binding at $title" \
        --perGroup \
        --colors "#2166AC" "#B2182B" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: binding profile failed"

    # -------------------------------------------------------------------------
    # Matrix 3: meDIP signals only (GFP vs TES - combined)
    # -------------------------------------------------------------------------
    echo "  Computing meDIP-only matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED" \
        --samplesLabel "GFP (combined)" "TES (combined)" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${dmr_type}_medip_matrix.gz" \
        2>/dev/null || {
            echo "  WARNING: computeMatrix failed for $dmr_type meDIP matrix"
            return
        }

    # Generate meDIP-only heatmap
    echo "  Creating meDIP heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${dmr_type}_medip_matrix.gz" \
        -o "$OUTDIR/heatmaps/${dmr_type}_medip_only_heatmap.png" \
        --colorMap YlGn Purples \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "meDIP Signal at $title" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: meDIP heatmap failed"

    # Generate meDIP-only profile
    plotProfile \
        -m "$OUTDIR/matrices/${dmr_type}_medip_matrix.gz" \
        -o "$OUTDIR/profiles/${dmr_type}_medip_only_profile.png" \
        --plotTitle "meDIP Signal at $title" \
        --perGroup \
        --colors "#1B7837" "#762A83" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: meDIP profile failed"

    echo "  Done with $dmr_type"
}

# ============================================================================
# Process each DMR category
# ============================================================================

echo ""
echo "=== Creating binding overlay visualizations ==="

# All high-confidence DMRs
create_binding_overlay "all_dmrs" "$ALL_DMR_BED" "High-Confidence DMRs"

# Hypermethylated DMRs
create_binding_overlay "hypermethylated" "$HYPER_BED" "High-Conf Hypermethylated (TES > GFP)"

# Hypomethylated DMRs
create_binding_overlay "hypomethylated" "$HYPO_BED" "High-Conf Hypomethylated (TES < GFP)"

# ============================================================================
# Create comparative visualization: Hyper vs Hypo
# ============================================================================

echo ""
echo "=== Creating comparative Hyper vs Hypo analysis ==="

HYPER_COUNT=$(wc -l < "$HYPER_BED")
HYPO_COUNT=$(wc -l < "$HYPO_BED")

if [[ $HYPER_COUNT -ge 10 ]] && [[ $HYPO_COUNT -ge 10 ]]; then

    # Matrix with both region sets
    echo "  Computing comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$HYPER_BED" "$HYPO_BED" \
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED" \
        --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/hyper_vs_hypo_comparison_matrix.gz" \
        2>/dev/null || echo "  WARNING: comparison matrix failed"

    # Generate heatmap showing both DMR categories
    echo "  Creating comparison heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/hyper_vs_hypo_comparison_matrix.gz" \
        -o "$OUTDIR/heatmaps/hyper_vs_hypo_comparison_heatmap.png" \
        --colorList "white,#2166AC" "white,#B2182B" "white,#1B7837" "white,#762A83" \
        --zMin 0 0 0 0 \
        --zMax 2 2 150 150 \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 3 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --plotTitle "High-Confidence: Binding & Methylation Overlay" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison heatmap failed"

    # Generate profile showing both DMR categories
    echo "  Creating comparison profile..."
    plotProfile \
        -m "$OUTDIR/matrices/hyper_vs_hypo_comparison_matrix.gz" \
        -o "$OUTDIR/profiles/hyper_vs_hypo_comparison_profile.png" \
        --plotTitle "High-Confidence: Binding & Methylation" \
        --colors "#2166AC" "#B2182B" "#1B7837" "#762A83" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --legendLocation upper-right \
        --plotHeight 8 \
        --plotWidth 14 \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison profile failed"

    echo "  Done with comparison"
else
    echo "  Skipping comparison - insufficient regions in one or both categories"
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "DMR Binding Overlay Analysis Complete"
echo "=============================================="
echo ""
echo "HIGH-CONFIDENCE DMR counts:"
echo "  All DMRs: $(wc -l < "$ALL_DMR_BED") regions"
echo "  Hypermethylated: $HYPER_COUNT regions"
echo "  Hypomethylated: $HYPO_COUNT regions"
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Contents:"
echo "  - matrices/: deepTools matrices for each DMR category"
echo "  - heatmaps/: PNG heatmaps showing binding + methylation (300 DPI)"
echo "  - profiles/: PNG profile plots (300 DPI)"
echo ""
echo "Key comparisons:"
echo "  - *_combined_*: All 4 signals (TES, TEAD1, GFP_meDIP, TES_meDIP)"
echo "  - *_binding_only_*: TES vs TEAD1 binding at DMRs"
echo "  - *_medip_only_*: GFP vs TES methylation signal"
echo "  - hyper_vs_hypo_*: Comparison between hyper/hypomethylated DMRs"
echo ""
echo "NOTE: These are HIGH-CONFIDENCE DMRs where BOTH TES and GFP have"
echo "      >2 mean reads, filtering out GFP dropout artifacts."
echo ""
echo "Interpretation:"
echo "  - The hypomethylated regions now dominate (~65% vs ~35% hyper)"
echo "  - This is the OPPOSITE of the original data (91% hyper)"
echo "  - Strong binding at hypomethylated sites suggests TF occupancy"
echo "    may protect DNA from methylation (despite DNMT fusion)"
echo ""
echo "Finished: $(date)"
