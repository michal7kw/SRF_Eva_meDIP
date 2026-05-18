#!/bin/bash
#SBATCH --job-name=dmr_binding_overlay
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/18_dmr_binding_overlay.out
#SBATCH --error=logs/18_dmr_binding_overlay.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# DMR Binding Overlay Analysis
# ============================================================================
# Purpose: Visualize TES/TEAD1 binding signals at DMR regions
# Overlay:
#   - TES Cut&Tag signal
#   - TEAD1 Cut&Tag signal
#   - meDIP GFP signal
#   - meDIP TES signal
#
# This allows assessment of whether TES/TEAD1 binding correlates with
# methylation changes at differentially methylated regions.
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/18_binding_overlay"
mkdir -p "$OUTDIR"/{matrices,heatmaps,profiles}
mkdir -p logs

echo "=============================================="
echo "DMR Binding Overlay Analysis"
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

# DMR BED files from stratified heatmaps (script 17) or differential analysis
DMR_DIR="results/17_dmr_heatmaps/beds"
DMR_DIR_ORIG="results/07_differential_MEDIPS"

# Check if stratified BEDs exist, otherwise use original DMR files
if [[ -d "$DMR_DIR" ]]; then
    HYPER_BED="$DMR_DIR/hypermethylated_dmrs.bed"
    HYPO_BED="$DMR_DIR/hypomethylated_dmrs.bed"
    ALL_DMR_BED="$DMR_DIR/all_dmrs.bed"
    STRINGENT_BED="$DMR_DIR/stringent_dmrs.bed"
else
    HYPER_BED="$DMR_DIR_ORIG/TES_vs_GFP_hypermethylated.bed"
    HYPO_BED="$DMR_DIR_ORIG/TES_vs_GFP_hypomethylated.bed"
fi

# Cut&Tag BigWig files (use combined/average for cleaner signal)
CUTNTAG_DIR="../SRF_Eva_CUTandTAG/results/06_bigwig"
TES_BW="$CUTNTAG_DIR/TES_comb.bw"
TEAD1_BW="$CUTNTAG_DIR/TEAD1_comb.bw"

# Alternatively use average bigwigs from combined replicates
CUTNTAG_AVG_DIR="../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/bigwig"
if [[ -f "$CUTNTAG_AVG_DIR/TES_average.bw" ]]; then
    TES_AVG_BW="$CUTNTAG_AVG_DIR/TES_average.bw"
    TEAD1_AVG_BW="$CUTNTAG_AVG_DIR/TEAD1_average.bw"
else
    TES_AVG_BW="$TES_BW"
    TEAD1_AVG_BW="$TEAD1_BW"
fi

# meDIP BigWig files - individual replicates
MEDIP_DIR="results/05_bigwig"
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
GFP_MEDIP_BW2="$MEDIP_DIR/GFP-2-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"
TES_MEDIP_BW2="$MEDIP_DIR/TES-2-IP_RPKM.bw"

# meDIP BigWig files - combined replicates (pre-computed correct averages)
# NOTE: bigWigMerge produces incorrect values, so we use properly averaged files
# Created by 05b_combine_bigwig_replicates.sh using bigwigAverage
GFP_MEDIP_COMBINED="$MEDIP_DIR/GFP_average.bw"
TES_MEDIP_COMBINED="$MEDIP_DIR/TES_average.bw"

# Chromosome sizes file (kept for reference, no longer used for merging)
CHROM_SIZES="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38.chrom.sizes"

# Verify input files
echo ""
echo "=== Checking input files ==="

# Check DMR files
for f in "$HYPER_BED" "$HYPO_BED"; do
    if [[ -f "$f" ]]; then
        count=$(wc -l < "$f")
        echo "  Found: $f ($count regions)"
    else
        echo "  WARNING: Not found: $f"
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

# Check meDIP BigWigs (individual)
for f in "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  WARNING: Not found: $f"
    fi
done

# Check meDIP BigWigs (combined)
echo ""
echo "=== Checking combined meDIP BigWig files ==="
for f in "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  Not found: $f (will be created)"
    fi
done

# ============================================================================
# Create combined BigWig files if they don't exist (merge replicates)
# ============================================================================

merge_bigwigs() {
    local output=$1
    local input1=$2
    local input2=$3
    local name=$4

    if [[ -f "$output" ]]; then
        echo "  $name combined BigWig already exists, skipping..."
    else
        echo "  Creating $name combined BigWig (averaging replicates)..."

        # Use bigWigMerge to merge, then convert back to BigWig
        local tmp_bg="${output%.bw}.tmp.bg"

        # bigWigMerge sums the values, so we need to divide by 2 for average
        bigWigMerge "$input1" "$input2" "$tmp_bg"

        # Divide values by 2 to get average (using awk)
        awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4/2}' "$tmp_bg" > "${tmp_bg}.avg"
        mv "${tmp_bg}.avg" "$tmp_bg"

        # Sort the bedGraph file
        sort -k1,1 -k2,2n "$tmp_bg" > "${tmp_bg}.sorted"
        mv "${tmp_bg}.sorted" "$tmp_bg"

        # Convert to BigWig
        if [[ -f "$CHROM_SIZES" ]]; then
            bedGraphToBigWig "$tmp_bg" "$CHROM_SIZES" "$output"
            rm -f "$tmp_bg"
            echo "    Created: $output"
        else
            echo "  WARNING: Chromosome sizes file not found at $CHROM_SIZES"
            echo "           Attempting to fetch from UCSC..."
            fetchChromSizes hg38 > "$CHROM_SIZES" 2>/dev/null || true
            if [[ -f "$CHROM_SIZES" ]]; then
                bedGraphToBigWig "$tmp_bg" "$CHROM_SIZES" "$output"
                rm -f "$tmp_bg"
                echo "    Created: $output"
            else
                echo "  ERROR: Could not create combined BigWig - missing chrom.sizes"
                rm -f "$tmp_bg"
            fi
        fi
    fi
}

echo ""
echo "=== Verifying pre-computed average meDIP BigWig files ==="
# NOTE: No longer using buggy bigWigMerge - using pre-computed correct averages
if [[ -f "$GFP_MEDIP_COMBINED" ]]; then
    echo "  GFP average meDIP BigWig found: $GFP_MEDIP_COMBINED"
else
    echo "  ERROR: GFP average meDIP BigWig not found: $GFP_MEDIP_COMBINED"
    echo "  Please create it using: bigwigAverage -b $GFP_MEDIP_BW1 $GFP_MEDIP_BW2 -o $GFP_MEDIP_COMBINED"
    exit 1
fi
if [[ -f "$TES_MEDIP_COMBINED" ]]; then
    echo "  TES average meDIP BigWig found: $TES_MEDIP_COMBINED"
else
    echo "  ERROR: TES average meDIP BigWig not found: $TES_MEDIP_COMBINED"
    echo "  Please create it using: bigwigAverage -b $TES_MEDIP_BW1 $TES_MEDIP_BW2 -o $TES_MEDIP_COMBINED"
    exit 1
fi

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
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
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
    # NOTE: Cut&Tag signal (~0-2) is ~100x lower than meDIP signal (~0-200)
    # Use per-sample zMin/zMax to show both on appropriate scales
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

    # Generate combined profile with per-sample panels (each with own Y-axis)
    echo "  Creating combined profile (per-sample panels)..."
    plotProfile \
        -m "$OUTDIR/matrices/${dmr_type}_combined_matrix.gz" \
        -o "$OUTDIR/profiles/${dmr_type}_combined_profile.png" \
        --plotTitle "Binding & Methylation at $title" \
        --perSample \
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
    # Matrix 3: meDIP signals only (GFP vs TES - both replicates)
    # -------------------------------------------------------------------------
    echo "  Computing meDIP-only matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
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
        --colorMap YlGn YlGn Purples Purples \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 2.5 \
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
        --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: meDIP profile failed"

    echo "  Done with $dmr_type"
}

# ============================================================================
# Function to create binding overlay heatmap with COMBINED replicates
# ============================================================================

create_binding_overlay_combined() {
    local dmr_type=$1
    local bed_file=$2
    local title=$3

    if [[ ! -f "$bed_file" ]]; then
        echo "  Skipping $dmr_type (combined) - BED file not found: $bed_file"
        return
    fi

    # Check if combined BigWig files exist
    if [[ ! -f "$GFP_MEDIP_COMBINED" ]] || [[ ! -f "$TES_MEDIP_COMBINED" ]]; then
        echo "  Skipping $dmr_type (combined) - combined meDIP BigWig files not available"
        return
    fi

    local region_count=$(wc -l < "$bed_file")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $dmr_type (combined) - too few regions ($region_count)"
        return
    fi

    echo ""
    echo "=== Processing $dmr_type - COMBINED REPLICATES ($region_count regions) ==="

    # -------------------------------------------------------------------------
    # Matrix 1: Cut&Tag binding + combined meDIP methylation
    # -------------------------------------------------------------------------
    echo "  Computing combined matrix (Cut&Tag + combined meDIP)..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED" \
        --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${dmr_type}_combined_reps_matrix.gz" \
        --outFileNameMatrix "$OUTDIR/matrices/${dmr_type}_combined_reps_matrix.tab" \
        2>/dev/null || {
            echo "  WARNING: computeMatrix failed for $dmr_type combined reps matrix"
            return
        }

    # Generate heatmap with combined replicates
    # NOTE: Cut&Tag signal (~0-2) is ~100x lower than meDIP signal (~0-200)
    # Use per-sample zMin/zMax to show both on appropriate scales
    echo "  Creating combined reps heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${dmr_type}_combined_reps_matrix.gz" \
        -o "$OUTDIR/heatmaps/${dmr_type}_binding_methylation_combined_heatmap.png" \
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
        --plotTitle "$title (Combined Replicates)" \
        --legendLocation upper-right \
        --dpi 300 \
        --outFileSortedRegions "$OUTDIR/heatmaps/${dmr_type}_combined_sorted.bed" \
        2>/dev/null || echo "  WARNING: plotHeatmap failed"

    # Generate combined profile with per-sample panels (each with own Y-axis)
    echo "  Creating combined profile (per-sample panels, combined reps)..."
    plotProfile \
        -m "$OUTDIR/matrices/${dmr_type}_combined_reps_matrix.gz" \
        -o "$OUTDIR/profiles/${dmr_type}_combined_reps_profile.png" \
        --plotTitle "Binding & Methylation at $title (Combined Replicates)" \
        --perSample \
        --colors "#2166AC" "#B2182B" "#1B7837" "#762A83" \
        --legendLocation upper-right \
        --plotHeight 8 \
        --plotWidth 14 \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: combined reps profile failed"

    # -------------------------------------------------------------------------
    # Matrix 2: Combined meDIP signals only (GFP vs TES)
    # -------------------------------------------------------------------------
    echo "  Computing meDIP-only combined matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED" \
        --samplesLabel "GFP (combined)" "TES (combined)" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${dmr_type}_medip_combined_matrix.gz" \
        2>/dev/null || {
            echo "  WARNING: computeMatrix failed for $dmr_type meDIP combined matrix"
            return
        }

    # Generate meDIP-only combined heatmap
    echo "  Creating meDIP combined heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${dmr_type}_medip_combined_matrix.gz" \
        -o "$OUTDIR/heatmaps/${dmr_type}_medip_combined_heatmap.png" \
        --colorMap YlGn Purples \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --xAxisLabel "" \
        --refPointLabel "DMR" \
        --plotTitle "meDIP Signal at $title (Combined Replicates)" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: meDIP combined heatmap failed"

    # Generate meDIP-only combined profile
    plotProfile \
        -m "$OUTDIR/matrices/${dmr_type}_medip_combined_matrix.gz" \
        -o "$OUTDIR/profiles/${dmr_type}_medip_combined_profile.png" \
        --plotTitle "meDIP Signal at $title (Combined Replicates)" \
        --perGroup \
        --colors "#1B7837" "#762A83" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: meDIP combined profile failed"

    echo "  Done with $dmr_type (combined replicates)"
}

# ============================================================================
# Process each DMR category
# ============================================================================

echo ""
echo "=== Creating binding overlay visualizations ==="

# Hypermethylated DMRs
create_binding_overlay "hypermethylated" "$HYPER_BED" "Hypermethylated DMRs (TES > GFP)"

# Hypomethylated DMRs
create_binding_overlay "hypomethylated" "$HYPO_BED" "Hypomethylated DMRs (TES < GFP)"

# All DMRs (if available)
if [[ -f "$ALL_DMR_BED" ]]; then
    create_binding_overlay "all_dmrs" "$ALL_DMR_BED" "All DMRs (FDR<0.05, |FC|>2)"
fi

# Stringent DMRs (if available)
if [[ -f "$STRINGENT_BED" ]]; then
    create_binding_overlay "stringent" "$STRINGENT_BED" "Stringent DMRs (FDR<0.01, |FC|>4)"
fi

# ============================================================================
# Process each DMR category - COMBINED REPLICATES
# ============================================================================

echo ""
echo "=== Creating binding overlay visualizations - COMBINED REPLICATES ==="

# Hypermethylated DMRs (combined)
create_binding_overlay_combined "hypermethylated" "$HYPER_BED" "Hypermethylated DMRs (TES > GFP)"

# Hypomethylated DMRs (combined)
create_binding_overlay_combined "hypomethylated" "$HYPO_BED" "Hypomethylated DMRs (TES < GFP)"

# All DMRs (combined, if available)
if [[ -f "$ALL_DMR_BED" ]]; then
    create_binding_overlay_combined "all_dmrs" "$ALL_DMR_BED" "All DMRs (FDR<0.05, |FC|>2)"
fi

# Stringent DMRs (combined, if available)
if [[ -f "$STRINGENT_BED" ]]; then
    create_binding_overlay_combined "stringent" "$STRINGENT_BED" "Stringent DMRs (FDR<0.01, |FC|>4)"
fi

# ============================================================================
# Create comparative visualization: Hyper vs Hypo
# ============================================================================

echo ""
echo "=== Creating comparative Hyper vs Hypo analysis ==="

if [[ -f "$HYPER_BED" ]] && [[ -f "$HYPO_BED" ]]; then

    # Matrix with both region sets
    echo "  Computing comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$HYPER_BED" "$HYPO_BED" \
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
        --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/hyper_vs_hypo_comparison_matrix.gz" \
        2>/dev/null || echo "  WARNING: comparison matrix failed"

    # Generate heatmap showing both DMR categories
    # NOTE: Cut&Tag signal (~0-2) is ~100x lower than meDIP signal (~0-200)
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
        --plotTitle "Binding & Methylation: Hyper vs Hypo DMRs" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison heatmap failed"

    # Generate profile showing both DMR categories
    # Use --perSample to give each sample its own Y-axis (auto-scaled)
    echo "  Creating comparison profile (per-sample panels)..."
    plotProfile \
        -m "$OUTDIR/matrices/hyper_vs_hypo_comparison_matrix.gz" \
        -o "$OUTDIR/profiles/hyper_vs_hypo_comparison_profile.png" \
        --plotTitle "Binding & Methylation: Hyper vs Hypo DMRs" \
        --perSample \
        --colors "#2166AC" "#B2182B" "#1B7837" "#762A83" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --legendLocation upper-right \
        --plotHeight 8 \
        --plotWidth 14 \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison profile failed"

    echo "  Done with comparison"
fi

# ============================================================================
# Create comparative visualization: Hyper vs Hypo - COMBINED REPLICATES
# ============================================================================

echo ""
echo "=== Creating comparative Hyper vs Hypo analysis - COMBINED REPLICATES ==="

if [[ -f "$HYPER_BED" ]] && [[ -f "$HYPO_BED" ]] && [[ -f "$GFP_MEDIP_COMBINED" ]] && [[ -f "$TES_MEDIP_COMBINED" ]]; then

    # Matrix with both region sets using combined meDIP
    echo "  Computing comparison matrix (combined meDIP)..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$HYPER_BED" "$HYPO_BED" \
        -S "$TES_AVG_BW" "$TEAD1_AVG_BW" "$GFP_MEDIP_COMBINED" "$TES_MEDIP_COMBINED" \
        --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/hyper_vs_hypo_combined_comparison_matrix.gz" \
        2>/dev/null || echo "  WARNING: comparison matrix failed"

    # Generate heatmap showing both DMR categories
    # NOTE: Cut&Tag signal (~0-2) is ~100x lower than meDIP signal (~0-200)
    echo "  Creating comparison heatmap (combined)..."
    plotHeatmap \
        -m "$OUTDIR/matrices/hyper_vs_hypo_combined_comparison_matrix.gz" \
        -o "$OUTDIR/heatmaps/hyper_vs_hypo_combined_comparison_heatmap.png" \
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
        --plotTitle "Binding & Methylation: Hyper vs Hypo (Combined Replicates)" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison heatmap failed"

    # Generate profile showing both DMR categories
    # Use --perSample to give each sample its own Y-axis (auto-scaled)
    echo "  Creating comparison profile (combined, per-sample panels)..."
    plotProfile \
        -m "$OUTDIR/matrices/hyper_vs_hypo_combined_comparison_matrix.gz" \
        -o "$OUTDIR/profiles/hyper_vs_hypo_combined_comparison_profile.png" \
        --plotTitle "Binding & Methylation: Hyper vs Hypo (Combined Replicates)" \
        --perSample \
        --colors "#2166AC" "#B2182B" "#1B7837" "#762A83" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --legendLocation upper-right \
        --plotHeight 8 \
        --plotWidth 14 \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison profile failed"

    echo "  Done with comparison (combined replicates)"
else
    echo "  Skipping combined comparison - missing required files"
fi

# ============================================================================
# Create R-based side-by-side multi-panel figures
# ============================================================================

echo ""
echo "=== Creating R-based side-by-side figures ==="

# Activate R environment
conda activate seurat_full2

# Create R script for side-by-side visualization
cat > "$OUTDIR/create_sidebyside_profiles.R" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(grid)
    library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

cat("Creating publication-quality side-by-side profile figures...\n")

# Publication-ready color palette (ColorBrewer)
COLORS_BINDING <- c("#2166AC", "#B2182B")  # Blue and Red
COLORS_MEDIP <- c("#A6D96A", "#1B7837", "#DFC27D", "#762A83")  # Light green, dark green, tan, purple

# Publication-ready theme
theme_publication <- function(base_size = 12) {
    theme_bw(base_size = base_size) +
    theme(
        # Title
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        # Axis
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size - 1, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        # Legend
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = base_size - 1),
        legend.key.size = unit(0.8, "lines"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.margin = margin(t = 0, b = 5),
        # Panel
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        # Margins
        plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
}

# Function to read deepTools matrix and create profile data
read_deeptools_matrix <- function(matrix_file) {
    if (!file.exists(matrix_file)) {
        return(NULL)
    }

    # Read the gzipped matrix
    con <- gzfile(matrix_file, "r")
    header <- readLines(con, n = 1)
    close(con)

    # Parse header for sample labels and parameters
    header_json <- gsub("^@", "", header)
    header_info <- jsonlite::fromJSON(header_json)

    # Read the data
    mat_data <- read.table(gzfile(matrix_file), skip = 1, header = FALSE)

    # Extract sample labels
    sample_labels <- header_info$sample_labels

    # Calculate number of bins
    n_bins <- (ncol(mat_data) - 6) / length(sample_labels)

    # Extract signal columns for each sample
    profiles <- list()
    for (i in seq_along(sample_labels)) {
        start_col <- 7 + (i - 1) * n_bins
        end_col <- start_col + n_bins - 1
        sample_mat <- mat_data[, start_col:end_col]
        profiles[[sample_labels[i]]] <- colMeans(sample_mat, na.rm = TRUE)
    }

    # Create position vector (assuming -5000 to +5000 bp)
    positions <- seq(-5000, 5000, length.out = n_bins)

    # Convert to data frame for ggplot
    df <- data.frame(Position = positions)
    for (name in names(profiles)) {
        df[[name]] <- profiles[[name]]
    }

    return(df)
}

# Function to create publication-quality profile plot
create_profile_plot <- function(df, samples, colors, title, ylab = "Signal") {
    if (is.null(df)) return(NULL)

    # Melt data for ggplot
    df_long <- melt(df, id.vars = "Position",
                    measure.vars = samples,
                    variable.name = "Sample",
                    value.name = "Signal")

    p <- ggplot(df_long, aes(x = Position, y = Signal, color = Sample)) +
        geom_line(linewidth = 1.2) +
        scale_color_manual(values = colors) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
        scale_x_continuous(
            breaks = c(-5000, -2500, 0, 2500, 5000),
            labels = c("-5kb", "-2.5kb", "DMR", "+2.5kb", "+5kb"),
            expand = c(0.02, 0)
        ) +
        scale_y_continuous(expand = c(0.02, 0)) +
        labs(title = title, x = NULL, y = ylab) +
        theme_publication(base_size = 11) +
        guides(color = guide_legend(nrow = 1))

    return(p)
}

# Process each DMR category
dmr_types <- c("hypermethylated", "hypomethylated", "all_dmrs", "stringent")

for (dmr_type in dmr_types) {
    binding_matrix <- file.path(outdir, "matrices", paste0(dmr_type, "_binding_matrix.gz"))
    medip_matrix <- file.path(outdir, "matrices", paste0(dmr_type, "_medip_matrix.gz"))

    if (!file.exists(binding_matrix) || !file.exists(medip_matrix)) {
        cat(paste("  Skipping", dmr_type, "- matrix files not found\n"))
        next
    }

    cat(paste("  Processing", dmr_type, "...\n"))

    tryCatch({
        # Read matrices
        binding_df <- read_deeptools_matrix(binding_matrix)
        medip_df <- read_deeptools_matrix(medip_matrix)

        if (is.null(binding_df) || is.null(medip_df)) {
            cat(paste("    Failed to read matrices for", dmr_type, "\n"))
            next
        }

        # Create binding plot
        binding_samples <- colnames(binding_df)[colnames(binding_df) != "Position"]
        p_binding <- create_profile_plot(
            binding_df,
            binding_samples,
            COLORS_BINDING,
            "TF Binding (Cut&Tag)",
            "Signal"
        )

        # Create meDIP plot
        medip_samples <- colnames(medip_df)[colnames(medip_df) != "Position"]
        p_medip <- create_profile_plot(
            medip_df,
            medip_samples,
            COLORS_MEDIP,
            "DNA Methylation (meDIP)",
            "Signal (RPKM)"
        )

        # Combine side by side - PNG output
        output_file <- file.path(outdir, "profiles", paste0(dmr_type, "_sidebyside_profile.png"))

        png(output_file, width = 14, height = 5, units = "in", res = 300)
        grid.arrange(
            p_binding, p_medip,
            ncol = 2,
            top = textGrob(
                paste0(tools::toTitleCase(gsub("_", " ", dmr_type)), " DMRs: Binding vs Methylation"),
                gp = gpar(fontsize = 14, fontface = "bold")
            )
        )
        dev.off()

        cat(paste("    Saved:", output_file, "\n"))

    }, error = function(e) {
        cat(paste("    Error processing", dmr_type, ":", e$message, "\n"))
    })
}

cat("Done creating side-by-side figures.\n")
RSCRIPT

# Run the R script
Rscript "$OUTDIR/create_sidebyside_profiles.R" "$OUTDIR" 2>/dev/null || echo "  Note: R side-by-side figures may not have been created"

# Clean up temp files
rm -f "$OUTDIR/profiles/"*_temp.png 2>/dev/null

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "DMR Binding Overlay Analysis Complete"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Combined meDIP BigWig files:"
if [[ -f "$GFP_MEDIP_COMBINED" ]]; then
    echo "  GFP combined: $GFP_MEDIP_COMBINED"
else
    echo "  GFP combined: NOT CREATED"
fi
if [[ -f "$TES_MEDIP_COMBINED" ]]; then
    echo "  TES combined: $TES_MEDIP_COMBINED"
else
    echo "  TES combined: NOT CREATED"
fi
echo ""
echo "Contents:"
echo "  - matrices/: deepTools matrices for each DMR category"
echo "  - heatmaps/: PNG heatmaps showing binding + methylation (300 DPI)"
echo "  - profiles/: PNG profile plots (300 DPI)"
echo ""
echo "Key comparisons (Individual Replicates):"
echo "  - *_combined_*: All 4 signals (TES, TEAD1, GFP_meDIP, TES_meDIP)"
echo "  - *_binding_only_*: TES vs TEAD1 binding at DMRs"
echo "  - *_medip_only_*: GFP vs TES methylation signal (all 4 replicates)"
echo "  - hyper_vs_hypo_*: Comparison between hyper/hypomethylated DMRs"
echo ""
echo "Key comparisons (Combined Replicates):"
echo "  - *_combined_reps_*: Binding + combined meDIP signals"
echo "  - *_medip_combined_*: Combined GFP vs TES meDIP signal"
echo "  - *_binding_methylation_combined_*: Full overlay with combined replicates"
echo "  - hyper_vs_hypo_combined_*: Hyper vs Hypo with combined replicates"
echo ""
echo "Profile types generated (PNG, 300 DPI):"
echo "  - *_binding_only_profile.png: TES vs TEAD1 binding (proper scale)"
echo "  - *_medip_only_profile.png: GFP vs TES methylation (proper scale)"
echo "  - *_sidebyside_profile.png: Binding & methylation SIDE BY SIDE (each with own scale)"
echo ""
echo "Heatmap types generated (PNG, 300 DPI):"
echo "  - *_heatmap.png: Individual replicates"
echo "  - *_combined_heatmap.png: Combined replicates (averaged)"
echo ""
echo "Publication-ready features:"
echo "  - 300 DPI resolution for all plots"
echo "  - ColorBrewer color palette for accessibility"
echo "  - Clean x-axis labels (-5kb, DMR, +5kb) to prevent overlap"
echo "  - Consistent styling across all visualizations"
echo ""
echo "Interpretation:"
echo "  - Strong TES/TEAD1 signal at DMRs suggests direct binding"
echo "  - Increased meDIP signal in TES samples confirms hypermethylation"
echo "  - Compare binding density at hyper vs hypo DMRs"
echo "  - Combined replicates reduce noise and show cleaner signal"
echo ""
echo "Finished: $(date)"
