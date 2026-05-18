#!/bin/bash
#SBATCH --job-name=dmr_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/17_dmr_stratified_heatmaps.out
#SBATCH --error=logs/17_dmr_stratified_heatmaps.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# DMR Stratified Heatmaps
# ============================================================================
# Purpose: Create heatmaps for different DMR categories:
#   1. All DMRs (FDR<0.05, |FC|>2)
#   2. Stringent DMRs (FDR<0.01, |FC|>4)
#   3. Hypermethylated only
#   4. Hypomethylated only
#
# Shows meDIP signal (GFP and TES) at DMR regions
# Generates both individual replicate and combined replicate heatmaps
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/17_dmr_heatmaps"
mkdir -p "$OUTDIR"/{matrices,heatmaps,beds}
mkdir -p logs

echo "=============================================="
echo "DMR Stratified Heatmaps"
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

# Define input files
DMR_CSV="results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"
HYPER_BED="results/07_differential_MEDIPS/TES_vs_GFP_hypermethylated.bed"
HYPO_BED="results/07_differential_MEDIPS/TES_vs_GFP_hypomethylated.bed"

# BigWig files - individual replicates
GFP_BW1="results/05_bigwig/GFP-1-IP_RPKM.bw"
GFP_BW2="results/05_bigwig/GFP-2-IP_RPKM.bw"
TES_BW1="results/05_bigwig/TES-1-IP_RPKM.bw"
TES_BW2="results/05_bigwig/TES-2-IP_RPKM.bw"

# Combined BigWig files - use pre-computed correct averages
# NOTE: bigWigMerge produces incorrect values, so we use the properly averaged files
# Created by 05b_combine_bigwig_replicates.sh using bigwigAverage
GFP_COMBINED="results/05_bigwig/GFP_average.bw"
TES_COMBINED="results/05_bigwig/TES_average.bw"

# Chromosome sizes file (kept for reference, no longer used for merging)
CHROM_SIZES="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38.chrom.sizes"

# Verify input files exist
# for f in "$DMR_CSV" "$HYPER_BED" "$HYPO_BED" "$GFP_BW1" "$GFP_BW2" "$TES_BW1" "$TES_BW2"; do
#     if [[ ! -f "$f" ]]; then
#         echo "ERROR: Required file not found: $f"
#         exit 1
#     fi
# done

# echo "Input files verified."

# ============================================================================
# Create combined BigWig files (merge replicates)
# ============================================================================
echo ""
echo "=== Creating combined replicate BigWig files ==="

# Function to merge BigWig files by averaging
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

# Verify pre-computed average BigWig files exist (no longer using buggy bigWigMerge)
if [[ -f "$GFP_COMBINED" ]]; then
    echo "  GFP average BigWig found: $GFP_COMBINED"
else
    echo "  ERROR: GFP average BigWig not found: $GFP_COMBINED"
    echo "  Please create it using: bigwigAverage -b $GFP_BW1 $GFP_BW2 -o $GFP_COMBINED"
    exit 1
fi
if [[ -f "$TES_COMBINED" ]]; then
    echo "  TES average BigWig found: $TES_COMBINED"
else
    echo "  ERROR: TES average BigWig not found: $TES_COMBINED"
    echo "  Please create it using: bigwigAverage -b $TES_BW1 $TES_BW2 -o $TES_COMBINED"
    exit 1
fi

# ============================================================================
# Create BED files for different DMR categories
# ============================================================================
echo ""
echo "=== Creating DMR category BED files ==="

# 1. All DMRs (already in CSV, convert to BED)
echo "Creating all DMRs BED..."
tail -n +2 "$DMR_CSV" | awk -F',' '{
    gsub(/"/, "", $1);  # Remove quotes from chr
    print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7"\t."
}' > "$OUTDIR/beds/all_dmrs.bed"
ALL_DMR_COUNT=$(wc -l < "$OUTDIR/beds/all_dmrs.bed")
echo "  All DMRs: $ALL_DMR_COUNT"

# 2. Stringent DMRs (FDR<0.01, |FC|>4 i.e. |logFC|>2)
echo "Creating stringent DMRs BED..."
tail -n +2 "$DMR_CSV" | awk -F',' '{
    gsub(/"/, "", $1);
    fdr = $10+0;
    logfc = $7+0;
    if (fdr < 0.01 && (logfc > 2 || logfc < -2)) {
        print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"logfc"\t."
    }
}' > "$OUTDIR/beds/stringent_dmrs.bed"
STRINGENT_COUNT=$(wc -l < "$OUTDIR/beds/stringent_dmrs.bed")
echo "  Stringent DMRs (FDR<0.01, |FC|>4): $STRINGENT_COUNT"

# 3. Copy hypermethylated and hypomethylated BEDs
cp "$HYPER_BED" "$OUTDIR/beds/hypermethylated_dmrs.bed"
cp "$HYPO_BED" "$OUTDIR/beds/hypomethylated_dmrs.bed"
HYPER_COUNT=$(wc -l < "$OUTDIR/beds/hypermethylated_dmrs.bed")
HYPO_COUNT=$(wc -l < "$OUTDIR/beds/hypomethylated_dmrs.bed")
echo "  Hypermethylated: $HYPER_COUNT"
echo "  Hypomethylated: $HYPO_COUNT"

# ============================================================================
# Compute matrices and generate heatmaps for each category
# ============================================================================

# Function to create heatmap for a DMR category
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

# Function to create heatmap for a DMR category using combined replicates
create_dmr_heatmap_combined() {
    local category=$1
    local bed_file=$2
    local title=$3

    local region_count=$(wc -l < "$bed_file")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $category (combined) - too few regions ($region_count)"
        return
    fi

    # Check if combined BigWig files exist
    if [[ ! -f "$GFP_COMBINED" ]] || [[ ! -f "$TES_COMBINED" ]]; then
        echo "  Skipping $category (combined) - combined BigWig files not available"
        return
    fi

    echo ""
    echo "Processing $category - Combined Replicates ($region_count regions)..."

    # Compute matrix - center on DMR, extend 5kb each side
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

    # Generate heatmap - sorted by mean signal
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

echo ""
echo "=== Generating heatmaps for each DMR category ==="

# Generate heatmaps for all categories - Individual Replicates
echo ""
echo "--- Individual Replicate Heatmaps ---"
create_dmr_heatmap "all_dmrs" "$OUTDIR/beds/all_dmrs.bed" "All DMRs (TES vs GFP, FDR<0.05, |FC|>2)"
create_dmr_heatmap "stringent_dmrs" "$OUTDIR/beds/stringent_dmrs.bed" "Stringent DMRs (FDR<0.01, |FC|>4)"
create_dmr_heatmap "hypermethylated" "$OUTDIR/beds/hypermethylated_dmrs.bed" "Hypermethylated DMRs (TES > GFP)"
create_dmr_heatmap "hypomethylated" "$OUTDIR/beds/hypomethylated_dmrs.bed" "Hypomethylated DMRs (TES < GFP)"

# Generate heatmaps for all categories - Combined Replicates
echo ""
echo "=== Generating heatmaps for COMBINED replicates ==="
echo "--- Combined Replicate Heatmaps ---"
create_dmr_heatmap_combined "all_dmrs" "$OUTDIR/beds/all_dmrs.bed" "All DMRs (TES vs GFP, FDR<0.05, |FC|>2)"
create_dmr_heatmap_combined "stringent_dmrs" "$OUTDIR/beds/stringent_dmrs.bed" "Stringent DMRs (FDR<0.01, |FC|>4)"
create_dmr_heatmap_combined "hypermethylated" "$OUTDIR/beds/hypermethylated_dmrs.bed" "Hypermethylated DMRs (TES > GFP)"
create_dmr_heatmap_combined "hypomethylated" "$OUTDIR/beds/hypomethylated_dmrs.bed" "Hypomethylated DMRs (TES < GFP)"

# ============================================================================
# Create combined side-by-side comparison
# ============================================================================
echo ""
echo "=== Creating side-by-side comparison ==="

# Only create if both hyper and hypo have enough regions
if [[ -f "$OUTDIR/matrices/hypermethylated_matrix.gz" ]] && [[ -f "$OUTDIR/matrices/hypomethylated_matrix.gz" ]]; then
    echo "Creating comparison plots (individual replicates)..."

    # Create overlaid profile comparison using individual replicates
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$OUTDIR/beds/hypermethylated_dmrs.bed" "$OUTDIR/beds/hypomethylated_dmrs.bed" \
        -S "$GFP_BW1" "$TES_BW1" \
        --samplesLabel "GFP" "TES" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/hyper_vs_hypo_matrix.gz" \
        2>/dev/null

    plotProfile \
        -m "$OUTDIR/matrices/hyper_vs_hypo_matrix.gz" \
        -o "$OUTDIR/heatmaps/hyper_vs_hypo_comparison.png" \
        --plotTitle "Methylation at Hyper vs Hypo DMRs" \
        --perGroup \
        --colors "#2166AC" "#B2182B" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --legendLocation upper-right \
        --refPointLabel "DMR" \
        --dpi 300 \
        2>/dev/null

    echo "  Done with individual replicate comparison"
fi

# Create combined replicate comparison
if [[ -f "$GFP_COMBINED" ]] && [[ -f "$TES_COMBINED" ]]; then
    HYPER_COUNT=$(wc -l < "$OUTDIR/beds/hypermethylated_dmrs.bed")
    HYPO_COUNT=$(wc -l < "$OUTDIR/beds/hypomethylated_dmrs.bed")

    if [[ $HYPER_COUNT -ge 10 ]] && [[ $HYPO_COUNT -ge 10 ]]; then
        echo "Creating comparison plots (combined replicates)..."

        # Create overlaid profile comparison using combined replicates
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
            --plotTitle "Methylation at Hyper vs Hypo DMRs (Combined Replicates)" \
            --perGroup \
            --colors "#2166AC" "#B2182B" \
            --regionsLabel "Hypermethylated" "Hypomethylated" \
            --legendLocation upper-right \
            --refPointLabel "DMR" \
            --dpi 300 \
            2>/dev/null

        # Also create a heatmap for the combined comparison
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
            --plotTitle "Hyper vs Hypo DMRs (Combined Replicates)" \
            --legendLocation best \
            --regionsLabel "Hypermethylated" "Hypomethylated" \
            --dpi 300 \
            2>/dev/null

        echo "  Done with combined replicate comparison"
    fi
fi

# ============================================================================
# Generate summary
# ============================================================================
echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="
echo ""
echo "DMR Categories:"
echo "  All DMRs: $ALL_DMR_COUNT regions"
echo "  Stringent DMRs: $STRINGENT_COUNT regions"
echo "  Hypermethylated: $HYPER_COUNT regions"
echo "  Hypomethylated: $HYPO_COUNT regions"
echo ""
echo "Combined BigWig files:"
if [[ -f "$GFP_COMBINED" ]]; then
    echo "  GFP combined: $GFP_COMBINED"
else
    echo "  GFP combined: NOT CREATED"
fi
if [[ -f "$TES_COMBINED" ]]; then
    echo "  TES combined: $TES_COMBINED"
else
    echo "  TES combined: NOT CREATED"
fi
echo ""
echo "Output files in: $OUTDIR/"
echo "  - beds/: BED files for each category"
echo "  - matrices/: deepTools matrices (individual + combined replicates)"
echo "  - heatmaps/: PNG heatmaps and profiles (300 DPI)"
echo ""
echo "Heatmap types generated (PNG, 300 DPI):"
echo "  - *_heatmap.png: Individual replicates (GFP-1, GFP-2, TES-1, TES-2)"
echo "  - *_combined_heatmap.png: Combined replicates (GFP avg, TES avg)"
echo "  - *_profile.png: Signal profiles (individual replicates)"
echo "  - *_combined_profile.png: Signal profiles (combined replicates)"
echo "  - hyper_vs_hypo_*.png: Comparison of hyper vs hypo DMRs"
echo ""
echo "Publication-ready features:"
echo "  - 300 DPI resolution for all plots"
echo "  - ColorBrewer color palette for accessibility"
echo "  - Clean x-axis labels (-5kb, DMR, +5kb) to prevent overlap"
echo ""
echo "Finished: $(date)"
