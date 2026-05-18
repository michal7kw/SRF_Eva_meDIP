#!/bin/bash
#SBATCH --job-name=binding_causes_methylation
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/19_binding_causes_methylation.out
#SBATCH --error=logs/19_binding_causes_methylation.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# Does TES Binding Cause Methylation?
# ============================================================================
#
# This script addresses the CAUSAL question: Does TES binding lead to DNA
# methylation at target sites?
#
# APPROACH: Start with TES/TEAD1 BINDING SITES (Cut&Tag peaks) and examine
# whether methylation (meDIP signal) is elevated at these locations.
#
# This is the OPPOSITE of script 18, which asked:
#   "At DMRs, is there TF binding?" (Region selection: meDIP)
#
# Here we ask:
#   "At TF binding sites, is there methylation change?" (Region selection: Cut&Tag)
#
# Key comparisons:
#   1. TES binding sites: Is meDIP signal higher in TES vs GFP?
#   2. TEAD1 binding sites: Is meDIP signal different? (endogenous control)
#   3. TES-unique vs TEAD1-unique vs Shared sites
#   4. Promoter vs non-promoter binding sites
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/19_binding_causes_methylation"
mkdir -p "$OUTDIR"/{beds,matrices,heatmaps,profiles,statistics}
mkdir -p logs

echo "=============================================="
echo "Does TES Binding Cause Methylation?"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment with deepTools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ============================================================================
# Define input files
# ============================================================================

# Cut&Tag peak files (BINDING SITES - these define the regions)
CUTNTAG_PEAKS="../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/consensus_peaks"
TES_PEAKS="$CUTNTAG_PEAKS/TES_consensus_peaks.bed"
TEAD1_PEAKS="$CUTNTAG_PEAKS/TEAD1_consensus_peaks.bed"

# Alternative: use individual replicate peaks if consensus not available
CUTNTAG_PEAKS_ALT="../SRF_Eva_CUTandTAG/results/05_peaks_narrow"

# Cut&Tag BigWig files (for binding signal visualization)
CUTNTAG_BW="../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/bigwig"
TES_BINDING_BW="$CUTNTAG_BW/TES_average.bw"
TEAD1_BINDING_BW="$CUTNTAG_BW/TEAD1_average.bw"

# Fallback BigWig locations
CUTNTAG_BW_ALT="../SRF_Eva_CUTandTAG/results/06_bigwig"

# meDIP BigWig files
MEDIP_DIR="results/05_bigwig"
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
GFP_MEDIP_BW2="$MEDIP_DIR/GFP-2-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"
TES_MEDIP_BW2="$MEDIP_DIR/TES-2-IP_RPKM.bw"

# Combined meDIP BigWigs - use pre-computed correct averages
# NOTE: bigWigMerge produces incorrect values, so we use properly averaged files
# Created by 05b_combine_bigwig_replicates.sh using bigwigAverage
GFP_MEDIP_COMBINED="$MEDIP_DIR/GFP_average.bw"
TES_MEDIP_COMBINED="$MEDIP_DIR/TES_average.bw"

# DMR file for intersection analysis
DMR_FILE="results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv"

# ============================================================================
# Verify and locate input files
# ============================================================================

echo ""
echo "=== Checking input files ==="

# Find TES peaks
if [[ -f "$TES_PEAKS" ]]; then
    echo "  Found TES consensus peaks: $TES_PEAKS"
elif [[ -f "$CUTNTAG_PEAKS_ALT/TES-1_peaks.narrowPeak" ]]; then
    echo "  Using individual TES peaks (will merge)"
    TES_PEAKS="$OUTDIR/beds/TES_merged_peaks.bed"
    cat "$CUTNTAG_PEAKS_ALT"/TES-*_peaks.narrowPeak | \
        cut -f1-3 | sort -k1,1 -k2,2n | \
        bedtools merge -i - > "$TES_PEAKS"
else
    echo "  ERROR: No TES peaks found"
    exit 1
fi

# Find TEAD1 peaks
if [[ -f "$TEAD1_PEAKS" ]]; then
    echo "  Found TEAD1 consensus peaks: $TEAD1_PEAKS"
elif [[ -f "$CUTNTAG_PEAKS_ALT/TEAD1-1_peaks.narrowPeak" ]]; then
    echo "  Using individual TEAD1 peaks (will merge)"
    TEAD1_PEAKS="$OUTDIR/beds/TEAD1_merged_peaks.bed"
    cat "$CUTNTAG_PEAKS_ALT"/TEAD1-*_peaks.narrowPeak | \
        cut -f1-3 | sort -k1,1 -k2,2n | \
        bedtools merge -i - > "$TEAD1_PEAKS"
else
    echo "  ERROR: No TEAD1 peaks found"
    exit 1
fi

# Find binding BigWigs
if [[ -f "$TES_BINDING_BW" ]]; then
    echo "  Found TES binding BigWig: $TES_BINDING_BW"
else
    TES_BINDING_BW="$CUTNTAG_BW_ALT/TES-1.bw"
    if [[ -f "$TES_BINDING_BW" ]]; then
        echo "  Using individual TES BigWig: $TES_BINDING_BW"
    else
        echo "  WARNING: No TES binding BigWig found"
    fi
fi

if [[ -f "$TEAD1_BINDING_BW" ]]; then
    echo "  Found TEAD1 binding BigWig: $TEAD1_BINDING_BW"
else
    TEAD1_BINDING_BW="$CUTNTAG_BW_ALT/TEAD1-1.bw"
    if [[ -f "$TEAD1_BINDING_BW" ]]; then
        echo "  Using individual TEAD1 BigWig: $TEAD1_BINDING_BW"
    else
        echo "  WARNING: No TEAD1 binding BigWig found"
    fi
fi

# Check meDIP files
for f in "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  ERROR: Not found: $f"
        exit 1
    fi
done

# ============================================================================
# Create peak categories
# ============================================================================

echo ""
echo "=== Creating peak categories ==="

# Count input peaks
TES_COUNT=$(wc -l < "$TES_PEAKS")
TEAD1_COUNT=$(wc -l < "$TEAD1_PEAKS")
echo "  TES peaks: $TES_COUNT"
echo "  TEAD1 peaks: $TEAD1_COUNT"

# Create standardized BED files (first 3 columns only)
cut -f1-3 "$TES_PEAKS" | sort -k1,1 -k2,2n > "$OUTDIR/beds/TES_all_peaks.bed"
cut -f1-3 "$TEAD1_PEAKS" | sort -k1,1 -k2,2n > "$OUTDIR/beds/TEAD1_all_peaks.bed"

# TES-unique peaks (TES but not TEAD1)
bedtools intersect -a "$OUTDIR/beds/TES_all_peaks.bed" \
    -b "$OUTDIR/beds/TEAD1_all_peaks.bed" -v > "$OUTDIR/beds/TES_unique_peaks.bed"
TES_UNIQUE=$(wc -l < "$OUTDIR/beds/TES_unique_peaks.bed")
echo "  TES-unique peaks: $TES_UNIQUE"

# TEAD1-unique peaks (TEAD1 but not TES)
bedtools intersect -a "$OUTDIR/beds/TEAD1_all_peaks.bed" \
    -b "$OUTDIR/beds/TES_all_peaks.bed" -v > "$OUTDIR/beds/TEAD1_unique_peaks.bed"
TEAD1_UNIQUE=$(wc -l < "$OUTDIR/beds/TEAD1_unique_peaks.bed")
echo "  TEAD1-unique peaks: $TEAD1_UNIQUE"

# Shared peaks (TES and TEAD1 overlap)
bedtools intersect -a "$OUTDIR/beds/TES_all_peaks.bed" \
    -b "$OUTDIR/beds/TEAD1_all_peaks.bed" -u > "$OUTDIR/beds/shared_peaks.bed"
SHARED=$(wc -l < "$OUTDIR/beds/shared_peaks.bed")
echo "  Shared peaks (TES ∩ TEAD1): $SHARED"

# ============================================================================
# Intersect binding sites with DMRs
# ============================================================================

echo ""
echo "=== Intersecting binding sites with DMRs ==="

# Convert DMR CSV to BED
if [[ -f "$DMR_FILE" ]]; then
    # Skip header, extract chr, start, end, logFC
    tail -n +2 "$DMR_FILE" | \
        awk -F',' '{gsub(/"/, "", $1); gsub(/"/, "", $2); gsub(/"/, "", $3); gsub(/"/, "", $7);
                    print $1"\t"$2"\t"$3"\t"$7}' | \
        sort -k1,1 -k2,2n > "$OUTDIR/beds/all_dmrs.bed"

    DMR_COUNT=$(wc -l < "$OUTDIR/beds/all_dmrs.bed")
    echo "  Total DMRs: $DMR_COUNT"

    # TES peaks overlapping DMRs
    bedtools intersect -a "$OUTDIR/beds/TES_all_peaks.bed" \
        -b "$OUTDIR/beds/all_dmrs.bed" -u > "$OUTDIR/beds/TES_peaks_with_DMR.bed"
    TES_DMR=$(wc -l < "$OUTDIR/beds/TES_peaks_with_DMR.bed")
    echo "  TES peaks overlapping DMRs: $TES_DMR ($(echo "scale=1; 100*$TES_DMR/$TES_COUNT" | bc)%)"

    # TES peaks NOT overlapping DMRs
    bedtools intersect -a "$OUTDIR/beds/TES_all_peaks.bed" \
        -b "$OUTDIR/beds/all_dmrs.bed" -v > "$OUTDIR/beds/TES_peaks_no_DMR.bed"
    TES_NO_DMR=$(wc -l < "$OUTDIR/beds/TES_peaks_no_DMR.bed")
    echo "  TES peaks without DMRs: $TES_NO_DMR ($(echo "scale=1; 100*$TES_NO_DMR/$TES_COUNT" | bc)%)"

    # TEAD1 peaks overlapping DMRs
    bedtools intersect -a "$OUTDIR/beds/TEAD1_all_peaks.bed" \
        -b "$OUTDIR/beds/all_dmrs.bed" -u > "$OUTDIR/beds/TEAD1_peaks_with_DMR.bed"
    TEAD1_DMR=$(wc -l < "$OUTDIR/beds/TEAD1_peaks_with_DMR.bed")
    echo "  TEAD1 peaks overlapping DMRs: $TEAD1_DMR ($(echo "scale=1; 100*$TEAD1_DMR/$TEAD1_COUNT" | bc)%)"

    # TEAD1 peaks NOT overlapping DMRs
    bedtools intersect -a "$OUTDIR/beds/TEAD1_all_peaks.bed" \
        -b "$OUTDIR/beds/all_dmrs.bed" -v > "$OUTDIR/beds/TEAD1_peaks_no_DMR.bed"

    # Hypermethylated DMRs only (logFC > 0)
    awk '$4 > 0' "$OUTDIR/beds/all_dmrs.bed" > "$OUTDIR/beds/hypermethylated_dmrs.bed"
    HYPER_COUNT=$(wc -l < "$OUTDIR/beds/hypermethylated_dmrs.bed")
    echo "  Hypermethylated DMRs: $HYPER_COUNT"

    # TES peaks overlapping HYPER DMRs
    bedtools intersect -a "$OUTDIR/beds/TES_all_peaks.bed" \
        -b "$OUTDIR/beds/hypermethylated_dmrs.bed" -u > "$OUTDIR/beds/TES_peaks_with_hyperDMR.bed"
    TES_HYPER=$(wc -l < "$OUTDIR/beds/TES_peaks_with_hyperDMR.bed")
    echo "  TES peaks overlapping hypermethylated DMRs: $TES_HYPER"

else
    echo "  WARNING: DMR file not found, skipping intersection analysis"
fi

# ============================================================================
# Function to create methylation profile at binding sites
# ============================================================================

create_methylation_at_binding() {
    local peak_type=$1
    local bed_file=$2
    local title=$3

    if [[ ! -f "$bed_file" ]]; then
        echo "  Skipping $peak_type - BED file not found"
        return
    fi

    local region_count=$(wc -l < "$bed_file")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $peak_type - too few regions ($region_count)"
        return
    fi

    echo ""
    echo "=== Processing $peak_type ($region_count regions) ==="

    # -------------------------------------------------------------------------
    # Matrix: meDIP signal at binding sites
    # -------------------------------------------------------------------------
    echo "  Computing meDIP matrix at $peak_type..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${peak_type}_medip_matrix.gz" \
        2>/dev/null || {
            echo "  WARNING: computeMatrix failed for $peak_type"
            return
        }

    # Generate meDIP heatmap
    echo "  Creating meDIP heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${peak_type}_medip_matrix.gz" \
        -o "$OUTDIR/heatmaps/${peak_type}_medip_heatmap.png" \
        --colorMap YlGn YlGn Purples Purples \
        --sortRegions descend \
        --sortUsing mean \
        --sortUsingSamples 3 4 \
        --heatmapHeight 15 \
        --heatmapWidth 2.5 \
        --xAxisLabel "" \
        --refPointLabel "Peak" \
        --plotTitle "meDIP at $title" \
        --legendLocation upper-right \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: heatmap failed"

    # Generate meDIP profile
    echo "  Creating meDIP profile..."
    plotProfile \
        -m "$OUTDIR/matrices/${peak_type}_medip_matrix.gz" \
        -o "$OUTDIR/profiles/${peak_type}_medip_profile.png" \
        --plotTitle "meDIP Signal at $title" \
        --perGroup \
        --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
        --legendLocation upper-right \
        --refPointLabel "Peak" \
        --plotHeight 6 \
        --plotWidth 10 \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: profile failed"

    # -------------------------------------------------------------------------
    # Matrix: Combined binding + meDIP (for context)
    # -------------------------------------------------------------------------
    if [[ -f "$TES_BINDING_BW" ]] && [[ -f "$TEAD1_BINDING_BW" ]]; then
        echo "  Computing combined matrix (binding + meDIP)..."
        computeMatrix reference-point \
            --referencePoint center \
            -b 5000 -a 5000 \
            -R "$bed_file" \
            -S "$TES_BINDING_BW" "$TEAD1_BINDING_BW" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
            --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
            --skipZeros \
            --binSize 50 \
            --numberOfProcessors 8 \
            -o "$OUTDIR/matrices/${peak_type}_combined_matrix.gz" \
            2>/dev/null || echo "  WARNING: combined matrix failed"

        # Generate combined heatmap
        echo "  Creating combined heatmap..."
        plotHeatmap \
            -m "$OUTDIR/matrices/${peak_type}_combined_matrix.gz" \
            -o "$OUTDIR/heatmaps/${peak_type}_combined_heatmap.png" \
            --colorList "white,#2166AC" "white,#B2182B" "white,#1B7837" "white,#762A83" \
            --zMin 0 0 0 0 \
            --zMax 3 3 150 150 \
            --sortRegions descend \
            --sortUsing mean \
            --sortUsingSamples 1 \
            --heatmapHeight 15 \
            --heatmapWidth 3 \
            --xAxisLabel "" \
            --refPointLabel "Peak" \
            --plotTitle "Binding & Methylation at $title" \
            --legendLocation upper-right \
            --dpi 300 \
            2>/dev/null || echo "  WARNING: combined heatmap failed"

        # Profile with per-sample panels
        plotProfile \
            -m "$OUTDIR/matrices/${peak_type}_combined_matrix.gz" \
            -o "$OUTDIR/profiles/${peak_type}_combined_profile.png" \
            --plotTitle "Binding & Methylation at $title" \
            --perSample \
            --colors "#2166AC" "#B2182B" "#1B7837" "#762A83" \
            --legendLocation upper-right \
            --refPointLabel "Peak" \
            --plotHeight 8 \
            --plotWidth 14 \
            --dpi 300 \
            2>/dev/null || echo "  WARNING: combined profile failed"
    fi

    echo "  Done with $peak_type"
}

# ============================================================================
# Process each binding site category
# ============================================================================

echo ""
echo "=== Creating methylation profiles at binding sites ==="

# All TES binding sites
create_methylation_at_binding "TES_all" "$OUTDIR/beds/TES_all_peaks.bed" "All TES Binding Sites"

# All TEAD1 binding sites
create_methylation_at_binding "TEAD1_all" "$OUTDIR/beds/TEAD1_all_peaks.bed" "All TEAD1 Binding Sites"

# TES-unique binding sites
create_methylation_at_binding "TES_unique" "$OUTDIR/beds/TES_unique_peaks.bed" "TES-Unique Binding Sites"

# TEAD1-unique binding sites
create_methylation_at_binding "TEAD1_unique" "$OUTDIR/beds/TEAD1_unique_peaks.bed" "TEAD1-Unique Binding Sites"

# Shared binding sites
create_methylation_at_binding "shared" "$OUTDIR/beds/shared_peaks.bed" "Shared (TES+TEAD1) Binding Sites"

# TES peaks with DMR overlap
if [[ -f "$OUTDIR/beds/TES_peaks_with_DMR.bed" ]]; then
    create_methylation_at_binding "TES_with_DMR" "$OUTDIR/beds/TES_peaks_with_DMR.bed" "TES Sites WITH DMR Overlap"
fi

# TES peaks without DMR overlap
if [[ -f "$OUTDIR/beds/TES_peaks_no_DMR.bed" ]]; then
    create_methylation_at_binding "TES_no_DMR" "$OUTDIR/beds/TES_peaks_no_DMR.bed" "TES Sites WITHOUT DMR Overlap"
fi

# TES peaks with hypermethylated DMR overlap
if [[ -f "$OUTDIR/beds/TES_peaks_with_hyperDMR.bed" ]]; then
    create_methylation_at_binding "TES_with_hyperDMR" "$OUTDIR/beds/TES_peaks_with_hyperDMR.bed" "TES Sites WITH Hypermethylated DMR"
fi

# ============================================================================
# Create comparative visualizations
# ============================================================================

echo ""
echo "=== Creating comparative visualizations ==="

# Compare TES vs TEAD1 binding sites
if [[ -f "$OUTDIR/beds/TES_all_peaks.bed" ]] && [[ -f "$OUTDIR/beds/TEAD1_all_peaks.bed" ]]; then
    echo "  Computing TES vs TEAD1 comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$OUTDIR/beds/TES_all_peaks.bed" "$OUTDIR/beds/TEAD1_all_peaks.bed" \
        -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
        --samplesLabel "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/TES_vs_TEAD1_sites_medip_matrix.gz" \
        2>/dev/null || echo "  WARNING: comparison matrix failed"

    # Heatmap
    plotHeatmap \
        -m "$OUTDIR/matrices/TES_vs_TEAD1_sites_medip_matrix.gz" \
        -o "$OUTDIR/heatmaps/TES_vs_TEAD1_sites_medip_heatmap.png" \
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
        2>/dev/null || echo "  WARNING: comparison heatmap failed"

    # Profile
    plotProfile \
        -m "$OUTDIR/matrices/TES_vs_TEAD1_sites_medip_matrix.gz" \
        -o "$OUTDIR/profiles/TES_vs_TEAD1_sites_medip_profile.png" \
        --plotTitle "meDIP at TES vs TEAD1 Binding Sites" \
        --perGroup \
        --colors "#1B7837" "#762A83" \
        --regionsLabel "TES sites" "TEAD1 sites" \
        --legendLocation upper-right \
        --refPointLabel "Peak" \
        --plotHeight 6 \
        --plotWidth 10 \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: comparison profile failed"
fi

# Compare TES-unique vs TEAD1-unique vs Shared
if [[ -f "$OUTDIR/beds/TES_unique_peaks.bed" ]] && \
   [[ -f "$OUTDIR/beds/TEAD1_unique_peaks.bed" ]] && \
   [[ -f "$OUTDIR/beds/shared_peaks.bed" ]]; then

    echo "  Computing 3-way comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$OUTDIR/beds/TES_unique_peaks.bed" \
           "$OUTDIR/beds/shared_peaks.bed" \
           "$OUTDIR/beds/TEAD1_unique_peaks.bed" \
        -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
        --samplesLabel "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/binding_categories_medip_matrix.gz" \
        2>/dev/null || echo "  WARNING: 3-way matrix failed"

    # Profile comparing 3 binding categories
    plotProfile \
        -m "$OUTDIR/matrices/binding_categories_medip_matrix.gz" \
        -o "$OUTDIR/profiles/binding_categories_medip_profile.png" \
        --plotTitle "meDIP at Different Binding Site Categories" \
        --perGroup \
        --colors "#1B7837" "#762A83" \
        --regionsLabel "TES-unique" "Shared" "TEAD1-unique" \
        --legendLocation upper-right \
        --refPointLabel "Peak" \
        --plotHeight 6 \
        --plotWidth 12 \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: 3-way profile failed"

    # Heatmap
    plotHeatmap \
        -m "$OUTDIR/matrices/binding_categories_medip_matrix.gz" \
        -o "$OUTDIR/heatmaps/binding_categories_medip_heatmap.png" \
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
        2>/dev/null || echo "  WARNING: 3-way heatmap failed"
fi

# Compare TES sites WITH vs WITHOUT DMR overlap
if [[ -f "$OUTDIR/beds/TES_peaks_with_DMR.bed" ]] && \
   [[ -f "$OUTDIR/beds/TES_peaks_no_DMR.bed" ]]; then

    echo "  Computing TES +/- DMR comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$OUTDIR/beds/TES_peaks_with_DMR.bed" "$OUTDIR/beds/TES_peaks_no_DMR.bed" \
        -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
        --samplesLabel "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/TES_DMR_overlap_matrix.gz" \
        2>/dev/null || echo "  WARNING: DMR overlap matrix failed"

    # Profile
    plotProfile \
        -m "$OUTDIR/matrices/TES_DMR_overlap_matrix.gz" \
        -o "$OUTDIR/profiles/TES_DMR_overlap_profile.png" \
        --plotTitle "meDIP at TES Sites: With vs Without DMR Overlap" \
        --perGroup \
        --colors "#1B7837" "#762A83" \
        --regionsLabel "TES+DMR" "TES only" \
        --legendLocation upper-right \
        --refPointLabel "Peak" \
        --plotHeight 6 \
        --plotWidth 10 \
        --dpi 300 \
        2>/dev/null || echo "  WARNING: DMR overlap profile failed"
fi

# ============================================================================
# Calculate quantitative statistics
# ============================================================================

echo ""
echo "=== Calculating quantitative statistics ==="

# Switch to R environment for statistics
conda activate seurat_full2

cat > "$OUTDIR/statistics/calculate_binding_methylation_stats.R" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

cat("Calculating methylation statistics at binding sites...\n\n")

# ============================================================================
# Load data
# ============================================================================

# Load DMRs
dmr_file <- "results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv"
if (file.exists(dmr_file)) {
    dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)
    dmrs <- dmrs %>%
        mutate(
            dmr_id = paste(chr, start, end, sep = "_"),
            direction = ifelse(logFC > 0, "hyper", "hypo")
        )
    cat(sprintf("Loaded %d DMRs\n", nrow(dmrs)))
} else {
    stop("DMR file not found")
}

# Function to load BED file
load_bed <- function(filepath) {
    if (!file.exists(filepath)) return(NULL)
    df <- read.table(filepath, stringsAsFactors = FALSE)
    colnames(df)[1:3] <- c("chr", "start", "end")
    df$peak_id <- paste(df$chr, df$start, df$end, sep = "_")
    return(df)
}

# Load binding site categories
beds_dir <- file.path(outdir, "beds")
tes_all <- load_bed(file.path(beds_dir, "TES_all_peaks.bed"))
tead1_all <- load_bed(file.path(beds_dir, "TEAD1_all_peaks.bed"))
tes_unique <- load_bed(file.path(beds_dir, "TES_unique_peaks.bed"))
tead1_unique <- load_bed(file.path(beds_dir, "TEAD1_unique_peaks.bed"))
shared <- load_bed(file.path(beds_dir, "shared_peaks.bed"))
tes_with_dmr <- load_bed(file.path(beds_dir, "TES_peaks_with_DMR.bed"))
tes_no_dmr <- load_bed(file.path(beds_dir, "TES_peaks_no_DMR.bed"))

# ============================================================================
# Calculate statistics
# ============================================================================

stats_list <- list()

add_stat <- function(category, metric, value) {
    stats_list[[length(stats_list) + 1]] <<- data.frame(
        Category = category,
        Metric = metric,
        Value = as.character(value),
        stringsAsFactors = FALSE
    )
}

# Overview
add_stat("Overview", "Total_TES_peaks", nrow(tes_all))
add_stat("Overview", "Total_TEAD1_peaks", nrow(tead1_all))
add_stat("Overview", "Total_DMRs", nrow(dmrs))
add_stat("Overview", "Hypermethylated_DMRs", sum(dmrs$direction == "hyper"))
add_stat("Overview", "Hypomethylated_DMRs", sum(dmrs$direction == "hypo"))

# Binding site categories
add_stat("Binding_Categories", "TES_unique_peaks", nrow(tes_unique))
add_stat("Binding_Categories", "TEAD1_unique_peaks", nrow(tead1_unique))
add_stat("Binding_Categories", "Shared_peaks", nrow(shared))

# DMR overlap
if (!is.null(tes_with_dmr) && !is.null(tes_no_dmr)) {
    n_tes_dmr <- nrow(tes_with_dmr)
    n_tes_no_dmr <- nrow(tes_no_dmr)
    pct_tes_dmr <- round(100 * n_tes_dmr / nrow(tes_all), 1)

    add_stat("DMR_Overlap", "TES_peaks_with_DMR", n_tes_dmr)
    add_stat("DMR_Overlap", "TES_peaks_without_DMR", n_tes_no_dmr)
    add_stat("DMR_Overlap", "Pct_TES_peaks_with_DMR", paste0(pct_tes_dmr, "%"))

    # TEAD1 DMR overlap (if available)
    tead1_with_dmr <- load_bed(file.path(beds_dir, "TEAD1_peaks_with_DMR.bed"))
    if (!is.null(tead1_with_dmr)) {
        n_tead1_dmr <- nrow(tead1_with_dmr)
        pct_tead1_dmr <- round(100 * n_tead1_dmr / nrow(tead1_all), 1)
        add_stat("DMR_Overlap", "TEAD1_peaks_with_DMR", n_tead1_dmr)
        add_stat("DMR_Overlap", "Pct_TEAD1_peaks_with_DMR", paste0(pct_tead1_dmr, "%"))
    }
}

# Fisher's exact test: Is TES binding enriched at DMRs compared to TEAD1?
# This tests whether TES specifically recruits methylation machinery

# We need genome-wide background for proper enrichment
# Simplified: compare TES vs TEAD1 DMR overlap rates
if (!is.null(tes_with_dmr) && !is.null(tead1_all)) {
    tead1_with_dmr <- load_bed(file.path(beds_dir, "TEAD1_peaks_with_DMR.bed"))
    if (!is.null(tead1_with_dmr)) {
        # Contingency table
        tes_dmr <- nrow(tes_with_dmr)
        tes_no <- nrow(tes_all) - tes_dmr
        tead1_dmr <- nrow(tead1_with_dmr)
        tead1_no <- nrow(tead1_all) - tead1_dmr

        contingency <- matrix(c(tes_dmr, tes_no, tead1_dmr, tead1_no), nrow = 2,
                              dimnames = list(DMR = c("Yes", "No"), TF = c("TES", "TEAD1")))

        fisher_result <- fisher.test(contingency)

        add_stat("Enrichment_Test", "Fisher_OR_TES_vs_TEAD1_at_DMRs", round(fisher_result$estimate, 3))
        add_stat("Enrichment_Test", "Fisher_pvalue", format(fisher_result$p.value, scientific = TRUE, digits = 3))
        add_stat("Enrichment_Test", "Fisher_95CI_low", round(fisher_result$conf.int[1], 3))
        add_stat("Enrichment_Test", "Fisher_95CI_high", round(fisher_result$conf.int[2], 3))

        cat("\n=== Fisher's Exact Test: TES vs TEAD1 DMR Overlap ===\n")
        cat(sprintf("  TES peaks with DMR: %d / %d (%.1f%%)\n", tes_dmr, nrow(tes_all), 100*tes_dmr/nrow(tes_all)))
        cat(sprintf("  TEAD1 peaks with DMR: %d / %d (%.1f%%)\n", tead1_dmr, nrow(tead1_all), 100*tead1_dmr/nrow(tead1_all)))
        cat(sprintf("  Odds Ratio: %.3f (95%% CI: %.3f - %.3f)\n",
                    fisher_result$estimate, fisher_result$conf.int[1], fisher_result$conf.int[2]))
        cat(sprintf("  p-value: %s\n", format(fisher_result$p.value, scientific = TRUE, digits = 3)))

        if (fisher_result$estimate > 1 && fisher_result$p.value < 0.05) {
            cat("  Interpretation: TES binding is ENRICHED at DMRs compared to TEAD1\n")
        } else if (fisher_result$estimate < 1 && fisher_result$p.value < 0.05) {
            cat("  Interpretation: TES binding is DEPLETED at DMRs compared to TEAD1\n")
        } else {
            cat("  Interpretation: No significant difference in DMR overlap\n")
        }
    }
}

# ============================================================================
# Write output
# ============================================================================

all_stats <- do.call(rbind, stats_list)
output_file <- file.path(outdir, "statistics", "binding_methylation_statistics.csv")
write.csv(all_stats, output_file, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", output_file))

# Summary
cat("\n=== SUMMARY ===\n")
cat(sprintf("TES binding sites: %d\n", nrow(tes_all)))
cat(sprintf("  - TES-unique: %d\n", nrow(tes_unique)))
cat(sprintf("  - Shared with TEAD1: %d\n", nrow(shared)))
if (!is.null(tes_with_dmr)) {
    cat(sprintf("  - Overlapping DMRs: %d (%.1f%%)\n",
                nrow(tes_with_dmr), 100*nrow(tes_with_dmr)/nrow(tes_all)))
}
cat(sprintf("\nTEAD1 binding sites: %d\n", nrow(tead1_all)))
cat(sprintf("  - TEAD1-unique: %d\n", nrow(tead1_unique)))

RSCRIPT

# Run statistics script
Rscript "$OUTDIR/statistics/calculate_binding_methylation_stats.R" "$OUTDIR" 2>&1 | \
    grep -v "^WARNING\|overwriting" || echo "  Statistics calculation completed"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Analysis Complete: Does TES Binding Cause Methylation?"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Key outputs:"
echo ""
echo "BED files (regions):"
echo "  beds/TES_all_peaks.bed          - All TES binding sites"
echo "  beds/TEAD1_all_peaks.bed        - All TEAD1 binding sites"
echo "  beds/TES_unique_peaks.bed       - TES-specific binding sites"
echo "  beds/TEAD1_unique_peaks.bed     - TEAD1-specific binding sites"
echo "  beds/shared_peaks.bed           - Shared TES+TEAD1 sites"
echo "  beds/TES_peaks_with_DMR.bed     - TES sites overlapping DMRs"
echo "  beds/TES_peaks_no_DMR.bed       - TES sites without DMR overlap"
echo ""
echo "Heatmaps (meDIP signal at binding sites):"
ls -1 "$OUTDIR/heatmaps/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Profiles (meDIP signal at binding sites):"
ls -1 "$OUTDIR/profiles/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Statistics:"
echo "  statistics/binding_methylation_statistics.csv"
echo ""
echo "=============================================="
echo "INTERPRETATION GUIDE"
echo "=============================================="
echo ""
echo "If TES binding CAUSES methylation, we expect:"
echo "  1. Higher meDIP signal in TES samples at TES binding sites"
echo "  2. Enrichment of DMRs at TES binding sites (vs TEAD1)"
echo "  3. Stronger methylation at TES-unique sites than TEAD1-unique"
echo ""
echo "Key comparisons:"
echo "  - TES_all_medip_profile.png: Is TES meDIP > GFP meDIP at TES sites?"
echo "  - TES_vs_TEAD1_sites_medip_profile.png: Compare methylation at TES vs TEAD1 sites"
echo "  - binding_categories_medip_profile.png: TES-unique vs Shared vs TEAD1-unique"
echo "  - TES_DMR_overlap_profile.png: TES sites with vs without DMR"
echo ""
echo "Finished: $(date)"
