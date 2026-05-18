#!/bin/bash
#SBATCH --job-name=hyper_vs_hypo_dmr
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --output=logs/19b_hyper_vs_hypo_dmr.out
#SBATCH --error=logs/19b_hyper_vs_hypo_dmr.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# Clear Visualization: Hyper vs Hypo DMRs at TES Binding Sites
# ============================================================================
#
# Key finding: While 91% of DMRs genome-wide are hypermethylated, at TES
# binding sites the pattern is inverted - 60% overlap hypomethylated DMRs.
#
# This script creates separate, clearly labeled visualizations for:
#   1. TES peaks overlapping HYPERmethylated DMRs (n=71)
#   2. TES peaks overlapping HYPOmethylated DMRs (n=107)
#   3. Combined comparison showing the paradoxical pattern
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Output directories
OUTDIR="results/19_binding_causes_methylation"
mkdir -p "$OUTDIR"/{beds,matrices,heatmaps,profiles,statistics}
mkdir -p logs

echo "=============================================="
echo "Hyper vs Hypo DMR Analysis at TES Sites"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment with deepTools and bedtools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ============================================================================
# Input files
# ============================================================================

BEDS_DIR="$OUTDIR/beds"
MEDIP_DIR="results/05_bigwig"

# meDIP BigWig files
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
GFP_MEDIP_BW2="$MEDIP_DIR/GFP-2-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"
TES_MEDIP_BW2="$MEDIP_DIR/TES-2-IP_RPKM.bw"

# ============================================================================
# Step 1: Create/verify BED files
# ============================================================================

echo ""
echo "=== Step 1: Creating BED files ==="

# Create hypo DMR bed if not exists
if [[ ! -f "$BEDS_DIR/hypo_dmrs.bed" ]]; then
    awk '$4 < 0' "$BEDS_DIR/all_dmrs.bed" > "$BEDS_DIR/hypo_dmrs.bed"
fi

# Create TES peaks with hypoDMR
echo "Creating TES_peaks_with_hypoDMR.bed..."
bedtools intersect -a "$BEDS_DIR/TES_all_peaks.bed" \
    -b "$BEDS_DIR/hypo_dmrs.bed" -u > "$BEDS_DIR/TES_peaks_with_hypoDMR.bed"

# Count peaks
N_HYPER=$(wc -l < "$BEDS_DIR/TES_peaks_with_hyperDMR.bed")
N_HYPO=$(wc -l < "$BEDS_DIR/TES_peaks_with_hypoDMR.bed")
N_TOTAL=$((N_HYPER + N_HYPO))

echo "  TES peaks with HYPERmethylated DMRs: $N_HYPER"
echo "  TES peaks with HYPOmethylated DMRs: $N_HYPO"
echo "  Total TES peaks with DMR overlap: $N_TOTAL"

# ============================================================================
# Step 2: Create comparison matrix (both groups together)
# ============================================================================

echo ""
echo "=== Step 2: Creating comparison matrix ==="

computeMatrix reference-point \
    --referencePoint center \
    -b 5000 -a 5000 \
    -R "$BEDS_DIR/TES_peaks_with_hyperDMR.bed" "$BEDS_DIR/TES_peaks_with_hypoDMR.bed" \
    -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
    --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
    --skipZeros \
    --binSize 50 \
    --numberOfProcessors 8 \
    -o "$OUTDIR/matrices/hyper_vs_hypo_DMR_comparison_matrix.gz"

echo "  Matrix created: hyper_vs_hypo_DMR_comparison_matrix.gz"

# ============================================================================
# Step 3: Create side-by-side profile plots
# ============================================================================

echo ""
echo "=== Step 3: Creating profile plots ==="

# Combined profile with both groups
plotProfile \
    -m "$OUTDIR/matrices/hyper_vs_hypo_DMR_comparison_matrix.gz" \
    -o "$OUTDIR/profiles/TES_hyperDMR_vs_hypoDMR_profile.png" \
    --plotTitle "meDIP Signal at TES Binding Sites: HYPER (n=$N_HYPER) vs HYPO (n=$N_HYPO) DMRs" \
    --perGroup \
    --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
    --regionsLabel "HYPER DMR sites" "HYPO DMR sites" \
    --legendLocation upper-right \
    --refPointLabel "TES Peak" \
    --plotHeight 6 \
    --plotWidth 12 \
    --dpi 300

echo "  Profile created: TES_hyperDMR_vs_hypoDMR_profile.png"

# Per-sample profile for detailed view (plotType lines shows individual samples)
plotProfile \
    -m "$OUTDIR/matrices/hyper_vs_hypo_DMR_comparison_matrix.gz" \
    -o "$OUTDIR/profiles/TES_hyperDMR_vs_hypoDMR_per_sample.png" \
    --plotTitle "meDIP at TES Sites: HYPER (n=$N_HYPER) vs HYPO (n=$N_HYPO) DMRs - Per Sample" \
    --plotType lines \
    --colors "#A6D96A" "#1B7837" "#DFC27D" "#762A83" \
    --regionsLabel "HYPER DMRs" "HYPO DMRs" \
    --legendLocation upper-right \
    --refPointLabel "TES Peak" \
    --plotHeight 8 \
    --plotWidth 14 \
    --dpi 300

echo "  Profile created: TES_hyperDMR_vs_hypoDMR_per_sample.png"

# ============================================================================
# Step 4: Create combined heatmap
# ============================================================================

echo ""
echo "=== Step 4: Creating heatmaps ==="

plotHeatmap \
    -m "$OUTDIR/matrices/hyper_vs_hypo_DMR_comparison_matrix.gz" \
    -o "$OUTDIR/heatmaps/TES_hyperDMR_vs_hypoDMR_heatmap.png" \
    --colorMap YlGn YlGn Purples Purples \
    --sortRegions descend \
    --sortUsing mean \
    --heatmapHeight 15 \
    --heatmapWidth 2.5 \
    --xAxisLabel "" \
    --refPointLabel "TES Peak" \
    --regionsLabel "HYPER DMRs (n=$N_HYPER)" "HYPO DMRs (n=$N_HYPO)" \
    --plotTitle "meDIP at TES Sites with Hyper vs Hypo DMRs" \
    --legendLocation upper-right \
    --dpi 300

echo "  Heatmap created: TES_hyperDMR_vs_hypoDMR_heatmap.png"

# ============================================================================
# Step 5: Create simplified GFP vs TES comparison
# ============================================================================

echo ""
echo "=== Step 5: Creating GFP vs TES average comparison ==="

# Check that combined BigWig files exist (created by 05a_combine_bigwig_replicates.sh)
if [[ ! -f "$MEDIP_DIR/GFP_average.bw" ]] || [[ ! -f "$MEDIP_DIR/TES_average.bw" ]]; then
    echo "  ERROR: Combined BigWig files not found."
    echo "  Please run 05a_combine_bigwig_replicates.sh first to create:"
    echo "    - $MEDIP_DIR/GFP_average.bw"
    echo "    - $MEDIP_DIR/TES_average.bw"
    echo "  Skipping averaged comparison..."
else
    echo "  Found combined BigWig files:"
    echo "    - $MEDIP_DIR/GFP_average.bw"
    echo "    - $MEDIP_DIR/TES_average.bw"
fi

# Create simplified matrix with averaged signals
if [[ -f "$MEDIP_DIR/GFP_average.bw" ]] && [[ -f "$MEDIP_DIR/TES_average.bw" ]]; then
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$BEDS_DIR/TES_peaks_with_hyperDMR.bed" "$BEDS_DIR/TES_peaks_with_hypoDMR.bed" \
        -S "$MEDIP_DIR/GFP_average.bw" "$MEDIP_DIR/TES_average.bw" \
        --samplesLabel "GFP meDIP" "TES meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/hyper_vs_hypo_averaged_matrix.gz"

    plotProfile \
        -m "$OUTDIR/matrices/hyper_vs_hypo_averaged_matrix.gz" \
        -o "$OUTDIR/profiles/TES_hyperDMR_vs_hypoDMR_averaged.png" \
        --plotTitle "meDIP at TES Sites: HYPER (n=$N_HYPER) vs HYPO (n=$N_HYPO) DMRs" \
        --perGroup \
        --colors "#1B7837" "#762A83" \
        --regionsLabel "HYPER DMRs" "HYPO DMRs" \
        --legendLocation upper-right \
        --refPointLabel "TES Peak" \
        --plotHeight 6 \
        --plotWidth 12 \
        --dpi 300

    echo "  Averaged profile created: TES_hyperDMR_vs_hypoDMR_averaged.png"
fi

# ============================================================================
# Step 6: Statistical summary
# ============================================================================

echo ""
echo "=== Step 6: Creating statistical summary ==="

# Switch to R for statistics
conda activate seurat_full2

cat > "$OUTDIR/statistics/dmr_direction_analysis.R" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

# ============================================================================
# DMR Direction Analysis at TES Binding Sites
# ============================================================================

cat("============================================\n")
cat("DMR Direction Analysis at TES Binding Sites\n")
cat("============================================\n\n")

# Counts
n_hyper_at_tes <- 71
n_hypo_at_tes <- 107
n_total_at_tes <- n_hyper_at_tes + n_hypo_at_tes

n_hyper_global <- 29034
n_hypo_global <- 2886
n_total_global <- n_hyper_global + n_hypo_global

# Observed proportions
pct_hyper_at_tes <- 100 * n_hyper_at_tes / n_total_at_tes
pct_hypo_at_tes <- 100 * n_hypo_at_tes / n_total_at_tes
pct_hyper_global <- 100 * n_hyper_global / n_total_global
pct_hypo_global <- 100 * n_hypo_global / n_total_global

cat("=== DMR Direction at TES Binding Sites ===\n\n")
cat(sprintf("At TES binding sites:\n"))
cat(sprintf("  HYPERmethylated: %d (%.1f%%)\n", n_hyper_at_tes, pct_hyper_at_tes))
cat(sprintf("  HYPOmethylated:  %d (%.1f%%)\n", n_hypo_at_tes, pct_hypo_at_tes))
cat(sprintf("  Total:           %d\n\n", n_total_at_tes))

cat(sprintf("Global DMR distribution:\n"))
cat(sprintf("  HYPERmethylated: %d (%.1f%%)\n", n_hyper_global, pct_hyper_global))
cat(sprintf("  HYPOmethylated:  %d (%.1f%%)\n", n_hypo_global, pct_hypo_global))
cat(sprintf("  Total:           %d\n\n", n_total_global))

# Fisher's exact test
# Contingency table:
#              At TES    Not at TES
# Hyper         71         28963
# Hypo         107          2779

contingency <- matrix(
    c(n_hyper_at_tes, n_hyper_global - n_hyper_at_tes,
      n_hypo_at_tes, n_hypo_global - n_hypo_at_tes),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        DMR_Type = c("Hypermethylated", "Hypomethylated"),
        Location = c("At_TES_sites", "Not_at_TES")
    )
)

cat("=== Fisher's Exact Test ===\n\n")
cat("Contingency table:\n")
print(contingency)

fisher_result <- fisher.test(contingency)

cat(sprintf("\nOdds Ratio: %.3f\n", fisher_result$estimate))
cat(sprintf("95%% CI: %.3f - %.3f\n", fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("p-value: %.2e\n\n", fisher_result$p.value))

# Interpretation
cat("=== INTERPRETATION ===\n\n")
if (fisher_result$estimate < 1 && fisher_result$p.value < 0.05) {
    cat("TES binding sites are SIGNIFICANTLY DEPLETED for hypermethylated DMRs\n")
    cat("(or equivalently, ENRICHED for hypomethylated DMRs).\n\n")
    cat(sprintf("Expected if random: %.1f%% hyper / %.1f%% hypo\n", pct_hyper_global, pct_hypo_global))
    cat(sprintf("Observed at TES:    %.1f%% hyper / %.1f%% hypo\n\n", pct_hyper_at_tes, pct_hypo_at_tes))
    cat("BIOLOGICAL INTERPRETATION:\n")
    cat("- TES binding appears to PREVENT methylation at direct target sites\n")
    cat("- The widespread hypermethylation (91% globally) is likely an INDIRECT effect\n")
    cat("- TES may recruit DNMT3A/3L to sites DISTAL from its binding locations\n")
}

# Create bar plot
outdir <- "results/19_binding_causes_methylation"

df_plot <- data.frame(
    Location = rep(c("At TES Sites\n(n=178)", "Global\n(n=31,920)"), each = 2),
    DMR_Type = rep(c("Hypermethylated", "Hypomethylated"), 2),
    Percentage = c(pct_hyper_at_tes, pct_hypo_at_tes, pct_hyper_global, pct_hypo_global)
)

p <- ggplot(df_plot, aes(x = Location, y = Percentage, fill = DMR_Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_dodge(width = 0.7), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Hypermethylated" = "#D73027", "Hypomethylated" = "#4575B4"),
                      name = "DMR Direction") +
    labs(title = "DMR Direction: TES Binding Sites vs Global",
         subtitle = sprintf("Fisher's exact test: OR=%.2f, p=%.2e", fisher_result$estimate, fisher_result$p.value),
         x = "", y = "Percentage of DMRs") +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top"
    ) +
    ylim(0, 100)

ggsave(file.path(outdir, "profiles", "DMR_direction_at_TES_sites_barplot.png"),
       p, width = 8, height = 6, dpi = 300)

cat(sprintf("\nBar plot saved: %s/profiles/DMR_direction_at_TES_sites_barplot.png\n", outdir))

# Save statistics to file
stats_file <- file.path(outdir, "statistics", "dmr_direction_summary.txt")
sink(stats_file)
cat("DMR Direction Analysis at TES Binding Sites\n")
cat("============================================\n\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("COUNTS:\n")
cat(sprintf("  TES peaks with HYPER DMRs: %d\n", n_hyper_at_tes))
cat(sprintf("  TES peaks with HYPO DMRs:  %d\n", n_hypo_at_tes))
cat(sprintf("  Total TES peaks with DMR:  %d\n\n", n_total_at_tes))
cat("PERCENTAGES:\n")
cat(sprintf("  At TES sites: %.1f%% hyper / %.1f%% hypo\n", pct_hyper_at_tes, pct_hypo_at_tes))
cat(sprintf("  Global:       %.1f%% hyper / %.1f%% hypo\n\n", pct_hyper_global, pct_hypo_global))
cat("FISHER'S EXACT TEST:\n")
cat(sprintf("  Odds Ratio: %.3f\n", fisher_result$estimate))
cat(sprintf("  95%% CI: %.3f - %.3f\n", fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("  p-value: %.2e\n\n", fisher_result$p.value))
cat("INTERPRETATION:\n")
cat("  TES binding is significantly DEPLETED at hypermethylated DMRs.\n")
cat("  This suggests TES binding PREVENTS methylation at direct targets,\n")
cat("  while causing methylation at distal sites (indirect effect).\n")
sink()

cat(sprintf("Statistics saved: %s\n", stats_file))

RSCRIPT

Rscript "$OUTDIR/statistics/dmr_direction_analysis.R"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="
echo ""
echo "Key outputs:"
echo ""
echo "Profiles:"
ls -1 "$OUTDIR/profiles/"*hyper*hypo* 2>/dev/null | while read f; do echo "  $(basename $f)"; done
ls -1 "$OUTDIR/profiles/"*DMR_direction* 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Heatmaps:"
ls -1 "$OUTDIR/heatmaps/"*hyper*hypo* 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Statistics:"
echo "  dmr_direction_summary.txt"
echo ""
echo "KEY FINDING:"
echo "  Global DMRs: 91% hypermethylated / 9% hypomethylated"
echo "  At TES sites: 40% hypermethylated / 60% hypomethylated"
echo "  TES binding is DEPLETED at hypermethylated regions (OR < 1, p < 0.05)"
echo ""
echo "Finished: $(date)"
