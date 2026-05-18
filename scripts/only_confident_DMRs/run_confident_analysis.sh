#!/bin/bash

# ============================================================================
# Master Script: High-Confidence DMR Analysis Pipeline
# ============================================================================
# Purpose: Run all high-confidence DMR analysis scripts in order
#
# This pipeline filters DMRs to retain only regions where BOTH TES and GFP
# samples have meaningful signal (>2 reads), addressing the GFP library
# quality issues that cause artifactual results.
#
# Usage:
#   ./run_confident_analysis.sh           # Submit all jobs with dependencies
#   ./run_confident_analysis.sh --dry-run # Show what would be submitted
#
# Scripts executed:
#   1. 00_filter_confident_dmrs.sh - Filter DMRs for high confidence
#   2. 17_dmr_stratified_heatmaps_confident.sh - Create heatmaps
#   3. 18_dmr_binding_overlay_confident.sh - Create binding overlays
#   4. 29_degs_down_binding_methylation_confident.sh - DEGs binding+methylation
#   5. 30_peaked_genes_metagene_confident.sh - Peaked genes metagene profiles
#   6. 31_enhancer_analysis_confident.sh - Enhancer analysis (DEGs DOWN)
#   7. 32_encode_enhancer_methylation_confident.sh - ENCODE enhancer methylation
#   8. 33_encode_enhancer_degs_down_confident.sh - ENCODE enhancer DEGs DOWN
# ============================================================================

set -e

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

SCRIPT_DIR="scripts/only_confident_DMRs"

echo "=============================================="
echo "High-Confidence DMR Analysis Pipeline"
echo "=============================================="
echo ""

# Check for dry-run mode
DRY_RUN=false
if [[ "$1" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "DRY RUN MODE - showing what would be submitted"
    echo ""
fi

# Create logs directory
mkdir -p "$SCRIPT_DIR/logs"

# Make scripts executable
chmod +x "$SCRIPT_DIR"/*.sh

echo "=== Pipeline Steps ==="
echo ""
echo "Step 1: Filter DMRs (00_filter_confident_dmrs.sh)"
echo "  - Filters for DMRs where BOTH samples have >2 reads"
echo "  - Expected: ~2,000 DMRs (from ~32,000 original)"
echo "  - Output: results/07_differential_MEDIPS_confident/"
echo ""
echo "Step 2: Create Heatmaps (17_dmr_stratified_heatmaps_confident.sh)"
echo "  - Creates heatmaps for high-confidence DMR categories"
echo "  - Output: results/17_dmr_heatmaps_confident/"
echo ""
echo "Step 3: Create Binding Overlays (18_dmr_binding_overlay_confident.sh)"
echo "  - Visualizes TES/TEAD1 binding at high-confidence DMRs"
echo "  - Output: results/18_binding_overlay_confident/"
echo ""
echo "Step 4: DEGs Binding+Methylation (29_degs_down_binding_methylation_confident.sh)"
echo "  - Integrative analysis: DEGs stratified by binding and DMR status"
echo "  - Output: scripts/only_confident_DMRs/output/29_degs_down_binding_methylation_confident/"
echo ""
echo "Step 5: Peaked Genes Metagene (30_peaked_genes_metagene_confident.sh)"
echo "  - Metagene profiles for genes with Cut&Tag peaks, stratified by DMR"
echo "  - Output: scripts/only_confident_DMRs/output/30_peaked_genes_metagene_confident/"
echo ""
echo "Step 6: Enhancer Analysis (31_enhancer_analysis_confident.sh)"
echo "  - Enhancers of DEGs DOWN stratified by DMR status"
echo "  - Output: scripts/only_confident_DMRs/output/31_enhancer_analysis_confident/"
echo ""
echo "Step 7: ENCODE Enhancer Methylation (32_encode_enhancer_methylation_confident.sh)"
echo "  - ENCODE enhancers stratified by TES binding and DMR status"
echo "  - Output: scripts/only_confident_DMRs/output/32_encode_enhancer_confident/"
echo ""
echo "Step 8: ENCODE Enhancer DEGs DOWN (33_encode_enhancer_degs_down_confident.sh)"
echo "  - ENCODE enhancers of DEGs DOWN stratified by binding and DMR"
echo "  - Output: scripts/only_confident_DMRs/output/33_encode_enhancer_degs_down_confident/"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    echo "=== DRY RUN: Commands that would be executed ==="
    echo ""
    echo "sbatch $SCRIPT_DIR/00_filter_confident_dmrs.sh"
    echo "sbatch --dependency=afterok:\$JOB1 $SCRIPT_DIR/17_dmr_stratified_heatmaps_confident.sh"
    echo "sbatch --dependency=afterok:\$JOB2 $SCRIPT_DIR/18_dmr_binding_overlay_confident.sh"
    echo "sbatch --dependency=afterok:\$JOB1 $SCRIPT_DIR/29_degs_down_binding_methylation_confident.sh"
    echo "sbatch --dependency=afterok:\$JOB1 $SCRIPT_DIR/30_peaked_genes_metagene_confident.sh"
    echo "sbatch --dependency=afterok:\$JOB1 $SCRIPT_DIR/31_enhancer_analysis_confident.sh"
    echo "sbatch --dependency=afterok:\$JOB1 $SCRIPT_DIR/32_encode_enhancer_methylation_confident.sh"
    echo "sbatch --dependency=afterok:\$JOB7 $SCRIPT_DIR/33_encode_enhancer_degs_down_confident.sh"
    echo ""
    echo "To actually run, execute without --dry-run flag"
    exit 0
fi

echo "=== Submitting Jobs ==="
echo ""

# Step 1: Filter DMRs
echo "Submitting Step 1: Filter confident DMRs..."
JOB1=$(sbatch --parsable "$SCRIPT_DIR/00_filter_confident_dmrs.sh")
echo "  Job ID: $JOB1"

# Step 2: Create heatmaps (depends on Step 1)
echo "Submitting Step 2: Create heatmaps (depends on $JOB1)..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 "$SCRIPT_DIR/17_dmr_stratified_heatmaps_confident.sh")
echo "  Job ID: $JOB2"

# Step 3: Create binding overlays (depends on Step 2)
echo "Submitting Step 3: Create binding overlays (depends on $JOB2)..."
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 "$SCRIPT_DIR/18_dmr_binding_overlay_confident.sh")
echo "  Job ID: $JOB3"

# Step 4: DEGs binding+methylation (depends on Step 1 only - can run in parallel with 2&3)
echo "Submitting Step 4: DEGs binding+methylation (depends on $JOB1)..."
JOB4=$(sbatch --parsable --dependency=afterok:$JOB1 "$SCRIPT_DIR/29_degs_down_binding_methylation_confident.sh")
echo "  Job ID: $JOB4"

# Step 5: Peaked genes metagene (depends on Step 1)
echo "Submitting Step 5: Peaked genes metagene (depends on $JOB1)..."
JOB5=$(sbatch --parsable --dependency=afterok:$JOB1 "$SCRIPT_DIR/30_peaked_genes_metagene_confident.sh")
echo "  Job ID: $JOB5"

# Step 6: Enhancer analysis (depends on Step 1)
echo "Submitting Step 6: Enhancer analysis (depends on $JOB1)..."
JOB6=$(sbatch --parsable --dependency=afterok:$JOB1 "$SCRIPT_DIR/31_enhancer_analysis_confident.sh")
echo "  Job ID: $JOB6"

# Step 7: ENCODE enhancer methylation (depends on Step 1)
echo "Submitting Step 7: ENCODE enhancer methylation (depends on $JOB1)..."
JOB7=$(sbatch --parsable --dependency=afterok:$JOB1 "$SCRIPT_DIR/32_encode_enhancer_methylation_confident.sh")
echo "  Job ID: $JOB7"

# Step 8: ENCODE enhancer DEGs DOWN (depends on Step 7 for ENCODE data)
echo "Submitting Step 8: ENCODE enhancer DEGs DOWN (depends on $JOB7)..."
JOB8=$(sbatch --parsable --dependency=afterok:$JOB7 "$SCRIPT_DIR/33_encode_enhancer_degs_down_confident.sh")
echo "  Job ID: $JOB8"

echo ""
echo "=== Pipeline Submitted ==="
echo ""
echo "Job IDs:"
echo "  Step 1 (Filter DMRs):           $JOB1"
echo "  Step 2 (Heatmaps):              $JOB2"
echo "  Step 3 (Binding Overlay):       $JOB3"
echo "  Step 4 (DEGs Bind+Meth):        $JOB4"
echo "  Step 5 (Peaked Genes):          $JOB5"
echo "  Step 6 (Enhancer Analysis):     $JOB6"
echo "  Step 7 (ENCODE Enhancer):       $JOB7"
echo "  Step 8 (ENCODE DEGs DOWN):      $JOB8"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo ""
echo "View logs:"
echo "  tail -f $SCRIPT_DIR/logs/00_filter_confident_dmrs.out"
echo "  tail -f $SCRIPT_DIR/logs/17_dmr_stratified_heatmaps_confident.out"
echo "  tail -f $SCRIPT_DIR/logs/18_dmr_binding_overlay_confident.out"
echo "  tail -f $SCRIPT_DIR/logs/29_degs_down_binding_methylation_confident.out"
echo "  tail -f $SCRIPT_DIR/logs/30_peaked_genes_metagene_confident.out"
echo "  tail -f $SCRIPT_DIR/logs/31_enhancer_analysis_confident.out"
echo "  tail -f $SCRIPT_DIR/logs/32_encode_enhancer_methylation_confident.out"
echo "  tail -f $SCRIPT_DIR/logs/33_encode_enhancer_degs_down_confident.out"
echo ""
echo "Expected results:"
echo "  results/07_differential_MEDIPS_confident/filtering_summary.txt"
echo "  results/17_dmr_heatmaps_confident/heatmaps/*.png"
echo "  results/18_binding_overlay_confident/heatmaps/*.png"
echo "  $SCRIPT_DIR/output/29_degs_down_binding_methylation_confident/*.png"
echo "  $SCRIPT_DIR/output/30_peaked_genes_metagene_confident/*.png"
echo "  $SCRIPT_DIR/output/31_enhancer_analysis_confident/*.png"
echo "  $SCRIPT_DIR/output/32_encode_enhancer_confident/*.png"
echo "  $SCRIPT_DIR/output/33_encode_enhancer_degs_down_confident/*.png"
echo ""
echo "Key finding to verify:"
echo "  Original: ~91% hypermethylated / ~9% hypomethylated"
echo "  Expected: ~35% hypermethylated / ~65% hypomethylated (REVERSED!)"
echo ""
