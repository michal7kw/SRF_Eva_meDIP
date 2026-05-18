#!/bin/bash

################################################################################
# Master Pipeline: meDIP Heatmap and Metaprofile Analysis
################################################################################
#
# Purpose:
#   Orchestrates the complete meDIP visualization pipeline from gene set
#   definition through advanced integrated visualizations.
#
# Pipeline Steps:
#   1. Define gene sets based on RNA-seq expression (12_define_gene_sets.sh)
#   2. Compute meDIP signal matrices at promoters (13_compute_matrices.sh)
#   3. Generate heatmaps with deepTools (14_plot_heatmaps.sh)
#   4. Generate metaprofiles with deepTools (15_plot_metaprofiles.sh)
#   5. Advanced integrated visualizations (16_advanced_visualization.sh)
#
# Usage:
#   ./master_heatmap_pipeline.sh [--dry-run] [--from-step N]
#
# Options:
#   --dry-run       Print commands without submitting jobs
#   --from-step N   Start from step N (skip earlier steps)
#
# Example:
#   # Run entire pipeline
#   ./master_heatmap_pipeline.sh
#
#   # Test without submitting
#   ./master_heatmap_pipeline.sh --dry-run
#
#   # Start from matrix computation (skip gene set definition)
#   ./master_heatmap_pipeline.sh --from-step 2
#
################################################################################

set -euo pipefail

# Define directories
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts"
LOG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/logs"

mkdir -p "$LOG_DIR"

# Parse command line arguments
DRY_RUN=false
START_STEP=1

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --from-step)
            START_STEP="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run] [--from-step N]"
            exit 1
            ;;
    esac
done

# Function to submit job with dependency
submit_job() {
    local step=$1
    local script=$2
    local dependency=$3
    local job_name=$4

    if [ $step -lt $START_STEP ]; then
        echo "  Skipping step $step (starting from step $START_STEP)"
        return
    fi

    echo ""
    echo "========================================"
    echo "Step $step: $job_name"
    echo "========================================"
    echo "Script: $script"

    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] Would submit: sbatch $script"
        if [ -n "$dependency" ]; then
            echo "[DRY RUN] With dependency: $dependency"
        fi
        return
    fi

    if [ -n "$dependency" ]; then
        JOB_ID=$(sbatch --dependency=afterok:$dependency --parsable "$script")
    else
        JOB_ID=$(sbatch --parsable "$script")
    fi

    echo "  Submitted job ID: $JOB_ID"
    echo "$JOB_ID"
}

################################################################################
# Pipeline Execution
################################################################################

echo "========================================"
echo "meDIP Heatmap & Metaprofile Pipeline"
echo "========================================"
echo "Start time: $(date)"
echo ""

if [ "$DRY_RUN" = true ]; then
    echo "*** DRY RUN MODE - No jobs will be submitted ***"
    echo ""
fi

if [ $START_STEP -gt 1 ]; then
    echo "*** Starting from step $START_STEP ***"
    echo ""
fi

# Track job IDs for dependencies
PREV_JOB_ID=""

################################################################################
# Step 1: Define Gene Sets
################################################################################

STEP1_SCRIPT="${SCRIPT_DIR}/12_define_gene_sets.sh"

if [ -f "$STEP1_SCRIPT" ]; then
    STEP1_JOB=$(submit_job 1 "$STEP1_SCRIPT" "" "Define Gene Sets")
    PREV_JOB_ID=$STEP1_JOB
else
    echo "ERROR: Step 1 script not found: $STEP1_SCRIPT"
    exit 1
fi

################################################################################
# Step 2: Compute Signal Matrices
################################################################################

STEP2_SCRIPT="${SCRIPT_DIR}/13_compute_matrices.sh"

if [ -f "$STEP2_SCRIPT" ]; then
    STEP2_JOB=$(submit_job 2 "$STEP2_SCRIPT" "$PREV_JOB_ID" "Compute Signal Matrices")
    PREV_JOB_ID=$STEP2_JOB
else
    echo "ERROR: Step 2 script not found: $STEP2_SCRIPT"
    exit 1
fi

################################################################################
# Step 3: Generate Heatmaps (depends on Step 2)
################################################################################

STEP3_SCRIPT="${SCRIPT_DIR}/14_plot_heatmaps.sh"

if [ -f "$STEP3_SCRIPT" ]; then
    STEP3_JOB=$(submit_job 3 "$STEP3_SCRIPT" "$PREV_JOB_ID" "Generate Heatmaps")
else
    echo "WARNING: Step 3 script not found: $STEP3_SCRIPT"
    STEP3_JOB=""
fi

################################################################################
# Step 4: Generate Metaprofiles (depends on Step 2, parallel with Step 3)
################################################################################

STEP4_SCRIPT="${SCRIPT_DIR}/15_plot_metaprofiles.sh"

if [ -f "$STEP4_SCRIPT" ]; then
    # Use STEP2_JOB as dependency (can run in parallel with heatmaps)
    STEP4_JOB=$(submit_job 4 "$STEP4_SCRIPT" "$STEP2_JOB" "Generate Metaprofiles")
else
    echo "WARNING: Step 4 script not found: $STEP4_SCRIPT"
    STEP4_JOB=""
fi

################################################################################
# Step 5: Advanced Visualization (depends on Step 2)
################################################################################

STEP5_SCRIPT="${SCRIPT_DIR}/16_advanced_visualization.sh"

if [ -f "$STEP5_SCRIPT" ]; then
    # Use STEP2_JOB as dependency (can run in parallel with heatmaps/metaprofiles)
    STEP5_JOB=$(submit_job 5 "$STEP5_SCRIPT" "$STEP2_JOB" "Advanced Visualization")
else
    echo "WARNING: Step 5 script not found: $STEP5_SCRIPT"
    STEP5_JOB=""
fi

################################################################################
# Summary
################################################################################

echo ""
echo "========================================"
echo "Pipeline Submission Complete"
echo "========================================"
echo "End time: $(date)"
echo ""

if [ "$DRY_RUN" = false ]; then
    echo "Submitted jobs:"
    echo "  Step 1 (Gene Sets):       ${STEP1_JOB:-N/A}"
    echo "  Step 2 (Matrices):        ${STEP2_JOB:-N/A}"
    echo "  Step 3 (Heatmaps):        ${STEP3_JOB:-N/A}"
    echo "  Step 4 (Metaprofiles):    ${STEP4_JOB:-N/A}"
    echo "  Step 5 (Advanced Viz):    ${STEP5_JOB:-N/A}"
    echo ""
    echo "Monitor job status:"
    echo "  squeue -u $USER"
    echo ""
    echo "View logs:"
    echo "  tail -f ${LOG_DIR}/*.out"
    echo ""
    echo "Cancel all jobs:"
    if [ -n "${STEP1_JOB:-}" ]; then
        echo "  scancel ${STEP1_JOB} ${STEP2_JOB:-} ${STEP3_JOB:-} ${STEP4_JOB:-} ${STEP5_JOB:-}"
    fi
else
    echo "*** DRY RUN COMPLETE - No jobs were submitted ***"
    echo ""
    echo "To run the pipeline, execute:"
    echo "  $0"
fi

echo ""
echo "========================================"

################################################################################
# Pipeline Completion Message
################################################################################

if [ "$DRY_RUN" = false ]; then
    cat << 'EOF'

Expected Runtime: 4-6 hours total
  - Step 1: ~30 minutes
  - Step 2: ~2-3 hours (most time-consuming)
  - Step 3: ~30-60 minutes
  - Step 4: ~20-30 minutes
  - Step 5: ~10-20 minutes

Expected Outputs:
  results/12_gene_sets/         - Gene lists and promoter BED files
  results/13_matrices/          - meDIP signal matrices (deepTools)
  results/14_heatmaps/          - Heatmap visualizations (PDF)
  results/15_metaprofiles/      - Metaprofile plots (PDF)
  results/16_advanced_visualization/ - Integrated analyses (PDF + CSV)

Key Visualizations to Check:
  1. 14_heatmaps/*_GFP_vs_TES_sidebyside.pdf - Side-by-side comparison
  2. 14_heatmaps/*_clustered.pdf - K-means clustering patterns
  3. 15_metaprofiles/all_gene_sets_*_metaprofile.pdf - Average signals
  4. 16_advanced_visualization/integrated_medip_rnaseq_heatmap.pdf - ComplexHeatmap

Biological Questions Addressed:
  ✓ Do highly expressed genes have different promoter methylation?
  ✓ Are upregulated genes hypomethylated? Downregulated hypermethylated?
  ✓ Are there distinct clusters of coordinated methylation + expression?
  ✓ Does TES alter methylation-expression relationships?
  (Note: TESmut excluded from analysis - failed sample)

EOF
fi
