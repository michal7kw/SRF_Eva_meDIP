#!/bin/bash
#SBATCH --job-name=dmr_feature_dist
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/08b_feature_distribution.out
#SBATCH --error=logs/08b_feature_distribution.err

################################################################################
# Script: 08b_genomic_feature_distribution.sh
# Purpose: Run enhanced genomic feature distribution visualization for DMRs
#
# Description:
#   Creates publication-quality multi-panel plots showing DMR distribution
#   across genomic features, with separate panels for hypermethylated and
#   hypomethylated DMRs.
#
# Prerequisites:
#   - Step 08 (08_annotation.R) must have completed successfully
#   - Annotated DMR files must exist in results/08_annotation/
#
# Output:
#   - Multi-panel PDFs and PNGs in results/08_annotation/feature_distribution/
#   - Summary statistics CSV
#
# Runtime: ~10-15 minutes
################################################################################

set -euo pipefail

echo "=============================================="
echo "DMR Genomic Feature Distribution Analysis"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node: ${SLURM_JOB_NODELIST:-local}"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate annotation_enrichment

# Set working directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
cd $BASE_DIR

# Verify input files exist
if [ ! -d "results/08_annotation" ] || [ -z "$(ls -A results/08_annotation/*_annotated.csv 2>/dev/null)" ]; then
    echo "ERROR: No annotated DMR files found in results/08_annotation/"
    echo "Please run 08_annotation.R first (step 08)"
    exit 1
fi

echo "Input files found:"
ls -la results/08_annotation/*_annotated.csv
echo ""

# Create output directory
mkdir -p results/08_annotation/feature_distribution

# Run R script
echo "Running visualization script..."
echo ""
Rscript "./scripts/08b_genomic_feature_distribution.R"

echo ""
echo "=============================================="
echo "Output files:"
echo "=============================================="
ls -la results/08_annotation/feature_distribution/

echo ""
echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="
