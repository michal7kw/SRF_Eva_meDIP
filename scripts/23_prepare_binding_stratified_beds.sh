#!/bin/bash
#SBATCH --job-name=23_prep_beds
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/23_prepare_binding_stratified_beds_%j.out
#SBATCH --error=logs/23_prepare_binding_stratified_beds_%j.err

# ============================================================================
# 23_prepare_binding_stratified_beds.sh
# ============================================================================
# Purpose: SLURM wrapper for R script that creates gene sets for
#          binding-stratified MeDIP metaprofile analysis
#
# Usage: sbatch 23_prepare_binding_stratified_beds.sh
# ============================================================================

echo "============================================"
echo "Binding-Stratified BED File Preparation"
echo "============================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "============================================"

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create logs directory if needed
mkdir -p logs

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

echo ""
echo "Using R version:"
R --version | head -1

echo ""
echo "Running R script..."
echo ""

# Run the R script
Rscript scripts/23_prepare_binding_stratified_beds.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "============================================"
    echo "BED file preparation completed successfully!"
    echo "============================================"
    echo ""
    echo "Output files:"
    ls -lh results/23_binding_stratified_metaprofiles/beds/
    echo ""
    echo "Next step: Run 23_compute_binding_matrices.sh"
else
    echo ""
    echo "ERROR: R script failed!"
    exit 1
fi

echo ""
echo "Finished: $(date)"
