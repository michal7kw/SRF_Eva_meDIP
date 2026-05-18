#!/bin/bash
#SBATCH --job-name=volcano_plots
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/volcano_plots_%j.out
#SBATCH --error=logs/volcano_plots_%j.err

################################################################################
# Create Volcano Plots for meDIP Differential Methylation
################################################################################

echo "=========================================="
echo "meDIP Volcano Plot Generation"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Setup
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP
mkdir -p logs

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

echo "Running R script..."
Rscript scripts/create_volcano_plots.R

echo ""
echo "=========================================="
echo "Job completed: $(date)"
echo "=========================================="
