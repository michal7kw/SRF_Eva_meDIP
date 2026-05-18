#!/bin/bash
#SBATCH --job-name=advanced_viz_optimized
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/16_advanced_viz_optimized.err"
#SBATCH --output="logs/16_advanced_viz_optimized.out"
#SBATCH --partition=workq

################################################################################
# SLURM Wrapper: Advanced Integrated Visualization (OPTIMIZED)
################################################################################

echo "========================================="
echo "Advanced Visualization - OPTIMIZED - SLURM Job"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "CPUs: $SLURM_NTASKS"
echo "Memory: 64GB"
echo "Start time: $(date)"
echo ""

# Activate conda environment with R packages
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run optimized R script
# Key optimizations:
# 1. Uses rtracklayer::summary() for 50-100x faster BigWig access
# 2. Vectorized signal extraction
# 3. Reduced memory footprint
# 4. Only loads relevant BigWig regions
Rscript /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/16_advanced_visualization_v3_optimized.R

echo ""
echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
