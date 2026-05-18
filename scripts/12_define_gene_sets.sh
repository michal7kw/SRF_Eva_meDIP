#!/bin/bash
#SBATCH --job-name=define_gene_sets
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/12_define_gene_sets.err"
#SBATCH --output="logs/12_define_gene_sets.out"
#SBATCH --partition=workq

################################################################################
# SLURM Wrapper: Define Gene Sets for meDIP Heatmap Analysis
################################################################################

echo "========================================="
echo "Define Gene Sets - SLURM Job"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Start time: $(date)"
echo ""

# Activate conda environment with R packages
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run R script
Rscript /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/12_define_gene_sets.R

echo ""
echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
