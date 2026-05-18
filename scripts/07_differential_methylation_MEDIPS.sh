#!/bin/bash
#SBATCH --job-name=medip_MEDIPS
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/07_differential_methylation_MEDIPS.out
#SBATCH --error=logs/07_differential_methylation_MEDIPS.err

################################################################################
# meDIP-seq Differential Methylation Analysis using MEDIPS
#
# This script runs biologically rigorous quantitative methylation analysis
# using MEDIPS instead of peak calling.
#
# Why MEDIPS?
# - Treats methylation as quantitative (not binary peaks)
# - Corrects for CpG density bias (critical for meDIP!)
# - Genome-wide coverage in uniform windows
# - No arbitrary thresholds
# - Proper statistical framework
#
# Runtime: ~2-4 hours (genome-wide analysis)
# Memory: 64GB recommended
################################################################################

echo "========================================"
echo "meDIP-seq: MEDIPS Analysis"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo ""

# Activate conda environment with R and MEDIPS
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-medips-analysis

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "R version: $(R --version | head -n1)"
echo ""

# Set working directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts"
cd $BASE_DIR

echo "Working directory: $(pwd)"
echo ""

# Run MEDIPS analysis
echo "Starting MEDIPS differential methylation analysis..."
echo "This will take 2-4 hours for genome-wide analysis"
echo ""

Rscript "07_differential_methylation_MEDIPS.R"

EXIT_CODE=$?

echo ""
echo "========================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "MEDIPS analysis completed successfully!"
else
    echo "ERROR: MEDIPS analysis failed with exit code $EXIT_CODE"
fi
echo "End time: $(date)"
echo "========================================"

exit $EXIT_CODE
