#!/bin/bash
#SBATCH --job-name=medip_MEDIPS_INPUT
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/07_differential_methylation_MEDIPS_with_INPUT.out
#SBATCH --error=logs/07_differential_methylation_MEDIPS_with_INPUT.err

################################################################################
# meDIP-seq Differential Methylation Analysis using MEDIPS WITH INPUT
#
# This script runs MEDIPS analysis WITH PROPER INPUT NORMALIZATION.
#
# Key Difference from Original:
#   - INPUT samples are used for normalization (ISet1, ISet2 parameters)
#   - Corrects for copy number variation and chromatin accessibility bias
#   - Results saved to a SEPARATE directory to preserve original results
#
# INPUT Normalization Benefits:
#   - Removes genomic copy number effects
#   - Corrects for chromatin accessibility differences
#   - Accounts for sequencing depth variations between IP and INPUT
#   - More rigorous for meDIP-seq data interpretation
#
# NOTE: Using TESmut-1-INPUT as COMMON INPUT for all samples because:
#   - TES-1-INPUT does not exist in raw data
#   - GFP-1-INPUT was not processed through the pipeline
#   - TESmut-1-INPUT is from the same cell line (SNB19)
#
# Output Directory:
#   results/07_differential_MEDIPS_INPUT_normalized/
#   (Does NOT overwrite original results in 07_differential_MEDIPS/)
#
# Runtime: ~2-4 hours (genome-wide analysis)
# Memory: 64GB recommended
################################################################################

echo "========================================"
echo "meDIP-seq: MEDIPS Analysis WITH INPUT"
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

echo "========================================"
echo "INPUT NORMALIZATION ENABLED"
echo "========================================"
echo "This analysis uses a COMMON INPUT sample for normalization:"
echo "  TES-1-IP, TES-2-IP  <->  TESmut-1-INPUT"
echo "  GFP-1-IP, GFP-2-IP  <->  TESmut-1-INPUT"
echo ""
echo "Rationale: TES-1-INPUT does not exist, GFP-1-INPUT not processed."
echo "           TESmut-1-INPUT used as common background (same cell line)."
echo ""
echo "Output will be saved to:"
echo "  results/07_differential_MEDIPS_INPUT_normalized/"
echo ""
echo "(Original results preserved in 07_differential_MEDIPS/)"
echo "========================================"
echo ""

# Run MEDIPS analysis WITH INPUT
echo "Starting MEDIPS differential methylation analysis WITH INPUT normalization..."
echo "This will take 2-4 hours for genome-wide analysis"
echo ""

Rscript "07_differential_methylation_MEDIPS_with_INPUT.R"

EXIT_CODE=$?

echo ""
echo "========================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "MEDIPS analysis WITH INPUT completed successfully!"
    echo ""
    echo "Results saved to:"
    echo "  results/07_differential_MEDIPS_INPUT_normalized/"
    echo ""
    echo "Key output files:"
    echo "  - TES_vs_GFP_DMRs_FDR05_FC2.csv  (significant DMRs)"
    echo "  - TES_vs_GFP_all_windows.csv     (all genome-wide results)"
    echo "  - PCA_plot.pdf                    (sample clustering)"
    echo "  - *_MA_plot.pdf                   (fold change vs expression)"
    echo "  - *_volcano_plot.pdf              (statistical significance)"
else
    echo "ERROR: MEDIPS analysis failed with exit code $EXIT_CODE"
fi
echo "End time: $(date)"
echo "========================================"

exit $EXIT_CODE
