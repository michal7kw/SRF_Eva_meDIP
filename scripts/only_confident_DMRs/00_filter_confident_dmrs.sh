#!/bin/bash
#SBATCH --job-name=filter_confident_dmrs
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=logs/00_filter_confident_dmrs.out
#SBATCH --error=logs/00_filter_confident_dmrs.err

# ============================================================================
# Filter DMRs to retain only high-confidence regions
# ============================================================================
# This script filters DMRs where BOTH TES and GFP samples have >2 mean reads,
# addressing the GFP library quality issues that cause artifactual
# "hypermethylation" calls.
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

# Create logs and output directories
mkdir -p logs
mkdir -p output/07_differential_MEDIPS_confident

echo "=============================================="
echo "High-Confidence DMR Filtering"
echo "Started: $(date)"
echo "=============================================="

# Activate R environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run R script
Rscript 00_filter_confident_dmrs.R

echo ""
echo "=============================================="
echo "Filtering complete!"
echo "Finished: $(date)"
echo "=============================================="
