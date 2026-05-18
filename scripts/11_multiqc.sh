#!/bin/bash
#SBATCH --job-name=medip_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/11_multiqc.out
#SBATCH --error=logs/11_multiqc.err

################################################################################
# Script: 11_multiqc.sh
# Purpose: Generate comprehensive QC report with MultiQC
#
# Description:
#   Aggregates QC metrics from all pipeline steps into a single interactive
#   HTML report for easy quality assessment and cross-sample comparison.
#
# MultiQC aggregates:
#   - FastQC: Raw and trimmed read quality
#   - Trim Galore: Trimming statistics
#   - Bowtie2: Alignment rates
#   - Samtools: BAM statistics
#   - Picard: Duplication metrics
#   - MACS2: Peak calling statistics
#
# Output: results/multiqc_report.html
################################################################################

set -euo pipefail
echo "=============================================="
echo "meDIP-seq Pipeline: MultiQC Summary"
echo "=============================================="
echo "Start: $(date)"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate multiqc

echo "MultiQC version:"
multiqc --version

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
RESULTS_DIR="${BASE_DIR}/results"
OUT_DIR="${BASE_DIR}/results"

cd "${BASE_DIR}"

echo ""
echo "Running MultiQC..."
echo "Searching for QC files in: ${RESULTS_DIR}"
echo ""

multiqc \
    "${RESULTS_DIR}" \
    --outdir "${OUT_DIR}" \
    --filename medip_multiqc_report \
    --title "meDIP-seq Quality Control Report" \
    --comment "Multi-omics epigenomics project: TES methylation analysis" \
    --force \
    --verbose

if [ $? -eq 0 ]; then
    echo ""
    echo "=============================================="
    echo "MultiQC completed successfully!"
    echo "=============================================="
    echo ""
    echo "Report: ${OUT_DIR}/medip_multiqc_report.html"
    echo ""
    echo "Open in browser to review:"
    echo "  - Sample quality metrics"
    echo "  - Alignment statistics"
    echo "  - Duplication rates"
    echo "  - Cross-sample comparisons"
    echo ""
else
    echo "ERROR: MultiQC failed"
    exit 1
fi

echo "End: $(date)"
