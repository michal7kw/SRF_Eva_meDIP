#!/bin/bash
#SBATCH --job-name=medip_fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/01_fastqc.out
#SBATCH --error=logs/01_fastqc.err

################################################################################
# Script: 01_fastqc.sh
# Purpose: Quality control assessment of raw meDIP-seq FASTQ files
#
# Description:
#   This script performs quality control analysis on raw sequencing data using
#   FastQC. It processes all paired-end FASTQ files (R1 and R2) and generates
#   HTML reports with quality metrics.
#
# Why FastQC for meDIP-seq:
#   - Identifies adapter contamination that can reduce alignment accuracy
#   - Detects low-quality bases that should be trimmed
#   - Checks for unexpected GC content biases
#   - Identifies over-represented sequences (PCR artifacts)
#   - Essential QC before proceeding with alignment
#
# Input:
#   - Raw FASTQ files in ../fastq/ directory
#   - Paired-end format: *_R1_001.fastq.gz, *_R2_001.fastq.gz
#
# Output:
#   - FastQC HTML reports in ../results/01_fastqc/
#   - Per-sample quality metrics
#   - Zipped data files with detailed statistics
#
# Expected Quality Metrics:
#   - Per base sequence quality: Phred score >30 for most positions
#   - Per sequence quality score: Peak at high quality (>30)
#   - Per base GC content: Relatively stable across read length
#   - Sequence length: 150 bp (typical for NovaSeq)
#   - Adapter content: Should be present (trimmed in next step)
#
# Runtime: ~2-4 hours for all samples
# Memory: ~2GB per sample (16GB total safe for parallel processing)
################################################################################

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

echo "========================================="
echo "meDIP-seq Pipeline - Step 01: FastQC"
echo "========================================="
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate quality

echo "Loaded environment: quality"
echo "FastQC version:"
fastqc --version
echo ""

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
FASTQ_DIR="${BASE_DIR}/fastq"
OUTPUT_DIR="${BASE_DIR}/results/01_fastqc"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if FASTQ directory exists
if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

# Count FASTQ files
NUM_FILES=$(find "${FASTQ_DIR}" -name "*.fastq.gz" | wc -l)
echo "Found ${NUM_FILES} FASTQ files to process"
echo ""

# Run FastQC on all FASTQ files
echo "Running FastQC..."
echo "Output directory: ${OUTPUT_DIR}"
echo ""

fastqc \
    --outdir "${OUTPUT_DIR}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --quiet \
    "${FASTQ_DIR}"/*.fastq.gz

# Check if FastQC completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "FastQC completed successfully!"

    # Count output files
    NUM_HTML=$(find "${OUTPUT_DIR}" -name "*_fastqc.html" | wc -l)
    echo "Generated ${NUM_HTML} HTML reports"

    # List output files
    echo ""
    echo "Generated reports:"
    ls -lh "${OUTPUT_DIR}"/*_fastqc.html

    echo ""
    echo "========================================="
    echo "Next steps:"
    echo "1. Review FastQC HTML reports in: ${OUTPUT_DIR}"
    echo "2. Check for adapter contamination and quality issues"
    echo "3. Run step 02: Adapter trimming (02_trim.sh)"
    echo "========================================="
else
    echo "ERROR: FastQC failed!"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "========================================="

# Key Quality Metrics to Review in FastQC Reports:
#
# 1. Per base sequence quality
#    - Green: Good quality (Phred >28)
#    - Orange/Red: Poor quality (Phred <20)
#    - Action: If many red regions, consider more aggressive trimming
#
# 2. Per sequence quality scores
#    - Peak should be at high Phred scores (>30)
#    - Broad or low-quality peak indicates problems
#
# 3. Per base sequence content
#    - First ~10-15 bp often show bias (random priming artifact)
#    - Normal for meDIP-seq, will be trimmed
#
# 4. Per sequence GC content
#    - Should show single peak around 45-50% (human genome)
#    - Multiple peaks may indicate contamination
#
# 5. Sequence Duplication Levels
#    - meDIP-seq expects moderate duplication (enrichment assay)
#    - 10-30% duplication is normal
#    - >50% suggests over-amplification or low complexity
#
# 6. Adapter Content
#    - Common adapters: Illumina Universal, Nextera
#    - Peaks at 3' end of reads
#    - Will be removed by Trim Galore in next step
#
# 7. Overrepresented sequences
#    - Check if they match adapters
#    - Genomic sequences may indicate PCR artifacts
#
# 8. Read length
#    - Should be consistent (150 bp typical)
#    - Will decrease after trimming
#
# Troubleshooting:
# - If fastqc fails: Check FASTQ file integrity (try unzipping one file)
# - If low quality throughout: May need to exclude samples from analysis
# - If high adapter content: Normal, will be trimmed in step 02
# - If very low GC content: Check for contamination or wrong reference genome
