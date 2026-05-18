#!/bin/bash
#SBATCH --job-name=medip_trim
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/02_trim.out
#SBATCH --error=logs/02_trim.err

################################################################################
# Script: 02_trim.sh
# Purpose: Adapter trimming and quality filtering of meDIP-seq reads
#
# Description:
#   This script uses Trim Galore (wrapper around Cutadapt + FastQC) to:
#   1. Remove Illumina adapter sequences from read ends
#   2. Trim low-quality bases (Phred < 20)
#   3. Remove very short reads after trimming (<20 bp)
#   4. Maintain paired-end read pairing
#   5. Generate FastQC reports on trimmed data
#
# Why Trim Galore for meDIP-seq:
#   - Adapter contamination reduces alignment accuracy
#   - Low-quality bases introduce sequencing errors
#   - Improves signal-to-noise ratio for peak calling
#   - Automatic adapter detection (no need to specify sequences)
#   - Built-in paired-end mode maintains read pairing
#
# Trimming Strategy:
#   - Quality threshold: Phred 20 (99% base call accuracy)
#   - Minimum length: 20 bp (too short reads are uninformative)
#   - Stringency: 3 (adapter overlap required for trimming)
#   - Both ends trimmed independently, then paired
#
# Input:
#   - Raw FASTQ files from ../fastq/
#   - Paired-end: *_R1_001.fastq.gz, *_R2_001.fastq.gz
#
# Output:
#   - Trimmed FASTQ files: ../results/02_trimmed/*_val_1.fq.gz, *_val_2.fq.gz
#   - Trimming reports: *_trimming_report.txt
#   - FastQC reports: *_fastqc.html
#
# Expected Results:
#   - Read length reduction: ~140-150 bp (from 150 bp)
#   - Bases trimmed: 5-15% typically
#   - Reads removed: <5% (most should pass filters)
#   - Improved per-base quality scores in FastQC
#
# Runtime: ~4-8 hours for all samples
# Memory: ~2-4 GB per sample (32 GB safe for parallel processing)
################################################################################

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

echo "========================================="
echo "meDIP-seq Pipeline - Step 02: Trim Galore"
echo "========================================="
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate tg

echo "Loaded environment: tg"
echo "Trim Galore version:"
trim_galore --version
echo ""

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
FASTQ_DIR="${BASE_DIR}/fastq"
OUTPUT_DIR="${BASE_DIR}/results/02_trimmed"
CONFIG_FILE="${BASE_DIR}/config/samples.txt"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if FASTQ directory exists
if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

# Check if config file exists
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Sample config file not found: ${CONFIG_FILE}"
    exit 1
fi

# Read sample names (skip comments and empty lines)
SAMPLES=$(grep -v "^#" "${CONFIG_FILE}" | grep -v "^$" | awk '{print $1}')

# Count samples
NUM_SAMPLES=$(echo "${SAMPLES}" | wc -w)
echo "Processing ${NUM_SAMPLES} samples"
echo ""

# Trim Galore parameters
QUALITY=20          # Phred quality threshold
STRINGENCY=3        # Adapter overlap required
MIN_LENGTH=20       # Minimum read length after trimming
CORES=2             # Cores per sample (total: NUM_SAMPLES * CORES should not exceed SLURM_CPUS_PER_TASK)

echo "Trim Galore parameters:"
echo "  Quality threshold: ${QUALITY}"
echo "  Adapter stringency: ${STRINGENCY}"
echo "  Minimum length: ${MIN_LENGTH} bp"
echo "  Cores per sample: ${CORES}"
echo ""

# Process each sample
for SAMPLE in ${SAMPLES}; do
    echo "========================================="
    echo "Processing: ${SAMPLE}"
    echo "Time: $(date)"
    echo ""

    # Define input files
    R1="${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

    # Check if input files exist
    if [ ! -f "${R1}" ]; then
        echo "WARNING: R1 file not found for ${SAMPLE}: ${R1}"
        echo "Skipping..."
        echo ""
        continue
    fi

    if [ ! -f "${R2}" ]; then
        echo "WARNING: R2 file not found for ${SAMPLE}: ${R2}"
        echo "Skipping..."
        echo ""
        continue
    fi

    echo "Input files:"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    echo ""

    # Run Trim Galore
    echo "Running Trim Galore..."

    trim_galore \
        --paired \
        --quality ${QUALITY} \
        --stringency ${STRINGENCY} \
        --length ${MIN_LENGTH} \
        --fastqc \
        --cores ${CORES} \
        --output_dir "${OUTPUT_DIR}" \
        "${R1}" "${R2}"

    # Check if trimming succeeded
    if [ $? -eq 0 ]; then
        echo "Trimming completed successfully for ${SAMPLE}"

        # Check output files
        VAL1="${OUTPUT_DIR}/${SAMPLE}_R1_001_val_1.fq.gz"
        VAL2="${OUTPUT_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"

        if [ -f "${VAL1}" ] && [ -f "${VAL2}" ]; then
            # Get file sizes
            SIZE_R1_ORIG=$(stat -c%s "${R1}" | numfmt --to=iec-i --suffix=B)
            SIZE_R1_TRIM=$(stat -c%s "${VAL1}" | numfmt --to=iec-i --suffix=B)
            SIZE_R2_ORIG=$(stat -c%s "${R2}" | numfmt --to=iec-i --suffix=B)
            SIZE_R2_TRIM=$(stat -c%s "${VAL2}" | numfmt --to=iec-i --suffix=B)

            echo ""
            echo "Output files created:"
            echo "  R1: ${SIZE_R1_ORIG} -> ${SIZE_R1_TRIM}"
            echo "  R2: ${SIZE_R2_ORIG} -> ${SIZE_R2_TRIM}"
        else
            echo "WARNING: Expected output files not found!"
        fi
    else
        echo "ERROR: Trimming failed for ${SAMPLE}"
        # Continue with other samples instead of exiting
    fi

    echo ""
done

# Summary
echo "========================================="
echo "Trimming Summary"
echo "========================================="

# Count output files
NUM_TRIMMED=$(find "${OUTPUT_DIR}" -name "*_val_*.fq.gz" | wc -l)
EXPECTED=$((NUM_SAMPLES * 2))  # Each sample has R1 and R2

echo "Expected output files: ${EXPECTED}"
echo "Generated output files: ${NUM_TRIMMED}"

if [ ${NUM_TRIMMED} -eq ${EXPECTED} ]; then
    echo "Status: All samples processed successfully!"
else
    echo "WARNING: Some samples may have failed!"
    echo "Check logs for details."
fi

# List all trimmed files
echo ""
echo "Trimmed files:"
ls -lh "${OUTPUT_DIR}"/*_val_*.fq.gz 2>/dev/null || echo "No trimmed files found!"

# List trimming reports
echo ""
echo "Trimming reports:"
ls -lh "${OUTPUT_DIR}"/*_trimming_report.txt 2>/dev/null || echo "No reports found!"

# List FastQC reports
echo ""
echo "FastQC reports on trimmed data:"
ls -lh "${OUTPUT_DIR}"/*_fastqc.html 2>/dev/null || echo "No FastQC reports found!"

echo ""
echo "========================================="
echo "Next steps:"
echo "1. Review trimming reports: ${OUTPUT_DIR}/*_trimming_report.txt"
echo "2. Compare FastQC reports before/after trimming"
echo "3. Check that most reads (>95%) passed filters"
echo "4. Run step 03: Alignment (03_align.sh)"
echo "========================================="

echo ""
echo "End time: $(date)"
echo "========================================="

# Understanding Trimming Reports:
#
# Key metrics in *_trimming_report.txt:
#
# 1. Total reads processed
#    - Should match FastQC read count
#
# 2. Reads with adapters
#    - Percentage with Illumina adapters detected
#    - Typical: 20-80% (normal for meDIP-seq)
#
# 3. Reads that were too short
#    - Discarded because <20 bp after trimming
#    - Should be <5% (if higher, check library quality)
#
# 4. Quality-trimmed bases
#    - Total bases removed due to low quality
#    - Typical: 5-15% of total bases
#
# 5. Adapter-trimmed bases
#    - Total bases removed as adapters
#    - Depends on insert size distribution
#
# 6. Total written (filtered)
#    - Reads passing all filters
#    - Should be >95% of input reads
#
# Interpreting Results:
#
# Good trimming:
#   - >95% reads retained
#   - 5-15% bases trimmed
#   - Improved quality scores in FastQC
#   - Adapter content reduced to near zero
#
# Problems to watch for:
#   - >20% reads discarded: Library quality issue
#   - <5% bases trimmed: May not be improving quality
#   - >50% reads with adapters: Short insert sizes
#   - Very uneven R1 vs R2 trimming: Asymmetric quality issues
#
# Troubleshooting:
#   - If many reads too short: Increase --length threshold or check library prep
#   - If little trimming: Data may already be high quality (good!)
#   - If R1 and R2 very different: Check sequencing quality per read
#   - If trim_galore fails: Check disk space and file permissions
