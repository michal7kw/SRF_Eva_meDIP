#!/bin/bash
#SBATCH --job-name=medip_align_array
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=1-5                # Process samples 1-5 in parallel (TESmut IP excluded, TESmut-1-INPUT kept as common control)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          # Reduced from 32 since we're parallel
#SBATCH --mem=64G                  # Reduced from 128G per job
#SBATCH --time=08:00:00            # Reduced time per sample
#SBATCH --output=logs/03_align_array_%a.out
#SBATCH --error=logs/03_align_array_%a.err

################################################################################
# Script: 03_align_array.sh
# Purpose: PARALLEL alignment of meDIP-seq reads using SLURM job arrays
#
# Description:
#   This is a parallelized version of 03_align.sh that processes each sample
#   as a separate SLURM array job. This dramatically reduces total runtime:
#
#   Sequential (03_align.sh):     3-5 hours for all 5 samples
#   Parallel (03_align_array.sh): ~30-45 minutes (all samples run simultaneously)
#
# How Job Arrays Work:
#   - SLURM launches 5 independent jobs simultaneously (array=1-5)
#   - Each job processes one sample (SLURM_ARRAY_TASK_ID selects which sample)
#   - Jobs run on different compute nodes (if available)
#   - %A = Job array ID (same for all), %a = Task ID (1-5)
#   - NOTE: TESmut IP samples excluded from analysis (failed sample)
#   - NOTE: TESmut-1-INPUT is retained as common INPUT control
#
# Resource Allocation:
#   - 8 CPUs per sample (vs 32 in sequential version)
#   - 32 GB RAM per sample (vs 128 GB total)
#   - Total cluster usage: 5 jobs × 8 CPUs = 40 CPUs (if nodes available)
#   - More efficient use of cluster resources
#
# Advantages over Sequential:
#   ✓ 8-12x faster total runtime (parallel processing)
#   ✓ Better cluster resource utilization
#   ✓ Isolated failures (one sample failing doesn't affect others)
#   ✓ Easier to rerun failed samples (just resubmit specific array indices)
#
# Disadvantages:
#   - Uses more total cluster resources simultaneously
#   - May wait in queue if cluster is busy
#   - Need to aggregate results after completion
#
# When to Use:
#   - Use this version when cluster has available nodes
#   - Use sequential version (03_align.sh) if cluster is very busy
#   - Use sequential for debugging (easier to follow single job log)
################################################################################

set -euo pipefail

echo "========================================="
echo "meDIP-seq ARRAY Alignment - Sample ${SLURM_ARRAY_TASK_ID}"
echo "========================================="
echo "Job Array ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment

echo "Loaded environment: alignment"
echo "Bowtie2 version: $(bowtie2 --version | head -1)"
echo "Samtools version: $(samtools --version | head -1)"
echo ""

# Define directories and files
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
INPUT_DIR="${BASE_DIR}/results/02_trimmed"
OUTPUT_DIR="${BASE_DIR}/results/03_aligned"
CONFIG_FILE="${BASE_DIR}/config/samples.txt"

# Reference genome (Bowtie2 index)
GENOME_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/GRCh38"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if genome index exists
if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    echo "ERROR: Bowtie2 index not found: ${GENOME_INDEX}"
    echo "Please create Bowtie2 index for hg38 genome"
    exit 1
fi

# Get sample name for this array task
# Extract line number matching SLURM_ARRAY_TASK_ID from config file
SAMPLE=$(grep -v "^#" "${CONFIG_FILE}" | grep -v "^$" | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "${SAMPLE}" ]; then
    echo "ERROR: Could not find sample for array task ${SLURM_ARRAY_TASK_ID}"
    echo "Check that config file has ${SLURM_ARRAY_TASK_ID} samples"
    exit 1
fi

echo "========================================="
echo "Processing: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}/5"
echo ""

# Bowtie2 alignment parameters
THREADS=${SLURM_CPUS_PER_TASK}
MAX_INSERT=700
MIN_INSERT=0

echo "Bowtie2 parameters:"
echo "  Preset: --very-sensitive"
echo "  Max insert size: ${MAX_INSERT} bp"
echo "  Min insert size: ${MIN_INSERT} bp"
echo "  Threads: ${THREADS}"
echo ""

# Define input files (trimmed FASTQ)
R1="${INPUT_DIR}/${SAMPLE}_R1_001_val_1.fq.gz"
R2="${INPUT_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"

# Check if input files exist
if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "ERROR: Trimmed files not found for ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    exit 1
fi

echo "Input files:"
echo "  R1: ${R1}"
echo "  R2: ${R2}"
echo ""

# Define output files
SAM="${OUTPUT_DIR}/${SAMPLE}.sam"
BAM_UNSORTED="${OUTPUT_DIR}/${SAMPLE}.bam"
BAM_SORTED="${OUTPUT_DIR}/${SAMPLE}.sorted.bam"
STATS="${OUTPUT_DIR}/${SAMPLE}.align_stats.txt"

# Run Bowtie2 alignment
echo "Running Bowtie2 alignment..."

bowtie2 \
    --very-sensitive \
    -X ${MAX_INSERT} \
    -I ${MIN_INSERT} \
    --no-mixed \
    --no-discordant \
    --threads ${THREADS} \
    -x "${GENOME_INDEX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    -S "${SAM}" \
    2> "${STATS}"

# Check if alignment succeeded
if [ $? -ne 0 ]; then
    echo "ERROR: Bowtie2 alignment failed for ${SAMPLE}"
    rm -f "${SAM}"
    exit 1
fi

echo "Alignment completed"
echo ""

# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view \
    -@ ${THREADS} \
    -bS \
    "${SAM}" \
    -o "${BAM_UNSORTED}"

if [ $? -ne 0 ]; then
    echo "ERROR: SAM to BAM conversion failed for ${SAMPLE}"
    rm -f "${SAM}" "${BAM_UNSORTED}"
    exit 1
fi

# Remove SAM file to save space
rm -f "${SAM}"
echo "SAM file removed to save space"
echo ""

# Sort BAM by coordinates
echo "Sorting BAM file by coordinates..."
samtools sort \
    -@ ${THREADS} \
    -m 3G \
    "${BAM_UNSORTED}" \
    -o "${BAM_SORTED}"

if [ $? -ne 0 ]; then
    echo "ERROR: BAM sorting failed for ${SAMPLE}"
    rm -f "${BAM_UNSORTED}" "${BAM_SORTED}"
    exit 1
fi

# Remove unsorted BAM to save space
rm -f "${BAM_UNSORTED}"
echo "Unsorted BAM removed to save space"
echo ""

# Index sorted BAM
echo "Indexing sorted BAM..."
samtools index "${BAM_SORTED}"

if [ $? -ne 0 ]; then
    echo "WARNING: BAM indexing failed for ${SAMPLE}"
fi

echo ""

# Print alignment statistics
echo "========================================="
echo "Alignment Statistics: ${SAMPLE}"
echo "========================================="
cat "${STATS}"
echo ""

# BAM statistics
echo "BAM File Statistics:"
samtools flagstat "${BAM_SORTED}" | tee "${OUTPUT_DIR}/${SAMPLE}.flagstat.txt"
echo ""

# Get file size
SIZE=$(stat -c%s "${BAM_SORTED}" | numfmt --to=iec-i --suffix=B)
echo "Output file size: ${SIZE}"
echo ""

# Extract alignment rate for summary
ALIGN_RATE=$(grep "overall alignment rate" "${STATS}" | awk '{print $1}')
echo "Overall alignment rate: ${ALIGN_RATE}"
echo ""

echo "========================================="
echo "Sample ${SAMPLE} completed successfully!"
echo "End time: $(date)"
echo "========================================="

################################################################################
# Usage Instructions:
#
# 1. Submit all samples as array job:
#    sbatch scripts/03_align_array.sh
#
# 2. Monitor all jobs:
#    squeue -u $USER
#
# 3. Check individual sample logs:
#    tail -f logs/03_align_array_JOBID_1.out  # Sample 1 (TES-1-IP)
#    tail -f logs/03_align_array_JOBID_2.out  # Sample 2 (TES-2-IP)
#    ...etc
#
# 4. Rerun specific failed samples (e.g., sample 3 and 5):
#    sbatch --array=3,5 scripts/03_align_array.sh
#
# 5. Cancel all array jobs:
#    scancel JOBID
#
# 6. Cancel specific array task (e.g., task 3):
#    scancel JOBID_3
#
# Output files:
#   - BAM files: results/03_aligned/*.sorted.bam
#   - Indexes: results/03_aligned/*.sorted.bam.bai
#   - Stats: results/03_aligned/*.align_stats.txt
#   - Flagstat: results/03_aligned/*.flagstat.txt
#
# Expected completion time:
#   - All 5 samples: 30-45 minutes (parallel)
#   - vs 3-5 hours (sequential)
#
################################################################################
