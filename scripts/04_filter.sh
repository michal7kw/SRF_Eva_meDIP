#!/bin/bash
#SBATCH --job-name=medip_filter
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/04_filter.out
#SBATCH --error=logs/04_filter.err

################################################################################
# Script: 04_filter.sh
# Purpose: Filter and deduplicate aligned meDIP-seq BAM files
#
# Description:
#   This script performs quality filtering and PCR duplicate removal on aligned
#   BAM files. It applies multiple filtering steps to ensure high-quality data
#   for peak calling and downstream analysis.
#
# Filtering Steps:
#   1. Remove unmapped reads
#   2. Remove improperly paired reads (wrong orientation/insert size)
#   3. Remove low-quality alignments (MAPQ < 30)
#   4. Remove PCR duplicates (Picard MarkDuplicates)
#
# Why Each Filter:
#
#   Unmapped reads: No genomic location, cannot contribute to peaks
#
#   Improperly paired: Indicates structural issues or mapping errors
#   - Discordant pairs (wrong orientation)
#   - Abnormal insert sizes (too large/small)
#
#   Low MAPQ (<30): Ambiguous alignments
#   - MAPQ 30 = 0.1% chance of incorrect mapping
#   - Removes multimappers and poor alignments
#   - Critical for accurate peak calling
#
#   PCR duplicates: Technical artifacts from library amplification
#   - Inflate read counts artificially
#   - Create false peaks or exaggerate real ones
#   - 10-30% duplication typical for meDIP-seq
#   - Marking (not removing) preserves counts for QC
#
# Input:
#   - Sorted BAM files from ../results/03_aligned/
#
# Output:
#   - Filtered & deduplicated BAM: ../results/04_filtered/*_filtered_dedup.bam
#   - BAM index: *_filtered_dedup.bam.bai
#   - Duplication metrics: *_dup_metrics.txt
#   - Filtering statistics: *_filter_stats.txt
#
# Expected Results:
#   - Duplication rate: 10-30% (higher suggests over-amplification)
#   - Reads retained: 60-85% of aligned reads
#   - MAPQ distribution: Peak at 42 (unique alignments)
#
# Runtime: ~4-6 hours for all samples
# Memory: ~8-16 GB per sample (32 GB safe)
################################################################################

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

echo "========================================="
echo "meDIP-seq Pipeline - Step 04: Filtering"
echo "========================================="
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

echo "Loaded environment: alignment"
echo "Samtools version:"
samtools --version | head -1
echo "Picard version:"
java -jar /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment/share/picard-*/picard.jar MarkDuplicates --version 2>&1 | head -1 || echo "Picard loaded"
echo ""

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
INPUT_DIR="${BASE_DIR}/results/03_aligned"
OUTPUT_DIR="${BASE_DIR}/results/04_filtered"
CONFIG_FILE="${BASE_DIR}/config/samples.txt"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if config file exists
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Sample config file not found: ${CONFIG_FILE}"
    exit 1
fi

# Read sample names
SAMPLES=$(grep -v "^#" "${CONFIG_FILE}" | grep -v "^$" | awk '{print $1}')
NUM_SAMPLES=$(echo "${SAMPLES}" | wc -w)

echo "Processing ${NUM_SAMPLES} samples"
echo "CPU cores: ${SLURM_CPUS_PER_TASK}"
echo ""

# Filtering parameters
MAPQ_THRESHOLD=30      # Minimum mapping quality (99.9% confidence)
THREADS=${SLURM_CPUS_PER_TASK}

echo "Filtering parameters:"
echo "  MAPQ threshold: ${MAPQ_THRESHOLD}"
echo "  Remove unmapped: Yes"
echo "  Remove improperly paired: Yes"
echo "  Remove duplicates: Yes (marked and removed)"
echo ""

# Find Picard JAR
PICARD_JAR=$(find /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment/share -name "picard.jar" | head -1)

if [ -z "${PICARD_JAR}" ]; then
    echo "ERROR: Picard JAR not found"
    exit 1
fi

echo "Picard JAR: ${PICARD_JAR}"
echo ""

# Process each sample
for SAMPLE in ${SAMPLES}; do
    echo "========================================="
    echo "Processing: ${SAMPLE}"
    echo "Time: $(date)"
    echo ""

    # Define input file
    INPUT_BAM="${INPUT_DIR}/${SAMPLE}.sorted.bam"

    # Check if input exists
    if [ ! -f "${INPUT_BAM}" ]; then
        echo "WARNING: Input BAM not found: ${INPUT_BAM}"
        echo "Skipping..."
        echo ""
        continue
    fi

    echo "Input: ${INPUT_BAM}"
    echo ""

    # Define intermediate and output files
    FILTERED_BAM="${OUTPUT_DIR}/${SAMPLE}_filtered.bam"
    RGTAG_BAM="${OUTPUT_DIR}/${SAMPLE}_filtered_rg.bam"
    MARKDUP_BAM="${OUTPUT_DIR}/${SAMPLE}_filtered_markdup.bam"
    FINAL_BAM="${OUTPUT_DIR}/${SAMPLE}_filtered_dedup.bam"
    DUP_METRICS="${OUTPUT_DIR}/${SAMPLE}_dup_metrics.txt"
    FILTER_STATS="${OUTPUT_DIR}/${SAMPLE}_filter_stats.txt"

    # Step 1: Filter for quality, proper pairing, and mapping quality
    echo "Step 1: Filtering for quality and proper pairing..."
    echo "  - Removing unmapped reads"
    echo "  - Removing improperly paired reads"
    echo "  - Filtering MAPQ < ${MAPQ_THRESHOLD}"
    echo ""

    samtools view \
        -@ ${THREADS} \
        -b \
        -f 2 \
        -F 1804 \
        -q ${MAPQ_THRESHOLD} \
        "${INPUT_BAM}" \
        -o "${FILTERED_BAM}"

    # SAM flags explained:
    # -f 2: Require properly paired (0x2)
    # -F 1804: Filter out:
    #   - 0x4: Read unmapped
    #   - 0x8: Mate unmapped
    #   - 0x400: Duplicate (from previous runs)
    #   - 0x800: Supplementary alignment
    # -q 30: Minimum MAPQ of 30

    if [ $? -ne 0 ]; then
        echo "ERROR: Filtering failed for ${SAMPLE}"
        rm -f "${FILTERED_BAM}"
        continue
    fi

    echo "Filtering completed"
    echo ""

    # Get initial read count
    READS_BEFORE=$(samtools view -c "${INPUT_BAM}")
    READS_AFTER_FILTER=$(samtools view -c "${FILTERED_BAM}")

    echo "Reads before filtering: ${READS_BEFORE}"
    echo "Reads after filtering: ${READS_AFTER_FILTER}"
    echo ""

    # Step 2: Add read groups (required by Picard)
    echo "Step 2: Adding read group information..."
    echo "  Sample: ${SAMPLE}"
    echo "  Library: ${SAMPLE}"
    echo ""

    java -Xmx8g -jar "${PICARD_JAR}" AddOrReplaceReadGroups \
        INPUT="${FILTERED_BAM}" \
        OUTPUT="${RGTAG_BAM}" \
        RGID="${SAMPLE}" \
        RGLB="${SAMPLE}" \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM="${SAMPLE}" \
        VALIDATION_STRINGENCY=LENIENT \
        SORT_ORDER=coordinate

    if [ $? -ne 0 ]; then
        echo "ERROR: AddOrReplaceReadGroups failed for ${SAMPLE}"
        rm -f "${FILTERED_BAM}" "${RGTAG_BAM}"
        continue
    fi

    echo "Read group addition completed"
    echo ""

    # Step 3: Mark PCR duplicates with Picard
    echo "Step 3: Marking PCR duplicates with Picard..."

    java -Xmx16g -jar "${PICARD_JAR}" MarkDuplicates \
        INPUT="${RGTAG_BAM}" \
        OUTPUT="${MARKDUP_BAM}" \
        METRICS_FILE="${DUP_METRICS}" \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=LENIENT \
        CREATE_INDEX=false \
        ASSUME_SORT_ORDER=coordinate

    if [ $? -ne 0 ]; then
        echo "ERROR: MarkDuplicates failed for ${SAMPLE}"
        rm -f "${FILTERED_BAM}" "${RGTAG_BAM}" "${MARKDUP_BAM}" "${DUP_METRICS}"
        continue
    fi

    echo "Duplicate marking completed"
    echo ""

    # Step 4: Remove marked duplicates
    echo "Step 4: Removing marked duplicates..."

    samtools view \
        -@ ${THREADS} \
        -b \
        -F 1024 \
        "${MARKDUP_BAM}" \
        -o "${FINAL_BAM}"

    # -F 1024: Filter out duplicate flag (0x400)

    if [ $? -ne 0 ]; then
        echo "ERROR: Duplicate removal failed for ${SAMPLE}"
        rm -f "${FILTERED_BAM}" "${RGTAG_BAM}" "${MARKDUP_BAM}" "${FINAL_BAM}"
        continue
    fi

    echo "Duplicate removal completed"
    echo ""

    # Step 5: Index final BAM
    echo "Step 5: Indexing final BAM..."
    samtools index "${FINAL_BAM}"

    if [ $? -ne 0 ]; then
        echo "WARNING: Indexing failed for ${SAMPLE}"
    fi

    echo ""

    # Clean up intermediate files
    rm -f "${FILTERED_BAM}" "${RGTAG_BAM}" "${MARKDUP_BAM}"
    echo "Intermediate files removed"
    echo ""

    # Generate statistics
    READS_FINAL=$(samtools view -c "${FINAL_BAM}")
    READS_REMOVED=$((READS_AFTER_FILTER - READS_FINAL))
    READS_DUPLICATES=$((READS_AFTER_FILTER - READS_FINAL))

    PCT_RETAINED=$(awk "BEGIN {printf \"%.2f\", (${READS_FINAL}/${READS_BEFORE})*100}")
    PCT_FILTERED=$(awk "BEGIN {printf \"%.2f\", ((${READS_BEFORE}-${READS_AFTER_FILTER})/${READS_BEFORE})*100}")
    PCT_DUPLICATES=$(awk "BEGIN {printf \"%.2f\", (${READS_DUPLICATES}/${READS_AFTER_FILTER})*100}")

    # Write statistics file
    cat > "${FILTER_STATS}" << EOF
Sample: ${SAMPLE}
Date: $(date)

Read Counts:
  Initial reads: ${READS_BEFORE}
  After quality filtering: ${READS_AFTER_FILTER}
  After deduplication: ${READS_FINAL}

Reads Removed:
  By quality filters: $((READS_BEFORE - READS_AFTER_FILTER))
  By deduplication: ${READS_DUPLICATES}
  Total removed: $((READS_BEFORE - READS_FINAL))

Percentages:
  Retained (final/initial): ${PCT_RETAINED}%
  Filtered by quality: ${PCT_FILTERED}%
  Duplicates: ${PCT_DUPLICATES}%

Filtering Criteria:
  MAPQ threshold: ${MAPQ_THRESHOLD}
  Properly paired: Required
  Unmapped: Removed
  Duplicates: Removed
EOF

    echo "Filtering statistics:"
    cat "${FILTER_STATS}"
    echo ""

    # Extract duplication rate from Picard metrics
    if [ -f "${DUP_METRICS}" ]; then
        echo "Picard duplication metrics:"
        grep -A 1 "^LIBRARY" "${DUP_METRICS}" || echo "Metrics extraction failed"
        echo ""
    fi

    # Get file size
    SIZE=$(stat -c%s "${FINAL_BAM}" | numfmt --to=iec-i --suffix=B)
    echo "Final BAM size: ${SIZE}"
    echo ""

    # Generate flagstat
    echo "Final BAM statistics:"
    samtools flagstat "${FINAL_BAM}" | tee "${OUTPUT_DIR}/${SAMPLE}_filtered_dedup.flagstat.txt"
    echo ""

done

# Summary
echo "========================================="
echo "Filtering Summary"
echo "========================================="

# Count successful outputs
NUM_FILTERED=$(find "${OUTPUT_DIR}" -name "*_filtered_dedup.bam" | wc -l)

echo "Expected samples: ${NUM_SAMPLES}"
echo "Filtered samples: ${NUM_FILTERED}"

if [ ${NUM_FILTERED} -eq ${NUM_SAMPLES} ]; then
    echo "Status: All samples filtered successfully!"
else
    echo "WARNING: Some samples may have failed!"
    echo "Check logs for details."
fi

echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/*_filtered_dedup.bam 2>/dev/null || echo "No BAM files found!"

echo ""
echo "========================================="
echo "Filtering Statistics Summary"
echo "========================================="

# Summarize retention rates
for STATS in "${OUTPUT_DIR}"/*_filter_stats.txt; do
    if [ -f "${STATS}" ]; then
        SAMPLE=$(basename "${STATS}" _filter_stats.txt)
        RETAINED=$(grep "Retained" "${STATS}" | awk '{print $NF}')
        DUPLICATES=$(grep "Duplicates:" "${STATS}" | awk '{print $NF}')
        echo "${SAMPLE}: ${RETAINED} retained, ${DUPLICATES} duplicates"
    fi
done

echo ""
echo "========================================="
echo "Next steps:"
echo "1. Review filtering statistics (*_filter_stats.txt)"
echo "2. Check duplication rates (*_dup_metrics.txt)"
echo "3. Expected: 60-85% reads retained, 10-30% duplicates"
echo "4. Load filtered BAM in IGV to inspect quality"
echo "5. Run step 05: BigWig generation (05_bigwig.sh)"
echo "========================================="

echo ""
echo "End time: $(date)"
echo "========================================="

# Understanding Filtering Results:
#
# 1. Quality Filtering (MAPQ < 30)
#    Removes:
#      - Multimapping reads (align to multiple locations)
#      - Low-confidence alignments
#      - Randomly placed reads
#    Typical removal: 10-25% of aligned reads
#
# 2. Proper Pairing Filter
#    Removes:
#      - Discordant pairs (wrong orientation)
#      - Abnormal insert sizes
#      - Orphan reads (mate unmapped)
#    Typical removal: <5% of aligned reads
#
# 3. PCR Duplicate Removal
#    Removes:
#      - Reads with identical start positions (likely PCR artifacts)
#      - Optical duplicates (sequencer artifacts)
#    Expected duplication rates:
#      - Good library: 10-20%
#      - Acceptable: 20-30%
#      - Concerning: >30% (over-amplification)
#      - Very low (<5%): May indicate low starting material
#
# Overall Retention Rate:
#    - Good: 70-85% of initial aligned reads
#    - Acceptable: 60-70%
#    - Concerning: <60% (quality issues)
#
# Red Flags:
#    - >40% duplicates: Over-amplification, consider deeper sequencing
#    - <50% retention: Poor library quality or alignment issues
#    - Very different rates between samples: Batch effects or prep issues
#
# What to Check in IGV:
#    - Coverage should be smoother after duplicate removal
#    - Peaks should be more distinct
#    - Suspicious "tower" peaks should be reduced
#
# Troubleshooting:
#    - High duplicate rate: May need more input DNA or fewer PCR cycles
#    - Low retention: Check alignment quality and MAPQ threshold
#    - Picard fails: Increase -Xmx memory allocation
#    - Very uneven coverage: May indicate technical artifacts
