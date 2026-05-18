#!/bin/bash
#SBATCH --job-name=medip_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_bigwig.out
#SBATCH --error=logs/05_bigwig.err

################################################################################
# Script: 05_bigwig.sh
# Purpose: Generate normalized BigWig coverage tracks for visualization
#
# Description:
#   This script converts filtered BAM files to BigWig format with RPKM
#   normalization for genome browser visualization (IGV, UCSC). BigWig files
#   provide continuous coverage values that show methylation enrichment across
#   the genome.
#
# Why BigWig Tracks:
#   - Compact format (~100x smaller than BAM)
#   - Fast random access for genome browsers
#   - Normalized coverage enables cross-sample comparison
#   - Essential for visual validation of peaks
#   - Shares tracks with collaborators easily
#
# Normalization Strategy: RPKM (Reads Per Kilobase per Million mapped reads)
#   Formula: (reads / (fragment_length/1000)) / (total_reads/1000000)
#
#   Why RPKM for meDIP-seq:
#   - Normalizes for sequencing depth (different samples have different read counts)
#   - Normalizes for region size (longer regions get more reads)
#   - Allows direct comparison of enrichment across samples
#   - Standard for epigenomics visualization
#
# Alternative Normalizations (not used here):
#   - TPM: Total Per Million (similar to RPKM)
#   - CPM: Counts Per Million (no length normalization)
#   - Raw counts: No normalization (not comparable across samples)
#
# Tool: deepTools bamCoverage
#   - Industry standard for ChIP-seq/meDIP-seq track generation
#   - Built-in normalization methods
#   - Efficient multi-threading
#   - Handles both single-end and paired-end data
#
# Input:
#   - Filtered BAM files from ../results/04_filtered/
#
# Output:
#   - BigWig tracks: ../results/05_bigwig/*_RPKM.bw
#   - Track statistics: *_bigwig_stats.txt
#
# Expected Results:
#   - File size: 50-200 MB per sample
#   - Smooth coverage across genome
#   - Visible enrichment at CpG islands and promoters (IP samples)
#   - INPUT samples show flat baseline coverage
#
# Runtime: ~2-4 hours for all samples
# Memory: ~4-8 GB per sample
################################################################################

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

echo "========================================="
echo "meDIP-seq Pipeline - Step 05: BigWig Generation"
echo "========================================="
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig-generation

echo "Loaded environment: seurat_full2"
echo "deepTools version:"
bamCoverage --version 2>&1 | head -1 || echo "bamCoverage loaded"
echo ""

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
INPUT_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/05_bigwig"
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

# bamCoverage parameters
THREADS=${SLURM_CPUS_PER_TASK}
BIN_SIZE=10                # Resolution in base pairs (10 bp for high resolution)
NORMALIZATION="RPKM"       # Normalization method
EFFECTIVE_GENOME_SIZE=2913022398  # hg38 effective genome size (non-N bases)

echo "bamCoverage parameters:"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Normalization: ${NORMALIZATION}"
echo "  Effective genome size: ${EFFECTIVE_GENOME_SIZE}"
echo "  Threads: ${THREADS}"
echo ""

# Process each sample
for SAMPLE in ${SAMPLES}; do
    echo "========================================="
    echo "Processing: ${SAMPLE}"
    echo "Time: $(date)"
    echo ""

    # Define input file
    INPUT_BAM="${INPUT_DIR}/${SAMPLE}_filtered_dedup.bam"

    # Check if input exists
    if [ ! -f "${INPUT_BAM}" ]; then
        echo "WARNING: Input BAM not found: ${INPUT_BAM}"
        echo "Skipping..."
        echo ""
        continue
    fi

    echo "Input: ${INPUT_BAM}"
    echo ""

    # Check if BAM is indexed
    if [ ! -f "${INPUT_BAM}.bai" ]; then
        echo "BAM index not found. Creating index..."
        samtools index "${INPUT_BAM}"
    fi

    # Define output file
    OUTPUT_BW="${OUTPUT_DIR}/${SAMPLE}_RPKM.bw"
    STATS_FILE="${OUTPUT_DIR}/${SAMPLE}_bigwig_stats.txt"

    # Run bamCoverage
    echo "Generating BigWig track..."

    bamCoverage \
        --bam "${INPUT_BAM}" \
        --outFileName "${OUTPUT_BW}" \
        --outFileFormat bigwig \
        --binSize ${BIN_SIZE} \
        --normalizeUsing ${NORMALIZATION} \
        --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
        --numberOfProcessors ${THREADS} \
        --extendReads \
        --verbose

    # bamCoverage options explained:
    #   --extendReads: Extends reads to fragment length (important for paired-end)
    #   --binSize 10: 10 bp resolution (balance between file size and detail)
    #   --normalizeUsing RPKM: Normalize by depth and region size
    #   --effectiveGenomeSize: Non-N bases in hg38 (for accurate normalization)

    if [ $? -ne 0 ]; then
        echo "ERROR: bamCoverage failed for ${SAMPLE}"
        rm -f "${OUTPUT_BW}"
        continue
    fi

    echo "BigWig generation completed"
    echo ""

    # Get file sizes for comparison
    BAM_SIZE=$(stat -c%s "${INPUT_BAM}" | numfmt --to=iec-i --suffix=B)
    BW_SIZE=$(stat -c%s "${OUTPUT_BW}" | numfmt --to=iec-i --suffix=B)

    echo "File sizes:"
    echo "  Input BAM: ${BAM_SIZE}"
    echo "  Output BigWig: ${BW_SIZE}"
    echo ""

    # Get read count from BAM
    TOTAL_READS=$(samtools view -c -F 0x4 "${INPUT_BAM}")

    # Write statistics
    cat > "${STATS_FILE}" << EOF
Sample: ${SAMPLE}
Date: $(date)

Input:
  BAM file: ${INPUT_BAM}
  BAM size: ${BAM_SIZE}
  Total reads: ${TOTAL_READS}

Output:
  BigWig file: ${OUTPUT_BW}
  BigWig size: ${BW_SIZE}

Parameters:
  Bin size: ${BIN_SIZE} bp
  Normalization: ${NORMALIZATION}
  Effective genome size: ${EFFECTIVE_GENOME_SIZE}
  Extend reads: Yes
EOF

    echo "BigWig statistics:"
    cat "${STATS_FILE}"
    echo ""

done

# Summary
echo "========================================="
echo "BigWig Generation Summary"
echo "========================================="

# Count successful outputs
NUM_BIGWIG=$(find "${OUTPUT_DIR}" -name "*_RPKM.bw" | wc -l)

echo "Expected BigWig files: ${NUM_SAMPLES}"
echo "Generated BigWig files: ${NUM_BIGWIG}"

if [ ${NUM_BIGWIG} -eq ${NUM_SAMPLES} ]; then
    echo "Status: All BigWig files generated successfully!"
else
    echo "WARNING: Some samples may have failed!"
    echo "Check logs for details."
fi

echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/*_RPKM.bw 2>/dev/null || echo "No BigWig files found!"

echo ""
echo "========================================="
echo "Visualization Instructions"
echo "========================================="
echo ""
echo "IGV (Integrative Genomics Viewer):"
echo "  1. Load hg38 genome"
echo "  2. File → Load from File → Select *_RPKM.bw files"
echo "  3. Compare IP vs INPUT samples side-by-side"
echo "  4. Navigate to known methylated regions (e.g., imprinted genes)"
echo "  5. Check for enrichment at CpG islands"
echo ""
echo "UCSC Genome Browser:"
echo "  1. Go to https://genome.ucsc.edu"
echo "  2. My Data → Custom Tracks"
echo "  3. Upload BigWig files"
echo "  4. Adjust track height and color for comparison"
echo ""
echo "Expected Patterns:"
echo "  - IP samples: Peaks at promoters and CpG islands"
echo "  - INPUT samples: Flat baseline coverage"
echo "  - TES vs GFP: Different methylation patterns"
echo "  - Good signal-to-noise: >5-fold enrichment at peaks"
echo ""

echo ""
echo "========================================="
echo "Next steps:"
echo "1. Load BigWig files in IGV or UCSC browser"
echo "2. Visually inspect enrichment patterns"
echo "3. Compare IP vs INPUT samples"
echo "4. Navigate to candidate genes from RNA-seq"
echo "5. Run step 06: Peak calling (06_peak_calling.sh)"
echo "========================================="

echo ""
echo "End time: $(date)"
echo "========================================="

# Understanding BigWig Tracks:
#
# What to Look For in IGV:
#
# 1. IP vs INPUT Comparison:
#    - IP should show clear peaks over INPUT baseline
#    - INPUT should have relatively flat coverage
#    - Good enrichment: 5-20 fold IP over INPUT at peaks
#
# 2. CpG Island Enrichment:
#    - Load CpG island track in IGV
#    - IP samples should show enrichment at CpG islands
#    - Particularly strong at gene promoters
#
# 3. Coverage Uniformity:
#    - Check multiple chromosomes
#    - Coverage should be relatively even (no huge spikes)
#    - Big spikes may indicate repetitive regions or artifacts
#
# 4. Signal-to-Noise Ratio:
#    - Background (intergenic regions): Low RPKM (0-5)
#    - Peaks (methylated regions): High RPKM (20-100+)
#    - Good SNR: >5-10 fold enrichment
#
# 5. Sample Comparison:
#    - Normalize track heights for visual comparison
#    - Look for condition-specific peaks (TES vs GFP)
#    - Check biological replicates for consistency
#
# Quality Issues to Watch For:
#
# 1. Very low signal:
#    - All RPKM values < 5
#    - May indicate failed IP or low sequencing depth
#    - Check total read counts
#
# 2. No peaks visible:
#    - IP looks like INPUT
#    - Failed antibody enrichment
#    - Wrong antibody or poor IP conditions
#
# 3. Extremely high signal:
#    - RPKM values > 500
#    - May indicate PCR artifacts or highly duplicated regions
#    - Check duplication rates from step 04
#
# 4. Uneven coverage:
#    - Some chromosomes have much higher signal
#    - May indicate copy number variation or contamination
#    - Compare with INPUT to distinguish
#
# 5. Replicates don't match:
#    - Different peak patterns between replicates
#    - Technical variability or batch effects
#    - May need to exclude bad replicates
#
# Troubleshooting:
#    - bamCoverage fails: Check BAM index exists
#    - File size too large: Increase --binSize (try 25 or 50 bp)
#    - Out of memory: Reduce --numberOfProcessors
#    - Tracks look weird: Check RPKM values are in reasonable range (1-100)
#
# Advanced Visualization:
#    - Create track hubs for sharing with collaborators
#    - Generate metaplots showing average enrichment across gene features
#    - Use deepTools computeMatrix + plotHeatmap for publication figures
