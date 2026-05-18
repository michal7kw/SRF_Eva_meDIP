#!/bin/bash
#SBATCH --job-name=medip_filter_bounds
#SBATCH --account=kubacki.michal
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/04b_filter_bounds.out
#SBATCH --error=logs/04b_filter_bounds.err

# Script to filter reads that are outside the bounds of standard GRCh38
# This fixes the R mismatch error: "supplied start is > refwidth"

set -euo pipefail

# Load environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
INPUT_DIR="${BASE_DIR}/results/04_filtered"
# We will overwrite the input files or create new ones? 
# Safer to create a backup first, then overwrite, so downstream scripts work unchanged.
BACKUP_DIR="${BASE_DIR}/results/04_filtered/backup_original"

mkdir -p "${BACKUP_DIR}"

# Define standard GRCh38 lengths (from BSgenome.Hsapiens.NCBI.GRCh38)
# We will use an awk script to filter BAMs based on these lengths
# Load chromosome sizes from file
CHROM_SIZES_FILE="${BASE_DIR}/scripts/GRCh38_EBV.chrom.sizes.tsv"

if [ ! -f "${CHROM_SIZES_FILE}" ]; then
    echo "ERROR: Chromosome sizes file not found: ${CHROM_SIZES_FILE}"
    exit 1
fi

# Create a clean version for awk (removing 'chr' prefix to match BAM if needed)
# Checking first few lines of BAM to see if they use 'chr' or just numbers
# The previous `samtools view` showed "SN:1", "SN:10", etc. (no 'chr' prefix)
# The provided file has "chr1", "chr2", etc.
# We need to strip 'chr' prefix from the sizes file for it to match the BAM SN names!

sed 's/^chr//' "${CHROM_SIZES_FILE}" > chrom.sizes

# Function to filter a single BAM
filter_bam() {
    local BAM="$1"
    local SAMPLE=$(basename "${BAM}" _filtered_dedup.bam)
    local TEMP_BAM="${BAM}.temp.bam"
    
    echo "Processing ${SAMPLE}..."
    
    # Backup if not exists
    if [ ! -f "${BACKUP_DIR}/${SAMPLE}_filtered_dedup.bam" ]; then
        cp "${BAM}" "${BACKUP_DIR}/"
        echo "  Backed up original to ${BACKUP_DIR}"
    fi

    # Filter using samtools and awk
    # Strategy: 
    # 1. Convert to SAM
    # 2. Filter lines where pos + length > chrom_length
    # 3. Convert back to BAM
    
    echo "  Filtering out-of-bounds reads..."
    
    samtools view -h "${BAM}" | \
    awk '
    BEGIN {
        # Load chromosome sizes
        while ((getline < "chrom.sizes") > 0) {
            sizes[$1] = $2
        }
    }
    {
        if (/^@/) { print; next } # Print header
        
        ref = $3
        pos = $4
        # Calculate read length from CIGAR or sequence
        # Simplified: use sequence length
        len = length($10)
        
        # Check bounds
        if (ref in sizes) {
            if (pos + len <= sizes[ref]) {
                print
            }
        } else {
            # Keep reads on other contigs (or drop them? Let is drop them to be safe for R)
            # print
        }
    }' | \
    samtools view -b - > "${TEMP_BAM}"
    
    # Replace original
    mv "${TEMP_BAM}" "${BAM}"
    samtools index "${BAM}"
    
    echo "  Done. Filtered BAM saved to ${BAM}"
}

# Run for all BAM files
for BAM in "${INPUT_DIR}"/*_filtered_dedup.bam; do
    filter_bam "${BAM}"
done

echo "Bound filtering complete. Ready for MEDIPS."
rm chrom.sizes
