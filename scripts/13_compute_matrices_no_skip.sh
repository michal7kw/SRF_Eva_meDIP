#!/bin/bash
#SBATCH --job-name=compute_matrices
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/13_compute_matrices.err"
#SBATCH --output="logs/13_compute_matrices.out"
#SBATCH --partition=workq

################################################################################
# Compute meDIP Signal Matrices at Promoters Using deepTools
################################################################################

echo "========================================="
echo "Compute Signal Matrices - deepTools"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Cores: $SLURM_NTASKS"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new  # Has deepTools installed

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
BIGWIG_DIR="${BASE_DIR}/results/05_bigwig"
BED_DIR="${BASE_DIR}/results/12_gene_sets"
OUT_DIR="${BASE_DIR}/results/13_matrices"

mkdir -p ${OUT_DIR}

# Define BigWig files
# NOTE: TESmut samples excluded from analysis (failed sample)

# GFP samples (control)
GFP_BIGWIGS=(
    "${BIGWIG_DIR}/GFP-1-IP_RPKM.bw"
    "${BIGWIG_DIR}/GFP-2-IP_RPKM.bw"
)

# TES samples (treated)
TES_BIGWIGS=(
    "${BIGWIG_DIR}/TES-1-IP_RPKM.bw"
    "${BIGWIG_DIR}/TES-2-IP_RPKM.bw"
)

# Check if BigWig files exist
echo "Checking BigWig files..."
for bw in "${GFP_BIGWIGS[@]}" "${TES_BIGWIGS[@]}"; do
    if [ ! -f "$bw" ]; then
        echo "  WARNING: Missing $bw"
    else
        echo "  Found: $(basename $bw)"
    fi
done
echo ""

# Define gene sets
GENE_SETS=(
    "all_genes"
    "highly_expressed"
    "lowly_expressed"
    "upregulated"
    "downregulated"
    "unchanged"
)

echo "========================================="
echo "Computing Matrices for Each Gene Set"
echo "========================================="
echo ""

# Parameters
UPSTREAM=5000      # 5kb upstream of TSS
DOWNSTREAM=5000    # 5kb downstream of TSS
BINSIZE=50         # 50bp bins
CORES=$SLURM_NTASKS

################################################################################
# Loop through gene sets and compute matrices
################################################################################

for GENE_SET in "${GENE_SETS[@]}"; do
    echo "----------------------------------------"
    echo "Processing: ${GENE_SET}"
    echo "----------------------------------------"

    BED_FILE="${BED_DIR}/${GENE_SET}_promoters.bed"

    # Check if BED file exists
    if [ ! -f "$BED_FILE" ]; then
        echo "  WARNING: BED file not found: $BED_FILE"
        echo "  Skipping ${GENE_SET}"
        echo ""
        continue
    fi

    # Count regions
    NUM_REGIONS=$(wc -l < "$BED_FILE")
    echo "  Regions: $NUM_REGIONS"
    echo ""

    ########################################################################
    # 1. GFP (Control) Matrix
    ########################################################################
    echo "  [1/2] Computing GFP matrix..."

    # NOTE: Using 'TSS' because BED files now contain just the TSS point (1bp),
    # so we center on the reference point (TSS)
    computeMatrix reference-point \
        --referencePoint TSS \
        --scoreFileName "${GFP_BIGWIGS[@]}" \
        --regionsFileName "$BED_FILE" \
        --outFileName "${OUT_DIR}/${GENE_SET}_GFP_matrix.gz" \
        --outFileNameMatrix "${OUT_DIR}/${GENE_SET}_GFP_matrix.tab" \
        --upstream $UPSTREAM \
        --downstream $DOWNSTREAM \
        --binSize $BINSIZE \
        \
        --numberOfProcessors $CORES \
        2>&1 | grep -v "^$"

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "  ✓ GFP matrix saved"
    else
        echo "  ✗ ERROR computing GFP matrix"
    fi
    echo ""

    ########################################################################
    # 2. TES (Treated) Matrix
    ########################################################################
    echo "  [2/2] Computing TES matrix..."

    computeMatrix reference-point \
        --referencePoint TSS \
        --scoreFileName "${TES_BIGWIGS[@]}" \
        --regionsFileName "$BED_FILE" \
        --outFileName "${OUT_DIR}/${GENE_SET}_TES_matrix.gz" \
        --outFileNameMatrix "${OUT_DIR}/${GENE_SET}_TES_matrix.tab" \
        --upstream $UPSTREAM \
        --downstream $DOWNSTREAM \
        --binSize $BINSIZE \
        \
        --numberOfProcessors $CORES \
        2>&1 | grep -v "^$"

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "  ✓ TES matrix saved"
    else
        echo "  ✗ ERROR computing TES matrix"
    fi
    echo ""

    # NOTE: TESmut matrix computation removed (TESmut samples excluded - failed sample)

    ########################################################################
    # 3. Combined Matrix (GFP + TES for side-by-side comparison)
    ########################################################################
    echo "  [3/3] Computing combined GFP+TES matrix..."

    computeMatrix reference-point \
        --referencePoint TSS \
        --scoreFileName "${GFP_BIGWIGS[@]}" "${TES_BIGWIGS[@]}" \
        --regionsFileName "$BED_FILE" \
        --outFileName "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz" \
        --outFileNameMatrix "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_matrix.tab" \
        --upstream $UPSTREAM \
        --downstream $DOWNSTREAM \
        --binSize $BINSIZE \
        \
        --numberOfProcessors $CORES \
        2>&1 | grep -v "^$"

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "  ✓ Combined matrix saved"
    else
        echo "  ✗ ERROR computing combined matrix"
    fi
    echo ""

done

################################################################################
# Summary
################################################################################

echo "========================================="
echo "Matrix Computation Complete"
echo "========================================="
echo ""

echo "Generated matrices in: ${OUT_DIR}/"
echo ""

# Count output files
NUM_MATRICES=$(ls -1 ${OUT_DIR}/*.gz 2>/dev/null | wc -l)
echo "Total matrix files created: $NUM_MATRICES"
echo ""

# List all matrices
echo "Matrix files:"
ls -lh ${OUT_DIR}/*.gz
echo ""

echo "Next steps:"
echo "  1. Run 14_plot_heatmaps.sh to generate heatmaps"
echo "  2. Run 15_plot_metaprofiles.sh to generate metaprofiles"
echo ""

echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
