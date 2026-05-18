#!/bin/bash
#SBATCH --job-name=23_compute_binding_matrices
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/23_compute_binding_matrices.err"
#SBATCH --output="logs/23_compute_binding_matrices.out"
#SBATCH --partition=workq

################################################################################
# 23_compute_binding_matrices.sh
################################################################################
# Purpose: Compute deepTools matrices for binding-stratified gene sets
#
# Gene Sets:
#   1. direct_bound_downregulated - DEGs down + TES/TEAD binding
#   2. indirect_unbound_downregulated - DEGs down + NO binding
#   3. random_unbound_unchanged - Control genes (unchanged + NO binding)
#
# Output: Matrices for metaprofile visualization
################################################################################

echo "========================================="
echo "Compute Binding-Stratified Matrices"
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
BED_DIR="${BASE_DIR}/results/23_binding_stratified_metaprofiles/beds"
OUT_DIR="${BASE_DIR}/results/23_binding_stratified_metaprofiles/matrices"

mkdir -p ${OUT_DIR}

# Define BigWig files
GFP_BIGWIGS=(
    "${BIGWIG_DIR}/GFP-1-IP_RPKM.bw"
    "${BIGWIG_DIR}/GFP-2-IP_RPKM.bw"
)

TES_BIGWIGS=(
    "${BIGWIG_DIR}/TES-1-IP_RPKM.bw"
    "${BIGWIG_DIR}/TES-2-IP_RPKM.bw"
)

# Check if BigWig files exist
echo "Checking BigWig files..."
for bw in "${GFP_BIGWIGS[@]}" "${TES_BIGWIGS[@]}"; do
    if [ ! -f "$bw" ]; then
        echo "  ERROR: Missing $bw"
        exit 1
    else
        echo "  Found: $(basename $bw)"
    fi
done
echo ""

# Define gene sets for binding-stratified analysis
GENE_SETS=(
    "direct_bound_downregulated"
    "indirect_unbound_downregulated"
    "random_unbound_unchanged"
)

# Check if BED files exist
echo "Checking BED files..."
for GENE_SET in "${GENE_SETS[@]}"; do
    BED_FILE="${BED_DIR}/${GENE_SET}_promoter.bed"
    if [ ! -f "$BED_FILE" ]; then
        echo "  ERROR: Missing $BED_FILE"
        echo "  Please run 23_prepare_binding_stratified_beds.sh first"
        exit 1
    else
        NUM_REGIONS=$(wc -l < "$BED_FILE")
        echo "  Found: ${GENE_SET}_promoter.bed ($NUM_REGIONS regions)"
    fi
done
echo ""

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
# Loop through gene sets and compute combined matrices (GFP + TES)
################################################################################

for GENE_SET in "${GENE_SETS[@]}"; do
    echo "----------------------------------------"
    echo "Processing: ${GENE_SET}"
    echo "----------------------------------------"

    BED_FILE="${BED_DIR}/${GENE_SET}_promoter.bed"
    NUM_REGIONS=$(wc -l < "$BED_FILE")
    echo "  Regions: $NUM_REGIONS"
    echo ""

    ########################################################################
    # Compute combined GFP + TES matrix for overlay visualization
    ########################################################################
    echo "  Computing GFP vs TES combined matrix..."

    computeMatrix reference-point \
        --referencePoint TSS \
        --scoreFileName "${GFP_BIGWIGS[@]}" "${TES_BIGWIGS[@]}" \
        --regionsFileName "$BED_FILE" \
        --outFileName "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz" \
        --outFileNameMatrix "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_matrix.tab" \
        --upstream $UPSTREAM \
        --downstream $DOWNSTREAM \
        --binSize $BINSIZE \
        --skipZeros \
        --numberOfProcessors $CORES

    if [ $? -eq 0 ]; then
        echo "  [OK] Matrix saved: ${OUT_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz"
    else
        echo "  [ERROR] Failed to compute matrix for ${GENE_SET}"
        exit 1
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

# List all matrices with sizes
echo "Matrix files:"
ls -lh ${OUT_DIR}/*.gz
echo ""

echo "Next step:"
echo "  Run 23_plot_binding_metaprofiles.sh to generate metaprofiles"
echo ""

echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
