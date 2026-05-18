#!/bin/bash
#SBATCH --job-name=23_binding_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/23_plot_binding_metaprofiles_%j.err"
#SBATCH --output="logs/23_plot_binding_metaprofiles_%j.out"
#SBATCH --partition=workq

################################################################################
# 23_plot_binding_metaprofiles.sh
################################################################################
# Purpose: Generate metaprofiles comparing MeDIP signals for binding-stratified
#          gene sets. Each plot shows GFP vs TES overlay to assess differential
#          methylation between conditions.
#
# Gene Sets:
#   1. Direct targets (DEGs down + TES/TEAD binding) - expect differential MeDIP
#   2. Indirect targets (DEGs down + NO binding) - expect similar MeDIP
#   3. Random control (unchanged + NO binding) - baseline control
#
# Output: PDF metaprofiles for each category
################################################################################

echo "========================================="
echo "Plot Binding-Stratified Metaprofiles"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new  # Has deepTools installed

# Define directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
MATRIX_DIR="${BASE_DIR}/results/23_binding_stratified_metaprofiles/matrices"
OUT_DIR="${BASE_DIR}/results/23_binding_stratified_metaprofiles/plots"
TEMP_DIR="${OUT_DIR}/temp_matrices"

mkdir -p ${OUT_DIR}
mkdir -p ${TEMP_DIR}

# Define gene sets with descriptive titles
declare -A TITLES
TITLES["direct_bound_downregulated"]="DEGs Down WITH TES/TEAD1 Binding (Direct Targets)"
TITLES["indirect_unbound_downregulated"]="DEGs Down WITHOUT TES/TEAD1 Binding (Indirect)"
TITLES["random_unbound_unchanged"]="Unchanged Genes WITHOUT Binding (Random Control)"

# Check if matrices exist
echo "Checking input matrices..."
for GENE_SET in "direct_bound_downregulated" "indirect_unbound_downregulated" "random_unbound_unchanged"; do
    MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz"
    if [ ! -f "$MATRIX" ]; then
        echo "  ERROR: Missing $MATRIX"
        echo "  Please run 23_compute_binding_matrices.sh first"
        exit 1
    else
        echo "  Found: $(basename $MATRIX)"
    fi
done
echo ""

echo "========================================="
echo "1. Individual Metaprofiles (GFP vs TES)"
echo "========================================="
echo ""

################################################################################
# Generate individual metaprofiles for each category
################################################################################

for GENE_SET in "direct_bound_downregulated" "indirect_unbound_downregulated" "random_unbound_unchanged"; do
    echo "Processing: ${GENE_SET}"
    echo "  Title: ${TITLES[$GENE_SET]}"

    MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz"
    TITLE="${TITLES[$GENE_SET]}"

    # Generate metaprofile with GFP and TES overlay
    plotProfile \
        --matrixFile "$MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_metaprofile.pdf" \
        --outFileNameData "${OUT_DIR}/${GENE_SET}_metaprofile.tab" \
        --perGroup \
        --plotTitle "$TITLE" \
        --refPointLabel "TSS" \
        --regionsLabel "Genes" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --colors darkblue blue darkred red \
        --plotHeight 10 \
        --plotWidth 12 \
        --yAxisLabel "Mean meDIP Signal (RPKM)" \
        --legendLocation upper-right \
        2>&1

    if [ -f "${OUT_DIR}/${GENE_SET}_metaprofile.pdf" ]; then
        echo "  [OK] Metaprofile saved: ${OUT_DIR}/${GENE_SET}_metaprofile.pdf"
    else
        echo "  [ERROR] Failed to generate metaprofile for ${GENE_SET}"
    fi
    echo ""
done

echo "========================================="
echo "2. Combined Comparative Metaprofile"
echo "========================================="
echo ""

################################################################################
# Merge all matrices to create a combined comparison plot
################################################################################

# Collect matrices in order
MATRICES=(
    "${MATRIX_DIR}/direct_bound_downregulated_GFP_vs_TES_matrix.gz"
    "${MATRIX_DIR}/indirect_unbound_downregulated_GFP_vs_TES_matrix.gz"
    "${MATRIX_DIR}/random_unbound_unchanged_GFP_vs_TES_matrix.gz"
)

echo "  Merging 3 matrices for combined plot..."

# Merge matrices using computeMatrixOperations rbind
computeMatrixOperations rbind \
    --matrixFile "${MATRICES[@]}" \
    --outFileName "${TEMP_DIR}/all_binding_categories_merged.gz" \
    2>&1

if [ -f "${TEMP_DIR}/all_binding_categories_merged.gz" ]; then
    plotProfile \
        --matrixFile "${TEMP_DIR}/all_binding_categories_merged.gz" \
        --outFileName "${OUT_DIR}/binding_comparison_metaprofile.pdf" \
        --outFileNameData "${OUT_DIR}/binding_comparison_metaprofile.tab" \
        --perGroup \
        --plotTitle "meDIP Metaprofile: Binding Status Comparison" \
        --refPointLabel "TSS" \
        --plotHeight 12 \
        --plotWidth 14 \
        --yAxisLabel "Mean meDIP Signal (RPKM)" \
        --legendLocation upper-right \
        2>&1

    if [ -f "${OUT_DIR}/binding_comparison_metaprofile.pdf" ]; then
        echo "  [OK] Combined metaprofile saved: ${OUT_DIR}/binding_comparison_metaprofile.pdf"
    else
        echo "  [ERROR] Failed to generate combined metaprofile"
    fi
else
    echo "  [ERROR] Failed to merge matrices"
fi
echo ""

################################################################################
# Generate heatmaps (additional visualization)
################################################################################

echo "========================================="
echo "3. Heatmaps for Each Category"
echo "========================================="
echo ""

for GENE_SET in "direct_bound_downregulated" "indirect_unbound_downregulated" "random_unbound_unchanged"; do
    echo "Processing: ${GENE_SET}"

    MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz"
    TITLE="${TITLES[$GENE_SET]}"

    # Generate heatmap
    plotHeatmap \
        --matrixFile "$MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_heatmap.pdf" \
        --sortUsing mean \
        --sortUsingSamples 1 2 \
        --averageTypeSummaryPlot mean \
        --colorMap YlOrRd \
        --plotTitle "$TITLE" \
        --refPointLabel "TSS" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --heatmapHeight 15 \
        --heatmapWidth 8 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        2>&1

    if [ -f "${OUT_DIR}/${GENE_SET}_heatmap.pdf" ]; then
        echo "  [OK] Heatmap saved: ${OUT_DIR}/${GENE_SET}_heatmap.pdf"
    else
        echo "  [ERROR] Failed to generate heatmap for ${GENE_SET}"
    fi
    echo ""
done

# Clean up temp directory
rm -rf ${TEMP_DIR}

################################################################################
# Summary
################################################################################

echo "========================================="
echo "Metaprofile Generation Complete"
echo "========================================="
echo ""

echo "Output directory: ${OUT_DIR}/"
echo ""

echo "Generated files:"
ls -lh ${OUT_DIR}/*.pdf 2>/dev/null
echo ""

echo "Data tables (for quantitative analysis):"
ls -lh ${OUT_DIR}/*.tab 2>/dev/null
echo ""

echo "Interpretation guide:"
echo "  - Direct targets: Expect DIFFERENTIAL methylation (GFP ≠ TES)"
echo "  - Indirect targets: Expect SIMILAR methylation (GFP ≈ TES)"
echo "  - Random control: Baseline - no expected difference"
echo ""

echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
