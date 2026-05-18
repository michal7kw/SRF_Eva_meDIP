#!/bin/bash
#SBATCH --job-name=plot_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/15_plot_metaprofiles.err"
#SBATCH --output="logs/15_plot_metaprofiles.out"
#SBATCH --partition=workq

################################################################################
# Generate meDIP Metaprofiles (Average Signal Plots) Using deepTools
################################################################################

echo "========================================="
echo "Plot Metaprofiles - deepTools"
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
MATRIX_DIR="${BASE_DIR}/results/13_matrices"
OUT_DIR="${BASE_DIR}/results/15_metaprofiles"
TEMP_DIR="${OUT_DIR}/temp_matrices"

mkdir -p ${OUT_DIR}
mkdir -p ${TEMP_DIR}

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
echo "Generating Metaprofiles"
echo "========================================="
echo ""

################################################################################
# 1. Individual metaprofiles for each gene set (GFP + TES overlay)
################################################################################

echo "----------------------------------------"
echo "1. Individual Metaprofiles (GFP + TES)"
echo "----------------------------------------"
echo ""

for GENE_SET in "${GENE_SETS[@]}"; do
    echo "Processing: ${GENE_SET}"

    COMBINED_MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz"

    if [ ! -f "$COMBINED_MATRIX" ]; then
        echo "  WARNING: Matrix not found for ${GENE_SET}"
        echo "  Skipping..."
        echo ""
        continue
    fi

    # Generate metaprofile with both GFP and TES
    plotProfile \
        --matrixFile "$COMBINED_MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_metaprofile.pdf" \
        --outFileNameData "${OUT_DIR}/${GENE_SET}_metaprofile.tab" \
        --perGroup \
        --plotTitle "meDIP Metaprofile: ${GENE_SET//_/ }" \
        --refPointLabel "TSS" \
        --regionsLabel "${GENE_SET//_/ }" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --colors darkblue blue darkred red \
        --plotHeight 10 \
        --plotWidth 12 \
        --yAxisLabel "Mean meDIP Signal (RPKM)" \
        2>&1

    if [ -f "${OUT_DIR}/${GENE_SET}_metaprofile.pdf" ]; then
        echo "  ✓ Metaprofile saved"
    else
        echo "  ✗ ERROR generating metaprofile"
    fi
    echo ""
done

################################################################################
# 2. Combined metaprofiles: All gene sets on one plot (GFP condition)
################################################################################

echo "----------------------------------------"
echo "2. Combined Metaprofile: All Gene Sets (GFP)"
echo "----------------------------------------"
echo ""

# Collect all GFP matrices and merge using computeMatrixOperations rbind
GFP_MATRICES=()
LABELS=()

for GENE_SET in "${GENE_SETS[@]}"; do
    MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_matrix.gz"
    if [ -f "$MATRIX" ]; then
        GFP_MATRICES+=("$MATRIX")
        LABELS+=("${GENE_SET//_/ }")
    fi
done

if [ ${#GFP_MATRICES[@]} -gt 1 ]; then
    echo "  Merging ${#GFP_MATRICES[@]} matrices..."

    # Merge matrices using computeMatrixOperations rbind
    computeMatrixOperations rbind \
        --matrixFile "${GFP_MATRICES[@]}" \
        --outFileName "${TEMP_DIR}/all_gene_sets_GFP_merged.gz" \
        2>&1

    if [ -f "${TEMP_DIR}/all_gene_sets_GFP_merged.gz" ]; then
        plotProfile \
            --matrixFile "${TEMP_DIR}/all_gene_sets_GFP_merged.gz" \
            --outFileName "${OUT_DIR}/all_gene_sets_GFP_metaprofile.pdf" \
            --outFileNameData "${OUT_DIR}/all_gene_sets_GFP_metaprofile.tab" \
            --perGroup \
            --plotTitle "meDIP Metaprofile: All Gene Sets (GFP Control)" \
            --refPointLabel "TSS" \
            --colors '#e41a1c' '#377eb8' '#4daf4a' '#984ea3' '#ff7f00' '#a65628' \
            --plotHeight 10 \
            --plotWidth 14 \
            --yAxisLabel "Mean meDIP Signal (RPKM)" \
            --legendLocation upper-right \
            2>&1

        if [ -f "${OUT_DIR}/all_gene_sets_GFP_metaprofile.pdf" ]; then
            echo "  ✓ Combined GFP metaprofile saved"
        else
            echo "  ✗ ERROR generating combined GFP metaprofile"
        fi
    else
        echo "  ✗ ERROR merging GFP matrices"
    fi
else
    echo "  ✗ Not enough GFP matrices found (need at least 2)"
fi
echo ""

################################################################################
# 3. Combined metaprofiles: All gene sets on one plot (TES condition)
################################################################################

echo "----------------------------------------"
echo "3. Combined Metaprofile: All Gene Sets (TES)"
echo "----------------------------------------"
echo ""

# Collect all TES matrices
TES_MATRICES=()

for GENE_SET in "${GENE_SETS[@]}"; do
    MATRIX="${MATRIX_DIR}/${GENE_SET}_TES_matrix.gz"
    if [ -f "$MATRIX" ]; then
        TES_MATRICES+=("$MATRIX")
    fi
done

if [ ${#TES_MATRICES[@]} -gt 1 ]; then
    echo "  Merging ${#TES_MATRICES[@]} matrices..."

    # Merge matrices using computeMatrixOperations rbind
    computeMatrixOperations rbind \
        --matrixFile "${TES_MATRICES[@]}" \
        --outFileName "${TEMP_DIR}/all_gene_sets_TES_merged.gz" \
        2>&1

    if [ -f "${TEMP_DIR}/all_gene_sets_TES_merged.gz" ]; then
        plotProfile \
            --matrixFile "${TEMP_DIR}/all_gene_sets_TES_merged.gz" \
            --outFileName "${OUT_DIR}/all_gene_sets_TES_metaprofile.pdf" \
            --outFileNameData "${OUT_DIR}/all_gene_sets_TES_metaprofile.tab" \
            --perGroup \
            --plotTitle "meDIP Metaprofile: All Gene Sets (TES Treated)" \
            --refPointLabel "TSS" \
            --colors '#e41a1c' '#377eb8' '#4daf4a' '#984ea3' '#ff7f00' '#a65628' \
            --plotHeight 10 \
            --plotWidth 14 \
            --yAxisLabel "Mean meDIP Signal (RPKM)" \
            --legendLocation upper-right \
            2>&1

        if [ -f "${OUT_DIR}/all_gene_sets_TES_metaprofile.pdf" ]; then
            echo "  ✓ Combined TES metaprofile saved"
        else
            echo "  ✗ ERROR generating combined TES metaprofile"
        fi
    else
        echo "  ✗ ERROR merging TES matrices"
    fi
else
    echo "  ✗ Not enough TES matrices found (need at least 2)"
fi
echo ""

################################################################################
# 4. Expression-stratified metaprofiles (Highly vs Lowly expressed)
################################################################################

echo "----------------------------------------"
echo "4. Expression-Stratified Metaprofile"
echo "----------------------------------------"
echo ""

HIGH_MATRIX="${MATRIX_DIR}/highly_expressed_GFP_vs_TES_matrix.gz"
LOW_MATRIX="${MATRIX_DIR}/lowly_expressed_GFP_vs_TES_matrix.gz"

if [ -f "$HIGH_MATRIX" ] && [ -f "$LOW_MATRIX" ]; then
    echo "  Merging high and low expression matrices..."

    # Merge matrices using computeMatrixOperations rbind
    computeMatrixOperations rbind \
        --matrixFile "$HIGH_MATRIX" "$LOW_MATRIX" \
        --outFileName "${TEMP_DIR}/expression_stratified_merged.gz" \
        2>&1

    if [ -f "${TEMP_DIR}/expression_stratified_merged.gz" ]; then
        plotProfile \
            --matrixFile "${TEMP_DIR}/expression_stratified_merged.gz" \
            --outFileName "${OUT_DIR}/expression_stratified_metaprofile.pdf" \
            --outFileNameData "${OUT_DIR}/expression_stratified_metaprofile.tab" \
            --perGroup \
            --plotTitle "meDIP Metaprofile: High vs Low Expression" \
            --refPointLabel "TSS" \
            --plotHeight 10 \
            --plotWidth 14 \
            --yAxisLabel "Mean meDIP Signal (RPKM)" \
            --legendLocation upper-right \
            2>&1

        if [ -f "${OUT_DIR}/expression_stratified_metaprofile.pdf" ]; then
            echo "  ✓ Expression-stratified metaprofile saved"
        else
            echo "  ✗ ERROR generating expression-stratified metaprofile"
        fi
    else
        echo "  ✗ ERROR merging expression matrices"
    fi
else
    echo "  ✗ Required matrices not found"
    [ ! -f "$HIGH_MATRIX" ] && echo "    Missing: $HIGH_MATRIX"
    [ ! -f "$LOW_MATRIX" ] && echo "    Missing: $LOW_MATRIX"
fi
echo ""

################################################################################
# 5. Regulation-stratified metaprofiles (Up vs Down vs Unchanged)
################################################################################

echo "----------------------------------------"
echo "5. Regulation-Stratified Metaprofile"
echo "----------------------------------------"
echo ""

UP_MATRIX="${MATRIX_DIR}/upregulated_GFP_vs_TES_matrix.gz"
DOWN_MATRIX="${MATRIX_DIR}/downregulated_GFP_vs_TES_matrix.gz"
UNCHANGED_MATRIX="${MATRIX_DIR}/unchanged_GFP_vs_TES_matrix.gz"

if [ -f "$UP_MATRIX" ] && [ -f "$DOWN_MATRIX" ] && [ -f "$UNCHANGED_MATRIX" ]; then
    echo "  Merging regulation category matrices..."

    # Merge matrices using computeMatrixOperations rbind
    computeMatrixOperations rbind \
        --matrixFile "$UP_MATRIX" "$DOWN_MATRIX" "$UNCHANGED_MATRIX" \
        --outFileName "${TEMP_DIR}/regulation_stratified_merged.gz" \
        2>&1

    if [ -f "${TEMP_DIR}/regulation_stratified_merged.gz" ]; then
        plotProfile \
            --matrixFile "${TEMP_DIR}/regulation_stratified_merged.gz" \
            --outFileName "${OUT_DIR}/regulation_stratified_metaprofile.pdf" \
            --outFileNameData "${OUT_DIR}/regulation_stratified_metaprofile.tab" \
            --perGroup \
            --plotTitle "meDIP Metaprofile: Gene Regulation Categories" \
            --refPointLabel "TSS" \
            --plotHeight 10 \
            --plotWidth 14 \
            --yAxisLabel "Mean meDIP Signal (RPKM)" \
            --legendLocation upper-right \
            2>&1

        if [ -f "${OUT_DIR}/regulation_stratified_metaprofile.pdf" ]; then
            echo "  ✓ Regulation-stratified metaprofile saved"
        else
            echo "  ✗ ERROR generating regulation-stratified metaprofile"
        fi
    else
        echo "  ✗ ERROR merging regulation matrices"
    fi
else
    echo "  ✗ Required matrices not found"
    [ ! -f "$UP_MATRIX" ] && echo "    Missing: $UP_MATRIX"
    [ ! -f "$DOWN_MATRIX" ] && echo "    Missing: $DOWN_MATRIX"
    [ ! -f "$UNCHANGED_MATRIX" ] && echo "    Missing: $UNCHANGED_MATRIX"
fi
echo ""

# NOTE: TES vs TESmut comparison removed (TESmut samples excluded - failed sample)

# Clean up temp directory
rm -rf ${TEMP_DIR}

################################################################################
# Summary
################################################################################

echo "========================================="
echo "Metaprofile Generation Complete"
echo "========================================="
echo ""

echo "Generated metaprofiles in: ${OUT_DIR}/"
echo ""

# Count output files
NUM_PROFILES=$(ls -1 ${OUT_DIR}/*.pdf 2>/dev/null | wc -l)
echo "Total metaprofile files created: $NUM_PROFILES"
echo ""

# List all metaprofiles
echo "Metaprofile files:"
ls -lh ${OUT_DIR}/*.pdf
echo ""

echo "Metaprofile types generated:"
echo "  1. Individual gene set metaprofiles (GFP + TES overlay)"
echo "  2. Combined all gene sets (GFP condition) - merged matrix"
echo "  3. Combined all gene sets (TES condition) - merged matrix"
echo "  4. Expression-stratified (high vs low) - merged matrix"
echo "  5. Regulation-stratified (up vs down vs unchanged) - merged matrix"
echo "  (Note: TESmut excluded from analysis - failed sample)"
echo ""

echo "Data tables (tab-delimited):"
ls -lh ${OUT_DIR}/*.tab 2>/dev/null || echo "  No tab files found"
echo ""

echo "Expected output files:"
echo "  Individual: {all_genes,highly_expressed,lowly_expressed,upregulated,downregulated,unchanged}_metaprofile.pdf"
echo "  Combined:   all_gene_sets_{GFP,TES}_metaprofile.pdf"
echo "  Stratified: expression_stratified_metaprofile.pdf, regulation_stratified_metaprofile.pdf"
echo ""

echo "Next steps:"
echo "  1. Review metaprofiles to understand global methylation patterns"
echo "  2. Run 16_advanced_visualization.R for integrated analysis"
echo "  3. Compare heatmaps with metaprofiles for consistency"
echo ""

echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
