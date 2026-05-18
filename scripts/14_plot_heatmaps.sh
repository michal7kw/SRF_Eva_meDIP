#!/bin/bash
#SBATCH --job-name=plot_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/14_plot_heatmaps.err"
#SBATCH --output="logs/14_plot_heatmaps.out"
#SBATCH --partition=workq

################################################################################
# Generate meDIP Heatmaps Using deepTools
################################################################################

echo "========================================="
echo "Plot Heatmaps - deepTools"
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
OUT_DIR="${BASE_DIR}/results/14_heatmaps"

mkdir -p ${OUT_DIR}

# Define gene sets
GENE_SETS=(
    "all_genes"
    "highly_expressed"
    "lowly_expressed"
    "upregulated"
    "downregulated"
    "unchanged"
)

# Color schemes
COLOR_MAP="RdYlBu_r"  # Red (high) - Yellow - Blue (low), reversed

echo "========================================="
echo "Generating Heatmaps for Each Gene Set"
echo "========================================="
echo ""

################################################################################
# Loop through gene sets and generate heatmaps
################################################################################

for GENE_SET in "${GENE_SETS[@]}"; do
    echo "----------------------------------------"
    echo "Processing: ${GENE_SET}"
    echo "----------------------------------------"

    # Check if matrices exist
    GFP_MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_matrix.gz"
    TES_MATRIX="${MATRIX_DIR}/${GENE_SET}_TES_matrix.gz"
    COMBINED_MATRIX="${MATRIX_DIR}/${GENE_SET}_GFP_vs_TES_matrix.gz"

    if [ ! -f "$GFP_MATRIX" ] || [ ! -f "$TES_MATRIX" ]; then
        echo "  WARNING: Matrix files not found for ${GENE_SET}"
        echo "  Skipping..."
        echo ""
        continue
    fi

    ########################################################################
    # 1. Side-by-side heatmap: GFP | TES (Combined Matrix)
    ########################################################################
    echo "  [1/5] Generating side-by-side heatmap (GFP | TES)..."

    plotHeatmap \
        --matrixFile "$COMBINED_MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_sidebyside.pdf" \
        --outFileSortedRegions "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_regions.bed" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap $COLOR_MAP \
        --zMin 0 \
        --whatToShow 'heatmap and colorbar' \
        --heatmapHeight 20 \
        --heatmapWidth 6 \
        --xAxisLabel "Distance from TSS (bp)" \
        --refPointLabel "TSS" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --plotTitle "meDIP Signal: ${GENE_SET//_/ } (GFP vs TES)" \
        --legendLocation upper-left \
        2>&1 | grep -v "^$"

    if [ $? -eq 0 ]; then
        echo "  ✓ Side-by-side heatmap saved"
    else
        echo "  ✗ ERROR generating side-by-side heatmap"
    fi
    echo ""

    ########################################################################
    # 2. GFP heatmap only (sorted by mean signal)
    ########################################################################
    echo "  [2/5] Generating GFP heatmap..."

    plotHeatmap \
        --matrixFile "$GFP_MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_GFP_heatmap.pdf" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap $COLOR_MAP \
        --zMin 0 \
        --whatToShow 'heatmap and colorbar' \
        --heatmapHeight 20 \
        --heatmapWidth 4 \
        --xAxisLabel "Distance from TSS (bp)" \
        --refPointLabel "TSS" \
        --samplesLabel "GFP-1" "GFP-2" \
        --plotTitle "meDIP Signal (GFP): ${GENE_SET//_/ }" \
        --legendLocation upper-left \
        2>&1 | grep -v "^$"

    if [ $? -eq 0 ]; then
        echo "  ✓ GFP heatmap saved"
    else
        echo "  ✗ ERROR generating GFP heatmap"
    fi
    echo ""

    ########################################################################
    # 3. TES heatmap only (sorted by mean signal)
    ########################################################################
    echo "  [3/5] Generating TES heatmap..."

    plotHeatmap \
        --matrixFile "$TES_MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_TES_heatmap.pdf" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap $COLOR_MAP \
        --zMin 0 \
        --whatToShow 'heatmap and colorbar' \
        --heatmapHeight 20 \
        --heatmapWidth 4 \
        --xAxisLabel "Distance from TSS (bp)" \
        --refPointLabel "TSS" \
        --samplesLabel "TES-1" "TES-2" \
        --plotTitle "meDIP Signal (TES): ${GENE_SET//_/ }" \
        --legendLocation upper-left \
        2>&1 | grep -v "^$"

    if [ $? -eq 0 ]; then
        echo "  ✓ TES heatmap saved"
    else
        echo "  ✗ ERROR generating TES heatmap"
    fi
    echo ""

    ########################################################################
    # 4. K-means clustered heatmap (GFP + TES combined, 4 clusters)
    ########################################################################
    echo "  [4/5] Generating k-means clustered heatmap..."

    plotHeatmap \
        --matrixFile "$COMBINED_MATRIX" \
        --outFileName "${OUT_DIR}/${GENE_SET}_GFP_vs_TES_clustered.pdf" \
        --kmeans 4 \
        --colorMap $COLOR_MAP \
        --zMin 0 \
        --whatToShow 'heatmap and colorbar' \
        --heatmapHeight 20 \
        --heatmapWidth 6 \
        --xAxisLabel "Distance from TSS (bp)" \
        --refPointLabel "TSS" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --plotTitle "meDIP Signal (K-means): ${GENE_SET//_/ }" \
        --legendLocation upper-left \
        2>&1 | grep -v "^$"

    if [ $? -eq 0 ]; then
        echo "  ✓ Clustered heatmap saved"
    else
        echo "  ✗ ERROR generating clustered heatmap"
    fi
    echo ""

    # NOTE: TESmut heatmap generation removed (TESmut samples excluded - failed sample)
    echo ""

done

################################################################################
# Generate comparative heatmaps across gene sets
################################################################################

echo "========================================="
echo "Generating Comparative Visualizations"
echo "========================================="
echo ""

# For upregulated vs downregulated comparison
UP_COMBINED="${MATRIX_DIR}/upregulated_GFP_vs_TES_matrix.gz"
DOWN_COMBINED="${MATRIX_DIR}/downregulated_GFP_vs_TES_matrix.gz"

if [ -f "$UP_COMBINED" ] && [ -f "$DOWN_COMBINED" ]; then
    echo "Creating upregulated vs downregulated comparison..."

    plotHeatmap \
        --matrixFile "$UP_COMBINED" "$DOWN_COMBINED" \
        --outFileName "${OUT_DIR}/upregulated_vs_downregulated_comparison.pdf" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap $COLOR_MAP \
        --zMin 0 \
        --whatToShow 'heatmap and colorbar' \
        --heatmapHeight 20 \
        --heatmapWidth 6 \
        --xAxisLabel "Distance from TSS (bp)" \
        --refPointLabel "TSS" \
        --regionsLabel "Upregulated genes" "Downregulated genes" \
        --plotTitle "meDIP Signal: DEGs Comparison" \
        --legendLocation upper-left \
        2>&1 | grep -v "^$"

    echo "  ✓ Comparative heatmap saved"
    echo ""
fi

################################################################################
# Summary
################################################################################

echo "========================================="
echo "Heatmap Generation Complete"
echo "========================================="
echo ""

echo "Generated heatmaps in: ${OUT_DIR}/"
echo ""

# Count output files
NUM_HEATMAPS=$(ls -1 ${OUT_DIR}/*.pdf 2>/dev/null | wc -l)
echo "Total heatmap files created: $NUM_HEATMAPS"
echo ""

# List all heatmaps
echo "Heatmap files:"
ls -lh ${OUT_DIR}/*.pdf
echo ""

echo "Heatmap types generated for each gene set:"
echo "  1. Side-by-side (GFP | TES) - for direct comparison"
echo "  2. GFP only - control condition"
echo "  3. TES only - treated condition"
echo "  4. K-means clustered - identifies methylation patterns"
echo "  (Note: TESmut excluded from analysis - failed sample)"
echo ""

echo "Next steps:"
echo "  1. Review heatmaps to identify methylation patterns"
echo "  2. Run 15_plot_metaprofiles.sh to generate average signal plots"
echo "  3. Run 16_advanced_visualization.R for integrated analysis"
echo ""

echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
