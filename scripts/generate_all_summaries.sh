#!/bin/bash
#SBATCH --job-name=generate_summaries
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/generate_summaries_%j.out
#SBATCH --error=logs/generate_summaries_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

################################################################################
# Master Script: Generate All Summary Statistics
#
# Runs both meDIP and integration summary scripts to produce:
# 1. meDIP summary statistics (DMRs, correlations)
# 2. Integration summary statistics (binding + methylation + expression)
#
# Usage:
#   sbatch scripts/generate_all_summaries.sh
#   OR
#   bash scripts/generate_all_summaries.sh
################################################################################

set -e  # Exit on error

# Define paths
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
MEDIP_DIR="${BASE_DIR}/meDIP"
INTEGRATION_DIR="${BASE_DIR}/SRF_Eva_integrated_analysis"

# Create log directory if running interactively
mkdir -p "${BASE_DIR}/logs"

echo "========================================"
echo "Multi-Omics Summary Generation Pipeline"
echo "========================================"
echo "Start time: $(date)"
echo "Working directory: ${BASE_DIR}"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

################################################################################
# Step 1: Generate meDIP Summary Statistics
################################################################################

echo "========================================"
echo "Step 1: meDIP Summary Statistics"
echo "========================================"

cd "${MEDIP_DIR}"

if [ -f "scripts/summary_medip_statistics.R" ]; then
    echo "Running: summary_medip_statistics.R"
    Rscript scripts/summary_medip_statistics.R

    if [ $? -eq 0 ]; then
        echo "SUCCESS: meDIP summary generated"
    else
        echo "ERROR: meDIP summary failed"
        exit 1
    fi
else
    echo "ERROR: summary_medip_statistics.R not found"
    exit 1
fi

echo ""

################################################################################
# Step 2: Generate Integration Summary Statistics
################################################################################

echo "========================================"
echo "Step 2: Integration Summary Statistics"
echo "========================================"

cd "${INTEGRATION_DIR}"

if [ -f "scripts/summary_integration_statistics.R" ]; then
    echo "Running: summary_integration_statistics.R"
    Rscript scripts/summary_integration_statistics.R

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Integration summary generated"
    else
        echo "ERROR: Integration summary failed"
        exit 1
    fi
else
    echo "ERROR: summary_integration_statistics.R not found"
    exit 1
fi

echo ""

################################################################################
# Step 3: Generate Markdown Summary
################################################################################

echo "========================================"
echo "Step 3: Generate Markdown Summary"
echo "========================================"

cd "${BASE_DIR}"

if [ -f "scripts/generate_markdown_summary.R" ]; then
    echo "Running: generate_markdown_summary.R"
    Rscript scripts/generate_markdown_summary.R

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Markdown summary generated"
    else
        echo "ERROR: Markdown summary failed"
        exit 1
    fi
else
    echo "WARNING: generate_markdown_summary.R not found, skipping"
fi

echo ""

################################################################################
# Step 4: Combine All Summaries
################################################################################

echo "========================================"
echo "Step 4: Combining Summaries"
echo "========================================"

cd "${BASE_DIR}"

# Create combined output directory
COMBINED_DIR="${BASE_DIR}/results_summary"
mkdir -p "${COMBINED_DIR}"

# Copy all summary files to combined directory
echo "Copying summary files to: ${COMBINED_DIR}"

# meDIP summaries
if [ -d "${MEDIP_DIR}/results/summary" ]; then
    cp -v "${MEDIP_DIR}/results/summary/"*.csv "${COMBINED_DIR}/" 2>/dev/null || true
fi

# Integration summaries
if [ -d "${INTEGRATION_DIR}/results/summary" ]; then
    cp -v "${INTEGRATION_DIR}/results/summary/"*.csv "${COMBINED_DIR}/" 2>/dev/null || true
fi

echo ""

################################################################################
# Step 5: Generate Text Report
################################################################################

echo "========================================"
echo "Step 5: Generating Text Report"
echo "========================================"

# Create a simple combined text report
REPORT_FILE="${COMBINED_DIR}/COMBINED_SUMMARY_REPORT.txt"

cat > "${REPORT_FILE}" << 'EOF'
================================================================================
MULTI-OMICS ANALYSIS SUMMARY REPORT
TES vs GFP Comparison (SNB19 Glioblastoma Cells)
================================================================================

Generated: DATE_PLACEHOLDER

This report summarizes the key findings from the integrated multi-omics analysis
combining Cut&Tag (binding), meDIP-seq (methylation), and RNA-seq (expression).

================================================================================
FILES GENERATED
================================================================================

meDIP Statistics:
  - medip_summary_statistics.csv      : All meDIP statistics (long format)
  - medip_summary_wide.csv            : Statistics in wide format
  - medip_top_dmrs_by_category.csv    : Top 50 DMRs per category

Integration Statistics:
  - integration_summary_statistics.csv : All integration statistics
  - integration_summary_wide.csv       : Statistics in wide format
  - integration_top_genes_by_category.csv : Top 50 genes per regulatory class

================================================================================
KEY FINDINGS (populated from CSV files)
================================================================================

See individual CSV files for complete statistics.

Key categories included:
  - DMR_Overview: Total DMRs, direction breakdown
  - DMR_Effect_Size: Fold change statistics
  - DMR_Genomic_Distribution: Promoter vs gene body vs intergenic
  - Direct_Targets: Binding + expression classification
  - Meth_Expr_Correlation: Pearson and Spearman correlations
  - DMR_DEG_Enrichment: Fisher's exact test for DMR enrichment at DEGs

================================================================================
METHODOLOGY
================================================================================

1. DIFFERENTIAL METHYLATION (meDIP-seq)
   - Method: MEDIPS with CpG density normalization
   - Windows: 500bp genome-wide
   - Statistics: edgeR negative binomial GLM
   - Thresholds: FDR < 0.05, |log2FC| > 1 (stringent)

2. DIFFERENTIAL EXPRESSION (RNA-seq)
   - Method: DESeq2
   - Thresholds: padj < 0.05, |log2FC| > 1

3. BINDING ANALYSIS (Cut&Tag)
   - Peaks: MACS2 narrow peaks
   - Annotation: ChIPseeker
   - Distance threshold: 50kb for enhancer assignment

4. INTEGRATION
   - Peak-to-gene mapping via promoter/enhancer annotation
   - Gene ID harmonization: Ensembl <-> Entrez <-> Symbol
   - DMR-gene mapping via ChIPseeker annotation

================================================================================
EOF

# Replace date placeholder
sed -i "s/DATE_PLACEHOLDER/$(date)/" "${REPORT_FILE}"

echo "Created: ${REPORT_FILE}"

################################################################################
# Summary
################################################################################

echo ""
echo "========================================"
echo "SUMMARY GENERATION COMPLETE"
echo "========================================"
echo ""
echo "Output files in: ${COMBINED_DIR}/"
echo ""
ls -lh "${COMBINED_DIR}/"
echo ""
echo "End time: $(date)"
echo "========================================"
