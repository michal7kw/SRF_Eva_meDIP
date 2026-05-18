#!/bin/bash

################################################################################
# Script: master_script.sh
# Purpose: Master pipeline orchestration for meDIP-seq analysis
#
# Description:
#   Sequential execution of all meDIP-seq analysis steps with dependency
#   management. Submits jobs to SLURM and waits for completion before
#   proceeding to dependent steps.
#
# Usage:
#   ./scripts/master_script.sh              # Run complete pipeline
#   ./scripts/master_script.sh --from 3     # Start from step 3
#   ./scripts/master_script.sh --dry-run    # Preview without submitting
#
# Pipeline Steps:
#   01. FastQC - Raw read quality control
#   02. Trim Galore - Adapter trimming
#   03. Bowtie2 - Genome alignment
#   04. Filtering - Quality filtering and deduplication
#   05. BigWig - Coverage track generation
#   06. MACS2 - Peak calling
#   07. DiffBind - Differential methylation
#   08. ChIPseeker - Genomic annotation
#   09. HOMER - Motif enrichment
#   10. Integration - Multi-omics analysis
#   11. MultiQC - QC report generation
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse command line arguments
START_STEP=1
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --from)
            START_STEP="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --from N      Start from step N (default: 1)"
            echo "  --dry-run     Show commands without executing"
            echo "  --help        Show this help message"
            echo ""
            echo "Steps:"
            echo "  1. FastQC"
            echo "  2. Trim Galore"
            echo "  3. Bowtie2 alignment"
            echo "  4. Filtering & deduplication"
            echo "  5. BigWig generation"
            echo "  6. Peak calling"
            echo "  7. Differential methylation"
            echo "  8. Genomic annotation"
            echo "  9. Motif enrichment"
            echo "  10. Multi-omics integration"
            echo "  11. MultiQC report"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname "${SCRIPT_DIR}")"

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}meDIP-seq Analysis Pipeline${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""
echo "Base directory: ${BASE_DIR}"
echo "Script directory: ${SCRIPT_DIR}"
echo "Start step: ${START_STEP}"
echo "Dry run: ${DRY_RUN}"
echo ""

# Function to submit job and wait for completion
submit_and_wait() {
    local STEP=$1
    local SCRIPT=$2
    local DESCRIPTION=$3

    if [ ${START_STEP} -gt ${STEP} ]; then
        echo -e "${YELLOW}[SKIP]${NC} Step ${STEP}: ${DESCRIPTION}"
        return 0
    fi

    echo -e "${GREEN}[RUN]${NC} Step ${STEP}: ${DESCRIPTION}"
    echo "  Script: ${SCRIPT}"

    if [ ${DRY_RUN} -eq 1 ]; then
        echo "  [DRY RUN] Would submit: sbatch ${SCRIPT}"
        return 0
    fi

    # Submit job
    cd "${BASE_DIR}"
    JOB_OUTPUT=$(sbatch "${SCRIPT}")
    JOB_ID=$(echo "${JOB_OUTPUT}" | awk '{print $NF}')

    if [ -z "${JOB_ID}" ]; then
        echo -e "${RED}[ERROR]${NC} Failed to submit job"
        exit 1
    fi

    echo "  Job ID: ${JOB_ID}"
    echo "  Waiting for completion..."

    # Wait for job to complete
    while true; do
        JOB_STATE=$(squeue -j ${JOB_ID} -h -o %T 2>/dev/null || echo "COMPLETED")

        if [ -z "${JOB_STATE}" ] || [ "${JOB_STATE}" = "COMPLETED" ]; then
            echo -e "  ${GREEN}[DONE]${NC} Job completed"
            break
        elif [ "${JOB_STATE}" = "FAILED" ] || [ "${JOB_STATE}" = "CANCELLED" ]; then
            echo -e "  ${RED}[FAILED]${NC} Job ${JOB_STATE}"
            echo "  Check logs: ${BASE_DIR}/logs/"
            exit 1
        else
            echo "  Status: ${JOB_STATE}"
            sleep 60
        fi
    done

    echo ""
}

# Function to run R script directly (for steps that don't need SLURM)
run_r_script() {
    local STEP=$1
    local SCRIPT=$2
    local DESCRIPTION=$3

    if [ ${START_STEP} -gt ${STEP} ]; then
        echo -e "${YELLOW}[SKIP]${NC} Step ${STEP}: ${DESCRIPTION}"
        return 0
    fi

    echo -e "${GREEN}[RUN]${NC} Step ${STEP}: ${DESCRIPTION}"
    echo "  Script: ${SCRIPT}"

    if [ ${DRY_RUN} -eq 1 ]; then
        echo "  [DRY RUN] Would run: Rscript ${SCRIPT}"
        return 0
    fi

    cd "${BASE_DIR}"
    source /opt/common/tools/ric.cosr/miniconda3/bin/activate
    conda activate r_chipseq_env

    Rscript "${SCRIPT}"

    if [ $? -eq 0 ]; then
        echo -e "  ${GREEN}[DONE]${NC} Script completed"
    else
        echo -e "  ${RED}[FAILED]${NC} Script failed"
        exit 1
    fi

    echo ""
}

# Run pipeline steps
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}Starting Pipeline Execution${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

submit_and_wait 1 "scripts/01_fastqc.sh" "FastQC quality control"
submit_and_wait 2 "scripts/02_trim.sh" "Trim Galore adapter trimming"
submit_and_wait 3 "scripts/03_align.sh" "Bowtie2 genome alignment"
submit_and_wait 4 "scripts/04_filter.sh" "Filtering and deduplication"
submit_and_wait 5 "scripts/05_bigwig.sh" "BigWig track generation"
submit_and_wait 6 "scripts/06_peak_calling.sh" "MACS2 peak calling"

# R-based analyses can be run with SLURM or directly
echo -e "${BLUE}Running R-based analyses...${NC}"
echo ""

if [ ${START_STEP} -le 7 ]; then
    run_r_script 7 "scripts/07_differential_methylation.R" "Differential methylation analysis"
fi

if [ ${START_STEP} -le 8 ]; then
    run_r_script 8 "scripts/08_annotation.R" "Genomic annotation"
fi

submit_and_wait 9 "scripts/09_motif_analysis.sh" "HOMER motif enrichment"

if [ ${START_STEP} -le 10 ]; then
    run_r_script 10 "scripts/integrative.R" "Multi-omics integration"
fi

submit_and_wait 11 "scripts/11_multiqc.sh" "MultiQC report generation"

# Pipeline completion
echo -e "${BLUE}=========================================${NC}"
echo -e "${GREEN}Pipeline Completed Successfully!${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""
echo "Results directory: ${BASE_DIR}/results/"
echo ""
echo "Key output files:"
echo "  - MultiQC report: results/medip_multiqc_report.html"
echo "  - BigWig tracks: results/05_bigwig/*_RPKM.bw"
echo "  - Peaks: results/06_peaks/*_peaks.broadPeak"
echo "  - DMRs: results/07_differential/*_DMRs_FDR05.csv"
echo "  - Annotations: results/08_annotation/*_annotated.csv"
echo "  - Integration: results/integrative/epigenetic_regulation.csv"
echo ""
echo "Next steps:"
echo "  1. Review MultiQC report for overall quality"
echo "  2. Load BigWig tracks in IGV"
echo "  3. Examine DMRs for biological insights"
echo "  4. Investigate direct epigenetic targets"
echo ""
echo -e "${GREEN}Analysis complete!${NC}"
