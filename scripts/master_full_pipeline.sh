#!/bin/bash

################################################################################
# Script: master_full_pipeline.sh
# Purpose: Complete meDIP-seq pipeline with SLURM dependency management
#
# Description:
#   Orchestrates all meDIP-seq analysis steps using SLURM job dependencies
#   to maximize parallelism. Jobs that don't depend on each other run
#   concurrently while respecting data dependencies.
#
# Usage:
#   ./scripts/master_full_pipeline.sh              # Run complete pipeline
#   ./scripts/master_full_pipeline.sh --from 5     # Start from step 5
#   ./scripts/master_full_pipeline.sh --dry-run    # Preview without submitting
#   ./scripts/master_full_pipeline.sh --skip-core  # Skip core (01-05), run analysis only
#
# Pipeline Steps and Dependencies:
#   ┌─ 01_fastqc.sh ─────────────────────────────────────────────────────────────┐
#   │                                                                             │
#   └─ 02_trim.sh ──┬─ 03_align_array.sh ─── 04_filter.sh ──┬─ 05_bigwig.sh ──┐  │
#                   │                                        │                 │  │
#                   │                                        └─ 07_diff*.sh ───┤  │
#                   │                                                          │  │
#                   │  ┌───────────────────────────────────────────────────────┘  │
#                   │  │                                                          │
#                   │  ├─ 08_annotation.sh                                        │
#                   │  │                                                          │
#                   │  └─ 12_define_gene_sets.sh ─── 13_compute_matrices.sh ──┬──┤
#                   │                                                          │  │
#                   │                                         14_heatmaps.sh ──┤  │
#                   │                                                          │  │
#                   │                                      15_metaprofiles.sh ─┤  │
#                   │                                                          │  │
#                   │                                   16_advanced_viz.sh ────┤  │
#                   │                                                          │  │
#                   │                                   17_dmr_heatmaps.sh ────┤  │
#                   │                                                          │  │
#                   │                                   18_dmr_overlay.sh ─────┘  │
#                   │                                                             │
#                   └─ 11_multiqc.sh (runs last) ─────────────────────────────────┘
#
# NOTE: TESmut IP samples excluded from analysis (failed sample)
#       TESmut-1-INPUT is retained as common INPUT control
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname "${SCRIPT_DIR}")"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "${LOG_DIR}"

# Parse command line arguments
START_STEP=1
DRY_RUN=0
SKIP_CORE=0

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
        --skip-core)
            SKIP_CORE=1
            START_STEP=7
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --from N        Start from step N (default: 1)"
            echo "  --dry-run       Show commands without executing"
            echo "  --skip-core     Skip core processing (01-05), assume data exists"
            echo "  --help          Show this help message"
            echo ""
            echo "Pipeline Steps:"
            echo "  Core Processing:"
            echo "    1. FastQC (raw read QC)"
            echo "    2. Trim Galore (adapter trimming)"
            echo "    3. Bowtie2 alignment (parallel array job)"
            echo "    4. Filtering & deduplication"
            echo "    5. BigWig generation"
            echo ""
            echo "  Differential Analysis:"
            echo "    7. Differential methylation (MEDIPS)"
            echo "    8. Genomic annotation"
            echo ""
            echo "  Visualization Pipeline:"
            echo "    12. Define gene sets"
            echo "    13. Compute matrices"
            echo "    14. Plot heatmaps"
            echo "    15. Plot metaprofiles"
            echo "    16. Advanced visualization"
            echo "    17. DMR stratified heatmaps"
            echo "    18. DMR binding overlay"
            echo ""
            echo "  QC Report:"
            echo "    11. MultiQC (runs last)"
            echo ""
            echo "Examples:"
            echo "  $0                    # Run complete pipeline"
            echo "  $0 --dry-run          # Preview all jobs"
            echo "  $0 --from 7           # Start from differential analysis"
            echo "  $0 --skip-core        # Skip steps 1-5, run analysis only"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}      meDIP-seq Complete Analysis Pipeline            ${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""
echo "Base directory: ${BASE_DIR}"
echo "Script directory: ${SCRIPT_DIR}"
echo "Start step: ${START_STEP}"
echo "Dry run: $([[ ${DRY_RUN} -eq 1 ]] && echo 'Yes' || echo 'No')"
echo "Skip core: $([[ ${SKIP_CORE} -eq 1 ]] && echo 'Yes' || echo 'No')"
echo ""
echo -e "${YELLOW}NOTE: TESmut IP samples excluded (failed sample)${NC}"
echo -e "${YELLOW}      TESmut-1-INPUT retained as common INPUT control${NC}"
echo ""

# Track job IDs for dependencies
declare -A JOB_IDS

# Function to submit job with optional dependencies
submit_job() {
    local step=$1
    local script=$2
    local job_name=$3
    local dependencies=$4  # Comma-separated list of job names to wait for

    # Skip if before start step
    if [ ${step} -lt ${START_STEP} ]; then
        echo -e "${YELLOW}[SKIP]${NC} Step ${step}: ${job_name}"
        return 0
    fi

    # Check if script exists
    local script_path="${SCRIPT_DIR}/${script}"
    if [ ! -f "${script_path}" ]; then
        echo -e "${RED}[MISSING]${NC} Step ${step}: ${job_name} - Script not found: ${script}"
        return 1
    fi

    echo -e "${GREEN}[SUBMIT]${NC} Step ${step}: ${job_name}"
    echo "         Script: ${script}"

    # Build dependency string
    local dep_string=""
    if [ -n "${dependencies}" ]; then
        local dep_ids=""
        IFS=',' read -ra DEP_NAMES <<< "$dependencies"
        for dep_name in "${DEP_NAMES[@]}"; do
            if [ -n "${JOB_IDS[$dep_name]:-}" ]; then
                if [ -n "$dep_ids" ]; then
                    dep_ids="${dep_ids}:${JOB_IDS[$dep_name]}"
                else
                    dep_ids="${JOB_IDS[$dep_name]}"
                fi
            fi
        done
        if [ -n "$dep_ids" ]; then
            dep_string="--dependency=afterok:${dep_ids}"
            echo "         Dependencies: ${dependencies} (${dep_ids})"
        fi
    fi

    if [ ${DRY_RUN} -eq 1 ]; then
        echo -e "         ${CYAN}[DRY RUN]${NC} sbatch ${dep_string} ${script_path}"
        JOB_IDS[$job_name]="DRY_${step}"
        return 0
    fi

    # Submit job
    cd "${BASE_DIR}"
    local job_output
    if [ -n "${dep_string}" ]; then
        job_output=$(sbatch ${dep_string} "${script_path}" 2>&1)
    else
        job_output=$(sbatch "${script_path}" 2>&1)
    fi

    local job_id=$(echo "${job_output}" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")

    if [ -z "${job_id}" ]; then
        echo -e "         ${RED}[ERROR]${NC} Failed to submit: ${job_output}"
        return 1
    fi

    JOB_IDS[$job_name]="${job_id}"
    echo -e "         ${GREEN}Job ID: ${job_id}${NC}"
    echo ""
}

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}                 Submitting Jobs                      ${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

################################################################################
# CORE PROCESSING PIPELINE (Steps 1-5)
################################################################################

echo -e "${CYAN}--- Core Processing Pipeline ---${NC}"
echo ""

# Step 1: FastQC (independent - can run in parallel with everything)
submit_job 1 "01_fastqc.sh" "fastqc" ""

# Step 2: Trimming (independent of FastQC)
submit_job 2 "02_trim.sh" "trim" ""

# Step 3: Alignment (depends on trimming)
submit_job 3 "03_align_array.sh" "align" "trim"

# Step 4: Filtering (depends on alignment)
submit_job 4 "04_filter.sh" "filter" "align"

# Step 5: BigWig generation (depends on filtering)
submit_job 5 "05_bigwig.sh" "bigwig" "filter"

################################################################################
# DIFFERENTIAL ANALYSIS (Steps 7-8)
################################################################################

echo -e "${CYAN}--- Differential Analysis ---${NC}"
echo ""

# Step 7: Differential methylation (depends on filtering, parallel with BigWig)
# Use the MEDIPS version (not peak-based)
submit_job 7 "07_differential_methylation_MEDIPS.sh" "diffmeth" "filter"

# Step 8: Annotation (depends on differential analysis)
submit_job 8 "08_annotation.sh" "annotation" "diffmeth"

################################################################################
# VISUALIZATION PIPELINE (Steps 12-18)
################################################################################

echo -e "${CYAN}--- Visualization Pipeline ---${NC}"
echo ""

# Step 12: Define gene sets (depends on differential analysis)
submit_job 12 "12_define_gene_sets.sh" "genesets" "diffmeth"

# Step 13: Compute matrices (depends on BigWig + gene sets)
submit_job 13 "13_compute_matrices.sh" "matrices" "bigwig,genesets"

# The following steps can run in parallel after matrices are computed:

# Step 14: Plot heatmaps (depends on matrices)
submit_job 14 "14_plot_heatmaps.sh" "heatmaps" "matrices"

# Step 15: Plot metaprofiles (depends on matrices, parallel with heatmaps)
submit_job 15 "15_plot_metaprofiles.sh" "metaprofiles" "matrices"

# Step 16: Advanced visualization (depends on BigWig + differential, parallel)
submit_job 16 "16_advanced_visualization_v3_optimized.sh" "advancedviz" "bigwig,diffmeth"

# Step 17: DMR stratified heatmaps (depends on BigWig + differential)
submit_job 17 "17_dmr_stratified_heatmaps.sh" "dmrheatmaps" "bigwig,diffmeth"

# Step 18: DMR binding overlay (depends on BigWig + differential)
submit_job 18 "18_dmr_binding_overlay.sh" "dmroverlay" "bigwig,diffmeth"

################################################################################
# QC REPORT (Step 11 - runs last)
################################################################################

echo -e "${CYAN}--- QC Report ---${NC}"
echo ""

# Step 11: MultiQC (depends on all core steps)
# Wait for key analysis steps to complete
submit_job 11 "11_multiqc.sh" "multiqc" "fastqc,trim,align,filter,bigwig"

################################################################################
# SUMMARY
################################################################################

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}              Pipeline Submission Summary             ${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

if [ ${DRY_RUN} -eq 1 ]; then
    echo -e "${CYAN}*** DRY RUN - No jobs were actually submitted ***${NC}"
    echo ""
fi

echo "Submitted Jobs:"
echo "---------------"
printf "  %-20s %s\n" "Step" "Job ID"
printf "  %-20s %s\n" "----" "------"

for job_name in fastqc trim align filter bigwig diffmeth annotation genesets matrices heatmaps metaprofiles advancedviz dmrheatmaps dmroverlay multiqc; do
    if [ -n "${JOB_IDS[$job_name]:-}" ]; then
        printf "  %-20s %s\n" "${job_name}" "${JOB_IDS[$job_name]}"
    fi
done

echo ""
echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}                  Parallelism Map                     ${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""
echo "Jobs running concurrently at each stage:"
echo ""
echo "Stage 1: [fastqc] + [trim]"
echo "    │"
echo "Stage 2: [align]"
echo "    │"
echo "Stage 3: [filter]"
echo "    │"
echo "Stage 4: [bigwig] + [diffmeth]"
echo "    │         │"
echo "    │         └─→ [annotation]"
echo "    │         └─→ [genesets]"
echo "    │                 │"
echo "    └────────────────→[matrices]"
echo "                       │"
echo "Stage 5: [heatmaps] + [metaprofiles] + [advancedviz] + [dmrheatmaps] + [dmroverlay]"
echo "    │"
echo "Stage 6: [multiqc]"
echo ""

if [ ${DRY_RUN} -eq 0 ]; then
    echo -e "${BLUE}======================================================${NC}"
    echo -e "${BLUE}                   Monitoring                         ${NC}"
    echo -e "${BLUE}======================================================${NC}"
    echo ""
    echo "Check job status:"
    echo "  squeue -u \$USER"
    echo ""
    echo "View specific job logs:"
    echo "  tail -f ${LOG_DIR}/<step>*.out"
    echo ""
    echo "Cancel all pipeline jobs:"
    all_jobs=""
    for job_name in fastqc trim align filter bigwig diffmeth annotation genesets matrices heatmaps metaprofiles advancedviz dmrheatmaps dmroverlay multiqc; do
        if [ -n "${JOB_IDS[$job_name]:-}" ] && [[ "${JOB_IDS[$job_name]}" != DRY_* ]]; then
            all_jobs="${all_jobs} ${JOB_IDS[$job_name]}"
        fi
    done
    if [ -n "${all_jobs}" ]; then
        echo "  scancel${all_jobs}"
    fi
    echo ""
fi

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}               Expected Runtime                       ${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""
echo "Estimated total time: 8-12 hours (with parallelism)"
echo ""
echo "  Core Processing (sequential):"
echo "    01_fastqc:      ~2-4 hours"
echo "    02_trim:        ~4-8 hours"
echo "    03_align:       ~30-45 min (parallel array)"
echo "    04_filter:      ~4-6 hours"
echo "    05_bigwig:      ~2-4 hours"
echo ""
echo "  Analysis (runs in parallel with bigwig):"
echo "    07_diffmeth:    ~2-4 hours"
echo "    08_annotation:  ~30 min"
echo ""
echo "  Visualization (runs in parallel):"
echo "    12_genesets:    ~30 min"
echo "    13_matrices:    ~2-3 hours"
echo "    14_heatmaps:    ~30-60 min"
echo "    15_metaprofiles: ~20-30 min"
echo "    16_advancedviz: ~10-20 min"
echo "    17_dmrheatmaps: ~30-60 min"
echo "    18_dmroverlay:  ~30-60 min"
echo ""
echo "  QC Report:"
echo "    11_multiqc:     ~30 min"
echo ""

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}                Expected Outputs                      ${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""
echo "results/"
echo "├── 01_fastqc/              # Raw read QC reports"
echo "├── 02_trimmed/             # Trimmed FASTQ files"
echo "├── 03_aligned/             # Sorted BAM files"
echo "├── 04_filtered/            # Filtered & deduplicated BAMs"
echo "├── 05_bigwig/              # RPKM-normalized coverage tracks"
echo "├── 07_differential_MEDIPS/ # DMRs and differential analysis"
echo "├── 08_annotation/          # Annotated DMRs"
echo "├── 12_gene_sets/           # Gene lists and promoter BEDs"
echo "├── 13_matrices/            # deepTools signal matrices"
echo "├── 14_heatmaps/            # Methylation heatmaps"
echo "├── 15_metaprofiles/        # Average signal profiles"
echo "├── 16_advanced_visualization/ # Integrated plots"
echo "├── 17_dmr_stratified/      # DMR-stratified heatmaps"
echo "├── 18_dmr_binding/         # DMR + binding overlay"
echo "└── medip_multiqc_report.html"
echo ""

echo -e "${GREEN}Pipeline submission complete!${NC}"
echo ""
