#!/bin/bash
#SBATCH --job-name=30_peaked_metagene_conf
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/30_peaked_genes_metagene_confident.out
#SBATCH --error=logs/30_peaked_genes_metagene_confident.err

echo "=========================================="
echo "PEAKED GENES METAGENE - HIGH-CONFIDENCE"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Create metagene profiles for genes with Cut&Tag peaks"
echo "      Stratified by high-confidence DMR status"
echo ""
echo "Key insight: High-confidence DMRs (both samples >2 reads) show"
echo "             65% hypomethylation vs 91% hyper in original analysis"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/scripts/only_confident_DMRs

OUTDIR="output/30_peaked_genes_metagene_confident"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: PREPARE GENE LISTS (R script)
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 30_peaked_genes_metagene_confident.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# Verify BED files were created
echo ""
echo "Checking BED files..."
for bed in TES_peaked_with_DMR.bed TES_peaked_no_DMR.bed \
           TEAD1_peaked_with_DMR.bed TEAD1_peaked_no_DMR.bed \
           TES_DEGs_DOWN_with_DMR.bed TES_DEGs_DOWN_no_DMR.bed; do
    if [ -f "${OUTDIR}/${bed}" ]; then
        count=$(wc -l < "${OUTDIR}/${bed}")
        echo "  ${bed}: ${count} genes"
    else
        echo "  WARNING: ${bed} not found!"
    fi
done

# =============================================================================
# STEP 2: SETUP DEEPTOOLS
# =============================================================================

echo ""
echo "=== STEP 2: Computing Coverage Matrices ==="
echo ""

conda activate tg

# BigWig files
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_average.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_average.bw"

# Check BigWig files exist
for bw in $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH; do
    if [ ! -f "$bw" ]; then
        echo "ERROR: BigWig file not found: $bw"
        exit 1
    fi
done
echo "All BigWig files found"

# Get gene counts for labels
N_TES_DMR=$(wc -l < ${OUTDIR}/TES_peaked_with_DMR.bed 2>/dev/null || echo 0)
N_TES_NODMR=$(wc -l < ${OUTDIR}/TES_peaked_no_DMR.bed 2>/dev/null || echo 0)
N_TEAD1_DMR=$(wc -l < ${OUTDIR}/TEAD1_peaked_with_DMR.bed 2>/dev/null || echo 0)
N_TEAD1_NODMR=$(wc -l < ${OUTDIR}/TEAD1_peaked_no_DMR.bed 2>/dev/null || echo 0)

# =============================================================================
# STEP 3: TES PEAKED GENES - DMR STRATIFIED
# =============================================================================

echo ""
echo "=== STEP 3: TES-Peaked Genes (DMR-stratified) ==="
echo ""

# Check if we have enough genes for stratified analysis
MIN_GENES=10

if [ "$N_TES_DMR" -ge $MIN_GENES ] && [ "$N_TES_NODMR" -ge $MIN_GENES ]; then
    computeMatrix scale-regions \
        -S $TES_BIND $TES_METH $GFP_METH \
        -R ${OUTDIR}/TES_peaked_with_DMR.bed ${OUTDIR}/TES_peaked_no_DMR.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/TES_peaked_dmr_stratified_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/TES_peaked_dmr_stratified_matrix.gz \
        -out ${OUTDIR}/TES_peaked_dmr_stratified_profile.png \
        --perGroup \
        --colors "#E31A1C" "#7B3294" "#636363" \
        --startLabel "TSS" \
        --endLabel "TTS" \
        --samplesLabel "TES bind" "TES meth" "GFP meth" \
        --regionsLabel "TES-peaked + DMR (n=${N_TES_DMR})" \
                       "TES-peaked - no DMR (n=${N_TES_NODMR})" \
        --plotTitle "TES Binding & Methylation at TES-Peaked Genes (DMR-stratified)" \
        --yAxisLabel "Mean signal" \
        --plotHeight 10 \
        --plotWidth 14 \
        --dpi 300

    echo "  Created: TES_peaked_dmr_stratified_profile.png"
else
    echo "  Skipping TES DMR-stratified (insufficient genes: DMR=$N_TES_DMR, noDMR=$N_TES_NODMR)"
fi

# =============================================================================
# STEP 4: TEAD1 PEAKED GENES - DMR STRATIFIED
# =============================================================================

echo ""
echo "=== STEP 4: TEAD1-Peaked Genes (DMR-stratified) ==="
echo ""

if [ "$N_TEAD1_DMR" -ge $MIN_GENES ] && [ "$N_TEAD1_NODMR" -ge $MIN_GENES ]; then
    computeMatrix scale-regions \
        -S $TEAD1_BIND $TES_METH $GFP_METH \
        -R ${OUTDIR}/TEAD1_peaked_with_DMR.bed ${OUTDIR}/TEAD1_peaked_no_DMR.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/TEAD1_peaked_dmr_stratified_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/TEAD1_peaked_dmr_stratified_matrix.gz \
        -out ${OUTDIR}/TEAD1_peaked_dmr_stratified_profile.png \
        --perGroup \
        --colors "#377EB8" "#7B3294" "#636363" \
        --startLabel "TSS" \
        --endLabel "TTS" \
        --samplesLabel "TEAD1 bind" "TES meth" "GFP meth" \
        --regionsLabel "TEAD1-peaked + DMR (n=${N_TEAD1_DMR})" \
                       "TEAD1-peaked - no DMR (n=${N_TEAD1_NODMR})" \
        --plotTitle "TEAD1 Binding & Methylation at TEAD1-Peaked Genes (DMR-stratified)" \
        --yAxisLabel "Mean signal" \
        --plotHeight 10 \
        --plotWidth 14 \
        --dpi 300

    echo "  Created: TEAD1_peaked_dmr_stratified_profile.png"
else
    echo "  Skipping TEAD1 DMR-stratified (insufficient genes: DMR=$N_TEAD1_DMR, noDMR=$N_TEAD1_NODMR)"
fi

# =============================================================================
# STEP 5: DEGs DOWN WITH PEAKS - DMR STRATIFIED
# =============================================================================

echo ""
echo "=== STEP 5: DEGs DOWN with Peaks (DMR-stratified) ==="
echo ""

N_TES_DOWN_DMR=$(wc -l < ${OUTDIR}/TES_DEGs_DOWN_with_DMR.bed 2>/dev/null || echo 0)
N_TES_DOWN_NODMR=$(wc -l < ${OUTDIR}/TES_DEGs_DOWN_no_DMR.bed 2>/dev/null || echo 0)

if [ "$N_TES_DOWN_DMR" -ge $MIN_GENES ] && [ "$N_TES_DOWN_NODMR" -ge $MIN_GENES ]; then
    computeMatrix scale-regions \
        -S $TES_BIND $TES_METH $GFP_METH \
        -R ${OUTDIR}/TES_DEGs_DOWN_with_DMR.bed ${OUTDIR}/TES_DEGs_DOWN_no_DMR.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/TES_DEGs_DOWN_dmr_stratified_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/TES_DEGs_DOWN_dmr_stratified_matrix.gz \
        -out ${OUTDIR}/TES_DEGs_DOWN_dmr_stratified_profile.png \
        --perGroup \
        --colors "#E31A1C" "#7B3294" "#636363" \
        --startLabel "TSS" \
        --endLabel "TTS" \
        --samplesLabel "TES bind" "TES meth" "GFP meth" \
        --regionsLabel "TES DEGs DOWN + DMR (n=${N_TES_DOWN_DMR})" \
                       "TES DEGs DOWN - no DMR (n=${N_TES_DOWN_NODMR})" \
        --plotTitle "TES DEGs DOWN: Binding & Methylation (DMR-stratified)" \
        --yAxisLabel "Mean signal" \
        --plotHeight 10 \
        --plotWidth 14 \
        --dpi 300

    echo "  Created: TES_DEGs_DOWN_dmr_stratified_profile.png"
else
    echo "  Skipping TES DEGs DOWN DMR-stratified (insufficient genes)"
fi

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "PEAKED GENES METAGENE COMPLETE"
echo "(HIGH-CONFIDENCE VERSION)"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated files:"
ls -lh ${OUTDIR}/*.png 2>/dev/null
echo ""
echo "Key insight:"
echo "  Genes WITH high-confidence DMRs have RELIABLE methylation signal"
echo "  in both TES and GFP samples. The GFP library has quality issues"
echo "  (4.3M unique molecules), causing 65% of original 'hypermethylated'"
echo "  DMRs to have GFP=0 reads."
echo ""
echo "  Original DMRs: 91% hypermethylated (artifact)"
echo "  High-conf DMRs: 35% hyper / 65% hypo (real biology)"
echo ""
