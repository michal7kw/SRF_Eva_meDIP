#!/bin/bash
#SBATCH --job-name=22_cpg_island_context
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=logs/22_cpg_island_context.out
#SBATCH --error=logs/22_cpg_island_context.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# CpG Island Context Analysis for DMRs
# ============================================================================
#
# Purpose: Stratify DMRs by CpG island geography (Island/Shore/Shelf/Open Sea)
#          to understand methylation context and avoid analysis bias
#
# CpG Geography Definitions (UCSC standard):
#   - CpG Island: High CpG density regions
#   - Shore: ±2kb flanking islands (tissue-specific methylation)
#   - Shelf: 2-4kb flanking islands
#   - Open Sea: >4kb from any island
#
# Best Practices Implemented:
#   - Proper genome background for enrichment analysis
#   - Account for region size differences
#   - Report both count and proportion statistics
#   - Consistent chromosome naming (UCSC format)
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/22_cpg_context"
mkdir -p "$OUTDIR"/{beds,plots}
mkdir -p logs

echo "=============================================="
echo "CpG Island Context Analysis for DMRs"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ============================================================================
# Define input files
# ============================================================================

# CpG islands from HOMER (UCSC format: chr1, chr2, etc.)
CPG_ISLANDS="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/annotation_enrichment/share/homer/data/genomes/hg38/annotations/basic/cpgIsland.ann.txt"

# Chromosome sizes (UCSC format)
CHROM_SIZES="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/annotation_enrichment/share/homer/data/genomes/hg38/chrom.sizes"

# DMR file (Ensembl format: 1, 2, etc. - needs conversion)
DMR_FILE="results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

echo ""
echo "=== Checking input files ==="

# Verify input files exist
for f in "$CPG_ISLANDS" "$CHROM_SIZES" "$DMR_FILE"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  ERROR: Not found: $f"
        exit 1
    fi
done

# ============================================================================
# Step 1: Convert CpG islands to BED format
# ============================================================================

echo ""
echo "=== Step 1: Converting CpG islands to BED format ==="

# HOMER format: ID chr start end strand annotation
# Convert to BED: chr start end name
awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $3, $4, $1}' "$CPG_ISLANDS" | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/cpg_islands.bed"

ISLAND_COUNT=$(wc -l < "$OUTDIR/beds/cpg_islands.bed")
echo "  CpG Islands: $ISLAND_COUNT"

# ============================================================================
# Step 2: Generate CpG shores (±2kb from islands)
# ============================================================================

echo ""
echo "=== Step 2: Generating CpG shores (±2kb from islands) ==="

# Create standard chromosome sizes (only main chromosomes for bedtools)
grep -E "^chr[0-9]+\s|^chrX\s|^chrY\s" "$CHROM_SIZES" > "$OUTDIR/beds/hg38_main_chroms.sizes"

# Left shore (upstream 2kb)
bedtools flank -i "$OUTDIR/beds/cpg_islands.bed" \
    -g "$OUTDIR/beds/hg38_main_chroms.sizes" \
    -l 2000 -r 0 | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/cpg_shores_left.bed"

# Right shore (downstream 2kb)
bedtools flank -i "$OUTDIR/beds/cpg_islands.bed" \
    -g "$OUTDIR/beds/hg38_main_chroms.sizes" \
    -l 0 -r 2000 | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/cpg_shores_right.bed"

# Combine and merge shores, then remove overlap with islands
cat "$OUTDIR/beds/cpg_shores_left.bed" "$OUTDIR/beds/cpg_shores_right.bed" | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - | \
    bedtools subtract -a - -b "$OUTDIR/beds/cpg_islands.bed" > "$OUTDIR/beds/cpg_shores.bed"

SHORE_COUNT=$(wc -l < "$OUTDIR/beds/cpg_shores.bed")
echo "  CpG Shores (±2kb, non-overlapping with islands): $SHORE_COUNT regions"

# ============================================================================
# Step 3: Generate CpG shelves (2-4kb from islands)
# ============================================================================

echo ""
echo "=== Step 3: Generating CpG shelves (2-4kb from islands) ==="

# Expand islands by 4kb on each side
bedtools slop -i "$OUTDIR/beds/cpg_islands.bed" \
    -g "$OUTDIR/beds/hg38_main_chroms.sizes" \
    -b 4000 | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > "$OUTDIR/beds/cpg_extended_4kb.bed"

# Expand islands by 2kb on each side (for subtraction)
bedtools slop -i "$OUTDIR/beds/cpg_islands.bed" \
    -g "$OUTDIR/beds/hg38_main_chroms.sizes" \
    -b 2000 | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > "$OUTDIR/beds/cpg_extended_2kb.bed"

# Shelves = 4kb extension minus 2kb extension
bedtools subtract -a "$OUTDIR/beds/cpg_extended_4kb.bed" \
    -b "$OUTDIR/beds/cpg_extended_2kb.bed" > "$OUTDIR/beds/cpg_shelves.bed"

SHELF_COUNT=$(wc -l < "$OUTDIR/beds/cpg_shelves.bed")
echo "  CpG Shelves (2-4kb from islands): $SHELF_COUNT regions"

# ============================================================================
# Step 4: Convert DMRs to BED format (Ensembl -> UCSC)
# ============================================================================

echo ""
echo "=== Step 4: Converting DMRs to BED format ==="

# DMR file is CSV with header: "chr","start","stop","CpG_count",...
# Chr column uses Ensembl format (1, 2, 3) - need to add "chr" prefix
tail -n +2 "$DMR_FILE" | \
    awk -F',' 'BEGIN{OFS="\t"} {
        gsub(/"/, "", $1);  # Remove quotes from chr
        gsub(/"/, "", $2);  # Remove quotes from start
        gsub(/"/, "", $3);  # Remove quotes from stop
        gsub(/"/, "", $7);  # Remove quotes from logFC
        # Add chr prefix if not present
        chr = $1;
        if (chr !~ /^chr/) chr = "chr" chr;
        print chr, $2, $3, "DMR_" NR, $7
    }' | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/dmrs.bed"

DMR_COUNT=$(wc -l < "$OUTDIR/beds/dmrs.bed")
echo "  Total DMRs: $DMR_COUNT"

# ============================================================================
# Step 5: Annotate DMRs with CpG context
# ============================================================================

echo ""
echo "=== Step 5: Annotating DMRs with CpG context ==="

# Intersect DMRs with each CpG region type
# -wa writes original DMR, -u reports DMR once even if multiple overlaps

# Islands
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" -b "$OUTDIR/beds/cpg_islands.bed" \
    -wa -u > "$OUTDIR/beds/dmrs_in_islands.bed" 2>/dev/null || touch "$OUTDIR/beds/dmrs_in_islands.bed"

# Shores
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" -b "$OUTDIR/beds/cpg_shores.bed" \
    -wa -u > "$OUTDIR/beds/dmrs_in_shores.bed" 2>/dev/null || touch "$OUTDIR/beds/dmrs_in_shores.bed"

# Shelves
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" -b "$OUTDIR/beds/cpg_shelves.bed" \
    -wa -u > "$OUTDIR/beds/dmrs_in_shelves.bed" 2>/dev/null || touch "$OUTDIR/beds/dmrs_in_shelves.bed"

# Combined islands + shores + shelves (for finding open sea)
cat "$OUTDIR/beds/cpg_islands.bed" "$OUTDIR/beds/cpg_shores.bed" "$OUTDIR/beds/cpg_shelves.bed" | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > "$OUTDIR/beds/cpg_all_regions.bed"

# Open Sea (DMRs not overlapping any CpG region)
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" -b "$OUTDIR/beds/cpg_all_regions.bed" \
    -v > "$OUTDIR/beds/dmrs_in_opensea.bed" 2>/dev/null || touch "$OUTDIR/beds/dmrs_in_opensea.bed"

# Count DMRs in each context
ISLAND_DMRS=$(wc -l < "$OUTDIR/beds/dmrs_in_islands.bed")
SHORE_DMRS=$(wc -l < "$OUTDIR/beds/dmrs_in_shores.bed")
SHELF_DMRS=$(wc -l < "$OUTDIR/beds/dmrs_in_shelves.bed")
OPENSEA_DMRS=$(wc -l < "$OUTDIR/beds/dmrs_in_opensea.bed")

echo "  DMRs in CpG Islands: $ISLAND_DMRS"
echo "  DMRs in CpG Shores: $SHORE_DMRS"
echo "  DMRs in CpG Shelves: $SHELF_DMRS"
echo "  DMRs in Open Sea: $OPENSEA_DMRS"

# ============================================================================
# Step 6: Create annotated DMR file with CpG context
# ============================================================================

echo ""
echo "=== Step 6: Creating annotated DMR file ==="

# Create a combined annotation file
{
    echo -e "chr\tstart\tend\tdmr_id\tlogFC\tcpg_context"

    # Islands
    awk 'BEGIN{OFS="\t"} {print $0, "Island"}' "$OUTDIR/beds/dmrs_in_islands.bed"

    # Shores (exclude those already in islands)
    bedtools intersect -a "$OUTDIR/beds/dmrs_in_shores.bed" -b "$OUTDIR/beds/dmrs_in_islands.bed" -v | \
        awk 'BEGIN{OFS="\t"} {print $0, "Shore"}'

    # Shelves (exclude those already in islands or shores)
    cat "$OUTDIR/beds/dmrs_in_islands.bed" "$OUTDIR/beds/dmrs_in_shores.bed" | \
        sort -k1,1 -k2,2n | bedtools merge -i - > "$OUTDIR/beds/temp_island_shore.bed"
    bedtools intersect -a "$OUTDIR/beds/dmrs_in_shelves.bed" -b "$OUTDIR/beds/temp_island_shore.bed" -v | \
        awk 'BEGIN{OFS="\t"} {print $0, "Shelf"}'

    # Open Sea
    awk 'BEGIN{OFS="\t"} {print $0, "OpenSea"}' "$OUTDIR/beds/dmrs_in_opensea.bed"

} | sort -k1,1 -k2,2n > "$OUTDIR/dmr_cpg_annotation.tsv"

rm -f "$OUTDIR/beds/temp_island_shore.bed"

echo "  Created: $OUTDIR/dmr_cpg_annotation.tsv"

# ============================================================================
# Step 7: Calculate genome background for enrichment
# ============================================================================

echo ""
echo "=== Step 7: Calculating genome background ==="

# Calculate total genome coverage for each CpG context
ISLAND_BP=$(awk '{sum += $3-$2} END {print sum}' "$OUTDIR/beds/cpg_islands.bed")
SHORE_BP=$(awk '{sum += $3-$2} END {print sum}' "$OUTDIR/beds/cpg_shores.bed")
SHELF_BP=$(awk '{sum += $3-$2} END {print sum}' "$OUTDIR/beds/cpg_shelves.bed")
ALL_CPG_BP=$(awk '{sum += $3-$2} END {print sum}' "$OUTDIR/beds/cpg_all_regions.bed")

# Total mappable genome (main chromosomes only)
GENOME_BP=$(awk '{sum += $2} END {print sum}' "$OUTDIR/beds/hg38_main_chroms.sizes")
OPENSEA_BP=$((GENOME_BP - ALL_CPG_BP))

echo "  Genome coverage:"
echo "    CpG Islands: $ISLAND_BP bp"
echo "    CpG Shores: $SHORE_BP bp"
echo "    CpG Shelves: $SHELF_BP bp"
echo "    Open Sea: $OPENSEA_BP bp"
echo "    Total genome: $GENOME_BP bp"

# Save background statistics
{
    echo -e "context\tregion_count\ttotal_bp\tpct_genome"
    echo -e "Island\t$ISLAND_COUNT\t$ISLAND_BP\t$(echo "scale=4; $ISLAND_BP * 100 / $GENOME_BP" | bc)"
    echo -e "Shore\t$SHORE_COUNT\t$SHORE_BP\t$(echo "scale=4; $SHORE_BP * 100 / $GENOME_BP" | bc)"
    echo -e "Shelf\t$SHELF_COUNT\t$SHELF_BP\t$(echo "scale=4; $SHELF_BP * 100 / $GENOME_BP" | bc)"
    echo -e "OpenSea\tNA\t$OPENSEA_BP\t$(echo "scale=4; $OPENSEA_BP * 100 / $GENOME_BP" | bc)"
} > "$OUTDIR/genome_background.tsv"

echo "  Created: $OUTDIR/genome_background.tsv"

# ============================================================================
# Step 8: Run R analysis for statistics and visualization
# ============================================================================

echo ""
echo "=== Step 8: Running R analysis ==="

# Switch to R environment
conda activate r_chipseq_env

Rscript scripts/22_cpg_island_context.R

echo ""
echo "=============================================="
echo "CpG Context Analysis Complete"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Key outputs:"
echo "  BED files:"
ls -1 "$OUTDIR/beds/"*.bed 2>/dev/null | while read f; do echo "    $(basename $f)"; done
echo ""
echo "  Analysis files:"
ls -1 "$OUTDIR/"*.tsv "$OUTDIR/"*.csv 2>/dev/null | while read f; do echo "    $(basename $f)"; done
echo ""
echo "  Plots:"
ls -1 "$OUTDIR/plots/"*.png 2>/dev/null | while read f; do echo "    $(basename $f)"; done
echo ""
echo "Finished: $(date)"
