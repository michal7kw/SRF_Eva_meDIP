#!/bin/bash
#SBATCH --job-name=23_chromatin_state_integration
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=logs/23_chromatin_state_integration.out
#SBATCH --error=logs/23_chromatin_state_integration.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# Chromatin State Integration with DMRs
# ============================================================================
#
# Purpose: Overlay DMRs with chromatin state maps from Roadmap Epigenomics
#          to understand epigenetic context of methylation changes
#
# ChromHMM 15-state model (Roadmap Epigenomics):
#   1_TssA       - Active TSS
#   2_TssAFlnk   - Flanking Active TSS
#   3_TxFlnk     - Transcr. at gene 5' and 3'
#   4_Tx         - Strong transcription
#   5_TxWk       - Weak transcription
#   6_EnhG       - Genic enhancers
#   7_Enh        - Enhancers
#   8_ZNF/Rpts   - ZNF genes & repeats
#   9_Het        - Heterochromatin
#   10_TssBiv    - Bivalent/Poised TSS
#   11_BivFlnk   - Flanking Bivalent TSS/Enh
#   12_EnhBiv    - Bivalent Enhancer
#   13_ReprPC    - Repressed PolyComb
#   14_ReprPCWk  - Weak Repressed PolyComb
#   15_Quies     - Quiescent/Low
#
# Reference cell types used:
#   - E081: Fetal Brain (most relevant for glioblastoma)
#   - E017: IMR90 (alternative reference)
#
# Best Practices:
#   - Multiple reference cell types for comparison
#   - Proper genome background for enrichment
#   - Clear documentation of reference limitations
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/23_chromatin_states"
mkdir -p "$OUTDIR"/{beds,plots,reference}
mkdir -p logs

echo "=============================================="
echo "Chromatin State Integration Analysis"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ============================================================================
# Define input files and reference URLs
# ============================================================================

# DMR file (Ensembl format)
DMR_FILE="results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

# Roadmap Epigenomics ChromHMM files (15-state model, hg38 liftover)
# Using multiple cell types for robustness
CHROMHMM_BASE="https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final"

# Primary: E081 - Fetal Brain (most relevant for glioblastoma context)
CHROMHMM_E081="E081_15_coreMarks_hg38lift_mnemonics.bed.gz"

# Secondary: E017 - IMR90 fetal lung fibroblasts (well-characterized)
CHROMHMM_E017="E017_15_coreMarks_hg38lift_mnemonics.bed.gz"

# Chromosome sizes
CHROM_SIZES="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/annotation_enrichment/share/homer/data/genomes/hg38/chrom.sizes"

echo ""
echo "=== Checking input files ==="

if [[ -f "$DMR_FILE" ]]; then
    echo "  Found DMR file: $DMR_FILE"
else
    echo "  ERROR: DMR file not found: $DMR_FILE"
    exit 1
fi

# ============================================================================
# Step 1: Download ChromHMM reference files
# ============================================================================

echo ""
echo "=== Step 1: Downloading ChromHMM reference files ==="

# Download E081 (Fetal Brain)
if [[ -f "$OUTDIR/reference/E081_chromhmm.bed" ]]; then
    echo "  E081 (Fetal Brain) already downloaded"
else
    echo "  Downloading E081 (Fetal Brain)..."
    wget -q -O "$OUTDIR/reference/$CHROMHMM_E081" "$CHROMHMM_BASE/$CHROMHMM_E081" 2>/dev/null || {
        echo "  WARNING: Could not download E081 from Roadmap"
        echo "  Trying alternative source..."
        # Alternative: try ENCODE
        echo "  If download fails, please manually download from:"
        echo "    $CHROMHMM_BASE/$CHROMHMM_E081"
    }

    if [[ -f "$OUTDIR/reference/$CHROMHMM_E081" ]]; then
        gunzip -c "$OUTDIR/reference/$CHROMHMM_E081" > "$OUTDIR/reference/E081_chromhmm.bed"
        echo "  Downloaded and extracted E081"
    fi
fi

# Download E017 (IMR90) as alternative
if [[ -f "$OUTDIR/reference/E017_chromhmm.bed" ]]; then
    echo "  E017 (IMR90) already downloaded"
else
    echo "  Downloading E017 (IMR90)..."
    wget -q -O "$OUTDIR/reference/$CHROMHMM_E017" "$CHROMHMM_BASE/$CHROMHMM_E017" 2>/dev/null || {
        echo "  WARNING: Could not download E017"
    }

    if [[ -f "$OUTDIR/reference/$CHROMHMM_E017" ]]; then
        gunzip -c "$OUTDIR/reference/$CHROMHMM_E017" > "$OUTDIR/reference/E017_chromhmm.bed"
        echo "  Downloaded and extracted E017"
    fi
fi

# Check if at least one reference is available
if [[ ! -f "$OUTDIR/reference/E081_chromhmm.bed" && ! -f "$OUTDIR/reference/E017_chromhmm.bed" ]]; then
    echo ""
    echo "ERROR: No ChromHMM reference files available."
    echo "Please manually download from Roadmap Epigenomics:"
    echo "  $CHROMHMM_BASE/$CHROMHMM_E081"
    echo ""
    echo "Then extract to: $OUTDIR/reference/E081_chromhmm.bed"
    exit 1
fi

# Use E081 if available, otherwise E017
if [[ -f "$OUTDIR/reference/E081_chromhmm.bed" ]]; then
    CHROMHMM_REF="$OUTDIR/reference/E081_chromhmm.bed"
    REF_NAME="E081_FetalBrain"
    echo "  Using E081 (Fetal Brain) as primary reference"
else
    CHROMHMM_REF="$OUTDIR/reference/E017_chromhmm.bed"
    REF_NAME="E017_IMR90"
    echo "  Using E017 (IMR90) as primary reference"
fi

# Show reference statistics
echo ""
echo "  Reference file statistics:"
wc -l "$CHROMHMM_REF" | awk '{print "    Total regions: " $1}'
cut -f4 "$CHROMHMM_REF" | sort | uniq -c | sort -rn | head -5 | while read count state; do
    echo "    $state: $count"
done

# ============================================================================
# Step 2: Convert DMRs to BED format (Ensembl -> UCSC)
# ============================================================================

echo ""
echo "=== Step 2: Converting DMRs to BED format ==="

# DMR file format: "chr","start","stop","CpG_count","mean_group1","mean_group2","logFC",...
tail -n +2 "$DMR_FILE" | \
    awk -F',' 'BEGIN{OFS="\t"} {
        gsub(/"/, "", $1);
        gsub(/"/, "", $2);
        gsub(/"/, "", $3);
        gsub(/"/, "", $7);
        chr = $1;
        if (chr !~ /^chr/) chr = "chr" chr;
        print chr, $2, $3, "DMR_" NR, $7
    }' | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/dmrs.bed"

DMR_COUNT=$(wc -l < "$OUTDIR/beds/dmrs.bed")
echo "  Total DMRs: $DMR_COUNT"

# ============================================================================
# Step 3: Intersect DMRs with chromatin states
# ============================================================================

echo ""
echo "=== Step 3: Intersecting DMRs with chromatin states ==="

# Use bedtools intersect with -loj to keep all DMRs even if no overlap
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" \
    -b "$CHROMHMM_REF" \
    -loj > "$OUTDIR/beds/dmr_chromatin_intersect.bed"

# Count DMRs per chromatin state
echo "  DMRs per chromatin state:"
cut -f9 "$OUTDIR/beds/dmr_chromatin_intersect.bed" | \
    sort | uniq -c | sort -rn | head -20

# ============================================================================
# Step 4: Create annotated DMR file with chromatin state
# ============================================================================

echo ""
echo "=== Step 4: Creating annotated DMR file ==="

# For DMRs overlapping multiple states, assign the state with maximum overlap
# Using bedtools intersect with -wo for overlap length
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" \
    -b "$CHROMHMM_REF" \
    -wo > "$OUTDIR/beds/dmr_chromatin_overlap.bed" 2>/dev/null || touch "$OUTDIR/beds/dmr_chromatin_overlap.bed"

# For DMRs with no overlap, label as "Unknown"
bedtools intersect -a "$OUTDIR/beds/dmrs.bed" \
    -b "$CHROMHMM_REF" \
    -v > "$OUTDIR/beds/dmrs_no_chromatin.bed" 2>/dev/null || touch "$OUTDIR/beds/dmrs_no_chromatin.bed"

# Create consolidated annotation (best overlap for each DMR)
{
    echo -e "chr\tstart\tend\tdmr_id\tlogFC\tchromatin_state\toverlap_bp"

    # DMRs with chromatin state overlap
    if [[ -s "$OUTDIR/beds/dmr_chromatin_overlap.bed" ]]; then
        # Sort by DMR and overlap length, keep best match
        sort -k4,4 -k10,10nr "$OUTDIR/beds/dmr_chromatin_overlap.bed" | \
            awk 'BEGIN{OFS="\t"} !seen[$4]++ {print $1, $2, $3, $4, $5, $9, $10}'
    fi

    # DMRs without overlap
    if [[ -s "$OUTDIR/beds/dmrs_no_chromatin.bed" ]]; then
        awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, "Unknown", 0}' "$OUTDIR/beds/dmrs_no_chromatin.bed"
    fi

} | sort -k1,1 -k2,2n > "$OUTDIR/dmr_chromatin_annotation.tsv"

echo "  Created: dmr_chromatin_annotation.tsv"

# ============================================================================
# Step 5: Calculate genome background
# ============================================================================

echo ""
echo "=== Step 5: Calculating genome background ==="

# Calculate total bp per chromatin state
{
    echo -e "chromatin_state\ttotal_bp\tregion_count"
    awk 'BEGIN{OFS="\t"} {
        state = $4;
        bp = $3 - $2;
        count[state]++;
        total[state] += bp;
    } END {
        for (s in total) {
            print s, total[s], count[s]
        }
    }' "$CHROMHMM_REF" | sort -k2,2nr
} > "$OUTDIR/genome_background_chromatin.tsv"

echo "  Created: genome_background_chromatin.tsv"

# Total genome coverage
TOTAL_GENOME_BP=$(awk 'NR>1 {sum += $2} END {print sum}' "$OUTDIR/genome_background_chromatin.tsv")
echo "  Total genome covered by ChromHMM: $TOTAL_GENOME_BP bp"

# ============================================================================
# Step 6: Run R analysis for statistics and visualization
# ============================================================================

echo ""
echo "=== Step 6: Running R analysis ==="

# Create R script for analysis
cat > "$OUTDIR/analyze_chromatin_states.R" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
})

cat("Chromatin State Analysis - R Script\n")
cat("====================================\n\n")

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP")
OUTDIR <- "results/23_chromatin_states"
PLOTDIR <- file.path(OUTDIR, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

# Load data
dmr_annot <- read.delim(file.path(OUTDIR, "dmr_chromatin_annotation.tsv"),
                        stringsAsFactors = FALSE)
background <- read.delim(file.path(OUTDIR, "genome_background_chromatin.tsv"),
                         stringsAsFactors = FALSE)

cat("Loaded", nrow(dmr_annot), "annotated DMRs\n")
cat("Loaded", nrow(background), "chromatin states\n\n")

# Define chromatin state categories and colors
state_categories <- data.frame(
    state = c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
              "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het",
              "10_TssBiv", "11_BivFlnk", "12_EnhBiv",
              "13_ReprPC", "14_ReprPCWk", "15_Quies", "Unknown"),
    category = c("Active_TSS", "Active_TSS", "Transcription", "Transcription", "Transcription",
                 "Enhancer", "Enhancer", "ZNF_Repeats", "Heterochromatin",
                 "Bivalent", "Bivalent", "Bivalent",
                 "Repressed", "Repressed", "Quiescent", "Unknown"),
    color = c("#FF0000", "#FF6969", "#00FF00", "#008000", "#006400",
              "#FFFF00", "#FFD700", "#66CDAA", "#8A91D0",
              "#CD5C5C", "#BDB76B", "#808000",
              "#808080", "#C0C0C0", "#FFFFFF", "#000000")
)

# Summary statistics by chromatin state
dmr_summary <- dmr_annot %>%
    group_by(chromatin_state) %>%
    summarise(
        n_dmrs = n(),
        mean_logFC = mean(logFC, na.rm = TRUE),
        median_logFC = median(logFC, na.rm = TRUE),
        n_hyper = sum(logFC > 0, na.rm = TRUE),
        n_hypo = sum(logFC < 0, na.rm = TRUE),
        pct_hyper = n_hyper / n() * 100,
        .groups = "drop"
    ) %>%
    mutate(pct_total = n_dmrs / sum(n_dmrs) * 100)

# Add genome background
total_genome_bp <- sum(background$total_bp)
background_pct <- background %>%
    mutate(pct_genome = total_bp / total_genome_bp * 100) %>%
    select(chromatin_state, total_bp, pct_genome)

dmr_summary <- dmr_summary %>%
    left_join(background_pct, by = "chromatin_state")

# Calculate enrichment
dmr_summary <- dmr_summary %>%
    mutate(
        expected_dmrs = sum(n_dmrs) * (pct_genome / 100),
        fold_enrichment = n_dmrs / expected_dmrs,
        log2_enrichment = log2(fold_enrichment)
    )

# Handle Inf values
dmr_summary$log2_enrichment[is.infinite(dmr_summary$log2_enrichment)] <- NA

print(head(dmr_summary, 20))

# Save summary
write.csv(dmr_summary, file.path(OUTDIR, "chromatin_state_summary.csv"),
          row.names = FALSE)

# ============================================================================
# Visualization 1: DMR distribution by chromatin state
# ============================================================================

# Order states by number of DMRs
state_order <- dmr_summary %>%
    arrange(desc(n_dmrs)) %>%
    pull(chromatin_state)

dmr_summary$chromatin_state <- factor(dmr_summary$chromatin_state,
                                       levels = state_order)

p1 <- ggplot(dmr_summary, aes(x = chromatin_state, y = n_dmrs)) +
    geom_bar(stat = "identity", fill = "#3182BD", color = "black", linewidth = 0.3) +
    geom_text(aes(label = n_dmrs), vjust = -0.3, size = 3) +
    labs(title = "DMR Distribution by Chromatin State",
         subtitle = paste0("Reference: Roadmap Epigenomics, Total DMRs: ", sum(dmr_summary$n_dmrs)),
         x = "Chromatin State", y = "Number of DMRs") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave(file.path(PLOTDIR, "dmr_by_chromatin_state.png"), p1,
       width = 12, height = 7, dpi = 300)
cat("Saved: dmr_by_chromatin_state.png\n")

# ============================================================================
# Visualization 2: Enrichment plot
# ============================================================================

p2 <- ggplot(dmr_summary %>% filter(!is.na(log2_enrichment)),
             aes(x = reorder(chromatin_state, -log2_enrichment), y = log2_enrichment)) +
    geom_bar(stat = "identity", aes(fill = log2_enrichment > 0),
             color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                      guide = "none") +
    labs(title = "DMR Enrichment by Chromatin State",
         subtitle = "log2(Observed/Expected based on genomic coverage)",
         x = "Chromatin State", y = "log2(Fold Enrichment)") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave(file.path(PLOTDIR, "dmr_enrichment_chromatin_state.png"), p2,
       width = 12, height = 7, dpi = 300)
cat("Saved: dmr_enrichment_chromatin_state.png\n")

# ============================================================================
# Visualization 3: logFC distribution by chromatin state
# ============================================================================

# Order states for boxplot
dmr_annot$chromatin_state <- factor(dmr_annot$chromatin_state,
                                     levels = state_order)

p3 <- ggplot(dmr_annot, aes(x = chromatin_state, y = logFC)) +
    geom_boxplot(fill = "#9ECAE1", outlier.size = 0.5, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Methylation Change by Chromatin State",
         x = "Chromatin State", y = "log2(Fold Change)") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(PLOTDIR, "logFC_by_chromatin_state.png"), p3,
       width = 12, height = 7, dpi = 300)
cat("Saved: logFC_by_chromatin_state.png\n")

# ============================================================================
# Visualization 4: Category-level summary (grouped states)
# ============================================================================

# Add category to DMR annotation
dmr_with_cat <- dmr_annot %>%
    left_join(state_categories %>% select(state, category),
              by = c("chromatin_state" = "state"))

# Handle unknowns
dmr_with_cat$category[is.na(dmr_with_cat$category)] <- "Unknown"

cat_summary <- dmr_with_cat %>%
    group_by(category) %>%
    summarise(
        n_dmrs = n(),
        mean_logFC = mean(logFC, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(pct = n_dmrs / sum(n_dmrs) * 100)

# Define category colors
cat_colors <- c("Active_TSS" = "#FF4136", "Transcription" = "#2ECC40",
                "Enhancer" = "#FFDC00", "ZNF_Repeats" = "#66CDAA",
                "Heterochromatin" = "#8A91D0", "Bivalent" = "#FF851B",
                "Repressed" = "#AAAAAA", "Quiescent" = "#DDDDDD",
                "Unknown" = "#000000")

p4 <- ggplot(cat_summary, aes(x = reorder(category, -n_dmrs), y = n_dmrs, fill = category)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = cat_colors, guide = "none") +
    geom_text(aes(label = paste0(n_dmrs, "\n(", round(pct, 1), "%)")),
              vjust = -0.2, size = 3) +
    labs(title = "DMR Distribution by Chromatin Category",
         x = "Chromatin Category", y = "Number of DMRs") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(PLOTDIR, "dmr_by_chromatin_category.png"), p4,
       width = 10, height = 7, dpi = 300)
cat("Saved: dmr_by_chromatin_category.png\n")

# ============================================================================
# Statistical tests
# ============================================================================

# Kruskal-Wallis test
kw_test <- kruskal.test(logFC ~ chromatin_state, data = dmr_annot)
cat("\nKruskal-Wallis test (logFC ~ chromatin_state):\n")
cat("  Chi-squared =", round(kw_test$statistic, 2), "\n")
cat("  p-value =", format(kw_test$p.value, digits = 3, scientific = TRUE), "\n")

# Save statistical results
sink(file.path(OUTDIR, "statistical_tests.txt"))
cat("Chromatin State Analysis - Statistical Tests\n")
cat("=============================================\n\n")
print(kw_test)
sink()

cat("\nAnalysis complete!\n")
cat("Output directory:", OUTDIR, "\n")
RSCRIPT

# Run R analysis
conda activate r_chipseq_env
Rscript "$OUTDIR/analyze_chromatin_states.R"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Chromatin State Analysis Complete"
echo "=============================================="
echo ""
echo "Reference used: $REF_NAME"
echo "Output directory: $OUTDIR/"
echo ""
echo "Key outputs:"
echo "  - dmr_chromatin_annotation.tsv: DMRs with chromatin state"
echo "  - chromatin_state_summary.csv: Summary statistics"
echo "  - statistical_tests.txt: Kruskal-Wallis test results"
echo ""
echo "Plots:"
ls -1 "$OUTDIR/plots/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "IMPORTANT CAVEAT:"
echo "  ChromHMM states are from Roadmap Epigenomics reference ($REF_NAME),"
echo "  not SNB19 glioblastoma cells. Results should be interpreted with"
echo "  this limitation in mind."
echo ""
echo "Finished: $(date)"
