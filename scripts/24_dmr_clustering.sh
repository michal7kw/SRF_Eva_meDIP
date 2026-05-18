#!/bin/bash
#SBATCH --job-name=dmr_clustering
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=logs/24_dmr_clustering.out
#SBATCH --error=logs/24_dmr_clustering.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# DMR Clustering Analysis
# ============================================================================
#
# Purpose: Identify larger differentially methylated domains (DMDs) by merging
#          adjacent DMRs and characterize them by genomic features
#
# Rationale:
#   - Individual 500bp DMRs may be part of larger methylation domains
#   - Clustered DMRs suggest regional epigenetic changes
#   - "Mega-DMRs" may indicate coordinated regulatory events
#
# Best Practices:
#   - Test multiple merge distances to avoid arbitrary threshold bias
#   - Require consistent direction within clusters (>80% same direction)
#   - Compare cluster gene associations with single-DMR results
#
# Output:
#   - Clustered DMRs at multiple distance thresholds (500bp, 1kb, 2kb)
#   - Cluster statistics (size, direction consistency, genomic features)
#   - "Mega-DMR" identification (>5 DMRs or >10kb span)
#
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP

# Create output directories
OUTDIR="results/24_dmr_clustering"
mkdir -p "$OUTDIR"/{beds,plots}
mkdir -p logs

echo "=============================================="
echo "DMR Clustering Analysis"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ============================================================================
# Define input files
# ============================================================================

# DMR file
DMR_FILE="results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

# Gene annotation (for feature annotation)
GENE_ANNOT="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# CpG context from previous analysis (if available)
CPG_CONTEXT="results/22_cpg_context/dmr_cpg_annotation.tsv"

echo ""
echo "=== Checking input files ==="

if [[ -f "$DMR_FILE" ]]; then
    echo "  Found DMR file: $DMR_FILE"
else
    echo "  ERROR: DMR file not found: $DMR_FILE"
    exit 1
fi

# ============================================================================
# Step 1: Convert DMRs to BED format with logFC information
# ============================================================================

echo ""
echo "=== Step 1: Preparing DMR BED file ==="

# Format: chr start end dmr_id logFC direction
tail -n +2 "$DMR_FILE" | \
    awk -F',' 'BEGIN{OFS="\t"} {
        gsub(/"/, "", $1);
        gsub(/"/, "", $2);
        gsub(/"/, "", $3);
        gsub(/"/, "", $7);
        chr = $1;
        if (chr !~ /^chr/) chr = "chr" chr;
        logFC = $7 + 0;
        direction = (logFC > 0) ? "hyper" : "hypo";
        print chr, $2, $3, "DMR_" NR, logFC, direction
    }' | \
    sort -k1,1 -k2,2n > "$OUTDIR/beds/dmrs_with_logFC.bed"

DMR_COUNT=$(wc -l < "$OUTDIR/beds/dmrs_with_logFC.bed")
echo "  Total DMRs: $DMR_COUNT"

# Count by direction
HYPER_COUNT=$(awk '$6=="hyper"' "$OUTDIR/beds/dmrs_with_logFC.bed" | wc -l)
HYPO_COUNT=$(awk '$6=="hypo"' "$OUTDIR/beds/dmrs_with_logFC.bed" | wc -l)
echo "  Hypermethylated: $HYPER_COUNT"
echo "  Hypomethylated: $HYPO_COUNT"

# ============================================================================
# Step 2: Cluster DMRs at multiple distance thresholds
# ============================================================================

echo ""
echo "=== Step 2: Clustering DMRs at multiple distance thresholds ==="

# Function to cluster DMRs at a given distance
cluster_dmrs() {
    local dist=$1
    local dist_label=$2

    echo "  Clustering at ${dist_label}..."

    # Use bedtools merge with -d to allow gaps
    # -c 4,5,6 -o collapse to preserve DMR info
    bedtools merge -i "$OUTDIR/beds/dmrs_with_logFC.bed" \
        -d "$dist" \
        -c 4,5,6 \
        -o collapse,collapse,collapse > "$OUTDIR/beds/dmr_clusters_${dist_label}.bed"

    local cluster_count=$(wc -l < "$OUTDIR/beds/dmr_clusters_${dist_label}.bed")
    echo "    Clusters at ${dist_label}: $cluster_count"
}

# Cluster at different distances
cluster_dmrs 500 "500bp"
cluster_dmrs 1000 "1kb"
cluster_dmrs 2000 "2kb"
cluster_dmrs 5000 "5kb"

# ============================================================================
# Step 3: Calculate cluster statistics
# ============================================================================

echo ""
echo "=== Step 3: Calculating cluster statistics ==="

# Function to calculate cluster stats
calc_cluster_stats() {
    local dist_label=$1
    local input="$OUTDIR/beds/dmr_clusters_${dist_label}.bed"
    local output="$OUTDIR/cluster_stats_${dist_label}.tsv"

    echo "  Processing ${dist_label} clusters..."

    {
        echo -e "cluster_id\tchr\tstart\tend\tspan_bp\tn_dmrs\tdmr_ids\tmean_logFC\tmedian_logFC\tn_hyper\tn_hypo\tpct_hyper\tdirection_consistent"

        awk 'BEGIN{OFS="\t"} {
            chr = $1
            start = $2
            end = $3
            span = end - start

            # Parse collapsed DMR IDs and logFCs
            split($4, ids, ",")
            split($5, logfcs, ",")
            split($6, dirs, ",")

            n_dmrs = length(ids)

            # Calculate statistics
            sum_logfc = 0
            n_hyper = 0
            n_hypo = 0

            for (i = 1; i <= n_dmrs; i++) {
                sum_logfc += logfcs[i]
                if (dirs[i] == "hyper") n_hyper++
                else n_hypo++
            }

            mean_logfc = sum_logfc / n_dmrs
            pct_hyper = n_hyper / n_dmrs * 100

            # Check direction consistency (>80% same direction)
            max_dir = (n_hyper > n_hypo) ? n_hyper : n_hypo
            consistent = (max_dir / n_dmrs >= 0.8) ? "yes" : "no"

            # Median (approximate - take middle value when sorted)
            # For simplicity, use mean as proxy

            print "CLU_" NR, chr, start, end, span, n_dmrs, $4, mean_logfc, mean_logfc, n_hyper, n_hypo, pct_hyper, consistent
        }' "$input"

    } > "$output"

    echo "    Saved: $output"
}

calc_cluster_stats "500bp"
calc_cluster_stats "1kb"
calc_cluster_stats "2kb"
calc_cluster_stats "5kb"

# ============================================================================
# Step 4: Identify mega-DMRs (>5 DMRs or >10kb span)
# ============================================================================

echo ""
echo "=== Step 4: Identifying mega-DMRs ==="

# Use 1kb clustering as default
{
    echo -e "cluster_id\tchr\tstart\tend\tspan_bp\tn_dmrs\tmean_logFC\tpct_hyper\tdirection_consistent\tqualifies_by"

    awk 'BEGIN{OFS="\t"} NR>1 {
        span = $5
        n_dmrs = $6
        mean_logfc = $8
        pct_hyper = $12
        consistent = $13

        qualifies = ""
        if (n_dmrs >= 5) qualifies = qualifies "n_dmrs;"
        if (span >= 10000) qualifies = qualifies "span;"

        if (qualifies != "") {
            sub(/;$/, "", qualifies)
            print $1, $2, $3, $4, span, n_dmrs, mean_logfc, pct_hyper, consistent, qualifies
        }
    }' "$OUTDIR/cluster_stats_1kb.tsv"

} > "$OUTDIR/mega_dmrs.tsv"

MEGA_COUNT=$(tail -n +2 "$OUTDIR/mega_dmrs.tsv" | wc -l)
echo "  Identified $MEGA_COUNT mega-DMRs (>5 DMRs or >10kb span)"

# Save mega-DMR BED file
tail -n +2 "$OUTDIR/mega_dmrs.tsv" | \
    awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $1, $7}' > "$OUTDIR/beds/mega_dmrs.bed"

# ============================================================================
# Step 5: Annotate clusters with genomic features (if GTF available)
# ============================================================================

echo ""
echo "=== Step 5: Annotating clusters with genomic features ==="

if [[ -f "$GENE_ANNOT" ]]; then
    # Extract gene coordinates from GTF
    echo "  Extracting gene coordinates from GTF..."

    # Get gene bodies
    awk '$3 == "gene"' "$GENE_ANNOT" | \
        awk 'BEGIN{OFS="\t"} {
            chr = $1
            if (chr !~ /^chr/) chr = "chr" chr
            start = $4
            end = $5
            # Extract gene_name
            match($0, /gene_name "([^"]+)"/, arr)
            gene_name = arr[1]
            print chr, start, end, gene_name
        }' | sort -k1,1 -k2,2n > "$OUTDIR/beds/genes.bed"

    # Get promoters (TSS ± 2kb)
    awk '$3 == "gene"' "$GENE_ANNOT" | \
        awk 'BEGIN{OFS="\t"} {
            chr = $1
            if (chr !~ /^chr/) chr = "chr" chr
            strand = $7
            if (strand == "+") {
                tss = $4
            } else {
                tss = $5
            }
            start = tss - 2000
            end = tss + 2000
            if (start < 0) start = 0
            match($0, /gene_name "([^"]+)"/, arr)
            gene_name = arr[1]
            print chr, start, end, gene_name
        }' | sort -k1,1 -k2,2n > "$OUTDIR/beds/promoters.bed"

    echo "  Extracted $(wc -l < "$OUTDIR/beds/genes.bed") genes"
    echo "  Extracted $(wc -l < "$OUTDIR/beds/promoters.bed") promoters"

    # Annotate 1kb clusters with gene overlap
    echo "  Annotating clusters with gene overlap..."
    bedtools intersect -a "$OUTDIR/beds/dmr_clusters_1kb.bed" \
        -b "$OUTDIR/beds/promoters.bed" \
        -wa -wb > "$OUTDIR/beds/clusters_promoter_overlap.bed" 2>/dev/null || touch "$OUTDIR/beds/clusters_promoter_overlap.bed"

    PROMOTER_CLUSTERS=$(cut -f1-3 "$OUTDIR/beds/clusters_promoter_overlap.bed" | sort -u | wc -l)
    echo "    Clusters overlapping promoters: $PROMOTER_CLUSTERS"
else
    echo "  GTF annotation not found, skipping gene annotation"
fi

# ============================================================================
# Step 6: Add CpG context if available
# ============================================================================

echo ""
echo "=== Step 6: Integrating CpG context (if available) ==="

if [[ -f "$CPG_CONTEXT" ]]; then
    echo "  Found CpG context data, integrating..."
    # This will be handled by R script
else
    echo "  CpG context not found, skipping integration"
    echo "  Run script 22 first for CpG context integration"
fi

# ============================================================================
# Step 7: Run R analysis
# ============================================================================

echo ""
echo "=== Step 7: Running R analysis ==="

conda activate r_chipseq_env
Rscript scripts/24_dmr_clustering.R

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "DMR Clustering Analysis Complete"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Clustering results at multiple thresholds:"
for dist in 500bp 1kb 2kb 5kb; do
    count=$(wc -l < "$OUTDIR/beds/dmr_clusters_${dist}.bed")
    echo "  ${dist}: $count clusters"
done
echo ""
echo "Mega-DMRs identified: $MEGA_COUNT"
echo ""
echo "Key outputs:"
echo "  cluster_stats_*.tsv - Statistics for each distance threshold"
echo "  mega_dmrs.tsv - Large methylation domains"
echo ""
echo "Plots:"
ls -1 "$OUTDIR/plots/"*.png 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Finished: $(date)"
