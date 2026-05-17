#!/usr/bin/env bash
# =============================================================================
# rerun_deseq_rep12.sh — Re-run DESeq2 analysis on replicates 1 & 2 only
# Reuses: dds.RData from full run (subsets in R), existing workflow scripts
# Produces: all QC plots, volcano, MA, DE tables in a separate output subdir
#
# Usage: bash rerun_deseq_rep12.sh <config.conf>
# =============================================================================
set -euo pipefail
[[ $# -ne 1 ]] && { echo "Usage: $0 <config.conf>" >&2; exit 1; }

CONFIG="$(realpath "$1")"
source "$CONFIG"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RSCRIPT="${RSCRIPT_BIN:-Rscript}"

# ── Paths ─────────────────────────────────────────────────────────────────────
FULL_COUNTS="${OUTDIR}/analysis/counts"          # from original run
FULL_DDS="${FULL_COUNTS}/dds.RData"
OUTDIR_R12="${OUTDIR}/analysis_rep12"            # new subdir — original untouched
COUNTS_R12="${OUTDIR_R12}/counts"
DE_R12="${OUTDIR_R12}/DE"
FIGURES_R12="${OUTDIR_R12}/figures"

[[ -f "$FULL_DDS" ]] || { echo "ERROR: dds.RData not found: $FULL_DDS" >&2; exit 1; }
mkdir -p "$COUNTS_R12" "$DE_R12" "$FIGURES_R12"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Step 1: Subset dds to rep 1 & 2, re-estimate size factors, write outputs ──
log "Step 1 — Subsetting dds.RData to replicates 1 & 2..."
"$RSCRIPT" - << REOF
suppressPackageStartupMessages({ library(DESeq2); library(data.table) })

load("${FULL_DDS}")   # loads 'dds'

# Subset to replicates 1 and 2
keep <- colData(dds)\$replicate %in% c(1, 2)
cat(sprintf("  Keeping %d / %d samples\n", sum(keep), ncol(dds)))
cat("  Samples kept:", paste(colnames(dds)[keep], collapse=", "), "\n")

dds <- dds[, keep]
dds <- estimateSizeFactors(dds)   # re-estimate on subset

# Write size factors
SF     <- sizeFactors(dds)
counts_mat <- counts(dds)
total_mapped <- colSums(counts_mat)
rpm    <- sweep(counts_mat, 2, total_mapped/1e6, "/")
mean_rpm <- rowMeans(rpm)
SF_rpm <- SF * mean(mean_rpm[mean_rpm > 1])

fwrite(data.table(sample_id=colnames(dds), size_factor=SF, sf_rpm=SF_rpm),
       "${COUNTS_R12}/size_factors.tsv", sep="\t")

# Write raw and normalized count tables
raw <- as.data.frame(counts(dds))
fwrite(cbind(gene=rownames(raw), raw),
       "${COUNTS_R12}/raw_counts.tsv", sep="\t")

norm <- as.data.frame(counts(dds, normalized=TRUE))
fwrite(cbind(gene=rownames(norm), norm),
       "${COUNTS_R12}/normalized_counts.tsv", sep="\t")

# Save subsetted dds
save(dds, file="${COUNTS_R12}/dds.RData")
cat("[Step 1] Done. dds saved to ${COUNTS_R12}/dds.RData\n")
REOF

# ── Step 2: QC plots (PCA, clustering, heatmap) ───────────────────────────────
log "Step 2 — QC plots (PCA, clustering, heatmaps)..."
"$RSCRIPT" "$REPO/scripts/Rscripts/deseq2_qc_plots.R" \
  --countsrdata "${COUNTS_R12}/dds.RData" \
  --outdir      "${FIGURES_R12}"

# ── Step 3: DE analysis (volcano, MA, tables) ─────────────────────────────────
log "Step 3 — Differential expression analysis..."
[[ "${SPECIES:-mouse}" == "human" ]] && _GTF="$GTF_HUMAN" || _GTF="$GTF_MOUSE"
GTF="$_GTF" \
"$RSCRIPT" "$REPO/scripts/Rscripts/deseq2_de.R" \
  --countsrdata "${COUNTS_R12}/dds.RData" \
  --contrasts   "$(realpath "${CONTRASTS}")" \
  --outdir      "${DE_R12}" \
  --padj        "${PADJ_THRESHOLD:-0.05}" \
  --lfc         "${LFC_THRESHOLD:-1}"

log "Done. All rep1+2 outputs in: ${OUTDIR_R12}/"
echo ""
echo "  QC figures : ${FIGURES_R12}/"
echo "  DE results : ${DE_R12}/"
echo "  Count tables: ${COUNTS_R12}/"
