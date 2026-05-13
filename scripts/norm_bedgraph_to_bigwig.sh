#!/usr/bin/env bash
# =============================================================================
# norm_bedgraph_to_bigwig.sh
# =============================================================================
# ORIGIN: ADAPTED — RNA-seq/scripts/ShellScripts/normBedGrpaphToBigWigM.sh
# v2 change: non-standard chromosomes/scaffolds filtered out before conversion
#   so BigWigs are clean for UCSC Genome Browser.
#   Keeps: chr1-chr22, chrX, chrY, chrM  (works for human and mouse)
#
# Usage: norm_bedgraph_to_bigwig.sh <bedGraph[.gz]> <chrom_sizes> <outdir> <kentutils_dir>
# =============================================================================
set -euo pipefail
BG="$1"; CHROM="$2"; OUTDIR="$3"; KENT="$4"
mkdir -p "$OUTDIR"
BASE="$(basename "$BG")"; STEM="${BASE%.gz}"; STEM="${STEM%.bedGraph}"
SORTED="$OUTDIR/${STEM}S.bedGraph"
trap 'rm -f "$SORTED"' EXIT

# Decompress, filter to standard chromosomes, sort
if [[ "$BG" == *.gz ]]; then gunzip -c "$BG"; else cat "$BG"; fi \
  | grep -E '^chr([0-9]+|X|Y|M)\t' \
  | LC_COLLATE=C sort -k1,1 -k2,2n > "$SORTED"

"${KENT}/bedGraphToBigWig" "$SORTED" "$CHROM" "$OUTDIR/${STEM}S.bw"
gzip "$SORTED"
trap - EXIT
echo "BigWig written: $OUTDIR/${STEM}S.bw"
