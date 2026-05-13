#!/usr/bin/env bash
# ORIGIN: ADAPTED — RNA-seq/scripts/ShellScripts/normBedGrpaphToBigWigM.sh
# Changes: all paths parameterized; accepts gzip; cleanup via trap.
# Usage: norm_bedgraph_to_bigwig.sh <bedGraph[.gz]> <chrom_sizes> <outdir> <kentutils_dir>
set -euo pipefail
mkdir -p "$3"
BASE="$(basename "$1")"; STEM="${BASE%.gz}"; STEM="${STEM%.bedGraph}"
SORTED="$3/${STEM}S.bedGraph"
trap 'rm -f "$SORTED"' EXIT
if [[ "$1" == *.gz ]]; then gunzip -c "$1"; else cat "$1"; fi \
  | LC_COLLATE=C sort -k1,1 -k2,2n > "$SORTED"
"$4/bedGraphToBigWig" "$SORTED" "$2" "$3/${STEM}S.bw"
gzip "$SORTED"; trap - EXIT
