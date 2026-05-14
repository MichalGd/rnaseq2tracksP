#!/usr/bin/env bash
# ORIGIN: NEW v1
# Usage: bam_sort_index.sh <input.bam> <outdir> <threads>
set -euo pipefail
BAM="$1"; OUTDIR="$2"; THR="${3:-4}"
BASE="$(basename "$BAM" _Aligned.out.bam)"
SORTED="${OUTDIR}/${BASE}_sortedS.bam"
samtools sort -@ "$THR" -o "$SORTED" "$BAM"
samtools index "$SORTED"
