#!/usr/bin/env bash
# ORIGIN: NEW
# Usage: bam_sort_index.sh <input_bam> <outdir> <threads>
set -euo pipefail
mkdir -p "$2"
STEM="$(basename "$1" .bam)"; STEM="${STEM%_Aligned.out}"
samtools sort -@ "${3:-4}" -o "$2/${STEM}_sortedS.bam" "$1"
samtools index "$2/${STEM}_sortedS.bam"
