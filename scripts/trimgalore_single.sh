#!/usr/bin/env bash
# ORIGIN: ADAPTED — fastq2tracks/scripts/trimgalore_batch.1.0.sh
# Changes: single-sample wrapper; SE/PE mode via parameter.
# Usage: trimgalore_single.sh <R1> <R2_or_empty> <outdir> <qual> <minlen> <SE|PE> <sample_id>
set -euo pipefail
mkdir -p "$3"
if [[ "$6" == "PE" ]]; then
  trim_galore --paired --quality "$4" --length "$5" --basename "$7" -o "$3" "$1" "$2"
else
  trim_galore --quality "$4" --length "$5" --basename "$7" -o "$3" "$1"
fi
