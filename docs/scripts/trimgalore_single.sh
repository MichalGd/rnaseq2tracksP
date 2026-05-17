#!/usr/bin/env bash
# ORIGIN: ADAPTED v1 — parameterised TrimGalore SE/PE wrapper
# Usage: trimgalore_single.sh <R1> <R2_or_empty> <outdir> <quality> <min_len> <layout> <sample_id>
set -euo pipefail
R1="$1"; R2="$2"; OUTDIR="$3"; QUAL="${4:-20}"; MINLEN="${5:-20}"
LAYOUT="$6"; SID="$7"
if [[ "$LAYOUT" == "PE" ]]; then
  trim_galore --paired --quality "$QUAL" --length "$MINLEN" \
    --fastqc --output_dir "$OUTDIR" --basename "$SID" "$R1" "$R2"
else
  trim_galore --quality "$QUAL" --length "$MINLEN" \
    --fastqc --output_dir "$OUTDIR" "$R1"
fi
