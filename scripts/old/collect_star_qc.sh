#!/usr/bin/env bash
# ORIGIN: NEW v4.0 — parse STAR Log.final.out into TSV + symlinks for MultiQC
# Usage: collect_star_qc.sh <star_logs_dir> <outdir>
set -euo pipefail
LOGDIR="$1"; OUTDIR="$2/star"
mkdir -p "$OUTDIR"
for f in "$LOGDIR"/*Log.final.out; do
  [[ -f "$f" ]] || continue
  ln -sf "$(realpath "$f")" "$OUTDIR/$(basename "$f")"
done
TSV="$OUTDIR/star_alignment_summary.tsv"
printf "sample\tuniquely_mapped_pct\tmulti_mapped_pct\ttoo_short_pct\tnum_input_reads\tuniquely_mapped_reads\tavg_mapped_length\n" > "$TSV"
for f in "$OUTDIR"/*Log.final.out; do
  [[ -f "$f" ]] || continue
  S="$(basename "$f" | sed 's/_Log\.final\.out$//')"
  ex() { grep "$1" "$f" | awk -F'|' '{gsub(/[ \t%]/,"",$2);print $2}'; }
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$S" "$(ex 'Uniquely mapped reads %')" \
    "$(ex '% of reads mapped to multiple loci')" \
    "$(ex '% of reads unmapped: too short')" \
    "$(ex 'Number of input reads')" \
    "$(ex 'Uniquely mapped reads number')" \
    "$(ex 'Average mapped length')" >> "$TSV"
done
echo "[collect_star_qc.sh] Written: $TSV"
