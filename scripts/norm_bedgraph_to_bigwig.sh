#!/usr/bin/env bash
# ORIGIN: ADAPTED v1→v3 — chr filter using SPECIES, CHROMOSOME_NAMING, REGULAR_CHROMS_ONLY
# Usage: norm_bedgraph_to_bigwig.sh <input.bedGraph.gz_or_.bedGraph> <chrom.sizes> <outdir> <kentutils_dir>
# Env vars (exported by master): SPECIES CHROMOSOME_NAMING REGULAR_CHROMS_ONLY
set -euo pipefail
INPUT="$1"; CS="$2"; OUTDIR="$3"; KENT="$4"
mkdir -p "$OUTDIR"
BASE="$(basename "$INPUT" .gz)"; BASE="${BASE%.bedGraph}"

# Decompress if needed
if [[ "$INPUT" == *.gz ]]; then
  DECOFILE="$OUTDIR/${BASE}.bedGraph"
  zcat "$INPUT" > "$DECOFILE"
else
  DECOFILE="$INPUT"
fi

SORTED="$OUTDIR/${BASE}.sorted.bedGraph"
sort -k1,1 -k2,2n "$DECOFILE" > "$SORTED"
[[ "$INPUT" == *.gz ]] && rm -f "$DECOFILE"

FINAL="$SORTED"
SP="${SPECIES:-mouse}"; NAMING="${CHROMOSOME_NAMING:-ucsc}"; FILTER="${REGULAR_CHROMS_ONLY:-true}"

chr_pattern() {
  local sp="$1" naming="$2"
  case "${sp}:${naming}" in
    human:ucsc)    echo '^(chr([1-9]|1[0-9]|2[0-2])|chrX|chrY|chrM)[[:space:]]' ;;
    human:ensembl) echo '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)[[:space:]]' ;;
    mouse:ucsc)    echo '^(chr([1-9]|1[0-9])|chrX|chrY|chrM)[[:space:]]' ;;
    mouse:ensembl) echo '^([1-9]|1[0-9]|X|Y|MT)[[:space:]]' ;;
    *) echo '' ;;
  esac
}

if [[ "$FILTER" == "true" ]]; then
  ALLCHR="$OUTDIR/${BASE}.all_chromosomes.bedGraph.gz"
  gzip -c "$SORTED" > "$ALLCHR"
  PAT="$(chr_pattern "$SP" "$NAMING")"
  if [[ -n "$PAT" ]]; then
    FILTERED="$OUTDIR/${BASE}.filtered.bedGraph"
    grep -E "$PAT" "$SORTED" > "$FILTERED"
    FINAL="$FILTERED"
  fi
fi

BW="$OUTDIR/${BASE}.bw"
"${KENT}/bedGraphToBigWig" "$FINAL" "$CS" "$BW"
gzip -f "$SORTED"
[[ -f "$OUTDIR/${BASE}.filtered.bedGraph" ]] && rm -f "$OUTDIR/${BASE}.filtered.bedGraph"
echo "[norm_bedgraph_to_bigwig.sh] Written: $BW"
