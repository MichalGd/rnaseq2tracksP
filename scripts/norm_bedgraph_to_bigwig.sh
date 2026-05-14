#!/usr/bin/env bash
# =============================================================================
# norm_bedgraph_to_bigwig.sh
# =============================================================================
# ORIGIN: ADAPTED — RNA-seq/scripts/ShellScripts/normBedGrpaphToBigWigM.sh
# v3 changes:
#   - CHROMOSOME_NAMING=ucsc|ensembl  (both human and mouse)
#   - REGULAR_CHROMS_ONLY=true|false  (escape hatch for custom genomes)
#   - keeps <stem>.all_chromosomes.bedGraph for debugging when filtering is on
#   - accepts SPECIES from environment (set by master script)
#
# Usage (called by master with env vars pre-set):
#   norm_bedgraph_to_bigwig.sh <bedGraph[.gz]> <chrom_sizes> <outdir> <kentutils_dir>
# =============================================================================
set -euo pipefail
BG="$1"; CHROM="$2"; OUTDIR="$3"; KENT="$4"
mkdir -p "$OUTDIR"

BASE="$(basename "$BG")"; STEM="${BASE%.gz}"; STEM="${STEM%.bedGraph}"
ALL_BG="$OUTDIR/${STEM}.all_chromosomes.bedGraph"
SORTED="$OUTDIR/${STEM}S.bedGraph"
trap 'rm -f "$SORTED"' EXIT

# Decompress
if [[ "$BG" == *.gz ]]; then gunzip -c "$BG" > "$ALL_BG"
else cp "$BG" "$ALL_BG"; fi

# Build chromosome filter pattern based on CHROMOSOME_NAMING and SPECIES
SPECIES="${SPECIES:-mouse}"
CHR_NAMING="${CHROMOSOME_NAMING:-ucsc}"
REGULAR="${REGULAR_CHROMS_ONLY:-true}"

chr_pattern() {
  local sp="$1" nm="$2"
  if   [[ "$nm" == "ucsc"    && "$sp" == "human" ]]; then echo '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\t'
  elif [[ "$nm" == "ucsc"    && "$sp" == "mouse" ]]; then echo '^chr([1-9]|1[0-9]|X|Y|M)\t'
  elif [[ "$nm" == "ensembl" && "$sp" == "human" ]]; then echo '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)\t'
  elif [[ "$nm" == "ensembl" && "$sp" == "mouse" ]]; then echo '^([1-9]|1[0-9]|X|Y|MT)\t'
  else echo "ERROR: unknown CHROMOSOME_NAMING=$nm / SPECIES=$sp" >&2; exit 1; fi
}

if [[ "$REGULAR" == "true" ]]; then
  PATTERN="$(chr_pattern "$SPECIES" "$CHR_NAMING")"
  grep -E "$PATTERN" "$ALL_BG" \
    | LC_COLLATE=C sort -k1,1 -k2,2n > "$SORTED"
  echo "  Filtered to canonical ${CHR_NAMING} ${SPECIES} chromosomes"
else
  LC_COLLATE=C sort -k1,1 -k2,2n "$ALL_BG" > "$SORTED"
  echo "  REGULAR_CHROMS_ONLY=false — all chromosomes retained"
fi

"${KENT}/bedGraphToBigWig" "$SORTED" "$CHROM" "$OUTDIR/${STEM}S.bw"
gzip -f "$SORTED"
trap - EXIT

# Keep the all-chromosomes bedGraph only when filtering was active (for debugging)
# Remove it if filtering is off (it is already the input = $SORTED)
[[ "$REGULAR" == "true" ]] && gzip -f "$ALL_BG" || rm -f "$ALL_BG"

echo "BigWig written: $OUTDIR/${STEM}S.bw"
