#!/usr/bin/env bash
# ORIGIN: NEW v3 / UPDATED v4 — 7 pre-run checks including RSeQC
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${1:-$REPO/config/config.conf}"
PASS=0; FAIL=0; WARN=0
ok()   { echo "  [PASS] $*"; PASS=$((PASS+1)); }
fail() { echo "  [FAIL] $*"; FAIL=$((FAIL+1)); }
warn() { echo "  [WARN] $*"; WARN=$((WARN+1)); }
section() { echo ""; echo "── $* ─────────────────────────────────────"; }

section "1. Bash syntax"
while IFS= read -r -d '' f; do
  bash -n "$f" 2>/dev/null && ok "$(basename "$f")" || fail "$(basename "$f")"
done < <(find "$REPO/scripts" -name "*.sh" -print0)

section "2. R packages"
for pkg in DESeq2 apeglm Rsamtools GenomicAlignments rtracklayer \
           GenomicFeatures GenomicRanges vsn pheatmap RColorBrewer \
           ggplot2 data.table optparse knitr kableExtra rmarkdown; do
  Rscript -e "library($pkg,quietly=TRUE)" 2>/dev/null && ok "R: $pkg" || fail "R: $pkg"
done

section "3. Core tools"
for t in STAR samtools bedtools fastqc trim_galore multiqc Rscript; do
  command -v "$t" &>/dev/null && ok "$t" || warn "$t not found"
done

section "4. Kent utils"
[[ -f "$CONFIG" ]] && source "$CONFIG" 2>/dev/null || true
[[ -x "${KENTUTILS_DIR:-}/bedGraphToBigWig" ]] \
  && ok "bedGraphToBigWig" || warn "bedGraphToBigWig not found"

section "5. RSeQC"
RSEQC_DIR="${RSEQC_BIN_DIR:-}"
for py in infer_experiment.py read_distribution.py geneBody_coverage.py \
          junction_annotation.py junction_saturation.py; do
  if [[ -n "$RSEQC_DIR" && -x "$RSEQC_DIR/$py" ]]; then ok "RSeQC: $py"
  elif command -v "$py" &>/dev/null; then ok "RSeQC: $py"
  else warn "RSeQC: $py not found"; fi
done

section "6. Config"
if [[ -f "$CONFIG" ]]; then
  for v in SPECIES LIBRARY_LAYOUT SAMPLESHEET OUTDIR RUN_RSEQC \
            REGULAR_CHROMS_ONLY CHROMOSOME_NAMING STRAND_TOLERANCE_PCT; do
    [[ -n "${!v:-}" ]] && ok "$v=${!v}" || warn "$v empty"
  done
  case "${SPECIES:-}" in
    human) for v in STAR_INDEX_HUMAN GTF_HUMAN CHROM_SIZES_HUMAN RSEQC_BED_HUMAN; do
      [[ -n "${!v:-}" ]] && ok "$v" || warn "$v empty"; done ;;
    mouse) for v in STAR_INDEX_MOUSE GTF_MOUSE CHROM_SIZES_MOUSE RSEQC_BED_MOUSE; do
      [[ -n "${!v:-}" ]] && ok "$v" || warn "$v empty"; done ;;
    *) warn "SPECIES not set" ;;
  esac
else warn "config.conf not found at $CONFIG"; fi

section "7. Samplesheet"
SS="${SAMPLESHEET:-$REPO/config/samplesheet.csv}"
if [[ -f "$SS" ]]; then
  N=$(grep -vc '^[[:space:]]*#\|^sample_id' "$SS" || true)
  ok "$N data rows in $SS"
else warn "Samplesheet not found: $SS"; fi

echo ""; echo "════════════════════════════════════════"
echo "Smoke test: $PASS passed  $FAIL failed  $WARN warnings"
echo "════════════════════════════════════════"
[[ $FAIL -eq 0 ]] || { echo "Fix FAIL items before running."; exit 1; }
