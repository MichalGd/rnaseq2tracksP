#!/usr/bin/env bash
# ORIGIN: NEW v4.0 — preflight dependency and file check
set -euo pipefail
CONFIG="$1"; source "$CONFIG"
FAIL=0; WARN=0
fail() { echo "  [FAIL] $*" >&2; FAIL=$((FAIL+1)); }
warn() { echo "  [WARN] $*"; WARN=$((WARN+1)); }
ok()   { echo "  [ OK ] $*"; }
section() { echo ""; echo "── $* ────────────────────────────────────────"; }

section "1. Core tools"
for t in STAR samtools bedtools fastqc trim_galore "${FASTQC_BIN:-fastqc}" \
         "${MULTIQC_BIN:-multiqc}" "${RSCRIPT_BIN:-Rscript}"; do
  command -v "$t" &>/dev/null && ok "$t" || fail "$t not in PATH"
done

section "2. bedGraphToBigWig"
[[ -x "${KENTUTILS_DIR:-}/bedGraphToBigWig" ]] \
  && ok "bedGraphToBigWig: ${KENTUTILS_DIR}/bedGraphToBigWig" \
  || fail "bedGraphToBigWig not found in KENTUTILS_DIR='${KENTUTILS_DIR:-}'"

section "3. R packages"
RSCRIPT="${RSCRIPT_BIN:-Rscript}"
for pkg in DESeq2 apeglm Rsamtools GenomicAlignments rtracklayer \
           GenomicFeatures txdbmaker GenomicRanges vsn pheatmap RColorBrewer \
           ggplot2 data.table optparse knitr kableExtra rmarkdown ashr \
           clusterProfiler enrichplot ReactomePA fgsea msigdbr; do
  "$RSCRIPT" -e "library($pkg,quietly=TRUE)" 2>/dev/null \
    && ok "R: $pkg" || fail "R pkg missing: $pkg"
done

section "4. RSeQC"
RSEQC_DIR="${RSEQC_BIN_DIR:-}"
for py in infer_experiment.py read_distribution.py geneBody_coverage.py \
          junction_annotation.py junction_saturation.py; do
  if [[ -n "$RSEQC_DIR" && -x "$RSEQC_DIR/$py" ]]; then ok "RSeQC: $py"
  elif command -v "$py" &>/dev/null; then ok "RSeQC: $py (PATH)"
  else warn "RSeQC: $py not found — module will be skipped"; fi
done

section "5. RSeQC BED"
RSEQC_BED="${RSEQC_BED_MOUSE:-}"
[[ "${SPECIES:-mouse}" == "human" ]] && RSEQC_BED="${RSEQC_BED_HUMAN:-}"
if [[ -n "$RSEQC_BED" && -f "$RSEQC_BED" ]]; then ok "RSEQC_BED: $RSEQC_BED"
elif [[ "${RUN_RSEQC:-true}" == "true" ]]; then
  fail "RSEQC_BED not set/missing (required when RUN_RSEQC=true)"
else warn "RSEQC_BED not set — RSeQC skipped"; fi

section "6. Genome files"
SP="${SPECIES:-mouse}"
case "$SP" in
  human) IDX="${STAR_INDEX_HUMAN:-}"; G="${GTF_HUMAN:-}"; CS="${CHROM_SIZES_HUMAN:-}" ;;
  mouse) IDX="${STAR_INDEX_MOUSE:-}"; G="${GTF_MOUSE:-}"; CS="${CHROM_SIZES_MOUSE:-}" ;;
  *) fail "SPECIES must be human|mouse"; IDX=""; G=""; CS="" ;;
esac
[[ -d "$IDX" ]] && ok "STAR index: $IDX" || fail "STAR index missing: $IDX"
[[ -f "$G"   ]] && ok "GTF: $G"          || fail "GTF missing: $G"
[[ -f "$CS"  ]] && ok "chrom.sizes: $CS" || fail "chrom.sizes missing: $CS"

section "7. Samplesheet"
SS="${SAMPLESHEET:-config/samplesheet.csv}"
if [[ -f "$SS" ]]; then
  N=$(grep -vc '^[[:space:]]*#\|^sample_id' "$SS" || true)
  ok "Samplesheet: $N samples"
else fail "Samplesheet not found: $SS"; fi

echo ""
echo "════════════════════════════════════════"
echo "Preflight: $FAIL FAIL  $WARN WARN"
echo "════════════════════════════════════════"
[[ $FAIL -eq 0 ]] || { echo "Fix FAIL items before running." >&2; exit 1; }
echo "Preflight passed."
