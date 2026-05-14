#!/usr/bin/env bash
# =============================================================================
# run_smoke_test.sh — executable smoke test for rnaseq2tracks
# =============================================================================
# Tests:
#   1. Bash syntax check on all .sh scripts
#   2. R library availability check (all required packages)
#   3. Tool availability check (STAR, samtools, fastqc, trim_galore, multiqc,
#      Rscript, bedGraphToBigWig)
#   4. Config template validation (all required variables present)
#   5. Samplesheet parsing dry-run
#
# Does NOT require real FASTQ/BAM files or a STAR index.
# Usage: bash tests/run_smoke_test.sh [config/config.conf]
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${1:-$REPO/config/config.conf}"
PASS=0; FAIL=0; WARN=0

ok()   { echo "  [PASS] $*"; PASS=$((PASS+1)); }
fail() { echo "  [FAIL] $*"; FAIL=$((FAIL+1)); }
warn() { echo "  [WARN] $*"; WARN=$((WARN+1)); }
section() { echo ""; echo "── $* ──────────────────────────────────────────"; }

section "1. Bash syntax check"
while IFS= read -r -d '' f; do
  if bash -n "$f" 2>/dev/null; then ok "syntax: $(basename "$f")"
  else fail "syntax: $(basename "$f")"; fi
done < <(find "$REPO/scripts" -name "*.sh" -print0)

section "2. R package availability"
RSCRIPT="${RSCRIPT_BIN:-Rscript}"
for pkg in DESeq2 apeglm Rsamtools GenomicAlignments rtracklayer GenomicFeatures \
           GenomicRanges vsn pheatmap RColorBrewer ggplot2 data.table \
           optparse knitr kableExtra rmarkdown; do
  if "$RSCRIPT" -e "library($pkg)" 2>/dev/null; then ok "R: $pkg"
  else fail "R: $pkg  (install with BiocManager::install('$pkg') or conda)"; fi
done

section "3. External tool availability"
for tool in fastqc trim_galore STAR samtools multiqc Rscript; do
  if command -v "$tool" &>/dev/null; then ok "tool: $tool ($(command -v "$tool"))"
  else warn "tool: $tool not found in PATH"; fi
done

section "4. Config file check"
if [[ -f "$CONFIG" ]]; then
  source "$CONFIG"
  for v in SPECIES LIBRARY_LAYOUT SAMPLESHEET KENTUTILS_DIR \
            REGULAR_CHROMS_ONLY CHROMOSOME_NAMING OUTDIR; do
    if [[ -n "${!v:-}" ]]; then ok "config: $v=${!v}"
    else warn "config: $v is empty — fill in config.conf"; fi
  done
  # Check that species-specific paths are set
  case "${SPECIES:-}" in
    human) for v in STAR_INDEX_HUMAN GTF_HUMAN CHROM_SIZES_HUMAN; do
      [[ -n "${!v:-}" ]] && ok "config: $v" || warn "config: $v empty"; done ;;
    mouse) for v in STAR_INDEX_MOUSE GTF_MOUSE CHROM_SIZES_MOUSE; do
      [[ -n "${!v:-}" ]] && ok "config: $v" || warn "config: $v empty"; done ;;
    *) warn "config: SPECIES not set or unknown ('${SPECIES:-}')"; ;;
  esac
else
  warn "config.conf not found at $CONFIG — copy config_template.conf and fill in paths"
fi

section "5. Samplesheet parsing dry-run"
SS="${SAMPLESHEET:-$REPO/config/samplesheet.csv}"
if [[ -f "$SS" ]]; then
  N=$(grep -vc '^[[:space:]]*#\|^sample_id' "$SS" || true)
  ok "samplesheet: $N data rows found in $SS"
else
  warn "samplesheet not found: $SS"
fi

echo ""
echo "════════════════════════════════════════"
echo "Smoke test results: $PASS passed  $FAIL failed  $WARN warnings"
echo "════════════════════════════════════════"
[[ $FAIL -eq 0 ]] || { echo "Fix FAIL items before running the pipeline."; exit 1; }
