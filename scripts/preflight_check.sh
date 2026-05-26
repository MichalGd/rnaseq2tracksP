#!/usr/bin/env bash
# ORIGIN: NEW v4.0 — preflight dependency and file check
set -euo pipefail
CONFIG="$1"

# ── Sanitize CRLF before sourcing ─────────────────────────────────────────────
_sanitize_crlf() {
    local file="$1" label="$2"
    [[ -f "$file" ]] || return
    if grep -qP '\r' "$file" 2>/dev/null; then
        echo "  [FIX ] CRLF detected in $label — stripping \\r: $file"
        sed -i 's/\r//' "$file"
    fi
}
_sanitize_crlf "$CONFIG" "config"
# SAMPLESHEET and CONTRASTS are sourced from config, so source first then sanitize
source "$CONFIG"
[[ -n "${SAMPLESHEET:-}" ]] && _sanitize_crlf "$SAMPLESHEET" "samplesheet"
[[ -n "${CONTRASTS:-}"   ]] && _sanitize_crlf "$CONTRASTS"   "contrasts"

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
[[ "${SAMPLESHEET:-}" != /* ]] && warn "SAMPLESHEET is a relative path — use absolute paths in config"
[[ "${CONTRASTS:-}"   != /* ]] && warn "CONTRASTS is a relative path — use absolute paths in config"
SS="${SAMPLESHEET:-config/samplesheet.csv}"
if [[ -f "$SS" ]]; then
  N=$(grep -vc '^[[:space:]]*#\|^sample_id' "$SS" || true)
  ok "Samplesheet: $N samples"
else fail "Samplesheet not found: $SS"; fi

section "8. Input file format"
# Required columns
_check_header() {
    local file="$1" label="$2"; shift 2
    [[ -f "$file" ]] || { fail "$label not found: $file"; return; }
    local hdr; hdr=$(head -1 "$file")
    local col; for col in "$@"; do
        echo "$hdr" | grep -qF "$col" || fail "$label missing column: '$col'"
    done
    echo "$hdr" | grep -qF "$1" && ok "$label columns OK"
}
_check_header "${SAMPLESHEET:-}" "samplesheet" \
    sample_id fastq_R1 condition replicate strandedness
[[ -n "${CONTRASTS:-}" ]] && _check_header "${CONTRASTS:-}" "contrasts" \
    contrast_id numerator denominator

# Contrast levels must exist in samplesheet conditions
if [[ -n "${CONTRASTS:-}" && -f "${CONTRASTS:-}" && -f "${SAMPLESHEET:-}" ]]; then
    python3 - <<PYCHECK
import csv, sys
conds = set()
with open("${SAMPLESHEET}") as f:
    for r in csv.DictReader(f): conds.add(r["condition"].strip())
errs = []
with open("${CONTRASTS}") as f:
    for r in csv.DictReader(f):
        for side in ("numerator","denominator"):
            v = r[side].strip()
            if v not in conds:
                errs.append(f"contrast '{r['contrast_id']}': {side} '{v}' not in samplesheet")
if errs:
    for e in errs: print(f"  [FAIL] {e}")
    sys.exit(1)
else:
    print("  [ OK ] All contrast levels match samplesheet conditions")
PYCHECK
    [[ $? -ne 0 ]] && FAIL=$((FAIL+1))
fi

echo ""
echo "════════════════════════════════════════"
echo "Preflight: $FAIL FAIL  $WARN WARN"
echo "════════════════════════════════════════"
[[ $FAIL -eq 0 ]] || { echo "Fix FAIL items before running." >&2; exit 1; }
echo "Preflight passed."

# ── FastQ Screen (Step 2b — optional) ────────────────────────────────────────
FASTQSCREEN_CONF="${FASTQSCREEN_CONF:-$REPO/config/fastq_screen.conf}"
if command -v fastq_screen &>/dev/null; then
    log "OK  fastq_screen: $(fastq_screen --version 2>&1 | head -1)"
    if [[ -f "$FASTQSCREEN_CONF" ]]; then
        log "OK  fastq_screen.conf: $FASTQSCREEN_CONF"
        # verify all DATABASE paths resolve
        while IFS= read -r line; do
            [[ "$line" =~ ^DATABASE ]] || continue
            db_path=$(echo "$line" | awk '{print $3}')
            if ls "${db_path}".1.bt2 &>/dev/null || ls "${db_path}".1.bt2l &>/dev/null; then
                log "OK  bowtie2 index: $db_path"
            else
                log "WARN missing bowtie2 index: $db_path (Step 2b will skip this database)"
            fi
        done < "$FASTQSCREEN_CONF"
    else
        log "WARN fastq_screen.conf not found: $FASTQSCREEN_CONF (Step 2b will be skipped)"
    fi
else
    log "WARN fastq_screen not in PATH (Step 2b will be skipped)"
fi
