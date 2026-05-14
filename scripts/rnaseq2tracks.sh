#!/usr/bin/env bash
# =============================================================================
# rnaseq2tracks.sh — master RNA-seq workflow orchestrator (v3.0)
# =============================================================================
# Language rationale
#   Bash : subprocess fan-out, job throttling, STAR/samtools/TrimGalore/FastQC/MultiQC
#   R    : ALL analytical steps (coverage, normalization, DE, QC plots, report)
#
# v3 changes vs v2:
#   - SPECIES-driven genome path selection (human/mouse, one config for both)
#   - CHROMOSOME_NAMING (ucsc|ensembl) + REGULAR_CHROMS_ONLY passed to BigWig step
#   - norm_bedgraph_to_bigwig.sh receives SPECIES/CHROMOSOME_NAMING/REGULAR_CHROMS_ONLY
#     as environment variables
#
# Usage: ./scripts/rnaseq2tracks.sh config/config.conf
# =============================================================================
set -euo pipefail
[[ $# -ne 1 ]] && { echo "Usage: $0 <config_file>" >&2; exit 1; }
CONFIG="$(realpath "$1")"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Job throttle ─────────────────────────────────────────────────────────────
declare -a _PIDS=()
submit() {
  while [[ ${#_PIDS[@]} -ge ${MAX_JOBS:-8} ]]; do
    local live=()
    for p in "${_PIDS[@]}"; do kill -0 "$p" 2>/dev/null && live+=("$p"); done
    _PIDS=("${live[@]+"${live[@]}"}"); [[ ${#_PIDS[@]} -ge ${MAX_JOBS:-8} ]] && sleep 2
  done
  eval "$@" &; _PIDS+=($!)
}
wait_all() {
  local ok=0
  for p in "${_PIDS[@]+"${_PIDS[@]}"}"; do wait "$p" || ok=1; done
  _PIDS=()
  [[ $ok -eq 0 ]] || { log "ERROR: a background job failed"; exit 1; }
}

# ── Species-driven genome path resolution ────────────────────────────────────
SPECIES="${SPECIES:-mouse}"
case "$SPECIES" in
  human)
    STAR_INDEX="${STAR_INDEX_HUMAN:?ERROR: STAR_INDEX_HUMAN not set in $CONFIG}"
    GTF="${GTF_HUMAN:?ERROR: GTF_HUMAN not set in $CONFIG}"
    CHROM_SIZES="${CHROM_SIZES_HUMAN:?ERROR: CHROM_SIZES_HUMAN not set in $CONFIG}"
    ;;
  mouse)
    STAR_INDEX="${STAR_INDEX_MOUSE:?ERROR: STAR_INDEX_MOUSE not set in $CONFIG}"
    GTF="${GTF_MOUSE:?ERROR: GTF_MOUSE not set in $CONFIG}"
    CHROM_SIZES="${CHROM_SIZES_MOUSE:?ERROR: CHROM_SIZES_MOUSE not set in $CONFIG}"
    ;;
  *)
    echo "ERROR: SPECIES must be 'human' or 'mouse', got '$SPECIES'" >&2; exit 1 ;;
esac
log "Species: $SPECIES  Naming: ${CHROMOSOME_NAMING:-ucsc}  Regular chrs only: ${REGULAR_CHROMS_ONLY:-true}"

# Export for child scripts (norm_bedgraph_to_bigwig.sh reads these)
export SPECIES CHROMOSOME_NAMING="${CHROMOSOME_NAMING:-ucsc}" \
       REGULAR_CHROMS_ONLY="${REGULAR_CHROMS_ONLY:-true}"

# ── Validate required config vars ────────────────────────────────────────────
for v in STAR_INDEX GTF CHROM_SIZES GENOME_ASSEMBLY SAMPLESHEET KENTUTILS_DIR \
          LIBRARY_LAYOUT MAX_JOBS STAR_THREADS; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: $v not resolved in $CONFIG" >&2; exit 1; }
done
[[ "$LIBRARY_LAYOUT" =~ ^(SE|PE)$ ]] || \
  { echo "ERROR: LIBRARY_LAYOUT must be SE or PE" >&2; exit 1; }

# ── Step 1: Output folders ────────────────────────────────────────────────────
OUTDIR="${OUTDIR:-$(pwd)/rnaseq2tracks_output}"
log "STEP 1 — Output: $OUTDIR"
mkdir -p \
  "$OUTDIR/fastQC/raw"          "$OUTDIR/fastQC/trimmed" \
  "$OUTDIR/multiQC/raw"         "$OUTDIR/multiQC/trimmed" \
  "$OUTDIR/multiQC/alignments"  "$OUTDIR/multiQC/final" \
  "$OUTDIR/trimmedFastq"        "$OUTDIR/STARalignments" \
  "$OUTDIR/STARlogs"            "$OUTDIR/STARgeneCounts" \
  "$OUTDIR/bams"                "$OUTDIR/bedGraph/raw" \
  "$OUTDIR/bedGraph/normalized" "$OUTDIR/bedGraph/merged" \
  "$OUTDIR/bigwig"              "$OUTDIR/analysis/counts" \
  "$OUTDIR/analysis/DE"         "$OUTDIR/analysis/figures" \
  "$OUTDIR/reports"

# ── Parse samplesheet ─────────────────────────────────────────────────────────
declare -a SID R1 R2 COND REP STRAND
while IFS=',' read -r f1 f2 f3 f4 f5 f6 _rest; do
  [[ "$f1" =~ ^[[:space:]]*# || "$f1" == "sample_id" ]] && continue
  if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    SID+=("$f1"); R1+=("$f2"); R2+=("$f3"); COND+=("$f4"); REP+=("$f5"); STRAND+=("$f6")
  else
    SID+=("$f1"); R1+=("$f2"); R2+=(""); COND+=("$f3"); REP+=("$f4"); STRAND+=("$f5")
  fi
done < <(grep -v '^[[:space:]]*#' "$SAMPLESHEET")
N=${#SID[@]}; log "Loaded $N samples  layout=$LIBRARY_LAYOUT"

# ── Steps 2–3: FastQC / MultiQC raw ──────────────────────────────────────────
log "STEP 2 — FastQC raw"
for ((i=0;i<N;i++)); do
  if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/raw' --threads ${FASTQC_THREADS:-4} '${R1[$i]}' '${R2[$i]}'"
  else
    submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/raw' --threads ${FASTQC_THREADS:-4} '${R1[$i]}'"
  fi
done; wait_all
log "STEP 3 — MultiQC raw"
"${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/raw" -n multiQC_raw \
  -o "$OUTDIR/multiQC/raw" --data-format tsv --export -q

# ── Step 4: TrimGalore ────────────────────────────────────────────────────────
log "STEP 4 — TrimGalore ($LIBRARY_LAYOUT)"
for ((i=0;i<N;i++)); do
  submit "$REPO_ROOT/scripts/trimgalore_single.sh \
    '${R1[$i]}' '${R2[$i]}' '$OUTDIR/trimmedFastq' \
    '${TRIM_QUALITY:-20}' '${TRIM_MIN_LENGTH:-20}' '$LIBRARY_LAYOUT' '${SID[$i]}'"
done; wait_all

# ── Steps 5–6: FastQC / MultiQC trimmed ──────────────────────────────────────
log "STEP 5 — FastQC trimmed"
while IFS= read -r -d '' fq; do
  submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/trimmed' --threads ${FASTQC_THREADS:-4} '$fq'"
done < <(find "$OUTDIR/trimmedFastq" -name "*.fq.gz" -print0 2>/dev/null); wait_all
log "STEP 6 — MultiQC trimmed"
"${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/trimmed" -n multiQC_trimmed \
  -o "$OUTDIR/multiQC/trimmed" --data-format tsv --export -q

# ── Step 7: STAR alignment ────────────────────────────────────────────────────
log "STEP 7 — STAR alignment ($LIBRARY_LAYOUT)"
for ((i=0;i<N;i++)); do
  if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    _r1="$OUTDIR/trimmedFastq/${SID[$i]}_val_1.fq.gz"
    _r2="$OUTDIR/trimmedFastq/${SID[$i]}_val_2.fq.gz"
    submit "$REPO_ROOT/scripts/star_PE_single.sh \
      '$STAR_INDEX' '$OUTDIR/STARalignments' '${SID[$i]}' '$_r1' '$_r2' \
      '${STAR_THREADS:-15}' '${TMPDIR:-/tmp}'"
  else
    _r1="$OUTDIR/trimmedFastq/${SID[$i]}_trimmed.fq.gz"
    submit "$REPO_ROOT/scripts/star_SE_single.sh \
      '$STAR_INDEX' '$OUTDIR/STARalignments' '${SID[$i]}' '$_r1' \
      '${STAR_THREADS:-15}' '${TMPDIR:-/tmp}'"
  fi
done; wait_all
mv "$OUTDIR/STARalignments/"*ReadsPerGene.out.tab "$OUTDIR/STARgeneCounts/" 2>/dev/null || true
mv "$OUTDIR/STARalignments/"*Log.final.out         "$OUTDIR/STARlogs/"       2>/dev/null || true

# ── Step 8: samtools sort + index ─────────────────────────────────────────────
log "STEP 8 — samtools sort + index"
for bam in "$OUTDIR/STARalignments/"*Aligned.out.bam; do
  [[ -f "$bam" ]] || continue
  submit "$REPO_ROOT/scripts/bam_sort_index.sh '$bam' '$OUTDIR/bams' '${SAMTOOLS_THREADS:-4}'"
done; wait_all

# ── Step 9: MultiQC alignments ────────────────────────────────────────────────
log "STEP 9 — MultiQC alignments"
"${MULTIQC_BIN:-multiqc}" "$OUTDIR/STARlogs" -n multiQC_alignments \
  -o "$OUTDIR/multiQC/alignments" --data-format tsv --export -q

# ── Step 10: Raw bedGraph (R) ─────────────────────────────────────────────────
log "STEP 10 — Raw bedGraph (R: Rsamtools + GenomicAlignments + rtracklayer)"
"${RSCRIPT_BIN:-Rscript}" "$REPO_ROOT/scripts/Rscripts/bam_to_bedgraph.R" \
  --samplesheet "$SAMPLESHEET" --bamdir "$OUTDIR/bams" \
  --outdir "$OUTDIR/bedGraph/raw" --layout "$LIBRARY_LAYOUT"

# ── Step 11: DESeq2 normalization (R) ────────────────────────────────────────
log "STEP 11 — DESeq2 normalization (R)"
"${RSCRIPT_BIN:-Rscript}" "$REPO_ROOT/scripts/Rscripts/deseq2_normalize.R" \
  --samplesheet "$SAMPLESHEET" --countdir "$OUTDIR/STARgeneCounts" \
  --gtf "$GTF" --layout "$LIBRARY_LAYOUT" \
  --outdir "$OUTDIR/analysis/counts" --design "${DESIGN_FORMULA:-~ condition}"

# ── Step 12: Normalized bedGraph (R) ─────────────────────────────────────────
log "STEP 12 — Normalized bedGraph (R)"
"${RSCRIPT_BIN:-Rscript}" "$REPO_ROOT/scripts/Rscripts/normalize_bedgraph.R" \
  --samplesheet "$SAMPLESHEET" \
  --sffile "$OUTDIR/analysis/counts/size_factors.tsv" \
  --rawbgdir "$OUTDIR/bedGraph/raw" \
  --outdir "$OUTDIR/bedGraph/normalized" --layout "$LIBRARY_LAYOUT"

# ── Step 13: BigWig (SPECIES/NAMING env vars inherited) ──────────────────────
log "STEP 13 — BigWig  [species=$SPECIES  naming=$CHROMOSOME_NAMING  filter=$REGULAR_CHROMS_ONLY]"
for bg in "$OUTDIR/bedGraph/normalized/"*.bedGraph.gz; do
  [[ -f "$bg" ]] || continue
  submit "$REPO_ROOT/scripts/norm_bedgraph_to_bigwig.sh \
    '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '$KENTUTILS_DIR'"
done; wait_all

# ── Steps 14–15: Merge replicates + BigWig ────────────────────────────────────
if [[ "${MERGE_REPLICATES:-true}" == "true" ]]; then
  log "STEP 14 — Merge replicates (R)"
  "${RSCRIPT_BIN:-Rscript}" "$REPO_ROOT/scripts/Rscripts/merge_bedgraph_replicates.R" \
    --samplesheet "$SAMPLESHEET" --bgdir "$OUTDIR/bedGraph/normalized" \
    --outdir "$OUTDIR/bedGraph/merged" --genome "$GENOME_ASSEMBLY" \
    --layout "$LIBRARY_LAYOUT"
  log "STEP 15 — BigWig (merged)"
  for bg in "$OUTDIR/bedGraph/merged/"*.bedGraph; do
    [[ -f "$bg" ]] || continue
    submit "$REPO_ROOT/scripts/norm_bedgraph_to_bigwig.sh \
      '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '$KENTUTILS_DIR'"
  done; wait_all
fi

# ── Step 16: DESeq2 DE (R) ────────────────────────────────────────────────────
if [[ "${RUN_DE:-true}" == "true" ]]; then
  log "STEP 16 — DESeq2 DE (R)"
  "${RSCRIPT_BIN:-Rscript}" "$REPO_ROOT/scripts/Rscripts/deseq2_de.R" \
    --countsrdata "$OUTDIR/analysis/counts/dds.RData" \
    --contrasts "${CONTRASTS:-config/contrasts.csv}" \
    --outdir "$OUTDIR/analysis/DE" \
    --padj "${PADJ_THRESHOLD:-0.05}" --lfc "${LFC_THRESHOLD:-1}"
fi

# ── Step 17: QC plots (R) ─────────────────────────────────────────────────────
log "STEP 17 — DESeq2 QC plots (R)"
"${RSCRIPT_BIN:-Rscript}" "$REPO_ROOT/scripts/Rscripts/deseq2_qc_plots.R" \
  --countsrdata "$OUTDIR/analysis/counts/dds.RData" \
  --outdir "$OUTDIR/analysis/figures"

# ── Step 18: UCSC tracks ──────────────────────────────────────────────────────
if [[ "${UCSC_TRACKS:-true}" == "true" ]]; then
  log "STEP 18 — UCSC tracks"
  "$REPO_ROOT/scripts/create_ucsc_tracks.sh" \
    "$OUTDIR/bigwig" "$OUTDIR/reports/ucsc_tracks.txt" \
    "$OUTDIR/reports/bigwig_summary.txt" "${UCSC_BASE_URL:-PLACEHOLDER_URL}"
fi

# ── Step 19: MultiQC final ────────────────────────────────────────────────────
log "STEP 19 — MultiQC final"
"${MULTIQC_BIN:-multiqc}" \
  "$OUTDIR/fastQC/raw" "$OUTDIR/fastQC/trimmed" "$OUTDIR/STARlogs" \
  -n multiQC_final -o "$OUTDIR/multiQC/final" --data-format tsv --export -q

# ── Step 20: R Markdown report ────────────────────────────────────────────────
log "STEP 20 — Pipeline report (R Markdown)"
"${RSCRIPT_BIN:-Rscript}" -e "
  rmarkdown::render(
    input       = '$REPO_ROOT/scripts/Rscripts/pipeline_report.Rmd',
    params      = list(outdir='$OUTDIR', config='$CONFIG',
                       samplesheet='$SAMPLESHEET', species='${SPECIES}',
                       layout='${LIBRARY_LAYOUT}'),
    output_file = '$OUTDIR/reports/pipeline_report.html',
    quiet       = TRUE)
" || log "WARNING: pipeline report failed (check pandoc) — continuing"

log "rnaseq2tracks v3 complete.  Results: $OUTDIR"
