#!/usr/bin/env bash
# =============================================================================
# rnaseq2tracks.sh — master orchestrator (v4.3 — enrichment analysis)
# =============================================================================
# v4.2 changes vs v4.1:
#   Step 10 — bam_to_bedgraph.R:          1 R job per sample (parallel)
#   Step 12 — normalize_bedgraph.R:        1 R job per sample (parallel)
#   Step 14 — merge_bedgraph_replicates.R: 1 R job per condition (parallel)
#   Per-sample resume in Steps 10 and 12 (finer granularity than before)
#
# Usage: ./scripts/rnaseq2tracks.sh config/config.conf
# =============================================================================
set -euo pipefail
[[ $# -ne 1 ]] && { echo "Usage: $0 <config>" >&2; exit 1; }
CONFIG="$(realpath "$1")"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"
SAMPLESHEET="$(realpath "${SAMPLESHEET}")"
[[ -n "${CONTRASTS:-}" ]] && CONTRASTS="$(realpath "${CONTRASTS}")"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
skip() { log "SKIP — $* (output exists; set FORCE_RERUN=1 to rerun)"; }

FORCE_RERUN="${FORCE_RERUN:-0}"
done_check() {
  [[ "$FORCE_RERUN" == "1" ]] && return 1
  [[ -e "$1" ]] && return 0 || return 1
}

# ── Job throttle ──────────────────────────────────────────────────────────────
declare -a _PIDS=()
submit() {
  while [[ ${#_PIDS[@]} -ge ${MAX_JOBS:-8} ]]; do
    local live=()
    for p in "${_PIDS[@]}"; do kill -0 "$p" 2>/dev/null && live+=("$p"); done
    _PIDS=("${live[@]+"${live[@]}"}"); [[ ${#_PIDS[@]} -ge ${MAX_JOBS:-8} ]] && sleep 2
  done
  eval "$@" &
  _PIDS+=($!)
}
wait_all() {
  local ok=0
  for p in "${_PIDS[@]+"${_PIDS[@]}"}"; do wait "$p" || ok=1; done
  _PIDS=()
  [[ $ok -eq 0 ]] || { log "ERROR: a background job failed"; exit 1; }
}

# ── Step 0: Preflight ─────────────────────────────────────────────────────────
log "STEP 0 — Preflight checks"
"$REPO/scripts/preflight_check.sh" "$CONFIG"

# ── Species path resolution ───────────────────────────────────────────────────
SPECIES="${SPECIES:-mouse}"
case "$SPECIES" in
  human)
    STAR_INDEX="${STAR_INDEX_HUMAN:?STAR_INDEX_HUMAN not set}"
    GTF="${GTF_HUMAN:?GTF_HUMAN not set}"
    CHROM_SIZES="${CHROM_SIZES_HUMAN:?CHROM_SIZES_HUMAN not set}"
    RSEQC_BED="${RSEQC_BED_HUMAN:-}" ;;
  mouse)
    STAR_INDEX="${STAR_INDEX_MOUSE:?STAR_INDEX_MOUSE not set}"
    GTF="${GTF_MOUSE:?GTF_MOUSE not set}"
    CHROM_SIZES="${CHROM_SIZES_MOUSE:?CHROM_SIZES_MOUSE not set}"
    RSEQC_BED="${RSEQC_BED_MOUSE:-}" ;;
  *) echo "ERROR: SPECIES must be human|mouse" >&2; exit 1 ;;
esac
export SPECIES CHROMOSOME_NAMING="${CHROMOSOME_NAMING:-ucsc}" \
       REGULAR_CHROMS_ONLY="${REGULAR_CHROMS_ONLY:-true}"

[[ "$LIBRARY_LAYOUT" =~ ^(SE|PE)$ ]] || \
  { echo "ERROR: LIBRARY_LAYOUT must be SE|PE" >&2; exit 1; }

# ── Step 1: Output tree ───────────────────────────────────────────────────────
OUTDIR="${OUTDIR:-$(pwd)/rnaseq2tracks_output}"
log "STEP 1 — Output: $OUTDIR"
mkdir -p \
  "$OUTDIR/fastQC/raw"          "$OUTDIR/fastQC/trimmed" \
  "$OUTDIR/multiQC/raw"         "$OUTDIR/multiQC/trimmed" \
  "$OUTDIR/multiQC/alignments"  "$OUTDIR/multiQC/final" \
  "$OUTDIR/trimmedFastq"        "$OUTDIR/STARalignments" \
  "$OUTDIR/STARlogs"            "$OUTDIR/STARgeneCounts" \
  "$OUTDIR/bams" \
  "$OUTDIR/07_qc/star"          "$OUTDIR/07_qc/rseqc" \
  "$OUTDIR/07_qc/multiqc" \
  "$OUTDIR/bedGraph/raw"        "$OUTDIR/bedGraph/normalized" \
  "$OUTDIR/bedGraph/merged" \
  "$OUTDIR/bigwig" \
  "$OUTDIR/analysis/counts"     "$OUTDIR/analysis/DE" \
  "$OUTDIR/analysis/figures"    "$OUTDIR/reports"

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
N=${#SID[@]}; log "Loaded $N samples  layout=$LIBRARY_LAYOUT  species=$SPECIES"

# ── Step 2–3: FastQC / MultiQC raw ───────────────────────────────────────────
_s2="$OUTDIR/fastQC/raw/$(basename "${R1[0]}" .fq.gz)_fastqc.html"
if done_check "$_s2"; then skip "STEP 2 — FastQC raw"
else
  log "STEP 2 — FastQC raw"
  for ((i=0;i<N;i++)); do
    if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
      submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/raw' \
        --threads ${FASTQC_THREADS:-4} '${R1[$i]}' '${R2[$i]}'"
    else
      submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/raw' \
        --threads ${FASTQC_THREADS:-4} '${R1[$i]}'"
    fi
  done; wait_all
fi
if done_check "$OUTDIR/multiQC/raw/multiQC_raw.html"; then skip "STEP 3 — MultiQC raw"
else
  log "STEP 3 — MultiQC raw"
  "${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/raw" -n multiQC_raw \
    -o "$OUTDIR/multiQC/raw" --data-format tsv --export -q
fi

# ── Step 4: TrimGalore ────────────────────────────────────────────────────────
_s4=$(if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
  echo "$OUTDIR/trimmedFastq/${SID[0]}_val_1.fq.gz"
else echo "$OUTDIR/trimmedFastq/${SID[0]}_trimmed.fq.gz"; fi)
if done_check "$_s4"; then skip "STEP 4 — TrimGalore"
else
  log "STEP 4 — TrimGalore ($LIBRARY_LAYOUT)"
  for ((i=0;i<N;i++)); do
    submit "$REPO/scripts/trimgalore_single.sh \
      '${R1[$i]}' '${R2[$i]}' '$OUTDIR/trimmedFastq' \
      '${TRIM_QUALITY:-20}' '${TRIM_MIN_LENGTH:-20}' '$LIBRARY_LAYOUT' '${SID[$i]}'"
  done; wait_all
fi

# ── Step 5–6: FastQC / MultiQC trimmed ───────────────────────────────────────
if done_check "$OUTDIR/multiQC/trimmed/multiQC_trimmed.html"; then
  skip "STEP 5+6 — FastQC/MultiQC trimmed"
else
  log "STEP 5 — FastQC trimmed"
  while IFS= read -r -d '' fq; do
    submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/trimmed' \
      --threads ${FASTQC_THREADS:-4} '$fq'"
  done < <(find "$OUTDIR/trimmedFastq" -name "*.fq.gz" -print0 2>/dev/null); wait_all
  log "STEP 6 — MultiQC trimmed"
  "${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/trimmed" -n multiQC_trimmed \
    -o "$OUTDIR/multiQC/trimmed" --data-format tsv --export -q
fi

# ── Step 7: STAR ──────────────────────────────────────────────────────────────
_s7="$OUTDIR/STARlogs/${SID[0]}_Log.final.out"
if done_check "$_s7"; then skip "STEP 7 — STAR alignment"
else
  log "STEP 7 — STAR alignment"
  for ((i=0;i<N;i++)); do
    if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
      _r1="$OUTDIR/trimmedFastq/${SID[$i]}_val_1.fq.gz"
      _r2="$OUTDIR/trimmedFastq/${SID[$i]}_val_2.fq.gz"
      submit "$REPO/scripts/star_PE_single.sh \
        '$STAR_INDEX' '$OUTDIR/STARalignments' '${SID[$i]}' '$_r1' '$_r2' \
        '${STAR_THREADS:-15}' '${TMPDIR:-/tmp}'"
    else
      _r1="$OUTDIR/trimmedFastq/${SID[$i]}_trimmed.fq.gz"
      submit "$REPO/scripts/star_SE_single.sh \
        '$STAR_INDEX' '$OUTDIR/STARalignments' '${SID[$i]}' '$_r1' \
        '${STAR_THREADS:-15}' '${TMPDIR:-/tmp}'"
    fi
  done; wait_all
  mv "$OUTDIR/STARalignments/"*ReadsPerGene.out.tab "$OUTDIR/STARgeneCounts/" 2>/dev/null || true
  mv "$OUTDIR/STARalignments/"*Log.final.out         "$OUTDIR/STARlogs/"       2>/dev/null || true
fi

# ── Step 8: samtools sort + index ─────────────────────────────────────────────
_s8="$OUTDIR/bams/${SID[0]}_sortedS.bam"
if done_check "$_s8"; then skip "STEP 8 — samtools sort+index"
else
  log "STEP 8 — samtools sort + index"
  for bam in "$OUTDIR/STARalignments/"*Aligned.out.bam; do
    [[ -f "$bam" ]] || continue
    submit "$REPO/scripts/bam_sort_index.sh '$bam' '$OUTDIR/bams' '${SAMTOOLS_THREADS:-4}'"
  done; wait_all
fi

# ── Step 9: MultiQC alignments ────────────────────────────────────────────────
if done_check "$OUTDIR/multiQC/alignments/multiQC_alignments.html"; then
  skip "STEP 9 — MultiQC alignments"
else
  log "STEP 9 — MultiQC alignments"
  "${MULTIQC_BIN:-multiqc}" "$OUTDIR/STARlogs" -n multiQC_alignments \
    -o "$OUTDIR/multiQC/alignments" --data-format tsv --export -q
fi

# ── Step 9b: STAR QC summary ──────────────────────────────────────────────────
if done_check "$OUTDIR/07_qc/star/star_alignment_summary.tsv"; then
  skip "STEP 9b — STAR alignment summary"
else
  log "STEP 9b — STAR alignment summary"
  "$REPO/scripts/collect_star_qc.sh" "$OUTDIR/STARlogs" "$OUTDIR/07_qc"
fi

# ── Step 10: bam_to_bedgraph.R — PARALLEL (1 job per sample) ─────────────────
log "STEP 10 — bam_to_bedgraph.R (parallel: 1 job per sample)"
_any_s10_missing=0
for ((i=0;i<N;i++)); do
  _fwd="$OUTDIR/bedGraph/raw/${SID[$i]}_FwdS.bedGraph.gz"
  _uns="$OUTDIR/bedGraph/raw/${SID[$i]}_unstranded.bedGraph.gz"
  if [[ -f "$_fwd" || -f "$_uns" ]] && [[ "$FORCE_RERUN" != "1" ]]; then
    log "  SKIP ${SID[$i]} (bedGraph exists)"
  else
    _any_s10_missing=1
    submit "${RSCRIPT_BIN:-Rscript} '$REPO/scripts/Rscripts/bam_to_bedgraph.R' \
      --sample_id '${SID[$i]}' \
      --bam '$OUTDIR/bams/${SID[$i]}_sortedS.bam' \
      --strandedness '${STRAND[$i]}' \
      --outdir '$OUTDIR/bedGraph/raw' \
      --layout '$LIBRARY_LAYOUT'"
  fi
done
wait_all
[[ $_any_s10_missing -eq 0 ]] && skip "STEP 10 — all bedGraphs already present" || true

# ── Step 10b: Strand consistency (always runs — fast safety check) ────────────
log "STEP 10b — Strand consistency check"
"$REPO/scripts/check_strand_consistency.sh" \
  "$SAMPLESHEET" "$OUTDIR/bams" "$LIBRARY_LAYOUT" "${STRAND_TOLERANCE_PCT:-5}" "${MAX_JOBS:-8}"

# ── Step 10c: RSeQC — background (steps 11-18 run in parallel) ─────────────
_s10c_sentinel="$OUTDIR/07_qc/rseqc/infer_experiment/${SID[0]}_infer_experiment.txt"
RSEQC_BG_PID=""
if [[ "${RUN_RSEQC:-true}" == "true" && -n "${RSEQC_BED:-}" && -f "${RSEQC_BED:-/dev/null}" ]]; then
  if done_check "$_s10c_sentinel" && done_check "$OUTDIR/07_qc/multiqc/multiQC_rseqc.html"; then
    skip "STEP 10c — RSeQC (already complete)"
  else
    log "STEP 10c — RSeQC launching in background (PID will follow)"
    (
      if ! done_check "$_s10c_sentinel"; then
        "$REPO/scripts/run_rnaseq_qc.sh"           "$SAMPLESHEET" "$OUTDIR/bams" "$OUTDIR/07_qc" "$RSEQC_BED"           "${RSEQC_BIN_DIR:-}" "$LIBRARY_LAYOUT" "${MAX_JOBS:-8}"
      fi
      if ! done_check "$OUTDIR/07_qc/multiqc/multiQC_rseqc.html"; then
        MQC_RSEQC=()
        for d in read_distribution junction_annotation junction_saturation genebody; do
          [[ -d "$OUTDIR/07_qc/rseqc/$d" ]] && MQC_RSEQC+=("$OUTDIR/07_qc/rseqc/$d")
        done
        [[ ${#MQC_RSEQC[@]} -gt 0 ]] &&         "${MULTIQC_BIN:-multiqc}" "${MQC_RSEQC[@]}"           -n multiQC_rseqc -o "$OUTDIR/07_qc/multiqc"           --data-format tsv --export -q || true
      fi
    ) &
    RSEQC_BG_PID=$!
    log "STEP 10c — RSeQC running in background PID=$RSEQC_BG_PID"
  fi
else
  log "STEP 10c — RSeQC SKIPPED (RUN_RSEQC=false or RSEQC_BED not found)"
fi

# ── Step 11: DESeq2 normalization (must be serial — needs all samples) ────────
if done_check "$OUTDIR/analysis/counts/dds.RData"; then
  skip "STEP 11 — DESeq2 normalization"
else
  log "STEP 11 — DESeq2 normalization"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_normalize.R" \
    --samplesheet "$SAMPLESHEET" --countdir "$OUTDIR/STARgeneCounts" \
    --gtf "$GTF" --layout "$LIBRARY_LAYOUT" \
    --outdir "$OUTDIR/analysis/counts" --design "${DESIGN_FORMULA:-~ condition}"
fi

# ── Step 12: normalize_bedgraph.R — PARALLEL (1 job per sample) ──────────────
log "STEP 12 — normalize_bedgraph.R (parallel: 1 job per sample)"
for ((i=0;i<N;i++)); do
  _nfwd="$OUTDIR/bedGraph/normalized/${SID[$i]}_FwdS_norm.bedGraph.gz"
  _nuns="$OUTDIR/bedGraph/normalized/${SID[$i]}_unstranded_norm.bedGraph.gz"
  if [[ -f "$_nfwd" || -f "$_nuns" ]] && [[ "$FORCE_RERUN" != "1" ]]; then
    log "  SKIP ${SID[$i]} (normalized bedGraph exists)"
  else
    submit "${RSCRIPT_BIN:-Rscript} '$REPO/scripts/Rscripts/normalize_bedgraph.R' \
      --sample_id '${SID[$i]}' \
      --strandedness '${STRAND[$i]}' \
      --sffile '$OUTDIR/analysis/counts/size_factors.tsv' \
      --rawbgdir '$OUTDIR/bedGraph/raw' \
      --outdir '$OUTDIR/bedGraph/normalized' \
      --layout '$LIBRARY_LAYOUT'"
  fi
done; wait_all

# ── Step 13: BigWig per sample ────────────────────────────────────────────────
_s13="$OUTDIR/bigwig/${SID[0]}_FwdS_norm.bw"
[[ ! -f "$_s13" ]] && _s13="$OUTDIR/bigwig/${SID[0]}_unstranded_norm.bw"
if done_check "$_s13"; then skip "STEP 13 — BigWig per sample"
else
  log "STEP 13 — BigWig [species=$SPECIES naming=$CHROMOSOME_NAMING filter=$REGULAR_CHROMS_ONLY]"
  for bg in "$OUTDIR/bedGraph/normalized/"*_norm.bedGraph.gz; do
    [[ -f "$bg" ]] || continue
    submit "$REPO/scripts/norm_bedgraph_to_bigwig.sh \
      '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '${KENTUTILS_DIR}'"
  done; wait_all
fi

# ── Step 14: merge_bedgraph_replicates.R — PARALLEL (1 job per condition) ─────
_n14=$(find "$OUTDIR/bedGraph/merged" -name "*_merged.bedGraph" 2>/dev/null | wc -l)
if [[ "$_n14" -gt 0 ]] && [[ "$FORCE_RERUN" != "1" ]]; then
  skip "STEP 14 — merge_bedgraph_replicates.R"
else
  log "STEP 14 — merge_bedgraph_replicates.R (parallel: 1 job per condition)"
  # Build unique conditions with their sample IDs and strandedness
  declare -A COND_SIDS COND_STRAND
  for ((i=0;i<N;i++)); do
    c="${COND[$i]}"
    COND_SIDS["$c"]="${COND_SIDS[$c]:-}${COND_SIDS[$c]:+,}${SID[$i]}"
    COND_STRAND["$c"]="${STRAND[$i]}"
  done
  for c in "${!COND_SIDS[@]}"; do
    submit "${RSCRIPT_BIN:-Rscript} '$REPO/scripts/Rscripts/merge_bedgraph_replicates.R' \
      --condition '$c' \
      --sample_ids '${COND_SIDS[$c]}' \
      --strandedness '${COND_STRAND[$c]}' \
      --bgdir '$OUTDIR/bedGraph/normalized' \
      --outdir '$OUTDIR/bedGraph/merged' \
      --layout '$LIBRARY_LAYOUT'"
  done; wait_all
fi

# ── Step 15: merged BigWigs ───────────────────────────────────────────────────
_n15=$(find "$OUTDIR/bigwig" -name "*_merged.bw" 2>/dev/null | wc -l)
if [[ "$_n15" -gt 0 ]] && [[ "$FORCE_RERUN" != "1" ]]; then
  skip "STEP 15 — merged BigWigs"
else
  log "STEP 15 — merged BigWigs"
  for bg in "$OUTDIR/bedGraph/merged/"*_merged.bedGraph; do
    [[ -f "$bg" ]] || continue
    submit "$REPO/scripts/norm_bedgraph_to_bigwig.sh \
      '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '${KENTUTILS_DIR}'"
  done; wait_all
fi

# ── Step 16: DESeq2 DE ────────────────────────────────────────────────────────
_n16=$(find "$OUTDIR/analysis/DE" -name "*_DE_results.tsv" 2>/dev/null | wc -l)
if [[ "$_n16" -gt 0 ]] && [[ "$FORCE_RERUN" != "1" ]]; then
  skip "STEP 16 — DESeq2 DE"
else
  if [[ -f "${CONTRASTS_FILE:-$REPO/config/contrasts.csv}" ]]; then
    export GTF
  log "STEP 16 — DESeq2 DE"
    "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_de.R" \
      --countsrdata "$OUTDIR/analysis/counts/dds.RData" \
      --contrasts "${CONTRASTS_FILE:-$REPO/config/contrasts.csv}" \
      --outdir "$OUTDIR/analysis/DE" \
      --padj "${DE_PADJ_THRESHOLD:-0.05}" --lfc "${DE_LFC_THRESHOLD:-1}"
  else
    log "STEP 16 — DESeq2 DE SKIPPED (no contrasts.csv)"
  fi
fi

# ── Step 17: DESeq2 QC plots ──────────────────────────────────────────────────
if done_check "$OUTDIR/analysis/figures/PCA.pdf"; then
  skip "STEP 17 — DESeq2 QC plots"
else
  log "STEP 17 — DESeq2 QC plots"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_qc_plots.R" \
    --countsrdata "$OUTDIR/analysis/counts/dds.RData" \
    --outdir "$OUTDIR/analysis/figures"
fi

# ── Step 18: UCSC tracks ──────────────────────────────────────────────────────
if done_check "$OUTDIR/reports/ucsc_tracks.txt"; then
  skip "STEP 18 — UCSC tracks"
else
  if [[ -n "${UCSC_BASE_URL:-}" ]]; then
    log "STEP 18 — UCSC tracks"
    "$REPO/scripts/create_ucsc_tracks.sh" \
      "$OUTDIR/bigwig" "$OUTDIR/reports/ucsc_tracks.txt" \
      "$OUTDIR/reports/bigwig_summary.txt" "$UCSC_BASE_URL"
  else
    log "STEP 18 — UCSC tracks SKIPPED (UCSC_BASE_URL not set)"
  fi
fi

# ── Wait for background RSeQC before final MultiQC ──────────────────────────
if [[ -n "${RSEQC_BG_PID:-}" ]]; then
  log "STEP 19 — waiting for background RSeQC (PID $RSEQC_BG_PID)..."
  wait "$RSEQC_BG_PID" || log "WARNING: RSeQC background job had errors — continuing"
fi
# ── Step 19: MultiQC final ────────────────────────────────────────────────────
if done_check "$OUTDIR/multiQC/final/multiQC_final.html"; then
  skip "STEP 19 — MultiQC final"
else
  log "STEP 19 — MultiQC final"
  MQC_SOURCES=("$OUTDIR/fastQC" "$OUTDIR/multiQC/raw" "$OUTDIR/multiQC/trimmed"
               "$OUTDIR/STARlogs" "$OUTDIR/multiQC/alignments")
  for d in read_distribution junction_annotation junction_saturation genebody; do
    [[ -d "$OUTDIR/07_qc/rseqc/$d" ]] && MQC_SOURCES+=("$OUTDIR/07_qc/rseqc/$d")
  done
  "${MULTIQC_BIN:-multiqc}" "${MQC_SOURCES[@]}" \
  --ignore "$OUTDIR/07_qc/multiqc" \
    -n multiQC_final -o "$OUTDIR/multiQC/final" \
    --data-format tsv --export -q
fi

# ── Step 20: Pipeline report ──────────────────────────────────────────────────
if done_check "$OUTDIR/reports/pipeline_report.html"; then
  skip "STEP 20 — Pipeline report"
else
  log "STEP 20 — Pipeline report"
  "${RSCRIPT_BIN:-Rscript}" -e "
    rmarkdown::render(
      input             = '$REPO/scripts/Rscripts/pipeline_report.Rmd',
      output_file       = 'pipeline_report.html',
      output_dir        = '$OUTDIR/reports',
      intermediates_dir = '$OUTDIR/reports',
      knit_root_dir     = '$OUTDIR',
      params = list(
        outdir      = '$OUTDIR',
        config      = '$CONFIG',
        samplesheet = '$SAMPLESHEET',
        species     = '$SPECIES',
        layout      = '$LIBRARY_LAYOUT'
      ),
      quiet = TRUE
    )
  " || log "WARNING: pipeline report failed — check pandoc"
fi


# ── Step 21: Gene enrichment analysis (ORA + GSEA) ───────────────────────────
if done_check "$OUTDIR/analysis/enrichment/.enrichment_done"; then
  skip "STEP 21 — Enrichment analysis"
else
  log "STEP 21 — Enrichment analysis (ORA + GSEA)"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_enrichment.R" \
    --dedir     "$OUTDIR/analysis/DE" \
    --contrasts "$(realpath "${CONTRASTS}")" \
    --outdir    "$OUTDIR/analysis/enrichment" \
    --species   "$SPECIES" \
    --padj      "${PADJ_THRESHOLD:-0.05}" \
    --lfc       "${LFC_THRESHOLD:-1}" \
    --minGS     "${ENRICHMENT_MINGS:-10}" \
    --maxGS     "${ENRICHMENT_MAXGS:-500}" \
    && touch "$OUTDIR/analysis/enrichment/.enrichment_done" \
    || log "WARNING: enrichment analysis failed — check deseq2_enrichment.R"
fi

log "rnaseq2tracks v4.3 complete.  Results: $OUTDIR"
