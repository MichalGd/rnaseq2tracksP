#!/usr/bin/env bash
# =============================================================================
# rnaseq2tracks.sh — master orchestrator (v4.1 — resume/checkpoint)
# =============================================================================
# Language division:
#   Bash : subprocess fan-out, job throttling, STAR, samtools, TrimGalore,
#          FastQC, MultiQC, RSeQC, BigWig conversion
#   R    : ALL analytical steps (coverage, normalization, DE, QC plots, report)
#
# Steps:
#   0    preflight_check.sh          tools + R packages + RSeQC + genome files
#   1    mkdir                        output tree incl. 07_qc/
#   2-3  FastQC / MultiQC raw
#   4    TrimGalore SE/PE
#   5-6  FastQC / MultiQC trimmed
#   7    STAR alignment (--quantMode GeneCounts)
#   8    samtools sort + index
#   9    MultiQC alignment logs
#   9b   collect_star_qc.sh          STAR summary TSV + MultiQC symlinks
#   10   bam_to_bedgraph.R           strand-aware coverage (R)
#   10b  check_strand_consistency.sh Fwd+Rev vs Total sanity check
#   10c  run_rnaseq_qc.sh            RSeQC 4 modules + geneBody
#        MultiQC RSeQC
#   11   deseq2_normalize.R          counts, SF, SF_rpm (R)
#   12   normalize_bedgraph.R        SF_rpm scaling (R)
#   13   norm_bedgraph_to_bigwig.sh  chr filter → BigWig
#   14   merge_bedgraph_replicates.R replicate averaging (R)
#   15   norm_bedgraph_to_bigwig.sh  merged BigWigs
#   16   deseq2_de.R                 Wald + apeglm LFC (R)
#   17   deseq2_qc_plots.R           PCA, clustering, heatmaps (R)
#   18   create_ucsc_tracks.sh
#   19   MultiQC final (all sources incl. RSeQC)
#   20   pipeline_report.Rmd → HTML (R)
#
# Resume behaviour:
#   Each step checks for a sentinel file/directory before running.
#   If the expected output already exists, the step is skipped with [SKIP].
#   To force re-run a step, delete its output or set FORCE_RERUN=1.
#
# Usage: ./scripts/rnaseq2tracks.sh config/config.conf
# =============================================================================
set -euo pipefail
[[ $# -ne 1 ]] && { echo "Usage: $0 <config>" >&2; exit 1; }
CONFIG="$(realpath "$1")"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
skip() { log "SKIP — $* (output exists; set FORCE_RERUN=1 to rerun)"; }

FORCE_RERUN="${FORCE_RERUN:-0}"
# done <sentinel> — returns 0 if step should be skipped, 1 if it must run
done_check() {
  [[ "$FORCE_RERUN" == "1" ]] && return 1
  local sentinel="$1"
  if [[ -e "$sentinel" ]]; then return 0; else return 1; fi
}

# ── Job throttle ─────────────────────────────────────────────────────────────
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
# Sentinel: fastQC/raw/<first_sample>_fastqc.html
_s2_sentinel="$OUTDIR/fastQC/raw/$(basename "${R1[0]}" .fq.gz)_fastqc.html"
if done_check "$_s2_sentinel"; then
  skip "STEP 2 — FastQC raw"
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

if done_check "$OUTDIR/multiQC/raw/multiQC_raw.html"; then
  skip "STEP 3 — MultiQC raw"
else
  log "STEP 3 — MultiQC raw"
  "${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/raw" -n multiQC_raw \
    -o "$OUTDIR/multiQC/raw" --data-format tsv --export -q
fi

# ── Step 4: TrimGalore ────────────────────────────────────────────────────────
# Sentinel: trimmedFastq/<first_sample>_val_1.fq.gz (PE) or _trimmed.fq.gz (SE)
if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
  _s4_sentinel="$OUTDIR/trimmedFastq/${SID[0]}_val_1.fq.gz"
else
  _s4_sentinel="$OUTDIR/trimmedFastq/${SID[0]}_trimmed.fq.gz"
fi
if done_check "$_s4_sentinel"; then
  skip "STEP 4 — TrimGalore"
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
# Sentinel: STARlogs/<first_sample>_Log.final.out  (logs are moved after alignment)
_s7_sentinel="$OUTDIR/STARlogs/${SID[0]}_Log.final.out"
if done_check "$_s7_sentinel"; then
  skip "STEP 7 — STAR alignment"
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
_s8_sentinel="$OUTDIR/bams/${SID[0]}_sortedS.bam"
if done_check "$_s8_sentinel"; then
  skip "STEP 8 — samtools sort+index"
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

# ── Step 9b: STAR alignment summary ──────────────────────────────────────────
if done_check "$OUTDIR/07_qc/star/star_alignment_summary.tsv"; then
  skip "STEP 9b — STAR alignment summary"
else
  log "STEP 9b — STAR alignment summary"
  "$REPO/scripts/collect_star_qc.sh" "$OUTDIR/STARlogs" "$OUTDIR/07_qc"
fi

# ── Step 10: bedGraph (R) ─────────────────────────────────────────────────────
# Sentinel: first sample's bedGraph (unstranded or FwdS)
_s10_sentinel_fwd="$OUTDIR/bedGraph/raw/${SID[0]}_FwdS.bedGraph.gz"
_s10_sentinel_uns="$OUTDIR/bedGraph/raw/${SID[0]}_unstranded.bedGraph.gz"
if done_check "$_s10_sentinel_fwd" || done_check "$_s10_sentinel_uns"; then
  skip "STEP 10 — bam_to_bedgraph.R"
else
  log "STEP 10 — bam_to_bedgraph.R"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/bam_to_bedgraph.R" \
    --samplesheet "$SAMPLESHEET" --bamdir "$OUTDIR/bams" \
    --outdir "$OUTDIR/bedGraph/raw" --layout "$LIBRARY_LAYOUT"
fi

# ── Step 10b: Strand consistency ──────────────────────────────────────────────
# Always runs (fast, no output file; guards data integrity)
log "STEP 10b — Strand consistency check"
"$REPO/scripts/check_strand_consistency.sh" \
  "$SAMPLESHEET" "$OUTDIR/bams" "$LIBRARY_LAYOUT" "${STRAND_TOLERANCE_PCT:-5}"

# ── Step 10c: RSeQC ───────────────────────────────────────────────────────────
_s10c_sentinel="$OUTDIR/07_qc/rseqc/infer_experiment/${SID[0]}_infer_experiment.txt"
if [[ "${RUN_RSEQC:-true}" == "true" && -n "${RSEQC_BED:-}" && -f "${RSEQC_BED:-/dev/null}" ]]; then
  if done_check "$_s10c_sentinel"; then
    skip "STEP 10c — RSeQC"
  else
    log "STEP 10c — RSeQC RNA-seq QC"
    "$REPO/scripts/run_rnaseq_qc.sh" \
      "$SAMPLESHEET" "$OUTDIR/bams" "$OUTDIR/07_qc" "$RSEQC_BED" \
      "${RSEQC_BIN_DIR:-}" "$LIBRARY_LAYOUT" "${MAX_JOBS:-8}"
  fi
  if done_check "$OUTDIR/07_qc/multiqc/multiQC_rseqc.html"; then
    skip "STEP 10c — MultiQC RSeQC"
  else
    log "STEP 10c — MultiQC RSeQC"
    MQC_RSEQC=()
    for d in read_distribution junction_annotation junction_saturation genebody; do
      [[ -d "$OUTDIR/07_qc/rseqc/$d" ]] && MQC_RSEQC+=("$OUTDIR/07_qc/rseqc/$d")
    done
    [[ ${#MQC_RSEQC[@]} -gt 0 ]] && \
    "${MULTIQC_BIN:-multiqc}" "${MQC_RSEQC[@]}" \
      -n multiQC_rseqc -o "$OUTDIR/07_qc/multiqc" \
      --data-format tsv --export -q || true
  fi
else
  log "STEP 10c — RSeQC SKIPPED (RUN_RSEQC=false or RSEQC_BED not found)"
fi

# ── Step 11: DESeq2 normalization ─────────────────────────────────────────────
if done_check "$OUTDIR/analysis/counts/dds.RData"; then
  skip "STEP 11 — DESeq2 normalization"
else
  log "STEP 11 — DESeq2 normalization"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_normalize.R" \
    --samplesheet "$SAMPLESHEET" --countdir "$OUTDIR/STARgeneCounts" \
    --gtf "$GTF" --layout "$LIBRARY_LAYOUT" \
    --outdir "$OUTDIR/analysis/counts" --design "${DESIGN_FORMULA:-~ condition}"
fi

# ── Step 12: Normalize bedGraph ───────────────────────────────────────────────
_s12_sentinel_fwd="$OUTDIR/bedGraph/normalized/${SID[0]}_FwdS_norm.bedGraph.gz"
_s12_sentinel_uns="$OUTDIR/bedGraph/normalized/${SID[0]}_unstranded_norm.bedGraph.gz"
if done_check "$_s12_sentinel_fwd" || done_check "$_s12_sentinel_uns"; then
  skip "STEP 12 — normalize_bedgraph.R"
else
  log "STEP 12 — normalize_bedgraph.R"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/normalize_bedgraph.R" \
    --samplesheet "$SAMPLESHEET" \
    --sffile "$OUTDIR/analysis/counts/size_factors.tsv" \
    --rawbgdir "$OUTDIR/bedGraph/raw" \
    --outdir "$OUTDIR/bedGraph/normalized" \
    --layout "$LIBRARY_LAYOUT"
fi

# ── Step 13: BigWig ───────────────────────────────────────────────────────────
_s13_sentinel="$OUTDIR/bigwig/${SID[0]}_FwdS_norm.bw"
[[ ! -f "$_s13_sentinel" ]] && _s13_sentinel="$OUTDIR/bigwig/${SID[0]}_unstranded_norm.bw"
if done_check "$_s13_sentinel"; then
  skip "STEP 13 — norm_bedgraph_to_bigwig (per-sample)"
else
  log "STEP 13 — BigWig [species=$SPECIES naming=$CHROMOSOME_NAMING filter=$REGULAR_CHROMS_ONLY]"
  for bg in "$OUTDIR/bedGraph/normalized/"*_norm.bedGraph.gz; do
    [[ -f "$bg" ]] || continue
    submit "$REPO/scripts/norm_bedgraph_to_bigwig.sh \
      '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '${KENTUTILS_DIR}'"
  done; wait_all
fi

# ── Step 14–15: Merge replicates + merged BigWigs ────────────────────────────
_s14_sentinel="$OUTDIR/bedGraph/merged"
_n_merged=$(find "$OUTDIR/bedGraph/merged" -name "*_merged.bedGraph" 2>/dev/null | wc -l)
if [[ "$_n_merged" -gt 0 ]] && [[ "$FORCE_RERUN" != "1" ]]; then
  skip "STEP 14 — merge_bedgraph_replicates.R"
else
  log "STEP 14 — merge_bedgraph_replicates.R"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/merge_bedgraph_replicates.R" \
    --samplesheet "$SAMPLESHEET" \
    --bgdir "$OUTDIR/bedGraph/normalized" \
    --outdir "$OUTDIR/bedGraph/merged" \
    --layout "$LIBRARY_LAYOUT"
fi

_n_merged_bw=$(find "$OUTDIR/bigwig" -name "*_merged.bw" 2>/dev/null | wc -l)
if [[ "$_n_merged_bw" -gt 0 ]] && [[ "$FORCE_RERUN" != "1" ]]; then
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
_n_de=$(find "$OUTDIR/analysis/DE" -name "*_DE_results.tsv" 2>/dev/null | wc -l)
if [[ "$_n_de" -gt 0 ]] && [[ "$FORCE_RERUN" != "1" ]]; then
  skip "STEP 16 — DESeq2 DE"
else
  if [[ -f "${CONTRASTS_FILE:-$REPO/config/contrasts.csv}" ]]; then
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
      '$REPO/scripts/Rscripts/pipeline_report.Rmd',
      output_file = '$OUTDIR/reports/pipeline_report.html',
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

log "rnaseq2tracks v4.1 complete.  Results: $OUTDIR"
