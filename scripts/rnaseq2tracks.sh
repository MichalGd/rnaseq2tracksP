#!/usr/bin/env bash
# =============================================================================
# rnaseq2tracks.sh — master orchestrator (v4.0)
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
#   11   deseq2_normalize.R          counts, SF, SF_rpm, FPKM, TPM (R)
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
# Usage: ./scripts/rnaseq2tracks.sh config/config.conf
# =============================================================================
set -euo pipefail
[[ $# -ne 1 ]] && { echo "Usage: $0 <config>" >&2; exit 1; }
CONFIG="$(realpath "$1")"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
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

[[ "$LIBRARY_LAYOUT" =~ ^(SE|PE)$ ]] || { echo "ERROR: LIBRARY_LAYOUT must be SE|PE" >&2; exit 1; }

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
log "STEP 3 — MultiQC raw"
"${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/raw" -n multiQC_raw \
  -o "$OUTDIR/multiQC/raw" --data-format tsv --export -q

# ── Step 4: TrimGalore ────────────────────────────────────────────────────────
log "STEP 4 — TrimGalore ($LIBRARY_LAYOUT)"
for ((i=0;i<N;i++)); do
  submit "$REPO/scripts/trimgalore_single.sh \
    '${R1[$i]}' '${R2[$i]}' '$OUTDIR/trimmedFastq' \
    '${TRIM_QUALITY:-20}' '${TRIM_MIN_LENGTH:-20}' '$LIBRARY_LAYOUT' '${SID[$i]}'"
done; wait_all

# ── Step 5–6: FastQC / MultiQC trimmed ───────────────────────────────────────
log "STEP 5 — FastQC trimmed"
while IFS= read -r -d '' fq; do
  submit "${FASTQC_BIN:-fastqc} --outdir '$OUTDIR/fastQC/trimmed' \
    --threads ${FASTQC_THREADS:-4} '$fq'"
done < <(find "$OUTDIR/trimmedFastq" -name "*.fq.gz" -print0 2>/dev/null); wait_all
log "STEP 6 — MultiQC trimmed"
"${MULTIQC_BIN:-multiqc}" "$OUTDIR/fastQC/trimmed" -n multiQC_trimmed \
  -o "$OUTDIR/multiQC/trimmed" --data-format tsv --export -q

# ── Step 7: STAR ──────────────────────────────────────────────────────────────
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

# ── Step 8: samtools sort + index ─────────────────────────────────────────────
log "STEP 8 — samtools sort + index"
for bam in "$OUTDIR/STARalignments/"*Aligned.out.bam; do
  [[ -f "$bam" ]] || continue
  submit "$REPO/scripts/bam_sort_index.sh '$bam' '$OUTDIR/bams' '${SAMTOOLS_THREADS:-4}'"
done; wait_all

# ── Step 9: MultiQC alignments ────────────────────────────────────────────────
log "STEP 9 — MultiQC alignments"
"${MULTIQC_BIN:-multiqc}" "$OUTDIR/STARlogs" -n multiQC_alignments \
  -o "$OUTDIR/multiQC/alignments" --data-format tsv --export -q

# ── Step 9b: STAR alignment summary ──────────────────────────────────────────
log "STEP 9b — STAR alignment summary"
"$REPO/scripts/collect_star_qc.sh" "$OUTDIR/STARlogs" "$OUTDIR/07_qc"

# ── Step 10: bedGraph (R) ─────────────────────────────────────────────────────
log "STEP 10 — bam_to_bedgraph.R"
"${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/bam_to_bedgraph.R" \
  --samplesheet "$SAMPLESHEET" --bamdir "$OUTDIR/bams" \
  --outdir "$OUTDIR/bedGraph/raw" --layout "$LIBRARY_LAYOUT"

# ── Step 10b: Strand consistency ──────────────────────────────────────────────
log "STEP 10b — Strand consistency check"
"$REPO/scripts/check_strand_consistency.sh" \
  "$SAMPLESHEET" "$OUTDIR/bams" "$LIBRARY_LAYOUT" "${STRAND_TOLERANCE_PCT:-5}"

# ── Step 10c: RSeQC ───────────────────────────────────────────────────────────
if [[ "${RUN_RSEQC:-true}" == "true" && -n "${RSEQC_BED:-}" && -f "${RSEQC_BED:-/dev/null}" ]]; then
  log "STEP 10c — RSeQC RNA-seq QC"
  "$REPO/scripts/run_rnaseq_qc.sh" \
    "$SAMPLESHEET" "$OUTDIR/bams" "$OUTDIR/07_qc" "$RSEQC_BED" \
    "${RSEQC_BIN_DIR:-}" "$LIBRARY_LAYOUT" "${MAX_JOBS:-8}"
  log "STEP 10c — MultiQC RSeQC"
  MQC_RSEQC=()
  for d in read_distribution junction_annotation junction_saturation genebody; do
    [[ -d "$OUTDIR/07_qc/rseqc/$d" ]] && MQC_RSEQC+=("$OUTDIR/07_qc/rseqc/$d")
  done
  [[ ${#MQC_RSEQC[@]} -gt 0 ]] && \
  "${MULTIQC_BIN:-multiqc}" "${MQC_RSEQC[@]}" \
    -n multiQC_rseqc -o "$OUTDIR/07_qc/multiqc" \
    --data-format tsv --export -q || true
else
  log "STEP 10c — RSeQC SKIPPED (RUN_RSEQC=false or RSEQC_BED not found)"
fi

# ── Step 11: DESeq2 normalization ─────────────────────────────────────────────
log "STEP 11 — DESeq2 normalization"
"${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_normalize.R" \
  --samplesheet "$SAMPLESHEET" --countdir "$OUTDIR/STARgeneCounts" \
  --gtf "$GTF" --layout "$LIBRARY_LAYOUT" \
  --outdir "$OUTDIR/analysis/counts" --design "${DESIGN_FORMULA:-~ condition}"

# ── Step 12: Normalize bedGraph ───────────────────────────────────────────────
log "STEP 12 — normalize_bedgraph.R"
"${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/normalize_bedgraph.R" \
  --samplesheet "$SAMPLESHEET" \
  --sffile "$OUTDIR/analysis/counts/size_factors.tsv" \
  --rawbgdir "$OUTDIR/bedGraph/raw" \
  --outdir "$OUTDIR/bedGraph/normalized" --layout "$LIBRARY_LAYOUT"

# ── Step 13: BigWig ───────────────────────────────────────────────────────────
log "STEP 13 — BigWig [species=$SPECIES naming=$CHROMOSOME_NAMING filter=$REGULAR_CHROMS_ONLY]"
for bg in "$OUTDIR/bedGraph/normalized/"*.bedGraph.gz; do
  [[ -f "$bg" ]] || continue
  submit "$REPO/scripts/norm_bedgraph_to_bigwig.sh \
    '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '$KENTUTILS_DIR'"
done; wait_all

# ── Step 14–15: Replicate merging ────────────────────────────────────────────
if [[ "${MERGE_REPLICATES:-true}" == "true" ]]; then
  log "STEP 14 — merge_bedgraph_replicates.R"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/merge_bedgraph_replicates.R" \
    --samplesheet "$SAMPLESHEET" --bgdir "$OUTDIR/bedGraph/normalized" \
    --outdir "$OUTDIR/bedGraph/merged" --genome "$GENOME_ASSEMBLY" \
    --layout "$LIBRARY_LAYOUT"
  log "STEP 15 — merged BigWigs"
  for bg in "$OUTDIR/bedGraph/merged/"*.bedGraph; do
    [[ -f "$bg" ]] || continue
    submit "$REPO/scripts/norm_bedgraph_to_bigwig.sh \
      '$bg' '$CHROM_SIZES' '$OUTDIR/bigwig' '$KENTUTILS_DIR'"
  done; wait_all
fi

# ── Step 16: DE ───────────────────────────────────────────────────────────────
if [[ "${RUN_DE:-true}" == "true" ]]; then
  log "STEP 16 — DESeq2 DE"
  "${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_de.R" \
    --countsrdata "$OUTDIR/analysis/counts/dds.RData" \
    --contrasts "${CONTRASTS:-config/contrasts.csv}" \
    --outdir "$OUTDIR/analysis/DE" \
    --padj "${PADJ_THRESHOLD:-0.05}" --lfc "${LFC_THRESHOLD:-1}"
fi

# ── Step 17: QC plots ─────────────────────────────────────────────────────────
log "STEP 17 — DESeq2 QC plots"
"${RSCRIPT_BIN:-Rscript}" "$REPO/scripts/Rscripts/deseq2_qc_plots.R" \
  --countsrdata "$OUTDIR/analysis/counts/dds.RData" \
  --outdir "$OUTDIR/analysis/figures"

# ── Step 18: UCSC tracks ──────────────────────────────────────────────────────
if [[ "${UCSC_TRACKS:-true}" == "true" ]]; then
  log "STEP 18 — UCSC tracks"
  "$REPO/scripts/create_ucsc_tracks.sh" \
    "$OUTDIR/bigwig" "$OUTDIR/reports/ucsc_tracks.txt" \
    "$OUTDIR/reports/bigwig_summary.txt" "${UCSC_BASE_URL:-PLACEHOLDER}"
fi

# ── Step 19: MultiQC final ────────────────────────────────────────────────────
log "STEP 19 — MultiQC final"
MQC_FINAL=("$OUTDIR/fastQC/raw" "$OUTDIR/fastQC/trimmed" \
           "$OUTDIR/STARlogs" "$OUTDIR/07_qc/star")
for d in read_distribution junction_annotation junction_saturation genebody; do
  [[ -d "$OUTDIR/07_qc/rseqc/$d" ]] && MQC_FINAL+=("$OUTDIR/07_qc/rseqc/$d")
done
"${MULTIQC_BIN:-multiqc}" "${MQC_FINAL[@]}" \
  -n multiQC_final -o "$OUTDIR/multiQC/final" --data-format tsv --export -q

# ── Step 20: HTML report ──────────────────────────────────────────────────────
log "STEP 20 — Pipeline report"
"${RSCRIPT_BIN:-Rscript}" -e "
  rmarkdown::render(
    input       = '$REPO/scripts/Rscripts/pipeline_report.Rmd',
    params      = list(outdir='$OUTDIR', config='$CONFIG',
                       samplesheet='$SAMPLESHEET', species='${SPECIES}',
                       layout='${LIBRARY_LAYOUT}'),
    output_file = '$OUTDIR/reports/pipeline_report.html',
    quiet=TRUE)
" || log "WARNING: pipeline report failed — check pandoc"

log "rnaseq2tracks v4 complete.  Results: $OUTDIR"
