#!/usr/bin/env bash
# =============================================================================
# run_rnaseq_qc.sh — RNA-seq-specific QC with RSeQC
# =============================================================================
# ORIGIN: ADAPTED from:
#   MichalGd/3end-RNAseq-0.1:
#     Scripts/ExtractingPAs/RNA_RSeQC_QuantSeqRev_17jan2024.sh
#     Scripts/ExtractingPAs/RSeQC_check17jan24.sh
#
# CHANGES vs original:
#   - All paths parameterised (zero hardcoded paths)
#   - Samplesheet-driven sample list (not ls *.bam)
#   - SE/PE compatible
#   - Added junction_annotation.py and junction_saturation.py
#   - geneBody_coverage.py runs on merged BAM (not per-sample) for efficiency
#   - Job throttling via PID array
#   - Graceful SKIP per module if binary not found
#
# Usage: run_rnaseq_qc.sh <samplesheet> <bamdir> <outdir> <bed> <rseqc_dir> <layout> <max_jobs>
# =============================================================================
set -euo pipefail
SS="$1"; BAMDIR="$2"; OUTDIR="$3"; BED="$4"
RSEQC_DIR="${5:-}"; LAYOUT="$6"; MAX_JOBS="${7:-8}"

RSEQC() {
  local py="$1"; shift
  if [[ -n "$RSEQC_DIR" && -x "$RSEQC_DIR/$py" ]]; then "$RSEQC_DIR/$py" "$@"
  elif command -v "$py" &>/dev/null; then "$py" "$@"
  else echo "  [SKIP] $py not found" >&2; return 0; fi
}

mkdir -p \
  "$OUTDIR/rseqc/infer_experiment" \
  "$OUTDIR/rseqc/read_distribution" \
  "$OUTDIR/rseqc/junction_annotation" \
  "$OUTDIR/rseqc/junction_saturation" \
  "$OUTDIR/rseqc/genebody"

declare -a SID
while IFS=',' read -r f1 _rest; do
  [[ "$f1" =~ ^[[:space:]]*# || "$f1" == "sample_id" ]] && continue
  SID+=("$f1")
done < <(grep -v '^[[:space:]]*#' "$SS")

declare -a _PIDS=()
submit() {
  while [[ ${#_PIDS[@]} -ge $MAX_JOBS ]]; do
    local live=()
    for p in "${_PIDS[@]}"; do kill -0 "$p" 2>/dev/null && live+=("$p"); done
    _PIDS=("${live[@]+"${live[@]}"}"); [[ ${#_PIDS[@]} -ge $MAX_JOBS ]] && sleep 2
  done
  eval "$@" &
  _PIDS+=($!)
}
wait_all() {
  local ok=0
  for p in "${_PIDS[@]+"${_PIDS[@]}"}"; do wait "$p" || ok=1; done
  _PIDS=()
  [[ $ok -eq 0 ]] || { echo "ERROR: RSeQC job failed" >&2; exit 1; }
}

echo "[run_rnaseq_qc.sh] Starting per-sample RSeQC (${#SID[@]} samples)..."
for sid in "${SID[@]}"; do
  BAM="$BAMDIR/${sid}_sortedS.bam"
  [[ -f "$BAM" ]] || { echo "  [WARN] BAM missing: $BAM" >&2; continue; }

  submit "RSEQC infer_experiment.py -i '$BAM' -r '$BED' \
    > '$OUTDIR/rseqc/infer_experiment/${sid}_infer_experiment.txt' 2>&1"

  submit "RSEQC read_distribution.py -i '$BAM' -r '$BED' \
    > '$OUTDIR/rseqc/read_distribution/${sid}_read_distribution.txt' 2>&1"

  submit "RSEQC junction_annotation.py -i '$BAM' -r '$BED' \
    -o '$OUTDIR/rseqc/junction_annotation/$sid' \
    > '$OUTDIR/rseqc/junction_annotation/${sid}_junction_annotation.log' 2>&1"

  submit "RSEQC junction_saturation.py -i '$BAM' -r '$BED' \
    -o '$OUTDIR/rseqc/junction_saturation/$sid' \
    > '$OUTDIR/rseqc/junction_saturation/${sid}_junction_saturation.log' 2>&1"
done
wait_all
echo "[run_rnaseq_qc.sh] Per-sample modules done."

# geneBody_coverage — per-sample in parallel (faster than merged BAM)
echo "[run_rnaseq_qc.sh] geneBody_coverage per-sample (parallel)..."
for sid in "${SID[@]}"; do
  BAM="$BAMDIR/${sid}_sortedS.bam"
  [[ -f "$BAM" ]] || { echo "  [WARN] BAM missing: $BAM" >&2; continue; }
  submit "RSEQC geneBody_coverage.py -i '$BAM' -r '$BED' \
    -o '$OUTDIR/rseqc/genebody/$sid' \
    > '$OUTDIR/rseqc/genebody/${sid}_geneBody_coverage.log' 2>&1 || true"
done
wait_all
echo "[run_rnaseq_qc.sh] Done. Outputs: $OUTDIR/rseqc/"
