#!/usr/bin/env bash
# =============================================================================
# check_strand_consistency.sh v4.1 — parallel per-sample strand check
# ORIGIN: NEW v4.0 — verify strand-split read counts vs total mapped
# CHANGES v4.1: parallel execution via submit/wait_all (same pattern as
#               run_rnaseq_qc.sh); results collected in samplesheet order
# Usage: check_strand_consistency.sh <samplesheet> <bamdir> <layout> \
#                                    [tolerance_pct=5] [max_jobs=8]
# =============================================================================
set -euo pipefail
SS="$1"; BAMDIR="$2"; LAYOUT="$3"; TOL="${4:-5}"; MAX_JOBS="${5:-8}"

declare -a SID STRAND
while IFS=',' read -r f1 f2 f3 f4 f5 f6 _rest; do
  [[ "$f1" =~ ^[[:space:]]*# || "$f1" == "sample_id" ]] && continue
  SID+=("$f1")
  [[ "$LAYOUT" == "PE" ]] && STRAND+=("$f6") || STRAND+=("$f5")
done < <(grep -v '^[[:space:]]*#' "$SS")

TMPDIR_RESULTS=$(mktemp -d)
trap 'rm -rf "$TMPDIR_RESULTS"' EXIT

# ── Job throttle ──────────────────────────────────────────────────────────────
declare -a _PIDS=()
submit() {
  while [[ ${#_PIDS[@]} -ge $MAX_JOBS ]]; do
    local live=()
    for p in "${_PIDS[@]}"; do kill -0 "$p" 2>/dev/null && live+=("$p"); done
    _PIDS=("${live[@]+"${live[@]}"}"); [[ ${#_PIDS[@]} -ge $MAX_JOBS ]] && sleep 1
  done
  eval "$@" &
  _PIDS+=($!)
}
wait_all() {
  local ok=0
  for p in "${_PIDS[@]+"${_PIDS[@]}"}"; do wait "$p" || ok=1; done
  _PIDS=()
  [[ $ok -eq 0 ]] || { echo "ERROR: a strand-check job failed" >&2; exit 1; }
}

# ── Per-sample function — writes one-line result to tmp file ──────────────────
check_sample() {
  local sid="$1" strand="$2" bam="$3" layout="$4" tol="$5" threads="$6" outf="$7"
  strand=$(echo "$strand" | tr -d '[:space:]\r')
  [[ "$strand" == "unstranded" ]] && { echo "SKIP $sid" > "$outf"; return 0; }
  [[ -f "$bam" ]] || { echo "WARN $sid BAM missing: $bam" > "$outf"; return 0; }

  TOTAL=$(samtools view -c -@ "$threads" -F 0x900 "$bam")
  if [[ "$layout" == "PE" ]]; then
    FWD=$(samtools  view -c -@ "$threads" -F 0x900 -f 0x50 "$bam")
    FWD2=$(samtools view -c -@ "$threads" -F 0x910 -f 0x80 "$bam")
    STRAND_FWD=$(( FWD + FWD2 ))
  else
    STRAND_FWD=$(samtools view -c -@ "$threads" -F 0x900 -f 0x10 "$bam")
  fi
  STRAND_REV=$(( TOTAL - STRAND_FWD ))
  SUM=$(( STRAND_FWD + STRAND_REV ))
  DIFF=$(( SUM - TOTAL )); [[ $DIFF -lt 0 ]] && DIFF=$(( -DIFF ))
  PCT=$(( TOTAL > 0 ? DIFF * 100 / TOTAL : 0 ))

  if [[ $PCT -gt $tol ]]; then
    echo "FAIL $sid Fwd=$STRAND_FWD Rev=$STRAND_REV Sum=$SUM Total=$TOTAL diff=${PCT}% > ${tol}%" > "$outf"
  else
    echo "OK   $sid Fwd=$STRAND_FWD Rev=$STRAND_REV Total=$TOTAL diff=${PCT}%" > "$outf"
  fi
}
export -f check_sample

# ── Submit all samples in parallel ───────────────────────────────────────────
# Divide samtools threads across parallel jobs to avoid CPU over-subscription
SAMTOOLS_T=$(( ${SAMTOOLS_THREADS:-4} > MAX_JOBS ? ${SAMTOOLS_THREADS:-4} / MAX_JOBS : 1 ))

echo "[check_strand_consistency.sh] Checking ${#SID[@]} samples (max_jobs=$MAX_JOBS, tol=${TOL}%)..."
for i in "${!SID[@]}"; do
  sid="${SID[$i]}"
  strand="${STRAND[$i]}"
  bam="$BAMDIR/${sid}_sortedS.bam"
  outf="$TMPDIR_RESULTS/${sid}.result"
  submit "check_sample '$sid' '$strand' '$bam' '$LAYOUT' '$TOL' '$SAMTOOLS_T' '$outf'"
done
wait_all

# ── Collect results in original samplesheet order ────────────────────────────
echo ""
FAIL=0
for i in "${!SID[@]}"; do
  sid="${SID[$i]}"
  outf="$TMPDIR_RESULTS/${sid}.result"
  [[ -f "$outf" ]] || { echo "  [WARN] No result file for $sid" >&2; continue; }
  result=$(cat "$outf")
  status="${result%% *}"
  msg="${result#* }"
  case "$status" in
    OK)   printf "  [ OK ] %s\n" "$msg" ;;
    FAIL) printf "  [FAIL] %s\n" "$msg" >&2; FAIL=$(( FAIL + 1 )) ;;
    WARN) printf "  [WARN] %s\n" "$msg" ;;
    SKIP) : ;;
  esac
done

echo ""
[[ $FAIL -eq 0 ]] || {
  echo "ERROR: Strand check failed for $FAIL sample(s). Check samplesheet strandedness column." >&2
  exit 1
}
echo "[check_strand_consistency.sh] All stranded samples passed."
