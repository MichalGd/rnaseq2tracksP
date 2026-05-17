#!/usr/bin/env bash
# ORIGIN: NEW v4.0 — verify strand-split read counts vs total mapped
# Usage: check_strand_consistency.sh <samplesheet> <bamdir> <layout> [tolerance_pct]
set -euo pipefail
SS="$1"; BAMDIR="$2"; LAYOUT="$3"; TOL="${4:-5}"
declare -a SID STRAND
while IFS=',' read -r f1 f2 f3 f4 f5 f6 _rest; do
  [[ "$f1" =~ ^[[:space:]]*# || "$f1" == "sample_id" ]] && continue
  SID+=("$f1")
  [[ "$LAYOUT" == "PE" ]] && STRAND+=("$f6") || STRAND+=("$f5")
done < <(grep -v '^[[:space:]]*#' "$SS")
FAIL=0
for i in "${!SID[@]}"; do
  sid="${SID[$i]}"; strand="${STRAND[$i]}"
  [[ "$strand" == "unstranded" ]] && continue
  BAM="$BAMDIR/${sid}_sortedS.bam"
  [[ -f "$BAM" ]] || { echo "  [WARN] BAM missing, skip: $BAM"; continue; }
  TOTAL=$(samtools view -c -@ ${SAMTOOLS_THREADS:-4} -F 0x900 "$BAM")
  if [[ "$LAYOUT" == "PE" ]]; then
    FWD=$(samtools view -c -@ ${SAMTOOLS_THREADS:-4} -F 0x900 -f 0x50 "$BAM")
    FWD2=$(samtools view -c -@ ${SAMTOOLS_THREADS:-4} -F 0x910 -f 0x80 "$BAM")
    STRAND_FWD=$(( FWD + FWD2 ))
  else
    STRAND_FWD=$(samtools view -c -@ ${SAMTOOLS_THREADS:-4} -F 0x900 -f 0x10 "$BAM")
  fi
  STRAND_REV=$(( TOTAL - STRAND_FWD ))
  SUM=$(( STRAND_FWD + STRAND_REV ))
  DIFF=$(( SUM - TOTAL )); [[ $DIFF -lt 0 ]] && DIFF=$(( -DIFF ))
  PCT=$(( TOTAL > 0 ? DIFF * 100 / TOTAL : 0 ))
  if [[ $PCT -gt $TOL ]]; then
    echo "  [FAIL] $sid: Fwd=$STRAND_FWD Rev=$STRAND_REV Sum=$SUM Total=$TOTAL diff=${PCT}% > ${TOL}%" >&2
    FAIL=$((FAIL+1))
  else
    printf "  [ OK ] %s  Fwd=%d  Rev=%d  Total=%d  diff=%d%%\n" \
      "$sid" "$STRAND_FWD" "$STRAND_REV" "$TOTAL" "$PCT"
  fi
done
echo ""
[[ $FAIL -eq 0 ]] || { echo "ERROR: Strand check failed for $FAIL sample(s). Check samplesheet strandedness column." >&2; exit 1; }
echo "[check_strand_consistency.sh] All stranded samples passed."
