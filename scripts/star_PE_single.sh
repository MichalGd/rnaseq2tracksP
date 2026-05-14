#!/usr/bin/env bash
# ORIGIN: NEW v1 — parameterised STAR PE wrapper
# Usage: star_PE_single.sh <index> <outdir> <sample_id> <R1> <R2> <threads> <tmpdir>
set -euo pipefail
INDEX="$1"; OUTDIR="$2"; SID="$3"; R1="$4"; R2="$5"; THR="${6:-15}"; TMP="${7:-/tmp}"
STAR --genomeDir "$INDEX" \
  --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --outSAMattributes NH HI NM MD \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "${OUTDIR}/${SID}_" \
  --outMultimapperOrder Random --outSAMmultNmax 1 \
  --chimSegmentMin 15 \
  --quantMode GeneCounts \
  --genomeLoad LoadAndKeep \
  --peOverlapNbasesMin 10 \
  --outTmpDir "${TMP}/STAR_tmp_${SID}_$$" \
  --runThreadN "$THR"
rm -rf "${TMP}/STAR_tmp_${SID}_$$"
