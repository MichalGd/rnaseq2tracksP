#!/usr/bin/env bash
# ORIGIN: ADAPTED v1 — parameterised STAR SE wrapper
# Usage: star_SE_single.sh <index> <outdir> <sample_id> <R1.fq.gz> <threads> <tmpdir>
set -euo pipefail
INDEX="$1"; OUTDIR="$2"; SID="$3"; R1="$4"; THR="${5:-15}"; TMP="${6:-/tmp}"
STAR --genomeDir "$INDEX" \
  --readFilesIn "$R1" --readFilesCommand zcat \
  --outSAMattributes NH HI NM MD \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "${OUTDIR}/${SID}_" \
  --outMultimapperOrder Random --outSAMmultNmax 1 \
  --chimSegmentMin 15 \
  --quantMode GeneCounts \
  --genomeLoad LoadAndKeep \
  --outTmpDir "${TMP}/STAR_tmp_${SID}_$$" \
  --runThreadN "$THR"
rm -rf "${TMP}/STAR_tmp_${SID}_$$"
