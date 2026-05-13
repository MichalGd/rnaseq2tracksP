#!/usr/bin/env bash
# ORIGIN: NEW (mirrors star_SE_single.sh; adds PE flags)
# Usage: star_PE_single.sh <index> <outdir> <sample_id> <R1> <R2> <threads> <tmpdir>
set -euo pipefail
STAR --genomeDir "$1" --readFilesIn "$4" "$5" --readFilesCommand zcat \
  --outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$2/$3_" --outMultimapperOrder Random \
  --outSAMmultNmax 1 --chimSegmentMin 15 --quantMode GeneCounts \
  --genomeLoad LoadAndKeep --outTmpDir "${7:-/tmp}/STAR_${3}_$$" \
  --runThreadN "${6:-15}" --peOverlapNbasesMin 10
rm -rf "${7:-/tmp}/STAR_${3}_$$"
