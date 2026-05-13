#!/usr/bin/env bash
# ORIGIN: ADAPTED — RNA-seq/scripts/ShellScripts/STAR_SE_singleFileM.sh
# Changes: all paths parameterized; STAR tmp includes PID (parallel-safe).
# Usage: star_SE_single.sh <index> <outdir> <sample_id> <R1> <threads> <tmpdir>
set -euo pipefail
STAR --genomeDir "$1" --readFilesIn "$4" --readFilesCommand zcat \
  --outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$2/$3_" --outMultimapperOrder Random \
  --outSAMmultNmax 1 --chimSegmentMin 15 --quantMode GeneCounts \
  --genomeLoad LoadAndKeep --outTmpDir "${6:-/tmp}/STAR_${3}_$$" \
  --runThreadN "${5:-15}"
rm -rf "${6:-/tmp}/STAR_${3}_$$"
