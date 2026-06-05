#!/usr/bin/env bash
# ORIGIN: ADAPTED v1 from fastq2tracks create_ucsc_tracks.sh pattern
# Usage: create_ucsc_tracks.sh <bigwig_dir> <trackfile_out> <summary_out> <base_url>
set -euo pipefail
BWDIR="$1"; TRACKFILE="$2"; SUMMARY="$3"; BASEURL="$4"
COLORS=("220,50,32" "32,50,220" "50,180,50" "180,50,180" "50,180,180" "180,180,50")
i=0
> "$TRACKFILE"; > "$SUMMARY"
for bw in "$BWDIR"/*.bw; do
  [[ -f "$bw" ]] || continue
  NAME="$(basename "$bw" .bw)"
  C="${COLORS[$((i % ${#COLORS[@]}))]}"
  [[ "$NAME" == *Rev* || "$NAME" == *rev* ]] && \
    { echo "track type=bigWig name=\"$NAME\" description=\"$NAME\" bigDataUrl=${BASEURL}/$(basename "$bw") visibility=full color=$C negateValues=on" >> "$TRACKFILE"; } || \
    { echo "track type=bigWig name=\"$NAME\" description=\"$NAME\" bigDataUrl=${BASEURL}/$(basename "$bw") visibility=full color=$C" >> "$TRACKFILE"; }
  echo "$NAME  $bw" >> "$SUMMARY"
  i=$((i+1))
done
echo "[create_ucsc_tracks.sh] $i tracks written to $TRACKFILE"
