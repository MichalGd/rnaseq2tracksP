#!/usr/bin/env bash
# ORIGIN: ADAPTED — fastq2tracks/scripts/create_ucsc_tracks.1.0.sh
# Changes: all paths parameterized; strand-aware colour coding.
# Usage: create_ucsc_tracks.sh <bigwig_dir> <tracks_out> <summary_out> <base_url>
set -euo pipefail
mkdir -p "$(dirname "$2")" "$(dirname "$3")"
> "$2"; printf "track_num\tfilename\tsample_name\tfile_size_MB\n" > "$3"
N=0
for bw in "$1"/*.bw; do
  [[ -f "$bw" ]] || continue; N=$((N+1))
  FNAME="$(basename "$bw")"; SAMP="${FNAME%.bw}"; MB=$(du -m "$bw" | cut -f1)
  [[ "$FNAME" == *_Fwd* ]] && COL="0,0,200" || ([[ "$FNAME" == *_Rev* ]] && COL="200,0,0" || COL="0,0,0")
  printf 'track type=bigWig name="%s" description="%s" bigDataUrl=%s/%s visibility=full color=%s autoScale=off alwaysZero=on gridDefault=on graphType=bar windowingFunction=mean viewLimits=0:5 maxHeightPixels=200:50:10\n' \
    "$SAMP" "$SAMP" "${4:-PLACEHOLDER_URL}" "$FNAME" "$COL" >> "$2"
  printf "%d\t%s\t%s\t%s\n" "$N" "$FNAME" "$SAMP" "$MB" >> "$3"
done
echo "Generated $N track definitions -> $2"
