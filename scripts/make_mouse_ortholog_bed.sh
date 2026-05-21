#!/usr/bin/env bash
# make_mouse_ortholog_bed.sh
# ──────────────────────────────────────────────────────────────────────────────
# Creates a mm39 BED12 annotation for RSeQC containing transcripts from mouse
# orthologs of the human top-1000 expressed protein-coding genes.
#
# Usage:
#   bash make_mouse_ortholog_bed.sh \
#     <mm39_full.bed> \
#     <mouse_ensg_list.txt> \
#     <mm39_annotation.gtf> \
#     <output.bed>
#
# Arguments:
#   mm39_full.bed        Full mm39 BED12 annotation (e.g. from GENCODE or RSeQC site)
#   mouse_ensg_list.txt  List of mouse ENSMUSG IDs to keep (one per line)
#   mm39_annotation.gtf  Matching mm39 GTF (same version as the BED)
#   output.bed           Output filtered BED12
#
# The script selects all transcripts belonging to the requested ENSMUSG IDs,
# preserving all isoforms (analogous to the human make_top2000_proteincoding_bed.sh).
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

BED_IN="${1:?ERROR: provide mm39 BED12 as arg 1}"
ENSG_LIST="${2:?ERROR: provide mouse ENSMUSG list as arg 2}"
GTF="${3:?ERROR: provide mm39 GTF as arg 3}"
BED_OUT="${4:?ERROR: provide output BED path as arg 4}"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "Inputs:"
echo "  mm39 BED12 : $BED_IN ($(wc -l < "$BED_IN") transcripts)"
echo "  ENSG list  : $ENSG_LIST ($(wc -l < "$ENSG_LIST") genes)"
echo "  GTF        : $GTF"
echo "  output     : $BED_OUT"

# ── Step 1: Build ENST→ENSG map from mm39 GTF ────────────────────────────────
log "Step 1 — building transcript→gene map from mm39 GTF..."
awk '
/^#/ { next }
$3 == "transcript" {
    match($0, /transcript_id "([^"]+)"/, t)
    match($0, /gene_id "([^"]+)"/, g)
    tx  = t[1]; sub(/\.[0-9]+$/, "", tx)
    gene = g[1]; sub(/\.[0-9]+$/, "", gene)
    print tx "\t" gene
}' "$GTF" | sort -u > "$TMP/tx2gene.tsv"
echo "  Transcripts mapped: $(wc -l < "$TMP/tx2gene.tsv")"

# ── Step 2: Select ENST IDs for target ENSG ──────────────────────────────────
log "Step 2 — selecting transcripts for target genes..."
sort -u "$ENSG_LIST" > "$TMP/ensg_target.txt"
awk 'NR==FNR{keep[$1]=1; next} ($2 in keep){print $1}' \
    "$TMP/ensg_target.txt" "$TMP/tx2gene.tsv" \
    | sort -u > "$TMP/enst_keep.txt"
echo "  Transcripts to keep: $(wc -l < "$TMP/enst_keep.txt")"
[[ $(wc -l < "$TMP/enst_keep.txt") -eq 0 ]] && \
    { echo "ERROR: no transcripts found — check GTF/ENSG list match mm39"; exit 1; }

# ── Step 3: Filter BED12 ─────────────────────────────────────────────────────
log "Step 3 — filtering BED12..."
mkdir -p "$(dirname "$BED_OUT")"
awk 'NR==FNR{keep[$1]=1; next} {tx=$4; sub(/\.[0-9]+$/, "", tx); if (tx in keep) print}' \
    "$TMP/enst_keep.txt" "$BED_IN" > "$BED_OUT"

# ── Step 4: Summary ──────────────────────────────────────────────────────────
N_TX=$(wc -l < "$BED_OUT")
N_GENES=$(awk 'NR==FNR{m[$1]=$2; next} {tx=$4; sub(/\.[0-9]+$/,"",tx); g=m[tx]; if(g!="") genes[g]=1} END{print length(genes)}' \
    "$TMP/tx2gene.tsv" "$BED_OUT")

echo ""
log "✓ Complete"
printf "  Target ENSG requested : %d\n" "$(wc -l < "$TMP/ensg_target.txt")"
printf "  Transcripts selected  : %d\n" "$(wc -l < "$TMP/enst_keep.txt")"
printf "  Output BED rows       : %d\n" "$N_TX"
printf "  Unique genes in BED   : %d\n" "$N_GENES"
printf "  Output file           : %s\n" "$BED_OUT"
[[ $N_TX -lt 100 ]] && echo "  [WARN] Very few rows — verify GTF and BED are same mm39 version (e.g. both GENCODE vM31)"
