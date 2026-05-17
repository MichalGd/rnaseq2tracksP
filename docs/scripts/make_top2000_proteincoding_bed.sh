#!/usr/bin/env bash
set -eo pipefail

SS="${1:?ERROR: provide samplesheet.csv as arg 1}"
COUNTDIR="${2:?ERROR: provide STAR counts directory as arg 2}"
GTF="${3:?ERROR: provide annotation.gtf as arg 3}"
BED="${4:?ERROR: provide annotation.bed as arg 4}"
OUTBED="${5:?ERROR: provide output.bed as arg 5}"
N="${6:-2000}"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "Inputs:"
echo "  samplesheet : $SS"
echo "  counts dir  : $COUNTDIR"
echo "  GTF         : $GTF"
echo "  BED         : $BED  ($(wc -l < "$BED") transcripts)"
echo "  output      : $OUTBED"
echo "  top-N       : $N"

# ── Step 1: Parse samplesheet ─────────────────────────────────────────────────
log "Step 1 — parsing samplesheet..."
FIRST_DATA=$(grep -v '^[[:space:]]*#\|^sample_id' "$SS" | head -1)
NFIELDS=$(echo "$FIRST_DATA" | awk -F',' '{print NF}')
if [[ $NFIELDS -ge 6 ]]; then
    STRAND_IDX=6; echo "  Layout detected: PE (strandedness in CSV column 6)"
else
    STRAND_IDX=5; echo "  Layout detected: SE (strandedness in CSV column 5)"
fi

grep -v '^[[:space:]]*#\|^sample_id' "$SS" | \
while IFS=',' read -r f1 f2 f3 f4 f5 f6 _rest; do
    SID="$f1"; [[ -z "$SID" ]] && continue
    STRAND=$(echo "$([[ $STRAND_IDX -eq 6 ]] && echo "$f6" || echo "$f5")" | tr -d '[:space:]\r')
    case "$STRAND" in
        unstranded) SCOL=2 ;; forward) SCOL=3 ;; reverse) SCOL=4 ;;
        *) SCOL=2; echo "  [WARN] Unknown strandedness '$STRAND' for $SID" ;;
    esac
    printf '%s\t%s\n' "$SID" "$SCOL"
done > "$TMP/sample_cols.tsv"

N_SAMPLES=$(wc -l < "$TMP/sample_cols.tsv")
echo "  Samples loaded: $N_SAMPLES"
awk '{printf "    %-30s strand_col=%s\n", $1, $2}' "$TMP/sample_cols.tsv"
[[ $N_SAMPLES -eq 0 ]] && { echo "ERROR: no samples parsed"; exit 1; }

# ── Step 2: Sum raw counts across all samples ─────────────────────────────────
log "Step 2 — summing raw counts across $N_SAMPLES samples..."
FIRST=1
while IFS=$'\t' read -r SID SCOL; do
    TABFILE="$COUNTDIR/${SID}_ReadsPerGene.out.tab"
    if [[ ! -f "$TABFILE" ]]; then echo "  [WARN] Missing: $TABFILE"; continue; fi
    awk -v col="$SCOL" 'NR>4 { print $1 "\t" $col }' "$TABFILE" \
        | sort -k1,1 > "$TMP/sample_current.tsv"
    if [[ $FIRST -eq 1 ]]; then
        cp "$TMP/sample_current.tsv" "$TMP/counts_sum.tsv"; FIRST=0
    else
        join -t$'\t' "$TMP/counts_sum.tsv" "$TMP/sample_current.tsv" \
            | awk -F'\t' '{print $1 "\t" $2+$3}' > "$TMP/counts_sum_new.tsv"
        mv "$TMP/counts_sum_new.tsv" "$TMP/counts_sum.tsv"
    fi
done < "$TMP/sample_cols.tsv"
echo "  Genes with summed counts: $(wc -l < "$TMP/counts_sum.tsv")"

# ── Step 3: Extract protein-coding gene IDs from GTF ─────────────────────────
log "Step 3 — extracting protein-coding gene IDs from GTF..."
awk '
/^#/ { next }
$3 == "gene" && /gene_type "protein_coding"/ {
    match($0, /gene_id "([^"]+)"/, arr)
    gene = arr[1]; sub(/\.[0-9]+$/, "", gene); print gene
}' "$GTF" | sort -u > "$TMP/pc_genes.txt"
echo "  Protein-coding genes: $(wc -l < "$TMP/pc_genes.txt")"

# ── Step 4: Rank and select top-N protein-coding genes ───────────────────────
log "Step 4 — ranking and selecting top-$N protein-coding genes..."
awk '
NR == FNR { pc[$1]=1; next }
{ gene=$1; sub(/\.[0-9]+$/, "", gene); if (gene in pc) print $2 "\t" gene }
' "$TMP/pc_genes.txt" "$TMP/counts_sum.tsv" \
    | sort -t$'\t' -k1,1rn \
    | awk -F'\t' '{print $2}' \
    | head -n "$N" > "$TMP/top_genes.txt" || true
echo "  Top-$N genes selected: $(wc -l < "$TMP/top_genes.txt")"
[[ $(wc -l < "$TMP/top_genes.txt") -eq 0 ]] && { echo "ERROR: no genes selected"; exit 1; }

# ── Step 5: Build ENST→ENSG map from GTF ─────────────────────────────────────
log "Step 5 — building transcript→gene map from GTF..."
awk '
/^#/ { next }
$3 == "transcript" && /gene_type "protein_coding"/ {
    match($0, /transcript_id "([^"]+)"/, t); match($0, /gene_id "([^"]+)"/, g)
    tx=t[1]; sub(/\.[0-9]+$/,"",tx); gene=g[1]; sub(/\.[0-9]+$/,"",gene)
    print tx "\t" gene
}' "$GTF" | sort -u > "$TMP/tx2gene.tsv"
echo "  Protein-coding transcripts mapped: $(wc -l < "$TMP/tx2gene.tsv")"

# ── Step 6: Select transcripts for top-N genes ───────────────────────────────
log "Step 6 — selecting transcripts for top genes..."
awk 'NR==FNR{keep[$1]=1;next} {if($2 in keep) print $1}' \
    "$TMP/top_genes.txt" "$TMP/tx2gene.tsv" | sort -u > "$TMP/tx_keep.txt"
echo "  Transcripts to keep: $(wc -l < "$TMP/tx_keep.txt")"

# ── Step 7: Filter BED12 ──────────────────────────────────────────────────────
log "Step 7 — filtering BED12..."
mkdir -p "$(dirname "$OUTBED")"
awk 'NR==FNR{keep[$1]=1;next}{tx=$4;sub(/\.[0-9]+$/,"",tx);if(tx in keep)print}' \
    "$TMP/tx_keep.txt" "$BED" > "$OUTBED"

# ── Step 8: Summary ───────────────────────────────────────────────────────────
N_ROWS=$(wc -l < "$OUTBED")
N_GENES_OUT=$(awk 'NR==FNR{m[$1]=$2;next}{tx=$4;sub(/\.[0-9]+$/,"",tx);g=m[tx];if(g!="")genes[g]=1}END{print length(genes)}' \
    "$TMP/tx2gene.tsv" "$OUTBED")
echo ""
log "✓ Complete"
printf "  Samples processed             : %d\n" "$N_SAMPLES"
printf "  Protein-coding genes ranked   : %d\n" "$(wc -l < "$TMP/pc_genes.txt")"
printf "  Top-%d genes selected         : %d\n" "$N" "$(wc -l < "$TMP/top_genes.txt")"
printf "  Output BED rows (transcripts) : %d\n" "$N_ROWS"
printf "  Unique genes in output BED    : %d\n" "$N_GENES_OUT"
printf "  Output file                   : %s\n" "$OUTBED"
[[ $N_ROWS -lt 100 ]] && echo "  [WARN] Very few rows — verify GTF/BED are the same genome version"
