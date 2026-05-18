# Usage

## 1. Copy and fill config

```bash
cp config/config_template.conf config/config.conf
# config.conf is in .gitignore — it contains your server paths
```

Required variables:
- `SPECIES`, `GENOME_ASSEMBLY`
- `STAR_INDEX_HUMAN` / `STAR_INDEX_MOUSE`
- `GTF_HUMAN` / `GTF_MOUSE`
- `CHROM_SIZES_HUMAN` / `CHROM_SIZES_MOUSE`
- `RSEQC_BED_HUMAN` / `RSEQC_BED_MOUSE`
- `KENTUTILS_DIR`
- `OUTDIR`
- `CONTRASTS` (path to contrasts.csv)

Enrichment-specific variables (optional, have defaults):
- `PADJ_THRESHOLD` — ORA gene list adjusted p-value cutoff (default `0.05`)
- `LFC_THRESHOLD` — ORA gene list |log2FC| cutoff (default `1`)
- `ENRICHMENT_MINGS` — minimum gene set size (default `10`)
- `ENRICHMENT_MAXGS` — maximum gene set size (default `500`)

## 2. Prepare samplesheet

Copy template and add one row per sample. For PE:
```
sample_id,fastq_R1,fastq_R2,condition,replicate,strandedness
KO_1,/data/KO_1_R1.fq.gz,/data/KO_1_R2.fq.gz,KO,1,reverse
WT_1,/data/WT_1_R1.fq.gz,/data/WT_1_R2.fq.gz,WT,1,reverse
```

See `examples/samplesheet_example_PE.csv` and `examples/samplesheet_example_SE.csv`.

## 3. Prepare contrasts

```
contrast_id,numerator,denominator
KO_vs_WT,KO,WT
```

See `examples/contrasts_example.csv`.

## 4. Smoke test (recommended)

```bash
bash tests/run_smoke_test.sh config/config.conf
```

This checks bash syntax, all R packages (including enrichment packages), core tools, UCSC kentutils, RSeQC binaries, config variables, samplesheet, and contrasts file.

## 5. Run

```bash
./scripts/rnaseq2tracks.sh config/config.conf
```

## Rerunning specific steps

Steps use sentinel files and output-existence checks to skip completed work.

```bash
# Rerun DE (Step 16) and enrichment (Step 21) only
rm -f $OUTDIR/analysis/DE/*_DE_results.tsv
rm -f $OUTDIR/analysis/enrichment/.enrichment_done
./scripts/rnaseq2tracks.sh config/config.conf

# Rerun enrichment only (Step 21)
rm -f $OUTDIR/analysis/enrichment/.enrichment_done
./scripts/rnaseq2tracks.sh config/config.conf

# Rerun with custom thresholds (LFC=0 for broader sensitivity)
DE_LFC_THRESHOLD=0 DE_PADJ_THRESHOLD=0.05 \
PADJ_THRESHOLD=0.05 LFC_THRESHOLD=0 \
./scripts/rnaseq2tracks.sh config/config.conf

# Force rerun of all steps
FORCE_RERUN=1 ./scripts/rnaseq2tracks.sh config/config.conf
```

## Replicate subset analysis

To re-run DE on a subset of replicates (e.g., replicates 1 and 2 only):

```bash
bash scripts/rerun_deseq_rep12.sh config/config.conf
```

Output goes to `analysis_rep12/` inside `OUTDIR`.

## Post-run

```bash
# Free STAR shared memory (if genomeLoad=LoadAndKeep)
STAR --genomeDir $STAR_INDEX --genomeLoad Remove

# Check HTML pipeline report
firefox $OUTDIR/reports/pipeline_report.html

# Check RSeQC MultiQC
firefox $OUTDIR/07_qc/multiqc/multiQC_rseqc.html

# Browse enrichment results
ls $OUTDIR/analysis/enrichment/

# Copy BigWigs to web server for UCSC
rsync -av $OUTDIR/bigwig/*.bw user@webserver:/path/to/bigwig/
```
