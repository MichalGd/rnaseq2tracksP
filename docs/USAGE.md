# Usage

## 1. Copy and fill config

```bash
cp config/config_template.conf config/config.conf
# Add config.conf to .gitignore — it contains your server paths
```

Required variables in config.conf:
- `SPECIES`, `GENOME_ASSEMBLY`
- `STAR_INDEX_HUMAN` / `STAR_INDEX_MOUSE`
- `GTF_HUMAN` / `GTF_MOUSE`
- `CHROM_SIZES_HUMAN` / `CHROM_SIZES_MOUSE`
- `RSEQC_BED_HUMAN` / `RSEQC_BED_MOUSE`
- `KENTUTILS_DIR`
- `OUTDIR`

## 2. Prepare samplesheet

Copy template and add one row per sample. For PE:
```
sample_id,fastq_R1,fastq_R2,condition,replicate,strandedness
KO_12_1,/data/KO_12_1_1__ERR14875937_1.fq.gz,/data/KO_12_1_2__ERR14875937_2.fq.gz,KO,1,reverse
```

## 3. Prepare contrasts

```
contrast_id,numerator,denominator
KO_vs_WT,KO,WT
```

## 4. Smoke test (recommended)

```bash
bash tests/run_smoke_test.sh config/config.conf
```

## 5. Run

```bash
./scripts/rnaseq2tracks.sh config/config.conf
```

## Post-run

```bash
# Free STAR shared memory (if genomeLoad=LoadAndKeep)
STAR --genomeDir $STAR_INDEX --genomeLoad Remove

# Check HTML report
firefox $OUTDIR/reports/pipeline_report.html

# Check RSeQC MultiQC
firefox $OUTDIR/07_qc/multiqc/multiQC_rseqc.html

# Copy BigWigs to web server for UCSC
rsync -av $OUTDIR/bigwig/*.bw user@webserver:/path/to/bigwig/
```
