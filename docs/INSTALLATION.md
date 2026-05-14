# Installation

## Conda (recommended)

```bash
conda env create -f environment.yml
conda activate rnaseq2tracks
```

## RSeQC

Included in `environment.yml`. Verify:
```bash
infer_experiment.py --version
read_distribution.py --version
```

If installed elsewhere, set `RSEQC_BIN_DIR` in config.

## RSeQC BED12 annotation

### Option 1: Download prebuilt

```bash
# hg38
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V42_Comprehensive.bed.gz
gunzip hg38_GENCODE_V42_Comprehensive.bed.gz

# mm39 — use mm10 file or generate (see option 2)
```

### Option 2: Generate from GTF

```bash
conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed
gtfToGenePred annotation.gtf annotation.genePred
genePredToBed annotation.genePred annotation.bed
```

## STAR index

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/star_index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/annotation.gtf.gz \
  --sjdbOverhang 149 \
  --runThreadN 16
```

## Chromosome sizes

```bash
fetchChromSizes mm39 > mm39.chrom.sizes
fetchChromSizes hg38 > hg38.chrom.sizes
```
