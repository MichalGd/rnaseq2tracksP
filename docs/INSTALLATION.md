# Installation

## Conda (recommended)

```bash
conda env create -f environment.yml
conda activate rnaseq2tracks
```

This installs all tools and R packages including the enrichment analysis dependencies (clusterProfiler, ReactomePA, fgsea, msigdbr, org.Hs.eg.db, org.Mm.eg.db).

## Verify installation

```bash
bash tests/run_smoke_test.sh config/config.conf
```

## RSeQC

Included in `environment.yml`. Verify:
```bash
infer_experiment.py --version
read_distribution.py --version
```

If installed elsewhere, set `RSEQC_BIN_DIR` in `config.conf`.

## RSeQC BED12 annotation

### Option 1: Download prebuilt

```bash
# hg38 (GENCODE v42)
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V42_Comprehensive.bed.gz
gunzip hg38_GENCODE_V42_Comprehensive.bed.gz

# mm39
wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm39_GENCODE_M31_Comprehensive.bed.gz
gunzip mm39_GENCODE_M31_Comprehensive.bed.gz
```

### Option 2: Generate from GTF

```bash
conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed
gtfToGenePred annotation.gtf annotation.genePred
genePredToBed annotation.genePred annotation.bed
```

Set in `config.conf`:
```bash
RSEQC_BED_HUMAN="/path/to/hg38_GENCODE_V42_Comprehensive.bed"
RSEQC_BED_MOUSE="/path/to/mm39_GENCODE_M31_Comprehensive.bed"
```

## STAR index

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/star_index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/annotation.gtf \
  --sjdbOverhang 149 \
  --runThreadN 16
```

`sjdbOverhang` should be set to read length − 1.

## Chromosome sizes

```bash
fetchChromSizes hg38 > hg38.chrom.sizes
fetchChromSizes mm39 > mm39.chrom.sizes
```

Or generate from genome FASTA:
```bash
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > genome.chrom.sizes
```

## UCSC kentutils

`bedGraphToBigWig` is not available via conda. Download from:
```
https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
```

```bash
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x bedGraphToBigWig
```

Set `KENTUTILS_DIR` in `config.conf` to the directory containing the binary.
