# Installation

```bash
conda env create -f environment.yml
conda activate rnaseq2tracks
```

## STAR index

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/star_index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/annotation.gtf.gz \
  --sjdbOverhang 149 --runThreadN 16
```

## Chrom sizes

```bash
fetchChromSizes mm39 > /path/to/mm39.chrom.sizes
```
