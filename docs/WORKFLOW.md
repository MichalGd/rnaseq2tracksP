# Workflow documentation

## Step table

| Step | Tool | Language | Notes |
|------|------|----------|-------|
| 1 | mkdir | Bash | Output tree |
| 2 | FastQC | Bash | Parallel raw reads |
| 3 | MultiQC | Bash | Raw QC |
| 4 | TrimGalore | Bash | SE or PE; `--basename` for PE renaming |
| 5 | FastQC | Bash | Trimmed reads |
| 6 | MultiQC | Bash | Trimmed QC |
| 7 | STAR | Bash | `--quantMode GeneCounts`; PE: `--peOverlapNbasesMin 10` |
| 8 | samtools | Bash | Sort + index BAM |
| 9 | MultiQC | Bash | Alignment logs |
| 10 | `bam_to_bedgraph.R` | R | Strand-aware coverage — Rsamtools + GenomicAlignments |
| 11 | `deseq2_normalize.R` | R | Counts, SF, SF_rpm, FPKM, TPM, dds.RData |
| 12 | `normalize_bedgraph.R` | R | rtracklayer SF_rpm scaling |
| 13 | `norm_bedgraph_to_bigwig.sh` | Bash | Chr filter → BigWig; keeps all_chrs intermediate |
| 14 | `merge_bedgraph_replicates.R` | R | GRanges disjoin mean (if MERGE_REPLICATES=true) |
| 15 | `norm_bedgraph_to_bigwig.sh` | Bash | Merged BigWigs |
| 16 | `deseq2_de.R` | R | Wald + apeglm (if RUN_DE=true) |
| 17 | `deseq2_qc_plots.R` | R | PCA, clustering, heatmaps |
| 18 | `create_ucsc_tracks.sh` | Bash | ucsc_tracks.txt |
| 19 | MultiQC | Bash | Final unified |
| 20 | `pipeline_report.Rmd` | R Markdown | HTML |

## Strandedness → STAR column

| Value | STAR col | Library type |
|-------|---------|-------------|
| `unstranded` | 2 | Standard non-stranded |
| `forward` | 3 | Read 1 on RNA strand |
| `reverse` | 4 | Read 2 on RNA strand (NEBNext Ultra II, TruSeq Stranded, dUTP) |

## Chromosome filter matrix (v3)

| SPECIES | CHROMOSOME_NAMING | Pattern retained |
|---------|-------------------|-----------------|
| human | ucsc | chr1–chr22, chrX, chrY, chrM |
| human | ensembl | 1–22, X, Y, MT |
| mouse | ucsc | chr1–chr19, chrX, chrY, chrM |
| mouse | ensembl | 1–19, X, Y, MT |
| any | (any) | All — when `REGULAR_CHROMS_ONLY=false` |
