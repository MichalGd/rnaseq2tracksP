# Workflow — step reference (v4)

## Step table

| Step | Script | Language | Origin | Notes |
|------|--------|----------|--------|-------|
| 0 | `preflight_check.sh` | Bash | NEW v4 | Tools, R pkgs, RSeQC, genome files |
| 1 | mkdir | Bash | — | Output tree incl. 07_qc/ |
| 2 | FastQC | Bash | — | Parallel raw reads |
| 3 | MultiQC | Bash | — | Raw QC |
| 4 | `trimgalore_single.sh` | Bash | ADAPTED v1 | SE or PE; `--basename` |
| 5 | FastQC | Bash | — | Trimmed |
| 6 | MultiQC | Bash | — | Trimmed QC |
| 7 | `star_SE_single.sh` / `star_PE_single.sh` | Bash | ADAPTED/NEW v1 | `--quantMode GeneCounts` |
| 8 | `bam_sort_index.sh` | Bash | NEW v1 | samtools sort + index |
| 9 | MultiQC | Bash | — | Alignment logs |
| 9b | `collect_star_qc.sh` | Bash | NEW v4 | STAR TSV + MultiQC symlinks |
| 10 | `bam_to_bedgraph.R` | R | NEW v2 | Rsamtools + GenomicAlignments |
| 10b | `check_strand_consistency.sh` | Bash | NEW v4 | Fwd+Rev vs Total sanity |
| 10c | `run_rnaseq_qc.sh` | Bash/RSeQC | ADAPTED v4 | 4 RSeQC modules + geneBody |
| 10c | MultiQC RSeQC | Bash | — | 07_qc/multiqc/ |
| 11 | `deseq2_normalize.R` | R | ADAPTED v1 | SF, SF_rpm, FPKM, TPM |
| 12 | `normalize_bedgraph.R` | R | ADAPTED v2 | rtracklayer SF_rpm |
| 13 | `norm_bedgraph_to_bigwig.sh` | Bash | ADAPTED v3 | Chr filter + all_chrs intermediate |
| 14 | `merge_bedgraph_replicates.R` | R | ADAPTED v1 | GRanges disjoin mean |
| 15 | `norm_bedgraph_to_bigwig.sh` | Bash | — | Merged BigWigs |
| 16 | `deseq2_de.R` | R | ADAPTED v1 | Wald + apeglm |
| 17 | `deseq2_qc_plots.R` | R | ADAPTED v1 | PCA, clustering, heatmaps |
| 18 | `create_ucsc_tracks.sh` | Bash | ADAPTED v1 | ucsc_tracks.txt |
| 19 | MultiQC final | Bash | — | All sources incl. RSeQC |
| 20 | `pipeline_report.Rmd` | R Markdown | NEW v1 | HTML report |

## Strandedness → STAR column

| Value | STAR col | Library |
|-------|---------|---------|
| `unstranded` | 2 | Non-stranded |
| `forward` | 3 | Read 1 on RNA strand |
| `reverse` | 4 | dUTP / NEBNext Ultra II / TruSeq Stranded |

## Chromosome filter

| SPECIES | CHROMOSOME_NAMING | Retained |
|---------|-------------------|---------|
| human | ucsc | chr1–22, chrX, chrY, chrM |
| human | ensembl | 1–22, X, Y, MT |
| mouse | ucsc | chr1–19, chrX, chrY, chrM |
| mouse | ensembl | 1–19, X, Y, MT |
