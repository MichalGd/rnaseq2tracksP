# Script reference — origin and language decisions

## File-level summary

| Script | Language | Origin | v3 change |
|--------|----------|--------|-----------|
| `rnaseq2tracks.sh` | Bash | NEW v1 | Species path resolution; env var export |
| `star_SE_single.sh` | Bash | ADAPTED | — |
| `star_PE_single.sh` | Bash | NEW v1 | — |
| `trimgalore_single.sh` | Bash | ADAPTED | — |
| `bam_sort_index.sh` | Bash | NEW v1 | — |
| `norm_bedgraph_to_bigwig.sh` | Bash | ADAPTED | CHROMOSOME_NAMING + REGULAR_CHROMS_ONLY + all_chrs intermediate |
| `create_ucsc_tracks.sh` | Bash | ADAPTED | — |
| `Rscripts/bam_to_bedgraph.R` | R | NEW v2 | — |
| `Rscripts/deseq2_normalize.R` | R | ADAPTED | — |
| `Rscripts/normalize_bedgraph.R` | R | ADAPTED | — |
| `Rscripts/deseq2_de.R` | R | ADAPTED | — |
| `Rscripts/deseq2_qc_plots.R` | R | ADAPTED | — |
| `Rscripts/merge_bedgraph_replicates.R` | R | ADAPTED | — |
| `Rscripts/pipeline_report.Rmd` | R Markdown | NEW v1 | Species param added |
| `tests/check_bash_syntax.sh` | Bash | VERBATIM | — |
| `tests/run_smoke_test.sh` | Bash | NEW v3 | 5-category executable test |
| `LICENSE` | — | NEW v3 | MIT |
| `CITATION.cff` | — | NEW v3 | CFF 1.2.0 |

## Why Bash stays for orchestration

PID-based parallel job throttling (`kill -0` + `wait`) is ~30 lines in Bash.
The equivalent using `processx::run()` in R would require ~200 lines of error
handling + parallel dispatch for zero functional gain. All external tools
(STAR, samtools, TrimGalore, FastQC, MultiQC, bedGraphToBigWig) are shell-native.

## Why R handles coverage (v2+)

Old Bash approach wrote 2 × split BAMs per stranded sample via `samtools view -F/-f 16`.
`bam_to_bedgraph.R` uses `readGAlignmentPairs()` → GRanges strand subsetting → `coverage()` →
`export.bedGraph()` — no split BAMs, no intermediate disk I/O, correct PE orientation.

## Chromosome filter (v3)

`norm_bedgraph_to_bigwig.sh` reads `SPECIES`, `CHROMOSOME_NAMING`, and
`REGULAR_CHROMS_ONLY` from the environment (set by `rnaseq2tracks.sh`).
Pattern generation is done in a `chr_pattern()` shell function supporting
all four combinations: `{ucsc,ensembl} × {human,mouse}`.
