# Changelog

## [1.0] — 2026-05-12

### Language architecture
- **Bash**: orchestration, STAR, samtools, TrimGalore, FastQC, MultiQC, Kent utils, UCSC
- **R**: all analytical work — coverage, normalization, DE, QC plots, report

### Key improvement over v1
- `strand_split_bedgraph.sh` (samtools view flag-filter → bedtools genomecov) replaced by
  `bam_to_bedgraph.R` (Rsamtools + GenomicAlignments + rtracklayer)
  — no intermediate strand-split BAMs; correct PE pair orientation
- `normalize_bedgraph.R`: data.table fread/fwrite → rtracklayer import/export
- All R scripts write `sessionInfo_*.txt`
- `pipeline_report.Rmd`: kableExtra + ggplot2 SF_rpm bar chart
