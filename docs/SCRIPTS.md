# Script reference

## Language decision per script

| Script | Language | Origin | Role |
|--------|----------|--------|------|
| `rnaseq2tracks.sh` | Bash | NEW | 20-step master orchestrator |
| `star_SE_single.sh` | Bash | ADAPTED | STAR SE alignment |
| `star_PE_single.sh` | Bash | NEW | STAR PE alignment |
| `trimgalore_single.sh` | Bash | ADAPTED | TrimGalore SE/PE |
| `bam_sort_index.sh` | Bash | NEW | samtools sort + index |
| `norm_bedgraph_to_bigwig.sh` | Bash | ADAPTED | Sort bedGraph → BigWig |
| `create_ucsc_tracks.sh` | Bash | ADAPTED | UCSC track file |
| `Rscripts/bam_to_bedgraph.R` | **R** | **NEW v2** | Strand-aware coverage (Rsamtools + GenomicAlignments + rtracklayer); replaces strand_split_bedgraph.sh |
| `Rscripts/deseq2_normalize.R` | **R** | ADAPTED | DESeq2 counts, SF, SF_rpm, FPKM, TPM |
| `Rscripts/normalize_bedgraph.R` | **R** | ADAPTED | rtracklayer SF_rpm scaling |
| `Rscripts/deseq2_de.R` | **R** | ADAPTED | DE + apeglm shrinkage |
| `Rscripts/deseq2_qc_plots.R` | **R** | ADAPTED | PCA, clustering, heatmaps |
| `Rscripts/merge_bedgraph_replicates.R` | **R** | ADAPTED | GRanges disjoin averaging |
| `Rscripts/pipeline_report.Rmd` | **R Markdown** | NEW | HTML report + SF_rpm chart |

## Why Bash stays for orchestration

Bash is idiomatic for:
- Background job fan-out with PID arrays and `kill -0` health checks
- Calling STAR, samtools, TrimGalore — tools that expect stdin/stdout/file paths
- Simple text file generation (UCSC tracks)

Using `processx::run()` or `system2()` inside an R master would add ~200 lines
of R boilerplate (error handling, stdout capture, parallel dispatch) to replace
~30 lines of clean Bash. No benefit.

## Why R handles coverage

Old approach (`strand_split_bedgraph.sh`):
1. `samtools view -b -F 16 bam > fwd.bam` (writes disk)
2. `samtools view -b -f 16 bam > rev.bam` (writes disk)
3. `bedtools genomecov -ibam fwd.bam -bg > fwd.bedGraph`
4. `bedtools genomecov -ibam rev.bam -bg > rev.bedGraph`

New approach (`bam_to_bedgraph.R`):
1. `readGAlignments/readGAlignmentPairs` → GAlignments object in memory
2. `coverage(ga[strand(ga)=="+"])` → Rle coverage
3. `export.bedGraph()` → file

No split BAMs. PE pair orientation handled correctly by `GAlignmentPairs`.
