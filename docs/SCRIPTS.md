# Script origin table

| Script | Language | Origin | v4 change |
|--------|----------|--------|-----------|
| `rnaseq2tracks.sh` | Bash | NEW v1 | Steps 0, 9b, 10b, 10c; RSeQC in Step 19 |
| `preflight_check.sh` | Bash | **NEW v4** | — |
| `run_rnaseq_qc.sh` | Bash | **ADAPTED v4** | From RNA_RSeQC_QuantSeqRev_17jan2024.sh + RSeQC_check17jan24.sh (MichalGd/3end-RNAseq-0.1) |
| `collect_star_qc.sh` | Bash | **NEW v4** | — |
| `check_strand_consistency.sh` | Bash | **NEW v4** | — |
| `star_SE_single.sh` | Bash | ADAPTED v1 | — |
| `star_PE_single.sh` | Bash | NEW v1 | — |
| `trimgalore_single.sh` | Bash | ADAPTED v1 | — |
| `bam_sort_index.sh` | Bash | NEW v1 | — |
| `norm_bedgraph_to_bigwig.sh` | Bash | ADAPTED v1→v3 | chr filter |
| `create_ucsc_tracks.sh` | Bash | ADAPTED v1 | — |
| `Rscripts/bam_to_bedgraph.R` | R | NEW v2 | — |
| `Rscripts/deseq2_normalize.R` | R | ADAPTED v1 | — |
| `Rscripts/normalize_bedgraph.R` | R | ADAPTED v2 | — |
| `Rscripts/deseq2_de.R` | R | ADAPTED v1 | — |
| `Rscripts/deseq2_qc_plots.R` | R | ADAPTED v1 | — |
| `Rscripts/merge_bedgraph_replicates.R` | R | ADAPTED v1 | — |
| `Rscripts/pipeline_report.Rmd` | R Markdown | NEW v1 | STAR + infer_experiment table |
| `tests/check_bash_syntax.sh` | Bash | VERBATIM | — |
| `tests/run_smoke_test.sh` | Bash | NEW v3→v4 | Added RSeQC check |

## run_rnaseq_qc.sh — adaptation notes

**Source scripts from `MichalGd/3end-RNAseq-0.1`:**
- `RNA_RSeQC_QuantSeqRev_17jan2024.sh`: hardcoded paths; `for` loop over `ls *.bam`; runs read_distribution + read_duplication only
- `RSeQC_check17jan24.sh`: single-BAM wrapper for read_distribution.py and read_duplication.py

**Changes in v4 run_rnaseq_qc.sh:**
1. All paths parameterised — zero hardcoded paths
2. Samplesheet-driven sample list
3. SE/PE aware
4. Replaced `read_duplication.py` with `junction_annotation.py` + `junction_saturation.py` (more informative for standard RNA-seq)
5. `geneBody_coverage.py` on merged BAM (original ran per-sample, very slow for large cohorts)
6. PID-array job throttling
7. Graceful SKIP per module when binary not found
