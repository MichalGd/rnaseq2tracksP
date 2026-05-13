# rnaseq2tracks 1.0

End-to-end RNA-seq workflow: raw FASTQ → count matrices → normalized BigWig tracks → differential expression.

## Language map

| Layer | Language | Rationale |
|-------|----------|-----------|
| Orchestration + external tools | **Bash** | STAR/samtools/TrimGalore/FastQC/MultiQC are shell-native; PID job throttling is idiomatic |
| Strand-aware coverage | **R** — Rsamtools + GenomicAlignments + rtracklayer | No intermediate split-BAMs; type-safe GRanges; correct PE orientation |
| bedGraph normalization | **R** — rtracklayer | Direct GRanges import/export; no awk/sed |
| Count matrix + normalization | **R** — DESeq2 | Native |
| Replicate merging | **R** — GenomicRanges disjoin | Exact replication of original mergeBedGraph logic |
| Differential expression | **R** — DESeq2 + apeglm | Native |
| QC plots | **R** — ggplot2, pheatmap | Native |
| Report | **R Markdown** | kableExtra + ggplot2 SF_rpm bar chart + sessionInfo |
| BigWig conversion | **Bash** — Kent utils | bedGraphToBigWig has no R equivalent |
| UCSC track file | **Bash** | Simple text generation |

## Quick start

```bash
git clone https://github.com/MichalGd/rnaseq2tracks.git && cd rnaseq2tracks
conda env create -f environment.yml && conda activate rnaseq2tracks
cp config/config_template.conf config/config.conf   # fill in paths
cp config/samplesheet_template_PE.csv config/samplesheet.csv
cp config/contrasts_template.csv config/contrasts.csv
bash tests/check_bash_syntax.sh
./scripts/rnaseq2tracks.sh config/config.conf
```

## Input FASTQ naming (PE example)
```
KO_12_1_1__ERR14875937_1.fq.gz   # R1
KO_12_1_2__ERR14875937_2.fq.gz   # R2 — set sample_id=KO_12_1 in samplesheet
```

## Samplesheet columns

| Column | PE | SE | Values |
|--------|----|----|--------|
| `sample_id` | ✓ | ✓ | unique, no spaces |
| `fastq_R1` | ✓ | ✓ | absolute path |
| `fastq_R2` | ✓ | — | absolute path |
| `condition` | ✓ | ✓ | e.g. KO, WT |
| `replicate` | ✓ | ✓ | 1, 2, 3 … |
| `strandedness` | ✓ | ✓ | `unstranded` / `forward` / `reverse` |

## Docs
[Workflow](docs/WORKFLOW.md) · [Scripts](docs/SCRIPTS.md) · [Installation](docs/INSTALLATION.md) · [Usage](docs/USAGE.md) · [Known Issues](docs/KNOWN_ISSUES.md) · [GitHub Upload](docs/GITHUB_UPLOAD.md)
