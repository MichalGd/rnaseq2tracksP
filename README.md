<p align="center">
  <h1 align="center">rnaseq2tracks</h1>
  <p align="center">
    End-to-end RNA-seq pipeline: raw FASTQ → count matrices → normalized BigWig tracks → differential expression
  </p>
  <p align="center">
    <img src="https://img.shields.io/badge/version-2.0-blue" alt="version"/>
    <img src="https://img.shields.io/badge/language-Bash%20%7C%20R-informational" alt="language"/>
    <img src="https://img.shields.io/badge/aligner-STAR-green" alt="STAR"/>
    <img src="https://img.shields.io/badge/DE-DESeq2-red" alt="DESeq2"/>
    <img src="https://img.shields.io/badge/tracks-UCSC%20compatible-blueviolet" alt="UCSC"/>
    <img src="https://img.shields.io/badge/layout-SE%20%7C%20PE-orange" alt="layout"/>
    <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="license"/>
  </p>
</p>

---

## Workflow overview

```mermaid
flowchart TD
    A([📁 Raw FASTQ\nSE or PE]) --> B[FastQC\nraw reads]
    B --> C[MultiQC\nraw QC report]
    A --> D[TrimGalore\nadapter & quality trim]
    D --> E[FastQC\ntrimmed reads]
    E --> F[MultiQC\ntrimmed QC report]
    D --> G[STAR alignment\n--quantMode GeneCounts]
    G --> H[samtools\nsort + index BAM]
    H --> I[MultiQC\nalignment QC]
    G --> J[⭐ DESeq2\ncount matrix\nSF · SF_rpm · FPKM · TPM]
    H --> K[⭐ bam_to_bedgraph.R\nRsamtools + GenomicAlignments\nstrand-aware coverage]
    J --> L[⭐ normalize_bedgraph.R\nSF_rpm scaling\nrtracklayer]
    K --> L
    L --> M[bedGraphToBigWig\ncanonical chrs only\n✅ UCSC-ready]
    L --> N[⭐ merge_bedgraph_replicates.R\nGenomicRanges disjoin mean]
    N --> O[bedGraphToBigWig\nmerged BigWigs]
    J --> P[⭐ DESeq2 DE\nWald test + apeglm LFC shrinkage]
    J --> Q[⭐ QC plots\nPCA · clustering · heatmaps]
    M --> R[UCSC custom\ntrack file]
    O --> R
    P --> S([📊 Results\ncounts · BigWigs · DE tables\nHTML report])
    Q --> S
    R --> S
    I --> S

    style A fill:#4a90d9,color:#fff,stroke:#2c5f8a
    style S fill:#27ae60,color:#fff,stroke:#1a7a40
    style J fill:#fff3cd,stroke:#f0a500
    style K fill:#fff3cd,stroke:#f0a500
    style L fill:#fff3cd,stroke:#f0a500
    style N fill:#fff3cd,stroke:#f0a500
    style P fill:#fff3cd,stroke:#f0a500
    style Q fill:#fff3cd,stroke:#f0a500
```

> ⭐ = implemented in **R** &nbsp;|&nbsp; All other steps: **Bash** shell

---

## Features

- **Single-end and paired-end** support — set one parameter in config
- **Strand-aware BigWig tracks** — forward and reverse strand per sample
- **UCSC-compatible BigWigs** — non-standard chromosomes and scaffolds filtered out automatically
- **DESeq2 SF_rpm normalization** — size-factor anchored to mean RPM; publication-ready scale
- **Differential expression** — Wald test + apeglm LFC shrinkage, volcano and MA plots per contrast
- **Replicate merging** — optional averaged BigWigs per condition for cleaner visualization
- **Full QC** — FastQC + MultiQC at three stages; PCA, sample clustering, top-50 heatmaps
- **Reproducibility** — `sessionInfo.txt` written by every R module
- **HTML pipeline report** — kableExtra tables + SF_rpm bar chart; self-contained

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/MichalGd/rnaseq2tracks.git
cd rnaseq2tracks

# 2. Environment
conda env create -f environment.yml
conda activate rnaseq2tracks

# 3. Configure
cp config/config_template.conf  config/config.conf    # fill in paths
cp config/samplesheet_template_PE.csv config/samplesheet.csv
cp config/contrasts_template.csv config/contrasts.csv

# 4. Run
./scripts/rnaseq2tracks.sh config/config.conf
```

---

## Input FASTQ naming

Paired-end (supported format):
```
KO_12_1_1__ERR14875937_1.fq.gz   ← R1
KO_12_1_2__ERR14875937_2.fq.gz   ← R2
```
Set `sample_id=KO_12_1` in the samplesheet. TrimGalore `--basename` handles renaming automatically.

---

## Samplesheet

| Column | PE | SE | Values |
|---|---|---|---|
| `sample_id` | ✓ | ✓ | unique, no spaces |
| `fastq_R1` | ✓ | ✓ | absolute path |
| `fastq_R2` | ✓ | — | absolute path |
| `condition` | ✓ | ✓ | e.g. `KO`, `WT` |
| `replicate` | ✓ | ✓ | `1`, `2`, `3` … |
| `strandedness` | ✓ | ✓ | `unstranded` / `forward` / `reverse` |

Mixed strandedness across samples is supported — the `strandedness` column is per-sample.

---

## Strandedness guide

| Value | STAR column | Typical library |
|---|---|---|
| `unstranded` | 2 | Standard non-stranded |
| `forward` | 3 | 1st read on RNA strand |
| `reverse` | 4 | 2nd read on RNA strand *(NEBNext Ultra II, TruSeq Stranded, dUTP)* |

---

## Outputs

```
<OUTDIR>/
├── analysis/
│   ├── counts/          raw_counts.tsv · normalized_counts.tsv · fpkm · tpm · size_factors.tsv · dds.RData
│   ├── DE/              <contrast>_DE_results.tsv · volcano.pdf · MA_plot.pdf
│   └── figures/         PCA.pdf · sample_clustering.pdf · top50_heatmap.pdf
├── bigwig/
│   ├── <sample>_FwdS.bw          per-sample forward strand  ┐
│   ├── <sample>_RevS.bw          per-sample reverse strand  ┘ always produced
│   ├── <condition>_Fwd_mergedS.bw  replicate-merged forward  ┐
│   └── <condition>_Rev_mergedS.bw  replicate-merged reverse  ┘ if MERGE_REPLICATES=true
├── multiQC/             raw · trimmed · alignments · final
├── STARlogs/
├── STARgeneCounts/
└── reports/
    ├── pipeline_report.html
    └── ucsc_tracks.txt
```

---

## Key configuration variables

```bash
STAR_INDEX=""          # STAR genome index directory
GTF=""                 # Genome GTF (gzip OK)
CHROM_SIZES=""         # chrom.sizes (fetchChromSizes)
KENTUTILS_DIR=""       # directory containing bedGraphToBigWig
LIBRARY_LAYOUT="PE"    # SE | PE
MERGE_REPLICATES="true"
RUN_DE="true"
```

Full template: [`config/config_template.conf`](config/config_template.conf)

---

## Language map

| Step | Language | Why |
|---|---|---|
| Orchestration, STAR, samtools, TrimGalore, FastQC, MultiQC | **Bash** | Shell-native tools; PID job throttling idiomatic in Bash |
| Strand-aware coverage | **R** — Rsamtools + GenomicAlignments + rtracklayer | No intermediate split-BAMs; correct PE pair orientation |
| bedGraph normalization | **R** — rtracklayer | Type-safe GRanges import/export |
| Count matrix + normalization | **R** — DESeq2 | Native |
| Replicate merging | **R** — GenomicRanges | Exact `disjoin` + mean score logic |
| Differential expression | **R** — DESeq2 + apeglm | Native |
| QC plots | **R** — ggplot2, pheatmap | Native |
| Report | **R Markdown** | kableExtra + ggplot2 |
| BigWig conversion | **Bash** — Kent utils | No R equivalent for bedGraphToBigWig |

---

## Documentation

| Doc | Contents |
|---|---|
| [WORKFLOW.md](docs/WORKFLOW.md) | Step-by-step table, strandedness guide, SF_rpm formula |
| [SCRIPTS.md](docs/SCRIPTS.md) | Per-script origin tags (ADAPTED / NEW), language rationale |
| [INSTALLATION.md](docs/INSTALLATION.md) | Conda setup, STAR index, chrom.sizes |
| [USAGE.md](docs/USAGE.md) | Config reference, post-run commands |
| [OUTPUTS.md](docs/OUTPUTS.md) | Full output directory tree with column descriptions |
| [KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md) | Memory, shared memory, apeglm, pandoc, Seqinfo |
| [GITHUB_UPLOAD.md](docs/GITHUB_UPLOAD.md) | Git Bash upload instructions |

---

## License

MIT © Michal Gdula
