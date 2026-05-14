<p align="center">
  <h1 align="center">rnaseq2tracks</h1>
  <p align="center">End-to-end RNA-seq: raw FASTQ → counts → normalized BigWigs → differential expression</p>
  <p align="center">
    <img src="https://img.shields.io/badge/version-4.0-blue"/>
    <img src="https://img.shields.io/badge/language-Bash%20%7C%20R-informational"/>
    <img src="https://img.shields.io/badge/aligner-STAR-green"/>
    <img src="https://img.shields.io/badge/QC-RSeQC-orange"/>
    <img src="https://img.shields.io/badge/DE-DESeq2-red"/>
    <img src="https://img.shields.io/badge/layout-SE%20%7C%20PE-orange"/>
    <img src="https://img.shields.io/badge/species-human%20%7C%20mouse-lightblue"/>
    <img src="https://img.shields.io/badge/license-MIT-lightgrey"/>
  </p>
</p>

---

## Workflow

```mermaid
flowchart TD
    A([📁 Raw FASTQ\nSE or PE]) --> B[Step 0: preflight_check.sh\ntools · R pkgs · RSeQC · genome files]
    B --> C[FastQC raw]
    C --> D[MultiQC raw]
    A --> E[TrimGalore]
    E --> F[FastQC trimmed]
    F --> G[MultiQC trimmed]
    E --> H[STAR\n--quantMode GeneCounts]
    H --> I[samtools\nsort + index]
    H --> J[MultiQC\nalignment logs]
    H --> K[collect_star_qc.sh\nSTAR summary TSV]

    I --> L[⭐ bam_to_bedgraph.R\nRsamtools + GenomicAlignments]
    I --> M[⚙️ check_strand_consistency.sh\nFwd+Rev vs Total]
    I --> N[🔬 run_rnaseq_qc.sh\ninfer_experiment · read_distribution\njunction_annotation · junction_saturation\ngeneBody_coverage]
    K --> O[07_qc/star/]
    N --> P[07_qc/rseqc/]
    O --> Q[MultiQC RSeQC\n07_qc/multiqc/]
    P --> Q

    L --> R[⭐ deseq2_normalize.R\ncounts · SF · SF_rpm · FPKM · TPM]
    L --> S[⭐ normalize_bedgraph.R\nSF_rpm · rtracklayer]
    R --> S
    S --> T[norm_bedgraph_to_bigwig.sh\nchr filter · BigWig]
    S --> U[⭐ merge_bedgraph_replicates.R]
    U --> V[Merged BigWigs]
    R --> W[⭐ deseq2_de.R\nWald + apeglm LFC]
    R --> X[⭐ deseq2_qc_plots.R\nPCA · heatmaps]
    T --> Y[UCSC tracks]
    W --> Z([📊 Reports\nHTML · BigWigs · DE tables])
    X --> Z; Y --> Z; Q --> Z

    style A fill:#4a90d9,color:#fff
    style Z fill:#27ae60,color:#fff
    style B fill:#ffe0e0,stroke:#cc0000
    style M fill:#ffe0e0,stroke:#cc0000
    style N fill:#e8f4ff,stroke:#0066cc
    style Q fill:#e8f4ff,stroke:#0066cc
    style L fill:#fff3cd,stroke:#f0a500
    style R fill:#fff3cd,stroke:#f0a500
    style S fill:#fff3cd,stroke:#f0a500
    style U fill:#fff3cd,stroke:#f0a500
    style W fill:#fff3cd,stroke:#f0a500
    style X fill:#fff3cd,stroke:#f0a500
```

> ⭐ R &nbsp;|&nbsp; ⚙️ sanity check &nbsp;|&nbsp; 🔬 RSeQC &nbsp;|&nbsp; other = Bash

---

## Features

- **SE and PE** support; choice set in `config.conf`
- **Human and mouse** in one config — switch with `SPECIES=`
- **Strand-aware BigWig tracks** (forward / reverse) or unstranded
- **UCSC-compatible BigWigs** — canonical chromosomes only (UCSC or Ensembl naming)
- **Preflight check** — validates tools, R packages, RSeQC binaries, BED file, genome paths
- **STAR alignment summary TSV** from Log.final.out — included in MultiQC
- **RSeQC QC module** — infer_experiment, read_distribution, junction_annotation, junction_saturation, geneBody_coverage; integrated into MultiQC
- **Strand consistency sanity check** — hard fail if Fwd+Rev diverges from Total > tolerance
- **DESeq2 SF_rpm normalization** (size factor × mean-RPM anchor)
- **DESeq2 DE** — Wald + apeglm LFC shrinkage per contrast
- **Replicate merging** — GRanges disjoin mean BigWigs
- **HTML pipeline report** — STAR summary, size factors, infer_experiment table, output index
- **Executable smoke test** — 7 checks before full run

---

## Quick start

```bash
git clone https://github.com/MichalGd/rnaseq2tracks.git && cd rnaseq2tracks
conda env create -f environment.yml && conda activate rnaseq2tracks

cp config/config_template.conf  config/config.conf   # fill all paths
cp config/samplesheet_template_PE.csv config/samplesheet.csv
cp config/contrasts_template.csv config/contrasts.csv

bash tests/run_smoke_test.sh config/config.conf      # pre-flight
./scripts/rnaseq2tracks.sh   config/config.conf
```

---

## Input FASTQ naming (PE)

```
KO_12_1_1__ERR14875937_1.fq.gz   ← R1
KO_12_1_2__ERR14875937_2.fq.gz   ← R2
```
`sample_id = KO_12_1` in samplesheet.

---

## Samplesheet columns

| Column | PE | SE | Values |
|--------|----|----|--------|
| `sample_id` | ✓ | ✓ | unique, no spaces |
| `fastq_R1` | ✓ | ✓ | absolute path |
| `fastq_R2` | ✓ | — | absolute path |
| `condition` | ✓ | ✓ | `KO`, `WT`, … |
| `replicate` | ✓ | ✓ | `1`, `2`, `3` … |
| `strandedness` | ✓ | ✓ | `unstranded` / `forward` / `reverse` |

---

## Strandedness

| Value | STAR col | Library type |
|-------|---------|-------------|
| `unstranded` | 2 | Non-stranded |
| `forward` | 3 | Read 1 on RNA strand |
| `reverse` | 4 | dUTP / NEBNext Ultra II / TruSeq Stranded |

Run `infer_experiment.py` output (in `07_qc/rseqc/infer_experiment/`) to validate.

---

## Output tree

```
<OUTDIR>/
├── 07_qc/
│   ├── star/               STAR Log.final.out symlinks + summary TSV
│   ├── rseqc/              infer_experiment · read_distribution
│   │                       junction_annotation · junction_saturation · genebody
│   └── multiqc/            multiQC_rseqc.html
├── analysis/
│   ├── counts/             raw_counts.tsv · normalized_counts.tsv · size_factors.tsv · dds.RData
│   ├── DE/                 DE tables · volcano plots
│   └── figures/            PCA · clustering · heatmaps
├── bigwig/                 per-sample Fwd/Rev + merged BigWigs
├── multiQC/                raw · trimmed · alignments · final
└── reports/
    ├── pipeline_report.html
    └── ucsc_tracks.txt
```

---

## Key new config variables (v4)

```bash
RSEQC_BED_MOUSE="/path/to/mm39_GENCODE_vM31.bed"
RSEQC_BED_HUMAN="/path/to/hg38_GENCODE_V45.bed"
RSEQC_BIN_DIR=""        # empty = use PATH
RUN_RSEQC="true"
STRAND_TOLERANCE_PCT="5"
```

---

## Documentation

| File | Contents |
|------|---------|
| [RSEQC.md](docs/RSEQC.md) | Module descriptions, BED files, metric interpretation |
| [WORKFLOW.md](docs/WORKFLOW.md) | Full step table |
| [SCRIPTS.md](docs/SCRIPTS.md) | Origin tags for every file |
| [INSTALLATION.md](docs/INSTALLATION.md) | Conda, STAR index, BED generation |
| [USAGE.md](docs/USAGE.md) | Config reference, post-run |
| [OUTPUTS.md](docs/OUTPUTS.md) | Full output tree with column docs |
| [KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md) | STAR shared memory, apeglm, gzip |
| [GITHUB_UPLOAD.md](docs/GITHUB_UPLOAD.md) | 3-command upload + tagging |

---

## Citation

See [`CITATION.cff`](CITATION.cff)

## License

MIT © Michal Gdula
