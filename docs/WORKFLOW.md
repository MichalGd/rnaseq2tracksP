# Workflow — step reference (v5.1)

## Step table

| Step | Script | Language | Origin | Notes |
|------|--------|----------|--------|-------|
| 0 | `preflight_check.sh` | Bash | NEW v4 | Tools, R pkgs, RSeQC, FastQ Screen conf + indexes, genome files |
| 1 | mkdir | Bash | — | Output tree including 07_qc/, fastQScreen/, and analysis/enrichment/ |
| 2 | FastQC | Bash | — | Parallel raw reads |
| 3 | MultiQC | Bash | — | Raw QC |
| 2b | `fastq_screen` | Bash | **NEW v5** | Species swap + mycoplasma screen on raw reads; bowtie2 against 5-species panel |
| 4 | `trimgalore_single.sh` | Bash | ADAPTED v1 | SE or PE; `--basename` |
| 5 | FastQC | Bash | — | Trimmed |
| 6 | MultiQC | Bash | — | Trimmed QC |
| 7 | `star_SE_single.sh` / `star_PE_single.sh` | Bash | ADAPTED/NEW v1 | `--quantMode GeneCounts` |
| 8 | `bam_sort_index.sh` | Bash | NEW v1 | samtools sort + index |
| 9 | MultiQC | Bash | — | Alignment logs |
| 9b | `collect_star_qc.sh` | Bash | NEW v4 | STAR TSV + MultiQC symlinks |
| 10 | `bam_to_bedgraph.R` | R | NEW v2 | Rsamtools + GenomicAlignments; parallel per sample |
| 10b | `check_strand_consistency.sh` | Bash | NEW v4 | Fwd+Rev vs Total sanity check |
| 10c | `run_rnaseq_qc.sh` | Bash/RSeQC | ADAPTED v4 | 4 RSeQC modules + geneBody coverage |
| 10c | MultiQC RSeQC | Bash | — | 07_qc/multiqc/ |
| 11 | `deseq2_normalize.R` | R | ADAPTED v1 | SF, SF_rpm, FPKM, TPM |
| 12 | `normalize_bedgraph.R` | R | ADAPTED v2 | rtracklayer SF_rpm; parallel per sample |
| 13 | `norm_bedgraph_to_bigwig.sh` | Bash | ADAPTED v3 | Chr filter + all_chrs intermediate |
| 14 | `merge_bedgraph_replicates.R` | R | ADAPTED v1 | GRanges disjoin mean; parallel per condition |
| 15 | `norm_bedgraph_to_bigwig.sh` | Bash | — | Merged BigWigs |
| 16 | `deseq2_de.R` | R | ADAPTED v4.2 | Wald + apeglm/ashr LFC shrinkage; unshrunken + shrunken volcano and MA plots; annotated count tables |
| 17 | `deseq2_qc_plots.R` | R | ADAPTED v1 | PCA, clustering, heatmaps |
| 18 | `create_ucsc_tracks.sh` | Bash | ADAPTED v1 | ucsc_tracks.txt |
| 19 | MultiQC final | Bash | — | All sources including RSeQC and FastQ Screen |
| 20 | `pipeline_report.Rmd` | R Markdown | NEW v1 | HTML report |
| 21 | `deseq2_enrichment.R` | R | **NEW v5** | ORA + GSEA: GO BP/MF/CC, KEGG, Reactome, MSigDB Hallmarks |
| 22 | `cleanup_intermediates()` (inline) + `cleanup_existing_run.sh` | Bash | **NEW v5.1** | Removes large regenerable intermediate files after confirmed pipeline success. Gated on 3 completion sentinels. Controlled by `CLEANUP_INTERMEDIATES`, `CLEANUP_DRYRUN`, `CLEANUP_ALLCHR_BEDGRAPH` in config. See [`docs/CLEANUP.md`](CLEANUP.md) |

---

## Step 2b — FastQ Screen detail

Step 2b runs immediately after raw FastQC (Step 2) and **before** TrimGalore (Step 4). It screens a random subset of raw reads against a multi-species bowtie2 panel.

**Purpose:** detect two classes of pre-alignment quality problems that are invisible to FastQC:
- **Species swap** — a sample mapping predominantly to the wrong organism (e.g. human reads in a mouse experiment)
- **Mycoplasma contamination** — a common cell culture contaminant that can account for up to a few percent of reads

**Reference panel:**

| Database | Genome | Version |
|---|---|---|
| Mouse | GRCm39 | GENCODE M31 |
| Human | GRCh38 | GENCODE v42 |
| Zebrafish | GRCz11 | Ensembl 112 |
| Drosophila | BDGP6.46 | Ensembl 112 |
| Mycoplasma | Combined 8-species | NCBI RefSeq |

**Config parameters:**

| Variable | Default | Description |
|---|---|---|
| `FASTQSCREEN_CONF` | `config/fastq_screen.conf` | Path to database config |
| `FASTQSCREEN_THREADS` | `4` | bowtie2 threads per sample |
| `FASTQSCREEN_SUBSET` | `200000` | Reads sampled per file (0 = all) |

**Step behaviour:** skips gracefully (warning logged) if `fastq_screen` is absent from PATH or the conf file is missing. The pipeline continues to Step 4 in that case.

**Action thresholds:**

| Observation | Threshold | Action |
|---|---|---|
| Non-expected species | > 5 % | Warning — check sample provenance |
| Non-expected species | > 20 % | Fail — exclude sample |
| Mycoplasma | > 0.5 % | Warning — notify cell culture team |
| Mycoplasma | > 2 % | Fail — confirmed contamination, exclude sample |

**Output:** `fastQScreen/<sample>_screen.txt|html|png` per sample. All results are included automatically in the final MultiQC report (Step 19) as a stacked bar chart panel.

See [`docs/FASTQSCREEN.md`](FASTQSCREEN.md) for full setup instructions including bowtie2 index building.

---

## Step 21 — Gene enrichment analysis detail

Step 21 runs after Step 16 (DE) and processes each contrast independently.

**Input:** `analysis/DE/<contrast_id>_DE_results.tsv`

**Gene ID conversion:** Ensembl IDs (version suffix stripped) → Entrez IDs via `bitr()` (org.Hs.eg.db / org.Mm.eg.db). Multi-mapping is deduplicated by keeping the first match per Ensembl ID.

**ORA gene list:** genes with `padj < PADJ_THRESHOLD` and `|log2FC| > LFC_THRESHOLD`. Background: all expressed genes with a valid Entrez mapping and non-NA padj.

**GSEA ranking metric:** `sign(log2FoldChange) × −log10(padj + 1e−300)` — robust to ashr LFC shrinkage, avoids tied pile-up at zero. Genes are deduplicated by keeping the highest `|rank_metric|` per Entrez ID.

**Databases and methods:**

| Database | ORA | GSEA | pvalueCutoff (ORA) | pvalueCutoff (GSEA) |
|----------|-----|------|--------------------|---------------------|
| GO Biological Process | ✓ | ✓ | 0.25 | 0.25 |
| GO Molecular Function | ✓ | ✓ | 0.05 | 0.25 |
| GO Cellular Component | ✓ | — | 0.05 | — |
| KEGG | ✓ | ✓ | 0.05 | 0.25 |
| Reactome | ✓ | ✓ | 0.05 | 0.25 |
| MSigDB Hallmarks | — | ✓ | — | 0.05 (plot only) |

**Output:** `analysis/enrichment/<contrast_id>/` — TSV result tables and PDF/PNG plots (dotplot, barplot, cnetplot for ORA; dotplot and barplot for GSEA; NES barplot for Hallmarks).

**Step caching:** Completion is tracked by `analysis/enrichment/.enrichment_done`. Delete this file to rerun Step 21 without rerunning the full pipeline.

---

## Step 22 — Post-run cleanup detail

Step 22 runs after Step 21 and removes large intermediate files that can be fully regenerated from raw FASTQs and pipeline scripts. It is **disabled by default** (`CLEANUP_INTERMEDIATES=0`).

**Sentinels — all must exist before any file is deleted:**

| Sentinel file | Created at step |
|---|---|
| `multiQC/final/multiQC_final.html` | Step 19 |
| `reports/pipeline_report.html` | Step 20 |
| `analysis/enrichment/.enrichment_done` | Step 21 |

At least one `bigwig/*.bw` file must also be present. If any sentinel is missing, cleanup is skipped with a warning.

**Config parameters:**

| Variable | Default | Description |
|---|---|---|
| `CLEANUP_INTERMEDIATES` | `0` | `1` = enable cleanup |
| `CLEANUP_DRYRUN` | `0` | `1` = print what would be deleted without deleting |
| `CLEANUP_ALLCHR_BEDGRAPH` | `0` | `1` = also remove `bigwig/*.all_chromosomes.bedGraph.gz` |

See [`docs/CLEANUP.md`](CLEANUP.md) for full details, kept/removed file lists, and regeneration instructions.

---

## Strandedness → STAR column

| Value | STAR col | Library |
|-------|---------|---------|
| `unstranded` | 2 | Non-stranded |
| `forward` | 3 | Read 1 on RNA strand |
| `reverse` | 4 | dUTP / NEBNext Ultra II / TruSeq Stranded |

---

## Chromosome filter

| SPECIES | CHROMOSOME_NAMING | Retained |
|---------|-------------------|---------| 
| human | ucsc | chr1–22, chrX, chrY, chrM |
| human | ensembl | 1–22, X, Y, MT |
| mouse | ucsc | chr1–19, chrX, chrY, chrM |
| mouse | ensembl | 1–19, X, Y, MT |

---

## Rerunning specific steps

Steps use sentinel files and output-existence checks to skip completed work.

```bash
# Rerun FastQ Screen only (Step 2b)
rm -f fastQScreen/*

# Rerun DE (Step 16) and enrichment (Step 21) only
rm -f analysis/DE/*_DE_results.tsv
rm -f analysis/enrichment/.enrichment_done

# Rerun enrichment only (Step 21)
rm -f analysis/enrichment/.enrichment_done

# Rerun with custom thresholds (LFC=0 for broader sensitivity)
DE_LFC_THRESHOLD=0 DE_PADJ_THRESHOLD=0.05 \
PADJ_THRESHOLD=0.05 LFC_THRESHOLD=0 \
./scripts/rnaseq2tracks.sh config/config.conf

# Force rerun of all steps
FORCE_RERUN=1 ./scripts/rnaseq2tracks.sh config/config.conf
```

## Cleanup — removing large intermediates

After a successful run, free disk space by removing regenerable files:

```bash
# Dry-run first — see what would be deleted
bash scripts/cleanup_existing_run.sh /path/to/output/ --dry-run

# Live cleanup
bash scripts/cleanup_existing_run.sh /path/to/output/

# Or enable automatic cleanup in config for future runs:
# CLEANUP_INTERMEDIATES=1 in config/config.conf
```

See [`docs/CLEANUP.md`](CLEANUP.md) for the full reference.
