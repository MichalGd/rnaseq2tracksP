# Storage Cleanup — Removing Large Intermediate Files

**Added in v5.1**

After a successful pipeline run, large intermediate files that can be fully regenerated from the original FASTQ files and pipeline scripts can be safely removed. Two mechanisms are provided: automatic cleanup (built into `rnaseq2tracks.sh` Step 22) and a standalone one-shot script for already-completed runs.

---

## What is removed vs. kept

### Files removed (all safely regenerable)

| Directory | Pattern | Created at step | Typical size |
|---|---|---|---|
| `trimmedFastq/` | `*_trimmed.fq.gz`, `*_val_1.fq.gz`, `*_val_2.fq.gz` | Step 4 (TrimGalore) | ~80% of raw FASTQ |
| `STARalignments/` | `*_Aligned.out.bam` | Step 7 (STAR) | 5–10× raw FASTQ |
| `STARalignments/` | `*_SJ.out.tab` | Step 7 (STAR) | Small–medium |
| `bams/` | `*_sortedS.bam`, `*_sortedS.bam.bai` | Step 8 (samtools) | ~same as unsorted BAM |
| `bedGraph/raw/` | `*.bedGraph.gz` (un-normalised per-sample coverage) | Step 10 (bam_to_bedgraph.R) | 200 MB – 1 GB per sample |
| `bedGraph/merged/` | `*.bedGraph` (uncompressed merged, no `.gz`) | Step 14 (merge_bedgraph_replicates.R) | 1–3 GB per condition |
| `bigwig/` *(optional)* | `*.all_chromosomes.bedGraph.gz` | Step 13/15 (norm_bedgraph_to_bigwig.sh) | 200 MB – 2 GB |

### Files kept

| Directory | Contents |
|---|---|
| Input path (user-defined) | Original raw FASTQ files — never touched by pipeline |
| `bedGraph/normalized/` | `*_norm.bedGraph.gz` — size-factor-normalised per-sample coverage |
| `bigwig/` | `*.sorted.bedGraph.gz`, `*.bw` — final genome browser tracks |
| `STARgeneCounts/` | `*_ReadsPerGene.out.tab` — STAR gene counts (input for DESeq2) |
| `STARlogs/` | `*_Log.final.out` — alignment statistics |
| `fastQC/`, `fastQScreen/` | QC HTML/zip reports |
| `multiQC/` | MultiQC HTML reports |
| `07_qc/` | RSeQC outputs |
| `analysis/counts/` | `dds.RData`, `raw_counts.tsv`, `normalized_counts.tsv`, `size_factors.tsv` |
| `analysis/DE/` | DE result tables, volcano and MA plots |
| `analysis/tables/` | `*_annotated_counts.csv` — merged annotation + counts + DE |
| `analysis/enrichment/` | ORA and GSEA results (TSV + PDF/PNG) |
| `analysis/figures/` | PCA, clustering, heatmap PDFs |
| `reports/` | `pipeline_report.html`, UCSC track definitions |

> **Note:** `bedGraph/raw/*.bedGraph.gz` files are already gzip-compressed when written by `bam_to_bedgraph.R`. The cleanup targets the entire `bedGraph/raw/` directory content. The `bedGraph/merged/` cleanup targets only uncompressed `.bedGraph` files (plain text), leaving any `.gz` files intact.

---

## Option 1 — Automatic cleanup on future runs (Step 22)

The cleanup function is built into `rnaseq2tracks.sh` as **Step 22**. It is disabled by default and runs only after all pipeline completion sentinels are confirmed.

### Enable in `config/config.conf`

```bash
# Remove large intermediate files after successful run (default: 0 = off)
CLEANUP_INTERMEDIATES=1

# Dry-run: print what would be deleted without actually deleting (default: 0)
CLEANUP_DRYRUN=0

# Also remove bigwig/*.all_chromosomes.bedGraph.gz (default: 0 = keep)
CLEANUP_ALLCHR_BEDGRAPH=0
```

### Completion sentinels checked before any deletion

Step 22 aborts if any of these files are missing:

| Sentinel file | Created at step |
|---|---|
| `multiQC/final/multiQC_final.html` | Step 19 |
| `reports/pipeline_report.html` | Step 20 |
| `analysis/enrichment/.enrichment_done` | Step 21 |

At least one `bigwig/*.bw` file must also exist. If any sentinel is missing, cleanup is skipped with a warning and the pipeline exits normally.

### Recommended workflow for first use

```bash
# 1. Enable cleanup in dry-run mode
CLEANUP_INTERMEDIATES=1
CLEANUP_DRYRUN=1

# 2. Run the pipeline (all done_check steps will skip; Step 22 will print what would be deleted)
./scripts/rnaseq2tracks.sh config/config.conf

# 3. Review the Step 22 output. If it looks correct, disable dry-run:
CLEANUP_DRYRUN=0
```

---

## Option 2 — One-shot cleanup for an existing completed run

Use `scripts/cleanup_existing_run.sh` to clean up a run that has already finished.

```bash
# Dry-run first — see every file and total space freed
bash scripts/cleanup_existing_run.sh /path/to/output/ --dry-run

# Live run
bash scripts/cleanup_existing_run.sh /path/to/output/

# Also remove all_chromosomes.bedGraph.gz files
bash scripts/cleanup_existing_run.sh /path/to/output/ --allchr
```

The script performs the same sentinel checks as Step 22 before deleting anything.

---

## Regenerating deleted files

All deleted files can be reconstructed from the original FASTQ files and the pipeline scripts:

| Deleted file type | How to regenerate |
|---|---|
| Trimmed FASTQs | Re-run Steps 4–6 (`trimgalore_single.sh`) |
| Unsorted BAMs | Re-run Step 7 (`star_SE/PE_single.sh`) from trimmed FASTQs |
| Sorted BAMs + indices | Re-run Step 8 (`bam_sort_index.sh`) from unsorted BAMs |
| Raw bedGraphs | Re-run Step 10 (`bam_to_bedgraph.R`) from sorted BAMs |
| Merged bedGraphs | Re-run Steps 12–14 (`normalize_bedgraph.R`, `merge_bedgraph_replicates.R`) |

Set `FORCE_RERUN=1` in config to force re-execution of already-completed steps.

---

## Scripts added

| Script | Location | Purpose |
|---|---|---|
| `cleanup_existing_run.sh` | `scripts/` | One-shot cleanup for completed runs |
| Cleanup function (Step 22) | `scripts/rnaseq2tracks.sh` | Automatic post-run cleanup |

---

## Config parameters added

| Parameter | Default | Description |
|---|---|---|
| `CLEANUP_INTERMEDIATES` | `0` | `1` = enable automatic cleanup after successful run |
| `CLEANUP_DRYRUN` | `0` | `1` = print files that would be deleted without deleting |
| `CLEANUP_ALLCHR_BEDGRAPH` | `0` | `1` = also remove `bigwig/*.all_chromosomes.bedGraph.gz` |
