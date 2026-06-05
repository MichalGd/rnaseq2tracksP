# Changelog

## [5.1] — 2026-06-05

### New feature: post-run storage cleanup (Step 22)

Adds automatic removal of large intermediate files after a successful pipeline run.
Cleanup is disabled by default and gated on three pipeline-completion sentinels.

### New scripts
| Script | Origin | Description |
|--------|--------|-------------|
| `scripts/cleanup_existing_run.sh` | NEW v5.1 | One-shot cleanup for already-completed runs. Checks completion sentinels, prints space freed per category. Options: `--dry-run`, `--allchr` |

### Modified scripts
| Script | Change |
|--------|--------|
| `scripts/rnaseq2tracks.sh` | Step 22 added: `cleanup_intermediates()` function with sentinel-gated deletion of trimmed FASTQs, unsorted BAMs, sorted BAMs+indices, raw bedGraphs, uncompressed merged bedGraphs |

### New config parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `CLEANUP_INTERMEDIATES` | `0` | `1` = enable automatic cleanup after successful run |
| `CLEANUP_DRYRUN` | `0` | `1` = dry-run mode (print only, no deletion) |
| `CLEANUP_ALLCHR_BEDGRAPH` | `0` | `1` = also remove `bigwig/*.all_chromosomes.bedGraph.gz` |

### New documentation
- `docs/CLEANUP.md` — full cleanup reference: what is removed, what is kept, sentinel logic, regeneration instructions

### Files removed by cleanup
| Directory | Pattern | Step produced |
|---|---|---|
| `trimmedFastq/` | `*_trimmed.fq.gz`, `*_val_[12].fq.gz` | Step 4 |
| `STARalignments/` | `*_Aligned.out.bam`, `*_SJ.out.tab` | Step 7 |
| `bams/` | `*_sortedS.bam`, `*_sortedS.bam.bai` | Step 8 |
| `bedGraph/raw/` | `*.bedGraph.gz` (un-normalised) | Step 10 |
| `bedGraph/merged/` | `*.bedGraph` (uncompressed only) | Step 14 |
| `bigwig/` *(optional)* | `*.all_chromosomes.bedGraph.gz` | Steps 13/15 |

---

## [5.0] — 2026-05-18

### New scripts
| Script | Origin | Description |
|--------|--------|-------------|
| `scripts/Rscripts/deseq2_enrichment.R` | NEW v5.0 | Gene set enrichment analysis (ORA + GSEA) per contrast: GO BP/MF/CC, KEGG, Reactome, MSigDB Hallmarks |
| `scripts/rerun_deseq_rep12.sh` | NEW v5.0 | Utility script to re-run DESeq2 DE on a replicate subset (Steps 1–3 equivalent) |

### Modified scripts
| Script | Change |
|--------|--------|
| `scripts/rnaseq2tracks.sh` | Step 21 added: gene enrichment analysis (ORA + GSEA) via `deseq2_enrichment.R`; version bumped to v4.3 |
| `scripts/Rscripts/deseq2_de.R` | v4.2: unshrunken (`res_raw`) and shrunken (`res_shrunk`) results produced separately; six volcano PDFs (_raw, _shrunk, _raw_clipped, _shrunk_clipped) and two MA PDFs per contrast; annotated count tables with GTF-derived gene annotation |

### New config parameters
| Parameter | Description |
|-----------|-------------|
| `PADJ_THRESHOLD` | adjusted p-value threshold for ORA gene list in enrichment (default 0.05) |
| `LFC_THRESHOLD` | \|log2FC\| threshold for ORA gene list in enrichment (default 1) |
| `ENRICHMENT_MINGS` | minimum gene set size (default 10) |
| `ENRICHMENT_MAXGS` | maximum gene set size (default 500) |

### New output
- `analysis/enrichment/<contrast_id>/` — per-contrast ORA and GSEA TSV tables and PDF/PNG plots
- `analysis/enrichment/.enrichment_done` — sentinel file for step caching
- `analysis/DE/<contrast_id>_volcano_raw.pdf` — unshrunken LFC volcano plot
- `analysis/DE/<contrast_id>_volcano_shrunk.pdf` — shrunken LFC volcano plot
- `analysis/DE/<contrast_id>_volcano_*_clipped.pdf` — axis-clipped versions
- `analysis/DE/<contrast_id>_MA_raw.pdf` and `_MA_shrunk.pdf` — MA plots
- `analysis/tables/<contrast_id>_annotated_counts.csv` — merged annotation + counts + DE per contrast

### Updated dependencies (`environment.yml`)
Added:
- `bioconductor-clusterprofiler>=4.10`
- `bioconductor-enrichplot>=1.22`
- `bioconductor-reactomepa>=1.46`
- `bioconductor-fgsea>=1.28`
- `bioconductor-org.mm.eg.db>=3.18`
- `bioconductor-org.hs.eg.db>=3.18`
- `r-ashr>=2.2`
- `r-msigdbr>=7.5`

---

## [4.0] — 2026-05-14

### New scripts
| Script | Origin | Description |
|--------|--------|-------------|
| `scripts/preflight_check.sh` | NEW | Validates tools, R packages, RSeQC, genome files before run |
| `scripts/run_rnaseq_qc.sh` | ADAPTED | infer_experiment, read_distribution, junction_annotation, junction_saturation, geneBody_coverage |
| `scripts/collect_star_qc.sh` | NEW | Parses STAR Log.final.out → TSV summary; symlinks for MultiQC |
| `scripts/check_strand_consistency.sh` | NEW | Verifies Fwd+Rev counts ≈ Total mapped; fails on divergence |

### Modified scripts
| Script | Change |
|--------|--------|
| `scripts/rnaseq2tracks.sh` | Step 0 preflight; Steps 9b, 10b, 10c added; Step 19 MultiQC includes 07_qc |
| `config/config_template.conf` | RSEQC_BED_HUMAN/MOUSE, RSEQC_BIN_DIR, RUN_RSEQC, STRAND_TOLERANCE_PCT |

### New docs
- `docs/RSEQC.md` — comprehensive module documentation with annotation file guide and metric interpretation

### New output directories
- `07_qc/star/` — STAR alignment summary and Log.final.out symlinks
- `07_qc/rseqc/` — all RSeQC per-sample outputs and geneBody merged BAM
- `07_qc/multiqc/` — dedicated RSeQC MultiQC report

---

## [3.0] — 2026-05-13

- Multi-species config (human and mouse in one file)
- Configurable chromosome filter (UCSC/Ensembl, human/mouse)
- LICENSE and CITATION.cff added
- Executable smoke test

---

## [2.0] — 2026-05-13

- `bam_to_bedgraph.R` (Rsamtools); `normalize_bedgraph.R` (rtracklayer)
- UCSC chromosome filter
- `docs/OUTPUTS.md`

---

## [1.0] — 2026-05-12

- Initial release
