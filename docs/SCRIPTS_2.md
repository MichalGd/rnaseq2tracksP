# Script reference (v5)

## Script origin table

| Script | Language | Origin | Change in v5 |
|--------|----------|--------|--------------|
| `rnaseq2tracks.sh` | Bash | NEW v1 | Step 21 (enrichment) added; version v4.3 |
| `preflight_check.sh` | Bash | NEW v4 | — |
| `run_rnaseq_qc.sh` | Bash | ADAPTED v4 | — |
| `collect_star_qc.sh` | Bash | NEW v4 | — |
| `check_strand_consistency.sh` | Bash | NEW v4 | — |
| `rerun_deseq_rep12.sh` | Bash | **NEW v5** | Utility: re-run DE on replicate subset |
| `star_SE_single.sh` | Bash | ADAPTED v1 | — |
| `star_PE_single.sh` | Bash | NEW v1 | — |
| `trimgalore_single.sh` | Bash | ADAPTED v1 | — |
| `bam_sort_index.sh` | Bash | NEW v1 | — |
| `norm_bedgraph_to_bigwig.sh` | Bash | ADAPTED v1→v3 | — |
| `create_ucsc_tracks.sh` | Bash | ADAPTED v1 | — |
| `Rscripts/bam_to_bedgraph.R` | R | NEW v2 | — |
| `Rscripts/deseq2_normalize.R` | R | ADAPTED v1 | — |
| `Rscripts/normalize_bedgraph.R` | R | ADAPTED v2 | — |
| `Rscripts/deseq2_de.R` | R | ADAPTED v4.2 | Unshrunken + shrunken results; 6 volcano + 2 MA plots; annotated count tables |
| `Rscripts/deseq2_enrichment.R` | R | **NEW v5** | ORA + GSEA per contrast: GO, KEGG, Reactome, MSigDB Hallmarks |
| `Rscripts/deseq2_qc_plots.R` | R | ADAPTED v1 | — |
| `Rscripts/merge_bedgraph_replicates.R` | R | ADAPTED v1 | — |
| `Rscripts/pipeline_report.Rmd` | R Markdown | NEW v1 | — |
| `tests/check_bash_syntax.sh` | Bash | VERBATIM | — |
| `tests/run_smoke_test.sh` | Bash | NEW v3→v5 | Added enrichment R package checks (section 3), contrasts file check (section 9) |

---

## deseq2_de.R — v4.2 changes

Both unshrunken (`res_raw`) and shrunken (`res_shrunk`) DESeq2 results are now produced per contrast:

- **`res_raw`** — unshrunken Wald test results; used for volcano and MA plot raw versions and for correct p-value distribution visualization
- **`res_shrunk`** — LFC-shrunken results via apeglm (ashr fallback); canonical reported values in `_DE_results.tsv`

Per contrast, six PDFs and PNGs are produced:
- `_volcano_raw` and `_volcano_shrunk` — full-range volcano plots
- `_volcano_raw_clipped` and `_volcano_shrunk_clipped` — axis-clipped versions (xlim ±4, ylim 0–25)
- `_MA_raw` and `_MA_shrunk` — MA plots

Annotated count tables (`analysis/tables/<contrast_id>_annotated_counts.csv`) merge GTF-derived gene annotation (gene_id, gene_name, gene_type, chr, strand, protein_coding flag), raw counts, normalized counts, and DE statistics for all expressed genes. These are suitable for downstream filtering and visualization.

---

## deseq2_enrichment.R — new in v5

**Usage:**
```bash
Rscript scripts/Rscripts/deseq2_enrichment.R \
  --dedir   <analysis/DE> \
  --contrasts <config/contrasts.csv> \
  --outdir  <analysis/enrichment> \
  --species human|mouse \
  --padj    0.05 \
  --lfc     1 \
  --minGS   10 \
  --maxGS   500
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--dedir` | — | Directory containing `*_DE_results.tsv` files |
| `--contrasts` | — | Contrasts CSV file |
| `--outdir` | — | Output directory (created if absent) |
| `--species` | `mouse` | `human` or `mouse` |
| `--padj` | `0.05` | Adjusted p-value threshold for ORA gene selection |
| `--lfc` | `0` | \|log2FC\| threshold for ORA gene selection |
| `--minGS` | `10` | Minimum gene set size |
| `--maxGS` | `500` | Maximum gene set size |

**GSEA ranking metric:**
The script uses `sign(log2FoldChange) × −log10(padj + 1e−300)` as the ranking metric. This approach is deliberately robust to ashr-shrunken LFC values, which tend to cluster near zero and can cause tied rankings if used directly. The sign-times-significance metric preserves directionality while using the more reliable p-value for magnitude.

**ORA vs GSEA:**
- ORA tests whether a specific subset of genes (significant DE genes) is over-represented in a gene set relative to a background. It is sensitive to the significance threshold used to define the gene list.
- GSEA uses the full ranked gene list and tests whether a gene set is consistently enriched at the top or bottom of the ranking. It does not require an arbitrary significance cutoff and is generally more sensitive when the number of significant genes is small.

---

## run_rnaseq_qc.sh — adaptation notes

**Source scripts from `MichalGd/3end-RNAseq-0.1`:**
- `RNA_RSeQC_QuantSeqRev_17jan2024.sh`: hardcoded paths; `for` loop over `ls *.bam`; runs read_distribution + read_duplication only
- `RSeQC_check17jan24.sh`: single-BAM wrapper for read_distribution.py and read_duplication.py

**Changes in v4 `run_rnaseq_qc.sh`:**
1. All paths parameterised — zero hardcoded paths
2. Samplesheet-driven sample list
3. SE/PE aware
4. Replaced `read_duplication.py` with `junction_annotation.py` + `junction_saturation.py` (more informative for standard RNA-seq)
5. `geneBody_coverage.py` runs on merged BAM (original ran per-sample, which is slow for large cohorts)
6. PID-array job throttling
7. Graceful SKIP per module when binary not found
