# Workflow — step reference (v5)

## Step table

| Step | Script | Language | Origin | Notes |
|------|--------|----------|--------|-------|
| 0 | `preflight_check.sh` | Bash | NEW v4 | Tools, R pkgs, RSeQC, genome files |
| 1 | mkdir | Bash | — | Output tree including 07_qc/ and analysis/enrichment/ |
| 2 | FastQC | Bash | — | Parallel raw reads |
| 3 | MultiQC | Bash | — | Raw QC |
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
| 19 | MultiQC final | Bash | — | All sources including RSeQC |
| 20 | `pipeline_report.Rmd` | R Markdown | NEW v1 | HTML report |
| 21 | `deseq2_enrichment.R` | R | **NEW v5** | ORA + GSEA: GO BP/MF/CC, KEGG, Reactome, MSigDB Hallmarks |

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

## Strandedness → STAR column

| Value | STAR col | Library |
|-------|---------|---------|
| `unstranded` | 2 | Non-stranded |
| `forward` | 3 | Read 1 on RNA strand |
| `reverse` | 4 | dUTP / NEBNext Ultra II / TruSeq Stranded |

## Chromosome filter

| SPECIES | CHROMOSOME_NAMING | Retained |
|---------|-------------------|---------|
| human | ucsc | chr1–22, chrX, chrY, chrM |
| human | ensembl | 1–22, X, Y, MT |
| mouse | ucsc | chr1–19, chrX, chrY, chrM |
| mouse | ensembl | 1–19, X, Y, MT |

## Rerunning specific steps

Steps use sentinel files and output-existence checks to skip completed work. To rerun only DE and enrichment:

```bash
# Remove sentinels for Steps 16 and 21
rm -f analysis/DE/*_DE_results.tsv
rm -f analysis/enrichment/.enrichment_done

# Run with custom thresholds
DE_LFC_THRESHOLD=0 DE_PADJ_THRESHOLD=0.05 \
PADJ_THRESHOLD=0.05 LFC_THRESHOLD=0 \
./scripts/rnaseq2tracks.sh config/config.conf
```

To force rerun of all steps:
```bash
FORCE_RERUN=1 ./scripts/rnaseq2tracks.sh config/config.conf
```
