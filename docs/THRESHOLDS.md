# Statistical Thresholds: Differential Expression and Gene Set Enrichment Analysis
This document explains the statistical thresholds used in the `rnaseq2tracksP` pipeline (v4.x) for differential gene expression (DEG) analysis (Step 16) and gene set enrichment analysis (GSEA/ORA, Step 21). It serves as the authoritative reference for configuring `config.conf` and understanding what each parameter controls.

***
## Overview of the Two-Stage Design
The pipeline separates statistical stringency into two conceptually distinct stages:

1. **DEG calling** (Step 16, `deseq2_de.R`): Identifies significantly differentially expressed genes per contrast. Uses strict thresholds appropriate for gene-level claims.
2. **Gene set enrichment** (Step 21, `deseq2_enrichment.R`): Tests whether biological pathways/gene sets are coordinately regulated. Uses deliberately relaxed input thresholds to maximize pathway-level power.

This two-stage design is intentional — pathway analysis benefits from an inclusive gene universe and uses its own internal multiple-testing correction to control false discovery.

***
## Stage 1 — Differential Gene Expression (Step 16)
### Parameters in `config.conf`
```bash
PADJ_THRESHOLD=0.05    # BH-adjusted p-value cutoff for DEG significance
LFC_THRESHOLD=0        # |log2 fold change| minimum (0 = any direction, no magnitude filter)
```

These values are passed to `deseq2_de.R` via `--padj` and `--lfc`.
### What Each Threshold Controls
| Parameter | Default | Controls | Applied to |
|-----------|---------|----------|------------|
| `PADJ_THRESHOLD` | 0.05 | BH-adjusted p-value cutoff | Significance annotation in volcano/MA plots; `*_significant.tsv` output |
| `LFC_THRESHOLD` | 0 | Minimum \|log2FC\| for significance tier | Same as above |
### Significance Tiers
`deseq2_de.R` classifies each gene into one of three tiers used for plot colouring and output filtering:

| Tier | Colour | Criteria |
|------|--------|----------|
| `ns` | grey | padj ≥ 0.10 or NA |
| `trend` | green | 0.05 ≤ padj < 0.10 |
| `sig` | red | padj < 0.05 |

The `LFC_THRESHOLD` additionally gates which genes appear in `*_significant.tsv` (genes with `padj < PADJ_THRESHOLD AND |log2FC| ≥ LFC_THRESHOLD`).
### Rationale for `LFC_THRESHOLD=0`
Setting `LFC_THRESHOLD=0` means **any statistically significant change is reported**, regardless of fold-change magnitude. This is appropriate when:
- Comparing conditions with modest but biologically meaningful regulation (e.g., transcription factor knockdowns, early time points)
- The downstream enrichment step is expected to provide biological context for small but consistent changes
- You do not want to pre-filter genes before pathway analysis

A more stringent alternative (`LFC_THRESHOLD=1`, i.e., ≥2-fold change) is standard for studies where biological effect size is of primary interest. The pipeline defaults to `0` to be inclusive and leave magnitude filtering to the enrichment stage.
### LFC Shrinkage
`deseq2_de.R` applies **apeglm shrinkage** (with ashr as fallback) to log2 fold changes before writing output TSVs. Shrunken LFCs are used for all reported values and plots; unshrunken Wald estimates are used only for p-value computation and volcano plot variants labelled `raw`. This is the current best practice recommended by the DESeq2 developers.

***
## Stage 2 — Gene Set Enrichment Analysis (Step 21)
Step 21 runs two complementary approaches in parallel for each contrast: **Over-Representation Analysis (ORA)** and **GSEA** (rank-based). These differ fundamentally in their input and the thresholds they apply.
### Input Gene Filtering (`config.conf` parameters)
```bash
PADJ_THRESHOLD=0.05    # Reused: filters the significant gene list passed to ORA
LFC_THRESHOLD=0        # Reused: |log2FC| minimum for ORA input genes
ENRICHMENT_MINGS=10    # Minimum gene set size (applied to all methods)
ENRICHMENT_MAXGS=500   # Maximum gene set size (applied to all methods)
```

> **Note:** `PADJ_THRESHOLD` and `LFC_THRESHOLD` are shared between Step 16 and Step 21 in `config.conf`. Step 21 uses them to define its ORA input gene list independently — changing them affects both stages simultaneously.
### ORA (Over-Representation Analysis)
ORA tests whether significantly DE genes are over-represented in a gene set relative to a background of all expressed genes (all genes with non-NA padj in the dataset).

#### ORA Input Gene Selection

```r
sig_genes <- de[!is.na(padj) & padj <= opt$padj & abs(log2FoldChange) >= opt$lfc, ENTREZID]
background <- de[!is.na(padj), ENTREZID]
```

The background is the full set of tested genes — this is the statistically correct approach and avoids the bias introduced by using genome-wide gene lists as background.

#### ORA p-value and q-value Cutoffs (Internal to `deseq2_enrichment.R`)

These cutoffs control which enriched gene sets appear in the output TSVs and plots. They are **not** exposed in `config.conf` and are hardcoded in the script:

| Database | `pvalueCutoff` | `qvalueCutoff` | BH correction |
|----------|---------------|---------------|---------------|
| GO Biological Process | 0.05 | 0.2 | Yes |
| GO Molecular Function | 0.05 | 0.2 | Yes |
| GO Cellular Component | 0.05 | 0.2 | Yes |
| KEGG | 0.05 | — | Yes |
| Reactome | 0.05 | — | Yes |

- `pvalueCutoff`: Raw hypergeometric p-value threshold before BH adjustment
- `qvalueCutoff`: BH-adjusted p-value (q-value) threshold; filters to gene sets that pass FDR correction
- The `qvalueCutoff=0.2` for GO terms is a deliberate relaxation relative to the raw p-value cutoff, acknowledging that GO term redundancy and correlation inflate BH-adjusted values
### GSEA (Gene Set Enrichment Analysis)
GSEA operates on the **entire ranked gene list**, not a binary significant/non-significant split. All genes with valid LFC and padj values contribute. This makes GSEA more powerful than ORA for detecting coordinated changes of modest magnitude.

#### Ranking Metric

```r
rank_metric = sign(log2FoldChange) * -log10(padj + 1e-300)
```

This metric combines directionality (sign of LFC) with statistical confidence (-log10 padj). It is robust to apeglm/ashr shrinkage (which compresses LFC magnitudes without altering sign) and avoids tied pile-ups at padj=1. Genes with the largest magnitude significant changes rank at the extremes; non-significant genes cluster near zero.

#### GSEA p-value Cutoffs (Internal to `deseq2_enrichment.R`)

| Database | `pvalueCutoff` | BH correction | Permutations |
|----------|---------------|---------------|-------------|
| GO Biological Process | 0.25 | Yes | 10,000 |
| GO Molecular Function | 0.25 | Yes | 10,000 |
| KEGG | 0.25 | Yes | 10,000 |
| Reactome | 0.25 | Yes | 10,000 |
| MSigDB Hallmarks (fgsea) | reported at all padj; plot filtered at padj < 0.05 | Yes | 10,000 |

The `pvalueCutoff=0.25` for GSEA is deliberately more permissive than ORA's `0.05`. This is **field-standard practice** for rank-based enrichment because:
- GSEA with BH correction is inherently conservative — the permutation-based null distribution is already stringent
- Many biologically meaningful pathway activations produce NES with adjusted p between 0.05–0.25, particularly with moderate sample sizes (n=2–3 replicates)
- The NES magnitude and leading-edge gene coherence provide additional interpretive filters beyond the p-value alone
### Gene Set Size Filtering
```bash
ENRICHMENT_MINGS=10    # Excludes very small gene sets (too few genes for reliable statistics)
ENRICHMENT_MAXGS=500   # Excludes very large gene sets (too broad to be interpretable)
```

Applied identically to all ORA and GSEA methods. The `10/500` window is the community consensus default in clusterProfiler documentation and widely adopted in the literature.

***
## Full Threshold Reference Table
| Parameter | Location | Stage | Default | Recommended range | Notes |
|-----------|----------|-------|---------|------------------|-------|
| `PADJ_THRESHOLD` | `config.conf` | DEG + ORA input | 0.05 | 0.05–0.1 | Shared between Step 16 and Step 21 |
| `LFC_THRESHOLD` | `config.conf` | DEG + ORA input | 0 | 0–1 | 0 = inclusive; 1 = ≥2-fold magnitude filter |
| `ENRICHMENT_MINGS` | `config.conf` | All GSEA/ORA | 10 | 10–15 | Minimum gene set size |
| `ENRICHMENT_MAXGS` | `config.conf` | All GSEA/ORA | 500 | 300–500 | Maximum gene set size |
| ORA `pvalueCutoff` (GO BP/MF/CC) | `deseq2_enrichment.R` | ORA output filter | 0.05 | 0.05 | Hardcoded; raw hypergeometric p |
| ORA `qvalueCutoff` (GO) | `deseq2_enrichment.R` | ORA output filter | 0.2 | 0.1–0.25 | BH-adjusted; relaxed to account for GO redundancy |
| ORA `pvalueCutoff` (KEGG, Reactome) | `deseq2_enrichment.R` | ORA output filter | 0.05 | 0.05 | Standard |
| GSEA `pvalueCutoff` (all databases) | `deseq2_enrichment.R` | GSEA output filter | 0.25 | 0.2–0.25 | Standard for rank-based GSEA with BH correction |
| Hallmarks plot filter | `deseq2_enrichment.R` | fgsea plot only | padj < 0.05 | 0.05–0.25 | Only affects the barplot, not TSV output |
| GSEA permutations | `deseq2_enrichment.R` | All GSEA | 10,000 | 1,000–10,000 | 10,000 recommended for final results |
