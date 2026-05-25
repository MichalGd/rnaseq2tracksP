# Gene Set Enrichment Analysis
This document describes the gene set enrichment analysis (GSEA) implemented in Step 21 of the `rnaseq2tracksP` pipeline, executed by [`scripts/Rscripts/deseq2_enrichment.R`](../scripts/Rscripts/deseq2_enrichment.R).

***
## Purpose
After differential expression analysis (Step 16), each contrast produces a ranked list of genes with associated fold-changes and adjusted p-values. Step 21 translates these gene-level statistics into biological pathway and process-level insights by asking two complementary questions:

- **ORA:** Are significantly DE genes over-represented in any known gene set?
- **GSEA:** Is a gene set coordinately shifted across the entire expression ranking, even if individual genes do not reach significance?

Both methods are run for every contrast defined in `contrasts.csv`. Results are written to `analysis/enrichment/trast_id>/`.

***
## Tools and Databases Used
### R Packages
| Package | Version used | Reference |
|---------|-------------|-----------|
| [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) | Bioconductor | ORA and GSEA for GO and KEGG |
| [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) | Bioconductor | ORA and GSEA for Reactome pathways |
| [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) | Bioconductor | Fast pre-ranked GSEA for MSigDB Hallmarks |
| [msigdbr](https://igordot.github.io/msigdbr/) | CRAN | MSigDB gene sets for multiple organisms |
| [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) | Bioconductor | Visualization (dotplots, cnetplots) |
### Databases
| Database | URL | Used for |
|----------|-----|---------|
| Gene Ontology (GO) | [geneontology.org](https://geneontology.org) | BP, MF, CC ORA and GSEA |
| KEGG Pathways | [kegg.jp](https://www.kegg.jp) | Pathway ORA and GSEA |
| Reactome | [reactome.org](https://reactome.org) | Pathway ORA and GSEA |
| MSigDB Hallmarks | [gsea-msigdb.org/gsea/msigdb](https://www.gsea-msigdb.org/gsea/msigdb) | Hallmark gene set GSEA (fgsea) |

***
## Methods Overview
### 1. ORA — Over-Representation Analysis
**Implemented with:** [`clusterProfiler::enrichGO`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [`enrichKEGG`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [`ReactomePA::enrichPathway`](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html)

**What it does:** Tests whether the list of significantly DE genes contains more members of a gene set than expected by chance, using a hypergeometric test.

**Input:**
- *Significant genes:* genes passing `PADJ_THRESHOLD` (default 0.05) **and** `|log2FC| ≥ LFC_THRESHOLD` (default 0) from the contrast DE results
- *Background:* all genes with a non-NA adjusted p-value in that contrast (i.e. all tested genes — not the whole genome)

**Databases run:** GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG, Reactome

**Key output files per contrast:**

| File | Contents |
|------|----------|
| `*_ORA_GOBP.tsv` | GO Biological Process enrichment table |
| `*_ORA_GOMF.tsv` | GO Molecular Function enrichment table |
| `*_ORA_GOCC.tsv` | GO Cellular Component enrichment table |
| `*_ORA_KEGG.tsv` | KEGG pathway enrichment table |
| `*_ORA_Reactome.tsv` | Reactome pathway enrichment table |
| `*_ORA_*_dotplot.{pdf,png}` | Dot plot (size = gene count, colour = adj. p-value) |
| `*_ORA_*_barplot.{pdf,png}` | Bar plot ranked by gene count |
| `*_ORA_*_cnetplot.{pdf,png}` | Network plot linking genes to enriched terms |

**Interpretation:** Each TSV row is one enriched gene set. The most important columns are `Description` (gene set name), `GeneRatio` (fraction of your sig genes in that set), `p.adjust` (BH-corrected p-value), and `geneID` (the overlapping gene symbols). Rows are sorted by ascending p-value. Focus on sets with `p.adjust < 0.05` and a `GeneRatio` that is biologically plausible (not a single-gene set hitting a very broad term).

***
### 2. GSEA — Gene Set Enrichment Analysis
**Implemented with:** [`clusterProfiler::gseGO`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [`gseKEGG`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [`ReactomePA::gsePathway`](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html)

**What it does:** Ranks *all* expressed genes by a scoring metric and tests whether gene set members cluster at the top or bottom of the ranking. Unlike ORA, it does not require a binary significant/non-significant cutoff and is more sensitive to coordinated changes of modest magnitude.

**Input:**
- *Ranked list:* all genes with valid `log2FoldChange` and `padj`, sorted by:

```
rank_metric = sign(log2FoldChange) × −log10(padj + 1e-300)
```

Genes most strongly upregulated with high confidence score highest; genes most strongly downregulated with high confidence score lowest. This metric is robust to apeglm/ashr LFC shrinkage and avoids tied pile-ups at padj = 1.

- *Gene ID:* Ensembl IDs are converted to Entrez IDs via `AnnotationDbi::bitr()`. Genes without a valid Entrez mapping are excluded.

**Databases run:** [GO](https://geneontology.org) Biological Process, GO Molecular Function, [KEGG](https://www.kegg.jp), [Reactome](https://reactome.org), [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) Hallmarks (via fgsea)

**Key output files per contrast:**

| File | Contents |
|------|----------|
| `*_GSEA_GOBP.tsv` | GO BP GSEA results table |
| `*_GSEA_GOMF.tsv` | GO MF GSEA results table |
| `*_GSEA_KEGG.tsv` | KEGG GSEA results table |
| `*_GSEA_Reactome.tsv` | Reactome GSEA results table |
| `*_GSEAHallmarks.tsv` | MSigDB Hallmarks results (fgsea) |
| `*_GSEA_*_dotplot.{pdf,png}` | Dot plot split by enrichment direction |
| `*_GSEA_*_barplot.{pdf,png}` | NES bar plot (red = up, blue = down) |
| `*_GSEAHallmarks_barplot.{pdf,png}` | Hallmarks NES bar plot (padj < 0.05 only) |

**Interpretation:** The key column is `NES` (Normalized Enrichment Score). Positive NES means the gene set is enriched among upregulated genes; negative NES means enrichment among downregulated genes. `p.adjust` is the BH-corrected permutation p-value. Sets with `|NES| > 1.5` and `p.adjust < 0.25` are generally considered noteworthy. Check `leading_edge` / `core_enrichment` to identify the specific genes driving the enrichment.

> **Note on GSEA with few DEGs:** GSEA requires a sufficient number of genes to produce a meaningful ranking. Contrasts with fewer than 100 genes in the ranked list are automatically skipped with a log message. This avoids a known [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) crash on near-empty ranked distributions.

***
### 3. MSigDB Hallmarks (fgsea)
**Implemented with:** [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html) + [`msigdbr`](https://igordot.github.io/msigdbr/)

**What it does:** Runs fast pre-ranked GSEA against the 50 [MSigDB Hallmark gene sets](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H) — curated, non-redundant gene sets representing well-defined biological states and processes. Gene sets are loaded via `msigdbr` with species-appropriate ortholog mapping (no internet connection required at runtime).

**Input:** Same ranked list as GSEA above (Entrez IDs, sorted by rank metric).

**Output:** `*_GSEAHallmarks.tsv` (all results) and `*_GSEAHallmarks_barplot.pdf/png` (only sets with `padj < 0.05`).

**Interpretation:** Hallmarks are particularly useful for an initial overview — 50 broad, well-characterised pathways. The barplot gives an immediate visual summary of which hallmark processes are up or down-regulated in each contrast.

***
## Thresholds and Settings
### Configurable via `config.conf`
| Parameter | Default | Effect |
|-----------|---------|--------|
| `PADJ_THRESHOLD` | `0.05` | BH-adjusted p-value cutoff for ORA input gene list |
| `LFC_THRESHOLD` | `0` | Minimum \|log2FC\| for ORA input (0 = any direction) |
| `ENRICHMENT_MINGS` | `10` | Minimum gene set size for all methods |
| `ENRICHMENT_MAXGS` | `500` | Maximum gene set size for all methods |
### Hardcoded in `deseq2_enrichment.R`
| Method | Parameter | Value | Rationale |
|--------|-----------|-------|-----------|
| ORA (all databases) | `pvalueCutoff` | `0.05` | Standard hypergeometric p-value threshold |
| ORA GO terms | `qvalueCutoff` | `0.2` | Relaxed BH threshold; [GO](https://geneontology.org) terms are correlated, which inflates FDR |
| ORA (all) | `pAdjustMethod` | `"BH"` | Benjamini-Hochberg correction |
| GSEA (all databases) | `pvalueCutoff` | `0.25` | Standard for rank-based GSEA with BH correction; permutation null is already stringent |
| GSEA | `nPermSimple` | `10000` | High permutation count for stable p-values |
| [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) Hallmarks | `nPermSimple` | `10000` | Same as above |
| Hallmarks barplot | padj filter | `< 0.05` | Only significantly enriched Hallmarks are plotted |
| GSEA skip threshold | `length(rankedList)` | `< 100` | Prevents crashes on near-empty contrasts |

> **Why is the GSEA `pvalueCutoff` set to 0.25?**  
> Rank-based GSEA with BH correction is inherently conservative. The permutation-based null already accounts for multiple testing at the gene level. Many biologically meaningful pathway enrichments produce adjusted p-values between 0.05 and 0.25, particularly with small sample sizes (n = 2 replicates). The NES magnitude and leading-edge gene coherence serve as additional interpretive filters beyond the p-value alone.
### Species and Database Settings
The organism database ([`org.Mm.eg.db`](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) for mouse, [`org.Hs.eg.db`](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) for human) is selected automatically based on the `SPECIES` variable in `config.conf`. [KEGG](https://www.kegg.jp) organism codes (`mmu`/`hsa`) and [Reactome](https://reactome.org) organism names (`mouse`/`human`) are set accordingly. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) Hallmarks are loaded with [`msigdbr`](https://igordot.github.io/msigdbr/) using species-specific ortholog mapping.

***
## Recommended Files to Inspect First
For a new contrast, check these files in order:

1. **`*_GSEAHallmarks_barplot.pdf`**  
   The fastest overview. 50 [MSigDB Hallmarks](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H) cover most major biological processes. The NES barplot immediately shows which broad pathways are up or down. Start here before diving into individual GO/KEGG results.

2. **`*_ORA_GOBP.tsv` + `*_ORA_GOBP_dotplot.pdf`**  
   [GO Biological Process](https://geneontology.org) is the richest database for mechanistic insight. The dotplot shows the top enriched terms with significance and overlap size at a glance. Filter to `p.adjust < 0.05` rows.

3. **`*_ORA_Reactome.tsv` + `*_ORA_Reactome_dotplot.pdf`**  
   [Reactome](https://reactome.org) pathways are more curated and less redundant than GO BP. Useful for identifying specific signalling cascades or metabolic routes.

4. **`*_GSEA_KEGG.tsv` + `*_GSEA_KEGG_barplot.pdf`**  
   [KEGG](https://www.kegg.jp) GSEA complements ORA by capturing coordinated but sub-threshold changes in metabolic and signalling pathways. The NES barplot separates up- and down-regulated pathways clearly.

5. **`*_ORA_*_cnetplot.pdf`** (any database)  
   The concept network plot links enriched gene sets to the individual genes driving them. Essential for identifying which specific genes are responsible for an enriched term — particularly useful when a term is unexpected.
