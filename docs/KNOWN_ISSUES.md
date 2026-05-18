# Known issues

1. **STAR shared memory**: If the pipeline is interrupted, shared memory may not be released.
   Fix: `STAR --genomeDir $STAR_INDEX --genomeLoad Remove`

2. **apeglm falls back to ashr**: For multi-group designs where contrasts are not simple
   model coefficients, apeglm fails with `coef %in% resultsNamesDDS is not TRUE`.
   The script automatically falls back to `ashr`, which is statistically valid. This
   warning is expected and does not affect result quality.

3. **geneBody_coverage memory**: Merging large BAMs (>30 samples, >30 M reads each)
   may require >64 GB RAM. Reduce `MAX_JOBS` or run geneBody_coverage.py separately.

4. **RSeQC BED file format**: Must be BED12. BED6 will cause RSeQC to report all reads
   as intergenic. Verify with `awk 'NF==12' annotation.bed | wc -l`.

5. **CHROMOSOME_NAMING mismatch**: If BigWig chromosomes do not match chrom.sizes,
   `bedGraphToBigWig` will exit. Verify: `head -1 annotation.gtf` vs `head -1 chrom.sizes`.

6. **TrimGalore PE output naming**: TrimGalore uses `--basename` to produce
   `${SID}_val_1.fq.gz` / `${SID}_val_2.fq.gz`. If R1/R2 filenames are
   non-standard, check `trimgalore_single.sh` `--basename` handling.

7. **junction_annotation.py on 3'-mRNA-seq**: Very few junctions are expected for
   QuantSeq / 3'-end libraries. A high "unknown" fraction and low overall junction
   count are normal for this library type.

8. **Strand check false positives with low coverage**: Very low-coverage samples
   (<1 M mapped reads) may have noisy strand ratios. Increase `STRAND_TOLERANCE_PCT`
   to 15 for pilot experiments.

9. **R.utils not installed**: `bam_to_bedgraph.R` uses R.utils for gzip output.
   Install with: `install.packages("R.utils")` or via conda (`r-r.utils` in `environment.yml`).

10. **MultiQC version**: RSeQC module support requires MultiQC >= 1.14.
    Check: `multiqc --version`

11. **MultiQC colour warnings**: MultiQC >= 1.20 expects hex colour codes. If your
    `multiqc_config.yaml` specifies RGB tuples (e.g., `77,175,74`), you will see
    `Error converting color` warnings. These are cosmetic — reports are generated
    correctly. Fix by converting to hex: `77,175,74` → `#4DAF4A`; `55,126,184` → `#377EB8`.

12. **Enrichment Step 21 — GSEA "no term enriched"**: With few DE genes or small contrasts,
    GSEA may find no enriched terms in some databases. This is expected and not an error.
    ORA results are typically available when GSEA returns nothing. Consider relaxing
    `DE_LFC_THRESHOLD` to 0 (catching all padj < 0.05 genes) to increase GSEA input size.

13. **Enrichment Entrez ID mapping**: Only ~60% of Ensembl gene IDs in the human genome
    map to Entrez IDs via `bitr()`. Unannotated and non-protein-coding genes are commonly
    not mapped. This is expected and does not indicate a problem. The ORA background is
    correctly restricted to genes with valid Entrez mappings.
