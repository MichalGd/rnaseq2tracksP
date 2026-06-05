# Outputs

## Full directory tree

```
<OUTDIR>/
├── fastQC/
│   ├── raw/                   FastQC HTML + zip per sample (raw)
│   └── trimmed/               FastQC HTML + zip (trimmed)
├── fastQScreen/               FastQ Screen per-sample reports (Step 2b)
│   ├── *_screen.txt           Tab-separated mapping % per reference database
│   ├── *_screen.html          Interactive per-sample bar chart
│   └── *_screen.png           Static plot
├── multiQC/
│   ├── raw/multiQC_raw.html
│   ├── trimmed/multiQC_trimmed.html
│   ├── alignments/multiQC_alignments.html
│   └── final/multiQC_final.html        ← includes RSeQC + FastQ Screen outputs  [SENTINEL]
├── trimmedFastq/              trimmed FASTQ (.fq.gz)                              [CLEANUP — Step 22]
├── STARalignments/            *Aligned.out.bam, *_SJ.out.tab                     [CLEANUP — Step 22]
├── STARlogs/                  *Log.final.out
├── STARgeneCounts/            *ReadsPerGene.out.tab                               ← KEPT (DESeq2 input)
├── bams/                      *_sortedS.bam + .bai                               [CLEANUP — Step 22]
├── 07_qc/
│   ├── star/
│   │   ├── *_Log.final.out    (symlinks)
│   │   └── star_alignment_summary.tsv
│   ├── rseqc/
│   │   ├── infer_experiment/  *_infer_experiment.txt
│   │   ├── read_distribution/ *_read_distribution.txt
│   │   ├── junction_annotation/ *.*
│   │   ├── junction_saturation/ *.*
│   │   └── genebody/          all_samples.geneBodyCoverage.*
│   └── multiqc/
│       └── multiQC_rseqc.html
├── bedGraph/
│   ├── raw/                   *_FwdS.bedGraph.gz  *_RevS.bedGraph.gz             [CLEANUP — Step 22]
│   ├── normalized/            *_FwdS_norm.bedGraph.gz  *_RevS_norm.bedGraph.gz   ← KEPT
│   └── merged/                <condition>_FwdS_norm_merged.bedGraph              [CLEANUP — uncompressed only]
├── bigwig/                    *_FwdS_norm.bw  *_RevS_norm.bw                     ← KEPT
│                              <condition>_FwdS_norm_merged.bw                    ← KEPT
│                              *.all_chromosomes.bedGraph.gz                      [CLEANUP — optional, --allchr flag]
├── analysis/
│   ├── counts/
│   │   ├── raw_counts.tsv
│   │   ├── normalized_counts.tsv
│   │   ├── size_factors.tsv        sample_id, size_factor, sf_rpm
│   │   └── dds.RData
│   ├── DE/
│   │   ├── *_DE_results.tsv        all genes, shrunken LFC, sorted by padj
│   │   ├── *_significant.tsv       filtered by PADJ_THRESHOLD and LFC_THRESHOLD
│   │   ├── *_volcano_raw.pdf/.png
│   │   ├── *_volcano_shrunk.pdf/.png
│   │   ├── *_volcano_raw_clipped.pdf/.png
│   │   ├── *_volcano_shrunk_clipped.pdf/.png
│   │   ├── *_MA_raw.pdf/.png
│   │   ├── *_MA_shrunk.pdf/.png
│   │   └── deseq2_de_sessionInfo.txt
│   ├── tables/
│   │   └── *_annotated_counts.csv  gene annotation + raw + norm counts + DE per contrast
│   ├── figures/
│   │   ├── PCA.pdf
│   │   ├── sample_clustering.pdf
│   │   └── heatmaps.pdf
│   └── enrichment/
│       ├── .enrichment_done        sentinel file (delete to rerun Step 21)       [SENTINEL]
│       ├── deseq2_enrichment_sessionInfo.txt
│       └── <contrast_id>/
│           ├── *_ORA_GOBP.tsv      ORA results table — GO Biological Process
│           ├── *_ORA_GOMF.tsv      ORA results table — GO Molecular Function
│           ├── *_ORA_GOCC.tsv      ORA results table — GO Cellular Component
│           ├── *_ORA_KEGG.tsv      ORA results table — KEGG
│           ├── *_ORA_Reactome.tsv  ORA results table — Reactome
│           ├── *_ORA_*_dotplot.pdf/.png
│           ├── *_ORA_*_barplot.pdf/.png
│           ├── *_ORA_*_cnetplot.pdf/.png
│           ├── *_GSEA_GOBP.tsv     GSEA results table — GO Biological Process
│           ├── *_GSEA_GOMF.tsv     GSEA results table — GO Molecular Function
│           ├── *_GSEA_KEGG.tsv     GSEA results table — KEGG
│           ├── *_GSEA_Reactome.tsv GSEA results table — Reactome
│           ├── *_GSEA_Hallmarks.tsv GSEA results table — MSigDB Hallmarks
│           ├── *_GSEA_*_dotplot.pdf/.png
│           ├── *_GSEA_*_barplot.pdf/.png
│           └── *_GSEA_Hallmarks_barplot.pdf/.png
└── reports/
    ├── pipeline_report.html                                                       [SENTINEL]
    └── ucsc_tracks.txt
```

### Storage cleanup (Step 22)

Directories and files marked `[CLEANUP]` above are large intermediate files that are fully regenerable from raw FASTQs and pipeline scripts. They are removed automatically by Step 22 when `CLEANUP_INTERMEDIATES=1` in `config/config.conf`, after verifying that all three `[SENTINEL]` files exist.

Use `cleanup_existing_run.sh` to clean up runs that have already completed:

```bash
bash scripts/cleanup_existing_run.sh /path/to/output/ --dry-run   # preview
bash scripts/cleanup_existing_run.sh /path/to/output/             # live
```

See [`docs/CLEANUP.md`](CLEANUP.md) for full details, kept/removed file lists, and regeneration instructions.

---

## Key output files

| File | Description |
|------|-------------|
| `fastQScreen/*_screen.txt` | Per-sample mapping % against Mouse, Human, Zebrafish, Drosophila, Mycoplasma panels |
| `star_alignment_summary.tsv` | Uniquely mapped %, multi-mapped %, input reads per sample |
| `size_factors.tsv` | DESeq2 size factor and SF_rpm per sample |
| `raw_counts.tsv` | Raw gene counts matrix |
| `normalized_counts.tsv` | DESeq2 size-factor normalized counts |
| `*_DE_results.tsv` | Full DE results with shrunken LFC for all expressed genes |
| `*_significant.tsv` | DE results filtered by PADJ_THRESHOLD and LFC_THRESHOLD |
| `*_annotated_counts.csv` | GTF annotation + raw/norm counts + DE statistics merged per gene |
| `*_ORA_*.tsv` | Over-representation analysis results (gene ID, description, p.adjust, gene ratio) |
| `*_GSEA_*.tsv` | GSEA results (NES, p.adjust, leading edge genes) |
| `*_infer_experiment.txt` | Strandedness fraction — use to validate samplesheet |
| `*_read_distribution.txt` | % reads in CDS, UTR, intron, intergenic |
| `multiQC_final.html` | Complete QC report including RSeQC and FastQ Screen |

---

## FastQ Screen output columns (`*_screen.txt`)

| Column | Description |
|--------|-------------|
| `Genome` | Database name as defined in `fastq_screen.conf` |
| `#Reads_processed` | Total reads in the screened subset |
| `#Unmapped` | Reads with no hit in this database |
| `%Unmapped` | Percentage unmapped |
| `#One_hit_one_library` | Reads mapping uniquely to this database only |
| `%One_hit_one_library` | Percentage uniquely mapped to this database only |
| `#Multiple_hits_one_library` | Reads with multiple hits in this database only |
| `#One_hit_multiple_libraries` | Reads mapping to this and at least one other database |
| `%One_hit_multiple_libraries` | Percentage of cross-mapping reads |

**Interpretation:** for a clean mouse sample, `%One_hit_one_library` for Mouse should be ~90–97 %. Mycoplasma `%One_hit_one_library` > 0.5 % warrants investigation; > 2 % is a confirmed contamination flag.

---

## Enrichment output columns

### ORA tables (`*_ORA_*.tsv`)
| Column | Description |
|--------|-------------|
| `ID` | Gene set identifier |
| `Description` | Gene set name |
| `GeneRatio` | Ratio of sig genes in set vs total sig genes |
| `BgRatio` | Ratio of background genes in set |
| `pvalue` | Fisher's exact test p-value |
| `p.adjust` | BH-adjusted p-value |
| `qvalue` | q-value |
| `geneID` | Gene symbols in the set (readable) |
| `Count` | Number of sig genes in the set |

### GSEA tables (`*_GSEA_*.tsv`)
| Column | Description |
|--------|-------------|
| `ID` | Gene set identifier |
| `Description` | Gene set name |
| `NES` | Normalized Enrichment Score (positive = upregulated, negative = downregulated) |
| `pvalue` | GSEA permutation p-value |
| `p.adjust` | BH-adjusted p-value |
| `core_enrichment` | Leading edge genes |
