# Output directory structure

```
<OUTDIR>/
├── fastQC/
│   ├── raw/                    FastQC HTML + zip per FASTQ file
│   └── trimmed/                FastQC HTML + zip per trimmed FASTQ
│
├── multiQC/
│   ├── raw/                    multiQC_raw.html
│   ├── trimmed/                multiQC_trimmed.html
│   ├── alignments/             multiQC_alignments.html
│   └── final/                  multiQC_final.html (all stages combined)
│
├── trimmedFastq/               TrimGalore output .fq.gz files
│
├── STARalignments/             Intermediate unsorted BAMs (can be deleted after step 8)
├── STARlogs/                   *_Log.final.out per sample
├── STARgeneCounts/             *_ReadsPerGene.out.tab per sample
│
├── bams/
│   └── <sample_id>_sortedS.bam      sorted, indexed BAM
│   └── <sample_id>_sortedS.bam.bai
│
├── bedGraph/
│   ├── raw/
│   │   ├── <sample>_Fwd.bedGraph.gz       raw forward coverage
│   │   ├── <sample>_Rev.bedGraph.gz       raw reverse coverage
│   │   └── <sample>_unstranded.bedGraph.gz  (if strandedness=unstranded)
│   ├── normalized/
│   │   ├── <sample>_Fwd_norm.bedGraph.gz  SF_rpm scaled
│   │   └── <sample>_Rev_norm.bedGraph.gz
│   └── merged/
│       ├── <condition>_Fwd_merged.bedGraph  replicate-averaged
│       └── <condition>_Rev_merged.bedGraph
│
├── bigwig/
│   ├── <sample>_FwdS.bw               per-sample forward  ┐ always produced
│   ├── <sample>_RevS.bw               per-sample reverse  ┘
│   ├── <condition>_Fwd_mergedS.bw     merged forward      ┐ if MERGE_REPLICATES=true
│   └── <condition>_Rev_mergedS.bw     merged reverse      ┘
│
├── analysis/
│   ├── counts/
│   │   ├── raw_counts.tsv             gene_id + annotation + raw integer counts
│   │   ├── normalized_counts.tsv      DESeq2 SF-normalized counts
│   │   ├── fpkm_counts.tsv            FPKM per gene
│   │   ├── tpm_counts.tsv             TPM per gene
│   │   ├── size_factors.tsv           sample_id · condition · SF · SF_rpm · exonic_reads
│   │   ├── library_stats.tsv          exonic reads, N_unmapped, N_multimapping, etc.
│   │   ├── dds.RData                  DESeqDataSet (dds + geneInfo + ss)
│   │   └── sessionInfo_normalize.txt
│   ├── DE/
│   │   ├── <contrast_id>_DE_results.tsv  gene_id · gene_name · log2FC · LFC_shrunken · padj
│   │   ├── <contrast_id>_MA_plot.pdf
│   │   ├── <contrast_id>_volcano.pdf
│   │   └── sessionInfo_DE.txt
│   └── figures/
│       ├── PCA.pdf
│       ├── sample_clustering.pdf
│       ├── top50_heatmap.pdf
│       └── sessionInfo_qc.txt
│
└── reports/
    ├── pipeline_report.html      self-contained HTML (kableExtra + ggplot2)
    ├── ucsc_tracks.txt           paste into UCSC Genome Browser > My Data > Custom Tracks
    └── bigwig_summary.txt        track_num · filename · sample_name · file_size_MB
```

## Count table columns

All count tables share a common left block of annotation columns:

| Column | Description |
|---|---|
| `gene_id` | Ensembl gene ID |
| `gene_name` | Gene symbol |
| `gene_type` | Biotype (protein_coding, lncRNA, …) |
| `seqnames` | Chromosome |
| `start` / `end` | Gene coordinates (GTF-based) |
| `medianTxLen` | Median transcript length (used for FPKM/TPM) |
| `<sample_id>` × N | Count value per sample |

## size_factors.tsv columns

| Column | Description |
|---|---|
| `SF` | DESeq2 geometric mean size factor |
| `SF_rpm` | SF anchored to mean RPM across all samples (used for bedGraph scaling) |
| `exonic_reads` | Total reads mapping to exonic features |
