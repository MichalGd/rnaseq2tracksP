# Outputs

## Full directory tree

```
<OUTDIR>/
├── fastQC/
│   ├── raw/                   FastQC HTML + zip per sample (raw)
│   └── trimmed/               FastQC HTML + zip (trimmed)
├── multiQC/
│   ├── raw/multiQC_raw.html
│   ├── trimmed/multiQC_trimmed.html
│   ├── alignments/multiQC_alignments.html
│   └── final/multiQC_final.html        ← includes RSeQC outputs
├── trimmedFastq/              trimmed FASTQ (.fq.gz)
├── STARalignments/            *Aligned.out.bam
├── STARlogs/                  *Log.final.out
├── STARgeneCounts/            *ReadsPerGene.out.tab
├── bams/                      *_sortedS.bam + .bai
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
│   ├── raw/                   *_FwdS.bedGraph.gz  *_RevS.bedGraph.gz
│   ├── normalized/            *_FwdS_norm.bedGraph.gz  *_all_chromosomes.bedGraph.gz
│   └── merged/                <condition>_FwdS_norm_merged.bedGraph
├── bigwig/                    *_FwdS_norm.bw  *_RevS_norm.bw
│                              <condition>_FwdS_norm_merged.bw
├── analysis/
│   ├── counts/
│   │   ├── raw_counts.tsv
│   │   ├── normalized_counts.tsv
│   │   ├── size_factors.tsv   sample_id, size_factor, sf_rpm
│   │   └── dds.RData
│   ├── DE/                    *_DE_results.tsv  *_significant.tsv  *_volcano.pdf
│   └── figures/               PCA.pdf  sample_clustering.pdf  heatmaps
└── reports/
    ├── pipeline_report.html
    └── ucsc_tracks.txt
```

## Key output files

| File | Description |
|------|-------------|
| `star_alignment_summary.tsv` | Uniquely mapped %, multi-mapped %, input reads per sample |
| `size_factors.tsv` | DESeq2 size factor + SF_rpm per sample |
| `raw_counts.tsv` | Raw gene counts matrix |
| `normalized_counts.tsv` | DESeq2 size-factor normalized counts |
| `*_infer_experiment.txt` | Strandedness fraction (validate samplesheet) |
| `*_read_distribution.txt` | % reads in CDS, UTR, intron, intergenic |
| `multiQC_final.html` | Complete QC report |
