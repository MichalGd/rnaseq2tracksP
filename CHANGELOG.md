# Changelog

## [2.0] — 2026-05-13

### Changes from v1.0

#### BigWig (fix)
- `norm_bedgraph_to_bigwig.sh`: non-standard chromosomes and scaffolds are now
  filtered out before `bedGraphToBigWig`. Only `chr1–chr22, chrX, chrY, chrM`
  are retained. BigWigs are clean for UCSC Genome Browser.

#### README
- Mermaid flowchart diagram showing all 20 pipeline steps with colour coding
  (R steps highlighted; Bash steps plain)
- Badges: version, language, aligner, DE tool, UCSC compatibility, layout, license
- Outputs directory tree
- Language map table
- Strandedness guide
- Configuration key variables section
- Documentation index table

## [1.0] — 2026-05-12

### Initial release
- Bash orchestration; R for all analytical steps
- `bam_to_bedgraph.R`: Rsamtools + GenomicAlignments + rtracklayer
  (replaced Bash strand-split approach)
- DESeq2 SF_rpm normalization; apeglm LFC shrinkage
- SE and PE support; per-sample strandedness
- Replicate merging via GRanges disjoin
- R Markdown HTML pipeline report
