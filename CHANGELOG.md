# Changelog

## [3.0] — 2026-05-13

### New features
- **Multi-species config** — single `config.conf` holds human and mouse genome paths;
  switch with `SPECIES=human|mouse`; master script resolves STAR_INDEX / GTF / CHROM_SIZES automatically
- **Configurable chromosome filter** — `CHROMOSOME_NAMING=ucsc|ensembl` and
  `REGULAR_CHROMS_ONLY=true|false`; all four human/mouse × UCSC/Ensembl combinations supported
- **Debugging intermediate** — `<stem>.all_chromosomes.bedGraph.gz` retained alongside BigWig
  when `REGULAR_CHROMS_ONLY=true` (removed when false)
- **LICENSE** (MIT) added
- **CITATION.cff** added
- **Executable smoke test** `tests/run_smoke_test.sh` — checks Bash syntax, R packages,
  tool availability, config completeness, and samplesheet parsing before a real run

### Files changed
| File | Change |
|---|---|
| `scripts/norm_bedgraph_to_bigwig.sh` | CHROMOSOME_NAMING + REGULAR_CHROMS_ONLY logic |
| `scripts/rnaseq2tracks.sh` | Species path resolution; env var export to BigWig step |
| `config/config_template.conf` | Human + mouse path pairs; CHROMOSOME_NAMING; REGULAR_CHROMS_ONLY |
| `README.md` | Species badge; chromosome filter section; smoke test in quick start |
| `tests/run_smoke_test.sh` | New executable test (5 checks) |
| `LICENSE` | New — MIT |
| `CITATION.cff` | New |

## [2.0] — 2026-05-13
- `bam_to_bedgraph.R` (Rsamtools + GenomicAlignments) replaces Bash strand-split approach
- `normalize_bedgraph.R` upgraded to rtracklayer import/export
- UCSC chromosome filter added to `norm_bedgraph_to_bigwig.sh`
- `docs/OUTPUTS.md` added
- All R scripts write `sessionInfo.txt`

## [1.0] — 2026-05-12
- Initial release: Bash orchestration + R analytical modules
- DESeq2 SF_rpm normalization; apeglm LFC shrinkage
- SE and PE support; per-sample strandedness
- Replicate merging via GRanges disjoin
- R Markdown HTML pipeline report
