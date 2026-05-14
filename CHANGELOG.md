# Changelog

## [4.0] — 2026-05-14

### New scripts
| Script | Origin | Description |
|--------|--------|-------------|
| `scripts/preflight_check.sh` | NEW | Validates tools, R packages, RSeQC, genome files before run |
| `scripts/run_rnaseq_qc.sh` | ADAPTED from RNA_RSeQC_QuantSeqRev_17jan2024.sh + RSeQC_check17jan24.sh | infer_experiment, read_distribution, junction_annotation, junction_saturation, geneBody_coverage |
| `scripts/collect_star_qc.sh` | NEW | Parses STAR Log.final.out → TSV summary; symlinks for MultiQC |
| `scripts/check_strand_consistency.sh` | NEW | Verifies Fwd+Rev counts ≈ Total mapped; fails on divergence |

### Modified scripts
| Script | Change |
|--------|--------|
| `scripts/rnaseq2tracks.sh` | Step 0 preflight; Steps 9b, 10b, 10c added; Step 19 MultiQC includes 07_qc |
| `config/config_template.conf` | RSEQC_BED_HUMAN/MOUSE, RSEQC_BIN_DIR, RUN_RSEQC, STRAND_TOLERANCE_PCT |

### New docs
- `docs/RSEQC.md` — comprehensive module documentation with annotation file guide and metric interpretation

### New output directory
- `07_qc/star/` — STAR alignment summary + Log.final.out symlinks
- `07_qc/rseqc/` — all RSeQC per-sample outputs + geneBody merged BAM
- `07_qc/multiqc/` — dedicated RSeQC MultiQC report

## [3.0] — 2026-05-13
- Multi-species config (human+mouse in one file)
- Configurable chromosome filter (UCSC/Ensembl, human/mouse)
- LICENSE + CITATION.cff
- Executable smoke test

## [2.0] — 2026-05-13
- bam_to_bedgraph.R (Rsamtools); normalize_bedgraph.R (rtracklayer)
- UCSC chromosome filter; docs/OUTPUTS.md

## [1.0] — 2026-05-12
- Initial release
