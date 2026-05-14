# RSeQC RNA-seq QC — documentation

## Required BED12 annotation

All RSeQC modules need a BED12 transcript annotation file.

### Download prebuilt files

| Genome | Source |
|--------|--------|
| hg38 (GENCODE v42) | https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/ |
| mm39 | https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/ |

### Generate from GTF

```bash
gtfToGenePred annotation.gtf annotation.genePred
genePredToBed annotation.genePred annotation.bed
```

Set in config:
```bash
RSEQC_BED_HUMAN="/path/to/hg38_GENCODE_V45_Comprehensive.bed"
RSEQC_BED_MOUSE="/path/to/mm39_GENCODE_vM31_Comprehensive.bed"
```

---

## Modules

### infer_experiment.py

**What it does:** Classifies read orientation relative to transcript annotation
to infer library strandedness.

**Output:** `07_qc/rseqc/infer_experiment/<sample>_infer_experiment.txt`

**Interpretation:**
```
Fraction of reads failed to determine: 0.02
Fraction of reads explained by "1++,1--,2+-,2-+": 0.03   ← forward
Fraction of reads explained by "1+-,1-+,2++,2--": 0.95   ← reverse
```
- Dominant fraction > 0.85 → stranded library
- Both ≈ 0.5 → unstranded
- Use to validate your samplesheet `strandedness` column

**MultiQC:** ✅ strandedness inference plot

---

### read_distribution.py

**What it does:** Assigns reads to genomic features (CDS exons, UTRs, introns, intergenic).

**Output:** `07_qc/rseqc/read_distribution/<sample>_read_distribution.txt`

**Interpretation:**
- Good RNA-seq: >60% reads in CDS + UTR exons
- High intronic (>15%): possible genomic DNA contamination or degradation
- For 3' mRNA-seq (QuantSeq REV): dominant 3'UTR fraction is expected

**MultiQC:** ✅ stacked bar chart

---

### junction_annotation.py

**What it does:** Classifies splice junctions as known (in BED annotation), novel, or partial.

**Output:** `07_qc/rseqc/junction_annotation/<sample>.*`

**Interpretation:**
- >80% known junctions → library accurately reflects annotated splicing
- High novel fraction → unannotated transcription or mapping artefacts

**MultiQC:** ✅ junction annotation donut charts

---

### junction_saturation.py

**What it does:** Estimates whether sequencing depth is sufficient to detect all splice junctions.

**Output:** `07_qc/rseqc/junction_saturation/<sample>.*`

**Interpretation:**
- Curve plateauing at full depth → sufficient; rising steeply → undersaturated
- Critical for alternative splicing analyses

**MultiQC:** ✅ saturation curves

---

### geneBody_coverage.py

**What it does:** Plots read coverage distribution from 5' to 3' across gene bodies.
Run once on a merged BAM from all samples.

**Output:** `07_qc/rseqc/genebody/all_samples.*`

**Interpretation:**
- Flat curve → even coverage (full-length library)
- 3' bias → poly-A capture or mRNA degradation; expected for QuantSeq REV
- 5' bias → 5'-enrichment protocol
- Irregular/bimodal → rRNA contamination or library quality issues

**MultiQC:** ✅ gene body coverage curves

---

## Output directory

```
<OUTDIR>/07_qc/
├── star/
│   ├── <sample>_Log.final.out   (symlinks)
│   └── star_alignment_summary.tsv
├── rseqc/
│   ├── infer_experiment/
│   ├── read_distribution/
│   ├── junction_annotation/
│   ├── junction_saturation/
│   └── genebody/
└── multiqc/
    └── multiQC_rseqc.html
```

---

## Strand consistency check

`check_strand_consistency.sh` verifies for each stranded sample:

```
|Fwd_reads + Rev_reads - Total_primary_mapped| / Total < STRAND_TOLERANCE_PCT
```

Default tolerance: 5%. Adjust with `STRAND_TOLERANCE_PCT="10"` in config.

**Typical fail cause:** samplesheet `strandedness` column does not match the
actual library chemistry. Fix by running `infer_experiment.py` first, then
correcting the samplesheet.

---

## Installation

```bash
conda install -c bioconda rseqc
# or
pip install RSeQC

# verify
infer_experiment.py --version
```

If installed to a non-PATH location:
```bash
RSEQC_BIN_DIR="/path/to/rseqc/bin"
```
