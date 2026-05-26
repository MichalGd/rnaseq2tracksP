# FastQ Screen — Species Swap & Mycoplasma Contamination Check

FastQ Screen is run as **Step 2b** — immediately after raw FastQC and before TrimGalore. It screens a subset of raw reads against a panel of reference genomes to detect sample mislabelling (species swap) and mycoplasma contamination.

---

## Why this check matters

| Problem | What it means | When it happens |
|---|---|---|
| **Species swap** | >10 % of reads map to the wrong organism | Mislabelled tubes, shared pipettes, LIMS entry errors |
| **Mycoplasma contamination** | >0.5 % of reads map to mycoplasma genomes | Common cell culture contaminant; affects gene expression globally |

Both problems are invisible in FastQC and become apparent only after alignment — when it is too late to rescue the experiment cleanly. Running FastQ Screen before trimming allows early detection and exclusion of affected samples before computation-heavy steps.

---

## Reference panel

The following five databases are screened:

| Database | Genome | Purpose |
|---|---|---|
| **Mouse** | GRCm39 (GENCODE M31) | Expected primary organism |
| **Human** | GRCh38 (GENCODE v42) | Species swap detection |
| **Zebrafish** | GRCz11 (Ensembl 112) | Species swap detection |
| **Drosophila** | BDGP6.46 (Ensembl 112) | Species swap detection |
| **Mycoplasma** | Combined 8-species genome | Contamination screen |

Bowtie2 indexes must be built once and paths set in `config/fastq_screen.conf`.

---

## Configuration

Three variables control Step 2b in `config/config.conf`:

```bash
FASTQSCREEN_CONF="config/fastq_screen.conf"   # path to database config
FASTQSCREEN_THREADS=4                          # bowtie2 threads per sample
FASTQSCREEN_SUBSET=200000                      # reads sampled per file (0 = all)
```

The step **skips gracefully** if `fastq_screen` is not in PATH or the conf file is missing — a warning is logged and the pipeline continues.

---

## Database config (`config/fastq_screen.conf`)

```
THREADS         4
BOWTIE2         bowtie2

DATABASE Mouse        /path/to/fastq_screen_db/mouse/mouse
DATABASE Human        /path/to/fastq_screen_db/human/human
DATABASE Zebrafish    /path/to/fastq_screen_db/zebrafish/zebrafish
DATABASE Drosophila   /path/to/fastq_screen_db/drosophila/drosophila
DATABASE Mycoplasma   /path/to/fastq_screen_db/mycoplasma/mycoplasma
```

### Building bowtie2 indexes (one-time setup)

```bash
# Example: mouse (GRCm39)
bowtie2-build --threads 8 GRCm39.primary_assembly.genome.fa \
  /path/to/fastq_screen_db/mouse/mouse

# Repeat for each DATABASE entry
```

---

## Output

FastQ Screen writes per-sample results to `<OUTDIR>/fastQScreen/`:

| File | Content |
|---|---|
| `<sample>_screen.txt` | Tab-separated mapping percentages per database |
| `<sample>_screen.html` | Interactive per-sample bar chart |
| `<sample>_screen.png` | Static plot |

All results are automatically included in the **final MultiQC report** (`multiQC/final/multiQC_final.html`) as a stacked bar chart panel.

---

## Interpreting results

### Expected clean mouse sample

| Database | % Mapped (one hit) |
|---|---|
| Mouse | ~90–97 % |
| Human | < 1 % |
| Zebrafish | < 0.5 % |
| Drosophila | < 0.5 % |
| Mycoplasma | **< 0.1 %** |

### Action thresholds

| Observation | Threshold | Action |
|---|---|---|
| Non-expected species > 5 % | Warning | Investigate sample provenance; check LIMS |
| Non-expected species > 20 % | Fail | Exclude sample; do not proceed to alignment |
| Mycoplasma > 0.5 % | Warning | Flag sample; notify cell culture team |
| Mycoplasma > 2 % | Fail | Exclude sample; confirmed contamination |

---

## Manual single-sample test

```bash
fastq_screen \
  --conf config/fastq_screen.conf \
  --outdir /tmp/fqs_test \
  --threads 4 \
  --subset 200000 \
  --aligner bowtie2 \
  /path/to/sample_R1.fq.gz

cat /tmp/fqs_test/*_screen.txt
```

---

## See also

- `config/fastq_screen.conf` — database paths and thread settings
- `environment.yml` — `fastq-screen` and `bowtie2` conda packages
- `docs/INSTALLATION.md` — full index building instructions
