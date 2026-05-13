# Usage

```bash
./scripts/rnaseq2tracks.sh config/config.conf
```

After run:
```bash
STAR --genomeDir $STAR_INDEX --genomeLoad Remove  # free shared memory
open $OUTDIR/reports/pipeline_report.html
```
