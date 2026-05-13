# Known issues

1. **bam_to_bedgraph.R memory** — `readGAlignmentPairs()` loads all reads. For >1B reads, set `yieldSize` on BamFile or fall back to bedtools genomecov.
2. **STAR shared memory** — After run: `STAR --genomeDir $STAR_INDEX --genomeLoad Remove`. Low-memory systems: replace `LoadAndKeep` with `NoSharedMemory`.
3. **apeglm coef naming** — Falls back to `type="normal"` if coefficient not found; informative message lists available coefficients.
4. **bedGraphToBigWig non-standard chrs** — Pre-filter to standard chromosomes if needed.
5. **UCSC tracks need HTTP URL** — Set `UCSC_BASE_URL` before run.
6. **Pandoc required** — Pipeline continues if report fails; install pandoc via conda.
7. **Seqinfo network call** — `merge_bedgraph_replicates.R` calls `Seqinfo(genome=)` on first use; install BSgenome package for offline use.
