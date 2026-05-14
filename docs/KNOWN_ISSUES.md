# Known issues

1. **STAR shared memory**: If pipeline is interrupted, shared memory may not be released.
   Fix: `STAR --genomeDir $STAR_INDEX --genomeLoad Remove`

2. **apeglm convergence warning**: For very small experiments (<3 replicates/condition),
   apeglm may not converge. Script falls back to `type="normal"` automatically.

3. **geneBody_coverage memory**: Merging large BAMs (>30 samples, >30M reads each)
   may require >64 GB RAM. Reduce `MAX_JOBS` or run separately.

4. **RSeQC BED file format**: Must be BED12. BED6 will cause RSeQC to report all reads
   as intergenic. Verify with `awk 'NF==12' annotation.bed | wc -l`.

5. **CHROMOSOME_NAMING mismatch**: If BigWig chromosomes don't match chrom.sizes,
   bedGraphToBigWig will exit. Verify alignment: `head -1 annotation.gtf` vs
   `head -1 chrom.sizes`.

6. **TrimGalore PE output naming**: TrimGalore uses `--basename` to produce
   `${SID}_val_1.fq.gz` / `${SID}_val_2.fq.gz`. If R1/R2 filenames are
   non-standard, check trimgalore_single.sh `--basename` handling.

7. **junction_annotation.py on 3' mRNA-seq**: Very few junctions expected for
   QuantSeq / 3'-end libraries. High "unknown" fraction is normal; low overall
   junction count is expected.

8. **strand_check false positives with low coverage**: Very low-coverage samples
   (<1M mapped reads) may have noisy strand ratios. Increase STRAND_TOLERANCE_PCT
   to 15 for pilot experiments.

9. **R.utils::gzip not installed**: bam_to_bedgraph.R uses R.utils for gzip.
   Install with: `install.packages("R.utils")`

10. **MultiQC version**: RSeQC module support requires MultiQC >= 1.14.
    Check: `multiqc --version`
