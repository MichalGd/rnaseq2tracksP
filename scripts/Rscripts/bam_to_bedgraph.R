#!/usr/bin/env Rscript
# =============================================================================
# bam_to_bedgraph.R
# =============================================================================
# ORIGIN: NEW (v2) — replaces strand_split_bedgraph.sh from v1.
#   Rsamtools + GenomicAlignments + rtracklayer.
#   No intermediate BAMs; correct PE pair orientation.
#
# Outputs:
#   <sample_id>_unstranded.bedGraph.gz    (strandedness == "unstranded")
#   <sample_id>_Fwd.bedGraph.gz           (stranded)
#   <sample_id>_Rev.bedGraph.gz           (stranded)
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(Rsamtools)
  library(GenomicAlignments); library(rtracklayer)
})
option_list <- list(
  make_option("--samplesheet", type="character"),
  make_option("--bamdir",      type="character"),
  make_option("--outdir",      type="character"),
  make_option("--layout",      type="character", default="PE")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

ss <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)
if (opt$layout == "PE") {
  colnames(ss)[1:6] <- c("sample_id","fastq_R1","fastq_R2","condition","replicate","strandedness")
} else {
  colnames(ss)[1:5] <- c("sample_id","fastq_R1","condition","replicate","strandedness")
}

export_gz <- function(gr, out_gz) {
  tmp <- sub("\\.gz$", "", out_gz)
  export.bedGraph(gr, tmp)
  system2("gzip", c("-f", tmp))
  message("  Written: ", basename(out_gz))
}

cov_to_gr <- function(cov) {
  gr <- GRanges(as(cov, "GRanges"))
  gr[score(gr) != 0]
}

for (i in seq_len(nrow(ss))) {
  sid    <- ss$sample_id[i]
  strand <- ss$strandedness[i]
  bam    <- file.path(opt$bamdir, paste0(sid, "_sortedS.bam"))
  if (!file.exists(bam)) { warning("BAM missing: ", bam); next }
  message("Coverage: ", sid, "  strandedness=", strand)

  bf     <- BamFile(bam, asMates=(opt$layout == "PE"))
  param  <- ScanBamParam(flag=scanBamFlag(
    isSecondaryAlignment=FALSE, isSupplementaryAlignment=FALSE,
    isProperPair=if(opt$layout == "PE") TRUE else NA))

  if (strand == "unstranded") {
    ga  <- if (opt$layout == "PE") as(readGAlignmentPairs(bf, param=param), "GAlignments") \
           else readGAlignments(bf, param=param)
    export_gz(cov_to_gr(coverage(ga)),
              file.path(opt$outdir, paste0(sid, "_unstranded.bedGraph.gz")))
  } else {
    ga <- if (opt$layout == "PE") as(readGAlignmentPairs(bf, param=param), "GAlignments") \
          else readGAlignments(bf, param=param)
    export_gz(cov_to_gr(coverage(ga[strand(ga) == "+"])),
              file.path(opt$outdir, paste0(sid, "_Fwd.bedGraph.gz")))
    export_gz(cov_to_gr(coverage(ga[strand(ga) == "-"])),
              file.path(opt$outdir, paste0(sid, "_Rev.bedGraph.gz")))
  }
}
message("bam_to_bedgraph.R done.  Outputs in: ", opt$outdir)
