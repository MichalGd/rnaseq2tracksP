# =============================================================================
# bam_to_bedgraph.R — strand-aware coverage from BAM to bedGraph
# ORIGIN: NEW v2 — replaces Bash strand-split approach
# Language decision: Rsamtools + GenomicAlignments provides strand-correct
#   PE handling without split BAMs; rtracklayer export.bedGraph for output
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(Rsamtools); library(GenomicAlignments)
  library(rtracklayer); library(GenomicRanges)
})
option_list <- list(
  make_option("--samplesheet"), make_option("--bamdir"),
  make_option("--outdir"),      make_option("--layout", default="PE")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

ss <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)
layout <- opt$layout

for (i in seq_len(nrow(ss))) {
  sid      <- ss$sample_id[i]
  strand_  <- if (layout=="PE") ss$strandedness[i] else ss$strandedness[i]
  bam_path <- file.path(opt$bamdir, paste0(sid, "_sortedS.bam"))
  if (!file.exists(bam_path)) { warning(sid,": BAM missing, skip"); next }

  message("[bam_to_bedgraph.R] ", sid, "  strand=", strand_)
  bf <- BamFile(bam_path)

  if (strand_ == "unstranded") {
    if (layout == "PE")
      reads <- readGAlignmentPairs(bf, param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)))
    else
      reads <- readGAlignments(bf, param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)))
    cov <- coverage(reads)
    bg <- rtracklayer::export(cov, connection=NULL, format="bedGraph")
    bg_gr <- as(bg, "GRanges")
    bg_df <- as.data.frame(bg_gr)[,c("seqnames","start","end","score")]
    bg_df$start <- bg_df$start - 1
    bg_df <- bg_df[bg_df$score > 0, ]
    out <- file.path(opt$outdir, paste0(sid, "_unstranded.bedGraph"))
    write.table(bg_df, out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    R.utils::gzip(out, overwrite=TRUE)

  } else {
    # Strand-specific: forward and reverse
    for (sg in c("Fwd","Rev")) {
      if (layout == "PE") {
        # dUTP / reverse library: R2 on RNA strand
        #   Fwd (+) strand: R2 NOT reverse  OR  R1 reverse
        #   Rev (-) strand: the complement
        if (sg == "Fwd") {
          flag_1 <- scanBamFlag(isSecondaryAlignment=FALSE, isFirstMateRead=TRUE,  isMinusStrand=TRUE)
          flag_2 <- scanBamFlag(isSecondaryAlignment=FALSE, isFirstMateRead=FALSE, isMinusStrand=FALSE)
        } else {
          flag_1 <- scanBamFlag(isSecondaryAlignment=FALSE, isFirstMateRead=TRUE,  isMinusStrand=FALSE)
          flag_2 <- scanBamFlag(isSecondaryAlignment=FALSE, isFirstMateRead=FALSE, isMinusStrand=TRUE)
        }
        r1 <- readGAlignments(bf, param=ScanBamParam(flag=flag_1))
        r2 <- readGAlignments(bf, param=ScanBamParam(flag=flag_2))
        reads_gr <- c(granges(r1), granges(r2))
      } else {
        if (sg == "Fwd")
          flag_ <- scanBamFlag(isSecondaryAlignment=FALSE, isMinusStrand=TRUE)
        else
          flag_ <- scanBamFlag(isSecondaryAlignment=FALSE, isMinusStrand=FALSE)
        reads_gr <- granges(readGAlignments(bf, param=ScanBamParam(flag=flag_)))
      }
      cov <- coverage(reads_gr)
      bg_df <- as.data.frame(GRanges(cov))
      bg_df <- bg_df[bg_df$score > 0, c("seqnames","start","end","score")]
      bg_df$start <- bg_df$start - 1
      if (sg == "Rev") bg_df$score <- -abs(bg_df$score)
      out <- file.path(opt$outdir, paste0(sid, "_", sg, "S.bedGraph"))
      write.table(bg_df, out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
      R.utils::gzip(out, overwrite=TRUE)
    }
  }
}
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir, "bam_to_bedgraph_sessionInfo.txt"))
message("[bam_to_bedgraph.R] Done.")
