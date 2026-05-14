# =============================================================================
# merge_bedgraph_replicates.R — average replicates per condition via GRanges
# ORIGIN: ADAPTED from RNA-seq repo replicate merge logic
# Language decision: GRanges disjoin + mean avoids sort/bedtools issues
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges); library(data.table)
})
option_list <- list(
  make_option("--samplesheet"),
  make_option("--bgdir"),
  make_option("--outdir"),
  make_option("--genome",  default="mm39"),
  make_option("--layout",  default="PE")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

ss <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)
conditions <- unique(ss$condition)

suffixes_for <- function(strand) {
  if (strand=="unstranded") "_unstranded_norm" else c("_FwdS_norm","_RevS_norm")
}

for (cond in conditions) {
  sids <- ss$sample_id[ss$condition == cond]
  strand_types <- unique(ss$strandedness[ss$condition == cond])

  for (sfx in suffixes_for(strand_types[1])) {
    bgs <- lapply(sids, function(s) {
      f <- file.path(opt$bgdir, paste0(s, sfx, ".bedGraph.gz"))
      if (!file.exists(f)) return(NULL)
      import(f, format="bedGraph")
    })
    bgs <- Filter(Negate(is.null), bgs)
    if (length(bgs) == 0) next

    # GRanges disjoin + mean
    all_gr <- do.call(c, bgs)
    dj     <- disjoin(all_gr)
    scores <- sapply(bgs, function(bg) {
      ov <- findOverlaps(dj, bg)
      s  <- rep(0, length(dj))
      s[queryHits(ov)] <- bg$score[subjectHits(ov)]
      s
    })
    dj$score <- if (is.matrix(scores)) rowMeans(scores) else mean(scores)

    out <- file.path(opt$outdir, paste0(cond, sfx, "_merged.bedGraph"))
    export(dj, con=out, format="bedGraph")
    message("[merge_bedgraph_replicates.R] ", basename(out))
  }
}
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,"merge_bedgraph_sessionInfo.txt"))
message("[merge_bedgraph_replicates.R] Done.")
