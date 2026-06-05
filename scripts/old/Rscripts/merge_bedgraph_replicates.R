# =============================================================================
# merge_bedgraph_replicates.R — average replicates: ONE condition per call
# ORIGIN: ADAPTED v4→v4.2 — condition loop moved to Bash
# Called by: rnaseq2tracks.sh Step 14 (parallel, 1 job per condition)
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
})

option_list <- list(
  make_option("--condition",    type="character", help="Condition label"),
  make_option("--sample_ids",   type="character",
              help="Comma-separated sample IDs belonging to this condition"),
  make_option("--strandedness", type="character", default="reverse",
              help="unstranded / forward / reverse (all replicates same)"),
  make_option("--bgdir",        type="character", help="Normalized bedGraph dir"),
  make_option("--outdir",       type="character", help="Output dir for merged bedGraphs"),
  make_option("--layout",       type="character", default="PE")
)
opt   <- parse_args(OptionParser(option_list = option_list))
sids  <- strsplit(opt$sample_ids, ",")[[1]]

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

suffixes <- if (opt$strandedness == "unstranded") "_unstranded_norm" else c("_FwdS_norm", "_RevS_norm")

for (sfx in suffixes) {
  message("[merge_bedgraph_replicates.R] condition=", opt$condition,
          "  suffix=", sfx, "  n_reps=", length(sids))
  bgs <- lapply(sids, function(s) {
    f <- file.path(opt$bgdir, paste0(s, sfx, ".bedGraph.gz"))
    if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
    import(f, format = "bedGraph")
  })
  bgs <- Filter(Negate(is.null), bgs)
  if (length(bgs) == 0) { warning("No bedGraphs for ", opt$condition, sfx); next }

  if (length(bgs) == 1) {
    merged <- bgs[[1]]
  } else {
    dj <- disjoin(do.call(c, bgs))
    scores <- sapply(bgs, function(bg) {
      ov <- findOverlaps(dj, bg)
      s  <- rep(0, length(dj))
      s[queryHits(ov)] <- bg$score[subjectHits(ov)]
      s
    })
    dj$score <- rowMeans(scores)
    merged   <- dj
  }

  out <- file.path(opt$outdir, paste0(opt$condition, sfx, "_merged.bedGraph"))
  export(merged, con = out, format = "bedGraph")
  message("  Written: ", basename(out))
}

writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,
                     paste0(opt$condition, "_merge_bedgraph_sessionInfo.txt")))
message("[merge_bedgraph_replicates.R] Done: ", opt$condition)
