# =============================================================================
# normalize_bedgraph.R — apply SF_rpm to bedGraph using rtracklayer
# ORIGIN: ADAPTED v2 — rtracklayer replaces Bash awk arithmetic
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(data.table)
})
option_list <- list(
  make_option("--samplesheet"), make_option("--sffile"),
  make_option("--rawbgdir"),    make_option("--outdir"),
  make_option("--layout", default="PE")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

sf_tbl <- fread(opt$sffile, sep="\t")
ss     <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)

for (sid in ss$sample_id) {
  sf_rpm <- sf_tbl$sf_rpm[sf_tbl$sample_id == sid]
  strand_ <- ss$strandedness[ss$sample_id == sid]
  suffixes <- if (strand_ == "unstranded") "_unstranded" else c("_FwdS","_RevS")

  for (sfx in suffixes) {
    bgz <- file.path(opt$rawbgdir, paste0(sid, sfx, ".bedGraph.gz"))
    if (!file.exists(bgz)) { warning(bgz," missing, skip"); next }
    bg  <- import(bgz, format="bedGraph")
    bg$score <- bg$score / sf_rpm
    out <- file.path(opt$outdir, paste0(sid, sfx, "_norm.bedGraph.gz"))
    export(bg, con=out, format="bedGraph")
    message("[normalize_bedgraph.R] ", basename(out))
  }
}
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,"normalize_bedgraph_sessionInfo.txt"))
message("[normalize_bedgraph.R] Done.")
