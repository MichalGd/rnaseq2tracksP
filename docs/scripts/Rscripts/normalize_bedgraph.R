# =============================================================================
# normalize_bedgraph.R — SF_rpm scaling: ONE sample per invocation
# ORIGIN: ADAPTED v4→v4.2 — loop moved to Bash; script now single-sample
# Called by: rnaseq2tracks.sh Step 12 (parallel, 1 job per sample)
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(data.table)
})

option_list <- list(
  make_option("--sample_id",    type="character", help="Sample ID"),
  make_option("--strandedness", type="character", default="reverse",
              help="unstranded / forward / reverse"),
  make_option("--sffile",       type="character", help="size_factors.tsv path"),
  make_option("--rawbgdir",     type="character", help="Raw bedGraph directory"),
  make_option("--outdir",       type="character", help="Output directory"),
  make_option("--layout",       type="character", default="PE")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

sf_tbl <- fread(opt$sffile, sep = "\t")
sf_rpm <- sf_tbl$sf_rpm[sf_tbl$sample_id == opt$sample_id]
if (length(sf_rpm) == 0 || is.na(sf_rpm))
  stop("[normalize_bedgraph.R] sf_rpm not found for sample: ", opt$sample_id)

message("[normalize_bedgraph.R] ", opt$sample_id,
        "  sf_rpm=", round(sf_rpm, 4),
        "  strand=", opt$strandedness)

suffixes <- if (opt$strandedness == "unstranded") "_unstranded" else c("_FwdS", "_RevS")

for (sfx in suffixes) {
  bgz <- file.path(opt$rawbgdir, paste0(opt$sample_id, sfx, ".bedGraph.gz"))
  if (!file.exists(bgz)) {
    warning("[normalize_bedgraph.R] missing: ", bgz, " — skipping")
    next
  }
  bg <- import(bgz, format = "bedGraph")
  bg$score <- bg$score / sf_rpm
  out <- file.path(opt$outdir, paste0(opt$sample_id, sfx, "_norm.bedGraph.gz"))
  export(bg, con = out, format = "bedGraph")
  message("  Written: ", basename(out))
}

writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,
                     paste0(opt$sample_id, "_normalize_bedgraph_sessionInfo.txt")))
message("[normalize_bedgraph.R] Done: ", opt$sample_id)
