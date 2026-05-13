#!/usr/bin/env Rscript
# =============================================================================
# normalize_bedgraph.R — SF_rpm scaling via rtracklayer import/export
# =============================================================================
# ORIGIN: ADAPTED — RNA-seq/workflows/GeneratingNormalizedBedgraphsDEseq2.sh
#   v2: switched from data.table fread/fwrite to rtracklayer import.bedGraph /
#   export.bedGraph for type-safe GRanges arithmetic.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(rtracklayer) })
option_list <- list(
  make_option("--samplesheet", type="character"),
  make_option("--sffile",      type="character"),
  make_option("--rawbgdir",    type="character"),
  make_option("--outdir",      type="character"),
  make_option("--layout",      type="character", default="PE")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

ss <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)
if (opt$layout=="PE") colnames(ss)[1:6] <- c("sample_id","fastq_R1","fastq_R2","condition","replicate","strandedness")
else                  colnames(ss)[1:5] <- c("sample_id","fastq_R1","condition","replicate","strandedness")

sf_df <- read.table(opt$sffile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(sf_df) <- sf_df$sample_id

norm_bg <- function(bg_gz, sf_rpm, out_gz) {
  gr <- import.bedGraph(bg_gz); score(gr) <- score(gr)/sf_rpm
  tmp <- sub("\\.gz$","",out_gz); export.bedGraph(gr, tmp)
  system2("gzip", c("-f", tmp)); message("  Normalized: ", basename(out_gz))
}

for (i in seq_len(nrow(ss))) {
  sid <- ss$sample_id[i]; strand <- ss$strandedness[i]
  sf_rpm <- sf_df[sid,"SF_rpm"]
  if (is.na(sf_rpm)||sf_rpm==0) { warning("SF_rpm missing for ",sid); next }
  sfx_list <- if (strand=="unstranded") c("_unstranded.bedGraph.gz") \
              else c("_Fwd.bedGraph.gz","_Rev.bedGraph.gz")
  for (sfx in sfx_list) {
    bg_in  <- file.path(opt$rawbgdir, paste0(sid, sfx))
    bg_out <- file.path(opt$outdir,   sub("\\.bedGraph\\.gz$","_norm.bedGraph.gz",basename(bg_in)))
    if (file.exists(bg_in)) norm_bg(bg_in, sf_rpm, bg_out)
    else warning("Not found: ", bg_in)
  }
}
message("normalize_bedgraph.R done.  Outputs in: ", opt$outdir)
