#!/usr/bin/env Rscript
# =============================================================================
# merge_bedgraph_replicates.R — average replicate normalized bedGraphs
# =============================================================================
# ORIGIN: ADAPTED — RNA-seq/workflows/mergingRNAseqBedraphReplicatesForVisualization.sh
#   (mergeBedGraph2/3/4 inline R functions)
#   v2: uses rtracklayer import.bedGraph/export.bedGraph; sessionInfo written.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(rtracklayer); library(GenomicRanges); library(Rsamtools) })
option_list <- list(
  make_option("--samplesheet",type="character"), make_option("--bgdir",type="character"),
  make_option("--outdir",type="character"), make_option("--genome",type="character",default="mm39"),
  make_option("--layout",type="character",default="PE")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir,recursive=TRUE,showWarnings=FALSE)
ss <- read.csv(opt$samplesheet,comment.char="#",stringsAsFactors=FALSE)
if (opt$layout=="PE") colnames(ss)[1:6] <- c("sample_id","fastq_R1","fastq_R2","condition","replicate","strandedness")
else colnames(ss)[1:5] <- c("sample_id","fastq_R1","condition","replicate","strandedness")
seqinf <- tryCatch(Seqinfo(genome=opt$genome), error=function(e){message("Seqinfo unavailable for ",opt$genome);NULL})
merge_bgs <- function(files, seqinf, out_name, out_dir) {
  files <- files[file.exists(files)]; n <- length(files)
  if (n<2){message("  Skip (<2 files): ",out_name);return(invisible(NULL))}
  message("  Merging ",n," -> ",out_name)
  grl <- lapply(files, function(f){ gr <- import.bedGraph(f); if(!is.null(seqinf)) seqlevels(gr) <- seqlevels(seqinf); gr })
  all_gr <- do.call(c,grl)
  r <- disjoin(all_gr,ignore.strand=TRUE,with.revmap=TRUE)
  r$score <- vapply(r$revmap, function(idx) mean(score(all_gr)[idx]), numeric(1)); r$revmap <- NULL
  out_file <- file.path(out_dir,paste0(out_name,".bedGraph"))
  export.bedGraph(r,out_file); message("    -> ",out_file)
}
groups <- unique(ss[,c("condition","strandedness")])
for (i in seq_len(nrow(groups))) {
  cond <- groups$condition[i]; strand <- groups$strandedness[i]
  sids <- ss$sample_id[ss$condition==cond & ss$strandedness==strand]
  sfx_map <- if(strand=="unstranded") list("unstranded"="_unstranded_norm.bedGraph.gz") \
             else list("Fwd"="_Fwd_norm.bedGraph.gz","Rev"="_Rev_norm.bedGraph.gz")
  for (sfx_name in names(sfx_map)) {
    merge_bgs(file.path(opt$bgdir,paste0(sids,sfx_map[[sfx_name]])),seqinf,
              paste0(cond,"_",sfx_name,"_merged"),opt$outdir)
  }
}
writeLines(capture.output(sessionInfo()), file.path(opt$outdir,"sessionInfo_merge.txt"))
message("merge_bedgraph_replicates.R done.  Outputs in: ", opt$outdir)
