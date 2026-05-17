# =============================================================================
# deseq2_normalize.R — DESeq2 size factor + SF_rpm normalization
# ORIGIN: ADAPTED from RNA-seq repo normalizedBedGraphs scripts
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(GenomicFeatures); library(txdbmaker)
  library(data.table)
})
option_list <- list(
  make_option("--samplesheet"), make_option("--countdir"),
  make_option("--gtf"),         make_option("--layout",   default="PE"),
  make_option("--outdir"),      make_option("--design",   default="~ condition")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

ss   <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)
sids <- ss$sample_id
col  <- if (opt$layout=="PE") 4L else 2L  # STAR unstranded col

# Load ReadsPerGene
count_list <- lapply(sids, function(s) {
  f <- file.path(opt$countdir, paste0(s,"_ReadsPerGene.out.tab"))
  d <- fread(f, skip=4, header=FALSE, col.names=c("gene","unstranded","fwd","rev"))
  # per-sample strand column selection
  strand_ <- ss$strandedness[ss$sample_id==s]
  sc <- switch(strand_, unstranded=2L, forward=3L, reverse=4L, 2L)
  out <- d[, .(gene, count=.SD[[sc]]), .SDcols=names(d)]
  setnames(out, "count", s)
  out
})
counts_mat <- Reduce(function(a,b) merge(a,b,by="gene"), count_list)
rn <- counts_mat$gene; counts_mat <- as.matrix(counts_mat[,-1])
rownames(counts_mat) <- rn; colnames(counts_mat) <- sids

coldata <- data.frame(
  row.names   = ss$sample_id,
  condition   = ss$condition,
  replicate   = ss$replicate,
  stringsAsFactors = FALSE
)
dds <- DESeqDataSetFromMatrix(counts_mat, coldata, as.formula(opt$design))
dds <- estimateSizeFactors(dds)
SF  <- sizeFactors(dds)

# Mean exonic RPM anchor for SF_rpm
txdb   <- txdbmaker::makeTxDbFromGFF(opt$gtf)
exons  <- exonsBy(txdb, by="gene")
glen   <- sum(width(reduce(exons))) / 1000
total_mapped <- colSums(counts_mat)
rpm    <- sweep(counts_mat, 2, total_mapped/1e6, "/")
mean_rpm <- rowMeans(rpm)
SF_rpm <- SF * mean(mean_rpm[mean_rpm > 1])   # anchor

# Outputs
fwrite(data.table(sample_id=sids, size_factor=SF, sf_rpm=SF_rpm),
       file.path(opt$outdir,"size_factors.tsv"), sep="\t")
raw <- as.data.frame(counts(dds))
fwrite(cbind(gene=rownames(raw), raw),
       file.path(opt$outdir,"raw_counts.tsv"), sep="\t")
norm <- as.data.frame(counts(dds, normalized=TRUE))
fwrite(cbind(gene=rownames(norm), norm),
       file.path(opt$outdir,"normalized_counts.tsv"), sep="\t")
save(dds, file=file.path(opt$outdir,"dds.RData"))
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,"deseq2_normalize_sessionInfo.txt"))
message("[deseq2_normalize.R] Done.")
