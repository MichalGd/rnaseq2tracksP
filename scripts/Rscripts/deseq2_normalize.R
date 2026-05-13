#!/usr/bin/env Rscript
# =============================================================================
# deseq2_normalize.R
# ORIGIN: ADAPTED — RNA-seq/workflows/GeneratingNormalizedBedgraphsDEseq2.sh
#   v2: sessionInfo written; rtracklayer/GenomicFeatures unchanged.
# Outputs: raw_counts.tsv, normalized_counts.tsv, fpkm_counts.tsv,
#          tpm_counts.tsv, size_factors.tsv, library_stats.tsv,
#          dds.RData, sessionInfo_normalize.txt
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(data.table)
  library(rtracklayer); library(GenomicFeatures)
})
option_list <- list(
  make_option("--samplesheet", type="character"),
  make_option("--countdir",    type="character"),
  make_option("--gtf",         type="character"),
  make_option("--layout",      type="character", default="PE"),
  make_option("--outdir",      type="character"),
  make_option("--design",      type="character", default="~ condition")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

ss <- read.csv(opt$samplesheet, comment.char="#", stringsAsFactors=FALSE)
if (opt$layout == "PE") {
  colnames(ss)[1:6] <- c("sample_id","fastq_R1","fastq_R2","condition","replicate","strandedness")
} else {
  colnames(ss)[1:5] <- c("sample_id","fastq_R1","condition","replicate","strandedness")
}
strand_col <- function(s) switch(s,"unstranded"=2L,"forward"=3L,"reverse"=4L,stop(paste("Unknown:",s)))

countData <- NULL
for (i in seq_len(nrow(ss))) {
  sid <- ss$sample_id[i]
  tab <- file.path(opt$countdir, paste0(sid,"_ReadsPerGene.out.tab"))
  if (!file.exists(tab)) stop("File not found: ", tab)
  col <- strand_col(ss$strandedness[i])
  dt  <- data.frame(fread(tab))[, c(1, col)]
  colnames(dt) <- c("gene_id", sid)
  countData <- if (is.null(countData)) dt else merge(countData, dt, by="gene_id", all=TRUE)
}
gids <- countData$gene_id; rownames(countData) <- gids; countData$gene_id <- NULL
countDataHead <- countData[1:4,,drop=FALSE]
countData <- countData[5:nrow(countData),,drop=FALSE]
countData <- apply(countData, 2, as.integer); rownames(countData) <- gids[5:length(gids)]

coldata <- data.frame(condition=factor(ss$condition), replicate=factor(ss$replicate), row.names=ss$sample_id)
dds <- DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=as.formula(opt$design))
dds <- estimateSizeFactors(dds); SF <- sizeFactors(dds)
rawCounts <- as.data.frame(counts(dds,normalized=FALSE))
normCounts <- as.data.frame(counts(dds,normalized=TRUE))

message("Loading GTF...")
gtf <- import(opt$gtf); gtf_genes <- gtf[gtf$type=="gene"]
geneInfo <- as.data.frame(gtf_genes[,c("gene_id","gene_name","gene_type")],stringsAsFactors=FALSE)
geneInfo$seqnames <- as.character(seqnames(gtf_genes))
geneInfo$start <- start(gtf_genes); geneInfo$end <- end(gtf_genes)
rownames(geneInfo) <- geneInfo$gene_id
txdb <- suppressMessages(makeTxDbFromGFF(opt$gtf))
tx_lens <- transcriptLengths(txdb)
med_tx <- tapply(tx_lens$tx_len, tx_lens$gene_id, median)
geneInfo$medianTxLen <- med_tx[geneInfo$gene_id]
geneInfo$medianTxLen[is.na(geneInfo$medianTxLen)] <- 1000

fpkmCounts <- sweep(normCounts, 1, geneInfo[rownames(normCounts),"medianTxLen"]/1000, "/")
rpk <- sweep(rawCounts, 1, geneInfo[rownames(rawCounts),"medianTxLen"]/1000, "/")
tpmCounts <- sweep(rpk, 2, colSums(rpk)/1e6, "/")

exonicReads <- colSums(rawCounts); featuredNorm <- exonicReads/SF
SF_rpm <- SF * mean(featuredNorm/(colSums(rawCounts)/1e6))

sf_table <- data.frame(sample_id=ss$sample_id, condition=ss$condition, replicate=ss$replicate,
  strandedness=ss$strandedness, SF=SF[ss$sample_id], SF_rpm=SF_rpm[ss$sample_id],
  exonic_reads=exonicReads[ss$sample_id], row.names=NULL)

statsTable <- rbind(exonicReads, apply(countDataHead,2,as.numeric), exonicReads/SF, SF, SF_rpm)
rownames(statsTable) <- c("exonic_reads","N_unmapped","N_multimapping","N_noFeature","N_ambiguous","exonicNorm","SF","SF_rpm")

write.table(cbind(geneInfo[rownames(rawCounts),],rawCounts),  file.path(opt$outdir,"raw_counts.tsv"),        sep="\t",quote=FALSE,row.names=FALSE)
write.table(cbind(geneInfo[rownames(normCounts),],normCounts),file.path(opt$outdir,"normalized_counts.tsv"),  sep="\t",quote=FALSE,row.names=FALSE)
write.table(cbind(geneInfo[rownames(fpkmCounts),],fpkmCounts),file.path(opt$outdir,"fpkm_counts.tsv"),        sep="\t",quote=FALSE,row.names=FALSE)
write.table(cbind(geneInfo[rownames(tpmCounts),],tpmCounts),  file.path(opt$outdir,"tpm_counts.tsv"),         sep="\t",quote=FALSE,row.names=FALSE)
write.table(sf_table,          file.path(opt$outdir,"size_factors.tsv"),  sep="\t",quote=FALSE,row.names=FALSE)
write.table(t(statsTable),     file.path(opt$outdir,"library_stats.tsv"), sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)
save(dds,geneInfo,ss, file=file.path(opt$outdir,"dds.RData"))
writeLines(capture.output(sessionInfo()), file.path(opt$outdir,"sessionInfo_normalize.txt"))
message("deseq2_normalize.R done.  Outputs in: ", opt$outdir)
