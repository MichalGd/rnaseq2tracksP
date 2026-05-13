#!/usr/bin/env Rscript
# =============================================================================
# deseq2_qc_plots.R — PCA, sample-distance heatmap, top-50 heatmap
# =============================================================================
# ORIGIN: ADAPTED — RNA-seq/workflows/SampleOverview_distances_PCA_heatmaps.sh
#   v2: sessionInfo written.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(DESeq2); library(vsn); library(pheatmap); library(RColorBrewer) })
option_list <- list(make_option("--countsrdata",type="character"), make_option("--outdir",type="character"))
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir,recursive=TRUE,showWarnings=FALSE)
load(opt$countsrdata); dds <- estimateSizeFactors(dds); vsd <- vst(dds,blind=FALSE)
pdf(file.path(opt$outdir,"PCA.pdf"),height=6,width=7)
print(plotPCA(vsd,intgroup=c("condition","replicate"))); dev.off()
sampleDists <- dist(t(assay(vsd))); sdm <- as.matrix(sampleDists)
rownames(sdm) <- paste(vsd$condition,vsd$replicate,sep="_"); colnames(sdm) <- NULL
pdf(file.path(opt$outdir,"sample_clustering.pdf"),height=6,width=8)
pheatmap(sdm,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,
  col=colorRampPalette(rev(brewer.pal(9,"Blues")))(255)); dev.off()
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
ann <- as.data.frame(colData(dds)[,"condition",drop=FALSE])
pdf(file.path(opt$outdir,"top50_heatmap.pdf"),height=10,width=8)
pheatmap(assay(ntd)[select,],cluster_rows=TRUE,show_rownames=TRUE,cluster_cols=FALSE,
  annotation_col=ann,labels_row=geneInfo[rownames(assay(ntd)[select,]),"gene_name"]); dev.off()
writeLines(capture.output(sessionInfo()), file.path(opt$outdir,"sessionInfo_qc.txt"))
message("deseq2_qc_plots.R done.  Outputs in: ", opt$outdir)
