# =============================================================================
# deseq2_qc_plots.R — PCA, sample clustering, heatmaps
# ORIGIN: ADAPTED from RNA-seq repo QC plot scripts
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(ggplot2)
  library(pheatmap); library(RColorBrewer); library(vsn)
})
option_list <- list(
  make_option("--countsrdata"),
  make_option("--outdir")
)
opt <- parse_args(OptionParser(option_list=option_list))

save_dual <- function(p, path_pdf, w=7, h=6) {
  ggsave(path_pdf, p, width=w, height=h)
  ggsave(sub("\\.pdf$", ".png", path_pdf), p, width=w, height=h, dpi=300)
}
load(opt$countsrdata)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)

# PCA
pca <- plotPCA(vsd, intgroup=c("condition","replicate"), returnData=TRUE)
pca$replicate <- as.factor(pca$replicate)
pct <- round(100*attr(pca,"percentVar"), 1)
p <- ggplot(pca, aes(PC1, PC2, color=condition, shape=replicate, label=name)) +
  geom_point(size=3) + ggrepel::geom_text_repel(size=3) +
  xlab(paste0("PC1 (",pct[1],"%)")) + ylab(paste0("PC2 (",pct[2],"%)")) +
  theme_bw()
ggsave(file.path(opt$outdir,"PCA.pdf"), p, width=7, height=5)

# Sample correlation heatmap
sampleDists <- dist(t(assay(vsd)))
mat <- as.matrix(sampleDists)
pheatmap(mat,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colorRampPalette(rev(brewer.pal(9,"Blues")))(255),
  filename=file.path(opt$outdir,"sample_clustering.pdf"))

# Top 500 variable genes heatmap
rv   <- rowVars(assay(vsd))
top  <- head(order(rv, decreasing=TRUE), 500)
pheatmap(assay(vsd)[top,],
  scale="row", show_rownames=FALSE,
  annotation_col=as.data.frame(colData(vsd)[,c("condition","replicate")]),
  filename=file.path(opt$outdir,"top500_variable_genes_heatmap.pdf"))

# Mean-SD plot
pdf(file.path(opt$outdir,"meanSD_plot.pdf"))
meanSdPlot(assay(vsd))
dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,"deseq2_qc_plots_sessionInfo.txt"))
message("[deseq2_qc_plots.R] Done.")
