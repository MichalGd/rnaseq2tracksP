#!/usr/bin/env Rscript
# =============================================================================
# deseq2_de.R — DESeq2 Wald test + apeglm LFC shrinkage per contrast
# =============================================================================
# ORIGIN: ADAPTED — RNA-seq/workflows/DifferentialGeneExpressionWithDEseq2.sh
#   v2: informative apeglm fallback message; sessionInfo written.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(DESeq2); library(apeglm); library(ggplot2) })
option_list <- list(
  make_option("--countsrdata", type="character"),
  make_option("--contrasts",   type="character"),
  make_option("--outdir",      type="character"),
  make_option("--padj",  type="double", default=0.05),
  make_option("--lfc",   type="double", default=1.0)
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
load(opt$countsrdata); dds <- DESeq(dds)
contrasts_df <- read.csv(opt$contrasts, comment.char="#", stringsAsFactors=FALSE)
colnames(contrasts_df)[1:3] <- c("contrast_id","condition_a","condition_b")
for (i in seq_len(nrow(contrasts_df))) {
  cid <- contrasts_df$contrast_id[i]; ca <- contrasts_df$condition_a[i]; cb <- contrasts_df$condition_b[i]
  message("Contrast: ",cid," (",ca," vs ",cb,")")
  res <- results(dds, contrast=c("condition",ca,cb), alpha=opt$padj)
  coef_name <- paste0("condition_",ca,"_vs_",cb)
  res_lfc <- tryCatch(lfcShrink(dds,coef=coef_name,type="apeglm",quiet=TRUE),
    error=function(e){
      message("  apeglm coef '",coef_name,"' not found. Available: ",paste(resultsNames(dds),collapse=", "),
              "\n  Using type='normal'.")
      lfcShrink(dds,contrast=c("condition",ca,cb),type="normal",quiet=TRUE)})
  res_df <- as.data.frame(res,stringsAsFactors=FALSE); res_df$gene_id <- rownames(res_df)
  res_df$LFC_shrunken <- res_lfc$log2FoldChange; res_df$LFC_padj <- res_lfc$padj
  res_df$padj[is.na(res_df$padj)] <- 1
  out_df <- merge(geneInfo[,c("gene_id","gene_name","gene_type","seqnames","start","end")],res_df,by="gene_id",all.y=TRUE)
  write.table(out_df, file.path(opt$outdir,paste0(cid,"_DE_results.tsv")), sep="\t",quote=FALSE,row.names=FALSE)
  pdf(file.path(opt$outdir,paste0(cid,"_MA_plot.pdf")),width=6,height=5)
  plotMA(res,alpha=opt$padj,main=paste("MA:",cid)); dev.off()
  plot_df <- out_df[!is.na(out_df$padj)&!is.na(out_df$log2FoldChange),]
  plot_df$sig <- ifelse(plot_df$padj<opt$padj&abs(plot_df$log2FoldChange)>=opt$lfc,"significant","ns")
  p <- ggplot(plot_df,aes(x=log2FoldChange,y=-log10(padj),color=sig))+
    geom_point(size=0.5,alpha=0.6)+scale_color_manual(values=c("significant"="red","ns"="grey60"))+
    geom_vline(xintercept=c(-opt$lfc,opt$lfc),linetype="dashed")+
    geom_hline(yintercept=-log10(opt$padj),linetype="dashed")+
    labs(title=paste("Volcano:",cid),x="log2FC",y="-log10(padj)")+theme_bw(base_size=12)
  ggsave(file.path(opt$outdir,paste0(cid,"_volcano.pdf")),p,width=6,height=5)
  message("  Significant: ",sum(plot_df$sig=="significant",na.rm=TRUE))
}
writeLines(capture.output(sessionInfo()), file.path(opt$outdir,"sessionInfo_DE.txt"))
message("deseq2_de.R done.  Outputs in: ", opt$outdir)
