# =============================================================================
# deseq2_de.R — Wald test + apeglm LFC shrinkage per contrast
# ORIGIN: ADAPTED — from RNA-seq repo DE analysis scripts
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(apeglm)
  library(ggplot2); library(data.table)
})
option_list <- list(
  make_option("--countsrdata"),
  make_option("--contrasts"),
  make_option("--outdir"),
  make_option("--padj",  default="0.05"),
  make_option("--lfc",   default="1")
)
opt  <- parse_args(OptionParser(option_list=option_list))
load(opt$countsrdata)   # loads 'dds'
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dds <- DESeq(dds)

contrasts <- fread(opt$contrasts, sep=",", comment.char="#",
                   col.names=c("contrast_id","numerator","denominator"))
padj_thr <- as.numeric(opt$padj)
lfc_thr  <- as.numeric(opt$lfc)

for (i in seq_len(nrow(contrasts))) {
  cid <- contrasts$contrast_id[i]
  num <- contrasts$numerator[i]
  den <- contrasts$denominator[i]
  message("[deseq2_de.R] contrast: ", cid, "  (", num, " vs ", den, ")")
  coef <- paste0("condition_", num, "_vs_", den)

  res <- tryCatch(
    lfcShrink(dds, coef=coef, type="apeglm"),
    error=function(e) {
      message("apeglm failed (", conditionMessage(e), "), using type='normal'")
      lfcShrink(dds, contrast=c("condition",num,den), type="normal")
    }
  )

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  fwrite(res_df[order(res_df$padj, na.last=TRUE),],
         file.path(opt$outdir, paste0(cid, "_DE_results.tsv")), sep="\t")

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_thr &
                  abs(res_df$log2FoldChange) > lfc_thr, ]
  fwrite(sig, file.path(opt$outdir, paste0(cid, "_significant.tsv")), sep="\t")

  p <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(color=(!is.na(padj) & padj<padj_thr & abs(log2FoldChange)>lfc_thr)),
               alpha=.4, size=.7) +
    scale_color_manual(values=c("grey60","red2")) +
    labs(title=cid, color="Significant") + theme_bw()
  ggsave(file.path(opt$outdir, paste0(cid,"_volcano.pdf")), p, width=6, height=5)
}
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,"deseq2_de_sessionInfo.txt"))
message("[deseq2_de.R] Done.")
