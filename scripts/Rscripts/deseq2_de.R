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

contrasts <- fread(opt$contrasts, sep=",",
                   col.names=c("contrast_id","numerator","denominator"))
padj_thr <- as.numeric(opt$padj)   # 0.05 — hard significant
lfc_thr  <- as.numeric(opt$lfc)

# Colour scheme shared by volcano and MA plots
# RED        : padj < 0.05  (significant)
# LIGHT GREEN: 0.05 <= padj < 0.10  (trend)
# GREY       : padj >= 0.10 or NA  (not significant)
sig_colors  <- c("ns"="#AAAAAA", "trend"="#90EE90", "sig"="#CC0000")
sig_labels  <- c("ns"="padj ≥ 0.10", "trend"="0.05 ≤ padj < 0.10", "sig"="padj < 0.05")

classify_sig <- function(padj_vec) {
  tier <- rep("ns", length(padj_vec))
  tier[!is.na(padj_vec) & padj_vec <  0.10] <- "trend"
  tier[!is.na(padj_vec) & padj_vec <  0.05] <- "sig"
  factor(tier, levels=c("ns","trend","sig"))
}

for (i in seq_len(nrow(contrasts))) {
  cid <- contrasts$contrast_id[i]
  num <- contrasts$numerator[i]
  den <- contrasts$denominator[i]
  message("[deseq2_de.R] contrast: ", cid, "  (", num, " vs ", den, ")")
  coef <- paste0("condition_", num, "_vs_", den)

  res <- tryCatch(
    lfcShrink(dds, coef=coef, type="apeglm"),
    error=function(e) {
      message("apeglm failed (", conditionMessage(e), "), using type='ashr'")
      lfcShrink(dds, contrast=c("condition",num,den), type="ashr")
    }
  )

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$sig_tier <- classify_sig(res_df$padj)

  fwrite(res_df[order(res_df$padj, na.last=TRUE),],
         file.path(opt$outdir, paste0(cid, "_DE_results.tsv")), sep="\t")

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_thr &
                  abs(res_df$log2FoldChange) > lfc_thr, ]
  fwrite(sig, file.path(opt$outdir, paste0(cid, "_significant.tsv")), sep="\t")

  # ── Volcano plot ────────────────────────────────────────────────────────────
  # Plot ns first so sig points are drawn on top
  res_plot <- res_df[order(res_df$sig_tier), ]
  p_vol <- ggplot(res_plot, aes(log2FoldChange, -log10(pvalue), color=sig_tier)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=sig_colors, labels=sig_labels, name="Significance") +
    geom_vline(xintercept=c(-lfc_thr, lfc_thr), linetype="dashed",
               color="grey40", linewidth=0.4) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               color="grey40", linewidth=0.4) +
    labs(title=cid, x="log2 Fold Change", y="-log10(p-value)") +
    theme_bw() +
    theme(legend.position="right")
  ggsave(file.path(opt$outdir, paste0(cid, "_volcano.pdf")), p_vol, width=7, height=5)

  # ── MA plot ─────────────────────────────────────────────────────────────────
  # baseMean on log10 x-axis; LFC on y-axis; same colour tiers
  res_ma <- res_df[res_df$baseMean > 0, ]
  res_ma <- res_ma[order(res_ma$sig_tier), ]
  p_ma <- ggplot(res_ma, aes(log10(baseMean), log2FoldChange, color=sig_tier)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=sig_colors, labels=sig_labels, name="Significance") +
    geom_hline(yintercept=c(-lfc_thr, lfc_thr), linetype="dashed",
               color="grey40", linewidth=0.4) +
    geom_hline(yintercept=0, color="black", linewidth=0.4) +
    labs(title=cid, x="log10(mean expression)", y="log2 Fold Change") +
    theme_bw() +
    theme(legend.position="right")
  ggsave(file.path(opt$outdir, paste0(cid, "_MA.pdf")), p_ma, width=7, height=5)
}

writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir,"deseq2_de_sessionInfo.txt"))

# =============================================================================
# Annotated count tables — one CSV per contrast
# Merges: GTF annotation + raw counts + normalized counts + DE results
# =============================================================================
message("[deseq2_de.R] Building annotated count tables...")

counts_dir <- dirname(opt$countsrdata)   # analysis/counts/
gtf_path   <- Sys.getenv("GTF")          # passed via environment from rnaseq2tracks.sh
tables_dir <- file.path(dirname(opt$outdir), "tables")
dir.create(tables_dir, recursive=TRUE, showWarnings=FALSE)

if (gtf_path == "" || !file.exists(gtf_path)) {
  message("[deseq2_de.R] [WARN] GTF not found via $GTF env var — skipping annotated tables")
} else {
  suppressPackageStartupMessages(library(data.table))

  # Gene annotation from GTF
  gtf_dt <- fread(cmd=paste("grep -v '^#'", gtf_path),
                  sep="\t", header=FALSE, quote="",
                  col.names=c("chr","source","feature","start","end",
                              "score","strand","frame","attributes"))
  gtf_dt <- gtf_dt[feature == "gene"]
  extract_attr <- function(attrs, key) {
    m <- regmatches(attrs, regexpr(paste0(key, ' "[^"]+'), attrs))
    sub(paste0(key, ' "'), "", m)
  }
  genes_ann <- data.table(
    gene_id      = sub("\\.[0-9]+$", "", extract_attr(gtf_dt$attributes, "gene_id")),
    gene_id_full = extract_attr(gtf_dt$attributes, "gene_id"),
    gene_name    = ifelse(grepl('gene_name "', gtf_dt$attributes),
                          extract_attr(gtf_dt$attributes, "gene_name"), NA_character_),
    gene_type    = ifelse(grepl('gene_type "', gtf_dt$attributes),
                          extract_attr(gtf_dt$attributes, "gene_type"), NA_character_),
    chr          = gtf_dt$chr,
    start        = gtf_dt$start,
    end          = gtf_dt$end,
    strand       = gtf_dt$strand
  )
  genes_ann[, protein_coding := gene_type == "protein_coding"]

  # Raw counts
  raw <- fread(file.path(counts_dir, "raw_counts.tsv"))
  setnames(raw, "gene", "gene_id_full")
  raw[, gene_id := sub("\\.[0-9]+$", "", gene_id_full)]
  raw_scols <- setdiff(names(raw), c("gene_id_full","gene_id"))
  setnames(raw, raw_scols, paste0("raw_", raw_scols))

  # Normalized counts
  norm <- fread(file.path(counts_dir, "normalized_counts.tsv"))
  setnames(norm, "gene", "gene_id_full")
  norm[, gene_id := sub("\\.[0-9]+$", "", gene_id_full)]
  norm_scols <- setdiff(names(norm), c("gene_id_full","gene_id"))
  setnames(norm, norm_scols, paste0("norm_", norm_scols))

  # Base table
  base_tbl <- merge(genes_ann,
                    raw[, c("gene_id", paste0("raw_", raw_scols)), with=FALSE],
                    by="gene_id", all.x=TRUE)
  base_tbl <- merge(base_tbl,
                    norm[, c("gene_id", paste0("norm_", norm_scols)), with=FALSE],
                    by="gene_id", all.x=TRUE)

  # One table per contrast
  for (i in seq_len(nrow(contrasts))) {
    cid     <- contrasts$contrast_id[i]
    de_file <- file.path(opt$outdir, paste0(cid, "_DE_results.tsv"))
    if (!file.exists(de_file)) next
    de <- fread(de_file)
    de[, gene_id := sub("\\.[0-9]+$", "", gene)]
    de_cols <- intersect(c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj"), names(de))
    de_sub  <- de[, de_cols, with=FALSE]
    setnames(de_sub, setdiff(de_cols,"gene_id"),
             paste0(setdiff(de_cols,"gene_id"), "_", cid))
    tbl <- merge(base_tbl, de_sub, by="gene_id", all.x=TRUE)
    padj_col <- paste0("padj_", cid)
    base_col <- paste0("baseMean_", cid)
    tbl <- tbl[order(get(padj_col), -get(base_col), na.last=TRUE)]
    outfile <- file.path(tables_dir, paste0(cid, "_annotated_counts.csv"))
    fwrite(tbl, outfile, sep=",", na="NA")
    message(sprintf("[deseq2_de.R] Table written: %s  (%d genes x %d cols)",
                    outfile, nrow(tbl), ncol(tbl)))
  }
}

message("[deseq2_de.R] Done.")
