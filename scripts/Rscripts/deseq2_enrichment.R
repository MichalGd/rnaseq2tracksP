# =============================================================================
# deseq2_enrichment.R - Gene enrichment analysis (ORA + GSEA) per contrast
# ORIGIN: NEW v4.3
#
# Input : _DE_results.tsv files written by deseq2_de.R (shrunken LFC, padj)
# Output: per contrast - per database:
#           TSV tables  : ORA and GSEA results
#           PDF plots   : dotplot, barplot, cnetplot (ORA); dotplot, barplot (GSEA)
#
# Databases covered (all offline, no internet required):
#   GO BP / MF / CC  - via clusterProfiler + org.Mm.eg.db / org.Hs.eg.db
#   KEGG             - via clusterProfiler enrichKEGG / gseKEGG
#   Reactome         - via ReactomePA
#   MSigDB Hallmarks - via fgsea + msigdbr
#
# Usage:
#   Rscript deseq2_enrichment.R \
#     --dedir     <path/to/DE/>           \
#     --contrasts <contrasts.csv>         \
#     --outdir    <path/to/enrichment/>   \
#     --species   mouse|human             \
#     --padj      0.05                    \
#     --lfc       1                       \
#     --minGS     10                      \
#     --maxGS     500
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ReactomePA)
  library(fgsea)
  library(msigdbr)
})

option_list <- list(
  make_option("--dedir",     type="character"),
  make_option("--contrasts", type="character"),
  make_option("--outdir",    type="character"),
  make_option("--species",   type="character", default="mouse"),
  make_option("--padj",      type="double",    default=0.05),
  make_option("--lfc",       type="double",    default=0),
  make_option("--minGS",     type="integer",   default=10L),
  make_option("--maxGS",     type="integer",   default=500L)
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# -- Species-specific settings -------------------------------------------------
if (opt$species == "human") {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  orgdb       <- org.Hs.eg.db
  kegg_org    <- "hsa"
  reactome_org <- "human"
  msig_species <- "Homo sapiens"
} else {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  orgdb       <- org.Mm.eg.db
  kegg_org    <- "mmu"
  reactome_org <- "mouse"
  msig_species <- "Mus musculus"
}

log <- function(...) message(sprintf("[deseq2_enrichment.R] %s", paste0(...)))

# -- ID conversion helper ------------------------------------------------------
# Strips Ensembl version suffixes and converts to Entrez IDs
to_entrez <- function(ensembl_ids) {
  ids_clean <- sub("[.][0-9]+$", "", ensembl_ids)
  map <- bitr(ids_clean, fromType="ENSEMBL", toType="ENTREZID",
              OrgDb=orgdb, drop=TRUE)
  map
}

# -- Plot helpers --------------------------------------------------------------
save_plot <- function(p, path, w=8, h=6) {
  tryCatch({
    ggsave(path, p, width=w, height=h, limitsize=FALSE)
    ggsave(sub("[.]pdf$", ".png", path), p, width=w, height=h, dpi=300, limitsize=FALSE)},
    error=function(e) message("  [WARN] plot failed: ", basename(path), " - ", conditionMessage(e))
  )
}

plot_ora <- function(res, title, outpfx, n=20) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible(NULL))
  tryCatch({
    save_plot(dotplot(res, showCategory=n, title=title),
              paste0(outpfx, "_dotplot.pdf"))
  }, error=function(e) message("  [WARN] dotplot: ", conditionMessage(e)))
  tryCatch({
    df <- as.data.frame(res)[seq_len(min(n, nrow(as.data.frame(res)))), ]
    df$Description <- factor(df$Description, levels=rev(df$Description))
    p <- ggplot(df, aes(Count, Description, fill=p.adjust)) +
      geom_col() +
      scale_fill_gradient(low="#CC0000", high="#AAAAAA", name="adj. p-value") +
      labs(title=title, x="Gene Count", y=NULL) +
      theme_bw()
    save_plot(p, paste0(outpfx, "_barplot.pdf"), h=max(4, nrow(df)*0.35+2))
  }, error=function(e) message("  [WARN] barplot: ", conditionMessage(e)))
  tryCatch({
    save_plot(cnetplot(res, showCategory=min(5,nrow(as.data.frame(res))),
                       node_label="gene"),
              paste0(outpfx, "_cnetplot.pdf"), w=10, h=8)
  }, error=function(e) message("  [WARN] cnetplot: ", conditionMessage(e)))
}

plot_gsea <- function(res, title, outpfx, n=20) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible(NULL))
  tryCatch({
    save_plot(dotplot(res, showCategory=n, title=title, split=".sign") +
                facet_grid(.~.sign),
              paste0(outpfx, "_dotplot.pdf"), w=10, h=6)
  }, error=function(e) message("  [WARN] GSEA dotplot: ", conditionMessage(e)))
  tryCatch({
    df <- as.data.frame(res)
    df <- df[order(df$NES), ][seq_len(min(n, nrow(df))), ]
    df$Description <- factor(df$Description, levels=df$Description)
    p <- ggplot(df, aes(NES, Description,
                        fill=ifelse(NES > 0, "up", "down"))) +
      geom_col() +
      scale_fill_manual(values=c("up"="#CC0000","down"="#2166AC"), name="Direction") +
      labs(title=title, x="NES", y=NULL) +
      theme_bw()
    save_plot(p, paste0(outpfx, "_barplot.pdf"), h=max(4, nrow(df)*0.35+2))
  }, error=function(e) message("  [WARN] GSEA barplot: ", conditionMessage(e)))
}

# -- MSigDB gene sets (loaded once) -------------------------------------------
log("Loading MSigDB Hallmarks for ", msig_species, "...")
msig_raw <- msigdbr(species=msig_species, collection="H")
# Column name changed across msigdbr versions
entrez_col <- if ("entrez_id" %in% names(msig_raw)) "entrez_id" else "entrez_gene"
msig_H <- split(as.character(msig_raw[[entrez_col]]), msig_raw$gs_name)

# -- Contrasts -----------------------------------------------------------------
contrasts <- fread(opt$contrasts, sep=",",
                   col.names=c("contrast_id","numerator","denominator"))

# -- Per-contrast loop ---------------------------------------------------------
for (i in seq_len(nrow(contrasts))) {
  cid <- contrasts$contrast_id[i]
  log("Processing contrast: ", cid)

  de_file <- file.path(opt$dedir, paste0(cid, "_DE_results.tsv"))
  if (!file.exists(de_file)) {
    log("  [SKIP] DE results not found: ", de_file); next
  }

  de <- fread(de_file)
  de[, gene_clean := sub("[.][0-9]+$", "", gene)]

  # -- ID conversion ------------------------------------------------------------
  id_map <- tryCatch(
    to_entrez(de$gene_clean),
    error=function(e) { log("  [WARN] ID conversion failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(id_map) || nrow(id_map) == 0) {
    log("  [SKIP] No Entrez IDs mapped for ", cid); next
  }
  de_mapped <- merge(de, id_map, by.x="gene_clean", by.y="ENSEMBL", all.x=FALSE)
  # Deduplicate: one Entrez ID per Ensembl gene (keep first match)
  de_mapped <- de_mapped[!duplicated(de_mapped$gene_clean), ]
  log(sprintf("  Mapped %d / %d genes to Entrez IDs (deduplicated)", nrow(de_mapped), nrow(de)))

  # -- Gene lists ---------------------------------------------------------------
  # ORA: significant genes (padj < threshold AND |LFC| > threshold)
  sig_genes <- de_mapped[!is.na(padj) & padj < opt$padj &
                           abs(log2FoldChange) > opt$lfc, ENTREZID]
  log(sprintf("  Sig genes for ORA: %d", length(sig_genes)))

  # Background: all genes with valid Entrez mapping and non-NA padj
  bg_genes <- de_mapped[!is.na(padj), ENTREZID]

  # GSEA: ranked list - stat or LFC - -log10(padj), sorted descending
  # Use log2FoldChange * -log10(padj+1e-300) as ranking metric
  gsea_df <- de_mapped[!is.na(log2FoldChange) & !is.na(padj)]
  # Prefer Wald stat if available (better GSEA ranking); fallback to LFC*-log10(padj)
  # Ranking metric: sign(LFC) * -log10(padj) - robust to ashr shrinkage, no tied pile-up
  gsea_df[, rank_metric := sign(log2FoldChange) * -log10(padj + 1e-300)]
  gsea_df <- gsea_df[order(-rank_metric)]
  # Deduplicate: keep highest |rank_metric| per Entrez ID
  gsea_df <- gsea_df[, .SD[which.max(abs(rank_metric))], by=ENTREZID]
  ranked_list <- sort(setNames(gsea_df$rank_metric, gsea_df$ENTREZID), decreasing=TRUE)

  # -- Output directory per contrast --------------------------------------------
  cdir <- file.path(opt$outdir, cid)
  dir.create(cdir, recursive=TRUE, showWarnings=FALSE)

  # ============================================================================
  # ORA
  # ============================================================================
  if (length(sig_genes) >= 5) {

    # -- GO BP ------------------------------------------------------------------
    log("  ORA: GO BP")
    ora_gobp <- tryCatch(
      enrichGO(gene=sig_genes, universe=bg_genes, OrgDb=orgdb,
               ont="BP", keyType="ENTREZID",
               minGSSize=opt$minGS, maxGSSize=opt$maxGS,
               pAdjustMethod="BH", pvalueCutoff=0.25, qvalueCutoff=0.2,
               readable=TRUE),
      error=function(e) { message("    GO BP ORA error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(ora_gobp)) {
      fwrite(as.data.frame(ora_gobp),
             file.path(cdir, paste0(cid, "_ORA_GOBP.tsv")), sep="\t")
      plot_ora(ora_gobp, paste0(cid, " ORA - GO:BP"),
               file.path(cdir, paste0(cid, "_ORA_GOBP")))
    }

    # -- GO MF ------------------------------------------------------------------
    log("  ORA: GO MF")
    ora_gomf <- tryCatch(
      enrichGO(gene=sig_genes, universe=bg_genes, OrgDb=orgdb,
               ont="MF", keyType="ENTREZID",
               minGSSize=opt$minGS, maxGSSize=opt$maxGS,
               pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
               readable=TRUE),
      error=function(e) { message("    GO MF ORA error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(ora_gomf)) {
      fwrite(as.data.frame(ora_gomf),
             file.path(cdir, paste0(cid, "_ORA_GOMF.tsv")), sep="\t")
      plot_ora(ora_gomf, paste0(cid, " ORA - GO:MF"),
               file.path(cdir, paste0(cid, "_ORA_GOMF")))
    }

    # -- GO CC ------------------------------------------------------------------
    log("  ORA: GO CC")
    ora_gocc <- tryCatch(
      enrichGO(gene=sig_genes, universe=bg_genes, OrgDb=orgdb,
               ont="CC", keyType="ENTREZID",
               minGSSize=opt$minGS, maxGSSize=opt$maxGS,
               pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
               readable=TRUE),
      error=function(e) { message("    GO CC ORA error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(ora_gocc)) {
      fwrite(as.data.frame(ora_gocc),
             file.path(cdir, paste0(cid, "_ORA_GOCC.tsv")), sep="\t")
      plot_ora(ora_gocc, paste0(cid, " ORA - GO:CC"),
               file.path(cdir, paste0(cid, "_ORA_GOCC")))
    }

    # -- KEGG -------------------------------------------------------------------
    log("  ORA: KEGG")
    ora_kegg <- tryCatch(
      enrichKEGG(gene=sig_genes, universe=bg_genes,
                 organism=kegg_org,
                 minGSSize=opt$minGS, maxGSSize=opt$maxGS,
                 pAdjustMethod="BH", pvalueCutoff=0.05,
                 ),
      error=function(e) { message("    KEGG ORA error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(ora_kegg)) {
      # Convert IDs to gene symbols for readability
      ora_kegg <- tryCatch(setReadable(ora_kegg, OrgDb=orgdb, keyType="ENTREZID"),
                           error=function(e) ora_kegg)
      fwrite(as.data.frame(ora_kegg),
             file.path(cdir, paste0(cid, "_ORA_KEGG.tsv")), sep="\t")
      plot_ora(ora_kegg, paste0(cid, " ORA - KEGG"),
               file.path(cdir, paste0(cid, "_ORA_KEGG")))
    }

    # -- Reactome ---------------------------------------------------------------
    log("  ORA: Reactome")
    ora_reactome <- tryCatch(
      enrichPathway(gene=sig_genes, universe=bg_genes,
                    organism=reactome_org,
                    minGSSize=opt$minGS, maxGSSize=opt$maxGS,
                    pAdjustMethod="BH", pvalueCutoff=0.05,
                    readable=TRUE),
      error=function(e) { message("    Reactome ORA error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(ora_reactome)) {
      fwrite(as.data.frame(ora_reactome),
             file.path(cdir, paste0(cid, "_ORA_Reactome.tsv")), sep="\t")
      plot_ora(ora_reactome, paste0(cid, " ORA - Reactome"),
               file.path(cdir, paste0(cid, "_ORA_Reactome")))
    }

  } else {
    log("  [SKIP] ORA - fewer than 5 sig genes")
  }

  # ============================================================================
  # GSEA
  # ============================================================================
  if (length(ranked_list) < 10) {
    log("  [SKIP] GSEA - fewer than 10 ranked genes"); next
  }

  # -- GSEA GO BP --------------------------------------------------------------
  log("  GSEA: GO BP")
  gsea_gobp <- tryCatch(
    gseGO(geneList=ranked_list, OrgDb=orgdb, ont="BP", keyType="ENTREZID",
          minGSSize=opt$minGS, maxGSSize=opt$maxGS,
          pAdjustMethod="BH", pvalueCutoff=0.25,
          nPermSimple=10000, verbose=FALSE),
    error=function(e) { message("    GSEA GO BP error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(gsea_gobp)) {
    gsea_gobp <- tryCatch(setReadable(gsea_gobp, OrgDb=orgdb, keyType="ENTREZID"),
                          error=function(e) gsea_gobp)
    fwrite(as.data.frame(gsea_gobp),
           file.path(cdir, paste0(cid, "_GSEA_GOBP.tsv")), sep="\t")
    plot_gsea(gsea_gobp, paste0(cid, " GSEA - GO:BP"),
              file.path(cdir, paste0(cid, "_GSEA_GOBP")))
  }

  # -- GSEA GO MF --------------------------------------------------------------
  log("  GSEA: GO MF")
  gsea_gomf <- tryCatch(
    gseGO(geneList=ranked_list, OrgDb=orgdb, ont="MF", keyType="ENTREZID",
          minGSSize=opt$minGS, maxGSSize=opt$maxGS,
          pAdjustMethod="BH", pvalueCutoff=0.25,
          nPermSimple=10000, verbose=FALSE),
    error=function(e) { message("    GSEA GO MF error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(gsea_gomf)) {
    gsea_gomf <- tryCatch(setReadable(gsea_gomf, OrgDb=orgdb, keyType="ENTREZID"),
                          error=function(e) gsea_gomf)
    fwrite(as.data.frame(gsea_gomf),
           file.path(cdir, paste0(cid, "_GSEA_GOMF.tsv")), sep="\t")
    plot_gsea(gsea_gomf, paste0(cid, " GSEA - GO:MF"),
              file.path(cdir, paste0(cid, "_GSEA_GOMF")))
  }

  # -- GSEA KEGG ---------------------------------------------------------------
  log("  GSEA: KEGG")
  gsea_kegg <- tryCatch(
    gseKEGG(geneList=ranked_list, organism=kegg_org,
            minGSSize=opt$minGS, maxGSSize=opt$maxGS,
            pAdjustMethod="BH", pvalueCutoff=0.25,
            verbose=FALSE),
    error=function(e) { message("    GSEA KEGG error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(gsea_kegg)) {
    gsea_kegg <- tryCatch(setReadable(gsea_kegg, OrgDb=orgdb, keyType="ENTREZID"),
                          error=function(e) gsea_kegg)
    fwrite(as.data.frame(gsea_kegg),
           file.path(cdir, paste0(cid, "_GSEA_KEGG.tsv")), sep="\t")
    plot_gsea(gsea_kegg, paste0(cid, " GSEA - KEGG"),
              file.path(cdir, paste0(cid, "_GSEA_KEGG")))
  }

  # -- GSEA Reactome ------------------------------------------------------------
  log("  GSEA: Reactome")
  gsea_reactome <- tryCatch(
    gsePathway(geneList=ranked_list, organism=reactome_org,
               minGSSize=opt$minGS, maxGSSize=opt$maxGS,
               pAdjustMethod="BH", pvalueCutoff=0.25, verbose=FALSE),
    error=function(e) { message("    GSEA Reactome error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(gsea_reactome)) {
    gsea_reactome <- tryCatch(setReadable(gsea_reactome, OrgDb=orgdb, keyType="ENTREZID"),
                              error=function(e) gsea_reactome)
    fwrite(as.data.frame(gsea_reactome),
           file.path(cdir, paste0(cid, "_GSEA_Reactome.tsv")), sep="\t")
    plot_gsea(gsea_reactome, paste0(cid, " GSEA - Reactome"),
              file.path(cdir, paste0(cid, "_GSEA_Reactome")))
  }

  # -- GSEA MSigDB Hallmarks (fgsea) --------------------------------------------
  log("  GSEA: MSigDB Hallmarks (fgsea)")
  fgsea_H <- tryCatch(
    fgsea(pathways=msig_H, stats=ranked_list,
          minSize=opt$minGS, maxSize=opt$maxGS,
          nPermSimple=10000),
    error=function(e) { message("    fgsea Hallmarks error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(fgsea_H) && nrow(fgsea_H) > 0) {
    fgsea_H <- fgsea_H[order(padj)]
    # Convert leadingEdge list to string for TSV output
    fgsea_out <- copy(fgsea_H)
    fgsea_out[, leadingEdge := sapply(leadingEdge, paste, collapse=",")]
    fwrite(fgsea_out,
           file.path(cdir, paste0(cid, "_GSEA_Hallmarks.tsv")), sep="\t")
    # Plot significant Hallmarks
    sig_H <- fgsea_H[padj < 0.05]
    if (nrow(sig_H) > 0) {
      sig_H[, pathway_short := sub("HALLMARK_", "", pathway)]
      p_H <- ggplot(sig_H[order(NES)],
                    aes(x=NES, y=reorder(pathway_short, NES),
                        fill=ifelse(NES > 0, "up", "down"))) +
        geom_col() +
        scale_fill_manual(values=c("up"="#CC0000","down"="#2166AC"),
                          name="Direction") +
        labs(title=paste0(cid, " GSEA - MSigDB Hallmarks (padj<0.05)"),
             x="Normalized Enrichment Score", y=NULL) +
        theme_bw() + theme(legend.position="right")
      save_plot(p_H, file.path(cdir, paste0(cid, "_GSEA_Hallmarks_barplot.pdf")),
                w=10, h=max(4, nrow(sig_H)*0.35 + 2))
    }
  }

  log(sprintf("  Done: %s - %s/", cid, cdir))
}

writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir, "deseq2_enrichment_sessionInfo.txt"))
message("[deseq2_enrichment.R] All contrasts complete.")
