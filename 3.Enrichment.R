## =========================================================
## GO + KEGG enrichment
## =========================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(stringr)
  library(scales)     # for alpha()
})

## -------------------------
## 0) User settings
## -------------------------
infile <- "mygenes.csv"       # input gene list file
show_genes <- TRUE             # show geneID text under term label
go_top_per_ont <- 4            # top terms per ontology (BP/CC/MF)
kegg_top_n <- 12               # top KEGG pathways
padj_method <- "BH"            # FDR adjustment
p_cut <- 1                     # keep all, filter later if needed
q_cut <- 1

## -------------------------
## 1) Read genes (robust column detection)
## -------------------------
genes_df <- read.csv(infile, check.names = FALSE)

pick_gene_col <- function(df){
  cand <- c("gene","Gene","symbol","SYMBOL","genes","Genes")
  hit <- intersect(cand, names(df))
  if (length(hit) > 0) return(hit[1])
  return(names(df)[1])  # fallback: first column
}

gene_col <- pick_gene_col(genes_df)

input_genes <- genes_df[[gene_col]] %>%
  as.character() %>%
  str_trim() %>%
  .[. != ""] %>%
  unique()

cat("Input genes:", length(input_genes), "\n")

## -------------------------
## 2) ID mapping (SYMBOL -> ENTREZID); report failures
## -------------------------
ID <- bitr(input_genes,
           fromType = "SYMBOL",
           toType   = c("ENTREZID","SYMBOL"),
           OrgDb    = org.Hs.eg.db)

mapped <- unique(ID$SYMBOL)
failed <- setdiff(input_genes, mapped)

cat("Mapped:", length(mapped), " | Failed:", length(failed), "\n")
write.csv(data.frame(failed_gene = failed), "failed_symbol_mapping.csv", row.names = FALSE)

if (nrow(ID) == 0) stop("No genes mapped. Check gene symbols / species / input column.")

## -------------------------
## 3) Enrichment: GO + KEGG
## -------------------------
ego <- enrichGO(
  gene          = ID$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = padj_method,
  pvalueCutoff  = p_cut,
  qvalueCutoff  = q_cut,
  readable      = TRUE
)

ekegg <- enrichKEGG(
  gene          = ID$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = padj_method,
  pvalueCutoff  = p_cut,
  qvalueCutoff  = q_cut
)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

go_df   <- as.data.frame(ego)
kegg_df <- as.data.frame(ekegg)

write.csv(go_df,   "GO_enrichment_results.csv",   row.names = FALSE)
write.csv(kegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
