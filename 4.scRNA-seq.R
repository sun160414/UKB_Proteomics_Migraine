## =========================================================
## scRNA-seq pipeline (Seurat + Harmony)
## - robust input (multiple .h5)
## - QC + SCT + Harmony + UMAP
## - annotation plots
## - composition plots (bar/alluvial/donut)
## - DE: Migraine vs HC (all cells + per celltype)
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
  library(harmony)
})

## -------------------------
## 0) Palette
## -------------------------
custom_palette <- function() {
  c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896",
    "#9467BD", "#C5B0D5", "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7",
    "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#7FC97F", "#BEAED4", "#FDC086", "#386CB0",
    "#F0027F", "#BF5B16")
}

## -------------------------
## 1) User params (EDIT HERE)
## -------------------------
data_dir <- "/home/ug1268u4/UKB/Project/Proteomics_Migraine/3.scRNA-seq/GSE269117"
out_dir  <- "/home/ug1268u4/UKB/Project/Proteomics_Migraine/3.scRNA-seq/_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define case/control by sample IDs (orig.ident)
case_samples <- c("GSM8306599", "GSM8306600", "GSM8306606", "GSM8306608", "GSM8306615")
case_label <- "Migraine"
ctrl_label <- "HC"

# QC thresholds (generic; tune per dataset)
qc <- list(
  min_features = 600,
  max_features = 13500,
  min_counts   = 500,
  max_counts   = 18000,
  max_mt       = 25,
  max_hb       = 10
)

# Harmony + clustering
pcs_use <- 1:30
resolution_use <- 0.4

# DE settings
de_min_pct <- 0.1
de_logfc_threshold <- 0
min_cells_celltype <- 50
min_cells_per_group <- 20

## -------------------------
## 2) Read multiple 10X .h5 and merge
## -------------------------
setwd(data_dir)
h5_files <- list.files(pattern = "\\.h5$")
stopifnot(length(h5_files) > 0)

sce_list <- lapply(h5_files, function(f){
  x <- Read10X_h5(f)
  sample_id <- str_split(f, "_", simplify = TRUE)[1,1]
  CreateSeuratObject(x, project = sample_id, min.features = 200, min.cells = 3)
})

sample_ids <- substr(h5_files, 1, 10)  # fallback add.cell.ids
pbmc <- merge(sce_list[[1]], y = sce_list[-1], add.cell.ids = sample_ids)

## If you have Seurat v5 layers, join them
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  pbmc <- JoinLayers(pbmc)
}

## -------------------------
## 3) QC metrics + filtering
## -------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

hb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
hb_genes <- CaseMatch(hb_genes, rownames(pbmc))
pbmc[["percent.HB"]] <- PercentageFeatureSet(pbmc, features = hb_genes)

# QC violin before
qc_feats <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")
theme_blankx <- theme(axis.title.x = element_blank())

p_before <- wrap_plots(lapply(qc_feats, function(f){
  VlnPlot(pbmc, group.by = "orig.ident", features = f, pt.size = 0) +
    theme_blankx + NoLegend() +
    scale_fill_manual(values = custom_palette())
}), nrow = 2)

ggsave(file.path(out_dir, "QC_violin_before.pdf"), p_before, width = 14, height = 8)

# Filter
pbmc <- subset(
  pbmc,
  subset =
    nFeature_RNA > qc$min_features & nFeature_RNA < qc$max_features &
    nCount_RNA   > qc$min_counts   & nCount_RNA   < qc$max_counts &
    percent.mt   < qc$max_mt &
    percent.HB   < qc$max_hb
)

p_after <- wrap_plots(lapply(qc_feats, function(f){
  VlnPlot(pbmc, group.by = "orig.ident", features = f, pt.size = 0) +
    theme_blankx + NoLegend() +
    scale_fill_manual(values = custom_palette())
}), nrow = 2)

ggsave(file.path(out_dir, "QC_violin_after.pdf"), p_after, width = 14, height = 8)

saveRDS(pbmc, file.path(out_dir, "pbmc_after_qc.rds"))

## -------------------------
## 4) Normalization (SCT) + PCA + Harmony + UMAP
## -------------------------
# For large data: avoid parallel memory blow-ups
if (requireNamespace("future", quietly = TRUE)) {
  library(future)
  plan("sequential")
  options(future.globals.maxSize = 5e9)
}

pbmc <- SCTransform(pbmc, verbose = FALSE)
DefaultAssay(pbmc) <- "SCT"

pbmc <- RunPCA(pbmc, verbose = FALSE)

# Harmony (batch = orig.ident)
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident", assay.use = "SCT", max.iter.harmony = 20)

pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs_use)
pbmc <- FindClusters(pbmc, resolution = resolution_use)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs_use)

saveRDS(pbmc, file.path(out_dir, "pbmc_harmony_umap.rds"))

## -------------------------
## 5) UMAP (square + optional legend)
## -------------------------
pal <- custom_palette()
p_umap <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", raster = TRUE, pt.size = 0.5) +
  scale_color_manual(values = pal) +
  coord_fixed() +
  theme(aspect.ratio = 1)

ggsave(file.path(out_dir, "UMAP_orig_no_legend.pdf"),
       p_umap + theme(legend.position = "none"), width = 4, height = 4)

ggsave(file.path(out_dir, "UMAP_orig_with_legend.pdf"),
       p_umap + theme(legend.position = "right"), width = 6, height = 4)

## -------------------------
## 6) Group label (case/control) from sample IDs
## -------------------------
pbmc$group <- ifelse(pbmc$orig.ident %in% case_samples, case_label, ctrl_label)
pbmc$group <- factor(pbmc$group, levels = c(ctrl_label, case_label))
write.csv(table(pbmc$orig.ident, pbmc$group), file.path(out_dir, "group_table.csv"))

## -------------------------
## 7) Optional: markers & DotPlot / FeaturePlot
## -------------------------
# Use your existing celltype column if present; otherwise use clusters
celltype_col <- if ("celltype.main" %in% colnames(pbmc@meta.data)) "celltype.main" else "seurat_clusters"

markers <- unique(c(
  "CD3D","CD3E","CD3G","CD4","CCR7","CD8A","CD8B","GZMK",
  "NKG7","GNLY","PRF1","KLRD1","KLRC1",
  "MS4A1","CD79A","CD79B","CD74","HLA-DRA",
  "LYZ","S100A8","S100A9","LST1","CTSD","FCGR3A"
))

p_dot <- DotPlot(pbmc, features = markers, group.by = celltype_col) +
  RotatedAxis() +
  scale_color_gradientn(colours = c("#1F77B4", "#f0f0f0", "#cc0000"))

ggsave(file.path(out_dir, "DotPlot_celltype.pdf"), p_dot, width = 10, height = 4)

# FeaturePlots (safe gradient)
feat_plots <- lapply(markers, function(g){
  FeaturePlot(pbmc, features = g, raster = TRUE, pt.size = 0.4) +
    coord_fixed() + theme(aspect.ratio = 1) +
    ggtitle(g) +
    scale_color_gradientn(colours = c("#f0f0f0", "#a00000"))
})

ggsave(file.path(out_dir, "FeaturePlot_markers.pdf"),
       wrap_plots(feat_plots, ncol = 4), width = 16, height = 24)

## -------------------------
## 8) Cell composition (per sample + per group)
## -------------------------
meta <- pbmc@meta.data %>%
  mutate(celltype = as.character(.data[[celltype_col]])) %>%
  mutate(celltype = ifelse(is.na(celltype) | celltype == "", "Unknown", celltype))

# proportions per sample
prop_sample <- meta %>%
  count(orig.ident, celltype, name = "count") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# proportions per group
prop_group <- meta %>%
  count(group, celltype, name = "count") %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

write.csv(prop_sample, file.path(out_dir, "celltype_proportion_per_sample.csv"), row.names = FALSE)
write.csv(prop_group,  file.path(out_dir, "celltype_proportion_per_group.csv"),  row.names = FALSE)

# color map per celltype (stable)
cell_levels <- sort(unique(prop_group$celltype))
cell_cols <- setNames(custom_palette()[seq_along(cell_levels)], cell_levels)

# stacked bar (group)
p_bar_group <- ggplot(prop_group, aes(x = group, y = proportion, fill = celltype)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = cell_cols) +
  theme_bw() +
  labs(x = "Group", y = "Proportion", fill = NULL)

ggsave(file.path(out_dir, "stacked_bar_group.pdf"), p_bar_group, width = 3.5, height = 4)

# alluvial (group)
p_alluvial <- ggplot(prop_group,
                     aes(x = group, y = proportion, fill = celltype,
                         stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.6) +
  geom_flow(width = 0.6, alpha = 0.3, color = "white") +
  scale_fill_manual(values = cell_cols) +
  theme_classic() +
  labs(x = "Group", y = "Proportion")

ggsave(file.path(out_dir, "alluvial_group.pdf"), p_alluvial, width = 4, height = 4.5)

# donut (group)
total_cells <- prop_group %>% group_by(group) %>% summarise(total_cells = sum(count), .groups="drop")

p_donut <- ggplot(prop_group, aes(x = 3, y = proportion, fill = celltype)) +
  geom_col(width = 1.5, color = "white") +
  facet_grid(. ~ group) +
  coord_polar(theta = "y") +
  xlim(c(0.2, 3.8)) +
  scale_fill_manual(values = cell_cols) +
  theme_void() +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3) +
  geom_text(data = total_cells,
            aes(x = 3, y = 0, label = paste0("Total\n", total_cells)),
            inherit.aes = FALSE, size = 5, fontface = "bold")

ggsave(file.path(out_dir, "donut_group.pdf"), p_donut, width = 8, height = 4)

## -------------------------
## 9) Differential expression: Migraine vs HC
## -------------------------
Idents(pbmc) <- pbmc$group

de_all <- FindMarkers(
  pbmc,
  ident.1 = case_label,
  ident.2 = ctrl_label,
  test.use = "wilcox",
  min.pct = de_min_pct,
  logfc.threshold = de_logfc_threshold,
  only.pos = FALSE
) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

write.csv(de_all, file.path(out_dir, "DE_case_vs_ctrl_allcells.csv"), row.names = FALSE)

# DE per celltype
celltypes <- sort(unique(as.character(meta$celltype)))

de_list <- lapply(celltypes, function(ct){
  cells_ct <- rownames(pbmc@meta.data)[meta$celltype == ct]
  if (length(cells_ct) < min_cells_celltype) return(NULL)
  
  obj_ct <- subset(pbmc, cells = cells_ct)
  if (!all(c(ctrl_label, case_label) %in% unique(obj_ct$group))) return(NULL)
  if (min(table(obj_ct$group)) < min_cells_per_group) return(NULL)
  
  Idents(obj_ct) <- obj_ct$group
  de <- FindMarkers(
    obj_ct,
    ident.1 = case_label,
    ident.2 = ctrl_label,
    test.use = "wilcox",
    min.pct = de_min_pct,
    logfc.threshold = de_logfc_threshold,
    only.pos = FALSE
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(celltype = ct) %>%
    arrange(p_val_adj, desc(avg_log2FC))
  
  de
})

de_by_celltype <- bind_rows(de_list)
write.csv(de_by_celltype, file.path(out_dir, "DE_case_vs_ctrl_byCelltype.csv"), row.names = FALSE)

cat("DONE. Outputs in:", out_dir, "\n")