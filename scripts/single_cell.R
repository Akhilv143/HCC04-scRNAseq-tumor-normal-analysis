## ============================================================
## HCC04 scRNA-seq: Tumor vs Normal Analysis
## Dataset: GEO GSE149614
## Patient: HCC04
## Platform: 10x Genomics Chromium Single Cell 3'
## ============================================================

## ============================================================
## STEP 0: Setup folders and packages
## ============================================================

dir.create("figures", showWarnings = FALSE)
dir.create("figures/clustering", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/markers", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/composition", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/differential_expression", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/pathway_enrichment", recursive = TRUE, showWarnings = FALSE)

dir.create("tables", showWarnings = FALSE)
dir.create("tables/markers", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/differential_expression", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/pathway_enrichment", recursive = TRUE, showWarnings = FALSE)

dir.create("objects", showWarnings = FALSE)

library(data.table)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

set.seed(123)

## ============================================================
## STEP 1: Load GEO cell metadata
## ============================================================

meta <- read.table(
  "metadata/GSE149614_HCC.metadata.updated.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

cat("Metadata dimensions:", dim(meta), "\n")
print(colnames(meta))
print(table(meta$site))
print(table(meta$patient, meta$site))

## ============================================================
## STEP 2: Select matched HCC04 Normal and Tumor cells
## ============================================================

pilot_patient <- "HCC04"

wanted_cells <- rownames(meta)[
  meta$patient == pilot_patient &
    meta$site %in% c("Normal", "Tumor")
]

cat("Selected cells:", length(wanted_cells), "\n")
print(table(meta[wanted_cells, "site"]))

## ============================================================
## STEP 3: Memory-safe selective count-matrix loading
## ============================================================

## The full count matrix is not included in GitHub.
## Download it from GEO and store it locally at:
## data/GSE149614_HCC.scRNAseq.S71915.count.txt

select_cols <- c("V1", wanted_cells)

counts_pilot <- fread(
  "data/GSE149614_HCC.scRNAseq.S71915.count.txt",
  select = select_cols,
  showProgress = TRUE
)

gene_names <- counts_pilot$V1
counts_pilot[, V1 := NULL]

counts_dense <- as.matrix(counts_pilot)

rm(counts_pilot)
gc()

counts_sparse <- Matrix(
  counts_dense,
  sparse = TRUE
)

rownames(counts_sparse) <- gene_names
colnames(counts_sparse) <- wanted_cells

rm(counts_dense, gene_names)
gc()

cat("Count-matrix dimensions:", dim(counts_sparse), "\n")

## ============================================================
## STEP 4: Create Seurat object and attach metadata
## ============================================================

hcc04 <- CreateSeuratObject(
  counts = counts_sparse,
  project = "GSE149614_HCC04",
  min.cells = 3,
  min.features = 200
)

rm(counts_sparse)
gc()

pilot_meta <- meta[colnames(hcc04), , drop = FALSE]

hcc04 <- AddMetaData(
  object = hcc04,
  metadata = pilot_meta
)

cat("Seurat object summary:\n")
print(hcc04)

cat("Cells by site:\n")
print(table(hcc04$site))

cat("Cells by GEO-published cell type:\n")
print(table(hcc04$celltype))

## ============================================================
## STEP 5: Quality control metrics and visualizations
## ============================================================

hcc04[["percent.mt"]] <- PercentageFeatureSet(
  hcc04,
  pattern = "^MT-"
)

hcc04[["percent.ribo"]] <- PercentageFeatureSet(
  hcc04,
  pattern = "^RP[SL]"
)

cat("QC summary: nFeature_RNA\n")
print(summary(hcc04$nFeature_RNA))

cat("QC summary: nCount_RNA\n")
print(summary(hcc04$nCount_RNA))

cat("QC summary: percent.mt\n")
print(summary(hcc04$percent.mt))

p_qc_violin <- VlnPlot(
  hcc04,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "site",
  ncol = 3,
  pt.size = 0.05
)

ggsave(
  "figures/HCC04_QC_violin.png",
  p_qc_violin,
  width = 12,
  height = 5,
  dpi = 300
)

p_qc_counts_features <- FeatureScatter(
  hcc04,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "site"
)

ggsave(
  "figures/HCC04_QC_counts_vs_features.png",
  p_qc_counts_features,
  width = 7,
  height = 6,
  dpi = 300
)

p_qc_counts_mt <- FeatureScatter(
  hcc04,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "site"
)

ggsave(
  "figures/HCC04_QC_counts_vs_mitochondrial.png",
  p_qc_counts_mt,
  width = 7,
  height = 6,
  dpi = 300
)

## ============================================================
## STEP 6: QC filtering
## ============================================================

cells_before_qc <- ncol(hcc04)

hcc04 <- subset(
  hcc04,
  subset =
    nFeature_RNA >= 500 &
    nFeature_RNA <= 6500 &
    nCount_RNA <= 50000 &
    percent.mt < 15
)

cells_after_qc <- ncol(hcc04)

cat("Cells before QC:", cells_before_qc, "\n")
cat("Cells after QC:", cells_after_qc, "\n")
cat("Cells removed:", cells_before_qc - cells_after_qc, "\n")

print(table(hcc04$site))
print(table(hcc04$celltype))

## ============================================================
## STEP 7: Normalization, variable features, PCA
## ============================================================

hcc04 <- NormalizeData(
  hcc04,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

hcc04 <- FindVariableFeatures(
  hcc04,
  selection.method = "vst",
  nfeatures = 2000
)

top10_variable_genes <- head(VariableFeatures(hcc04), 10)

p_variable_features <- VariableFeaturePlot(hcc04)

p_variable_features <- LabelPoints(
  plot = p_variable_features,
  points = top10_variable_genes,
  repel = TRUE
)

ggsave(
  "figures/HCC04_variable_features.png",
  p_variable_features,
  width = 8,
  height = 6,
  dpi = 300
)

hcc04 <- ScaleData(
  hcc04,
  features = VariableFeatures(hcc04),
  vars.to.regress = "percent.mt"
)

hcc04 <- RunPCA(
  hcc04,
  features = VariableFeatures(hcc04),
  npcs = 30
)

p_elbow <- ElbowPlot(hcc04, ndims = 30)

ggsave(
  "figures/HCC04_elbow_plot.png",
  p_elbow,
  width = 7,
  height = 5,
  dpi = 300
)

## ============================================================
## STEP 8: Neighbours, Louvain clusters and UMAP
## ============================================================

n_pcs <- 15

hcc04 <- FindNeighbors(
  hcc04,
  dims = 1:n_pcs
)

hcc04 <- FindClusters(
  hcc04,
  resolution = 0.5
)

hcc04 <- RunUMAP(
  hcc04,
  dims = 1:n_pcs,
  seed.use = 123
)

cat("Cells per Seurat cluster:\n")
print(table(Idents(hcc04)))

p_clusters <- DimPlot(
  hcc04,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("HCC04: Seurat clusters")

p_site <- DimPlot(
  hcc04,
  reduction = "umap",
  group.by = "site"
) +
  ggtitle("HCC04: Normal versus Tumor")

p_celltype <- DimPlot(
  hcc04,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("HCC04: GEO-published cell types")

ggsave(
  "figures/clustering/HCC04_UMAP_clusters.png",
  p_clusters,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/clustering/HCC04_UMAP_tumor_normal.png",
  p_site,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/clustering/HCC04_UMAP_published_celltypes.png",
  p_celltype,
  width = 9,
  height = 7,
  dpi = 300
)

## ============================================================
## STEP 9: Cluster markers and annotation validation
## ============================================================

all_markers <- FindAllMarkers(
  hcc04,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  all_markers,
  "tables/markers/HCC04_all_Seurat_cluster_markers.csv",
  row.names = FALSE
)

top5_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  )

write.csv(
  top5_markers,
  "tables/markers/HCC04_top5_marker_genes_per_cluster.csv",
  row.names = FALSE
)

cluster_celltype_table <- table(
  Seurat_cluster = Idents(hcc04),
  Published_celltype = hcc04$celltype
)

write.csv(
  as.data.frame(cluster_celltype_table),
  "tables/markers/HCC04_cluster_vs_published_celltype.csv",
  row.names = FALSE
)

print(cluster_celltype_table)

heatmap_genes <- unique(top5_markers$gene)

p_heatmap <- DoHeatmap(
  hcc04,
  features = heatmap_genes,
  size = 3
) +
  NoLegend()

ggsave(
  "figures/markers/HCC04_cluster_marker_heatmap.png",
  p_heatmap,
  width = 12,
  height = 10,
  dpi = 300
)

canonical_markers <- list(
  Hepatocyte = c("ALB", "APOA1", "TTR", "CYP3A4", "KRT8", "KRT18"),
  T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"),
  Myeloid = c("LYZ", "TYROBP", "C1QA", "C1QB", "CD68", "FCER1G"),
  B_cell = c("MS4A1", "CD79A", "CD74", "HLA-DRA"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1")
)

canonical_markers <- lapply(
  canonical_markers,
  function(x) x[x %in% rownames(hcc04)]
)

p_dotplot <- DotPlot(
  hcc04,
  features = canonical_markers
) +
  RotatedAxis()

ggsave(
  "figures/markers/HCC04_canonical_marker_dotplot.png",
  p_dotplot,
  width = 13,
  height = 7,
  dpi = 300
)

## ============================================================
## STEP 10: Normal vs Tumor cell-type composition
## ============================================================

celltype_comp <- hcc04@meta.data %>%
  count(site, celltype, name = "n_cells") %>%
  group_by(site) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  celltype_comp,
  "tables/HCC04_celltype_composition.csv",
  row.names = FALSE
)

p_composition <- ggplot(
  celltype_comp,
  aes(x = site, y = proportion, fill = celltype)
) +
  geom_col(color = "white") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "HCC04 cell-type composition: Normal versus Tumor",
    x = NULL,
    y = "Cell proportion",
    fill = "Published cell type"
  ) +
  theme_classic()

ggsave(
  "figures/composition/HCC04_celltype_composition.png",
  p_composition,
  width = 8,
  height = 5,
  dpi = 300
)

## ============================================================
## STEP 11: Strict hepatocyte-like cell selection
## ============================================================

hepatocyte_strict <- subset(
  hcc04,
  subset =
    celltype == "Hepatocyte" &
    ALB > 0 &
    APOA1 > 0 &
    PTPRC == 0 &
    MS4A1 == 0 &
    CD79A == 0
)

cat("Strict hepatocyte-like cells by tissue:\n")
print(table(hepatocyte_strict$site))

saveRDS(
  hepatocyte_strict,
  "objects/HCC04_strict_hepatocyte_subset.rds"
)

## ============================================================
## STEP 12: Tumor vs Normal differential expression
## ============================================================

Idents(hepatocyte_strict) <- "site"

de_hepatocyte_strict <- FindMarkers(
  hepatocyte_strict,
  ident.1 = "Tumor",
  ident.2 = "Normal",
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

de_hepatocyte_strict$gene <- rownames(de_hepatocyte_strict)

de_hepatocyte_strict <- de_hepatocyte_strict %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

write.csv(
  de_hepatocyte_strict,
  "tables/differential_expression/HCC04_strict_hepatocyte_Tumor_vs_Normal_DE.csv",
  row.names = FALSE
)

tumor_up_strict <- de_hepatocyte_strict %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC >= 0.5
  ) %>%
  arrange(desc(avg_log2FC))

normal_up_strict <- de_hepatocyte_strict %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC <= -0.5
  ) %>%
  arrange(avg_log2FC)

write.csv(
  tumor_up_strict,
  "tables/differential_expression/HCC04_strict_hepatocyte_Tumor_upregulated_genes.csv",
  row.names = FALSE
)

write.csv(
  normal_up_strict,
  "tables/differential_expression/HCC04_strict_hepatocyte_Normal_upregulated_genes.csv",
  row.names = FALSE
)

cat("Tumor-upregulated genes:", nrow(tumor_up_strict), "\n")
cat("Normal-upregulated genes:", nrow(normal_up_strict), "\n")

print(head(tumor_up_strict, 25))

## ============================================================
## STEP 13: Volcano plot
## ============================================================

volcano_strict <- de_hepatocyte_strict %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC >= 0.5 ~ "Tumor-up",
      p_val_adj < 0.05 & avg_log2FC <= -0.5 ~ "Normal-up",
      TRUE ~ "Not significant"
    ),
    neg_log10_fdr = -log10(pmax(p_val_adj, 1e-300))
  )

label_strict <- tumor_up_strict %>%
  slice_max(
    avg_log2FC,
    n = 15,
    with_ties = FALSE
  ) %>%
  pull(gene)

p_volcano_strict <- ggplot(
  volcano_strict,
  aes(
    x = avg_log2FC,
    y = neg_log10_fdr,
    color = significance
  )
) +
  geom_point(alpha = 0.65, size = 1.1) +
  scale_color_manual(
    values = c(
      "Tumor-up" = "#D73027",
      "Normal-up" = "#4575B4",
      "Not significant" = "grey75"
    )
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40"
  ) +
  ggrepel::geom_text_repel(
    data = subset(volcano_strict, gene %in% label_strict),
    aes(label = gene),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.4,
    point.padding = 0.2
  ) +
  labs(
    title = "HCC04 hepatocyte-like cells: Tumor vs Normal",
    subtitle = "Strict marker filtering; single-patient exploratory analysis",
    x = "Average log2 fold change (Tumor / Normal)",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_classic()

ggsave(
  "figures/differential_expression/HCC04_strict_hepatocyte_Tumor_vs_Normal_volcano.png",
  p_volcano_strict,
  width = 9,
  height = 7,
  dpi = 300
)

## ============================================================
## STEP 14: GO Biological Process enrichment
## ============================================================

go_tumor_strict <- enrichGO(
  gene = tumor_up_strict$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

go_tumor_strict_df <- as.data.frame(go_tumor_strict)

write.csv(
  go_tumor_strict_df,
  "tables/pathway_enrichment/HCC04_strict_hepatocyte_GO_Tumor_upregulated.csv",
  row.names = FALSE
)

cat("Significant GO Biological Process terms:", nrow(go_tumor_strict_df), "\n")

print(head(go_tumor_strict_df, 15))

p_go_tumor_strict <- dotplot(
  go_tumor_strict,
  showCategory = 15,
  title = "HCC04: GO BP of Tumor-upregulated hepatocyte genes"
)

ggsave(
  "figures/pathway_enrichment/HCC04_strict_hepatocyte_GO_Tumor_upregulated.png",
  p_go_tumor_strict,
  width = 10,
  height = 8,
  dpi = 300
)

## ============================================================
## STEP 15: Save final objects and session information
## ============================================================

saveRDS(
  hcc04,
  "objects/HCC04_GSE149614_full_processed_seurat.rds"
)

saveRDS(
  hepatocyte_strict,
  "objects/HCC04_GSE149614_strict_hepatocyte_subset.rds"
)

writeLines(
  capture.output(sessionInfo()),
  "tables/sessionInfo.txt"
)

cat("Analysis completed successfully.\n")
