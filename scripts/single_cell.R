# STEP 1: Create output folders

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

# STEP 2: Load required packages

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

# STEP 3: Load GEO cell metadata

meta <- read.table(
  "metadata/GSE149614_HCC.metadata.updated.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# STEP 4: Inspect metadata and patient tissue distribution

dim(meta)
colnames(meta)
table(meta$site)
table(meta$patient, meta$site)
head(meta)

# STEP 5: Select HCC04 Normal and Tumor cells

pilot_patient <- "HCC04"

wanted_cells <- rownames(meta)[
  meta$patient == pilot_patient &
    meta$site %in% c("Normal", "Tumor")
]

length(wanted_cells)
table(meta[wanted_cells, "site"])

# STEP 6: Selectively load HCC04 count-matrix columns

select_cols <- c("V1", wanted_cells)

counts_pilot <- fread(
  "data/GSE149614_HCC.scRNAseq.S71915.count.txt",
  select = select_cols,
  showProgress = TRUE
)

gene_names <- counts_pilot$V1
counts_pilot[, V1 := NULL]

# STEP 7: Convert selected count matrix to sparse format

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

dim(counts_sparse)

# STEP 8: Create Seurat object

hcc04 <- CreateSeuratObject(
  counts = counts_sparse,
  project = "GSE149614_HCC04",
  min.cells = 3,
  min.features = 200
)

rm(counts_sparse)
gc()

hcc04

# STEP 9: Attach GEO cell metadata

pilot_meta <- meta[colnames(hcc04), , drop = FALSE]

hcc04 <- AddMetaData(
  object = hcc04,
  metadata = pilot_meta
)

table(hcc04$site)
table(hcc04$sample)
table(hcc04$celltype)

# STEP 10: Calculate QC metrics

hcc04[["percent.mt"]] <- PercentageFeatureSet(
  hcc04,
  pattern = "^MT-"
)

hcc04[["percent.ribo"]] <- PercentageFeatureSet(
  hcc04,
  pattern = "^RP[SL]"
)

summary(hcc04$nFeature_RNA)
summary(hcc04$nCount_RNA)
summary(hcc04$percent.mt)
summary(hcc04$percent.ribo)

# STEP 11: Generate QC plots

p_qc_violin <- VlnPlot(
  hcc04,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "site",
  ncol = 3,
  pt.size = 0.05
)

p_qc_counts_features <- FeatureScatter(
  hcc04,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "site"
)

p_qc_counts_mt <- FeatureScatter(
  hcc04,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "site"
)

p_qc_violin
p_qc_counts_features
p_qc_counts_mt

ggsave(
  "figures/HCC04_QC_violin.png",
  p_qc_violin,
  width = 12,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/HCC04_QC_counts_vs_features.png",
  p_qc_counts_features,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/HCC04_QC_counts_vs_mitochondrial.png",
  p_qc_counts_mt,
  width = 7,
  height = 6,
  dpi = 300
)

# STEP 12: Filter low-quality cells

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

table(hcc04$site)
table(hcc04$celltype)

# STEP 13: Normalize expression data

hcc04 <- NormalizeData(
  hcc04,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# STEP 14: Identify variable genes

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

p_variable_features

ggsave(
  "figures/HCC04_variable_features.png",
  p_variable_features,
  width = 8,
  height = 6,
  dpi = 300
)

# STEP 15: Scale data and run PCA

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

p_elbow

ggsave(
  "figures/HCC04_elbow_plot.png",
  p_elbow,
  width = 7,
  height = 5,
  dpi = 300
)

# STEP 16: Build neighbour graph

n_pcs <- 15

hcc04 <- FindNeighbors(
  hcc04,
  dims = 1:n_pcs
)

# STEP 17: Cluster cells and calculate UMAP

hcc04 <- FindClusters(
  hcc04,
  resolution = 0.5
)

hcc04 <- RunUMAP(
  hcc04,
  dims = 1:n_pcs,
  seed.use = 123
)

table(Idents(hcc04))

# STEP 18: Save UMAP visualizations

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
  ggtitle("HCC04: Tumor versus Normal")

p_celltype <- DimPlot(
  hcc04,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("HCC04: Published cell-type annotations")

p_clusters
p_site
p_celltype

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

# STEP 19: Identify Seurat cluster markers

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

top5_markers

# STEP 20: Validate clusters against GEO cell labels

cluster_celltype_table <- table(
  Seurat_cluster = Idents(hcc04),
  Published_celltype = hcc04$celltype
)

cluster_celltype_table

cluster_celltype_proportion <- round(
  prop.table(cluster_celltype_table, margin = 1),
  2
)

cluster_celltype_proportion

write.csv(
  as.data.frame(cluster_celltype_table),
  "tables/markers/HCC04_cluster_vs_published_celltype.csv",
  row.names = FALSE
)

# STEP 21: Create marker heatmap and canonical marker DotPlot

heatmap_genes <- unique(top5_markers$gene)

p_heatmap <- DoHeatmap(
  hcc04,
  features = heatmap_genes,
  size = 3
) +
  NoLegend()

p_heatmap

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

p_dotplot

ggsave(
  "figures/markers/HCC04_canonical_marker_dotplot.png",
  p_dotplot,
  width = 13,
  height = 7,
  dpi = 300
)

# STEP 22: Compare Normal and Tumor cell-type composition

celltype_comp <- hcc04@meta.data %>%
  count(site, celltype, name = "n_cells") %>%
  group_by(site) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

celltype_comp

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
    title = "HCC04 cell-type composition: Normal vs Tumor",
    x = NULL,
    y = "Cell proportion",
    fill = "Published cell type"
  ) +
  theme_classic()

p_composition

ggsave(
  "figures/composition/HCC04_celltype_composition.png",
  p_composition,
  width = 8,
  height = 5,
  dpi = 300
)

# STEP 23: Select published hepatocytes

hepatocytes <- subset(
  hcc04,
  subset = celltype == "Hepatocyte"
)

table(hepatocytes$site)
ncol(hepatocytes)

# STEP 24: Apply strict hepatocyte-like filtering

p_hepatocyte_markers <- FeaturePlot(
  hcc04,
  features = c("ALB", "APOA1", "TTR", "PTPRC", "MS4A1", "CD79A"),
  ncol = 3
)

p_hepatocyte_markers

ggsave(
  "figures/markers/HCC04_hepatocyte_filter_markers.png",
  p_hepatocyte_markers,
  width = 12,
  height = 8,
  dpi = 300
)

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

table(hepatocyte_strict$site)
ncol(hepatocyte_strict)

saveRDS(
  hepatocyte_strict,
  "objects/HCC04_strict_hepatocyte_subset.rds"
)

# STEP 25: Run Tumor vs Normal differential expression

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

# STEP 26: Export tumor- and normal-upregulated genes

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

head(tumor_up_strict, 25)
head(normal_up_strict, 25)

# STEP 27: Create strict hepatocyte volcano plot

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

p_volcano_strict

ggsave(
  "figures/differential_expression/HCC04_strict_hepatocyte_Tumor_vs_Normal_volcano.png",
  p_volcano_strict,
  width = 9,
  height = 7,
  dpi = 300
)

# STEP 28: GO Biological Process enrichment and save final outputs

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

nrow(go_tumor_strict_df)
head(go_tumor_strict_df, 15)

write.csv(
  go_tumor_strict_df,
  "tables/pathway_enrichment/HCC04_strict_hepatocyte_GO_Tumor_upregulated.csv",
  row.names = FALSE
)

p_go_tumor_strict <- dotplot(
  go_tumor_strict,
  showCategory = 15,
  title = "HCC04: GO BP of Tumor-upregulated hepatocyte genes"
)

p_go_tumor_strict

ggsave(
  "figures/pathway_enrichment/HCC04_strict_hepatocyte_GO_Tumor_upregulated.png",
  p_go_tumor_strict,
  width = 10,
  height = 8,
  dpi = 300
)

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

cat("HCC04 scRNA-seq analysis completed successfully.\n")
