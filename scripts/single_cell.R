## ===== STEP 1: Load libraries =====
library(data.table)
library(Seurat)
library(dplyr)
library(ggplot2)

## ===== STEP 2: Load metadata (small file, has existing annotations) =====
meta <- read.table("GSE149614_HCC.metadata.updated.txt",
                   header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
dim(meta)                # 71915 cells x 7 metadata columns
colnames(meta)            # sample, res.3, site, patient, stage, virus, celltype
table(meta$site)          # Lymph / Normal / PVTT / Tumor counts
table(meta$patient, meta$site)   # which patients have matched Tumor + Normal

## ===== STEP 3: Pick a pilot patient with matched Tumor + Normal =====
pilot_patient <- "HCC04"
wanted_cells <- rownames(meta)[meta$patient == pilot_patient & meta$site %in% c("Tumor", "Normal")]
length(wanted_cells)                  # 6897
table(meta[wanted_cells, "site"])     # Normal 3396, Tumor 3501

## ===== STEP 4: Load ONLY those cells from the 3.5GB count matrix (memory-safe) =====
select_cols <- c("V1", wanted_cells)
counts_pilot <- fread("GSE149614_HCC.scRNAseq.S71915.count.txt", select = select_cols)
counts_pilot <- as.data.frame(counts_pilot)
rownames(counts_pilot) <- counts_pilot$V1
counts_pilot$V1 <- NULL
dim(counts_pilot)          # 25712 genes x 6897 cells

## ===== STEP 5: Build Seurat object + attach metadata =====
hcc04 <- CreateSeuratObject(counts = as.matrix(counts_pilot), project = "HCC04_pilot",
                            min.cells = 3, min.features = 200)

pilot_meta <- meta[colnames(hcc04), ]
hcc04 <- AddMetaData(hcc04, metadata = pilot_meta)

table(hcc04$site)         # Normal 3396, Tumor 3501
table(hcc04$celltype)     # B, Endothelial, Fibroblast, Hepatocyte, Myeloid, T/NK