# HCC04 scRNA-seq: Tumor vs Normal Hepatocyte Analysis

## Overview
Single-cell RNA-seq analysis of matched tumor and non-tumor liver samples from HCC04 in GEO accession GSE149614.

## Dataset
- GEO accession: GSE149614
- Study: A Single-Cell Atlas of the Multicellular Ecosystem of Primary and Metastatic Hepatocellular Carcinoma
- Technology: 10x Genomics Chromium Single Cell 3'
- Patient: HCC04
- Initial cells: 6,897 (3,396 Normal; 3,501 Tumor)
- Strict hepatocyte-like subset: 2,289 cells (211 Normal; 2,078 Tumor)

## Workflow
1. Memory-safe loading of the GEO count matrix using `data.table::fread()`
2. Metadata-based selection of matched HCC04 tumor and normal cells
3. QC, normalization, variable-feature selection, PCA, graph clustering, and UMAP
4. Cluster-marker discovery and validation against GEO-provided cell labels
5. Cell-type composition comparison between Normal and Tumor
6. Strict hepatocyte-like cell filtering
7. Tumor-vs-normal differential expression and GO Biological Process enrichment

## Main Results
- Recovered hepatocyte, T/NK, myeloid, endothelial, B-cell, and fibroblast populations.
- Tumor-associated hepatocyte-like cells showed increased biosynthetic and mitochondrial programs.
- Tumor-upregulated genes included ALDH3A1, CHI3L1, CYP17A1, GSTA3, CDC20, UBE2C, CENPM, GAGE1, and GAGE13.
- Enriched pathways included ribosome biogenesis, rRNA processing, translation, RNA splicing, mitochondrial translation, and mitochondrial organization.

## Limitation
This is an exploratory cell-level analysis from one matched patient. It is not a patient-level statistical analysis and does not establish universal HCC biomarkers.

## Repository Structure
- `scripts/`: R workflow
- `figures/`: QC, clustering, markers, composition, DE, and enrichment plots
- `tables/`: marker, DE, composition, and enrichment output tables

## Citation
Lu Y, Yang A, Quan C, Pan Y, et al. A single-cell atlas of the multicellular ecosystem of primary and metastatic hepatocellular carcinoma. Nature Communications. 2022;13:4594. PMID: 35933472.
