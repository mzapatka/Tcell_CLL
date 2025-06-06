# scRNA/VDJ Analysis for CLL and rLN

This repository houses the code and scripts used in "T cell landscape definition by multi-omics identifies galectin-9 as a novel immunotherapy target in chronic lymphocytic leukemia (CLL)".

## Overview

Chronic lymphocytic leukemia (CLL) is a form of leukemia characterized by the presence of abnormal B lymphocytes. This study employs multi-omics approaches to define the T cell landscape in CLL, with a focus on identifying potential therapeutic targets. Among the findings, galectin-9 is highlighted as a promising target for immunotherapy.

## This repository provides:

R scripts for data analysis and visualization. Shell scripts for preprocessing and automation tasks. Reproducible workflows for multi-omics data integration and downstream analysis.

### Key Features

Comprehensive multi-omics analysis of T cell landscapes in CLL. Identification of galectin-9 as a potential immunotherapy target. Scripts for data processing, statistical analysis, and visualization.

### Repository Structure

/R: Contains R scripts for statistical analysis and visualization.


/Shell: Contains shell scripts for preprocessing and workflow automation.


This repository contains the analysis pipeline and code for single-cell RNA-seq (scRNA) and Variable Domain of the Immune Receptor Sequencing (VDJ) data from Chronic Lymphocytic Leukemia (CLL) and reactive lymph node (rLN).

## 1\. Data Preprocessing

* **Cell Ranger**: Scripts for preprocessing scRNA and VDJ data.
  * `scRNA_cellranger.sh`
  * `scRNA_cellrange_vdj.sh`

## 2\. Seurat QC and Normalization

### a. Doublet Detection by DoubletFinder and ambient RNA removal

`scRNA_doubletFinder_RemoveBackground.R`

### c. Merging Datasets (rLN and CLL LN)

 Merge datasets using Seurat's `merge` function.

### d. Clustering and Annotation

#### 1\. Creation of CLL Object Only

*  Implement clustering and annotation for CLL dataset.
  * Supplementary Figure 6A: `umap_harmony_by_cluster.pdf`
  * Supplementary Figures 6B-D: Same as provided in `Fig4A-B_SupplFig5B,D_scRNAseq_PBandLN.R`

#### 2\. Creation of rLN Object Only

*  Implement clustering and annotation for rLN dataset.
  * Supplementary Figure 9B: `umap_harmony_by_cluster.pdf`
  * Supplementary Figure 9C: `umap_by_pbmc3kII.pdf`
  * Supplementary Figure 9E: `SupplFig9E_ViolinPlot_TIM3_GAL9.R`

#### 3\. Creation of CLL LN and PB Object

`Fig4A-B_SupplFig5B,D_scRNAseq_PBandLN.R`

### e. Export Count Matrix with Cell Type and VDJ Annotation

`.qs` object `scDataMerged_harmony2_23102024.qs` (PB and LN samples together)

## 3\. TCR Analysis

`Fig5A-F_TCR.R` and `SupplFig7A-H.R`

## 4\. Pseudotime Analysis

`Fig4F-J_Pseudotime.R`
  * Figures:
    * Fig 4F: `destiny_CD8_diff_map2D_orig.pdf` 
    * Fig 4G-H: `cheatmap_genes_with_4kmanual.pdf`
    * Fig 4I: `destiny_CD4_diff_map12.pdf`, `destiny_CD4_diff_map12zoomB.pdf` (zoom figure)
    * Fig 4J: `cheatmap_genes_with_4kmanual.pdf`
    * Supplementary Figures:
      * Suppl Fig 6F: `destiny_CD8_diff_steps.pdf`
      * Suppl Fig 6G: `destiny_CD4_diff_steps.pdf`

## 5\. Signatures Analysis

`Fig4D-E_SupplFig6E_Exh_Signatures.R`
  * Figures:
    * Fig 4D: `vln_EX_sig_1_V3.pdf`
    * Fig 4E: `vln_TPEX_sig_2_V3.pdf`
    * Supplementary Figure 6E: `umap_TPEX_sig_2.pdf`

## 6\. CellChat Analysis

### a. CLL Object Only

`Fig6A-C_SupplFig9A_CellChat_CLLonly.R`

### b. Comparison of CLL vs rLN
`Fig6D-F_SupplFig9D.R`
*  Original figure file names:
  * Fig 6D: `circos_D_interactions.pdf`
  * Fig 6E: `n_interaction_HM_2.pdf`
  * Fig 6F: `pairwise_dot_diff_interactions_candidates_bidirectional_ordered_0.01.pdf`
  * Supplementary Figure 9D: `compare_n_interactions.pdf`

## 7\. CyTOF Processing 

### a. Data Processing

`Cytof_DataProcessing.R` 

### b. Figures Generation

* **Figures Code **:
  * Fig 1B-E: `Fig1B_CYTOF_UMAP.R`, `Fig1C_SupplFig2B_CyTOF_UMAP_perMarker.R`, `Fig1D_CyTOF_Heatmapt_Markers.R`, `Fig1E_ClusterCorrHeatmapt.R`
    * File name Fig 1B: `UMAP_Seurat.pdf_noLegend.pdf` and `UMAPannottaion_TcellsWithoutLegend.pdf`
    * File name Fig 1C: `cytofALLMarkers_TcellsDaf2-001.png`
  * Fig 3A-B: `Fig3A_CYTOF_PCA.R`, `Fig3B_UMAP_perCondition.R`
  * Fig 3D: `Fig3D_CyTOFabundanceLimma.R`
  * Supplementary Figures:
    * Suppl Fig 1B: 
      * File name: `DafN_Tcells_Table_perSampleClusterCells.xlsx` and `Counts_TcellsDaf22.pdf`
    * Suppl Fig 1C: `SupplFig1C_UMAP_persample.R`
    * Suppl Fig 2A-B: `Fig1C_SupplFig2B_CyTOF_UMAP_perMarker.R`
    * Suppl Fig 4A: `SupplFig4A_SampleClustering.R`
