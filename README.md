# CyteTypeR Tutorial

A learning repository for **CyteTypeR**, an R package for automated cell type annotation of single-cell RNA-seq data using the CyteType API.

## Overview

This tutorial demonstrates the CyteTypeR workflow for annotating scRNA-seq clusters with predicted cell types. The package submits cluster markers to the CyteType API and returns annotations with an interactive HTML report.

## Requirements

```r
# Core packages
install.packages(c("Seurat", "tidyverse"))

# CyteTypeR
devtools::install_github("NygenAnalytics/CyteTypeR")

# Additional dependencies
install.packages(c("SCpubr", "glmGamPoi", "tictoc"))
remotes::install_github("satijalab/seurat-data")
```

## Quick Start

```r
library(Seurat)
library(CyteTypeR)

# 1. Load and preprocess data
obj <- LoadData("ifnb") %>%
  SCTransform() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.4)

# 2. Find markers
markers <- FindAllMarkers(obj, group.by = "seurat_clusters")

# 3. Prepare for CyteTypeR
prepped_data <- PrepareCyteTypeR(obj, markers, n_top_genes = 10)

# 4. Annotate
annotated_obj <- CyteTypeR(obj, prepped_data, study_context = "PBMC IFNB stimulation")
```

## Resources

- [CyteTypeR GitHub](https://github.com/NygenAnalytics/CyteTypeR)
- [CyteType Python](https://github.com/NygenAnalytics/CyteType)
- [Example Report](https://nygen-labs-prod--cytetype-api.modal.run/report/34fac9e9-3c43-4c46-95f4-6b2994e57ada)

## Files

- `cytetyper_guide.qmd` - Full tutorial (Quarto document)
- `query.json` - Example API request structure
