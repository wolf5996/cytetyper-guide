# 🧬 CyteTypeR Tutorial

A hands-on learning guide for **CyteTypeR** — an R package that automatically annotates your scRNA-seq clusters using the [CyteType API](https://github.com/NygenAnalytics/CyteTypeR).

> **TL;DR** — Feed it a Seurat object with cluster markers, get back cell type annotations + a beautiful interactive report. ✨

## 🔬 What It Does

```
Seurat Object ──→ PrepareCyteTypeR() ──→ CyteTypeR API ──→ Annotated Object + HTML Report
     │                    │                                        │
  clusters            extracts                              predicted cell
  markers            markers,                               types added to
  UMAP              metadata,                                 metadata
                   coordinates
```

## 📦 Requirements

```r
# Core
install.packages(c("Seurat", "tidyverse"))

# CyteTypeR
devtools::install_github("NygenAnalytics/CyteTypeR")

# Dependencies
install.packages(c("SCpubr", "glmGamPoi", "tictoc"))
remotes::install_github("satijalab/seurat-data")
```

## 🚀 Quick Start

```r
library(Seurat)
library(CyteTypeR)

# 1. Load and preprocess
obj <- LoadData("ifnb")
obj <- SCTransform(obj, vst.flavor = "v2")
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:20)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)

# 2. Find markers
markers <- FindAllMarkers(obj, group.by = "seurat_clusters")

# 3. Prepare for CyteTypeR
prepped_data <- PrepareCyteTypeR(obj, markers, n_top_genes = 10)

# 4. Annotate 🎯
annotated_obj <- CyteTypeR(obj, prepped_data, study_context = "PBMC IFNB stimulation")
```

## 🚀 How to Use

1. Clone this repo
2. Open `cytetype.Rproj` in RStudio (this sets your working directory and ensures reproducibility)
3. Open and run through `cytetyper_guide.qmd` step by step

## 📂 Files

| File | Description |
|------|-------------|
| `cytetype.Rproj` | 🔧 RStudio project file — open this first |
| `cytetyper_guide.qmd` | 📖 Full step-by-step tutorial (Quarto) |
| `query.json` | 🔍 Example API request/response structure |

## 🔗 Resources

| | Link |
|---|------|
| 🧪 CyteTypeR (R) | [github.com/NygenAnalytics/CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR) |
| 🐍 CyteType (Python) | [github.com/NygenAnalytics/CyteType](https://github.com/NygenAnalytics/CyteType) |
| 📊 Example Report | [Interactive HTML Report](https://nygen-labs-prod--cytetype-api.modal.run/report/34fac9e9-3c43-4c46-95f4-6b2994e57ada) |

---

Made with ☕ while learning single-cell bioinformatics.
