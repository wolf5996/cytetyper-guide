# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a learning repository for **CyteTypeR**, an R package that provides automated cell type annotation for single-cell RNA-seq data using the CyteType API.

- **CyteTypeR GitHub**: https://github.com/NygenAnalytics/CyteTypeR
- **Python version**: https://github.com/NygenAnalytics/CyteType
- **Purpose**: Tutorial and practice with CyteTypeR API-driven annotation workflow

## Repository Structure

- `CyteTypeR_guide.qmd` - Main tutorial guide (Quarto document)
- `query.json` - Example API request/response data showing cluster metadata and marker structure
- `cytetype.Rproj` - RStudio project configuration

## Required R Packages

The CyteTypeR workflow depends on:
- `Seurat` - Core scRNA-seq analysis
- `SeuratData` - Example datasets (e.g., IFNB)
- `SCpubr` - Visualization
- `tidyseurat` - Tidy data manipulation for Seurat objects
- `CyteTypeR` - Main package (install via `devtools::install_github("NygenAnalytics/CyteTypeR")`)
- `glmGamPoi` - Differential expression backend
- `tictoc` - Timing utilities

## Typical Workflow

The CyteTypeR annotation pipeline follows this sequence:

1. **Load Data**: Load scRNA-seq data (e.g., `SeuratData::LoadData(ds = "ifnb")`)

2. **Pre-processing**: SCTransform normalization, PCA, UMAP
   ```r
   obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dims = 1:20)
   ```

3. **Clustering**: Find neighbors and clusters
   ```r
   obj %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.4)
   ```

4. **Differential Expression**: Identify cluster markers
   ```r
   markers <- FindAllMarkers(obj, group.by = "seurat_clusters")
   ```

5. **Prepare Data**: Format for CyteTypeR API
   ```r
   prepped_data <- PrepareCyteTypeR(
     obj, markers,
     n_top_genes = 10,
     group_key = "seurat_clusters",
     aggregate_metadata = TRUE,
     coordinates_key = "umap"
   )
   ```

6. **Annotate**: Submit to CyteType API
   ```r
   annotated_obj <- CyteTypeR(
     obj = obj,
     prepped_data = prepped_data,
     study_context = "description of experiment",
     metadata = list(title = "...", run_label = "...", experiment_name = "...")
   )
   ```

## Output

CyteTypeR returns:
- Updated Seurat object with predicted cell type annotations in metadata
- Interactive HTML report with cluster-level visualizations and marker analysis
- Example report: https://nygen-labs-prod--cytetype-api.modal.run/report/34fac9e9-3c43-4c46-95f4-6b2994e57ada?cluster=7&tab=markers

## Development Commands

**Render Quarto document**:
```bash
quarto render CyteTypeR_guide.qmd
```

**Run R in project context**:
```bash
R
# or open cytetype.Rproj in RStudio
```

## Key Architecture Notes

- CyteTypeR acts as an API client that formats Seurat objects into the expected JSON structure (see `query.json` for example)
- The API expects: cluster labels, aggregated metadata percentages, UMAP coordinates, and top marker genes per cluster
- The `PrepareCyteTypeR()` function handles the transformation from Seurat to API format
- Annotations are added back to the Seurat object's metadata upon return
