# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tutorial repository for **CyteTypeR**, an R package for automated cell type annotation of scRNA-seq data via the CyteType API.

- **CyteTypeR GitHub**: https://github.com/NygenAnalytics/CyteTypeR
- **Python version**: https://github.com/NygenAnalytics/CyteType

## Development Commands

**Render the tutorial:**
```bash
quarto render cytetyper_guide.qmd
```

**Install CyteTypeR (from GitHub):**
```r
devtools::install_github("NygenAnalytics/CyteTypeR")
```

## Key Files

- `cytetyper_guide.qmd` — Main tutorial (Quarto). Walks through the full pipeline: load IFNB data → SCTransform → PCA/UMAP → clustering → FindAllMarkers → PrepareCyteTypeR → CyteTypeR API call.
- `query.json` — Real API request/response example (~4.7MB). Shows the JSON structure `PrepareCyteTypeR()` produces: `clusterLabels`, `clusterMetadata` (aggregated metadata percentages), `coordinates` (UMAP embeddings), and `markers` (top genes per cluster with log2FC/p-values).

## Architecture

CyteTypeR is an API client. The workflow has two key transformation boundaries:

1. **Seurat → JSON** (`PrepareCyteTypeR()`): Extracts cluster labels, aggregates categorical metadata into per-cluster percentages, pulls UMAP coordinates, and selects top N marker genes per cluster. The output JSON structure matches `query.json`.

2. **API → Seurat** (`CyteTypeR()`): Sends the prepared JSON to the CyteType API with a `study_context` string and metadata. Returns predicted cell type annotations that get added to the Seurat object's metadata columns, plus an interactive HTML report URL.

The API expects these fields in `input_data`: `clusterLabels` (map of index→label), `clusterMetadata` (nested per-cluster percentage breakdowns of each metadata column), `coordinates` (UMAP x/y per cell), and `markers` (per-cluster gene lists with stats).
