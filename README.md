# ARV Cell Transcriptomics Analysis Toolkit (2024)

This repository contains a collection of R scripts developed for analyzing single-cell transcriptomics data. It emphasizes identifying cell-type-specific gene expression patterns, conducting dimensionality reduction analyses, and visualizing complex intersections of differentially expressed genes (DEGs).

## Features

- **Differential Expression Analysis:** Identify and visualize DEGs across multiple cell types.
- **Principal Component Analysis (PCA):** Reduce data dimensionality and explore gene expression variability across cell populations.
- **Gene Ontology (GO) Enrichment:** Discover biological processes associated with cell-type-specific DEGs.
- **Intersection Visualization:** Generate UpSet plots to visualize shared and unique gene expression patterns among cell types.

## Repository Structure

| Script/File                  | Description                                                  |
|------------------------------|--------------------------------------------------------------|
| `all_DEGs_plots.R`           | Generates plots of DEGs across various cell types.           |
| `pca_cell_type_only.R`       | Performs PCA specifically focused on cell-type differences.  |
| `topGo_unique_DEGs_Celltype.R`| Conducts GO enrichment analysis using topGO.                 |
| `upset_overall.R`            | Creates UpSet plots for DEG intersections among cell types.  |
| `ASC_scripts.zip`            | Additional scripts/data related to adipose-derived stem cells (ASCs). |

## Requirements

- R (version 4.0 or later recommended)
- Required packages: `topGO`, `ggplot2`, `UpSetR`, `Seurat`, and other dependencies as specified within scripts.

## Usage

Clone the repository and ensure all dependencies are installed. Execute scripts within an R environment or IDE (e.g., RStudio).

```bash
git clone https://github.com/Zubair2021/ARV_Cell_Transcriptomics_2024.git
cd ARV_Cell_Transcriptomics_2024
```

Run the desired scripts directly within R:

```R
source("all_DEGs_plots.R")
source("pca_cell_type_only.R")
source("topGo_unique_DEGs_Celltype.R")
source("upset_overall.R")
```

## Contributions

Contributions to improve the analysis scripts and add additional features are welcome. Please fork the repository and submit a pull request.
