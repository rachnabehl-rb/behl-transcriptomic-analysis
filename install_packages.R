# Install the R packages used across the analysis scripts.
# Usage: source("install_packages.R")

cran_pkgs <- c(
  "dplyr", "tidyr", "tibble", "stringr", "readr", "readxl", "purrr", "tidyverse",
  "data.table", "glue", "ggplot2", "ggrepel", "ggnewscale", "ggalluvial",
  "patchwork", "gridExtra", "scales", "scattermore", "RColorBrewer", "pheatmap",
  "circlize", "dendextend", "VennDiagram", "UpSetR", "adegenet", "WGCNA",
  "car", "emmeans", "multcomp"
)

bioc_pkgs <- c(
  "DESeq2", "BiocParallel", "genefilter", "limma",
  "AnnotationDbi", "org.At.tair.db", "clusterProfiler",
  # WGCNA dependencies that are distributed through Bioconductor
  "impute", "preprocessCore", "GO.db"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

still_missing <- c(bioc_pkgs, cran_pkgs)[
  !vapply(c(bioc_pkgs, cran_pkgs), requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing)) {
  warning("Could not install: ", paste(still_missing, collapse = ", "))
} else {
  message("All packages installed.")
}
