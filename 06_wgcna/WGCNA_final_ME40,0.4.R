## ============================================================
## FRESH WGCNA PIPELINE USING COLLEAGUE SETTINGS
## Starting from:
##   1. vsd.csv
##   2. Metadata_Sample_Area_medians_WITH_GROUPS.csv
##
## Main outputs:
##   S23A_sample_dendrogram_with_trait_bars.pdf
##   S23BC_soft_threshold_Rsquare_meanConnectivity.pdf
##   S23D_Clustering_of_module_eigengenes_before_merge.pdf
##   S24A_Dynamic_and_Merged_tree_cut.pdf
##   S25_ME_vs_Combined_groups.pdf
##   ME_vs_Treatment.pdf
##   ME_vs_Genotype.pdf
##   ME_vs_Generation.pdf
##   ME_vs_phenotypes.pdf
## ============================================================

rm(list = ls())
set.seed(65)

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(tibble)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

## ============================================================
## 0. PATHS AND SETTINGS
## ============================================================

wd <- "~/Documents/DNA_seq/RNA_seq/WGCNA"
setwd(wd)

outdir <- file.path(wd, "WGCNA_final")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

vsd_file  <- "vsd.csv"
meta_file <- "Metadata_Sample_Area_medians_WITH_GROUPS.csv"

trait_cols <- c("Area week 3", "Area week 4", "Weight.of.seeds")

## Colleague-style WGCNA settings
networkType_use <- "unsigned"
hclust_method <- "average"

minModuleSize <- 40
deepSplit_use <- 2
cutHeight_use <- 0.995
pamRespectsDendro_use <- FALSE

MEdissThres <- 0.40

## ============================================================
## 1. LOAD VSD FILE
## ============================================================

vsd1 <- read.csv(
  vsd_file,
  sep = ",",
  header = TRUE,
  check.names = FALSE
)

cat("\nOriginal vsd.csv dimensions:\n")
print(dim(vsd1))

if (ncol(vsd1) >= 2) {
  vsd1[, 2:ncol(vsd1)] <- lapply(vsd1[, 2:ncol(vsd1)], function(x) {
    suppressWarnings(as.numeric(x))
  })
}

## Transpose: samples as rows, genes as columns
vsd <- t(vsd1)

## First row after transpose contains gene IDs
colnames(vsd) <- vsd[1, ]
vsd <- vsd[-1, , drop = FALSE]

vsd <- as.data.frame(vsd)
vsd$SampleName <- rownames(vsd)

vsd <- vsd %>%
  mutate(SampleName = trimws(SampleName))

cat("\nExpression table after transpose:\n")
print(dim(vsd))

## ============================================================
## 2. LOAD METADATA
## ============================================================

meta <- read.csv(
  meta_file,
  sep = ",",
  header = TRUE,
  check.names = FALSE
) %>%
  mutate(SampleName = trimws(SampleName))

cat("\nMetadata dimensions:\n")
print(dim(meta))

cat("\nMetadata columns:\n")
print(colnames(meta))

required_cols <- c(
  "SampleName",
  "Genotype",
  "Treatment",
  "Generation",
  trait_cols
)

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required metadata columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (!"Combined" %in% colnames(meta)) {
  meta <- meta %>%
    mutate(Combined = paste(Genotype, Treatment, Generation, sep = "_"))
}

meta[trait_cols] <- lapply(meta[trait_cols], function(x) {
  suppressWarnings(as.numeric(x))
})

meta <- meta %>%
  mutate(
    Genotype   = str_squish(as.character(Genotype)),
    Treatment  = str_squish(as.character(Treatment)),
    Generation = str_squish(as.character(Generation)),
    Combined   = str_squish(as.character(Combined))
  )

## Optional but helpful ordering
meta$Genotype <- factor(meta$Genotype, levels = c("WT", "gsnor"))

meta$Treatment <- factor(
  meta$Treatment,
  levels = c("Control", "Drought", "CO2", "Ozone", "Heat", "Combination")
)

meta$Generation <- factor(meta$Generation, levels = c("1", "3", "5"))

## ============================================================
## 3. JOIN EXPRESSION AND METADATA
## ============================================================

vsd_final <- dplyr::inner_join(vsd, meta, by = "SampleName")

cat("\nNumber of samples after joining expression + metadata:\n")
print(nrow(vsd_final))

if (nrow(vsd_final) == 0) {
  stop("No matching SampleName values between vsd.csv and metadata.")
}

non_gene_cols <- c(
  "SampleName",
  "Genotype",
  "Treatment",
  "Generation",
  "Combined",
  trait_cols
)

gene_cols <- setdiff(colnames(vsd_final), non_gene_cols)

cat("\nNumber of gene columns:\n")
print(length(gene_cols))

vsd_final[gene_cols] <- lapply(vsd_final[gene_cols], function(x) {
  suppressWarnings(as.numeric(x))
})

cat("\nAny NA in expression matrix after numeric conversion?\n")
print(any(is.na(vsd_final[, gene_cols])))

## ============================================================
## 4. COLLAPSE TO GROUP MEDIANS
## ============================================================

m_vsd <- vsd_final %>%
  group_by(Genotype, Treatment, Generation) %>%
  summarise(
    across(all_of(gene_cols), ~ median(.x, na.rm = TRUE)),
    across(all_of(trait_cols), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

m_vsd$SampleName <- with(
  m_vsd,
  paste(Genotype, Treatment, Generation, sep = "_")
)

m_vsd <- as.data.frame(m_vsd)

cat("\nGroup-median table dimensions:\n")
print(dim(m_vsd))

cat("\nGroup-level samples:\n")
print(m_vsd$SampleName)

## ============================================================
## 5. BUILD WGCNA MATRICES
## ============================================================

datTraits <- m_vsd %>%
  select(SampleName, all_of(trait_cols)) %>%
  column_to_rownames("SampleName") %>%
  as.data.frame()

datExpr <- m_vsd %>%
  select(SampleName, all_of(gene_cols)) %>%
  column_to_rownames("SampleName") %>%
  as.data.frame()

datExpr[] <- lapply(datExpr, function(x) suppressWarnings(as.numeric(x)))
datTraits[] <- lapply(datTraits, function(x) suppressWarnings(as.numeric(x)))

stopifnot(all(rownames(datExpr) == rownames(datTraits)))

cat("\ndatExpr dimensions before quality check:\n")
print(dim(datExpr))

cat("\ndatTraits dimensions:\n")
print(dim(datTraits))

## ============================================================
## 6. SAMPLE DENDROGRAM WITH TRAIT BARS
## ============================================================

sampleTree <- hclust(
  dist(datExpr),
  method = hclust_method
)

traitColors <- WGCNA::numbers2colors(
  datTraits,
  signed = FALSE
)

pdf(
  file.path(outdir, "S23A_sample_dendrogram_with_trait_bars.pdf"),
  width = 14,
  height = 7,
  useDingbats = FALSE
)

par(mar = c(5, 40, 3, 0))

WGCNA::plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = colnames(datTraits),
  main = "Sample dendrogram and trait heatmap"
)

dev.off()

png(
  file.path(outdir, "S23A_sample_dendrogram_with_trait_bars.png"),
  width = 2400,
  height = 1200,
  res = 180
)

par(mar = c(5, 40, 3, 0))

WGCNA::plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = colnames(datTraits),
  main = "Sample dendrogram and trait heatmap"
)

dev.off()

## ============================================================
## 7. QUALITY CHECK
## ============================================================

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(
      paste(
        "Removing genes:",
        paste(colnames(datExpr)[!gsg$goodGenes], collapse = ", ")
      )
    )
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(
      paste(
        "Removing samples:",
        paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")
      )
    )
  }
  
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  datTraits <- datTraits[rownames(datExpr), , drop = FALSE]
}

cat("\ndatExpr dimensions after goodSamplesGenes:\n")
print(dim(datExpr))

## Important:
## No top-50% variance filtering here.
## This follows the colleague-style approach.

## ============================================================
## 8. SOFT-THRESHOLDING: AUTOMATIC SELECTION
## ============================================================

powers <- c(1:14, seq(from = 15, to = 21, by = 2))

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = networkType_use
)

pdf(
  file.path(outdir, "S23BC_soft_threshold_Rsquare_meanConnectivity.pdf"),
  width = 14,
  height = 6,
  useDingbats = FALSE
)

par(mfrow = c(1, 2))
cex1 <- 0.9

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, unsigned R²",
  type = "n",
  main = "Scale independence"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

abline(h = 0.87, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)

dev.off()

write.csv(
  sft$fitIndices,
  file.path(outdir, "soft_threshold_fitIndices.csv"),
  row.names = FALSE
)

softPower <- sft$powerEstimate

cat("\nAutomatic softPower selected by pickSoftThreshold:\n")
print(softPower)

## Safety fallback only if WGCNA cannot estimate power
if (is.na(softPower)) {
  
  fi <- sft$fitIndices
  r2 <- -sign(fi[, 3]) * fi[, 2]
  candidate_powers <- fi[r2 >= 0.87, 1]
  
  if (length(candidate_powers) > 0) {
    softPower <- min(candidate_powers)
  } else {
    softPower <- 12
  }
  
  cat("\nNo automatic powerEstimate. Using fallback softPower:\n")
  print(softPower)
}

cat("\nFinal softPower used for adjacency/TOM:\n")
print(softPower)

## ============================================================
## 9. ADJACENCY AND TOM
## ============================================================

adjacency <- WGCNA::adjacency(
  datExpr,
  power = softPower,
  type = networkType_use
)

TOM <- WGCNA::TOMsimilarity(adjacency)

colnames(TOM) <- colnames(datExpr)
rownames(TOM) <- colnames(datExpr)

dissTOM <- 1 - TOM

WGCNA::collectGarbage()

## ============================================================
## 10. GENE TREE
## ============================================================

geneTree <- hclust(
  as.dist(dissTOM),
  method = hclust_method
)

pdf(
  file.path(outdir, "S24A_gene_cluster_dendrogram_only.pdf"),
  width = 14,
  height = 7,
  useDingbats = FALSE
)

plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Cluster Dendrogram",
  labels = FALSE,
  hang = 0.04
)

dev.off()

## ============================================================
## 11. DYNAMIC TREE CUT
## ============================================================

unmergedLabels <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplit_use,
  cutHeight = cutHeight_use,
  minClusterSize = minModuleSize,
  pamRespectsDendro = pamRespectsDendro_use
)

unmergedColors <- WGCNA::labels2colors(unmergedLabels)

cat("\nUnmerged module sizes:\n")
print(sort(table(unmergedColors), decreasing = TRUE))

cat("\nGrey genes before merging:\n")
print(sum(unmergedColors == "grey"))

pdf(
  file.path(outdir, "S24A_Dynamic_tree_cut.pdf"),
  width = 14,
  height = 8,
  useDingbats = FALSE
)

WGCNA::plotDendroAndColors(
  geneTree,
  unmergedColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Cluster Dendrogram"
)

dev.off()

## ============================================================
## 12. MODULE EIGENGENE CLUSTERING BEFORE MERGE
## This matches the colleague-style plot.
## ============================================================

MElist_unmerged <- WGCNA::moduleEigengenes(
  datExpr,
  colors = unmergedColors
)

MEs_unmerged <- MElist_unmerged$eigengenes

MEdiss_unmerged <- 1 - cor(MEs_unmerged, use = "p")

MEtree_unmerged <- hclust(
  as.dist(MEdiss_unmerged),
  method = hclust_method
)

pdf(
  file.path(outdir, "S23D_Clustering_of_module_eigengenes_before_merge.pdf"),
  width = 14,
  height = 8,
  useDingbats = FALSE
)

par(mar = c(8, 4, 4, 2))

plot(
  MEtree_unmerged,
  main = "Clustering of module eigengenes",
  xlab = "",
  sub = "",
  cex = 0.8
)

abline(h = 0.40, col = "red", lwd = 2)
#abline(h = 0.30, col = "green", lwd = 2)

dev.off()

## ============================================================
## 13. MERGE CLOSE MODULES
## ============================================================

merge <- WGCNA::mergeCloseModules(
  datExpr,
  unmergedColors,
  cutHeight = MEdissThres,
  verbose = 3
)

mergedColors <- merge$colors
moduleColors <- mergedColors
MEs <- merge$newMEs
MEs0 <- WGCNA::orderMEs(MEs)

cat("\nMerged module sizes:\n")
print(sort(table(mergedColors), decreasing = TRUE))

cat("\nGrey genes after merging:\n")
print(sum(mergedColors == "grey"))

cat("\nFinal modules including grey:\n")
print(length(unique(mergedColors)))

cat("\nFinal modules excluding grey:\n")
print(length(setdiff(unique(mergedColors), "grey")))

## ============================================================
## 14. DYNAMIC AND MERGED TREE CUT
## ============================================================

pdf(
  file.path(outdir, "S24A_Dynamic_and_Merged_tree_cut.pdf"),
  width = 14,
  height = 8,
  useDingbats = FALSE
)

WGCNA::plotDendroAndColors(
  geneTree,
  cbind(unmergedColors, mergedColors),
  c("DynamicTree cut", "MergedTree cut"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = FALSE,
  guideHang = 0.05,
  main = "Cluster Dendrogram"
)

dev.off()

png(
  file.path(outdir, "S24A_Dynamic_and_Merged_tree_cut.png"),
  width = 2400,
  height = 1400,
  res = 180
)

WGCNA::plotDendroAndColors(
  geneTree,
  cbind(unmergedColors, mergedColors),
  c("DynamicTree cut", "MergedTree cut"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = FALSE,
  guideHang = 0.05,
  main = "Cluster Dendrogram"
)

dev.off()

## ============================================================
## 15. MODULE EIGENGENE CLUSTERING AFTER MERGE
## ============================================================

MEdiss_merged <- 1 - cor(MEs0, use = "p")

MEtree_merged <- hclust(
  as.dist(MEdiss_merged),
  method = hclust_method
)

pdf(
  file.path(outdir, "S23D_Clustering_of_module_eigengenes_after_merge.pdf"),
  width = 14,
  height = 8,
  useDingbats = FALSE
)

par(mar = c(8, 4, 4, 2))

plot(
  MEtree_merged,
  main = "Clustering of merged module eigengenes",
  xlab = "",
  sub = "",
  cex = 0.8
)

abline(h = MEdissThres, col = "red", lwd = 2)

dev.off()

## ============================================================
## 16. PREPARE GROUP METADATA FOR MODULE-TRAIT HEATMAPS
## ============================================================

meta2 <- m_vsd %>%
  dplyr::select(
    Genotype,
    Treatment,
    Generation,
    SampleName,
    all_of(trait_cols)
  )

rownames(meta2) <- meta2$SampleName

meta2 <- meta2[rownames(datExpr), , drop = FALSE]

stopifnot(all(rownames(meta2) == rownames(datExpr)))
stopifnot(all(rownames(MEs0) == rownames(datExpr)))

meta2$Combined_like <- with(
  meta2,
  paste(Genotype, Treatment, Generation, sep = "_")
)

combined_levels <- expand.grid(
  Genotype = c("WT", "gsnor"),
  Treatment = c("Control", "Drought", "CO2", "Ozone", "Heat", "Combination"),
  Generation = c("1", "3", "5"),
  stringsAsFactors = FALSE
) %>%
  mutate(Combined_like = paste(Genotype, Treatment, Generation, sep = "_")) %>%
  pull(Combined_like)

combined_levels <- combined_levels[combined_levels %in% meta2$Combined_like]

## ============================================================
## 17. HELPER FUNCTION FOR BINARIZED GROUP HEATMAPS
## Safer version: avoids empty bin_t problem
## ============================================================

plot_ME_vs_factor <- function(fvec,
                              factor_name,
                              outfile_prefix,
                              x_order = NULL,
                              width = 20,
                              height = 14,
                              x_angle = 45,
                              show_text = TRUE) {
  
  ## -----------------------------
  ## Clean factor vector
  ## -----------------------------
  fvec <- as.character(fvec)
  fvec <- stringr::str_squish(fvec)
  
  names(fvec) <- rownames(datExpr)
  
  ## Remove NA / empty groups if present
  if (any(is.na(fvec) | fvec == "")) {
    stop("Some values in ", factor_name, " are NA or empty. Please check fvec.")
  }
  
  ## -----------------------------
  ## Decide group order
  ## -----------------------------
  if (is.null(x_order)) {
    x_order_use <- unique(fvec)
  } else {
    x_order <- as.character(x_order)
    x_order <- stringr::str_squish(x_order)
    
    ## Keep only groups actually present in fvec
    x_order_use <- x_order[x_order %in% unique(fvec)]
    
    ## If x_order gives zero match, stop with useful message
    if (length(x_order_use) == 0) {
      cat("\nValues present in fvec:\n")
      print(sort(unique(fvec)))
      
      cat("\nValues present in x_order:\n")
      print(x_order)
      
      stop(
        "No overlap between x_order and actual values in ",
        factor_name,
        ". This caused the 0-column heatmap error."
      )
    }
  }
  
  ## -----------------------------
  ## Make binary group matrix
  ## Each column = one group
  ## Each row = one sample/group-median sample
  ## -----------------------------
  f <- factor(fvec, levels = x_order_use)
  
  bin_t <- model.matrix(~ 0 + f)
  colnames(bin_t) <- levels(f)
  rownames(bin_t) <- rownames(datExpr)
  
  ## Remove columns with all zero, just in case
  bin_t <- bin_t[, colSums(bin_t, na.rm = TRUE) > 0, drop = FALSE]
  
  if (ncol(bin_t) == 0) {
    stop("bin_t has 0 columns after filtering. Check group names/order.")
  }
  
  ## -----------------------------
  ## Correlate module eigengenes with binary groups
  ## -----------------------------
  moduleTraitCor <- cor(
    MEs0,
    bin_t,
    use = "p",
    method = "pearson"
  )
  
  moduleTraitPvalue <- WGCNA::corPvalueFisher(
    moduleTraitCor,
    nrow(datExpr)
  )
  
  ## -----------------------------
  ## Row labels
  ## -----------------------------
  MEColors <- substring(names(MEs0), 3)
  MEColorNames <- paste(MEColors, sep = "")
  
  ## -----------------------------
  ## Text matrix
  ## -----------------------------
  textMatrix <- matrix(
    paste(
      signif(moduleTraitCor, 2),
      "\n(",
      signif(moduleTraitPvalue, 1),
      ")",
      sep = ""
    ),
    nrow = nrow(moduleTraitCor),
    ncol = ncol(moduleTraitCor),
    dimnames = dimnames(moduleTraitCor)
  )
  
  ## -----------------------------
  ## Save correlation and p-value tables
  ## -----------------------------
  write.csv(
    moduleTraitCor,
    file.path(outdir, paste0(outfile_prefix, "_cor.csv"))
  )
  
  write.csv(
    moduleTraitPvalue,
    file.path(outdir, paste0(outfile_prefix, "_pvalue.csv"))
  )
  
  ## -----------------------------
  ## Plot heatmap
  ## -----------------------------
  pdf(
    file.path(outdir, paste0(outfile_prefix, ".pdf")),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  par(mar = c(14, 14, 4, 6))
  
  WGCNA::labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = stringr::str_replace_all(colnames(bin_t), "_", "\n"),
    yLabels = MEColorNames,
    xLabelsAngle = x_angle,
    ySymbols = MEColorNames,
    colorLabels = TRUE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = if (show_text) textMatrix else NULL,
    setStdMargins = FALSE,
    cex.text = 0.45,
    cex.lab.x = 0.9,
    cex.lab.y = 1.0,
    cex.legendLabel = 1,
    legendLabel = "Pearson's correlation",
    zlim = c(-1, 1),
    main = paste("Module Eigengene -", factor_name, "relationships")
  )
  
  dev.off()
  
  cat("\nFinished:", outfile_prefix, "\n")
  cat("Number of groups plotted:", ncol(bin_t), "\n")
}
## ============================================================
## 18. S25-STYLE HEATMAP: MODULES VS COMBINED GROUPS
## ============================================================
## ============================================================
## Correct old-style x-axis order for S25 heatmap
## Order:
## WT first, then gsnor
## Within each genotype:
## Treatment order: Control, Drought, CO2, Ozone, Heat, Combination
## Within each treatment:
## Generation order: 1, 3, 5
## ============================================================

genotype_order <- c("WT", "gsnor")

treatment_order <- c(
  "Control",
  "Drought",
  "CO2",
  "Ozone",
  "Heat",
  "Combination"
)

generation_order <- c("1", "3", "5")

combined_levels <- unlist(
  lapply(genotype_order, function(g) {
    unlist(
      lapply(treatment_order, function(t) {
        paste(g, t, generation_order, sep = "_")
      })
    )
  })
)

## Keep only groups actually present in your data
combined_levels <- combined_levels[combined_levels %in% meta2$Combined_like]

cat("\nFinal old-style x-axis order:\n")
print(combined_levels)

plot_ME_vs_factor(
  fvec = meta2$Combined_like,
  factor_name = "Combined groups",
  outfile_prefix = "S25_ME_vs_Combined_groups",
  x_order = combined_levels,
  width = 15,
  height = 16,
  x_angle = 45,
  show_text = TRUE
)

## ============================================================
## 19. ADDITIONAL HEATMAPS
## ============================================================

plot_ME_vs_factor(
  fvec = meta2$Treatment,
  factor_name = "Treatment",
  outfile_prefix = "ME_vs_Treatment",
  x_order = c("Control", "Drought", "CO2", "Ozone", "Heat", "Combination"),
  width = 10,
  height = 14,
  x_angle = 0,
  show_text = TRUE
)

plot_ME_vs_factor(
  fvec = meta2$Genotype,
  factor_name = "Genotype",
  outfile_prefix = "ME_vs_Genotype",
  x_order = c("WT", "gsnor"),
  width = 8,
  height = 14,
  x_angle = 0,
  show_text = TRUE
)

plot_ME_vs_factor(
  fvec = meta2$Generation,
  factor_name = "Generation",
  outfile_prefix = "ME_vs_Generation",
  x_order = c("1", "3", "5"),
  width = 8,
  height = 14,
  x_angle = 0,
  show_text = TRUE
)

## ============================================================
## 20. PHENOTYPE HEATMAP
## ============================================================

meta3 <- meta2[, trait_cols, drop = FALSE]

moduleTraitCor <- cor(
  MEs0,
  meta3,
  use = "p",
  method = "pearson"
)

moduleTraitPvalue <- WGCNA::corPvalueFisher(
  moduleTraitCor,
  nrow(datExpr)
)

MEColors <- substring(names(MEs0), 3)
MEColorNames <- paste(MEColors, sep = "")

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")",
  sep = ""
)

dim(textMatrix) <- dim(moduleTraitCor)

write.csv(
  moduleTraitCor,
  file.path(outdir, "ME_vs_phenotypes_cor.csv")
)

write.csv(
  moduleTraitPvalue,
  file.path(outdir, "ME_vs_phenotypes_pvalue.csv")
)

pdf(
  file.path(outdir, "ME_vs_phenotypes.pdf"),
  width = 10,
  height = 14,
  useDingbats = FALSE
)

par(mar = c(10, 14, 4, 6))

WGCNA::labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(meta3),
  yLabels = MEColorNames,
  xLabelsAngle = 45,
  ySymbols = MEColorNames,
  colorLabels = TRUE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  cex.lab.x = 1.0,
  cex.lab.y = 1.0,
  cex.legendLabel = 1,
  legendLabel = "Pearson's correlation",
  zlim = c(-1, 1),
  main = "Module Eigengene - phenotype relationships"
)

dev.off()

## ============================================================
## 21. HUB GENES
## ============================================================

hubs <- WGCNA::chooseTopHubInEachModule(
  datExpr,
  mergedColors
)

hubs <- data.frame(
  gene_id = hubs,
  modules = names(hubs),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  hubs,
  file.path(outdir, "Hub_genes.csv"),
  row.names = FALSE
)

## ============================================================
## 22. GENE MODULE INFORMATION
## ============================================================

gene_mod <- data.frame(
  ID = colnames(datExpr),
  module_color = mergedColors,
  stringsAsFactors = FALSE
)

gene_mod$module_label <- match(
  gene_mod$module_color,
  sort(unique(gene_mod$module_color))
)

gene_mod <- gene_mod %>%
  select(ID, module_label, module_color)

write.csv(
  gene_mod,
  file.path(outdir, "genes_modules_info.csv"),
  row.names = FALSE
)

## ============================================================
## 23. MODULE SIZES
## ============================================================

module_sizes <- table(mergedColors)

module_size_df <- data.frame(
  module_color = names(module_sizes),
  size = as.integer(module_sizes)
) %>%
  arrange(desc(size))

write.csv(
  module_size_df,
  file.path(outdir, "ordered_modules_sizes.csv"),
  row.names = FALSE
)

## ============================================================
## 24. SAVE OBJECTS
## ============================================================

saveRDS(datExpr, file.path(outdir, "datExpr_final.rds"))
saveRDS(datTraits, file.path(outdir, "datTraits_final.rds"))
saveRDS(vsd_final, file.path(outdir, "vsd_final_joined.rds"))
saveRDS(m_vsd, file.path(outdir, "m_vsd_group_median.rds"))
saveRDS(gene_cols, file.path(outdir, "gene_cols.rds"))

saveRDS(adjacency, file.path(outdir, "adjacency_final.rds"))
saveRDS(TOM, file.path(outdir, "TOM_final.rds"))
saveRDS(geneTree, file.path(outdir, "geneTree_final.rds"))

saveRDS(unmergedColors, file.path(outdir, "unmergedColors_final.rds"))
saveRDS(mergedColors, file.path(outdir, "mergedColors_final.rds"))
saveRDS(moduleColors, file.path(outdir, "moduleColors_final.rds"))
saveRDS(MEs0, file.path(outdir, "MEs0_final.rds"))

## ============================================================
## 25. FINAL SUMMARY
## ============================================================

cat("\n============================================================\n")
cat("WGCNA finished successfully.\n")
cat("Output directory:\n")
print(outdir)

cat("\nFinal settings:\n")
cat("softPower:", softPower, "\n")
cat("networkType:", networkType_use, "\n")
cat("hclust method:", hclust_method, "\n")
cat("deepSplit:", deepSplit_use, "\n")
cat("cutHeight:", cutHeight_use, "\n")
cat("minModuleSize:", minModuleSize, "\n")
cat("pamRespectsDendro:", pamRespectsDendro_use, "\n")
cat("MEdissThres:", MEdissThres, "\n")

cat("\nFinal module numbers:\n")
cat("Modules including grey:", length(unique(mergedColors)), "\n")
cat("Modules excluding grey:", length(setdiff(unique(mergedColors), "grey")), "\n")
cat("Grey genes:", sum(mergedColors == "grey"), "\n")

cat("\nModule size table:\n")
print(module_size_df)

cat("============================================================\n")


##leftover heatmaps
## ============================================================
## EXTRA HELPERS FOR ADDITIONAL HEATMAPS
## ============================================================

plot_labeled_heatmap_custom <- function(cor_mat,
                                        p_mat,
                                        outfile_prefix,
                                        main_title,
                                        x_labels = colnames(cor_mat),
                                        y_labels = rownames(cor_mat),
                                        width = 12,
                                        height = 10,
                                        x_angle = 45,
                                        cex_text = 0.55,
                                        cex_lab_x = 1,
                                        cex_lab_y = 1,
                                        sep_after = NULL,
                                        bottom_mar = 12,
                                        left_mar = 12) {
  
  stopifnot(all(dim(cor_mat) == dim(p_mat)))
  
  textMatrix <- matrix(
    "",
    nrow = nrow(cor_mat),
    ncol = ncol(cor_mat),
    dimnames = dimnames(cor_mat)
  )
  
  for (i in seq_len(nrow(cor_mat))) {
    for (j in seq_len(ncol(cor_mat))) {
      if (!is.na(cor_mat[i, j]) && !is.na(p_mat[i, j])) {
        textMatrix[i, j] <- paste0(
          signif(cor_mat[i, j], 2),
          "\n(",
          signif(p_mat[i, j], 1),
          ")"
        )
      } else {
        textMatrix[i, j] <- ""
      }
    }
  }
  
  write.csv(cor_mat, file.path(outdir, paste0(outfile_prefix, "_cor.csv")))
  write.csv(p_mat,   file.path(outdir, paste0(outfile_prefix, "_pvalue.csv")))
  
  pdf(
    file.path(outdir, paste0(outfile_prefix, ".pdf")),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  par(mar = c(bottom_mar, left_mar, 4, 6))
  
  WGCNA::labeledHeatmap(
    Matrix = cor_mat,
    xLabels = x_labels,
    yLabels = y_labels,
    ySymbols = y_labels,
    colorLabels = FALSE,
    xLabelsAngle = x_angle,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = cex_text,
    cex.lab.x = cex_lab_x,
    cex.lab.y = cex_lab_y,
    cex.legendLabel = 1,
    legendLabel = "Pearson's correlation",
    zlim = c(-1, 1),
    main = main_title
  )
  
  if (!is.null(sep_after)) {
    for (s in sep_after) abline(v = s + 0.5, lwd = 2)
  }
  
  dev.off()
}


make_binary_design <- function(group_vec, group_order = NULL) {
  
  group_vec <- as.character(group_vec)
  group_vec <- stringr::str_squish(group_vec)
  
  if (is.null(group_order)) {
    group_order <- unique(group_vec)
  } else {
    group_order <- group_order[group_order %in% unique(group_vec)]
  }
  
  f <- factor(group_vec, levels = group_order)
  mm <- model.matrix(~ 0 + f)
  colnames(mm) <- levels(f)
  rownames(mm) <- names(group_vec)
  
  mm
}


cor_by_group_with_trait <- function(MEs_mat,
                                    meta_df,
                                    trait_col,
                                    group_col = "Combined_like",
                                    group_order = NULL,
                                    min_n = 2) {
  
  stopifnot(all(rownames(MEs_mat) == rownames(meta_df)))
  
  groups_present <- as.character(meta_df[[group_col]])
  
  if (is.null(group_order)) {
    group_order <- unique(groups_present)
  } else {
    group_order <- group_order[group_order %in% unique(groups_present)]
  }
  
  module_names <- colnames(MEs_mat)
  
  cor_mat <- matrix(
    NA_real_,
    nrow = length(module_names),
    ncol = length(group_order),
    dimnames = list(sub("^ME", "", module_names), group_order)
  )
  
  p_mat <- matrix(
    NA_real_,
    nrow = length(module_names),
    ncol = length(group_order),
    dimnames = list(sub("^ME", "", module_names), group_order)
  )
  
  n_mat <- matrix(
    0,
    nrow = length(module_names),
    ncol = length(group_order),
    dimnames = list(sub("^ME", "", module_names), group_order)
  )
  
  for (g in group_order) {
    idx <- which(as.character(meta_df[[group_col]]) == g)
    
    y_all <- meta_df[[trait_col]][idx]
    
    for (m in seq_along(module_names)) {
      x_all <- MEs_mat[idx, module_names[m]]
      ok <- complete.cases(x_all, y_all)
      
      x <- x_all[ok]
      y <- y_all[ok]
      n_mat[m, g] <- sum(ok)
      
      if (length(x) >= min_n &&
          stats::sd(x) > 0 &&
          stats::sd(y) > 0) {
        
        cor_mat[m, g] <- suppressWarnings(cor(x, y, method = "pearson"))
        p_mat[m, g]   <- suppressWarnings(cor.test(x, y, method = "pearson")$p.value)
      }
    }
  }
  
  list(cor = cor_mat, p = p_mat, n = n_mat)
}

## ============================================================
## OLD-STYLE COMBINED ORDER
## WT first, then gsnor
## Within genotype: treatment order
## Within treatment: generation order 1,3,5
## ============================================================

genotype_order   <- c("WT", "gsnor")
treatment_order  <- c("Control", "Drought", "CO2", "Ozone", "Heat", "Combination")
generation_order <- c("1", "3", "5")

combined_levels <- unlist(
  lapply(genotype_order, function(g) {
    unlist(
      lapply(treatment_order, function(t) {
        paste(g, t, generation_order, sep = "_")
      })
    )
  })
)

combined_levels <- combined_levels[combined_levels %in% meta2$Combined_like]


## ============================================================
## S26A: Module eigengenes vs Genotype x Generation
## Uses group-median eigengenes (MEs0)
## ============================================================

meta2$Genotype_Generation <- paste(meta2$Genotype, meta2$Generation, sep = "_")

gg_order <- c("WT_1", "WT_3", "WT_5", "gsnor_1", "gsnor_3", "gsnor_5")
gg_order <- gg_order[gg_order %in% meta2$Genotype_Generation]

bin_gg <- make_binary_design(
  group_vec = setNames(meta2$Genotype_Generation, rownames(meta2)),
  group_order = gg_order
)

cor_gg <- cor(MEs0, bin_gg, use = "p", method = "pearson")
p_gg   <- WGCNA::corPvalueFisher(cor_gg, nrow(datExpr))

plot_labeled_heatmap_custom(
  cor_mat = cor_gg,
  p_mat = p_gg,
  outfile_prefix = "S26A_ME_vs_Genotype_x_Generation",
  main_title = "Module Eigengenes – Genotype × Generation (WT first; generations ordered)",
  x_labels = c("WT 1", "WT 3", "WT 5", "gsnor 1", "gsnor 3", "gsnor 5")[seq_len(ncol(cor_gg))],
  y_labels = sub("^ME", "", rownames(cor_gg)),
  width = 11,
  height = 13,
  x_angle = 45,
  cex_text = 0.7,
  cex_lab_x = 1.0,
  cex_lab_y = 1.0,
  sep_after = 3,
  bottom_mar = 8,
  left_mar = 10
)

## ============================================================
## S26B: Module eigengenes vs Genotype x Treatment
## Uses group-median eigengenes (MEs0)
## ============================================================

meta2$Genotype_Treatment <- paste(meta2$Genotype, meta2$Treatment, sep = "_")

gt_order <- c(
  paste("WT", treatment_order, sep = "_"),
  paste("gsnor", treatment_order, sep = "_")
)

gt_order <- gt_order[gt_order %in% meta2$Genotype_Treatment]

bin_gt <- make_binary_design(
  group_vec = setNames(meta2$Genotype_Treatment, rownames(meta2)),
  group_order = gt_order
)

cor_gt <- cor(MEs0, bin_gt, use = "p", method = "pearson")
p_gt   <- WGCNA::corPvalueFisher(cor_gt, nrow(datExpr))

xlabs_gt <- c(
  paste("WT", treatment_order),
  paste("gsnor", treatment_order)
)
xlabs_gt <- xlabs_gt[seq_len(ncol(cor_gt))]

plot_labeled_heatmap_custom(
  cor_mat = cor_gt,
  p_mat = p_gt,
  outfile_prefix = "S26B_ME_vs_Genotype_x_Treatment",
  main_title = "Module Eigengenes – Genotype × Treatment (binarized)",
  x_labels = xlabs_gt,
  y_labels = sub("^ME", "", rownames(cor_gt)),
  width = 14,
  height = 12,
  x_angle = 45,
  cex_text = 0.55,
  cex_lab_x = 1.0,
  cex_lab_y = 1.0,
  sep_after = 6,
  bottom_mar = 8,
  left_mar = 10
)

## ============================================================
## SAMPLE-LEVEL MEs FOR WITHIN-GROUP PHENOTYPE CORRELATIONS
## ============================================================

meta_samples <- vsd_final %>%
  dplyr::select(
    SampleName,
    Genotype,
    Treatment,
    Generation,
    all_of(trait_cols)
  ) %>%
  mutate(
    Genotype   = as.character(Genotype),
    Treatment  = as.character(Treatment),
    Generation = as.character(Generation),
    Combined_like = paste(Genotype, Treatment, Generation, sep = "_")
  )

rownames(meta_samples) <- meta_samples$SampleName

datExpr_samples <- vsd_final %>%
  dplyr::select(SampleName, all_of(gene_cols)) %>%
  column_to_rownames("SampleName") %>%
  as.data.frame()

datExpr_samples[] <- lapply(datExpr_samples, function(x) suppressWarnings(as.numeric(x)))

meta_samples <- meta_samples[rownames(datExpr_samples), , drop = FALSE]

stopifnot(all(rownames(meta_samples) == rownames(datExpr_samples)))

## Align module colors to sample-level expression columns
## ============================================================
## FIX: keep only final WGCNA genes in sample-level expression
## ============================================================

## mergedColors corresponds to final filtered datExpr genes
names(mergedColors) <- colnames(datExpr)

## Check how many sample-level genes are not in final WGCNA datExpr
missing_genes <- setdiff(colnames(datExpr_samples), colnames(datExpr))

cat("\nGenes present in datExpr_samples but removed from final WGCNA datExpr:\n")
cat(length(missing_genes), "\n")

## Keep only genes used in the final WGCNA network
datExpr_samples <- datExpr_samples[, colnames(datExpr), drop = FALSE]

## Now module colors will match perfectly
moduleColors_for_samples <- mergedColors[colnames(datExpr_samples)]

stopifnot(length(moduleColors_for_samples) == ncol(datExpr_samples))
stopifnot(!any(is.na(moduleColors_for_samples)))
stopifnot(all(names(moduleColors_for_samples) == colnames(datExpr_samples)))

cat("\nSample-level expression matrix after matching final WGCNA genes:\n")
print(dim(datExpr_samples))

cat("\nModule colors matched successfully.\n")
print(table(moduleColors_for_samples))

MEs_samples <- WGCNA::moduleEigengenes(
  datExpr_samples,
  colors = moduleColors_for_samples
)$eigengenes

MEs_samples <- WGCNA::orderMEs(MEs_samples)

stopifnot(all(rownames(MEs_samples) == rownames(meta_samples)))

## ============================================================
## S27A: Module eigengenes vs Area week 3 within groups
## ============================================================

res_area3 <- cor_by_group_with_trait(
  MEs_mat   = MEs_samples,
  meta_df   = meta_samples,
  trait_col = "Area week 3",
  group_col = "Combined_like",
  group_order = combined_levels,
  min_n = 2
)

plot_labeled_heatmap_custom(
  cor_mat = res_area3$cor,
  p_mat   = res_area3$p,
  outfile_prefix = "S27A_ME_vs_Area_week_3_within_groups",
  main_title = "Module eigengenes vs Area week 3 (within groups)",
  x_labels = gsub("_", " ", colnames(res_area3$cor)),
  y_labels = rownames(res_area3$cor),
  width = 15,
  height = 16,
  x_angle = 45,
  cex_text = 0.38,
  cex_lab_x = 0.9,
  cex_lab_y = 0.95,
  sep_after = 18,
  bottom_mar = 12,
  left_mar = 10
)

## ============================================================
## S27B: Module eigengenes vs Area week 4 within groups
## ============================================================

res_area4 <- cor_by_group_with_trait(
  MEs_mat   = MEs_samples,
  meta_df   = meta_samples,
  trait_col = "Area week 4",
  group_col = "Combined_like",
  group_order = combined_levels,
  min_n = 2
)

plot_labeled_heatmap_custom(
  cor_mat = res_area4$cor,
  p_mat   = res_area4$p,
  outfile_prefix = "S27B_ME_vs_Area_week_4_within_groups",
  main_title = "Module eigengenes vs Area week 4 (within groups)",
  x_labels = gsub("_", " ", colnames(res_area4$cor)),
  y_labels = rownames(res_area4$cor),
  width = 15,
  height = 16,
  x_angle = 45,
  cex_text = 0.38,
  cex_lab_x = 0.9,
  cex_lab_y = 0.95,
  sep_after = 18,
  bottom_mar = 12,
  left_mar = 10
)

## ============================================================
## S28A: Module eigengenes vs Weight.of.seeds within groups
## Ignore S28B as requested
## ============================================================

res_seed <- cor_by_group_with_trait(
  MEs_mat   = MEs_samples,
  meta_df   = meta_samples,
  trait_col = "Weight.of.seeds",
  group_col = "Combined_like",
  group_order = combined_levels,
  min_n = 2
)

plot_labeled_heatmap_custom(
  cor_mat = res_seed$cor,
  p_mat   = res_seed$p,
  outfile_prefix = "S28A_ME_vs_Weight_of_seeds_within_groups",
  main_title = "Module eigengenes vs Weight.of.seeds (within groups)",
  x_labels = gsub("_", " ", colnames(res_seed$cor)),
  y_labels = rownames(res_seed$cor),
  width = 15,
  height = 16,
  x_angle = 45,
  cex_text = 0.38,
  cex_lab_x = 0.9,
  cex_lab_y = 0.95,
  sep_after = 18,
  bottom_mar = 12,
  left_mar = 10
)

write.csv(res_area3$n, file.path(outdir, "S27A_ME_vs_Area_week_3_within_groups_n.csv"))
write.csv(res_area4$n, file.path(outdir, "S27B_ME_vs_Area_week_4_within_groups_n.csv"))
write.csv(res_seed$n,  file.path(outdir, "S28A_ME_vs_Weight_of_seeds_within_groups_n.csv"))
