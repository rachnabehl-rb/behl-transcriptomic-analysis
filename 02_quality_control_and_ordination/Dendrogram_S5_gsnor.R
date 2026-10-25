rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(genefilter)
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(stringr)
  library(tibble)
})

register(MulticoreParam(8))

# ============================================================
# 1. Input
# ============================================================

count_file <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Count_matrix.txt"

# Same directory as your WT code
out_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/WT_only_PCA_DAPC_exact_old_input"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ntop <- 2000

# ============================================================
# 2. Read count matrix
# ============================================================

df <- read.table(
  count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# ============================================================
# 3. Exact gsnor sample order
# ============================================================

gsnor_gen1_samples <- c(
  "gs1C", "gs3C", "gs6C", "gs7C", "gs8C",
  "gs1D", "gs2D", "gs4D", "gs5D",
  "gs3CO", "gs4CO", "gs8CO", "gs9CO",
  "gs1O", "gs3O", "gs5O", "gs6O",
  "gs5H", "gs7H", "gs8H", "gs10H",
  "gs2T", "gs6T", "gs7T", "gs10T"
)

gsnor_gen3_samples <- c(
  "gs1C4_2", "gs3C4_1", "gs6C2_1", "gs6C2_3", "gs7C3_2",
  "gs1D3_1", "gs2D1_2", "gs4D3_2", "gs5D3_3",
  "gs3CO3_2", "gs4CO4_2", "gs8CO1_3", "gs9CO1_2",
  "gs1O2_2", "gs3O3_3", "gs5O4_2", "gs6O2_2",
  "gs5H3_2", "gs7H4_4", "gs8H1_2", "gs10H3_4",
  "gs2T4_2", "gs6T3_1", "gs7T2_4", "gs10T2_5"
)

gsnor_gen5_samples <- c(
  "gs1C4_212", "gs3C4_111", "gs6C2_311", "gs6C2_312", "gs7C3_211",
  "gs1D3_112", "gs2D1_211", "gs4D3_221", "gs5D3_322",
  "gs3CO3_212", "gs4CO4_222", "gs8CO1_322", "gs9CO1_212",
  "gs1O2_211", "gs3O3_322", "gs5O4_212", "gs6O2_221",
  "gs5H3_221", "gs7H4_422", "gs8H1_221", "gs10H3_412",
  "gs2T4_211", "gs6T3_111", "gs7T2_412", "gs10T2_522"
)

gsnor_samples <- c(
  gsnor_gen1_samples,
  gsnor_gen3_samples,
  gsnor_gen5_samples
)

missing_samples <- setdiff(gsnor_samples, colnames(df))

if (length(missing_samples) > 0) {
  stop(
    "These gsnor samples are missing from the count matrix:\n",
    paste(missing_samples, collapse = ", ")
  )
}

counts <- dplyr::select(df, all_of(gsnor_samples))

# ============================================================
# 4. Condition factor
# ============================================================

condition <- factor(
  c(
    rep("gsnorgen1_Control", 5),
    rep("gsnorgen1_Drought", 4),
    rep("gsnorgen1_CO2", 4),
    rep("gsnorgen1_Ozone", 4),
    rep("gsnorgen1_Heat", 4),
    rep("gsnorgen1_Combination", 4),
    
    rep("gsnorgen3_Control", 5),
    rep("gsnorgen3_Drought", 4),
    rep("gsnorgen3_CO2", 4),
    rep("gsnorgen3_Ozone", 4),
    rep("gsnorgen3_Heat", 4),
    rep("gsnorgen3_Combination", 4),
    
    rep("gsnorgen5_Control", 5),
    rep("gsnorgen5_Drought", 4),
    rep("gsnorgen5_CO2", 4),
    rep("gsnorgen5_Ozone", 4),
    rep("gsnorgen5_Heat", 4),
    rep("gsnorgen5_Combination", 4)
  ),
  levels = c(
    "gsnorgen1_Control", "gsnorgen1_Drought", "gsnorgen1_CO2",
    "gsnorgen1_Ozone", "gsnorgen1_Heat", "gsnorgen1_Combination",
    "gsnorgen3_Control", "gsnorgen3_Drought", "gsnorgen3_CO2",
    "gsnorgen3_Ozone", "gsnorgen3_Heat", "gsnorgen3_Combination",
    "gsnorgen5_Control", "gsnorgen5_Drought", "gsnorgen5_CO2",
    "gsnorgen5_Ozone", "gsnorgen5_Heat", "gsnorgen5_Combination"
  )
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# ============================================================
# 5. Batch assignment
# Batch1 = gsnor G1 + G3
# Batch2 = gsnor G5
# ============================================================

coldata$batch <- factor(
  ifelse(
    rownames(coldata) %in% c(gsnor_gen1_samples, gsnor_gen3_samples),
    "Batch1",
    "Batch2"
  ),
  levels = c("Batch1", "Batch2")
)

# ============================================================
# 6. Group / generation / treatment for plotting
# ============================================================

coldata$generation <- factor(
  ifelse(
    rownames(coldata) %in% gsnor_gen1_samples, "Gen1",
    ifelse(
      rownames(coldata) %in% gsnor_gen3_samples,
      "Gen3",
      "Gen5"
    )
  ),
  levels = c("Gen1", "Gen3", "Gen5")
)

coldata$group <- factor(
  ifelse(
    rownames(coldata) %in% gsnor_gen1_samples, "gsnor Gen 1",
    ifelse(
      rownames(coldata) %in% gsnor_gen3_samples,
      "gsnor Gen 3",
      "gsnor Gen5"
    )
  ),
  levels = c("gsnor Gen 1", "gsnor Gen 3", "gsnor Gen5")
)

control_samples <- c(
  "gs1C", "gs3C", "gs6C", "gs7C", "gs8C",
  "gs1C4_2", "gs3C4_1", "gs6C2_1", "gs6C2_3", "gs7C3_2",
  "gs1C4_212", "gs3C4_111", "gs6C2_311", "gs6C2_312", "gs7C3_211"
)

drought_samples <- c(
  "gs1D", "gs2D", "gs4D", "gs5D",
  "gs1D3_1", "gs2D1_2", "gs4D3_2", "gs5D3_3",
  "gs1D3_112", "gs2D1_211", "gs4D3_221", "gs5D3_322"
)

CO2_samples <- c(
  "gs3CO", "gs4CO", "gs8CO", "gs9CO",
  "gs3CO3_2", "gs4CO4_2", "gs8CO1_3", "gs9CO1_2",
  "gs3CO3_212", "gs4CO4_222", "gs8CO1_322", "gs9CO1_212"
)

ozone_samples <- c(
  "gs1O", "gs3O", "gs5O", "gs6O",
  "gs1O2_2", "gs3O3_3", "gs5O4_2", "gs6O2_2",
  "gs1O2_211", "gs3O3_322", "gs5O4_212", "gs6O2_221"
)

heat_samples <- c(
  "gs5H", "gs7H", "gs8H", "gs10H",
  "gs5H3_2", "gs7H4_4", "gs8H1_2", "gs10H3_4",
  "gs5H3_221", "gs7H4_422", "gs8H1_221", "gs10H3_412"
)

combination_samples <- c(
  "gs2T", "gs6T", "gs7T", "gs10T",
  "gs2T4_2", "gs6T3_1", "gs7T2_4", "gs10T2_5",
  "gs2T4_211", "gs6T3_111", "gs7T2_412", "gs10T2_522"
)

coldata$treatment <- factor(
  ifelse(
    rownames(coldata) %in% control_samples, "Control",
    ifelse(
      rownames(coldata) %in% drought_samples, "Drought",
      ifelse(
        rownames(coldata) %in% CO2_samples, "CO2",
        ifelse(
          rownames(coldata) %in% ozone_samples, "Ozone",
          ifelse(
            rownames(coldata) %in% heat_samples, "Heat",
            "Combination"
          )
        )
      )
    )
  )
)

cat("\ngsnor batch table:\n")
print(table(coldata$batch, coldata$group))

cat("\ngsnor treatment table:\n")
print(table(coldata$generation, coldata$treatment))

write.csv(
  coldata,
  file.path(out_dir, "gsnor_coldata_PCA.csv")
)

# ============================================================
# 7. DESeq2 VST
# ============================================================

counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- estimateSizeFactors(dds)

idx <- rowSums(counts(dds, normalized = TRUE) >= 3) >= 2
dds <- dds[idx, ]

vsd <- vst(dds, blind = FALSE)

# ============================================================
# 8. Batch correction for PCA visualization
# ============================================================

vsd_bc <- vsd

if ("batch" %in% colnames(colData(vsd_bc))) {
  assay(vsd_bc) <- limma::removeBatchEffect(
    assay(vsd_bc),
    batch = vsd_bc$batch
  )
}

obj <- vsd_bc
# ============================================================
# WT normal dendrogram + circular dendrogram
# Use directly after:
# obj <- vsd_bc
# ============================================================

suppressPackageStartupMessages({
  library(dendextend)
  library(circlize)
})

# ------------------------------------------------------------
# 1. Colour labels by sample name / treatment
# CO must be checked before C and O
# ------------------------------------------------------------

assign_color <- function(sample_names) {
  
  sapply(sample_names, function(sample) {
    
    if (grepl("CO", sample)) {
      "darkgreen"       # CO2
    } else if (grepl("C", sample)) {
      "black"           # Control
    } else if (grepl("D", sample)) {
      "blue"            # Drought
    } else if (grepl("O", sample)) {
      "red"             # Ozone
    } else if (grepl("H", sample)) {
      "magenta"         # Heat
    } else if (grepl("T", sample)) {
      "brown"           # Combination
    } else {
      "black"
    }
  })
}

# ------------------------------------------------------------
# 2. Use full VST batch-corrected matrix directly
# No ntop selection
# ------------------------------------------------------------

mat <- assay(obj)

# samples x genes distance
sampleDists <- dist(
  t(mat),
  method = "euclidean"
)

# ------------------------------------------------------------
# 3. Normal dendrogram first
# ------------------------------------------------------------

sampleTree <- as.dendrogram(
  hclust(sampleDists, method = "complete")
)

# Colour branches by condition if available
if ("condition" %in% colnames(colData(obj))) {
  k_branch <- length(levels(droplevels(colData(obj)$condition)))
} else {
  k_branch <- 18
}

dend_colored <- color_branches(
  sampleTree,
  k = k_branch
)

# Colour sample labels
label_colors <- assign_color(labels(dend_colored))

labels_colors(dend_colored) <- label_colors
labels_cex(dend_colored) <- 0.55

# ------------------------------------------------------------
# 4. Plot normal dendrogram
# ------------------------------------------------------------

par(mar = c(8, 4, 4, 2))

plot(
  dend_colored,
  main = "Dendrogram of Sample Clustering",
  ylab = "Distance"
)

legend(
  "topright",
  legend = c("Control", "Drought", "CO2", "Ozone", "Heat", "Combination"),
  col = c("black", "blue", "darkgreen", "red", "magenta", "brown"),
  pch = 19,
  bty = "n",
  cex = 0.8,
  title = "Treatment"
)

# ------------------------------------------------------------
# 5. Circular dendrogram from same dendrogram object
# ------------------------------------------------------------

circos.clear()

circos.par(
  start.degree = 0,
  gap.after = 0,
  track.margin = c(0, 0.02),
  cell.padding = c(0, 0, 0, 0)
)

par(mar = rep(0, 4))

circlize_dendrogram(
  dend_colored,
  labels_track_height = 0.2,
  dend_track_height = 0.4
)

circos.clear()

