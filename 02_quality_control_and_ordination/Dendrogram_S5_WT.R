rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(genefilter)
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(adegenet)
  library(stringr)
  library(tibble)
})

register(MulticoreParam(8))

# ============================================================
# 1. Input
# ============================================================

count_file <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Count_matrix.txt"

out_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Dendrogram"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ntop <- 2000
seed_use <- 1
n_rep_xval <- 50
n_rep_cv <- 50

# Optional: only use this if you know the old WT n.pca value
force_n_pca_WT <- NULL

# Samples to exclude if present
bad_samples <- c(
  "WT_Pf_2-1-4",
  "WT7O_new",
  "WT8C_new"
)

# ============================================================
# 2. Read count matrix
# ============================================================

df <- read.table(
  count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

df <- df[, !colnames(df) %in% bad_samples, drop = FALSE]

cat("\nExcluded samples still present?\n")
print(intersect(bad_samples, colnames(df)))

# ============================================================
# 3. Exact WT sample order from old PCA script
# ============================================================

WT_samples <- c(
  # Gen1
  "WT1C", "WT3C", "WT5C", "WT6C", "WT8C", "WT10C",
  "WT3D", "WT4D", "WT6D", "WT10D",
  "WT3CO", "WT4CO", "WT5CO", "WT7CO",
  "WT1O", "WT2O", "WT5O", "WT7O",
  "WT1H", "WT5H", "WT6H", "WT8H",
  "WT1T", "WT2T", "WT3T", "WT6T",
  
  # Gen3
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4",
  
  # Gen5
  "WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111",
  "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211",
  "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212",
  "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112",
  "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221",
  "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422"
)

missing_samples <- setdiff(WT_samples, colnames(df))

if (length(missing_samples) > 0) {
  stop(
    "These WT samples are missing from the count matrix:\n",
    paste(missing_samples, collapse = ", ")
  )
}

counts <- dplyr::select(df, all_of(WT_samples))

# ============================================================
# 4. Condition exactly as old WT PCA script
# ============================================================

condition <- factor(
  c(
    rep("WTgen1_Control", 6),
    rep("WTgen1_Drought", 4),
    rep("WTgen1_CO2", 4),
    rep("WTgen1_Ozone", 4),
    rep("WTgen1_Heat", 4),
    rep("WTgen1_Combination", 4),
    
    rep("WTgen3_Control", 5),
    rep("WTgen3_Drought", 4),
    rep("WTgen3_CO2", 4),
    rep("WTgen3_Ozone", 4),
    rep("WTgen3_Heat", 4),
    rep("WTgen3_Combination", 4),
    
    rep("WTgen5_Control", 5),
    rep("WTgen5_Drought", 4),
    rep("WTgen5_CO2", 4),
    rep("WTgen5_Ozone", 4),
    rep("WTgen5_Heat", 4),
    rep("WTgen5_Combination", 4)
  ),
  levels = c(
    "WTgen1_Control", "WTgen1_Drought", "WTgen1_CO2",
    "WTgen1_Ozone", "WTgen1_Heat", "WTgen1_Combination",
    "WTgen3_Control", "WTgen3_Drought", "WTgen3_CO2",
    "WTgen3_Ozone", "WTgen3_Heat", "WTgen3_Combination",
    "WTgen5_Control", "WTgen5_Drought", "WTgen5_CO2",
    "WTgen5_Ozone", "WTgen5_Heat", "WTgen5_Combination"
  )
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# ============================================================
# 5. Exact four WT batches
# ============================================================

batch1_samples <- c(
  "WT1C", "WT6C", "WT8C",
  "WT3D", "WT4D",
  "WT3CO", "WT7CO",
  "WT5O", "WT7O",
  "WT1H", "WT5H",
  "WT2T", "WT3T"
)

batch2_samples <- c(
  "WT3C", "WT5C", "WT10C",
  "WT6D", "WT10D",
  "WT4CO", "WT5CO",
  "WT1O", "WT2O",
  "WT6H", "WT8H",
  "WT1T", "WT6T"
)

batch3_samples <- c(
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4"
)

batch4_samples <- c(
  "WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111",
  "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211",
  "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212",
  "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112",
  "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221",
  "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422"
)

coldata$batch <- factor(
  ifelse(
    rownames(coldata) %in% batch1_samples, "Batch1",
    ifelse(
      rownames(coldata) %in% batch2_samples, "Batch2",
      ifelse(
        rownames(coldata) %in% batch3_samples, "Batch3",
        ifelse(
          rownames(coldata) %in% batch4_samples, "Batch4",
          NA
        )
      )
    )
  ),
  levels = c("Batch1", "Batch2", "Batch3", "Batch4")
)

if (any(is.na(coldata$batch))) {
  stop(
    "These samples were not assigned to a batch:\n",
    paste(rownames(coldata)[is.na(coldata$batch)], collapse = ", ")
  )
}

# ============================================================
# 6. Group / generation / treatment exactly for plotting
# ============================================================

WT_gen1_samples <- c(
  "WT1C", "WT3C", "WT5C", "WT6C", "WT8C", "WT10C",
  "WT3D", "WT4D", "WT6D", "WT10D",
  "WT3CO", "WT4CO", "WT5CO", "WT7CO",
  "WT1O", "WT2O", "WT5O", "WT7O",
  "WT1H", "WT5H", "WT6H", "WT8H",
  "WT1T", "WT2T", "WT3T", "WT6T"
)

WT_gen3_samples <- batch3_samples
WT_gen5_samples <- batch4_samples

coldata$group <- factor(
  ifelse(
    rownames(coldata) %in% WT_gen1_samples, "WT Gen 1",
    ifelse(
      rownames(coldata) %in% WT_gen3_samples, "WT Gen 3",
      "WT Gen5"
    )
  ),
  levels = c("WT Gen 1", "WT Gen 3", "WT Gen5")
)

coldata$generation <- factor(
  ifelse(
    rownames(coldata) %in% WT_gen1_samples, "Gen1",
    ifelse(
      rownames(coldata) %in% WT_gen3_samples, "Gen3",
      "Gen5"
    )
  ),
  levels = c("Gen1", "Gen3", "Gen5")
)

control_samples <- c(
  "WT1C", "WT3C", "WT5C", "WT6C", "WT8C", "WT10C",
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111"
)

drought_samples <- c(
  "WT3D", "WT4D", "WT6D", "WT10D",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211"
)

CO2_samples <- c(
  "WT3CO", "WT4CO", "WT5CO", "WT7CO",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212"
)

ozone_samples <- c(
  "WT1O", "WT2O", "WT5O", "WT7O",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112"
)

heat_samples <- c(
  "WT1H", "WT5H", "WT6H", "WT8H",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221"
)

combination_samples <- c(
  "WT1T", "WT2T", "WT3T", "WT6T",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4",
  "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422"
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

cat("\nWT batch table:\n")
print(table(coldata$batch, coldata$group))

cat("\nWT treatment table:\n")
print(table(coldata$generation, coldata$treatment))

write.csv(
  coldata,
  file.path(out_dir, "WT_coldata_exact_old_input.csv")
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
# 8. Batch correction exactly like old PCA code
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
#saved as 10X 10 PDF
