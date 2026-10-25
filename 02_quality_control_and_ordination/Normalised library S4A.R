rm(list = ls())

library(DESeq2)
library(limma)

count_file <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Count_matrix.txt"

df <- read.delim(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(df) <- make.unique(df[[1]])

# Use all WT and gsnor/gs samples in the count matrix.
# WT_Pf_2-1-4 is excluded.
sample_cols <- grep("^(WT|gs)[0-9]", colnames(df), value = TRUE)
sample_cols <- setdiff(sample_cols, "WT_Pf_2-1-4")

# Optional: remove these if they are not part of the final experiment
sample_cols <- setdiff(sample_cols, c("WT7O_new", "WT8C_new"))

sample_info <- data.frame(
  sample = sample_cols,
  stringsAsFactors = FALSE
)

sample_info$genotype <- ifelse(grepl("^WT", sample_info$sample), "WT", "gsnor-ko")

# G1 = original names without "_"
# G3 = names ending like _1, _2, _3, _4, _5
# G5 = duplicated/extra columns ending like _111, _212, _322, etc.
sample_info$generation <- ifelse(
  !grepl("_", sample_info$sample),
  "G1",
  ifelse(grepl("_[0-9]+$", sample_info$sample) &
           nchar(sub(".*_", "", sample_info$sample)) == 1,
         "G3",
         "G5")
)

sample_info$stress_code <- sub(
  "^(WT|gs)[0-9]+(CO|C|D|O|H|T).*",
  "\\2",
  sample_info$sample
)

sample_info$stress <- c(
  C  = "Control",
  D  = "Drought",
  CO = "CO2",
  O  = "Ozone",
  H  = "Heat",
  T  = "Combination"
)[sample_info$stress_code]

group_levels <- c(
  "WT G1", "gsnor-ko G1",
  "WT G3", "gsnor-ko G3",
  "WT G5", "gsnor-ko G5"
)

stress_levels <- c(
  "Control", "Drought", "CO2",
  "Ozone", "Heat", "Combination"
)

sample_info$group <- paste(sample_info$genotype, sample_info$generation)
sample_info$group <- factor(sample_info$group, levels = group_levels)
sample_info$stress <- factor(sample_info$stress, levels = stress_levels)

# Order inside each group:
# Control, Drought, CO2, Ozone, Heat, Combination
sample_info <- sample_info[order(
  sample_info$group,
  sample_info$stress,
  sample_info$sample
), ]

counts <- as.matrix(df[, sample_info$sample])
storage.mode(counts) <- "integer"

condition <- factor(
  paste(sample_info$group, sample_info$stress, sep = "_")
)

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
  "gs1C", "gs3C", "gs6C", "gs7C", "gs8C",
  "gs1D", "gs2D", "gs4D", "gs5D",
  "gs3CO", "gs4CO", "gs8CO", "gs9CO",
  "gs1O", "gs3O", "gs5O", "gs6O",
  "gs5H", "gs7H", "gs8H", "gs10H",
  "gs2T", "gs6T", "gs7T", "gs10T",
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4",
  "gs1C4_2", "gs3C4_1", "gs6C2_1", "gs6C2_3", "gs7C3_2",
  "gs1D3_1", "gs2D1_2", "gs4D3_2", "gs5D3_3",
  "gs3CO3_2", "gs4CO4_2", "gs8CO1_3", "gs9CO1_2",
  "gs1O2_2", "gs3O3_3", "gs5O4_2", "gs6O2_2",
  "gs5H3_2", "gs7H4_4", "gs8H1_2", "gs10H3_4",
  "gs2T4_2", "gs6T3_1", "gs7T2_4", "gs10T2_5"
)

# Batch 4: all generation 5 samples, both WT and gsnor-ko
batch4_samples <- sample_info$sample[sample_info$generation == "G5"]

batch <- ifelse(
  sample_info$sample %in% batch4_samples,
  "Batch4",
  ifelse(
    sample_info$sample %in% batch1_samples,
    "Batch1",
    ifelse(
      sample_info$sample %in% batch2_samples,
      "Batch2",
      ifelse(
        sample_info$sample %in% batch3_samples,
        "Batch3",
        NA
      )
    )
  )
)

if (any(is.na(batch))) {
  stop(
    "These samples were not assigned to a batch: ",
    paste(sample_info$sample[is.na(batch)], collapse = ", ")
  )
}

coldata <- data.frame(
  row.names = sample_info$sample,
  condition = condition,
  batch = factor(batch, levels = c("Batch1", "Batch2", "Batch3", "Batch4"))
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Same filtering logic from your old code
dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- estimateSizeFactors(dds)

idx <- rowSums(counts(dds, normalized = TRUE) >= 3) >= 2
dds <- dds[idx, ]

# VST normalization
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Same batch correction step from your old code
assay(vsd) <- limma::removeBatchEffect(
  assay(vsd),
  batch = vsd$batch
)

statusCol <- c(
  Control     = "#66C2A5",
  Drought     = "#FFD92F",
  CO2         = "#8DA0CB",
  Ozone       = "#E78AC3",
  Heat        = "#A6D854",
  Combination = "grey"
)

box_cols <- statusCol[as.character(sample_info$stress)]

setwd("~/Documents/DNA_seq/RNA_seq/figures_raw")

pdf("Figure_A_normalised_VST_distributions.pdf", width = 30, height = 15)

par(mar = c(12, 5, 4, 2) + 0.1)

boxplot(
  assay(vsd)[, sample_info$sample],
  xlab = "",
  ylab = "Log2 counts per million",
  las = 2,
  main = "Normalised Distributions",
  col = box_cols,
  border = "black",
  outline = TRUE,
  pch = 21,
  cex = 0.55,
  outcol = "black",
  outbg = "black",
  cex.axis = 1.0,
  cex.lab = 1.5
)

abline(
  h = median(assay(vsd)),
  col = "blue",
  lwd = 2
)

dev.off()
