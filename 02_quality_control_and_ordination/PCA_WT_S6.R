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

out_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/WT_only_PCA_DAPC_exact_old_input"
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

saveRDS(obj, file.path(out_dir, "WT_vsd_bc_object.rds"))

write.csv(
  assay(obj),
  file.path(out_dir, "WT_vst_batch_corrected_matrix.csv")
)

# ============================================================
# 9. PCA PC1 vs PC2
# ============================================================

plot_pcs <- function(vsd_obj,
                     color_by = "treatment",
                     shape_by = "group",
                     ntop = 2000,
                     pcx = 1,
                     pcy = 2,
                     label_points = FALSE) {
  
  cd <- as.data.frame(colData(vsd_obj))
  
  mat <- assay(vsd_obj)
  
  rv <- genefilter::rowVars(mat)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca <- prcomp(
    t(mat[sel, ]),
    center = TRUE,
    scale. = FALSE
  )
  
  percentVar <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  df_pca <- as.data.frame(pca$x)
  df_pca$sample <- rownames(df_pca)
  df_pca <- cbind(df_pca, cd)
  
  xcol <- paste0("PC", pcx)
  ycol <- paste0("PC", pcy)
  
  p <- ggplot(df_pca, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(
      aes(
        color = .data[[color_by]],
        shape = .data[[shape_by]]
      ),
      size = 3,
      alpha = 0.9
    ) +
    facet_grid(. ~ group, scales = "fixed") +
    xlab(sprintf("PC%d: %.1f%% variance", pcx, percentVar[pcx])) +
    ylab(sprintf("PC%d: %.1f%% variance", pcy, percentVar[pcy])) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black")
    ) +
    labs(
      title = "WT PCA",
      color = "Treatment",
      shape = "Group"
    )
  
  if (label_points) {
    p <- p +
      geom_text(
        aes(label = sample),
        vjust = -0.7,
        size = 3,
        show.legend = FALSE
      )
  }
  
  list(
    plot = p,
    data = df_pca,
    pca = pca,
    percentVar = percentVar,
    selected_genes = rownames(mat)[sel]
  )
}

pca_out <- plot_pcs(
  obj,
  color_by = "treatment",
  shape_by = "group",
  ntop = ntop,
  pcx = 1,
  pcy = 2
)

pca_out$plot

##saved as 5X7 for publication

ggsave(
  file.path(out_dir, "WT_PCA_PC1_PC2.pdf"),
  pca_out$plot,
  width = 8,
  height = 4
)

ggsave(
  file.path(out_dir, "WT_PCA_PC1_PC2.png"),
  pca_out$plot,
  width = 8,
  height = 4,
  dpi = 300
)

write.csv(
  pca_out$data,
  file.path(out_dir, "WT_PCA_PC1_PC2_coordinates.csv"),
  row.names = FALSE
)

