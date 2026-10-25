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

saveRDS(
  obj,
  file.path(out_dir, "gsnor_vsd_bc_object_for_PCA.rds")
)

write.csv(
  assay(obj),
  file.path(out_dir, "gsnor_vst_batch_corrected_matrix_for_PCA.csv")
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
  
  p <- ggplot(
    df_pca,
    aes(x = .data[[xcol]], y = .data[[ycol]])
  ) +
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
      title = "gsnor-ko PCA",
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

ggsave(
  file.path(out_dir, "gsnor_PCA_PC1_PC2.pdf"),
  pca_out$plot,
  width = 8,
  height = 4
)

ggsave(
  file.path(out_dir, "gsnor_PCA_PC1_PC2.png"),
  pca_out$plot,
  width = 8,
  height = 4,
  dpi = 300
)

write.csv(
  pca_out$data,
  file.path(out_dir, "gsnor_PCA_PC1_PC2_coordinates.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    PC = paste0("PC", seq_along(pca_out$percentVar)),
    percent_variance = pca_out$percentVar
  ),
  file.path(out_dir, "gsnor_PCA_percent_variance.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = pca_out$selected_genes),
  file.path(out_dir, paste0("gsnor_PCA_top_", ntop, "_genes.csv")),
  row.names = FALSE
)

cat("\nDone. gsnor PCA results saved in:\n")
cat(out_dir, "\n")
