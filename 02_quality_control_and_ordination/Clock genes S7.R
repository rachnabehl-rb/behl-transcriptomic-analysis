rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# 1. Input paths and settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"
count_file <- file.path(root_dir, "Count_matrix.txt")

out_dir <- file.path(root_dir, "Clock_gene_figure_VST_all_then_filter")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SAVE_PLOTS <- TRUE

padj_cutoff <- 0.05
lfc_cutoff  <- 0

# PCA settings
LABEL_PCA_POINTS <- FALSE
LOCK_PCA_X_AXIS <- TRUE
PCA_X_LIMITS <- c(-7, 5)

# ============================================================
# 2. Clock gene lists
# ============================================================

# Old PCA / delta-ME clock list, including PRR9 and LNK1
clock_gene_table_all <- tibble::tribble(
  ~gene_id,     ~symbol, ~phase,
  "AT2G46830",  "CCA1",  "Morning",
  "AT1G01060",  "LHY",   "Morning",
  "AT2G46790",  "PRR9",  "Morning",
  "AT5G02810",  "PRR7",  "Morning",
  "AT3G09600",  "RVE8",  "Morning",
  "AT5G02100",  "RVE4",  "Morning",
  "AT5G64170",  "LNK1",  "Morning",
  "AT3G54500",  "LNK2",  "Morning",
  
  "AT5G61380",  "PRR1",  "Evening",
  "AT5G24470",  "PRR5",  "Evening",
  "AT1G22770",  "GI",    "Evening",
  "AT2G25930",  "ELF3",  "Evening",
  "AT2G40080",  "ELF4",  "Evening",
  "AT3G46640",  "LUX",   "Evening",
  "AT5G57360",  "ZTL",   "Evening",
  "AT1G68050",  "FKF1",  "Evening",
  "AT5G08330",  "CHE",   "Evening"
)

# Display order requested for panels C and D
gene_order_symbols_display <- c(
  "CCA1", "LHY", "PRR7", "RVE8", "RVE4", "LNK2",
  "PRR1", "PRR5", "GI", "ELF3", "ELF4", "LUX", "ZTL", "FKF1", "CHE"
)

clock_gene_table_display <- clock_gene_table_all %>%
  filter(symbol %in% gene_order_symbols_display) %>%
  mutate(
    symbol = factor(symbol, levels = gene_order_symbols_display),
    phase = factor(phase, levels = c("Morning", "Evening")),
    gene_label = paste0(gene_id, " (", symbol, ")")
  ) %>%
  arrange(symbol)

clock_gene_table_all <- clock_gene_table_all %>%
  mutate(
    phase = factor(phase, levels = c("Morning", "Evening")),
    gene_label = paste0(gene_id, " (", symbol, ")")
  )

clock_genes_pca <- unique(clock_gene_table_all$gene_id)

# ============================================================
# 3. Read count matrix
# ============================================================

counts_df <- read.table(
  count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

bad_samples <- c(
  "WT_Pf_2-1-4",
  "WT7O_new",
  "WT8C_new"
)

sample_cols <- grep("^(WT|gs)[0-9]", colnames(counts_df), value = TRUE)
sample_cols <- setdiff(sample_cols, bad_samples)

# ============================================================
# 4. Metadata
# ============================================================

stress_map <- c(
  "C"  = "Control",
  "D"  = "Drought",
  "CO" = "CO2",
  "O"  = "Ozone",
  "H"  = "Heat",
  "T"  = "Combination"
)

treatment_order <- c(
  "Control", "Drought", "CO2", "Ozone", "Heat", "Combination"
)

generation_order <- c("G1", "G3", "G5")

sample_info <- data.frame(
  sample = sample_cols,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Genotype = ifelse(str_detect(sample, "^WT"), "WT", "gsnor-ko"),
    
    Generation = case_when(
      !str_detect(sample, "_") ~ "G1",
      nchar(sub(".*_", "", sample)) == 1 ~ "G3",
      TRUE ~ "G5"
    ),
    
    stress_code = str_match(sample, "^(WT|gs)[0-9]+(CO|C|D|O|H|T).*")[, 3],
    Treatment = stress_map[stress_code],
    Treatment_label = ifelse(Treatment == "Combination", "Combi", Treatment),
    condition = paste(Genotype, Generation, Treatment, sep = "_")
  ) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "gsnor-ko")),
    Generation = factor(Generation, levels = generation_order),
    Treatment = factor(Treatment, levels = treatment_order),
    condition = factor(condition)
  ) %>%
  arrange(Genotype, Generation, Treatment, sample)

if (any(is.na(sample_info$Treatment))) {
  stop(
    "Treatment could not be assigned for these samples:\n",
    paste(sample_info$sample[is.na(sample_info$Treatment)], collapse = ", ")
  )
}

write.csv(
  sample_info,
  file.path(out_dir, "sample_metadata_clock_figure.csv"),
  row.names = FALSE
)

counts <- as.matrix(counts_df[, sample_info$sample, drop = FALSE])
counts[is.na(counts)] <- 0
counts <- round(counts)
storage.mode(counts) <- "integer"

coldata <- data.frame(
  row.names = sample_info$sample,
  sample_info
)

# ============================================================
# 5. DESeq2 object and VST all genes first
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- estimateSizeFactors(dds)

idx <- rowSums(counts(dds, normalized = TRUE) >= 3) >= 2
dds <- dds[idx, ]

# This is the key requested step:
# normalise all genes first using VST, blind = FALSE
vst_obj <- vst(dds, blind = FALSE)

vst_mat_all <- assay(vst_obj)

cat("\nVST object created using all genes first.\n")
cat("VST matrix dimensions:\n")
print(dim(vst_mat_all))

# ============================================================
# 6. Match clock genes after VST
# ============================================================

match_clock_rows <- function(clock_table, mat) {
  
  rn_clean <- toupper(sub("\\..*$", "", rownames(mat)))
  
  clock_table %>%
    mutate(
      row_id = sapply(seq_len(n()), function(i) {
        
        gid <- toupper(gene_id[i])
        sym <- toupper(as.character(symbol[i]))
        
        idx_gene <- which(rn_clean == gid)
        idx_sym  <- which(rn_clean == sym)
        
        if (length(idx_gene) > 0) {
          rownames(mat)[idx_gene[1]]
        } else if (length(idx_sym) > 0) {
          rownames(mat)[idx_sym[1]]
        } else {
          NA_character_
        }
      })
    )
}

clock_all_matched <- match_clock_rows(clock_gene_table_all, vst_mat_all)

clock_display_matched <- match_clock_rows(clock_gene_table_display, vst_mat_all) %>%
  mutate(
    symbol = factor(symbol, levels = gene_order_symbols_display),
    gene_label = paste0(gene_id, " (", symbol, ")")
  ) %>%
  arrange(symbol)

clock_present_pca <- clock_all_matched %>%
  filter(!is.na(row_id))

clock_missing_pca <- clock_all_matched %>%
  filter(is.na(row_id))

clock_present_display <- clock_display_matched %>%
  filter(!is.na(row_id)) %>%
  arrange(symbol)

clock_missing_display <- clock_display_matched %>%
  filter(is.na(row_id))

cat("\nClock genes present for PCA:", nrow(clock_present_pca), "/", nrow(clock_gene_table_all), "\n")
print(clock_present_pca %>% select(gene_id, symbol, phase, row_id))

cat("\nClock genes missing for PCA:\n")
print(clock_missing_pca %>% select(gene_id, symbol, phase))

cat("\nClock genes present for C/D display:", nrow(clock_present_display), "/", nrow(clock_gene_table_display), "\n")
print(clock_present_display %>% select(gene_id, symbol, phase, row_id))

if (nrow(clock_present_pca) < 3) {
  stop("Too few clock genes present for PCA.")
}

# Matrix for B using old full clock list
clock_mat_all <- vst_mat_all[clock_present_pca$row_id, , drop = FALSE]
rownames(clock_mat_all) <- clock_present_pca$gene_label

# Matrix for C/D display order
clock_mat_display <- vst_mat_all[clock_present_display$row_id, , drop = FALSE]
rownames(clock_mat_display) <- clock_present_display$gene_label

gene_order_labels_display <- clock_present_display %>%
  arrange(symbol) %>%
  pull(gene_label)

# ============================================================
# 7. Common colours
# ============================================================

treatment_cols <- c(
  "Control" = "#00BA38",
  "Drought" = "#00BFC4",
  "CO2" = "#F8766D",
  "Ozone" = "#F564E3",
  "Heat" = "#619CFF",
  "Combination" = "#B79F00"
)

genotype_cols <- c(
  "WT" = "#F8766D",
  "gsnor-ko" = "#00BFC4"
)

# ============================================================
# A. PCA using DESeq2::plotPCA
# VST all genes first, then filter clock genes
# ============================================================

vst_clock <- vst_obj[clock_present_pca$row_id, ]

colData(vst_clock)$Genotype <- factor(
  colData(vst_clock)$Genotype,
  levels = c("WT", "gsnor-ko")
)

colData(vst_clock)$Treatment <- factor(
  colData(vst_clock)$Treatment,
  levels = treatment_order
)

colData(vst_clock)$Generation <- factor(
  colData(vst_clock)$Generation,
  levels = generation_order
)

pca_data <- plotPCA(
  vst_clock,
  intgroup = c("Genotype", "Treatment", "Generation"),
  #ntop = nrow(vst_clock),
  returnData = TRUE
)

percentVar <- round(100 * attr(pca_data, "percentVar"), 1)

cat("\nPanel A PCA variance:\n")
cat("PC1 =", percentVar[1], "%\n")
cat("PC2 =", percentVar[2], "%\n")
cat("PC1 range:\n")
print(range(pca_data$PC1, na.rm = TRUE))
cat("PC2 range:\n")
print(range(pca_data$PC2, na.rm = TRUE))

pA <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Treatment, shape = Genotype)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  facet_grid(. ~ Generation, scales = "fixed") +
  scale_color_manual(
    values = treatment_cols,
    labels = c(
      "Control" = "Control",
      "Drought" = "Drought",
      "CO2" = expression(CO[2]),
      "Ozone" = expression(O[3]),
      "Heat" = "Heat",
      "Combination" = "Combi"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "WT" = 17,
      "gsnor-ko" = 16
    )
  ) +
  labs(
    title = "(A)",
    x = paste0("PC1 (", percentVar[1], "% variance)"),
    y = paste0("PC2 (", percentVar[2], "% variance)"),
    color = "Treatments",
    shape = "Genotype"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

if (LOCK_PCA_X_AXIS) {
  pA <- pA + coord_cartesian(xlim = PCA_X_LIMITS)
}

if (LABEL_PCA_POINTS) {
  pA <- pA +
    ggrepel::geom_text_repel(
      aes(label = name),
      size = 2.5,
      max.overlaps = 100,
      show.legend = FALSE
    )
}

pA
##saved as 7X7 PDF
# ============================================================
# B. Delta ME using raw VST values
# Morning mean minus evening mean
# Uses the old full clock list including PRR9 and LNK1
# ============================================================

morning_genes <- clock_present_pca %>%
  filter(phase == "Morning") %>%
  pull(gene_label)

evening_genes <- clock_present_pca %>%
  filter(phase == "Evening") %>%
  pull(gene_label)

x_group_levels_B <- unlist(
  lapply(treatment_order, function(trt) {
    paste(trt, generation_order)
  })
)

x_group_labels_B <- gsub("Combination", "Combi", x_group_levels_B)

delta_df <- data.frame(
  sample = colnames(clock_mat_all),
  morning_score = colMeans(clock_mat_all[morning_genes, , drop = FALSE], na.rm = TRUE),
  evening_score = colMeans(clock_mat_all[evening_genes, , drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
) %>%
  mutate(
    delta_ME = morning_score - evening_score
  ) %>%
  left_join(
    sample_info %>%
      select(sample, Genotype, Treatment, Generation),
    by = "sample"
  ) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "gsnor-ko")),
    Generation = factor(Generation, levels = generation_order),
    Treatment = factor(Treatment, levels = treatment_order),
    x_group = paste(Treatment, Generation),
    x_group = factor(x_group, levels = x_group_levels_B)
  )

pB <- ggplot(
  delta_df,
  aes(x = x_group, y = delta_ME, color = Genotype)
) +
  geom_boxplot(
    fill = NA,
    outlier.shape = NA,
    width = 0.65,
    linewidth = 0.8,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.75
    ),
    size = 1.8,
    alpha = 0.85
  ) +
  scale_color_manual(values = genotype_cols) +
  scale_x_discrete(labels = x_group_labels_B) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(
    title = "(B)",
    x = NULL,
    y = expression(Delta*"ME (Morning - Evening)"),
    color = "Genotype"
  )

pB

# ============================================================
# C. Raw VST heatmap for control samples only
# Gene names/IDs on right
# Requested gene order
# ============================================================

control_samples <- sample_info %>%
  filter(Treatment == "Control") %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "gsnor-ko")),
    Generation = factor(Generation, levels = generation_order)
  ) %>%
  arrange(Generation, Genotype, sample)

control_mat <- clock_mat_display[, control_samples$sample, drop = FALSE]

control_heat_df <- as.data.frame(control_mat) %>%
  rownames_to_column("gene_label") %>%
  pivot_longer(
    cols = -gene_label,
    names_to = "sample",
    values_to = "vst_value"
  ) %>%
  left_join(
    clock_present_display %>%
      select(gene_label, phase, symbol, gene_id),
    by = "gene_label"
  ) %>%
  left_join(
    control_samples %>%
      select(sample, Genotype, Generation, Treatment),
    by = "sample"
  ) %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(gene_order_labels_display)),
    sample = factor(sample, levels = control_samples$sample),
    phase = factor(phase, levels = c("Morning", "Evening"))
  )

pC <- ggplot(
  control_heat_df,
  aes(x = sample, y = gene_label, fill = vst_value)
) +
  geom_tile(color = "grey70", linewidth = 0.25) +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#67A9CF", "#FFFFBF", "#F46D43", "#D73027"),
    name = "VST"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y.right = element_text(size = 8),
    axis.text.y.left = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "right"
  ) +
  labs(
    title = "(C)",
    x = NULL,
    y = NULL
  )

pC

# ============================================================
# D. pheatmap-style log2FC heatmap
# Same style as attached figure:
# one heatmap matrix, columns ordered as genotype-treatment blocks,
# each block has G1, G3, G5
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

# ------------------------------------------------------------
# 1. Define column order
# ------------------------------------------------------------

genotype_order_D <- c("WT", "gsnor-ko")
treatment_order_D <- c("Drought", "CO2", "Ozone", "Heat", "Combination")
generation_order <- c("G1", "G3", "G5")

comparison_order_D <- unlist(
  lapply(genotype_order_D, function(gt) {
    unlist(
      lapply(treatment_order_D, function(trt) {
        paste(gt, trt, generation_order)
      })
    )
  })
)

# Column labels should only show G1/G3/G5 repeatedly
column_labels_D <- rep(generation_order, times = length(genotype_order_D) * length(treatment_order_D))

# Gaps after every 3 columns
gaps_col_D <- seq(3, length(comparison_order_D) - 3, by = 3)

# ------------------------------------------------------------
# 2. Build full grid
# ------------------------------------------------------------

plot_grid_D <- tidyr::expand_grid(
  gene_label = gene_order_labels_display,
  genotype = genotype_order_D,
  treatment = treatment_order_D,
  generation = generation_order
) %>%
  mutate(
    comparison = paste(genotype, treatment, generation)
  )

gene_info <- clock_present_display %>%
  select(gene_label, gene_id, symbol, phase)

lfc_df_plot <- plot_grid_D %>%
  left_join(
    lfc_df_raw,
    by = c(
      "gene_label",
      "genotype",
      "treatment",
      "generation",
      "comparison"
    )
  ) %>%
  left_join(
    gene_info,
    by = "gene_label",
    suffix = c("", ".lookup")
  ) %>%
  mutate(
    gene_id = coalesce(gene_id, gene_id.lookup),
    symbol = coalesce(symbol, as.character(symbol.lookup)),
    phase = coalesce(phase, as.character(phase.lookup))
  ) %>%
  select(
    gene_id, symbol, gene_label, phase,
    genotype, treatment, generation,
    comparison, log2FoldChange, padj, log2FC_sig
  )

# ------------------------------------------------------------
# 3. Remove genes with no significant log2FC anywhere
# ------------------------------------------------------------

lfc_df_plot <- lfc_df_plot %>%
  group_by(gene_label) %>%
  filter(any(!is.na(log2FC_sig))) %>%
  ungroup()

gene_order_D_present <- gene_order_labels_display[
  gene_order_labels_display %in% unique(as.character(lfc_df_plot$gene_label))
]

cat("\nGenes retained in Panel D:\n")
print(gene_order_D_present)

# ------------------------------------------------------------
# 4. Convert long data to matrix
# ------------------------------------------------------------

lfc_df_plot <- lfc_df_plot %>%
  mutate(
    gene_label = factor(gene_label, levels = gene_order_D_present),
    comparison = factor(comparison, levels = comparison_order_D)
  ) %>%
  arrange(gene_label, comparison)

mat_D <- lfc_df_plot %>%
  select(gene_label, comparison, log2FC_sig) %>%
  pivot_wider(
    names_from = comparison,
    values_from = log2FC_sig
  ) %>%
  as.data.frame()

rownames(mat_D) <- mat_D$gene_label
mat_D$gene_label <- NULL

mat_D <- mat_D[gene_order_D_present, comparison_order_D, drop = FALSE]
mat_D <- as.matrix(mat_D)

# ------------------------------------------------------------
# 5. Optional row annotation for morning/evening genes
# ------------------------------------------------------------

row_anno_D <- lfc_df_plot %>%
  distinct(gene_label, phase) %>%
  mutate(gene_label = as.character(gene_label)) %>%
  filter(gene_label %in% rownames(mat_D)) %>%
  arrange(match(gene_label, rownames(mat_D))) %>%
  column_to_rownames("gene_label")

row_anno_D$phase <- factor(row_anno_D$phase, levels = c("Morning", "Evening"))

annotation_colors_D <- list(
  phase = c(
    "Morning" = "blue",
    "Evening" = "red"
  )
)

# ------------------------------------------------------------
# 6. Colour scale
# ------------------------------------------------------------

max_abs_D <- max(abs(mat_D), na.rm = TRUE)

if (!is.finite(max_abs_D)) {
  max_abs_D <- 1
}

breaks_D <- seq(-max_abs_D, max_abs_D, length.out = 101)

cols_D <- colorRampPalette(
  c("#313695", "white", "#A50026")
)(100)

# ------------------------------------------------------------
# 7. Draw pheatmap
# ------------------------------------------------------------

pD_heatmap <- pheatmap(
  mat_D,
  color = cols_D,
  breaks = breaks_D,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = column_labels_D,
  gaps_col = gaps_col_D,
  annotation_row = row_anno_D,
  annotation_colors = annotation_colors_D,
  border_color = "black",
  na_col = "grey92",
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 7,
  angle_col = 90,
  legend = TRUE,
  annotation_names_row = FALSE,
  main = "",
  silent = TRUE
)

# View in RStudio plot window
grid::grid.newpage()
grid::grid.draw(pD_heatmap$gtable)

# Save directly
pdf(
  file.path(out_dir, "Panel_D_pheatmap_exact_style.pdf"),
  width = 12,
  height = 3.5
)
grid::grid.draw(pD_heatmap$gtable)
dev.off()

png(
  file.path(out_dir, "Panel_D_pheatmap_exact_style.png"),
  width = 12,
  height = 3.5,
  units = "in",
  res = 300
)
grid::grid.draw(pD_heatmap$gtable)
dev.off()

# For combining with patchwork later
pD <- patchwork::wrap_elements(full = pD_heatmap$gtable)
# ============================================================
# Combined figure
# ============================================================

combined_clock_figure <- (pA | pB) / pC / pD +
  plot_layout(
    heights = c(1.15, 1, 1)
  )

combined_clock_figure

# ============================================================
# Save outputs
# ============================================================

if (SAVE_PLOTS) {
  
  ggsave(
    file.path(out_dir, "Panel_A_clock_gene_PCA_plotPCA_VST_all_then_filter.pdf"),
    pA,
    width = 7,
    height = 4
  )
  
  ggsave(
    file.path(out_dir, "Panel_A_clock_gene_PCA_plotPCA_VST_all_then_filter.png"),
    pA,
    width = 7,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "Panel_B_delta_ME_raw_VST.pdf"),
    pB,
    width = 7,
    height = 4
  )
  
  ggsave(
    file.path(out_dir, "Panel_B_delta_ME_raw_VST.png"),
    pB,
    width = 7,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "Panel_C_control_clock_gene_raw_VST_heatmap.pdf"),
    pC,
    width = 11,
    height = 4
  )
  
  ggsave(
    file.path(out_dir, "Panel_C_control_clock_gene_raw_VST_heatmap.png"),
    pC,
    width = 11,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "Panel_D_clock_gene_significant_log2FC_heatmap_filtered.pdf"),
    pD,
    width = 11,
    height = 4
  )
  
  ggsave(
    file.path(out_dir, "Panel_D_clock_gene_significant_log2FC_heatmap_filtered.png"),
    pD,
    width = 11,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "Combined_clock_gene_figure_final.pdf"),
    combined_clock_figure,
    width = 12,
    height = 12
  )
  
  ggsave(
    file.path(out_dir, "Combined_clock_gene_figure_final.png"),
    combined_clock_figure,
    width = 12,
    height = 12,
    dpi = 300
  )
  
  write.csv(
    pca_data,
    file.path(out_dir, "Panel_A_PCA_coordinates.csv"),
    row.names = FALSE
  )
  
  write.csv(
    delta_df,
    file.path(out_dir, "Panel_B_delta_ME_scores.csv"),
    row.names = FALSE
  )
  
  write.csv(
    control_heat_df,
    file.path(out_dir, "Panel_C_control_raw_VST_heatmap_data.csv"),
    row.names = FALSE
  )
  
  write.csv(
    lfc_df_plot,
    file.path(out_dir, "Panel_D_clock_gene_log2FC_heatmap_data_filtered.csv"),
    row.names = FALSE
  )
  
  write.csv(
    clock_present_pca,
    file.path(out_dir, "clock_genes_present_for_PCA.csv"),
    row.names = FALSE
  )
  
  write.csv(
    clock_missing_pca,
    file.path(out_dir, "clock_genes_missing_for_PCA.csv"),
    row.names = FALSE
  )
}

cat("\nDone. Final clock-gene figure saved in:\n")
cat(out_dir, "\n")