rm(list = ls())

# ============================================================
# 0. Packages
# ============================================================

if (!requireNamespace("scattermore", quietly = TRUE)) {
  install.packages("scattermore")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(scattermore)
})

# ============================================================
# 1. Input settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "Volcano_plots_no_grid_LIGHT_PDF_no_Cairo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 0

padj_floor <- 1e-250
x_lim <- 10

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

treatment_order <- c(
  "Drought", "CO2", "Ozone", "Heat", "Combination"
)

genotype_order <- c("WT", "gsnor-ko")
generation_order <- c("1", "3", "5")

# ============================================================
# 2. Light PDF device
# No Cairo
# ============================================================

pdf_light <- function(filename, width, height, ...) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE,
    compress = TRUE,
    onefile = FALSE
  )
}

# ============================================================
# 3. Read DESeq2 result file
# ============================================================

read_deseq_file <- function(file) {
  
  dat <- tryCatch(
    read.table(
      file,
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE
    ),
    error = function(e1) {
      tryCatch(
        read.table(
          file,
          header = TRUE,
          row.names = 1,
          check.names = FALSE,
          stringsAsFactors = FALSE
        ),
        error = function(e2) {
          read.csv(
            file,
            header = TRUE,
            row.names = 1,
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }
      )
    }
  )
  
  dat <- as.data.frame(dat)
  
  if (!"gene_id" %in% colnames(dat)) {
    dat$gene_id <- rownames(dat)
  }
  
  dat
}

# ============================================================
# 4. Detect treatment from full DESeq2 files only
# ============================================================

get_treatment_exact <- function(file_name) {
  
  fname <- basename(file_name)
  fname_no_ext <- tools::file_path_sans_ext(fname)
  fname_low <- tolower(trimws(fname_no_ext))
  
  case_when(
    fname_low == "drought" ~ "Drought",
    fname_low == "co2" ~ "CO2",
    fname_low == "ozone" ~ "Ozone",
    fname_low == "heat" ~ "Heat",
    fname_low == "combination" ~ "Combination",
    TRUE ~ NA_character_
  )
}

# ============================================================
# 5. Read all full DESeq2 result files
# ============================================================

volcano_list <- list()
file_map <- list()

for (folder in folders) {
  
  folder_path <- file.path(root_dir, folder)
  
  if (!dir.exists(folder_path)) {
    warning("Folder not found: ", folder_path)
    next
  }
  
  all_files <- list.files(
    folder_path,
    pattern = "\\.(txt|tsv|csv)$",
    full.names = TRUE
  )
  
  files <- all_files[
    !str_detect(
      tolower(basename(all_files)),
      "pval|pvalue|pvales|upregulated|downregulated"
    )
  ]
  
  genotype <- ifelse(str_detect(folder, "^WT"), "WT", "gsnor-ko")
  
  generation <- str_extract(folder, "gen[0-9]+") %>%
    str_remove("gen")
  
  for (file in files) {
    
    treatment <- get_treatment_exact(file)
    
    file_map[[length(file_map) + 1]] <- data.frame(
      folder = folder,
      genotype = genotype,
      generation = generation,
      treatment_detected = treatment,
      file = basename(file),
      stringsAsFactors = FALSE
    )
    
    if (is.na(treatment)) {
      message("Skipping non-treatment file: ", basename(file))
      next
    }
    
    dat <- read_deseq_file(file)
    
    required_cols <- c("log2FoldChange", "padj")
    
    if (!all(required_cols %in% colnames(dat))) {
      warning("Skipping file because required columns are missing: ", basename(file))
      next
    }
    
    dat2 <- dat %>%
      mutate(
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj),
        genotype = genotype,
        generation = generation,
        treatment = treatment,
        file = basename(file)
      ) %>%
      select(
        gene_id,
        genotype,
        generation,
        treatment,
        file,
        log2FoldChange,
        padj
      )
    
    volcano_list[[length(volcano_list) + 1]] <- dat2
  }
}

volcano_df <- bind_rows(volcano_list)
file_map_df <- bind_rows(file_map)

write.csv(
  file_map_df,
  file.path(out_dir, "Volcano_file_detection_map.csv"),
  row.names = FALSE
)

if (nrow(volcano_df) == 0) {
  stop("No volcano data generated. Check folder path and DESeq2 result files.")
}

# ============================================================
# 6. Prepare volcano data
# ============================================================

volcano_df <- volcano_df %>%
  mutate(
    genotype = factor(genotype, levels = genotype_order),
    generation = factor(generation, levels = generation_order),
    treatment = factor(treatment, levels = treatment_order),
    
    neg_log10_padj = -log10(pmax(padj, padj_floor)),
    
    direction = case_when(
      !is.na(padj) &
        !is.na(log2FoldChange) &
        padj < padj_cutoff &
        log2FoldChange > lfc_cutoff ~ "Up",
      
      !is.na(padj) &
        !is.na(log2FoldChange) &
        padj < padj_cutoff &
        log2FoldChange < -lfc_cutoff ~ "Down",
      
      TRUE ~ "Not significant"
    ),
    
    direction = factor(
      direction,
      levels = c("Down", "Up", "Not significant")
    ),
    
    generation_label = case_when(
      generation == "1" ~ "Generation 1\n(Stress)",
      generation == "3" ~ "Generation 3\n(Stress)",
      generation == "5" ~ "Generation 5\n(No Stress)"
    ),
    
    generation_label = factor(
      generation_label,
      levels = c(
        "Generation 1\n(Stress)",
        "Generation 3\n(Stress)",
        "Generation 5\n(No Stress)"
      )
    )
  )

write.csv(
  volcano_df,
  file.path(out_dir, "Volcano_all_gene_data.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Count up/down genes for labels
# ============================================================

count_df <- volcano_df %>%
  group_by(genotype, generation, generation_label, treatment) %>%
  summarise(
    n_down = sum(direction == "Down", na.rm = TRUE),
    n_up = sum(direction == "Up", na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  count_df,
  file.path(out_dir, "Volcano_up_down_counts.csv"),
  row.names = FALSE
)

cat("\nUp/down counts for volcano labels:\n")
print(count_df)

# ============================================================
# 8. Volcano plotting function
# No Cairo, no ggrastr
# Points rasterized using scattermore
# Grid lines removed
# ============================================================

make_volcano_panel <- function(genotype_name, panel_title, show_legend = TRUE) {
  
  plot_data <- volcano_df %>%
    filter(genotype == genotype_name)
  
  label_data <- count_df %>%
    filter(genotype == genotype_name)
  
  label_down <- label_data %>%
    mutate(
      x = -x_lim * 0.78,
      y = Inf,
      label = n_down
    )
  
  label_up <- label_data %>%
    mutate(
      x = x_lim * 0.78,
      y = Inf,
      label = n_up
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      color = direction
    )
  ) +
    scattermore::geom_scattermore(
      pixels = c(700, 500),
      pointsize = 5,
      alpha = 0.65,
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "orange"
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "orange"
    ) +
    geom_text(
      data = label_down,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      color = "#1F4AFF",
      hjust = 0,
      vjust = 1.4,
      size = 5
    ) +
    geom_text(
      data = label_up,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      color = "#FF2B1A",
      hjust = 1,
      vjust = 1.4,
      size = 5
    ) +
    facet_grid(
      generation_label ~ treatment,
      scales = "free_y",
      switch = "y",
      drop = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Down" = "#1F4AFF",
        "Up" = "#FF2B1A",
        "Not significant" = "black"
      ),
      breaks = c("Down", "Up"),
      labels = c(
        "Downregulated\ngenes",
        "Upregulated genes"
      )
    ) +
    coord_cartesian(
      xlim = c(-x_lim, x_lim)
    ) +
    labs(
      title = panel_title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted p-value",
      color = NULL
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = 4,
          alpha = 1
        )
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0
      ),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = 13,
        color = "black"
      ),
      strip.text.y.left = element_text(
        size = 11,
        color = "#B00000",
        face = "bold",
        angle = 90
      ),
      strip.placement = "outside",
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        colour = "grey45",
        fill = NA,
        linewidth = 0.5
      ),
      
      axis.text = element_text(
        size = 7,
        color = "black"
      ),
      axis.title = element_text(
        size = 8,
        color = "black"
      ),
      
      legend.position = ifelse(show_legend, "top", "none"),
      legend.justification = "right",
      legend.text = element_text(size = 10),
      legend.margin = margin(0, 0, 0, 0),
      
      panel.spacing = unit(0.10, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  p
}

# ============================================================
# 9. Make WT and gsnor-ko volcano panels
# ============================================================

p_wt <- make_volcano_panel(
  genotype_name = "WT",
  panel_title = "(A) Volcano plots (WT)",
  show_legend = TRUE
)

p_gsnor <- make_volcano_panel(
  genotype_name = "gsnor-ko",
  panel_title = expression("(B) Volcano plots (" * italic(gsnor) * "-ko)"),
  show_legend = FALSE
)

print(p_wt)
print(p_gsnor)

# ============================================================
# 10. Combined figure
# ============================================================

combined_volcano <- p_wt / p_gsnor +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Supplementary Figure S14: RNA seq",
    theme = theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0,
        margin = margin(0, 0, 8, 0)
      )
    )
  )

print(combined_volcano)

# ============================================================
# 11. Save light PDFs and PNGs
# ============================================================

ggsave(
  file.path(out_dir, "Volcano_plots_WT_no_grid_LIGHT_no_Cairo.pdf"),
  p_wt,
  device = pdf_light,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_dir, "Volcano_plots_WT_no_grid_LIGHT_no_Cairo.png"),
  p_wt,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(out_dir, "Volcano_plots_gsnor_ko_no_grid_LIGHT_no_Cairo.pdf"),
  p_gsnor,
  device = pdf_light,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_dir, "Volcano_plots_gsnor_ko_no_grid_LIGHT_no_Cairo.png"),
  p_gsnor,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(out_dir, "Supplementary_Figure_S14_RNAseq_volcano_plots_no_grid_LIGHT_no_Cairo.pdf"),
  combined_volcano,
  device = pdf_light,
  width = 11,
  height = 14
)

ggsave(
  file.path(out_dir, "Supplementary_Figure_S14_RNAseq_volcano_plots_no_grid_LIGHT_no_Cairo.png"),
  combined_volcano,
  width = 11,
  height = 14,
  dpi = 300
)

cat("\nDone. Light volcano PDFs without Cairo saved in:\n")
cat(out_dir, "\n")

