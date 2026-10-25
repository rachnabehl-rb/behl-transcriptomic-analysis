############################################################
## Epigenetic regulator DEG heatmap
## Exact gene IDs provided manually
## No Excel matching, no org.At.tair.db, no automatic mapping
############################################################

rm(list = ls())

## =========================================================
## 0. Packages
## =========================================================

needed_cran <- c(
  "dplyr", "tidyr", "stringr", "ggplot2",
  "patchwork", "scales", "readr", "tibble"
)

for (pkg in needed_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
  library(tibble)
})

## =========================================================
## 1. Paths and settings
## =========================================================

ROOT_DIR <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(
  ROOT_DIR,
  "Epigenetic_regulator_heatmap_MANUAL_AGI_IDS"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

genotype_order <- c("WT", "gsnor-ko")

treatment_order <- c(
  "Drought", "CO2", "Ozone", "Heat", "Combination"
)

generation_order <- c("G1", "G3", "G5")

treatment_label_map <- c(
  "Drought"     = "Drought",
  "CO2"         = "CO2",
  "Ozone"       = "O3",
  "Heat"        = "Heat",
  "Combination" = "Combi"
)

## DEG threshold
padj_cutoff <- 0.05

## TRUE = colour only padj < 0.05 cells.
## FALSE = show raw log2FoldChange for all cells.
USE_ONLY_SIGNIFICANT_DEG <- TRUE

## Colour scale
lfc_limit <- 6

## =========================================================
## 2. Manual gene table: exact old-plot genes only
## =========================================================

target_genes <- tibble::tribble(
  ~panel, ~gene_label,    ~gene_id,     ~old_order,
  
  ## A: DNA methylation / RdDM / silencing
  "A", "RDR2",        "AT4G11130",  1,
  "A", "AGO6",        "AT2G32940",  2,
  "A", "DMS3",        "AT3G49250",  3,
  "A", "CMT3",        "AT1G69770",  4,
  "A", "MORC3",       "AT4G36270",  5,
  "A", "MORC6",       "AT1G19100",  6,
  "A", "ADCP1",       "AT1G09320",  7,
  "A", "DDM1",        "AT5G66750",  8,
  "A", "H2A.W.6",     "AT5G59870",  9,
  "A", "DRM2",        "AT5G14620", 10,
  "A", "RDM1",        "AT3G22680", 11,
  
  ## B: Polycomb / PcG-related
  "B", "BPC1",        "AT2G01930",  1,
  "B", "BPC2",        "AT2G35550",  2,
  "B", "BPC3",        "AT1G68120",  3,
  "B", "CLF",         "AT2G23380",  4,
  "B", "FIE",         "AT3G20740",  5,
  "B", "MSI1",        "AT5G58230",  6,
  
  ## C: histone / chromatin
  "C", "UBP22",       "AT5G10790",  1,
  "C", "DET1",        "AT4G10180",  2,
  "C", "HUB1",        "AT2G44950",  3,
  "C", "HUB2",        "AT1G55250",  4,
  "C", "HTB1",        "AT1G07790",  5,
  "C", "HTB2",        "AT5G22880",  6,
  "C", "HTB3",        "AT2G28720",  7,
  "C", "HTB4",        "AT5G59910",  8,
  "C", "HTB5",        "AT2G37470",  9,
  "C", "HTB6",        "AT3G53650", 10,
  "C", "HTB7",        "AT3G09480", 11,
  "C", "HTB8",        "AT1G08170", 12,
  "C", "HTB9",        "AT3G45980", 13,
  "C", "HTB11",       "AT3G46030", 14,
  
  ## D: RNA methylation
  "D", "ALKBH10B",    "AT4G02940",  1,
  "D", "FIP37",       "AT3G54170",  2,
  "D", "MTA",         "AT4G10760",  3,
  
  ## E: methyl cycle / one-carbon / SAM AT3G17390
  "E", "SAMS3",       "AT2G36880",  1,
  "E", "SAMS4",       "AT3G17390",  2,
  "E", "SAHH2",       "AT3G23810",  3,
  "E", "MTHFD1",      "AT3G12290",  4,
  "E", "FPGS2",       "AT3G10160",  5,
  "E", "MTHFR2",      "AT2G44160",  6,
  "E", "CGS1/MTO1",   "AT3G01120",  7
) %>%
  mutate(
    gene_id = toupper(gene_id)
  ) %>%
  arrange(
    factor(panel, levels = c("A", "B", "C", "D", "E")),
    old_order
  ) %>%
  group_by(panel) %>%
  mutate(
    gene_order_in_panel = row_number()
  ) %>%
  ungroup()

write_csv(
  target_genes,
  file.path(out_dir, "MANUAL_gene_IDs_used_for_heatmap.csv")
)

cat("\nManual genes used:\n")
print(target_genes)

## =========================================================
## 3. Helper functions
## =========================================================

clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- gsub("^transcript:", "", x, ignore.case = TRUE)
  x <- gsub("^gene:", "", x, ignore.case = TRUE)
  x <- sub("\\..*$", "", x)
  x <- trimws(x)
  x <- toupper(x)
  stringr::str_extract(x, "AT([1-5]|C|M)G\\d{5}")
}

is_agi_like <- function(x) {
  x <- toupper(as.character(x))
  grepl("^AT([1-5]|C|M)G\\d{5}", x)
}

get_treatment_exact <- function(file_name) {
  
  fname <- basename(file_name)
  fname_no_ext <- tools::file_path_sans_ext(fname)
  fname_low <- tolower(trimws(fname_no_ext))
  
  dplyr::case_when(
    fname_low == "drought"     ~ "Drought",
    fname_low == "co2"         ~ "CO2",
    fname_low == "ozone"       ~ "Ozone",
    fname_low == "heat"        ~ "Heat",
    fname_low == "combination" ~ "Combination",
    TRUE ~ NA_character_
  )
}

read_deseq_file_robust <- function(file) {
  
  ## Main case: DESeq2 table with first unnamed column as rowname/gene ID.
  dat <- tryCatch(
    read.table(
      file,
      header = TRUE,
      row.names = 1,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  
  if (!is.null(dat) &&
      all(c("log2FoldChange", "padj") %in% colnames(dat))) {
    
    dat <- as.data.frame(dat)
    dat$gene_id <- rownames(dat)
    return(dat)
  }
  
  ## Whitespace fallback
  dat <- tryCatch(
    read.table(
      file,
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  
  if (!is.null(dat) &&
      all(c("log2FoldChange", "padj") %in% colnames(dat))) {
    
    dat <- as.data.frame(dat)
    dat$gene_id <- rownames(dat)
    return(dat)
  }
  
  ## Repair shifted-header case
  dat <- tryCatch(
    read.table(
      file,
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat)) {
    dat <- tryCatch(
      read.csv(
        file,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )
  }
  
  if (is.null(dat)) {
    stop("Could not read file: ", file)
  }
  
  dat <- as.data.frame(dat)
  
  first_col_gene_fraction <- mean(is_agi_like(dat[[1]]), na.rm = TRUE)
  
  if (is.finite(first_col_gene_fraction) && first_col_gene_fraction > 0.5) {
    
    old_names <- colnames(dat)
    
    colnames(dat) <- make.unique(
      c("gene_id", old_names[seq_len(ncol(dat) - 1)])
    )
    
  } else {
    
    gene_candidates <- c("gene_id", "gene", "Gene", "AGI", "TAIR", "id", "ID")
    
    gene_col <- gene_candidates[gene_candidates %in% colnames(dat)][1]
    
    if (is.na(gene_col)) {
      gene_col <- colnames(dat)[1]
    }
    
    dat$gene_id <- dat[[gene_col]]
  }
  
  if (!all(c("log2FoldChange", "padj") %in% colnames(dat))) {
    stop(
      "log2FoldChange/padj missing in file:\n",
      file,
      "\nDetected columns:\n",
      paste(colnames(dat), collapse = ", ")
    )
  }
  
  dat
}

## =========================================================
## 4. Read DESeq2 result files
## =========================================================

all_res_list <- list()
file_map <- list()

for (folder in folders) {
  
  folder_path <- file.path(ROOT_DIR, folder)
  
  if (!dir.exists(folder_path)) {
    warning("Folder not found: ", folder_path)
    next
  }
  
  genotype <- ifelse(str_detect(folder, "^WT"), "WT", "gsnor-ko")
  
  generation_num <- str_extract(folder, "gen[0-9]+") %>%
    str_remove("gen")
  
  generation <- paste0("G", generation_num)
  
  all_files <- list.files(
    folder_path,
    pattern = "\\.(txt|csv|tsv)$",
    full.names = TRUE
  )
  
  files <- all_files[
    !str_detect(
      tolower(basename(all_files)),
      "pval|pvalue|pvales|upregulated|downregulated"
    )
  ]
  
  for (file in files) {
    
    treatment <- get_treatment_exact(file)
    
    file_map[[length(file_map) + 1]] <- data.frame(
      folder = folder,
      genotype = genotype,
      generation = generation,
      treatment = treatment,
      file = basename(file),
      stringsAsFactors = FALSE
    )
    
    if (is.na(treatment)) {
      next
    }
    
    dat <- read_deseq_file_robust(file)
    
    dat2 <- dat %>%
      mutate(
        gene_id = clean_gene_id(gene_id),
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj),
        genotype = genotype,
        generation = generation,
        treatment = treatment,
        significant = !is.na(padj) & padj < padj_cutoff
      ) %>%
      filter(!is.na(gene_id), gene_id != "") %>%
      dplyr::select(
        gene_id,
        genotype,
        generation,
        treatment,
        log2FoldChange,
        padj,
        significant
      )
    
    all_res_list[[length(all_res_list) + 1]] <- dat2
  }
}

all_res <- bind_rows(all_res_list)
file_map_df <- bind_rows(file_map)

write_csv(
  file_map_df,
  file.path(out_dir, "DESeq2_file_detection_map.csv")
)

if (nrow(all_res) == 0) {
  stop("No DESeq2 result data were read.")
}

## =========================================================
## 5. Check missing expected files
## =========================================================

expected_files <- expand_grid(
  genotype = genotype_order,
  treatment = treatment_order,
  generation = generation_order
)

available_files <- file_map_df %>%
  filter(!is.na(treatment)) %>%
  distinct(genotype, treatment, generation)

missing_expected_files <- expected_files %>%
  anti_join(
    available_files,
    by = c("genotype", "treatment", "generation")
  )

write_csv(
  missing_expected_files,
  file.path(out_dir, "MISSING_expected_genotype_treatment_generation_files.csv")
)

if (nrow(missing_expected_files) > 0) {
  warning(
    "Some expected files are missing. Check:\n",
    file.path(out_dir, "MISSING_expected_genotype_treatment_generation_files.csv")
  )
}

## =========================================================
## 6. X-axis layout
## =========================================================

x_lookup <- expand_grid(
  genotype = genotype_order,
  treatment = treatment_order,
  generation = generation_order
) %>%
  mutate(
    genotype = factor(genotype, levels = genotype_order),
    treatment = factor(treatment, levels = treatment_order),
    generation = factor(generation, levels = generation_order),
    
    treatment_gap = 0.35,
    genotype_gap = 1.30,
    
    treatment_i = as.integer(treatment),
    generation_i = as.integer(generation),
    
    wt_width =
      length(treatment_order) * length(generation_order) +
      (length(treatment_order) - 1) * treatment_gap,
    
    x_base = case_when(
      genotype == "WT" ~ 0,
      genotype == "gsnor-ko" ~ wt_width + genotype_gap
    ),
    
    x_pos =
      x_base +
      (treatment_i - 1) * (length(generation_order) + treatment_gap) +
      generation_i,
    
    x_label = paste(
      treatment_label_map[as.character(treatment)],
      generation
    )
  ) %>%
  dplyr::select(
    genotype,
    treatment,
    generation,
    x_pos,
    x_label
  )

write_csv(
  x_lookup,
  file.path(out_dir, "x_axis_layout_check.csv")
)

wt_range <- x_lookup %>%
  filter(genotype == "WT") %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    xmid = mean(range(x_pos))
  )

gs_range <- x_lookup %>%
  filter(genotype == "gsnor-ko") %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    xmid = mean(range(x_pos))
  )

separator_x <- mean(c(wt_range$xmax, gs_range$xmin))

## =========================================================
## 7. Prepare heatmap data
## =========================================================

plot_grid <- expand_grid(
  gene_id = target_genes$gene_id,
  genotype = genotype_order,
  treatment = treatment_order,
  generation = generation_order
) %>%
  left_join(
    target_genes,
    by = "gene_id"
  ) %>%
  left_join(
    x_lookup,
    by = c("genotype", "treatment", "generation")
  )

plot_data <- plot_grid %>%
  left_join(
    all_res,
    by = c("gene_id", "genotype", "treatment", "generation")
  ) %>%
  group_by(
    panel,
    gene_label,
    gene_id,
    old_order,
    gene_order_in_panel,
    genotype,
    treatment,
    generation,
    x_pos,
    x_label
  ) %>%
  summarise(
    log2FoldChange_raw = dplyr::first(log2FoldChange),
    padj = dplyr::first(padj),
    significant = dplyr::first(significant),
    .groups = "drop"
  ) %>%
  mutate(
    significant = ifelse(is.na(significant), FALSE, significant),
    
    plot_lfc = case_when(
      is.na(log2FoldChange_raw) ~ 0,
      USE_ONLY_SIGNIFICANT_DEG & !significant ~ 0,
      TRUE ~ log2FoldChange_raw
    ),
    
    plot_lfc = pmax(pmin(plot_lfc, lfc_limit), -lfc_limit)
  )

write_csv(
  plot_data,
  file.path(out_dir, "heatmap_plot_data_MANUAL_AGI_IDS.csv")
)

sanity_check <- plot_data %>%
  group_by(genotype, treatment, generation) %>%
  summarise(
    n_genes = n_distinct(gene_label),
    n_unique_raw_lfc = n_distinct(round(log2FoldChange_raw, 5), na.rm = TRUE),
    n_unique_plot_lfc = n_distinct(round(plot_lfc, 5), na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  sanity_check,
  file.path(out_dir, "SANITY_logFC_uniformity_check.csv")
)

missing_values <- plot_data %>%
  filter(is.na(log2FoldChange_raw))

write_csv(
  missing_values,
  file.path(out_dir, "gene_cells_with_no_DESeq2_value_found.csv")
)

## =========================================================
## 8. Plot function
## =========================================================

make_one_panel <- function(panel_letter,
                           show_x_axis = TRUE,
                           show_headers = FALSE,
                           show_legend = TRUE) {
  
  dat <- plot_data %>%
    filter(panel == panel_letter)
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  gene_levels <- target_genes %>%
    filter(panel == panel_letter) %>%
    arrange(old_order) %>%
    pull(gene_label)
  
  dat <- dat %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_levels))
    )
  
  min_x <- min(x_lookup$x_pos)
  max_x <- max(x_lookup$x_pos)
  
  p <- ggplot(
    dat,
    aes(
      x = x_pos,
      y = gene_label,
      fill = plot_lfc
    )
  ) +
    geom_tile(
      width = 0.95,
      height = 0.92,
      color = "black",
      linewidth = 0.25
    ) +
    geom_vline(
      xintercept = separator_x,
      color = "#B22222",
      linewidth = 0.55,
      linetype = "dashed"
    ) +
    scale_fill_gradient2(
      low = "#00008B",
      mid = "grey90",
      high = "#D7191C",
      midpoint = 0,
      limits = c(-lfc_limit, lfc_limit),
      oob = scales::squish,
      na.value = "grey90",
      name = expression(log[2]~"fold change")
    ) +
    scale_x_continuous(
      breaks = x_lookup$x_pos,
      labels = x_lookup$x_label,
      limits = c(min_x - 0.55, max_x + 0.85),
      expand = expansion(mult = c(0.002, 0.002))
    ) +
    scale_y_discrete(
      position = "right",
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      
      axis.text.y.right = element_text(
        color = "black",
        size = 9.8,
        hjust = 0
      ),
      axis.text.y.left = element_blank(),
      axis.ticks.y = element_blank(),
      
      axis.text.x = if (show_x_axis) {
        element_text(
          angle = 55,
          hjust = 1,
          vjust = 1,
          color = "black",
          size = 9
        )
      } else {
        element_blank()
      },
      
      axis.ticks.x = if (show_x_axis) {
        element_line(color = "black", linewidth = 0.25)
      } else {
        element_blank()
      },
      
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      
      plot.margin = margin(6, 95, 6, 28)
    )
  
  p <- p +
    annotate(
      "text",
      x = min_x - 1.55,
      y = Inf,
      label = paste0("(", panel_letter, ")"),
      hjust = 0,
      vjust = 1.15,
      size = 6,
      fontface = "bold"
    )
  
  if (show_headers) {
    p <- p +
      annotate(
        "text",
        x = wt_range$xmid,
        y = Inf,
        label = "WT",
        vjust = -0.65,
        size = 5,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = gs_range$xmid,
        y = Inf,
        label = "italic(gsnor)*'-ko'",
        parse = TRUE,
        vjust = -0.65,
        size = 5,
        fontface = "bold"
      )
  }
  
  p
}

## =========================================================
## 9. Combined figure
## =========================================================

panel_plots <- list()

for (pp in c("A", "B", "C", "D", "E")) {
  
  panel_plots[[pp]] <- make_one_panel(
    panel_letter = pp,
    show_x_axis = pp == "E",
    show_headers = pp == "A",
    show_legend = TRUE
  )
}

panel_plots <- panel_plots[!sapply(panel_plots, is.null)]

panel_heights <- target_genes %>%
  count(panel, name = "n_genes") %>%
  arrange(panel)

heights_use <- panel_heights$n_genes[
  panel_heights$panel %in% names(panel_plots)
]

combined_plot <- wrap_plots(panel_plots, ncol = 1) +
  plot_layout(
    heights = heights_use,
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

print(combined_plot)

fig_height <- max(7.5, min(16, 2.8 + 0.20 * nrow(target_genes)))

ggsave(
  file.path(out_dir, "Epigenetic_regulator_heatmap_MANUAL_AGI_IDS.pdf"),
  combined_plot,
  width = 13.8,
  height = fig_height,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(out_dir, "Epigenetic_regulator_heatmap_MANUAL_AGI_IDS.png"),
  combined_plot,
  width = 13.8,
  height = fig_height,
  dpi = 300,
  limitsize = FALSE
)

## =========================================================
## 10. Individual panels
## =========================================================

for (pp in names(panel_plots)) {
  
  p_ind <- make_one_panel(
    panel_letter = pp,
    show_x_axis = TRUE,
    show_headers = TRUE,
    show_legend = TRUE
  )
  
  n_genes_panel <- target_genes %>%
    filter(panel == pp) %>%
    nrow()
  
  ggsave(
    file.path(out_dir, paste0("Panel_", pp, "_MANUAL_AGI_IDS.pdf")),
    p_ind,
    width = 13.8,
    height = max(2.0, 1.2 + 0.20 * n_genes_panel),
    useDingbats = FALSE,
    limitsize = FALSE
  )
}

cat("\nDone.\n")
cat("Output saved in:\n")
cat(out_dir, "\n\n")

cat("Important check files:\n")
cat(file.path(out_dir, "MANUAL_gene_IDs_used_for_heatmap.csv"), "\n")
cat(file.path(out_dir, "DESeq2_file_detection_map.csv"), "\n")
cat(file.path(out_dir, "MISSING_expected_genotype_treatment_generation_files.csv"), "\n")
cat(file.path(out_dir, "SANITY_logFC_uniformity_check.csv"), "\n")
cat(file.path(out_dir, "heatmap_plot_data_MANUAL_AGI_IDS.csv"), "\n")
cat(file.path(out_dir, "gene_cells_with_no_DESeq2_value_found.csv"), "\n")

