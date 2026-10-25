rm(list = ls())

# ============================================================
# 0. Packages
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
}

if (!requireNamespace("org.At.tair.db", quietly = TRUE)) {
  BiocManager::install("org.At.tair.db", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(AnnotationDbi)
  library(org.At.tair.db)
})

# ============================================================
# 1. Input settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "Persistent_Switcher_heatmaps_TAIR_annotated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 0

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

treatment_order <- c(
  "Drought", "CO2", "Ozone", "Heat", "Combination"
)

treatment_labels <- c(
  "Drought"     = "Drought",
  "CO2"         = "CO\u2082",
  "Ozone"       = "O\u2083",
  "Heat"        = "Heat",
  "Combination" = "Combi"
)

generation_order <- c("G1", "G3", "G5")

low_col  <- "#00008B"
mid_col  <- "grey94"
high_col <- "#D7191C"

fill_limit <- 10

# ============================================================
# 2. Exact plotting rules
# ============================================================
# required_n = number of significant generations required.
#
# Corrected rule:
# WT persistent Ozone = 2 generations
#
# gsnor-ko Drought = 2 generations and G5 forced to None
# because there are no significant genes in G5 drought.

rule_table <- tibble::tribble(
  ~genotype,   ~class,        ~treatment,     ~required_n, ~force_G5_None,
  
  # WT persistent
  "WT",        "Persistent",  "Drought",      3,           FALSE,
  "WT",        "Persistent",  "CO2",          2,           FALSE,
  "WT",        "Persistent",  "Ozone",        2,           FALSE,
  "WT",        "Persistent",  "Heat",         3,           FALSE,
  "WT",        "Persistent",  "Combination",  3,           FALSE,
  
  # WT switcher
  "WT",        "Switcher",    "Drought",      3,           FALSE,
  "WT",        "Switcher",    "CO2",          2,           FALSE,
  "WT",        "Switcher",    "Ozone",        2,           FALSE,
  "WT",        "Switcher",    "Heat",         3,           FALSE,
  "WT",        "Switcher",    "Combination",  3,           FALSE,
  
  # gsnor-ko persistent
  "gsnor-ko",  "Persistent",  "Drought",      2,           TRUE,
  "gsnor-ko",  "Persistent",  "CO2",          3,           FALSE,
  "gsnor-ko",  "Persistent",  "Ozone",        3,           FALSE,
  "gsnor-ko",  "Persistent",  "Heat",         3,           FALSE,
  "gsnor-ko",  "Persistent",  "Combination",  3,           FALSE,
  
  # gsnor-ko switcher
  "gsnor-ko",  "Switcher",    "Drought",      2,           TRUE,
  "gsnor-ko",  "Switcher",    "CO2",          3,           FALSE,
  "gsnor-ko",  "Switcher",    "Ozone",        3,           FALSE,
  "gsnor-ko",  "Switcher",    "Heat",         3,           FALSE,
  "gsnor-ko",  "Switcher",    "Combination",  3,           FALSE
)

# ============================================================
# 3. Helper functions
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
  dat$gene_id <- rownames(dat)
  dat
}

get_treatment_exact <- function(file_name) {
  
  fname <- basename(file_name)
  fname_no_ext <- tools::file_path_sans_ext(fname)
  fname_low <- tolower(trimws(fname_no_ext))
  
  case_when(
    fname_low == "drought"     ~ "Drought",
    fname_low == "co2"         ~ "CO2",
    fname_low == "ozone"       ~ "Ozone",
    fname_low == "heat"        ~ "Heat",
    fname_low == "combination" ~ "Combination",
    TRUE ~ NA_character_
  )
}

get_symbol_column <- function(dat) {
  
  cand <- c(
    "symbol", "SYMBOL",
    "gene_symbol", "GeneSymbol",
    "gene_name", "GeneName",
    "name", "Name",
    "external_gene_name",
    "Description", "description"
  )
  
  hit <- cand[cand %in% colnames(dat)]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}

classify_gene <- function(s1, s3, s5) {
  
  sig_status <- c(s1, s3, s5)
  sig_status <- sig_status[sig_status != "None"]
  
  if (length(sig_status) < 2) {
    return(NA_character_)
  }
  
  if (length(unique(sig_status)) == 1) {
    return("Persistent")
  }
  
  if (all(c("Up", "Down") %in% sig_status)) {
    return("Switcher")
  }
  
  NA_character_
}

# ============================================================
# 4. TAIR annotation using org.At.tair.db
# ============================================================

get_tair_annotation_from_bioc <- function(gene_ids) {
  
  gene_ids <- unique(toupper(sub("\\..*$", "", gene_ids)))
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  
  available_keytypes <- AnnotationDbi::keytypes(org.At.tair.db)
  available_columns  <- AnnotationDbi::columns(org.At.tair.db)
  
  keytype_use <- if ("TAIR" %in% available_keytypes) {
    "TAIR"
  } else if ("ALIAS" %in% available_keytypes) {
    "ALIAS"
  } else {
    stop("Could not find TAIR or ALIAS keytype in org.At.tair.db.")
  }
  
  columns_use <- intersect(
    c("TAIR", "SYMBOL", "GENENAME"),
    available_columns
  )
  
  ann <- AnnotationDbi::select(
    org.At.tair.db,
    keys = gene_ids,
    keytype = keytype_use,
    columns = columns_use
  )
  
  ann <- as.data.frame(ann)
  
  if ("TAIR" %in% colnames(ann)) {
    ann$gene_id <- ann$TAIR
  } else {
    ann$gene_id <- ann[[keytype_use]]
  }
  
  if (!"SYMBOL" %in% colnames(ann)) {
    ann$SYMBOL <- NA_character_
  }
  
  if (!"GENENAME" %in% colnames(ann)) {
    ann$GENENAME <- NA_character_
  }
  
  ann_out <- ann %>%
    mutate(
      gene_id = toupper(sub("\\..*$", "", gene_id)),
      tair_symbol = as.character(SYMBOL),
      tair_name = as.character(GENENAME)
    ) %>%
    group_by(gene_id) %>%
    summarise(
      tair_symbol = {
        x <- tair_symbol[!is.na(tair_symbol) & tair_symbol != ""]
        if (length(x) > 0) x[1] else NA_character_
      },
      tair_name = {
        x <- tair_name[!is.na(tair_name) & tair_name != ""]
        if (length(x) > 0) x[1] else NA_character_
      },
      .groups = "drop"
    )
  
  ann_out
}

make_gene_label <- function(gene_id,
                            tair_symbol,
                            tair_name,
                            file_symbol = NA_character_) {
  
  symbol <- dplyr::coalesce(tair_symbol, file_symbol)
  symbol <- ifelse(is.na(symbol), "", symbol)
  symbol <- trimws(symbol)
  
  tair_name <- ifelse(is.na(tair_name), "", tair_name)
  tair_name <- trimws(tair_name)
  
  label <- ifelse(
    symbol != "" & symbol != gene_id,
    paste0(gene_id, " (", symbol, ")"),
    ifelse(
      tair_name != "" & tair_name != gene_id,
      paste0(gene_id, " (", tair_name, ")"),
      gene_id
    )
  )
  
  label
}

# ============================================================
# 5. Read all full DESeq2 result files
# ============================================================

all_res_list <- list()
file_map <- list()

for (folder in folders) {
  
  folder_path <- file.path(root_dir, folder)
  
  if (!dir.exists(folder_path)) {
    warning("Folder not found: ", folder_path)
    next
  }
  
  genotype <- ifelse(grepl("^WT", folder), "WT", "gsnor-ko")
  
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
    
    if (is.na(treatment)) next
    
    dat <- read_deseq_file(file)
    
    if (!all(c("log2FoldChange", "padj") %in% colnames(dat))) {
      warning("Skipping file because log2FoldChange/padj missing: ", file)
      next
    }
    
    sym_col <- get_symbol_column(dat)
    
    if (!is.na(sym_col)) {
      dat$file_symbol <- as.character(dat[[sym_col]])
    } else {
      dat$file_symbol <- NA_character_
    }
    
    dat2 <- dat %>%
      mutate(
        gene_id = toupper(sub("\\..*$", "", gene_id)),
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj),
        genotype = genotype,
        generation = generation,
        treatment = treatment,
        status = case_when(
          !is.na(padj) &
            !is.na(log2FoldChange) &
            padj < padj_cutoff &
            log2FoldChange > lfc_cutoff ~ "Up",
          
          !is.na(padj) &
            !is.na(log2FoldChange) &
            padj < padj_cutoff &
            log2FoldChange < -lfc_cutoff ~ "Down",
          
          TRUE ~ "None"
        )
      ) %>%
      dplyr::select(
        gene_id,
        file_symbol,
        genotype,
        generation,
        treatment,
        log2FoldChange,
        padj,
        status
      )
    
    all_res_list[[length(all_res_list) + 1]] <- dat2
  }
}

all_res <- bind_rows(all_res_list)
file_map_df <- bind_rows(file_map)

write.csv(
  file_map_df,
  file.path(out_dir, "heatmap_file_detection_map.csv"),
  row.names = FALSE
)

if (nrow(all_res) == 0) {
  stop("No DESeq2 result data were read.")
}

# ============================================================
# 6. Get TAIR annotation
# ============================================================

tair_ann <- get_tair_annotation_from_bioc(all_res$gene_id)

write.csv(
  tair_ann,
  file.path(out_dir, "TAIR_annotation_from_org_At_tair_db.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Collapse possible duplicates and make wide table
# ============================================================

all_res_clean <- all_res %>%
  group_by(gene_id, genotype, treatment, generation) %>%
  summarise(
    file_symbol = {
      x <- file_symbol[!is.na(file_symbol) & file_symbol != ""]
      if (length(x) > 0) x[1] else NA_character_
    },
    log2FoldChange = dplyr::first(log2FoldChange),
    padj = dplyr::first(padj),
    status = dplyr::first(status),
    .groups = "drop"
  )

wide_lfc <- all_res_clean %>%
  dplyr::select(gene_id, file_symbol, genotype, treatment, generation, log2FoldChange) %>%
  pivot_wider(
    names_from = generation,
    values_from = log2FoldChange,
    names_prefix = "lfc_"
  )

wide_status <- all_res_clean %>%
  dplyr::select(gene_id, genotype, treatment, generation, status) %>%
  pivot_wider(
    names_from = generation,
    values_from = status,
    names_prefix = "status_",
    values_fill = "None"
  )

wide_padj <- all_res_clean %>%
  dplyr::select(gene_id, genotype, treatment, generation, padj) %>%
  pivot_wider(
    names_from = generation,
    values_from = padj,
    names_prefix = "padj_"
  )

wide_symbol <- all_res_clean %>%
  dplyr::select(gene_id, file_symbol) %>%
  filter(!is.na(file_symbol), file_symbol != "") %>%
  distinct() %>%
  group_by(gene_id) %>%
  summarise(
    file_symbol_final = dplyr::first(file_symbol),
    .groups = "drop"
  )

gene_wide <- wide_lfc %>%
  left_join(wide_status, by = c("gene_id", "genotype", "treatment")) %>%
  left_join(wide_padj, by = c("gene_id", "genotype", "treatment")) %>%
  left_join(wide_symbol, by = "gene_id") %>%
  left_join(tair_ann, by = "gene_id") %>%
  mutate(
    file_symbol = coalesce(file_symbol, file_symbol_final),
    gene_label = make_gene_label(
      gene_id = gene_id,
      tair_symbol = tair_symbol,
      tair_name = tair_name,
      file_symbol = file_symbol
    )
  ) %>%
  dplyr::select(-file_symbol_final)

needed_cols <- c(
  "lfc_G1", "lfc_G3", "lfc_G5",
  "status_G1", "status_G3", "status_G5",
  "padj_G1", "padj_G3", "padj_G5"
)

for (cc in needed_cols) {
  if (!cc %in% colnames(gene_wide)) {
    if (str_detect(cc, "^status")) {
      gene_wide[[cc]] <- "None"
    } else {
      gene_wide[[cc]] <- NA_real_
    }
  }
}

# ============================================================
# 8. Classify persistent/switcher and apply corrected rules
# ============================================================

# ============================================================
# 8. Classify persistent/switcher and apply corrected rules
# FIXED: remove None/None/None genes before max/mean calculation
# ============================================================

gene_classified <- gene_wide %>%
  rowwise() %>%
  mutate(
    n_sig = sum(c(status_G1, status_G3, status_G5) != "None"),
    
    class = classify_gene(status_G1, status_G3, status_G5),
    
    sig_pattern = paste(status_G1, status_G3, status_G5, sep = "_"),
    
    sig_status_string = paste(
      c(status_G1, status_G3, status_G5)[c(status_G1, status_G3, status_G5) != "None"],
      collapse = "_"
    ),
    
    dominant_direction = case_when(
      class == "Persistent" &
        sig_status_string != "" &
        length(unique(strsplit(sig_status_string, "_")[[1]])) == 1 &
        strsplit(sig_status_string, "_")[[1]][1] == "Up" ~ "Up",
      
      class == "Persistent" &
        sig_status_string != "" &
        length(unique(strsplit(sig_status_string, "_")[[1]])) == 1 &
        strsplit(sig_status_string, "_")[[1]][1] == "Down" ~ "Down",
      
      TRUE ~ NA_character_
    ),
    
    # Use log2FC only where the gene is significant.
    # Non-significant generations become NA and plot white.
    plot_lfc_G1 = ifelse(status_G1 == "None", NA_real_, lfc_G1),
    plot_lfc_G3 = ifelse(status_G3 == "None", NA_real_, lfc_G3),
    plot_lfc_G5 = ifelse(status_G5 == "None", NA_real_, lfc_G5)
  ) %>%
  ungroup() %>%
  
  # Remove genes that are not persistent/switcher before calculating max/mean
  filter(
    !is.na(class),
    n_sig >= 2
  ) %>%
  
  rowwise() %>%
  mutate(
    mean_lfc_sig = mean(
      c(plot_lfc_G1, plot_lfc_G3, plot_lfc_G5),
      na.rm = TRUE
    ),
    
    max_abs_lfc_sig = {
      x <- abs(c(plot_lfc_G1, plot_lfc_G3, plot_lfc_G5))
      x <- x[!is.na(x)]
      if (length(x) == 0) NA_real_ else max(x)
    }
  ) %>%
  ungroup() %>%
  
  left_join(
    rule_table,
    by = c("genotype", "class", "treatment")
  ) %>%
  filter(
    !is.na(required_n),
    n_sig == required_n
  ) %>%
  filter(
    ifelse(force_G5_None, status_G5 == "None", TRUE)
  ) %>%
  mutate(
    treatment = factor(treatment, levels = treatment_order),
    class = factor(class, levels = c("Persistent", "Switcher"))
  )

write.csv(
  gene_classified,
  file.path(out_dir, "persistent_switcher_gene_table_TAIR_corrected.csv"),
  row.names = FALSE
)

summary_counts <- gene_classified %>%
  count(genotype, treatment, class, required_n, force_G5_None, name = "n_genes") %>%
  arrange(genotype, class, treatment)

write.csv(
  summary_counts,
  file.path(out_dir, "persistent_switcher_counts_TAIR_corrected.csv"),
  row.names = FALSE
)

cat("\nFinal counts after corrected rules:\n")
print(summary_counts)

# ============================================================
# 9. Prepare heatmap data
# ============================================================

prepare_heatmap_data <- function(genotype_name, class_name) {
  
  dat <- gene_classified %>%
    filter(
      genotype == genotype_name,
      class == class_name
    ) %>%
    mutate(
      treatment = factor(treatment, levels = treatment_order),
      treatment_label = treatment_labels[as.character(treatment)]
    )
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  if (class_name == "Persistent") {
    dat <- dat %>%
      arrange(
        treatment,
        desc(required_n),
        dominant_direction,
        desc(abs(mean_lfc_sig)),
        desc(max_abs_lfc_sig),
        gene_label
      )
  } else {
    dat <- dat %>%
      arrange(
        treatment,
        desc(required_n),
        sig_pattern,
        desc(max_abs_lfc_sig),
        gene_label
      )
  }
  
  dat_long <- dat %>%
    dplyr::select(
      gene_id,
      gene_label,
      genotype,
      treatment,
      treatment_label,
      class,
      required_n,
      n_sig,
      sig_pattern,
      status_G1,
      status_G3,
      status_G5,
      plot_lfc_G1,
      plot_lfc_G3,
      plot_lfc_G5
    ) %>%
    pivot_longer(
      cols = starts_with("plot_lfc_"),
      names_to = "Generation",
      values_to = "log2FC"
    ) %>%
    mutate(
      Generation = recode(
        Generation,
        "plot_lfc_G1" = "G1",
        "plot_lfc_G3" = "G3",
        "plot_lfc_G5" = "G5"
      ),
      Generation = factor(Generation, levels = generation_order)
    )
  
  dat_long
}

# ============================================================
# 10. Make one treatment heatmap
# ============================================================

make_single_treatment_plot <- function(dat_long_treat,
                                       class_name,
                                       show_legend = FALSE) {
  
  trt_label <- unique(dat_long_treat$treatment_label)
  req_n <- unique(dat_long_treat$required_n)
  
  dat_gene_order <- dat_long_treat %>%
    distinct(gene_label, sig_pattern, required_n)
  
  gene_levels <- rev(dat_gene_order$gene_label)
  
  n_genes <- length(gene_levels)
  
  y_font <- case_when(
    n_genes <= 20  ~ 6.5,
    n_genes <= 60  ~ 5.2,
    n_genes <= 120 ~ 4.0,
    TRUE ~ 2.4
  )
  
  p <- ggplot(
    dat_long_treat,
    aes(
      x = Generation,
      y = factor(gene_label, levels = gene_levels),
      fill = log2FC
    )
  ) +
    geom_tile(
      color = "grey25",
      linewidth = 0.22
    ) +
    scale_fill_gradient2(
      low = low_col,
      mid = mid_col,
      high = high_col,
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      oob = scales::squish,
      na.value = "white",
      name = "log2FC"
    ) +
    scale_y_discrete(position = "right") +
    labs(
      title = paste0(trt_label, " (", tolower(class_name), ", ", req_n, " gen)"),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0
      ),
      
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      
      axis.text.x = element_text(
        color = "black",
        size = 8
      ),
      axis.text.y.right = element_text(
        color = "black",
        size = y_font
      ),
      axis.text.y.left = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(3, 3, 3, 3)
    )
  
  p
}

# ============================================================
# 11. Make full figure for one genotype/class
# ============================================================

make_full_figure <- function(genotype_name, class_name) {
  
  dat_long <- prepare_heatmap_data(genotype_name, class_name)
  
  if (is.null(dat_long) || nrow(dat_long) == 0) {
    message("No data for ", genotype_name, " ", class_name)
    return(NULL)
  }
  
  trt_counts <- dat_long %>%
    distinct(treatment, treatment_label, gene_label) %>%
    count(treatment, treatment_label, name = "n_genes") %>%
    mutate(
      treatment = factor(treatment, levels = treatment_order)
    ) %>%
    arrange(desc(n_genes))
  
  largest_treatment <- as.character(trt_counts$treatment[1])
  
  other_treatments <- dat_long %>%
    distinct(treatment) %>%
    mutate(treatment = factor(treatment, levels = treatment_order)) %>%
    arrange(treatment) %>%
    filter(treatment != largest_treatment) %>%
    pull(treatment) %>%
    as.character()
  
  dat_big <- dat_long %>%
    filter(treatment == largest_treatment)
  
  p_big <- make_single_treatment_plot(
    dat_big,
    class_name = class_name,
    show_legend = TRUE
  )
  
  if (length(other_treatments) > 0) {
    
    left_plots <- lapply(other_treatments, function(trt) {
      
      dat_t <- dat_long %>%
        filter(treatment == trt)
      
      make_single_treatment_plot(
        dat_t,
        class_name = class_name,
        show_legend = FALSE
      )
    })
    
    left_heights <- sapply(other_treatments, function(trt) {
      dat_long %>%
        filter(treatment == trt) %>%
        distinct(gene_label) %>%
        nrow()
    })
    
    left_heights <- pmax(left_heights, 8)
    
    left_combined <- wrap_plots(
      left_plots,
      ncol = 1
    ) +
      plot_layout(heights = left_heights)
    
    final_plot <- left_combined | p_big
    
    final_plot <- final_plot +
      plot_layout(
        widths = c(1, 1.45),
        guides = "collect"
      )
    
  } else {
    
    final_plot <- p_big
  }
  
  final_plot <- final_plot &
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    )
  
  final_plot
}

# ============================================================
# 12. Make four figures
# ============================================================

p_wt_persistent <- make_full_figure("WT", "Persistent")
p_wt_switcher   <- make_full_figure("WT", "Switcher")
p_gs_persistent <- make_full_figure("gsnor-ko", "Persistent")
p_gs_switcher   <- make_full_figure("gsnor-ko", "Switcher")

if (!is.null(p_wt_persistent)) print(p_wt_persistent)
if (!is.null(p_wt_switcher))   print(p_wt_switcher)
if (!is.null(p_gs_persistent)) print(p_gs_persistent)
if (!is.null(p_gs_switcher))   print(p_gs_switcher)

# ============================================================
# 13. Save four figures
# ============================================================

save_heatmap_figure <- function(plot_obj, genotype_name, class_name, file_stub) {
  
  if (is.null(plot_obj)) return(NULL)
  
  dat_long <- prepare_heatmap_data(genotype_name, class_name)
  
  n_genes_total <- dat_long %>%
    distinct(treatment, gene_label) %>%
    nrow()
  
  fig_height <- max(7, min(22, 3.5 + 0.08 * n_genes_total))
  
  ggsave(
    file.path(out_dir, paste0(file_stub, ".pdf")),
    plot_obj,
    width = 9.5,
    height = fig_height,
    useDingbats = FALSE
  )
  
  ggsave(
    file.path(out_dir, paste0(file_stub, ".png")),
    plot_obj,
    width = 9.5,
    height = fig_height,
    dpi = 300
  )
}

save_heatmap_figure(
  p_wt_persistent,
  "WT",
  "Persistent",
  "WT_persistent_heatmap_TAIR_corrected"
)

save_heatmap_figure(
  p_wt_switcher,
  "WT",
  "Switcher",
  "WT_switcher_heatmap_TAIR_corrected"
)

save_heatmap_figure(
  p_gs_persistent,
  "gsnor-ko",
  "Persistent",
  "gsnor_ko_persistent_heatmap_TAIR_corrected"
)

save_heatmap_figure(
  p_gs_switcher,
  "gsnor-ko",
  "Switcher",
  "gsnor_ko_switcher_heatmap_TAIR_corrected"
)

cat("\nDone. TAIR-annotated corrected heatmaps saved in:\n")
cat(out_dir, "\n")
