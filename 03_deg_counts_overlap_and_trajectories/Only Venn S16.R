rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(VennDiagram)
  library(grid)
})

# ============================================================
# 1. Input paths and settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "DEG_UpSet_Venn_individual_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 0

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

treatments <- c(
  "Drought", "CO2", "Ozone", "Heat", "Combination"
)

genotypes <- c("WT", "gsnor-ko")
generations <- c("G1", "G3", "G5")

# More visible intersections for Heat and Combination
top_n_upset <- c(
  "Drought" = 40,
  "CO2" = 40,
  "Ozone" = 40,
  "Heat" = 70,
  "Combination" = 70
)

# Original Venn colours
venn_cols <- c(
  "Drought" = "red",
  "CO2" = "blue",
  "Ozone" = "green",
  "Heat" = "yellow",
  "Combination" = "purple"
)

# ============================================================
# 2. Helper functions
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

get_genotype_generation <- function(folder) {
  
  genotype <- ifelse(str_detect(folder, "^WT"), "WT", "gsnor-ko")
  
  generation <- str_extract(folder, "gen[0-9]+") %>%
    str_remove("gen") %>%
    paste0("G", .)
  
  data.frame(
    folder = folder,
    genotype = genotype,
    generation = generation,
    stringsAsFactors = FALSE
  )
}

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

get_deg_genes <- function(dat, padj_cutoff = 0.05, lfc_cutoff = 0) {
  
  if (!all(c("log2FoldChange", "padj") %in% colnames(dat))) {
    stop("Required columns log2FoldChange and padj are missing.")
  }
  
  dat %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < padj_cutoff,
      abs(log2FoldChange) >= lfc_cutoff
    ) %>%
    pull(gene_id) %>%
    unique()
}

draw_empty_page <- function(label) {
  grid::grid.newpage()
  grid::grid.text(
    label,
    x = 0.5,
    y = 0.5,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
}

# ============================================================
# 3. Read all DESeq2 result files and make DEG lists
# ============================================================

deg_gene_lists <- list()
deg_records <- list()
file_map <- list()

for (folder in folders) {
  
  folder_path <- file.path(root_dir, folder)
  
  if (!dir.exists(folder_path)) {
    warning("Folder not found: ", folder_path)
    next
  }
  
  fg <- get_genotype_generation(folder)
  
  all_files <- list.files(
    folder_path,
    pattern = "\\.(txt|csv|tsv)$",
    full.names = TRUE
  )
  
  # Use only full DESeq2 result files
  files <- all_files[
    !str_detect(
      tolower(basename(all_files)),
      "upregulated|downregulated|pval|pvalue|pvales"
    )
  ]
  
  for (file in files) {
    
    treatment <- get_treatment_exact(file)
    
    file_map[[length(file_map) + 1]] <- data.frame(
      folder = folder,
      genotype = fg$genotype,
      generation = fg$generation,
      treatment = treatment,
      file = basename(file),
      stringsAsFactors = FALSE
    )
    
    if (is.na(treatment)) next
    
    dat <- read_deseq_file(file)
    
    if (!all(c("log2FoldChange", "padj") %in% colnames(dat))) {
      warning("Skipping file because required columns are missing: ", file)
      next
    }
    
    deg_genes <- get_deg_genes(
      dat,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    )
    
    key <- paste(fg$genotype, fg$generation, treatment, sep = "__")
    deg_gene_lists[[key]] <- deg_genes
    
    if (length(deg_genes) > 0) {
      deg_records[[length(deg_records) + 1]] <- data.frame(
        genotype = rep(fg$genotype, length(deg_genes)),
        generation = rep(fg$generation, length(deg_genes)),
        treatment = rep(treatment, length(deg_genes)),
        gene_id = deg_genes,
        stringsAsFactors = FALSE
      )
    } else {
      message("No DEGs found for: ", key)
    }
  }
}

file_map_df <- bind_rows(file_map)

write.csv(
  file_map_df,
  file.path(out_dir, "DESeq2_file_detection_map.csv"),
  row.names = FALSE
)

if (length(deg_records) > 0) {
  deg_df <- bind_rows(deg_records)
} else {
  deg_df <- data.frame(
    genotype = character(),
    generation = character(),
    treatment = character(),
    gene_id = character(),
    stringsAsFactors = FALSE
  )
}

write.csv(
  deg_df,
  file.path(out_dir, "all_DEG_gene_lists_long.csv"),
  row.names = FALSE
)

deg_counts <- deg_df %>%
  count(genotype, generation, treatment, name = "n_DEGs")

all_combinations <- tidyr::expand_grid(
  genotype = genotypes,
  generation = generations,
  treatment = treatments
)

deg_counts <- all_combinations %>%
  left_join(
    deg_counts,
    by = c("genotype", "generation", "treatment")
  ) %>%
  mutate(
    n_DEGs = ifelse(is.na(n_DEGs), 0, n_DEGs)
  ) %>%
  arrange(genotype, generation, treatment)

write.csv(
  deg_counts,
  file.path(out_dir, "DEG_counts_summary.csv"),
  row.names = FALSE
)

cat("\nDEG counts summary:\n")
print(deg_counts)

# ============================================================
# 4. Individual manual UpSet plots
# Single-dot intersections first, then combined intersections
# No combined figures
# ============================================================

make_upset_sets_for_treatment <- function(trt) {
  
  set_order <- c(
    "WT_gen1",
    "gsnor_gen1",
    "WT_gen3",
    "gsnor_gen3",
    "WT_gen5",
    "gsnor_gen5"
  )
  
  set_keys <- c(
    "WT_gen1"     = paste("WT", "G1", trt, sep = "__"),
    "gsnor_gen1"  = paste("gsnor-ko", "G1", trt, sep = "__"),
    "WT_gen3"     = paste("WT", "G3", trt, sep = "__"),
    "gsnor_gen3"  = paste("gsnor-ko", "G3", trt, sep = "__"),
    "WT_gen5"     = paste("WT", "G5", trt, sep = "__"),
    "gsnor_gen5"  = paste("gsnor-ko", "G5", trt, sep = "__")
  )
  
  upset_sets <- list()
  
  for (nm in set_order) {
    key <- set_keys[[nm]]
    
    if (key %in% names(deg_gene_lists)) {
      upset_sets[[nm]] <- deg_gene_lists[[key]]
    } else {
      warning("Missing DEG list for: ", key)
      upset_sets[[nm]] <- character(0)
    }
  }
  
  upset_sets
}

make_membership_df <- function(set_list) {
  
  all_genes <- sort(unique(unlist(set_list)))
  
  if (length(all_genes) == 0) {
    return(NULL)
  }
  
  mem <- data.frame(
    gene_id = all_genes,
    stringsAsFactors = FALSE
  )
  
  for (set_name in names(set_list)) {
    mem[[set_name]] <- all_genes %in% set_list[[set_name]]
  }
  
  mem
}

make_manual_upset_plot <- function(set_list, trt, top_n = 40) {
  
  set_order <- names(set_list)
  mem <- make_membership_df(set_list)
  
  if (is.null(mem)) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = paste0("UpSet plot - ", trt, "\nNo DEGs found"),
          size = 6,
          fontface = "bold"
        ) +
        theme_void()
    )
  }
  
  # ------------------------------------------------------------
  # Set sizes
  # ------------------------------------------------------------
  
  set_size_df <- data.frame(
    set = set_order,
    set_size = lengths(set_list),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      set = factor(set, levels = rev(set_order))
    )
  
  # ------------------------------------------------------------
  # Exact intersections
  # ------------------------------------------------------------
  
  mem2 <- mem %>%
    rowwise() %>%
    mutate(
      intersection = paste(
        set_order[c_across(all_of(set_order))],
        collapse = "|"
      )
    ) %>%
    ungroup()
  
  inter_df <- mem2 %>%
    count(intersection, name = "intersection_size") %>%
    filter(intersection != "") %>%
    mutate(
      degree = stringr::str_count(intersection, "\\|") + 1
    ) %>%
    arrange(
      degree,
      desc(intersection_size)
    ) %>%
    slice_head(n = top_n) %>%
    mutate(
      intersection_id = paste0("I", row_number()),
      intersection_id = factor(intersection_id, levels = intersection_id)
    )
  
  if (nrow(inter_df) == 0) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = paste0("UpSet plot - ", trt, "\nNo non-empty intersections"),
          size = 6,
          fontface = "bold"
        ) +
        theme_void()
    )
  }
  
  # ------------------------------------------------------------
  # Dot matrix data
  # ------------------------------------------------------------
  
  matrix_df <- inter_df %>%
    select(intersection_id, intersection) %>%
    tidyr::separate_rows(intersection, sep = "\\|") %>%
    rename(set = intersection) %>%
    mutate(
      set = factor(set, levels = rev(set_order))
    )
  
  all_matrix_df <- tidyr::expand_grid(
    intersection_id = inter_df$intersection_id,
    set = factor(rev(set_order), levels = rev(set_order))
  ) %>%
    left_join(
      matrix_df %>% mutate(active = TRUE),
      by = c("intersection_id", "set")
    ) %>%
    mutate(
      active = ifelse(is.na(active), FALSE, active)
    )
  
  line_df <- matrix_df %>%
    group_by(intersection_id) %>%
    summarise(
      ymin = min(as.numeric(set)),
      ymax = max(as.numeric(set)),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------
  # Main bar plot
  # ------------------------------------------------------------
  
  p_bar <- ggplot(
    inter_df,
    aes(x = intersection_id, y = intersection_size)
  ) +
    geom_col(fill = "grey25", width = 0.75) +
    geom_text(
      aes(label = intersection_size),
      vjust = -0.25,
      size = 3
    ) +
    labs(
      title = paste0("UpSet plot- ", trt),
      y = "Intersection Size",
      x = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(color = "black")
    )
  
  # ------------------------------------------------------------
  # Dot matrix plot
  # ------------------------------------------------------------
  
  p_matrix <- ggplot(
    all_matrix_df,
    aes(x = intersection_id, y = set)
  ) +
    geom_point(
      aes(color = active),
      size = 2.8
    ) +
    geom_segment(
      data = line_df,
      aes(
        x = intersection_id,
        xend = intersection_id,
        y = ymin,
        yend = ymax
      ),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    scale_color_manual(
      values = c(
        "TRUE" = "grey25",
        "FALSE" = "grey85"
      ),
      guide = "none"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.ticks.y = element_blank(),
      axis.line = element_blank()
    )
  
  # ------------------------------------------------------------
  # Set-size bar plot
  # ------------------------------------------------------------
  
  p_setsize <- ggplot(
    set_size_df,
    aes(x = set_size, y = set)
  ) +
    geom_col(fill = "grey25", width = 0.7) +
    geom_text(
      aes(label = set_size),
      hjust = -0.15,
      size = 3
    ) +
    labs(x = "Set Size", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(color = "black")
    ) +
    expand_limits(x = max(set_size_df$set_size, na.rm = TRUE) * 1.15)
  
  # ------------------------------------------------------------
  # Combine components of one individual UpSet plot
  # ------------------------------------------------------------
  
  p_right <- p_bar / p_matrix +
    plot_layout(heights = c(2.2, 1))
  
  p_final <- p_setsize | p_right
  
  p_final <- p_final +
    plot_layout(widths = c(1.1, 4))
  
  return(p_final)
}

save_manual_upset_for_treatment <- function(trt, out_dir) {
  
  upset_sets <- make_upset_sets_for_treatment(trt)
  
  set_sizes <- data.frame(
    treatment = trt,
    set = names(upset_sets),
    n_genes = lengths(upset_sets),
    stringsAsFactors = FALSE
  )
  
  print(set_sizes)
  
  write.csv(
    set_sizes,
    file.path(out_dir, paste0("UpSet_", trt, "_set_sizes.csv")),
    row.names = FALSE
  )
  
  p <- make_manual_upset_plot(
    set_list = upset_sets,
    trt = trt,
    top_n = top_n_upset[[trt]]
  )
  
  ggsave(
    file.path(out_dir, paste0("UpSet_", trt, "_manual.pdf")),
    p,
    width = 12,
    height = 6,
    useDingbats = FALSE
  )
  
  ggsave(
    file.path(out_dir, paste0("UpSet_", trt, "_manual.png")),
    p,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

for (trt in treatments) {
  
  cat("\nMaking individual UpSet plot for:", trt, "\n")
  
  save_manual_upset_for_treatment(
    trt = trt,
    out_dir = out_dir
  )
}

# ============================================================
# 5. Individual Venn diagrams
# Original colours
# Special gsnor-ko G5 layout like attached
# No Cairo
# ============================================================

venn_order_default <- c(
  "Drought",
  "CO2",
  "Ozone",
  "Heat",
  "Combination"
)

# For gsnor-ko G5 only:
# CO2 left, Heat top, Combination top/right, Ozone right, Drought bottom
venn_order_gsnor_G5 <- c(
  "CO2",
  "Heat",
  "Combination",
  "Ozone",
  "Drought"
)

make_venn_sets_for_group <- function(gt, gen) {
  
  order_use <- if (gt == "gsnor-ko" && gen == "G5") {
    venn_order_gsnor_G5
  } else {
    venn_order_default
  }
  
  venn_sets <- list()
  
  for (trt in order_use) {
    
    key <- paste(gt, gen, trt, sep = "__")
    
    if (key %in% names(deg_gene_lists)) {
      venn_sets[[trt]] <- deg_gene_lists[[key]]
    } else {
      warning("Missing DEG list for: ", key)
      venn_sets[[trt]] <- character(0)
    }
  }
  
  venn_sets
}

draw_venn_plot <- function(venn_sets, gt, gen) {
  
  is_special <- gt == "gsnor-ko" && gen == "G5"
  
  if (is_special) {
    
    venn_grob <- VennDiagram::venn.diagram(
      x = venn_sets,
      filename = NULL,
      category.names = names(venn_sets),
      fill = venn_cols[names(venn_sets)],
      alpha = 0.45,
      col = "black",
      lwd = 1.2,
      cex = 0.80,
      cat.cex = 0.95,
      cat.col = "black",
      cat.fontface = "plain",
      fontface = "plain",
      margin = 0.01,
      scaled = FALSE,
      euler.d = FALSE,
      
      # Label positions for attached-like gsnor-ko G5 layout
      # CO2 left, Heat top, Combination top/right, Ozone right, Drought bottom
      cat.pos = c(180, 90, 55, 0, -90),
      cat.dist = c(0.05, 0.04, 0.05, 0.05, 0.04)
    )
    
  } else {
    
    venn_grob <- VennDiagram::venn.diagram(
      x = venn_sets,
      filename = NULL,
      category.names = names(venn_sets),
      fill = venn_cols[names(venn_sets)],
      alpha = 0.45,
      col = "black",
      lwd = 1.2,
      cex = 0.85,
      cat.cex = 0.90,
      cat.col = "black",
      cat.fontface = "plain",
      fontface = "plain",
      margin = 0.04,
      scaled = FALSE,
      euler.d = FALSE
    )
  }
  
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  
  grid::grid.text(
    paste(gt, gen),
    x = 0.5,
    y = 0.98,
    gp = grid::gpar(fontsize = 10, fontface = "bold")
  )
}

save_venn_plot <- function(gt, gen, out_dir) {
  
  venn_sets <- make_venn_sets_for_group(gt, gen)
  
  set_sizes <- data.frame(
    genotype = gt,
    generation = gen,
    treatment = names(venn_sets),
    n_genes = lengths(venn_sets),
    stringsAsFactors = FALSE
  )
  
  print(set_sizes)
  
  safe_gt <- gsub("[^A-Za-z0-9]+", "_", gt)
  
  write.csv(
    set_sizes,
    file.path(out_dir, paste0("Venn_", safe_gt, "_", gen, "_set_sizes.csv")),
    row.names = FALSE
  )
  
  pdf_file <- file.path(out_dir, paste0("Venn_", safe_gt, "_", gen, ".pdf"))
  png_file <- file.path(out_dir, paste0("Venn_", safe_gt, "_", gen, ".png"))
  
  if (all(lengths(venn_sets) == 0)) {
    
    pdf(pdf_file, width = 4.5, height = 4.5, useDingbats = FALSE)
    draw_empty_page(paste(gt, gen, "\nNo DEGs found"))
    dev.off()
    
    png(png_file, width = 4.5, height = 4.5, units = "in", res = 300)
    draw_empty_page(paste(gt, gen, "\nNo DEGs found"))
    dev.off()
    
    return(invisible(NULL))
  }
  
  pdf(pdf_file, width = 4.5, height = 4.5, useDingbats = FALSE)
  tryCatch(
    draw_venn_plot(venn_sets, gt, gen),
    error = function(e) {
      draw_empty_page(paste0(gt, " ", gen, "\nVenn failed"))
      warning(e$message)
    }
  )
  dev.off()
  
  png(png_file, width = 4.5, height = 4.5, units = "in", res = 300)
  tryCatch(
    draw_venn_plot(venn_sets, gt, gen),
    error = function(e) {
      draw_empty_page(paste0(gt, " ", gen, "\nVenn failed"))
      warning(e$message)
    }
  )
  dev.off()
}

for (gen in generations) {
  for (gt in genotypes) {
    
    cat("\nSaving individual Venn plot for:", gt, gen, "\n")
    
    save_venn_plot(
      gt = gt,
      gen = gen,
      out_dir = out_dir
    )
  }
}

cat("\nDone. Individual UpSet and Venn plots saved in:\n")
cat(out_dir, "\n")



###venn diagrma (wqorking)
# ============================================================
# 5. Individual Venn diagrams
# No Cairo
# Original colours
# gsnor-ko G5 becomes 4-set Venn because Drought = 0
# ============================================================

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
})

venn_cols <- c(
  "Drought"     = "red",
  "CO2"         = "blue",
  "Ozone"       = "green",
  "Heat"        = "yellow",
  "Combination" = "purple"
)

venn_order_default <- c(
  "Drought",
  "CO2",
  "Ozone",
  "Heat",
  "Combination"
)

# For gsnor-ko G5, Drought has 0 DEGs, so make it a 4-set Venn
# This matches your attached example.
venn_order_gsnor_G5 <- c(
  "CO2",
  "Heat",
  "Combination",
  "Ozone"
)

make_venn_sets_for_group <- function(gt, gen) {
  
  order_use <- if (gt == "gsnor-ko" && gen == "G5") {
    venn_order_gsnor_G5
  } else {
    venn_order_default
  }
  
  venn_sets <- list()
  
  for (trt in order_use) {
    
    key <- paste(gt, gen, trt, sep = "__")
    
    if (key %in% names(deg_gene_lists)) {
      venn_sets[[trt]] <- deg_gene_lists[[key]]
    } else {
      warning("Missing DEG list for: ", key)
      venn_sets[[trt]] <- character(0)
    }
  }
  
  # Remove empty sets only for Venn plotting
  # This prevents distorted zero-set circles.
  venn_sets <- venn_sets[lengths(venn_sets) > 0]
  
  venn_sets
}

draw_empty_page <- function(label) {
  grid::grid.newpage()
  grid::grid.text(
    label,
    x = 0.5,
    y = 0.5,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
}

draw_venn_plot <- function(venn_sets, gt, gen) {
  
  n_sets <- length(venn_sets)
  
  if (n_sets < 2) {
    draw_empty_page(paste(gt, gen, "\nLess than 2 non-empty sets"))
    return(invisible(NULL))
  }
  
  is_special <- gt == "gsnor-ko" && gen == "G5"
  
  if (is_special) {
    
    venn_grob <- VennDiagram::venn.diagram(
      x = venn_sets,
      filename = NULL,
      category.names = names(venn_sets),
      fill = venn_cols[names(venn_sets)],
      alpha = 0.45,
      col = "black",
      lwd = 1.2,
      cex = 0.85,
      cat.cex = 0.90,
      cat.col = "black",
      cat.fontface = "plain",
      fontface = "plain",
      margin = 0.04,
      scaled = FALSE,
      euler.d = FALSE,
      
      # 4-set label positions:
      # CO2 left, Heat top, Combination top/right, Ozone right
      cat.pos = c(180, 90, 35, 0),
      cat.dist = c(0.05, 0.05, 0.05, 0.05)
    )
    
  } else {
    
    venn_grob <- VennDiagram::venn.diagram(
      x = venn_sets,
      filename = NULL,
      category.names = names(venn_sets),
      fill = venn_cols[names(venn_sets)],
      alpha = 0.45,
      col = "black",
      lwd = 1.2,
      cex = 0.78,
      cat.cex = 0.80,
      cat.col = "black",
      cat.fontface = "plain",
      fontface = "plain",
      margin = 0.05,
      scaled = FALSE,
      euler.d = FALSE
    )
  }
  
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  
  invisible(venn_grob)
}

save_venn_plot <- function(gt, gen, out_dir) {
  
  venn_sets <- make_venn_sets_for_group(gt, gen)
  
  set_sizes <- data.frame(
    genotype = gt,
    generation = gen,
    treatment = names(venn_sets),
    n_genes = lengths(venn_sets),
    stringsAsFactors = FALSE
  )
  
  print(set_sizes)
  
  safe_gt <- gsub("[^A-Za-z0-9]+", "_", gt)
  
  write.csv(
    set_sizes,
    file.path(out_dir, paste0("Venn_", safe_gt, "_", gen, "_set_sizes.csv")),
    row.names = FALSE
  )
  
  pdf_file <- file.path(out_dir, paste0("Venn_", safe_gt, "_", gen, ".pdf"))
  png_file <- file.path(out_dir, paste0("Venn_", safe_gt, "_", gen, ".png"))
  
  pdf(pdf_file, width = 4, height = 4, useDingbats = FALSE)
  tryCatch(
    draw_venn_plot(venn_sets, gt, gen),
    error = function(e) {
      draw_empty_page(paste0(gt, " ", gen, "\nVenn failed"))
      warning(e$message)
    }
  )
  dev.off()
  
  png(png_file, width = 4, height = 4, units = "in", res = 300)
  tryCatch(
    draw_venn_plot(venn_sets, gt, gen),
    error = function(e) {
      draw_empty_page(paste0(gt, " ", gen, "\nVenn failed"))
      warning(e$message)
    }
  )
  dev.off()
}

# ============================================================
# Save individual Venn plots only
# ============================================================

for (gen in generations) {
  for (gt in genotypes) {
    
    cat("\nSaving individual Venn plot for:", gt, gen, "\n")
    
    save_venn_plot(
      gt = gt,
      gen = gen,
      out_dir = out_dir
    )
  }
}



