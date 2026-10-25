rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(VennDiagram)
  library(grid)
  library(gridExtra)
})

# ============================================================
# 1. Input
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "DEG_Venn_FigureC_colored_FIXED")
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

genotype_order <- c("WT", "gsnor-ko")
generation_order <- c("1", "3", "5")

generation_colors <- c(
  "G1" = "#9EE8F5",   # light blue
  "G3" = "#9BF29B",   # light green
  "G5" = "#E7A7EA"    # light purple
)

title_map <- c(
  "Drought" = "Drought",
  "CO2" = "CO\u2082",
  "Ozone" = "O\u2083",
  "Heat" = "Heat",
  "Combination" = "Combi"
)

# ============================================================
# 2. Read DESeq2 result file
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

# ============================================================
# 3. Detect exact full DESeq2 result files only
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
# 4. Collect DEG gene lists
# ============================================================

deg_gene_lists <- list()
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
  
  # Keep only full DESeq2 result files:
  # Drought.txt, CO2.txt, Ozone.txt, Heat.txt, Combination.txt
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
      treatment = treatment,
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
    
    dat <- dat %>%
      mutate(
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj)
      )
    
    deg_genes <- dat %>%
      filter(
        !is.na(padj),
        !is.na(log2FoldChange),
        padj < padj_cutoff,
        abs(log2FoldChange) >= lfc_cutoff
      ) %>%
      pull(gene_id) %>%
      unique()
    
    key <- paste(genotype, treatment, paste0("G", generation), sep = "__")
    
    deg_gene_lists[[key]] <- deg_genes
  }
}

file_map_df <- bind_rows(file_map)

cat("\nDetected files used:\n")
print(file_map_df)

write.csv(
  file_map_df,
  file.path(out_dir, "Venn_file_detection_map.csv"),
  row.names = FALSE
)

if (length(deg_gene_lists) == 0) {
  stop("No DEG gene lists were created. Please check file names and folder paths.")
}

# ============================================================
# 5. Get G1/G3/G5 DEG sets
# ============================================================

get_venn_sets <- function(genotype_name, treatment_name) {
  
  keys <- paste(
    genotype_name,
    treatment_name,
    paste0("G", generation_order),
    sep = "__"
  )
  
  sets <- lapply(keys, function(k) {
    if (k %in% names(deg_gene_lists)) {
      deg_gene_lists[[k]]
    } else {
      character(0)
    }
  })
  
  names(sets) <- c("G1", "G3", "G5")
  
  sets
}

# ============================================================
# 6. Compute Venn counts
# ============================================================

get_venn_counts <- function(sets) {
  
  G1 <- unique(sets$G1)
  G3 <- unique(sets$G3)
  G5 <- unique(sets$G5)
  
  list(
    area1 = length(G1),
    area2 = length(G3),
    area3 = length(G5),
    n12 = length(intersect(G1, G3)),
    n13 = length(intersect(G1, G5)),
    n23 = length(intersect(G3, G5)),
    n123 = length(Reduce(intersect, list(G1, G3, G5))),
    
    G1_only = length(setdiff(G1, union(G3, G5))),
    G3_only = length(setdiff(G3, union(G1, G5))),
    G5_only = length(setdiff(G5, union(G1, G3))),
    G1_G3_only = length(setdiff(intersect(G1, G3), G5)),
    G1_G5_only = length(setdiff(intersect(G1, G5), G3)),
    G3_G5_only = length(setdiff(intersect(G3, G5), G1)),
    G1_G3_G5 = length(Reduce(intersect, list(G1, G3, G5)))
  )
}

# ============================================================
# 7. Helper to draw VennDiagram grob/list properly
# ============================================================

draw_vd <- function(vd_obj) {
  if (inherits(vd_obj, "gList") || is.list(vd_obj)) {
    for (i in seq_along(vd_obj)) {
      grid.draw(vd_obj[[i]])
    }
  } else {
    grid.draw(vd_obj)
  }
}

# ============================================================
# 8. Draw standard triple Venn
# ============================================================

draw_triple_grob <- function(vc,
                             title_name,
                             scaled_use = FALSE,
                             euler_use = FALSE) {
  
  g <- grid.grabExpr({
    
    grid.newpage()
    
    vd <- draw.triple.venn(
      area1 = vc$area1,
      area2 = vc$area2,
      area3 = vc$area3,
      
      n12 = vc$n12,
      n13 = vc$n13,
      n23 = vc$n23,
      n123 = vc$n123,
      
      category = c("G1", "G3", "G5"),
      
      fill = c(
        generation_colors["G1"],
        generation_colors["G3"],
        generation_colors["G5"]
      ),
      alpha = c(0.55, 0.55, 0.55),
      
      col = c("black", "black", "black"),
      lwd = 1.5,
      
      cex = 1.1,
      cat.cex = 1.1,
      cat.fontface = "plain",
      fontface = "plain",
      cat.col = c("black", "black", "black"),
      
      euler.d = euler_use,
      scaled = scaled_use,
      
      margin = 0.08,
      ind = FALSE
    )
    
    draw_vd(vd)
    
    grid.text(
      title_name,
      x = 0.5,
      y = 0.96,
      gp = gpar(fontsize = 16, fontface = "plain")
    )
  })
  
  return(g)
}

# ============================================================
# 9. Draw pairwise Venn for cases with only 2 non-empty generations
# ============================================================

draw_pairwise_grob <- function(non_empty_sets,
                               title_name,
                               scaled_use = TRUE) {
  
  set_names <- names(non_empty_sets)
  A <- unique(non_empty_sets[[1]])
  B <- unique(non_empty_sets[[2]])
  
  g <- grid.grabExpr({
    
    grid.newpage()
    
    vd <- draw.pairwise.venn(
      area1 = length(A),
      area2 = length(B),
      cross.area = length(intersect(A, B)),
      
      category = set_names,
      
      fill = generation_colors[set_names],
      alpha = rep(0.55, 2),
      
      col = c("black", "black"),
      lwd = 1.5,
      
      cex = 1.1,
      cat.cex = 1.1,
      cat.fontface = "plain",
      fontface = "plain",
      cat.col = c("black", "black"),
      
      scaled = scaled_use,
      euler.d = TRUE,
      
      margin = 0.08,
      ind = FALSE
    )
    
    draw_vd(vd)
    
    grid.text(
      title_name,
      x = 0.5,
      y = 0.96,
      gp = gpar(fontsize = 16, fontface = "plain")
    )
  })
  
  return(g)
}

# ============================================================
# 10. Draw one-circle plot for cases with only 1 non-empty generation
# ============================================================

draw_single_grob <- function(non_empty_sets,
                             title_name) {
  
  set_name <- names(non_empty_sets)[1]
  n_genes <- length(non_empty_sets[[1]])
  
  g <- grid.grabExpr({
    
    grid.newpage()
    
    grid.circle(
      x = 0.5,
      y = 0.5,
      r = 0.30,
      gp = gpar(
        fill = generation_colors[set_name],
        alpha = 0.55,
        col = "black",
        lwd = 1.5
      )
    )
    
    grid.text(
      n_genes,
      x = 0.5,
      y = 0.5,
      gp = gpar(fontsize = 14)
    )
    
    grid.text(
      set_name,
      x = 0.5,
      y = 0.17,
      gp = gpar(fontsize = 12)
    )
    
    grid.text(
      title_name,
      x = 0.5,
      y = 0.96,
      gp = gpar(fontsize = 16)
    )
  })
  
  return(g)
}

# ============================================================
# 11. Main Venn function
# ============================================================

make_venn_grob <- function(genotype_name,
                           treatment_name,
                           title_name = treatment_name) {
  
  sets <- get_venn_sets(genotype_name, treatment_name)
  vc <- get_venn_counts(sets)
  
  set_sizes <- sapply(sets, length)
  non_empty_sets <- sets[set_sizes > 0]
  
  # No DEG case
  if (length(non_empty_sets) == 0) {
    
    g <- grid.grabExpr({
      grid.newpage()
      grid.text(
        paste0(title_name, "\nNo DEGs"),
        gp = gpar(fontsize = 14)
      )
    })
    
    return(g)
  }
  
  # One generation only
  if (length(non_empty_sets) == 1) {
    return(draw_single_grob(non_empty_sets, title_name))
  }
  
  # Two generations only
  # This helps gsnor-ko Drought look like your screenshot if G5 is empty.
  if (length(non_empty_sets) == 2) {
    
    g <- tryCatch(
      draw_pairwise_grob(
        non_empty_sets = non_empty_sets,
        title_name = title_name,
        scaled_use = TRUE
      ),
      error = function(e) {
        draw_pairwise_grob(
          non_empty_sets = non_empty_sets,
          title_name = title_name,
          scaled_use = FALSE
        )
      }
    )
    
    return(g)
  }
  
  # Special screenshot-style cases:
  # WT CO2 and gsnor-ko Drought are attempted as scaled/euler plots.
  special_scaled <- (
    genotype_name == "WT" && treatment_name == "CO2"
  ) || (
    genotype_name == "gsnor-ko" && treatment_name == "Drought"
  )
  
  if (special_scaled) {
    
    g_scaled <- tryCatch(
      draw_triple_grob(
        vc = vc,
        title_name = title_name,
        scaled_use = TRUE,
        euler_use = TRUE
      ),
      error = function(e) {
        message("Scaled/euler failed for ", genotype_name, " ", treatment_name,
                ". Using standard Venn instead.")
        NULL
      }
    )
    
    if (!is.null(g_scaled)) {
      return(g_scaled)
    }
  }
  
  # Default standard 3-circle Venn
  g <- draw_triple_grob(
    vc = vc,
    title_name = title_name,
    scaled_use = FALSE,
    euler_use = FALSE
  )
  
  return(g)
}

# ============================================================
# 12. Generate individual Venn diagrams
# ============================================================

venn_grobs <- list()

for (genotype_name in genotype_order) {
  
  for (treatment_name in treatment_order) {
    
    title_name <- title_map[[treatment_name]]
    
    g <- make_venn_grob(
      genotype_name = genotype_name,
      treatment_name = treatment_name,
      title_name = title_name
    )
    
    plot_name <- paste(genotype_name, treatment_name, sep = "_")
    plot_name <- str_replace_all(plot_name, "[- ]", "_")
    
    venn_grobs[[plot_name]] <- g
    
    pdf(
      file.path(out_dir, paste0("Venn_", plot_name, ".pdf")),
      width = 4,
      height = 4
    )
    grid.newpage()
    grid.draw(g)
    dev.off()
    
    png(
      file.path(out_dir, paste0("Venn_", plot_name, ".png")),
      width = 1200,
      height = 1200,
      res = 300
    )
    grid.newpage()
    grid.draw(g)
    dev.off()
  }
}

# ============================================================
# 13. Combined Figure C-style panel
# ============================================================

get_grob_name <- function(genotype_name, treatment_name) {
  paste(genotype_name, treatment_name, sep = "_") %>%
    str_replace_all("[- ]", "_")
}

WT_grobs <- lapply(treatment_order, function(trt) {
  venn_grobs[[get_grob_name("WT", trt)]]
})

gsnor_grobs <- lapply(treatment_order, function(trt) {
  venn_grobs[[get_grob_name("gsnor-ko", trt)]]
})

blank <- nullGrob()

WT_panel <- arrangeGrob(
  grobs = list(
    WT_grobs[[1]], WT_grobs[[2]], WT_grobs[[3]],
    blank,        WT_grobs[[4]], WT_grobs[[5]]
  ),
  layout_matrix = rbind(
    c(1, 2, 3),
    c(4, 5, 6)
  )
)

gsnor_panel <- arrangeGrob(
  grobs = list(
    gsnor_grobs[[1]], gsnor_grobs[[2]], gsnor_grobs[[3]],
    blank,           gsnor_grobs[[4]], gsnor_grobs[[5]]
  ),
  layout_matrix = rbind(
    c(1, 2, 3),
    c(4, 5, 6)
  )
)

WT_title <- textGrob(
  "WT",
  gp = gpar(fontsize = 18, fontface = "bold")
)

gsnor_title <- textGrob(
  "gsnor-ko",
  gp = gpar(fontsize = 18, fontface = "bold.italic")
)

combined_panel <- arrangeGrob(
  WT_title,
  WT_panel,
  gsnor_title,
  gsnor_panel,
  ncol = 1,
  heights = c(0.08, 1, 0.08, 1)
)

pdf(
  file.path(out_dir, "FigureC_DEG_Venn_colored_combined.pdf"),
  width = 10,
  height = 12
)
grid.newpage()
grid.draw(combined_panel)
dev.off()

png(
  file.path(out_dir, "FigureC_DEG_Venn_colored_combined.png"),
  width = 3000,
  height = 3600,
  res = 300
)
grid.newpage()
grid.draw(combined_panel)
dev.off()

# ============================================================
# 14. Save overlap counts
# ============================================================

overlap_summary <- list()

for (genotype_name in genotype_order) {
  
  for (treatment_name in treatment_order) {
    
    sets <- get_venn_sets(genotype_name, treatment_name)
    vc <- get_venn_counts(sets)
    
    overlap_summary[[length(overlap_summary) + 1]] <- data.frame(
      genotype = genotype_name,
      treatment = treatment_name,
      G1_total = vc$area1,
      G3_total = vc$area2,
      G5_total = vc$area3,
      G1_only = vc$G1_only,
      G3_only = vc$G3_only,
      G5_only = vc$G5_only,
      G1_G3_only = vc$G1_G3_only,
      G1_G5_only = vc$G1_G5_only,
      G3_G5_only = vc$G3_G5_only,
      G1_G3_G5 = vc$G1_G3_G5
    )
  }
}

overlap_summary <- bind_rows(overlap_summary)

write.csv(
  overlap_summary,
  file.path(out_dir, "FigureC_Venn_overlap_counts.csv"),
  row.names = FALSE
)

cat("\nOverlap summary:\n")
print(overlap_summary)

cat("\nFixed coloured Venn diagrams saved in:\n")
cat(out_dir, "\n")

