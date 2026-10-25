## ============================================================
## Mercator / MapMan enrichment for selected WGCNA modules
## Reference-style plot
##
## Modules:
## magenta, lightcyan1, darkgreen, ivory
##
## Special rule:
## If lightcyan1 has no significant terms, plot top 5 terms
## ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

## ============================================================
## 1. Paths
## ============================================================

wd <- "~/Documents/DNA_seq/RNA_seq/WGCNA"

wgcna_outdir <- file.path(wd, "WGCNA_final")

outdir <- file.path(wgcna_outdir, "Mercator_module_enrichment_selected_modules")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

gene_module_file <- file.path(wgcna_outdir, "genes_modules_info.csv")

mercator_file <- "~/Documents/DNA_seq/RNA_seq/Mercator/Mercator_Annot_Arabidopsis_singlecopy_geneid.csv"

## ============================================================
## 2. User settings
## ============================================================

target_modules <- c(
  "magenta",
  "lightcyan1",
  "darkgreen",
  "ivory"
)

## Significance cutoffs
PV_CUT <- 0.05
Q_CUT  <- 0.20

## Mercator term size filters
MIN_GENES_IN_TERM <- 3
MAX_GENES_IN_TERM <- 500

## Special fallback
FALLBACK_MODULE <- "lightcyan1"
FALLBACK_TOP_N <- 5

## ============================================================
## 3. Read WGCNA gene-module table
## ============================================================

gene_mod <- data.table::fread(gene_module_file, data.table = FALSE)

cat("\nColumns in gene-module file:\n")
print(colnames(gene_mod))

if (!"ID" %in% colnames(gene_mod)) {
  stop("genes_modules_info.csv must contain a gene ID column named 'ID'.")
}

if (!"module_color" %in% colnames(gene_mod)) {
  stop("genes_modules_info.csv must contain a module column named 'module_color'.")
}

gene_mod <- gene_mod %>%
  dplyr::mutate(
    gene_id = toupper(as.character(ID)),
    gene_id = stringr::str_replace(gene_id, "\\.\\d+$", ""),
    module_color = as.character(module_color)
  ) %>%
  dplyr::filter(!is.na(gene_id), gene_id != "") %>%
  dplyr::distinct(gene_id, module_color)

cat("\nAvailable WGCNA modules:\n")
print(sort(unique(gene_mod$module_color)))

missing_modules <- setdiff(target_modules, unique(gene_mod$module_color))

if (length(missing_modules) > 0) {
  stop(
    paste0(
      "These target modules were not found in genes_modules_info.csv: ",
      paste(missing_modules, collapse = ", "),
      "\n\nAvailable modules are:\n",
      paste(sort(unique(gene_mod$module_color)), collapse = ", ")
    )
  )
}

cat("\nTarget module sizes in WGCNA table:\n")
print(
  gene_mod %>%
    dplyr::filter(module_color %in% target_modules) %>%
    dplyr::group_by(module_color) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(gene_id),
      .groups = "drop"
    )
)

## ============================================================
## 4. Read Mercator annotation file
## Expected structure:
## Arabidopsis_gene_id, level1, level2, ..., level8
## ============================================================

merc_raw <- data.table::fread(
  mercator_file,
  sep = ",",
  header = TRUE,
  data.table = FALSE,
  encoding = "UTF-8"
)

colnames(merc_raw) <- trimws(colnames(merc_raw))
colnames(merc_raw) <- gsub("^\ufeff", "", colnames(merc_raw))

cat("\nColumns in Mercator file:\n")
print(colnames(merc_raw))

possible_gene_cols <- c(
  "Arabidopsis_gene_id",
  "arabidopsis_gene_id",
  "gene_id",
  "GeneID",
  "Geneid",
  "Gene",
  "gene",
  "ID",
  "Identifier",
  "IDENTIFIER",
  "AGI",
  "Locus",
  "locus"
)

gene_col <- possible_gene_cols[possible_gene_cols %in% colnames(merc_raw)][1]

if (is.na(gene_col)) {
  stop(
    paste0(
      "Could not detect gene ID column in Mercator file.\n",
      "Available columns are:\n",
      paste(colnames(merc_raw), collapse = ", ")
    )
  )
}

cat("\nDetected Mercator gene column:\n")
print(gene_col)

level_cols <- grep(
  "^level[0-9]+$",
  colnames(merc_raw),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nDetected Mercator level columns:\n")
print(level_cols)

if (length(level_cols) == 0) {
  stop("No columns named level1, level2, ... were detected in Mercator file.")
}

## ============================================================
## 5. Build TERM2GENE table
## ============================================================

term2gene <- merc_raw %>%
  dplyr::mutate(
    gene_id = toupper(as.character(.data[[gene_col]])),
    gene_id = stringr::str_replace(gene_id, "\\.\\d+$", "")
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(level_cols),
    names_to = "mercator_level",
    values_to = "mercator_term"
  ) %>%
  dplyr::mutate(
    mercator_level = tolower(as.character(mercator_level)),
    mercator_term = as.character(mercator_term),
    mercator_term = stringr::str_squish(mercator_term)
  ) %>%
  dplyr::filter(
    !is.na(gene_id),
    gene_id != "",
    !is.na(mercator_term),
    mercator_term != "",
    mercator_term != "NA",
    mercator_term != "na",
    mercator_term != "-",
    !stringr::str_detect(tolower(mercator_term), "^not assigned$"),
    !stringr::str_detect(tolower(mercator_term), "^unassigned$")
  ) %>%
  dplyr::distinct(gene_id, mercator_level, mercator_term) %>%
  dplyr::mutate(
    term_id = paste(mercator_level, mercator_term, sep = "__")
  )

term2gene$mercator_level <- factor(
  term2gene$mercator_level,
  levels = paste0("level", 1:8)
)

cat("\nTERM2GENE dimensions:\n")
print(dim(term2gene))

cat("\nExample TERM2GENE rows:\n")
print(head(term2gene, 20))

write.csv(
  term2gene,
  file.path(outdir, "TERM2GENE_Mercator_long_format.csv"),
  row.names = FALSE
)

## ============================================================
## 6. Define universe/background
## Universe = WGCNA genes with Mercator annotation
## ============================================================

universe_genes <- intersect(
  unique(gene_mod$gene_id),
  unique(term2gene$gene_id)
)

cat("\nUniverse genes present in both WGCNA and Mercator:\n")
print(length(universe_genes))

if (length(universe_genes) < 100) {
  warning("Very few genes overlap between WGCNA and Mercator annotation. Check gene IDs.")
}

gene_mod_bg <- gene_mod %>%
  dplyr::filter(gene_id %in% universe_genes)

term2gene_bg <- term2gene %>%
  dplyr::filter(gene_id %in% universe_genes)

cat("\nTarget module sizes after matching Mercator universe:\n")
print(
  gene_mod_bg %>%
    dplyr::filter(module_color %in% target_modules) %>%
    dplyr::group_by(module_color) %>%
    dplyr::summarise(
      n_genes_in_universe = dplyr::n_distinct(gene_id),
      .groups = "drop"
    )
)

## Filter terms by size
term_sizes <- term2gene_bg %>%
  dplyr::distinct(term_id, gene_id) %>%
  dplyr::group_by(term_id) %>%
  dplyr::summarise(
    term_size = dplyr::n_distinct(gene_id),
    .groups = "drop"
  )

valid_terms <- term_sizes %>%
  dplyr::filter(
    term_size >= MIN_GENES_IN_TERM,
    term_size <= MAX_GENES_IN_TERM
  ) %>%
  dplyr::pull(term_id)

term2gene_bg <- term2gene_bg %>%
  dplyr::filter(term_id %in% valid_terms)

cat("\nNumber of Mercator terms after size filtering:\n")
print(length(unique(term2gene_bg$term_id)))

## ============================================================
## 7. Hypergeometric enrichment function
## ============================================================

run_hypergeom_enrichment <- function(module_name, gene_mod_bg, term2gene_bg, universe_genes) {
  
  module_genes <- gene_mod_bg %>%
    dplyr::filter(module_color == module_name) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  module_genes <- intersect(module_genes, universe_genes)
  
  N <- length(universe_genes)
  n <- length(module_genes)
  
  cat("\nRunning enrichment for module:", module_name, "\n")
  cat("Module genes in universe:", n, "\n")
  
  if (n < 3) {
    warning(paste("Skipping module", module_name, "- fewer than 3 genes in universe."))
    return(NULL)
  }
  
  term_list <- term2gene_bg %>%
    dplyr::distinct(term_id, mercator_level, mercator_term, gene_id) %>%
    dplyr::group_by(term_id, mercator_level, mercator_term) %>%
    dplyr::summarise(
      term_genes = list(unique(gene_id)),
      term_size = length(unique(gene_id)),
      .groups = "drop"
    )
  
  res <- term_list %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      overlap_genes = list(intersect(module_genes, term_genes)),
      overlap_count = length(overlap_genes),
      p_value = phyper(
        q = overlap_count - 1,
        m = term_size,
        n = N - term_size,
        k = n,
        lower.tail = FALSE
      ),
      gene_ratio = overlap_count / n,
      bg_ratio = term_size / N,
      enrichment_ratio = gene_ratio / bg_ratio,
      overlap_gene_ids = paste(overlap_genes, collapse = ";")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(overlap_count > 0) %>%
    dplyr::mutate(
      module = module_name,
      p_adjust = p.adjust(p_value, method = "BH")
    ) %>%
    dplyr::select(
      module,
      term_id,
      mercator_level,
      mercator_term,
      overlap_count,
      term_size,
      p_value,
      p_adjust,
      gene_ratio,
      bg_ratio,
      enrichment_ratio,
      overlap_gene_ids
    ) %>%
    dplyr::arrange(p_adjust, p_value, dplyr::desc(overlap_count))
  
  return(res)
}

## ============================================================
## 8. Run enrichment for selected modules
## ============================================================

enrichment_selected <- purrr::map_dfr(
  target_modules,
  run_hypergeom_enrichment,
  gene_mod_bg = gene_mod_bg,
  term2gene_bg = term2gene_bg,
  universe_genes = universe_genes
)

if (nrow(enrichment_selected) == 0) {
  stop("No enrichment results produced. Check module names and annotation.")
}

write.csv(
  enrichment_selected,
  file.path(outdir, "Selected_modules_Mercator_enrichment_all_results.csv"),
  row.names = FALSE
)

cat("\nTop enrichment results:\n")
print(head(enrichment_selected, 30))

## ============================================================
## 9. Choose terms for plotting
## Max 7 significant terms per module
## If lightcyan1 has no significant terms, plot top 5 terms
## ============================================================

MAX_SIG_TERMS_PER_MODULE <- 7
FALLBACK_MODULE <- "lightcyan1"
FALLBACK_TOP_N <- 5

## IMPORTANT: new module order
module_plot_order <- c("ivory", "magenta", "lightcyan1", "darkgreen")

## IMPORTANT: desired Mercator level order inside each module
level_plot_order <- c("level1", "level2", "level3", "level4", "level5")

sig_selected <- enrichment_selected %>%
  dplyr::filter(
    p_value <= PV_CUT,
    p_adjust <= Q_CUT
  )

plot_df_list <- list()

for (m in module_plot_order) {
  
  sig_m <- sig_selected %>%
    dplyr::filter(module == m) %>%
    dplyr::arrange(p_adjust, p_value, dplyr::desc(overlap_count)) %>%
    dplyr::slice_head(n = MAX_SIG_TERMS_PER_MODULE)
  
  if (nrow(sig_m) > 0) {
    
    plot_df_list[[m]] <- sig_m %>%
      dplyr::mutate(plot_reason = "significant_top7")
    
  } else if (m == FALLBACK_MODULE) {
    
    message(
      paste0(
        "No significant terms found for ", FALLBACK_MODULE,
        ". Plotting top ", FALLBACK_TOP_N, " ranked terms instead."
      )
    )
    
    plot_df_list[[m]] <- enrichment_selected %>%
      dplyr::filter(module == m) %>%
      dplyr::arrange(p_adjust, p_value, dplyr::desc(overlap_count)) %>%
      dplyr::slice_head(n = FALLBACK_TOP_N) %>%
      dplyr::mutate(plot_reason = "top5_fallback_not_significant")
    
  } else {
    
    warning(
      paste0(
        "No significant terms found for module ", m,
        ". This module will not appear in the plot."
      )
    )
  }
}

plot_df <- dplyr::bind_rows(plot_df_list)

if (nrow(plot_df) == 0) {
  stop("No terms selected for plotting.")
}

## ============================================================
## Reorder terms exactly how you want
## ============================================================

plot_df <- plot_df %>%
  dplyr::mutate(
    module = factor(module, levels = module_plot_order),
    mercator_level = as.character(mercator_level),
    mercator_level = factor(mercator_level, levels = level_plot_order),
    level_num = as.numeric(mercator_level),
    p_adjust_safe = pmax(p_adjust, .Machine$double.xmin),
    neg_log10_padj = -log10(p_adjust_safe),
    mercator_term_clean = stringr::str_squish(mercator_term),
    mercator_term_wrapped = stringr::str_wrap(mercator_term_clean, width = 34)
  ) %>%
  dplyr::arrange(
    module,                 # first arrange by module
    level_num,              # then level1 -> level5
    p_adjust,               # then more significant first
    p_value,
    dplyr::desc(overlap_count),
    mercator_term_clean
  ) %>%
  dplyr::mutate(
    row_id = dplyr::row_number()
  )

## reverse y so first rows appear at the top
plot_df <- plot_df %>%
  dplyr::mutate(
    y_id = factor(row_id, levels = rev(row_id))
  )

cat("\nNumber of plotted terms per module:\n")
print(
  plot_df %>%
    dplyr::group_by(module, mercator_level, plot_reason) %>%
    dplyr::summarise(
      n_terms = dplyr::n(),
      min_p_adjust = min(p_adjust, na.rm = TRUE),
      max_p_adjust = max(p_adjust, na.rm = TRUE),
      .groups = "drop"
    )
)

write.csv(
  plot_df,
  file.path(outdir, "Selected_modules_terms_used_for_reference_style_plot.csv"),
  row.names = FALSE
)
## ============================================================
## 10. Reference-style bubble plot
## ============================================================

graphics.off()

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

## module colors in the SAME order
module_colors <- c(
  ivory      = "ivory3",
  magenta    = "magenta",
  lightcyan1 = "lightcyan1",
  darkgreen  = "darkgreen"
)

## level colors
level_colors <- c(
  level1 = "#00A65A",
  level2 = "#E69F00",
  level3 = "#7570B3",
  level4 = "#E7298A",
  level5 = "#66A61E"
)

x_max <- ceiling(max(plot_df$neg_log10_padj, na.rm = TRUE))
x_break_max <- ceiling(x_max / 4) * 4

if (x_break_max < 4) {
  x_break_max <- 4
}

x_breaks <- seq(0, x_break_max, by = 4)

bar_module_x <- x_break_max + 0.50
bar_level_x  <- x_break_max + 0.95
x_limit_max  <- x_break_max + 1.30

y_labels <- setNames(plot_df$mercator_term_wrapped, plot_df$y_id)

p_reference <- ggplot(plot_df, aes(y = y_id)) +
  
  geom_point(
    aes(
      x = neg_log10_padj,
      size = overlap_count,
      color = p_adjust
    ),
    alpha = 0.9
  ) +
  
  scale_color_gradient(
    name = "adjusted p",
    low = "#56B4E9",
    high = "#132B43",
    guide = guide_colorbar(reverse = TRUE)
  ) +
  
  scale_size_continuous(
    name = "Count",
    range = c(2, 7)
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    aes(
      x = bar_module_x,
      fill = module
    ),
    width = 0.28,
    height = 0.92
  ) +
  
  scale_fill_manual(
    name = "Module",
    values = module_colors,
    breaks = module_plot_order,
    drop = FALSE
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    aes(
      x = bar_level_x,
      fill = mercator_level
    ),
    width = 0.28,
    height = 0.92
  ) +
  
  scale_fill_manual(
    name = "Mercator level",
    values = level_colors,
    breaks = level_plot_order,
    drop = FALSE,
    na.value = "white"
  ) +
  
  annotate(
    "text",
    x = bar_module_x,
    y = nrow(plot_df) + 1.0,
    label = "Module",
    angle = 90,
    size = 3
  ) +
  
  annotate(
    "text",
    x = bar_level_x,
    y = nrow(plot_df) + 1.0,
    label = "Level",
    angle = 90,
    size = 3
  ) +
  
  scale_y_discrete(labels = y_labels) +
  
  scale_x_continuous(
    breaks = x_breaks,
    limits = c(0, x_limit_max),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  labs(
    x = expression(-log[10]("adjusted p")),
    y = "Mercator Terms",
    title = "Mercator enrichment of selected WGCNA modules"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5.5, 25, 5.5, 5.5)
  )

plot_height <- max(6, 0.32 * nrow(plot_df) + 2)
plot_height <- min(plot_height, 18)

ggsave(
  filename = file.path(outdir, "Selected_modules_Mercator_reference_style_plot_reordered.pdf"),
  plot = p_reference,
  width = 13,
  height = plot_height,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "Selected_modules_Mercator_reference_style_plot_reordered.png"),
  plot = p_reference,
  width = 13,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)

print(p_reference)
## ============================================================
## 11. Save top terms table
## ============================================================

top_terms_table <- plot_df %>%
  dplyr::select(
    module,
    plot_reason,
    mercator_level,
    mercator_term,
    overlap_count,
    term_size,
    p_value,
    p_adjust,
    neg_log10_padj,
    enrichment_ratio,
    overlap_gene_ids
  )

write.csv(
  top_terms_table,
  file.path(outdir, "Selected_modules_plotted_Mercator_terms.csv"),
  row.names = FALSE
)

cat("\nPlotted terms:\n")
print(top_terms_table, n = 100)

## ============================================================
## 12. Optional: check immunity/defense terms in selected modules
## ============================================================

immune_keywords <- c(
  "immun",
  "immune",
  "defen",
  "defense",
  "defence",
  "pathogen",
  "biotic",
  "disease",
  "resistance",
  "resistant",
  "NLR",
  "TIR",
  "NB-LRR",
  "LRR",
  "receptor",
  "effector",
  "salicylic",
  "jasmon",
  "ethylene",
  "WRKY",
  "pathogenesis",
  "PR protein",
  "hypersensitive",
  "oxidative burst",
  "ROS",
  "response to stress",
  "external stimuli",
  "stimulus"
)

immune_pattern <- paste(immune_keywords, collapse = "|")

immune_selected <- enrichment_selected %>%
  dplyr::filter(
    stringr::str_detect(
      mercator_term,
      stringr::regex(immune_pattern, ignore_case = TRUE)
    )
  ) %>%
  dplyr::arrange(p_adjust, p_value, dplyr::desc(overlap_count))

immune_selected_sig <- immune_selected %>%
  dplyr::filter(
    p_value <= PV_CUT,
    p_adjust <= Q_CUT
  )

write.csv(
  immune_selected,
  file.path(outdir, "Selected_modules_immunity_defense_like_all_results.csv"),
  row.names = FALSE
)

write.csv(
  immune_selected_sig,
  file.path(outdir, "Selected_modules_immunity_defense_like_SIGNIFICANT_results.csv"),
  row.names = FALSE
)

cat("\nSignificant immunity/defense-like terms in selected modules:\n")
print(
  immune_selected_sig %>%
    dplyr::select(
      module,
      mercator_level,
      mercator_term,
      overlap_count,
      term_size,
      p_value,
      p_adjust,
      enrichment_ratio,
      overlap_gene_ids
    ),
  n = 100
)

cat("\nDone. Outputs saved in:\n")
print(outdir)
