############################################################
## Epigenetic + shared GO/Mercator plot
## Corrected remembered logic:
## 1. Term must be present in an individual genotype-treatment
##    for at least 2 generations.
## 2. From those recurrent terms, keep only terms shared
##    between WT and gsnor, irrespective of treatment.
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(scales)
})

set.seed(42)

## =========================================================
## 1. Paths
## =========================================================

ROOT_DIR <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

combined_csv <- file.path(
  ROOT_DIR,
  "DESeq_split_up_down",
  "Mercator_result_standardized",
  "ALL_samples_Mercator_clusterprof_combined.csv"
)

out_dir <- file.path(ROOT_DIR, "GO_Mercator_epigenetic_shared_recurrent_logic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(combined_csv)) {
  stop("Combined enrichment file not found:\n", combined_csv)
}

all_results <- read_csv(combined_csv, show_col_types = FALSE)

## =========================================================
## 2. Settings
## =========================================================

# Same filtering as your old plotting code
PLOT_Q_CUT <- 0.20
PLOT_P_CUT <- 0.05

genotype_order <- c("WT", "gsnor")

treatment_order <- c(
  "Drought", "CO2", "Ozone", "Heat", "Combination"
)

generation_order <- c("Gen1", "Gen3", "Gen5")
direction_order <- c("up", "down")
levels_vec <- paste0("level", 1:8)

level_cols <- c(
  "level1" = "#1B9E77",
  "level2" = "#E66101",
  "level3" = "#7570B3",
  "level4" = "#E7298A",
  "level5" = "#66A61E",
  "level6" = "#E6AB02",
  "level7" = "#A6761D",
  "level8" = "#666666"
)

## =========================================================
## 3. Epigenetic terms from old plot
## =========================================================

epigenetic_terms <- c(
  "Coenzyme metabolism",
  "Plant organogenesis",
  "Amino acid metabolism",
  "RNA modification",
  "chromatin structure",
  "SnRK1-kinase regulatory system",
  "RNA transport",
  "DNA wrapping",
  "helix-turn-helix DNA-binding domain",
  "aspartate-derived amino acids",
  "histone H3-H4 tetramer",
  "histone H2A-H2B dimer",
  "CULIN-based E3 ubiquitin ligase activities",
  "Ring E3 ubiquitin ligase activities",
  "histone (H2B/HTB)",
  "histone (H3)"
)

# Only used for ordering if those terms appear.
# Extra shared terms found by the logic are added below.
shared_terms_preferred_order <- c(
  "Photosynthesis",
  "Nutrient uptake",
  "photophosphorylation",
  "phenolics",
  "S-glutathionylation",
  "photosystem II",
  "mitochondrial ribosome biogenesis",
  "beta-hairpin exposed by alpha/beta-scaffold structure",
  "photosystem I",
  "flavonoid biosynthesis",
  "glutaredoxin",
  "regulatory protein (FLZ) of SnRK1 complex",
  "WRKY transcription factor activity",
  "small ribosomal subunit proteome",
  "LHC-II complex",
  "PS-I complex",
  "AP2/ERF family",
  "transcription factor (WRKY)",
  "ERF subfamily",
  "subgroup ERF-IX transcription factor"
)

## =========================================================
## 4. Helper functions
## =========================================================

term_norm <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u2013|\u2014|\u2212", "-")
  x <- str_replace_all(x, "−", "-")
  x <- tolower(x)
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  x <- str_squish(x)
  x
}

safe_neglog10 <- function(x) {
  -log10(pmax(as.numeric(x), 1e-300))
}

display_treatment <- function(x) {
  case_when(
    x == "CO2" ~ "CO2",
    x == "Ozone" ~ "O3",
    x == "Combination" ~ "Combi",
    TRUE ~ x
  )
}

generation_short <- function(x) {
  case_when(
    x == "Gen1" ~ "G1",
    x == "Gen3" ~ "G3",
    x == "Gen5" ~ "G5",
    TRUE ~ x
  )
}

make_size_breaks <- function(max_count) {
  if (is.na(max_count) || max_count <= 10) return(c(5, 10))
  if (max_count <= 30) return(c(10, 20, 30))
  if (max_count <= 80) return(c(25, 50, 75))
  if (max_count <= 120) return(c(25, 50, 75, 100))
  c(50, 100, 150, 200)
}

## =========================================================
## 5. Clean enrichment table
## =========================================================

needed_cols <- c(
  "Description", "p.adjust", "pvalue", "Count",
  "sample_id", "genotype", "generation", "treatment",
  "Sets", "Level"
)

missing_cols <- setdiff(needed_cols, colnames(all_results))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns:\n",
    paste(missing_cols, collapse = ", ")
  )
}

sig_results <- all_results %>%
  mutate(
    Description = as.character(Description),
    term_norm = term_norm(Description),
    p.adjust = as.numeric(p.adjust),
    pvalue = as.numeric(pvalue),
    Count = as.numeric(Count),
    sample_id = as.character(sample_id),
    genotype = as.character(genotype),
    generation = as.character(generation),
    treatment = as.character(treatment),
    Sets = as.character(Sets),
    Level = as.character(Level),
    direction = case_when(
      str_detect(Sets, "_up$") ~ "up",
      str_detect(Sets, "_down$") ~ "down",
      TRUE ~ NA_character_
    ),
    neg_log10_padj = safe_neglog10(p.adjust)
  ) %>%
  filter(
    !is.na(Description),
    !is.na(term_norm),
    !is.na(direction),
    !is.na(p.adjust),
    !is.na(pvalue),
    p.adjust <= PLOT_Q_CUT,
    pvalue <= PLOT_P_CUT,
    genotype %in% genotype_order,
    treatment %in% treatment_order,
    generation %in% generation_order,
    direction %in% direction_order
  )

## =========================================================
## 6. Epigenetic terms
## =========================================================

epigenetic_norm <- term_norm(epigenetic_terms)

epigenetic_hits <- sig_results %>%
  filter(term_norm %in% epigenetic_norm) %>%
  mutate(term_group = "Epigenetic")

## =========================================================
## 7. Corrected shared-term logic
## =========================================================
# Step 1:
# For each GO term + genotype + treatment,
# check if it is present in at least 2 generations.
#
# Step 2:
# From those recurrent genotype-treatment terms,
# keep terms that are recurrent in both WT and gsnor.
#
# Treatment does not need to be the same between genotypes.

recurrent_by_genotype_treatment <- sig_results %>%
  filter(!term_norm %in% epigenetic_norm) %>%
  distinct(
    term_norm,
    Description,
    Level,
    genotype,
    treatment,
    generation
  ) %>%
  group_by(term_norm, genotype, treatment) %>%
  summarise(
    n_generations_present = n_distinct(generation),
    generations_present = paste(sort(unique(generation)), collapse = ";"),
    descriptions_seen = paste(sort(unique(Description)), collapse = " | "),
    levels_seen = paste(sort(unique(Level)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_generations_present >= 2)

shared_recurrent_terms <- recurrent_by_genotype_treatment %>%
  group_by(term_norm) %>%
  summarise(
    n_genotypes_recurrent = n_distinct(genotype),
    genotypes_recurrent = paste(sort(unique(genotype)), collapse = ";"),
    recurrent_contexts = paste(
      paste(genotype, treatment, generations_present, sep = ":"),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_genotypes_recurrent == 2,
    str_detect(genotypes_recurrent, "WT"),
    str_detect(genotypes_recurrent, "gsnor")
  )

shared_terms_norm <- unique(shared_recurrent_terms$term_norm)

# Plot all significant occurrences of selected shared terms,
# not only the recurrent treatment where they were selected.
shared_hits <- sig_results %>%
  filter(term_norm %in% shared_terms_norm) %>%
  filter(!term_norm %in% epigenetic_norm) %>%
  mutate(term_group = "Shared")

## =========================================================
## 8. Save selection summaries
## =========================================================

write_csv(
  recurrent_by_genotype_treatment,
  file.path(out_dir, "Step1_terms_recurrent_in_individual_genotype_treatment.csv")
)

write_csv(
  shared_recurrent_terms,
  file.path(out_dir, "Step2_recurrent_terms_shared_between_WT_and_gsnor.csv")
)

write_csv(
  shared_hits,
  file.path(out_dir, "Shared_hits_all_occurrences_plotted.csv")
)

## =========================================================
## 9. Combine plot terms
## =========================================================

plot_hits <- bind_rows(epigenetic_hits, shared_hits)

if (nrow(plot_hits) == 0) {
  stop("No epigenetic/shared terms found after filtering.")
}

write_csv(
  plot_hits,
  file.path(out_dir, "Epigenetic_and_shared_terms_used_for_plot.csv")
)

cat("\nNumber of selected terms:\n")
cat("Epigenetic terms found:", n_distinct(epigenetic_hits$term_norm), "\n")
cat("Shared recurrent terms found:", n_distinct(shared_hits$term_norm), "\n")

## =========================================================
## 10. X-axis layout
## =========================================================

x_lookup <- expand_grid(
  genotype = genotype_order,
  treatment = treatment_order,
  generation = generation_order
) %>%
  mutate(
    genotype = factor(genotype, levels = genotype_order),
    treatment = factor(treatment, levels = treatment_order),
    generation = factor(generation, levels = generation_order)
  ) %>%
  arrange(genotype, treatment, generation) %>%
  mutate(
    sample_id = paste(genotype, treatment, generation, sep = "_"),
    x_pos = row_number(),
    x_label = paste(
      display_treatment(as.character(treatment)),
      generation_short(as.character(generation))
    )
  )

wt_range <- x_lookup %>%
  filter(genotype == "WT") %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    xmid = mean(range(x_pos))
  )

gs_range <- x_lookup %>%
  filter(genotype == "gsnor") %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    xmid = mean(range(x_pos))
  )

separator_x <- mean(c(wt_range$xmax, gs_range$xmin))

## =========================================================
## 11. Term order
## =========================================================

epigenetic_order_norm <- term_norm(epigenetic_terms)
shared_order_norm <- term_norm(shared_terms_preferred_order)

term_order_table <- bind_rows(
  tibble(
    term_group = "Epigenetic",
    term_norm = epigenetic_order_norm,
    preferred_order = seq_along(epigenetic_order_norm)
  ),
  tibble(
    term_group = "Shared",
    term_norm = shared_order_norm,
    preferred_order = seq_along(shared_order_norm)
  )
)

term_summary <- plot_hits %>%
  mutate(
    term_key = paste(term_group, Level, term_norm, sep = "__")
  ) %>%
  group_by(term_group, term_key, term_norm, Description, Level) %>%
  summarise(
    best_score = max(neg_log10_padj, na.rm = TRUE),
    max_count = max(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    term_order_table,
    by = c("term_group", "term_norm")
  ) %>%
  mutate(
    term_group = factor(term_group, levels = c("Epigenetic", "Shared")),
    level_num = as.integer(str_remove(Level, "level"))
  ) %>%
  arrange(
    term_group,
    is.na(preferred_order),
    preferred_order,
    level_num,
    desc(best_score),
    desc(max_count),
    Description
  ) %>%
  mutate(
    term_label = Description
  )

dup_labels <- term_summary$term_label[
  duplicated(term_summary$term_label) |
    duplicated(term_summary$term_label, fromLast = TRUE)
]

term_summary <- term_summary %>%
  mutate(
    term_label = ifelse(
      term_label %in% dup_labels,
      paste0(term_label, " (", Level, ")"),
      term_label
    )
  )

term_levels <- rev(term_summary$term_key)

term_label_lookup <- term_summary$term_label
names(term_label_lookup) <- term_summary$term_key

## =========================================================
## 12. Plot data
## =========================================================

actual_points <- plot_hits %>%
  mutate(
    term_key = paste(term_group, Level, term_norm, sep = "__")
  ) %>%
  filter(term_key %in% term_summary$term_key) %>%
  group_by(
    term_group,
    term_key,
    Description,
    Level,
    sample_id,
    genotype,
    treatment,
    generation,
    direction
  ) %>%
  summarise(
    Count = max(Count, na.rm = TRUE),
    p.adjust = min(p.adjust, na.rm = TRUE),
    pvalue = min(pvalue, na.rm = TRUE),
    neg_log10_padj = max(neg_log10_padj, na.rm = TRUE),
    .groups = "drop"
  )

plot_points <- expand_grid(
  term_key = term_summary$term_key,
  sample_id = x_lookup$sample_id,
  direction = direction_order
) %>%
  left_join(
    term_summary %>%
      select(term_group, term_key, Description, Level),
    by = "term_key"
  ) %>%
  left_join(
    actual_points %>%
      select(
        term_key,
        sample_id,
        direction,
        Count,
        p.adjust,
        pvalue,
        neg_log10_padj
      ),
    by = c("term_key", "sample_id", "direction")
  ) %>%
  left_join(
    x_lookup %>%
      select(sample_id, x_pos, x_label, genotype, treatment, generation),
    by = "sample_id"
  ) %>%
  mutate(
    Count = as.numeric(Count),
    Count_real = Count,
    Count_for_size = ifelse(is.na(Count), 1, Count),
    neg_log10_padj = ifelse(is.na(neg_log10_padj), 0, neg_log10_padj),
    term_group = factor(term_group, levels = c("Epigenetic", "Shared")),
    term_key = factor(term_key, levels = term_levels),
    direction = factor(direction, levels = direction_order)
  )

strip_df <- term_summary %>%
  mutate(
    term_group = factor(term_group, levels = c("Epigenetic", "Shared")),
    term_key = factor(term_key, levels = term_levels),
    Level = factor(Level, levels = levels_vec),
    x_pos = max(x_lookup$x_pos) + 0.70
  )

header_df <- tibble(
  x_pos = c(wt_range$xmid, gs_range$xmid),
  label = c("WT", "italic(gsnor)*'-ko'")
)

max_col <- max(plot_points$neg_log10_padj, na.rm = TRUE)
if (!is.finite(max_col) || max_col == 0) max_col <- 1

max_count <- max(actual_points$Count, na.rm = TRUE)
size_breaks <- make_size_breaks(max_count)

## =========================================================
## 13. Plot
## =========================================================

p_recurrent_shared <- ggplot() +
  
  geom_vline(
    xintercept = separator_x,
    color = "#B22222",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = plot_points,
    aes(
      x = x_pos,
      y = term_key,
      color = neg_log10_padj,
      size = Count_for_size,
      shape = direction
    ),
    alpha = 0.95,
    position = position_dodge(width = 0.35)
  ) +
  
  geom_tile(
    data = strip_df,
    aes(
      x = x_pos,
      y = term_key,
      fill = Level
    ),
    width = 0.22,
    height = 0.85
  ) +
  
  geom_text(
    data = header_df,
    aes(
      x = x_pos,
      y = Inf,
      label = label
    ),
    parse = TRUE,
    inherit.aes = FALSE,
    vjust = -0.70,
    fontface = "bold",
    size = 5.3
  ) +
  
  annotate(
    "text",
    x = max(x_lookup$x_pos) + 1.18,
    y = -Inf,
    label = "Mercator levels",
    angle = 55,
    hjust = 0,
    vjust = -0.2,
    color = "red3",
    size = 3.5
  ) +
  
  facet_grid(
    term_group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_y_discrete(
    labels = term_label_lookup
  ) +
  
  scale_x_continuous(
    breaks = x_lookup$x_pos,
    labels = x_lookup$x_label,
    limits = c(0.5, max(x_lookup$x_pos) + 1.55),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  
  scale_color_gradientn(
    colors = c("#1f0088", "#8c1ab1", "#e66f61", "#ffe51f"),
    limits = c(0, max_col),
    name = expression(-log[10]("adjusted p"))
  ) +
  
  scale_size_continuous(
    range = c(0.25, 7.5),
    breaks = size_breaks,
    name = "Gene Count"
  ) +
  
  scale_shape_manual(
    values = c(
      "up" = 16,
      "down" = 17
    ),
    labels = c(
      "up" = "Up-regulated",
      "down" = "Down-regulated"
    ),
    name = "Direction"
  ) +
  
  scale_fill_manual(
    values = level_cols,
    breaks = levels_vec[levels_vec %in% unique(as.character(strip_df$Level))],
    drop = FALSE,
    name = "Mercator level"
  ) +
  
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 3),
    fill = guide_legend(order = 4)
  ) +
  
  labs(
    x = "Treatment and Generation",
    y = "Mercator Terms"
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    panel.grid.major.x = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey35",
      fill = NA,
      linewidth = 0.55
    ),
    
    strip.background = element_rect(
      fill = "grey85",
      color = "grey55",
      linewidth = 0.5
    ),
    strip.text.y.left = element_text(
      size = 13,
      face = "bold",
      angle = 90
    ),
    
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 8.4
    ),
    axis.text.y = element_text(
      color = "black",
      size = 9.2
    ),
    axis.title.x = element_text(
      size = 14,
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 14,
      face = "bold",
      margin = margin(r = 8)
    ),
    
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    
    plot.margin = margin(35, 35, 15, 15)
  )

print(p_recurrent_shared)

## =========================================================
## 14. Save
## =========================================================

n_terms_total <- nrow(term_summary)
fig_height <- max(8, min(18, 4 + 0.22 * n_terms_total))

ggsave(
  file.path(out_dir, "Epigenetic_shared_recurrent_GO_Mercator_plot.pdf"),
  p_recurrent_shared,
  width = 13,
  height = fig_height,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(out_dir, "Epigenetic_shared_recurrent_GO_Mercator_plot.png"),
  p_recurrent_shared,
  width = 13,
  height = fig_height,
  dpi = 300,
  limitsize = FALSE
)

cat("\nDone.\n")
cat("Plot saved in:\n")
cat(out_dir, "\n\n")

cat("Number of plotted terms:\n")
print(
  term_summary %>%
    count(term_group, name = "n_terms")
)

cat("\nShared recurrent terms selected:\n")
print(shared_recurrent_terms)

cat("\nRecurrent genotype-treatment contexts:\n")
print(recurrent_by_genotype_treatment)

