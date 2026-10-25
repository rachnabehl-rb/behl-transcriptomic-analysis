## ============================================================
## CFU Panel B plot from Excel sheets G1, G3, G5
## Same logic as old code:
## - log10 transformation
## - Gaussian GLM on log10_value ~ group
## - Tukey compact letter display
## - ordered groups:
##   WT Control, WT Combination, gsnor Control, gsnor Combination
## ============================================================

rm(list = ls())
set.seed(798)

## ============================================================
## Load libraries
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(car)
  library(emmeans)
  library(multcomp)
  library(stringr)
  library(patchwork)
})

## ============================================================
## Input and output
## ============================================================

## Change this path if needed
excel_file <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Supplementary File S22.xlsx"

out_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/CFU_panel_B_old_logic"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Expected sheet names
generation_order <- c("G1", "G3", "G5")

## ============================================================
## Read Excel sheets
## ============================================================

available_sheets <- readxl::excel_sheets(excel_file)

cat("Sheets found:\n")
print(available_sheets)

sheets_use <- intersect(generation_order, available_sheets)

if (length(sheets_use) != 3) {
  stop(
    "Expected sheets named G1, G3, G5. Found matching sheets: ",
    paste(sheets_use, collapse = ", ")
  )
}

all_data <- map_dfr(sheets_use, function(sh) {
  
  dat <- readxl::read_excel(excel_file, sheet = sh)
  
  ## Make column names easier to handle
  colnames(dat) <- tolower(str_squish(colnames(dat)))
  
  ## Required columns
  required_cols <- c("genotype", "treatment", "value")
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "Sheet ", sh, " is missing columns: ",
      paste(missing_cols, collapse = ", "),
      "\nDetected columns are: ",
      paste(colnames(dat), collapse = ", ")
    )
  }
  
  dat %>%
    mutate(generation = sh)
})

## ============================================================
## Clean data
## ============================================================

final.data <- all_data %>%
  mutate(
    genotype  = str_squish(as.character(genotype)),
    treatment = str_squish(as.character(treatment)),
    value     = as.numeric(value),
    generation = as.character(generation),
    
    genotype = case_when(
      str_detect(tolower(genotype), "wt") ~ "WT",
      str_detect(tolower(genotype), "gsnor") ~ "gsnor",
      TRUE ~ genotype
    ),
    
    treatment = case_when(
      str_detect(tolower(treatment), "control") ~ "Control",
      str_detect(tolower(treatment), "combi|combination") ~ "Combination",
      TRUE ~ treatment
    ),
    
    genotype = factor(genotype, levels = c("WT", "gsnor")),
    treatment = factor(treatment, levels = c("Control", "Combination")),
    generation = factor(generation, levels = generation_order)
  ) %>%
  filter(
    !is.na(genotype),
    !is.na(treatment),
    !is.na(value),
    !is.na(generation)
  )

## Check cleaned names
cat("\nCleaned genotype values:\n")
print(unique(final.data$genotype))

cat("\nCleaned treatment values:\n")
print(unique(final.data$treatment))

cat("\nCounts per group:\n")
print(
  final.data %>%
    count(generation, genotype, treatment)
)

## ============================================================
## Log10 transformation
## Same rule as old code:
## if zero values exist, use log10(value + 1)
## ============================================================

if (any(final.data$value <= 0, na.rm = TRUE)) {
  
  final.data <- final.data %>%
    mutate(log10_value = log10(value + 1))
  
  y_lab <- expression(log[10](CFU + 1) / leaf~disk)
  
} else {
  
  final.data <- final.data %>%
    mutate(log10_value = log10(value))
  
  y_lab <- expression(log[10](CFU) / leaf~disk)
}

## ============================================================
## Create ordered group variable
## ============================================================

final.data <- final.data %>%
  mutate(
    group = paste(genotype, treatment),
    group = factor(
      group,
      levels = c(
        "WT Control",
        "WT Combination",
        "gsnor Control",
        "gsnor Combination"
      )
    )
  )

cat("\nGroup order:\n")
print(levels(final.data$group))

write.csv(
  final.data,
  file.path(out_dir, "CFU_cleaned_log10_data.csv"),
  row.names = FALSE
)

## ============================================================
## GLM + Tukey CLD per generation
## Same old logic:
## Gaussian GLM on log10-transformed values
## ============================================================

run_glm_one_generation <- function(dat_gen) {
  
  gen_name <- unique(as.character(dat_gen$generation))
  
  dat_gen <- dat_gen %>%
    mutate(
      group = factor(
        group,
        levels = c(
          "WT Control",
          "WT Combination",
          "gsnor Control",
          "gsnor Combination"
        )
      )
    )
  
  ## Optional assumption checks
  shapiro_overall <- tryCatch(
    shapiro.test(dat_gen$log10_value)$p.value,
    error = function(e) NA_real_
  )
  
  shapiro_group <- dat_gen %>%
    group_by(group) %>%
    summarise(
      n = n(),
      shapiro_p = tryCatch(
        shapiro.test(log10_value)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  levene_p <- tryCatch(
    car::leveneTest(log10_value ~ group, data = dat_gen)[["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  
  ## GLM on log10 values, same as old code
  glm_fit <- glm(
    log10_value ~ group,
    family = gaussian(link = "identity"),
    data = dat_gen
  )
  
  anova_tbl <- car::Anova(
    glm_fit,
    type = 2,
    test.statistic = "F"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(generation = gen_name)
  
  emm <- emmeans(glm_fit, ~ group)
  
  pairwise_tbl <- pairs(emm, adjust = "tukey") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(generation = gen_name)
  
  letters_df <- multcomp::cld(
    emm,
    adjust = "tukey",
    Letters = letters,
    reverse = TRUE,
    alpha = 0.05
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      .group = str_remove_all(.group, "\\s+"),
      group = factor(
        group,
        levels = c(
          "WT Control",
          "WT Combination",
          "gsnor Control",
          "gsnor Combination"
        )
      ),
      generation = gen_name
    ) %>%
    dplyr::select(generation, group, .group, everything())
  
  assumption_tbl <- tibble(
    generation = gen_name,
    shapiro_overall_p = shapiro_overall,
    levene_p = levene_p
  )
  
  list(
    glm_fit = glm_fit,
    anova = anova_tbl,
    pairwise = pairwise_tbl,
    letters = letters_df,
    assumptions = assumption_tbl,
    shapiro_group = shapiro_group %>% mutate(generation = gen_name)
  )
}

glm_results <- final.data %>%
  group_split(generation) %>%
  map(run_glm_one_generation)

letters_df <- map_dfr(glm_results, "letters")
anova_df <- map_dfr(glm_results, "anova")
pairwise_df <- map_dfr(glm_results, "pairwise")
assumption_df <- map_dfr(glm_results, "assumptions")
shapiro_group_df <- map_dfr(glm_results, "shapiro_group")

write.csv(
  letters_df,
  file.path(out_dir, "GLM_Tukey_CLD_letters_by_generation.csv"),
  row.names = FALSE
)

write.csv(
  anova_df,
  file.path(out_dir, "GLM_Anova_by_generation.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_df,
  file.path(out_dir, "GLM_Tukey_pairwise_by_generation.csv"),
  row.names = FALSE
)

write.csv(
  assumption_df,
  file.path(out_dir, "GLM_assumption_summary_by_generation.csv"),
  row.names = FALSE
)

write.csv(
  shapiro_group_df,
  file.path(out_dir, "Shapiro_by_group_generation.csv"),
  row.names = FALSE
)

cat("\nCLD letters:\n")
print(
  letters_df %>%
    dplyr::select(generation, group, .group)
)

## ============================================================
## Label positions for letters
## ============================================================

label_df <- final.data %>%
  group_by(generation, group) %>%
  summarise(
    y_pos = max(log10_value, na.rm = TRUE) + 0.25,
    .groups = "drop"
  ) %>%
  left_join(
    letters_df %>% dplyr::select(generation, group, .group),
    by = c("generation", "group")
  ) %>%
  mutate(
    generation = factor(generation, levels = generation_order),
    group = factor(
      group,
      levels = c(
        "WT Control",
        "WT Combination",
        "gsnor Control",
        "gsnor Combination"
      )
    )
  )

write.csv(
  label_df,
  file.path(out_dir, "CLD_label_positions.csv"),
  row.names = FALSE
)

## ============================================================
## Plot
## Similar to old figure:
## - 3 panels G1, G3, G5
## - x-axis: Control Combi Control Combi
## - colour by genotype
## ============================================================

my_fill_colors <- c(
  "WT" = "cyan2",
  "gsnor" = "#FFA07F"
)

my_point_colors <- c(
  "WT" = "cyan3",
  "gsnor" = "#F08A6A"
)

x_labels <- c(
  "WT Control" = "Control",
  "WT Combination" = "Combi",
  "gsnor Control" = "Control",
  "gsnor Combination" = "Combi"
)

## Use same y-limits as your old code, but adapt if needed
y_lower <- 2.5
y_upper <- max(9.5, max(label_df$y_pos, na.rm = TRUE) + 0.3)

p <- ggplot(
  final.data,
  aes(
    x = group,
    y = log10_value,
    fill = genotype
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.75,
    color = "black",
    linewidth = 0.65
  ) +
  geom_jitter(
    aes(color = genotype),
    width = 0.12,
    height = 0,
    size = 2.2,
    alpha = 0.75
  ) +
  geom_text(
    data = label_df,
    aes(
      x = group,
      y = y_pos,
      label = .group
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "red"
  ) +
  facet_wrap(
    ~ generation,
    nrow = 1
  ) +
  scale_fill_manual(
    values = my_fill_colors,
    name = "Genotype",
    labels = c(
      "WT" = "WT",
      "gsnor" = expression(italic(gsnor) * "-ko")
    )
  ) +
  scale_color_manual(
    values = my_point_colors,
    guide = "none"
  ) +
  scale_x_discrete(
    labels = x_labels
  ) +
  coord_cartesian(
    ylim = c(y_lower, y_upper)
  ) +
  theme_classic(base_size = 14) +
  labs(
    tag = "(B)",
    x = "Treatments",
    y = y_lab,
    fill = "Genotype"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 12,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 14,
      color = "black"
    ),
    
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    
    plot.tag = element_text(
      size = 22,
      face = "bold"
    ),
    plot.tag.position = c(0.01, 0.98),
    
    plot.margin = margin(10, 15, 10, 10)
  )

print(p)

## ============================================================
## Save plot
## ============================================================

ggsave(
  filename = file.path(out_dir, "Panel_B_GLM_CFU_boxplot_letters_old_logic.pdf"),
  plot = p,
  width = 12,
  height = 5,
  useDingbats = FALSE
)

ggsave(
  filename = file.path(out_dir, "Panel_B_GLM_CFU_boxplot_letters_old_logic.png"),
  plot = p,
  width = 12,
  height = 5,
  dpi = 300
)

cat("\nDone.\n")
cat("Output saved in:\n")
cat(out_dir, "\n")

