rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(scales)
})

# ============================================================
# 1. Input
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

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
  
  if (!"gene_id" %in% colnames(dat)) {
    dat$gene_id <- rownames(dat)
  }
  
  dat
}

# ============================================================
# 3. Detect treatment ONLY from full result files
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
# 4. Count up/down DEGs from full DESeq2 files only
# ============================================================

deg_counts <- list()
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
  
  # ----------------------------------------------------------
  # Keep ONLY full DESeq2 result files:
  # Drought.txt, CO2.txt, Ozone.txt, Heat.txt, Combination.txt
  # Exclude pvalues, upregulated, downregulated files
  # ----------------------------------------------------------
  
  files <- all_files[
    !str_detect(
      tolower(basename(all_files)),
      "pval|pvalue|upregulated|downregulated"
    )
  ]
  
  if (length(files) == 0) {
    warning("No full DESeq2 result files found in: ", folder_path)
    next
  }
  
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
      warning(
        "Skipping file because required columns are missing: ",
        basename(file)
      )
      next
    }
    
    dat <- dat %>%
      mutate(
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj)
      )
    
    sig <- dat %>%
      filter(
        !is.na(padj),
        !is.na(log2FoldChange),
        padj < padj_cutoff,
        abs(log2FoldChange) >= lfc_cutoff
      )
    
    n_up <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
    n_down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
    
    deg_counts[[length(deg_counts) + 1]] <- data.frame(
      genotype = genotype,
      generation = generation,
      treatment = treatment,
      file = basename(file),
      n_up = n_up,
      n_down = n_down,
      stringsAsFactors = FALSE
    )
  }
}

deg_summary_raw <- bind_rows(deg_counts)
file_map_df <- bind_rows(file_map)

cat("\nDetected full DESeq2 files:\n")
print(file_map_df)

if (nrow(deg_summary_raw) == 0) {
  stop("No DEG count data were generated. Please check folder paths and file names.")
}

# ============================================================
# 5. Duplicate check
# ============================================================

duplicates <- deg_summary_raw %>%
  count(genotype, generation, treatment) %>%
  filter(n > 1)

if (nrow(duplicates) > 0) {
  cat("\nDuplicate full treatment files detected:\n")
  print(duplicates)
  
  cat("\nFiles causing duplicates:\n")
  print(
    deg_summary_raw %>%
      semi_join(
        duplicates,
        by = c("genotype", "generation", "treatment")
      ) %>%
      arrange(genotype, generation, treatment)
  )
  
  stop("Duplicate full treatment files still detected. Check file names.")
}

# ============================================================
# 6. Complete missing combinations with zeros
# ============================================================

deg_summary <- deg_summary_raw %>%
  mutate(
    genotype = factor(genotype, levels = genotype_order),
    generation = factor(generation, levels = generation_order),
    treatment = factor(treatment, levels = treatment_order)
  ) %>%
  complete(
    genotype,
    generation,
    treatment,
    fill = list(
      n_up = 0,
      n_down = 0,
      file = NA_character_
    )
  ) %>%
  arrange(
    genotype,
    treatment,
    generation
  )

cat("\nDEG summary:\n")
print(deg_summary)

cat("\nCombination rows check:\n")
print(
  deg_summary %>%
    filter(treatment == "Combination")
)

# ============================================================
# 7. Prepare plotting dataframe
# ============================================================

plot_df <- deg_summary %>%
  mutate(
    sample_label = case_when(
      genotype == "WT" ~ paste0("WT", treatment, "_", generation),
      genotype == "gsnor-ko" ~ paste0("gs", treatment, "_", generation)
    )
  ) %>%
  select(
    genotype,
    generation,
    treatment,
    sample_label,
    n_up,
    n_down
  ) %>%
  pivot_longer(
    cols = c(n_down, n_up),
    names_to = "direction",
    values_to = "count"
  ) %>%
  mutate(
    direction = recode(
      direction,
      "n_down" = "Down",
      "n_up" = "Up"
    ),
    signed_count = ifelse(direction == "Down", -count, count)
  )

label_order <- deg_summary %>%
  mutate(
    sample_label = case_when(
      genotype == "WT" ~ paste0("WT", treatment, "_", generation),
      genotype == "gsnor-ko" ~ paste0("gs", treatment, "_", generation)
    )
  ) %>%
  arrange(
    genotype,
    treatment,
    generation
  ) %>%
  pull(sample_label)

plot_df$sample_label <- factor(
  plot_df$sample_label,
  levels = rev(unique(label_order))
)

max_abs <- max(abs(plot_df$signed_count), na.rm = TRUE)
label_offset <- max_abs * 0.035

plot_df <- plot_df %>%
  mutate(
    label_x = case_when(
      direction == "Up" & count > 0 ~ signed_count + label_offset,
      direction == "Down" & count > 0 ~ signed_count - label_offset,
      TRUE ~ signed_count
    ),
    label_hjust = case_when(
      direction == "Up" ~ 0,
      direction == "Down" ~ 1,
      TRUE ~ 0.5
    ),
    label_text = ifelse(count > 0, as.character(count), "")
  )

# ============================================================
# 8. Plot
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = signed_count,
    y = sample_label,
    fill = direction
  )
) +
  geom_col(width = 0.75) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    aes(
      x = label_x,
      label = label_text,
      hjust = label_hjust
    ),
    size = 3.2,
    colour = "black"
  ) +
  facet_grid(
    genotype ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Down" = "#2C7FB8",
      "Up" = "#E41A1C"
    ),
    breaks = c("Down", "Up")
  ) +
  scale_x_continuous(
    limits = c(-max_abs * 1.25, max_abs * 1.25),
    breaks = scales::pretty_breaks(n = 5),
    labels = abs
  ) +
  labs(
    x = "Gene number (|log2FC| ≥ 0)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.y = element_text(
      angle = 270,
      face = "bold",
      size = 11
    ),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 10)
  )

p

#saved 10X6 portrait

# ============================================================
# 9. Save outputs
# ============================================================

write.csv(
  file_map_df,
  file.path(root_dir, "DEG_file_detection_map_FULL_FILES_ONLY.csv"),
  row.names = FALSE
)

write.csv(
  deg_summary,
  file.path(root_dir, "DEG_up_down_counts_summary.csv"),
  row.names = FALSE
)

ggsave(
  file.path(root_dir, "DEG_up_down_diverging_barplot.pdf"),
  p,
  width = 7,
  height = 9
)

ggsave(
  file.path(root_dir, "DEG_up_down_diverging_barplot.png"),
  p,
  width = 7,
  height = 9,
  dpi = 300
)