rm(list = ls())

# ============================================================
# 0. Packages
# ============================================================

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  install.packages("ggalluvial")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ============================================================
# 1. Input settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "Gene_fate_alluvial_plots_light")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 0

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

treatment_order <- c("Drought", "CO2", "Ozone", "Heat", "Combination")
generation_order <- c("1", "3", "5")
genotype_order   <- c("WT", "gsnor-ko")

# Colours same style as your figure
fate_cols <- c(
  "Persistent" = "#1b9e77",  # teal
  "Fading"     = "#7570b3",  # purple
  "Emergent"   = "#66a61e",  # green
  "Switcher"   = "#d95f02"   # orange
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
  
  if (!"gene_id" %in% colnames(dat)) {
    dat$gene_id <- rownames(dat)
  }
  
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

# ============================================================
# 3. Read all full DESeq2 result files
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
  generation <- stringr::str_extract(folder, "gen[0-9]+") %>% stringr::str_remove("gen")
  
  all_files <- list.files(
    folder_path,
    pattern = "\\.(txt|csv|tsv)$",
    full.names = TRUE
  )
  
  # use only full DESeq2 files
  files <- all_files[
    !stringr::str_detect(
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
      warning("Skipping file because log2FoldChange/padj missing: ", basename(file))
      next
    }
    
    dat2 <- dat %>%
      mutate(
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj),
        genotype = genotype,
        generation = generation,
        treatment = treatment,
        status = case_when(
          !is.na(padj) & !is.na(log2FoldChange) &
            padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ "Up",
          
          !is.na(padj) & !is.na(log2FoldChange) &
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
          
          TRUE ~ "None"
        )
      ) %>%
      select(gene_id, genotype, generation, treatment, log2FoldChange, padj, status)
    
    all_res_list[[length(all_res_list) + 1]] <- dat2
  }
}

all_res <- bind_rows(all_res_list)
file_map_df <- bind_rows(file_map)

write.csv(
  file_map_df,
  file.path(out_dir, "alluvial_file_detection_map.csv"),
  row.names = FALSE
)

if (nrow(all_res) == 0) {
  stop("No DESeq2 data were read. Check folder paths and file names.")
}

# ============================================================
# 4. Create wide table of G1/G3/G5 status
# ============================================================

status_wide <- all_res %>%
  mutate(
    genotype   = factor(genotype, levels = genotype_order),
    generation = factor(generation, levels = generation_order),
    treatment  = factor(treatment, levels = treatment_order),
    status     = factor(status, levels = c("Down", "None", "Up"))
  ) %>%
  select(gene_id, genotype, treatment, generation, status) %>%
  distinct() %>%
  pivot_wider(
    names_from = generation,
    values_from = status,
    names_prefix = "G",
    values_fill = "None"
  ) %>%
  mutate(
    G1 = as.character(G1),
    G3 = as.character(G3),
    G5 = as.character(G5)
  )

# remove genes that are None in all generations
status_wide <- status_wide %>%
  filter(!(G1 == "None" & G3 == "None" & G5 == "None"))

# ============================================================
# 5. Classify fate
# ============================================================

status_wide <- status_wide %>%
  rowwise() %>%
  mutate(
    has_up   = any(c(G1, G3, G5) == "Up"),
    has_down = any(c(G1, G3, G5) == "Down"),
    
    fate = case_when(
      has_up & has_down ~ "Switcher",
      G1 != "None" & G5 == G1 ~ "Persistent",
      G5 == "None" & (G1 != "None" | G3 != "None") ~ "Fading",
      G1 == "None" & (G3 != "None" | G5 != "None") ~ "Emergent",
      TRUE ~ "Fading"
    )
  ) %>%
  ungroup() %>%
  select(-has_up, -has_down) %>%
  mutate(
    fate = factor(fate, levels = c("Persistent", "Fading", "Emergent", "Switcher"))
  )

# ============================================================
# 6. Labels and plotting order
# ============================================================

status_wide <- status_wide %>%
  mutate(
    genotype_plot = case_when(
      genotype == "WT"       ~ "WT",
      genotype == "gsnor-ko" ~ "gsnor-ko"
    ),
    genotype_plot = factor(genotype_plot, levels = c("WT", "gsnor-ko")),
    
    treatment_plot = case_when(
      treatment == "Drought"     ~ "Drought",
      treatment == "CO2"         ~ "CO2",
      treatment == "Ozone"       ~ "O3",
      treatment == "Heat"        ~ "Heat",
      treatment == "Combination" ~ "Combi"
    ),
    treatment_plot = factor(
      treatment_plot,
      levels = c("Drought", "CO2", "O3", "Heat", "Combi")
    ),
    
    alluvium_id = paste(genotype, treatment, gene_id, sep = "__")
  )

write.csv(
  status_wide,
  file.path(out_dir, "gene_fate_status_wide.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Aggregate trajectories for alluvial plot
# IMPORTANT: use wide alluvia format, not lodes format
# ============================================================

alluvial_plot_df <- status_wide %>%
  count(
    genotype_plot,
    treatment_plot,
    fate,
    G1,
    G3,
    G5,
    name = "n_genes"
  ) %>%
  mutate(
    G1 = factor(G1, levels = c("Up", "None", "Down")),
    G3 = factor(G3, levels = c("Up", "None", "Down")),
    G5 = factor(G5, levels = c("Up", "None", "Down")),
    
    fate = factor(
      fate,
      levels = c("Persistent", "Fading", "Emergent", "Switcher")
    ),
    
    genotype_plot = factor(
      genotype_plot,
      levels = c("WT", "gsnor-ko")
    ),
    
    treatment_plot = factor(
      treatment_plot,
      levels = c("Drought", "CO2", "O3", "Heat", "Combi")
    )
  )

write.csv(
  alluvial_plot_df,
  file.path(out_dir, "gene_fate_alluvial_aggregated_counts.csv"),
  row.names = FALSE
)

cat("\nAggregated alluvial input:\n")
print(alluvial_plot_df)

cat("\nTotal genes plotted:\n")
print(sum(alluvial_plot_df$n_genes))

# ============================================================
# 8. Plot alluvial figure - light background
# ============================================================

p_alluvial <- ggplot(
  alluvial_plot_df,
  aes(
    y = n_genes,
    axis1 = G1,
    axis2 = G3,
    axis3 = G5
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = fate),
    width = 0.16,
    alpha = 0.95,
    color = NA
  ) +
  ggalluvial::geom_stratum(
    width = 0.16,
    fill = "white",
    color = "grey40",
    linewidth = 0.35
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3.6,
    color = "black"
  ) +
  scale_x_discrete(
    limits = c("G1", "G3", "G5"),
    expand = c(0.05, 0.05)
  ) +
  facet_grid(
    treatment_plot ~ genotype_plot,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = fate_cols,
    breaks = c("Persistent", "Fading", "Emergent", "Switcher"),
    drop = FALSE
  ) +
  labs(
    x = "Generations",
    y = "Number of genes",
    fill = "Category:"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "grey95", color = NA),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey60",
      fill = NA,
      linewidth = 0.5
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "grey60",
      linewidth = 0.5
    ),
    strip.text.x = element_text(
      size = 12,
      color = "black"
    ),
    strip.text.y.left = element_text(
      size = 12,
      color = "black",
      angle = 90
    ),
    
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 12,
      color = "black"
    ),
    
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    
    panel.spacing = unit(0.12, "lines"),
    plot.margin = margin(6, 6, 6, 6)
  )

print(p_alluvial)

# ============================================================
# 9. Save plot
# ============================================================

ggsave(
  file.path(out_dir, "Gene_fate_alluvial_light_FIXED.pdf"),
  p_alluvial,
  width = 9,
  height = 8,
  useDingbats = FALSE
)

ggsave(
  file.path(out_dir, "Gene_fate_alluvial_light_FIXED.png"),
  p_alluvial,
  width = 9,
  height = 8,
  dpi = 300
)

cat("\nDone. Fixed alluvial plot saved in:\n")
cat(out_dir, "\n")
