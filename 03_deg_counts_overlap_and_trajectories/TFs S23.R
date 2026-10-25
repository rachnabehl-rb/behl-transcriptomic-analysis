rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ============================================================
# 1. Input settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "TF_heatmap_like_original")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 0

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

genotype_order <- c("WT", "gsnor-ko")

treatment_order <- c(
  "Drought", "CO2", "Ozone", "Heat", "Combination"
)

treatment_label_map <- c(
  "Drought"     = "Drought",
  "CO2"         = "CO\u2082",
  "Ozone"       = "O\u2083",
  "Heat"        = "Heat",
  "Combination" = "Combi"
)

generation_order <- c("G1", "G3", "G5")

# Colour range
fill_limit <- 8

# ============================================================
# 2. TF genes and families from your attached figure
# ============================================================

tf_table <- tibble::tribble(
  ~gene_id,    ~family,
  "AT4G36780", "BES1",
  "AT1G69010", "bHLH",
  "AT3G56970", "bHLH",
  "AT3G56880", "bHLH",
  "AT4G09180", "bHLH",
  "AT4G09820", "bHLH",
  "AT4G16430", "bHLH",
  "AT4G17880", "bHLH",
  "AT1G77920", "bZIP",
  "AT3G10800", "bZIP",
  "AT2G41835", "C2H2",
  "AT3G57000", "C2H2",
  "AT1G25440", "CO-like",
  "AT1G68520", "CO-like",
  "AT1G73870", "CO-like",
  "AT2G24790", "DBB",
  "AT3G15210", "ERF",
  "AT4G31060", "ERF",
  "AT5G61590", "ERF",
  "AT5G61600", "ERF",
  "AT3G10760", "G2-like",
  "AT1G61730", "GeBP",
  "AT4G25210", "GeBP",
  "AT1G14920", "GRAS",
  "AT1G50420", "GRAS",
  "AT3G50650", "GRAS",
  "AT5G66770", "GRAS",
  "AT2G22840", "GRF",
  "AT5G53660", "GRF",
  "AT4G03250", "HB-other",
  "AT2G22800", "HD-ZIP",
  "AT3G61180", "HD-ZIP",
  "AT3G49940", "LBD",
  "AT1G49010", "MYB",
  "AT1G66230", "MYB",
  "AT4G38620", "MYB",
  "AT5G60890", "MYB",
  "AT2G46830", "MYB related",
  "AT5G47390", "MYB related",
  "AT3G01000", "NAC",
  "AT5G05660", "NF-X1",
  "AT4G14540", "NF-YB",
  "AT5G38140", "NF-YC",
  "AT5G08330", "TCP",
  "AT2G24570", "WRKY",
  "AT5G26170", "WRKY",
  "AT5G39760", "ZF-HD"
) %>%
  mutate(
    gene_id = toupper(gene_id),
    gene_order = row_number()
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
      sep = "\t",
      check.names = FALSE,
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

# ============================================================
# 4. Read DESeq2 result files
# ============================================================

all_res_list <- list()
file_map <- list()

for (folder in folders) {
  
  folder_path <- file.path(root_dir, folder)
  
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
  
  # keep only full DESeq2 result files
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
  file.path(out_dir, "TF_heatmap_file_detection_map.csv"),
  row.names = FALSE
)

if (nrow(all_res) == 0) {
  stop("No DESeq2 result data were read.")
}

# ============================================================
# 5. Prepare full plotting grid
# ============================================================

plot_grid <- tidyr::expand_grid(
  gene_id = tf_table$gene_id,
  genotype = genotype_order,
  treatment = treatment_order,
  generation = generation_order
) %>%
  left_join(tf_table, by = "gene_id")

plot_data <- plot_grid %>%
  left_join(
    all_res,
    by = c("gene_id", "genotype", "treatment", "generation")
  ) %>%
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    genotype = factor(genotype, levels = genotype_order),
    treatment = factor(treatment, levels = treatment_order),
    generation = factor(generation, levels = generation_order),
    treatment_label = treatment_label_map[as.character(treatment)]
  )

write.csv(
  plot_data,
  file.path(out_dir, "TF_heatmap_plot_data.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Make manual x positions with gaps like original figure
# ============================================================

x_lookup <- plot_data %>%
  distinct(genotype, treatment, generation, treatment_label) %>%
  mutate(
    genotype_i = as.integer(genotype),
    treatment_i = as.integer(treatment),
    generation_i = as.integer(generation),
    
    # gap between treatment blocks
    treatment_gap = 0.35,
    
    # larger gap between WT and gsnor-ko
    genotype_gap = 1.35,
    
    x_base = case_when(
      genotype == "WT" ~ 0,
      genotype == "gsnor-ko" ~ 
        length(treatment_order) * length(generation_order) +
        (length(treatment_order) - 1) * treatment_gap +
        genotype_gap
    ),
    
    x_pos = x_base +
      (treatment_i - 1) * (length(generation_order) + treatment_gap) +
      generation_i
  ) %>%
  dplyr::select(genotype, treatment, generation, treatment_label, x_pos)

plot_data <- plot_data %>%
  left_join(
    x_lookup,
    by = c("genotype", "treatment", "generation", "treatment_label")
  )

# ============================================================
# 7. Make manual y positions with family gaps
# ============================================================

row_lookup <- tf_table %>%
  mutate(
    family_change = family != lag(family, default = dplyr::first(family)),
    family_block = cumsum(family_change),
    family_gap = family_block * 0.35,
    y_from_top = gene_order + family_gap,
    y_pos = max(y_from_top) - y_from_top + 1
  )

plot_data <- plot_data %>%
  left_join(
    row_lookup %>% dplyr::select(gene_id, y_pos),
    by = "gene_id"
  )

family_labels <- row_lookup %>%
  group_by(family) %>%
  summarise(
    y_pos = mean(y_pos),
    .groups = "drop"
  )

gene_labels <- row_lookup %>%
  dplyr::select(gene_id, family, y_pos)

# ============================================================
# 8. Axis labels and separators
# ============================================================

x_axis_df <- x_lookup %>%
  distinct(genotype, treatment, generation, treatment_label, x_pos) %>%
  arrange(x_pos) %>%
  mutate(
    x_label = paste(treatment_label, generation)
  )

wt_x <- x_lookup %>%
  filter(genotype == "WT") %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    xmid = mean(range(x_pos))
  )

gs_x <- x_lookup %>%
  filter(genotype == "gsnor-ko") %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    xmid = mean(range(x_pos))
  )

separator_x <- mean(c(wt_x$xmax, gs_x$xmin))

right_gene_x   <- max(x_lookup$x_pos) + 0.55
right_family_x <- max(x_lookup$x_pos) + 4.20
right_header_x <- max(x_lookup$x_pos) + 2.60

# ============================================================
# 9. Heatmap plot
# ============================================================

p_tf_heatmap <- ggplot() +
  geom_tile(
    data = plot_data,
    aes(
      x = x_pos,
      y = y_pos,
      fill = log2FoldChange
    ),
    width = 0.96,
    height = 0.88,
    color = "grey25",
    linewidth = 0.25
  ) +
  scale_fill_gradient2(
    low = "#00008B",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-fill_limit, fill_limit),
    oob = scales::squish,
    na.value = "grey90",
    name = NULL
  ) +
  
  # red separator between WT and gsnor-ko
  geom_vline(
    xintercept = separator_x,
    color = "#B22222",
    linewidth = 1.2
  ) +
  
  # gene labels
  geom_text(
    data = gene_labels,
    aes(
      x = right_gene_x,
      y = y_pos,
      label = gene_id
    ),
    hjust = 0,
    size = 3.2,
    color = "black"
  ) +
  
  # family labels
  geom_text(
    data = family_labels,
    aes(
      x = right_family_x,
      y = y_pos,
      label = family
    ),
    hjust = 0,
    size = 3.4,
    color = "black"
  ) +
  
  # thin line between gene IDs and family names
  geom_segment(
    aes(
      x = right_family_x - 0.28,
      xend = right_family_x - 0.28,
      y = min(row_lookup$y_pos) - 0.4,
      yend = max(row_lookup$y_pos) + 0.6
    ),
    linewidth = 0.45,
    color = "black"
  ) +
  
  # WT and gsnor-ko headings
  annotate(
    "text",
    x = wt_x$xmid,
    y = max(row_lookup$y_pos) + 1.6,
    label = "WT",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = gs_x$xmid,
    y = max(row_lookup$y_pos) + 1.6,
    label = "italic(gsnor)*'-ko'",
    parse = TRUE,
    fontface = "bold",
    size = 5
  ) +
  
  # right-side header
  annotate(
    "text",
    x = right_header_x,
    y = max(row_lookup$y_pos) + 1.6,
    label = "Transcription factors | Families",
    fontface = "bold",
    size = 4.0,
    hjust = 0.5
  ) +
  
  scale_x_continuous(
    breaks = x_axis_df$x_pos,
    labels = x_axis_df$x_label,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Treatment and Generation",
    y = NULL
  ) +
  coord_cartesian(
    xlim = c(min(x_lookup$x_pos) - 0.7, right_family_x + 2.4),
    ylim = c(min(row_lookup$y_pos) - 0.5, max(row_lookup$y_pos) + 2.4),
    clip = "off"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = "black",
      size = 10
    ),
    axis.title.x = element_text(
      color = "black",
      size = 13,
      face = "bold",
      margin = margin(t = 22)
    ),
    
    legend.position = c(0.86, 0.93),
    legend.direction = "vertical",
    legend.key.height = unit(1.7, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.text = element_text(size = 10),
    
    plot.margin = margin(30, 120, 65, 20)
  )

print(p_tf_heatmap)

# ============================================================
# 10. Save
# ============================================================

ggsave(
  file.path(out_dir, "TF_heatmap_like_original.pdf"),
  p_tf_heatmap,
  width = 12,
  height = 10,
  useDingbats = FALSE
)

ggsave(
  file.path(out_dir, "TF_heatmap_like_original.png"),
  p_tf_heatmap,
  width = 12,
  height = 10,
  dpi = 300
)

cat("\nDone. TF heatmap saved in:\n")
cat(out_dir, "\n")
