# Revised copy: uses the eligible genes from each individual DESeq2 contrast
# as the background for its up/down Mercator enrichment tests. Re-run this
# script, then re-run "Go term_logic revived_Fig 4.R" to refresh Figure 4.
################### Standardized Mercator/GO analysis + attached-style plotting
################### clusterProfiler enrichment, old logic, new figure style

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(stringr)
  library(glue)
  library(purrr)
  library(readr)
  library(patchwork)
  library(scales)
})

set.seed(42)

## =========================
## CONFIG
## =========================

ROOT_DIR <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

SPLIT_DIR <- file.path(ROOT_DIR, "DESeq_split_up_down")
MERCATOR_DIR <- file.path(ROOT_DIR, "Mercator")

OUT_ROOT <- file.path(SPLIT_DIR, "Mercator_result_standardized_updated")
PLOT_DIR <- file.path(ROOT_DIR, "GO_Mercator_plots_attached_style")

dir.create(SPLIT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# DEG thresholds for split gene lists
PADJ_CUT <- 0.05
LFC_CUT  <- 0

# Enrichment thresholds: same as your old code
PV_CUT <- 0.05
Q_CUT  <- 0.20
MIN_GS <- 10
MAX_GS <- 500

# Plot filtering: same as your old plotting code
PLOT_Q_CUT <- 0.20
PLOT_P_CUT <- 0.05

levels_vec <- paste0("level", 1:8)

folders <- c(
  "WT_gen1", "WT_gen3", "WT_gen5",
  "gsnor_gen1", "gsnor_gen3", "gsnor_gen5"
)

treatments <- c("Drought", "CO2", "Ozone", "Heat", "Combination")

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

## =========================
## Helpers
## =========================

normalize_agi <- function(x) {
  x <- trimws(as.character(x))
  x <- str_replace(x, regex("^transcript:", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^gene:", ignore_case = TRUE), "")
  x <- sub("\\..*$", "", x)
  x <- toupper(x)
  str_extract(x, "AT([1-5]|C|M)G\\d{5}")
}

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
  dat$gene <- rownames(dat)
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

read_mercator_file <- function(fp) {
  
  ln <- readLines(fp, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  
  pos <- regexpr(",", ln, fixed = TRUE)
  keep <- pos > 0
  
  if (!all(keep)) {
    message(basename(fp), ": skipping ", sum(!keep), " line(s) without a comma.")
  }
  
  ln <- ln[keep]
  pos <- pos[keep]
  
  Term <- substr(ln, 1, pos - 1)
  Gene <- substr(ln, pos + 1, nchar(ln))
  
  strip <- function(z) {
    gsub('^["\\\']|["\\\']$', "", trimws(z))
  }
  
  tibble(
    Term = strip(Term),
    Gene = strip(Gene)
  )
}

prep_TERM2GENE <- function(df, minGSSize = MIN_GS, maxGSSize = MAX_GS) {
  
  stopifnot(all(c("ID", "gene") %in% colnames(df)))
  
  df <- df %>%
    mutate(
      ID = trimws(as.character(ID)),
      gene = trimws(as.character(gene))
    ) %>%
    filter(
      ID != "",
      gene != "",
      !is.na(ID),
      !is.na(gene)
    ) %>%
    distinct(ID, gene)
  
  sz <- df %>%
    count(ID, name = "set_size")
  
  keep_ids <- sz %>%
    filter(
      set_size >= minGSSize,
      set_size <= maxGSSize
    ) %>%
    pull(ID)
  
  df <- df %>%
    filter(ID %in% keep_ids)
  
  data.frame(
    ID = df$ID,
    gene = df$gene,
    stringsAsFactors = FALSE
  )
}

clean_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

infer_set_full_from_filename <- function(fp) {
  
  b <- tools::file_path_sans_ext(basename(fp))
  s <- clean_key(b)
  
  direction <- if (str_detect(s, "(^|_)up($|_)")) {
    "up"
  } else if (str_detect(s, "(^|_)down($|_)")) {
    "down"
  } else {
    NA_character_
  }
  
  geno <- if (str_detect(s, "(^|_)wt($|_)")) {
    "wt"
  } else if (str_detect(s, "gsnor")) {
    "gsnor"
  } else {
    NA_character_
  }
  
  gen <- str_match(s, "gen_?(1|3|5)")[, 2]
  
  if (is.na(gen)) {
    gen <- str_match(s, "(^|_)g(1|3|5)($|_)")[, 3]
  }
  
  trt <- str_match(s, "(drought|co2|ozone|heat|combination)")[, 2]
  
  if (any(is.na(c(geno, gen, trt, direction)))) {
    return(NA_character_)
  }
  
  paste0(geno, "_gen", gen, "_", trt, "_", direction)
}

parse_set_full_vec <- function(set_full_vec) {
  
  s <- clean_key(set_full_vec)
  
  mA <- str_match(
    s,
    "^(wt|gsnor)_gen(1|3|5)_(drought|co2|ozone|heat|combination)_(up|down)$"
  )
  
  mB <- str_match(
    s,
    "^(wt|gsnor)_(drought|co2|ozone|heat|combination)_gen(1|3|5)_(up|down)$"
  )
  
  useA <- !is.na(mA[, 2])
  
  geno_raw <- ifelse(useA, mA[, 2], mB[, 2])
  gen_raw  <- ifelse(useA, mA[, 3], mB[, 4])
  trt_raw  <- ifelse(useA, mA[, 4], mB[, 3])
  dir_raw  <- ifelse(useA, mA[, 5], mB[, 5])
  
  genotype <- case_when(
    is.na(geno_raw) ~ NA_character_,
    geno_raw == "wt" ~ "WT",
    TRUE ~ "gsnor"
  )
  
  generation <- ifelse(
    is.na(gen_raw),
    NA_character_,
    paste0("Gen", gen_raw)
  )
  
  treatment <- case_when(
    is.na(trt_raw) ~ NA_character_,
    trt_raw == "co2" ~ "CO2",
    trt_raw == "combination" ~ "Combination",
    TRUE ~ str_to_title(trt_raw)
  )
  
  sets_key <- case_when(
    is.na(trt_raw) ~ NA_character_,
    trt_raw == "combination" ~ "combi",
    TRUE ~ trt_raw
  )
  
  Sets <- ifelse(
    is.na(sets_key) | is.na(dir_raw),
    NA_character_,
    paste0(sets_key, "_", dir_raw)
  )
  
  sample_id <- ifelse(
    is.na(genotype) | is.na(treatment) | is.na(generation),
    NA_character_,
    paste0(genotype, "_", treatment, "_", generation)
  )
  
  tibble(
    set_full   = set_full_vec,
    genotype   = genotype,
    generation = generation,
    treatment  = treatment,
    Sets       = Sets,
    sample_id  = sample_id,
    direction  = dir_raw
  )
}

read_split_file <- function(fp, use_filename = TRUE) {
  
  df <- tryCatch(
    readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(df)) {
    df <- tryCatch(
      readr::read_csv(fp, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
  }
  
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  
  gene_candidates <- c(
    "gene", "Gene", "gene_id", "Gene_id", "GENE_ID",
    "AGI", "agi", "TAIR", "tair", "id", "ID"
  )
  
  gene_col <- gene_candidates[gene_candidates %in% names(df)][1]
  
  if (is.na(gene_col)) {
    gene_col <- names(df)[1]
  }
  
  set_full_one <- NA_character_
  
  if (use_filename) {
    set_full_one <- infer_set_full_from_filename(fp)
  }
  
  if (is.na(set_full_one)) {
    
    sets_col <- if ("Sets" %in% names(df)) {
      "Sets"
    } else if ("sets" %in% names(df)) {
      "sets"
    } else {
      NA_character_
    }
    
    if (!is.na(sets_col)) {
      vals <- df[[sets_col]] %>%
        as.character() %>%
        na.omit()
      
      vals <- vals[nzchar(trimws(vals))]
      
      if (length(vals)) {
        set_full_one <- vals[which.max(table(vals))]
      }
    }
    
    if (is.na(set_full_one)) {
      set_full_one <- tools::file_path_sans_ext(basename(fp))
    }
  }
  
  tibble(
    gene_id  = normalize_agi(df[[gene_col]]),
    set_full = rep(set_full_one, nrow(df))
  ) %>%
    filter(
      !is.na(gene_id),
      gene_id != "",
      !is.na(set_full),
      set_full != ""
    )
}

## =========================
## STEP 1: Create split up/down files from full DESeq2 results
## =========================

cat("\nCreating split up/down files in:\n")
cat(SPLIT_DIR, "\n")

split_summary <- list()
eligible_background <- list()

for (folder in folders) {
  
  folder_path <- file.path(ROOT_DIR, folder)
  
  if (!dir.exists(folder_path)) {
    warning("Folder not found: ", folder_path)
    next
  }
  
  genotype_name <- ifelse(str_detect(folder, "^WT"), "WT", "gsnor")
  
  generation_num <- str_extract(folder, "gen[0-9]+") %>%
    str_remove("gen")
  
  generation_name <- paste0("Gen", generation_num)
  
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
    
    treatment_name <- get_treatment_exact(file)
    
    if (is.na(treatment_name)) {
      next
    }
    
    dat <- read_deseq_file(file)
    
    if (!all(c("log2FoldChange", "padj") %in% colnames(dat))) {
      warning("Missing log2FoldChange/padj in: ", file)
      next
    }
    
    dat2 <- dat %>%
      mutate(
        gene = normalize_agi(gene),
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj)
      ) %>%
      filter(
        !is.na(gene),
        gene != "",
        !is.na(log2FoldChange),
        !is.na(padj)
      )

    # Contrast-specific background: genes eligible for a DEG call in this
    # genotype, generation and treatment, including non-DEGs. Keep one
    # background for both up- and down-regulated gene lists.
    contrast_id <- paste(genotype_name, treatment_name, generation_name, sep = "_")
    if (contrast_id %in% names(eligible_background)) {
      stop("Duplicate DESeq2 result for contrast: ", contrast_id)
    }
    eligible_background[[contrast_id]] <- unique(dat2$gene)
    
    up_genes <- dat2 %>%
      filter(
        padj < PADJ_CUT,
        log2FoldChange > LFC_CUT
      ) %>%
      distinct(gene) %>%
      arrange(gene)
    
    down_genes <- dat2 %>%
      filter(
        padj < PADJ_CUT,
        log2FoldChange < -LFC_CUT
      ) %>%
      distinct(gene) %>%
      arrange(gene)
    
    treatment_key <- tolower(treatment_name)
    
    up_file <- file.path(
      SPLIT_DIR,
      paste0(genotype_name, "_", generation_name, "_", treatment_key, "_up.txt")
    )
    
    down_file <- file.path(
      SPLIT_DIR,
      paste0(genotype_name, "_", generation_name, "_", treatment_key, "_down.txt")
    )
    
    write_tsv(up_genes, up_file)
    write_tsv(down_genes, down_file)
    
    split_summary[[length(split_summary) + 1]] <- tibble(
      genotype = genotype_name,
      generation = generation_name,
      treatment = treatment_name,
      up_n = nrow(up_genes),
      down_n = nrow(down_genes),
      up_file = basename(up_file),
      down_file = basename(down_file)
    )
  }
}

split_summary_df <- bind_rows(split_summary)

write_csv(
  split_summary_df,
  file.path(SPLIT_DIR, "split_DEG_file_summary.csv")
)

cat("\nSplit file summary:\n")
print(split_summary_df)

## =========================
## STEP 2: Build TERM2GENE
## =========================

merc_files <- list.files(
  MERCATOR_DIR,
  pattern = "_Arabidopsis_mercator\\.csv$",
  full.names = TRUE
)

if (!length(merc_files)) {
  stop("No *_Arabidopsis_mercator.csv in: ", MERCATOR_DIR)
}

term2gene <- merc_files %>%
  set_names(basename(.)) %>%
  imap_dfr(function(fp, nm) {
    read_mercator_file(fp) %>%
      mutate(Level = str_extract(nm, "level[1-8]"))
  }) %>%
  filter(Level %in% levels_vec) %>%
  mutate(
    Term = trimws(Term),
    Gene = normalize_agi(Gene)
  ) %>%
  filter(
    !is.na(Term),
    Term != "",
    !is.na(Gene),
    Gene != ""
  )

T2G_by_level <- map(levels_vec, ~{
  df <- term2gene %>%
    filter(Level == .x) %>%
    select(Term, Gene) %>%
    distinct()
  
  data.frame(
    ID = df$Term,
    gene = df$Gene,
    stringsAsFactors = FALSE
  )
}) %>%
  set_names(levels_vec)

message("\nTERM2GENE built for levels 1-8.")
print(table(term2gene$Level))

## =========================
## STEP 3: Load split files
## =========================

split_files <- list.files(
  SPLIT_DIR,
  pattern = "(_up|_down)\\.txt$",
  full.names = TRUE
)

if (!length(split_files)) {
  stop("No *_up.txt/_down.txt files found in: ", SPLIT_DIR)
}

mod2gene <- map_dfr(split_files, ~read_split_file(.x, use_filename = TRUE))

if (!nrow(mod2gene)) {
  stop("No genes loaded from split files.")
}

sets_meta <- parse_set_full_vec(unique(mod2gene$set_full))

bad <- sets_meta %>%
  filter(
    is.na(sample_id) |
      is.na(Sets) |
      is.na(treatment)
  )

if (nrow(bad)) {
  message("\nWARNING: some set labels could not be parsed. These will be dropped.")
  print(bad$set_full)
}

mod2gene2 <- mod2gene %>%
  left_join(sets_meta, by = "set_full") %>%
  filter(
    !is.na(sample_id),
    !is.na(Sets),
    !is.na(treatment),
    !is.na(genotype),
    !is.na(generation)
  )

if (!nrow(mod2gene2)) {
  stop("After parsing, no usable rows remain.")
}

message(
  "\nLoaded split files: ", length(split_files),
  " | Unique sets: ", n_distinct(mod2gene2$set_full),
  " | Genes total: ", nrow(mod2gene2)
)

message("\nParsed set summary:")
print(
  mod2gene2 %>%
    distinct(set_full, sample_id, genotype, generation, treatment, Sets) %>%
    arrange(sample_id) %>%
    head(30)
)

write_csv(
  mod2gene2,
  file.path(SPLIT_DIR, "mod2gene_loaded_for_Mercator.csv")
)

## =========================
## STEP 4: Run clusterProfiler enrichment
## =========================

all_results_list <- list()
background_audit <- list()

sample_dirs <- mod2gene2 %>%
  distinct(sample_id, genotype, generation, treatment) %>%
  mutate(out_dir = file.path(OUT_ROOT, sample_id))

walk(sample_dirs$out_dir, ~dir.create(.x, showWarnings = FALSE, recursive = TRUE))

for (sf in unique(mod2gene2$set_full)) {
  
  meta <- sets_meta %>%
    filter(set_full == sf)
  
  if (!nrow(meta)) {
    next
  }
  
  sample_id <- meta$sample_id[1]
  out_dir   <- file.path(OUT_ROOT, sample_id)
  set_label <- meta$Sets[1]
  
  genelist <- mod2gene2 %>%
    filter(set_full == sf) %>%
    pull(gene_id) %>%
    unique()
  
  genelist <- genelist[!is.na(genelist)]
  
  if (!length(genelist)) {
    next
  }

  universe_raw <- eligible_background[[sample_id]]
  if (is.null(universe_raw) || !length(universe_raw)) {
    stop("No contrast-specific eligible-gene background for: ", sample_id)
  }
  
  message("\n===== ", sample_id, " | ", set_label, " | genes=", length(genelist), " =====")
  
  for (lv in levels_vec) {
    
    t2g_raw <- T2G_by_level[[lv]]
    
    if (is.null(t2g_raw) || !nrow(t2g_raw)) {
      next
    }
    
    if (sum(genelist %in% t2g_raw$gene) == 0) {
      next
    }
    
    # Restrict the Mercator annotation to genes eligible in this contrast
    # BEFORE applying the 10-500 annotated-gene term-size filter.
    t2g <- prep_TERM2GENE(
      t2g_raw[t2g_raw$gene %in% universe_raw, , drop = FALSE],
      MIN_GS, MAX_GS
    )
    
    if (!nrow(t2g)) {
      next
    }
    
    universe_level <- intersect(universe_raw, unique(t2g$gene))
    genelist_level <- intersect(genelist, universe_level)

    background_audit[[length(background_audit) + 1L]] <- tibble(
      sample_id = sample_id,
      set_label = set_label,
      level = lv,
      eligible_genes_contrast = length(universe_raw),
      genes_in_level_background = length(universe_level),
      input_DEGs = length(genelist),
      DEGs_in_level_background = length(genelist_level),
      tested_terms = n_distinct(t2g$ID)
    )

    if (!length(genelist_level)) {
      next
    }
    
    enrich <- tryCatch(
      enricher(
        gene = genelist_level,
        universe = universe_level,
        TERM2GENE = t2g,
        pvalueCutoff = PV_CUT,
        qvalueCutoff = Q_CUT,
        pAdjustMethod = "BH",
        minGSSize = MIN_GS,
        maxGSSize = MAX_GS
      ),
      error = function(e) {
        message("  enricher error: ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(enrich)) {
      next
    }
    
    res <- as.data.frame(enrich)
    
    if (!nrow(res)) {
      next
    }
    
    csv_fp <- file.path(
      out_dir,
      paste0(lv, "_", set_label, "_Mercator_clusterprof.csv")
    )
    
    write.csv(res, csv_fp, row.names = FALSE)
    
    res2 <- res %>%
      mutate(
        sample_id  = meta$sample_id[1],
        genotype   = meta$genotype[1],
        generation = meta$generation[1],
        treatment  = meta$treatment[1],
        Sets       = set_label,
        Level      = lv
      )
    
    all_results_list[[length(all_results_list) + 1]] <- res2
  }
}

## =========================
## STEP 5: Combine results
## =========================

if (length(background_audit)) {
  write_csv(
    bind_rows(background_audit),
    file.path(OUT_ROOT, "Mercator_contrast_background_audit.csv")
  )
}

if (length(all_results_list)) {
  
  all_results <- bind_rows(all_results_list)
  
  if ("GeneRatio" %in% names(all_results)) {
    parts <- str_split(
      ifelse(is.na(all_results$GeneRatio), "", all_results$GeneRatio),
      "/",
      simplify = TRUE
    )
    
    all_results$GeneRatio_num <- suppressWarnings(
      as.numeric(parts[, 1]) / as.numeric(parts[, 2])
    )
  } else {
    all_results$GeneRatio_num <- NA_real_
  }
  
  comb_fp <- file.path(
    OUT_ROOT,
    "ALL_samples_Mercator_clusterprof_combined.csv"
  )
  
  write_csv(all_results, comb_fp)
  
  message("\nCombined CSV written: ", comb_fp)
  
  message("\nCounts by genotype/treatment/Sets/sample_id:")
  print(
    all_results %>%
      count(genotype, treatment, Sets, sample_id) %>%
      arrange(genotype, treatment, Sets, sample_id)
  )
  
} else {
  
  all_results <- tibble()
  message("\nNo results collected.")
}

if (!nrow(all_results)) {
  stop("No enrichment results found. Check split lists and Mercator annotation overlap.")
}

## =========================
## STEP 6: Attached-style plotting functions
## =========================

wrap_terms <- function(x, width = 42) {
  stringr::str_wrap(x, width = width)
}

safe_neglog10 <- function(x) {
  -log10(pmax(as.numeric(x), 1e-300))
}

make_size_breaks <- function(max_count) {
  if (is.na(max_count) || max_count <= 10) return(c(5, 10))
  if (max_count <= 30) return(c(10, 20, 30))
  if (max_count <= 60) return(c(10, 20, 40, 60))
  if (max_count <= 100) return(c(20, 40, 60, 80, 100))
  c(50, 100, 150, 200)
}

display_genotype <- function(x) {
  ifelse(x == "gsnor", "gsnor-ko", x)
}

direction_title <- function(x) {
  ifelse(x == "up", "upregulated", "downregulated")
}

make_attached_style_panel <- function(df,
                                      genotype_keep,
                                      treatment_keep,
                                      direction_keep,
                                      samples_keep,
                                      panel_letter = NULL,
                                      q_cut = PLOT_Q_CUT,
                                      p_cut = PLOT_P_CUT,
                                      term_wrap = 42,
                                      save_individual = TRUE) {
  
  trt_key <- tolower(treatment_keep)
  trt_key <- ifelse(trt_key == "combination", "combi", trt_key)
  set_keep <- paste0(trt_key, "_", direction_keep)
  
  dd <- df %>%
    filter(
      Sets == set_keep,
      treatment == treatment_keep,
      genotype == genotype_keep,
      sample_id %in% samples_keep,
      !is.na(p.adjust),
      !is.na(pvalue),
      p.adjust <= q_cut,
      pvalue <= p_cut
    )
  
  title_text <- paste0(
    ifelse(is.null(panel_letter), "", paste0("(", panel_letter, ") ")),
    display_genotype(genotype_keep),
    " ",
    treatment_keep,
    " (",
    direction_title(direction_keep),
    ")"
  )
  
  if (!nrow(dd)) {
    
    p_empty <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = paste0(title_text, "\nNo enriched terms"),
        size = 5
      ) +
      theme_void() +
      labs(title = title_text) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0)
      )
    
    return(list(plot = p_empty, n_terms = 0))
  }
  
  dd_sample <- dd %>%
    mutate(
      term_key = paste(Level, ID, sep = "__"),
      Description = as.character(Description),
      Level = as.character(Level),
      neg_log10_padj = safe_neglog10(p.adjust)
    ) %>%
    group_by(term_key, ID, Description, Level, sample_id) %>%
    summarise(
      pvalue = min(pvalue, na.rm = TRUE),
      p.adjust = min(p.adjust, na.rm = TRUE),
      Count = max(Count, na.rm = TRUE),
      GeneRatio = dplyr::first(GeneRatio),
      geneID = paste(unique(geneID), collapse = "/"),
      neg_log10_padj = max(neg_log10_padj, na.rm = TRUE),
      .groups = "drop"
    )
  
  term_summary <- dd_sample %>%
    group_by(term_key, ID, Description, Level) %>%
    summarise(
      best_padj = min(p.adjust, na.rm = TRUE),
      max_count = max(Count, na.rm = TRUE),
      best_score = max(neg_log10_padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      level_num = as.integer(str_remove(Level, "level")),
      term_label = wrap_terms(Description, width = term_wrap)
    ) %>%
    arrange(
      level_num,
      desc(best_score),
      desc(max_count),
      Description
    )
  
  term_key_levels <- rev(term_summary$term_key)
  term_label_lookup <- term_summary$term_label
  names(term_label_lookup) <- term_summary$term_key
  
  x_df <- tibble(
    sample_id = factor(samples_keep, levels = samples_keep),
    x_pos = seq_along(samples_keep),
    x_label = samples_keep
  )
  
  plot_df <- expand_grid(
    term_key = term_summary$term_key,
    sample_id = factor(samples_keep, levels = samples_keep)
  ) %>%
    left_join(
      term_summary %>%
        select(term_key, ID, Description, Level, term_label),
      by = "term_key"
    ) %>%
    left_join(
      dd_sample %>%
        select(
          term_key,
          sample_id,
          pvalue,
          p.adjust,
          Count,
          GeneRatio,
          geneID,
          neg_log10_padj
        ) %>%
        mutate(sample_id = factor(sample_id, levels = samples_keep)),
      by = c("term_key", "sample_id")
    ) %>%
    left_join(x_df, by = "sample_id") %>%
    mutate(
      Count = as.numeric(Count),
      Count = ifelse(is.na(Count), 0, Count),
      
      # This makes tiny dark dots where a term is absent in that generation,
      # matching your attached plot style.
      Count_for_size = ifelse(Count > 0, Count, 1),
      neg_log10_padj = ifelse(is.na(p.adjust), 0, safe_neglog10(p.adjust)),
      term_key = factor(term_key, levels = term_key_levels)
    )
  
  strip_df <- term_summary %>%
    mutate(
      term_key = factor(term_key, levels = term_key_levels),
      x_pos = max(x_df$x_pos) + 0.45,
      Level = factor(Level, levels = levels_vec)
    )
  
  max_col <- max(plot_df$neg_log10_padj, na.rm = TRUE)
  if (!is.finite(max_col) || max_col == 0) max_col <- 1
  
  max_count <- max(plot_df$Count, na.rm = TRUE)
  size_breaks <- make_size_breaks(max_count)
  
  level_breaks_present <- levels_vec[levels_vec %in% unique(as.character(strip_df$Level))]
  
  p <- ggplot() +
    geom_point(
      data = plot_df,
      aes(
        x = x_pos,
        y = term_key,
        color = neg_log10_padj,
        size = Count_for_size
      ),
      alpha = 0.95
    ) +
    geom_tile(
      data = strip_df,
      aes(
        x = x_pos,
        y = term_key,
        fill = Level
      ),
      width = 0.18,
      height = 0.88
    ) +
    scale_y_discrete(
      labels = term_label_lookup
    ) +
    scale_x_continuous(
      breaks = x_df$x_pos,
      labels = x_df$x_label,
      limits = c(0.5, max(x_df$x_pos) + 0.85),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_color_gradientn(
      colors = c("#1f0088", "#8c1ab1", "#e66f61", "#ffe51f"),
      limits = c(0, max_col),
      name = expression(-log[10]("adj p"))
    ) +
    scale_size_continuous(
      range = c(0.5, 6.2),
      breaks = size_breaks,
      name = "Gene Count"
    ) +
    scale_fill_manual(
      values = level_cols,
      breaks = level_breaks_present,
      drop = FALSE,
      name = "Mercator level"
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2),
      fill = guide_legend(order = 3)
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = "GO / Mercator term"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0,
        margin = margin(b = 8)
      ),
      panel.grid.major.x = element_line(
        color = "grey88",
        linewidth = 0.35
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "grey35",
        fill = NA,
        linewidth = 0.6
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black",
        size = 8
      ),
      axis.text.y = element_text(
        color = "black",
        size = 8
      ),
      axis.title.y = element_text(
        color = "black",
        size = 10
      ),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 10, 8, 8)
    )
  
  if (save_individual) {
    
    n_terms <- nrow(term_summary)
    fig_height <- max(5, min(25, 2.8 + 0.18 * n_terms))
    
    file_stub <- paste0(
      genotype_keep,
      "_",
      treatment_keep,
      "_",
      direction_keep,
      "_attached_style"
    ) %>%
      str_replace_all("[^A-Za-z0-9_\\-]", "_")
    
    ggsave(
      file.path(PLOT_DIR, paste0(file_stub, ".pdf")),
      p,
      width = 5,
      height = fig_height,
      useDingbats = FALSE,
      limitsize = FALSE
    )
    
    ggsave(
      file.path(PLOT_DIR, paste0(file_stub, ".png")),
      p,
      width = 5,
      height = fig_height,
      dpi = 300,
      limitsize = FALSE
    )
  }
  
  list(plot = p, n_terms = nrow(term_summary))
}

## =========================
## STEP 7: Panel letters
## =========================

panel_letter_table <- tibble::tribble(
  ~genotype, ~treatment,     ~direction, ~letter,
  "WT",      "Drought",      "up",       "A",
  "WT",      "Drought",      "down",     "B",
  "WT",      "CO2",          "up",       "C",
  "WT",      "CO2",          "down",     "D",
  "WT",      "Ozone",        "up",       "E",
  "WT",      "Ozone",        "down",     "F",
  "WT",      "Heat",         "up",       "G",
  "WT",      "Heat",         "down",     "H",
  "WT",      "Combination",  "up",       "I",
  "WT",      "Combination",  "down",     "J",
  "gsnor",   "Drought",      "up",       "K",
  "gsnor",   "Drought",      "down",     "L",
  "gsnor",   "CO2",          "up",       "M",
  "gsnor",   "CO2",          "down",     "N",
  "gsnor",   "Ozone",        "up",       "O",
  "gsnor",   "Ozone",        "down",     "P",
  "gsnor",   "Heat",         "up",       "Q",
  "gsnor",   "Heat",         "down",     "R",
  "gsnor",   "Combination",  "up",       "S",
  "gsnor",   "Combination",  "down",     "T"
)

get_panel_letter <- function(genotype_keep, treatment_keep, direction_keep) {
  panel_letter_table %>%
    filter(
      genotype == genotype_keep,
      treatment == treatment_keep,
      direction == direction_keep
    ) %>%
    pull(letter) %>%
    dplyr::first()
}

## =========================
## STEP 8: Make attached-style paired plots
## =========================

all_attached_panels <- list()

for (geno in c("WT", "gsnor")) {
  
  for (trt in treatments) {
    
    samples_keep <- paste0(
      geno,
      "_",
      trt,
      "_",
      c("Gen1", "Gen3", "Gen5")
    )
    
    up_letter <- get_panel_letter(geno, trt, "up")
    down_letter <- get_panel_letter(geno, trt, "down")
    
    cat("\nMaking attached-style plot: ", geno, " ", trt, " up/down\n")
    
    up_obj <- make_attached_style_panel(
      df = all_results,
      genotype_keep = geno,
      treatment_keep = trt,
      direction_keep = "up",
      samples_keep = samples_keep,
      panel_letter = up_letter,
      save_individual = TRUE
    )
    
    down_obj <- make_attached_style_panel(
      df = all_results,
      genotype_keep = geno,
      treatment_keep = trt,
      direction_keep = "down",
      samples_keep = samples_keep,
      panel_letter = down_letter,
      save_individual = TRUE
    )
    
    p_up <- up_obj$plot
    p_down <- down_obj$plot
    
    paired_plot <- p_up + p_down +
      plot_layout(widths = c(1, 1), guides = "keep") +
      plot_annotation(
        title = "Supplementary Figure S16: Functional enrichment",
        theme = theme(
          plot.title = element_text(
            size = 19,
            face = "bold",
            hjust = 0,
            margin = margin(b = 10)
          )
        )
      )
    
    max_terms <- max(up_obj$n_terms, down_obj$n_terms, na.rm = TRUE)
    fig_height <- max(8, min(26, 3.5 + 0.18 * max_terms))
    
    file_stub <- paste0(
      geno,
      "_",
      trt,
      "_up_down_attached_style"
    ) %>%
      str_replace_all("[^A-Za-z0-9_\\-]", "_")
    
    ggsave(
      file.path(PLOT_DIR, paste0(file_stub, ".pdf")),
      paired_plot,
      width = 14,
      height = fig_height,
      useDingbats = FALSE,
      limitsize = FALSE
    )
    
    ggsave(
      file.path(PLOT_DIR, paste0(file_stub, ".png")),
      paired_plot,
      width = 14,
      height = fig_height,
      dpi = 300,
      limitsize = FALSE
    )
    
    all_attached_panels[[file_stub]] <- paired_plot
  }
}

cat("\nDone.\n")
cat("Split DEG files saved in:\n")
cat(SPLIT_DIR, "\n\n")
cat("clusterProfiler enrichment results saved in:\n")
cat(OUT_ROOT, "\n\n")
cat("Attached-style plots saved in:\n")
cat(PLOT_DIR, "\n")
