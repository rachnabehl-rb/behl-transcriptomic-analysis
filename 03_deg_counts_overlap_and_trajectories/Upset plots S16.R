rm(list = ls())

suppressPackageStartupMessages({
  library(UpSetR)
  library(dplyr)
})

# ============================================================
# Input settings
# ============================================================

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"

out_dir <- file.path(root_dir, "Simple_UpSet_all_treatments_FINAL_Yaxis_fixed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 0

# ============================================================
# Helper functions
# ============================================================

read_deseq_file <- function(file) {
  
  dat <- read.table(
    file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  dat <- as.data.frame(dat)
  dat$gene_id <- rownames(dat)
  dat
}

get_deg_genes <- function(file) {
  
  dat <- read_deseq_file(file)
  
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

# ============================================================
# Drought
# Visual y-axis order wanted:
# WT_gen5, gsnor_gen1, WT_gen3, gsnor_gen3, WT_gen1
# gsnor_gen5 is removed automatically if zero
# ============================================================

WT_gen1_Drought <- get_deg_genes(file.path(root_dir, "WT_gen1", "Drought.txt"))
gsnor_gen1_Drought <- get_deg_genes(file.path(root_dir, "gsnor_gen1", "Drought.txt"))
WT_gen3_Drought <- get_deg_genes(file.path(root_dir, "WT_gen3", "Drought.txt"))
gsnor_gen3_Drought <- get_deg_genes(file.path(root_dir, "gsnor_gen3", "Drought.txt"))
WT_gen5_Drought <- get_deg_genes(file.path(root_dir, "WT_gen5", "Drought.txt"))
gsnor_gen5_Drought <- get_deg_genes(file.path(root_dir, "gsnor_gen5", "Drought.txt"))

drought_sets <- list(
  WT_gen1 = WT_gen1_Drought,
  gsnor_gen3 = gsnor_gen3_Drought,
  WT_gen3 = WT_gen3_Drought,
  gsnor_gen1 = gsnor_gen1_Drought,
  WT_gen5 = WT_gen5_Drought,
  gsnor_gen5 = gsnor_gen5_Drought
)

drought_sets <- drought_sets[lengths(drought_sets) > 0]
drought_input <- UpSetR::fromList(drought_sets)

drought_set_order <- c(
  "WT_gen1",
  "gsnor_gen3",
  "WT_gen3",
  "gsnor_gen1",
  "WT_gen5",
  "gsnor_gen5"
)
drought_set_order <- drought_set_order[drought_set_order %in% names(drought_sets)]

pdf(file.path(out_dir, "UpSet_Drought_FINAL.pdf"), width = 10, height = 6, useDingbats = FALSE)
suppressWarnings(
  UpSetR::upset(
    drought_input,
    sets = drought_set_order,
    keep.order = TRUE,
    nsets = length(drought_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.8,
    line.size = 0.8,
    text.scale = c(1.3, 1.2, 1.1, 1.1, 1.2, 1.1)
  )
)
dev.off()

png(file.path(out_dir, "UpSet_Drought_FINAL.png"), width = 10, height = 6, units = "in", res = 300)
suppressWarnings(
  UpSetR::upset(
    drought_input,
    sets = drought_set_order,
    keep.order = TRUE,
    nsets = length(drought_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.8,
    line.size = 0.8,
    text.scale = c(1.3, 1.2, 1.1, 1.1, 1.2, 1.1)
  )
)
dev.off()

# ============================================================
# CO2
# Visual y-axis order wanted:
# WT_gen1, gsnor_gen1, gsnor_gen5, WT_gen3, WT_gen5, gsnor_gen3
# Therefore sets must be reverse:
# gsnor_gen3, WT_gen5, WT_gen3, gsnor_gen5, gsnor_gen1, WT_gen1
# ============================================================

WT_gen1_CO2 <- get_deg_genes(file.path(root_dir, "WT_gen1", "CO2.txt"))
gsnor_gen1_CO2 <- get_deg_genes(file.path(root_dir, "gsnor_gen1", "CO2.txt"))
WT_gen3_CO2 <- get_deg_genes(file.path(root_dir, "WT_gen3", "CO2.txt"))
gsnor_gen3_CO2 <- get_deg_genes(file.path(root_dir, "gsnor_gen3", "CO2.txt"))
WT_gen5_CO2 <- get_deg_genes(file.path(root_dir, "WT_gen5", "CO2.txt"))
gsnor_gen5_CO2 <- get_deg_genes(file.path(root_dir, "gsnor_gen5", "CO2.txt"))

co2_sets <- list(
  gsnor_gen3 = gsnor_gen3_CO2,
  WT_gen5 = WT_gen5_CO2,
  WT_gen3 = WT_gen3_CO2,
  gsnor_gen5 = gsnor_gen5_CO2,
  gsnor_gen1 = gsnor_gen1_CO2,
  WT_gen1 = WT_gen1_CO2
)

co2_sets <- co2_sets[lengths(co2_sets) > 0]
co2_input <- UpSetR::fromList(co2_sets)

co2_set_order <- c(
  "gsnor_gen3",
  "WT_gen5",
  "WT_gen3",
  "gsnor_gen5",
  "gsnor_gen1",
  "WT_gen1"
)
co2_set_order <- co2_set_order[co2_set_order %in% names(co2_sets)]

pdf(file.path(out_dir, "UpSet_CO2_FINAL.pdf"), width = 10, height = 6, useDingbats = FALSE)
suppressWarnings(
  UpSetR::upset(
    co2_input,
    sets = co2_set_order,
    keep.order = TRUE,
    nsets = length(co2_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.8,
    line.size = 0.8,
    text.scale = c(1.3, 1.2, 1.1, 1.1, 1.2, 1.1)
  )
)
dev.off()

png(file.path(out_dir, "UpSet_CO2_FINAL.png"), width = 10, height = 6, units = "in", res = 300)
suppressWarnings(
  UpSetR::upset(
    co2_input,
    sets = co2_set_order,
    keep.order = TRUE,
    nsets = length(co2_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.8,
    line.size = 0.8,
    text.scale = c(1.3, 1.2, 1.1, 1.1, 1.2, 1.1)
  )
)
dev.off()

# ============================================================
# Ozone
# Visual y-axis order wanted:
# WT_gen1, WT_gen5, WT_gen3, gsnor_gen5, gsnor_gen1, gsnor_gen3
# Reverse for UpSetR:
# gsnor_gen3, gsnor_gen1, gsnor_gen5, WT_gen3, WT_gen5, WT_gen1
# ============================================================

WT_gen1_Ozone <- get_deg_genes(file.path(root_dir, "WT_gen1", "Ozone.txt"))
gsnor_gen1_Ozone <- get_deg_genes(file.path(root_dir, "gsnor_gen1", "Ozone.txt"))
WT_gen3_Ozone <- get_deg_genes(file.path(root_dir, "WT_gen3", "Ozone.txt"))
gsnor_gen3_Ozone <- get_deg_genes(file.path(root_dir, "gsnor_gen3", "Ozone.txt"))
WT_gen5_Ozone <- get_deg_genes(file.path(root_dir, "WT_gen5", "Ozone.txt"))
gsnor_gen5_Ozone <- get_deg_genes(file.path(root_dir, "gsnor_gen5", "Ozone.txt"))

ozone_sets <- list(
  gsnor_gen3 = gsnor_gen3_Ozone,
  gsnor_gen1 = gsnor_gen1_Ozone,
  gsnor_gen5 = gsnor_gen5_Ozone,
  WT_gen3 = WT_gen3_Ozone,
  WT_gen5 = WT_gen5_Ozone,
  WT_gen1 = WT_gen1_Ozone
)

ozone_sets <- ozone_sets[lengths(ozone_sets) > 0]
ozone_input <- UpSetR::fromList(ozone_sets)

ozone_set_order <- c(
  "gsnor_gen3",
  "gsnor_gen1",
  "gsnor_gen5",
  "WT_gen3",
  "WT_gen5",
  "WT_gen1"
)
ozone_set_order <- ozone_set_order[ozone_set_order %in% names(ozone_sets)]

pdf(file.path(out_dir, "UpSet_Ozone_FINAL.pdf"), width = 10, height = 6, useDingbats = FALSE)
suppressWarnings(
  UpSetR::upset(
    ozone_input,
    sets = ozone_set_order,
    keep.order = TRUE,
    nsets = length(ozone_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.8,
    line.size = 0.8,
    text.scale = c(1.3, 1.2, 1.1, 1.1, 1.2, 1.1)
  )
)
dev.off()

png(file.path(out_dir, "UpSet_Ozone_FINAL.png"), width = 10, height = 6, units = "in", res = 300)
suppressWarnings(
  UpSetR::upset(
    ozone_input,
    sets = ozone_set_order,
    keep.order = TRUE,
    nsets = length(ozone_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.8,
    line.size = 0.8,
    text.scale = c(1.3, 1.2, 1.1, 1.1, 1.2, 1.1)
  )
)
dev.off()

# ============================================================
# Heat
# Visual y-axis order wanted:
# WT_gen5, gsnor_gen3, gsnor_gen5, WT_gen1, WT_gen3, gsnor_gen1
# Reverse for UpSetR:
# gsnor_gen1, WT_gen3, WT_gen1, gsnor_gen5, gsnor_gen3, WT_gen5
# ============================================================

WT_gen1_Heat <- get_deg_genes(file.path(root_dir, "WT_gen1", "Heat.txt"))
gsnor_gen1_Heat <- get_deg_genes(file.path(root_dir, "gsnor_gen1", "Heat.txt"))
WT_gen3_Heat <- get_deg_genes(file.path(root_dir, "WT_gen3", "Heat.txt"))
gsnor_gen3_Heat <- get_deg_genes(file.path(root_dir, "gsnor_gen3", "Heat.txt"))
WT_gen5_Heat <- get_deg_genes(file.path(root_dir, "WT_gen5", "Heat.txt"))
gsnor_gen5_Heat <- get_deg_genes(file.path(root_dir, "gsnor_gen5", "Heat.txt"))

heat_sets <- list(
  gsnor_gen1 = gsnor_gen1_Heat,
  WT_gen3 = WT_gen3_Heat,
  WT_gen1 = WT_gen1_Heat,
  gsnor_gen5 = gsnor_gen5_Heat,
  gsnor_gen3 = gsnor_gen3_Heat,
  WT_gen5 = WT_gen5_Heat
)

heat_sets <- heat_sets[lengths(heat_sets) > 0]
heat_input <- UpSetR::fromList(heat_sets)

heat_set_order <- c(
  "gsnor_gen1",
  "WT_gen3",
  "WT_gen1",
  "gsnor_gen5",
  "gsnor_gen3",
  "WT_gen5"
)
heat_set_order <- heat_set_order[heat_set_order %in% names(heat_sets)]

pdf(file.path(out_dir, "UpSet_Heat_FINAL.pdf"), width = 12, height = 6, useDingbats = FALSE)
suppressWarnings(
  UpSetR::upset(
    heat_input,
    sets = heat_set_order,
    keep.order = TRUE,
    nsets = length(heat_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.5,
    line.size = 0.8,
    text.scale = c(1.2, 1.1, 1.0, 1.0, 1.1, 1.0)
  )
)
dev.off()

png(file.path(out_dir, "UpSet_Heat_FINAL.png"), width = 12, height = 6, units = "in", res = 300)
suppressWarnings(
  UpSetR::upset(
    heat_input,
    sets = heat_set_order,
    keep.order = TRUE,
    nsets = length(heat_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.5,
    line.size = 0.8,
    text.scale = c(1.2, 1.1, 1.0, 1.0, 1.1, 1.0)
  )
)
dev.off()

# ============================================================
# Combination
# Visual y-axis order wanted:
# WT_gen5, WT_gen3, gsnor_gen3, gsnor_gen5, WT_gen1, gsnor_gen1
# Reverse for UpSetR:
# gsnor_gen1, WT_gen1, gsnor_gen5, gsnor_gen3, WT_gen3, WT_gen5
# ============================================================

WT_gen1_Combination <- get_deg_genes(file.path(root_dir, "WT_gen1", "Combination.txt"))
gsnor_gen1_Combination <- get_deg_genes(file.path(root_dir, "gsnor_gen1", "Combination.txt"))
WT_gen3_Combination <- get_deg_genes(file.path(root_dir, "WT_gen3", "Combination.txt"))
gsnor_gen3_Combination <- get_deg_genes(file.path(root_dir, "gsnor_gen3", "Combination.txt"))
WT_gen5_Combination <- get_deg_genes(file.path(root_dir, "WT_gen5", "Combination.txt"))
gsnor_gen5_Combination <- get_deg_genes(file.path(root_dir, "gsnor_gen5", "Combination.txt"))

combination_sets <- list(
  gsnor_gen1 = gsnor_gen1_Combination,
  WT_gen1 = WT_gen1_Combination,
  gsnor_gen5 = gsnor_gen5_Combination,
  gsnor_gen3 = gsnor_gen3_Combination,
  WT_gen3 = WT_gen3_Combination,
  WT_gen5 = WT_gen5_Combination
)

combination_sets <- combination_sets[lengths(combination_sets) > 0]
combination_input <- UpSetR::fromList(combination_sets)

combination_set_order <- c(
  "gsnor_gen1",
  "WT_gen1",
  "gsnor_gen5",
  "gsnor_gen3",
  "WT_gen3",
  "WT_gen5"
)
combination_set_order <- combination_set_order[combination_set_order %in% names(combination_sets)]

pdf(file.path(out_dir, "UpSet_Combination_FINAL.pdf"), width = 12, height = 6, useDingbats = FALSE)
suppressWarnings(
  UpSetR::upset(
    combination_input,
    sets = combination_set_order,
    keep.order = TRUE,
    nsets = length(combination_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.5,
    line.size = 0.8,
    text.scale = c(1.2, 1.1, 1.0, 1.0, 1.1, 1.0)
  )
)
dev.off()

png(file.path(out_dir, "UpSet_Combination_FINAL.png"), width = 12, height = 6, units = "in", res = 300)
suppressWarnings(
  UpSetR::upset(
    combination_input,
    sets = combination_set_order,
    keep.order = TRUE,
    nsets = length(combination_set_order),
    nintersects = NA,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    main.bar.color = "grey25",
    sets.bar.color = "grey25",
    matrix.color = "grey25",
    point.size = 2.5,
    line.size = 0.8,
    text.scale = c(1.2, 1.1, 1.0, 1.0, 1.1, 1.0)
  )
)
dev.off()

cat("\nSaved all final UpSet plots in:\n")
cat(out_dir, "\n")

