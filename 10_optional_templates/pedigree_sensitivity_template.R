# Pedigree sensitivity checks for phenotyping and RNA-seq.
# Fill the three input paths below and verify the pedigree/parent IDs against
# propagation records. Sample-name prefixes alone are insufficient evidence.
# This script DOES NOT modify your original analyses or source data.

PHENOTYPE_FILE <- "phenotype_week4.csv"
METADATA_FILE  <- "rna_sample_metadata_verified.csv"
COUNT_FILE     <- "Count_matrix.txt"
OUTPUT_DIR     <- "pedigree_sensitivity_results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Required phenotype.csv columns:
# genotype, treatment, generation, pedigree_id, plant_id, rosette_area
# One row per measured plant (e.g. week 4 only); pedigree_id is the VERIFIED
# ancestor/parent family for generations 2-5. G1 pedigree_id may be blank.
# If you have multiple repeated measurements per plant, choose ONE occasion
# before running this script; repeated measurements are not new plants.

# Required RNA metadata.csv columns:
# sample_id, genotype, treatment, generation, pedigree_id
# One row per sequencing library. Verify G3/G5 parent IDs from propagation
# records. Sample IDs must match columns of the raw integer count matrix.
# Optional: RNA_CONTROL_SUBSET_ONLY <- TRUE to avoid all-contrast checks.

## 1. Phenotype: plant-level means versus equal-weight pedigree means.
if (file.exists(PHENOTYPE_FILE)) {
  ph <- read.csv(PHENOTYPE_FILE, stringsAsFactors = FALSE,
                 check.names = FALSE)
  required <- c("genotype", "treatment", "generation", "pedigree_id",
                "plant_id", "rosette_area")
  stopifnot(all(required %in% names(ph)))
  ph$generation <- as.integer(ph$generation)
  ph$rosette_area <- as.numeric(ph$rosette_area)
  stopifnot(all(!is.na(ph$generation)), all(is.finite(ph$rosette_area)))
  ph <- subset(ph, generation >= 2L)
  stopifnot(nrow(ph) > 0L, all(!is.na(ph$pedigree_id) & nzchar(ph$pedigree_id)))
  plant_key <- paste(ph$genotype, ph$treatment, ph$generation, ph$plant_id,
                     sep = "|")
  if (anyDuplicated(plant_key)) {
    stop("Phenotype file has repeated plant IDs within a group. Select one measurement per plant.")
  }
  library(dplyr)
  group_vars <- c("genotype", "treatment", "generation")
  plant_summary <- ph %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(plant_mean = mean(rosette_area), n_plants = n(), .groups = "drop")
  pedigree_means <- ph %>%
    group_by(across(all_of(c(group_vars, "pedigree_id")))) %>%
    summarise(pedigree_mean = mean(rosette_area),
              plants_per_pedigree = n(), .groups = "drop")
  pedigree_summary <- pedigree_means %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(equal_weight_pedigree_mean = mean(pedigree_mean),
              n_pedigrees = n(),
              min_plants_per_pedigree = min(plants_per_pedigree),
              max_plants_per_pedigree = max(plants_per_pedigree),
              .groups = "drop")
  comparison <- left_join(plant_summary, pedigree_summary, by = group_vars) %>%
    mutate(change_when_equal_weighted = equal_weight_pedigree_mean - plant_mean)
  write.csv(pedigree_means, file.path(OUTPUT_DIR, "phenotype_pedigree_means.csv"),
            row.names = FALSE)
  write.csv(comparison, file.path(OUTPUT_DIR, "phenotype_group_comparison.csv"),
            row.names = FALSE)

  # Compare each treatment with its control within genotype and generation;
  # this is descriptive, not a chamber-replicated treatment test.
  control <- comparison %>% filter(treatment == "Control") %>%
    select(genotype, generation,
           control_plant_mean = plant_mean,
           control_pedigree_mean = equal_weight_pedigree_mean)
  effects <- comparison %>% filter(treatment != "Control") %>%
    inner_join(control, by = c("genotype", "generation")) %>%
    mutate(plant_treatment_minus_control = plant_mean - control_plant_mean,
           pedigree_treatment_minus_control =
             equal_weight_pedigree_mean - control_pedigree_mean,
           same_direction = sign(plant_treatment_minus_control) ==
                            sign(pedigree_treatment_minus_control))
  write.csv(effects, file.path(OUTPUT_DIR, "phenotype_contrast_directions.csv"),
            row.names = FALSE)
  message("Phenotype summaries complete; inspect direction and magnitude of contrasts.")
} else {
  message("Phenotype file not found; phenotype analysis skipped: ", PHENOTYPE_FILE)
}

## 2. RNA-seq: keep one library per verified pedigree in G3/G5 controls.
# Each 5-library control should have one pedigree represented twice. Run BOTH
# possible subsets (remove either of those two libraries), not a selected best
# result. Treatment groups should already contain one library per pedigree.
if (file.exists(METADATA_FILE) && file.exists(COUNT_FILE)) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 first.")
  meta <- read.csv(METADATA_FILE, stringsAsFactors = FALSE,
                   check.names = FALSE)
  needed <- c("sample_id", "genotype", "treatment", "generation", "pedigree_id")
  stopifnot(all(needed %in% names(meta)))
  meta$generation <- as.integer(meta$generation)
  meta <- subset(meta, generation %in% c(3L, 5L))
  if (anyDuplicated(meta$sample_id)) stop("Duplicate RNA sample_id in metadata.")
  if (anyNA(meta$pedigree_id) || any(!nzchar(meta$pedigree_id))) {
    stop("Supply verified pedigree IDs for all G3/G5 RNA libraries.")
  }
  counts_all <- as.matrix(read.table(COUNT_FILE, header = TRUE,
                                      row.names = 1, check.names = FALSE,
                                      sep = "", quote = "", comment.char = ""))
  if (anyDuplicated(rownames(counts_all))) stop("Duplicate gene IDs in count matrix.")
  if (!all(meta$sample_id %in% colnames(counts_all))) {
    stop("RNA metadata sample IDs not found in raw count matrix: ",
         paste(setdiff(meta$sample_id, colnames(counts_all)), collapse = ", "))
  }
  if (anyNA(counts_all) || any(counts_all < 0) ||
      any(abs(counts_all - round(counts_all)) > 1e-8)) {
    stop("The input must be untransformed, nonnegative integer gene counts.")
  }
  storage.mode(counts_all) <- "integer"

  fit_pair <- function(pair_metadata) {
    pair_counts <- counts_all[, pair_metadata$sample_id, drop = FALSE]
    condition <- factor(ifelse(pair_metadata$treatment == "Control", "C", "S"),
                        levels = c("C", "S"))
    coldata <- data.frame(row.names = pair_metadata$sample_id, condition)
    dds <- DESeq2::DESeqDataSetFromMatrix(pair_counts, coldata, design = ~condition)
    # Matches the filters in the supplied G3/G5 pairwise DESeq2 scripts.
    dds <- dds[rowSums(DESeq2::counts(dds)) > 10, ]
    dds <- DESeq2::estimateSizeFactors(dds)
    keep <- rowSums(DESeq2::counts(dds, normalized = TRUE) >= 3) >= 2
    dds <- dds[keep, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    as.data.frame(DESeq2::results(dds, contrast = c("condition", "S", "C")))
  }

  output <- list()
  for (gt in unique(meta$genotype)) for (gn in sort(unique(meta$generation))) {
    strata <- subset(meta, genotype == gt & generation == gn)
    ctrl <- subset(strata, treatment == "Control")
    duplicated_families <- names(which(table(ctrl$pedigree_id) > 1L))
    if (nrow(ctrl) != 5L || length(unique(ctrl$pedigree_id)) != 4L ||
        length(duplicated_families) != 1L ||
        sum(ctrl$pedigree_id == duplicated_families) != 2L) {
      stop("Expected 5 controls from 4 verified pedigrees (one pair): ", gt, " G", gn)
    }
    twin_ids <- ctrl$sample_id[ctrl$pedigree_id == duplicated_families]
    for (tr in setdiff(unique(strata$treatment), "Control")) {
      treated <- subset(strata, treatment == tr)
      if (nrow(treated) != 4L || anyDuplicated(treated$pedigree_id)) {
        stop("Expected 4 treated libraries from distinct pedigrees: ",
             gt, " G", gn, " ", tr)
      }
      full <- fit_pair(rbind(ctrl, treated))
      full_sig <- rownames(full)[!is.na(full$padj) & full$padj < .05]
      for (drop_id in twin_ids) {
        reduced <- fit_pair(rbind(ctrl[ctrl$sample_id != drop_id, ], treated))
        common <- intersect(rownames(full), rownames(reduced))
        ok <- common[is.finite(full[common, "log2FoldChange"]) &
                     is.finite(reduced[common, "log2FoldChange"])]
        reduced_sig <- rownames(reduced)[!is.na(reduced$padj) & reduced$padj < .05]
        output[[length(output) + 1L]] <- data.frame(
          genotype = gt, generation = gn, treatment = tr,
          excluded_control_sample = drop_id,
          n_control_full = nrow(ctrl), n_control_sensitivity = nrow(ctrl) - 1L,
          n_treated = nrow(treated), n_genes_compared = length(ok),
          DEGs_full = length(full_sig), DEGs_sensitivity = length(reduced_sig),
          fraction_full_DEGs_retained =
            if (length(full_sig)) length(intersect(full_sig, reduced_sig)) / length(full_sig)
            else NA_real_,
          fraction_common_genes_same_LFC_sign =
            if (length(ok)) mean(sign(full[ok, "log2FoldChange"]) ==
                                  sign(reduced[ok, "log2FoldChange"])) else NA_real_,
          spearman_common_gene_LFC =
            if (length(ok) > 2L) suppressWarnings(cor(
              full[ok, "log2FoldChange"], reduced[ok, "log2FoldChange"],
              method = "spearman")) else NA_real_,
          stringsAsFactors = FALSE
        )
        message("Finished: ", gt, " G", gn, " ", tr, " excluding ", drop_id)
      }
    }
  }
  write.csv(do.call(rbind, output),
            file.path(OUTPUT_DIR, "rna_pedigree_control_sensitivity.csv"),
            row.names = FALSE)
} else {
  message("RNA count matrix or verified RNA metadata missing; RNA analysis skipped.")
}
