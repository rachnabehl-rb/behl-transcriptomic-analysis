## ============================================================
## MODULE RANKING AFTER WGCNA
## Use after fresh WGCNA pipeline
##
## Required saved objects:
##   datExpr_final.rds
##   MEs0_final.rds
##   m_vsd_group_median.rds
##   mergedColors_final.rds
##
## Main outputs:
##   MODULE_RANKING_ANOVA_long.csv
##   MODULE_RANKING_ANOVA_wide_summary.csv
##   MODULE_RANKING_Treatment_summary.csv
##   MODULE_RANKING_Genotype_summary.csv
##   MODULE_RANKING_Generation_summary.csv
##   MODULE_RANKING_FINAL_combined.csv
##   MODULES_RECOMMENDED_FOR_FOCUS.csv
##   HEATMAP_ranked_modules_Treatment_Genotype_Generation.pdf
## ============================================================

rm(list = ls())
set.seed(65)

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(tibble)
  library(stringr)
})

options(stringsAsFactors = FALSE)

## ============================================================
## 1. PATHS
## ============================================================

wd <- "~/Documents/DNA_seq/RNA_seq/WGCNA"

## Change this if your WGCNA output folder has another name
wgcna_outdir <- file.path(wd, "WGCNA_final")

ranking_outdir <- file.path(wgcna_outdir, "MODULE_RANKING")
dir.create(ranking_outdir, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## 2. LOAD WGCNA OBJECTS
## ============================================================

datExpr <- readRDS(file.path(wgcna_outdir, "datExpr_final.rds"))
MEs0 <- readRDS(file.path(wgcna_outdir, "MEs0_final.rds"))
m_vsd <- readRDS(file.path(wgcna_outdir, "m_vsd_group_median.rds"))
mergedColors <- readRDS(file.path(wgcna_outdir, "mergedColors_final.rds"))

trait_cols <- c("Area week 3", "Area week 4", "Weight.of.seeds")

cat("\nLoaded objects:\n")
cat("datExpr dimensions:\n")
print(dim(datExpr))

cat("MEs0 dimensions:\n")
print(dim(MEs0))

cat("m_vsd dimensions:\n")
print(dim(m_vsd))

cat("Number of module colors:\n")
print(length(mergedColors))

## ============================================================
## 3. PREPARE METADATA ALIGNED WITH MODULE EIGENGENES
## ============================================================

meta_rank <- m_vsd %>%
  dplyr::select(
    Genotype,
    Treatment,
    Generation,
    SampleName,
    dplyr::all_of(trait_cols)
  )

rownames(meta_rank) <- meta_rank$SampleName

meta_rank <- meta_rank[rownames(datExpr), , drop = FALSE]

stopifnot(all(rownames(meta_rank) == rownames(datExpr)))
stopifnot(all(rownames(MEs0) == rownames(datExpr)))

meta_rank <- meta_rank %>%
  mutate(
    Genotype = str_squish(as.character(Genotype)),
    Treatment = str_squish(as.character(Treatment)),
    Generation = str_squish(as.character(Generation))
  )

meta_rank$Genotype <- factor(
  meta_rank$Genotype,
  levels = c("WT", "gsnor")
)

meta_rank$Treatment <- factor(
  meta_rank$Treatment,
  levels = c("Control", "Drought", "CO2", "Ozone", "Heat", "Combination")
)

meta_rank$Generation <- factor(
  meta_rank$Generation,
  levels = c("1", "3", "5")
)

## Remove grey module from ranking
MEs_rank <- MEs0[, colnames(MEs0) != "MEgrey", drop = FALSE]

cat("\nSamples used for ranking:\n")
print(nrow(meta_rank))

cat("\nModules used for ranking, excluding grey:\n")
print(ncol(MEs_rank))

cat("\nTreatment counts:\n")
print(table(meta_rank$Treatment))

cat("\nGenotype counts:\n")
print(table(meta_rank$Genotype))

cat("\nGeneration counts:\n")
print(table(meta_rank$Generation))

## ============================================================
## 4. ANOVA RANKING
## Model: ME ~ Genotype + Treatment + Generation
## ============================================================

run_module_anova <- function(module_name, MEs_rank, meta_rank) {
  
  df <- data.frame(
    ME = as.numeric(MEs_rank[, module_name]),
    Genotype = meta_rank$Genotype,
    Treatment = meta_rank$Treatment,
    Generation = meta_rank$Generation
  )
  
  fit <- lm(ME ~ Genotype + Treatment + Generation, data = df)
  
  a <- anova(fit)
  
  total_ss <- sum(a$`Sum Sq`, na.rm = TRUE)
  residual_ss <- a["Residuals", "Sum Sq"]
  
  out <- a %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term != "Residuals") %>%
    mutate(
      module = sub("^ME", "", module_name),
      sum_sq = `Sum Sq`,
      mean_sq = `Mean Sq`,
      F_value = `F value`,
      p_value = `Pr(>F)`,
      eta2 = sum_sq / total_ss,
      partial_eta2 = sum_sq / (sum_sq + residual_ss)
    ) %>%
    dplyr::select(
      module,
      term,
      sum_sq,
      eta2,
      partial_eta2,
      F_value,
      p_value
    )
  
  return(out)
}

anova_rank_long <- purrr::map_dfr(
  colnames(MEs_rank),
  run_module_anova,
  MEs_rank = MEs_rank,
  meta_rank = meta_rank
)

anova_rank_long <- anova_rank_long %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(term, q_value, desc(partial_eta2))

write.csv(
  anova_rank_long,
  file.path(ranking_outdir, "MODULE_RANKING_ANOVA_long.csv"),
  row.names = FALSE
)

cat("\nTop Treatment modules by ANOVA:\n")
anova_rank_long %>%
  filter(term == "Treatment") %>%
  arrange(q_value, desc(partial_eta2)) %>%
  print(n = 20)

cat("\nTop Genotype modules by ANOVA:\n")
anova_rank_long %>%
  filter(term == "Genotype") %>%
  arrange(q_value, desc(partial_eta2)) %>%
  print(n = 20)

cat("\nTop Generation modules by ANOVA:\n")
anova_rank_long %>%
  filter(term == "Generation") %>%
  arrange(q_value, desc(partial_eta2)) %>%
  print(n = 20)

## ============================================================
## 5. WIDE ANOVA SUMMARY TABLE
## ============================================================

anova_rank_wide <- anova_rank_long %>%
  dplyr::select(module, term, partial_eta2, p_value, q_value) %>%
  pivot_wider(
    names_from = term,
    values_from = c(partial_eta2, p_value, q_value),
    names_sep = "_"
  )

anova_rank_wide <- anova_rank_wide %>%
  rowwise() %>%
  mutate(
    dominant_factor = c("Genotype", "Treatment", "Generation")[
      which.max(c(
        partial_eta2_Genotype,
        partial_eta2_Treatment,
        partial_eta2_Generation
      ))
    ],
    dominant_eta2 = max(
      c(
        partial_eta2_Genotype,
        partial_eta2_Treatment,
        partial_eta2_Generation
      ),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  arrange(desc(dominant_eta2))

write.csv(
  anova_rank_wide,
  file.path(ranking_outdir, "MODULE_RANKING_ANOVA_wide_summary.csv"),
  row.names = FALSE
)

cat("\nOverall ANOVA module ranking:\n")
print(anova_rank_wide, n = 30)

## ============================================================
## 6. BINARIZED FACTOR RANKING
## This ranks modules by strongest correlation with each level
## ============================================================

## ============================================================
## FIXED BINARIZED FACTOR RANKING FUNCTION
## Robust version: does not depend on Var1 / Var2 / Freq names
## ============================================================

make_factor_ranking <- function(fvec, factor_name, MEs_rank, ranking_outdir) {
  
  fvec <- as.character(fvec)
  
  ## Binarize factor levels
  bin_mat <- WGCNA::binarizeCategoricalVariable(
    fvec,
    includePairwise = FALSE,
    includeLevelVsAll = TRUE,
    minCount = 1,
    dropFirstLevelVsAll = FALSE
  )
  
  rownames(bin_mat) <- rownames(MEs_rank)
  
  stopifnot(nrow(bin_mat) == nrow(MEs_rank))
  stopifnot(all(rownames(bin_mat) == rownames(MEs_rank)))
  
  ## Correlation matrix: modules x factor levels
  cor_mat <- cor(
    MEs_rank,
    bin_mat,
    use = "p",
    method = "pearson"
  )
  
  p_mat <- WGCNA::corPvalueFisher(
    cor_mat,
    nrow(MEs_rank)
  )
  
  q_mat <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )
  
  ## Robust long-format conversion
  long_df <- expand.grid(
    module = rownames(cor_mat),
    level = colnames(cor_mat),
    stringsAsFactors = FALSE
  )
  
  long_df$cor <- as.vector(cor_mat)
  long_df$p_value <- as.vector(p_mat)
  long_df$q_value <- as.vector(q_mat)
  
  long_df <- long_df %>%
    mutate(
      module = sub("^ME", "", module),
      factor = factor_name,
      abs_cor = abs(cor)
    ) %>%
    arrange(q_value, desc(abs_cor))
  
  ## One summary row per module
  summary_df <- long_df %>%
    group_by(module, factor) %>%
    summarise(
      strongest_level = level[which.max(abs_cor)],
      max_abs_cor = max(abs_cor, na.rm = TRUE),
      signed_cor_at_strongest = cor[which.max(abs_cor)],
      min_p_value = min(p_value, na.rm = TRUE),
      min_q_value = min(q_value, na.rm = TRUE),
      n_levels_q_0.05 = sum(q_value <= 0.05, na.rm = TRUE),
      n_levels_q_0.10 = sum(q_value <= 0.10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(min_q_value, desc(max_abs_cor))
  
  ## Save outputs
  write.csv(
    long_df,
    file.path(ranking_outdir, paste0("MODULE_RANKING_", factor_name, "_long.csv")),
    row.names = FALSE
  )
  
  write.csv(
    summary_df,
    file.path(ranking_outdir, paste0("MODULE_RANKING_", factor_name, "_summary.csv")),
    row.names = FALSE
  )
  
  return(summary_df)
}

treatment_rank <- make_factor_ranking(
  fvec = meta_rank$Treatment,
  factor_name = "Treatment",
  MEs_rank = MEs_rank,
  ranking_outdir = ranking_outdir
)

genotype_rank <- make_factor_ranking(
  fvec = meta_rank$Genotype,
  factor_name = "Genotype",
  MEs_rank = MEs_rank,
  ranking_outdir = ranking_outdir
)

generation_rank <- make_factor_ranking(
  fvec = meta_rank$Generation,
  factor_name = "Generation",
  MEs_rank = MEs_rank,
  ranking_outdir = ranking_outdir
)

cat("\nTop treatment-associated modules by binarized correlation:\n")
print(treatment_rank, n = 30)

cat("\nTop genotype-associated modules by binarized correlation:\n")
print(genotype_rank, n = 30)

cat("\nTop generation-associated modules by binarized correlation:\n")
print(generation_rank, n = 30)

## ============================================================
## 7. PHENOTYPE RANKING - FIXED VERSION
## Robust version: does not depend on Var1 / Var2 / Freq names
## ============================================================

trait_mat <- meta_rank[, trait_cols, drop = FALSE]
trait_mat[] <- lapply(trait_mat, function(x) suppressWarnings(as.numeric(x)))

pheno_cor <- cor(
  MEs_rank,
  trait_mat,
  use = "p",
  method = "pearson"
)

pheno_p <- WGCNA::corPvalueFisher(
  pheno_cor,
  nrow(MEs_rank)
)

pheno_q <- matrix(
  p.adjust(as.vector(pheno_p), method = "BH"),
  nrow = nrow(pheno_p),
  ncol = ncol(pheno_p),
  dimnames = dimnames(pheno_p)
)

## Robust long-format conversion
pheno_rank_long <- expand.grid(
  module = rownames(pheno_cor),
  phenotype = colnames(pheno_cor),
  stringsAsFactors = FALSE
)

pheno_rank_long$cor <- as.vector(pheno_cor)
pheno_rank_long$p_value <- as.vector(pheno_p)
pheno_rank_long$q_value <- as.vector(pheno_q)

pheno_rank_long <- pheno_rank_long %>%
  mutate(
    module = sub("^ME", "", module),
    abs_cor = abs(cor)
  ) %>%
  arrange(q_value, desc(abs_cor))

pheno_rank_summary <- pheno_rank_long %>%
  group_by(module) %>%
  summarise(
    strongest_phenotype = phenotype[which.max(abs_cor)],
    max_abs_pheno_cor = max(abs_cor, na.rm = TRUE),
    signed_cor_at_strongest_phenotype = cor[which.max(abs_cor)],
    min_pheno_p = min(p_value, na.rm = TRUE),
    min_pheno_q = min(q_value, na.rm = TRUE),
    n_phenotypes_q_0.05 = sum(q_value <= 0.05, na.rm = TRUE),
    n_phenotypes_q_0.10 = sum(q_value <= 0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_pheno_q, desc(max_abs_pheno_cor))

write.csv(
  pheno_rank_long,
  file.path(ranking_outdir, "MODULE_RANKING_Phenotypes_long.csv"),
  row.names = FALSE
)

write.csv(
  pheno_rank_summary,
  file.path(ranking_outdir, "MODULE_RANKING_Phenotypes_summary.csv"),
  row.names = FALSE
)

cat("\nTop phenotype-associated modules:\n")
print(pheno_rank_summary, n = 30)

## ============================================================
## 8. FINAL COMBINED RANKING TABLE
## ============================================================

combined_rank <- anova_rank_wide %>%
  left_join(
    treatment_rank %>%
      dplyr::select(
        module,
        treatment_strongest_level = strongest_level,
        treatment_max_abs_cor = max_abs_cor,
        treatment_signed_cor = signed_cor_at_strongest,
        treatment_min_p = min_p_value,
        treatment_min_q = min_q_value
      ),
    by = "module"
  ) %>%
  left_join(
    genotype_rank %>%
      dplyr::select(
        module,
        genotype_strongest_level = strongest_level,
        genotype_max_abs_cor = max_abs_cor,
        genotype_signed_cor = signed_cor_at_strongest,
        genotype_min_p = min_p_value,
        genotype_min_q = min_q_value
      ),
    by = "module"
  ) %>%
  left_join(
    generation_rank %>%
      dplyr::select(
        module,
        generation_strongest_level = strongest_level,
        generation_max_abs_cor = max_abs_cor,
        generation_signed_cor = signed_cor_at_strongest,
        generation_min_p = min_p_value,
        generation_min_q = min_q_value
      ),
    by = "module"
  ) %>%
  left_join(
    pheno_rank_summary,
    by = "module"
  ) %>%
  arrange(desc(dominant_eta2))

write.csv(
  combined_rank,
  file.path(ranking_outdir, "MODULE_RANKING_FINAL_combined.csv"),
  row.names = FALSE
)

cat("\nFinal combined module ranking:\n")
print(combined_rank, n = 50)

## ============================================================
## 9. PRINT FOCUS MODULES
## ============================================================

cat("\n============================================================\n")
cat("Treatment modules passing q <= 0.05 by ANOVA:\n")
cat("============================================================\n")

treat_sig_005 <- anova_rank_long %>%
  filter(term == "Treatment", q_value <= 0.05) %>%
  arrange(q_value, desc(partial_eta2))

print(treat_sig_005, n = 50)

cat("\n============================================================\n")
cat("Treatment modules passing q <= 0.10 by ANOVA:\n")
cat("============================================================\n")

treat_sig_010 <- anova_rank_long %>%
  filter(term == "Treatment", q_value <= 0.10) %>%
  arrange(q_value, desc(partial_eta2))

print(treat_sig_010, n = 50)

cat("\n============================================================\n")
cat("Top 10 treatment modules by effect size, regardless of FDR:\n")
cat("============================================================\n")

top_treatment_by_effect <- anova_rank_long %>%
  filter(term == "Treatment") %>%
  arrange(desc(partial_eta2), p_value) %>%
  slice_head(n = 10)

print(top_treatment_by_effect, n = 10)

cat("\n============================================================\n")
cat("Top 10 treatment modules by strongest treatment-level correlation:\n")
cat("============================================================\n")

top_treatment_by_cor <- treatment_rank %>%
  arrange(desc(max_abs_cor), min_q_value) %>%
  slice_head(n = 10)

print(top_treatment_by_cor, n = 10)

cat("\n============================================================\n")
cat("Top 10 genotype modules by effect size:\n")
cat("============================================================\n")

top_genotype_by_effect <- anova_rank_long %>%
  filter(term == "Genotype") %>%
  arrange(desc(partial_eta2), p_value) %>%
  slice_head(n = 10)

print(top_genotype_by_effect, n = 10)

cat("\n============================================================\n")
cat("Top 10 generation modules by effect size:\n")
cat("============================================================\n")

top_generation_by_effect <- anova_rank_long %>%
  filter(term == "Generation") %>%
  arrange(desc(partial_eta2), p_value) %>%
  slice_head(n = 10)

print(top_generation_by_effect, n = 10)

recommended_modules <- unique(c(
  top_treatment_by_effect$module,
  top_treatment_by_cor$module,
  top_genotype_by_effect$module[1:min(5, nrow(top_genotype_by_effect))],
  top_generation_by_effect$module[1:min(5, nrow(top_generation_by_effect))],
  pheno_rank_summary$module[1:min(5, nrow(pheno_rank_summary))]
))

recommended_modules <- recommended_modules[!is.na(recommended_modules)]

cat("\n============================================================\n")
cat("Recommended modules to focus on:\n")
cat("============================================================\n")
print(recommended_modules)

write.csv(
  data.frame(module = recommended_modules),
  file.path(ranking_outdir, "MODULES_RECOMMENDED_FOR_FOCUS.csv"),
  row.names = FALSE
)

## ============================================================
## 10. HEATMAP OF RANKING SCORES
## ============================================================

ranking_heatmap_df <- combined_rank %>%
  dplyr::select(
    module,
    partial_eta2_Genotype,
    partial_eta2_Treatment,
    partial_eta2_Generation,
    genotype_max_abs_cor,
    treatment_max_abs_cor,
    generation_max_abs_cor,
    max_abs_pheno_cor
  ) %>%
  column_to_rownames("module")

ranking_heatmap_mat <- as.matrix(ranking_heatmap_df)

pdf(
  file.path(ranking_outdir, "HEATMAP_ranked_modules_Treatment_Genotype_Generation.pdf"),
  width = 10,
  height = 12,
  useDingbats = FALSE
)

par(mar = c(10, 12, 4, 4))

WGCNA::labeledHeatmap(
  Matrix = ranking_heatmap_mat,
  xLabels = colnames(ranking_heatmap_mat),
  yLabels = rownames(ranking_heatmap_mat),
  ySymbols = rownames(ranking_heatmap_mat),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = signif(ranking_heatmap_mat, 2),
  cex.text = 0.7,
  cex.lab.x = 0.9,
  cex.lab.y = 0.9,
  xLabelsAngle = 45,
  setStdMargins = FALSE,
  main = "Module ranking scores"
)

dev.off()

## ============================================================
## 11. FINAL MESSAGE
## ============================================================

cat("\n============================================================\n")
cat("MODULE RANKING FINISHED SUCCESSFULLY\n")
cat("Output directory:\n")
print(ranking_outdir)

cat("\nImportant output files:\n")
cat("1. MODULE_RANKING_FINAL_combined.csv\n")
cat("2. MODULE_RANKING_ANOVA_long.csv\n")
cat("3. MODULE_RANKING_Treatment_summary.csv\n")
cat("4. MODULE_RANKING_Genotype_summary.csv\n")
cat("5. MODULE_RANKING_Generation_summary.csv\n")
cat("6. MODULES_RECOMMENDED_FOR_FOCUS.csv\n")
cat("7. HEATMAP_ranked_modules_Treatment_Genotype_Generation.pdf\n")
cat("============================================================\n")

