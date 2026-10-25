# Run the original "Venn Fig3C.R" in R first, then run/source this file
# in the SAME R session. Alternatively, paste this block after section 14
# of the original script. It does not redraw or modify the Venn diagrams.
#
# Null hypothesis: within each genotype/treatment, each generation's DEG
# set is sampled at its observed size from genes eligible for a DEG call
# IN THAT GENERATION. DESeq2's independent filtering can produce padj=NA
# for a gene in one generation but a valid DEG in another, so do not force
# a common padj-defined background or alter the published Venn counts.
# This tests enrichment of the THREE-WAY intersection, not WT vs gsnor-ko.
# A Monte Carlo test with B=20000 resolves p-values down to 1/(B+1).
set.seed(20260925)
B <- 20000L

needed_objects <- c(
  "file_map_df", "deg_gene_lists", "root_dir", "out_dir",
  "genotype_order", "treatment_order", "generation_order",
  "read_deseq_file", "get_venn_sets"
)
missing_objects <- needed_objects[
  !vapply(needed_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop("Run Venn Fig3C.R first. Missing: ", paste(missing_objects, collapse = ", "))
}

# One-sided permutation test with the three original eligible-gene pools.
three_way_p <- function(pools, sizes, observed, B) {
  if (observed == 0L) return(1)
  at_least_observed <- 0L
  for (i in seq_len(B)) {
    a <- sample(pools[["G1"]], sizes[["G1"]], replace = FALSE)
    b <- sample(pools[["G3"]], sizes[["G3"]], replace = FALSE)
    c <- sample(pools[["G5"]], sizes[["G5"]], replace = FALSE)
    simulated <- length(intersect(intersect(a, b), c))
    if (simulated >= observed) at_least_observed <- at_least_observed + 1L
  }
  (at_least_observed + 1) / (B + 1)
}

read_eligible_ids <- function(genotype, treatment, generation) {
  folder <- if (genotype == "WT") {
    paste0("WT_gen", generation)
  } else {
    paste0("gsnor_gen", generation)
  }
  matches <- file_map_df[
    !is.na(file_map_df$treatment) &
      file_map_df$folder == folder &
      file_map_df$treatment == treatment,
    , drop = FALSE
  ]
  if (nrow(matches) != 1L) {
    stop("Expected exactly one full DESeq2 result for ", folder, "/",
         treatment, "; found ", nrow(matches))
  }
  file <- file.path(root_dir, folder, matches$file[[1L]])
  if (!file.exists(file)) stop("Missing full result file: ", file)
  dat <- read_deseq_file(file)
  if (!all(c("padj", "log2FoldChange", "gene_id") %in% names(dat))) {
    stop("Missing DESeq2 columns in ", file)
  }
  ids <- dat$gene_id[
    !is.na(suppressWarnings(as.numeric(dat$padj))) &
      !is.na(suppressWarnings(as.numeric(dat$log2FoldChange)))
  ]
  unique(as.character(ids[!is.na(ids) & nzchar(ids)]))
}

rows <- list()
for (genotype in genotype_order) {
  for (treatment in treatment_order) {
    eligible_by_generation <- lapply(generation_order, function(g) {
      read_eligible_ids(genotype, treatment, g)
    })
    names(eligible_by_generation) <- c("G1", "G3", "G5")
    common_eligible <- Reduce(intersect, eligible_by_generation)
    N_common <- length(common_eligible)
    N_by_generation <- vapply(eligible_by_generation, length, integer(1))
    if (any(N_by_generation == 0L)) {
      stop("No eligible genes in a generation for ", genotype, "/", treatment)
    }

    sets <- get_venn_sets(genotype, treatment)
    # Check each Venn set against ITS OWN generation's eligible genes.
    outside <- vapply(names(sets), function(g) {
      length(setdiff(sets[[g]], eligible_by_generation[[g]]))
    }, integer(1))
    if (any(outside > 0L)) {
      stop(
        "Figure 3C has DEGs outside their own tested-gene background for ",
        genotype, "/", treatment, ": ",
        paste(names(outside), outside, sep = "=", collapse = ", "),
        ". Check the DEG filter and result files."
      )
    }

    sizes <- vapply(sets, function(s) length(unique(s)), integer(1))
    x <- length(Reduce(intersect, lapply(sets, unique)))
    expected <- N_common * prod(as.double(sizes) / N_by_generation)
    p <- three_way_p(eligible_by_generation, sizes, x, B)

    rows[[length(rows) + 1L]] <- data.frame(
      genotype = genotype,
      treatment = treatment,
      eligible_G1_N = N_by_generation[["G1"]],
      eligible_G3_N = N_by_generation[["G3"]],
      eligible_G5_N = N_by_generation[["G5"]],
      eligible_common_N = N_common,
      G1_total = sizes[["G1"]],
      G3_total = sizes[["G3"]],
      G5_total = sizes[["G5"]],
      shared_G1_G3_G5 = x,
      pct_of_G1 = if (sizes[["G1"]] > 0) 100 * x / sizes[["G1"]] else NA_real_,
      pct_of_G3 = if (sizes[["G3"]] > 0) 100 * x / sizes[["G3"]] else NA_real_,
      pct_of_G5 = if (sizes[["G5"]] > 0) 100 * x / sizes[["G5"]] else NA_real_,
      expected_shared = expected,
      observed_over_expected = if (expected > 0) x / expected else NA_real_,
      p_enrichment = p,
      permutations = B,
      stringsAsFactors = FALSE
    )
  }
}

overlap_stats <- do.call(rbind, rows)
overlap_stats$padj_BH <- p.adjust(overlap_stats$p_enrichment, method = "BH")
output_file <- file.path(out_dir, "Figure3C_three_way_overlap_statistics.csv")
write.csv(overlap_stats, output_file, row.names = FALSE, na = "")

cat("\nAll comparisons (percentages use each generation's own DEG total):\n")
print(overlap_stats, row.names = FALSE, digits = 4)
cat("\nExamples for Results text:\n")
print(subset(overlap_stats,
             (genotype == "WT" & treatment == "Heat") |
               (genotype == "gsnor-ko" & treatment == "Combination")),
      row.names = FALSE, digits = 4)
cat("\nSaved: ", output_file, "\n", sep = "")
cat("Minimum resolvable Monte Carlo p-value: ", 1 / (B + 1), "\n", sep = "")
