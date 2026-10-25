# Sensitivity of the exact S18-S22 Persistent/Switcher classifications to
# post hoc |estimated log2FoldChange| > 0, 0.2, and 0.5.
# These are filters on estimates AFTER padj < 0.05, not DESeq2 lfcThreshold tests.
# Uses the FULL original DESeq2 per-gene results in the six folders below.
# Does not make heatmaps or install annotation/plotting packages.

root_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq"
out_dir <- file.path(root_dir, "Trajectory_LFC_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05
cutoffs <- c(0, 0.2, 0.5)
treatments <- c("Drought", "CO2", "Ozone", "Heat", "Combination")
generations <- c("G1", "G3", "G5")
genotypes <- c("WT", "gsnor-ko")
categories <- c("Persistent", "Switcher", "Unclassified")

read_full_result <- function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE)
  if (!length(first_line)) stop("Empty result file: ", path)
  sep <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else
    if (grepl(",", first_line, fixed = TRUE)) "," else ""
  x <- read.table(path, header = TRUE, row.names = 1, sep = sep,
                  check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "\"", comment.char = "")
  needed <- c("padj", "log2FoldChange")
  if (!all(needed %in% names(x)))
    stop("Expected full DESeq2 table (padj/log2FoldChange): ", path)
  id <- toupper(sub("\\..*$", "", rownames(x))) # same ID cleanup as heatmap script
  if (anyNA(id) || any(id == "") || anyDuplicated(id))
    stop("Missing/duplicated gene IDs after version stripping in ", path,
         "; resolve IDs before counting trajectories.")
  lfc <- suppressWarnings(as.numeric(as.character(x$log2FoldChange)))
  padj <- suppressWarnings(as.numeric(as.character(x$padj)))
  if (any(is.na(lfc) & !is.na(x$log2FoldChange)) ||
      any(is.na(padj) & !is.na(x$padj))) {
    stop("Nonnumeric log2FoldChange/padj in ", path)
  }
  data.frame(gene_id = id, lfc = lfc, padj = padj,
             stringsAsFactors = FALSE)
}

# The source script chooses only filenames exactly matching a treatment.
# Require one result file for each genotype x generation x treatment.
get_result_path <- function(genotype, generation, treatment) {
  folder <- paste0(if (genotype == "WT") "WT" else "gsnor",
                   "_gen", substring(generation, 2L))
  location <- file.path(root_dir, folder)
  if (!dir.exists(location)) stop("Missing folder: ", location)
  paths <- list.files(location, pattern = "\\.(txt|csv|tsv)$",
                      full.names = TRUE, ignore.case = TRUE)
  stems <- tools::file_path_sans_ext(basename(paths))
  hits <- paths[tolower(stems) == tolower(treatment)]
  if (length(hits) != 1L)
    stop("Expected exactly one full result file for ", genotype, " / ",
         generation, " / ", treatment, " in ", location,
         "; found ", length(hits), ".")
  hits
}

call_state <- function(lfc, padj, cutoff) {
  result <- rep("None", length(lfc))
  ok <- is.finite(lfc) & !is.na(padj) & padj < padj_cutoff
  result[ok & lfc > cutoff] <- "Up"
  result[ok & lfc < -cutoff] <- "Down"
  result
}

# Mirrors classify_gene() and subsequent n_sig == required_n filtering in
# Heatmaps S18-S22(1).R. Under the original rules WT CO2/ozone and gsnor-ko
# Drought require TWO non-None states; others require THREE.
trajectory_class <- function(g1, g3, g5, required_n, force_G5_None) {
  states <- cbind(g1, g3, g5)
  n_sig <- rowSums(states != "None")
  same_up <- rowSums(states == "Up") == n_sig
  same_down <- rowSums(states == "Down") == n_sig
  mixed <- rowSums(states == "Up") > 0L & rowSums(states == "Down") > 0L
  eligible <- n_sig == required_n
  if (force_G5_None) eligible <- eligible & g5 == "None"
  ans <- rep("Unclassified", length(g1))
  ans[eligible & (same_up | same_down)] <- "Persistent"
  ans[eligible & mixed] <- "Switcher"
  ans
}

input_audit <- list()
transition_rows <- list()
summary_rows <- list()
overview_rows <- list()
gene_rows <- list()
for (gt in genotypes) for (tr in treatments) {
  per_generation <- setNames(vector("list", length(generations)), generations)
  for (g in generations) {
    path <- get_result_path(gt, g, tr)
    per_generation[[g]] <- read_full_result(path)
    input_audit[[length(input_audit) + 1L]] <- data.frame(
      genotype = gt, treatment = tr, generation = g,
      source_file = path, n_gene_ids = nrow(per_generation[[g]]),
      n_testable_padj = sum(!is.na(per_generation[[g]]$padj)))
  }
  # Use the union, as the heatmap script does when widening its per-gene data.
  ids <- sort(unique(unlist(lapply(per_generation, `[[`, "gene_id"),
                            use.names = FALSE)))
  wide <- data.frame(gene_id = ids, stringsAsFactors = FALSE)
  for (g in generations) {
    z <- per_generation[[g]]
    ix <- match(ids, z$gene_id)
    wide[[paste0("lfc_", g)]] <- z$lfc[ix]
    wide[[paste0("padj_", g)]] <- z$padj[ix]
  }
  required_n <- if (gt == "WT" && tr %in% c("CO2", "Ozone")) 2L else
    if (gt == "gsnor-ko" && tr == "Drought") 2L else 3L
  force_none <- gt == "gsnor-ko" && tr == "Drought"
  for (cut in cutoffs) {
    label <- if (cut == 0) "0" else as.character(cut)
    for (g in generations)
      wide[[paste0("state_", g, "_", label)]] <-
        call_state(wide[[paste0("lfc_", g)]],
                   wide[[paste0("padj_", g)]], cut)
    s <- paste0("state_", generations, "_", label)
    wide[[paste0("pattern_", label)]] <-
      paste(wide[[s[1L]]], wide[[s[2L]]], wide[[s[3L]]], sep = "_")
    wide[[paste0("class_", label)]] <- trajectory_class(
      wide[[s[1L]]], wide[[s[2L]]], wide[[s[3L]]],
      required_n, force_none)
  }
  for (cut in cutoffs[-1L]) {
    label <- as.character(cut)
    original <- wide$class_0
    updated <- wide[[paste0("class_", label)]]
    patterns_changed <- wide$pattern_0 != wide[[paste0("pattern_", label)]]
    overview_rows[[length(overview_rows) + 1L]] <- data.frame(
      genotype = gt, treatment = tr, cutoff = cut,
      required_n = required_n,
      baseline_classified_n = sum(original != "Unclassified"),
      cutoff_classified_n = sum(updated != "Unclassified"),
      baseline_genes_changed_class_n =
        sum(original != "Unclassified" & original != updated),
      baseline_genes_changed_pattern_n =
        sum(original != "Unclassified" & patterns_changed),
      newly_classified_n =
        sum(original == "Unclassified" & updated != "Unclassified"),
      stringsAsFactors = FALSE)
    for (old_class in categories) for (new_class in categories) {
      transition_rows[[length(transition_rows) + 1L]] <- data.frame(
        genotype = gt, treatment = tr, cutoff = cut,
        required_n = required_n, force_G5_None = force_none,
        baseline_class = old_class, new_class = new_class,
        n_genes = sum(original == old_class & updated == new_class),
        stringsAsFactors = FALSE)
    }
    for (cl in categories[1:2]) {
      was <- original == cl
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        genotype = gt, treatment = tr, cutoff = cut,
        class = cl, required_n = required_n,
        baseline_n = sum(was),
        same_class_n = sum(was & updated == cl),
        same_class_same_pattern_n = sum(was & updated == cl & !patterns_changed),
        same_class_changed_pattern_n = sum(was & updated == cl & patterns_changed),
        switched_to_other_class_n = sum(was & updated %in% setdiff(categories[1:2], cl)),
        lost_classification_n = sum(was & updated == "Unclassified"),
        cutoff_class_n = sum(updated == cl),
        newly_classified_from_unclassified_n =
          sum(original == "Unclassified" & updated == cl),
        stringsAsFactors = FALSE)
    }
    # One row per gene originally or newly classified, at each new cutoff.
    keep <- original != "Unclassified" | updated != "Unclassified"
    if (any(keep)) {
      v <- wide[keep, c("gene_id", paste0("lfc_", generations),
                       paste0("padj_", generations),
                       paste0("state_", generations, "_0"),
                       paste0("state_", generations, "_", label),
                       "pattern_0", paste0("pattern_", label),
                       "class_0", paste0("class_", label))]
      names(v) <- sub(paste0("_", label, "$"), "_cutoff", names(v))
      names(v) <- sub("_0$", "_baseline", names(v))
      gene_rows[[length(gene_rows) + 1L]] <- data.frame(
        genotype = gt, treatment = tr, cutoff = cut,
        required_n = required_n, force_G5_None = force_none,
        v, changed_class = original[keep] != updated[keep],
        changed_pattern = patterns_changed[keep],
        stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  message("Processed ", gt, " / ", tr, "; ", length(ids), " gene IDs")
}

audit <- do.call(rbind, input_audit)
transitions <- do.call(rbind, transition_rows)
summary <- do.call(rbind, summary_rows)
overview <- do.call(rbind, overview_rows)
per_gene <- do.call(rbind, gene_rows)
write.csv(audit, file.path(out_dir, "trajectory_input_audit.csv"), row.names = FALSE)
write.csv(transitions, file.path(out_dir, "trajectory_transition_counts.csv"),
          row.names = FALSE)
write.csv(summary, file.path(out_dir, "trajectory_sensitivity_by_class.csv"),
          row.names = FALSE)
write.csv(overview, file.path(out_dir, "trajectory_sensitivity_overview.csv"),
          row.names = FALSE)
write.csv(per_gene, file.path(out_dir, "trajectory_gene_changes.csv"),
          row.names = FALSE)

# Consistency checks; changed state within the same class is counted separately.
stopifnot(all(summary$baseline_n == summary$same_class_n +
              summary$switched_to_other_class_n + summary$lost_classification_n))
writeLines(c(
  "Source: full per-gene DESeq2 results used by Heatmaps S18-S22(1).R.",
  "Statuses: Up/Down if padj < 0.05 and estimated log2FC > cutoff / < -cutoff; otherwise None.",
  "Persistent means same direction among required significant sampled points; Switcher means both directions.",
  "Original figure exception: WT CO2/Ozone and gsnor-ko Drought require two sampled points; all others require G1, G3, G5.",
  "For gsnor-ko Drought, G5 must be None, exactly as in the supplied script.",
  "These are post hoc estimated effect-size filters, NOT tests of |true log2FC| > cutoff.",
  "Same-class changed-pattern genes are included in same_class_n, and listed separately.",
  "trajectory_sensitivity_overview.csv answers how many baseline-classified genes changed class or exact three-point pattern.",
  "Input files must match those used for the published heatmaps; compare baseline counts against the original figure table."
), file.path(out_dir, "README_trajectory_sensitivity.txt"))
print(summary)
message("Wrote sensitivity tables to ", normalizePath(out_dir))
