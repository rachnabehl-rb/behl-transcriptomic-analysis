# Build the supplementary 3 x 5 interaction table from the DESeq2 summary.
# Run from the folder containing 'genotype_treatment_interaction_results'
# or 'upload', or set INPUT_FILE to its full path below.
setwd("~/Documents/DNA_seq/RNA_seq/genotype_treatment_interaction_results")
candidates <- c(
  "genotype_treatment_interaction_results/interaction_summary_all_generations.csv",
  "upload/interaction_summary_all_generations.csv",
  "interaction_summary_all_generations.csv"
)
available <- candidates[file.exists(candidates)]
if (!length(available)) {
  stop("Cannot find interaction_summary_all_generations.csv. Edit candidates with its location.")
}
INPUT_FILE <- available[1L]
OUTPUT_DIR <- "supplementary_interaction_table"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("generation", "treatment", "genes_with_testable_padj",
              "interaction_DEGs_padj_0.05")
if (!all(required %in% names(d))) {
  stop("Missing columns: ", paste(setdiff(required, names(d)), collapse = ", "))
}

d$generation <- sub("^[Gg]", "", as.character(d$generation))
gens <- c("1", "3", "5")
trts <- c("Drought", "CO2", "Ozone", "Heat", "Combination")
expected_keys <- as.vector(outer(gens, trts, paste, sep = "|"))
observed_keys <- paste(d$generation, d$treatment, sep = "|")
if (anyDuplicated(observed_keys) ||
    !setequal(observed_keys, expected_keys)) {
  stop("Expected exactly one summary row per generation and treatment (15 rows).")
}
if (anyNA(d[, required]) ||
    any(d$interaction_DEGs_padj_0.05 < 0) ||
    any(d$genes_with_testable_padj <= 0) ||
    any(d$interaction_DEGs_padj_0.05 > d$genes_with_testable_padj)) {
  stop("Missing or inconsistent counts in the summary.")
}

# Each cell: significant interaction genes / genes with testable adjusted P.
table_out <- data.frame(Generation = paste0("G", gens),
                        stringsAsFactors = FALSE, check.names = FALSE)
for (tr in trts) {
  cell <- vapply(gens, function(g) {
    r <- d[d$generation == g & d$treatment == tr, , drop = FALSE]
    paste0(format(r$interaction_DEGs_padj_0.05,
                  big.mark = ",", scientific = FALSE, trim = TRUE),
           " / ",
           format(r$genes_with_testable_padj,
                  big.mark = ",", scientific = FALSE, trim = TRUE))
  }, character(1))
  table_out[[tr]] <- cell
}

write.csv(table_out,
  file.path(OUTPUT_DIR, "Supplementary_genotype_treatment_interaction_table.csv"),
  row.names = FALSE, na = "")

# Markdown can be pasted into a response letter or converted to a Word table.
markdown <- c(
  paste0("| ", paste(names(table_out), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(table_out)), collapse = " | "), " |"),
  apply(table_out, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
  "",
  "Entries: genes with significant genotype x treatment interaction (BH-adjusted P < 0.05) / genes with a testable adjusted P value for that contrast.",
  "G5 treatment labels refer to prior treatment history; G4 and G5 plants were grown under control conditions."
)
writeLines(markdown,
  file.path(OUTPUT_DIR, "Supplementary_genotype_treatment_interaction_table.md"))
cat(paste(markdown, collapse = "\n"), "\n")
message("Saved CSV and Markdown table in: ", normalizePath(OUTPUT_DIR))
