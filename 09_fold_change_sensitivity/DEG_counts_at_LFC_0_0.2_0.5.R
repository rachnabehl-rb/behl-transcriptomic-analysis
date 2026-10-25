# Refit the original within-genotype treatment-versus-control DESeq2 contrasts
# and count significant genes at three post-hoc absolute log2FC cutoffs.
# Run from the directory containing upload/; edit INPUT_DIR if needed.
# This preserves the original adjusted-P values; it is NOT a formal test that
# the true fold change exceeds 0.2 or 0.5 and does not use shrunken log2FCs.

INPUT_DIR <- "~/Documents/DNA_seq/RNA_seq"
COUNT_FILE <- file.path(INPUT_DIR, "genes_gtf_merged_gene_counts.txt")
META_FILE <- file.path(INPUT_DIR, "Metadata_Sample_Area_medians_WITH_GROUPS.csv")
OUT_DIR <- "DEG_LFC_cutoff_sensitivity"
GENS <- c(1L, 3L, 5L)
GENOTYPES <- c("WT", "gsnor")
TREATMENTS <- c("Drought", "CO2", "Ozone", "Heat", "Combination")
LFC_CUTOFFS <- c(0, 0.2, 0.5)

if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 first.")
if (!all(file.exists(c(COUNT_FILE, META_FILE)))) stop("Input count/metadata file missing.")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

m <- read.csv(META_FILE, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("SampleName", "Genotype", "Treatment", "Generation")
if (!all(required %in% names(m))) stop("Missing metadata columns.")
m$Generation <- as.integer(sub("^[Gg]", "", as.character(m$Generation)))
if (anyDuplicated(m$SampleName)) stop("Duplicate sample names in metadata.")

raw <- as.matrix(read.table(COUNT_FILE, header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE, quote = "", comment.char = ""))
if (!all(m$SampleName %in% colnames(raw))) stop("Samples absent from count matrix.")
if (anyDuplicated(rownames(raw)) || !is.numeric(raw) || anyNA(raw) ||
    any(raw < 0) || any(abs(raw - round(raw)) > 1e-8)) {
  stop("Expected unique gene IDs and nonnegative integer counts.")
}
storage.mode(raw) <- "integer"

# Copied from the supplied DESeq2_WT_gen1 script: only WT G1 has documented
# within-genotype batch-2 sample assignments in the supplied analysis files.
wt_g1_batch2 <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D",
  "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")
wt_g1_names <- m$SampleName[m$Generation == 1L & m$Genotype == "WT"]
if (length(wt_g1_names) != 26L || !all(wt_g1_batch2 %in% wt_g1_names)) {
  stop("The WT G1 batch map does not match the supplied metadata.")
}

results_rows <- list()
for (g in GENS) for (gt in GENOTYPES) for (tr in TREATMENTS) {
  sub <- m[m$Generation == g & m$Genotype == gt &
             m$Treatment %in% c("Control", tr), , drop = FALSE]
  if (anyNA(sub$SampleName) || nrow(sub) < 6L ||
      !all(sub$SampleName %in% colnames(raw))) stop("Invalid sample subset.")
  sub$condition <- factor(ifelse(sub$Treatment == "Control", "C", "S"),
                          levels = c("C", "S"))
  if (any(table(sub$condition) < 2L)) stop("Too few replicates: G",g," ",gt," ",tr)
  if (gt == "WT" && g == 1L) {
    sub$batch <- factor(ifelse(sub$SampleName %in% wt_g1_batch2,
                               "batch2", "batch1"))
    design_formula <- ~ batch + condition
  } else {
    design_formula <- ~ condition
  }
  rownames(sub) <- sub$SampleName
  y <- raw[, sub$SampleName, drop = FALSE]
  if (!identical(colnames(y), rownames(sub))) stop("Count/metadata order mismatch.")
  if (qr(model.matrix(design_formula, sub))$rank !=
      ncol(model.matrix(design_formula, sub))) {
    stop("Design not estimable: G",g," ",gt," ",tr)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(y, sub, design = design_formula)
  dds <- dds[rowSums(DESeq2::counts(dds)) > 10, ]
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- dds[rowSums(DESeq2::counts(dds, normalized = TRUE) >= 3) >= 2, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("condition", "S", "C"),
                          alpha = 0.05)
  fc <- res$log2FoldChange
  sig <- !is.na(res$padj) & res$padj < 0.05 & !is.na(fc)

  for (cut in LFC_CUTOFFS) {
    # At zero, count all padj-significant genes; never count both Up and Down.
    up <- sig & fc > cut
    down <- sig & fc < -cut
    n_total <- if (cut == 0) sum(sig) else sum(up | down)
    n_zero_fc <- if (cut == 0) sum(sig & fc == 0) else 0L
    results_rows[[length(results_rows) + 1L]] <- data.frame(
      generation = paste0("G",g), genotype = gt, treatment = tr,
      abs_log2FC_cutoff = cut, padj_cutoff = 0.05,
      genes_with_testable_padj = sum(!is.na(res$padj)),
      DEG_total = n_total, DEG_up = sum(up), DEG_down = sum(down),
      DEG_exactly_zero_log2FC = n_zero_fc,
      n_control_libraries = sum(sub$condition == "C"),
      n_treated_libraries = sum(sub$condition == "S"),
      stringsAsFactors = FALSE)
  }
  message("Finished G",g," / ",gt," / ",tr)
}
tbl <- do.call(rbind, results_rows)
if (nrow(tbl) != 90L) stop("Expected 30 contrasts x 3 cutoffs = 90 rows.")
write.csv(tbl, file.path(OUT_DIR, "DEG_counts_by_LFC_threshold_long.csv"),
          row.names = FALSE)

# A convenient wide table with all three thresholds on one line per contrast.
wide <- reshape(tbl[, c("generation", "genotype", "treatment",
  "abs_log2FC_cutoff", "genes_with_testable_padj", "DEG_total", "DEG_up",
  "DEG_down")], idvar = c("generation", "genotype", "treatment"),
  timevar = "abs_log2FC_cutoff", direction = "wide")
wide <- wide[order(match(wide$generation, paste0("G", GENS)),
                   match(wide$genotype, GENOTYPES),
                   match(wide$treatment, TREATMENTS)), ]
rownames(wide) <- NULL
write.csv(wide, file.path(OUT_DIR, "DEG_counts_by_LFC_threshold_wide.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_sessionInfo.txt"))
message("Saved DEG threshold tables in: ", normalizePath(OUT_DIR))
