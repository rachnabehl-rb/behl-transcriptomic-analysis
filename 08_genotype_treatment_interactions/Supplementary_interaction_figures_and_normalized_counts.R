# Supplementary figures for the genotype x treatment RNA-seq analysis.
# Run from the directory containing the 'upload' folder, or change INPUT_DIR.
# Requires DESeq2 and ggplot2. Does not modify the manuscript or source files.

INPUT_DIR <- "~/Documents/DNA_seq/RNA_seq"
OUTPUT_DIR <- "supplementary_interaction_figures"
COUNT_FILE <- file.path(INPUT_DIR, "genes_gtf_merged_gene_counts.txt")
SAMPLE_FILE <- file.path(INPUT_DIR, "Metadata_Sample_Area_medians_WITH_GROUPS.csv")

# Edit this table to show other example genes. Use exact TAIR gene IDs.
examples <- data.frame(
  treatment = c("CO2", "Combination"),
  gene_id = c("AT5G65630", "AT4G09020"),
  stringsAsFactors = FALSE
)

SCATTER_LIMIT <- 6  # visual zoom in log2FC units; does not affect statistical tests
generations <- c(1L, 3L, 5L)
treatments <- c("Drought", "CO2", "Ozone", "Heat", "Combination")

for (pkg in c("ggplot2", "DESeq2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install the R package: ", pkg)
  }
}
for (f in c(COUNT_FILE, SAMPLE_FILE)) {
  if (!file.exists(f)) stop("Input file missing: ", f)
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# Figure S1: interaction count overview (G1, G3, G5 x five treatments).
overview <- list()
g5_tables <- list()
for (g in generations) {
  for (tr in treatments) {
    file <- file.path(INPUT_DIR,
      paste0("G", g, "_", tr, "_three_effects_for_gene_review.csv"))
    if (!file.exists(file)) stop("Missing interaction results: ", file)
    d <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    needed <- c("gene_id", "WT_log2FC", "WT_padj", "gsnor_log2FC",
                "gsnor_padj", "interaction_log2FC", "interaction_padj")
    if (!all(needed %in% names(d))) {
      stop("Missing columns in ", file, ": ",
           paste(setdiff(needed, names(d)), collapse = ", "))
    }
    if (anyDuplicated(d$gene_id)) stop("Duplicate gene IDs: ", file)
    d$generation <- paste0("G", g)
    d$treatment <- tr
    d$interaction_significant <- !is.na(d$interaction_padj) &
      d$interaction_padj < 0.05
    overview[[length(overview) + 1L]] <- data.frame(
      generation = paste0("G", g), treatment = tr,
      significant_interactions = sum(d$interaction_significant),
      testable_interactions = sum(!is.na(d$interaction_padj)),
      pct_of_testable = 100 * mean(d$interaction_significant[
        !is.na(d$interaction_padj)]),
      stringsAsFactors = FALSE
    )
    if (g == 5L && tr %in% c("CO2", "Combination", "Ozone")) {
      g5_tables[[tr]] <- d
    }
  }
}
overview <- do.call(rbind, overview)
overview$treatment <- factor(overview$treatment, levels = treatments)
overview$generation <- factor(overview$generation, levels = paste0("G", generations))
write.csv(overview, file.path(OUTPUT_DIR,
  "S1_interaction_count_overview_and_testable_genes.csv"), row.names = FALSE)

p1 <- ggplot2::ggplot(overview, ggplot2::aes(treatment, generation,
                                             fill = significant_interactions)) +
  ggplot2::geom_tile(color = "white", linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = format(significant_interactions,
    big.mark = ",", scientific = FALSE)), size = 4.5) +
  ggplot2::scale_fill_gradient(low = "#edf4f9", high = "#196b9a",
                               trans = "log1p", name = "Genes") +
  ggplot2::scale_y_discrete(limits = rev(paste0("G", generations))) +
  ggplot2::labs(x = "Treatment/history", y = "Generation",
    title = "Genes with genotype × treatment interactions (adjusted P < 0.05)",
    caption = "Counts are descriptive; testable-gene denominators differ by contrast (see companion CSV).") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
ggplot2::ggsave(file.path(OUTPUT_DIR, "SuppFig_S1_interaction_overview.pdf"),
                p1, width = 10.5, height = 5.2)

# Figure S2: G5 WT versus gsnor treatment-history effects.
scatter <- do.call(rbind, g5_tables)
scatter <- scatter[is.finite(scatter$WT_log2FC) &
                     is.finite(scatter$gsnor_log2FC), , drop = FALSE]
scatter$treatment <- factor(scatter$treatment,
                             levels = c("CO2", "Combination", "Ozone"))
sig <- scatter[scatter$interaction_significant, , drop = FALSE]

p2 <- ggplot2::ggplot(scatter,
                      ggplot2::aes(WT_log2FC, gsnor_log2FC)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, color = "grey82", linewidth = 0.3) +
  ggplot2::geom_abline(slope = 1, intercept = 0,
                       linetype = "dashed", color = "grey45") +
  ggplot2::geom_point(color = "grey55", alpha = 0.12, size = 0.45) +
  ggplot2::geom_point(data = sig, color = "#d94727",
                      alpha = 0.7, size = 0.7) +
  ggplot2::facet_wrap(~ treatment, nrow = 1) +
  ggplot2::coord_cartesian(xlim = c(-SCATTER_LIMIT, SCATTER_LIMIT),
                           ylim = c(-SCATTER_LIMIT, SCATTER_LIMIT)) +
  ggplot2::labs(x = "WT: treatment/history vs control (log2FC)",
    y = "gsnor: treatment/history vs control (log2FC)",
    title = "G5 genotype-dependent treatment-history effects",
    subtitle = "Orange: interaction adjusted P < 0.05; dashed line: equal effects",
    caption = paste("Axes zoomed to ±", SCATTER_LIMIT,
                    "log2FC; points outside the view are omitted from display only.")) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(OUTPUT_DIR, "SuppFig_S2_G5_WT_vs_gsnor_effects.pdf"),
                p2, width = 12, height = 5)

# Figure S3 and audit CSV: normalized counts for the example genes.
# Re-create the same G5 size-factor normalization and filtering as the
# interaction analysis. A Wald re-fit is unnecessary to recover these counts.
meta <- read.csv(SAMPLE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
required_meta <- c("SampleName", "Genotype", "Treatment", "Generation")
if (!all(required_meta %in% names(meta))) stop("Metadata columns missing.")
meta$Generation <- as.integer(sub("^[Gg]", "", meta$Generation))
m <- meta[meta$Generation == 5L, , drop = FALSE]
if (anyDuplicated(m$SampleName)) stop("Duplicate G5 sample names.")
m$genotype <- factor(m$Genotype, levels = c("WT", "gsnor"))
m$treatment <- factor(m$Treatment,
  levels = c("Control", treatments))
rownames(m) <- m$SampleName

raw <- as.matrix(read.table(COUNT_FILE, header = TRUE, row.names = 1,
                            sep = "\t", quote = "", comment.char = "",
                            check.names = FALSE))
if (!all(m$SampleName %in% colnames(raw))) {
  stop("Some G5 metadata sample names are absent from the raw count matrix.")
}
raw <- raw[, m$SampleName, drop = FALSE]
if (!is.numeric(raw) || anyNA(raw) || any(raw < 0) ||
    any(abs(raw - round(raw)) > 1e-8)) {
  stop("Input must be raw, nonnegative integer gene counts.")
}
storage.mode(raw) <- "integer"
dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw,
  colData = m, design = ~ genotype * treatment)
dds <- dds[rowSums(DESeq2::counts(dds)) > 10, ]
dds <- DESeq2::estimateSizeFactors(dds)
dds <- dds[rowSums(DESeq2::counts(dds, normalized = TRUE) >= 3) >= 2, ]
norm <- DESeq2::counts(dds, normalized = TRUE)

if (!all(examples$treatment %in% treatments)) {
  stop("An example treatment is absent from the treatment list.")
}
if (!all(examples$gene_id %in% rownames(norm))) {
  stop("These example genes were not retained: ",
       paste(setdiff(examples$gene_id, rownames(norm)), collapse = ", "))
}
if (anyDuplicated(paste(examples$treatment, examples$gene_id))) {
  stop("Repeated gene/treatment combinations in examples.")
}

long <- vector("list", nrow(examples))
for (i in seq_len(nrow(examples))) {
  tr <- examples$treatment[i]
  gene <- examples$gene_id[i]
  sel <- m$Treatment %in% c("Control", tr)
  sample_ids <- m$SampleName[sel]
  long[[i]] <- data.frame(
    gene_id = gene,
    history_treatment = tr,
    sample_id = sample_ids,
    genotype = m$Genotype[sel],
    observed_condition = m$Treatment[sel],
    raw_count = as.numeric(DESeq2::counts(dds)[gene, sample_ids]),
    size_factor = as.numeric(DESeq2::sizeFactors(dds)[sample_ids]),
    normalized_count = as.numeric(norm[gene, sample_ids]),
    stringsAsFactors = FALSE
  )
}
long <- do.call(rbind, long)
long$log2_normalized_count_plus1 <- log2(long$normalized_count + 1)
# Keep the four groups in biological order; controls repeat across panels.
long$group <- paste(long$genotype, long$observed_condition)
long$group <- factor(long$group, levels = unique(unlist(lapply(
  unique(examples$treatment), function(tr)
    paste(rep(c("WT", "gsnor"), each = 2), rep(c("Control", tr), 2))))))
long$panel <- paste0(long$gene_id, " (", long$history_treatment, ")")

write.csv(long[, c("gene_id", "history_treatment", "sample_id", "genotype",
  "observed_condition", "raw_count", "size_factor", "normalized_count",
  "log2_normalized_count_plus1")],
  file.path(OUTPUT_DIR, "G5_example_genes_normalized_counts_by_sample.csv"),
  row.names = FALSE)

p3 <- ggplot2::ggplot(long,
                      ggplot2::aes(group, log2_normalized_count_plus1,
                                   color = genotype)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12,
                                                          height = 0, seed = 21),
                      size = 2.4, alpha = 0.9) +
  ggplot2::stat_summary(fun = median, geom = "point", shape = 95,
                        size = 9, color = "black") +
  ggplot2::facet_wrap(~ panel, scales = "free", nrow = 1) +
  ggplot2::scale_color_manual(values = c(WT = "#3273a6", gsnor = "#c65636")) +
  ggplot2::labs(x = "Genotype and G5 history", y = "log2(normalized count + 1)",
    title = "Individual G5 libraries for selected interaction genes",
    subtitle = "Each dot is one RNA-seq library; black horizontal marker is the median") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(OUTPUT_DIR, "SuppFig_S3_G5_example_gene_counts.pdf"),
                p3, width = 11.5, height = 5)

writeLines(capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "R_sessionInfo.txt"))
message("Finished. Files are in: ", normalizePath(OUTPUT_DIR))
