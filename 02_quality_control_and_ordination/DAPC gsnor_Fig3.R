rm(list = ls())

## ============================================================
## DAPC for gsnor-ko RNA-seq: combined Generation × Treatment classes
## Batch1 = gsnor G1 + G3
## Batch2 = gsnor G5
## ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(adegenet)
  library(genefilter)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(tidyr)
})

## ------------------------------------------------------------
## 0. Input
## ------------------------------------------------------------

count_file <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Count_matrix.txt"

out_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/WT_only_PCA_DAPC_exact_old_input"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

## ------------------------------------------------------------
## 1. Read count matrix
## ------------------------------------------------------------

df <- read.table(
  count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------
## 2. Define gsnor sample lists
## ------------------------------------------------------------

gsnor_gen1_samples <- c(
  "gs1C", "gs3C", "gs6C", "gs7C", "gs8C",
  "gs1D", "gs2D", "gs4D", "gs5D",
  "gs3CO", "gs4CO", "gs8CO", "gs9CO",
  "gs1O", "gs3O", "gs5O", "gs6O",
  "gs5H", "gs7H", "gs8H", "gs10H",
  "gs2T", "gs6T", "gs7T", "gs10T"
)

gsnor_gen3_samples <- c(
  "gs1C4_2", "gs3C4_1", "gs6C2_1", "gs6C2_3", "gs7C3_2",
  "gs1D3_1", "gs2D1_2", "gs4D3_2", "gs5D3_3",
  "gs3CO3_2", "gs4CO4_2", "gs8CO1_3", "gs9CO1_2",
  "gs1O2_2", "gs3O3_3", "gs5O4_2", "gs6O2_2",
  "gs5H3_2", "gs7H4_4", "gs8H1_2", "gs10H3_4",
  "gs2T4_2", "gs6T3_1", "gs7T2_4", "gs10T2_5"
)

gsnor_gen5_samples <- c(
  "gs1C4_212", "gs3C4_111", "gs6C2_311", "gs6C2_312", "gs7C3_211",
  "gs1D3_112", "gs2D1_211", "gs4D3_221", "gs5D3_322",
  "gs3CO3_212", "gs4CO4_222", "gs8CO1_322", "gs9CO1_212",
  "gs1O2_211", "gs3O3_322", "gs5O4_212", "gs6O2_221",
  "gs5H3_221", "gs7H4_422", "gs8H1_221", "gs10H3_412",
  "gs2T4_211", "gs6T3_111", "gs7T2_412", "gs10T2_522"
)

gsnor_all_samples <- c(
  gsnor_gen1_samples,
  gsnor_gen3_samples,
  gsnor_gen5_samples
)

missing_samples <- setdiff(gsnor_all_samples, colnames(df))

if (length(missing_samples) > 0) {
  stop(
    "These gsnor samples are missing from the count matrix:\n",
    paste(missing_samples, collapse = ", ")
  )
}

counts <- df[, gsnor_all_samples, drop = FALSE]

## ------------------------------------------------------------
## 3. Define condition factor
## ------------------------------------------------------------

condition <- factor(
  c(
    rep("gsnorgen1_Control", 5),
    rep("gsnorgen1_Drought", 4),
    rep("gsnorgen1_CO2", 4),
    rep("gsnorgen1_Ozone", 4),
    rep("gsnorgen1_Heat", 4),
    rep("gsnorgen1_Combination", 4),
    
    rep("gsnorgen3_Control", 5),
    rep("gsnorgen3_Drought", 4),
    rep("gsnorgen3_CO2", 4),
    rep("gsnorgen3_Ozone", 4),
    rep("gsnorgen3_Heat", 4),
    rep("gsnorgen3_Combination", 4),
    
    rep("gsnorgen5_Control", 5),
    rep("gsnorgen5_Drought", 4),
    rep("gsnorgen5_CO2", 4),
    rep("gsnorgen5_Ozone", 4),
    rep("gsnorgen5_Heat", 4),
    rep("gsnorgen5_Combination", 4)
  ),
  levels = c(
    "gsnorgen1_Control", "gsnorgen1_Drought", "gsnorgen1_CO2",
    "gsnorgen1_Ozone", "gsnorgen1_Heat", "gsnorgen1_Combination",
    "gsnorgen3_Control", "gsnorgen3_Drought", "gsnorgen3_CO2",
    "gsnorgen3_Ozone", "gsnorgen3_Heat", "gsnorgen3_Combination",
    "gsnorgen5_Control", "gsnorgen5_Drought", "gsnorgen5_CO2",
    "gsnorgen5_Ozone", "gsnorgen5_Heat", "gsnorgen5_Combination"
  )
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

## ------------------------------------------------------------
## 4. Assign generation, treatment, group, batch
## ------------------------------------------------------------

coldata$generation <- factor(
  ifelse(
    rownames(coldata) %in% gsnor_gen1_samples, "Gen1",
    ifelse(
      rownames(coldata) %in% gsnor_gen3_samples, "Gen3",
      "Gen5"
    )
  ),
  levels = c("Gen1", "Gen3", "Gen5")
)

coldata$group <- factor(
  ifelse(
    rownames(coldata) %in% gsnor_gen1_samples, "gsnor Gen 1",
    ifelse(
      rownames(coldata) %in% gsnor_gen3_samples, "gsnor Gen 3",
      "gsnor Gen5"
    )
  ),
  levels = c("gsnor Gen 1", "gsnor Gen 3", "gsnor Gen5")
)

## Batch design requested:
## Batch1 = G1 + G3
## Batch2 = G5
coldata$batch <- factor(
  ifelse(
    rownames(coldata) %in% c(gsnor_gen1_samples, gsnor_gen3_samples),
    "Batch1",
    "Batch2"
  ),
  levels = c("Batch1", "Batch2")
)

control_samples <- c(
  "gs1C", "gs3C", "gs6C", "gs7C", "gs8C",
  "gs1C4_2", "gs3C4_1", "gs6C2_1", "gs6C2_3", "gs7C3_2",
  "gs1C4_212", "gs3C4_111", "gs6C2_311", "gs6C2_312", "gs7C3_211"
)

drought_samples <- c(
  "gs1D", "gs2D", "gs4D", "gs5D",
  "gs1D3_1", "gs2D1_2", "gs4D3_2", "gs5D3_3",
  "gs1D3_112", "gs2D1_211", "gs4D3_221", "gs5D3_322"
)

CO2_samples <- c(
  "gs3CO", "gs4CO", "gs8CO", "gs9CO",
  "gs3CO3_2", "gs4CO4_2", "gs8CO1_3", "gs9CO1_2",
  "gs3CO3_212", "gs4CO4_222", "gs8CO1_322", "gs9CO1_212"
)

ozone_samples <- c(
  "gs1O", "gs3O", "gs5O", "gs6O",
  "gs1O2_2", "gs3O3_3", "gs5O4_2", "gs6O2_2",
  "gs1O2_211", "gs3O3_322", "gs5O4_212", "gs6O2_221"
)

heat_samples <- c(
  "gs5H", "gs7H", "gs8H", "gs10H",
  "gs5H3_2", "gs7H4_4", "gs8H1_2", "gs10H3_4",
  "gs5H3_221", "gs7H4_422", "gs8H1_221", "gs10H3_412"
)

combination_samples <- c(
  "gs2T", "gs6T", "gs7T", "gs10T",
  "gs2T4_2", "gs6T3_1", "gs7T2_4", "gs10T2_5",
  "gs2T4_211", "gs6T3_111", "gs7T2_412", "gs10T2_522"
)

coldata$treatment <- factor(
  ifelse(
    rownames(coldata) %in% control_samples, "Control",
    ifelse(
      rownames(coldata) %in% drought_samples, "Drought",
      ifelse(
        rownames(coldata) %in% CO2_samples, "CO2",
        ifelse(
          rownames(coldata) %in% ozone_samples, "Ozone",
          ifelse(
            rownames(coldata) %in% heat_samples, "Heat",
            "Combination"
          )
        )
      )
    )
  ),
  levels = c("CO2", "Combination", "Control", "Drought", "Heat", "Ozone")
)

cat("\ngsnor generation × treatment table:\n")
print(table(coldata$generation, coldata$treatment))

cat("\ngsnor batch table:\n")
print(table(coldata$batch, coldata$generation))

write.csv(
  coldata,
  "gsnor_coldata_DAPC.csv"
)

## ------------------------------------------------------------
## 5. DESeq2 object and filtering
## ------------------------------------------------------------

counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- estimateSizeFactors(dds)

idx <- rowSums(counts(dds, normalized = TRUE) >= 3) >= 2

dds <- dds[idx, ]

## ------------------------------------------------------------
## 6. VST transformation
## ------------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

## Batch correction ONLY for visualization/QC
## Do not use this for DESeq2 differential expression testing
vsd_bc <- vsd

if ("batch" %in% colnames(colData(vsd_bc))) {
  assay(vsd_bc) <- limma::removeBatchEffect(
    assay(vsd_bc),
    batch = colData(vsd_bc)$batch
  )
}

## Choose object for plotting
obj_plot <- vsd_bc

## For CV statistics, raw VST is safer
obj_cv <- vsd

saveRDS(obj_plot, "gsnor_vsd_batch_corrected_for_DAPC_plot.rds")
saveRDS(obj_cv, "gsnor_vsd_raw_for_DAPC_CV.rds")

## ------------------------------------------------------------
## 7. Add generation column
## ------------------------------------------------------------

add_generation <- function(vsd_obj) {
  
  cd <- as.data.frame(colData(vsd_obj))
  
  if (!"generation" %in% colnames(cd)) {
    
    if ("group" %in% colnames(cd)) {
      
      generation <- dplyr::recode(
        as.character(cd$group),
        "gsnor Gen 1" = "Gen1",
        "gsnor Gen 3" = "Gen3",
        "gsnor Gen5"  = "Gen5",
        .default      = as.character(cd$group)
      )
      
    } else if ("condition" %in% colnames(cd)) {
      
      generation_num <- stringr::str_match(as.character(cd$condition), "gen(\\d+)")[, 2]
      generation <- paste0("Gen", generation_num)
      
    } else {
      stop("No 'group', 'condition', or 'generation' column found in colData.")
    }
    
    colData(vsd_obj)$generation <- factor(
      generation,
      levels = c("Gen1", "Gen3", "Gen5")
    )
  }
  
  vsd_obj
}

obj_plot <- add_generation(obj_plot)
obj_cv   <- add_generation(obj_cv)

## ------------------------------------------------------------
## 8. Create combined Gen × Treatment class
## ------------------------------------------------------------

colData(obj_plot)$treat_gen <- interaction(
  colData(obj_plot)$generation,
  colData(obj_plot)$treatment,
  sep = "_",
  drop = TRUE
)

colData(obj_cv)$treat_gen <- interaction(
  colData(obj_cv)$generation,
  colData(obj_cv)$treatment,
  sep = "_",
  drop = TRUE
)

cat("\nCombined Generation × Treatment classes:\n")
print(table(colData(obj_plot)$treat_gen))

## ------------------------------------------------------------
## 9. Prepare expression matrix for DAPC
## ------------------------------------------------------------

get_X <- function(vsd_obj,
                  ntop = 2000,
                  scale_features = TRUE) {
  
  mat <- assay(vsd_obj)
  
  ## Select top variable genes
  rv <- genefilter::rowVars(mat)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  X <- t(mat[sel, , drop = FALSE])   # samples × genes
  
  if (scale_features) {
    X <- scale(X)
    X[is.na(X)] <- 0
  }
  
  attr(X, "selected_genes") <- rownames(mat)[sel]
  
  X
}

X_plot <- get_X(obj_plot, ntop = 2000, scale_features = TRUE)
X_cv   <- get_X(obj_cv,   ntop = 2000, scale_features = TRUE)

write.csv(
  data.frame(gene_id = attr(X_plot, "selected_genes")),
  "gsnor_DAPC_top_2000_genes_plot.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = attr(X_cv, "selected_genes")),
  "gsnor_DAPC_top_2000_genes_CV.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 10. Fit DAPC with cross-validation-selected PCs
## ------------------------------------------------------------

fit_dapc_xval <- function(X,
                          class,
                          seed = 1,
                          n_pca_max = 60,
                          n_rep = 50) {
  
  class <- droplevels(factor(class))
  
  ## Safe maximum number of PCs
  max_pcs <- min(
    n_pca_max,
    nrow(X) - nlevels(class) - 1,
    ncol(X) - 1
  )
  
  if (max_pcs < 2) {
    stop("Too few samples/classes to run DAPC safely.")
  }
  
  set.seed(seed)
  
  xval <- adegenet::xvalDapc(
    X,
    class,
    n.pca.max = max_pcs,
    training.set = 0.9,
    result = "groupMean",
    center = FALSE,
    scale = FALSE,
    n.rep = 50,
    xval.plot = FALSE
  )
  
  ## Robust extraction of best PC number
  best_n_pca <- suppressWarnings(
    as.numeric(xval$`Number of PCs Achieving Lowest MSE`)
  )
  
  if (is.na(best_n_pca)) {
    best_n_pca <- suppressWarnings(
      as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
    )
  }
  
  if (is.na(best_n_pca)) {
    rmse <- xval$`Root Mean Squared Error by Number of PCs of PCA`
    best_n_pca <- as.numeric(names(rmse)[which.min(rmse)])
  }
  
  n_da <- min(nlevels(class) - 1, 10)
  
  dapc_fit <- adegenet::dapc(
    X,
    class,
    n.pca = best_n_pca,
    n.da = n_da
  )
  
  list(
    dapc = dapc_fit,
    xval = xval,
    n_pca = best_n_pca,
    n_da = n_da,
    class = class
  )
}

res_comb <- fit_dapc_xval(
  X_plot,
  class = colData(obj_plot)$treat_gen,
  seed = 1,
  n_pca_max = 60,
  n_rep = 50
)

cat("\ngsnor selected n.pca:\n")
print(res_comb$n_pca)

cat("\ngsnor selected n.da:\n")
print(res_comb$n_da)

saveRDS(res_comb$xval, "gsnor_xvalDapc_combined_Generation_Treatment.rds")
saveRDS(res_comb$dapc, "gsnor_DAPC_combined_Generation_Treatment.rds")

## ------------------------------------------------------------
## 11. Convert DAPC output to ggplot data frame
## ------------------------------------------------------------

dapc_to_df <- function(dapc_fit, vsd_obj) {
  
  coords <- as.data.frame(dapc_fit$ind.coord)
  coords$sample <- rownames(coords)
  
  cd <- as.data.frame(colData(vsd_obj)) %>%
    tibble::rownames_to_column("sample")
  
  df <- coords %>%
    left_join(cd, by = "sample")
  
  df
}

df_comb <- dapc_to_df(res_comb$dapc, obj_plot)

write.csv(
  df_comb,
  "gsnor_DAPC_coordinates_combined_Generation_Treatment.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 12. Nice combined DAPC plot
## ------------------------------------------------------------

p_comb <- ggplot(df_comb, aes(x = LD1, y = LD2)) +
  geom_point(
    aes(color = treatment, shape = generation),
    size = 3,
    alpha = 0.9
  ) +
  stat_ellipse(
    aes(group = treat_gen, color = treatment),
    type = "norm",
    level = 0.68,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  facet_wrap(~ generation) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "DAPC of gsnor-ko RNA-seq samples",
    subtitle = paste0(
      "Classes = Generation × Treatment; n.pca = ",
      res_comb$n_pca,
      ", n.da = ",
      res_comb$n_da
    ),
    x = "LD1",
    y = "LD2",
    color = "Treatment",
    shape = "Generation"
  )

p_comb
##saved pdf as 5X7
ggsave(
  filename = "gsnor_DAPC_combined_Generation_Treatment.pdf",
  plot = p_comb,
  width = 11,
  height = 5.5
)

ggsave(
  filename = "gsnor_DAPC_combined_Generation_Treatment.png",
  plot = p_comb,
  width = 11,
  height = 5.5,
  dpi = 300
)

## ------------------------------------------------------------
## 13. Training-set assignment accuracy
##    This is optimistic, but useful as a first check
## ------------------------------------------------------------

pred_train <- predict(res_comb$dapc)

true_train <- colData(obj_plot)[rownames(pred_train$posterior), "treat_gen"]
true_train <- droplevels(factor(true_train))

train_confusion <- table(
  True = true_train,
  Pred = pred_train$assign
)

train_confusion

train_accuracy <- mean(pred_train$assign == true_train)

cat("\ngsnor training-set accuracy:\n")
print(train_accuracy)

write.csv(
  as.data.frame(train_confusion),
  "gsnor_DAPC_training_confusion_counts.csv",
  row.names = FALSE
)

write.csv(
  data.frame(training_accuracy = train_accuracy),
  "gsnor_DAPC_training_accuracy.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 14. Repeated CV accuracy
## ------------------------------------------------------------

dapc_cv_accuracy <- function(X,
                             class,
                             n_pca,
                             n_da,
                             train_prop = 0.8,
                             n_rep = 100,
                             seed = 1) {
  
  class <- droplevels(factor(class))
  lvls <- levels(class)
  
  acc <- numeric(n_rep)
  
  set.seed(seed)
  
  for (r in seq_len(n_rep)) {
    
    train_idx <- unlist(lapply(lvls, function(g) {
      idx <- which(class == g)
      n <- length(idx)
      
      ntr <- floor(n * train_prop)
      ntr <- max(ntr, 2)
      
      if ((n - ntr) < 1) {
        ntr <- n - 1
      }
      
      sample(idx, ntr)
    }))
    
    test_idx <- setdiff(seq_len(nrow(X)), train_idx)
    
    fit <- adegenet::dapc(
      X[train_idx, , drop = FALSE],
      class[train_idx],
      n.pca = min(n_pca, nrow(X[train_idx, , drop = FALSE]) - 1),
      n.da  = min(n_da, nlevels(class[train_idx]) - 1)
    )
    
    pred <- predict(
      fit,
      newdata = X[test_idx, , drop = FALSE]
    )$assign
    
    acc[r] <- mean(pred == class[test_idx])
  }
  
  data.frame(
    mean_accuracy = mean(acc),
    sd_accuracy = sd(acc),
    min_accuracy = min(acc),
    max_accuracy = max(acc)
  )
}

cv_summary <- dapc_cv_accuracy(
  X_cv,
  class = colData(obj_cv)$treat_gen,
  n_pca = res_comb$n_pca,
  n_da  = res_comb$n_da,
  train_prop = 0.8,
  n_rep = 100,
  seed = 1
)

cat("\ngsnor repeated CV summary:\n")
print(cv_summary)

write.csv(
  cv_summary,
  "gsnor_DAPC_CV_accuracy_summary.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 15. CV confusion matrix and top misclassifications
## ------------------------------------------------------------

cv_dapc_confusion <- function(X,
                              class,
                              n_pca,
                              n_da,
                              train_prop = 0.8,
                              n_rep = 200,
                              seed = 1) {
  
  class <- droplevels(factor(class))
  lvls <- levels(class)
  
  true_all <- character(0)
  pred_all <- character(0)
  
  set.seed(seed)
  
  for (r in seq_len(n_rep)) {
    
    train_idx <- unlist(lapply(lvls, function(g) {
      idx <- which(class == g)
      n <- length(idx)
      
      ntr <- floor(n * train_prop)
      ntr <- max(ntr, 2)
      
      if ((n - ntr) < 1) {
        ntr <- n - 1
      }
      
      sample(idx, ntr)
    }))
    
    test_idx <- setdiff(seq_len(nrow(X)), train_idx)
    
    fit <- tryCatch(
      adegenet::dapc(
        X[train_idx, , drop = FALSE],
        class[train_idx],
        n.pca = min(n_pca, nrow(X[train_idx, , drop = FALSE]) - 1),
        n.da  = min(n_da, nlevels(class[train_idx]) - 1)
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    pred <- tryCatch(
      predict(fit, newdata = X[test_idx, , drop = FALSE])$assign,
      error = function(e) NULL
    )
    
    if (is.null(pred)) next
    
    true_all <- c(true_all, as.character(class[test_idx]))
    pred_all <- c(pred_all, as.character(pred))
  }
  
  cm <- table(
    True = factor(true_all, levels = lvls),
    Pred = factor(pred_all, levels = lvls)
  )
  
  cm_pct <- prop.table(cm, margin = 1) * 100
  
  list(
    cm = cm,
    cm_pct = cm_pct,
    true_all = true_all,
    pred_all = pred_all
  )
}

cv_out <- cv_dapc_confusion(
  X_cv,
  class = colData(obj_cv)$treat_gen,
  n_pca = res_comb$n_pca,
  n_da  = res_comb$n_da,
  train_prop = 0.8,
  n_rep = 200,
  seed = 1
)

cat("\ngsnor CV confusion matrix counts:\n")
print(cv_out$cm)

cat("\ngsnor CV confusion matrix percentages:\n")
print(round(cv_out$cm_pct, 1))

## Per-class CV accuracy
per_class_acc <- data.frame(
  Class = rownames(cv_out$cm),
  n_test = rowSums(cv_out$cm),
  accuracy = diag(cv_out$cm) / pmax(1, rowSums(cv_out$cm))
) %>%
  arrange(accuracy)

per_class_acc

## Top misclassification pairs
top_pairs <- as.data.frame(cv_out$cm) %>%
  rename(n = Freq) %>%
  filter(True != Pred, n > 0) %>%
  group_by(True) %>%
  mutate(
    total_wrong_for_true_class = sum(n),
    mis_rate_within_true = n / total_wrong_for_true_class
  ) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  mutate(pair = paste0(True, " \u2192 ", Pred)) %>%
  select(pair, True, Pred, n, mis_rate_within_true)

head(top_pairs, 20)

## ------------------------------------------------------------
## 16. CV confusion heatmap
## ------------------------------------------------------------

cm_pct_df <- as.data.frame(cv_out$cm_pct) %>%
  rename(pct = Freq)

p_conf <- ggplot(cm_pct_df, aes(x = Pred, y = True, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.2) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "gsnor-ko DAPC CV confusion matrix",
    subtitle = "Row-wise percentage of held-out predictions",
    x = "Predicted class",
    y = "True class",
    fill = "%"
  )

p_conf

ggsave(
  filename = "gsnor_DAPC_CV_confusion_matrix.pdf",
  plot = p_conf,
  width = 8,
  height = 7
)

ggsave(
  filename = "gsnor_DAPC_CV_confusion_matrix.png",
  plot = p_conf,
  width = 8,
  height = 7,
  dpi = 300
)

## ------------------------------------------------------------
## 17. Top genes contributing to LD1 / LD2
## ------------------------------------------------------------

get_top_dapc_genes <- function(dapc_fit,
                               X,
                               axis = "LD1",
                               n = 30) {
  
  vc <- dapc_fit$var.contr
  
  if (is.null(dim(vc))) {
    score <- abs(vc)
    names(score) <- colnames(X)
  } else {
    
    if (is.null(rownames(vc))) {
      rownames(vc) <- colnames(X)
    }
    
    if (!axis %in% colnames(vc)) {
      stop(
        "Axis not found. Available axes are: ",
        paste(colnames(vc), collapse = ", ")
      )
    }
    
    score <- abs(vc[, axis])
  }
  
  data.frame(
    gene_id = names(score),
    contribution = as.numeric(score)
  ) %>%
    arrange(desc(contribution)) %>%
    slice_head(n = n)
}

top_LD1_genes <- get_top_dapc_genes(
  dapc_fit = res_comb$dapc,
  X = X_plot,
  axis = "LD1",
  n = 50
)

top_LD2_genes <- get_top_dapc_genes(
  dapc_fit = res_comb$dapc,
  X = X_plot,
  axis = "LD2",
  n = 50
)

top_LD1_genes
top_LD2_genes

write.csv(top_LD1_genes, "gsnor_DAPC_top_LD1_genes.csv", row.names = FALSE)
write.csv(top_LD2_genes, "gsnor_DAPC_top_LD2_genes.csv", row.names = FALSE)

## ------------------------------------------------------------
## 18. Save important CV tables
## ------------------------------------------------------------

write.csv(
  as.data.frame(cv_out$cm),
  "gsnor_DAPC_CV_confusion_counts.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(round(cv_out$cm_pct, 2)),
  "gsnor_DAPC_CV_confusion_percent.csv",
  row.names = FALSE
)

write.csv(
  per_class_acc,
  "gsnor_DAPC_per_class_CV_accuracy.csv",
  row.names = FALSE
)

write.csv(
  top_pairs,
  "gsnor_DAPC_top_misclassification_pairs.csv",
  row.names = FALSE
)

cat("\nDone. gsnor-ko DAPC results saved in:\n")
cat(out_dir, "\n")
