rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(genefilter)
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(adegenet)
  library(stringr)
  library(tibble)
})

register(MulticoreParam(8))

# ============================================================
# 1. Input
# ============================================================

count_file <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Count_matrix.txt"

out_dir <- "/Users/rachnabehl/Documents/DNA_seq/RNA_seq/WT_only_PCA_DAPC_exact_old_input"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ntop <- 2000
seed_use <- 1
n_rep_xval <- 50
n_rep_cv <- 50

# Optional: only use this if you know the old WT n.pca value
force_n_pca_WT <- NULL

# Samples to exclude if present
bad_samples <- c(
  "WT_Pf_2-1-4",
  "WT7O_new",
  "WT8C_new"
)

# ============================================================
# 2. Read count matrix
# ============================================================

df <- read.table(
  count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

df <- df[, !colnames(df) %in% bad_samples, drop = FALSE]

cat("\nExcluded samples still present?\n")
print(intersect(bad_samples, colnames(df)))

# ============================================================
# 3. Exact WT sample order from old PCA script
# ============================================================

WT_samples <- c(
  # Gen1
  "WT1C", "WT3C", "WT5C", "WT6C", "WT8C", "WT10C",
  "WT3D", "WT4D", "WT6D", "WT10D",
  "WT3CO", "WT4CO", "WT5CO", "WT7CO",
  "WT1O", "WT2O", "WT5O", "WT7O",
  "WT1H", "WT5H", "WT6H", "WT8H",
  "WT1T", "WT2T", "WT3T", "WT6T",
  
  # Gen3
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4",
  
  # Gen5
  "WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111",
  "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211",
  "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212",
  "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112",
  "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221",
  "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422"
)

missing_samples <- setdiff(WT_samples, colnames(df))

if (length(missing_samples) > 0) {
  stop(
    "These WT samples are missing from the count matrix:\n",
    paste(missing_samples, collapse = ", ")
  )
}

counts <- dplyr::select(df, all_of(WT_samples))

# ============================================================
# 4. Condition exactly as old WT PCA script
# ============================================================

condition <- factor(
  c(
    rep("WTgen1_Control", 6),
    rep("WTgen1_Drought", 4),
    rep("WTgen1_CO2", 4),
    rep("WTgen1_Ozone", 4),
    rep("WTgen1_Heat", 4),
    rep("WTgen1_Combination", 4),
    
    rep("WTgen3_Control", 5),
    rep("WTgen3_Drought", 4),
    rep("WTgen3_CO2", 4),
    rep("WTgen3_Ozone", 4),
    rep("WTgen3_Heat", 4),
    rep("WTgen3_Combination", 4),
    
    rep("WTgen5_Control", 5),
    rep("WTgen5_Drought", 4),
    rep("WTgen5_CO2", 4),
    rep("WTgen5_Ozone", 4),
    rep("WTgen5_Heat", 4),
    rep("WTgen5_Combination", 4)
  ),
  levels = c(
    "WTgen1_Control", "WTgen1_Drought", "WTgen1_CO2",
    "WTgen1_Ozone", "WTgen1_Heat", "WTgen1_Combination",
    "WTgen3_Control", "WTgen3_Drought", "WTgen3_CO2",
    "WTgen3_Ozone", "WTgen3_Heat", "WTgen3_Combination",
    "WTgen5_Control", "WTgen5_Drought", "WTgen5_CO2",
    "WTgen5_Ozone", "WTgen5_Heat", "WTgen5_Combination"
  )
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# ============================================================
# 5. Exact four WT batches
# ============================================================

batch1_samples <- c(
  "WT1C", "WT6C", "WT8C",
  "WT3D", "WT4D",
  "WT3CO", "WT7CO",
  "WT5O", "WT7O",
  "WT1H", "WT5H",
  "WT2T", "WT3T"
)

batch2_samples <- c(
  "WT3C", "WT5C", "WT10C",
  "WT6D", "WT10D",
  "WT4CO", "WT5CO",
  "WT1O", "WT2O",
  "WT6H", "WT8H",
  "WT1T", "WT6T"
)

batch3_samples <- c(
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4"
)

batch4_samples <- c(
  "WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111",
  "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211",
  "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212",
  "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112",
  "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221",
  "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422"
)

coldata$batch <- factor(
  ifelse(
    rownames(coldata) %in% batch1_samples, "Batch1",
    ifelse(
      rownames(coldata) %in% batch2_samples, "Batch2",
      ifelse(
        rownames(coldata) %in% batch3_samples, "Batch3",
        ifelse(
          rownames(coldata) %in% batch4_samples, "Batch4",
          NA
        )
      )
    )
  ),
  levels = c("Batch1", "Batch2", "Batch3", "Batch4")
)

if (any(is.na(coldata$batch))) {
  stop(
    "These samples were not assigned to a batch:\n",
    paste(rownames(coldata)[is.na(coldata$batch)], collapse = ", ")
  )
}

# ============================================================
# 6. Group / generation / treatment exactly for plotting
# ============================================================

WT_gen1_samples <- c(
  "WT1C", "WT3C", "WT5C", "WT6C", "WT8C", "WT10C",
  "WT3D", "WT4D", "WT6D", "WT10D",
  "WT3CO", "WT4CO", "WT5CO", "WT7CO",
  "WT1O", "WT2O", "WT5O", "WT7O",
  "WT1H", "WT5H", "WT6H", "WT8H",
  "WT1T", "WT2T", "WT3T", "WT6T"
)

WT_gen3_samples <- batch3_samples
WT_gen5_samples <- batch4_samples

coldata$group <- factor(
  ifelse(
    rownames(coldata) %in% WT_gen1_samples, "WT Gen 1",
    ifelse(
      rownames(coldata) %in% WT_gen3_samples, "WT Gen 3",
      "WT Gen5"
    )
  ),
  levels = c("WT Gen 1", "WT Gen 3", "WT Gen5")
)

coldata$generation <- factor(
  ifelse(
    rownames(coldata) %in% WT_gen1_samples, "Gen1",
    ifelse(
      rownames(coldata) %in% WT_gen3_samples, "Gen3",
      "Gen5"
    )
  ),
  levels = c("Gen1", "Gen3", "Gen5")
)

control_samples <- c(
  "WT1C", "WT3C", "WT5C", "WT6C", "WT8C", "WT10C",
  "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3",
  "WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111"
)

drought_samples <- c(
  "WT3D", "WT4D", "WT6D", "WT10D",
  "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2",
  "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211"
)

CO2_samples <- c(
  "WT3CO", "WT4CO", "WT5CO", "WT7CO",
  "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
  "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212"
)

ozone_samples <- c(
  "WT1O", "WT2O", "WT5O", "WT7O",
  "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1",
  "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112"
)

heat_samples <- c(
  "WT1H", "WT5H", "WT6H", "WT8H",
  "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2",
  "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221"
)

combination_samples <- c(
  "WT1T", "WT2T", "WT3T", "WT6T",
  "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4",
  "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422"
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
  )
)

cat("\nWT batch table:\n")
print(table(coldata$batch, coldata$group))

cat("\nWT treatment table:\n")
print(table(coldata$generation, coldata$treatment))

write.csv(
  coldata,
  file.path(out_dir, "WT_coldata_exact_old_input.csv")
)

# ============================================================
# 7. DESeq2 VST
# ============================================================

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



## ============================================================
## DAPC for WT RNA-seq: combined Generation × Treatment classes
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
## 1. VST transformation
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


## ------------------------------------------------------------
## 2. Add generation column
## ------------------------------------------------------------

add_generation <- function(vsd_obj) {
  
  cd <- as.data.frame(colData(vsd_obj))
  
  if (!"generation" %in% colnames(cd)) {
    
    if ("group" %in% colnames(cd)) {
      
      generation <- dplyr::recode(
        as.character(cd$group),
        "WT Gen 1" = "Gen1",
        "WT Gen 3" = "Gen3",
        "WT Gen 5" = "Gen5",
        "WT Gen5"  = "Gen5",
        .default   = as.character(cd$group)
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
## 3. Create combined Gen × Treatment class
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

table(colData(obj_plot)$treat_gen)


## ------------------------------------------------------------
## 4. Prepare expression matrix for DAPC
## ------------------------------------------------------------

get_X <- function(vsd_obj,
                  ntop = 2000,
                  scale_features = TRUE) {
  
  mat <- assay(vsd_obj)
  
  ## select top variable genes
  rv <- genefilter::rowVars(mat)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  X <- t(mat[sel, , drop = FALSE])   # samples × genes
  
  if (scale_features) {
    X <- scale(X)
    X[is.na(X)] <- 0
  }
  
  X
}

X_plot <- get_X(obj_plot, ntop = 2000, scale_features = TRUE)
X_cv   <- get_X(obj_cv,   ntop = 2000, scale_features = TRUE)


## ------------------------------------------------------------
## 5. Fit DAPC with cross-validation-selected PCs
## ------------------------------------------------------------

fit_dapc_xval <- function(X,
                          class,
                          seed = 1,
                          n_pca_max = 60,
                          n_rep = 50) {
  
  class <- droplevels(factor(class))
  
  ## safe maximum number of PCs
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
    n.pca.max  = max_pcs,
    training.set = 0.9,
    result = "groupMean",
    center = FALSE,
    scale = FALSE,
    n.rep = n_rep,
    xval.plot = FALSE
  )
  
  ## robust extraction of best PC number
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
    n.da  = n_da
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

res_comb$n_pca
res_comb$n_da


## ------------------------------------------------------------
## 6. Convert DAPC output to ggplot data frame
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

head(df_comb)


## ------------------------------------------------------------
## 7. Nice combined DAPC plot
## ------------------------------------------------------------

p_comb <- ggplot(df_comb, aes(x = LD1, y = LD2)) +
  geom_point(
    aes(color = treatment, shape = generation),
    size = 4,
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
    title = "DAPC of WT RNA-seq samples",
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

## ------------------------------------------------------------
## 8. Save DAPC plot
## ------------------------------------------------------------

ggsave(
  filename = "WT_DAPC_combined_Generation_Treatment.pdf",
  plot = p_comb,
  width = 11,
  height = 5.5
)

ggsave(
  filename = "WT_DAPC_combined_Generation_Treatment.png",
  plot = p_comb,
  width = 11,
  height = 5.5,
  dpi = 300
)


## ------------------------------------------------------------
## 9. Training-set assignment accuracy
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
train_accuracy


## ------------------------------------------------------------
## 10. Repeated CV accuracy
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

cv_summary


## ------------------------------------------------------------
## 11. CV confusion matrix and top misclassifications
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

## Confusion matrix counts
cv_out$cm

## Confusion matrix percentages
round(cv_out$cm_pct, 1)


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
## 12. CV confusion heatmap
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
    title = "DAPC CV confusion matrix",
    subtitle = "Row-wise percentage of held-out predictions",
    x = "Predicted class",
    y = "True class",
    fill = "%"
  )

p_conf

ggsave(
  filename = "WT_DAPC_CV_confusion_matrix.pdf",
  plot = p_conf,
  width = 8,
  height = 7
)

ggsave(
  filename = "WT_DAPC_CV_confusion_matrix.png",
  plot = p_conf,
  width = 8,
  height = 7,
  dpi = 300
)


## ------------------------------------------------------------
## 13. Top genes contributing to LD1 / LD2
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
      stop("Axis not found. Available axes are: ", paste(colnames(vc), collapse = ", "))
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

write.csv(top_LD1_genes, "WT_DAPC_top_LD1_genes.csv", row.names = FALSE)
write.csv(top_LD2_genes, "WT_DAPC_top_LD2_genes.csv", row.names = FALSE)


## ------------------------------------------------------------
## 14. Save important tables
## ------------------------------------------------------------

write.csv(as.data.frame(cv_out$cm),
          "WT_DAPC_CV_confusion_counts.csv",
          row.names = FALSE)

write.csv(as.data.frame(round(cv_out$cm_pct, 2)),
          "WT_DAPC_CV_confusion_percent.csv",
          row.names = FALSE)

write.csv(per_class_acc,
          "WT_DAPC_per_class_CV_accuracy.csv",
          row.names = FALSE)

write.csv(top_pairs,
          "WT_DAPC_top_misclassification_pairs.csv",
          row.names = FALSE)

write.csv(cv_summary,
          "WT_DAPC_CV_accuracy_summary.csv",
          row.names = FALSE)

