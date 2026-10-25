###preprocessing to generate vsd file
rm(list = ls())
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))
library("genefilter")
library(dplyr)
library(ggplot2)

# Here, I assess leaf responses to the enviromental stresses 
# Input

#with batch effect
df <- read.table('/Users/rachnabehl/Documents/DNA_seq/RNA_seq/Count_matrix.txt', header=TRUE, row.names=1)

### "C"  includes 2 leaf control samples
### "S" includes 2 leaf samples affected by stresses
counts <- dplyr::select(df, WT1C, WT3C, WT5C, WT6C, WT8C, WT10C, WT3D, WT4D, WT6D, WT10D, WT3CO, WT4CO, WT5CO, WT7CO, 
                 WT1O, WT2O, WT5O, WT7O, WT1H, WT5H, WT6H, WT8H, WT1T, WT2T, WT3T, WT6T, 
                 WT1C3_1, WT4C3_2, WT6C3_2, WT8C2_1, WT8C2_3, WT3D2_2, WT4D4_1, WT6D4_1, WT7D3_2, WT3CO4_1, WT4CO3_2, WT5CO1_3, WT7CO1_2,
                 WT1O4_1, WT2O3_2, WT5O3_2, WT7O1_1, WT1H2_2, WT5H3_3, WT6H2_1, WT8H2_2, WT1T4_1, WT2T1_2, WT3T3_4, WT6T2_4, 
                 WT1C3_121, WT4C3_211, WT4C3_212, WT6C3_212, WT8C2_111, WT3D2_211, WT4D4_111, WT6D4_122, WT7D3_211, WT3CO4_112, WT4CO3_211, WT5CO1_321, WT7CO1_212,
                 WT1O4_111, WT2O3_212, WT5O3_221, WT7O1_112, WT1H2_211, WT5H3_312, WT6H2_122, WT8H2_221, WT1T4_111, WT2T1_212, WT3T3_411, WT6T2_422,
                 gs1C, gs3C, gs6C, gs7C, gs8C, gs1D, gs2D, gs4D, gs5D, gs3CO, gs4CO, gs8CO, gs9CO, 
                 gs1O, gs3O, gs5O, gs6O, gs5H, gs7H, gs8H, gs10H, gs2T, gs6T, gs7T, gs10T, 
                 gs1C4_2, gs3C4_1, gs6C2_1, gs6C2_3, gs7C3_2, gs1D3_1, gs2D1_2, gs4D3_2, gs5D3_3, gs3CO3_2, gs4CO4_2, gs8CO1_3, gs9CO1_2,
                 gs1O2_2, gs3O3_3, gs5O4_2, gs6O2_2, gs5H3_2, gs7H4_4, gs8H1_2, gs10H3_4, gs2T4_2, gs6T3_1, gs7T2_4, gs10T2_5,
                 gs1C4_212, gs3C4_111, gs6C2_311, gs6C2_312, gs7C3_211, gs1D3_112, gs2D1_211, gs4D3_221, gs5D3_322, gs3CO3_212, gs4CO4_222, gs8CO1_322, gs9CO1_212, 
                 gs1O2_211, gs3O3_322, gs5O4_212, gs6O2_221, gs5H3_221, gs7H4_422, gs8H1_221, gs10H3_412, gs2T4_211, gs6T3_111, gs7T2_412, gs10T2_522) # , gs2T4_2


# Update the condition factor to include all stress conditions
condition <- factor(c(rep("WTgen1_Control", 6), 
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
                      rep("WTgen5_Combination", 4),
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
                      rep("gsnorgen5_Combination", 4)), 
                    level = c("WTgen1_Control", "WTgen1_Drought", "WTgen1_CO2", "WTgen1_Ozone", "WTgen1_Heat", "WTgen1_Combination", 
                              "WTgen3_Control", "WTgen3_Drought", "WTgen3_CO2", "WTgen3_Ozone", "WTgen3_Heat", "WTgen3_Combination",
                              "WTgen5_Control", "WTgen5_Drought", "WTgen5_CO2", "WTgen5_Ozone", "WTgen5_Heat", "WTgen5_Combination",
                              "gsnorgen1_Control", "gsnorgen1_Drought", "gsnorgen1_CO2", "gsnorgen1_Ozone", "gsnorgen1_Heat", "gsnorgen1_Combination",
                              "gsnorgen3_Control", "gsnorgen3_Drought", "gsnorgen3_CO2", "gsnorgen3_Ozone", "gsnorgen3_Heat", "gsnorgen3_Combination",
                              "gsnorgen5_Control", "gsnorgen5_Drought", "gsnorgen5_CO2", "gsnorgen5_Ozone", "gsnorgen5_Heat", "gsnorgen5_Combination"))

# Update the coldata with the new condition
coldata <- data.frame(row.names = colnames(counts), condition)

# Define sample names for Batch 1 and Batch 2 (only for WT gen1)
batch1_samples <- c("WT1C", "WT6C", "WT8C", "WT3D", "WT4D", "WT3CO", "WT7CO", "WT5O", "WT7O", "WT1H", "WT5H", "WT2T", "WT3T")
batch2_samples <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D", "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")
batch3_samples <- c("gs1C", "gs3C", "gs6C", "gs7C", "gs8C", "gs1D", "gs2D", "gs4D", "gs5D", "gs3CO", "gs4CO", "gs8CO", "gs9CO", 
                    "gs1O", "gs3O", "gs5O", "gs6O", "gs5H", "gs7H", "gs8H", "gs10H", "gs2T", "gs6T", "gs7T", "gs10T",
                    "WT1C3_1", "WT4C3_2", "WT6C3_2", "WT8C2_1", "WT8C2_3", "WT3D2_2", "WT4D4_1", "WT6D4_1", "WT7D3_2", "WT3CO4_1", "WT4CO3_2", "WT5CO1_3", "WT7CO1_2",
                    "WT1O4_1", "WT2O3_2", "WT5O3_2", "WT7O1_1", "WT1H2_2", "WT5H3_3", "WT6H2_1", "WT8H2_2", "WT1T4_1", "WT2T1_2", "WT3T3_4", "WT6T2_4",
                    "gs1C4_2", "gs3C4_1", "gs6C2_1", "gs6C2_3", "gs7C3_2", "gs1D3_1", "gs2D1_2", "gs4D3_2", "gs5D3_3", "gs3CO3_2", "gs4CO4_2", "gs8CO1_3", "gs9CO1_2", 
                    "gs1O2_2", "gs3O3_3", "gs5O4_2", "gs6O2_2", "gs5H3_2", "gs7H4_4", "gs8H1_2", "gs10H3_4", "gs2T4_2", "gs6T3_1", "gs7T2_4", "gs10T2_5")
batch4_samples <- c("WT1C3_121", "WT4C3_211", "WT4C3_212", "WT6C3_212", "WT8C2_111", "WT3D2_211", "WT4D4_111", "WT6D4_122", "WT7D3_211", "WT3CO4_112", "WT4CO3_211", "WT5CO1_321", "WT7CO1_212",
                    "WT1O4_111", "WT2O3_212", "WT5O3_221", "WT7O1_112", "WT1H2_211", "WT5H3_312", "WT6H2_122", "WT8H2_221", "WT1T4_111", "WT2T1_212", "WT3T3_411", "WT6T2_422",
                    "gs1C4_212", "gs3C4_111", "gs6C2_311", "gs6C2_312", "gs7C3_211", "gs1D3_112", "gs2D1_211", "gs4D3_221", "gs5D3_322", "gs3CO3_212", "gs4CO4_222", "gs8CO1_322", "gs9CO1_212", 
                    "gs1O2_211", "gs3O3_322", "gs5O4_212", "gs6O2_221", "gs5H3_221", "gs7H4_422", "gs8H1_221", "gs10H3_412", "gs2T4_211", "gs6T3_111", "gs7T2_412", "gs10T2_522")
#batch3_samples <- c("gs1C", "gs3C", "gs6C", "gs7C", "gs8C", "gs1D", "gs2D", "gs4D", "gs5D", "gs3CO", "gs4CO", "gs8CO", "gs9CO", "gs1O", "gs3O", "gs5O", "gs6O", "gs5H", "gs7H", "gs8H", "gs10H", "gs2T", "gs6T", "gs7T", "gs10T")
# Update the coldata with the new condition and batch
coldata <- data.frame(row.names = colnames(counts), condition)


# Assign batches based on sample names (only for WT gen1) #for 3 batches
coldata$batch <- factor(ifelse(rownames(coldata) %in% batch1_samples, "Batch1",
                               ifelse(rownames(coldata) %in% batch2_samples, "Batch2", 
                                      ifelse(rownames(coldata) %in% batch3_samples, "Batch3", "Batch4"))))
# Check the new coldata to make sure the groups are assigned correctly
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)


# Minimum pre-filtering without count normalization
## Filter out all rows that have any counts in less than 1 library
#dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- dds[rowSums(counts(dds)) > 10, ]


# Normalize counts and do Pre-filtering: keep only those genes that have at least 3 size-normalized fragment counts in at least 2 libraries
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized = TRUE) >= 3) >= 2
##add less than 10 counts removal in pre-filtering
dds <- dds[idx, ]


# Run DESeq2
vsd <- vst(dds, blind=FALSE)

head(vsd)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)

out <- data.frame(Gene = rownames(assay(vsd)), assay(vsd), check.names = FALSE)
setwd("~/Documents/DNA_seq/RNA_seq/WGCNA")
write.csv(out, "vsd_batch.csv", row.names = FALSE)
getwd()
