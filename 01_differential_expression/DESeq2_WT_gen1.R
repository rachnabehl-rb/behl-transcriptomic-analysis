##Drought
rm(list = ls())
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))
library("genefilter")
library(dplyr)

# Here, I asses Arabidopsis leaf responses to the enviromental stresses 
# Input
setwd("~/Documents/DNA_seq/RNA_seq/WT_gen1")
df <- read.table('~/Documents/DNA_seq/RNA_seq/Count_matrix.txt', header=TRUE, row.names=1)

batch1_samples <- c("WT1C", "WT6C", "WT8C", "WT3D", "WT4D", "WT3CO", "WT7CO", "WT5O", "WT7O", "WT1H", "WT5H", "WT2T", "WT3T")
batch2_samples <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D", "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")

##Same count and condition but run it only for one pair at a time!!
### "C"  includes 6 leaf control samples
### "S" includes 4 leaf samples affected by stresses
counts <- select(df, WT1C, WT3C, WT5C, WT6C, WT8C, WT10C, WT3D, WT4D, WT6D, WT10D)
condition <- factor(c("C","C","C","C","C","C","S","S","S","S"),level=c("C","S"))
batch <- factor(ifelse(colnames(counts) %in% batch1_samples, "batch1", "batch2"),level=c("batch1","batch2"))

coldata <- data.frame(row.names=colnames(counts),condition,batch)
dds <- DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~batch+condition)

# Minimum pre-filtering without count normalization
## Filter out all rows that have any counts in less than 1 library

dds <- dds[ rowSums(counts(dds)) > 10, ]

# Normalize counts and do Pre-filtering: keep only those genes that have at least 3 size-normalized fragment counts in at least 2 libraries

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 3, ) >= 2
dds <- dds[idx, ]

# Run DESeq2
dds <- DESeq(dds, parallel=TRUE)

resultsNames(dds)

# Find all genes that respond to stress
res <- results(dds, contrast=c("condition","S","C"),parallel=TRUE)
write.table(res, file="Drought.txt",quote=F,sep="\t")

## Exporting only the results which pass an adjusted p value threshold 0.05

resOrdered <- res[order(res$padj),]
resSig05 <- subset(resOrdered, padj < 0.05)
write.table(resSig05, file="Drought pvales.txt",quote=F,sep="\t")
nrow(resSig05)
## Subset genes that passed p < 0.05 and FC >= 2 thresholds (upregulated genes)
upregulated_genes <- resSig05[resSig05$log2FoldChange >= 0, ]

write.table(upregulated_genes, file = "Drought Upregulated Genes.txt", quote = FALSE, sep = "\t")

## Subset genes that passed p < 0.05 and FC >= 2 thresholds (downregulated genes)
downregulated_genes <- resSig05[resSig05$log2FoldChange <= 0, ]
write.table(downregulated_genes, file = "Drought Downregulated Genes.txt", quote = FALSE, sep = "\t")


##CO2
rm(list = ls())
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))
library("genefilter")
library(dplyr)

# Here, I asses Arabidopsis leaf responses to the enviromental stresses 
# Input
setwd("~/Documents/DNA_seq/RNA_seq/WT_gen1")
df <- read.table('~/Documents/DNA_seq/RNA_seq/Count_matrix.txt', header=TRUE, row.names=1)

batch1_samples <- c("WT1C", "WT6C", "WT8C", "WT3D", "WT4D", "WT3CO", "WT7CO", "WT5O", "WT7O", "WT1H", "WT5H", "WT2T", "WT3T")
batch2_samples <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D", "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")

##Same count and condition but run it only for one pair at a time!!
### "C"  includes 6 leaf control samples
### "S" includes 4 leaf samples affected by stresses
counts <- select(df, WT1C, WT3C, WT5C, WT6C, WT8C, WT10C, WT3CO, WT4CO, WT5CO, WT7CO)
condition <- factor(c("C","C","C","C","C","C","S","S","S","S"),level=c("C","S"))
batch <- factor(ifelse(colnames(counts) %in% batch1_samples, "batch1", "batch2"),level=c("batch1","batch2"))

coldata <- data.frame(row.names=colnames(counts),condition,batch)
dds <- DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~batch+condition)

# Minimum pre-filtering without count normalization
## Filter out all rows that have any counts in less than 1 library

dds <- dds[ rowSums(counts(dds)) > 10, ]

# Normalize counts and do Pre-filtering: keep only those genes that have at least 3 size-normalized fragment counts in at least 2 libraries

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 3, ) >= 2
dds <- dds[idx, ]

# Run DESeq2
dds <- DESeq(dds, parallel=TRUE)

resultsNames(dds)

# Find all genes that respond to stress
res <- results(dds, contrast=c("condition","S","C"),parallel=TRUE)
write.table(res, file="CO2.txt",quote=F,sep="\t")

## Exporting only the results which pass an adjusted p value threshold 0.05

resOrdered <- res[order(res$padj),]
resSig05 <- subset(resOrdered, padj < 0.05)
write.table(resSig05, file="CO2 pvales.txt",quote=F,sep="\t")

nrow(resSig05)
## Subset genes that passed p < 0.05 and FC >= 2 thresholds (upregulated genes)
upregulated_genes <- resSig05[resSig05$log2FoldChange >= 0, ]

write.table(upregulated_genes, file = "CO2 Upregulated Genes.txt", quote = FALSE, sep = "\t")

## Subset genes that passed p < 0.05 and FC >= 2 thresholds (downregulated genes)
downregulated_genes <- resSig05[resSig05$log2FoldChange <= 0, ]
write.table(downregulated_genes, file = "CO2 Downregulated Genes.txt", quote = FALSE, sep = "\t")


##Ozone
rm(list = ls())
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))
library("genefilter")
library(dplyr)

# Here, I asses Arabidopsis leaf responses to the enviromental stresses 
# Input
setwd("~/Documents/DNA_seq/RNA_seq/WT_gen1")
df <- read.table('~/Documents/DNA_seq/RNA_seq/Count_matrix.txt', header=TRUE, row.names=1)

batch1_samples <- c("WT1C", "WT6C", "WT8C", "WT3D", "WT4D", "WT3CO", "WT7CO", "WT5O", "WT7O", "WT1H", "WT5H", "WT2T", "WT3T")
batch2_samples <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D", "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")

##Same count and condition but run it only for one pair at a time!!
### "C"  includes 6 leaf control samples
### "S" includes 4 leaf samples affected by stresses
counts <- select(df, WT1C, WT3C, WT5C, WT6C, WT8C, WT10C, WT1O, WT2O, WT5O, WT7O)
condition <- factor(c("C","C","C","C","C","C","S","S","S","S"),level=c("C","S"))
batch <- factor(ifelse(colnames(counts) %in% batch1_samples, "batch1", "batch2"),level=c("batch1","batch2"))

coldata <- data.frame(row.names=colnames(counts),condition,batch)
dds <- DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~batch+condition)

# Minimum pre-filtering without count normalization
## Filter out all rows that have any counts in less than 1 library

dds <- dds[ rowSums(counts(dds)) > 10, ]

# Normalize counts and do Pre-filtering: keep only those genes that have at least 3 size-normalized fragment counts in at least 2 libraries

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 3, ) >= 2
dds <- dds[idx, ]

# Run DESeq2
dds <- DESeq(dds, parallel=TRUE)

resultsNames(dds)

# Find all genes that respond to stress
res <- results(dds, contrast=c("condition","S","C"),parallel=TRUE)
write.table(res, file="Ozone.txt",quote=F,sep="\t")

## Exporting only the results which pass an adjusted p value threshold 0.05

resOrdered <- res[order(res$padj),]
resSig05 <- subset(resOrdered, padj < 0.05)
write.table(resSig05, file="Ozone pvales.txt",quote=F,sep="\t")


## Subset genes that passed p < 0.05 and FC >= 2 thresholds (upregulated genes)
upregulated_genes <- resSig05[resSig05$log2FoldChange >= 0, ]

write.table(upregulated_genes, file = "Ozone Upregulated Genes.txt", quote = FALSE, sep = "\t")

## Subset genes that passed p < 0.05 and FC >= 2 thresholds (downregulated genes)
downregulated_genes <- resSig05[resSig05$log2FoldChange <= 0, ]
write.table(downregulated_genes, file = "Ozone Downregulated Genes.txt", quote = FALSE, sep = "\t")


##Heat
rm(list = ls())
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))
library("genefilter")
library(dplyr)

# Here, I asses Arabidopsis leaf responses to the enviromental stresses 
# Input
setwd("~/Documents/DNA_seq/RNA_seq/WT_gen1")
df <- read.table('~/Documents/DNA_seq/RNA_seq/Count_matrix.txt', header=TRUE, row.names=1)

batch1_samples <- c("WT1C", "WT6C", "WT8C", "WT3D", "WT4D", "WT3CO", "WT7CO", "WT5O", "WT7O", "WT1H", "WT5H", "WT2T", "WT3T")
batch2_samples <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D", "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")

##Same count and condition but run it only for one pair at a time!!
### "C"  includes 6 leaf control samples
### "S" includes 4 leaf samples affected by stresses
counts <- select(df, WT1C, WT3C, WT5C, WT6C, WT8C, WT10C, WT1H, WT5H, WT6H, WT8H)
condition <- factor(c("C","C","C","C","C","C","S","S","S","S"),level=c("C","S"))
batch <- factor(ifelse(colnames(counts) %in% batch1_samples, "batch1", "batch2"),level=c("batch1","batch2"))

coldata <- data.frame(row.names=colnames(counts),condition,batch)
dds <- DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~batch+condition)

# Minimum pre-filtering without count normalization
## Filter out all rows that have any counts in less than 1 library

dds <- dds[ rowSums(counts(dds)) > 10, ]

# Normalize counts and do Pre-filtering: keep only those genes that have at least 3 size-normalized fragment counts in at least 2 libraries

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 3, ) >= 2
dds <- dds[idx, ]

# Run DESeq2
dds <- DESeq(dds, parallel=TRUE)

resultsNames(dds)

# Find all genes that respond to stress
res <- results(dds, contrast=c("condition","S","C"),parallel=TRUE)
write.table(res, file="Heat.txt",quote=F,sep="\t")

## Exporting only the results which pass an adjusted p value threshold 0.05

resOrdered <- res[order(res$padj),]
resSig05 <- subset(resOrdered, padj < 0.05)
write.table(resSig05, file="Heat pvales.txt",quote=F,sep="\t")


## Subset genes that passed p < 0.05 and FC >= 2 thresholds (upregulated genes)
upregulated_genes <- resSig05[resSig05$log2FoldChange >= 0, ]

write.table(upregulated_genes, file = "Heat Upregulated Genes.txt", quote = FALSE, sep = "\t")

## Subset genes that passed p < 0.05 and FC >= 2 thresholds (downregulated genes)
downregulated_genes <- resSig05[resSig05$log2FoldChange <= 0, ]
write.table(downregulated_genes, file = "Heat Downregulated Genes.txt", quote = FALSE, sep = "\t")


##Combination
rm(list = ls())
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))
library("genefilter")
library(dplyr)

# Here, I asses Arabidopsis leaf responses to the enviromental stresses 
# Input
setwd("~/Documents/DNA_seq/RNA_seq/WT_gen1")
df <- read.table('~/Documents/DNA_seq/RNA_seq/Count_matrix.txt', header=TRUE, row.names=1)

batch1_samples <- c("WT1C", "WT6C", "WT8C", "WT3D", "WT4D", "WT3CO", "WT7CO", "WT5O", "WT7O", "WT1H", "WT5H", "WT2T", "WT3T")
batch2_samples <- c("WT3C", "WT5C", "WT10C", "WT6D", "WT10D", "WT4CO", "WT5CO", "WT1O", "WT2O", "WT6H", "WT8H", "WT1T", "WT6T")

##Same count and condition but run it only for one pair at a time!!
### "C"  includes 6 leaf control samples
### "S" includes 4 leaf samples affected by stresses
counts <- select(df, WT1C, WT3C, WT5C, WT6C, WT8C, WT10C, WT1T, WT2T, WT3T, WT6T)
condition <- factor(c("C","C","C","C","C","C","S","S","S","S"),level=c("C","S"))
batch <- factor(ifelse(colnames(counts) %in% batch1_samples, "batch1", "batch2"),level=c("batch1","batch2"))

coldata <- data.frame(row.names=colnames(counts),condition,batch)
dds <- DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~batch+condition)

# Minimum pre-filtering without count normalization
## Filter out all rows that have any counts in less than 1 library

dds <- dds[ rowSums(counts(dds)) > 10, ]

# Normalize counts and do Pre-filtering: keep only those genes that have at least 3 size-normalized fragment counts in at least 2 libraries

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 3, ) >= 2
dds <- dds[idx, ]

# Run DESeq2
dds <- DESeq(dds, parallel=TRUE)

resultsNames(dds)

# Find all genes that respond to stress
res <- results(dds, contrast=c("condition","S","C"),parallel=TRUE)
write.table(res, file="Combination.txt",quote=F,sep="\t")

## Exporting only the results which pass an adjusted p value threshold 0.05

resOrdered <- res[order(res$padj),]
resSig05 <- subset(resOrdered, padj < 0.05)
write.table(resSig05, file="Combination pvales.txt",quote=F,sep="\t")


## Subset genes that passed p < 0.05 and FC >= 2 thresholds (upregulated genes)
upregulated_genes <- resSig05[resSig05$log2FoldChange >= 0, ]

write.table(upregulated_genes, file = "Combination Upregulated Genes.txt", quote = FALSE, sep = "\t")

## Subset genes that passed p < 0.05 and FC >= 2 thresholds (downregulated genes)
downregulated_genes <- resSig05[resSig05$log2FoldChange <= 0, ]
write.table(downregulated_genes, file = "Combination Downregulated Genes.txt", quote = FALSE, sep = "\t")
