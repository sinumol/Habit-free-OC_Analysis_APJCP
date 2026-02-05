####################################################
## Differential Expression Analysis using DESeq2
## Tumor vs Normal (Habit-free samples from TCGA)
### Created on 22-05-2024
### Sinumol George
####################################################
##Dataset used for teh anlysis is avaible in TCGABiolink and Cbioportal
library(DESeq2)

## Input data
count_matrix <- combined_normal_tumor_hf
sample_info <- data.frame(
  condition = factor(condition_1, levels = c("Normal", "Tumor"))
)
rownames(sample_info) <- colnames(count_matrix)

## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ condition
)

## Run DESeq2
dds <- DESeq(dds)

## Extract results (Tumor vs Normal)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

## Order by adjusted p-value
res <- res[order(res$padj), ]

####################################################
PCA Analysis
####################################################

vsd <- vst(dds, blind = FALSE)

pdf("results/PCA_tumor_vs_normal.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

####################################################
## Significant Differentially Expressed Genes
####################################################

sig_deg <- subset(
  res,
  padj < 0.05 & abs(log2FoldChange) > 1
)

####################################################
