library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

# L4i_counts <- read.csv("/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i_counts_gene_symbol.csv") #MAC
L4i_counts <- read.csv("~/Dr. Z lab/L4i RNA seq/L4i-analysis/L4i_counts_gene_symbol.csv") # Windows
rownames(L4i_counts) <- L4i_counts$X
L4i_counts$X <- NULL
L4i_counts <- L4i_counts[,c( "CB62_E8", "E32C1_E8", "E32C4_E8",  "E32C6_E8", "E5C3_E8", "H9_E8", "RUES01_E8", "RUES02_E8",
                             "CB62_L4i", "E32C1_L4i", "E32C4_L4i", "E32C6_L4i", "E5C3_L4i", "H9_L4i", "RUES01_L4i","RUES02_L4i")]
samples <- colnames(L4i_counts)
cellline <- sub("_(E8|L4i)$", "", samples)
condition <- sub("^.*_", "", samples)

coldata <- data.frame(
  cellline = factor(cellline),
  condition = factor(condition, levels = c("E8", "L4i")),
  row.names = samples
)

dds <-DESeqDataSetFromMatrix(
  countData = L4i_counts,
  colData = coldata,
  design = ~ cellline + condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "L4i", "E8"))
res <- res[!is.na(res$padj) & res$padj < 0.05, ] # & abs(res$log2FoldChange) >= 1, ]
L4i_mat <- counts(dds, normalized = TRUE)
L4i_mat <- L4i_mat[rownames(res), , drop = FALSE]
L4i_mat <- t(scale(t(L4i_mat)))
L4i_mat <- L4i_mat[complete.cases(L4i_mat), , drop = FALSE]


# ff <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/clusterekmeanZGA_Riboseq_withoocyte_Vover0_k6.csv') #MAC
ff <- read.csv("~/Dr. Z lab/L4i RNA seq/L4i-analysis/clusterekmeanZGA_Riboseq_withoocyte_Vover0_k6.csv") # Windows
### https://humantfs.ccbr.utoronto.ca/download.php 
# tf <- read.delim('/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i-analysis/TF_names_v_1.01.txt') #MAC
tf <- read.delim('~/Dr. Z lab/L4i RNA seq/L4i-analysis/TF_names_v_1.01.txt')
tfs_to_add <- c('TPRXL', 'TPRX2', 'CPHXL', 'CPHXL2', 'DUXB',
                "FOXQ1","RUNX1","TBXT","MIXL1","HOXA1","EBF1","EBF2","HOXD8","HOXA7",
                "HOXD10","HOXA2","FOXD4","GBX2","GATA6","HOXC8","HOXB8","HOXA5","HOXB2",
                "DLX4","GATA5","FOXA2","CDX2","GATA4","SP5","SOX17","FOXE1","HOXB3",
                "HOXB4","GATA3","HAND1","HOXB13")
colnames(tf)[1] <- "tf"
tf <- rbind(tf, data.frame(tf = tfs_to_add))
tf <- unique(tf)

top_by_cluster <- data.frame()

zou_annot <- unique(ff$cluster_label)

for (k in zou_annot) {
  genes_in_cluster <- ff$gene[ff$cluster_label == k]
  genes_in_cluster <- intersect(genes_in_cluster, tf$tf)
  genes_in_cluster <- intersect(genes_in_cluster, rownames(L4i_mat)) 
  if (length(genes_in_cluster) < 2) next
  
  sub_L4i <- L4i_mat[genes_in_cluster, , drop = FALSE]
  sub_L4i <- sub_L4i[complete.cases(sub_L4i), , drop = FALSE]
  
  E8_cols  <- grepl("_E8$", colnames(sub_L4i))
  L4i_cols <- grepl("_L4i$", colnames(sub_L4i))
  
  L4i_vs_E8 <- data.frame(
    E8  = rowMeans(sub_L4i[, E8_cols]),
    L4i = rowMeans(sub_L4i[, L4i_cols])
  )
  L4i_vs_E8 <- L4i_vs_E8[L4i_vs_E8[, "L4i"] > 0, , drop = FALSE]
  L4i_vs_E8 <- L4i_vs_E8[order(L4i_vs_E8[, "L4i"], decreasing = TRUE), ]
  L4i_vs_E8$zou_annotation <- k
  top_by_cluster <- rbind(top_by_cluster, L4i_vs_E8)
  # 
  # top_by_cluster[[paste0(k)]] <- L4i_vs_E8
  
  print(pheatmap(as.matrix(L4i_vs_E8[, c("E8","L4i")]),
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 main = paste0("Zou: ", k)))
  
}



