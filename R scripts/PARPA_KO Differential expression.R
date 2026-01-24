library(DESeq2)

counts_matrix <- read.csv("C:/Users/willllllli/Documents/Dr. Z lab/L4i RNA seq/counts_matrix_gene_symbol.csv")\
rownames(counts_matrix) <- counts_matrix$X
counts_matrix$X <- NULL

ff <- read.csv("C:/Users/willllllli/Documents/Dr. Z lab/L4i RNA seq/hmESC_KO_genes.csv")

counts_matrix <- counts_matrix[rownames(counts_matrix) %in% ff$human_gene_symbol, ]



dds <-DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ cellline + condition
)

dds <- DESeq(dds)
result <- results(dds, contrast = c("condition", "L4i", "E8"))
result <- result[order(result$log2FoldChange), ]
head(result)