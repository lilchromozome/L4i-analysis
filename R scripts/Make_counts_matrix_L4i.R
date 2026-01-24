library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

counts <- read.table("/Users/zambidislab/Desktop/L4i_counts/count.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

counts_matrix <- counts[, 7:ncol(counts)]
rownames(counts_matrix) <- counts$Geneid
colnames(counts_matrix) = c("CB62_E8", "CB62_L4i", "E32C1_E8", "E32C1_L4i", "E32C4_E8", "E32C4_L4i", 
              "E32C6_E8", "E32C6_L4i", "E5C3_E8", "E5C3_L4i", "H9_E8", "H9_L4i",
              "RUES01_E8", "RUES01_L4i", "RUES02_E8", "RUES02_L4i")

write.csv(counts_matrix, "Desktop/L4i_counts/counts_matrix.csv")

ens <- rownames(counts_matrix)
symbols <- mapIds(org.Hs.eg.db,
                  keys = ens,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")
counts_matrix$gene_symbol <- symbols
rownames(counts_matrix) <- make.unique(ifelse(is.na(symbols), ens, symbols))
counts_matrix$gene_symbol <- NULL

write.csv(counts_matrix, "Desktop/L4i_counts/counts_matrix_gene_symbol.csv")


