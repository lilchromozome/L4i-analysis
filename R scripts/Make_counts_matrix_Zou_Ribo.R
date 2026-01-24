library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# featureCounts -a Homo_sapiens.GRCh38.109.gtf -o count.out -T 8 *.bam

counts <- read.table("/Users/willli/Documents/Zambidis lab/L4i RNAseq/zou bam/count.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

counts_matrix <- counts[, 7:ncol(counts)]
rownames(counts_matrix) <- counts$Geneid
colnames(counts_matrix) = c("zygote_r1", "zygote_r2", "twoC_r1", "twoC_r2", "fourC_r1", "fourC_r2", 
                            "eightC_r1", "eightC_r2", "ICM_r1", "ICM_r2")

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/L4i RNAseq/zou bam/zou_counts.csv")

ens <- rownames(counts_matrix)
symbols <- mapIds(org.Hs.eg.db,
                  keys = ens,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")
counts_matrix$gene_symbol <- symbols
rownames(counts_matrix) <- make.unique(ifelse(is.na(symbols), ens, symbols))
counts_matrix$gene_symbol <- NULL

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/L4i RNAseq/zou bam/zou_counts_gene_symbol.csv")
