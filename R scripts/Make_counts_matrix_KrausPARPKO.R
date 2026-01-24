library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)

# featureCounts -a Mus_musculus.GRCm39.115.gtf -o count.out -T 8 *.bam

counts <- read.table("/Users/willli/Documents/Zambidis lab/L4i RNAseq/Kraus BAM/count.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

counts_matrix <- counts[, 7:ncol(counts)]
rownames(counts_matrix) <- counts$Geneid
colnames(counts_matrix) = c("PARPKO_r1", "PARPKO_r2",  "PARPKO_r3", 
                            "WT_r1", "WT_r2", "WT_r3")

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/L4i RNAseq/Kraus BAM/mESC_PARP_KO.csv")

ens <- rownames(counts_matrix)
symbols <- mapIds(org.Mm.eg.db,
                  keys = ens,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")
counts_matrix$gene_symbol <- symbols
rownames(counts_matrix) <- make.unique(ifelse(is.na(symbols), ens, symbols))
counts_matrix$gene_symbol <- NULL

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/L4i RNAseq/Kraus BAM/mESC_PARP_KO_gene_symbol.csv")
