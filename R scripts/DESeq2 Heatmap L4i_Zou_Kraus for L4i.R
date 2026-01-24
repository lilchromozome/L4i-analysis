library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

L4i_counts <- read.csv("/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i_counts.csv")
rownames(L4i_counts) <- L4i_counts$X
L4i_counts$X <- NULL
L4i_counts <- L4i_counts[,c( "CB62_E8", "E32C1_E8", "E32C4_E8",  "E32C6_E8", "E5C3_E8", "H9_E8", "RUES01_E8", "RUES02_E8",
                             "CB62_L4i", "E32C1_L4i", "E32C4_L4i", "E32C6_L4i", "E5C3_L4i", "H9_L4i", "RUES01_L4i","RUES02_L4i")]

zou_embryo <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/zou_counts.csv')
rownames(zou_embryo) <- zou_embryo$X
zou_embryo$X <- NULL

mESC_counts <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/mESC_PARP_KO.csv')
rownames(mESC_counts) <- mESC_counts$X
mESC_counts$X <- NULL
mESC_counts <- mESC_counts[, c("WT_r1", "WT_r2", "WT_r3", "PARPKO_r1", "PARPKO_r2", "PARPKO_r3")]

ff <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/mESC PARP1 KO/clusterekmeanZGA_Riboseq_withoocyte_Vover0_k6.csv')
# map <- c(
#   "1" = "ICM_RNA",
#   "2" = "maternal_oocyte_RNA",
#   "3" = "maternal_1C_2C_RNA",
#   "4" = "hESC_RNA",
#   "5" = "four_cell_RNA",
#   "6" = "eight_cell_RNA"
# )
# ff$cluster_label <- map[as.character(ff$ph.kmeans.cluster)]
# write.csv(ff, '/Users/willli/Documents/Zambidis lab/L4i RNAseq/clusterekmeanZGA_Riboseq_withoocyte_Vover0_k6.csv')


### L4i processing ----------------------------------
samples <- colnames(L4i_counts)
cellline <- sub("_(E8|L4i)$", "", samples)
condition <- sub("^.*_", "", samples)

coldata <- data.frame(
  cellline = factor(cellline),
  condition = factor(condition, levels = c("E8", "L4i")),
  row.names = samples
)

all(colnames(L4i_counts) == rownames(coldata))

dds <-DESeqDataSetFromMatrix(
  countData = L4i_counts,
  colData = coldata,
  design = ~ cellline + condition
)
dds <- DESeq(dds)
L4i_mat <- counts(dds, normalized = TRUE)
L4i_mat <- t(scale(t(L4i_mat)))
L4i_mat <- L4i_mat[complete.cases(L4i_mat), , drop = FALSE]

### Zou processing ----------------------------------
embryo_sample <- colnames(zou_embryo)
embryo_type <- sub("_[^_]+$", "", embryo_sample)
replicate   <- sub("^.*_", "", embryo_sample)

embryo_coldata <- data.frame(
  embryo_type = factor(embryo_type),
  replicate   = factor(replicate),
  row.names   = embryo_sample
)

zou_dds <-DESeqDataSetFromMatrix(
  countData = zou_embryo,
  colData = embryo_coldata,
  design = ~ embryo_type + replicate
)
zou_dds <- DESeq(zou_dds)
zou_mat <- counts(zou_dds, normalized = TRUE)
zou_mat <- t(scale(t(zou_mat)))
zou_mat <- zou_mat[complete.cases(zou_mat), , drop = FALSE]

### mESC processing -------------
mESC_sample <- colnames(mESC_counts)
m_celltype <- sub("_[^_]+$", "", mESC_sample)
m_replicate <- sub("^.*_", "", mESC_sample)

mESC_coldata <- data.frame(
  m_celltype = factor(m_celltype),
  m_replicate   = factor(m_replicate),
  row.names   = mESC_sample
)

mESC_dds <-DESeqDataSetFromMatrix(
  countData = mESC_counts,
  colData = mESC_coldata,
  design = ~ m_celltype + m_replicate
)
mESC_dds <- DESeq(mESC_dds)
mESC_mat <- counts(mESC_dds, normalized = TRUE)
mESC_mat <- t(scale(t(mESC_mat)))
mESC_mat <- mESC_mat[complete.cases(mESC_mat), , drop = FALSE]

### Heatmaps ----------
top_genes <- read_csv("/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i_Zou_top_by_cluster.csv")
zou_annot <- unique(top_genes$zou_annotation)

for (z in zou_annot) {
  
  genes_z <- ff$human_gene_symbol[ff$Zou_annotation == z]
  genes_zm <- ff$mouse_gene_symbol[ff$Zou_annotation == z]
  
  genes_z <- intersect(genes_z, rownames(mat_scaled))  # L4i
  genes_z <- intersect(genes_z, rownames(zou_mat_scaled))  # Zou
  genes_zm <- intersect(genes_zm, rownames(mESC_mat_scaled)) # Kraus
  
  # skip tiny/empty groups
  if (length(genes_z) < 2) next
  
  sub_L4i <- mat_scaled[genes_z, , drop = FALSE]
  sub_zou <- zou_mat_scaled[genes_z, , drop = FALSE]
  sub_mESC <- mESC_mat_scaled[genes_zm, , drop = FALSE]
  
  # row-scale (relative per gene)
  sub_L4i_scaled <- t(scale(t(sub_L4i)))
  sub_L4i_scaled <- sub_L4i_scaled[complete.cases(sub_L4i_scaled), , drop = FALSE]
  sub_zou_scaled <- t(scale(t(sub_zou)))
  sub_zou_scaled <- sub_zou_scaled[complete.cases(sub_zou_scaled), , drop = FALSE]
  sub_mESC_scaled <- t(scale(t(sub_mESC)))
  sub_mESC_scaled <- sub_mESC_scaled[complete.cases(sub_mESC_scaled), , drop = FALSE]
  
  
  E8_cols  <- grepl("_E8$", colnames(sub_L4i_scaled))
  L4i_cols <- grepl("_L4i$", colnames(sub_L4i_scaled))
  
  zygote_cols  <- grepl("^zygote", colnames(sub_zou_scaled))
  twoC_cols <- grepl("^twoC", colnames(sub_zou_scaled))
  fourC_cols <- grepl("^fourC", colnames(sub_zou_scaled))
  eightC_cols <- grepl("^eightC", colnames(sub_zou_scaled))
  ICM_cols <- grepl("^ICM", colnames(sub_zou_scaled))
  
  PARPKO_cols <- grepl("^PARPKO", colnames(sub_mESC_scaled))
  WT_cols <- grepl("WT", colnames(sub_mESC_scaled))
  
  heatmap_mat_L4i <- cbind(
    E8  = rowMeans(sub_L4i_scaled[, E8_cols]),
    L4i = rowMeans(sub_L4i_scaled[, L4i_cols])
  )
  
  heatmap_mat_zou <- cbind(
    zygote  = rowMeans(sub_zou_scaled[, zygote_cols]),
    twoC = rowMeans(sub_zou_scaled[, twoC_cols]),
    fourC = rowMeans(sub_zou_scaled[, fourC_cols]),
    eightC = rowMeans(sub_zou_scaled[, eightC_cols]),
    ICM = rowMeans(sub_zou_scaled[, ICM_cols])
  )
  
  heatmap_mat_mESC <- cbind(
    PARPKO = rowMeans(sub_mESC_scaled[, PARPKO_cols]),
    WT = rowMeans(sub_mESC_scaled[, WT_cols])
  )
  
  heatmap_mat_L4i <- heatmap_mat_L4i[order(heatmap_mat_L4i[, "L4i"], decreasing = TRUE), ]
  heatmap_mat_zou <- heatmap_mat_zou[order(heatmap_mat_L4i[, "L4i"], decreasing = TRUE), ]
  heatmap_mat_mESC <- heatmap_mat_mESC[order(heatmap_mat_L4i[, "L4i"], decreasing = TRUE),]
  
  pheatmap(heatmap_mat_L4i,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste0("Zou: ", z))
  
  pheatmap(heatmap_mat_zou,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste0("Zou: ", z))
  
  pheatmap(heatmap_mat_mESC,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste0("Zou: ", z))
  
  
}


