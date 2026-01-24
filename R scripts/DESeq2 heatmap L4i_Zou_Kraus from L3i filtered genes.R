library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(grid)


L4i_counts <- read.csv("/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i_counts_gene_symbol.csv")
rownames(L4i_counts) <- L4i_counts$X
L4i_counts$X <- NULL
L4i_counts <- L4i_counts[,c( "CB62_E8", "E32C1_E8", "E32C4_E8",  "E32C6_E8", "E5C3_E8", "H9_E8", "RUES01_E8", "RUES02_E8",
                                   "CB62_L4i", "E32C1_L4i", "E32C4_L4i", "E32C6_L4i", "E5C3_L4i", "H9_L4i", "RUES01_L4i","RUES02_L4i")]

zou_embryo <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/zou_counts_gene_symbol.csv')
rownames(zou_embryo) <- zou_embryo$X
zou_embryo$X <- NULL

mESC_counts <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/mESC_PARPKO_gene_symbol.csv')
rownames(mESC_counts) <- mESC_counts$X
mESC_counts$X <- NULL
mESC_counts <- mESC_counts[, c('WT_r1', 'WT_r2', 'WT_r3', 'PARPKO_r1', 'PARPKO_r2', 'PARPKO_r3')]

ff <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/mESC PARP1 KO/hmESC_KO_genes.csv')

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
mat <- counts(dds, normalized = TRUE)

mat <- mat[rownames(mat) %in% ff$human_gene_symbol, , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]


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

zou_mat <- zou_mat[rownames(zou_mat) %in% ff$human_gene_symbol, , drop = FALSE]
zou_mat_scaled <- t(scale(t(zou_mat)))
zou_mat_scaled <- zou_mat_scaled[complete.cases(zou_mat_scaled), , drop = FALSE]

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

mESC_mat <- mESC_mat[rownames(mESC_mat) %in% ff$mouse_gene_symbol, , drop = FALSE]
mESC_mat_scaled <- t(scale(t(mESC_mat)))
mESC_mat_scaled <- mESC_mat_scaled[complete.cases(mESC_mat_scaled), , drop = FALSE]


### Plot Heatmaps ------------------
ht_list <- NULL
row_blocks <- list()
output_table <- data.frame()

interesting_genes <- c(
  "FOXQ1","RUNX1","TBXT","MIXL1","HOXA1","EBF1","EBF2","HOXD8","HOXA7",
  "HOXD10","HOXA2","FOXD4","GBX2","GATA6","HOXC8","HOXB8","HOXA5","HOXB2",
  "DLX4","GATA5","FOXA2","CDX2","GATA4","SP5","SOX17","FOXE1","HOXB3",
  "HOXB4","GATA3","HAND1","HOXB13"
)


label_fun <- function(labs) {
  labs <- as.character(labs)
  col <- ifelse(labs %in% interesting_genes, "black", "transparent")
  gpar(fontsize = 8, col = col)
}

Zou_annot <- c('ICM', 'maternal', 'eight cell', 'preZGA')

for (z in Zou_annot) {
  
  df_pairs <- unique(ff[ff$Zou_annotation == z, c("human_gene_symbol", "mouse_gene_symbol")])
  df_pairs <- df_pairs[!is.na(df_pairs$human_gene_symbol) & !is.na(df_pairs$mouse_gene_symbol), ]
  
  keep <- df_pairs$human_gene_symbol %in% rownames(mat_scaled) &  # L4i
          df_pairs$human_gene_symbol %in% rownames(zou_mat_scaled) &  # Zou
          df_pairs$mouse_gene_symbol %in% rownames(mESC_mat_scaled) # Kraus
  df_pairs <- df_pairs[keep, ]
  
  # skip tiny/empty groups
  if (length(df_pairs) < 2) next
  
  sub_L4i <- mat_scaled[df_pairs$human_gene_symbol, , drop = FALSE]
  sub_zou <- zou_mat_scaled[df_pairs$human_gene_symbol, , drop = FALSE]
  sub_mESC <- mESC_mat_scaled[df_pairs$mouse_gene_symbol, , drop = FALSE]
  
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
  
  #maintain same order
  ord_genes_h <- rownames(heatmap_mat_L4i[order(heatmap_mat_L4i[, "L4i"], decreasing = TRUE), , drop=FALSE])
  heatmap_mat_L4i <- heatmap_mat_L4i[ord_genes_h, , drop = FALSE]
  heatmap_mat_zou <- heatmap_mat_zou[ord_genes_h, , drop = FALSE]
  
  map_h2m <- setNames(ff$mouse_gene_symbol, ff$human_gene_symbol)
  ord_genes_m <- unname(map_h2m[ord_genes_h])
  ord_genes_m <- ord_genes_m[ord_genes_m %in% rownames(sub_mESC)]
  heatmap_mat_mESC <- cbind(
    WT     = rowMeans(sub_mESC[, WT_cols,     drop = FALSE]),
    PARPKO = rowMeans(sub_mESC[, PARPKO_cols, drop = FALSE])
  )
  heatmap_mat_mESC <- heatmap_mat_mESC[ord_genes_m, , drop = FALSE]
  
  ### build per-gene output table (same ordering as the heatmaps) -------
  df_out <- data.frame(
    human_gene_name = ord_genes_h,
    mouse_gene_name = ord_genes_m,
    zou_annotation  = z,
    zygote  = heatmap_mat_zou[ord_genes_h, "zygote"],
    twoC    = heatmap_mat_zou[ord_genes_h, "twoC"],
    fourC   = heatmap_mat_zou[ord_genes_h, "fourC"],
    eightC  = heatmap_mat_zou[ord_genes_h, "eightC"],
    ICM     = heatmap_mat_zou[ord_genes_h, "ICM"],
    E8      = heatmap_mat_L4i[ord_genes_h, "E8"],
    L4i     = heatmap_mat_L4i[ord_genes_h, "L4i"],
    WT_mESC     = heatmap_mat_mESC[ord_genes_m, "WT"],
    PARPKO_mESC = heatmap_mat_mESC[ord_genes_m, "PARPKO"],
    row.names = NULL,
    check.names = FALSE
  )
  output_table <- rbind(output_table, df_out)

  ### ComplexHeatmap ---------
  rowlabs_h <- rownames(heatmap_mat_L4i)
  rowlabs_h <- ifelse(rowlabs_h %in% interesting_genes, rowlabs_h, "")
  
  ht_zou <- Heatmap(heatmap_mat_zou, name = "Zou",
                 show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE,
                 row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, 
                 column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
  
  ht_mESC <- Heatmap(heatmap_mat_mESC, name = "mESC",
                    show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE,
                    row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, 
                    column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
  
  ht_L4i <- Heatmap(heatmap_mat_L4i, name = "L4i",
                    show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE,
                    row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, 
                    column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
  
  # labels: only interesting genes, others blank
  rowlabs_h <- rownames(heatmap_mat_L4i)
  rowlabs_h <- ifelse(rowlabs_h %in% interesting_genes, rowlabs_h, "")
  
  idx <- which(rowlabs_h != "")
  right_anno <- rowAnnotation(
    Genes = anno_mark(
      at = idx,
      labels = rowlabs_h[idx],
      side = "right",
      labels_gp = gpar(fontsize = 9)
    ),
    width = unit(2, "cm")
  )
  
  row_block <- ht_zou + ht_mESC + ht_L4i + right_anno
  draw(row_block, ht_gap = unit(2, "mm"))
}
write.csv(output_table, '/Users/willli/Documents/Zambidis lab/L4i RNAseq/mESC PARP1 KO/output_table.csv')
