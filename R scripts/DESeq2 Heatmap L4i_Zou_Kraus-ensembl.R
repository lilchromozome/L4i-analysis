library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(grid)
library(biomaRt)


L4i_counts <- read.csv("/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i analysis/L4i_counts.csv") #MAC
# L4i_counts <- read.csv("~/Dr. Z lab/L4i RNA seq/L4i-analysis/L4i_counts.csv")             #Windows
rownames(L4i_counts) <- L4i_counts$X
L4i_counts$X <- NULL
L4i_counts <- L4i_counts[,c( "CB62_E8", "E32C1_E8", "E32C4_E8",  "E32C6_E8", "E5C3_E8", "H9_E8", "RUES01_E8", "RUES02_E8",
                             "CB62_L4i", "E32C1_L4i", "E32C4_L4i", "E32C6_L4i", "E5C3_L4i", "H9_L4i", "RUES01_L4i","RUES02_L4i")]
# L4i_counts <- L4i_counts[,c( "RUES02_E8", "E5C3_E8",  "E32C1_E8",
#                              "RUES02_L4i",  "E5C3_L4i", "E32C1_L4i")]


zou_embryo <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i analysis/zou_counts.csv') #MAC
# zou_embryo <- read.csv("~/Dr. Z lab/L4i RNA seq/L4i-analysis/zou_counts.csv")               #Windows
rownames(zou_embryo) <- zou_embryo$X
zou_embryo$X <- NULL

mESC_counts <- read.csv('/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i analysis/mESC_PARP_KO.csv') #MAC
# mESC_counts <- read.csv("~/Dr. Z lab/L4i RNA seq/L4i-analysis/mESC_PARP_KO.csv")                          #Windows
rownames(mESC_counts) <- mESC_counts$X
mESC_counts$X <- NULL
mESC_counts <- mESC_counts[, c('WT_r1', 'WT_r2', 'WT_r3', 'PARPKO_r1', 'PARPKO_r2', 'PARPKO_r3')]

ff <- read.delim('/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i analysis/clusterekmeanZGA_RNAseq_withoocyte_k6.txt')   #MAC
# ff <-  read.delim("~/Dr. Z lab/L4i RNA seq/L4i-analysis/clusterekmeanZGA_RNAseq_withoocyte_k6.txt") # Windows
ff <- ff %>% rename(ICM = cluster1, maternal_oocyte = cluster2, maternal_1C_2C= cluster3, hESC=cluster4, four_cell=cluster5, eight_cell=cluster6)

# PULLED FROM https://resources.aertslab.org/cistarget/tf_lists/ 
tf <- read.delim('/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i analysis/allTFs_hg38.txt')  #MAC
# tf <- read.delim('~/Dr. Z lab/L4i RNA seq/L4i-analysis/TF_names_v_1.01.txt')  #Windows
tfs_to_add <- c('TPRXL', 'TPRX1', 'TPRX2', 'CPHXL', 'CPHXL2', 'DUXB',
                "FOXQ1","RUNX1","TBXT","MIXL1","HOXA1","EBF1","EBF2","HOXD8","HOXA7",
                "HOXD10","HOXA2","FOXD4","GBX2","GATA6","HOXC8","HOXB8","HOXA5","HOXB2",
                "DLX4","GATA5","FOXA2","CDX2","GATA4","SP5","SOX17","FOXE1","HOXB3",
                "HOXB4","GATA3","HAND1","HOXB13")
colnames(tf)[1] <- "tf"
tf <- rbind(tf, data.frame(tf = tfs_to_add))
tf <- unique(tf)

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
res <- results(dds, contrast = c("condition", "L4i", "E8"))
res <- res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange>0, ]
mat <- counts(dds, normalized = TRUE)
mat <- mat[rownames(res), , drop = FALSE]
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

m_res <- results(mESC_dds, contrast = c("m_celltype", "PARPKO", "WT"))
m_res <- m_res[!is.na(m_res$padj) & m_res$padj < 0.05 & m_res$log2FoldChange>0, ]
mESC_mat <- mESC_mat[rownames(m_res), , drop = FALSE]

mESC_mat_scaled <- t(scale(t(mESC_mat)))
mESC_mat_scaled <- mESC_mat_scaled[complete.cases(mESC_mat_scaled), , drop = FALSE]


### Plot Heatmaps ------------------
ht_list <- NULL
row_blocks <- list()
n_gene <- numeric()
output_table <- data.frame()

interesting_genes <- c(
  "FOXQ1","RUNX1","TBXT","MIXL1","HOXA1","EBF1","EBF2","HOXD8","HOXA7",
  "HOXD10","HOXA2","FOXD4","GBX2","GATA6","HOXC8","HOXB8","HOXA5","HOXB2",
  "DLX4","GATA5","FOXA2","CDX2","GATA4","SP5","SOX17","FOXE1","HOXB3",
  "HOXB4","GATA3","HAND1","HOXB13", tf$tf
)

label_fun <- function(labs) {
  labs <- as.character(labs)
  col <- ifelse(labs %in% interesting_genes, "black", "transparent")
  gpar(fontsize = 8, col = col)
}

Zou_annot <- colnames(ff)

L4i_ids <- rownames(mat_scaled)
zou_ids <- rownames(zou_mat_scaled)
mESC_ids  <- rownames(mESC_mat_scaled)

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# mouse ensemble -> human ensembl
map_m2h <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  values  = mESC_ids,
  mart    = mouse_mart
)

map_m2h <- map_m2h[
  !is.na(map_m2h$ensembl_gene_id) &
    !is.na(map_m2h$hsapiens_homolog_ensembl_gene) &
    !is.na(map_m2h$hsapiens_homolog_associated_gene_name)&
    map_m2h$ensembl_gene_id != "" &
    map_m2h$hsapiens_homolog_ensembl_gene != "" &
    map_m2h$hsapiens_homolog_associated_gene_name != "",
]



Zou_annot <- c('maternal_1C_2C', 'four_cell', 'eight_cell', 'ICM')
row_titles <- c('Zygote-2C', '4C', '8C-morula', 'Naive ICM')
for (i in seq_along(Zou_annot)) {
  z = Zou_annot[i]
  ## filter genes ------------------------
  cluster_genes <- ff[[z]]
  gene_set <- map_m2h[map_m2h$hsapiens_homolog_associated_gene_name %in% cluster_genes, ]
  # gene_set <- gene_set[gene_set$hsapiens_homolog_associated_gene_name %in% tf$tf, ]
  
  gene_set <- gene_set[gene_set$hsapiens_homolog_ensembl_gene %in% rownames(mat_scaled) &
              gene_set$hsapiens_homolog_ensembl_gene %in% rownames(zou_mat_scaled) &
              gene_set$ensembl_gene_id %in% rownames(mESC_mat_scaled), ]
  
  
  sub_L4i <- mat_scaled[gene_set$hsapiens_homolog_ensembl_gene, , drop = FALSE]
  sub_zou <- zou_mat_scaled[gene_set$hsapiens_homolog_ensembl_gene, , drop = FALSE]
  sub_mESC <- mESC_mat_scaled[gene_set$ensembl_gene_id, , drop = FALSE]
  
  
  
  
  ## Make heatmap mats -----------------
  E8_cols  <- grepl("_E8$", colnames(sub_L4i))
  L4i_cols <- grepl("_L4i$", colnames(sub_L4i))
  
  zygote_cols  <- grepl("^zygote", colnames(sub_zou))
  twoC_cols <- grepl("^twoC", colnames(sub_zou))
  fourC_cols <- grepl("^fourC", colnames(sub_zou))
  eightC_cols <- grepl("^eightC", colnames(sub_zou))
  ICM_cols <- grepl("^ICM", colnames(sub_zou))
  
  PARPKO_cols <- grepl("^PARPKO", colnames(sub_mESC))
  WT_cols <- grepl("WT", colnames(sub_mESC))
  
  heatmap_mat_L4i <- cbind(
    E8  = rowMeans(sub_L4i[, E8_cols]),
    L4i = rowMeans(sub_L4i[, L4i_cols])
  )
  
  heatmap_mat_zou <- cbind(
    zygote  = rowMeans(sub_zou[, zygote_cols]),
    twoC = rowMeans(sub_zou[, twoC_cols]),
    fourC = rowMeans(sub_zou[, fourC_cols]),
    eightC = rowMeans(sub_zou[, eightC_cols]),
    ICM = rowMeans(sub_zou[, ICM_cols])
  )
  
  heatmap_mat_mESC <- cbind(
    WT = rowMeans(sub_mESC[, WT_cols]),
    PARPKO = rowMeans(sub_mESC[, PARPKO_cols])
  )
  

  ## maintain same order ----------------
  target_order <- rownames(heatmap_mat_L4i)[order(heatmap_mat_L4i[, "L4i"], decreasing = TRUE)]
  ordered_gene_set <- gene_set[match(target_order, gene_set$hsapiens_homolog_ensembl_gene), ]
  heatmap_mat_L4i <- heatmap_mat_L4i[ordered_gene_set$hsapiens_homolog_ensembl_gene, , drop = FALSE]
  heatmap_mat_zou <- heatmap_mat_zou[ordered_gene_set$hsapiens_homolog_ensembl_gene, , drop = FALSE]
  heatmap_mat_mESC <- heatmap_mat_mESC[ordered_gene_set$ensembl_gene_id, , drop = FALSE]
  
  ### build per-gene output table (same ordering as the heatmaps) -------
  df_out <- data.frame(
    human_gene_name = ordered_gene_set$hsapiens_homolog_associated_gene_name,
    mouse_gene_name = ordered_gene_set$external_gene_name,
    zou_annotation  = z,
    zygote  = heatmap_mat_zou[ordered_gene_set$hsapiens_homolog_ensembl_gene, "zygote"],
    twoC    = heatmap_mat_zou[ordered_gene_set$hsapiens_homolog_ensembl_gene, "twoC"],
    fourC   = heatmap_mat_zou[ordered_gene_set$hsapiens_homolog_ensembl_gene, "fourC"],
    eightC  = heatmap_mat_zou[ordered_gene_set$hsapiens_homolog_ensembl_gene, "eightC"],
    ICM     = heatmap_mat_zou[ordered_gene_set$hsapiens_homolog_ensembl_gene, "ICM"],
    E8      = heatmap_mat_L4i[ordered_gene_set$hsapiens_homolog_ensembl_gene, "E8"],
    L4i     = heatmap_mat_L4i[ordered_gene_set$hsapiens_homolog_ensembl_gene, "L4i"],
    WT_mESC     = heatmap_mat_mESC[ordered_gene_set$ensembl_gene_id, "WT"],
    PARPKO_mESC = heatmap_mat_mESC[ordered_gene_set$ensembl_gene_id, "PARPKO"],
    row.names = NULL,
    check.names = FALSE
  )
  output_table <- rbind(output_table, df_out)
  
  ### ComplexHeatmap ---------
  rowlabs_h <- gene_set$hsapiens_homolog_associated_gene_name
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
  
  
  n_gene[i] <- nrow(sub_L4i)
  block_h <- unit(0.75 * n_gene[i] , "mm")
  
  ht_zou <- Heatmap(heatmap_mat_zou, name = "Zou",
                    show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE,
                    row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 90, heatmap_height = block_h,
                    column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
                    # , row_title = row_titles[i])
  
  ht_mESC <- Heatmap(heatmap_mat_mESC, name = "mESC",
                     show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE,
                     row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, heatmap_height = block_h,
                     column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
  
  ht_L4i <- Heatmap(heatmap_mat_L4i, name = "L4i",
                    show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE,
                    row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, heatmap_height = block_h,
                    column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15),
                    right_annotation = right_anno)
  
  
  # labels: only interesting genes, others blank
  rowlabs_h <- ordered_gene_set$hsapiens_homolog_associated_gene_name
  rowlabs_h <- ifelse(rowlabs_h %in% interesting_genes, rowlabs_h, "")


  row_blocks[[z]] <- ht_zou + ht_mESC + ht_L4i
  draw(
    row_blocks[[z]], 
    ht_gap = unit(2, "mm"),
    column_title = paste0(z, ":L4i and mESC upregulated"),
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
    )
}


nr <- length(Zou_annot)

grid.newpage()
hts <- unit(n_gene*1.1, "null")
pushViewport(viewport(layout = grid.layout(nrow = nr, ncol = 1, heights = hts)))

for (i in seq_along(Zou_annot)) {
  z <- Zou_annot[i]

  pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))

  draw(
    row_blocks[[z]],
    newpage = FALSE,
    show_heatmap_legend = F,
    show_annotation_legend = F,
  )

  popViewport()
}

popViewport()

# write_xlsx(output_table, '/Users/willli/Documents/Zambidis lab/L4i RNAseq/L4i analysis/TF_output_table.xlsx')
output_table <- data.frame()

