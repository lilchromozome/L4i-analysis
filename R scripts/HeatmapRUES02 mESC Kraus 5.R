library("DESeq2")
library(AnnotationHub)
library('ComplexHeatmap')
library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
library(writexl)
library(readxl)
library(dplyr)
library(circlize)

#Create ensemblID annotations
ah <- AnnotationHub()
query(AnnotationHub(), c("EnsDb", "Homo sapiens"))
EnsDb.Hsapiens.v109 <- ah[["AH109606"]]
EnsDb.Hsapiens.v109
edb <- EnsDb.Hsapiens.v109
seqlevelsStyle(edb) <- "UCSC"

#annotation transcription factors http://humantfs.ccbr.utoronto.ca/download.php
human_TF <- read.csv("~/../../Volumes/SSD_1/local\ data/genesets/human_TF_and_ETZ_added_PF_EnsDb_v109.csv")
rownames(human_TF) <- human_TF$gene_symbol
human_TF_counts <- as.data.frame(table(human_TF$annotation_TF2))
genes_humanTF <- rownames(human_TF)

human_TF$ensb_to_gene <- mapIds(x = EnsDb.Hsapiens.v109, keys = human_TF$Ensembl.ID,
                                column = "SYMBOL",
                                keytype = "GENEID",
                                multiVals = "first")
human_TF$ensb_to_gene[human_TF$ensb_to_gene==""]<-NA
human_TF$gene_symbol <- ifelse(is.na(human_TF$ensb_to_gene), human_TF$Ensembl.ID, human_TF$ensb_to_gene)
# human_TF$gene_symbol <- make.unique(human_TF$gene_symbol)
rownames(human_TF) <- human_TF$gene_symbol
human_TF_counts <- as.data.frame(table(human_TF$DBD))
human_TF_validated <- human_TF[human_TF$Is.TF.=="Yes", ]

#simplified annotation
human_TF$TF_class2 <- ifelse(grepl("ARID", human_TF$DBD), "ARID", 
                             ifelse(grepl("C2H2", human_TF$DBD),"C2H2_ZF", 
                                    ifelse(human_TF$DBD == "CCCH ZF", "CCCH_ZF",
                                           ifelse(grepl("ZF", human_TF$DBD),"other_ZF",
                                                  ifelse(grepl("Homeodomain", human_TF$DBD), "HOX",
                                                         ifelse(grepl("HMG", human_TF$DBD), "HMG/SOX",
                                                                ifelse(grepl("MBD", human_TF$DBD), "MBD",
                                                                       ifelse(grepl("MYB", human_TF$DBD), "MYB",
                                                                              ifelse(human_TF$DBD == "bHLH", "bHLH",
                                                                                     ifelse(human_TF$DBD == "bZIP", "bZIP",
                                                                                            ifelse(human_TF$DBD == "E2F", "E2F",
                                                                                                   ifelse(human_TF$DBD == "EBF1", "EBF1", ifelse(human_TF$DBD == "Forkhead", "FOX", ifelse(human_TF$DBD == "GATA", "GATA",
                                                                                                                                                                                           ifelse(human_TF$DBD == "mTERF", "mTERF", ifelse(human_TF$DBD == "Nuclear receptor", "Nuclear_receptor", ifelse(human_TF$DBD == "p53", "p53",
                                                                                                                                                                                                                                                                                                          ifelse(human_TF$DBD == "SMAD", "SMAD", ifelse(human_TF$DBD == "STAT", "STAT", ifelse(human_TF$DBD == "T-box", "T-box", ifelse(human_TF$DBD == "Unknown", "unknown_motif",
                                                                                                                                                                                                                                                                                                                                                                                                                                        "other motif")))))))))))))))))))))



#additional pioneer factors requested by ETZ
# CERS5, PLEKHA4, TPRXL, TPRX2, CPHXL, CPHXL2, DUXB
# ENSG00000139624, ENSG00000105559, ENSG00000180438, ENSG00000259009, ENSG00000283755, ENSG00000284484, ENSG00000282757
additional_PF <- read.csv("~/../../Volumes/SSD_1/local data/genesets/overexpressed pioneer factors.csv")
# additional_PF$gene_symbol <- make.unique(additional_PF$gene_symbol)
rownames(additional_PF) <- additional_PF$gene_symbol

transcription_regulators <- union(genes_humanTF, rownames(additional_PF))
transcription_regulators_ens <- union(human_TF$Ensembl.ID, additional_PF$EnsemblID)

# date format gene in Zou to be excluded
Zou_dates <- c("1-Sep","2-Mar.1",
               "10-Sep","2-Mar","4-Mar","5-Mar","7-Mar","8-Mar",
               "1-Mar.1","11-Sep","6-Sep","8-Sep",
               "10-Mar","11-Mar","12-Sep","14-Sep","9-Mar",
               "1-Mar","2-Sep","3-Mar","3-Sep","4-Sep","5-Sep","6-Mar","7-Sep","9-Sep")

#subsets Zou Translatome Riboseq
clusterssubsets <- read.delim("~/../../Volumes/SSD_1/local data/bulk_RNAseq/clusterekmeanZGA_RNAseq_withoocyte_k6.txt")
cluster1 <- list(clusterssubsets$cluster1)
cluster2 <- list(clusterssubsets$cluster2)
cluster3 <- list(clusterssubsets$cluster3)
cluster4 <- list(clusterssubsets$cluster4)
cluster5 <- list(clusterssubsets$cluster5)
cluster6 <- list(clusterssubsets$cluster6)

cluster1 <- lapply(cluster1, function(z){ z[!is.na(z) & z != "" & !(z %in% Zou_dates)]})
cluster2 <- lapply(cluster2, function(z){ z[!is.na(z) & z != "" & !(z %in% Zou_dates)]})
cluster3 <- lapply(cluster3, function(z){ z[!is.na(z) & z != "" & !(z %in% Zou_dates)]})
cluster4 <- lapply(cluster4, function(z){ z[!is.na(z) & z != "" & !(z %in% Zou_dates)]})
cluster5 <- lapply(cluster5, function(z){ z[!is.na(z) & z != "" & !(z %in% Zou_dates)]})
cluster6 <- lapply(cluster6, function(z){ z[!is.na(z) & z != "" & !(z %in% Zou_dates)]})

cluster1<- unlist(cluster1)
cluster2<- unlist(cluster2)
cluster3<- unlist(cluster3)
cluster4<- unlist(cluster4)
cluster5 <- unlist(cluster5)
cluster6 <- unlist(cluster6)

maternal_oocyte_RNA <- cluster2
maternal_1C_2C_RNA <- cluster3
four_cell_RNA <- cluster5
eight_cell_RNA <- cluster6
ICM_RNA <- cluster1
hESC_RNA <- cluster4

Ribo_seq_data <- read_xlsx(path = "~/../../Volumes/SSD_1/local data/bulk_RNAseq/clusterekmeanZGA_Riboseq_withoocyte_Vover0_k6.xlsx", col_names = T)
#subsets
Ribo_cluster1 <- Ribo_seq_data %>% subset(`ph$kmeans$cluster` == "1") %>% pull(gene)
Ribo_cluster2 <- Ribo_seq_data %>% subset(`ph$kmeans$cluster` == "2") %>% pull(gene)
Ribo_cluster3 <- Ribo_seq_data %>% subset(`ph$kmeans$cluster` == "3") %>% pull(gene)
Ribo_cluster4 <- Ribo_seq_data %>% subset(`ph$kmeans$cluster` == "4") %>% pull(gene)
Ribo_cluster5 <- Ribo_seq_data %>% subset(`ph$kmeans$cluster` == "5") %>% pull(gene)
Ribo_cluster6 <- Ribo_seq_data %>% subset(`ph$kmeans$cluster` == "6") %>% pull(gene)

maternal_1C2C_Ribo <- Ribo_cluster1
ICM_Ribo <- Ribo_cluster2
four_cell_Ribo <- Ribo_cluster3
hESC_Ribo <- Ribo_cluster4
eight_cell_Ribo <- Ribo_cluster5
maternal_oocyte_Ribo <- Ribo_cluster6

# eight cell added by ETZ
eight_cell_added <- c("TPRX2", "CPHXL2", "FOXD4L3", "FOXD4L6", "FOXD4L4", "BSX", "GSC2", "SHOX", "ONECUT3")
eight_cell_RNA <- union(eight_cell_RNA, eight_cell_added)
eight_cell_Ribo <- union(eight_cell_Ribo, eight_cell_added)
# combine maternal_1C_2C and 4C as preZGA
preZGA_RNA <- union(maternal_1C_2C_RNA, four_cell_RNA)
preZGA_Ribo <- union(maternal_1C2C_Ribo, four_cell_Ribo)

# interested clusters
preZGA_8C_ICM_RNA <- union(preZGA_RNA, union(eight_cell_RNA, ICM_RNA))
preZGA_8C_ICM_hESC_RNA <- union(preZGA_8C_ICM_RNA, hESC_RNA)
preZGA_8C_ICM_Ribo <- union(preZGA_Ribo, union(eight_cell_Ribo, ICM_Ribo))
preZGA_8C_ICM_hESC_Ribo <- union(preZGA_8C_ICM_Ribo, hESC_Ribo)


genes_DF <- read.delim("/Volumes/SSD_1/local data/ChIPseq/DB_peaks/Homo_sapiens.GRCh38.110.gene.bed", header=FALSE)

#DEGs RUES02
# Load summarized Experiment object
load("~/../../Volumes/SSD_1/local data/bulk_RNAseq/summarizeoverlapsobjectRUES02.Rda")
# View SumExp
se
# Load CSV
csvfile <- file.path("~/../../Volumes/SSD_1/local data/bulk_RNAseq/csvfileRUES02.csv")

#make it into an R table
(sampleTable <- read.csv(csvfile,row.names=1))
#create DEseq readable dataframe
(colData(se) <- DataFrame(sampleTable))
#create DEseq object, may change level/design
dds <- DESeqDataSet(se, design = ~ condition)


##Pre-filtering (not done)
keep <- rowSums(counts(dds)) >= 1
dds2 <- dds[keep,]

#define levels with either one
dds$condition <- relevel(dds$condition, ref = "primed")
dds2$condition <- relevel(dds$condition, ref = "primed")
#or
#dds$Status <- factor(dds$Status, levels=c("dux4","control"))

#filter lines with 0 info (not done)
#nrow(dds)
#dds <- dds[ rowSums(counts(dds)) > 1, ]
#nrow(dds)
# Filter: dds.sub removes X and Y chromosomes (not done)
#dds.sub <- dds[all(!seqnames(dds) %in% c("X", "Y")), ]
#nrow(dds.sub)

#differential expression analysis
dds <- DESeq(dds)
#or
#dds <- DESeq(dds.sub)

resultsNames(dds)
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
TIRN <- colMeans(mod_mat[dds$condition =="TIRN", ])
primed <- colMeans(mod_mat[dds$condition =="primed", ])
TIRN_Dux4 <- colMeans(mod_mat[dds$condition =="TIRN_Dux4", ])
primed_Dux4 <- colMeans(mod_mat[dds$condition =="primed_Dux4", ])
res_primed_DUX4_primed <- results(dds, contrast = primed_Dux4 - primed)
res_TIRN_primed <- results(dds, contrast = TIRN - primed)
res_TIRN_DUX4_primed <- results(dds, contrast = TIRN_Dux4 - primed)
res_TIRN_DUX4_TIRN <- results(dds, contrast = TIRN_Dux4 - TIRN)
res_TIRN_DUX4_primed_DUX4 <- results(dds, contrast = TIRN_Dux4 - primed_Dux4)
res <- cbind(res_primed_DUX4_primed, res_TIRN_primed, res_TIRN_DUX4_primed, res_TIRN_DUX4_TIRN, res_TIRN_DUX4_primed_DUX4)
res_df <- data.frame(res)

#Rename columns
colnames(res_df) <- c("baseMean_PD_P", "log2FoldChange_PD_P", "lfcSE_PD_P", "stat_PD_P", "pvalue_PD_P", "padj_PD_P", 
                      "baseMean_N_P", "log2FoldChange_N_P", "lfcSE_N_P", "stat_N_P", "pvalue_N_P", "padj_N_P", 
                      "baseMean_ND_P", "log2FoldChange_ND_P", "lfcSE_ND_P", "stat_ND_P", "pvalue_ND_P", "padj_ND_P",
                      "baseMean_ND_N", "log2FoldChange_ND_N", "lfcSE_ND_N", "stat_ND_N", "pvalue_ND_N", "padj_ND_N",
                      "baseMean_ND_PD", "log2FoldChange_ND_PD", "lfcSE_ND_PD", "stat_ND_PD", "pvalue_ND_PD", "padj_ND_PD")

#Create ensemblID annotations
ah <- AnnotationHub()
query(AnnotationHub(), c("EnsDb", "Homo sapiens"))
EnsDb.Hsapiens.v110 <- ah[["AH113665"]]
EnsDb.Hsapiens.v110

#Translate EnsemblIDs
res_df$EnsemblID <- rownames(res_df)
res_df$gene_symbol <- mapIds(x = EnsDb.Hsapiens.v110, keys = res_df$EnsemblID,
                             column = "SYMBOL",
                             keytype = "GENEID",
                             multiVals = "first")

res_df$gene_symbol2 <- res_df$gene_symbol
res_df$gene_symbol2[res_df$gene_symbol2==""]<-NA
res_df$gene_symbol2 <- ifelse(is.na(res_df$gene_symbol2), res_df$EnsemblID, res_df$gene_symbol2)
res_df$gene_symbol <- make.unique(res_df$gene_symbol2)
rownames(res_df) <- res_df$gene_symbol
res_df <- subset(res_df, select = -c(gene_symbol2))

res_df_annotated <- merge(res_df, human_TF, by = "gene_symbol", all.x = TRUE)
rownames(res_df_annotated) <- res_df_annotated$gene_symbol

# DEG
res_df_PD_P_UP <- subset(res_df, padj_PD_P < 0.05 & log2FoldChange_PD_P > 0.584962500721156)
genes_PD_P_UP <- unlist(res_df_PD_P_UP$gene_symbol)
res_df_PD_P_DOWN <- subset(res_df, padj_PD_P < 0.05 & log2FoldChange_PD_P < -0.584962500721156)
genes_PD_P_DOWN <- unlist(res_df_PD_P_DOWN$gene_symbol)

res_df_N_P_UP <- subset(res_df, padj_N_P < 0.05 & log2FoldChange_N_P > 0.584962500721156)
genes_N_P_UP <- unlist(res_df_N_P_UP$gene_symbol)
genes_N_P_UP_Ens <- unlist(res_df_N_P_UP$EnsemblID)
res_df_N_P_DOWN <- subset(res_df, padj_N_P < 0.05 & log2FoldChange_N_P < -0.584962500721156)
genes_N_P_DOWN <- unlist(res_df_N_P_DOWN$gene_symbol)

res_df_ND_P_UP <- subset(res_df, padj_ND_P < 0.05 & log2FoldChange_ND_P > 0.584962500721156)
genes_ND_P_UP <- unlist(res_df_ND_P_UP$gene_symbol)
res_df_ND_P_DOWN <- subset(res_df, padj_ND_P < 0.05 & log2FoldChange_ND_P < -0.584962500721156)
genes_ND_P_DOWN <- unlist(res_df_ND_P_DOWN$gene_symbol)

genes_UP <- union(genes_PD_P_UP, union(genes_N_P_UP, genes_ND_P_UP))
genes_DOWN <- union(genes_PD_P_DOWN, union(genes_N_P_DOWN, genes_ND_P_DOWN))
genes_DEG <- union(genes_UP, genes_DOWN)


genes_humanTF <- rownames(human_TF)

genes_humanTF_DEG <- intersect(genes_humanTF, genes_DEG)
genes_humanTF_UP <- intersect(genes_humanTF, genes_UP)

# heatmap
###########
rld_symbols <- rlog(dds)
rld_symbols <- rld_symbols[!is.na(rownames(rld_symbols))]
mat1 <- assay(rld_symbols)
#batch correction limma
mm <- model.matrix(~condition, colData(rld_symbols))
mat1 <- limma::removeBatchEffect(mat1, batch=rld_symbols$batch, design=mm)

#mean substraction
mat1 <- mat1 - rowMeans(mat1)
#Translate EnsemblIDs
mat1 <- data.frame(mat1)
mat1$EnsemblID <- rownames(mat1)
mat1$gene_symbol <- mapIds(x = EnsDb.Hsapiens.v110, keys = mat1$EnsemblID,
                           column = "SYMBOL",
                           keytype = "GENEID",
                           multiVals = "first")

mat1$gene_symbol2 <- mat1$gene_symbol
mat1$gene_symbol2[mat1$gene_symbol2==""]<-NA
mat1$gene_symbol2 <- ifelse(is.na(mat1$gene_symbol2), mat1$EnsemblID, mat1$gene_symbol2)
mat1$gene_symbol <- make.unique(mat1$gene_symbol2)
rownames(mat1) <- mat1$gene_symbol
mat1 <- subset(mat1, select = -c(gene_symbol2))
mat1 <- subset(mat1, select = -c(EnsemblID, gene_symbol))
var <- apply(mat1, 1, var)
selectedGenesV <- names(var[var>0])

#filter matrix with genes_TF_DEG and Zou
filteredTIRN_UP_V <- intersect(genes_N_P_UP, selectedGenesV)

filteredTF_V <- intersect(intersect(genes_humanTF, selectedGenesV), preZGA_8C_ICM)
filteredTF_V_TF_DEG <- intersect(intersect(genes_humanTF_DEG, selectedGenesV), preZGA_8C_ICM)
filteredTF_V_TF_UP <- intersect(intersect(genes_humanTF_UP, selectedGenesV), preZGA_8C_ICM)

res_df_TF_UP$Zou_annotation <- ifelse(res_df_TF_UP$gene_symbol %in% preZGA, "preZGA", 
                                      ifelse (res_df_TF_UP$gene_symbol %in% eight_cell, "eight cell", 
                                              ifelse (res_df_TF_UP$gene_symbol %in% ICM, "ICM", "NA")))
#write_xlsx(res_df_TF_UP, path = "~/../../Volumes/Extreme SSD/local data/bulk_RNAseq/bulkRNAseq_RUES02_FC_TF_UP_in_heatmap_with_Zou_annotations.xlsx", col_names = T)


HOX <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF3 == "HOX") %>% 
  pull(gene_symbol)
HOX_ANTP <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "HOX_ANTP") %>% 
  pull(gene_symbol)
HOX_PRD <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "HOX_PRD") %>% 
  pull(gene_symbol)
FOX <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "FOX") %>% 
  pull(gene_symbol)
EBF1 <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "EBF1") %>% 
  pull(gene_symbol)
C2H2_ZH_KLF <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "C2H2 ZH; KLF") %>% 
  pull(gene_symbol)
GATA <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "GATA") %>% 
  pull(gene_symbol)
TBX <- res_df_annotated %>%
  subset(res_df_annotated$annotation_TF_class2 == "T-box") %>% 
  pull(gene_symbol)

selectedTF <- c(HOX_ANTP, HOX_PRD, FOX, EBF1, C2H2_ZH_KLF, GATA)
selectedTF_UP_cat <- c("HOX_ANTP_UP", "HOX_PRD_UP", "FOX_UP", "EBF1UP", "C2H2_ZH_KLF_UP", "GATA_UP")
selectedTF_UP2 <- c(HOX_ANTP_UP, HOX_PRD_UP, FOX_UP, EBF1UP, C2H2_ZH_KLF_UP, GATA_UP, TBX_UP)
selectedTF_UP_cat2 <- c("HOX_ANTP_UP", "HOX_PRD_UP", "FOX_UP", "EBF1UP", "C2H2_ZH_KLF_UP", "GATA_UP", "TBX_UP")


#Kraus data

# Load summarized Experiment object
load("~/../../Volumes/SSD_2/aligned_star/summarizeoverlapsobject_mESCParp1KO_updated.Rda")
# View SumExp
se_Kraus
# Load CSV
csvfile <- file.path("~/../../Volumes/SSD_2/aligned_star/csv_mESC_Parp1KO_2.csv")

#make it into an R table
(sampleTable <- read.csv(csvfile,row.names=1))
#create DEseq readable dataframe
(colData(se_Kraus) <- DataFrame(sampleTable))
#create DEseq object, may change level/design
dds <- DESeqDataSet(se_Kraus, design = ~ condition2)

##Pre-filtering (not done)
keep <- rowSums(counts(dds)) >= 1
dds2 <- dds[keep,]

#define levels with either one
dds$condition2 <- relevel(dds$condition2, ref = "WT")

#differential expression analysis
dds <- DESeq(dds)
#or
#dds <- DESeq(dds.sub)

resultsNames(dds)
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
Parp1_KO <- colMeans(mod_mat[dds$condition2 =="Parp1_KO", ])
WT <- colMeans(mod_mat[dds$condition2 =="WT", ])

res_Parp1_KO <- results(dds, contrast = Parp1_KO - WT)

res_df_Kraus <- data.frame(res_Parp1_KO)

#Create ensemblID annotations
ah <- AnnotationHub()
query(AnnotationHub(), c("EnsDb", "Homo sapiens"))
EnsDb.Hsapiens.v110 <- ah[["AH113665"]]
EnsDb.Hsapiens.v110
query(AnnotationHub(), c("EnsDb", "Mus musculus"))
EnsDb.Mmusc.v111 <- ah[["AH116340"]]
EnsDb.Mmusc.v111


#Translate EnsemblIDs
res_df_Kraus$EnsemblID <- rownames(res_df_Kraus)
res_df_Kraus$gene_symbol <- mapIds(x = EnsDb.Mmusc.v111, keys = res_df_Kraus$EnsemblID,
                             column = "SYMBOL",
                             keytype = "GENEID",
                             multiVals = "first")

mouse_to_human_ID <- read.delim("/Volumes/SSD_1/local data/bulk_RNAseq/mart_export_mouseEnsembl_to_humanEnsembl.txt")

mouse_to_human_ID <- mouse_to_human_ID %>% dplyr::rename('EnsemblID' = Gene.stable.ID, 'human_EnsemblID' = Human.gene.stable.ID, 'human_gene_symbol' = Human.gene.name)


res_df_Kraus2 <- merge(res_df_Kraus, mouse_to_human_ID[,c("EnsemblID", "human_EnsemblID", "human_gene_symbol")], by = "EnsemblID", all.x = TRUE)

genes_mESCP1_up <- res_df_Kraus2 %>% subset(log2FoldChange > log2(1.5)) %>% pull(human_gene_symbol) %>% na.omit()
genes_mESCP1_up_Ens <- res_df_Kraus2 %>% subset(log2FoldChange > log2(1.5)) %>% pull(human_EnsemblID) %>% na.omit()

# heatmap
###########
rld_symbols <- rlog(dds)
rld_symbols <- rld_symbols[!is.na(rownames(rld_symbols))]
mat2 <- assay(rld_symbols)
#batch correction limma
mm <- model.matrix(~condition2, colData(rld_symbols))
#mat2 <- limma::removeBatchEffect(mat2, batch=rld_symbols$batch, design=mm)

#mean substraction
mat2 <- mat2 - rowMeans(mat2)
#Translate EnsemblIDs
mat2 <- data.frame(mat2)
mat2$EnsemblID <- rownames(mat2)
mat2$gene_symbol <- mapIds(x = EnsDb.Mmusc.v111, keys = mat2$EnsemblID,
                           column = "SYMBOL",
                           keytype = "GENEID",
                           multiVals = "first")

mat2 <- merge(mat2, mouse_to_human_ID[,c("EnsemblID", "human_EnsemblID", "human_gene_symbol")], by = "EnsemblID", all.x = TRUE)

mat2$gene_symbol2 <- mat2$human_gene_symbol
mat2$gene_symbol2[mat2$gene_symbol2 ==""] <- NA
mat2$gene_symbol2 <- ifelse(is.na(mat2$gene_symbol2), mat2$human_EnsemblID, mat2$gene_symbol2)
mat2$gene_symbol2 <- ifelse(is.na(mat2$gene_symbol2), mat2$EnsemblID, mat2$gene_symbol2)
mat2$gene_symbol3 <- make.unique(mat2$gene_symbol2)
rownames(mat2) <- mat2$gene_symbol3
mat2 <- subset(mat2, select = -c(gene_symbol, gene_symbol2, gene_symbol3, EnsemblID, human_EnsemblID, human_gene_symbol))

var <- apply(mat2, 1, var)
selectedGenesV2 <- names(var[var>0])


#Gambini data

# Load summarized Experiment object
load("~/../../Volumes/SSD_2/aligned_star/summarizeoverlapsobject_IWR_2C_embryo.Rda")
# View SumExp
se_Gambini
# Load CSV
csvfile <- file.path("~/../../Volumes/SSD_2/aligned_star/csv_IWR_mouse_embryo.csv")

#make it into an R table
(sampleTable <- read.csv(csvfile,row.names=1))
#create DEseq readable dataframe
(colData(se_Gambini) <- DataFrame(sampleTable))
#create DEseq object, may change level/design
dds <- DESeqDataSet(se_Gambini, design = ~ condition2)

##Pre-filtering (not done)
keep <- rowSums(counts(dds)) >= 1
dds2 <- dds[keep,]

#define levels with either one
dds$condition2 <- relevel(dds$condition2, ref = "DMSO")

#differential expression analysis
dds <- DESeq(dds)
#or
#dds <- DESeq(dds.sub)

resultsNames(dds)
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
IWR_2C <- colMeans(mod_mat[dds$condition2 =="IWR", ])
DMSO_2C <- colMeans(mod_mat[dds$condition2 =="DMSO", ])

res_IWR_2C <- results(dds, contrast = IWR_2C - DMSO_2C)

res_df_Gambini <- data.frame(res_IWR_2C)

#Create ensemblID annotations
ah <- AnnotationHub()
query(AnnotationHub(), c("EnsDb", "Homo sapiens"))
EnsDb.Hsapiens.v110 <- ah[["AH113665"]]
EnsDb.Hsapiens.v110
query(AnnotationHub(), c("EnsDb", "Mus musculus"))
EnsDb.Mmusc.v111 <- ah[["AH116340"]]
EnsDb.Mmusc.v111


#Translate EnsemblIDs
res_df_Gambini$EnsemblID <- rownames(res_df_Gambini)
res_df_Gambini$gene_symbol <- mapIds(x = EnsDb.Mmusc.v111, keys = res_df_Gambini$EnsemblID,
                                   column = "SYMBOL",
                                   keytype = "GENEID",
                                   multiVals = "first")

res_df_Gambini2 <- merge(res_df_Gambini, mouse_to_human_ID[,c("EnsemblID", "human_EnsemblID", "human_gene_symbol")], by = "EnsemblID", all.x = TRUE)

genes_IWR_2C_up <- res_df_Gambini2 %>% subset(log2FoldChange > log2(1.5)) %>% pull(human_gene_symbol) %>% na.omit()
genes_IWR_2C_up_Ens <- res_df_Gambini2 %>% subset(log2FoldChange > log2(1.5)) %>% pull(human_EnsemblID) %>% na.omit()

# heatmap
###########
rld_symbols <- rlog(dds)
rld_symbols <- rld_symbols[!is.na(rownames(rld_symbols))]
mat3 <- assay(rld_symbols)
#batch correction limma
mm <- model.matrix(~condition2, colData(rld_symbols))
#mat3 <- limma::removeBatchEffect(mat3, batch=rld_symbols$batch, design=mm)

#mean substraction
mat3 <- mat3 - rowMeans(mat3)
#Translate EnsemblIDs
mat3 <- data.frame(mat3)
mat3$EnsemblID <- rownames(mat3)
mat3$gene_symbol <- mapIds(x = EnsDb.Mmusc.v111, keys = mat3$EnsemblID,
                           column = "SYMBOL",
                           keytype = "GENEID",
                           multiVals = "first")

mat3 <- merge(mat3, mouse_to_human_ID[,c("EnsemblID", "human_EnsemblID", "human_gene_symbol")], by = "EnsemblID", all.x = TRUE)

mat3$gene_symbol2 <- mat3$human_gene_symbol
mat3$gene_symbol2[mat3$gene_symbol2 ==""] <- NA
mat3$gene_symbol2 <- ifelse(is.na(mat3$gene_symbol2), mat3$human_EnsemblID, mat3$gene_symbol2)
mat3$gene_symbol2 <- ifelse(is.na(mat3$gene_symbol2), mat3$EnsemblID, mat3$gene_symbol2)
mat3$gene_symbol3 <- make.unique(mat3$gene_symbol2)
rownames(mat3) <- mat3$gene_symbol3
mat3 <- subset(mat3, select = -c(gene_symbol, gene_symbol2, gene_symbol3, EnsemblID, human_EnsemblID, human_gene_symbol))

var <- apply(mat3, 1, var)
selectedGenesV3 <- names(var[var>0])



#filter matrix with genes_TF_DEG and Zou
genes_TIRN_up_mESP1KO_up <- intersect(genes_N_P_UP, genes_mESCP1_up)
genes_TIRN_up_mESP1KO_up_Ens <- intersect(genes_N_P_UP_Ens, genes_mESCP1_up_Ens)
genes_TIRN_up_mESP1KO_up_IWR2C_up <- intersect(intersect(genes_N_P_UP, genes_mESCP1_up), genes_IWR_2C_up)
genes_TIRN_up_mESP1KO_up_IWR2C_up_Ens <- intersect(intersect(genes_N_P_UP_Ens, genes_mESCP1_up_Ens), genes_IWR_2C_up_Ens)

genes_Kraus_TIRN_DF <- genes_DF %>% subset(V4 %in% genes_TIRN_up_mESP1KO_up_Ens)
#write.table(genes_Kraus_TIRN_DF, file="~/../../Volumes/Extreme SSD/local data/ChIPseq/DB_peaks/Homo_sapiens.GRCh38.110.genes_Parp1KO_TIRNup.txt", sep="\t", quote=F, row.names=F)

filteredTIRN_V <- intersect(intersect(intersect(genes_N_P_UP, selectedGenesV), selectedGenesV2), selectedGenesV3)
filteredTIRN_mESCP1_V <- intersect(intersect(intersect(genes_TIRN_up_mESP1KO_up, selectedGenesV), selectedGenesV2), selectedGenesV3)
filteredTIRN_mESCP1_IWR2C_V <- intersect(intersect(intersect(genes_TIRN_up_mESP1KO_up_IWR2C_up, selectedGenesV), selectedGenesV2), selectedGenesV3)

#add row means scaled data
filteredmat1 <- mat1[filteredTIRN_V, ]
filteredmat2 <- mat2[filteredTIRN_V, ]
filteredmat3 <- mat3[filteredTIRN_V, ]


filteredmat1 <- mat1[filteredTIRN_mESCP1_V, ]
filteredmat2 <- mat2[filteredTIRN_mESCP1_V, ]
filteredmat3 <- mat3[filteredTIRN_mESCP1_V, ]

filteredmat1 <- mat1[filteredTIRN_mESCP1_IWR2C_V, ]
filteredmat2 <- mat2[filteredTIRN_mESCP1_IWR2C_V, ]
filteredmat3 <- mat3[filteredTIRN_mESCP1_IWR2C_V, ]

filteredmat1 <- as.matrix(filteredmat1)
filteredmat1_scale <- t(scale(t(filteredmat1)))
filteredmat1_scale <- as.data.frame(filteredmat1_scale)
filteredmat1_scale$primed <- rowMeans(filteredmat1_scale[,c("RUES02E8_1_001", "RUES02E8_1_002", "RUES02E8_1_003", "RUES02E8_001", "RUES02E8_002", "RUES02E8_003")])
filteredmat1_scale$PD <- rowMeans(filteredmat1_scale[,c("RUES02E8Dux4_001", "RUES02E8Dux4_002", "RUES02E8Dux4_003")])
filteredmat1_scale$TIRN <- rowMeans(filteredmat1_scale[,c("RUES02L3i_1_001", "RUES02L3i_1_002", "RUES02L3i_1_003", "RUES02L3i_001", "RUES02L3i_002", "RUES02L3i_003")])
filteredmat1_scale$TD <- rowMeans(filteredmat1_scale[,c("RUES02L3iDux4_1_001", "RUES02L3iDux4_1_002", "RUES02L3iDux4_1_003", "RUES02L3iDux4_001", "RUES02L3iDux4_002", "RUES02L3iDux4_003")])
labels <- c("P", "PD", "T", "TD")
ordercol_mean <- c("primed", "PD", "TIRN", "TD")

filteredmat2 <- as.matrix(filteredmat2)
filteredmat2_scale <- t(scale(t(filteredmat2)))
filteredmat2_scale <- as.data.frame(filteredmat2_scale)
filteredmat2_scale$mESC <- rowMeans(filteredmat2_scale[,c("WT_mESC_001", "WT_mESC_002", "WT_mESC_003")])
filteredmat2_scale$Parp1KO_mESC <- rowMeans(filteredmat2_scale[,c("Parp1KO_mESC_001", "Parp1KO_mESC_002", "Parp1KO_mESC_003")])
labels2 <- c("mESC", "Parp1 KO mESC")
ordercol_mean2 <- c("mESC", "Parp1KO_mESC")

filteredmat3 <- as.matrix(filteredmat3)
filteredmat3_scale <- t(scale(t(filteredmat3)))
filteredmat3_scale <- as.data.frame(filteredmat3_scale)
filteredmat3_scale$IWR_2C <- rowMeans(filteredmat3_scale[,c("IWR_2C_1", "IWR_2C_2", "IWR_2C_3", "IWR_2C_4")])
filteredmat3_scale$DMSO_2C <- rowMeans(filteredmat3_scale[,c("DMSO_2C_1", "DMSO_2C_2", "DMSO_2C_3", "DMSO_2C_4")])
labels3 <- c("IWR-treated 2C", "DMSO control 2C")
ordercol_mean3 <- c("DMSO_2C", "IWR_2C")

filteredmat1_scale <- filteredmat1_scale %>% arrange(desc(RUES02L3i))
filteredmat2_scale <- filteredmat2_scale %>% arrange(desc(Parp1KO_mESC))
ht1 <- Heatmap(filteredmat1_scale[, ordercol_mean], name = "rlog",
        show_row_names = FALSE, col=rev(morecols(50)), cluster_columns = FALSE, cluster_row_slices = FALSE,
        show_row_dend = T, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
        row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht2 <- Heatmap(filteredmat2_scale[, ordercol_mean2], name = "rlog",
        show_row_names = FALSE, col=rev(morecols(50)), cluster_columns = FALSE, cluster_row_slices = FALSE,
        show_row_dend = FALSE, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
        row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht3 <- Heatmap(filteredmat3_scale[, ordercol_mean3], name = "rlog",
               show_row_names = FALSE, col=rev(morecols(50)),  cluster_columns = FALSE, cluster_row_slices = FALSE,
               show_row_dend = FALSE, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
               row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht3+ht2+ht1


###############Zou heatmap
# Zou paper FPKM
data <- read.csv("~/../../Volumes/SSD_1/local data/bulk_RNAseq/Riboseq_humanembryoFPKM_average.csv")
data <- unique(data)
data$gene <- make.unique(data$gene)
rownames(data) <- data$gene
data <- subset(data, select = -c(gene))
# data2 <- subset(data, select = -c(GV_RNA, MI_RNA, MII_RNA))
data$preZGA_Ribo <- rowMeans(data[,c("X1C_Ribo", "X2C_Ribo", "X4C_Ribo")])
data2 <- subset(data, select = c(preZGA_Ribo, X8C_Ribo, ICM_Ribo))

datamatrix <- as.matrix(data2)
var <- apply(datamatrix, 1, var)
selectedGenesV <- names(var[var>0])

log2datamatrix <- log2(datamatrix + 1)
log2var <- apply(log2datamatrix, 1, var)
log2selectedGenesV <- names(log2var[log2var>0])

# substract interested genes with Zou replicate log2


genes_TIRN_up_mESP1KO_up_in_Zou <- intersect(genes_TIRN_up_mESP1KO_up, log2selectedGenesV)
filteredTIRN_mESCP1_V4 <- intersect(intersect(intersect(intersect(genes_TIRN_up_mESP1KO_up, selectedGenesV), selectedGenesV2), selectedGenesV3), log2selectedGenesV)

log2filtereddatamatrix <- log2datamatrix[genes_TIRN_up_mESP1KO_up_in_Zou, ]

log2filtereddatamatrix_scale <- t(scale(t(log2filtereddatamatrix)))

#log2row_order <- rownames(filteredmat1_scale[Zou_GenesV,] %>% arrange(desc(TIRN)))

filteredmat1 <- mat1[filteredTIRN_mESCP1_V4, ]
filteredmat2 <- mat2[filteredTIRN_mESCP1_V4, ]
filteredmat3 <- mat3[filteredTIRN_mESCP1_V4, ]
filteredmat4 <- log2datamatrix[filteredTIRN_mESCP1_V4, ]

filteredmat1 <- as.matrix(filteredmat1)
filteredmat1_scale <- t(scale(t(filteredmat1)))
filteredmat1_scale <- as.data.frame(filteredmat1_scale)
filteredmat1_scale$primed <- rowMeans(filteredmat1_scale[,c("RUES02E8_1_001", "RUES02E8_1_002", "RUES02E8_1_003", "RUES02E8_001", "RUES02E8_002", "RUES02E8_003")])
filteredmat1_scale$PD <- rowMeans(filteredmat1_scale[,c("RUES02E8Dux4_001", "RUES02E8Dux4_002", "RUES02E8Dux4_003")])
filteredmat1_scale$TIRN <- rowMeans(filteredmat1_scale[,c("RUES02L3i_1_001", "RUES02L3i_1_002", "RUES02L3i_1_003", "RUES02L3i_001", "RUES02L3i_002", "RUES02L3i_003")])
filteredmat1_scale$TD <- rowMeans(filteredmat1_scale[,c("RUES02L3iDux4_1_001", "RUES02L3iDux4_1_002", "RUES02L3iDux4_1_003", "RUES02L3iDux4_001", "RUES02L3iDux4_002", "RUES02L3iDux4_003")])
labels <- c("P", "PD", "T", "TD")
ordercol_mean <- c("primed", "PD", "TIRN", "TD")

filteredmat2 <- as.matrix(filteredmat2)
filteredmat2_scale <- t(scale(t(filteredmat2)))
filteredmat2_scale <- as.data.frame(filteredmat2_scale)
filteredmat2_scale$mESC <- rowMeans(filteredmat2_scale[,c("WT_mESC_001", "WT_mESC_002", "WT_mESC_003")])
filteredmat2_scale$Parp1KO_mESC <- rowMeans(filteredmat2_scale[,c("Parp1KO_mESC_001", "Parp1KO_mESC_002", "Parp1KO_mESC_003")])
labels2 <- c("mESC", "Parp1 KO mESC")
ordercol_mean2 <- c("mESC", "Parp1KO_mESC")

filteredmat3 <- as.matrix(filteredmat3)
filteredmat3_scale <- t(scale(t(filteredmat3)))
filteredmat3_scale <- as.data.frame(filteredmat3_scale)
filteredmat3_scale$IWR_2C <- rowMeans(filteredmat3_scale[,c("IWR_2C_1", "IWR_2C_2", "IWR_2C_3", "IWR_2C_4")])
filteredmat3_scale$DMSO_2C <- rowMeans(filteredmat3_scale[,c("DMSO_2C_1", "DMSO_2C_2", "DMSO_2C_3", "DMSO_2C_4")])
labels3 <- c("IWR-treated 2C", "DMSO control 2C")
ordercol_mean3 <- c("DMSO_2C", "IWR_2C")

filteredmat4 <- as.matrix(filteredmat4)
filteredmat4_scale <- t(scale(t(filteredmat4)))
filteredmat4_scale <- as.data.frame(filteredmat4_scale)

#log2row_order <- rownames(filteredmat1_scale[Zou_GenesV,] %>% arrange(desc(RUES02L3i)))


ht1 <- Heatmap(filteredmat1_scale[, ordercol_mean], name = "rlog",
               show_row_names = FALSE, col=rev(morecols(50)), cluster_columns = FALSE, cluster_row_slices = FALSE,
               show_row_dend = T, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
               row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht2 <- Heatmap(filteredmat2_scale[, ordercol_mean2], name = "rlog",
               show_row_names = FALSE, col=rev(morecols(50)), cluster_columns = FALSE, cluster_row_slices = FALSE,
               show_row_dend = FALSE, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
               row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht3 <- Heatmap(filteredmat3_scale[, ordercol_mean3], name = "rlog",
               show_row_names = FALSE, col=rev(morecols(50)),  cluster_columns = FALSE, cluster_row_slices = FALSE,
               show_row_dend = FALSE, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
               row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht4 <- Heatmap(filteredmat4_scale, name = "log2(FPKM+1)",
               show_row_names = FALSE, col=rev(morecols(50)), cluster_columns = FALSE, cluster_row_slices = FALSE,
               show_row_dend = FALSE, clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", 
               row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))


ht1 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat1_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat1_scale)[rownames(filteredmat1_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))


ht2 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat2_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat2_scale)[rownames(filteredmat2_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))


ht3 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat3_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat3_scale)[rownames(filteredmat3_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

ht4 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat3_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat3_scale)[rownames(filteredmat3_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

ht4+ht3+ht2+ht1
ht4+ht2+ht3+ht1
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat1_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat1_scale)[rownames(filteredmat1_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))




annotatedgenes <- intersect(intersect(genes_humanTF, preZGA_8C_ICM_hESC_Ribo), genes_TIRN_up_mESP1KO_up)
annotatedgenes2 <- c("TFAP2C", "HAND1", "HAND2", "KLF17", "KLF2", "KLF4", "KLF5", "KLF17", "PRDM1", "SALL3", "SP5",
                    "YY2", "FOXA1", "FOXA2", "FOXB1", "FOXC1", "FOXO4", "GATA2", "GATA3", "GATA4", "GATA6", "TFCP2L1", "SOX17", "CDX2",
                    "DUXA", "HOXA1", "HOXA2", "HOXB5", "HOXB9", "LEUTX", "MIXL1", "PITX1", "PITX2", "MBD2", "SOX6",
                    "NR5A2", "RUNX1", "SMAD7", "STAT4", "EOMES", "TBXT", "TBX2", "TBX3", "DNMT3L", "BMP4", "TPRXL",
                    "TPRX2", "CPHXL", "CPHXL2", "DUXB", "ZSCAN2", "TET1", "EOMES", "EBF1", "EBF2", "WNT5A")
annotatedgenes <- intersect(intersect(selectedTF, preZGA_8C_ICM_Ribo), genes_TIRN_up_mESP1KO_up)
annotatedgenes <- union(annotatedgenes, annotatedgenes2)
annotatedgenes3 <- setdiff(annotatedgenes, c("EN1", "PRRX1", "LBX2", "DLX2"))

row_order <- rownames(filteredmat1_scale %>% arrange(desc(TIRN)))
filteredmat1_scale <- filteredmat1_scale %>% arrange(desc(TIRN))

filteredmat2_scale <-  filteredmat2_scale[match(rownames(filteredmat1_scale),rownames(filteredmat2_scale)), ]
filteredmat3_scale <-  filteredmat3_scale[match(rownames(filteredmat1_scale),rownames(filteredmat3_scale)), ]
filteredmat4_scale <-  filteredmat4_scale[match(rownames(filteredmat1_scale),rownames(filteredmat4_scale)), ]

ht1 <- Heatmap(filteredmat1_scale[, ordercol_mean], name = "rlog",
        show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE, row_order = row_order,
       cluster_rows = FALSE, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15)) +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat1_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat1_scale)[rownames(filteredmat1_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

ht1

ht2 <- Heatmap(filteredmat2_scale[, ordercol_mean2], name = "rlog",row_order = row_order,
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE, row_dend_reorder = -filteredmat1_scale$TIRN,
               cluster_rows = FALSE, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
ht2 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat2_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat2_scale)[rownames(filteredmat2_scale)%in%annotatedgenes], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
ht3+ht2+ht1

ht3 <- Heatmap(filteredmat3_scale, name = "log2(FPKM+1)",row_order = row_order,
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE, row_dend_reorder = -filteredmat1_scale$RUES02L3i,
               cluster_rows = FALSE, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht3+ht2+ht1
ht1
ht2  +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat2_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat2_scale)[rownames(filteredmat2_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

ht3 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat3_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat3_scale)[rownames(filteredmat3_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))


##########add split
clusterDF <- rownames(filteredmat1_scale)
clusterDF <- as.data.frame(clusterDF)
rownames(clusterDF) <- clusterDF$clusterDF
clusterDF$cluster <- NA
clusterDF$cluster <- ifelse(clusterDF$clusterDF %in% preZGA_Ribo, "pre ZGA",
                            ifelse(clusterDF$clusterDF %in% eight_cell_Ribo, "eight cell",
                                   ifelse(clusterDF$clusterDF %in% ICM_Ribo, "ICM",
                                          ifelse(clusterDF$clusterDF %in% hESC_Ribo, "hESC", NA))))

mylevelsZ <- c("pre ZGA", "eight cell", "ICM", "hESC")
split <- factor (x = clusterDF$cluster, levels = mylevelsZ)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)


ht1 <- Heatmap(filteredmat1_scale[, ordercol_mean], name = "rlog", split = split,
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE, row_order = row_order,
               cluster_rows = FALSE, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15)) +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat1_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat1_scale)[rownames(filteredmat1_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

ht1

ht2 <- Heatmap(filteredmat2_scale[, ordercol_mean2], name = "rlog", split = split, row_order = row_order,
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE,
               cluster_rows = FALSE, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
ht2 

ht3 <- Heatmap(filteredmat3_scale[, ordercol_mean3], name = "rlog", split = split,row_order = row_order,
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE,
               cluster_rows = T, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
ht3

ht4 <- Heatmap(filteredmat4_scale, name = "log2(FPKM+1)", split = split, row_order = row_order,
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE,
               cluster_rows = FALSE, row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))
ht4
ht3+ht2+ht1+ht4
ht1
ht2  +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat2_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat2_scale)[rownames(filteredmat2_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

ht3 +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat3_scale) %in% annotatedgenes3), 
                                 labels = rownames(filteredmat3_scale)[rownames(filteredmat3_scale)%in%annotatedgenes3], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))


ht3+ht2+ht1 +ht4+ rowAnnotation(link = anno_mark(at = which(rownames(filteredmat4_scale) %in% annotatedgenes3), 
                                                labels = rownames(filteredmat4_scale)[rownames(filteredmat4_scale)%in%annotatedgenes3], 
                                                labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))





####################
genes_TIRN_TIRNDup_mESP1KO_up <- intersect(union(genes_N_P_UP, genes_ND_P_UP), genes_mESCP1_up)

filteredTIRN_TIRND_mESCP1_V <- intersect(intersect(genes_TIRN_TIRNDup_mESP1KO_up, selectedGenesV), selectedGenesV2)


filteredmat1 <- mat1[filteredTIRN_mESCP1_V, ]
filteredmat1 <- as.matrix(filteredmat1)
filteredmat1_scale <- t(scale(t(filteredmat1)))
filteredmat1_scale <- as.data.frame(filteredmat1_scale)
filteredmat1_scale$RUES02E8 <- rowMeans(filteredmat1_scale[,c("RUES02E8_1_001", "RUES02E8_1_002", "RUES02E8_1_003", "RUES02E8_001", "RUES02E8_002", "RUES02E8_003")])
filteredmat1_scale$RUES02E8Dux4 <- rowMeans(filteredmat1_scale[,c("RUES02E8Dux4_001", "RUES02E8Dux4_002", "RUES02E8Dux4_003")])
filteredmat1_scale$RUES02L3i <- rowMeans(filteredmat1_scale[,c("RUES02L3i_1_001", "RUES02L3i_1_002", "RUES02L3i_1_003", "RUES02L3i_001", "RUES02L3i_002", "RUES02L3i_003")])
filteredmat1_scale$RUES02L3iDux4 <- rowMeans(filteredmat1_scale[,c("RUES02L3iDux4_1_001", "RUES02L3iDux4_1_002", "RUES02L3iDux4_1_003", "RUES02L3iDux4_001", "RUES02L3iDux4_002", "RUES02L3iDux4_003")])
labels <- c("P", "PD", "T", "TD")
ordercol_mean <- c("RUES02E8", "RUES02L3i")
#ordercol_mean <- c("RUES02E8","RUES02E8Dux4", "RUES02L3i", "RUES02L3iDux4")
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

filteredmat2 <- mat2[filteredTIRN_TIRND_mESCP1_V, ]
filteredmat2 <- as.matrix(filteredmat2)
filteredmat2_scale <- t(scale(t(filteredmat2)))
filteredmat2_scale <- as.data.frame(filteredmat2_scale)
filteredmat2_scale$mESC <- rowMeans(filteredmat2_scale[,c("WT_mESC_001", "WT_mESC_002", "WT_mESC_003")])
filteredmat2_scale$Parp1KO_mESC <- rowMeans(filteredmat2_scale[,c("Parp1KO_mESC_001", "Parp1KO_mESC_002", "Parp1KO_mESC_003")])
labels2 <- c("mESC", "Parp1 KO mESC")

ordercol_mean2 <- c("mESC", "Parp1KO_mESC")

ht1 <- Heatmap(filteredmat1_scale[, ordercol_mean], name = "rlog",
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE, row_dend_reorder = -filteredmat1_scale$RUES02L3i,
               show_row_dend = FALSE,clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15)) +
  rowAnnotation(link = anno_mark(at = which(rownames(filteredmat1_scale) %in% annotatedgenes), 
                                 labels = rownames(filteredmat1_scale)[rownames(filteredmat1_scale)%in%annotatedgenes], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
ht2 <- Heatmap(filteredmat2_scale[, ordercol_mean2], name = "rlog",
               show_row_names = F, col=rev(morecols(50)) , cluster_columns = FALSE, cluster_row_slices = FALSE, row_dend_reorder = -filteredmat1_scale$RUES02L3i,
               show_row_dend = FALSE,,clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", row_names_gp = gpar(fontsize = 7), row_gap = unit(1.5, "mm"), row_title_rot = 0, column_names_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15))

ht2+ht1

