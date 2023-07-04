rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(pheatmap)
load(file = "data/rdata/matrix_GSE205431_raw.Rdata")
hub_genes <- c("EGR1", "FOS", "FOSB", "KLF2", "JUNB")
Expr <- count_anno  %>%
    na.omit()
Expr <- Expr[hub_genes, ]
data_anno <- Expr
# data_anno[is.na(data_anno)] <- 1
annotation_col <- Targets %>%
    dplyr::select(group)
ann_colors <- list(
    group = c(SX = "#af561e", CTRL = "#1b809e")
)
pdf(
    file = "res/pic/final/Heatmap_hub_genes_gse205431.pdf",
    width = 60,
    height = 4
)
pheatmap(
    data_anno,
    annotation = annotation_col,
    annotation_colors = ann_colors,
    cluster_cols = TRUE,
    cluster_row = FALSE,
    color = colorRampPalette(
        c("green", "black", "red"),
        bias = 1
    )(512),
    show_colnames = TRUE,
    scale = "row", # 矫正
    border_color = "black",
    fontsize = 24,
    fontsize_row = 20,
    fontsize_col = 24,
    angle_col = "45",
    cellwidth = 95,
    cellheight = 25,
)
dev.off()
 hub_genes <- read.table(
    file = "res/data/genes/Targets_TF_genes_OA.txt",
    sep = "\t", check.names = TRUE, header = TRUE
)
hub_genes <- hub_genes$SYMBOL %>%
    tolower() %>%
    capitalize()
Expr <- fpkm_anno %>%
    column_to_rownames(var = "gene_name")
Expr <- Expr[hub_genes, ]
genes2show <- hub_genes
i <- which(v$genes$SYMBOL %in% genes2show)
data <- lcpm_normalised[i, ] %>%
    as.data.frame()
# 基因信息
geneid <- rownames(data)
genes <- select(Mus.musculus,
    keys = geneid,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
)
data_anno <- data %>%
    rownames_to_column(var = "ENSEMBL") %>%
    left_join(genes, by = "ENSEMBL") %>%
    relocate(SYMBOL) %>%
    dplyr::select(-ENSEMBL) %>%
    column_to_rownames(var = "SYMBOL")
annotation_col <- as.data.frame(v$targets) %>%
    dplyr::select(group)
ann_colors <- list(
    group = c(SX = "#ad2938", CTRL = "#1b809e")
)
pdf(
    file = "res/pic/final/Heatmap_com_genes_gse205431.pdf",
    width = 10,
    height = 40
)
pheatmap(
    data_anno,
    annotation = annotation_col,
    annotation_colors = ann_colors,
    cluster_cols = FALSE,
    cluster_row = TRUE,
    color = colorRampPalette(
        c("green", "black", "red"),
        bias = 1
    )(512),
    show_colnames = TRUE,
    scale = "row", # 矫正
    border_color = "black",
    fontsize = 24,
    fontsize_row = 20,
    fontsize_col = 24,
    angle_col = "45",
    cellwidth = 95,
    cellheight = 25,
)
dev.off()
