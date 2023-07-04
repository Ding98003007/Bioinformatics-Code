rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(pheatmap)
load(file = "data/rdata/matrix_GSE82107_raw.Rdata")
hub_genes <- c("EGR1", "FOS", "FOSB", "KLF2", "JUNB")
Expr <- Expr_GSE82107_anno
Expr <- Expr[hub_genes, ]
data_anno <- Expr
annotation_col <- Targets %>%
    dplyr::select(group)
ann_colors <- list(
    group = c(OA = "#af561e", CTRL = "#1b809e")
)
pdf(
    file = "res/pic/final/Heatmap_hub_genes_gse82107.pdf",
    width = 30,
    height = 6
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
hub_genes <- read.table(
    file = "res/data/genes/Targets_TF_genes_OA.txt",
    sep = "\t", check.names = TRUE, header = TRUE
)
hub_genes <- hub_genes$SYMBOL
Expr <- Expr_GSE82107_anno
Expr <- Expr[hub_genes, ]
data_anno <- Expr
annotation_col <- Targets %>%
    dplyr::select(group)
ann_colors <- list(
    group = c(OA = "#af561e", CTRL = "#1b809e")
)
pdf(
    file = "res/pic/final/Heatmap_com_genes_gse82107.pdf",
    width = 30,
    height = 35
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
