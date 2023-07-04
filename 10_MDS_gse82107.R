rm(list = ls())
library(tidyverse)
library(limma)
library(RColorBrewer)
load(file = "data/rdata/matrix_GSE82107_raw.Rdata")
pdf(
    file = "res/pic/final/MDSplots_gse82107.pdf", # nolint
    width = 10,
    height = 8
)
par(mfrow = c(1, 1))
group <- factor(Targets$group,
    levels = c("OA", "CTRL")
)
col_group <- group
levels(col_group) <- brewer.pal(
    nlevels(col_group),
    "Set1"
)
col_group <- as.character(col_group)
plotMDS(
    Expr_GSE82107_anno,
    labels = group,
    col = col_group,
)
title(main = "Sample groups MDS")
dev.off()
