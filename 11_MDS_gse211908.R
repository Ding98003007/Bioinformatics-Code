rm(list = ls())
library(tidyverse)
library(limma)
library(RColorBrewer)
load(file = "data/rdata/matrix_GSE205431_raw.Rdata")
group <- factor(Targets$group,
    levels = c("SX", "CTRL")
)
pdf(
    file = "res/pic/final/MDSplots_gse205431.pdf", # nolint
    width = 10,
    height = 8
)
par(mfrow = c(1, 1))
col_group <- group
levels(col_group) <- brewer.pal(
    nlevels(col_group),
    "Set1"
)
col_group <- as.character(col_group)
plotMDS(
    count,
    labels = group,
    col = col_group,
)
title(main = "Sample groups MDS")
dev.off()
