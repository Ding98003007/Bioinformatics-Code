rm(list = ls())
library(tidyverse)
library(RColorBrewer)
load(file = "data/rdata/matrix_GSE82107_raw.Rdata")
Expr_gse82107 <- Expr_GSE82107_anno
colors <- c(rep("#ff4400d1", 2), rep("#6200ff", 2))
pdf(
    file = "res/pic/final/Boxplots_expression_distributions_gse82107.pdf", # nolint
    width = 16,
    height = 10
)
par(mfrow = c(1, 1))
boxplot(Expr_gse82107,
    las = 2,
    col = colors,
    ylab = "Transformed Expression",
    xlab = "Distribution of Transformed Data"
)
title(
    main = "Distribution of Transformed Data of GSE82107",
    ylab = "Transformed Expression"
)
dev.off()
load(file = "data/rdata/matrix_GSE205431_raw.Rdata")
Expr_gse205431 <- tpm_anno
colors <- c(rep("#ff4400d1", 3), rep("#6200ff", 3))
pdf(
    file = "res/pic/final/Boxplots_expression_distributions_gse205431.pdf", # nolint
    width = 25,
    height = 12
)
par(mfrow = c(1, 1))
boxplot(Expr_gse205431,
    las = 2,
    col = colors,
    ylab = "Transformed Expression",
    xlab = "Distribution of Transformed Data"
)
title(
    main = "Distribution of Transformed Data of GSE205431",
    ylab = "Transformed Expression"
)
dev.off()
