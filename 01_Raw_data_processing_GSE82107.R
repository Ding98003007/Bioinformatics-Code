rm(list = ls())
library(tidyverse)
library(limma)
library(GEOquery)
library(AnnoProbe)
gset <- getGEO(
    filename = "data/raw/GSE82107_series_matrix.txt.gz",
    getGPL = FALSE
)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex)
}
Targets_raw <- read.table(
    file = "data/raw/GSE82107_group.txt",
    sep = "\t", check.names = TRUE, header = TRUE
) %>%
    arrange(sample_id) %>%
    column_to_rownames(var = "sample_id")
Targets <- Targets_raw %>%
    select(group)
probes_expr <- exprs(gset)
dim(probes_expr)
head(probes_expr[, 1:4])
gpl <- "GPL96"
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene <- idmap(
    gpl,
    type = "bioc"
)
head(probe2gene)
Expr_GSE82107_anno <- filterEM(
    probes_expr,
    probe2gene
)
head(Expr_GSE82107_anno)
Expr_GSE82107_anno <- Expr_GSE82107_anno[
    ,
    rownames(Targets)
]
identical(
    colnames(Expr_GSE82107_anno),
    rownames(Targets)
)
save(
    Expr_GSE82107_anno,
    Targets,
    file = "data/rdata/matrix_GSE82107_raw.Rdata"
)
write.table(
    Expr_GSE82107_anno,
    file = "res/data/matrix/matrix_GSE82107.txt",
    sep = "\t", quote = FALSE, row.names = TRUE
)
