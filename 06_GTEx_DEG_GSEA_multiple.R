rm(list = ls())
source("code/functions/correlation.R")
library(clusterProfiler)
library(DESeq2)
library(tidyverse)
library(future.apply)
library(ggplot2)
library(GSEABase)
library(GseaVis)
library(pheatmap)
plan(multicore)
Targets_gene <- c("EGR1")
exprSet <- read.table(
    file = "data/raw/GTEx_muscle_count.txt",
    sep = "\t", check.names = T, header = T
) %>%
    mutate(sum = rowSums(.[, -1])) %>%
    arrange(desc(sum)) %>%
    dplyr::select(-sum) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    column_to_rownames(var = "SYMBOL")
med <- median(as.numeric(exprSet[Targets_gene, ]))
group_list_expr <- as.data.frame(t(
    ifelse(exprSet[Targets_gene, ] > med, "High", "Low")
)) %>%
    rownames_to_column(var = "sample_id") %>%
    mutate(EXPLEVEL = EGR1) %>%
    dplyr::select(-EGR1)
group <- factor(group_list_expr$EXPLEVEL,
    levels = c("High", "Low")
)
conditions <- data.frame(
    sample = colnames(exprSet),
    group = group
) %>%
    column_to_rownames("sample")
dds <- DESeqDataSetFromMatrix(
    countData = exprSet,
    colData = conditions,
    design = ~group
)
dds_deg <- DESeq(dds, parallel = TRUE)
resultsNames(dds_deg)
DEGs_ebayes <- results(dds_deg,
    contrast = c("group", "High", "Low"),
    name = "group_High_vs_Low",
    tidy = TRUE
)
res_all <- DEGs_ebayes %>%
    # filter(padj < 0.05) %>%
    # filter(abs(log2FoldChange) > 0.5) %>%
    dplyr::select(row, log2FoldChange) %>%
    arrange(desc(log2FoldChange)) %>%
    na.omit()
geneList <- res_all$log2FoldChange # nolint
names(geneList) <- res_all$row
geneList <- sort(geneList,
    decreasing = TRUE
) %>%
    na.omit()
geneList <- geneList
d <- "data/MsigDB/symbols/h.all.gmt"
geneset <- read.gmt(d)
egmt <- GSEA(geneList,
    eps = 0,
    TERM2GENE = geneset,
    minGSSize = 100,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE,
)
res_hallmark <- as.data.frame(egmt)
pdf(
    file = "res/pic/final/dot_plot_target_gene_deg_gtex.pdf", # nolint
    width = 14,
    height = 8
)
dotplotGsea(
    data = res_hallmark,
    topn = 10,
    # order.by = "NES",
    add.seg = TRUE,
    scales = 'free'
)
dev.off()
phenotypes <- c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
    # "HALLMARK_ANGIOGENESIS"
)
pdf(
    file = "res/pic/final/gsea_plot_phenotypes_gtex.pdf", # nolint
    width = 10,
    height = 8
)
gseaNb(
    object = egmt,
    geneSetID = phenotypes,
    subPlot = 2,
    termWidth = 35,
    legend.position = c(0.8, 0.8),
    # addGene = genes,
    addPval = TRUE,
    pvalX = 0.05,
    pvalY = 0.05,
    # rmHt = FALSE,
    curveCol = jjAnno::useMyCol("stallion2", 3)
)
dev.off()
save(
    egmt,
    res_hallmark,
    phenotypes,
    file = "data/rdata/gtex_GSEA_res.Rdata"
)
