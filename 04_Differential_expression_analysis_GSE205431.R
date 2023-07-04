rm(list = ls())
library(tidyverse)
library(DESeq2)
library(future.apply)
library(clusterProfiler)
library(GseaVis)
plan(multicore)
load(file = "data/rdata/matrix_GSE205431_raw.Rdata")
exprSet <- count_anno %>%
    rownames_to_column(var = "gene_id") %>%
    mutate(sum = rowSums(.[, -1])) %>%
    arrange(desc(sum)) %>%
    dplyr::select(-sum) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")
group <- factor(Targets$group,
    levels = c("SX", "CTRL")
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
    contrast = c("group", "SX", "CTRL"),
    name = "group_SX_vs_CTRL",
    tidy = TRUE
)
geneList <- DEGs_ebayes$log2FoldChange 
names(geneList) <- DEGs_ebayes$row
geneList <- sort(geneList,
    decreasing = TRUE
) %>%
    na.omit()
d <- "data/MsigDB/symbols/h.all.gmt"
geneset <- read.gmt(d)
egmt <- GSEA(geneList,
    eps = 0,
    TERM2GENE = geneset,
    minGSSize = 10,
    maxGSSize = 2000,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE,
)
res_hallmark <- as.data.frame(egmt)
pdf(
    file = "res/pic/final/dot_plot_class_gse205431.pdf",
    width = 14,
    height = 8
)
dotplotGsea(
    data = res_hallmark,
    topn = 10,
    # order.by = "p.adjust",
    # add.seg = TRUE,
    # scales = "free"
)
dev.off()
phenotypes <- c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
    # "HALLMARK_ANGIOGENESIS"
)
pdf(
    file = "res/pic/final/gsea_plot_phenotypes_gse205431.pdf", # nolint
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
    file = "data/rdata/gse205431_GSEA_res.Rdata"
)
save(
    DEGs_ebayes,
    file = "data/rdata/deg_GSE205431_deg_oavsctrl.Rdata"
)
write.table(
    DEGs_ebayes,
    file = "res/data/genes/gse205431_deg_sig.txt",
    sep = "\t", quote = FALSE, row.names = TRUE
)
