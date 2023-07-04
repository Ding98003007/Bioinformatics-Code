rm(list = ls())
library(tidyverse)
library(limma)
library(future.apply)
library(clusterProfiler)
library(GseaVis)
plan(multicore)
load(file = "data/rdata/matrix_GSE82107_raw.Rdata")
group <- factor(Targets$group,
    levels = c("OA", "CTRL")
)
Expr <- Expr_GSE82107_anno
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(Expr)
design
contrast_matrix <- makeContrasts(
    OAvsCTRL = OA - CTRL,
    levels = colnames(design)
)
contrast_matrix
vfit <- lmFit(
    Expr,
    design
)
vfit <- contrasts.fit(
    vfit,
    contrasts = contrast_matrix
)
efit <- eBayes(vfit)
fit_ebayes <- decideTests(efit)
summary(fit_ebayes)
DEGs_ebayes <- topTable(
    efit,
    coef = 1,
    n = Inf
)
tfit <- treat(vfit)
fit_treat <- decideTests(tfit)
summary(fit_treat)
DEGs_treat <- topTreat(
    tfit,
    coef = 1,
    n = Inf
)
DEGs_ebayes_sig <- as.data.frame(DEGs_ebayes) %>%
    filter(adj.P.Val < 0.05)
DEGs_treat_sig <- as.data.frame(DEGs_treat) %>%
    filter(adj.P.Val < 0.05)

geneList <- DEGs_ebayes$logFC # nolint
names(geneList) <- rownames(DEGs_ebayes)
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
    file = "res/pic/final/dot_plot_class_GSE82107.pdf", # nolint
    width = 14,
    height = 10
)
dotplotGsea(
    data = res_hallmark,
    topn = 20,
    # order.by = "NES",
    # add.seg = TRUE,
    scales = "free"
)
dev.off()
phenotypes <- c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
    # "HALLMARK_ANGIOGENESIS"
)
pdf(
    file = "res/pic/final/gsea_plot_phenotypes_GSE82107.pdf", # nolint
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
    file = "data/rdata/GSE82107_GSEA_res.Rdata"
)
save(
    DEGs_ebayes,
    DEGs_ebayes_sig,
    DEGs_treat,
    DEGs_treat_sig,
    file = "data/rdata/deg_GSE82107_deg_oavsctrl.Rdata"
)
save(
    Expr,
    group,
    file = "data/tmp/heatmap_data_GSE82107_deg_oavsctrl.Rdata"
)
