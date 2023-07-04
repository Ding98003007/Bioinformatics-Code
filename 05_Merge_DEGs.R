rm(list = ls())
library(tidyverse)
load("data/rdata/deg_GSE82107_deg_oavsctrl.Rdata")
DEG_gse82107 <- DEGs_ebayes
load("data/rdata/deg_GSE205431_deg_oavsctrl.Rdata")
DEG_gse205431 <- DEGs_ebayes
p_filter <- 0.05
logFC_filter <- 0.5
DEG_gse82107_sig <- DEG_gse82107 %>%
    filter(P.Value < p_filter) %>%
    filter(abs(logFC) > logFC_filter) %>%
    # filter(logFC > logFC_filter) %>%
    rownames_to_column(var = "SYMBOL") %>%
    dplyr::select(SYMBOL, logFC) %>%
    arrange(desc(logFC)) %>%
    na.omit()
write.table(
    DEG_gse82107_sig,
    file = "res/data/genes/gse82107_deg_sig.txt",
    sep = "\t", quote = FALSE, row.names = TRUE
)
DEG_gse205431_sig <- DEG_gse205431 %>%
    filter(pvalue < p_filter) %>%
    filter(abs(log2FoldChange) > logFC_filter) %>%
    # filter(logFC > logFC_filter) %>%
    mutate(
        SYMBOL = row,
        logFC = log2FoldChange
    ) %>%
    dplyr::select(SYMBOL, logFC) %>%
    arrange(desc(logFC)) %>%
    na.omit()
write.table(
    DEG_gse205431_sig,
    file = "res/data/genes/gse205431_deg_sig.txt",
    sep = "\t", quote = FALSE, row.names = TRUE
)
COM_GENES <- intersect(DEG_gse82107_sig$SYMBOL, DEG_gse205431_sig$SYMBOL)
OA_DEG_COM <- DEG_gse82107_sig[
    which(DEG_gse82107_sig$SYMBOL %in% COM_GENES),
] %>%
    arrange(desc(logFC)) %>%
    mutate(logFC_OA = logFC) %>%
    dplyr::select(SYMBOL, logFC_OA) %>%
    left_join(
        DEG_gse205431_sig,
        by = "SYMBOL"
    ) %>%
    mutate(logFC_SX = logFC) %>%
    dplyr::select(SYMBOL, logFC_OA, logFC_SX)
write.table(OA_DEG_COM,
    file = "res/data/genes/com_genes_for_tf.txt",
    sep = "\t", row.names = FALSE
)
TF <- read.table(
    file = "res/data/genes/KnockTF-Analysis.txt",
    sep = "\t", check.names = TRUE, header = TRUE
) %>%
    arrange(Rank)
TF_OA <- TF[
    which(TF$TF %in% DEG_gse82107_sig$SYMBOL),
]
Targets_TF <- "EGR1"
Targets_TF_logfc <- DEG_gse82107_sig[
    which(DEG_gse82107_sig$SYMBOL == Targets_TF),
] %>%
    arrange(desc(logFC)) %>%
    mutate(logFC_OA = logFC) %>%
    dplyr::select(SYMBOL, logFC_OA) %>%
    left_join(
        DEG_gse205431_sig,
        by = "SYMBOL"
    ) %>%
    mutate(logFC_SX = logFC) %>%
    dplyr::select(SYMBOL, logFC_OA, logFC_SX)
write.table(Targets_TF_logfc,
    file = "res/data/genes/Targets_TF_logfc.txt",
    sep = "\t", row.names = FALSE
)
Targets_TF_genes <- as.character(str_split(
    TF[
        which(TF$TF == Targets_TF),
    ]$Overlapping_Genes,
    ",",
    simplify = T
) %>%
    t())
Targets_TF_genes_OA <- DEG_gse82107_sig[
    which(DEG_gse82107_sig$SYMBOL %in% Targets_TF_genes),
] %>%
    arrange(desc(logFC)) %>%
    mutate(logFC_OA = logFC) %>%
    dplyr::select(SYMBOL, logFC_OA) %>%
    left_join(
        DEG_gse205431_sig,
        by = "SYMBOL"
    ) %>%
    mutate(logFC_SX = logFC) %>%
    dplyr::select(SYMBOL, logFC_OA, logFC_SX) %>%
    mutate(
        sig = ifelse(logFC_OA > 0 & logFC_SX < 0,
            "diff_direction",
            ifelse(logFC_OA < 0 & logFC_SX > 0,
                "diff_direction",
                "same_direction"
            )
        )
    )
write.table(Targets_TF_genes_OA,
    file = "res/data/genes/Targets_TF_genes_OA.txt",
    sep = "\t", row.names = FALSE
)
