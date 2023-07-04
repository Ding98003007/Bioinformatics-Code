rm(list = ls())
library(tidyverse)
library(parallel)
source("code/functions/fpkm2tpm.R")
Targets <- read.table(
    file = "data/raw/GSE205431_group.txt",
    sep = "\t", check.names = T, header = T
) %>%
    arrange(sample)
gsm_id <- Targets %>%
    select(sample, sample_id)
count <- read.table(
    file = "data/raw/GSE205431_count.txt",
    sep = "\t", check.names = T, header = T, row.names = 1
) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    arrange(sample) %>%
    mutate(sample = str_remove(
        sample,
        pattern = "X"
    )) %>%
    left_join(gsm_id, by = "sample") %>%
    relocate(sample_id) %>%
    select(-sample) %>%
    column_to_rownames(var = "sample_id") %>%
    t() %>%
    as.data.frame()
identical(
    colnames(count),
    Targets$sample_id
)
genome_ref <- rtracklayer::import(
    "data/reference/Homo_sapiens.GRCh38.109.gtf.gz"
) %>%
    as.data.frame()
exon <- genome_ref[
    genome_ref$type == "gene",
    c("start", "end", "gene_id")
]
exon_by_geneid <- split(exon, exon$gene_id)
cl <- makeCluster(0.75 * detectCores())
efflen <- parLapply(cl, exon_by_geneid, function(x) {
    tmp <- apply(x, 1, function(y) {
        y[1]:y[2]
    }) 
    length(unique(unlist(tmp)))
})
gene_length <- data.frame(
    geneid = names(efflen),
    efflen = as.numeric(efflen)
)
expr1 <- count / gene_length$efflen
fpkm <- as.data.frame(t(t(expr1) / colSums(count)) * 10^9)
tpm <- as.data.frame(apply(fpkm, 2, fpkmToTpm)) 
colSums(tpm)
id2symbol <- genome_ref[
    genome_ref$type == "gene",
    c("gene_id", "gene_name")
] %>%
    distinct(gene_id, .keep_all = TRUE)
gene_type <- genome_ref[
    genome_ref$type == "gene",
    c("gene_id", "gene_biotype")
] %>%
    distinct(gene_id, .keep_all = TRUE)
mRNAs_tmp <- gene_type[gene_type$gene_biotype == "protein_coding", ]
mRNAs <- mRNAs_tmp$gene_id
count_anno <- count[mRNAs, ] %>%
    rownames_to_column(var = "gene_id") %>%
    mutate(gene_id = str_remove(
        gene_id,
        pattern = "\\..*"
    )) %>%
    left_join(id2symbol, by = "gene_id") %>%
    relocate(gene_name) %>%
    select(-gene_id) %>%
    mutate(sum = rowSums(.[, -1])) %>%
    arrange(desc(sum)) %>%
    select(-sum) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    na.omit()
write.table(count_anno,
    file = "data/tmp/tmp.txt",
    sep = "\t", row.names = FALSE, quote = F
)
count_anno <- read.table(
    file = "data/tmp/tmp.txt",
    sep = "\t", check.names = T, header = T, row.names = 1
)
fpkm_anno <- fpkm[mRNAs, ] %>%
    rownames_to_column(var = "gene_id") %>%
    mutate(gene_id = str_remove(
        gene_id,
        pattern = "\\..*"
    )) %>%
    left_join(id2symbol, by = "gene_id") %>%
    relocate(gene_name) %>%
    select(-gene_id) %>%
    mutate(sum = rowSums(.[, -1])) %>%
    arrange(desc(sum)) %>%
    select(-sum) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    na.omit()
write.table(fpkm_anno,
    file = "data/tmp/tmp.txt",
    sep = "\t", row.names = FALSE, quote = F
)
fpkm_anno <- read.table(
    file = "data/tmp/tmp.txt",
    sep = "\t", check.names = T, header = T, row.names = 1
)
tpm_anno <- tpm[mRNAs, ] %>%
    rownames_to_column(var = "gene_id") %>%
    mutate(gene_id = str_remove(
        gene_id,
        pattern = "\\..*"
    )) %>%
    left_join(id2symbol, by = "gene_id") %>%
    relocate(gene_name) %>%
    select(-gene_id) %>%
    mutate(sum = rowSums(.[, -1])) %>%
    arrange(desc(sum)) %>%
    select(-sum) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    na.omit()
write.table(tpm_anno,
    file = "data/tmp/tmp.txt",
    sep = "\t", row.names = FALSE, quote = F
)
tpm_anno <- read.table(
    file = "data/tmp/tmp.txt",
    sep = "\t", check.names = T, header = T, row.names = 1
)
Targets <- Targets %>%
    column_to_rownames(var = "sample_id") %>%
    select(group)
save(
    count_anno,
    fpkm_anno,
    tpm_anno,
    count,
    fpkm,
    tpm,
    Targets,
    file = "data/rdata/matrix_GSE205431_raw.Rdata"
)
write.table(
    count_anno,
    file = "res/data/matrix_GSE205431_count.txt",
    sep = "\t", quote = FALSE, row.names = TRUE
)
write.table(
    fpkm_anno,
    file = "res/data/matrix_GSE205431_fpkm.txt",
    sep = "\t", quote = FALSE, row.names = TRUE
)
