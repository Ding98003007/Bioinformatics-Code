batch_cor_spearman <- function(gene) {
  y <- as.numeric(exprSet[gene, ]) # nolint
  rownames <- rownames(exprSet) # nolint
  do.call(rbind, future_lapply(rownames, function(x) { # nolint
    dd <- cor.test(as.numeric(exprSet[x, ]), y, type = "spearman") # nolint
    data.frame(gene = gene, mRNAs = x, cor = dd$estimate, p.value = dd$p.value)
  }))
}
batch_cor_pearson <- function(gene) {
  y <- as.numeric(exprSet[gene, ]) # nolint
  rownames <- rownames(exprSet) # nolint
  do.call(rbind, future_lapply(rownames, function(x) { # nolint
    dd <- cor.test(as.numeric(exprSet[x, ]), y, type = "pearson") # nolint
    data.frame(gene = gene, mRNAs = x, cor = dd$estimate, p.value = dd$p.value)
  }))
}