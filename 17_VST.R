vst <- function(object, blind = TRUE, nsub = 1000, fitType = "parametric") {
    if (nrow(object) < nsub) {
        stop("less than 'nsub' rows,
  it is recommended to use varianceStabilizingTransformation directly")
    }
    if (is.null(colnames(object))) {
        colnames(object) <- seq_len(ncol(object))
    }
    if (is.matrix(object)) {
        matrixIn <- TRUE
        object <- DESeqDataSetFromMatrix(object, DataFrame(row.names = colnames(object)), ~1)
    } else {
        if (blind) {
            design(object) <- ~1
        }
        matrixIn <- FALSE
    }
    if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
        object <- estimateSizeFactors(object)
    }
    baseMean <- rowMeans(counts(object, normalized = TRUE))
    if (sum(baseMean > 5) < nsub) {
        stop("less than 'nsub' rows with mean normalized count > 5,
  it is recommended to use varianceStabilizingTransformation directly")
    }

    # subset to a specified number of genes with mean normalized count > 5
    object.sub <- object[baseMean > 5, ]
    baseMean <- baseMean[baseMean > 5]
    o <- order(baseMean)
    idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
    object.sub <- object.sub[idx, ]

    # estimate dispersion trend
    object.sub <- estimateDispersionsGeneEst(object.sub, quiet = TRUE)
    object.sub <- estimateDispersionsFit(object.sub, fitType = fitType, quiet = TRUE)

    # assign to the full object
    suppressMessages({
        dispersionFunction(object) <- dispersionFunction(object.sub)
    })

    # calculate and apply the VST
    vsd <- varianceStabilizingTransformation(object, blind = FALSE)
    if (matrixIn) {
        return(assay(vsd))
    } else {
        return(vsd)
    }
}
