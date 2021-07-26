#' An approximate implementation of the \code{--soloCellFilter  EmptyDrops_CR} filtering approach to identify empty droplets.
#'
#' An approximate implementation of the \code{--soloCellFilter  EmptyDrops_CR} filtering approach, 
#' which, itself, was reverse-engineered from the behavior of CellRanger 3+.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for barcodes.
#' @param expected A numeric scalar specifying the expected number of barcodes in this sample It is same as \code{nExpectedCells} in STARsolo.
#' @param max.percentile A numeric scalar specifying the percentile used in simple filtering, barcodes selected by simple filtering 
#' These will be regarded as real barcodes regardless of the \code{emptyDrops} result. 
#' It is same as \code{maxPercentile} in STARsolo.
#' @param max.min.ratio A numeric scalar specifying the maximum ratio of maximum UMI count and minimum UMI count used in simple filtering, 
#' maximum UMI count used in simple filtering is determined first by \code{expected*(1-max.percentile)}, minimum UMI count used in
#'  simple filtering is then determined by this ratio. It is same as \code{maxMinRatio} in STARsolo.
#' @param umi.min A numeric scalar specifying the minimum UMI count above which a sample will be included in ambient profiles. It is same as \code{umiMin} in STARsolo.
#' @param umi.min.frac.median A numeric scalar between 0 and 1 specifying that only the barcodes whose UMI count are above this number 
#' fraction of the median UMI count of the top \code{expected} barcodes will be included in the ambient profile. 
#' It is same as \code{umiMinFracMedian} in STARsolo.
#' @param cand.max.n An integer specifying the maximum number of ambient barcodes that are possible to be regarded as real barcodes. It is same as \code{canMaxN} in STARsolo.
#' @param ind.min An integer specifying the lowest UMI count ranking of the ambient pool, barcodes with UMI count ranking below
#' this number will not be included in the ambient pool. It is same as \code{indMin} in STARsolo. It is also same as \code{by.rank} in \code{emptyDrops}.
#' @param ind.max An integer specifying the highest UMI count ranking of the ambient pool, barcodes with UMI count ranking above
#' this number will not be included in the ambient pool. It is same as \code{indMax} in STARsolo.
#' @param fdr A numeric scalar specifying the FDR threshold to filter barcodes. Barcodes whose FDR returned by emptyDrops
#' is above this threshold will not be regarded as real barcodes. It is same as \code{FDR} in STARsolo.
#' @param round Logical scalar indicating whether to check for non-integer values in \code{m} and, if present, round them for ambient profile estimation (see \code{?\link{ambientProfileEmpty}}) and the multinomial simulations.
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be used.
#'
#' @details
#' This function is an approximate implementation of the  \code{--soloCellFilter  EmptyDrops_CR} filtering approach of STARsolo 2.7.9a
#' , which, itself, was reverse-engineered from the behavior of  CellRanger 3+. 
#' All parameters are default set as the default value used in starSolo 2.7.9a.
#' In most cases, users just need to specify the raw and unfiltered count matrix, \code{m}.
#' See \code{?\link{emptyDrops}} for an alternative approach for cell calling.
#' 
#' The main differences between \code{emptyDropsCellRanger} and \code{emptyDrops} are 
#' 1. \code{emptyDropsCellRanger} first applies a simple filtering strategy to identify 
#' retained cells according to the ranking of the total count of barcodes. This process is based on
#'  \code{expected}, \code{max.percentile}, \code{max.min.ratio}, \code{umi.min}, and \code{umi.min.frac.median}. 
#' 2. \code{emptyDropsCellRanger} takes barcodes whose total count rank within a certain range 
#' (by default, (45,000, 90,000]) as the input of SimpleGoodTuring. So in addition to the lower 
#' limit \code{lower}, it also has an bottom limit \code{bottom}.
#' 3. When computing ambient profile, \code{emptyDropsCellRanger} defines a candidate pool. 
#' Only the barcodes in the pool are involved in ambient profile computation and are assigned a p-value. 
#' By default, the pool include the 20,000 barcodes whose total count rank right after the barcodes 
#' selected by the simple filtering strategy described above.
#'  
#' @return
#' A DataFrame like \code{\link{emptyDrops}}, with an additional binary \code{is.cell} field demonstrating whether
#' barcodes are estimated as real barcodes.
#' 
#' @author
#' Dongze He, Rob Patro
#' 
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts()
#' 
#' # Identify likely cell-containing droplets.
#' e.out <- emptyDropsCellRanger(my.counts)
#' e.out
#' 
#' # Get matrix of estimated barcodes.
#' cell.counts <- my.counts[, e.out$is.cell]
#' 
#' @references
#' Kaminow et al. (2021).
#' STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data
#' \url{https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1}
#' 
#' @seealso
#' \code{\link{emptyDrops}}, for another method for calling barcodes.
#'
#' @name emptyDropsCellRanger
NULL

#' @export
#' @rdname emptyDropsCellRanger

#' @importFrom stats p.adjust
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom BiocParallel SerialParam
.empty_drops_cell_ranger  <- function(m, 
    # STARsolo arguments
    ## simple filtering
    expected=3000,            # nExpectedCells
    max.percentile=0.99,      # maxPercentile
    max.min.ratio=10,         # maxMinRatio

    ## emptyDrops_CR
    umi.min=500,              # umiMin
    umi.min.frac.median=0.01, # umiMinFracMedian
    can.max.n=20000,          # candMaxN
    ind.min=45000,            # indMin
    ind.max=90000,            # indMax
    fdr=0.01,                 # FDR

    # emptyDrops arguments
    round=TRUE,
    niters=10000,
    BPPARAM=SerialParam())
{  
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }
    
    # This function is an approximate implementation of the 
    # `--soloCellFilter  EmptyDrops_CR` filtering approach 
    # of STARsolo 2.7.9a (https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1),
    # which, itself, was reverse engineered from the behavior of 
    # CellRanger 3+. The original C++ code on which this 
    # function is based can be found at (https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
    #  and https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_emptyDrops_CR.cpp) 
    
    ambfun <- function(mat, totals) {
        # Simple Filtering
        # https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
        # line 36-61
        max.ind <- round(expected * (1 - max.percentile)) # maxind
        n.umi.max <- totals[order(totals, decreasing = TRUE)[min(length(totals), max.ind)]] # nUMImax

        # Barcodes with UMI count higher than retain will be regarded as real
        # barcodes without any further tests
        retain <- max(round(n.umi.max/max.min.ratio),1) # nUMImin

        # select barcodes to use as ambient solution
        # SoloFeature_emptyDrops_CR.cpp line 117-134
        ncells.simple <- sum(totals >= retain)
        min.umi  <- max(umi.min, round(umi.min.frac.median * totals[ncells.simple/2]))
        i.cand.last <- min(ncells.simple + can.max.n, sum(totals > min.umi))
        max.rank <- min(ind.max + 1, length(totals), i.cand.last)

        # NOTE: parallelization handled by setAutoBPPARAM above.
        discard <- rowSums(m) == 0
        if (any(discard)) {
            m <- m[!discard,,drop=FALSE]
        }

        o <- order(totals, decreasing=TRUE)
        ambient <- logical(length(totals))
        ambient[o[ind.min:max.rank]] <- TRUE
        ambient.m <- m[,ambient,drop=FALSE]
        ambient.prof <- rowSums(ambient.m)

        if (sum(ambient.prof)==0) {
            stop("no counts available to estimate the ambient profile")
        }
        ambient.prop <- .safe_good_turing(ambient.prof)

        # Barcodes to keep are not just !ambient, as we want to exclude the
        # barcodes that are too low to be even considered as ambient.
        keep <- (o < ind.min) 

        list(
            m=m, # this MUST have the same number of columns as input.
            discard=discard,
            ambient=ambient,
            ambient.m=ambient.m,
            ambient.prop=ambient.prop,
            keep=keep,
            metadata=list(lower = totals[o[ind.min]], bottom=totals[o[max.rank]], retain=retain)
        )
    }

    stats <- .test_empty_drops(m=m, ambient.FUN=ambfun, niters=niters, test.ambient=FALSE, ignore=NULL, alpha=Inf, round=round, BPPARAM=BPPARAM)

    tmp <- stats$PValue
    retain <- metadata(stats)$retain
    always <- stats$Total >= retain
    tmp[always] <- 0
    
    stats$FDR <- p.adjust(tmp, method="BH")
    stats$is.cell <- stats$FDR <= fdr
    stats$is.cell[is.na(stats$is.cell)] <- FALSE
    stats
}

#' @export
#' @rdname emptyDropsCellRanger
setGeneric("emptyDropsCellRanger", function(m, ...) standardGeneric("emptyDropsCellRanger"))

#' @export
#' @rdname emptyDropsCellRanger
setMethod("emptyDropsCellRanger", "ANY", .empty_drops_cell_ranger)

#' @export
#' @rdname emptyDropsCellRanger
#' @importFrom SummarizedExperiment assay
setMethod("emptyDropsCellRanger", "SummarizedExperiment", function(m, ..., assay.type="counts") {
    .empty_drops_cell_ranger(assay(m, assay.type), ...)
})
