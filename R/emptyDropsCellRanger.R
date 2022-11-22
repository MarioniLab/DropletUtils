#' CellRanger's emptyDrops variant 
#'
#' An approximate implementation of the \code{--soloCellFilter EmptyDrops_CR} filtering approach, 
#' which itself was reverse-engineered from the behavior of CellRanger 3.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for barcodes.
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such an object.
#' @param n.expected.cells An integer scalar specifying the number of expected cells in a sample. 
#' Corresponds to the \code{nExpectedCells} argument in \pkg{STARsolo}. 
#' @param max.percentile A numeric scalar between 0 and 1 used to define the maximum UMI count in the simple filtering algorithm. 
#' Corresponds to the \code{maxPercentile} argument in \pkg{STARsolo}. 
#' @param max.min.ratio An integer scalar specifying the ratio of the maximum and minimum UMI count in the simple filtering algorithm. 
#' Corresponds to the \code{maxMinRatio} argument in \pkg{STARsolo}.
#' @param umi.min An integer scalar specifying the minimum UMI count for inclusion of a barcode in the cell candidate pool. 
#' Corresponds to the \code{umiMin} argument in \pkg{STARsolo}.
#' @param umi.min.frac.median A numeric scalar between 0 and 1 used to define the minimum UMI count for inclusion of a barcode in the cell candidate pool.
#' Specifically, the minimum is defined as \code{umi.min.frac.median} times the median UMI count of the real cells assigned by the simple filtering algorithm. 
#' Corresponds to the \code{umiMinFracMedian} argument in \pkg{STARsolo}.
#' @param cand.max.n An integer scalar specifying the maximum number of barcodes that can be included in the cell candidate pool. 
#' In effect, this applies a minimum threshold that is defined as the \code{cand.max.n}-th largest UMI count among all cells that are \emph{not} selected by the simple filtering algorithm. 
#' Corresponds to the \code{candMaxN} in \pkg{STARsolo}.
#' @param ind.min An integer scalar specifying the lowest UMI count ranking for inclusion of a barcode in the ambient profile. 
#' Corresponds to the \code{indMin} argument in \pkg{STARsolo}.
#' @param ind.max An integer scalar specifying the highest UMI count ranking for inclusion of a barcode in the ambient profile. 
#' Corresponds to the \code{indMin} argument in \pkg{STARsolo}.
#' @param round A logical scalar indicating whether to check for non-integer values in \code{m} and, if present, round them for ambient profile estimation (see \code{?\link{ambientProfileEmpty}}) and the multinomial simulations.
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be used.
#' @param ... Further arguments to pass to individual methods.
#' Specifically, for the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.type String or integer specifying the assay of interest.
#'
#' @details
#' \code{emptyDropsCellRanger} splits each sample's barcodes into three subsets.
#' \enumerate{
#' \item The first subset contains barcodes that are selected by the \dQuote{simple filtering algorithm}, which are regarded as high quality cells without any further filtering.
#' The minimum threshold \eqn{T} for this subset is defined by taking the \code{max.percentile} percentile of the top \code{n.expected.cells} barcodes,
#' and then dividing by the \code{max.min.ratio} to obtain a minimum UMI count.
#' All barcodes here will have an FDR of zero.
#' \item The second subset contains the ambient pool and is defined as all barcodes with rankings between \code{ind.min} and \code{ind.max}. 
#' The barcodes that fall in this category will be used to compute the ambient profile.
#' None of these barcodes are considered to be potential cells.
#' \item The third subset contains the pool of barcodes that are potential cells, i.e., cell candidates.
#' This is defined as all barcodes with total counts below \eqn{T} and higher than all of the thresholds defined by \code{umi.min}, \code{umi.min.frac.median} and \code{cand.max.n}.
#' Only the barcodes within this subset will be tested for signficant deviations from the ambient profile, i.e., FDR is not \code{NaN}.
#' }
#'
#' As of time of writing, the arguments in \pkg{STARsolo} have a one-to-one correspondence with the arguments in \code{emptyDropsCellRanger}. 
#' All parameter defaults are set as the same as those used in STARsolo 2.7.9a.
#' 
#' The main differences between \code{emptyDropsCellRanger} and \code{emptyDrops} are:
#' \itemize{
#' \item \code{emptyDropsCellRanger} applies a simple filtering strategy to identify some real cells before any further investigation.
#' \item \code{emptyDropsCellRanger} takes barcodes whose total count ranks within a certain range - by default, \eqn{(45,000, 90,000]} - to compute the ambient profile.
#' In contrast, \code{emptyDrops} only defines the upper bound using \code{lower} or \code{by.rank}.
#' \item \code{emptyDropsCellRanger} defines a cell candidate pool according to three parameters, \code{umi.min},\code{umi.min.frac.median} and \code{cand.max.n}.
#' In \code{emptyDrops}, this is only defined by \code{lower}.
#' }
#' 
#' @return
#' A \linkS4class{DataFrame} with the same fields as that returned by \code{\link{emptyDrops}}.
#' 
#' @author
#' Dongze He, Rob Patro
#' 
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts(nempty=100000, nlarge=2000, nsmall=1000)
#' 
#' # Identify likely cell-containing droplets.
#' out <- emptyDropsCellRanger(my.counts)
#' out
#'
#' is.cell <- out$FDR <= 0.01
#' sum(is.cell, na.rm=TRUE)
#'
#' # Subsetting the matrix to the cell-containing droplets.
#' # (using 'which()' to handle NAs smoothly).
#' cell.counts <- my.counts[,which(is.cell),drop=FALSE]
#' dim(cell.counts)
#' 
#' @references
#' Kaminow B, Yunusov D, Dobin A (2021).
#' STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data.
#' \url{https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1}
#' 
#' @seealso
#' \code{\link{emptyDrops}}, for the original implementation.
#'
#' @name emptyDropsCellRanger
NULL

#' @importFrom stats p.adjust
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom BiocParallel SerialParam
.empty_drops_cell_ranger  <- function(m, 
                                      # STARsolo arguments
                                      ## simple filtering
                                      n.expected.cells=3000,            # nExpectedCells
                                      max.percentile=0.99,      # maxPercentile
                                      max.min.ratio=10,         # maxMinRatio
                                      
                                      ## emptyDrops_CR
                                      umi.min=500,              # umiMin
                                      umi.min.frac.median=0.01, # umiMinFracMedian
                                      cand.max.n=20000,         # candMaxN
                                      ind.min=45000,            # indMin
                                      ind.max=90000,            # indMax

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
        o <- order(totals, decreasing=TRUE)

        # Simple Filtering
        # https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
        # line 36-61
        max.ind <- round(n.expected.cells * (1 - max.percentile)) # maxind
        n.umi.max <- totals[o[min(length(totals), max.ind)]] # nUMImax
        
        # Barcodes with UMI count higher than retain will be regarded as real
        # barcodes without any further tests
        retain <- max(round(n.umi.max/max.min.ratio),1) # nUMImin
        
        # select barcodes to use as ambient solution
        # SoloFeature_emptyDrops_CR.cpp line 117-134
        ncells.simple <- sum(totals >= retain)
        min.umi  <- max(umi.min, round(umi.min.frac.median * totals[o[ncells.simple/2]]))
        i.cand.last <- min(ncells.simple + cand.max.n, sum(totals > min.umi))

        # NOTE: parallelization handled by setAutoBPPARAM above.
        discard <- rowSums(mat) == 0
        if (any(discard)) {
            mat <- mat[!discard,,drop=FALSE]
        }
        
        ambient <- logical(length(totals))
        ambient[o[min(ind.min, length(totals)):min(ind.max, length(totals))]] <- TRUE
        ambient.m <- mat[,ambient,drop=FALSE]
        ambient.prof <- rowSums(ambient.m)
        
        if (sum(ambient.prof)==0) {
            stop("no counts available to estimate the ambient profile")
        }
        ambient.prop <- .safe_good_turing(ambient.prof)
        
        # Barcodes to keep are not just !ambient, as we want to exclude the
        # barcodes that are too low to be even considered as ambient.
        # keep <- (o < (ncells.simple+1))
        keep <- logical(length(totals)) 
        keep[o[(ncells.simple+1):i.cand.last]] = TRUE
        
        list(
            m=mat, # this MUST have the same number of columns as input.
            discard=discard,
            ambient=ambient,
            ambient.m=ambient.m,
            ambient.prop=ambient.prop,
            keep=keep,
            metadata=list(lower = totals[o[min(ind.min, length(totals))]], bottom=totals[o[min(ind.max, length(totals))]], retain=retain)
        )
    }
    
    stats <- .test_empty_drops(m=m, ambient.FUN=ambfun, niters=niters, test.ambient=FALSE, ignore=NULL, alpha=Inf, round=round, BPPARAM=BPPARAM)
    
    tmp <- stats$PValue
    retain <- metadata(stats)$retain
    always <- stats$Total >= retain
    tmp[always] <- 0
    
    stats$FDR <- p.adjust(tmp, method="BH")
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
