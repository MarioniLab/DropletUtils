#' CellRanger's emptyDrops variant 
#'
#' An approximate implementation of the \code{--soloCellFilter EmptyDrops_CR} filtering approach, 
#' which itself was reverse-engineered from the behavior of CellRanger 3.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for barcodes.
#' @param n.expected.cells An integer scalar specifying the number of expected cells in a sample. 
#' @param max.percentile A numeric scalar between 0 and 1 used to define the maximum UMI count used in the simple filtering algorithm. 
#' @param max.min.ratio An integer scalar specifying the ratio of the maximum and minimum UMI count used in the simple filtering algorithm. 
#' @param umi.min An integer scalar specifying the minimum UMI count for inclusion of a barcode in the cell candidate pool. 
#' @param umi.min.frac.median A numeric scalar between 0 and 1 used to define the minimum UMI count for inclusion of a barcode in the cell candidate pool.
#' @param cand.max.n An integer scalar specifying the maximum number of barcodes that can be included in the cell candidate pool. 
#' @param ind.min An integer scalar specifying the lowest UMI count ranking for inclusion of a barcode in the ambient profile. 
#' @param ind.max An integer scalar specifying the highest UMI count ranking for inclusion of a barcode in the ambient profile. 
#' @param round A logical scalar indicating whether to check for non-integer values in \code{m} and, if present, round them for ambient profile estimation (see \code{?\link{ambientProfileEmpty}}) and the multinomial simulations.
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be used.
#'
#' @section Details about \code{emptyDropsCellRanger} arguments:
#' The arguments in \pkg{STARsolo} are in one-to-one correspondence with the arguments in \code{emptyDropsCellRanger}. All parameters defaults are set as the same as those used in STARsolo 2.7.9a.
#' \itemize{
#' \item \code{n.expected.cells}: This argument is used to define the number of expected cells in a sample. 
#' It is the same as \code{nExpectedCells} in \pkg{STARsolo}. 
#' \item \code{max.percentile}: Together with \code{n.expected.cells}, \code{max.percentile} defines the maximum UMI count that an expected cell can have. 
#' This argument is the same as \code{maxPercentile} in \pkg{STARsolo}. 
#' \item \code{max.min.ratio}: Based on the maximum UMI count computed using \code{n.expected.cells} and \code{max.percentile}, \code{max.min.ratio} defines the threshold of the minimum UMI count that an expected cell can have according to their ratio. 
#' The barcodes whose UMI count is above this threshold are assigned as real cells by the simple filtering algorithm (\code{nCellsSimple} in \pkg{STARsolo}). 
#' It is same as \code{maxMinRatio} in \pkg{STARsolo}.
#' \item \code{umi.min}: This argument is one of the three thresholds used to define the lower bound of the cell candidate pool.
#' \code{umi.min} specifies the minimum UMI count for inclusion of a barcode in the cell candidate pool. 
#' This argument is the same as \code{umiMin} in \pkg{STARsolo}.
#' \item \code{umi.min.frac.median}: This argument is one of the three thresholds used to define the lower bound of the cell candidate pool. 
#' This argument defines the minimum UMI as \code{umi.min.frac.median} times the median UMI count of the real cells assigned by the simple filtering algorithm. 
#' This argument is the same as \code{umiMinFracMedian} in \pkg{STARsolo}.
#' \item \code{cand.max.n}: This argument is one of the three thresholds used to define the lower bound of the cell candidate pool. 
#' This argument defines the minimum UMI as the UMI count of the barcodes whose UMI count ranks as the \code{cand.max.n}-th behind the real cells defined by the simple filtering algorithm. 
#' This argument is the same as \code{umiMinFracMedian} in \pkg{STARsolo}.
#' \item \code{ind.min}: This argument specifies the lowest UMI count ranking of the barcodes in the ambient pool. 
#' This argument is same as \code{indMin} in \pkg{STARsolo}.
#' \item \code{ind.max}: This argument specifies the highest UMI count ranking of the barcodes in the ambient pool. 
#' This argument is same as \code{indMin} in \pkg{STARsolo}.
#' }
#' 
#' @section Details about \code{emptyDropsCellRanger}:
#' \code{emptyDropsCellRanger} works on only three subsets of the barcodes in a sample. 
#' The first subset includes barcodes that are selected by the simple filtering pipeline. 
#' These barcodes will be regarded as high quality cells without any further filtering.
#' The second subset defines the ambient pool. 
#' This subset is defined by \code{ind.min} and \code{ind.max}. 
#' The barcodes that fall in this category will be used to compute the ambient profile. 
#' Real cells are expected to have significantly non-ambient profile.
#' The third subset defines the cell candidate pool. 
#' This subset includes the barcodes whose number of UMI count is lower than that of the cells selected by the simple filtering algorithm and higher than the threshold defined by \code{umi.min}, \code{umi.min.frac.median} and \code{cand.max.n} together.
#' Only the barcodes within this subset will be considered further (as potential non-empty barcodes) to compute the deviations from the ambient profile, and then the p-value and FDR. 
#' In other words, besides barcodes that are assigned as cells directly by the simple filtering algorithm, only the barcodes that fall in this category have the chance be selected as real cells (FDR is not \code{NAN}.).
#' 
#' @section Differences between \code{emptyDropsCellRanger} and \code{emptyDrops}:
#' The main differences between \code{emptyDropsCellRanger} and \code{emptyDrops} are:
#' \itemize{
#' \item \code{emptyDropsCellRanger} applies a simple filtering strategy to identify some real cells before any further investigation.
#' \item \code{emptyDropsCellRanger} takes barcodes whose total count ranks within a certain range - by default, (45,000, 90,000] - to compute the ambient profile.
#' While \code{emptyDrops} only defines the upper bound using \code{lower} or \code{by.rank}.
#' \item \code{emptyDropsCellRanger} defines a cell candidate pool according to three parameters, \code{umi.min},\code{umi.min.frac.median} and \code{cand.max.n}. While in \code{emptyDrops}, this is defined by \code{lower}.
#' }
#' 
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
#' my.counts <- DropletUtils:::simCounts()
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
#' Kaminow et al. (2021).
#' STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data
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
                                      can.max.n=20000,          # candMaxN
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
        # Simple Filtering
        # https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
        # line 36-61
        max.ind <- round(n.expected.cells * (1 - max.percentile)) # maxind
        n.umi.max <- totals[order(totals, decreasing = TRUE)[min(length(totals), max.ind)]] # nUMImax
        
        # Barcodes with UMI count higher than retain will be regarded as real
        # barcodes without any further tests
        retain <- max(round(n.umi.max/max.min.ratio),1) # nUMImin
        
        # select barcodes to use as ambient solution
        # SoloFeature_emptyDrops_CR.cpp line 117-134
        ncells.simple <- sum(totals >= retain)
        min.umi  <- max(umi.min, round(umi.min.frac.median * totals[ncells.simple/2]))
        i.cand.last <- min(ncells.simple + can.max.n, sum(totals > min.umi))

        # NOTE: parallelization handled by setAutoBPPARAM above.
        discard <- rowSums(mat) == 0
        if (any(discard)) {
            mat <- mat[!discard,,drop=FALSE]
        }
        
        o <- order(totals, decreasing=TRUE)
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