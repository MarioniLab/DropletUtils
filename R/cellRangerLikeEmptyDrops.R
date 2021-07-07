#' An approximate implementation of the `--soloCellFilter  EmptyDrops_CR` filtering approach to identify empty droplets.
#'
#' An approximate implementation of the `--soloCellFilter  EmptyDrops_CR` filtering approach, 
#' which, itself, was reverse-engineered from the behavior of CellRanger 3+.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for cells.
#' 
#' @param umiMin A numeric scalar specifying the  minimum UMI count above which a sample will be included in ambient profiles, 
#' as specified in the call to CellRanger.
#' 
#' @param umiMinFracMedian A numeric scalar between 0 and 1 specifying that only the samples whose UMI count are above \code{umiMinFracMedian} 
#' fraction of the median UMI count#'  of the top \code{nExpectedCells} samples will be included in the ambient profile. 
#' as specified in the call to CellRanger.
#' 
#' @param candMaxN An integer specifying the maximum number of ambient samples that are possible to be regarded as real cells, 
#' as specified in the call to CellRanger.
#' 
#' @param indMax An integer specifying the highest UMI count ranking of the ambient pool, cells with UMI count ranking above
#' this number will not be included in the ambient pool, as specified in the call to CellRanger.
#' 
#' @param indMin An integer specifying the lowest UMI count ranking of the ambient pool, cells with UMI count ranking below
#' this number will not be included in the ambient pool, as specified in the call to CellRanger.
#' 
#' @param fdr_thresh A numeric scalar specifying the FDR threshold to filter samples. Samples whose FDR returned by emptyDrops
#' is above this threshold will not be regarded as real cells, as specified in the call to CellRanger.
#' 
#' @param maxPercentile A numeric scalar specifying the percentile used in simple filtering, samples selected by simple filtering 
#' will be regarded as real cells regardless of the \code{emptyDrops} result, as specified in the call to CellRanger.
#' 
#' @param nExpectedCells A numeric scalar specifying the expected number of cells in this sample, as specified in the call to CellRanger.
#' 
#' @param maxMinRatio A numeric scalar specifying the maximum ratio of maximum UMI count and minimum UMI count used in simple filtering, 
#' maximum UMI count used in simple filtering is determined first by \code{nExpectedCells*(1-maxPercentile)}, minimum UMI count used in
#'  simple filtering is then determined by this ratio, as specified in the call to CellRanger..
#'  
#' @param seed Integer specifying the seed that will be used to run \code{emptyDrops}
#' @param ... For the generic, further arguments to pass to \code{emptyDrops}.
#'
#' @details
#' This function is an approximate implementation of the  `--soloCellFilter  EmptyDrops_CR` filtering approach of STARsolo 
#' (\url{https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1}), which, itself, was reverse engineered from the behavior of  CellRanger 3+. 
#' The original C++ code on which this function is based can be found at 
#' (\url{https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp}) 
#' All parameters are defaulty set as the default value used in starSolo and Cellranger.
#' In the most cases, users just need to specify the raw and unfiltered count matrix, \code{m}.
#' See \code{?\link{emptyDrops}} for an alternative approach for cell calling.
#' 
#' @return
#' A DataFrame like \code{\link{emptyDrops}}, with an additional binary \code{is.cell} field demonstrating whether
#' samples are estimated as real cells.
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
#' e.out <- cellRangerLikeEmptyDrops(my.counts)
#' e.out
#' 
#' # Get matrix of estimated cells.
#' cell.counts <- my.counts[, e.out$is.cell]
#' 
#' @references
#' Kaminow et al. (2021).
#' STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data
#' \url{https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1}
#' 
#' @seealso
#' \code{\link{emptyDrops}}, for another method for calling cells.
#'
#' @name cellRangerLikeEmptyDrops
NULL

# Authors: Dongze He, Rob Patro
# Center of Bioinformatics and Computational Biology, University of Maryland, College Park, Maryland, 20740

.cellRangerLikeEmptyDrops  <- function(m, 
                                      umiMin=500,
                                      umiMinFracMedian=0.01, 
                                      candMaxN=20000, 
                                      indMax=90000, 
                                      indMin=45000, 
                                      fdr_thresh=0.01, 
                                      maxPercentile=0.99, 
                                      nExpectedCells=3000, 
                                      maxMinRatio=10,
                                      seed=2718,
                                      ...
) {
  
  # This function is an approximate implementation of the 
  # `--soloCellFilter  EmptyDrops_CR` filtering approach 
  # of STARsolo (https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1),
  # which, itself, was reverse engineered from the behavior of 
  # CellRanger 3+. The original C++ code on which this 
  # function is based can be found at (https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp) 
  
  ###################################################################################################################    
  # get the sorted nUMI vector of cells 
  csums <- colSums2(m)
  indCount <- as.data.frame(cbind(1:length(csums), csums))
  colnames(indCount) <- c("index", "count")
  indCount <- indCount[order(indCount$count,decreasing = TRUE),]
  
  # Simple Filtering
  maxind <- round(nExpectedCells * (1 - maxPercentile))
  nUMImax <- indCount$count[min(ncol(m), maxind)]
  nUMImin <- round(nUMImax/maxMinRatio)
  ncellsSimple <- sum(indCount$count>=nUMImin)
  
  # set lower bound
  minUMI    <- max(umiMin, round(umiMinFracMedian * indCount$count[ncellsSimple/2]))
  
  ## we at most assign candMaxN samples in the ambient pool as real cells
  minUMI <- max(minUMI, indCount$count[min(ncellsSimple+candMaxN,nrow(indCount))])
  
  # emptyDrops
  ## ignore: the lower bound of UMI count, samples with UMI count less than ignore
  ## will not be considered as ambient cells.
  ignore_index <- min(ncol(m), indMax)
  ignore <- indCount$count[ignore_index]
  
  ## by.rank: cells with UMI count ranking lower than by.rank will be considered as 
  ## ambient cells
  by.rank <- indMin
  
  ## retain: samples with UMI count higher than retain will be regarded as cells 
  retain <- indCount$count[ncellsSimple]
  
  ## the cells with total UMI count between ignore and lower will be considered as ambient
  set.seed(seed)
  e.out <- emptyDrops(m, by.rank=by.rank, ignore=ignore, retain=retain, alpha=Inf)
  e.out$is.cell <- e.out$FDR < fdr_thresh
  e.out$is.cell[is.na(e.out$is.cell)] <- FALSE
  
  # further filter cells by minUMI
  e.out$is.cell[indCount[indCount$count<minUMI, "index"]] <- FALSE
  e.out
}

#' @export
#' @rdname cellRangerLikeEmptyDrops
setGeneric("cellRangerLikeEmptyDrops", function(m, ...) standardGeneric("cellRangerLikeEmptyDrops"))

#' @export
#' @rdname cellRangerLikeEmptyDrops
setMethod("cellRangerLikeEmptyDrops", "ANY", .cellRangerLikeEmptyDrops)

#' @export
#' @rdname cellRangerLikeEmptyDrops
#' @importFrom SummarizedExperiment assay
setMethod("cellRangerLikeEmptyDrops", "SummarizedExperiment", function(m, ..., assay.type="counts") {
  .cellRangerLikeEmptyDrops(assay(m, assay.type), ...)
})
