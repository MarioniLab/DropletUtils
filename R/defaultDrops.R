#' Call cells from their total number of UMIs
#'
#' Call cells according to the number of UMIs associated with each barcode, as implemented in CellRanger version 2.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param expected A numeric scalar specifying the expected number of cells in this sample, as specified in the call to CellRanger.
#' @param upper.quant A numeric scalar between 0 and 1 specifying the quantile of the top \code{expected} barcodes to consider for the first step of the algorithm.
#' @param lower.prop A numeric scalar between 0 and 1 specifying the fraction of molecules of the \code{upper.quant} quantile result that a barcode must exceed to be called as a cell.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' 
#' @details
#' The \code{defaultDrops} function will call cells based on library size similarly to the CellRanger software suite from 10X Genomics.
#' Default arguments correspond to an exact reproduction of CellRanger's algorithm, where the number of expected cells was also left at CellRanger default value.
#' 
#' The method computes the \code{upper.quant} quantile of the top \code{expected} barcodes, ordered by decreasing number of UMIs.
#' Any barcodes containing more molecules than \code{lower.prop} times this quantile is considered to be a cell, and is retained for further analysis.
#' 
#' This method may be vulnerable to calling very well-captured background RNA as cells, or missing real cells with low RNA content.
#' See \code{?\link{emptyDrops}} for an alternative approach for cell calling.
#' 
#' @return
#' A logical vector of length \code{ncol(m)}, indicating whether each column of \code{m} was called as a cell.
#' 
#' @author
#' Jonathan Griffiths
#' 
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts()
#' 
#' # Identify likely cell-containing droplets.
#' called <- defaultDrops(my.counts)
#' table(called)
#' 
#' # Get matrix of called cells.
#' cell.counts <- my.counts[, called]
#' 
#' @references
#' 10X Genomics (2017).
#' Cell Ranger Algorithms Overview.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview}
#' 
#' @seealso
#' \code{\link{emptyDrops}}, for another method for calling cells.
#'
#' @name defaultDrops
NULL

#' @importFrom stats quantile
#' @importFrom utils head
.default_drops <- function(m, expected = 3000, upper.quant = 0.99, lower.prop = 0.1) {
    if(upper.quant > 1 | upper.quant < 0){
        stop("'upper.quant' should be a numeric value between 0 and 1")
    }
    
    if(lower.prop > 1 | lower.prop < 0){
        stop("'lower.prop' should be a numeric value between 0 and 1")
    }
    
    libs <- colSums(m)
    o <- order(libs, decreasing = TRUE)
    top <- libs[head(o, n = expected)]
    
    threshold <- quantile(top, upper.quant)*lower.prop
    libs > threshold
}

#' @export
#' @rdname defaultDrops
setGeneric("defaultDrops", function(m, ...) standardGeneric("defaultDrops"))

#' @export
#' @rdname defaultDrops
setMethod("defaultDrops", "ANY", .default_drops)

#' @export
#' @rdname defaultDrops
#' @importFrom SummarizedExperiment assay
setMethod("defaultDrops", "SummarizedExperiment", function(m, ..., assay.type="counts") {
    .default_drops(assay(m, assay.type), ...)
})
