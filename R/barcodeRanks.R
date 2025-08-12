#' Calculate barcode ranks
#'
#' Compute barcode rank statistics and identify the knee and inflection points on the total count curve.
#' 
#' @param m A numeric matrix-like object containing UMI counts, where columns represent barcoded droplets and rows represent genes.
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param lower A numeric scalar specifying the lower bound on the total UMI count, 
#' at or below which all barcodes are assumed to correspond to empty droplets and excluded from knee/inflection point identification.
#' @param exclude.from An integer scalar specifying the number of highest ranking barcodes to exclude from knee/inflection point identification.
#' @param fit.bounds,df Deprecated and ignored.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' 
#' @details
#' Analyses of droplet-based scRNA-seq data often show a plot of the log-total count against the log-rank of each barcode
#' where the highest ranks have the largest totals.
#' This is equivalent to a transposed empirical cumulative density plot with log-transformed axes, 
#' which focuses on the barcodes with the largest counts.
#' To help create this plot, the \code{barcodeRanks} function will compute these ranks for all barcodes in \code{m}.
#' Barcodes with the same total count receive the same average rank to avoid problems with discrete runs of the same total.
#' 
#' The function will also identify the inflection and knee points on the curve for downstream use.
#' Both of these points correspond to a sharp transition between two components of the total count distribution, 
#' presumably reflecting the difference between empty droplets with little RNA and cell-containing droplets with much more RNA.
#' \itemize{
#' \item The inflection point is defined as the point on the log-rank/log-total curve where the first derivative is minimized.
#' If multiple inflection points are present, we choose the point that immediately follows the knee point.
#' \item To find the knee point, we draw a diagonal line that passes through the inflection point in the log-rank/log-total curve.
#' The knee point is defined as the location on the curve that is above and most distant from this line.
#' }
#' Only points with total counts above \code{lower} will be considered for knee/inflection point identification.
#' Similarly, the first \code{exclude.from} points will be ignored to avoid instability at the start of the curve.
#' 
#' @return
#' A \linkS4class{DataFrame} where each row corresponds to a column of \code{m}, and containing the following fields:
#' \describe{
#' \item{\code{rank}:}{Numeric, the rank of each barcode (averaged across ties).}
#' \item{\code{total}:}{Numeric, the total counts for each barcode.}
#' }
#' 
#' The metadata contains \code{knee}, a numeric scalar containing the total count at the knee point;
#' and \code{inflection}, a numeric scalar containing the total count at the inflection point.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' # Mocking up some data: 
#' set.seed(2000)
#' my.counts <- DropletUtils:::simCounts()
#' 
#' # Computing barcode rank statistics:
#' br.out <- barcodeRanks(my.counts)
#' names(br.out)
#' 
#' # Making a plot.
#' plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
#' o <- order(br.out$rank)
#' abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
#' abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
#' legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#'     legend=c("knee", "inflection"))
#' 
#' @seealso
#' \code{\link{emptyDrops}}, where this function is used.
#'
#' @export
#' @name barcodeRanks
NULL

#' @importFrom utils head
#' @importFrom Matrix colSums
#' @importFrom S4Vectors DataFrame metadata<-
.barcode_ranks <- function(m, lower=100, exclude.from=50, fit.bounds=NULL, df=20, ..., BPPARAM=SerialParam()) {
    old <- .parallelize(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    totals <- unname(.intColSums(m))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    keep[run.rank <= exclude.from] <- FALSE
    if (sum(keep) < 2L) { 
        stop("insufficient unique points for computing knee/inflection points")
    } 

    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])
    deriv <- diff(y) / diff(x)

    # Initial inflection point is defined as the minima in the first derivative.
    infl.index <- which.min(deriv)

    # Heuristically drawing a diagonal line (gradient -1) from the initial inflection point.
    # The knee is defined as the point on the curve with the maximum distance from that line.
    # The -1 is more or less pulled out of thin air based on what most curves look like;
    if (infl.index > 1) {
        infl.x <- x[infl.index]
        infl.y <- y[infl.index] 
        left.of.infl.x <- head(x, infl.index) # only considering points to the left of the inflection.
        left.of.infl.y <- head(y, infl.index)

        .find_knee <- function(gradient) {
            intercept <- infl.y - gradient * infl.x
            relative.dist <- left.of.infl.y - gradient * left.of.infl.x - intercept # vertical vs perpendicular distance is the same, relatively.
            knee.index <- which.max(relative.dist)
            if (relative.dist[knee.index] <= 0) { # if it's not above the line, we failed to find the knee.
                NULL
            } else {
                knee.index
            }
        }

        knee.index <- .find_knee(-1)

        # If there's nothing above the line with a fixed gradient, we fall back to an empirical gradient from the start of the curve.
        # This is more sensitive to the number of real cells, which stretches out the plateau and causes a leftward shift in the knee point.
        # But, at least we'll get something approximating a knee point.
        if (is.null(knee.index)) {
            gradient <- (infl.y - y[1]) / (infl.x - x[1])
            knee.index <- .find_knee(gradient)

            # If there's still nothing, we just set the knee index to the inflection point.
            if (is.null(knee.index)) {
                knee.index <- infl.index
            }
        }
    }

    # Refining the inflection point to the interval immediately following the knee point. 
    # This aims to protect against curves with multiple inflection points.
    up.to <- findInterval(x[knee.index] + 1, x)
    new.infl.index <- knee.index + which.min(deriv[knee.index:up.to]) - 1L
    infl.index <- new.infl.index

    knee <- 10^y[knee.index]
    inflection <- 10^y[infl.index]

    out <- DataFrame(
        rank=.reorder(run.rank, stuff$lengths, o), 
        total=.reorder(run.totals, stuff$lengths, o)
    )
    rownames(out) <- colnames(m)
    metadata(out) <- list(knee=knee, inflection=inflection)
    out
}

.reorder <- function(vals, lens, o) {
    out <- rep(vals, lens)
    out[o] <- out
    return(out)
}

#' @export
#' @rdname barcodeRanks
setGeneric("barcodeRanks", function(m, ...) standardGeneric("barcodeRanks"))

#' @export
#' @rdname barcodeRanks
setMethod("barcodeRanks", "ANY", .barcode_ranks)

#' @export
#' @rdname barcodeRanks
#' @importFrom SummarizedExperiment assay
setMethod("barcodeRanks", "SummarizedExperiment", function(m, ..., assay.type="counts") {
    .barcode_ranks(assay(m, assay.type), ...)
})
