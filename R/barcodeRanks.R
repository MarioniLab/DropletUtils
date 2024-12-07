#' Calculate barcode ranks
#'
#' Compute barcode rank statistics and identify the knee and inflection points on the total count curve.
#' 
#' @param m A numeric matrix-like object containing UMI counts, where columns represent barcoded droplets and rows represent genes.
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param lower A numeric scalar specifying the lower bound on the total UMI count, 
#' at or below which all barcodes are assumed to correspond to empty droplets.
#' @param fit.bounds A numeric vector of length 2, specifying the lower and upper bounds on the total UMI count
#' from which to obtain a section of the curve for spline fitting.
#' @param exclude.from An integer scalar specifying the number of highest ranking barcodes to exclude from spline fitting.
#' Ignored if \code{fit.bounds} is specified.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param df Deprecated and ignored.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the ANY method, further arguments to pass to \code{\link{smooth.spline}}.
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
#' The function will also identify the inflection and knee points on the curve for downstream use, 
#' Both of these points correspond to a sharp transition between two components of the total count distribution, 
#' presumably reflecting the difference between empty droplets with little RNA and cell-containing droplets with much more RNA.
#' \itemize{
#' \item The inflection point is computed as the point on the rank/total curve where the first derivative is minimized.
#' The derivative is computed directly from all points on the curve with total counts greater than \code{lower}.
#' This avoids issues with erratic behaviour of the curve at lower totals.
#' \item The knee point is defined as the point on the curve that is furthest from the straight line drawn between the \code{fit.bounds} locations on the curve.
#' We used to minimize the signed curvature to identify the knee point but this relies on the second derivative,
#' which was too unstable even after smoothing.
#' }
#'
#' If \code{fit.bounds} is not specified, the lower bound is automatically set to the inflection point
#' as this should lie below the knee point on typical curves.
#' The upper bound is set to the point at which the first derivative is closest to zero, 
#' i.e., the \dQuote{plateau} region before the knee point.
#' The first \code{exclude.from} barcodes with the highest totals are ignored in this process 
#' to avoid spuriously large numerical derivatives from unstable parts of the curve with low point density.
#'
#' Note that only points with total counts above \code{lower} will be considered for curve fitting,
#' regardless of how \code{fit.bounds} is defined.
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

#' @importFrom utils tail
#' @importFrom Matrix colSums
#' @importFrom S4Vectors DataFrame metadata<-
.barcode_ranks <- function(m, lower=100, fit.bounds=NULL, exclude.from=50, df=20, ..., BPPARAM=SerialParam()) {
    old <- .parallelize(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    totals <- unname(.intColSums(m))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    if (sum(keep)<3) { 
        stop("insufficient unique points for computing knee/inflection points")
    } 
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])
    
    # Numerical differentiation to identify bounds for spline fitting.
    edge.out <- .find_curve_bounds(x=x, y=y, exclude.from=exclude.from) 
    left.edge <- edge.out["left"]
    right.edge <- edge.out["right"]

    # As an aside: taking the right edge to get the total for the inflection point.
    # We use the numerical derivative as the spline is optimized for the knee.
    inflection <- 10^(y[right.edge])

    # We restrict curve fitting to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    if (is.null(fit.bounds)) {
        new.keep <- left.edge:right.edge
    } else {
        new.keep <- which(y > log10(fit.bounds[1]) & y < log10(fit.bounds[2]))
    }

    # Using the maximum distance to identify the knee point.
    if (length(new.keep) >= 4) {
        curx <- x[new.keep]
        cury <- y[new.keep]
        xbounds <- curx[c(1L, length(new.keep))]
        ybounds <- cury[c(1L, length(new.keep))]
        gradient <- diff(ybounds)/diff(xbounds)
        intercept <- ybounds[1] - xbounds[1] * gradient
        above <- which(cury >= curx * gradient + intercept)
        dist <- abs(gradient * curx[above] - cury[above] + intercept)/sqrt(gradient^2 + 1)
        knee <- 10^(cury[above[which.max(dist)]])
    } else {
        # Sane fallback upon overly aggressive filtering by 'exclude.from', 'lower'.
        knee <- 10^(y[new.keep[1]]) 
    }

    # Returning a whole stack of useful stats.
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

.find_curve_bounds <- function(x, y, exclude.from) 
# The upper/lower bounds are defined at the plateau and inflection, respectively.
# Some exclusion of the LHS points avoids problems with discreteness.
{
    d1n <- diff(y)/diff(x)

    skip <- min(length(d1n) - 1, sum(x <= log10(exclude.from)))
    d1n <- tail(d1n, length(d1n) - skip)

    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])

    c(left=left.edge, right=right.edge) + skip
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
