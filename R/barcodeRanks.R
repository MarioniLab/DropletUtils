#' Calculate barcode ranks
#'
#' Compute barcode rank statistics and identify the knee and inflection points on the total count curve.
#' 
#' @param m A numeric matrix-like object containing UMI counts, where columns represent barcoded droplets and rows represent genes.
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param lower A numeric scalar specifying the lower bound on the total UMI count, 
#' at or below which all barcodes are assumed to correspond to empty droplets and excluded from knee/inflection point identification.
#' @param exclude.from An integer scalar specifying the number of highest ranking barcodes to exclude from knee/inflection point identification.
#' @param window Numeric scalar specifying the length of the window (in log10 units) for knee/inflection point identification.
#' Larger values improve stability of estimates at the cost of sensitivity to changes in the curve.
#' @param gradient.threshold Numeric scalar specifying the maximum threshold on the gradient for identifying potential elbow points.
#' Lower values increase the stringency of elbow point identification.
#' @param fit.bounds,df Deprecated and ignored.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' 
#' @details
#' Analyses of droplet-based scRNA-seq data often show a plot of the log-total count against the log-rank of each barcode where the highest ranks have the largest totals.
#' This is equivalent to a transposed empirical cumulative density plot with log-transformed axes, which focuses on the barcodes with the largest counts.
#' To create this plot, the \code{barcodeRanks} function will compute these ranks for all barcodes in \code{m}.
#' Barcodes with the same total count receive the same average rank to avoid problems with discrete runs of the same total.
#' 
#' The function will also identify the inflection and knee points on the curve for downstream use.
#' Both of these points correspond to a sharp transition between two components of the total count distribution, 
#' presumably reflecting the difference between empty droplets with little RNA and cell-containing droplets with much more RNA.
#' Only points with total counts above \code{lower} will be considered for knee/inflection point identification.
#' Similarly, the first \code{exclude.from} points will be ignored to avoid instability at the start of the curve.
#'
#' The actual identification of the knee/inflection points is based on a simple curve-tracing algorithm.
#' We trace a window of fixed length \code{window} through the curve, and for each window, we consider the straight line connecting its ends: 
#' \itemize{
#' \item To find the knee, we filter for windows where the midpoint of the window lies above the end-connecting line.
#' Of these, we select the window with the shortest end-connecting line, i.e., the strongest curvature.
#' The midpoint of that window is defined as the knee.
#' \item To find the inflection, we pick the window with the lowest (i.e., most negative) gradient of the end-connecting line.
#' The midpoint of that window is defined as the inflection.
#' \item In cases with multiple knee/inflection points, we aim to report the earlier values, i.e., those with higher log-totals.
#' This is achieved by ignoring all windows after the first one that contains an \dQuote{elbow} point in the curve.
#' A window contains an elbow if its midpoint lies below the end-connecting line and the gradient is less than \code{gradient.threshold}.
#' }
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
.barcode_ranks <- function(m, lower=100, exclude.from=50, window=1, gradient.threshold=-1, fit.bounds=NULL, df=20, ..., BPPARAM=SerialParam()) {
    old <- .parallelize(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    totals <- unname(.intColSums(m))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    keep[run.rank <= exclude.from] <- FALSE
    if (sum(keep) < 2) {
        stop("insufficient unique points for computing knee/inflection points")
    }
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])

    # Scanning a window of length 'window' along the curve.
    # The length is calculated along the curve itself, not along the axes.
    dist.along.curve <- c(0, sqrt(diff(x)^2 + diff(y)^2))
    cumdist <- cumsum(dist.along.curve)
    rhs.loc <- cumdist + window
    to.scan <- rhs.loc <= cumdist[length(cumdist)]

    if (!any(to.scan)) {
        knee <- inflection <- 10^y[length(y)] # just pick the last point, whatever.
    } else {
        left.x <- x[to.scan]
        left.y <- y[to.scan]

        right.info <- .interpolate_on_curve(rhs.loc[to.scan], cumdist, dist.along.curve, x, y)
        right.x <- right.info$x
        right.y <- right.info$y

        window.gap <- sqrt((left.x - right.x)^2 + (left.y - right.y)^2) # distance in 2D space between the ends of the window.
        window.gradient <- (right.y - left.y) / (right.x - left.x) # gradient of the line between ends of the window
        window.intercept <- right.y - window.gradient * right.x # intercept of the line between ends of the window

        mid.info <- .interpolate_on_curve(cumdist[to.scan] + window/2, cumdist, dist.along.curve, x, y)
        mid.x <- mid.info$x
        mid.y <- mid.info$y

        mid.above <- mid.y > window.gradient * mid.x + window.intercept

        # We define the knee point at the window in 'mid.above' with the smallest 'window.gap', and the inflection point at the window with the most negative 'window.gradient'.
        # In practice, some curves contain multiple knees and inflections, and we would like to restrict ourselves to the first "meaningful" minima.
        # We do so by ignoring all points after the first "elbow point", i.e., window gradient below some threshold and the curve is below the line.
        elbow.index <- which(!mid.above & window.gradient < gradient.threshold)
        if (length(elbow.index) == 0) {
            infl.window <- which.min(window.gradient)
            maybe.knee <- which(mid.above)
        } else {
            first.elbow <- elbow.index[1]
            infl.window <- which.min(head(window.gradient, first.elbow))
            maybe.knee <- which(head(mid.above, first.elbow))
        }

        knee.window <- maybe.knee[which.min(window.gap[maybe.knee])]
        if (length(knee.window) == 0) {
            # Fallback if the curve is so weird that the midpoint is never above the window line.
            knee.window <- infl.window
        }

        # Picking an actual inflection/knee point based on the midpoint of the window.
        knee <- 10^mid.y[knee.window]
        inflection <- 10^mid.y[infl.window]
    }

    out <- DataFrame(
        rank=.reorder(run.rank, stuff$lengths, o), 
        total=.reorder(run.totals, stuff$lengths, o)
    )
    rownames(out) <- colnames(m)
    metadata(out) <- list(knee=knee, inflection=inflection)
    out
}

.interpolate_on_curve <- function(target.distance, cumulative.distance, stepwise.distance, x, y) {
    index <- findInterval(target.distance, cumulative.distance, checkSorted=FALSE, checkNA=FALSE, left.open=TRUE)
    index.p1 <- index + 1L
    prop <- (target.distance - cumulative.distance[index]) / stepwise.distance[index.p1]
    list(
        x=.interpolate_simple(prop, x[index], x[index.p1]),
        y=.interpolate_simple(prop, y[index], y[index.p1])
    )
}

.interpolate_simple <- function(prop, left.val, right.val) {
    left.val + prop * (right.val - left.val)
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
