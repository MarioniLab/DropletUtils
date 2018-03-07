#' @export
#' @importFrom stats smooth.spline predict fitted
barcodeRanks <- function(m, lower=100, fit.bounds=NULL, df=20, ...) 
# Returning statistics to construct a barcode-rank plot. Also calculates
# the knee and inflection points for further use.
#
# written by Aaron Lun
# created 22 December 2017    
{
    totals <- colSums(m)
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
    # The upper/lower bounds are defined at the plateau and inflection, respectively.
    d1n <- diff(y)/diff(x)
    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])

    # We restrict to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    if (is.null(fit.bounds)) {
        new.keep <- left.edge:right.edge
    } else {
        new.keep <- y > log10(fit.bounds[1]) & y < log10(fit.bounds[2])
    }

    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee point.
    fit <- smooth.spline(x[new.keep], y[new.keep], df=df, ...)
    d1 <- predict(fit, deriv=1)$y
    d2 <- predict(fit, deriv=2)$y
    curvature <- d2/(1 + d1^2)^1.5
    knee <- 10^(y[which.min(curvature)])

    # Taking the right edge to get the total for the inflection point.
    # We use the numerical derivative as the spline is optimized for the knee.
    inflection <- 10^(y[right.edge])

    # Returning a whole stack of useful stats.
    fitted.vals <- rep(NA_real_, length(keep))
    fitted.vals[keep][new.keep] <- 10^fitted(fit)
    return(list(rank=.reorder(run.rank, stuff$lengths, o), 
                total=.reorder(run.totals, stuff$lengths, o),
                fitted=.reorder(fitted.vals, stuff$lengths, o),
                knee=knee, inflection=inflection))
}

.reorder <- function(vals, lens, o) {
    out <- rep(vals, lens)
    out[o] <- out
    return(out)
}

