barcodeRanks <- function(m, lower=100, ...) 
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
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep]) 
    
    # Smoothing to avoid error multiplication upon differentiation.
    fit <- smooth.spline(x, y, ...)
    d1 <- predict(fit, deriv=1)$y
    d2 <- predict(fit, deriv=2)$y

    # Maximizing the curvature and returning the total for the knee point.
    curvature <- abs(d2)/(1 + d1^2)^1.5
    knee <- 10^(y[which.max(curvature)])

    # Maximizing the second derivative to get the total for the inflection point.
    inflection <- 10^(y[which.max(d2)])

    # Returning a whole stack of useful stats.
    fitted.vals <- rep(NA_real_, length(keep))
    fitted.vals[keep] <- 10^fitted(fit)
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

