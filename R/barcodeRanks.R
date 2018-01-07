barcodeRanks <- function(m, lower=200, df=50, ...) 
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
    fit <- smooth.spline(x, y, df=df, ...)
    d1 <- predict(fit, deriv=1)$y
    d2 <- predict(fit, deriv=2)$y

    # Minimizing the signed curvature and returning the total for the knee point.
    curvature <- d2/(1 + d1^2)^1.5
    knee <- 10^(y[which.min(curvature)])

    # Minimizing the first derivative to get the total for the inflection point.
    inflection <- 10^(y[which.min(d1)])

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

