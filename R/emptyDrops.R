testEmptyDrops <- function(m, lower=100, span=sqrt(2), npts=10000, test.ambient=FALSE, BPPARAM=SerialParam()) 
# A function to compute a non-ambient p-value for each barcode.
# 
# written by Aaron Lun
# created 1 August 2017
{
    discard <- rowSums(m) == 0
    m <- m[!discard,,drop=FALSE]
    ncells <- ncol(m)

    # Computing the average profile from the ambient cells.
    umi.sum <- colSums(m)
    ambient <- umi.sum <= lower # lower => "T" in the text.
    ambient.cells <- m[,ambient]
    ambient.prof <- rowSums(ambient.cells)
    ambient.prop <- goodTuringProportions(ambient.prof)

    # Removing supposed ambient cells from the matrix.
    if (!test.ambient) {
        keep <- !ambient
    } else {
        keep <- umi.sum > 0
    }
    obs <- m[,keep,drop=FALSE]
    obs.totals <- colSums(obs)

    # Calculating the likelihood ratio.
    if (is(obs, "dgCMatrix")) {
        i <- obs@i + 1L
        j <- rep(seq_len(ncol(obs)), diff(obs@p)) 
        x <- obs@x
    } else if (is(obs, "dgTMatrix")) {
        i <- obs@i + 1L
        j <- obs@j + 1L
        x <- obs@x
    } else {
        stop("unsupported matrix type")
    }
    
    p.n0 <- x * log(x/(ambient.prop[i]*obs.totals[j])) 
    by.col <- aggregate(p.n0, list(Col=j), sum)
    obs.LR <- numeric(length(obs.totals))
    obs.LR[by.col$Col] <- by.col$x

    # Computing a simulation with "npts" entries for any "span"-fold interval around the interrogation point..
    span <- log2(span)
    lower.pt <- log2(min(obs.totals))-span
    upper.pt <- log2(max(obs.totals))+span
    S <- round(npts * (upper.pt - lower.pt)/(2*span)) # npts => R in the text.
    sim.totals <- 2^seq(from=lower.pt, to=upper.pt, length.out=S)

    # Computing the deviance estimate for simulated runs.
    sim.LR <- bpvec(sim.totals, FUN=.simulate_dev, prop=ambient.prop, BPPARAM=BPPARAM)

    # Modelling the total-dependent trend in the simulated LR (and the variance around the trend).
    log.totals <- log2(sim.totals)
    trend.fit <- lowess(x=log.totals, y=sim.LR, f=0.2)
    sim.spread <- sim.LR/trend.fit$y
    expected.LR <- spline(trend.fit$x, trend.fit$y, xout=log2(obs.totals))$y  
    obs.spread <- obs.LR/expected.LR

    # Computing a p-value for each observed value.
    stats <- .compute_P(obs.totals, obs.spread, sim.totals, sim.spread, span)
    limited <- stats$limited
    p <- stats$p.value

    # Collating into some sensible output.
    all.p <- all.lr <- all.exp <- rep(NA_real_, ncells)
    all.lim <- rep(NA, ncells)
    all.p[keep] <- p
    all.lr[keep] <- obs.LR
    all.exp[keep] <- expected.LR
    all.lim[keep] <- limited
    return(DataFrame(Total=umi.sum, Deviance=all.lr, Expected=all.exp, PValue=all.p, Limited=all.lim, row.names=colnames(m)))
}

.simulate_dev <- function(totals, prop) {
    .Call(cxx_calculate_random_dev, totals, prop)
}

.compute_P <- function(obs.totals, obs.spread, sim.totals, sim.spread, span) {
    o <- order(obs.totals, obs.spread)
    perm.stats <- .Call(cxx_calculate_pval, obs.totals[o], obs.spread[o], sim.totals, sim.spread, 2^-span)
    perm.stats[[1]][o] <- perm.stats[[1]]
    perm.stats[[2]][o] <- perm.stats[[2]]
    limited <- perm.stats[[1]]==0L
    p <- (perm.stats[[1]]+1)/(perm.stats[[2]]+1)
    return(list(limited=limited, p.value=p))
}

findKneePoint <- function(m, lower=100) 
# A function to identify the knee point from a reverse cumulative distribution.
# 
# written by Aaron Lun
# created 7 August 2017    
{
    totals <- colSums(m)
    totals <- totals[totals >= lower]
    totals <- sort(totals, decreasing=TRUE)

    stuff <- rle(totals)
    y <- log(stuff$values)
    x <- log(cumsum(stuff$lengths) - (stuff$lengths-1)/2) # Get mid-rank of each run.

    # Smoothing to avoid error multiplication upon differentiation.
    fit <- smooth.spline(x, y)
    d1 <- predict(fit, deriv=1)$y
    d2 <- predict(fit, deriv=2)$y

    # Maximizing the curvature and returning the total at which this occurs.
    curvature <- abs(d2)/(1 + d1^2)^1.5
    return(exp(y[which.max(curvature)]))
}

emptyDrops <- function(m, lower=100, scale=1, ...) 
# Combined function that puts these all together, always keeping cells above the inflection
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    stats <- testEmptyDrops(m, lower=lower, ...)
    kneept <- findKneePoint(m, lower=lower)
    always <- stats$Total >= kneept*scale
    tmp <- stats$PValue
    tmp[always] <- 0
    stats$FDR <- p.adjust(tmp, method="BH")
    return(stats)
}
