computePValue <- function(m, lower=250, tol=0.5, npts=10000, BPPARAM=SerialParam()) 
# A function to compute a non-ambient p-value for each barcode.
# 
# written by Aaron Lun
# created 1 August 2017
# last modified 12 August 2017
{
    discard <- rowSums(m) == 0
    m <- m[!discard,]

    # Computing the average profile.
    umi.sum <- colSums(m)
    ambient <- (umi.sum <= lower) # lower => "T" in the text.
    ambient.cells <- m[,ambient]
    ambient.prof <- rowSums(ambient.cells)
    ambient.prop <- goodTuringProportions(ambient.prof)

    # Removing cells from the ambient.
    obs <- m[,!ambient]
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

    # Computing a simulation with "npts" entries for any "tol"-fold interval around the interrogation point..
    lower.pt <- log2(min(obs.totals))-tol
    upper.pt <- log2(max(obs.totals))+tol
    S <- round(npts * (upper.pt - lower.pt)/(2*tol)) # npts => R in the text.
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
    stats <- .compute_P(obs.totals, obs.spread, sim.totals, sim.spread, tol)
    limited <- stats$limited
    p <- stats$p.value

    # Collating into some sensible output.
    all.p <- all.lr <- all.exp <- rep(NA_real_, length(ambient))
    all.lim <- rep(NA, length(ambient))
    all.p[!ambient] <- p
    all.lr[!ambient] <- obs.LR
    all.exp[!ambient] <- expected.LR
    all.lim[!ambient] <- limited
    return(data.frame(Total=umi.sum, LR=all.lr, Expected=all.exp, PValue=all.p, 
                      Limited=all.lim, row.names=colnames(m)))
}

.simulate_dev <- function(totals, prop) {
    .Call("calculate_random_dev", totals, prop)
}

.compute_P <- function(obs.totals, obs.spread, sim.totals, sim.spread, tol) {
    o <- order(obs.totals, obs.spread)
    perm.stats <- .Call("calculate_pval", obs.totals[o], obs.spread[o], sim.totals, sim.spread, 2^-tol)
    perm.stats[[1]][o] <- perm.stats[[1]]
    perm.stats[[2]][o] <- perm.stats[[2]]
    limited <- perm.stats[[1]]==0L
    p <- (perm.stats[[1]]+1)/(perm.stats[[2]]+1)
    return(list(limited=limited, p.value=p))
}

findInflectionPoint <- function(m, lower=250) 
# A function to identify the inflection point from a reverse cumulative distribution.
# 
# written by Aaron Lun
# created 7 August 2017    
{
    totals <- colSums(m)
    totals <- totals[totals >= lower]
    
    stuff <- rle(sort(totals, decreasing=TRUE))
    y <- log(stuff$values)
    x <- log(cumsum(stuff$lengths))
    
    grad <- diff(y)/diff(x)
    threshold <- y[which.min(grad)+1]
    return(unname(exp(threshold)))
}

detectCells <- function(m, lower=250, scale=2, ...) 
# Combined function that puts these all together, always keeping cells above the inflection
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    stats <- computePValue(m, lower=lower, ...)
    inflection <- findInflectionPoint(m, lower=lower)
    always <- stats$Total >= inflection*scale
    tmp <- stats$PValue
    tmp[always] <- 0
    stats$FDR <- p.adjust(tmp, method="BH")
    return(stats)
}
