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
    umi.sum <- as.integer(round(colSums(m)))
    ambient <- umi.sum <= lower # lower => "T" in the text.
    ambient.cells <- m[,ambient]
    ambient.prof <- rowSums(ambient.cells)
    ambient.prop <- goodTuringProportions(ambient.prof)

    # Removing supposed ambient cells from the matrix.
    if (!test.ambient) {
        keep <- !ambient
    } else {
        keep <- umi.sum > 0L
    }
    obs <- m[,keep,drop=FALSE]
    obs.totals <- umi.sum[keep]

    # Calculating the log-multinomial probability for each cell.
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
   
    p.n0 <- x * log(ambient.prop[i]) - lfactorial(x)
    by.col <- aggregate(p.n0, list(Col=j), sum)
    obs.P <- numeric(length(obs.totals))
    obs.P[by.col$Col] <- by.col$x

    # Calculating the p-values using a Monte Carlo approach.
    o <- order(obs.totals, obs.P)
    re.P <- obs.P[o]
    re.totals <- rle(obs.totals[o]) # Ensure identity upon comparison.

    nworkers <- bpworkers(BPPARAM)
    per.core <- rep(ceiling(npts/nworkers), nworkers)
    per.core[1] <- npts - sum(per.core[-1]) # Making sure that we get exactly 'npts' iterations.

    out.values <- bplapply(per.core, FUN=.monte_carlo_pval, total.val=re.totals$values,
                           total.len=re.totals$lengths, P=re.P, ambient=ambient.prop, BPPARAM=BPPARAM)
    n.above <- Reduce("+", out.values)
    n.above[o] <- n.above

    # Computing a p-value for each observed probability 
    limited <- n.above==0L
    pval <- (n.above+1)/(npts+1)

    # Collating into some sensible output.
    all.p <- all.lr <- all.exp <- rep(NA_real_, ncells)
    all.lim <- rep(NA, ncells)
    all.p[keep] <- pval
    all.lr[keep] <- obs.P + lfactorial(obs.totals)
    all.lim[keep] <- limited
    return(DataFrame(Total=umi.sum, Probability=all.lr, PValue=all.p, Limited=all.lim, row.names=colnames(m)))
}

.monte_carlo_pval <- function(total.val, total.len, P, ambient, iterations) { 
    .Call(cxx_calculate_pval, total.val, total.len, P, ambient, iterations) 
}

emptyDrops <- function(m, lower=100, scale=1, test.args=list(), barcode.args=list()) 
# Combined function that puts these all together, always keeping cells above the inflection
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    stats <- do.call(testEmptyDrops, c(list(m, lower=lower), test.args))
    kneept <- do.call(barcodeRanks, c(list(m, lower=lower), barcode.args))$knee
    always <- stats$Total >= kneept*scale
    tmp <- stats$PValue
    tmp[always] <- 0
    stats$FDR <- p.adjust(tmp, method="BH")
    return(stats)
}
