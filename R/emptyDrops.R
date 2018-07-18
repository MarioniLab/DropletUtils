#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame
#' @importFrom edgeR goodTuringProportions
testEmptyDrops <- function(m, lower=100, niters=10000, test.ambient=FALSE, ignore=NULL, alpha=Inf, BPPARAM=SerialParam()) 
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
    # Also removing additional cells that don't pass some total count threshold, if required.
    if (!test.ambient) {
        keep <- !ambient
    } else {
        keep <- umi.sum > 0L
    }
    if (!is.null(ignore)) { 
        keep <- keep & umi.sum > ignore
    }
    obs <- m[,keep,drop=FALSE]
    obs.totals <- umi.sum[keep]

    # Estimating the alpha from the discarded ambient droplets, if desired.
    if (is.null(alpha)) {
        alpha <- .estimate_alpha(m[,ambient,drop=FALSE], ambient.prop, umi.sum[ambient]) 
    }

    # Calculating the log-multinomial probability for each cell.
    obs.P <- .compute_multinom_prob_data(obs, ambient.prop, alpha=alpha)
    rest.P <- .compute_multinom_prob_rest(obs.totals, alpha=alpha)

    # Computing the p-value for each observed probability.
    n.above <- .permute_counter(totals=obs.totals, probs=obs.P, ambient=ambient.prop, iter=niters, BPPARAM=BPPARAM, alpha=alpha)
    limited <- n.above==0L
    pval <- (n.above+1)/(niters+1)

    # Collating into some sensible output.
    all.p <- all.lr <- all.exp <- rep(NA_real_, ncells)
    all.lim <- rep(NA, ncells)
    all.p[keep] <- pval
    all.lr[keep] <- obs.P + rest.P 
    all.lim[keep] <- limited
    return(DataFrame(Total=umi.sum, LogProb=all.lr, PValue=all.p, Limited=all.lim, row.names=colnames(m)))
}

#' @importFrom BiocParallel bpnworkers SerialParam bplapply
.permute_counter <- function(totals, probs, ambient, iter, alpha=Inf, BPPARAM=SerialParam()) 
# Calculating the p-values using a Monte Carlo approach.
{
    o <- order(totals, probs)
    re.P <- probs[o]
    re.totals <- rle(totals[o]) # Ensure identity upon comparison.

    nworkers <- bpnworkers(BPPARAM)
    per.core <- rep(ceiling(iter/nworkers), nworkers)
    per.core[1] <- iter - sum(per.core[-1]) # Making sure that we get the exact number of iterations.

    out.values <- bplapply(per.core, FUN=.monte_carlo_pval, total.val=re.totals$values,
        total.len=re.totals$lengths, P=re.P, ambient=ambient, alpha=alpha, BPPARAM=BPPARAM)
    n.above <- Reduce("+", out.values)
    n.above[o] <- n.above
    return(n.above)
}

.monte_carlo_pval <- function(total.val, total.len, P, ambient, iterations, alpha) { 
    if (is.infinite(alpha)) {
        out <- .Call(cxx_montecarlo_pval, total.val, total.len, P, ambient, iterations) 
    } else {
        out <- .Call(cxx_montecarlo_pval_alpha, total.val, total.len, P, ambient, iterations, alpha) 
    }
    return(out)
}

#' @importFrom methods is
#' @importClassesFrom Matrix dgTMatrix
#' @importFrom stats aggregate
.compute_multinom_prob_data <- function(mat, prop, alpha=Inf)
# Efficiently calculates the data-dependent component of the log-multinomial probability.
# Also does so for the Dirichlet-multinomial log-probability for a given 'alpha'.
{
    if (is(mat, "dgTMatrix")) {
        i <- mat@i + 1L
        j <- mat@j + 1L
        x <- mat@x

        if (is.infinite(alpha)) {
            p.n0 <- x * log(prop[i]) - lfactorial(x)
        } else {
            p.n0 <- lgamma(alpha * prop[i] + x) - lfactorial(x) - lgamma(alpha * prop[i])
        }

        by.col <- aggregate(p.n0, list(Col=j), sum)
        obs.P <- numeric(ncol(mat))
        obs.P[by.col$Col] <- by.col$x
    } else {
        obs.P <- .Call(cxx_compute_multinom, mat, prop, alpha)
    }
    return(obs.P)
}

.compute_multinom_prob_rest <- function(totals, alpha=Inf) 
# Efficiently calculates the total-dependent component of the multinomial log-probability.
{
    if (is.infinite(alpha)) { 
        return(lfactorial(totals))
    } else {
        return(lfactorial(totals) + lgamma(alpha) - lgamma(totals + alpha))
    }
}

#' @importFrom methods is
#' @importClassesFrom Matrix dgTMatrix
#' @importMethodsFrom Matrix which
#' @importFrom stats optimize 
.estimate_alpha <- function(mat, prop, totals, interval=c(0.01, 10000))
# Efficiently finds the MLE for the overdispersion parameter of a Dirichlet-multinomial distribution.
{
    if (is(mat, "dgTMatrix")) {
        i <- mat@i + 1L
        j <- mat@j + 1L
        x <- mat@x
    } else {
        keep <- which(mat > 0, arr.ind=TRUE)
        i <- keep[,1]
        j <- keep[,2]
        x <- mat[keep]
    }

    per.prop <- prop[i] 
    LOGLIK <- function(alpha) {
        output <- numeric(length(alpha))
        for (adx in seq_along(alpha)) {
            cur.alpha <- alpha[adx]
            output[adx] <- lgamma(cur.alpha) * length(totals) - 
                sum(lgamma(totals + cur.alpha)) + 
                sum(lgamma(x + per.prop * cur.alpha)) - 
                sum(lgamma(per.prop * cur.alpha))
        }
        return(output)
    }

    optimize(LOGLIK, interval=interval, maximum=TRUE)$maximum
}

#' @export
#' @importFrom stats p.adjust
emptyDrops <- function(m, lower=100, retain=NULL, barcode.args=list(), ...) 
# Combined function that puts these all together, always keeping cells above the inflection
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    stats <- testEmptyDrops(m, lower=lower, ...)
    tmp <- stats$PValue
    
    if (is.null(retain)) {
        retain <- do.call(barcodeRanks, c(list(m, lower=lower), barcode.args))$knee
    }
    always <- stats$Total >= retain
    tmp[always] <- 0

    stats$FDR <- p.adjust(tmp, method="BH")
    return(stats)
}
