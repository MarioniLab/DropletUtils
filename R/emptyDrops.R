#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom edgeR goodTuringProportions
testEmptyDrops <- function(m, lower=100, niters=10000, test.ambient=FALSE, ignore=NULL, alpha=NULL, BPPARAM=SerialParam()) 
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

    output <- DataFrame(Total=umi.sum, LogProb=all.lr, PValue=all.p, Limited=all.lim, row.names=colnames(m))
    metadata(output) <- list(lower=lower, niters=niters, ambient=ambient.prop, alpha=alpha)
    output
}

#' @importFrom BiocParallel bpnworkers SerialParam bpmapply
#' @importFrom stats runif
.permute_counter <- function(totals, probs, ambient, iter, alpha=Inf, BPPARAM=SerialParam()) 
# Calculating the p-values using a Monte Carlo approach.
{
    o <- order(totals, probs)
    re.P <- probs[o]
    re.totals <- rle(totals[o]) # Ensure identity upon comparison.

    nworkers <- bpnworkers(BPPARAM)
    if (iter < nworkers) {
        per.core <- integer(nworkers)
        per.core[seq_len(iter)] <- 1L
    } else {
        per.core <- rep(ceiling(iter/nworkers), nworkers)
        per.core[1] <- iter - sum(per.core[-1]) # Making sure that we get the exact number of iterations.
    }

    # Creating seeds for the C++ PRNG to avoid disrupting the R seed in multi-core execution.
    seeds.per.core <- streams.per.core <- vector("list", nworkers)
    last <- 0L
    for (i in seq_len(nworkers)) {
        n <- per.core[i]
        seeds.per.core[[i]] <- runif(n, 0, 2^32)
        streams.per.core[[i]] <- last + seq_len(n)
        last <- last + n
    }

    out.values <- bpmapply(iterations=per.core, seeds=seeds.per.core, streams=streams.per.core,
        FUN=.monte_carlo_pval, 
        MoreArgs=list(
            total.val=re.totals$values, 
            total.len=re.totals$lengths, 
            P=re.P, 
            ambient=ambient, 
            alpha=alpha
        ), BPPARAM=BPPARAM, SIMPLIFY=FALSE)

    n.above <- Reduce("+", out.values)
    n.above[o] <- n.above
    return(n.above)
}

.monte_carlo_pval <- function(total.val, total.len, P, ambient, iterations, alpha, seeds, streams) 
# Wrapper function to preserve NAMESPACE in bpmapply.
{ 
    .Call(cxx_montecarlo_pval, total.val, total.len, P, ambient, iterations, alpha, seeds, streams) 
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
            prop.alpha <- per.prop * cur.alpha
            output[adx] <- lgamma(cur.alpha) * length(totals) - 
                sum(lgamma(totals + cur.alpha)) + 
                sum(lgamma(x + prop.alpha)) - 
                sum(lgamma(prop.alpha))
        }
        return(output)
    }

    optimize(LOGLIK, interval=interval, maximum=TRUE)$maximum
}

#' @export
#' @importFrom stats p.adjust
#' @importFrom S4Vectors metadata<- metadata
emptyDrops <- function(m, lower=100, retain=NULL, barcode.args=list(), ...) 
# Combined function that puts these all together, always keeping cells above the knee 
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    stats <- testEmptyDrops(m, lower=lower, ...)
    tmp <- stats$PValue
    
    if (is.null(retain)) {
        br.out <- do.call(barcodeRanks, c(list(m, lower=lower), barcode.args))
        retain <- metadata(br.out)$knee
    }
    always <- stats$Total >= retain
    tmp[always] <- 0

    metadata(stats)$retain <- retain
    stats$FDR <- p.adjust(tmp, method="BH")
    return(stats)
}
