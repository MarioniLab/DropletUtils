#' Identify empty droplets
#'
#' Distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment.
#' 
#' @param m A numeric matrix object - usually a \linkS4class{dgTMatrix} or \linkS4class{dgCMatrix} - 
#' containing droplet data \emph{prior to any filtering or cell calling}.
#' Columns represent barcoded droplets, rows represent genes.
#' @param lower A numeric scalar specifying the lower bound on the total UMI count, 
#' at or below which all barcodes are assumed to correspond to empty droplets.
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
#' @param test.ambient A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to \code{lower}.
#' @param ignore A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored (see Details for how this differs from \code{lower}).
#' @param alpha A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme.
#' @param BPPARAM A BiocParallelParam object indicating whether parallelization should be used to compute p-values.
#' @param retain A numeric scalar specifying the threshold for the total UMI count above which all barcodes are assumed to contain cells.
#' @param barcode.args Further arguments to pass to \code{\link{barcodeRanks}}.
#' @param round Logical scalar indicating whether to check for non-integer values in \code{m} and, if present, round them for ambient profile estimation (see \code{?\link{estimateAmbience}}) and the multinomial simulations.
#' @param ... Further arguments to pass to \code{testEmptyDrops}.
#' 
#' @section Details about \code{testEmptyDrops}:
#' The \code{testEmptyDrops} function first obtains an estimate of the composition of the ambient pool of RNA based on the barcodes with total UMI counts less than or equal to \code{lower} (see \code{?\link{estimateAmbience}} for details).
#' This assumes that a cell-containing droplet would generally have higher total counts than empty droplets containing RNA from the ambient pool.
#' Counts for the low-count barcodes are pooled together, and an estimate of the proportion vector for the ambient pool is calculated using \code{\link{goodTuringProportions}}.
#' The count vector for each barcode above \code{lower} is then tested for a significant deviation from these proportions.
#' 
#' Then, \code{testEmptyDrops} will test each barcode for significant deviations from the ambient profile.
#' The null hypothesis is that transcript molecules are included into droplets by multinomial sampling from the ambient profile.
#' For each barcode, the probability of obtaining its count vector based on the null model is computed.
#' Then, \code{niters} count vectors are simulated from the null model.
#' The proportion of simulated vectors with probabilities lower than the observed multinomial probability for that barcode is used to calculate the p-value.
#' 
#' We use this Monte Carlo approach as an exact multinomial p-value is difficult to calculate.
#' However, the p-value is lower-bounded by the value of \code{niters} (Phipson and Smyth, 2010), which can result in loss of power if \code{niters} is too small.
#' Users can check whether this loss of power has any practical consequence by checking the \code{Limited} field in the output.
#' If any barcodes have \code{Limited=TRUE} but does \emph{not} reject the null hypothesis, it suggests that \code{niters} should be increased.
#' 
#' The stability of the Monte Carlo $p$-values depends on \code{niters}, which is only set to a default of 10000 for speed.
#' Larger values improve stability with the only cost being that of time, so users should set \code{niters} to the largest value they are willing to wait for.
#' 
#' The \code{ignore} argument can also be set to ignore barcodes with total counts less than or equal to \code{ignore}.
#' This differs from the \code{lower} argument in that the ignored barcodes are not necessarily used to compute the ambient profile.
#' Users can interpret \code{ignore} as the minimum total count required for a barcode to be considered as a potential cell.
#' In contrast, \code{lower} is the maximum total count below which all barcodes are assumed to be empty droplets.
#' 
#' @section Details about \code{emptyDrops}:
#' The \code{emptyDrops} function identifies droplets that are likely to contain cells by calling \code{testEmptyDrops}.
#' The Benjamini-Hochberg correction is applied to the Monte Carlo p-values to correct for multiple testing.
#' Cells can then be defined by taking all barcodes with significantly non-ambient profiles, e.g., at a false discovery rate of 0.1\%.
#' 
#' Barcodes that contain more than \code{retain} total counts are always retained.
#' This ensures that large cells with profiles that are very similar to the ambient pool are not inadvertently discarded.
#' If \code{retain} is not specified, it is set to the total count at the knee point detected by \code{\link{barcodeRanks}}.
#' Manual specification of \code{retain} may be useful if the knee point was not correctly identified in complex log-rank curves.
#' Users can also set \code{retain=Inf} to disable automatic retention of barcodes with large totals.
#' 
#' All barcodes with total counts above \code{retain} are assigned p-values of zero \emph{during correction}, reflecting our assumption that they are true positives.
#' This ensures that their Monte Carlo p-values do not affect the correction of other genes, and also means that they will have FDR values of zero.
#' However, their original Monte Carlo p-values are still reported in the output, as these may be useful for diagnostic purposes.
#' 
#' In general, users should call \code{emptyDrops} rather than \code{testEmptyDrops}.
#' The latter is a \dQuote{no frills} version that is largely intended for use within other functions.
#' 
#' @section Handling overdispersion:
#' If \code{alpha} is set to a positive number, sampling is assumed to follow a Dirichlet-multinomial (DM) distribution.
#' The parameter vector of the DM distribution is defined as the estimated ambient profile scaled by \code{alpha}.
#' Smaller values of \code{alpha} model overdispersion in the counts, due to dependencies in sampling between molecules.
#' If \code{alpha=NULL}, a maximum likelihood estimate is obtained from the count profiles for all barcodes with totals less than or equal to \code{lower}.
#' If \code{alpha=Inf}, the sampling of molecules is modelled with a multinomial distribution.
#' 
#' Users can check whether the model is suitable by extracting the p-values for all barcodes with \code{test.ambient=TRUE}.
#' Under the null hypothesis, the p-values for presumed ambient barcodes (i.e., with total counts below \code{lower}) should be uniformly distributed.
#' Skews in the p-value distribution are indicative of an inaccuracy in the model and/or its estimates (of \code{alpha} or the ambient profile).
#' 
#' @section \code{NA} values in the results:
#' We assume that barcodes with total UMI counts below \code{lower} correspond to empty droplets.
#' These are used to estimate the ambient expression profile against which the remaining barcodes are tested.
#' Under this definition, these low-count barcodes cannot be cell-containing droplets and are excluded from the hypothesis testing.
#' 
#' However, it is still desirable for the number of rows of the output DataFrame to be the same as \code{ncol(m)}.
#' This allows easy subsetting of \code{m} based on a logical vector constructed from the output (e.g., to retain all FDR values below a threshold).
#' To satisfy this requirement, the rows for the excluded barcodes are filled in with \code{NA} values for all fields in the output.
#' We suggest using \code{\link{which}} to pick barcodes below a FDR threshold, see the Examples.
#' 
#' If \code{test.ambient=FALSE}, non-\code{NA} statistics will be reported for all barcodes.
#' This is occasionally useful for diagnostics to ensure that the p-values are well-calibrated for barcodes below \code{lower}.
#' Specifically, if the null hypothesis were true, p-values for low-count barcodes should have a uniform distribution.
#' Any strong peaks in the p-values near zero indicate that \code{emptyDrops} is not controlling the FDR correctly.
#'
#' @section Non-empty droplets versus cells:
#' Technically speaking, \code{emptyDrops} is designed to identify barcodes that correspond to non-empty droplets.
#' This is close to but not quite the same as identifying cells,
#' as droplets containing cell fragments, stripped nuclei and damaged cells will still be significantly non-empty.
#' As such, it may often be necessary to perform additional quality control on the significant barcodes;
#' we suggest doing so using methods from the \pkg{scater} package.
#'
#' On occasion, \code{emptyDrops} may identify many more non-empty droplets than the expected number of cells.
#' This is probably due to the generation of multiple cell fragments when a single cell is extensively damaged.
#' In such cases, it is informative to construct a MA plot comparing the average expression between retained low-count barcodes and discarded barcodes to see which genes are driving the differences (and thus contributing to the larger number of non-empty calls).
#' Mitochondrial and ribosomal genes are typical offenders; the former can be either up or down in the ambient solution, depending on whether the damage was severe enough to dissociate mitochondria from the cell fragments, while the latter is usually down in low-count barcodes due to loss of cytoplasmic RNA in cell fragments.
#'
#' To mitigate this effect, we can filtering out the problematic genes from the matrix provided to \code{emptyDrops}.
#' This eliminates their effect on the significance calculations and reduces the number of uninteresting non-empty calls,
#' see \url{https://github.com/MarioniLab/DropletUtils/issues/36} for an example.
#' Of course, the full set of genes can still be retained for downstream analysis.
#' 
#' @return
#' \code{testEmptyDrops} will return a DataFrame with the following components:
#' \describe{
#' \item{\code{Total}:}{Integer, the total UMI count for each barcode.}
#' \item{\code{LogProb}:}{Numeric, the log-probability of observing the barcode's count vector under the null model.}
#' \item{\code{PValue}:}{Numeric, the Monte Carlo p-value against the null model.}
#' \item{\code{Limited}:}{Logical, indicating whether a lower p-value could be obtained by increasing \code{niters}.}
#' }
#' 
#' \code{emptyDrops} will return a DataFrame like \code{testEmptyDrops}, with an additional \code{FDR} field.
#' 
#' The metadata of the output DataFrame will contains the ambient profile in \code{ambient}, the estimated/specified value of \code{alpha}, the specified value of \code{lower} and the number of iterations in \code{niters}.
#' For \code{emptyDrops}, the metadata will also contain the retention threshold in \code{retain}.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts()
#' 
#' # Identify likely cell-containing droplets.
#' out <- emptyDrops(my.counts)
#' out
#' 
#' is.cell <- out$FDR <= 0.01
#' sum(is.cell, na.rm=TRUE)
#'
#' # Subsetting the matrix to the cell-containing droplets.
#' # (using 'which()' to handle NAs smoothly).
#' cell.counts <- my.counts[,which(is.cell),drop=FALSE]
#' dim(cell.counts)
#' 
#' # Check if p-values are lower-bounded by 'niters'
#' # (increase 'niters' if any Limited==TRUE and Sig==FALSE)
#' table(Sig=is.cell, Limited=out$Limited)
#' 
#' @references
#' Lun A, Riesenfeld S, Andrews T, Dao TP, Gomes T, participants in the 1st Human Cell Atlas Jamboree, Marioni JC (2019).
#' Distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data.
#' \emph{Genome Biol.} 20, 63.
#' 
#' Phipson B, Smyth GK (2010).
#' Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
#' \emph{Stat. Appl. Genet. Mol. Biol.} 9:Article 39.
#' 
#' @seealso
#' \code{\link{barcodeRanks}}, for choosing the knee point.
#'
#' \code{\link{defaultDrops}}, for an implementation of the cell-calling method used by CellRanger version 2.
#'
#' \code{\link{estimateAmbience}}, for more details on estimation of the ambient profile.
#'
#' \code{\link{maximumAmbience}}, for estimating the maximum possible contribution of the ambient solution to a count profile.
#' 
#' @name emptyDrops
#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom Matrix rowSums colSums
testEmptyDrops <- function(m, lower=100, niters=10000, test.ambient=FALSE, ignore=NULL, alpha=NULL, 
    round=TRUE, BPPARAM=SerialParam()) 
{
    m <- .rounded_to_integer(m, round)
    astats <- .compute_ambient_stats(m, lower=lower)
    m <- astats$m
    umi.sum <- astats$umi.sum
    ambient <- astats$ambient
    ambient.prop <- astats$ambient.prop

    # Estimating the alpha from the discarded ambient droplets, if desired.
    if (is.null(alpha)) {
        ambient.m <- astats$ambient.m
        ambient.totals <- umi.sum[ambient]
        alpha <- .estimate_alpha(ambient.m, ambient.prop, ambient.totals) 
    }

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
    obs.m <- m[,keep,drop=FALSE]
    obs.totals <- umi.sum[keep]

    # Calculating the log-multinomial probability for each cell.
    obs.P <- .compute_multinom_prob_data(obs.m, ambient.prop, alpha=alpha)
    rest.P <- .compute_multinom_prob_rest(obs.totals, alpha=alpha)

    # Computing the p-value for each observed probability.
    n.above <- .permute_counter(totals=obs.totals, probs=obs.P, 
        ambient=ambient.prop, iter=niters, BPPARAM=BPPARAM, alpha=alpha)
    limited <- n.above==0L
    pval <- (n.above+1)/(niters+1)

    # Collating into some sensible output.
    ncells <- ncol(m)
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

    pcg.state <- .setup_pcg_state(per.core)

    out.values <- bpmapply(iterations=per.core, seeds=pcg.state$seeds, streams=pcg.state$streams,
        FUN=montecarlo_pval, 
        MoreArgs=list(
            totalval=re.totals$values, 
            totallen=re.totals$lengths, 
            prob=re.P, 
            ambient=ambient, 
            alpha=alpha
        ), BPPARAM=BPPARAM, SIMPLIFY=FALSE)

    n.above <- Reduce("+", out.values)
    n.above[o] <- n.above
    return(n.above)
}

#' @importFrom dqrng generateSeedVectors
.setup_pcg_state <- function(per.core) 
# Creating seeds for the C++ PRNG to avoid disrupting the R seed in multi-core execution.
{
    seeds.per.core <- streams.per.core <- vector("list", length(per.core))
    last <- 0L
    for (i in seq_along(per.core)) {
        N <- per.core[i]
        seeds.per.core[[i]] <- generateSeedVectors(N, nwords=2)
        streams.per.core[[i]] <- last + seq_len(N)
        last <- last + N
    }
    list(seeds=seeds.per.core, streams=streams.per.core)
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
        obs.P <- compute_multinom(mat, prop, alpha)
    }
    obs.P
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
#' @rdname emptyDrops
#' @importFrom stats p.adjust
#' @importFrom S4Vectors metadata<- metadata
emptyDrops <- function(m, lower=100, retain=NULL, barcode.args=list(), round=TRUE, ...) 
# Combined function that puts these all together, always keeping cells above the knee 
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    m <- .rounded_to_integer(m, round)
    stats <- testEmptyDrops(m, lower=lower, round=FALSE, ...)
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

#' @importFrom Matrix colSums rowSums
.rounded_to_integer <- function(m, round=TRUE) {
    if (round) {
        cs <- colSums(m)
        rs <- rowSums(m)
        if (!isTRUE(all.equal(cs, round(cs))) ||
            !isTRUE(all.equal(rs, round(rs)))) 
        {
            m <- round(m)
        }
    }
    m
}
