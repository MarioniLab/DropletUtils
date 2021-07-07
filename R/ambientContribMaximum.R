#' Ambient contribution by maximum scaling
#'
#' Compute the maximum contribution of the ambient solution to an expression profile for a group of droplets,
#' by scaling the ambient profile and testing for significant deviations in the count profile.
#'
#' @param y A numeric matrix-like object containing counts, where each row represents a gene and each column represents a cluster of cells (see Caveats).
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' 
#' \code{y} can also be a numeric vector of counts; this is coerced into a one-column matrix.
#' @param ambient A numeric vector of length equal to \code{nrow(y)},
#' containing the proportions of transcripts for each gene in the ambient solution.
#' Alternatively, a matrix where each row corresponds to a row of \code{y}
#' and each column contains a specific ambient profile for the corresponding column of \code{y}.
#' @param threshold Numeric scalar specifying the p-value threshold to use, see Details.
#' @param dispersion Numeric scalar specifying the dispersion to use in the negative binomial model.
#' Defaults to zero, i.e., a Poisson model.
#' @param num.points Integer scalar specifying the number of points to use for the grid search.
#' @param num.iter Integer scalar specifying the number of iterations to use for the grid search.
#' @param mode String indicating the output to return, see Value.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{controlAmbience}, arguments to pass to \code{ambientContribMaximum}.
#' 
#' @return 
#' If \code{mode="scale"},
#' a numeric vector is returned quantifying the maximum \dQuote{contribution} of the ambient solution to each column of \code{y}.
#' Scaling \code{ambient} by each entry yields the maximum ambient profile for the corresponding column of \code{y}.
#'
#' If \code{mode="profile"}, a numeric matrix is returned containing the maximum ambient profile for each column of \code{y}.
#' This is computed by scaling as described above; if \code{ambient} is a matrix, each column is scaled by the corresponding entry of the scaling vector.
#'
#' If \code{mode="proportion"}, a numeric matrix is returned containing the maximum proportion of counts in \code{y} that are attributable to ambient contamination.
#' This is computed by simply dividing the output of \code{mode="profile"} by \code{y} and capping all values at 1.
#'
#' @details
#' On occasion, it is useful to estimate the maximum possible contribution of the ambient solution to a count profile.
#' This represents the most pessimistic explanation of a particular expression pattern
#' and can be used to identify and discard suspect genes or clusters prior to downstream analyses.
#'
#' This function implements the following algorithm:
#' \enumerate{ 
#' \item We compute the mean ambient contribution for each gene by scaling \code{ambient} by some factor.
#' \code{ambient} itself is usually derived by summing counts across barcodes with low total counts,
#' see the output of \code{\link{emptyDrops}} for an example.
#' \item We compute a p-value for each gene based on the probability of observing a count equal to or below that in \code{y}, using the lower tail of a negative binomial (or Poisson) distribution with mean set to the ambient contribution.
#' The per-gene null hypothesis is that the expected count in \code{y} is equal to the sum of the scaled ambient proportion and some (non-negative) contribution from actual intra-cellular transcripts.
#' \item We combine p-values across all genes using Simes' method.
#' This represents the evidence against the joint null hypothesis (that all of the per-gene nulls are true).
#' \item We find the largest scaling factor that fails to reject this joint null at the specified \code{threshold}.
#' If \code{sum(ambient)} is equal to unity, this scaling factor can be interpreted as the maximum number of transcript molecules contributed to \code{y} by the ambient solution.
#' }
#' 
#' The process of going from a scaling factor to a combined p-value has no clean analytical solution,
#' so we use an iterative grid search to identify to largest possible scaling factor at a decent resolution.
#' \code{num.points} and \code{num.iter} control the resolution of the grid search,
#' and generally do not need to be changed.
#'
#' \code{maximumAmbience} is soft-deprecated; use \code{ambientContribMaximum} instead.
#'
#' @section Caveats:
#' The above algorithm is rather \emph{ad hoc} and offers little in the way of theoretical guarantees.
#' The p-value is used as a score rather than providing any meaningful error control.
#' Empirically, increasing \code{threshold} will return a higher scaling factor by making the estimation more robust to drop-outs in \code{y}, at the cost of increasing the risk of over-estimation of the ambient contribution.
#'
#' Our abuse of the p-value machinery means that the reported scaling often exceeds the actual contribution, especially at low counts where the reduced power fails to penalize overly large scaling factors.
#' Hence, the function works best when \code{y} contains aggregated counts for one or more groups of droplets with the same expected expression profile, e.g., clusters of related cells.
#' Higher counts provide more power to detect deviations, hopefully leading to a more accurate estimate of the scaling factor.
#' (On a practical note, this function is rather slow so it is more feasible to calculate on cluster-level profiles rather than per cell.)
#' 
#' Note that this function returns the \emph{maximum possible} contribution of the ambient solution to \code{y}, not the actual contribution.
#' In the most extreme case, if the ambient profile is similar to the expectation of \code{y} (e.g., due to sequencing a relatively homogeneous cell population), the maximum possible contribution of the ambient solution would be 100\% of \code{y}, and subtraction would yield an empty count vector!
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{ambientProfileEmpty}} and \code{\link{ambientProfileBimodal}}, to estimate the ambient profile.
#'
#' \code{\link{ambientContribSparse}} and \code{\link{ambientContribNegative}}, for other methods to estimate the ambient contribution.
#'
#' \code{\link{emptyDrops}}, which uses the ambient profile to call cells.
#'
#' @examples
#' # Making up some data for, e.g., a single cluster.
#' ambient <- c(runif(900, 0, 0.1), runif(100))
#' y <- rpois(1000, ambient * 100)
#' y[1:100] <- y[1:100] + rpois(100, 20) # actual biology.
#'
#' # Estimating the maximum possible scaling factor:
#' scaling <- ambientContribMaximum(y, ambient)
#' scaling
#'
#' # Estimating the maximum contribution to 'y' by 'ambient'.
#' contribution <- ambientContribMaximum(y, ambient, mode="profile")
#' DataFrame(ambient=drop(contribution), total=y)
#' 
#' @seealso 
#' \code{\link{ambientProfileEmpty}} or \code{\link{ambientProfileBimodal}}, to obtain an estimate to use in \code{ambient}.
#'
#' \code{\link{ambientContribNegative}} or \code{\link{ambientContribSparse}}, for other methods of estimating the contribution.
#'
#' @name ambientContribMaximum
NULL

#' @importFrom stats p.adjust ppois pnbinom
#' @importFrom BiocParallel SerialParam
.ambient_contrib_maximum <- function(y, ambient, threshold=0.1, dispersion=0, 
    num.points=100, num.iter=5, mode=c("scale", "profile", "proportion"), BPPARAM=SerialParam()) 
{
    mode <- match.arg(mode)

    if (is.null(dim(y))) {
        y <- cbind(y)
        colnames(y) <- NULL
    }

    scaling <- colBlockApply(y, BPPARAM=BPPARAM, .maximum_ambience_loop, ambient=ambient, 
        threshold=threshold, dispersion=dispersion, num.points=num.points, num.iter=num.iter)
    scaling <- unlist(scaling)

    .report_ambient_profile(scaling, ambient=ambient, y=y, mode=match.arg(mode))
}

#' @importFrom DelayedArray currentViewport makeNindexFromArrayViewport
.maximum_ambience_loop <- function(y, ambient, ...) {
    vp <- currentViewport()
    yidx <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)[[2]]
    output <- numeric(ncol(y))

    for (i in seq_along(output)) {
        if (is.null(dim(ambient))) {
            A <- ambient
        } else {
            if (is.null(yidx)) {
                idx <- i
            } else {
                idx <- yidx[i]
            }
            A <- ambient[,idx]
        }
        output[i] <- .maximum_ambience(y=y[,i], ambient=A, ...)
    }

    output
}

.maximum_ambience <- function(y, ambient, threshold, dispersion, num.points, num.iter) {
    if (dispersion==0) {
        FUN <- function(y, mu) {
            ppois(y, lambda=mu)
        }
    } else {
        FUN <- function(y, mu) {
            pnbinom(y, mu=mu, size=1/dispersion)
        }
    }

    if (!identical(names(y), names(ambient))) {
        warning("'y' and 'ambient' do not have the same feature names")
    }

    # Removing all-zero genes in the ambient profile.
    strip <- ambient==0
    y <- y[!strip]
    ambient <- ambient[!strip]

    # Executing a grid search.
    lower <- 0
    upper <- sum(y)/sum(ambient)
    iter <- 0

    while (iter <= num.iter) {
        scaling <- seq(from=lower, to=upper, length.out=num.points)
        triggered <- FALSE 

        for (i in 2:num.points) {
            p <- FUN(y, ambient * scaling[i])
            p <- p.adjust(p, method="BH")

            if (any(p <= threshold)) {
                lower <- scaling[i-1]
                upper <- scaling[i]
                triggered <- TRUE
                break
            }
        }

        if (!triggered) {
            lower <- scaling[num.points-1]
            upper <- scaling[num.points]
        }

        iter <- iter+1L
    }

    (lower+upper)/2
}

.report_ambient_profile <- function(scaling, ambient, y, mode) {
    if (mode=="scale") {
        names(scaling) <- colnames(y)
        scaling
    } else {
        if (is.null(dim(ambient))) {
            scaled.ambient <- outer(ambient, scaling)
        } else {
            scaled.ambient <- t(t(ambient) * scaling)
        }
        dimnames(scaled.ambient) <- dimnames(y)

        if (mode=="profile") {
            scaled.ambient
        } else {
            p <- scaled.ambient/y
            p[p > 1] <- 1
            p[y == 0] <- NaN # for fairness's sake.
            p
        }
    }
}

#' @export
#' @rdname ambientContribMaximum
maximumAmbience <- function(...) {
    ambientContribMaximum(...)
}

#' @export
#' @rdname ambientContribMaximum
setGeneric("ambientContribMaximum", function(y, ...) standardGeneric("ambientContribMaximum"))

#' @export 
#' @rdname ambientContribMaximum
setMethod("ambientContribMaximum", "ANY", .ambient_contrib_maximum)

#' @export
#' @rdname ambientContribMaximum
#' @importFrom SummarizedExperiment assay
setMethod("ambientContribMaximum", "SummarizedExperiment", function(y, ..., assay.type="counts") {
    .ambient_contrib_maximum(assay(y, assay.type), ...)
})
