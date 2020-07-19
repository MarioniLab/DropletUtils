#' Maximum ambient contribution
#'
#' Compute the maximum contribution of the ambient solution to a particular expression profile.
#'
#' @param y A numeric count matrix where each row represents a gene and each column represents an expression profile
#' (usually aggregated across multiple droplets in a sample).
#' This can also be a vector, in which case it is converted into a one-column matrix.
#' @param ambient A numeric vector of length equal to \code{nrow(y)},
#' containing the proportions of transcripts for each gene in the ambient solution.
#' Alternatively, a matrix where each row corresponds to a row of \code{y}
#' and each column contains a specific ambient profile for the corresponding column of \code{y}.
#' @param threshold Numeric scalar specifying the p-value threshold to use, see Details.
#' @param dispersion Numeric scalar specifying the dispersion to use in the negative binomial model.
#' Defaults to zero, i.e., a Poisson model.
#' @param num.points Integer scalar specifying the number of points to use for the grid search.
#' @param num.iter Integer scalar specifying the number of iterations to use for the grid search.
#' @param mode String indicating the output to return - the scaling factor, the ambient profile or the proportion of each gene's counts in \code{y} that is attributable to ambient contamination.
#' 
#' @return 
#' If \code{mode="scale"},
#' a numeric vector is returned quantifying the maximum \dQuote{contribution} of the ambient solution to each column of \code{y}.
#' Scaling columns of \code{ambient} by this vector yields the maximum ambient profile for each column of \code{y},
#' which can also be obtained by setting \code{mode="profile"}.
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
#' \item It computes the mean ambient contribution for each gene by scaling \code{ambient} by some factor.
#' \code{ambient} itself is usually derived by summing counts across barcodes with low total counts,
#' see the output of \code{\link{emptyDrops}} for an example.
#' \item It computes a p-value for each gene based on the probability of observing a count equal to or below that in \code{y},
#' using the lower tail of a negative binomial (or Poisson) distribution with mean set to the ambient contribution.
#' \item It combines p-values across genes using Simes' method.
#' The joint null hypothesis is that the expectation of \code{y} is equal to the sum of the scaled ambient proportions 
#' and some (non-negative) contribution from actual intra-cellular transcripts.
#' \item It finds the largest scaling factor that fails to reject this joint null at the specified \code{threshold}.
#' If \code{sum(ambient)} is equal to unity, this scaling factor can be interpreted as the maximum number of transcript molecules contributed to \code{y} by the ambient solution.
#' }
#' 
#' The process of going from a scaling factor to a combined p-value has no clean analytical solution,
#' so we use an iterative grid search to identify to largest possible scaling factor at a decent resolution.
#' \code{num.points} and \code{num.iter} control the resolution of the grid search,
#' and generally do not need to be changed.
#'
#' @section Caveats:
#' The algorithm implemented in this function is, admittedly, rather \emph{ad hoc} and offers little in the way of theoretical guarantees.
#' The reported scaling often exceeds the actual contribution, especially at low counts where the reduced power fails to penalize overly large scaling factors.
#' The p-value is largely used as a score rather than providing any meaningful error control.
#' Empirically, decreasing \code{threshold} will return a higher scaling factor by making the estimation more robust to drop-outs in \code{y}, at the cost of increasing the risk of over-estimation of the ambient contribution.
#'
#' It is also important to note that this function returns the \emph{maximum possible} contribution of the ambient solution to \code{y}, not the actual contribution.
#' It is probably unwise to attempt to obtain a \dQuote{cleaned} expression profile by subtracting the scaled ambient proportions from \code{y}.
#' In the most extreme case, if the ambient profile is similar to the expectation of \code{y} (e.g., due to sequencing a relatively homogeneous cell population), the maximum possible contribution of the ambient solution would be 100\% of \code{y}, and subtraction would yield an empty count vector!
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{emptyDrops}}, which uses the ambient profile to call cells.
#'
#' @examples
#' # Making up some data.
#' ambient <- c(runif(900, 0, 0.1), runif(100))
#' y <- rpois(1000, ambient * 50)
#' y <- y + rpois(1000, 5) # actual biology.
#'
#' # Estimating the maximum possible scaling factor:
#' scaling <- maximumAmbience(y, ambient)
#' scaling
#'
#' # Estimating the maximum contribution to 'y' by 'ambient'.
#' contribution <- maximumAmbience(y, ambient, mode="profile")
#' DataFrame(ambient=drop(contribution), total=y)
#' 
#' @seealso 
#' \code{\link{estimateAmbience}}, to obtain an estimate to use in \code{ambient}.
#'
#' \code{\link{controlAmbience}}, for a more accurate estimate when control features are available.
#' @export
#' @importFrom stats p.adjust ppois pnbinom
maximumAmbience <- function(y, ambient, threshold=0.1, dispersion=0, 
    num.points=100, num.iter=5, mode=c("scale", "profile", "proportion")) 
{
    mode <- match.arg(mode)
    args <- list(threshold=threshold, dispersion=dispersion, num.points=num.points, num.iter=num.iter, mode=mode)

    if (is.null(dim(y))) {
        y <- cbind(y)
        colnames(y) <- NULL
    }

    collated <- vector("list", ncol(y))
    names(collated) <- colnames(y)

    for (i in seq_along(collated)) {
        if (is.null(dim(ambient))) {
            A <- ambient
        } else {
            A <- ambient[,i]
        }
        collated[[i]] <- do.call(.maximum_ambience, c(list(y=y[,i], ambient=A), args))
    }

    if (mode=="scale") {
        unlist(collated)
    } else {
        do.call(cbind, collated)
    }
}

.maximum_ambience <- function(y, ambient, threshold, dispersion, num.points, num.iter, mode) {
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
    original.y <- y 
    original.ambient <- ambient
    
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

    scale <- (lower+upper)/2

    if (mode=="scale") {
        scale
    } else {
        profile <- scale * original.ambient
        if (mode=="profile") {
            profile
        } else {
            prop <- .clean_amb_props(profile, original.y)
            names(prop) <- names(original.ambient)
            prop
        }
    }
}

.clean_amb_props <- function(s, y) {
    p <- s/y
    p[p > 1] <- 1
    p[y == 0] <- NaN # for fairness's sake.
    p
}
