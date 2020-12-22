#' Ambient profile from bimodality
#'
#' Estimate the concentrations of tags in the ambient solution from a filtered HTO count matrix containing only counts for cells,
#' by assuming that each HTO has a bimodal abundance distribution with ambient and high-expressing components.
#'
#' @param x A numeric matrix-like object containing counts for each HTO (row) and cell (column).
#' @param min.prop Numeric scalar in (0, 1) specifying the expected minimum proportion of barcodes contributed by each sample.
#' @param ... Arguments to pass to \code{ambientProfileBimodal}.
#' 
#' @return A numeric vector of length equal to \code{nrow(x)}, containing the estimated ambient proportions for each HTO.
#'
#' @details
#' In some cases, we want to know the ambient profile but we only have the HTO count matrix for the cell-containing libraries.
#' This can be useful in functions such as \code{\link{hashedDrops}} or as a reference profile in \code{\link{medianSizeFactors}}.
#' However, as we only have the cell-containing libraries, we cannot use \code{\link{ambientProfileEmpty}}.
#' 
#' This function estimates the ambient profile by assuming that each HTO only labels a minority of the cells.
#' Under this assumption, each HTO's log-count distribution has a bimodal distribution where the lower mode represents ambient contamination.
#' We fit a two-component mixture model and identify all barcodes assigned to the lower component;
#' and the mean of those counts is used as an estimate of the ambient contribution.
#'
#' The initialization of the mixture model is controlled by \code{min.prop}, 
#' which starts the means of the lower and upper components at the \code{min.prop} and \code{1-min.prop} quantiles, respectively.
#' This means that each sample is expected to contribute somewhere between \code{[min.prop, 1-min.prop]} barcodes.
#' Larger values improve convergence but require stronger assumptions about the relative proportions of multiplexed samples.
#'
#' \code{inferAmbience} is soft-deprecated; use \code{ambientProfileBimodal} instead.
#'
#' @seealso
#' \code{\link{hashedDrops}}, where this function is used in the absence of an ambient profile.
#'
#' \code{\link{ambientProfileEmpty}}, which should be used when the raw matrix (prior to filtering for cells) is available.
#'
#' \code{\link{ambientContribSparse}} and related functions, to estimate the contribution of ambient contamination in each library.
#'
#' @author Aaron Lun
#' 
#' @examples
#' x <- rbind(
#'     rpois(1000, rep(c(100, 1), c(100, 900))),
#'     rpois(1000, rep(c(2, 100, 2), c(100, 100, 800))),
#'     rpois(1000, rep(c(3, 100, 3), c(200, 700, 100)))
#' )
#'
#' # Should be close to 1, 2, 3
#' ambientProfileBimodal(x)
#'
#' @export
ambientProfileBimodal <- function(x, min.prop=0.05) {
    ambient <- numeric(nrow(x))
    names(ambient) <- rownames(x)

    for (i in seq_along(ambient)) {
        current <- x[i,]
        chosen <- .get_lower_dist(log1p(current), min.prop)
        ambient[i] <- mean(current[chosen])
    }

    ambient
}

#' @export
#' @rdname ambientProfileBimodal
inferAmbience <- function(...) {
    ambientProfileBimodal(...)
}

#' @importFrom stats kmeans quantile
.get_lower_dist <- function(x, p)
# Effectively using Lloyd's algorithm as a special case of mixture models,
# to (i) avoid code dependencies and (ii) avoid problems with non-normal data.
{
    q <- quantile(x, c(p, 1-p))
    if (q[1]==q[2]) {
        # kmeans() fails in this case, so we just set everything less than or
        # equal to q[2] as the 'lower'. This should really only happen when 
        # almost all counts are zero, at which point this may as well be zero.
        x <= q[2]
    } else {
        out <- kmeans(x, centers=q)
        out$cluster == which.min(out$centers)
    }
}
