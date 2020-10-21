#' Infer the ambient profile
#'
#' Infer the ambient profile from a filtered HTO count matrix containing only counts for cells.
#'
#' @param x A numeric matrix-like object containing counts for each HTO (row) and cell (column).
#' @param min.prop Numeric scalar in (0, 1) specifying the expected minimum proportion of barcodes contributed by each sample.
#' 
#' @return A numeric vector of length equal to \code{nrow(x)}, containing the estimated ambient proportions for each HTO.
#'
#' @details
#' In some cases, we want to know the ambient profile but we only have the HTO count matrix for the cell-containing libraries.
#' This can be useful in functions such as \code{\link{hashedDrops}} or as a reference profile in \code{\link{medianSizeFactors}}.
#' However, as we only have the cell-containing libraries, we cannot use \code{\link{estimateAmbience}}.
#' 
#' This function allows us to obtain the ambient profile under the assumption that each HTO only labels a minority of the cells.
#' Specifically, it will fit a two-component mixture model to each HTO's count distribution.
#' All barcodes assigned to the lower component are considered to have background counts for that HTO,
#' and the mean of those counts is used as an estimate of the ambient contribution.
#'
#' The initialization of the mixture model is controlled by \code{min.prop}, 
#' which starts the means of the lower and upper components at the \code{min.prop} and \code{1-min.prop} quantiles, respectively.
#' This means that each sample is expected to contribute somewhere between \code{[min.prop, 1-min.prop]} barcodes.
#' Larger values improve convergence but require stronger assumptions about the relative proportions of multiplexed samples.
#'
#' @seealso
#' \code{\link{hashedDrops}}, where this function is used in the absence of an ambient profile.
#'
#' \code{\link{estimateAmbience}}, which should be used when the raw matrix (prior to filtering for cells) is available.
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
#' inferAmbience(x)
#'
#' @export
inferAmbience <- function(x, min.prop=0.05) {
    ambient <- numeric(nrow(x))
    names(ambient) <- rownames(x)

    for (i in seq_along(ambient)) {
        current <- x[i,]
        chosen <- .get_lower_dist(log1p(current), min.prop)
        ambient[i] <- mean(current[chosen])
    }

    ambient
}

#' @importFrom stats kmeans quantile
.get_lower_dist <- function(x, p)
# Effectively using Lloyd's algorithm as a special case of mixture models,
# to (i) avoid dependencies and (ii) avoid problems with non-normal data.
{
    out <- kmeans(x, centers=quantile(x, c(p, 1-p)))
    out$cluster == which.min(out$centers)
}
