#' Clean a HTO dataset
#'
#' Clean a HTO-containing count matrix where each column corresponds to a cell-containing barcode.
#'
#' @param x A numeric matrix-like object containing counts for each HTO (row) in each cell (column).
#' @param ambient A numeric vector of length equal to \code{nrow(x)}, containing the relative concentration of each HTO in the ambient solution.
#' Defaults to \code{\link{ambientProfileBimodal}(x)} if not explicitly provided.
#' @param controls A vector specifying the rows of \code{x} corresponding to control tags.
#' @param ... Further arguments to pass to \code{\link{isOutlier}}.
#' 
#' @author Aaron Lun
#'
#' @return A \linkS4class{DataFrame} with one row per column of \code{x}, containing the following fields:
#' \itemize{
#' \item \code{ambient.scale}, a numeric field containing the scaling factor for ambient contamination.
#' \item \code{zero.ambient}, a logical field indicating whether each cell has zero ambient contamination.
#' \item \code{high.ambient}, a logical field indicating whether each cell has unusually high levels of ambient contamination.
#' }
#'
#' If \code{controls} is supplied, the DataFrame also contains:
#' \itemize{
#' \item \code{sum.controls}, a numeric field containing the sum of counts for all control features.
#' \item \code{low.controls}, a logical field indicating whether each cell has unusually low levels of controls.
#' }
#'
#' Finally, the DataFrame contains \code{discard}, a logical field indicating whether a column in \code{x} should be discarded.
#'
#' @details
#' We remove cells for which there is no detectable ambient contamination, i.e., \code{ambient.scale} is zero.
#' We expect non-zero counts for most tags due to the deeply sequenced nature of HTO data.
#' If 50\% or more tags have zero counts, this is indicative of a failure in library preparation for that cell.
#'
#' We also remove cells for which the ambient contamination is unusually high, defined with \code{\link{isOutlier}} on the log-transformed \code{ambient.scale}.
#' These correspond to droplets that are contaminated by antibody conjugates, which cause all tags to have large counts.
#' The assumption here is that fewer than 50\% of tags geniunely have non-zero abundance in each cell.
#' 
#' Finally, if \code{controls} is specified, we remove cells with unusually low coverage of the controls, as defined with \code{\link{isOutlier}} on the log-sum count across control features.
#' Low control coverage is assumed to represent a failure of library preparation. 
#' 
#' @examples
#' x <- rbind(
#'     rpois(1000, rep(c(100, 10), c(100, 900))),
#'     rpois(1000, rep(c(20, 100, 20), c(100, 100, 800))),
#'     rpois(1000, rep(c(30, 100, 30), c(200, 700, 100)))
#' )
#'
#' # Adding a zero-ambient column plus a high-ambient column.
#' x <- cbind(0, x, 1000)
#' 
#' df <- cleanHashDrops(x)
#' df
#'
#' @seealso
#' \code{\link{ambientProfileBimodal}}, to infer the ambient profile.
#' 
#' \code{\link{ambientContribSparse}}, to estimate the ambient scaling for each droplet.
#'
#' \code{\link{isOutlier}}, to identify the outliers in a distribution of values.
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle medianSizeFactors isOutlier
cleanHashDrops <- function(x, ambient=NULL, controls=NULL, ...) {
    if (!is.null(ambient)) {
        ambient <- ambientProfileBimodal(x)
    }

    scale <- ambientContribSparse(x, ambient=ambient)
    non.zero <- scale > 0
    fail <- isOutlier(scale, type="higher", log=TRUE, subset=non.zero, ...)

    df <- DataFrame(
        ambient.scale=scale,
        zero.ambient=!non.zero, 
        high.ambient=fail,
        row.names=colnames(x)
    )
    discard <- !non.zero | fail

    if (!is.null(controls)) {
        df$sum.controls <- colSums(x[controls,,drop=FALSE])
        df$low.controls <- isOutlier(df$sum.controls, log=TRUE, type="lower", ...)
        discard <- discard | df$low.controls
    }

    df$discard <- discard
    df
}
