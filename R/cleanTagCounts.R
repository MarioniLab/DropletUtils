#' Clean a tag-based dataset
#'
#' Remove low-quality libraries from a count matrix where each row is a tag and each column corresponds to a cell-containing barcode.
#'
#' @param x A numeric matrix-like object containing counts for each tag (row) in each cell (column).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param ambient A numeric vector of length equal to \code{nrow(x)}, containing the relative concentration of each tag in the ambient solution.
#' Defaults to \code{\link{ambientProfileBimodal}(x)} if not explicitly provided.
#' @param controls A vector specifying the rows of \code{x} corresponding to control tags.
#' These are expected to be isotype controls that should not exhibit any real binding.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment, further arguments to pass to the ANY method.
#'
#' For the ANY method, further arguments to pass to \code{\link{isOutlier}}.
#' This includes \code{batch} to account for multi-batch experiments, and \code{nmads} to specify the stringency of the outlier-based filter.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param exclusive A character vector of names of mutually exclusive tags that should never be expressed on the same cell.
#' Alternatively, a list of vectors of mutually exclusive sets of tags - see \code{\link{ambientContribNegative}} for details.
#' @param sparse.prop Numeric scalar specifying the minimum proportion of tags that should be present per cell.
#' 
#' @author Aaron Lun
#'
#' @return A \linkS4class{DataFrame} with one row per column of \code{x}, containing the following fields:
#' \itemize{
#' \item \code{zero.ambient}, a logical field indicating whether each cell has zero ambient contamination.
#' \item \code{sum.controls}, a numeric field containing the sum of counts for all control features.
#' Only present if \code{controls} is supplied.
#' \item \code{high.controls}, a logical field indicating whether each cell has unusually high control total.
#' Only present if \code{controls} is supplied.
#' \item \code{ambient.scale}, a numeric field specifying the relative amount of ambient contamination.
#' Only present if \code{controls} is \emph{not} supplied.
#' \item \code{high.ambient}, a numeric field indicating whether each cell has unusually high ambient contamination.
#' Only present if \code{controls} is \emph{not} supplied.
#' \item \code{discard}, a logical field indicating whether a column in \code{x} should be discarded.
#' }
#'
#' @details
#' We remove cells for which there is no detectable ambient contamination.
#' Specifically, we expect non-zero counts for most tags due to the deeply sequenced nature of tag-based data.
#' If \code{sparse.prop} or more tags have zero counts, this is indicative of a failure in library preparation for that cell.
#'
#' We also remove cells for which the total control count is unusually high.
#' The control coverage is used as a proxy for non-specific binding, most notably from contamination of droplets by protein aggregates.
#' High levels of non-specific activity are undesirable as this masks the actual marker profile of affected cells.
#' The upper threshold is defined with \code{\link{isOutlier}} on the log-total control count.
#'
#' If \code{controls} is missing, we instead compute the ambient scaling factor for each cell.
#' This represents the amount of ambient contamination - see \code{?\link{ambientContribSparse}} for more details -
#' and cells with unusually high values are assumed to be affected by protein aggregates.
#' High outliers are again identified and removed based on the log-ambient scale.
#'
#' If \code{controls} is missing and \code{exclusive} is specified, the ambient scaling factor is computed by \code{\link{ambientContribNegative}} instead.
#' This can be helpful for explicitly removing cells with impossible marker combinations,
#' though it is only as comprehensive as the knowledge of mutually exclusive marker sets.
#'
#' % In principle, it would be better to filter on the ratio of \code{ambient.scale} as a fraction of \code{colSums(x)}.
#' % This would better distinguish between antibody conjugates and droplets that just happen to be deeply sequenced,
#' % as the latter would have a low percentage of counts attributable to ambient contamination.
#' % In practice, this would also remove real cells that have none of the tags, which may actually be biologically interesting.
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
#' df <- cleanTagCounts(x)
#' df
#'
#' @seealso
#' \code{\link{ambientContribSparse}}, to estimate the ambient contamination for each droplet.
#'
#' \code{\link{isOutlier}}, to identify the outliers in a distribution of values.
#'
#' @name cleanTagCounts
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle isOutlier
#' @importFrom Matrix rowMeans colSums
.clean_tag_counts <- function(x, controls, ..., ambient=NULL, exclusive=NULL, sparse.prop=0.5) {
    if (is.null(ambient)) {
       ambient <- ambientProfileBimodal(x) 
    }

    scale <- ambientContribSparse(x, ambient=ambient, prop=sparse.prop)
    df <- DataFrame(
        zero.ambient=scale==0,
        row.names=colnames(x)
    )

    if (!missing(controls)) {
        df$sum.controls <- colSums(x[controls,,drop=FALSE])
        df$high.controls <- extra <- isOutlier(df$sum.controls, log=TRUE, type="higher", ...)
    } else {
        if (!is.null(exclusive)) {
            if (!is.list(exclusive)) {
                exclusive <- as.list(exclusive)
            }
            scale <- ambientContribNegative(x, ambient=ambient, features=exclusive)
        }
        df$ambient.scale <- scale
        df$high.ambient <- extra <- isOutlier(scale, log=TRUE, type="higher", ...)
    }

    df$discard <- df$zero.ambient | extra
    df
}

#' @export
#' @rdname cleanTagCounts
setGeneric("cleanTagCounts", function(x, ...) standardGeneric("cleanTagCounts"))

#' @export
#' @rdname cleanTagCounts
setMethod("cleanTagCounts", "ANY", .clean_tag_counts)

#' @export
#' @rdname cleanTagCounts
#' @importFrom SummarizedExperiment assay
setMethod("cleanTagCounts", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .clean_tag_counts(assay(x, assay.type), ...)
})
