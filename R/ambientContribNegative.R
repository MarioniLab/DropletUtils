#' Ambient contribution from negative controls
#'
#' Estimate the contribution of the ambient solution to a particular expression profile,
#' based on the abundance of negative control features that should not be expressed in the latter.
#'
#' @param y A numeric matrix-like object containing counts, where each row represents a feature (e.g., a gene or a conjugated tag)
#' and each column represents either a cell or group of cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' 
#' \code{y} can also be a numeric vector of counts; this is coerced into a one-column matrix.
#' @param features A logical, integer or character vector specifying the negative control features in \code{y} and \code{ambient}.
#' 
#' Alternatively, a list of vectors specifying mutually exclusive sets of features.
#' @inheritParams ambientContribMaximum
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{controlAmbience}, arguments to pass to \code{ambientContribNegative}.
#'
#' @details
#' Negative control features should be those that cannot be expressed and thus fully attributable to ambient contamination.
#' This is most commonly determined \emph{a priori} from the biological context and experimental system.
#' For example, if spike-ins were introduced into the solution prior to cell capture,
#' these would serve as a gold standard for ambient contamination in \code{y}.
#' For single-nuclei sequencing, mitochondrial transcripts can serve a similar role 
#' under the assumption that all high-quality libraries are stripped nuclei.
#'
#' If \code{features} is a list, it is expected to contain multiple sets of mutually exclusive features.
#' Each cell should only express features in at most one set; no cell should express features in different sets.
#' The expression of multiple sets can thus be attributed to ambient contamination.
#' For this mode, an archetypal pairing is that of hemoglobins with immunoglobulins (Young and Behjati, 2018), 
#' which should not be co-expressed in any (known) cell type.
#'
#' \code{controlAmbience} is soft-deprecated; use \code{ambientContribNegative} instead.
#'
#' @return 
#' If \code{mode="scale"},
#' a numeric vector is returned quantifying the estimated \dQuote{contribution} of the ambient solution to each column of \code{y}.
#' Scaling \code{ambient} by each entry yields the maximum ambient profile for the corresponding column of \code{y}.
#'
#' If \code{mode="profile"}, a numeric matrix is returned containing the estimated ambient profile for each column of \code{y}.
#' This is computed by scaling as described above; if \code{ambient} is a matrix, each column is scaled by the corresponding entry of the scaling vector.
#'
#' If \code{mode="proportion"}, a numeric matrix is returned containing the estimated proportion of counts in \code{y} that are attributable to ambient contamination.
#' This is computed by simply dividing the output of \code{mode="profile"} by \code{y} and capping all values at 1.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Making up some data.
#' ambient <- c(runif(900, 0, 0.1), runif(100))
#' y <- rpois(1000, ambient * 50)
#' y <- y + c(integer(100), rpois(900, 5)) # actual biology, but first 100 genes silent.
#'
#' # Using the first 100 genes as negative controls:
#' scaling <- ambientContribNegative(y, ambient, features=1:100)
#' scaling
#'
#' # Estimating the negative control contribution to 'y' by 'ambient'.
#' contribution <- ambientContribNegative(y, ambient, features=1:100, mode="profile")
#' DataFrame(ambient=drop(contribution), total=y)
#'
#' @seealso 
#' \code{\link{ambientProfileEmpty}} or \code{\link{ambientProfileBimodal}}, to obtain a profile estimate to use in \code{ambient}.
#'
#' \code{\link{ambientContribMaximum}} or \code{\link{ambientContribSparse}},
#' for other methods of estimating contribution when negative control features are not available.
#'
#' @references
#' Young MD and Behjati S (2018).
#' SoupX removes ambient RNA contamination from droplet based single-cell RNA sequencing data.
#' \emph{biorXiv}.
#'
#' @name ambientContribNegative
NULL

#' @importFrom Matrix t
#' @importFrom scuttle sumCountsAcrossFeatures
#' @importFrom DelayedMatrixStats colMins
.ambient_contrib_negative <- function(y, ambient, features, mode=c("scale", "profile", "proportion")) {
    if (is.null(dim(y))) {
         y <- cbind(y)
    }

    if (!is.list(features)) {
        features <- list(features)
    } 
    sum.y <- sumCountsAcrossFeatures(y, features)

    if (is.null(dim(ambient))) {
        if (!identical(rownames(y), names(ambient))) {
            warning("'rownames(y)' and 'names(ambient)' are not the same")
        }
        sum.a <- sumCountsAcrossFeatures(cbind(ambient), features)
        sum.a <- drop(sum.a)
    } else {
        if (!identical(rownames(y), rownames(ambient))) {
            warning("'y' and 'ambient' do not have the same row names")
        }
        sum.a <- sumCountsAcrossFeatures(ambient, features)
    }

    props <- sum.y/sum.a
    if (nrow(props) == 1) {
        scaling <- drop(props)
    } else if (nrow(props) == 2) {
        # Special-cased for modest improvement to efficiency in the common
        # case of having two mutually exclusive sets.
        scaling <- colMins(props)
    } else {
        # Only one of these sets can be expressed, so the second-largest 
        # ratio must be attributable to ambient contamination (plus a
        # little bias from the ranking itself). 
        scaling <- apply(props, 2, function(x) sort(x, decreasing=TRUE)[2])
    }

    .report_ambient_profile(scaling, ambient=ambient, y=y, mode=match.arg(mode))
}

#' @export
#' @rdname ambientContribNegative
controlAmbience <- function(...) {
    ambientContribNegative(...)
}

#' @export
#' @rdname ambientContribNegative
setGeneric("ambientContribNegative", function(y, ...) standardGeneric("ambientContribNegative"))

#' @export 
#' @rdname ambientContribNegative
setMethod("ambientContribNegative", "ANY", .ambient_contrib_negative)

#' @export
#' @rdname ambientContribNegative
#' @importFrom SummarizedExperiment assay
setMethod("ambientContribNegative", "SummarizedExperiment", function(y, ..., assay.type="counts") {
    .ambient_contrib_negative(assay(y, assay.type), ...)
})

