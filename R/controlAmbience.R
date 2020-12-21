#' Ambient contribution from controls
#'
#' Estimate the contribution of the ambient solution to a particular expression profile,
#' based on the abundance of control features that should not be expressed in the latter.
#'
#' @inheritParams maximumAmbience
#' @param features A logical, integer or character vector specifying the control features in \code{y} and \code{ambient}.
#' 
#' Alternatively, a list of vectors specifying mutually exclusive sets of features.
#'
#' @details
#' Control features should be those that cannot be expressed and thus fully attributable to ambient contamination.
#' This is most commonly determined \emph{a priori} from the biological context and experimental system.
#' For example, if spike-ins were introduced into the solution prior to cell capture,
#' these would serve as a gold standard for ambient contamination in \code{y}.
#' For single-nuclei sequencing, mitochondrial transcripts can serve a similar role 
#' under the assumption that all high-quality libraries are stripped nuclei.
#'
#' If \code{features} is a list, it is expected to contain multiple sets of mutually exclusive features.
#' These features need not be controls but each cell should only express features in one set (or no sets).
#' The expression of multiple sets can thus be attributed to ambient contamination.
#' For this mode, an archetypal pairing is that of hemoglobins with immunoglobulins (Young and Behjati, 2018), 
#' which should not be co-expressed in any (known) cell type.
#'
#' @return 
#' If \code{mode="scale"},
#' a numeric vector is returned quantifying the estimated \dQuote{contribution} of the ambient solution to each column of \code{y}.
#' Scaling columns of \code{ambient} by this vector yields the estimated ambient profile for each column of \code{y},
#' which can also be obtained by setting \code{mode="profile"}.
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
#' # Using the first 100 genes as a control:
#' scaling <- controlAmbience(y, ambient, features=1:100)
#' scaling
#'
#' # Estimating the control contribution to 'y' by 'ambient'.
#' contribution <- controlAmbience(y, ambient, features=1:100, mode="profile")
#' DataFrame(ambient=drop(contribution), total=y)
#'
#' @seealso 
#' \code{\link{estimateAmbience}}, to obtain an estimate to use in \code{ambient}.
#'
#' \code{\link{maximumAmbience}}, when control features are not available.
#'
#' @references
#' Young MD and Behjati S (2018).
#' SoupX removes ambient RNA contamination from droplet based single-cell RNA sequencing data.
#' \emph{biorXiv}.
#'
#' @export
#' @importFrom Matrix t
#' @importFrom scuttle sumCountsAcrossFeatures
controlAmbience <- function(y, ambient, features, mode=c("scale", "profile", "proportion")) {
    if (is.null(dim(y))) {
         y <- cbind(y)
    }
    if (is.null(dim(ambient))) {
        ambient <- matrix(ambient, nrow(y), ncol(y), dimnames=list(names(ambient), NULL))
    }
    if (!identical(rownames(y), rownames(ambient))) {
        warning("'y' and 'ambient' do not have the same row names")
    }

    if (!is.list(features)) {
        features <- list(features)
    } 
    sum.y <- sumCountsAcrossFeatures(y, features)
    sum.a <- sumCountsAcrossFeatures(ambient, features)
    props <- sum.y/sum.a
    scaling <- apply(props, 2, min)

    .report_ambient_profile(scaling, ambient=ambient, y=y, mode=match.arg(mode))
}
