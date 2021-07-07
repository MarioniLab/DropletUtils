#' Ambient contribution by assuming sparsity 
#'
#' Estimate the contribution of the ambient solution to each droplet by assuming that no more than a certain percentage of features are actually present/expressed in the cell.
#'
#' @param y A numeric matrix-like object containing counts, where each row represents a feature (usually a conjugated tag) and each column represents a cell or group of cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' 
#' \code{y} can also be a numeric vector of counts; this is coerced into a one-column matrix.
#' @param ambient A numeric vector of length equal to \code{nrow(y)},
#' containing the proportions of transcripts for each feature in the ambient solution.
#' @param prop Numeric scalar specifying the maximum proportion of features that are expected to be present for any cell.
#' @inheritParams ambientContribMaximum
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' 
#' @return
#' If \code{mode="scale"}, a numeric vector is returned quantifying the estimated \dQuote{contribution} of the ambient solution to each column of \code{y}.
#' Scaling \code{ambient} by each entry yields the estimated ambient profile for the corresponding column of \code{y}.
#'
#' If \code{mode="profile"}, a numeric matrix is returned containing the estimated ambient profile for each column of \code{y}.
#' This is computed by scaling as described above; if \code{ambient} is a matrix, each column is scaled by the corresponding entry of the scaling vector.
#'
#' If \code{mode="proportion"}, a numeric matrix is returned containing the proportion of counts in \code{y} that are attributable to ambient contamination.
#' This is computed by simply dividing the output of \code{mode="profile"} by \code{y} and capping all values at 1.
#'
#' @details
#' The assumption here is that of sparsity, i.e., no more than \code{prop * nrow(y)} features should be actually present in each cell with a non-zero number of molecules.
#' This is reasonable for most tag-based applications where we would expect only 1-2 tags (for cell hashing) or a minority of tags (for general CITE-seq) to be present per cell.
#' Thus, counts for all other features must be driven by ambient contamination, allowing us to estimate a scaling factor for each cell based on the ratio to the ambient profile.
#'
#' For gene expression, the sparsity assumption is less justifiable.
#' Each cell could feasibly express a majority of the transcriptome (once we ignore constitutively silent features in the annotation, like pseudogenes).
#' The sparsity of gene expression data also yields less precise scale estimates, reducing their utility in downstream applications.
#' See \code{\link{ambientContribNegative}} or \code{\link{ambientContribMaximum}} instead, which operate from different assumptions.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{ambientProfileBimodal}}, to estimate the ambient profile for use in \code{ambient}.
#'
#' \code{\link{cleanTagCounts}}, where this function is used to estimate ambient scaling factors.
#' 
#' @examples
#' amb <- 1:10
#' y <- matrix(rpois(10000, lambda=amb), nrow=10)
#' y[sample(length(y), 1000, replace=TRUE)] <- 1000
#' 
#' scaling <- ambientContribSparse(y, ambient=amb)
#' hist(scaling)
#' 
#' @name ambientContribSparse
NULL

#' @importFrom DelayedMatrixStats colQuantiles
#' @importFrom DelayedArray setAutoBPPARAM
#' @importFrom BiocParallel SerialParam
.ambient_contrib_sparse <- function(y, ambient, prop=0.5, 
    mode=c("scale", "profile", "proportion"), BPPARAM=SerialParam())
{
    if (is.null(dim(y))) {
        y <- cbind(y)
    }

    old <- .parallelize(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    ambient[ambient==0] <- NA_real_ # trigger the formation of NAs.
    scaling <- colQuantiles(DelayedArray(y/ambient), na.rm=TRUE, probs=1-prop)
    .report_ambient_profile(scaling, ambient=ambient, y=y, mode=match.arg(mode))
}

#' @export
#' @rdname ambientContribSparse
setGeneric("ambientContribSparse", function(y, ...) standardGeneric("ambientContribSparse"))

#' @export 
#' @rdname ambientContribSparse
setMethod("ambientContribSparse", "ANY", .ambient_contrib_sparse)

#' @export
#' @rdname ambientContribSparse
#' @importFrom SummarizedExperiment assay
setMethod("ambientContribSparse", "SummarizedExperiment", function(y, ..., assay.type="counts") {
    .ambient_contrib_sparse(assay(y, assay.type), ...)
})
