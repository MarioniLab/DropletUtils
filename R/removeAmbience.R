#' Remove the ambient profile
#' 
#' Estimate and remove the ambient profile from a count matrix, given pre-existing groupings of similar cells.
#' This function is largely intended for plot beautification rather than real analysis.
#'
#' @param y A numeric matrix-like object containing counts for each gene (row) and cell or group of cells (column).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param ambient A numeric vector of length equal to \code{nrow(y)},
#' containing the proportions of transcripts for each gene in the ambient solution.
#' @param groups A vector of length equal to \code{ncol(y)}, specifying the assigned group for each cell.
#' This can also be a \linkS4class{DataFrame}, see \code{?\link{sumCountsAcrossCells}}.
#' @param features A vector of control features or a list of mutually exclusive feature sets, 
#' see \code{?\link{ambientContribNegative}} for more details.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the ANY method, Further arguments to pass to \code{\link{ambientContribMaximum}}.
#' @param assay.type Integer or string specifying the assay containing the count matrix.
#' @param dispersion Numeric scalar specifying the dispersion to use in the quantile-quantile mapping.
#' @param size.factors Numeric scalar specifying the size factors for each column of \code{y},
#' defaults to library size-derived size factors.
#' @param sink An optional \linkS4class{RealizationSink} object of the same dimensions as \code{y}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#'
#' @return
#' A numeric matrix-like object of the same dimensions as \code{y},
#' containing the counts after removing the ambient contamination.
#' The exact representation of the output will depend on the class of \code{y} and whether \code{sink} was used.
#'
#' @details
#' This function will aggregate counts from each group of related cells into an average profile.
#' For each group, we estimate the contribution of the ambient profile and subtract it from the average.
#' By default, this is done with \code{\link{ambientContribMaximum}}, but if enough is known about the biological system, users can specify \code{feaures} to use \code{\link{ambientContribNegative}} instead.
#'
#' We then perform quantile-quantile mapping of counts in \code{y} from the old to new averages.
#' This approach preserves the mean-variance relationship and improves the precision of estimate of the ambient contribution, but relies on a sensible grouping of similar cells, e.g., unsupervised clusters or cell type annotations.
#' As such, this function is best used at the end of the analysis to clean up expression matrices prior to visualization.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Making up some data.
#' ngenes <- 1000
#' ambient <- runif(ngenes, 0, 0.1)
#' cells <- c(runif(100) * 10, integer(900))
#' y <- matrix(rpois(ngenes * 100, cells + ambient), nrow=ngenes)
#'
#' # Pretending that all cells are in one group, in this example.
#' removed <- removeAmbience(y, ambient, groups=rep(1, ncol(y)))
#' summary(rowMeans(removed[1:100,]))
#' summary(rowMeans(removed[101:1000,]))
#'
#' @seealso
#' \code{\link{ambientContribMaximum}} and \code{\link{ambientContribNegative}}, to estimate the ambient contribution.
#'
#' \code{\link{estimateAmbience}}, to estimate the ambient profile.
#'
#' The \pkg{SoupX} package, which provides another implementation of the same general approach.
#' 
#' @name removeAmbience
NULL

#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment assay colData
#' @importFrom scuttle sumCountsAcrossCells librarySizeFactors
#' @importFrom BiocGenerics match
#' @importFrom DelayedArray blockApply colAutoGrid
.remove_ambience <- function(y, ambient, groups, features=NULL, ..., 
    size.factors=librarySizeFactors(y), dispersion=0.1, 
    sink=NULL, BPPARAM=SerialParam()) 
{
    summed <- sumCountsAcrossCells(y, groups, store.number=NULL)
    old.means <- assay(summed)
    sf.summed <- sumCountsAcrossCells(rbind(size.factors), groups)
    group.sf <- assay(sf.summed)[1,]

    if (is.null(features)) {
        profile <- ambientContribMaximum(old.means, ambient=ambient, ..., mode="profile")
    } else {
        profile <- ambientContribNegative(old.means, ambient=ambient, features=features, mode="profile")
    }
    new.means <- old.means - profile
    new.means[new.means < 0] <- 0

    if (is(groups, "DataFrame")) {
        matches <- match(groups, colData(summed))
    } else {
        matches <- match(groups, colData(summed)[,1])
    }

    if (!is.null(sink)) {
        # Force serial execution to avoid multi-threaded writing to sink.
        BPPARAM <- NULL
    }

    out <- blockApply(y, FUN=.remap_quantiles, 
         matches=matches, old.means=old.means, new.means=new.means, 
         dispersion=dispersion, size.factors=size.factors, group.sf=group.sf,
         sink=sink, 
         BPPARAM=BPPARAM, as.sparse=NA, grid=colAutoGrid(y))

    if (!is.null(sink)) {
        as(sink, "DelayedArray")
    } else {
        do.call(cbind, out)
    }
}

#' @importFrom methods is
#' @importFrom stats pnbinom qnbinom
#' @importFrom edgeR q2qnbinom
#' @importFrom DelayedArray makeNindexFromArrayViewport write_block currentViewport
.remap_quantiles <- function(block, matches, old.means, new.means, size.factors, group.sf, dispersion, sink) {
    vp <- currentViewport()
    idx <- makeNindexFromArrayViewport(vp, expand.RangeNSBS = TRUE)
    if (!is.null(idx[[2]])) {
        matches <- matches[idx[[2]]]
        size.factors <- size.factors[idx[[2]]]
    }

    old.means <- old.means[, matches, drop=FALSE]
    new.means <- new.means[, matches, drop=FALSE]
    rescale <- size.factors / group.sf[matches]
    old.means <- t(t(old.means) * rescale)
    new.means <- t(t(new.means) * rescale)

    originals <- as.matrix(block)
    p <- pnbinom(originals, mu=old.means, size=1/dispersion, log.p=TRUE)
    q <- qnbinom(p, mu=new.means, size=1/dispersion, log.p=TRUE)

    failed <- !is.finite(q)
    fallback <- q2qnbinom(originals[failed], input.mean=old.means[failed], output.mean=new.means[failed], dispersion=dispersion)
    fallback <- pmax(0, round(fallback))
    q[failed] <- fallback

    dim(q) <- dim(originals)

    if (!is.null(sink)) {
        write_block(sink, vp, q)
        NULL
    } else if (is(block, "SparseArraySeed")) {
        as(q, "CsparseMatrix")
    } else {
        q 
    }
}

#' @export
#' @rdname removeAmbience
setGeneric("removeAmbience", function(y, ...) standardGeneric("removeAmbience"))

#' @export
#' @rdname removeAmbience
setMethod("removeAmbience", "ANY", .remove_ambience)

#' @export
#' @rdname removeAmbience
#' @importFrom SummarizedExperiment assay
setMethod("removeAmbience", "SummarizedExperiment", function(y, ..., assay.type="counts") {
    .remove_ambience(assay(y, assay.type), ...)
})
