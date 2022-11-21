#' Estimate the ambient profile from empty droplets
#'
#' Estimate the transcript proportions in the ambient solution from an unfiltered count matrix,
#' assuming that low-count barcodes correspond to \dQuote{known empty} droplets.
#' Zeroes are filled in using the Good-Turing method.
#'
#' @inheritParams emptyDrops
#' @param by.rank An optional integer scalar, used as an alternative to \code{lower} to identifying assumed empty droplets - see Details.
#' @param good.turing Logical scalar indicating whether to perform Good-Turing estimation of the proportions.
#' @param known.empty an optional integer vector indexing barcodes that will be assumed to be empty, over-riding \code{lower} and \code{by.rank}.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{estimateAmbience}, arguments to pass to \code{ambientProfileEmpty}.
#'
#' @details
#' This function obtains an estimate of the composition of the ambient pool of RNA based on the barcodes with total UMI counts less than or equal to \code{lower}.
#' For each gene, counts are summed across all low-count barcodes and the resulting count vector is used for Good-Turing estimation of the proportions for each transcript.
#' The aim here is to obtain reasonable proportions for genes with counts of zero in low-count barcodes but non-zero counts for other barcodes (thus avoiding likelihoods of zero when modelling the latter with the proportions).
#'
#' This function will also attempt to detect whether \code{m} contains non-integer values by seeing if the column and row sums are discrete.
#' If such values are present, \code{m} is first \code{\link{round}}ed to the nearest integer value before proceeding.
#' This may be relevant when the count matrix is generated from pseudo-alignment methods like Alevin (see the \pkg{tximeta} package for details).
#' Rounding is performed by default as discrete count values are necessary for the Good-Turing algorithm, but if \code{m} is known to be discrete, setting \code{round=FALSE} can provide a small efficiency improvement.
#' 
#' Setting \code{good.turing=FALSE} may be convenient to obtain raw counts for use in further modelling.
#'
#' \code{estimateAmbience} is soft-deprecated; use \code{ambientProfileEmpty} instead.
#'
#' @section Behavior at zero counts:
#' Good-Turing returns zero probabilities for zero counts if none of the summed counts are equal to 1.
#' This is technically correct but not helpful, so we protect against this by adding a single \dQuote{pseudo-feature} with a count of 1 to the profile.
#' The modified profile is used to calculate a Good-Turing estimate of observing any feature that has zero counts, which is then divided to get the per-feature probability. 
#' We scale down all the other probabilities to make space for this new pseudo-probability, which has some properties of unclear utility (see \url{https://github.com/MarioniLab/DropletUtils/issues/39}). 
#' 
#' Note that genes with counts of zero across all barcodes in \code{m} automatically have proportions of zero.
#' This ensures that the estimation is not affected by the presence/absence of non-expressed genes in different annotations.
#' In any case, such genes are likely to be completely irrelevant to downstream steps and can be safely ignored.
#'
#' @section Finding the assumed empty droplets:
#' The default approach is to assume that all barcodes with total counts less than or equal to \code{lower} are empty.
#' This is generally effective but may not be adequate for datasets with unusually low or high sequencing depth, such that all or none of the barcodes are detected as empty respectively.
#' For example, there is no obvious choice for \code{lower} in CITE-seq data given that the coverage can be highly variable. 
#'
#' In such cases, an alternative approach can be used by passing an integer to the \code{by.rank} argument.
#' This specifies the number of barcodes with the highest total counts to ignore; the remaining barcodes are assumed to be ambient.
#' The idea is that, even if the exact threshold is unknown, we can be certain that a given experiment does not contain more than a particular number of genuine cell-containing barcodes based on the number of cells that were loaded into the machine.
#' By setting \code{by.rank} to something greater than this \emph{a priori} known number, we exclude the most likely candidates and use the remaining barcodes to compute the ambient profile.
#' 
#' Another alternative when working with some multimodal data, such as CITE-seq, could be to use statistics from one modality (e.g. mRNA counts) to define empty droplets in the other modality (e.g. CITE-seq) or combining CITE-seq with mRNA-counts.
#' In this case, one may set \code{known.empty} to an integer vector indexing barcodes in columns of `m` to mark cells for the ambient pool.
#' For the purpose of retaining cells, if \code{lower} is set, it will be used to define the pool of ambient RNA in \code{barcodeRanks}.
#' Otherwise the median of total counts of barcodes that have \code{known.empty} set will be used in its place.
#' 
#' @return
#' A numeric vector of length equal to \code{nrow(m)},
#' containing the estimated proportion of each transcript in the ambient solution.
#'
#' If \code{good.turing=FALSE}, the vector instead contains the sum of counts for each gene across all low-count barcodes.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts()
#' 
#' ambience <- ambientProfileEmpty(my.counts)
#' head(ambience)
#'
#' @seealso
#' \code{\link{emptyDrops}} and \code{\link{hashedDrops}}, where the ambient profile estimates are used for testing.
#'
#' \code{\link{ambientContribMaximum}} and related functions, to estimate the contribution of ambient contamination in each library.
#'
#' @name ambientProfileEmpty
NULL

#' @importFrom Matrix rowSums colSums
#' @importFrom BiocParallel bpstart bpstop SerialParam
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom DelayedArray setAutoBPPARAM
.ambient_profile_empty <- function(m, lower=100, by.rank=NULL, known.empty=NULL, round=TRUE, good.turing=TRUE, BPPARAM=SerialParam()) {
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    old <- .parallelize(BPPARAM)
    on.exit(setAutoBPPARAM(old), add=TRUE)

    m <- .rounded_to_integer(m, round)
    totals <- .intColSums(m)
    
    assumed.empty <- .get_putative_empty(totals, lower, by.rank, known.empty)

    if (good.turing) {
        a <- .compute_ambient_stats(m, totals, assumed.empty)
        output <- numeric(nrow(m))
        output[!a$discard] <- a$ambient.prop
        names(output) <- rownames(m)
    } else {
        output <- rowSums(m[,assumed.empty,drop=FALSE])
    }

    output
}

#' @export
#' @rdname ambientProfileEmpty
estimateAmbience <- function(...) {
    ambientProfileEmpty(...)
}

#' @importFrom Matrix rowSums
.compute_ambient_stats <- function(m, totals, assumed.empty) {
    # This doesn't invalidate 'totals', by definition.
    # NOTE: parallelization handled by setAutoBPPARAM above.
    discard <- rowSums(m) == 0
    if (any(discard)) {
        m <- m[!discard,,drop=FALSE]
    }

    # Computing the average profile from the ambient cells.
    ambient <- assumed.empty
    ambient.m <- m[,ambient,drop=FALSE]
    ambient.prof <- rowSums(ambient.m)

    if (sum(ambient.prof)==0) {
        stop("no counts available to estimate the ambient profile")
    }
    ambient.prop <- .safe_good_turing(ambient.prof)

    list(
        m=m, # this MUST have the same number of columns as input.
        discard=discard,
        ambient=ambient,
        ambient.m=ambient.m,
        ambient.prop=ambient.prop
    )
}

# inverse of `which(x)`
inv.which <- function(i, len, useNames = TRUE)
  stats::setNames('[<-'(logical(len), i, TRUE), if (useNames)
    names(i)
    else
      NULL)

.get_putative_empty <- function(totals, lower, by.rank, known.empty) {
    if (!is.null(known.empty)) {
      if (!is.integer(known.empty) || any(known.empty < 1) || any(known.empty > length(totals)) || any(duplicated(known.empty)))
        stop("If specified, 'known.empty' must be distinct, positive integers that index the barcodes.")
      if (!is.null(by.rank))
        warning("Ignoring 'by.rank' when 'known.empty' is set.")
      if (!is.null(lower))
        message(
          "When 'known.empty' is set, 'lower' is only used to mark droplets for unconditional retention via `barcodeRanks(...)`."
        )
    }
    if (is.null(by.rank) && is.null(known.empty)) {
      assumed.empty <- totals <= lower
    } else if (!is.null(known.empty)) {
      assumed.empty <- inv.which(known.empty, length(totals))
      if (is.null(lower))
        lower <- median(totals[assumed.empty])
    } else if (by.rank >= length(totals)) {
      stop("not have enough columns for supplied 'by.rank'")
    } else {
      #by.rank non null
      rt <- rank(-totals, ties.method = "first")
      assumed.empty <- rt > by.rank
      lower <- max(totals[rt > by.rank])
    }
    return(structure(assumed.empty, lower = lower))
}

#' @importFrom edgeR goodTuringProportions
.safe_good_turing <- function(ambient.prof) {
    ambient.prob <- goodTuringProportions(ambient.prof)

    still.zero <- ambient.prob<=0
    if (any(still.zero)) {
        pseudo.prob <- 1/sum(ambient.prof)
        ambient.prob[still.zero] <- pseudo.prob/sum(still.zero)
        ambient.prob[!still.zero] <- ambient.prob[!still.zero] * (1 - pseudo.prob)
    }

    ambient.prob
}

#' @export
#' @rdname ambientProfileEmpty
setGeneric("ambientProfileEmpty", function(m, ...) standardGeneric("ambientProfileEmpty"))

#' @export
#' @rdname ambientProfileEmpty
setMethod("ambientProfileEmpty", "ANY", .ambient_profile_empty)

#' @export
#' @rdname ambientProfileEmpty
#' @importFrom SummarizedExperiment assay
setMethod("ambientProfileEmpty", "SummarizedExperiment", function(m, ..., assay.type="counts") {
    .ambient_profile_empty(assay(m, assay.type), ...)    
})
