#' Estimate the ambient profile
#'
#' Estimate the transcript proportions in the ambient solution using the Good-Turing method.
#'
#' @inheritParams emptyDrops
#' @param good.turing Logical scalar indicating whether to perform Good-Turing estimation of the proportions.
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
#' Good-Turing returns zero probabilities for zero counts if none of the summed counts are equal to 1.
#' This is technically correct but not helpful, so we protect against this by adding a single \dQuote{pseudo-feature} with a count of 1 to the profile.
#' The modified profile is used to calculate a Good-Turing estimate of observing any feature that has zero counts, which is then divided to get the per-feature probability. 
#' We scale down all the other probabilities to make space for this new pseudo-probability, which has some properties of unclear utility (see \url{https://github.com/MarioniLab/DropletUtils/issues/39}). 
#' 
#' Note that genes with counts of zero across all barcodes in \code{m} automatically have proportions of zero.
#' This ensures that the estimation is not affected by the presence/absence of non-expressed genes in different annotations.
#' In any case, such genes are likely to be completely irrelevant to downstream steps and can be safely ignored.
#'
#' Setting \code{good.turing=FALSE} may be convenient to obtain raw counts for use in further modelling.
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
#' ambience <- estimateAmbience(my.counts)
#' head(ambience)
#'
#' @seealso
#' \code{\link{emptyDrops}} and \code{\link{hashedDrops}}, where the ambient profile estimates are used for testing.
#'
#' @export
#' @importFrom Matrix rowSums colSums
estimateAmbience <- function(m, lower=100, round=TRUE, good.turing=TRUE) {
    m <- .rounded_to_integer(m, round)

    if (good.turing) {
        a <- .compute_ambient_stats(m, lower=lower)
        output <- numeric(nrow(m))
        output[!a$discard] <- a$ambient.prop
        names(output) <- rownames(m)
    } else {
        output <- rowSums(m[,colSums(m) <= lower,drop=FALSE])
    }

    output
}

#' @importFrom Matrix colSums rowSums
.compute_ambient_stats <- function(m, lower) {
    discard <- rowSums(m) == 0
    if (any(discard)) {
        m <- m[!discard,,drop=FALSE]
    }
    ncells <- ncol(m)

    # Computing the average profile from the ambient cells.
    # Enforcing discreteness mainly for emptyDrops()'s Monte Carlo step.
    umi.sum <- as.integer(round(colSums(m)))

    ambient <- umi.sum <= lower # lower => "T" in the text.
    ambient.m <- m[,ambient,drop=FALSE]
    ambient.prof <- rowSums(ambient.m)

    if (sum(ambient.prof)==0) {
        stop("no counts available to estimate the ambient profile")
    }
    ambient.prop <- .safe_good_turing(ambient.prof)

    list(
        m=m, # this MUST have the same number of columns as input.
        discard=discard,
        umi.sum=umi.sum,
        ambient=ambient,
        ambient.m=ambient.m,
        ambient.prop=ambient.prop
    )
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
