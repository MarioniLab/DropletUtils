#' An approximate implementation of the \code{--soloCellFilter  EmptyDrops_CR} filtering approach to identify empty droplets.
#'
#' An approximate implementation of the \code{--soloCellFilter  EmptyDrops_CR} filtering approach, 
#' which, itself, was reverse-engineered from the behavior of CellRanger 3+.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for barcodes.
#' 
#' @param expected A numeric scalar specifying the expected number of barcodes in this sample It is same as \code{nExpectedCells} in STARsolo.
#' 
#' @param max.percentile A numeric scalar specifying the percentile used in simple filtering, barcodes selected by simple filtering 
#' will be regarded as real barcodes regardless of the \code{emptyDrops} result. It is same as \code{maxPercentile} in STARsolo.
#' 
#' @param max.min.ratio A numeric scalar specifying the maximum ratio of maximum UMI count and minimum UMI count used in simple filtering, 
#' maximum UMI count used in simple filtering is determined first by \code{expected*(1-max.percentile)}, minimum UMI count used in
#'  simple filtering is then determined by this ratio. It is same as \code{maxMinRatio} in STARsolo.
#'  
#' @param umi.min A numeric scalar specifying the minimum UMI count above which a sample will be included in ambient profiles. It is same as \code{umiMin} in STARsolo.
#' 
#' @param umi.min.frac.median A numeric scalar between 0 and 1 specifying that only the barcodes whose UMI count are above this number 
#' fraction of the median UMI count of the top \code{expected} barcodes will be included in the ambient profile. 
#' It is same as \code{umiMinFracMedian} in STARsolo.
#' 
#' @param cand.max.n An integer specifying the maximum number of ambient barcodes that are possible to be regarded as real barcodes. It is same as \code{canMaxN} in STARsolo.
#' 
#' @param by.rank.upper An integer specifying the lowest UMI count ranking of the ambient pool, barcodes with UMI count ranking below
#' this number will not be included in the ambient pool. It is same as \code{indMin} in STARsolo.

#' @param by.rank.lower An integer specifying the highest UMI count ranking of the ambient pool, barcodes with UMI count ranking above
#' this number will not be included in the ambient pool. It is same as \code{indMax} in STARsolo.
#' 
#' @param fdr.thresh A numeric scalar specifying the FDR threshold to filter barcodes. Barcodes whose FDR returned by emptyDrops
#' is above this threshold will not be regarded as real barcodes. It is same as \code{FDR} in STARsolo.
#'
#' @param barcode.args Further arguments to pass to \code{\link{barcodeRanks}}.
#' @param round Logical scalar indicating whether to check for non-integer values in \code{m} and, if present, round them for ambient profile estimation (see \code{?\link{ambientProfileEmpty}}) and the multinomial simulations.
#' @param alpha A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme (see \code{?\link{emptyDrops}}).
#' @param ... For the generic, further arguments to pass to to individual methods.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be used.

#' @details
#' This function is an approximate implementation of the  \code{--soloCellFilter  EmptyDrops_CR} filtering approach of STARsolo 2.7.9a
#' , which, itself, was reverse engineered from the behavior of  CellRanger 3+. 
#' All parameters are defaulty set as the default value used in starSolo 2.7.9a.
#' In the most cases, users just need to specify the raw and unfiltered count matrix, \code{m}.
#' See \code{?\link{emptyDrops}} for an alternative approach for cell calling.
#' 
#' @return
#' A DataFrame like \code{\link{emptyDrops}}, with an additional binary \code{is.cell} field demonstrating whether
#' barcodes are estimated as real barcodes.
#' 
#' @author
#' Dongze He, Rob Patro
#' 
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts()
#' 
#' # Identify likely cell-containing droplets.
#' e.out <- cellRangerLikeEmptyDrops(my.counts)
#' e.out
#' 
#' # Get matrix of estimated barcodes.
#' cell.counts <- my.counts[, e.out$is.cell]
#' 
#' @references
#' Kaminow et al. (2021).
#' STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data
#' \url{https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1}
#' 
#' @seealso
#' \code{\link{emptyDrops}}, for another method for calling barcodes.
#'
#' @name cellRangerLikeEmptyDrops
NULL

# Authors: Dongze He, Rob Patro
# Center of Bioinformatics and Computational Biology, University of Maryland, College Park, Maryland, 20740

#' @export
#' @rdname cellRangerLikeEmptyDrops
#' @importFrom BiocParallel bpstart bpstop SerialParam
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom Matrix colSums
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom beachmat colBlockApply
.test_empty_drops_cr <- function(m, lower=NULL, upper=NULL, niters=10000, ignore=NULL, retain=NULL, alpha=Inf, 
                           round=TRUE, by.rank.upper=45000, by.rank.lower=90000, BPPARAM=SerialParam()) 
{
  if (!.bpNotSharedOrUp(BPPARAM)) {
    bpstart(BPPARAM)
    on.exit(bpstop(BPPARAM))
  }
  
  old <- .parallelize(BPPARAM)
  on.exit(setAutoBPPARAM(old), add=TRUE)
  
  # NOTE: this is probably fine, as you don't have that many cell-containing
  # droplets per sample, so the sparse matrix will generally be small.
  m <- .realize_DA_to_memory(m, BPPARAM)
  
  m <- .rounded_to_integer(m, round)
  totals <- .intColSums(m)
  lower <- .get_lower_cr(totals, lower, by.rank.lower=by.rank.lower)
  upper <- .get_upper_cr(totals, upper, by.rank.upper=by.rank.upper)
  
  astats <- .compute_ambient_stats_cr(m, totals, lower=lower, upper=upper)

  m <- astats$m
  ambient <- astats$ambient
  ambient.prop <- astats$ambient.prop
  
  # Estimating the alpha from the discarded ambient droplets, if desired.
  if (is.null(alpha)) {
    ambient.m <- astats$ambient.m
    ambient.totals <- totals[ambient]
    alpha <- .estimate_alpha(ambient.m, ambient.prop, ambient.totals) 
  }
  
  # Removing supposed ambient barcodes from the matrix.
  # Also removing additional barcodes that don't pass some total count threshold.
  keep <- totals > 0L
  
  if (!is.null(ignore)) { 
    keep <- keep & totals >= ignore
  }
  if (!is.null(retain)) {
    keep <- keep & totals < retain
  }
  obs.m <- m[,keep,drop=FALSE]
  obs.totals <- totals[keep]
  
  # Calculating the log-multinomial probability for each cell.
  obs.P <- colBlockApply(obs.m, prop=ambient.prop, alpha=alpha, BPPARAM=BPPARAM, FUN=.compute_multinom_prob_data)
  obs.P <- unlist(obs.P, use.names=FALSE)
  rest.P <- .compute_multinom_prob_rest(obs.totals, alpha=alpha)
  
  # Computing the p-value for each observed probability.
  n.above <- .permute_counter(totals=obs.totals, probs=obs.P, 
                              ambient=ambient.prop, iter=niters, BPPARAM=BPPARAM, alpha=alpha)
  limited <- n.above==0L
  pval <- (n.above+1)/(niters+1)
  
  # Collating into some sensible output.
  ncells <- ncol(m)
  all.p <- all.lr <- all.exp <- rep(NA_real_, ncells)
  all.lim <- rep(NA, ncells)
  
  all.p[keep] <- pval
  all.lr[keep] <- obs.P + rest.P 
  all.lim[keep] <- limited
  
  output <- DataFrame(Total=totals, LogProb=all.lr, PValue=all.p, Limited=all.lim, row.names=colnames(m))
  metadata(output) <- list(lower=lower, upper=upper, niters=niters, ambient=ambient.prop, alpha=alpha)
  output
}

.get_lower_cr <- function(totals, lower, by.rank.lower) {
  if (!is.null(lower)) {
    lower
  } else if (by.rank.lower >= length(totals)) {
    stop("not have enough columns for supplied 'by.rank.upper'")
  } else {
    totals[order(totals, decreasing=TRUE)[by.rank.lower+1]]
  }
}

.get_upper_cr <- function(totals, upper, by.rank.upper) {
  if (!is.null(upper)) {
    upper
  } else if (by.rank.upper >= length(totals)) {
    totals[order(totals, decreasing=TRUE)[length(totals)]]
  } else {
    totals[order(totals, decreasing=TRUE)[by.rank.upper]]
  }
}

#' @importFrom Matrix rowSums
.compute_ambient_stats_cr <- function(m, totals, lower, upper) {
  # This doesn't invalidate 'totals', by definition.
  # NOTE: parallelization handled by setAutoBPPARAM above.
  discard <- rowSums(m) == 0
  if (any(discard)) {
    m <- m[!discard,,drop=FALSE]
  }
  
  # Computing the average profile from the ambient barcodes.
  ambient <- totals >= lower & totals <= upper # lower => "T" in the text.

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


#' @importFrom stats p.adjust
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom BiocParallel SerialParam
.cell_ranger_like_empty_drops  <- function(m, 
                                           # STARsolo arguments
                                           ## simple filtering
                                           expected=3000,            # nExpectedCells
                                           max.percentile=0.99,      # maxPercentile
                                           max.min.ratio=10,         # maxMinRatio
                                           ## emptyDrops_CR
                                           umi.min=500,              # umiMin
                                           umi.min.frac.median=0.01, # umiMinFracMedian
                                           can.max.n=20000,          # candMaxN
                                           by.rank.upper=45000,      # indMin
                                           by.rank.lower=90000,      # indMax
                                           fdr.threshold=0.01,       # FDR
                                           # emptyDrops arguments
                                           barcode.args=list(),
                                           round=TRUE,
                                           alpha=Inf, 
                                           ...,
                                           BPPARAM=SerialParam()
) {
  
  # This function is an approximate implementation of the 
  # `--soloCellFilter  EmptyDrops_CR` filtering approach 
  # of STARsolo 2.7.9a (https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1),
  # which, itself, was reverse engineered from the behavior of 
  # CellRanger 3+. The original C++ code on which this 
  # function is based can be found at (https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp) 
  
  ####################################################################################################################
  # setup emptyDrops
  if (!.bpNotSharedOrUp(BPPARAM)) {
    bpstart(BPPARAM)
    on.exit(bpstop(BPPARAM))
  }
  
  m <- .realize_DA_to_memory(m, BPPARAM)
  m <- .rounded_to_integer(m, round)

  ###################################################################################################################    
  # Simple Filtering
  # https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
  # line 36-61
  totals <- .intColSums(m)
  
  max.ind <- round(expected * (1 - max.percentile))
  n.umi.max <- totals[order(totals, decreasing = TRUE)[min(length(totals), max.ind)]]
  ## retain: barcodes with UMI count higher than retain will be regarded as real barcodes without any further tests
  retain <- max(round(n.umi.max/max.min.ratio),1) # n.umi.min
  ncells.simple <- sum(totals>=retain)
  
  ###################################################################################################################    
  # select candidate barcodes
  # SoloFeature_emptyDrops_CR.cpp line 117-134
  
  i.cand.first = ncells.simple+1
  min.umi  <- max(umi.min, round(umi.min.frac.median * totals[ncells.simple/2]))
  i.cand.last = min(ncells.simple + can.max.n, sum(totals > min.umi))

  ## ignore: barcodes with UMI count higher than retain will be regarded as emptyDrops,
  ## but barcodes with UMI counts between [ignore, lower] will be used in SimpleGoodTuring
  ignore = totals[order(totals, decreasing = TRUE)[i.cand.last]]
  
  stats <- .test_empty_drops_cr(m, 
                            round=FALSE, 
                            ignore=ignore,
                            retain=retain, 
                            by.rank.lower=by.rank.lower, 
                            by.rank.upper=by.rank.upper, 
                            ...,
                            BPPARAM=BPPARAM) 

  tmp <- stats$PValue
  
  # Possibly redefine 'lower' based on 'by.rank=' passed to testEmptyDrops.
  lower <- metadata(stats)$lower
  upper <- metadata(stats)$upper
  
  metadata(stats)$retain <- retain
  always <- stats$Total >= retain
  tmp[always] <- 0
  
  discard <- (stats$Total <= ignore)
  tmp[discard] <- NA_real_
  
  stats$FDR <- p.adjust(tmp, method="BH")
  stats$is.cell = stats$FDR <= fdr.threshold
  stats$is.cell[is.na(stats$is.cell)] <- FALSE
  stats
}

#' @export
#' @rdname cellRangerLikeEmptyDrops
setGeneric("cellRangerLikeEmptyDrops", function(m, ...) standardGeneric("cellRangerLikeEmptyDrops"))

#' @export
#' @rdname cellRangerLikeEmptyDrops
setMethod("cellRangerLikeEmptyDrops", "ANY", .cell_ranger_like_empty_drops)

#' @export
#' @rdname cellRangerLikeEmptyDrops
#' @importFrom SummarizedExperiment assay
setMethod("cellRangerLikeEmptyDrops", "SummarizedExperiment", function(m, ..., assay.type="counts") {
  .cell_ranger_like_empty_drops(assay(m, assay.type), ...)
})
