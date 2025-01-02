#' CellRanger's emptyDrops variant 
#'
#' An approximate implementation of the \code{--soloCellFilter EmptyDrops_CR} filtering approach, 
#' which itself was reverse-engineered from the behavior of CellRanger 3.
#' 
#' @param m A numeric matrix-like object containing counts, where columns represent barcoded droplets and rows represent features.
#' The matrix should only contain barcodes for an individual sample, prior to any filtering for barcodes.
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such an object.
#' @param n.expected.cells An integer scalar specifying the number of expected cells in a sample. 
#' If missing, will try to estimate this from the data using the order of magnitude algorithm from CellRanger.
#' Corresponds to the \code{nExpectedCells} argument in \pkg{STARsolo}. 
#' @param max.percentile A numeric scalar between 0 and 1 used to define the maximum UMI count in the simple filtering algorithm. 
#' Corresponds to the \code{maxPercentile} argument in \pkg{STARsolo}. 
#' @param max.min.ratio An integer scalar specifying the ratio of the maximum and minimum UMI count in the simple filtering algorithm. 
#' Corresponds to the \code{maxMinRatio} argument in \pkg{STARsolo}.
#' @param umi.min An integer scalar specifying the minimum UMI count for inclusion of a barcode in the cell candidate pool. 
#' Corresponds to the \code{umiMin} argument in \pkg{STARsolo}.
#' @param umi.min.frac.median A numeric scalar between 0 and 1 used to define the minimum UMI count for inclusion of a barcode in the cell candidate pool.
#' Specifically, the minimum is defined as \code{umi.min.frac.median} times the median UMI count of the real cells assigned by the simple filtering algorithm. 
#' Corresponds to the \code{umiMinFracMedian} argument in \pkg{STARsolo}.
#' @param cand.max.n An integer scalar specifying the maximum number of barcodes that can be included in the cell candidate pool. 
#' In effect, this applies a minimum threshold that is defined as the \code{cand.max.n}-th largest UMI count among all cells that are \emph{not} selected by the simple filtering algorithm. 
#' Corresponds to the \code{candMaxN} in \pkg{STARsolo}.
#' @param ind.min An integer scalar specifying the lowest UMI count ranking for inclusion of a barcode in the ambient profile. 
#' Corresponds to the \code{indMin} argument in \pkg{STARsolo}.
#' @param ind.max An integer scalar specifying the highest UMI count ranking for inclusion of a barcode in the ambient profile. 
#' Corresponds to the \code{indMin} argument in \pkg{STARsolo}.
#' @param round A logical scalar indicating whether to check for non-integer values in \code{m} and, if present, round them for ambient profile estimation (see \code{?\link{ambientProfileEmpty}}) and the multinomial simulations.
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be used.
#' @param ... Further arguments to pass to individual methods.
#' Specifically, for the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.type String or integer specifying the assay of interest.
#'
#' @details
#' \code{emptyDropsCellRanger} splits each sample's barcodes into three subsets.
#' \enumerate{
#' \item The first subset contains barcodes that are selected by the \dQuote{simple filtering algorithm}, which are regarded as high quality cells without any further filtering.
#' The minimum threshold \eqn{T} for this subset is defined by taking the \code{max.percentile} percentile of the top \code{n.expected.cells} barcodes,
#' and then dividing by the \code{max.min.ratio} to obtain a minimum UMI count.
#' (This is closely related to the algorithm used by \code{\link{defaultDrops}}.)
#' All barcodes identified in this manner will have an FDR of zero.
#' \item The second subset contains the ambient pool and is defined as all barcodes with rankings between \code{ind.min} and \code{ind.max}. 
#' The barcodes that fall in this category will be used to compute the ambient profile.
#' None of these barcodes are considered to be potential cells.
#' \item The third subset contains the pool of barcodes that are potential cells, i.e., cell candidates.
#' This is defined as all barcodes with total counts below \eqn{T} and higher than all of the thresholds defined by \code{umi.min}, \code{umi.min.frac.median} and \code{cand.max.n}.
#' Only the barcodes within this subset will be tested for signficant deviations from the ambient profile, i.e., FDR is not \code{NaN}.
#' }
#'
#' As of time of writing, the arguments in \pkg{STARsolo} have a one-to-one correspondence with the arguments in \code{emptyDropsCellRanger}. 
#' All parameter defaults are set as the same as those used in STARsolo 2.7.9a.
#' 
#' The main differences between \code{emptyDropsCellRanger} and \code{emptyDrops} are:
#' \itemize{
#' \item \code{emptyDropsCellRanger} does not use the knee point to identify \dQuote{presumed real} cells,
#' instead relying on a threshold based on the expected number of cells.
#' \item \code{emptyDropsCellRanger} takes barcodes whose total count ranks within a certain range - by default, \eqn{(45,000, 90,000]} - to compute the ambient profile.
#' In contrast, \code{emptyDrops} only defines the upper bound using \code{lower} or \code{by.rank}.
#' \item \code{emptyDropsCellRanger} defines a cell candidate pool according to three parameters, \code{umi.min}, \code{umi.min.frac.median} and \code{cand.max.n}.
#' In \code{emptyDrops}, this is only defined by \code{lower}.
#' }
#' 
#' @return
#' A \linkS4class{DataFrame} with the same fields as that returned by \code{\link{emptyDrops}}.
#' 
#' @author
#' Dongze He, Rob Patro
#' 
#' @examples
#' # Mocking up some data:
#' set.seed(0)
#' my.counts <- DropletUtils:::simCounts(nempty=100000, nlarge=2000, nsmall=1000)
#' 
#' # Identify likely cell-containing droplets.
#' out <- emptyDropsCellRanger(my.counts)
#' out
#'
#' is.cell <- out$FDR <= 0.01
#' sum(is.cell, na.rm=TRUE)
#'
#' # Subsetting the matrix to the cell-containing droplets.
#' # (using 'which()' to handle NAs smoothly).
#' cell.counts <- my.counts[,which(is.cell),drop=FALSE]
#' dim(cell.counts)
#' 
#' @references
#' Kaminow B, Yunusov D, Dobin A (2021).
#' STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data.
#' \url{https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1}
#' 
#' @seealso
#' \code{\link{emptyDrops}}, for the original implementation.
#'
#' @name emptyDropsCellRanger
NULL

#' @importFrom stats p.adjust
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom BiocParallel SerialParam
.empty_drops_cell_ranger  <- function(m, 
                                      # STARsolo arguments
                                      ## simple filtering
                                      n.expected.cells=NULL,            # nExpectedCells
                                      max.percentile=0.99,      # maxPercentile
                                      max.min.ratio=10,         # maxMinRatio
                                      
                                      ## emptyDrops_CR
                                      umi.min=500,              # umiMin
                                      umi.min.frac.median=0.01, # umiMinFracMedian
                                      cand.max.n=20000,         # candMaxN
                                      ind.min=45000,            # indMin
                                      ind.max=90000,            # indMax

                                      # emptyDrops arguments
                                      round=TRUE,
                                      niters=10000,
                                      BPPARAM=SerialParam())
{  
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # This function is an approximate implementation of the 
    # `--soloCellFilter  EmptyDrops_CR` filtering approach 
    # of STARsolo 2.7.9a (https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1),
    # which, itself, was reverse engineered from the behavior of 
    # CellRanger 3+. The original C++ code on which this 
    # function is based can be found at (https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
    #  and https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_emptyDrops_CR.cpp) 
    
    ambfun <- function(mat, totals) {
        o <- order(totals, decreasing=TRUE)

        # Simple Filtering
        # https://github.com/alexdobin/STAR/blob/master/source/SoloFeature_cellFiltering.cpp
        # line 36-61
        max.ind <- round(n.expected.cells * (1 - max.percentile)) # maxind
        n.umi.max <- totals[o[min(length(totals), max.ind)]] # nUMImax
        
        # Barcodes with UMI count higher than retain will be regarded as real
        # barcodes without any further tests
        retain <- max(round(n.umi.max/max.min.ratio),1) # nUMImin
        
        # select barcodes to use as ambient solution
        # SoloFeature_emptyDrops_CR.cpp line 117-134
        ncells.simple <- sum(totals >= retain)
        min.umi  <- max(umi.min, round(umi.min.frac.median * totals[o[ncells.simple/2]]))
        i.cand.last <- min(ncells.simple + cand.max.n, sum(totals > min.umi))

        # NOTE: parallelization handled by setAutoBPPARAM above.
        discard <- rowSums(mat) == 0
        if (any(discard)) {
            mat <- mat[!discard,,drop=FALSE]
        }
        
        ambient <- logical(length(totals))
        final.ind.min = min(ind.min, length(totals))
        final.ind.max = min(ind.max, length(totals))
        if ((final.ind.max - final.ind.min) < 100) {
            stop(paste0(
                "The ambient pool size (",
                final.ind.max - final.ind.min + 1,
                ") is too small; cannot proceed.",
                "Please adjust `ind.min` and `ind.max` to increase the size. ",
                "One suggestion is to set `ind.min = sum(colSums(m) > N)` and `ind.max = ncol(m)`, ",
                "where N is a threshold of UMI count, such as 100, and m is the count matrix."))
        }
        ambient[o[final.ind.min:final.ind.max]] <- TRUE
        ambient.m <- mat[,ambient,drop=FALSE]
        ambient.prof <- rowSums(ambient.m)
        
        if (sum(ambient.prof)==0) {
            stop("no counts available to estimate the ambient profile")
        }
        ambient.prop <- .safe_good_turing(ambient.prof)
        
        # Barcodes to keep are not just !ambient, as we want to exclude the
        # barcodes that are too low to be even considered as ambient.
        # keep <- (o < (ncells.simple+1))
        keep <- logical(length(totals)) 
        keep[o[(ncells.simple+1):i.cand.last]] = TRUE
        
        list(
            m=mat, # this MUST have the same number of columns as input.
            discard=discard,
            ambient=ambient,
            ambient.m=ambient.m,
            ambient.prop=ambient.prop,
            keep=keep,
            metadata=list(lower = totals[o[min(ind.min, length(totals))]], bottom=totals[o[min(ind.max, length(totals))]], retain=retain)
        )
    }
    
    if (is.null(n.expected.cells)) {
        n.expected.cells <- filter_cellular_barcodes_ordmag(colSums(m))$filtered_bcs
        if (n.expected.cells < 1) {
            warning("Could not estimate the number of expected cells; using 3000 instead.")
            n.expected.cells <- 3000
        }
    }

    stats <- .test_empty_drops(m=m, ambient.FUN=ambfun, niters=niters, test.ambient=FALSE, ignore=NULL, alpha=Inf, round=round, BPPARAM=BPPARAM)
    
    tmp <- stats$PValue
    retain <- metadata(stats)$retain
    always <- stats$Total >= retain
    tmp[always] <- 0
    
    stats$FDR <- p.adjust(tmp, method="BH")
    stats
}

#' @export
#' @rdname emptyDropsCellRanger
setGeneric("emptyDropsCellRanger", function(m, ...) standardGeneric("emptyDropsCellRanger"))

#' @export
#' @rdname emptyDropsCellRanger
setMethod("emptyDropsCellRanger", "ANY", .empty_drops_cell_ranger)

#' @export
#' @rdname emptyDropsCellRanger
#' @importFrom SummarizedExperiment assay
setMethod("emptyDropsCellRanger", "SummarizedExperiment", function(m, ..., assay.type="counts") {
    .empty_drops_cell_ranger(assay(m, assay.type), ...)
})

# The following functions are adopted from CellRanger
# https://github.com/10XGenomics/cellranger/blob/fe0616967771151209c3a6f97f4de92bf55ac0bd/lib/python/cellranger/cell_calling_helpers.py#L841
## Cell calling constants
ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100
ORDMAG_RECOVERED_CELLS_QUANTILE = 0.99
MIN_RECOVERED_CELLS_PER_GEM_GROUP = 50
MAX_RECOVERED_CELLS_PER_GEM_GROUP = 262144 # 1 << 18

find_within_ordmag <- function(x, baseline_idx) {
    x_ascending = sort(x)
    baseline = x_ascending[length(x_ascending) - baseline_idx]
    cutoff = round(0.1 * baseline)
    cutoff[cutoff <= 1] = 1
    # Return the index corresponding to the cutoff in descending order
    length(x) - findInterval(cutoff, x_ascending)
}

# Estimate the number of recovered cells by trying to find ordmag(recovered) =~ filtered.
# - Search for a result such that some loss(recovered_cells, filtered_cells) is minimized.
# - Here I'm using (obs - exp)**2 / exp, which approximates a proportion for small differences
#   but blows up for large differences.
# - Test over a log2-spaced range of values from 1..262_144
estimate_recovered_cells_ordmag <- function(nonzero_bc_counts, max_expected_cells) {
    recovered_cells = seq(1, log2(max_expected_cells), length.out = 2000)
    recovered_cells = unique(round(2^recovered_cells))

    baseline_bc_idx = round(recovered_cells * (1 - ORDMAG_RECOVERED_CELLS_QUANTILE))
    baseline_bc_idx[baseline_bc_idx > length(nonzero_bc_counts)] = length(nonzero_bc_counts) 
    
    filtered_cells = find_within_ordmag(nonzero_bc_counts, baseline_bc_idx)
    loss = (filtered_cells - recovered_cells)^2 / recovered_cells
    idx = which.min(loss)
    c(recovered_cells[idx], loss[idx])
}

# All barcodes that are close to within an order of magnitude of a top barcode.
# Takes all barcodes that are close to within an order of magnitude of a
# top barcode that likely represents a cell.
filter_cellular_barcodes_ordmag <- function(bc_counts) {
    nonzero_bc_counts = bc_counts[bc_counts > 0]
    if (length(nonzero_bc_counts) == 0) {
        stop("All barcodes have zero count; cannot proceed.")
    }
    
    set.seed(0)
    # Set the most cells to examine based on the empty drops range for this chemistry
    max_expected_cells = MAX_RECOVERED_CELLS_PER_GEM_GROUP
    recovered_cells_loss = 
        sapply(seq(1, ORDMAG_NUM_BOOTSTRAP_SAMPLES), function(x) {
        estimate_recovered_cells_ordmag(
            sample(nonzero_bc_counts, length(nonzero_bc_counts),replace = TRUE), max_expected_cells
        )
        })
    recovered_cells = mean(recovered_cells_loss[1,])
    loss = mean(recovered_cells_loss[2,])
    
    recovered_cells = max(round(recovered_cells), MIN_RECOVERED_CELLS_PER_GEM_GROUP)
    
    message(paste0("Found recovered_cells = ", recovered_cells, " with loss = ", loss))
    
    baseline_bc_idx = round(recovered_cells * (1 - ORDMAG_RECOVERED_CELLS_QUANTILE))
    baseline_bc_idx = min(baseline_bc_idx, length(nonzero_bc_counts))
    
    # Bootstrap sampling; run algo with many random samples of the data
    top_n_boot = sapply(seq(1, ORDMAG_NUM_BOOTSTRAP_SAMPLES), function(x) {
        find_within_ordmag(
        sample(nonzero_bc_counts, length(nonzero_bc_counts), replace = TRUE), baseline_bc_idx
        )
    })

    metrics = summarize_bootstrapped_top_n(top_n_boot, nonzero_bc_counts)
    
    # Get the filtered barcodes
    top_n = metrics.filtered_bcs

    if (top_n > len(nonzero_bc_counts)) {
        stop("Invalid selection of 0-count barcodes!")
    }
    return(metrics)
}

summarize_bootstrapped_top_n <- function(top_n_boot, nonzero_counts) {
    top_n_bcs_mean = mean(top_n_boot)
    # mimick cellranger, but not sure if we want to use sample var or population var
    top_n_bcs_var = sum((top_n_boot - top_n_bcs_mean)^2)/length(top_n_boot) 
    top_n_bcs_sd = sqrt(top_n_bcs_var)
    result = list(
        filtered_bcs = 0,
        filtered_bcs_lb = 0,
        filtered_bcs_ub = 0,
        filtered_bcs_var = 0,
        filtered_bcs_cv = 0,
        filtered_bcs_cutoff = 0
    )
    
    result$filtered_bcs_var = top_n_bcs_var
    result$filtered_bcs_cv = top_n_bcs_sd/top_n_bcs_mean
    result$filtered_bcs_lb = round(qnorm(0.025, mean = top_n_bcs_mean, sd = top_n_bcs_sd, lower.tail = TRUE), 0)
    result$filtered_bcs_ub = round(qnorm(0.975, mean = top_n_bcs_mean, sd = top_n_bcs_sd, lower.tail = TRUE), 0)
    
    nbcs = round(top_n_bcs_mean)
    result$filtered_bcs = nbcs
    
    # make sure that if a barcode with count x is selected, we select all barcodes with count >= x
    # this is true for each bootstrap sample, but is not true when we take the mean
    
    if (nbcs > 0) {
        sorted_counts = nonzero_counts[order(nonzero_counts, decreasing = TRUE)]
        
        cutoff = sorted_counts[nbcs]
        index = nbcs
        while (((index) < length(sorted_counts)) && sorted_counts[index] == cutoff) {
        index = index + 1
        }
        # if we end up grabbing too many barcodes, revert to initial estimate
        if ((index - nbcs) > 0.20 * nbcs) {
        return(result)
        }
        result$filtered_bcs = index
        result$filtered_bcs_cutoff = cutoff
    }
    return(result)
}
