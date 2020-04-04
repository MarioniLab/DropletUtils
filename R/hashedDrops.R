#' Demultiplex cell hashing data
#'
#' Demultiplex cell barcodes into their samples of origin based on the most significant hash tag oligo (HTO).
#' Also identify potential doublets based on the presence of multiple significant HTOs.
#'
#' Note that this function is still experimental; feedback is welcome.
#'
#' @param x A numeric/integer matrix-like object containing UMI counts.
#' Rows correspond to HTOs and columns correspond to cell barcodes.
#' Each barcode is assumed to correspond to a cell, i.e., cell calling is assumed to have already been performed.
#' @param prop A numeric vector of length equal to \code{nrow(x)},
#' specifying the relative abundance of each HTO across all samples.
#' @param pseudo.count A numeric scalar specifying the pseudo-count to add when computing log-fold changes.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per column of \code{x}, containing the following fields:
#' \itemize{
#' \item \code{Total}, integer specifying the total count across all HTOs for each barcode.
#' \item \code{Best}, integer specifying the HTO with the highest deviation from the expected proportions for each barcode.
#' \item \code{Second}, integer specifying the HTO with the second-highest deviation from the expected deviation.
#' \item \code{LogFC.1to2}, numeric containing the log-fold change between the counts of the best and second-best HTO.
#' \item \code{NMAD.1to2}, numeric specifying the number of MADs from the median for each value of \code{LogFC.1to2}.
#' \item \code{LogFC.2to3}, numeric containing the log-fold change between the counts of the second-best and third-best HTO.
#' \item \code{NMAD.2to3}, numeric specifying the number of MADs from the median for each value of \code{LogFC.2to3}.
#' }
#'
#' @details
#' The idea behind cell hashing is that cells from the same sample are stained with reagent conjugated with a single HTO.
#' Cells are mixed across multiple samples and subjected to droplet-based single-cell sequencing.
#' Cell barcode libraries can then be demultiplexed into individual samples based on whether their unique HTO is detected.
#'
#' We identify the sample of origin for each cell barcode as that corresponding to the HTO with the greatest deviation from its expected proportion in the multiplexed pool.
#' This deviation is quantified by the ratio of the count for each HTO to the corresponding value of \code{prop}.
#' The log-fold change between the largest and second-largest ratios is also reported for each barcode, with large log-fold changes representing confident assignment to a single sample.
#'
#' We identify potential doublets by looking at the log-fold change in the second- and third-largest ratios.
#' Large log-fold changes indicate that the second HTO is of higher abundance than the third HTO, consistent with a doublet. 
#' Note that we assume that the counts for the third and remaining HTOs are due to background contamination;
#' this ignores the (hopefully negligible) probability of higher-order multiplets.
#'
#' For both log-fold changes, we compute the number of median absolute deviations (MADs) from the median across all cells.
#' This allows users to easily identify doublets as large outliers, e.g., with \code{NMAD.2to3} above some threshold like 3 or 5.
#' More sophisticated approaches may also consider \emph{small} outliers with negative values in \code{NMAD.1to2},
#' where doublets will reduce the confidence of assignment to any single sample.
#'
#' We quantify deviation against \code{prop} to account for imbalances in the HTO distribution in the multiplexed sample.
#' Ideally, \code{prop} would be set to the HTO profile in the ambient solution in order to capture deviations due to the presence of labelled cells.
#' However, such information may not be available in many datasets where non-cell droplets are filtered out.
#' As such, the default \code{prop} is instead defined from the row sums of \code{x} under the assumption that the rate of ambient to non-ambient transcripts is the same for each sample.
#'
#' @section Use only on non-empty droplets:
#' This function assumes that cell calling has already been performed, e.g., with \code{\link{emptyDrops}}.
#' Specifically, \code{x} should only contain columns corresponding to non-empty droplets.
#' If empty droplets are included, their log-fold changes will simply reflect stochastic sampling in the ambient solution
#' and violate the assumptions involved in interpreting the \code{NMAD} fields.
#'
#' If \code{x} contains columns for both empty and non-empty droplets,
#' it is straightforward to simply run \code{\link{emptyDrops}} on the HTO count matrix to identify the latter.
#' Note that some fiddling with the \code{lower=} argument may be required,
#' depending on the sequencing depth of the HTO libraries.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up an example dataset with 10 HTOs and 10% doublets.
#' ncells <- 1000
#' nhto <- 10
#' y <- matrix(rpois(ncells*nhto, 50), nrow=nhto)
#' true.sample <- sample(nhto, ncells, replace=TRUE)
#' y[cbind(true.sample, seq_len(ncells))] <- 1000
#'
#' ndoub <- ncells/10
#' next.sample <- (true.sample[1:ndoub]  + 1) %% nrow(y)
#' next.sample[next.sample==0] <- nrow(y)
#' y[cbind(next.sample, seq_len(ndoub))] <- 500
#'
#' # Computing hashing statistics.
#' stats <- hashedDrops(y)
#'
#' # Doublets show up in the top-left,
#' # singlets in the bottom right.
#' plot(stats$LogFC.1to2, stats$LogFC.2to3)
#'
#' # Most cells should be singlets with low NMAD.
#' hist(stats$NMAD, breaks=50)
#'
#' # Identify doublets at a given NMAD threshold.
#' is.doublet <- stats$NMAD >= 5
#' summary(is.doublet)
#' 
#' # 'Best' contains the putative sample of origin.
#' table(stats$Best, true.sample) 
#'
#' @references
#' Stoeckius M, Zheng S, Houck-Loomis B et al. (2018)
#' Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics.
#' \emph{Genome Biol.} 19, 1:224
#'
#' 
#' @export
#' @importFrom Matrix t colSums
#' @importFrom S4Vectors DataFrame
hashedDrops <- function(x, prop=NULL, pseudo.count=3) {
    totals <- colSums(x)
    cell.names <- colnames(x)

    if (is.null(prop)) {
        prop <- rowMeans(x)
    }

    prop <- prop/mean(prop)
    discard <- prop == 0
    x <- x[!discard,,drop=FALSE]
    prop <- prop[!discard]

    output <- hashed_deltas(x, prop, pseudo.count)
    DataFrame(
        row.names=cell.names,
        Total=totals,
        Best=output$Best+1L,
        Second=output$Second+1L,
        LogFC.1to2=log2(output$FC),
        NMAD.1to2=.nmads(output$FC),
        LogFC.2to3=log2(output$FC2),
        NMAD.2to3=.nmads(output$FC2)
    )
}

#' @importFrom stats median mad
.nmads <- function(fc) {
    med <- median(fc)
    mad <- mad(fc)
    (fc - med)/mad
}
