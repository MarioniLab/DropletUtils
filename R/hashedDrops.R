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
#' @param ambient A numeric vector of length equal to \code{nrow(x)},
#' specifying the relative abundance of each HTO in the ambient solution.
#' See below for details.
#' @param pseudo.scale A numeric scalar specifying the scaling of the pseudo-count when computing log-fold changes.
#' Also serves as the minimum pseudo-count.
#' @param nmads A numeric scalar specifying the number of median absolute deviations (MADs) to use for calling outliers.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per column of \code{x}, containing the following fields:
#' \itemize{
#' \item \code{Total}, integer specifying the total count across all HTOs for each barcode.
#' \item \code{Best}, integer specifying the HTO with the highest abundance for each barcode.
#' \item \code{Second}, integer specifying the HTO with the second-highest abundance.
#' \item \code{LogFC.1to2}, numeric containing the log-fold change between the abundances of the best and second-best HTO.
#' \item \code{LogFC.2to3}, numeric containing the log-fold change between the abundances of the second-best and third-best HTO.
#' \item \code{Doublet}, logical specifying whether a barcode is a doublet.
#' \item \code{Confident}, logical specifying whether a barcode is a confidently assigned singlet.
#' }
#'
#' @details
#' The idea behind cell hashing is that cells from the same sample are stained with reagent conjugated with a single HTO.
#' Cells are mixed across multiple samples and subjected to droplet-based single-cell sequencing.
#' Cell barcode libraries can then be demultiplexed into individual samples based on whether their unique HTO is detected.
#'
#' We identify the sample of origin for each cell barcode as that corresponding to the most abundant HTO.
#' (See below for more details on exactly how \dQuote{most abundant} is defined.)
#' The log-fold change between the largest and second-largest ratios is also reported for each barcode, with large log-fold changes representing confident assignment to a single sample.
#' Libraries with low log-fold changes are then identified using an outlier detection strategy for potential removal.
#'
#' We also report the log-fold change between the second- and third-most abundant for each cell.
#' Large log-fold changes indicate that the second HTO is of much higher abundance than the third HTO, consistent with a doublet.
#' This reasoning assumes that the counts for the third and remaining HTOs are due to background contamination;
#' we ignore the (hopefully negligible) probability of higher-order multiplets.
#' Libraries with large log-fold changes are then identified as outliers for quality control.
#'
#' @section Outlier-based filtering:
#' We first identify putative doublets as those with \code{LogFC.2to3} values that are 3 MADs above the median.
#' This is done with a MAD computed from only the observations above the median,
#' to account for differences in the heaviness of the upper tail.
#' Doublet identities are reported in the \code{Doublet} field of the output. 
#' 
#' Of the non-doublet libraries, we consider them to be confidently assigned to a single sample if their \code{LogFC.1to2} values are \emph{not} less than 3 MADs below the median.
#' In this case, the MAD is computed using only observations below the median,
#' again to account for differences in the heaviness of the lower tail.
#' Confident assignments are marked in the \code{Confident} field of the output.
#' Note that any library for which \code{Confident} is {TRUE} will have \code{Doublet} set to \code{FALSE}, and vice versa.
#'
#' Besides the use of directionality,
#' another difference from \code{\link{mad}} is that we do not multiply by the normality constant.
#' This is primarily motivated by the admission that the distributions are not normal anyway.
#' In practice, this means that \code{nmad} is effectively lower than what it would otherwise be,
#' but that hardly matters all that much for an arbitrary threshold.
#'
#' @section Computing abundances:
#' HTO abundances require some care to compute due to the presence of ambient contamination in each library.
#' Ideally, the experiment would be performed in such a manner that the concentration of each HTO is the same.
#' However, if one HTO is present at higher concentration in the ambient solution,
#' this might incorrectly cause us to assign all barcodes to the corresponding sample.
#' 
#' To adjust for ambient contamination, we make the following series of assumptions:
#' \enumerate{
#' \item A minority of HTOs in a library are actually due to the presence of cell(s), the rest coming from the ambient solution.
#' \item The number of transcipt molecules (not necessarily UMIs) from the ambient solution is similar in each droplet.
#' \item Ambient contamination in each library follows the same profile as \code{ambient}.
#' }
#' 
#' We estimate the level of ambient contamination in each barcode using a \pkg{DESeq}-like normalization algorithm to compute the scaling factor to be applied to \code{ambient}.
#' (The requisite assumption of a non-DE majority is motivated by assumptions 1 and 3 above.)
#' Some care is taken when the number of HTOs is below 4 to ensure that scaling is not distorted by doublets.
#' We then subtract the scaled ambient proportions from the count profile of each cell,
#' yielding adjusted abundances that are (nominally) free from the effects of contamination.
#' Abundances that would otherwise be negative are set to zero.
#' 
#' Prior to computing the log-fold changes, we add a library-specific pseudo-count to all abundances.
#' This is defined for each library by taking the average ambient HTO count (i.e., the average of the scaled proportions) and scaling it by \code{pseudo.scale}.
#' In effect, the ambient contamination provides a natural pseudo-count that scales with capture efficiency and sequencing depth;
#' indeed, under assumption 2, any differences in ambient coverage between libraries can be solely attributed to such technical factors.
#'
#' This scaling property is useful as it ensures that shrinkage of the log-fold changes is not more severe for libraries that have not been sequenced as deeply.
#' We thus avoid excessive variability in the log-fold change distribution (and reduction in precision of outlier detection).
#' Another nice aspect of this adjustment is that it collapses to a no-op if the experiment is well-executed with identical concentrations of all HTOs in the ambient solution.
#' 
#' @section Getting the ambient proportions:
#' Ideally, \code{ambient} would be obtained from libraries that do not correspond to cell-containing droplets.
#' For example, we could get this information from the \code{\link{metadata}} of the \code{\link{emptyDrops}} output,
#' had we run \code{\link{emptyDrops}} on the HTO count matrix (see below).
#' 
#' In some cases (e.g., public data), counts are provided for only the cell-containing barcodes.
#' To handle this, we compute the median of each HTO across all barcodes to obtain a rough proxy for the ambient profile.
#' However, this does assume that there are at least 3 HTOs with evenly distributed numbers of cells in each sample.
#'
#' @section Use only on non-empty droplets:
#' This function assumes that cell calling has already been performed, e.g., with \code{\link{emptyDrops}}.
#' Specifically, \code{x} should only contain columns corresponding to non-empty droplets.
#' If empty droplets are included, their log-fold changes will simply reflect stochastic sampling in the ambient solution
#' and violate the assumptions involved in outlier detection.
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
#' hist(stats$NMAD.2to3, breaks=50)
#'
#' # Identify doublets at a given NMAD threshold.
#' is.doublet <- stats$NMAD >= 5
#' summary(is.doublet)
#' 
#' # Chcecking against the known truth, in this case
#' # 'Best' contains the putative sample of origin.
#' table(stats$Best, true.sample) 
#'
#' @references
#' Stoeckius M, Zheng S, Houck-Loomis B et al. (2018)
#' Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics.
#' \emph{Genome Biol.} 19, 1:224
#'
#' @seealso
#' \code{\link{emptyDrops}}, to identify which barcodes are likely to contain cells.
#' 
#' @export
#' @importFrom Matrix t colSums
#' @importFrom S4Vectors DataFrame
#' @importFrom stats median
hashedDrops <- function(x, ambient=NULL, pseudo.scale=1, nmads=3) {
    totals <- colSums(x)
    cell.names <- colnames(x)

    if (is.null(ambient)) {
        ambient <- vapply(seq_len(nrow(x)), function(i) median(x[i,]), 0)
    }

    discard <- ambient == 0
    x <- x[!discard,,drop=FALSE]
    ambient <- ambient[!discard]

    output <- hashed_deltas(x, ambient, pseudo.scale)
    lfc.1to2 <- log2(output$FC)
    lfc.2to3 <- log2(output$FC2)

    # Using *directional* MADs to avoid 
    med.2to3 <- median(lfc.2to3)
    delta.2to3 <- lfc.2to3[lfc.2to3 > med.2to3] - med.2to3
    upper.threshold <- med.2to3 + nmads * median(delta.2to3) 
    is.doublet <- lfc.2to3 > upper.threshold 

    lfc.singlet <- lfc.1to2[!is.doublet]
    med.singlet <- median(lfc.singlet)
    delta.singlet <- med.singlet - lfc.singlet[lfc.singlet < med.singlet]
    lower.threshold <- med.singlet - nmads * median(delta.singlet) 
    confident.singlet <- lfc.1to2 > lower.threshold & !is.doublet

    DataFrame(
        row.names=cell.names,
        Total=totals,
        Best=output$Best+1L,
        Second=output$Second+1L,
        LogFC.1to2=lfc.1to2,
        LogFC.2to3=lfc.2to3,
        Doublet=is.doublet,
        Confident=confident.singlet
    )
}
