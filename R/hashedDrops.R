#' Demultiplex cell hashing data
#'
#' Demultiplex cell barcodes into their samples of origin based on the most abundant hash tag oligo (HTO).
#' Also identify potential doublets based on the presence of multiple significant HTOs.
#'
#' Note that this function is still experimental; feedback is welcome.
#'
#' @param x A numeric/integer matrix-like object containing UMI counts.
#' Rows correspond to HTOs and columns correspond to cell barcodes.
#' Each barcode is assumed to correspond to a cell, i.e., cell calling is assumed to have already been performed.
#' @param ambient A numeric vector of length equal to \code{nrow(x)},
#' specifying the relative abundance of each HTO in the ambient solution - see Details.
#' @param min.prop Numeric scalar to be used to infer the ambient profile when \code{ambient=NULL}, see \code{\link{inferAmbience}}.
#' @param pseudo.count A numeric scalar specifying the minimum pseudo-count when computing log-fold changes.
#' @param doublet.nmads A numeric scalar specifying the number of median absolute deviations (MADs) to use to identify doublets.
#' @param doublet.min A numeric scalar specifying the minimum threshold on the log-fold change to use to identify doublets.
#' @param doublet.mixture Logical scalar indicating whether to use a 2-component mixture model to identify doublets.
#' @param confident.nmads A numeric scalar specifying the number of MADs to use to identify confidently assigned singlets.
#' @param confident.min A numeric scalar specifying the minimum threshold on the log-fold change to use to identify singlets.
#' @param combinations An integer matrix specifying valid \emph{combinations} of HTOs.
#' Each row corresponds to a single sample and specifies the indices of rows in \code{x} corresponding to the HTOs used to label that sample.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per column of \code{x}, containing the following fields:
#' \itemize{
#' \item \code{Total}, integer specifying the total count across all HTOs for each barcode.
#' \item \code{Best}, integer specifying the HTO with the highest abundance for each barcode.
#' \item \code{Second}, integer specifying the HTO with the second-highest abundance.
#' \item \code{LogFC}, numeric containing the log-fold change between the abundances of the best and second-best HTO.
#' \item \code{LogFC2}, numeric containing the log-fold change in the second-best HTO over the ambient contamination.
#' \item \code{Doublet}, logical specifying whether a barcode is a doublet.
#' \item \code{Confident}, logical specifying whether a barcode is a confidently assigned singlet.
#' }
#' In addition, the metadata contains \code{ambient}, a numeric vector containing the (estimate of the) ambient profile;
#' \code{doublet.threshold}, the threshold applied to \code{LogFC2} to identify doublets;
#' and \code{confident.threshold}, the threshold applied to non-doublet \code{LogFC} values to identify confident singlets.
#'
#' If \code{combinations} is specified, \code{Best} instead specifies the sample (i.e., row index of \code{combinations}).
#' The interpretation of \code{LogFC} and \code{LogFC2} are slightly different, and \code{Second} is not reported - see details.
#'
#' @details
#' The idea behind cell hashing is that cells from the same sample are stained with reagent conjugated with a single HTO.
#' Cells are mixed across multiple samples and subjected to droplet-based single-cell sequencing.
#' Cell barcode libraries can then be demultiplexed into individual samples based on whether their unique HTO is detected.
#'
#' We identify the sample of origin for each cell barcode as that corresponding to the most abundant HTO.
#' (See below for more details on exactly how \dQuote{most abundant} is defined.)
#' The log-fold change between the largest and second-largest abundances is also reported for each barcode, with large log-fold changes representing confident assignment to a single sample.
#' We also report the log-fold change of the second-most abundant HTO over the estimated level of ambient contamination.
#' Large log-fold changes indicate that the second HTO has greater abundance than expected, consistent with a doublet.
#'
#' To facilitate quality control, we explicitly identify problematic barcodes as outliers on the relevant metrics.
#' \itemize{
#' \item By default, we identify putative doublets as those with \code{LogFC2} values that are
#' (i) \code{doublet.nmads} MADs above the median and (ii) greater than \code{doublet.min}.
#' The hard threshold is more-or-less arbitrary and aims to avoid overly aggressive detection 
#' of large outliers in a naturally right-skewed distribution 
#' (given that the log-fold changes are positive by definition, and most of the distribution is located near zero).
#' \item Alternatively, if \code{doublet.mixture=TRUE}, we fit a two-component mixture model to the \code{LogFC2} distribution.
#' Doublets are identified as all members of the component with the larger mean.
#' This avoids the need for the arbitrary parameters mentioned above but only works when there are many doublets,
#' otherwise both components will be fitted to the non-doublet values.
#' (Initialization of the model assumes at least 5\% doublets.)
#' }
#' Of the non-doublet libraries, we consider them to be confidently assigned to a single sample if their \code{LogFC} values are (i) \emph{not} less than \code{confident.nmads} MADs below the median and (ii) greater than \code{confident.min}.
#' The hard threshold is again arbitrary, but this time it aims to avoid insufficiently aggressive outlier detection - 
#' typically from an inflation of the MAD when the \code{LogFC} values are large, positive and highly variable.
#'
#' @section Adjusting abundances for ambient contamination:
#' HTO abundances require some care to compute due to the presence of ambient contamination in each library.
#' Ideally, the experiment would be performed in such a manner that the concentration of each HTO is the same.
#' However, if one HTO is present at higher concentration in the ambient solution,
#' this might incorrectly cause us to assign all barcodes to the corresponding sample.
#' 
#' To adjust for ambient contamination, we assume that the ambient contamination in each library follows the same profile as \code{ambient}.
#' We further assume that a minority of HTOs in a library are actually driven by the presence of cell(s), the rest coming from the ambient solution.
#' We estimate the level of ambient contamination in each barcode by scaling \code{ambient}, using a \pkg{DESeq}-like normalization algorithm to compute the scaling factor.
#' (The requisite assumption of a non-DE majority follows from the two assumptions above.)
#' We then subtract the scaled ambient proportions from the HTO count profile to remove the effect of contamination.
#' Abundances that would otherwise be negative are set to zero.
#'
#' The scaling factor for each cell is defined by computing ratios between the HTO counts and \code{ambient}, and taking the median across all HTOs.
#' However, this strict definition is only used when there are at least 5 HTOs being considered.
#' For experiments with 3-4 HTOs, we assume that higher-order multiplets are negligible and define the scaling factor as the third-largest ratio.
#' For experiments with only 2 HTOs, the second-most abundant HTO is always used to estimate the ambient contamination.
#' 
#' Ideally, \code{ambient} would be obtained from libraries that do not correspond to cell-containing droplets.
#' For example, we could get this information from the \code{\link{metadata}} of the \code{\link{emptyDrops}} output,
#' had we run \code{\link{emptyDrops}} on the HTO count matrix (see below).
#' Unfortunately, in some cases (e.g., public data), counts are provided for only the cell-containing barcodes.
#' If \code{ambient=NULL}, the profile is inferred from \code{x} using \code{\link{inferAmbience}}.
#'
#' @section Computing the log-fold changes:
#' After subtraction of the ambient noise but before calculation of the log-fold changes,
#' we need to add a pseudo-count to ensure that the log-fold changes are well-defined.
#' We set the pseudo-count to the average ambient HTO count (i.e., the average of the scaled \code{ambient}), effectively exploiting the ambient contamination as a natural pseudo-count that scales with barcode-specific capture efficiency and sequencing depth.
#' (In libraries with low sequencing depth, we still enforce a minimum pseudo-count of \code{pseudo.count}.)
#'
#' This scaling behavior is useful as it ensures that shrinkage of the log-fold changes is not more severe for libraries that have not been sequenced as deeply.
#' We thus avoid excessive variability in the log-fold change distribution that would otherwise reduce the precision of outlier detection.
#' The implicit assumption is that the number of contaminating transcript molecules is roughly the same in each droplet, meaning that any differences in ambient coverage between libraries reflect technical biases that would also affect cell-derived HTO counts. 
#'
#' Another nice aspect of this entire procedure (subtraction and re-addition) is that it collapses to a no-op if the experiment is well-executed with identical concentrations of all HTOs in the ambient solution.
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
#' @section Handling 2 or fewer samples:
#' If \code{x} has no more than two rows, \code{Doublet}, \code{LogFC2} and \code{doublet.threshold} are set to \code{NA}.
#' Strictly speaking, doublet detection is not possible as the second HTO is always used to estimate the ambient scaling and thus \code{LogFC2} is always zero.
#' However, \code{Confident} calls are still available in the output of this function so assignment to the individual samples can still be performed.
#' In this scenario, the non-confident assignments are probably also doublets, though this cannot be said with much certainty.
#' 
#' If \code{x} has no more than one row, \code{Confident}, \code{LogFC} and \code{confident.threshold} are set to \code{NA}.
#' Obviously, if there is only one HTO, the identity of the assigned sample is a foregone conclusion.
#'
#' @section Resolving combinatorial hashes:
#' In some applications, samples are labelled with a combination of HTOs to enable achieve greater multiplexing throughput.
#' This is accommodated by passing \code{combinations} to specify the valid HTO combinations that were used for sample labelling.
#' Each row of \code{combinations} corresponds to a sample and should contain non-duplicated row indices of \code{x} corresponding to the HTOs used in that sample.
#'
#' The calculation for the single-HTO case is then generalized for HTO combinations.
#' The most important differences are that:
#' \itemize{
#' \item The reported \code{LogFC} is now the log-fold change between the \eqn{n}th most abundant HTO and the \eqn{n+1}th HTO,
#' where \eqn{n} is the number of HTOs in a valid combination. 
#' This captures the drop-off in abundance beyond the expected number of HTOs.
#' \item The reported \code{LogFC2} is now the log-fold change of the \eqn{n+1}th HTO over the ambient solution.
#' This captures the high abundance of the more-than-expected number of HTOs when doublets are present.
#' \item \code{Best} no longer refers to the row index of \code{x}, but instead to the row index of \code{combinations}.
#' This may contain \code{NA} values if a particular combination of HTOs is observed but not present in the expected set.
#' \item \code{Second} is no longer reported as we cannot conveniently determine the identity of the second sample.
#' }
#'
#' We also generalize the edge-case behavior when there are not enough HTOs to support doublet detection. 
#' Consider that an inter-sample doublet may result in either \eqn{n + 1} to \eqn{2n} abundant HTOs.
#' Estimation of the scaling factor will attempt to avoid using the top \eqn{2n} ratios.
#' If \code{nrow(x)} is equal to or less than \eqn{n + 1}, doublet statistics will not be reported at all, 
#' i.e., \code{Doublet} and \code{LogFC2} are set to \code{NA}.
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
#' # Doublets show up in the top-left, singlets in the bottom right.
#' plot(stats$LogFC, stats$LogFC2)
#'
#' # Most cells should be singlets with low second log-fold changes.
#' hist(stats$LogFC2, breaks=50)
#'
#' # Identify confident singlets or doublets at the given threshold.
#' summary(stats$Confident)
#' summary(stats$Doublet)
#' 
#' # Checking against the known truth, in this case
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
#' @importFrom stats median mad
#' @importFrom beachmat colBlockApply
hashedDrops <- function(x, ambient=NULL, min.prop=0.05, pseudo.count=5, 
    doublet.nmads=3, doublet.min=2, doublet.mixture=FALSE, confident.nmads=3, confident.min=2, combinations=NULL)
{ 
    totals <- colSums(x)
    cell.names <- colnames(x)

    if (is.null(ambient)) {
        ambient <- inferAmbience(x, min.prop)
    }

    discard <- ambient == 0
    x <- x[!discard,,drop=FALSE]
    ambient <- ambient[!discard]
    
    if (is.null(combinations)) {
        n.expected <- 1L
    } else {
        n.expected <- ncol(combinations)
    }

    output <- colBlockApply(x, FUN=hashed_deltas, prop=ambient, pseudo_count=pseudo.count, n_expected=n.expected)

    lfc <- log2(unlist(lapply(output, "[[", i="FC"), use.names=FALSE))
    lfc2 <- log2(unlist(lapply(output, "[[", i="FC2"), use.names=FALSE))
    best.sample <- do.call(cbind, lapply(output, "[[", i="Best"))
    second.sample <- unlist(lapply(output, "[[", i="Second"), use.names=FALSE)

    no.lfc2 <- all(is.na(lfc2))
    if (!no.lfc2) {
        if (!doublet.mixture) {
            # Using outlier detection to prune out doublets.
            med2 <- median(lfc2)
            mad2 <- mad(lfc2, center=med2)
            upper.threshold <- med2 + doublet.nmads * mad2
            upper.threshold <- max(upper.threshold, doublet.min)
            is.doublet <- lfc2 > upper.threshold 
        } else {
            is.doublet <- !.get_lower_dist(lfc2, p=0.05)
            upper.threshold <- max(lfc2[!is.doublet])
        }
    } else {
        # To keep the next step happy; we will flip it back to all-NA's soon enough.
        is.doublet <- logical(ncol(x))
        upper.threshold <- NA_real_
    }

    # Using outlier detection to identify confident singlets.
    if (!all(is.na(lfc))) {
        lfc.singlet <- lfc[!is.doublet]
        med.singlet <- median(lfc.singlet)
        mad.singlet <- mad(lfc.singlet, center=med.singlet)

        lower.threshold <- med.singlet - confident.nmads * mad.singlet
        lower.threshold <- max(lower.threshold, confident.min)
        confident.singlet <- lfc > lower.threshold & !is.doublet
    } else {
        confident.singlet <- rep(NA, ncol(x))
        lower.threshold <- NA_real_
    }

    if (no.lfc2) {
        is.doublet <- rep(NA, ncol(x))
    }

    # Translating from indices to combinations, if requested.
    best.sample <- best.sample + 1L
    second.sample <- second.sample + 1L

    if (is.null(combinations)) {
        best.sample <- drop(best.sample)
    } else {
        best.sample <- t(best.sample)
        combinations <- t(apply(combinations, 1, sort))
        colnames(combinations) <- colnames(best.sample) <- seq_len(n.expected)
        best.sample <- match(DataFrame(best.sample), DataFrame(combinations))
    }

    output <- DataFrame(
        row.names=cell.names,
        Total=totals,
        Best=best.sample,
        Second=second.sample,
        LogFC=lfc,
        LogFC2=lfc2,
        Doublet=is.doublet,
        Confident=confident.singlet
    )
   
    if (!is.null(combinations)) {
        output$Second <- NULL
    }

    metadata(output) <- list(
        ambient=ambient,
        confident.threshold=lower.threshold,
        doublet.threshold=upper.threshold
    )

    output
}
