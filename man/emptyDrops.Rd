\name{emptyDrops}
\alias{emptyDrops}
\alias{testEmptyDrops}

\title{Identify empty droplets}
\description{Distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment.}

\usage{
testEmptyDrops(m, lower=100, niters=10000, test.ambient=FALSE,
    ignore=NULL, alpha=NULL, BPPARAM=SerialParam())

emptyDrops(m, lower=100, retain=NULL, barcode.args=list(), ...)
}

\arguments{
\item{m}{A numeric matrix object - usually a \linkS4class{dgTMatrix} or \linkS4class{dgCMatrix} - containing droplet data \emph{prior to any filtering or cell calling}.
Columns represent barcoded droplets, rows represent genes.
}
\item{lower}{A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.}
\item{niters}{An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.}
\item{test.ambient}{A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to \code{lower}.}
\item{ignore}{A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored (see Details for how this differs from \code{lower}).}
\item{alpha}{A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme.}
\item{BPPARAM}{A BiocParallelParam object indicating whether parallelization should be used to compute p-values.}
\item{retain}{A numeric scalar specifying the threshold for the total UMI count above which all barcodes are assumed to contain cells.}
\item{barcode.args}{Further arguments to pass to \code{\link{barcodeRanks}}.}
\item{...}{Further arguments to pass to \code{testEmptyDrops}.}
}

\section{Details about \code{testEmptyDrops}}{
The \code{testEmptyDrops} function will obtain an estimate of the composition of the ambient pool of RNA based on the barcodes with total UMI counts less than or equal to \code{lower}.
This assumes that a cell-containing droplet would generally have higher total counts than empty droplets containing RNA from the ambient pool.
Counts for the low-count barcodes are pooled together, and an estimate of the proportion vector for the ambient pool is calculated using \code{\link{goodTuringProportions}}.
The count vector for each barcode above \code{lower} is then tested for a significant deviation from these proportions.

The null hypothesis is that transcript molecules are included into droplets by multinomial sampling from the ambient profile.
For each barcode, the probability of obtaining its count vector based on the null model is computed.
Then, \code{niters} count vectors are simulated from the null model.
The proportion of simulated vectors with probabilities lower than the observed multinomial probability for that barcode is used to calculate the p-value.
We use this Monte Carlo approach as an exact multinomial p-value is difficult to calculate.

The \code{ignore} argument can also be set to ignore barcodes with total counts less than or equal to \code{ignore}.
This differs from the \code{lower} argument in that the ignored barcodes are not necessarily used to compute the ambient profile.
Users can interpret \code{ignore} as the minimum total count required for a barcode to be considered as a potential cell.
In contrast, \code{lower} is the maximum total count below which all barcodes are assumed to be empty droplets.
}

\section{Details about \code{emptyDrops}}{
The \code{emptyDrops} function combines the results of \code{testEmptyDrops} with \code{\link{barcodeRanks}} to identify droplets that are likely to contain cells.
Barcodes that contain more than \code{retain} total counts are always retained.
This ensures that large cells with profiles that are very similar to the ambient pool are not inadvertently discarded.
If \code{retain} is not specified, it is set to the total count at the knee point detected by \code{\link{barcodeRanks}}.
Manual specification of \code{retain} may be useful if the knee point was not correctly identified in complex log-rank curves.
Users can also set \code{retain=Inf} to disable automatic retention of barcodes with large totals.

The Benjamini-Hochberg correction is also applied to the Monte Carlo p-values to correct for multiple testing.
Cells can then be defined by taking all barcodes with significantly non-ambient profiles, e.g., at a false discovery rate of 0.1\%.
All barcodes with total counts above \code{K} (or \code{retain}) are assigned p-values of zero \emph{during correction}, reflecting our assumption that they are true positives.
This ensures that their Monte Carlo p-values do not affect the correction of other genes, and also means that they will have FDR values of zero.
Nonetheless, their original Monte Carlo p-values are still reported in the output.

This function will also attempt to detect whether \code{m} contains non-integer values by seeing if the column and row sums are discrete.
If such values are present, \code{m} is first \code{\link{round}}ed to the nearest integer value before proceeding.
This may be relevant when the count matrix is generated from pseudo-alignment methods like Alevin (see the \pkg{tximeta} package for details).
Rounding is performed as discrete count values are necessary for the Good-Turing algorithm and for the multinomial simulations.
}

\section{Handling overdispersion}{
If \code{alpha} is set to a positive number, sampling is assumed to follow a Dirichlet-multinomial (DM) distribution.
The parameter vector of the DM distribution is defined as the estimated ambient profile scaled by \code{alpha}.
Smaller values of \code{alpha} model overdispersion in the counts, due to dependencies in sampling between molecules.
If \code{alpha=NULL}, a maximum likelihood estimate is obtained from the count profiles for all barcodes with totals less than or equal to \code{lower}.
If \code{alpha=Inf}, the sampling of molecules is modelled with a multinomial distribution.

Users can check whether the model is suitable by extracting the p-values for all barcodes with \code{test.ambient=TRUE}.
Under the null hypothesis, the p-values for presumed ambient barcodes (i.e., with total counts below \code{lower}) should be uniformly distributed.
Skews in the p-value distribution are indicative of an inaccuracy in the model and/or its estimates (of \code{alpha} or the ambient profile).
}

\value{
\code{testEmptyDrops} will return a DataFrame with the following components:
\describe{
\item{\code{Total}:}{Integer, the total UMI count for each barcode.}
\item{\code{LogProb}:}{Numeric, the log-probability of observing the barcode's count vector under the null model.}
\item{\code{PValue}:}{Numeric, the Monte Carlo p-value against the null model.}
\item{\code{Limited}:}{Logical, indicating whether a lower p-value could be obtained by increasing \code{niters}.}
}
For barcodes with counts below \code{lower}, \code{NA} values are returned for all fields if \code{test.ambient=FALSE}.
This is to ensure that the number of rows in the output DataFrame is identical to \code{ncol(m)}.

\code{emptyDrops} will return a DataFrame like \code{testEmptyDrops}, with an additional \code{FDR} field.

The metadata of the output DataFrame will contains the ambient profile in \code{ambient}, the estimated/specified value of \code{alpha}, the specified value of \code{lower} and the number of iterations in \code{niters}.
For \code{emptyDrops}, the metadata will also contain the retention threshold in \code{retain}.
}

\author{
Aaron Lun
}

\examples{
# Mocking up some data:
set.seed(0)
my.counts <- DropletUtils:::simCounts()

# Identify likely cell-containing droplets.
out <- emptyDrops(my.counts)
out

is.cell <- out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)

# Check if p-values are lower-bounded by 'niters'
# (increase 'niters' if any Limited==TRUE and Sig==FALSE)
table(Sig=is.cell, Limited=out$Limited)
}

\references{
Lun A, Riesenfeld S, Andrews T, Dao TP, Gomes T, participants in the 1st Human Cell Atlas Jamboree, Marioni JC (2018).
Distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data.
\emph{biorXiv}.

Phipson B, Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Stat. Appl. Genet. Mol. Biol.} 9:Article 39.
}

\seealso{
\code{\link{barcodeRanks}},
\code{\link{defaultDrops}}
}
