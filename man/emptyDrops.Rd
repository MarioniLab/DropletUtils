\name{emptyDrops}
\alias{emptyDrops}
\alias{testEmptyDrops}

\title{Identify empty droplets}
\description{Distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment.}

\usage{
testEmptyDrops(m, lower=100, span=sqrt(2), npts=10000, 
    test.ambient=FALSE, BPPARAM=SerialParam()) 

emptyDrops(m, lower=100, scale=1, test.args=list(),
    barcode.args=list()) 
}

\arguments{
\item{m}{A real sparse matrix object, either a dgTMatrix or dgCMatrix.
Columns represent barcoded droplets, rows represent cells.
}
\item{lower}{A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.}
\item{span, npts}{Numeric scalars controlling the accuracy of the p-value calculation.}
\item{test.ambient}{A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to \code{lower}.}
\item{BPPARAM}{A BiocParallelParam object indicating whether parallelization should be used to compute p-values.}
\item{scale}{A numeric scalar specifying the fold-change about the knee-point for the total UMI count, above which all barcodes are assumed to contain cells.}
\item{test.args}{Further arguments to pass to \code{testEmptyDrops}.}
\item{barcode.args}{Further arguments to pass to \code{\link{barcodeRanks}}.}
}

\details{
The \code{testEmptyDrops} function will obtain an estimate of the composition of the ambient pool of RNA based on the droplets with total UMI counts below \code{lower}.
This assumes that a cell-containing droplet would generally be able to achieve some minimum count.
The low-count droplets are pooled together, and an estimate of the proportion vector for the ambient pool is calculated using \code{\link[edgeR]{goodTuringProportions}}. 
The count vector for each droplet above \code{lower} is then tested for a significant deviation from these proportions,
under a Poisson sampling scheme for inclusion of transcript molecules into a droplet.

The distribution of the deviances under the null Poisson model follows a chi-squared distribution based on the saddlepoint approximation.
In practice, this approximation fails at low counts and a parametric simulation is instead required to obtain a null distribution.
This is done by simulating Poisson counts for droplets with the same composition as the ambient pool but at a range of total UMI counts.
The p-value for each real droplet is calculated by taking all simulated droplets within \code{span}-fold of the total UMI count,
and calculating the proportion of those simulated droplets that have deviances greater than the observed value.

For this calculation, a \code{span} value is used as the variance of the deviance is dependent on the total count.
Reducing the \code{span} will yield a more accurate p-value for each droplet.
The number of simulated droplets is chosen such that there will always be \code{npts} values within \code{span} at any total UMI count above \code{lower}.
Increasing \code{npts} will improve the precision and reproducibility of the p-values, at the cost of increasing the run-time.
It will also reduce the lower bound on the permutation p-values (see Phipson and Smyth, 2010).

The \code{emptyDrops} function combines the test results with distributional information from \code{\link{barcodeRanks}} to identify droplets that are likely to be cells.
Specifically, the total count at the knee point \code{K} is determined, and droplets that contain more than \code{scale*K} total counts are always retained.
This ensures that cells with profiles that are very similar to the ambient pool are not inadvertently discarded.
The Benjamini-Hochberg correction is also applied to the p-values from \code{testEmptyDrops} to correct for multiple testing.
Cells can then be defined by taking all droplets with significantly non-ambient profiles, e.g., at a false discovery rate of 1\%.
}

\value{
\code{testEmptyDrops} will return a DataFrame with the following components:
\describe{
\item{\code{Total}:}{Integer, the total UMI count for each droplet.}

\item{\code{Deviance}:}{Numeric, the observed deviance of the droplet count vector from the ambient profile.}
\item{\code{Expected}:}{Numeric, the expected deviance under a null model of Poisson sampling.}
\item{\code{PValue}:}{Numeric, the p-value computed using the simulated droplets.}
\item{\code{Limited}:}{Logical, indicating whether a lower p-value could be obtained by increasing \code{npts}.}
}
For droplets with counts below \code{lower}, \code{NA} values are returned for all fields if \code{test.ambient=FALSE}.
This is to ensure that the number of rows in the output DataFrame is identical to \code{ncol(m)}. 

\code{emptyDrops} will return a DataFrame like \code{testEmptyDrops}, with an additional \code{FDR} field.
Any droplets with counts larger than the product of \code{scale} and the knee point total count are assigned FDR values of zero.
}

\author{
Aaron Lun
}

\examples{
# Mocking up some data: generates a count matrix named 'my.counts'.
source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
class(my.counts)

# Identify likely cell-containing droplets. 
out <- emptyDrops(my.counts)
out

is.cell <- out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)

# Check if p-values are lower-bounded by 'npts'
# (increase 'npts' if any Limited==TRUE and Sig==FALSE)
table(Sig=is.cell, Limited=out$Limited)
}

\references{
Phipson B, Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Stat. Appl. Genet. Mol. Biol.} 9:Article 39.
}

\seealso{
\code{\link{barcodeRanks}}
}
