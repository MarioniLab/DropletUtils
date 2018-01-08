\name{barcodeRanks}
\alias{barcodeRanks}

\title{Calculate barcode ranks}
\description{Compute barcode rank statistics and identifry the knee and inflection points on the total count curve.}

\usage{
barcodeRanks(m, lower=100, fit.bounds=NULL, df=20, ...) 
}

\arguments{
\item{m}{A real sparse matrix object, either a dgTMatrix or dgCMatrix.
Columns represent barcoded droplets, rows represent cells.
}
\item{lower}{A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.}
\item{fit.bounds}{A numeric vector of length 2, specifying the lower and upper bouunds on the total UMI count for spline fitting.}
\item{df, ...}{Further arguments to pass to \code{\link{smooth.spline}}.}
}

\details{
Analyses of droplet-based scRNA-seq data often show a plot of the log-total count against the log-rank of each barcode, where the highest ranks have the largest totals.
This is equivalent to a transposed empirical cumulative density plot with log-transformed axes, which focuses on the barcodes with the largest counts.
The \code{barcodeRanks} function will compute these ranks for all barcodes.
Barcodes with the same total count receive the same average rank to avoid problems with discrete runs of the same total.

The function will also identify a number of interesting points on the curve for downstream use, namely the inflection and knee points.
Both of these points correspond to a sharp transition between two components of the total count distribution, 
presumably reflecting the difference between empty droplets with little RNA and cell-containing droplets with much more RNA.
\itemize{
\item The inflection point is computed as the point on the rank/total curve where the first derivative is minimized.
The derivative is computed directly from all points on the curve with total counts greater than \code{lower}.
This avoids issues with erratic behaviour of the curve at lower totals.
\item The knee point is defined as the point on the curve where the signed curvature is minimized.
This requires calculation of the second derivative, which is much more sensitive to noise in the curve.
To overcome this, a smooth spline is fitted to the log-total counts against the log-rank using the \code{\link{smooth.spline}} function.
Derivatives are then calculated from the fitted spline using \code{\link{predict}}.
}
We supply a default setting of \code{df} to avoid overfitting the spline, which results in unstability in the higher derivatives (and thus the curvature).
This and other arguments to \code{\link{smooth.spline}} can be tuned if the estimated knee point is not at an appropriate location.
We also restrict the fit to lie within the bounds defined by \code{fit.bounds} to focus on the region containing the knee point.
This allows us to obtain an accurate fit with low \code{df}, rather than attempting to model the entire curve.

If \code{fit.bounds} is not specified, the upper bound is automatically set to the inflection point, which should lie after the knee point.
The lower bound is set to the point at which the first derivative is closest to zero, i.e., the \dQuote{plateau} region before the knee point.
Note that only points with total counts above \code{lower} will be considered, regardless of how \code{fit.bounds} is defined.
}

\value{
A list with the following elements:
\describe{
\item{\code{rank}:}{A numeric vector of average ranks for each column of \code{m}.}
\item{\code{total}:}{A numeric vector of total counts for each column of \code{m}.}
\item{\code{fitted}:}{A numeric vector of fitted total counts from the spline for each column of \code{m}.
This is \code{NA} for points with \code{x} outside of \code{fit.bounds}.
}
\item{\code{knee}:}{A numeric scalar containing the total count at the knee point.}
\item{\code{inflection}:}{A numeric scalar containing the total count at the inflection point.}
}
}

\author{
Aaron Lun
}

\examples{
# Mocking up some data: 
set.seed(2000)
my.counts <- DropletUtils:::simCounts()

# Computing barcode rank statistics:
br.out <- barcodeRanks(my.counts)
names(br.out)

# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))
}

\seealso{
\code{\link{emptyDrops}}
}
