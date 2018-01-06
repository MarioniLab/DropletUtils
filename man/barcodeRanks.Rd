\name{barcodeRanks}
\alias{barcodeRanks}

\title{Calculate barcode ranks}
\description{Compute barcode rank statistics and identifry the knee and inflection points on the total count curve.}

\usage{
barcodeRanks(m, lower=100, ...) 
}

\arguments{
\item{m}{A real sparse matrix object, either a dgTMatrix or dgCMatrix.
Columns represent barcoded droplets, rows represent cells.
}
\item{lower}{A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.}
\item{...}{Further arguments to pass to \code{\link{smooth.spline}}.}
}

\details{
Analyses of droplet-based scRNA-seq data often show a plot of the log-total count against the log-rank of each barcode, where the highest ranks have the largest totals.
This is equivalent to a transposed empirical cumulative density plot with log-transformed axes.
It allows examination of the distribution of total counts across barcodes, focusing on the barcodes with the largest counts.

Here, ranks are computed by averaging the ranks for all droplets with the same total count.
The inflection point is defined as the point on the curve where the second derivative is maximized,
while the knee point is defined as the point on the curve where the curvature is maximized.
Both of these points correspond to a sharp transition between two components of the total count distribution, 
presumably reflecting the difference between empty droplets with little RNA and cell-containing droplets with much more RNA.

Identification of the knee and inflection points requires calculation of the numerical derivatives.
However, doing this directly on the log-totals yields unstable estimates of the derivatives, due to error propagation in the presence of noise in the curve.
To overcome this, a smooth spline is fitted to the log-total counts against the log-rank.
Derivatives are then calculated from the fitted spline using \code{\link{predict}}.
The curve becomes more erratic at lower counts, so only barcodes with counts above \code{lower} are used to fit the spline. 
}

\value{
A list with the following elements:
\describe{
\item{rank}{A numeric vector of average ranks for each column of \code{m}.}
\item{total}{A numeric vector of total counts for each column of \code{m}.}
\item{fitted}{A numeric vector of fitted total counts from the spline for each column of \code{m}.}
\item{knee}{A numeric scalar containing the total count at the knee point.}
\item{inflection}{A numeric scalar containing the total count at the inflection point.}
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
