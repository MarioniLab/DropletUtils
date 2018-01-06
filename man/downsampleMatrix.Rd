\name{downsampleMatrix}
\alias{downsampleMatrix}

\title{Downsample a count matrix}
\description{
Downsample a count matrix to a desired proportion for each cell.
}

\usage{
downsampleMatrix(x, prop, bycol=TRUE)
}

\arguments{
\item{x}{A numeric matrix of counts.}
\item{prop}{A numeric scalar or, if \code{bycol=TRUE}, a vector of length \code{ncol(x)}.
All values should lie in [0, 1] specifying the downsampling proportion for the matrix or for each cell.}
\item{bycol}{A logical scalar indicating whether downsampling should be performed on a column-by-column basis.}
}

\value{
A numeric matrix of downsampled counts, of the same type as \code{x}.
}

\details{
Given multiple batches of very different sequencing depths, it can be beneficial to downsample the deepest batches to match the coverage of the shallowest batches. 
This avoids differences in technical noise that can drive clustering by batch.

If \code{bycol=TRUE}, sampling without replacement is performed on the count vector for each cell.
This yields a new count vector where the total is equal to \code{prop} times the original total count. 
Each count in the returned matrix is guaranteed to be smaller than the original value in \code{x}.
Different proportions can be specified for different cells by setting \code{prop} to a vector.

If \code{bycol=FALSE}, downsampling without replacement is performed on the entire matrix.
This yields a new matrix where the total count across all cells is equal to \code{prop} times the original total.
The new total count for each cell may not be exactly equal to \code{prop} times the original value,
which may or may not be more appropriate than \code{bycol=TRUE} for particular applications.

Technically, downsampling on the reads with \code{\link{downsampleReads}} is more appropriate as it recapitulates the effect of differences in sequencing depth per cell.
However, in practice, the aim is to obtain cells that have similar total counts across batches, for which downsampling on the UMI counts is a more direct approach.

Note that this function was originally implemented in the \pkg{scater} package as \code{downsampleCounts}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{downsampleReads}}
}

\examples{
example(read10xCounts)
downsampled <- downsampleMatrix(counts(sce10x), prop = 0.5)
}
