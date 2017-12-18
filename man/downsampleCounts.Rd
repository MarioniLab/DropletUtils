\name{downsampleCounts}
\alias{downsampleCounts}

\title{Downsample a count matrix}
\description{
Downsample a count matrix to a desired proportion for each cell.
}

\usage{
downsampleCounts(x, prop)
}

\arguments{
\item{x}{A numeric matrix of counts.}
\item{prop}{A numeric scalar or vector of length \code{ncol(x)}.
All values should lie in [0, 1] indicating the downsampling proportion for each cell.}
}

\value{
An integer matrix of downsampled counts.
}

\details{
Given multiple 10X batches of very different sequencing depths, it can be beneficial to downsample the deepest batches to match the coverage of the shallowest batches. 
This avoids differences in technical noise that can drive clustering by batch. 
Downsampling without replacement is performed on the count vector for each cell,
yielding a new count vector where the total is equal to \code{prop} times the original total count. 
Each count in the returned matrix is guaranteed to be smaller than the original value in \code{x}.

This function provides an alternative to downsampling with the CellRanger \code{aggr} function or subsampling with the 10X Genomics R kit.
Notably, this function downsamples the UMI count matrix, rather than the reads used to obtain the UMI counts in the first place.
Technically, downsampling on the reads is more appropriate as it recapitulates the effect of differences in sequencing depth per cell.
However, in practice, the aim is to obtain cells that have similar total counts across batches, for which downsampling on the UMI counts is a more direct approach.
}

\author{
Aaron Lun
}

\examples{
example(read10xResults)
downsampled <- downsampleCounts(counts(sce10x), prop = 0.5)
}
