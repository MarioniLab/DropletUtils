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

Downsampling without replacement is performed on the counts in each cell to generate the output matrix. 
Each count in the returned matrix is guaranteed to be smaller than the original value in \code{x}. 
This provides an alternative to downsampling with the CellRanger \code{aggr} function.
}

\examples{
example(read10xResults)
downsampled <- downsampleCounts(counts(sce10x), prop = 0.5)
}
