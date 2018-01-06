\name{downsampleReads}
\alias{downsampleReads}

\title{Downsample reads in a 10X Genomics dataset}
\description{Generate a UMI count matrix after downsampling reads from the molecule information file produced by CellRanger for 10X Genomics data.}

\usage{
downsampleReads(sample, prop, barcode.length=NULL, bycol=FALSE)
}

\arguments{
\item{sample}{A string containing the path to the molecule information HDF5 file.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.}
\item{prop}{A numeric scalar or, if \code{bycol=TRUE}, a vector of length \code{ncol(x)}.
All values should lie in [0, 1] specifying the downsampling proportion for the matrix or for each cell.}
\item{bycol}{A logical scalar indicating whether downsampling should be performed on a column-by-column basis.}
}

\details{
This function downsamples the reads for each molecule by the specified \code{prop}, using the information in \code{sample}.
It then constructs a UMI count matrix based on the molecules with non-zero read counts.
The aim is to eliminate differences in technical noise that can drive clustering by batch, as described in \code{\link{downsampleMatrix}}.

Subsampling the reads with \code{downsampleReads} recapitulates the effect of differences in sequencing depth per cell.
This provides an alternative to downsampling with the CellRanger \code{aggr} function or subsampling with the 10X Genomics R kit.
Note that this differs from directly subsampling the UMI count matrix with \code{\link{downsampleMatrix}}.

If \code{bycol=TRUE}, sampling without replacement is performed on the reads for each cell.
The total number of reads for each cell after downsampling is guaranteed to be \code{prop} times the original total.
Different proportions can be specified for different cells by setting \code{prop} to a vector.

If \code{bycol=FALSE}, downsampling without replacement is performed on all reads from the entire dataset.
The total number of reads for each cell after downsampling may not be exactly equal to \code{prop} times the original value.
Note that this is the more natural approach and is the default, which differs from the default used in \code{\link{downsampleMatrix}}.
}

\value{
A numeric sparse matrix containing the downsampled UMI counts for each gene (row) and barcode (column).
}

\seealso{
\code{\link{downsampleMatrix}},
\code{\link{read10xMolInfo}}
}

\author{
Aaron Lun
}

\examples{
# Mocking up some 10X HDF5-formatted data.
out <- DropletUtils:::sim10xMolInfo(tempfile(), nsamples=1)

# Downsampling by the reads.
downsampleReads(out, barcode.length=4, prop=0.5)
}
