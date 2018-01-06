\name{defaultDrops}
\alias{defaultDrops}

\title{Call cells from number of UMIs}
\description{Call cells according to the number of UMIs associated with each barcode, as implemented in CellRanger.}

\usage{
defaultDrops(m, expected=3000, upper.quant=0.99, lower.prop=0.1)
}

\arguments{
\item{m}{A real sparse matrix object, either a dgTMatrix or dgCMatrix.
Columns represent barcoded droplets, rows represent cells.
The matrix should correspond to an individual sample.
}
\item{expected}{A numeric scalar specifying the expected number of cells in this sample, as specified in the call to CellRanger.}
\item{upper.quant}{A numeric scalar between 0 and 1 specifying the quantile of the top \code{expected} barcodes to consider for the first step of the algorithm}
\item{lower.prop}{A numeric scalar between 0 and 1 specifying the fraction of molecules of the \code{upper.quant} quantile result that a barcode must exceed to be called as a cell}
}

\details{
The \code{defaultDrops} function will call cells based on library size similarly to the CellRanger software suite from 10X Genomics.
Default arguments correspond to an exact reproduction of CellRanger's algorithm, where the number of expected cells was also left at CellRanger default value.

The method considers the \code{upper.quant} quantile of the top \code{expected} barcodes, ordered by decreasing number of UMIs, as a threshold. Any barcodes containing more molecules than \code{lower.prop} times this threshold is considered to be a cell, and is retained for further analysis.

This method may be vulnerable to calling very well-captured background RNA, or missing real cells that were poorly captured. See \code{?\link{emptyDrops}} for an alternative approach.
}

\value{
\code{defaultDrops} will return a logical vector of length \code{ncol(m)}.
Each element of the vector reports whether each column of \code{m} was called as a cell.
}

\author{
Jonathan Griffiths
}

\examples{
# Mocking up some data: 
set.seed(0)
my.counts <- DropletUtils:::simCounts()

# Identify likely cell-containing droplets. 
called <- defaultDrops(my.counts)
table(called)

# Get matrix of called cells.
cell.counts <- my.counts[, called]
}

\references{
10X Genomics (2017).
Cell Ranger Algorithms Overview. 
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview}
}

\seealso{
\code{\link{emptyDrops}}
}
