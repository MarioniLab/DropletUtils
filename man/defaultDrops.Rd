\name{cellrangerCall}
\alias{cellrangerCall}

\title{Call cells from number of UMIs}
\description{Call cells according to the number of UMIs associated wiht each barcode, like CellRanger's approach}

\usage{
cellrangerCall(m, exp.cells = 3000, upper = 0.99, lower = 0.1)
}

\arguments{
\item{m}{A real sparse matrix object, either a dgTMatrix or dgCMatrix.
Columns represent barcoded droplets, rows represent cells.
The matrix should correspond to an individual sample.
}
\item{exp.cells}{A numeric scalar specifying the expected number of cells in this sample.}
\item{upper}{A numeric scalar between 0 and 1 specifying the quantile of the top \code{exp.cells} to consider for the first step of the algorithm}
\item{lower}{A numeric scalar between 0 and 1 specifying the fraction of molecules of the \code{upper} quantile result for a barcode to be called as a cell}
}

\details{
The \code{testEmptyDrops} function will call cells based on library size similarly to CellRanger.
Default arguments correspond to an exact reproduction of CellRanger's algorithm, where the number of expected cells was also left at CellRanger default value.

The method considers the \code{upper} quantile top \code{exp.cells} barcodes, ordered by decreasing number of UMIs, as a threshold. Any barcodes containing more molecules than \code{lower} times this threshold is considered to be a cell, and is retained for further analysis.

This method may be vulnerable to calling very well-captured background RNA, or missing real cells that were poorly captured. See \code{?emptyDrops} for an alternative approach.
}

\value{
\code{cellrangerCall} will return a logical vector of length \code{ncol(m)}.
Each element of the vector reports whether each column of \code{m} was called as a cell.
}

\author{
Jonathan Griffiths
}

\examples{
# Mocking up some data: generates a count matrix named 'my.counts'.
source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
class(my.counts)

# Identify likely cell-containing droplets. 
called <- cellrangerCall(my.counts)
table(called)

# Get matrix of called cells
cell.counts <- my.counts[, called]
}

\references{
10X Genomics. 2018. 
Cell Ranger Algorithms Overview. Available at: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview. [Accessed 4 January 2018].
}

\seealso{
\code{\link{emptyDrops}}
}
