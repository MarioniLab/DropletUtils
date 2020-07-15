\name{get10xMolInfoStats}
\alias{get10xMolInfoStats}

\title{Get 10x cell statistics}
\description{Compute some basic per-cell statistics from the 10x molecule information file.}

\usage{
get10xMolInfoStats(sample, barcode.length=NULL, use.library=NULL)
}

\arguments{
\item{sample}{A string containing the path to the molecule information HDF5 file.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.}
\item{use.library}{An integer vector specifying the library indices for which to extract molecules from \code{sample}.
Alternatively, a character vector specifying the library type(s), e.g., \code{"Gene expression"}.}
}

\value{
A DataFrame containing one row per cell library, with the fields:
\describe{
    \item{\code{cell}:}{Character, the cell barcode.}
    \item{\code{gem_group}:}{Integer, the GEM group.}
    \item{\code{num.umis}:}{Integer, the number of UMIs assigned to this cell barcode/GEM group combination.}
    \item{\code{num.reads}:}{Integer, the number of reads for this combination.}
    \item{\code{num.genes}:}{Integer, the number of detected genes.}
}
}

\author{
    Aaron Lun
}

\seealso{
    \code{\link{read10xMolInfo}}
}

\examples{
# Mocking up some 10X HDF5-formatted data.
out <- DropletUtils:::simBasicMolInfo(tempfile())

get10xMolInfoStats(out)    
}
