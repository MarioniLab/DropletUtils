\name{read10xMatrix}
\alias{read10xMatrix}

\title{Read in the 10x count matrix}

\description{
Creates a sparse or HDF5-backed count matrix from the MatrixMarket file produced by CellRanger.
}

\usage{
read10xMatrix(file, hdf5.out=FALSE, chunk.size)
}

\arguments{
\item{file}{String containing the path to a MatrixMarket file, usually named \code{"matrix.mtx"}.}
\item{hdf5.out}{A logical scalar indicating whether a \code{\link[HDF5Array:HDF5Array-class]{HDF5Matrix}} object should be produced.}
\item{chunk.size}{An integer scalar specifying the chunk size when reading in records from \code{file}.}
}

\value{
A \code{dgCMatrix} object (or a \code{HDF5Matrix} object, if \code{hdf5.out=TRUE}) containing the counts for each gene (row) and cell barcode (column).
}

\details{
When \code{hdf5.out=FALSE}, \code{\link{readMM}} is used directly.
However, for very large 10x experiments with more than \code{.Machine$integer.max} non-zero entries, \code{dgCMatrix} may encounter integer overflows.
In such cases, setting \code{hdf5.out=TRUE} will produce a \code{HDF5Matrix} object instead.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{readMM}},
\code{\link{read10xCounts}}
}

\examples{
# Mocking up some 10X genomics output.
example(write10xCounts)

mm.path <- file.path(tmpdir, "matrix.mtx")
X <- read10xMatrix(mm.path)
altX <- read10xMatrix(mm.path, chunk.size=10, hdf5.out=TRUE)
}

\references{
10X Genomics (2017).
Gene-Barcode Matrices.
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices}
}
