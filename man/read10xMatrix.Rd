\name{read10xMatrix}
\alias{read10xMatrix}

\title{Read in the 10x count matrix}

\description{
Creates a sparse or HDF5-backed count matrix from the MatrixMarket file produced by CellRanger.
}

\usage{
read10xMatrix(file, chunk.size=NA, hdf5.out=FALSE)
}

\arguments{
\item{file}{String containing the path to a MatrixMarket file, usually named \code{"matrix.mtx"}.}
\item{chunk.size}{An integer scalar specifying the chunk size when reading in records from \code{file}.}
\item{hdf5.out}{A logical scalar indicating whether a \code{\link[HDF5Array:HDF5Array-class]{HDF5Matrix}} object should be produced.}
}

\value{
A \code{dgCMatrix} object (or a \code{HDF5Matrix} object, if \code{hdf5.out=TRUE}) containing the counts for each gene (row) and cell barcode (column).
}

\details{
When \code{chunk.size=NA}, \code{\link{readMM}} is used directly.
However, for large 10x experiments, chunking may be necessary to avoid problems with limited string lengths when \code{\link{scan}}ning from \code{file}.
This function will read in chunks from \code{file} and collate them to construct the final matrix.

By default, the function will produce a \code{dgCMatrix} object.
If \code{hdf5.out=TRUE} and \code{chunk.size} is not \code{NA}, a \code{HDF5Matrix} object will be produced instead.
This uses file-backed storage and may be more suitable than a sparse representation, especially for very large matrices where \code{dgCMatrix} encounters integer overflows.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{read10xCounts}}
}

\examples{
# Mocking up some 10X genomics output.
example(write10xCounts)

mm.path <- file.path(tmpdir, "matrix.mtx")
X <- read10xMatrix(mm.path)
Y <- read10xMatrix(mm.path, chunk.size=10)
Z <- read10xMatrix(mm.path, chunk.size=10, hdf5.out=TRUE)
}

\references{
10X Genomics (2017).
Gene-Barcode Matrices.
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices}
}
