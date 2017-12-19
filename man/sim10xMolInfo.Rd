\name{sim10xMolInfo}
\alias{sim10xMolInfo}

\title{Simulate 10X molecule information}
\description{Create a HDF5 file that mimics the molecule information file produced by CellRanger for 10X Genomics data.}

\usage{
sim10xMolInfo(prefix, umi.len=10, barcode.len=4, nmolecules=10000, swap.frac=0.2, 
    ngenes=20, nsamples=3, lambda=10, return.tab=FALSE)
}

\arguments{
\item{prefix}{A string containing the prefix to the file path.}
\item{umi.len}{An integer scalar specifying the length of the UMI.}
\item{barcode.len}{An integer scalar specifying the length of the cell barcode.}
\item{nmolecules}{An integer scalar specifying the number of transcript molecules to simulate.
Each molecule gets its own entry in the HDF5 file.}
\item{swap.frac}{A numeric scalar specifying the proportion of molecules that result from barcode swapping.}
\item{ngenes}{An integer scalar specifying the number of genes in the annotation.}
\item{nsamples}{An integer scalar specifying the number of samples to simulate.
Each sample has its own HDF5 file, named using the \code{prefix} followed by a number.}
\item{lambda}{A numeric scalar specifying the mean number of reads per molecule, under a Poisson sampling scheme.}
\item{return.tab}{A logical scalar indicating whether tables of known results should be returned.}
}

\details{
Many of the arguments are designed to test the correct identification of barcode-swapped molecules in \code{\link{swappedDrops}}.
Any molecules that arise due to barcode swapping have the number of reads set to 1.
Otherwise, molecules that are correctly assigned to their sample of origin have number of reads sampled from a Poisson distribution with mean \code{lambda}.

The number of possible cells is calculated as \code{4^barcode.len}.
Note that the cell barcode is saved as a 32-bit integer rather than as a 64-bit integer.
The latter is the actual specification of the \code{barcode} field in the 10X molecule information file.
Unfortunately, this is not easily achieved as R does not natively support 64-bit integers.
}

\author{
Jonathan Griffiths,
with modifications by Aaron Lun
}

\value{
HDF5 files are produced with the specified \code{prefix}.
If \code{return.tab=FALSE}, a character vector is also returned containing the file names.

If \code{return.tab=TRUE}, a named list is returned with the elements \code{files}, the character vector described above;
\code{original}, a data frame containing statistics for each unswapped molecule;
and \code{swapped}, a data frame containing statistics for each swapped molecule.
}

\seealso{
\code{\link{swappedDrops}}
}

\examples{
tmpdir <- tempfile()
dir.create(tmpdir)
sim10xMolInfo(tmpdir)
}
