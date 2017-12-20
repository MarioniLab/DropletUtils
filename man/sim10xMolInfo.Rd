\name{sim10xMolInfo}
\alias{sim10xMolInfo}

\title{Simulate 10X molecule information}
\description{Create a HDF5 file that mimics the molecule information file produced by CellRanger for 10X Genomics data.}

\usage{
sim10xMolInfo(prefix, nsamples=1, umi.length=10, barcode.length=4, 
    ngenes=20, nmolecules=10000, swap.frac=0.2, 
    ave.read=10, return.tab=FALSE)
}

\arguments{
\item{prefix}{A string containing the prefix to the file path.}
\item{nsamples}{An integer scalar specifying the number of samples to simulate.
Each sample has its own HDF5 file, named using the \code{prefix} followed by a number.}
\item{umi.length}{An integer scalar specifying the length of the UMI.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode.}
\item{ngenes}{An integer scalar specifying the number of genes in the annotation.}
\item{nmolecules}{An integer scalar specifying the number of transcript molecules to simulate.
Each molecule gets its own entry in the HDF5 file.}
\item{swap.frac}{A numeric scalar specifying the proportion of molecules that result from barcode swapping.}
\item{ave.read}{A numeric scalar specifying the mean number of reads per molecule, under a Poisson sampling scheme.}
\item{return.tab}{A logical scalar indicating whether tables of known results should be returned.}
}

\details{
Many of the arguments are designed to test the correct identification of barcode-swapped molecules in \code{\link{swappedDrops}}.
Any molecules that arise due to barcode swapping have the number of reads set to 1.
Molecules correctly assigned to their sample of origin have their number of reads set to 1 plus a sample from a Poisson distribution with mean \code{ave.reads}.
If \code{nsamples=1}, \code{swap.frac} is ignored and no swapped molecules will be generated.

The number of possible cells is calculated as \code{4^barcode.length}, with a similar calculation for the number of possible UMI barcodes.
These are represented in R as 32-bit integers with 2-bit encoding.
As such, the specified lengths for the UMI and cell barcode cannot be greater than 15 nt.

Note that the the \code{barcode} field in the 10X molecule information file uses 64-bit integers.
Unfortunately, this is not easily achieved as R does not natively support 64-bit integers, so we approximate this with a 32-bit integer instead.
}

\author{
Jonathan Griffiths,
with modifications by Aaron Lun
}

\value{
HDF5 files are produced with the specified \code{prefix}.
If \code{return.tab=FALSE}, a character vector is also returned containing the file names.

If \code{return.tab=TRUE}, a named list is returned with the elements \code{files}, the character vector described above;
\code{original}, a DataFrame containing statistics for each unswapped molecule (i.e., 2-bit encoded \code{cell} and \code{umi} barcodes, plus the \code{gene} and \code{sample} IDs);
and \code{swapped}, a DataFrame containing statistics for each swapped molecule.
}

\seealso{
\code{\link{read10xMolInfo}}
}

\examples{
sim10xMolInfo(tempfile(), nsamples=3)

sim10xMolInfo(tempfile(), nsamples=3, return.tab=TRUE)
}
