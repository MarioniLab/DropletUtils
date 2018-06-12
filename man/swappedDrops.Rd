\name{swappedDrops}
\alias{swappedDrops}

\title{Clean barcode-swapped droplet data}

\description{Remove the effects of barcode swapping on droplet-based single-cell RNA-seq data, specifically 10X Genomics datasets.}

\usage{
swappedDrops(samples, barcode.length=NULL, min.frac=0.8, get.swapped=FALSE, 
    get.diagnostics=FALSE, hdf5.out=TRUE) 
}

\arguments{
\item{samples}{A character vector containing paths to the molecule information HDF5 files, produced by CellRanger for 10X Genomics data.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.}
\item{min.frac}{A numeric scalar specifying the minimum fraction of reads required for a swapped molecule to be assigned to a sample.}
\item{get.swapped}{A logical scalar indicating whether the UMI counts corresponding to swapped molecules should also be returned.}
\item{get.diagnostics}{A logical scalar indicating whether to return the number of reads for each swapped molecule in each sample.}
\item{hdf5.out}{A logical scalar indicating whether the diagnotic matrix should be returned as a \linkS4class{HDF5Matrix}.}
}

\details{
Barcode swapping on the Illumina sequencer occurs when multiplexed samples undergo PCR re-amplification on the flow cell by excess primer with different barcodes.
This results in sequencing of the wrong barcode and molecules being assigned to incorrect samples after debarcoding.
With droplet data, there is the opportunity to remove such effects based on the combination of gene, UMI and cell barcode for each observed transcript molecule.
It is very unlikely that the same combination will arise from different molecules in multiple samples.
Thus, observation of the same combination across multiple samples is indicative of barcode swapping.

For 10X Genomics data, the number of reads corresponding to each gene-UMI-barcode combination is extracted from the molecule information file.
From the total number of reads assigned to that combination, the fraction of reads in each sample is calculated.
The sample with the largest fraction that is greater than \code{min.frac} is defined as the putative sample of origin to which the molecule is assigned.
This assumes that the swapping rate is low, so the sample of origin for a molecule should contain the majority of the reads.
In other all samples, reads for the combination are assumed to derive from swapping and do not contribute to the count matrix.
Setting \code{min.frac=1} will effectively remove all molecules that appear in multiple samples.
We do not recommend setting \code{min.frac} lower than 0.5.

If \code{diagnostics=TRUE}, a diagnostics matrix is returned containing the number of reads per gene-UMI-barcode combination in each sample.
Each row corresponds to a combination and each column corresponds to a sample.
This can be useful for examining the level of swapping across samples on a molecule-by-molecule basis, 
though for the sake of memory, the actual identity of the molecules is not returned.
By default, the matrix is returned as a \linkS4class{HDF5Matrix}, which reduces memory usage and avoids potential issues with integer overflow.
If \code{hdf5.out=FALSE}, a sparse matrix is returned instead, which is faster but uses more memory.
} 

\value{
A list is returned with the \code{cleaned} entry, itself a list of sparse matrices is returned.
Each matrix corresponds to a sample and contains the UMI count for each gene (row) and cell barcode (column) after removing swapped molecules.
All cell barcodes that were originally observed are reported as columns, though note that it is possible for some barcodes to contain no counts.

If \code{get.swapped=TRUE}, a \code{swapped} entry is returned in the top-level list.
This is a list containing sample-specific sparse matrices of UMI counts corresponding to the swapped molecules.
Adding the cleaned and swapped matrices for each sample should yield the total UMI count prior to removal of swapped molecules. 

If \code{get.diagnostics=TRUE}, the top-level list will also contain an additional \code{diagnostics} matrix.
}

\author{
Jonathan Griffiths,
with modifications by Aaron Lun
}

\seealso{
\code{\link{read10xMolInfo}}
}

\examples{
# Mocking up some 10x HDF5-formatted data, with swapping.
curfiles <- DropletUtils:::sim10xMolInfo(tempfile(), nsamples=3)

# Obtaining count matrices with swapping removed.
out <- swappedDrops(curfiles)
lapply(out$cleaned, dim)

out <- swappedDrops(curfiles, get.swapped=TRUE, get.diagnostics=TRUE)
names(out)
}

\references{
Griffiths JA, Lun ATL, Richard AC, Bach K, Marioni JC (2017).
Detection and removal of barcode swapping in single-cell RNA-seq data.
\emph{biorXiv} 177048.
}
