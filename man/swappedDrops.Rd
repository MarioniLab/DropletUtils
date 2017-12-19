\name{swappedDrops}
\alias{swappedDrops}

\title{Clean barcode-swapped droplet data}

\description{Remove the effects of barcode swapping on droplet-based single-cell RNA-seq data, specifically 10X Genomics datasets.}

\usage{
swappedDrops(samples, barcode.length, get.swapped=FALSE, min.frac=0.8)
}

\arguments{
\item{samples}{A character vector containing paths to the molecule information HDF5 files, produced by CellRanger for 10X Genomics data.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode.}
\item{get.swapped}{A logical scalar indicating whether the counts for swapped molecules should also be returned.}
\item{min.frac}{A numeric scalar specifying the minimum fraction of reads required for a swapped molecule to be assigned to a sample.}
}

\details{
Barcode swapping on the Illumina sequencer occurs when multiplexed samples undergo PCR re-amplification on the flow cell by excess primer with different barcodes.
This results in sequencing of the wrong barcode and molecules being assigned to incorrect samples after debarcoding.
With droplet data, there is the opportunity to remove such effects based on the combination of gene, UMI and cell barcode for each observed transcript molecule.
It is very unlikely that the same combination will arise from different molecules in multiple samples.
Thus, observation of the same combination across multiple samples is indicative of barcode swapping.

For each potentially swapped molecule in 10X Genomics data, the number of reads corresponding to that molecule is extracted from the molecule information file.
The fraction of reads in each sample is calculated, and the molecule is assigned to the sample with the largest fraction if it is greater than \code{min.frac}.
This assumes that the swapping rate is low, so the true sample of origin for a molecule should contain the majority of the reads.
Setting \code{min.frac=1} will effectively remove all molecules that appear in multiple samples.
We do not recommend setting \code{min.frac} lower than 0.5.

The \code{barcode.length} parameter needs to be set manually as the length of the cell barcode changes between versions of the 10X Genomics chemistry.
} 

\value{
If \code{get.swapped=FALSE}, a list of sparse matrices.
Each matrix corresponds to a sample and contains the UMI counts for each gene and barcode after removing swapped molecules.
All barcodes that were originally observed are reported as columns, though note that it is possible for some barcodes to contain no counts.

If \code{get.swapped=TRUE}, a list of two lists are returned.
The \code{cleaned} list is the same as described above. 
The \code{swapped} list contains sample-specific sparse matrices of UMI counts corresponding to the swapped molecules for each gene and barcode.
Adding these two matrices for each sample should yield the total UMI count.
}

\author{
Jonathan Griffiths,
with modifications by Aaron Lun
}

\seealso{
\code{\link{sim10xMolInfo}}
}

\examples{
# Setting up some temporary molInfo files.
tmpdir <- tempfile()
dir.create(tmpdir)
curfiles <- sim10xMolInfo(tmpdir)

# Getting swapped matrices, barcode length of 4 bp.
out <- swappedDrops(curfiles, 4)
lapply(out, dim)

out <- swappedDrops(curfiles, 4, get.swapped=TRUE)
names(out)
}
