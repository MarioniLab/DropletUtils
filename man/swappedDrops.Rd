\name{swappedDrops}
\alias{swappedDrops}

\title{Clean barcode-swapped droplet data}

\description{Remove the effects of barcode swapping on droplet-based single-cell RNA-seq data, specifically 10X Genomics datasets.}

\usage{
swappedDrops(samples, barcode.length=NULL, get.swapped=FALSE, 
    get.diagnostics=FALSE, min.frac=0.8)
}

\arguments{
\item{samples}{A character vector containing paths to the molecule information HDF5 files, produced by CellRanger for 10X Genomics data.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.}
\item{get.swapped}{A logical scalar indicating whether the UMI counts corresponding to swapped molecules should also be returned.}
\item{get.diagnostics}{A logical scalar indicating whether to return the number of reads for each swapped molecule in each sample.}
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
} 

\value{
A list is returned with the \code{cleaned} element, itself a list of sparse matrices is returned.
Each matrix corresponds to a sample and contains the UMI count for each gene (row) and barcode (column) after removing swapped molecules.
All barcodes that were originally observed are reported as columns, though note that it is possible for some barcodes to contain no counts.

If \code{get.swapped=TRUE}, a \code{swapped} element is returned in the top-level list.
This contains sample-specific sparse matrices of UMI counts corresponding to the swapped molecules.
Adding the cleaned and swapped matrices for each sample should yield the total UMI count prior to removal of swapped molecules. 

If \code{get.diagnostics=TRUE}, the top-level list will also contain an additional \code{diagnostics} matrix.
Each entry of this matrix contains the number of reads for each molecule (row) in each sample (column).
This can be useful for examining the level of swapping across samples on a molecule-by-molecule basis.
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
lapply(out, dim)

out <- swappedDrops(curfiles, get.swapped=TRUE, get.diagnostics=TRUE)
names(out)
}

\references{
Griffiths JA, Lun ATL, Richard AC, Bach K, Marioni JC (2017).
Detection and removal of barcode swapping in single-cell RNA-seq data.
\emph{biorXiv} 177048.
}
