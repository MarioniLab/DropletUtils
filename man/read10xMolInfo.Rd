\name{read10xMolInfo}
\alias{read10xMolInfo}

\title{Read the 10X molecule information file}
\description{Extract relevant fields from the molecule information HDF5 file, produced by CellRanger for 10X Genomics data.}

\usage{
read10xMolInfo(sample, barcode.length=NULL, keep.unmapped=FALSE, get.cell=TRUE,
    get.umi=TRUE, get.gem=TRUE, get.gene=TRUE, get.reads=TRUE,
    version=c("auto", "2", "3"))
}

\arguments{
\item{sample}{A string containing the path to the molecule information HDF5 file.}
\item{barcode.length}{An integer scalar specifying the length of the cell barcode.
Only relevant when \code{version="2"}.}
\item{keep.unmapped}{A logical scalar indicating whether unmapped molecules should be reported.}
\item{get.cell, get.umi, get.gem, get.gene, get.reads}{Logical scalar indicating whether the corresponding field should be extracted.}
\item{version}{String specifying the version of the 10X molecule information format to read data from.}
}

\value{
A list is returned containing two elements.
The first element is named \code{data} and is a DataFrame where each row corresponds to a single transcript molecule.
The fields are as follows:
\describe{
\item{\code{barcode}:}{Character, the cell barcode for each molecule.}
\item{\code{umi}:}{Integer, the processed UMI barcode in 2-bit encoding.} 
\item{\code{gem_group}:}{Integer, the GEM group.}
\item{\code{gene}:}{Integer, the index of the gene to which the molecule was assigned.
This refers to an entry in the \code{genes} vector, see below.}
\item{\code{reads}:}{Integer, the number of reads mapped to this molecule.}
}

The field will not be present in the DataFrame if the corresponding \code{get.*} argument is \code{FALSE}, 

The second element of the list is named \code{genes} and is a character vector containing the names of all genes in the annotation.
This contains the names of the various entries of \code{gene} for the individual molecules.
}

\details{
Molecules that were not assigned to any gene have \code{gene} set to \code{length(genes)+1}.
By default, these are removed when \code{keep.unmapped=FALSE}.

CellRanger 3.0 introduced a major change in the format of the molecule information files.
When \code{version="auto"}, the function will attempt to determine the version format of the file.
This can also be user-specified by setting \code{version} explicitly.

For files produced by version 2.2 of the CellRanger software, the length of the cell barcode is not given.
Instead, the barcode length is automatically inferred if \code{barcode.length=NULL} and \code{version="2"}.
Currently, version 1 of the 10X chemistry uses 14 nt barcodes, while version 2 uses 16 nt barcodes.

Setting any of the \code{get.*} arguments will (generally) avoid extraction of the corresponding field.
This can improve efficiency if that field is not necessary for further analysis.
Aside from the missing field, the results are guaranteed to be identical, i.e., same order and number of rows.
Note that some fields must be loaded to yield similar results, e.g., gene IDs will always be loaded to remove unmapped reads if \code{keep.unmapped=FALSE}.
}

\author{
Aaron Lun,
based on code by Jonathan Griffiths
}

\seealso{
\code{\link{makeCountMatrix}}
}

\examples{
# Mocking up some 10X HDF5-formatted data.
out <- DropletUtils:::sim10xMolInfo(tempfile())

# Reading the resulting file.
read10xMolInfo(out)
}

\references{
Zheng GX, Terry JM, Belgrader P, and others (2017).
Massively parallel digital transcriptional profiling of single cells. 
\emph{Nat Commun} 8:14049.

10X Genomics (2017).
Molecule info.
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/output/molecule_info}

10X Genomics (2018).
Molecule info.
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info}
}
