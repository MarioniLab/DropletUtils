\name{makeCountMatrix}
\alias{makeCountMatrix}

\title{Make a count matrix}
\description{Construct a count matrix from per-molecule information, typically the cell and gene of origin.}

\usage{
makeCountMatrix(gene, cell, all.genes=NULL, all.cells=NULL, value=NULL) 
}

\arguments{
\item{gene}{An integer or character vector specifying the gene to which each molecule was assigned.}
\item{cell}{An integer or character vector specifying the cell to which each molecule was assigned.}
\item{all.genes}{A character vector containing the names of all genes in the dataset.}
\item{all.cells}{A character vector containing the names of all cells in the dataset.}
\item{value}{A numeric vector containing values for each molecule.}
}

\details{
Each element of the vectors \code{gene}, \code{cell} and (if specified) \code{value} contain information for a single transcript molecule.
Each entry of the output matrix corresponds to a single \code{gene} and \code{cell} combination.
If multiple molecules are present with the same combination, their values in \code{value} are summed together, and the sum is used as the entry of the output matrix.

If \code{value=NULL}, it will default to a vector of all 1's.
Each entry of the output matrix represents the number of molecules with the corresponding combination, i.e., UMI counts.
Users can pass other metrics such as the number of reads covering each molecule.
This would yield a read count matrix rather than a UMI count matrix.

If \code{all.genes} is not specified, it is kept as \code{NULL} for integer \code{gene}.
Otherwise, it is defined as the sorted unique values of character \code{gene}.
The same occurs for \code{cell} and \code{all.cells}.

If \code{gene} is integer, its values should be positive and no greater than \code{length(all.genes)} if \code{all.genes!=NULL}.
If \code{gene} is character, its values should be a subset of those in \code{all.genes}.
The same requirements apply to \code{cell} and \code{all.cells}.
}

\value{
A sparse matrix where rows are genes, columns are cells and entries are the sum of \code{value} for each gene/cell combination.
Rows and columns are named if the \code{gene} or \code{cell} are character or if \code{all.genes} or \code{all.cells} are specified.
}

\author{
Aaron Lun
}

\examples{
nmolecules <- 100
gene.id <- sample(LETTERS, nmolecules, replace=TRUE)
cell.id <- sample(20, nmolecules, replace=TRUE)
makeCountMatrix(gene.id, cell.id)
}

\seealso{
\code{\link{read10xMolInfo}}
}

