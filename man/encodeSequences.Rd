\name{encodeSequences}
\alias{encodeSequences}

\title{Encode nucleotide sequences}
\description{Encode short nucleotide sequences into integers with a 2-bit encoding.}

\usage{
encodeSequences(sequences)
}

\arguments{
\item{sequences}{A character vector of short nucleotide sequences, e.g., UMIs or cell barcodes.}
}

\details{
Each pair of bits encodes a nucleotide - 00 is A, 01 is C, 10 is G and 11 is T.
The least significant byte contains the 3'-most nucleotides, and the remaining bits are set to zero.
Thus, the sequence \dQuote{CGGACT} is converted to the binary form:
\preformatted{    01 10 10 00 01 11
}
... which corresponds to the integer 1671.

A consequence of R's use of 32-bit integers means that no element of \code{sequences} can be more than 15 nt long.
Otherwise, integer overflow will occur.
}

\value{
An integer vector containing the encoded sequences.
}

\author{
Aaron Lun
}

\references{
10X Genomics (2017).
Molecule info.
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info}
}

\examples{
encodeSequences("CGGACT")
}
