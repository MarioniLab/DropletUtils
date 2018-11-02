\name{DropletUtilsnews}
\title{DropletUtils News}
\encoding{UTF-8}

\section{Version 1.4.0}{\itemize{
\item Removed read10xMatrix().
}}

\section{Version 1.2.0}{\itemize{
\item Added removeSwappedDrops() for removing swapping in other types of droplet-based data.

\item Added alpha= argument to testEmptyDrops() to support overdispersion during sampling.
Returned arguments and estimates in metadata of testEmptyDrops(), emptyDrops().

\item Added encodeSequences() for convenient 2-bit encoding of sequences.

\item Added get10xMolInfoStats() function to compute per-cell statistics from a molecule info file.

\item Deprecated read10xMatrix(), as it does not add much practical value over Matrix::readMM().

\item Support the 10X sparse HDF5 format in read10xCounts().

\item Support the 10X sparse HDF5 format in write10xCounts().
}}

\section{Version 1.0.0}{\itemize{
\item
New package DropletUtils, for handling droplet-based single-cell RNA sequencing data.
}}

