#' Downsample reads in a 10X Genomics dataset
#'
#' Generate a UMI count matrix after downsampling reads from the molecule information file 
#' produced by CellRanger for 10X Genomics data.
#' 
#' @param sample A string containing the path to the molecule information HDF5 file.
#' @param barcode.length An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.
#' @param prop A numeric scalar or, if \code{bycol=TRUE}, a vector of length \code{ncol(x)}.
#' All values should lie in [0, 1] specifying the downsampling proportion for the matrix or for each cell.
#' @param bycol A logical scalar indicating whether downsampling should be performed on a column-by-column basis.
#' @param features A character vector containing the names of the features on which to perform downsampling.
#' @param use.library An integer vector specifying the library indices for which to extract molecules from \code{sample}.
#' Alternatively, a character vector specifying the library type(s), e.g., \code{"Gene expression"}.
#' 
#' @details
#' This function downsamples the reads for each molecule by the specified \code{prop}, using the information in \code{sample}.
#' It then constructs a UMI count matrix based on the molecules with non-zero read counts.
#' The aim is to eliminate differences in technical noise that can drive clustering by batch, as described in \code{\link{downsampleMatrix}}.
#' 
#' Subsampling the reads with \code{downsampleReads} recapitulates the effect of differences in sequencing depth per cell.
#' This provides an alternative to downsampling with the CellRanger \code{aggr} function or subsampling with the 10X Genomics R kit.
#' Note that this differs from directly subsampling the UMI count matrix with \code{\link{downsampleMatrix}}.
#' 
#' If \code{bycol=FALSE}, downsampling without replacement is performed on all reads from the entire dataset.
#' The total number of reads for each cell after downsampling may not be exactly equal to \code{prop} times the original value.
#' Note that this is the more natural approach and is the default, which differs from the default used in \code{\link{downsampleMatrix}}.
#' 
#' If \code{bycol=TRUE}, sampling without replacement is performed on the reads for each cell.
#' The total number of reads for each cell after downsampling is guaranteed to be \code{prop} times the original total (rounded to the nearest integer).
#' Different proportions can be specified for different cells by setting \code{prop} to a vector, 
#' where each proportion corresponds to a cell/GEM combination in the order returned by \code{\link{get10xMolInfoStats}}.
#'
#' The \code{use.library} argument is intended for studies with multiple feature types, e.g., antibody capture or CRISPR tags.
#' As the reads for each feature type are generated in a separate sequencing library, it is generally most appropriate to downsample reads for each feature type separately.
#' This can be achieved by setting \code{use.library} to the name or index of the desired feature set.
#' The features of interest can also be directly specified with \code{features}.
#' (This will be intersected with any \code{use.library} choice if both are specified.)
#' 
#' @return
#' A numeric sparse matrix containing the downsampled UMI counts for each gene (row) and barcode (column).
#' If \code{features} is set, only the rows with names in \code{features} are returned.
#' 
#' @seealso
#' \code{\link{downsampleMatrix}}, for more general downsampling of the count matrix.
#'
#' \code{\link{read10xMolInfo}}, to read the contents of the molecule information file.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' # Mocking up some 10X HDF5-formatted data.
#' out <- DropletUtils:::simBasicMolInfo(tempfile())
#' 
#' # Downsampling by the reads.
#' downsampleReads(out, barcode.length=4, prop=0.5)
#' 
#'
#' @export
downsampleReads <- function(sample, prop, barcode.length=NULL, bycol=FALSE, features=NULL, use.library=NULL) 
# Downsamples the reads to produce a matrix of UMI counts, given
# a HDF5 file containing molecule information.
#
# written by Aaron Lun
# created 19 December 2017    
{
    incoming <- .extract_mol_info(sample, barcode.length=barcode.length, get.umi=FALSE, 
        use.library=use.library, subset.library.features=TRUE)

    if (!is.null(features)) {
        incoming <- .reindex_mol_info_features(incoming, incoming$genes %in% features)
    }
    
    out <- .get_cell_ordering(incoming$data$cell, incoming$data$gem_group)
    o <- out$order
    cell.id <- out$id
    unique.cells <- out$cell
    unique.gems <- out$gem
    run.length <- out$length

    if (bycol) {
        prop <- rep(prop, length.out=length(run.length))
        new.read.counts <- downsample_run_per_cell(run.length, incoming$data$reads[o], prop)
    } else {
        new.read.counts <- downsample_run(incoming$data$reads[o], prop)
    }

    # Assembling into a sparse matrix.
    keep <- new.read.counts>0L
    all.cells <- sprintf("%s-%i", unique.cells, unique.gems)
    makeCountMatrix(incoming$data$gene[o[keep]], cell.id[keep], 
                    all.genes=incoming$genes, all.cells=all.cells)
}
