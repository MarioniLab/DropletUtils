#' @export
downsampleReads <- function(sample, prop, barcode.length=NULL, bycol=FALSE) 
# Downsamples the reads to produce a matrix of UMI counts, given
# a HDF5 file containing molecule information.
#
# written by Aaron Lun
# created 19 December 2017    
{
    incoming <- read10xMolInfo(sample, barcode.length, get.umi=FALSE)
    
    out <- .get_cell_ordering(incoming$data$cell, incoming$data$gem_group)
    o <- out$order
    cell.id <- out$id
    unique.cells <- out$cell
    unique.gems <- out$gem
    run.length <- out$length

    if (bycol) {
        prop <- rep(prop, length.out=length(run.length))
    }
    new.read.counts <- downsample_runs(run.length, incoming$data$reads[o], prop, bycol)

    # Assembling into a sparse matrix.
    keep <- new.read.counts>0L
    all.cells <- sprintf("%s-%i", unique.cells, unique.gems)
    makeCountMatrix(incoming$data$gene[o[keep]], cell.id[keep], 
                    all.genes=incoming$genes, all.cells=all.cells)
}
