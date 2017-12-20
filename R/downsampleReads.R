downsampleReads <- function(sample, barcode.length, prop, bycol=FALSE) 
# Downsamples the reads to produce a matrix of UMI counts, given
# a HDF5 file containing molecule information.
#
# written by Aaron Lun
# created 19 December 2017    
{
    incoming <- read10xMolInfo(sample, barcode.length)
    if (length(incoming$data$cell)) { 
        cell.id <- paste0(incoming$data$cell, "-", incoming$data$gem_group)
    } else {
        cell.id <- character(0)
    }
    all.cells <- sort(unique(cell.id))
    
    # Obtaining new read counts for each molecule.
    o <- order(cell.id)
    cell.id <- cell.id[o]
    rout <- rle(cell.id)

    if (bycol) {
        prop <- rep(prop, length.out=length(all.cells))
    }
    new.read.counts <- .Call(cxx_downsample_runs, rout$lengths, incoming$data$reads[o], prop, bycol)

    # Assembling into a sparse matrix.
    keep <- new.read.counts>0L
    rid <- incoming$data$gene[o[keep]]
    cid <- match(cell.id[keep], all.cells)

    out <- sparseMatrix(i=rid, j=cid, x=rep(1, length(rid)),
                        dims=c(length(incoming$genes), length(all.cells)),
                        use.last.ij=FALSE, giveCsparse=TRUE)
    rownames(out) <- incoming$genes
    colnames(out) <- all.cells
    return(out)
}
