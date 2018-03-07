#' @export
#' @importFrom rhdf5 h5read
#' @importFrom S4Vectors DataFrame
read10xMolInfo <- function(sample, barcode.length=NULL, keep.unmapped=FALSE) 
# Utility function to read useful information from a 10X molecule information file.
#
# written by Aaron Lun
# based on code from Jonathan Griffiths
# created 20 December 2017    
{
    cell <- .Call(cxx_get_cell_barcodes, sample, "barcode", barcode.length)
    umi <- as.vector(h5read(sample, "/umi"))
    gem_group <- as.vector(h5read(sample, "/gem_group"))
    gene <- as.vector(h5read(sample, "/gene")) + 1L # Zero-indexed by default
    reads <- as.vector(h5read(sample, "/reads")) # Maybe useful for selective exclusion
    data <- DataFrame(cell=cell, umi=umi, gene=gene, reads=reads, gem_group=gem_group)

    # Defining the set of all barcodes, and of all genes.
    all.barcodes <- sort(unique(cell))
    gene.ids <- h5read(sample, "/gene_ids") 

    # Remove the unassigned gene entries.
    if (!keep.unmapped) { 
        keep <- gene <= length(gene.ids)
        data <- data[keep,]
    }

    # Don't define the total cell pool here, as higher level functions may want to use gem_group.
    return(list(data=data, genes=gene.ids))
}

