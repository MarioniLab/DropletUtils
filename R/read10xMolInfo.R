#' @export
#' @importFrom rhdf5 h5read
#' @importFrom S4Vectors DataFrame
read10xMolInfo <- function(sample, barcode.length=NULL, keep.unmapped=FALSE, 
    get.cell=TRUE, get.umi=TRUE, get.gem=TRUE, get.gene=TRUE, get.reads=TRUE)
# Utility function to read useful information from a 10X molecule information file.
#
# written by Aaron Lun
# based on code from Jonathan Griffiths
# created 20 December 2017    
{
    data <- list()

    if (get.cell) {
        data$cell <- .Call(cxx_get_cell_barcodes, sample, "barcode", barcode.length)
    }
    if (get.umi) {
        data$umi <- as.vector(h5read(sample, "/umi"))
    }
    if (get.gem) {
        data$gem_group <- as.vector(h5read(sample, "/gem_group"))
    }
    if (get.gene || !keep.unmapped) {
        data$gene <- as.vector(h5read(sample, "/gene")) + 1L # Zero-indexed by default
    }
    if (get.reads || length(data)==0) {
        nreads <- as.vector(h5read(sample, "/reads")) # Maybe useful for selective exclusion
        if (get.reads) {
            data$reads <- nreads
        } else {
            data$reads <- matrix(0L, length(nreads), 0) # Just to ensure we get the right number of rows. 
        }
    }

    data <- do.call(DataFrame, data)

    # Defining the set of all genes, removing unassigned gene entries.
    gene.ids <- h5read(sample, "/gene_ids") 
    if (!keep.unmapped) {
        keep <- data$gene <= length(gene.ids)
        if (!get.gene) {
            data$gene <- NULL
        }
        data <- data[keep,]
    }

    # Don't define the total cell pool here, as higher level functions may want to use gem_group.
    return(list(data=data, genes=gene.ids))
}

