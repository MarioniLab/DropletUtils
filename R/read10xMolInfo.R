#' @export
#' @importFrom rhdf5 h5read
#' @importFrom S4Vectors DataFrame
read10xMolInfo <- function(sample, barcode.length=NULL, keep.unmapped=FALSE, 
    get.cell=TRUE, get.umi=TRUE, get.gem=TRUE, get.gene=TRUE, get.reads=TRUE, 
    version=c("auto", "2", "3"))
# Utility function to read useful information from a 10X molecule information file.
#
# written by Aaron Lun
# based on code from Jonathan Griffiths
# created 20 December 2017    
{
    version <- match.arg(version)
    if (version=="auto") {
        available <- h5ls(sample, recursive=FALSE)
        version <- if ("barcode_idx" %in% available$name) "3" else "2"
    }

    data <- list()

    if (get.cell) {
        if (version=="3") {
            all.barcodes <- as.vector(h5read(sample, "/barcodes"))
            all.barcodes <- sub("-[0-9]+", "", all.barcodes) # removing GEM group.
            data$cell <- all.barcodes[as.vector(h5read(sample, "/barcode_idx")) + 1L]
        } else {
            data$cell <- .Call(cxx_get_cell_barcodes, sample, "barcode", barcode.length)
        }
    }

    if (get.umi) {
        data$umi <- as.vector(h5read(sample, "/umi"))
    }

    if (get.gem) {
        data$gem_group <- as.vector(h5read(sample, "/gem_group"))
    }

    if (get.gene || !keep.unmapped) {
        # Both of these are zero-indexed by default
        if (version=="3") {
            data$gene <- as.vector(h5read(sample, "/feature_idx")) + 1L 
        } else {
            data$gene <- as.vector(h5read(sample, "/gene")) + 1L
        }
    }

    if (get.reads || length(data)==0) {
        if (version=="3") {
            nreads <- as.vector(h5read(sample, "/count"))
        } else {
            nreads <- as.vector(h5read(sample, "/reads")) 
        }

        if (get.reads) {
            data$reads <- nreads
        } else {
            # Just to ensure we get the right number of rows,
            # if there were no other fields requested.
            data$reads <- matrix(0L, length(nreads), 0) 
        }
    }

    data <- do.call(DataFrame, data)

    # Defining the set of all genes, removing unassigned gene entries.
    if (version=="3") {
        gene.ids <- h5read(sample, "/features/id") 
    } else {
        gene.ids <- h5read(sample, "/gene_ids") 
    }
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

