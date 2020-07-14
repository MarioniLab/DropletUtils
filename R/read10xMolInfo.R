#' @export
#' @importFrom rhdf5 h5read H5Fopen H5Fclose H5Dopen H5Dclose
#' H5Dget_space H5Sget_simple_extent_dims H5Sclose
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
read10xMolInfo <- function(sample, barcode.length=NULL, keep.unmapped=FALSE, 
    get.cell=TRUE, get.umi=TRUE, get.gem=TRUE, get.gene=TRUE, get.reads=TRUE,
    get.library=TRUE, extract.library.info=FALSE, version=c("auto", "2", "3"))
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
            data$cell <- get_cell_barcodes(sample, "barcode", barcode.length)
        }
    }

    if (get.umi) {
        data$umi <- as.vector(h5read(sample, "/umi"))
    }

    if (get.gem) {
        data$gem_group <- as.vector(h5read(sample, "/gem_group"))
    }

    if (get.gene || !keep.unmapped) {
        # Both of these are zero-indexed by default.
        path <- if (version=="3") "/feature_idx" else "/gene"
        data$gene <- as.vector(h5read(sample, path)) + 1L 
    }

    if (get.reads) {
        path <- if (version=="3") "/count" else "/reads"
        data$reads <- as.vector(h5read(sample, path))
    }

    if (version=="3" && get.library) {
        data$library <- as.vector(h5read(sample, "/library_idx")) + 1L
    }

    if (length(data)==0) {
        # Just to ensure we get the right number of rows,
        # if there were no other fields requested.
        fhandle <- H5Fopen(sample)
        on.exit(H5Fclose(fhandle))
        dhandle <- H5Dopen(fhandle, "/umi")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        space <- H5Dget_space(dhandle)
        on.exit(H5Sclose(space), add=TRUE, after=FALSE)

        N <- H5Sget_simple_extent_dims(space)
        stopifnot(N$rank==1L)
        data <- make_zero_col_DFrame(N$size)
    } else {
        data <- do.call(DataFrame, data)
    }

    # Defining the set of all genes, removing unassigned gene entries.
    path <- if (version=="3") "/features/id" else "/gene_ids"
    gene.ids <- as.vector(h5read(sample, path))

    if (!keep.unmapped) {
        keep <- data$gene <= length(gene.ids)
        if (!get.gene) {
            data$gene <- NULL
        }
        data <- data[keep,,drop=FALSE]
    }

    # Don't define the total cell pool here, as higher level functions may want to use gem_group.
    output <- list(data=data, genes=gene.ids)

    if (version=='3' && extract.library.info) {
        lib.info <- h5read(sample, "/library_info")
        output$library.info <- jsonlite::fromJSON(lib.info, simplifyVector=FALSE)
    }

    output
}

