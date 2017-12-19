swappedDrops <- function(samples, barcode.length, get.swapped=FALSE, get.diagnostics=FALSE, min.frac=0.8)
# This function removes swapped reads between samples in 10X Genomics data.
#
# written by Jonathan Griffiths
# with modifications from Aaron Lun
# created 18 December 2017
{
    tabs <- lapply(samples, .readHDF5Data, barcode_len=barcode.length)

    # Warning if multiple GEM groups are observed, as this interferes with deswapping.
    for (i in seq_along(tabs)) {
        curgems <- tabs[[i]]$data$gem_group
        if (min(curgems)!=max(curgems)) { 
            warning("each sample should not contain multiple 10X runs")
            break
        }
        tabs[[i]]$data$gem_group <- NULL
    }
        
    # Identifying swapped molecules.
    swap.marks <- nreads <- vector("list", length(samples))
    for (i in seq_along(tabs)) {
        current <- tabs[[i]]$data
        swap.marks[[i]] <- paste(current$cell, current$umi, current$gene)
        nreads[[i]] <- current$reads
    }
    
    swap.marks <- unlist(swap.marks)
    nreads <- unlist(nreads)
    swap <- .findSwapped(swap.marks, nreads, min.frac)

    # Printing out per-molecule read counts. 
    if (get.diagnostics) {
        sample.ids <- rep(seq_along(tabs), unlist(lapply(tabs, FUN=function(x) { length(x$data[[1]]) })))
        all.molecules <- sort(unique(swap.marks))
        m <- match(swap.marks, all.molecules)
        diagnostics <- sparseMatrix(i=m, j=sample.ids, x=nreads, 
                                    dims=c(length(all.molecules), length(samples)))
        rownames(diagnostics) <- all.molecules
        colnames(diagnostics) <- names(samples)
    }

    # Producing an output sparseMatrix.
    cleaned <- swapped <- vector("list", length(samples))
    names(cleaned) <- names(swapped) <- names(samples)
    last <- 0L
    for (i in seq_along(tabs)) { 
        current <- tabs[[i]]$data
        anno <- tabs[[i]]$annotation
        curswap <- swap[last + seq_along(current[[1]])]

        mat <- sparseMatrix(i=current$gene[!curswap], 
                            j=match(current$cell[!curswap], anno$cells),
                            x=rep(1, sum(!curswap)),
                            dims=c(length(anno$genes), length(anno$cells)),
                            use.last.ij=FALSE, # Adds up the duplicates.
                            giveCsparse=TRUE)
        rownames(mat) <- anno$genes
        colnames(mat) <- anno$cells
        cleaned[[i]] <- mat

        # Adding up the reads to discard if requested.
        if (get.swapped) {
            mat <- sparseMatrix(i=current$gene[curswap], 
                                j=match(current$cell[curswap], anno$cells),
                                x=rep(1, sum(curswap)),
                                dims=c(length(anno$genes), length(anno$cells)),
                                use.last.ij=FALSE, # Adds up the duplicates.
                                giveCsparse=TRUE)
            rownames(mat) <- anno$genes
            colnames(mat) <- anno$cells
            swapped[[i]] <- mat
        }

        last <- last + length(current[[1]])
    }

    # Figuring out what kind of output to return.
    if (!get.swapped && !get.diagnostics) {
        return(cleaned)
    } else {
        output <- list(cleaned=cleaned)
        if (get.swapped){  
            output$swapped <- swapped
        } 
        if (get.diagnostics) { 
            output$diagnostics <- diagnostics
        }
        return(output)
    }
}

.readHDF5Data <- function(h5_loc, barcode_len) {
    cell <- .Call(cxx_get_cell_barcodes, h5_loc, "barcode", barcode_len)
    umi <- h5read(h5_loc, "/umi")
    gem_group <- h5read(h5_loc, "/gem_group")
    gene <- h5read(h5_loc, "/gene") + 1L #zero-indexed by default
    reads <- h5read(h5_loc, "/reads") #maybe useful for selective exclusion
  
    # Defining the set of all barcodes, and of all genes.
    all.barcodes <- sort(unique(cell))
    gene.ids <- h5read(h5_loc, "/gene_ids") 

    # Remove the unassigned gene entries.
    keep <- gene <= length(gene.ids)
    cell <- cell[keep]
    umi <- umi[keep]
    gem_group <- gem_group[keep]
    gene <- gene[keep]
    reads <- reads[keep]

    return(list(data=list(cell=cell, umi=umi, gene=gene, reads=reads, gem_group=gem_group), 
                annotation=list(genes=gene.ids, cells=all.barcodes)))
}

.findSwapped <- function(swap.marks, reads, min.frac=0.8, get.diagnostics=FALSE) { 
    # Identifies the molecules that are swapped or not. Technically we could 
    # use sparse matrices and max.col() to do this, but max.col() isn't supported
    # natively for sparse matrices and would call as.matrix() instead.
    o <- order(swap.marks)
    rout <- rle(swap.marks[o])
    is.swap <- .Call(cxx_find_swapped, rout$lengths, reads[o], min.frac)
    is.swap[o] <- is.swap
    return(is.swap)
}
