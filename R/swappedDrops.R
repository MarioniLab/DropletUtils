swappedDrops <- function(samples, barcode.length=NULL, get.swapped=FALSE, get.diagnostics=FALSE, min.frac=0.8)
# This function removes swapped reads between samples in 10X Genomics data.
#
# written by Jonathan Griffiths
# with modifications from Aaron Lun
# created 18 December 2017
{
    tabs <- lapply(samples, read10xMolInfo, barcode.length=barcode.length)
    for (i in seq_along(tabs)) {
        if (!identical(tabs[[1]]$genes, tabs[[i]]$genes)) {
            stop("gene information differs between samples")
        }

        curgems <- tabs[[i]]$data$gem_group
        if (length(curgems) && min(curgems)!=max(curgems)) { 
            warning("each sample should only contain GEM codes from one 10X run")
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
        all.genes <- tabs[[i]]$genes
        all.cells <- sort(unique(current$cell))
        curswap <- swap[last + seq_along(current[[1]])]

        mat <- sparseMatrix(i=current$gene[!curswap], 
                            j=match(current$cell[!curswap], all.cells),
                            x=rep(1, sum(!curswap)),
                            dims=c(length(all.genes), length(all.cells)),
                            use.last.ij=FALSE, # Adds up the duplicates.
                            giveCsparse=TRUE)
        rownames(mat) <- all.genes
        colnames(mat) <- all.cells
        cleaned[[i]] <- mat

        # Adding up the reads to discard if requested.
        if (get.swapped) {
            mat <- sparseMatrix(i=current$gene[curswap], 
                                j=match(current$cell[curswap], all.cells),
                                x=rep(1, sum(curswap)),
                                dims=c(length(all.genes), length(all.cells)),
                                use.last.ij=FALSE, # Adds up the duplicates.
                                giveCsparse=TRUE)
            rownames(mat) <- all.genes
            colnames(mat) <- all.cells
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
