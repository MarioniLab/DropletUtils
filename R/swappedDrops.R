#' @export
#' @importFrom Matrix sparseMatrix
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
    cells <- umis <- genes <- nreads <- vector("list", length(samples))
    for (i in seq_along(tabs)) {
        current <- tabs[[i]]$data
        cells[[i]] <- current$cell
        umis[[i]] <- current$umi
        genes[[i]] <- current$gene
        nreads[[i]] <- current$reads
    }

    cells <- unlist(cells)
    umis <- unlist(umis)
    genes <- unlist(genes)
    nreads <- unlist(nreads)
    swap.out <- .findSwapped(cells, umis, genes, nreads, min.frac, get.group=get.diagnostics)
    is.swap <- swap.out$swapped

    # Printing out per-molecule read counts. 
    if (get.diagnostics) {
        sample.ids <- rep(seq_along(tabs), unlist(lapply(tabs, FUN=function(x) { nrow(x$data) })))
        nmolecules <- if (length(swap.out$group)) max(swap.out$group) else 0
        diagnostics <- sparseMatrix(i=swap.out$group, j=sample.ids, x=nreads, 
                                    dims=c(nmolecules, length(samples)))
        colnames(diagnostics) <- names(samples)
    }

    # Iterating through the samples.
    cleaned <- swapped <- vector("list", length(samples))
    names(cleaned) <- names(swapped) <- names(samples)
    last <- 0L
    for (i in seq_along(tabs)) { 
        current <- tabs[[i]]$data
        all.genes <- tabs[[i]]$genes
        all.cells <- sort(unique(current$cell))

        # Forming count matrices from unswapped and swapped molecules. 
        curswap <- is.swap[last + seq_along(current[[1]])]
        cleaned[[i]] <- makeCountMatrix(current$gene[!curswap], current$cell[!curswap], 
                                        all.genes=all.genes, all.cells=all.cells)
        if (get.swapped) {
            swapped[[i]] <- makeCountMatrix(current$gene[curswap], current$cell[curswap], 
                                            all.genes=all.genes, all.cells=all.cells)
        }

        last <- last + length(current[[1]])
    }

    # Figuring out what kind of output to return.
    output <- list(cleaned=cleaned)
    if (get.swapped){  
        output$swapped <- swapped
    } 
    if (get.diagnostics) { 
        output$diagnostics <- diagnostics
    }
    return(output)
}

.findSwapped <- function(cells, umis, genes, reads, min.frac=0.8, get.group=FALSE)  
# Identifies the molecules that are swapped or not. Technically we could 
# use sparse matrices and max.col() to do this, but max.col() isn't supported
# natively for sparse matrices and would call as.matrix() instead.
{
    o <- order(cells, umis, genes)
    cells <- cells[o]
    umis <- umis[o]
    genes <- genes[o]
    reads <- reads[o]

    # Figures out the runs of the same type of cell/umi/gene combination.
    N <- length(o)
    if (N) { 
        is.diff <- cells[-N]!=cells[-1] | umis[-N]!=umis[-1] | genes[-N]!=genes[-1]
        runs <- diff(c(0L, which(is.diff), N))
    } else {
        runs <- integer(0)
    }

    # Identifying putative swapped reads.
    is.swap <- .Call(cxx_find_swapped, runs, reads, min.frac)
    is.swap[o] <- is.swap
    output <- list(swapped=is.swap) 

    if (get.group) {
        grouping <- rep(seq_along(runs), runs)
        grouping[o] <- grouping
        output$group <- grouping
    }
    return(output)
}
