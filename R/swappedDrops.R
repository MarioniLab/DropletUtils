#' @export
swappedDrops <- function(samples, barcode.length=NULL, min.frac=0.8, get.swapped=FALSE, get.diagnostics=FALSE, hdf5.out=TRUE)
# Removes swapped reads between samples in 10X Genomics data.
#
# written by Jonathan Griffiths
# with modifications from Aaron Lun
# created 18 December 2017
{
    ref.genes <- NULL
    cells <- umis <- genes <- nreads <- vector("list", length(samples))

    for (i in seq_along(samples)) {
        mol.info <- read10xMolInfo(samples[i], barcode.length=barcode.length)

        if (is.null(ref.genes)) {
            ref.genes <- mol.info$genes            
        } else if (!identical(ref.genes, mol.info$genes)) {
            stop("gene information differs between samples")
        }

        curgems <- mol.info$data$gem_group
        if (length(curgems) && min(curgems)!=max(curgems)) { 
            warning(paste0("sample '", samples[i], "' contains multiple GEM groups"))
        }

        current <- mol.info$data
        cells[[i]] <- current$cell
        umis[[i]] <- current$umi
        genes[[i]] <- current$gene
        nreads[[i]] <- current$reads
    }
        
    # Identifying swapped molecules.
    swap.out <- .Call(cxx_find_swapped_ultra, cells, umis, genes, nreads, min.frac, as.integer(get.diagnostics + hdf5.out))
    unswapped <- swap.out[[1]]

    cleaned <- swapped <- vector("list", length(samples))
    names(cleaned) <- names(swapped) <- names(samples)
    for (i in seq_along(samples)) { 
        cur.cells <- cells[[i]]
        all.cells <- sort(unique(cur.cells))
        cur.genes <- genes[[i]]

        # Forming count matrices from unswapped and swapped molecules. 
        notswap <- unswapped[[i]]
        cleaned[[i]] <- makeCountMatrix(cur.genes[notswap], cur.cells[notswap], all.genes=ref.genes, all.cells=all.cells)
        if (get.swapped) {
            curswap <- !notswap
            swapped[[i]] <- makeCountMatrix(cur.genes[curswap], cur.cells[curswap], all.genes=ref.genes, all.cells=all.cells)
        }
    }

    # Figuring out what kind of output to return.
    output <- list(cleaned=cleaned)
    if (get.swapped){  
        output$swapped <- swapped
    } 
    if (get.diagnostics) { 
        output$diagnostics <- swap.out[[2]]
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
