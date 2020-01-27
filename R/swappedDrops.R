#' @export
swappedDrops <- function(samples, barcode.length=NULL, ...)
# Removes swapped reads between samples in 10X Genomics data.
#
# written by Jonathan Griffiths
# with modifications from Aaron Lun
# created 18 December 2017
{
    ref.genes <- NULL
    cells <- umis <- genes <- nreads <- vector("list", length(samples))
    names(cells) <- names(samples)

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

    removeSwappedDrops(cells=cells, umis=umis, genes=genes, nreads=nreads, ref.genes=ref.genes, ...)
}

#' @export
removeSwappedDrops <- function(cells, umis, genes, nreads, ref.genes, min.frac=0.8,
    get.swapped=FALSE, get.diagnostics=FALSE, hdf5.out=TRUE)
# Core function to swappedDrops(), split off to accommodate non-HDF5 inputs.
# 
# written by Aaron Lun
# created 17 July 2018
{
    diag.code <- ifelse(get.diagnostics, 1L + as.integer(hdf5.out), 0L)
    swap.out <- find_swapped(cells, genes, umis, nreads, min.frac, diag.code)
    unswapped <- swap.out[[1]]

    nsamples <- length(cells)
    cleaned <- swapped <- vector("list", nsamples)
    names(cleaned) <- names(swapped) <- names(cells)
    for (i in seq_len(nsamples)) { 
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

