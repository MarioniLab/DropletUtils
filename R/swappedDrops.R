swappedDrops <- function(samples, barcode.length, output = c("cleaned", "total", "swapped"), use.fracs = FALSE, min.frac = 0.8)
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
        
    # Identifying swapped molecules, if requested.
    output <- match.arg(output)
    if (output=="total") { 
        dup <- logical(sum(lengths(nreads)))
    } else {
        dup.marks <- nreads <- vector("list", length(samples))
        for (i in seq_along(tabs)) {
            current <- tabs[[i]]$data
            dup.marks[[i]] <- paste(current$cell, current$umi, current$gene)
            nreads[[i]] <- current$reads
        }

        dup.marks <- unlist(dup.marks)
        nreads <- unlist(nreads)
        dup <- .findDuplicated(dup.marks, nreads, min.frac)
    }

    # Producing an output sparseMatrix.
    counts <- vector("list", length(samples))
    last <- 0L
    for (i in seq_along(tabs)) { 
        current <- tabs[[i]]$data
        anno <- tabs[[i]]$annotation
        curdup <- dup[last + seq_along(current[[1]])]

        mat <- sparseMatrix(i=current$gene[!curdup], 
                            j=match(current$cell[!curdup], anno$cells),
                            x=rep(1, sum(!curdup)),
                            dims=c(length(anno$genes), length(anno$cells)),
                            use.last.ij=FALSE, # Adds up the duplicates.
                            giveCsparse=TRUE)
        rownames(mat) <- anno$genes
        colnames(mat) <- anno$cells                

        counts[[i]] <- mat
        last <- last + length(current[[1]])
    }
    return(counts)
}

.readHDF5Data <- function(h5_loc, barcode_len) {
    cell <- .Call(cxx_get_cell_barcodes, h5_loc, "barcode", barcode_len)
    umi <- h5read(h5_loc, "/umi")
    gem_group <- h5read(h5_loc, "/gem_group")
    gene <- h5read(h5_loc, "/gene") + 1L #zero-indexed by default
    reads <- h5read(h5_loc, "/reads") #maybe useful for selective exclusion
  
    # Defining the set of all barcodes.
    all.barcodes <- unique(cell)
    gene.ids <- h5read(h5_loc, "/gene_ids") 

    # Remove the unassigned gene entries.
    keep <- gene <= length(gene.ids)
    cell <- cell[keep]
    umi <- umi[keep]
    gene <- gene[keep]
    reads <- reads[keep]

    return(list(data=list(cell=cell, umi=umi, gene=gene, reads=reads, gem_group=gem_group), 
                annotation=list(genes=gene.ids, cells=all.barcodes)))
}

.findDuplicated <- function(dup.marks, reads, min.frac=0.8) { 
    o <- order(dup.marks)
    rout <- rle(dup.marks[o])
    is.dup <- .Call(cxx_find_swapped, rout$lengths, reads[o], min.frac)
    is.dup[o] <- is.dup
    return(is.dup)
}
