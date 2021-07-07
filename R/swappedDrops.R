#' Clean barcode-swapped droplet data
#' 
#' Remove the effects of barcode swapping on droplet-based single-cell RNA-seq data, specifically 10X Genomics datasets.
#' 
#' @param samples A character vector containing paths to the molecule information HDF5 files, 
#' produced by CellRanger for 10X Genomics data.
#' Each file corresponds to one sample in a multiplexed pool.
#' @param barcode.length An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.
#' @param use.library An integer scalar specifying the library index for which to extract molecules from \code{sample}.
#' Alternatively, a string specifying the library type, e.g., \code{"Gene expression"}.
#' @param ... Further arguments to be passed to \code{removeSwappedDrops}.
#' @param cells A list of character vectors containing cell barcodes.
#' Each vector corresponds to one sample in a multiplexed pool, and each entry of the vector corresponds to one molecule.
#' @param umis A list of integer vectors containing encoded UMI sequences, organized as described for \code{cells}.
#' See \code{?\link{encodeSequences}} to convert sequences to integers.
#' @param genes A list of integer vectors specifying the gene indices, organized as described for \code{cells}.
#' Each index should refer to an element of \code{ref.genes}.
#' @param nreads A list of integer vectors containing the number of reads per molecule, organized as described for \code{cells}.
#' @param ref.genes A character vector containing the names or symbols of all genes.
#' @param min.frac A numeric scalar specifying the minimum fraction of reads required for a swapped molecule to be assigned to a sample.
#' @param get.swapped A logical scalar indicating whether the UMI counts corresponding to swapped molecules should be returned.
#' @param get.diagnostics A logical scalar indicating whether to return the number of reads for each molecule in each sample.
#' @param hdf5.out Deprecated and ignored.
#' 
#' @details
#' Barcode swapping on the Illumina sequencer occurs when multiplexed samples undergo PCR re-amplification on the flow cell by excess primer with different barcodes.
#' This results in sequencing of the wrong sample barcode and molecules being assigned to incorrect samples after debarcoding.
#' With droplet data, there is the opportunity to remove such effects based on the combination of gene, UMI and cell barcode for each observed transcript molecule.
#' It is very unlikely that the same combination will arise from different molecules in multiple samples.
#' Thus, observation of the same combination across multiple samples is indicative of barcode swapping.
#' 
#' We can remove swapped molecules based on the number of reads assigned to each gene-UMI-barcode combination.
#' From the total number of reads assigned to that combination, the fraction of reads in each sample is calculated.
#' The sample with the largest fraction that is greater than \code{min.frac} is defined as the putative sample of origin to which the molecule is assigned.
#' This assumes that the swapping rate is low, so the sample of origin for a molecule should contain the majority of the reads.
#' In other all samples, reads for the combination are assumed to derive from swapping and do not contribute to the count matrix.
#' Setting \code{min.frac=1} will effectively remove all molecules that appear in multiple samples.
#' We do not recommend setting \code{min.frac} lower than 0.5.
#' 
#' If \code{diagnostics=TRUE}, a diagnostics matrix is returned containing the number of reads per gene-UMI-barcode combination in each sample.
#' Each row corresponds to a combination and each column corresponds to a sample.
#' This can be useful for examining the level of swapping across samples on a molecule-by-molecule basis, 
#' though for the sake of memory, the actual identity of the molecules is not returned.
#' By default, the matrix is returned as a \linkS4class{HDF5Matrix}, which reduces memory usage and avoids potential issues with integer overflow.
#' If \code{hdf5.out=FALSE}, a sparse matrix is returned instead, which is faster but uses more memory.
#' 
#' \code{swappedDrops} is a wrapper around \code{removeSwappedDrops} that extracts the relevant data from the 10X Genomics molecule information file.
#' For other types of droplet-based data, it may be more convenient to call \code{removeSwappedDrops} directly.
#' 
#' @return 
#' A list is returned with the \code{cleaned} entry, itself a list of sparse matrices.
#' Each matrix corresponds to a sample and contains the UMI count for each gene (row) and cell barcode (column) after removing swapped molecules.
#' All cell barcodes that were originally observed are reported as columns, though note that it is possible for some barcodes to contain no counts.
#' 
#' If \code{get.swapped=TRUE}, a \code{swapped} entry is returned in the top-level list.
#' This is a list containing sample-specific sparse matrices of UMI counts corresponding to the swapped molecules.
#' Adding the cleaned and swapped matrices for each sample should yield the total UMI count prior to removal of swapped molecules. 
#' 
#' If \code{get.diagnostics=TRUE}, the top-level list will also contain an additional \code{diagnostics} matrix.
#' 
#' @section Format of the molecule information file:
#' \code{swappedDrops} makes a few assumptions about the nature of the data in each molecule information file.
#' These are necessary to simplify downstream processing and are generally acceptable in most cases.
#' 
#' Each molecule information file should contain data from only a single 10X run.
#' Users should \emph{not} combine multiple samples into a single molecule information file.
#' The function will emit a warning upon detecting multiple GEM groups from any molecule information file.
#' Molecules with different GEM groups will not be recognised as coming from a different sample, though they will be recognised as being derived from different cell-level libraries.
#' 
#' In files produced by CellRanger version 3.0, an additional per-molecule field is present indicating the (c)DNA library from which the molecule was derived.
#' Library preparation can be performed separately for different features (e.g., antibodies, CRISPR tags) such that one 10X run can contain data from multiple libraries.
#' This allows for arbitrarily complicated multiplexing schemes - for example, gene expression libraries might be multiplexed together across one set of samples,
#' while the antibody-derived libraries might be multiplexed across another \emph{different} set of samples.
#' For simplicity, we assume that multiplexing was performed across the same set of \code{samples} for all libraries therein.
#'
#' If a different multiplexing scheme was applied for each library type, users can set \code{use.library} to only check for swapping within a given library type(s).
#' For example, if the multiplexed set of samples for the gene expression libraries is different from the multiplexed set for the CRISPR libraries,
#' one could run \code{swappedDrops} separately on each set of samples with \code{use.library} set to the corresponding type. 
#' This avoids having to take the union of both sets of samples for a single \code{swappedDrops} run,
#' which could detect spurious swaps between samples that were never multiplexed together for the same library type.
#'
#' @author
#' Jonathan Griffiths,
#' with modifications by Aaron Lun
#' 
#' @seealso
#' \code{\link{read10xMolInfo}},
#' \code{\link{encodeSequences}}
#' 
#' @examples
#' # Mocking up some 10x HDF5-formatted data, with swapping.
#' curfiles <- DropletUtils:::simSwappedMolInfo(tempfile(), nsamples=3)
#' 
#' # Obtaining count matrices with swapping removed.
#' out <- swappedDrops(curfiles)
#' lapply(out$cleaned, dim)
#' 
#' out <- swappedDrops(curfiles, get.swapped=TRUE, get.diagnostics=TRUE)
#' names(out)
#' 
#' @references
#' Griffiths JA, Lun ATL, Richard AC, Bach K, Marioni JC (2018).
#' Detection and removal of barcode swapping in single-cell RNA-seq data.
#' \emph{Nat. Commun.} 9, 1:2667.
#'
#' @export
swappedDrops <- function(samples, barcode.length=NULL, use.library=NULL, ...) {
    ref.genes <- NULL
    cells <- umis <- genes <- nreads <- vector("list", length(samples))
    names(cells) <- names(samples)

    for (i in seq_along(samples)) {
        mol.info <- .extract_mol_info(samples[i], barcode.length=barcode.length, 
            use.library=use.library, subset.library.features=TRUE)

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
#' @rdname swappedDrops
removeSwappedDrops <- function(cells, umis, genes, nreads, ref.genes, min.frac=0.8, get.swapped=FALSE, 
    get.diagnostics=FALSE, hdf5.out=FALSE) 
# Core function to swappedDrops(), split off to accommodate non-HDF5 inputs.
# 
# written by Aaron Lun
# created 17 July 2018
{
    swap.out <- find_swapped(cells, genes, umis, nreads, min.frac, get.diagnostics)
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
        colnames(output$diagnostics) <- names(cleaned)
    }

    output
}

