#' Remove chimeric molecules 
#'
#' Remove chimeric molecules within each cell barcode's library in a droplet experiment.
#'
#' @param sample String containing paths to the molecule information HDF5 files, 
#' produced by CellRanger for 10X Genomics data.
#' @param barcode.length An integer scalar specifying the length of the cell barcode, see \code{\link{read10xMolInfo}}.
#' @param use.library An integer vector specifying the library indices for which to extract molecules from \code{sample}.
#' Alternatively, a character string specifying one or more library types, e.g., \code{"Gene expression"}.
#' @param ... Further arguments to be passed to \code{removeChimericDrops}.
#' @param cells Character vector containing cell barcodes, where each entry corresponds to one molecule.
#' @param umis Integer vector containing encoded UMI sequences, see \code{?\link{encodeSequences}} for details.
#' @param genes Integer vector specifying the gene indices.
#' Each index should refer to an element of \code{ref.genes}.
#' @param nreads Integer vector containing the number of reads per molecule.
#' @param ref.genes A character vector containing the names or symbols of all genes.
#' @param min.frac A numeric scalar specifying the minimum fraction of reads required for a chimeric molecule to be retained.
#' @param get.chimeric A logical scalar indicating whether the UMI counts corresponding to chimeric molecules should be returned.
#' @param get.diagnostics A logical scalar indicating whether to return statistics for each molecule grouping.
#'
#' @return 
#' A list is returned with the \code{cleaned} entry, 
#' a sparse matrix containing the UMI count for each gene (row) and cell barcode (column) after removing chimeric molecules.
#' All cell barcodes that were originally observed are reported as columns, 
#' though note that it is theoretically possible for some barcodes to contain no counts.
#' 
#' If \code{get.chimeric=TRUE}, a \code{chimeric} entry is returned in the list.
#' This is a sparse matrix of UMI counts corresponding to the chimeric molecules.
#' Adding the cleaned and chimeric matrices should yield the total UMI count prior to removal of swapped molecules. 
#' 
#' If \code{get.diagnostics=TRUE}, the top-level list will also contain an additional \code{diagnostics} \linkS4class{DataFrame}.
#' Each row corresponds to a group of molecules in the same cell with the same UMI.
#' The DataFrame holds the number of molecules in the group,
#' the sum of reads across all molecules in the group,
#' and the proportion of reads assigned to the most sequenced molecule.
#'
#' @details
#' Chimeric molecules are occasionally generated during library preparation for highly multiplexed droplet experiments.
#' Here, incomplete PCR products from one molecule hybridise to another molecule for extension using shared sequences like the poly-A tail for 3' protocols.
#' This produces an amplicon where the UMI and cell barcode originate from one transcript molecule but the gene sequence is from another.
#' If the second template is from another cell, this effect results in contamination of one cell's profile by another, similar to the contamination between samples discussed in \code{?\link{swappedDrops}}.
#'
#' Chimerism manifests as molecules that have the same UMI sequence and cell barcode but are assigned to different genes.
#' To remove them, this function will simply discard all molecules within the same cell that share UMI sequences.
#' Of course, this may also remove non-chimeric molecules that have the same UMI by chance,
#' but for typical UMI lengths (10-12 bp for 10X protocols) we expect UMI collisions to be very rare between molecules from the same cell.
#' 
#' Nonetheless, to mitigate losses due to collisions,
#' we retain any molecule that has a much greater number of reads compared to all other molecules with the same UMI in the same cell.
#' This is based on the expectation that the original non-chimeric molecule will have undergone more rounds of PCR amplification compared to its chimeric offspring, and thus will have higher read coverage.
#' For all molecules with the same UMI within a given cell, we compute the proportion of reads assigned to each molecule and we keep the molecule with a proportion above \code{min.frac}.
#' If no molecule passes this threshold, the entire set is discarded.
#'
#' The \code{use.library} argument can be used to only check for chimeras within a given feature type, e.g., CRISPR tags.
#' This is most relevant in situations where \code{sample} contains multiple libraries that involve different sets of shared sequences,
#' such that chimeras are unlikely to form between molecules from different libraries.
#' Analysis of just one library can be achieved by setting \code{use.library} to the name or index of the desired feature set.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Mocking up some 10x HDF5-formatted data.
#' curfile <- DropletUtils:::simBasicMolInfo(tempfile())
#' 
#' out <- chimericDrops(curfile)
#' dim(out$cleaned)
#' 
#' out2 <- chimericDrops(curfile, get.diagnostics=TRUE)
#' out2$diagnostics
#'
#' @references
#' Dixit A. (2016).
#' Correcting chimeric crosstalk in single cell RNA-seq experiments.
#' \emph{biorXiv}, \url{https://doi.org/10.1101/093237}
#'
#' @export
chimericDrops <- function(sample, barcode.length=NULL, use.library=NULL, ...) {
    mol.info <- .extract_mol_info(sample, barcode.length=barcode.length, 
        use.library=use.library, subset.library.features=TRUE)
    df <- mol.info$data
    removeChimericDrops(df$cell, df$umi, df$gene, df$reads, ref.genes=mol.info$genes, ...)
}

#' @export
#' @rdname chimericDrops
#' @importFrom S4Vectors DataFrame
removeChimericDrops <- function(cells, umis, genes, nreads, ref.genes, min.frac=0.8,
    get.chimeric=FALSE, get.diagnostics=FALSE)
{
    chim.out <- find_chimeric(cells, umis, nreads, min.frac, get.diagnostics)
    notswap <- chim.out[[1]]

    all.cells <- sort(unique(cells))
    output <- list(
        cleaned=makeCountMatrix(genes[notswap], cells[notswap], all.genes=ref.genes, all.cells=all.cells)
    )

    if (get.chimeric) {
        curswap <- !notswap
        output$chimeric <- makeCountMatrix(genes[curswap], cells[curswap], all.genes=ref.genes, all.cells=all.cells)
    }

    if (get.diagnostics) { 
        names(chim.out[[2]]) <- c("number", "reads", "prop")
        output$diagnostics <- DataFrame(chim.out[[2]])
    }

    output
}

