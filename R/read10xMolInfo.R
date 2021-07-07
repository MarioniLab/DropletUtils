#' Read the 10X molecule information file
#' 
#' Extract relevant fields from the molecule information HDF5 file, produced by CellRanger for 10X Genomics data.
#' 
#' @param sample A string containing the path to the molecule information HDF5 file.
#' @param barcode.length An integer scalar specifying the length of the cell barcode.
#' Only relevant when \code{version="2"}.
#' @param keep.unmapped A logical scalar indicating whether unmapped molecules should be reported.
#' @param get.cell,get.umi,get.gem,get.gene,get.reads,get.library 
#' Logical scalar indicating whether the corresponding field should be extracted for each molecule.
#' @param version String specifying the version of the 10X molecule information format to read data from.
#' @param extract.library.info Logical scalar indicating whether the library information should be extracted.
#' Only relevant when \code{version="3"}.
#' 
#' @return
#' A named list is returned containing \code{data}, 
#' a \linkS4class{DataFrame} where each row corresponds to a single transcript molecule.
#' This contains the following fields:
#' \describe{
#' \item{\code{barcode}:}{Character, the cell barcode for each molecule.}
#' \item{\code{umi}:}{Integer, the processed UMI barcode in 2-bit encoding.} 
#' \item{\code{gem_group}:}{Integer, the GEM group.}
#' \item{\code{gene}:}{Integer, the index of the gene to which the molecule was assigned.
#' This refers to an entry in the \code{genes} vector, see below.}
#' \item{\code{reads}:}{Integer, the number of reads mapped to this molecule.}
#' \item{\code{reads}:}{Integer, the number of reads mapped to this molecule.}
#' \item{\code{library}:}{Integer, the library index in cases where multiple libraries are present in the same file.
#' Only reported when \code{version="3"}.}
#' }
#' A field will not be present in the DataFrame if the corresponding \code{get.*} argument is \code{FALSE}, 
#' 
#' The second element of the list is \code{genes}, a character vector containing the names of all genes in the annotation.
#' This is indexed by the \code{gene} field in the \code{data} DataFrame.
#'
#' If \code{version="3"}, a \code{feature.type} entry is added to the list.
#' This is a character vector of the same length as \code{genes}, containing the feature type for each gene.
#'
#' If \code{extract.library.info=TRUE}, an additional element named \code{library.info} is returned.
#' This is a list of lists containing per-library information such as the \code{"library_type"}.
#' The \code{library} field in the \code{data} DataFrame indexes this list.
#' 
#' @details
#' Molecules that were not assigned to any gene have \code{gene} set to \code{length(genes)+1}.
#' By default, these are removed when \code{keep.unmapped=FALSE}.
#' 
#' CellRanger 3.0 introduced a major change in the format of the molecule information files.
#' When \code{version="auto"}, the function will attempt to determine the version format of the file.
#' This can also be user-specified by setting \code{version} explicitly.
#' 
#' For files produced by version 2.2 of the CellRanger software, the length of the cell barcode is not given.
#' Instead, the barcode length is automatically inferred if \code{barcode.length=NULL} and \code{version="2"}.
#' Currently, version 1 of the 10X chemistry uses 14 nt barcodes, while version 2 uses 16 nt barcodes.
#' 
#' Setting any of the \code{get.*} arguments will (generally) avoid extraction of the corresponding field.
#' This can improve efficiency if that field is not necessary for further analysis.
#' Aside from the missing field, the results are guaranteed to be identical, i.e., same order and number of rows.
#' 
#' @author
#' Aaron Lun,
#' based on code by Jonathan Griffiths
#' 
#' @seealso
#' \code{\link{makeCountMatrix}}, which creates a count matrix from this information. 
#' 
#' @examples
#' # Mocking up some 10X HDF5-formatted data.
#' out <- DropletUtils:::simBasicMolInfo(tempfile())
#' 
#' # Reading the resulting file.
#' read10xMolInfo(out)
#' 
#' @references
#' Zheng GX, Terry JM, Belgrader P, and others (2017).
#' Massively parallel digital transcriptional profiling of single cells. 
#' \emph{Nat Commun} 8:14049.
#' 
#' 10X Genomics (2017).
#' Molecule info.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/output/molecule_info}
#' 
#' 10X Genomics (2018).
#' Molecule info.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info}
#'
#' @export
#' @importFrom rhdf5 h5read H5Fopen H5Fclose H5Dopen H5Dclose
#' H5Dget_space H5Sget_simple_extent_dims H5Sclose
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
read10xMolInfo <- function(sample, barcode.length=NULL, keep.unmapped=FALSE, 
    get.cell=TRUE, get.umi=TRUE, get.gem=TRUE, get.gene=TRUE, get.reads=TRUE,
    get.library=TRUE, extract.library.info=FALSE, version=c("auto", "2", "3"))
{
    sample <- path.expand(sample) # protect against unexpanded tilde's.

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

    if (version=='3') {
        output$feature.type <- as.vector(h5read(sample, "/features/feature_type"))

        if (extract.library.info) {
            lib.info <- h5read(sample, "/library_info")
            output$library.info <- jsonlite::fromJSON(lib.info, simplifyVector=FALSE)
        }
    }

    output
}

.extract_mol_info <- function(sample, ..., use.library=NULL, get.library=FALSE, 
    extract.library.info=FALSE, subset.library.features=FALSE) 
{
    original.get <- get.library
    original.extract <- extract.library.info

    if (!is.null(use.library)) {
        get.library <- TRUE
        needs.lib.info <- is.character(use.library) || subset.library.features
        if (needs.lib.info) {
            extract.library.info <- TRUE
        }
    }

    output <- read10xMolInfo(sample, ..., get.library=get.library, extract.library.info=extract.library.info)

    has.library <- !is.null(output$data$library)

    if (!is.null(use.library) && has.library) {
        if (needs.lib.info) {
            available.lib <- vapply(output$library.info, FUN="[[", i="library_type", FUN.VALUE="")
            if (is.character(use.library)) {
                use.library <- which(available.lib %in% use.library) 
            }
        }

        output$data <- output$data[output$data$library %in% use.library,,drop=FALSE]

        if (subset.library.features) {
            all.types <- output$feature.type
            output <- .reindex_mol_info_features(output, all.types %in% available.lib[use.library])
        }
    }

    if (!original.get) {
        output$data$library <- NULL
    }
    if (!original.extract) {
        output$library.info <- NULL
    }
    output
}

.reindex_mol_info_features <- function(mol.info, keep) {
    present <- which(keep)
    mol.info$data$gene <- m <- match(mol.info$data$gene, present)
    mol.info$data <- mol.info$data[!is.na(m),,drop=FALSE]
    mol.info$genes <- mol.info$genes[present]
    mol.info$feature.type <- mol.info$feature.type[present]
    mol.info
}
