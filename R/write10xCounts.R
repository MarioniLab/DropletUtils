#' Write count data in the 10x format
#'
#' Create a directory containing the count matrix and cell/gene annotation from a sparse matrix of UMI counts, 
#' in the format produced by the CellRanger software suite.
#' 
#' @param x A sparse numeric matrix of UMI counts.
#' @param path A string containing the path to the output directory (for \code{type="sparse"}) or file (for \code{type="HDF5"}).
#' @param barcodes A character vector of cell barcodes, one per column of \code{x}.
#' @param gene.id A character vector of gene identifiers, one per row of \code{x}.
#' @param gene.symbol A character vector of gene symbols, one per row of \code{x}.
#' @param gene.type A character vector of gene types, expanded to one per row of \code{x}.
#' Only used when \code{version="3"}.
#' @param overwrite A logical scalar specifying whether \code{path} should be overwritten if it already exists.
#' @param type String specifying the type of 10X format to save \code{x} to.
#' This is either a directory containing a sparse matrix with row/column annotation (\code{"sparse"})
#' or a HDF5 file containing the same information (\code{"HDF5"}).
#' @param genome String specifying the genome for storage when \code{type="HDF5"}.
#' This can be a character vector with one genome per feature if \code{version="3"}.
#' @param version String specifying the version of the CellRanger format to produce.
#' @param chemistry,original.gem.groups,library.ids 
#' Strings containing metadata attributes to be added to the HDF5 file for \code{type="HDF5"}.
#' Their interpretation is not formally documented and is left to the user's imagination.
#' 
#' @details
#' This function will try to automatically detect the desired format based on whether \code{path} ends with \code{".h5"}.
#' If so, it assumes that \code{path} specifies a HDF5 file path and sets \code{type="HDF5"}.
#' Otherwise it will set \code{type="sparse"} under the assumption that \code{path} specifies a path to a directory.
#' 
#' Note that there were major changes in the output format for CellRanger version 3.0 to account for non-gene features such as antibody or CRISPR tags. 
#' Users can switch to this new format using \code{version="3"}.
#' See the documentation for \dQuote{latest} for this new format, otherwise see \dQuote{2.2} or earlier.
#'
#' The primary purpose of this function is to create files to use for testing \code{\link{read10xCounts}}.
#' In principle, it is possible to re-use the HDF5 matrices in \code{cellranger reanalyze}.
#' We recommend against doing so routinely due to CellRanger's dependence on undocumented metadata attributes that may change without notice.
#' 
#' @return 
#' For \code{type="sparse"}, a directory is produced at \code{path}.
#' If \code{version="2"}, this will contain the files \code{"matrix.mtx"}, \code{"barcodes.tsv"} and \code{"genes.tsv"}.
#' If \code{version="3"}, it will instead contain \code{"matrix.mtx.gz"}, \code{"barcodes.tsv.gz"} and \code{"features.tsv.gz"}.
#' 
#' For \code{type="HDF5"}, a HDF5 file is produced at \code{path} containing data in column-sparse format.
#' If \code{version="2"}, data are stored in the HDF5 group named \code{genome}.
#' If \code{version="3"}, data are stored in the group \code{"matrix"}.
#' 
#' A \code{TRUE} value is invisibly returned.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{read10xCounts}}, to read CellRanger matrices into R.
#' 
#' @examples
#' # Mocking up some count data.
#' library(Matrix)
#' my.counts <- matrix(rpois(1000, lambda=5), ncol=10, nrow=100)
#' my.counts <- as(my.counts, "dgCMatrix")
#' cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))
#' 
#' ngenes <- nrow(my.counts)
#' gene.ids <- paste0("ENSG0000", seq_len(ngenes))
#' gene.symb <- paste0("GENE", seq_len(ngenes))
#' 
#' # Writing this to file:
#' tmpdir <- tempfile()
#' write10xCounts(tmpdir, my.counts, gene.id=gene.ids, 
#'     gene.symbol=gene.symb, barcodes=cell.ids)
#' list.files(tmpdir)
#'
#' # Creating a version 3 HDF5 file:
#' tmph5 <- tempfile(fileext=".h5")
#' write10xCounts(tmph5, my.counts, gene.id=gene.ids, 
#'     gene.symbol=gene.symb, barcodes=cell.ids, version='3')
#' 
#' @references
#' 10X Genomics (2017).
#' Gene-Barcode Matrices.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/output/matrices}
#' 
#' 10X Genomics (2018).
#' Feature-Barcode Matrices.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices}
#' 
#' 10X Genomics (2018).
#' HDF5 Gene-Barcode Matrix Format.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/advanced/h5_matrices}
#' 
#' 10X Genomics (2018).
#' HDF5 Feature Barcode Matrix Format.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices}
#' 
#' @export
write10xCounts <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, gene.type="Gene Expression",
    overwrite=FALSE, type=c("auto", "sparse", "HDF5"), genome="unknown", version=c("2", "3"),
    chemistry="Single Cell 3' v3", original.gem.groups=1L, library.ids="custom")
{
    # Doing all the work on a temporary location next to 'path', as we have permissions there.
    # This avoids problems with 'path' already existing.
    temp.path <- tempfile(tmpdir=dirname(path)) 
    on.exit({ 
        if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
    })

    # Checking the values.
    if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
        stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
    }
    if (ncol(x)!=length(barcodes)) { 
        stop("'barcodes' must of of the same length as 'ncol(x)'")
    }

    # Determining what format to save in.
    version <- match.arg(version)
    type <- .type_chooser(path, match.arg(type))
    if (type=="sparse") {
        .write_sparse(temp.path, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
    } else {
        .write_hdf5(temp.path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
    }

    # We don't put this at the top as the write functions might fail; 
    # in which case, we would have deleted the existing 'path' for nothing.
    if (overwrite) {
        unlink(path, recursive=TRUE)
    } else if (file.exists(path)) { 
        stop("specified 'path' already exists")
    }
    file.rename(temp.path, path)
    return(invisible(TRUE))
}

#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, gene.type, version="2") {
    dir.create(path, showWarnings=FALSE)
    gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

    if (version=="3") {
        gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
        mhandle <- file.path(path, "matrix.mtx")
        bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
        fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
        on.exit({
            close(bhandle)
            close(fhandle)
        })
    } else {
        mhandle <- file.path(path, "matrix.mtx")
        bhandle <- file.path(path, "barcodes.tsv")
        fhandle <- file.path(path, "genes.tsv")
    }

    writeMM(x, file=mhandle)
    write(barcodes, file=bhandle)
    write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    if (version=="3") {
        # Annoyingly, writeMM doesn't take connection objects.
        gzip(mhandle)
    }

    return(NULL)
}

#' @importFrom rhdf5 h5createFile h5createGroup h5write h5writeAttribute H5Gopen H5Fopen H5Gclose H5Fclose
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.write_hdf5 <- function(path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version="3",
    chemistry="Single Cell 3' v3", original.gem.groups=1L, library.ids="custom")
{
    path <- path.expand(path) # protect against tilde's.
    h5createFile(path)

    if (version=="3") {
        group <- "matrix"
    } else {
        group <- genome
    }
    h5createGroup(path, group)

    h5write(barcodes, file=path, name=paste0(group, "/barcodes"))

    # Saving feature information.
    if (version=="3") {
        h5createGroup(path, file.path(group, "features"))

        h5write(gene.id, file=path, name=paste0(group, "/features/id"))
        h5write(gene.symbol, file=path, name=paste0(group, "/features/name"))
        h5write(rep(gene.type, length.out=length(gene.id)),
            file=path, name=paste0(group, "/features/feature_type"))

        h5write("genome", file=path, name=paste0(group, "/features/_all_tag_keys"))
        h5write(rep(genome, length.out=length(gene.id)),
            file=path, name=paste0(group, "/features/genome"))

        # Writing attributes.
        h5f <- H5Fopen(path)
        h5g <- H5Gopen(h5f, "/")
        h5writeAttribute(chemistry, h5obj=h5g, name="chemistry_description")
        h5writeAttribute("matrix", h5obj=h5g, name="filetype")
        h5writeAttribute(library.ids, h5obj=h5g, name="library_ids")
        h5writeAttribute(original.gem.groups, h5obj=h5g, name="original_gem_groups")
        h5writeAttribute(as.integer(version) - 1L, h5obj=h5g, name="version") # this is probably correct.
        H5Gclose(h5g)
        H5Fclose(h5f)

    } else {
        h5write(gene.id, file=path, name=paste0(group, "/genes"))
        h5write(gene.symbol, file=path, name=paste0(group, "/gene_names"))
    }

    # Saving matrix information.
    x <- as(x, "dgCMatrix")
    h5write(x@x, file=path, name=paste0(group, "/data"))
    h5write(dim(x), file=path, name=paste0(group, "/shape"))
    h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
    h5write(x@p, file=path, name=paste0(group, "/indptr"))

    return(NULL)
}

.type_chooser <- function(path, type) {
    if (type=="auto") {
        type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
    }
    type
}
