#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
read10xCounts <- function(samples, col.names=FALSE, type=c("auto", "sparse", "HDF5"), group=NULL)
# Reads in one or more 10X directories in 'samples', and produces
# a SingleCellExperiment object as the output.
#
# written by Davis McCarthy
# modifications by Aaron Lun
# some time ago.
{
    nsets <- length(samples)
    full_data <- vector("list", nsets)
    gene_info_list <- vector("list", nsets)
    cell_info_list <- vector("list", nsets)
    type <- match.arg(type)

    for (i in seq_len(nsets)) {
        run <- samples[i]
        type <- .type_chooser(run, type)

        if (type=="sparse") {
            info <- .read_from_sparse(run)
        } else {
            info <- .read_from_hdf5(run, group=group)
        }

        full_data[[i]] <- info$mat
        gene_info_list[[i]] <- info$gene.info
        cell.names <- info$cell.names
        cell_info_list[[i]] <- DataFrame(
          Sample = Rle(run, length(cell.names)),
          Barcode = cell.names, row.names=NULL)
    }

    # Checking gene uniqueness.
    if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
        stop("gene information differs between runs")
    }
    gene_info <- gene_info_list[[1]]
    colnames(gene_info) <- c("ID", "Symbol")
    rownames(gene_info) <- gene_info$ID

    # Forming the full data matrix.
    full_data <- do.call(cbind, full_data)

    # Adding the cell data (only using as colnames if there is only 1 set - guaranteed unique).
    cell_info <- do.call(rbind, cell_info_list)
    if (col.names && nsets == 1L) {
        colnames(full_data) <- cell_info$Barcode
    }

    SingleCellExperiment(list(counts = full_data), rowData = gene_info, colData = cell_info)
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
.read_from_sparse <- function(path) {
    barcode.loc <- file.path(path, "barcodes.tsv")
    gene.loc <- file.path(path, "genes.tsv")
    matrix.loc <- file.path(path, "matrix.mtx")

    list(
        mat=as(readMM(matrix.loc), "dgCMatrix"),
        cell.names=readLines(barcode.loc),
        gene.info=read.delim(gene.loc, header=FALSE, colClasses="character", stringsAsFactors=FALSE, quote="", comment.char="")
    )
}

#' @importFrom rhdf5 h5ls h5read
#' @importFrom HDF5Array TENxMatrix
#' @importFrom utils head
.read_from_hdf5 <- function(path, group=NULL) {
    if (is.null(group)) {
        available <- h5ls(path, recursive=FALSE)
        available <- available[available$otype=="H5I_GROUP",]

        if (nrow(available) > 1L) {
            to.see <- head(available$name, 3)
            if (length(to.see)==3L) {
                to.see[3] <- "..."
            }
            stop("more than one available group (", paste(to.see, collapse=", "), ")")
        } else if (nrow(available) == 0L) {
            stop("no available groups")
        }
        group <- available$name
    }

    mat <- TENxMatrix(path, group)
    dimnames(mat) <- NULL # for consistency.
    list(
        mat=mat,
        cell.names=as.character(h5read(path, paste0(group, "/barcodes"))),
        gene.info=data.frame(
            as.character(h5read(path, paste0(group, "/genes"))),
            as.character(h5read(path, paste0(group, "/gene_names"))),
            stringsAsFactors=FALSE
        )
    )
}
