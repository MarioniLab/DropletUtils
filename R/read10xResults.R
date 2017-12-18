read10xResults <- function(samples) 
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
    
    for (i in seq_len(nsets)) { 
        run <- samples[i]
        barcode.loc <- file.path(run, "barcodes.tsv")
        gene.loc <- file.path(run, "genes.tsv")
        matrix.loc <- file.path(run, "matrix.mtx")
        
        ## read sparse count matrix and cell barcodes.
        data_mat <- readMM(matrix.loc)
        data_mat <- as(data_mat, "dgCMatrix")
        cell.names <- read.table(barcode.loc, header = FALSE, colClasses = "character")[[1]]

        full_data[[i]] <- data_mat
        gene_info_list[[i]] <- read.table(gene.loc, header = FALSE, colClasses = "character")
        cell_info_list[[i]] <- DataFrame(Sample = run, Barcode = cell.names)
    }

    # Checking gene uniqueness. 
    if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
        stop("gene information differs between runs")
    }
    gene_info <- gene_info_list[[1]]
    colnames(gene_info) <- c("ID", "Symbol")
    rownames(gene_info) <- gene_info$id
    
    # Forming the full data matrix.
    full_data <- do.call(cbind, full_data)
    rownames(full_data) <- gene_info$ID

    # Adding the cell data (only using as colnames if there is only 1 set - guaranteed unique).
    cell_info <- do.call(rbind, cell_info_list)
    if (nsets == 1L) {
        colnames(full_data) <- cell_info$Barcode
    }

    SingleCellExperiment(list(counts = full_data), rowData = gene_info, colData = cell_info)
}

