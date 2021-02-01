#' @importFrom S4Vectors DataFrame
#' @importFrom stats rpois
simSwappedMolInfo <- function(prefix, nsamples=2, umi.length=10, barcode.length=4, 
    ngenes=20, nmolecules=10000, swap.frac=0.2, ave.read=10, 
    version=c("2", "3"), return.tab=FALSE)
# A function that creates a HDF5 file mimicking the molecule information from CellRanger.
# Used for testing the correctness of the swapping removal algorithm.   
#
# written by Jonathan Griffiths
# with modifications by Aaron Lun
# created 18 December 2017    
{
    # Generating original molecules (barcode.length <= 15, otherwise 32-bit integer will overflow). 
    noriginal <- round(nmolecules * (1-swap.frac))
    ncells <- 4L^as.integer(barcode.length)
    cell <- sample(ncells, noriginal, replace = TRUE)

    # UMIs are sampled without replacement, to guarantee uniqueness across all samples (w/o swapping).
    umi <- sample(4L^as.integer(umi.length), noriginal, replace = FALSE) 

    # Assigning each molecule to a gene and sample.
    gene <- sample(ngenes+1L, noriginal, replace = TRUE)
    sample <- sample(nsamples, noriginal, replace = TRUE)
    original <- DataFrame(cell = cell, umi = umi, gene = gene, sample = sample, gem_group=rep(1L, noriginal))

    version <- match.arg(version)
    if (version=="3") {
        original$library <- rep(1L, nrow(original))
    }

    # Creating swapped molecules.
    swapped <- original[sample(nrow(original), nmolecules - noriginal),]
    samp.vec <- seq_len(nsamples)
    new.sample <- swapped$sample
    for (x in samp.vec) {
        current <- swapped$sample==x
        new.sample[current] <- sample(samp.vec[-x], sum(current), replace=TRUE)
    }
    swapped$sample <- new.sample
    
    # Simulating the number of reads (swapped reads only get 1).
    original$reads <- rpois(nrow(original), lambda = ave.read) + 1L
    swapped$reads <- rep(1L, nrow(swapped))
    fulltab <- rbind(original, swapped)

    # Writing them to 10X-like HDF5 files.
    out.files <- paste0(prefix, ".", seq_len(nsamples), ".h5")

    for (i in seq_along(out.files)){
        sample <- seq_len(nsamples)[i]
        out.file <- out.files[i]
        if(file.exists(out.file)) {
            file.remove(out.file)
        }
        
        current <- fulltab[fulltab$sample==sample,]
        .write_stripped_mol_info(out.file, current, 
            barcode.length=barcode.length,
            gene.names=sprintf("ENSG%i", seq_len(ngenes)),
            feature.types=rep("A", ngenes),
            library.info=list(list(library_type="A", library_id=0, gem_group=1)),
            version=version)
    }

    if (return.tab) { 
        list(files=out.files, original=original, swapped=swapped)
    } else {
        out.files
    }
}

#' @importFrom S4Vectors DataFrame
simBasicMolInfo <- function(out.file, ngems=1, umi.length=10, barcode.length=4, 
    ngenes=20, nmolecules=10000, ave.read=10, version=c("2", "3"), return.tab=FALSE)
{
    ncells <- 4L^as.integer(barcode.length)
    cell <- sample(ncells, nmolecules, replace = TRUE)
    umi <- sample(4L^as.integer(umi.length), nmolecules)
    gene <- sample(ngenes+1L, nmolecules, replace = TRUE)
    reads <- pmax(1L, rpois(nmolecules, lambda = ave.read))
    gem_group <- sample(ngems, nmolecules, replace=TRUE)
    current <- DataFrame(cell = cell, umi = umi, gene = gene, reads=reads, gem_group=gem_group)

    version <- match.arg(version)
    if (version=="3") {
        library.info <- list(
            list(library_type="A", library_id=0L, gem_group=1L),
            list(library_type="B", library_id=1L, gem_group=1L),
            list(library_type="C", library_id=2L, gem_group=1L)
         )
         choices <- vapply(library.info, "[[", i="library_type", FUN.VALUE="")
         feature.types <- sample(choices, ngenes, replace=TRUE)

         m <- match(feature.types[current$gene], choices)
         m[is.na(m)] <- length(choices) + 1L
         current$library <- m
    } else {
        library.info <- NULL
        feature.types <- NULL 
    }

    .write_stripped_mol_info(out.file, current, 
        barcode.length=barcode.length,
        gene.names=sprintf("ENSG%i", seq_len(ngenes)),
        feature.types=feature.types,
        library.info=library.info,
        version=match.arg(version))

    if (return.tab) {
        list(files=out.file, original=current)
    } else {
        out.file
    }
}

#' @importFrom rhdf5 h5write h5createGroup h5createFile
.write_stripped_mol_info <- function(out.file, current, barcode.length, 
    gene.names, feature.types, library.info, version="2")
{
    out.file <- path.expand(out.file) # protect against tilde's.

    unlink(out.file)
    h5 <- h5createFile(out.file)

    # Technically these should be saved as 64-bit, but not possible here.
    if (version=="2") {
        h5write(current$cell - 1L, out.file, "barcode") 
    } else {
        actual.barcodes <- factor(.unmask_barcode(current$cell - 1L, barcode.length))
        h5write(as.integer(actual.barcodes) - 1L, out.file, "barcode_idx") 
        h5write(levels(actual.barcodes), out.file, "barcodes")
    }

    gene.field <- if (version=="2") "gene" else "feature_idx"
    h5write(current$gene - 1L, out.file, gene.field)

    read.field <- if (version=="2") "reads" else "count"
    h5write(current$reads, out.file, read.field)

    if (version=="2") {
        h5write(gene.names, out.file, "gene_ids")
    } else {
        h5createGroup(out.file, "features")
        h5write(gene.names, out.file, "features/id")
        h5write(feature.types, out.file, "features/feature_type")
    }

    if (version=="3") {
        h5write(current$library - 1L, out.file, "library_idx")
        h5write(as.character(jsonlite::toJSON(library.info, auto_unbox=TRUE)), out.file, "library_info")
    }

    h5write(current$umi, out.file, "umi")
    h5write(current$gem_group, out.file, "gem_group")
}

.unmask_barcode <- function(idx, blen) {
    seqs <- vector("list", blen)
    for (i in seq_len(blen)) {
        seqs[[i]] <- c("A", "C", "G", "T")[idx %% 4 + 1L]
        idx <- floor(idx/4)
    }
    do.call(paste0, rev(seqs))
}
