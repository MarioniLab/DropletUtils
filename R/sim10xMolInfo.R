sim10xMolInfo <- function(prefix, umi.len=10, barcode.len=4, nmolecules=10000, swap.frac=0.2, 
    ngenes=20, nsamples=3, lambda=10, return.tab=FALSE)
# A function that creates a HDF5 file mimicking the molecule information from CellRanger.
# Used for testing the correctness of the swapping removal algorithm.    
{
    # Generating original reads.
    noriginal <- round(nmolecules * (1-swap.frac))
    ncells <- 4L^barcode.len
    cell <- sample(ncells, noriginal, replace = TRUE) - 1L
    umi <- sample(4L^as.integer(umi.len), noriginal, replace = FALSE)

    # Assigning each read to a gene and sample.
    gene <- sample(c(0L, seq_len(ngenes)), noriginal, replace = TRUE)
    sample <- sample(nsamples, noriginal, replace = TRUE)
   
    # Creating swapped reads.
    original <- data.frame(cell = cell, umi = umi, gene = gene, sample = sample) 
    swapped <- original[sample(nrow(original), nmolecules - noriginal),]

    samp.vec <- seq_len(nsamples)
    for (x in samp.vec) {
        current <- swapped$sample==x
        new.sample <- sample(samp.vec[-x], sum(current), replace=TRUE)
        swapped$sample[current] <- new.sample
    }
    
    # Simulating the number of reads.
    original$reads <- rpois(nrow(original), lambda = 10) + 1L
    swapped$reads <- 1L
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
        h5 <- h5createFile(out.file)
        h5write(current$cell, out.file, "barcode")
        h5write(current$umi, out.file, "umi")
        h5write(current$gene, out.file, "gene")
        h5write(rep(1, nrow(current)), out.file, "gem_group")
        h5write(current$reads, out.file, "reads")
        h5write(array(paste0("ENSG", seq_len(ngenes))), out.file, "gene_ids")
    }

    if (return.tab) { 
        return(list(files=out.files, original=original, swapped=swapped))
    } else {
        return(out.files)
    }
}



