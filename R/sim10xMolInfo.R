#' @importFrom rhdf5 h5write h5createFile h5write.default
#' @importFrom stats rpois
sim10xMolInfo <- function(prefix, nsamples=1, umi.length=10, barcode.length=4, 
    ngenes=20, nmolecules=10000, swap.frac=0.2, ave.read=10, return.tab=FALSE)
# A function that creates a HDF5 file mimicking the molecule information from CellRanger.
# Used for testing the correctness of the swapping removal algorithm.   
#
# written by Jonathan Griffiths
# with modifications by Aaron Lun
# created 18 December 2017    
{
    if (nsamples==1L) {
        # Swapping disabled if there's only one sample being simulated.
        swap.frac <- 0
    }

    # Generating original molecules (barcode.length <= 15, otherwise 32-bit integer will overflow). 
    noriginal <- round(nmolecules * (1-swap.frac))
    ncells <- 4L^as.integer(barcode.length)
    cell <- sample(ncells, noriginal, replace = TRUE) - 1L

    # UMIs are sampled without replacement, to guarantee uniqueness across all samples (w/o swapping).
    umi <- sample(4L^as.integer(umi.length), noriginal, replace = FALSE) 

    # Assigning each molecule to a gene and sample.
    gene <- sample(c(0L, seq_len(ngenes)), noriginal, replace = TRUE)
    sample <- sample(nsamples, noriginal, replace = TRUE)
    original <- DataFrame(cell = cell, umi = umi, gene = gene, sample = sample) 
   
    # Creating swapped molecules (unless nsamples==1L, in which case we skip this).
    if (nsamples > 1L) { 
        swapped <- original[sample(nrow(original), nmolecules - noriginal),]
        samp.vec <- seq_len(nsamples)
        new.sample <- swapped$sample
        for (x in samp.vec) {
            current <- swapped$sample==x
            new.sample[current] <- sample(samp.vec[-x], sum(current), replace=TRUE)
        }
        swapped$sample <- new.sample
    } else {
        swapped <- original[integer(0),]
    }
    
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
        h5 <- h5createFile(out.file)
        h5write(current$cell, out.file, "barcode") # technically should be saved as 64-bit, but not possible here.
        h5write(current$umi, out.file, "umi")
        h5write(current$gene, out.file, "gene")
        h5write(rep(1, nrow(current)), out.file, "gem_group")
        h5write(current$reads, out.file, "reads")
        h5write(array(sprintf("ENSG%i", seq_len(ngenes))), out.file, "gene_ids")
    }

    if (return.tab) { 
        return(list(files=out.files, original=original, swapped=swapped))
    } else {
        return(out.files)
    }
}
