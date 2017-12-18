sim10xMolInfo <- function(prefix, num_mols=10000, dup_frac=0.2, ncells=100, ngenes=20, nsamples=3, lambda=10)
# An internal function that creates a HDF5 file mimicking the output from CellRanger.
# Used only for testing the correctness of the swapping removal algorithm.    
{
    # Generating original reads.
    base <- round(num_mols * (1-dup_frac) )
    cell <- sample(ncells, base, replace = TRUE)
    umi <- sample(4^10, base, replace = FALSE)

    # Assigning each read to a gene and sample.
    gene <- sample(c(0, seq_len(ngenes)), base, replace = TRUE)
    sample <- sample(nsamples, base, replace = TRUE)
   
    # Creating swapped reads.
    tab <- data.frame(cell = cell, umi = umi, gene = gene, sample = sample) 
    dups <- tab[sample(nrow(tab), num_mols - base),]

    samp_vec <- seq_len(nsamples)
    for (x in samp_vec) {
        current <- dups$sample==x
        new.sample <- sample(samp_vec[-x], sum(current), replace=TRUE)
        dups$sample[current] <- new.sample
    }
    
    # Simulating the number of reads.
    tab <- rbind(tab, dups)
    tab$reads <- rpois(nrow(tab), lambda = 10) + 1

    # Writing them to 10X-like HDF5 files.
    out_files <- paste0(prefix, "_", seq_len(nsamples), ".h5")
    for (i in seq_along(out_files)){
        sample <- seq_len(nsamples)[i]
        out_file <- out_files[i]
        if(file.exists(out_file)) {
            file.remove(out_file)
        }
        
        chosen <- tab$sample==sample
        h5 <- H5Fcreate(out_file)
        h5write(tab$cell[chosen], out_file, "barcode")
        h5write(tab$umi[chosen], out_file, "umi")
        h5write(tab$gene[chosen], out_file, "gene")
        h5write(tab$reads[chosen], out_file, "reads")
        h5write(array(paste0("ENSG", seq_len(ngenes))), out_file, "gene_ids")
    }

    return(out_files)
}



