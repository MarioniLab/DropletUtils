# This tests that swappedDrops works correctly.
# library(DropletUtils); library(testthat); source("test-swapping.R")

tmpdir <- tempfile()
dir.create(tmpdir)
ncells <- 100L
ngenes <- 20L
output <- sim10xMolInfo(tmpdir, return.tab=TRUE, ncells=ncells)
barcode <- ceiling(logb(ncells, 4))

test_that("Extraction of molecule information fields works correctly", {
    for (i in seq_along(output$files)) {
        ref.original <- output$original[output$original$sample==i,]        
        ref.swapped <- output$swapped[output$swapped$sample==i,]
        combined <- rbind(ref.original, ref.swapped)
        combined <- combined[combined$gene<ngenes,]

        current <- DropletUtils:::.readHDF5Data(output$files[i], barcode_len=barcode)
        expect_identical(as.integer(current$data$umi), combined$umi)
        expect_identical(as.integer(current$data$gene), combined$gene+1L)
        expect_identical(as.integer(current$data$reads), combined$reads)
        expect_identical(as.integer(current$data$gem_group), rep(1L, nrow(combined)))
        expect_identical(length(current$anno$cells), ncells)
        expect_identical(length(current$anno$genes), ngenes)
       
        # Checking that there is a 1:1 relationship between the cell barcodes and cell IDs.

        # Checking that using too little barcode length underestimates the number of cells.
        current2 <- DropletUtils:::.readHDF5Data(output$files[i], barcode_len=barcode - 1L)
        expect_true(length(current2$anno$cells) < ncells)
        current3 <- DropletUtils:::.readHDF5Data(output$files[i], barcode_len=barcode + 1L)
        expect_true(length(current3$anno$cells)==ncells)
    }
})

###########################
###########################

REFFUN <- function(original, swapped, min.frac) {
    combined <- rbind(original, swapped)
    combined <- combined[combined$gene<ngenes,] # Removing "unmapped" reads
    marking <- paste(combined$umi, combined$gene, combined$cell)    
    ref <- split(seq_len(nrow(combined)), marking)

    nsamples <- length(unique(original$sample))
    all.counts <- vector("list", nsamples)
    for (i in seq_len(nsamples)) { 
        all.counts[[i]] <- matrix(0, ngenes, ncells)
    }

    is.swapped <- !logical(nrow(combined))
    for (mol in seq_along(ref)) {
        current <- ref[[mol]]
        cur.reads <- combined$reads[current]
        all.props <- cur.reads/sum(cur.reads)
        chosen <- which.max(all.props)

        if (all.props[chosen] >= min.frac) {
            s <- combined$sample[current[chosen]]
            cur.gene <- combined$gene[current[chosen]] + 1L
            cur.cell <- combined$cell[current[chosen]]
            all.counts[[s]][cur.gene, cur.cell] <- all.counts[[s]][cur.gene, cur.cell] + 1
            is.swapped[current[chosen]] <- FALSE
        }
    }

    # Checking that the logical vector is correct.
    obs.swapped <- DropletUtils:::.findSwapped(marking, combined$reads, min.frac)
    expect_identical(obs.swapped, is.swapped)
    return(all.counts)
}

test_that("Removal of swapped drops works correctly", {
    # Figuring out the correspondence between cell ID and the reported barcode.
    cell.barcode <- .Call(DropletUtils:::cxx_get_cell_barcodes, output$files[1], "barcode", barcode)
    cell.id <- rhdf5::h5read(output$files[1], "barcode")
    barcodes.by.id <- cell.barcode[match(seq_len(ncells), cell.id)]

    # Constructing total matrices:
    combined <- rbind(output$original, output$swapped)
    combined <- combined[combined$gene<ngenes,] # Removing "unmapped" reads
    total.mat <- vector("list", length(output$files))
    for (s in seq_along(total.mat)) { 
        current.tab <- combined[combined$sample==s,]
        total.mat[[s]] <- Matrix::sparseMatrix(i=current.tab$gene+1, j=current.tab$cell,
                                               x=rep(1, nrow(current.tab)), 
                                               dims=c(ngenes, ncells))
    }

    # Matching them up for a specified min.frac of varying stringency.
    for (min.frac in c(0.5, 0.7, 1)) { 
        observed <- swappedDrops(output$files, barcode, get.swapped=TRUE, min.frac=min.frac)
        reference <- REFFUN(output$original, output$swapped, min.frac)
    
        # Checking that the cleaned object is correct.
        for (s in seq_along(reference)) {
            obs.mat <- as.matrix(observed$cleaned[[s]][,barcodes.by.id])
            ref.mat <- reference[[s]]
            rownames(ref.mat) <- rownames(obs.mat)
            colnames(ref.mat) <- barcodes.by.id
            expect_equal(obs.mat, ref.mat)  
        }

        # Checking that everything adds up to the total.
        for (s in seq_along(reference)) { 
            total <- (observed$cleaned[[s]] + observed$swapped[[s]])[,barcodes.by.id]
            ref.total <- total.mat[[s]]
            dimnames(ref.total) <- dimnames(total)
            expect_equal(ref.total, total)
        }
    }
})
