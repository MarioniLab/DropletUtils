# This tests that swappedDrops works correctly.
# library(DropletUtils); library(testthat); source("test-swapping.R")

tmpdir <- tempfile()
dir.create(tmpdir)
ngenes <- 20L
barcode <- 4L
output <- sim10xMolInfo(tmpdir, return.tab=TRUE, barcode=barcode)

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

        # Checking that there is a 1:1 relationship between the cell barcodes and cell IDs.
        by.barcode <- split(combined$cell, current$data$cell)
        expect_true(all(lengths(lapply(by.barcode, unique))==1L))
        by.cell.id <- split(current$data$cell, combined$cell)
        expect_true(all(lengths(lapply(by.cell.id, unique))==1L))

        # Checking that using too little barcode length underestimates the number of cells.
        current2 <- DropletUtils:::.readHDF5Data(output$files[i], barcode_len=barcode - 1L)
        expect_true(length(current2$anno$cells) < length(current$anno$cells))
        current3 <- DropletUtils:::.readHDF5Data(output$files[i], barcode_len=barcode + 1L)
        expect_true(length(current3$anno$cells) == length(current$anno$cells))

        # Checking annotation.
        expect_identical(current$anno$cells, sort(unique(current$data$cell)))
        expect_identical(length(current$anno$genes), ngenes)
    }
})

###########################
###########################

ncells <- 4L^barcode
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
            cur.cell <- combined$cell[current[chosen]] + 1L
            all.counts[[s]][cur.gene, cur.cell] <- all.counts[[s]][cur.gene, cur.cell] + 1
            is.swapped[current[chosen]] <- FALSE
        }
    }

    # Checking that the logical vector is correct.
    obs.swapped <- DropletUtils:::.findSwapped(marking, combined$reads, min.frac)
    expect_identical(obs.swapped, is.swapped)
    return(all.counts)
}

library(Matrix)
test_that("Removal of swapped drops works correctly", {
    for (nmolecules in c(10, 100, 1000, 10000)) { 
        output <- sim10xMolInfo(tmpdir, return.tab=TRUE, barcode=barcode, nmolecules=nmolecules)

        # Figuring out the correspondence between cell ID and the reported barcode.
        # This involves a bit of work when not all cells are available.
        retainer <- vector("list", length(output$files))
        for (i in seq_along(retainer)) { 
            cell.barcode <- .Call(DropletUtils:::cxx_get_cell_barcodes, output$files[i], "barcode", barcode)
            cell.id <- rhdf5::h5read(output$files[i], "barcode")
            m <- match(seq_len(ncells)-1, cell.id)
            keep <- which(!is.na(m))
            retainer[[i]] <- keep[order(cell.barcode[m[keep]])]
        }
    
        # Constructing total matrices:
        combined <- rbind(output$original, output$swapped)
        combined <- combined[combined$gene<ngenes,] # Removing "unmapped" reads
        total.mat <- vector("list", length(output$files))
        for (s in seq_along(total.mat)) { 
            current.tab <- combined[combined$sample==s,]
            total.mat[[s]] <- sparseMatrix(i=current.tab$gene+1, 
                                           j=current.tab$cell+1,
                                           x=rep(1, nrow(current.tab)), 
                                           dims=c(ngenes, ncells))
        }
    
        # Matching them up for a specified min.frac of varying stringency.
        for (min.frac in c(0.5, 0.7, 1)) { 
            observed <- swappedDrops(output$files, barcode, get.swapped=TRUE, min.frac=min.frac)

            # Checking that the cleaned object is correct.
            reference <- REFFUN(output$original, output$swapped, min.frac)
            for (s in seq_along(reference)) {
                obs.mat <- as.matrix(observed$cleaned[[s]])
                ref.mat <- reference[[s]][,retainer[[s]]]
                dimnames(ref.mat) <- dimnames(obs.mat)
                expect_equal(obs.mat, ref.mat)  
            }
    
            # Checking that everything adds up to the total.
            for (s in seq_along(reference)) { 
                total <- (observed$cleaned[[s]] + observed$swapped[[s]])
                ref.total <- total.mat[[s]][,retainer[[s]]]
                dimnames(ref.total) <- dimnames(total) 
                expect_equal(ref.total, total)
            }
        }

        # Further input/output tests.
        min.frac <- 0.9
        observed <- swappedDrops(output$files, barcode, get.swapped=TRUE, min.frac=min.frac)
        observed2 <- swappedDrops(output$files, barcode, min.frac=min.frac)
        expect_equal(observed$cleaned, observed2)
        
        observed3 <- swappedDrops(output$files, barcode, get.swapped=TRUE, get.diagnostics=TRUE, min.frac=min.frac)
        expect_equal(observed$cleaned, observed3$cleaned)
        expect_equal(observed$swapped, observed3$swapped)

        # Checking that the diagnostic field is consistent with the total.
        top.prop <- as.matrix(observed3$diagnostics)/rowSums(observed3$diagnostics)
        best.in.class <- max.col(top.prop)
        best.prop <- top.prop[(best.in.class - 1L) * nrow(top.prop) + seq_along(best.in.class)]
        for (s in seq_along(reference)) {
            expect_equal(sum(observed2[[s]]), sum(best.in.class==s & best.prop >= min.frac))
        }
    }
})
