# This tests that swappedDrops works correctly.
# library(DropletUtils); library(testthat); source("test-swapping.R")

tmpdir <- tempfile()
dir.create(tmpdir)
ngenes <- 20L
barcode <- 4L
ncells <- 4L^barcode

##########################################

# Defining a reference function to compare the results.
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
    obs.swapped <- DropletUtils:::.findSwapped(combined$cell, combined$umi, combined$gene, combined$reads, min.frac)$swapped
    expect_identical(obs.swapped, is.swapped)
    return(all.counts)
}

##########################################

library(Matrix)
set.seed(5717)
test_that("Removal of swapped drops works correctly", {
    for (nmolecules in c(10, 100, 1000, 10000)) { 
        output <- sim10xMolInfo(tmpdir, return.tab=TRUE, barcode.length=barcode, nsamples=3,
                                ngenes=ngenes, nmolecules=nmolecules)

        # Figuring out the correspondence between cell ID and the reported barcode.
        # This involves a bit of work to decode the barcode.
        retainer <- vector("list", length(output$files))
        for (i in seq_along(retainer)) { 
            info <- read10xMolInfo(output$files[i], barcode)
            all.cells <- sort(unique(info$data$cell))
            encoded <- unlist(lapply(strsplit(all.cells, ""), FUN=function(i) { sum(4^(rev(seq_along(i))-1)*c(A=0, C=1, G=2, T=3)[i]) }))
            retainer[[i]] <- encoded + 1L # to get back to 1-based indexing.
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
        expect_equal(observed$cleaned, observed2$cleaned)
        
        observed3 <- swappedDrops(output$files, barcode, get.swapped=TRUE, get.diagnostics=TRUE, min.frac=min.frac)
        expect_equal(observed$cleaned, observed3$cleaned)
        expect_equal(observed$swapped, observed3$swapped)

        # Checking that the diagnostic field is consistent with the total.
        top.prop <- as.matrix(observed3$diagnostics)/rowSums(observed3$diagnostics)
        best.in.class <- max.col(top.prop)
        best.prop <- top.prop[(best.in.class - 1L) * nrow(top.prop) + seq_along(best.in.class)]
        for (s in seq_along(reference)) {
            expect_equal(sum(observed2$cleaned[[s]]), sum(best.in.class==s & best.prop >= min.frac))
        }
    }
})

test_that("swappedDrops functions correctly for silly inputs", {
    output <- sim10xMolInfo(tmpdir, barcode.length=barcode, nsamples=3, ngenes=ngenes, nmolecules=0)
    deswapped <- swappedDrops(output)
    for (ref in deswapped$cleaned) {
        expect_identical(dim(ref), c(ngenes, 0L))
    }   

    output <- sim10xMolInfo(tmpdir, barcode.length=barcode, nsamples=3, ngenes=0, nmolecules=0)
    deswapped <- swappedDrops(output)
    for (ref in deswapped$cleaned) {
        expect_identical(dim(ref), c(0L, 0L))
    }   
})


