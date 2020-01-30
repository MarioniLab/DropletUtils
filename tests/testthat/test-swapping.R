# This tests that swappedDrops works correctly.
# library(DropletUtils); library(testthat); source("test-swapping.R")

tmpdir <- tempfile()
dir.create(tmpdir)
ngenes <- 20L
barcode <- 4L
ncells <- 4L^barcode

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
    return(all.counts)
}

library(Matrix)
set.seed(5717)
test_that("Removal of swapped drops works correctly", {
    for (nmolecules in c(10, 100, 1000, 10000)) { 
        output <- DropletUtils:::sim10xMolInfo(tmpdir, return.tab=TRUE, barcode.length=barcode, 
            nsamples=3, ngenes=ngenes, nmolecules=nmolecules)

        # Figuring out the correspondence between cell ID and the sorted barcode.
        was.mapped <- c(output$original$gene, output$swapped$gene)!=ngenes
        f <- factor(c(output$original$sample, output$swapped$sample), levels=seq_along(output$files))
        retainer <- split(c(output$original$cell, output$swapped$cell)[was.mapped], f[was.mapped], drop=FALSE)
        retainer <- lapply(retainer, FUN=function(i) {
            i <- unique(i)
            collected <- DropletUtils:::.unmask_barcode(i, barcode)
            stopifnot(!anyDuplicated(collected))
            i[order(collected)] + 1L # get back to 1-based indexing.
        })
   
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
                expect_true(all(reference[[s]][,-retainer[[s]]]==0))
            }
    
            # Checking that everything adds up to the total.
            for (s in seq_along(reference)) { 
                total <- observed$cleaned[[s]] + observed$swapped[[s]]
                ref.total <- total.mat[[s]][,retainer[[s]]]
                dimnames(ref.total) <- dimnames(total) 
                expect_equal(ref.total, total)
                expect_true(all(total.mat[[s]][,-retainer[[s]]]==0))
            }
        }
    }
})

test_that("Alternative input/output parameters work correctly", {
    for (nmolecules in c(10, 100, 1000, 10000)) { 
        output <- DropletUtils:::sim10xMolInfo(tmpdir, return.tab=TRUE, barcode.length=barcode, nsamples=3, ngenes=ngenes, nmolecules=nmolecules)

        # Further input/output tests.
        min.frac <- 0.9001
        observed <- swappedDrops(output$files, barcode, min.frac=min.frac)
        expect_equal(observed$swapped, NULL)
        expect_equal(observed$diagnostics, NULL)

        observed2 <- swappedDrops(output$files, barcode, get.swapped=TRUE, min.frac=min.frac)
        expect_equal(observed2$cleaned, observed$cleaned)
        expect_identical(lapply(observed2$swapped, dim), lapply(observed2$cleaned, dim))
        expect_equal(observed2$diagnostics, NULL)
        
        observed3 <- swappedDrops(output$files, barcode, get.swapped=TRUE, get.diagnostics=TRUE, min.frac=min.frac, hdf5.out=FALSE)
        expect_s4_class(observed3$diagnostics, "dgCMatrix")
        expect_equal(observed2$cleaned, observed3$cleaned)
        expect_equal(observed2$swapped, observed3$swapped)

        # Checking that the diagnostic field is consistent with the total.
        top.prop <- as.matrix(observed3$diagnostics)/rowSums(observed3$diagnostics)
        best.in.class <- max.col(top.prop)
        best.prop <- top.prop[(best.in.class - 1L) * nrow(top.prop) + seq_along(best.in.class)]
        for (s in seq_along(observed$cleaned)) {
            expect_equal(sum(observed$cleaned[[s]]), sum(best.in.class==s & best.prop >= min.frac))
        }

        # Checking that the HDF5 and sparse results are the same.
        observed4 <- swappedDrops(output$files, barcode, get.swapped=FALSE, get.diagnostics=TRUE, min.frac=min.frac, hdf5.out=TRUE)
#        expect_s4_class(observed4$diagnostics, "HDF5Array")
        expect_equivalent(as.matrix(observed3$diagnostics), as.matrix(observed4$diagnostics))
    }
})

test_that("swappedDrops functions correctly for silly inputs", {
    output <- DropletUtils:::sim10xMolInfo(tmpdir, barcode.length=barcode, nsamples=3, ngenes=ngenes, nmolecules=0)
    deswapped <- swappedDrops(output)
    for (ref in deswapped$cleaned) {
        expect_identical(dim(ref), c(ngenes, 0L))
    }   

    output <- DropletUtils:::sim10xMolInfo(tmpdir, barcode.length=barcode, nsamples=3, ngenes=0, nmolecules=0)
    deswapped <- swappedDrops(output)
    for (ref in deswapped$cleaned) {
        expect_identical(dim(ref), c(0L, 0L))
    }

    # Fails if you give it samples with different gene sets.
    tmpdir2 <- tempfile()
    dir.create(tmpdir2)

    o1 <- DropletUtils:::sim10xMolInfo(tmpdir, barcode.length=barcode, nsamples=3, ngenes=100, nmolecules=0)
    o2 <- DropletUtils:::sim10xMolInfo(tmpdir2, barcode.length=barcode, nsamples=3, ngenes=10, nmolecules=0)
    expect_error(swappedDrops(c(o1, o2)), "gene information differs")

    # Spits out a warning if you have multiple GEM groups.
    o1 <- DropletUtils:::sim10xMolInfo(tmpdir, barcode.length=barcode, nsamples=1, ngenes=100, nmolecules=100)
    rhdf5::h5write(sample(3L, 100, replace=TRUE), o1, "gem_group")
    expect_warning(swappedDrops(o1), "contains multiple GEM groups")

    # Responds to names in the sample paths.
    new.names <- LETTERS[seq_along(output)]
    deswapped2 <- swappedDrops(setNames(output, new.names), get.swapped=TRUE)
    expect_identical(names(deswapped2$cleaned), new.names)
    expect_identical(names(deswapped2$swapped), new.names)

    # removeSwappedDrops is not happy if manually specified lists don't match up.
    expect_error(removeSwappedDrops(cells=list(), umis=list(1L), genes=list(1L), nreads=list(1L), ref.genes=1:10), "lists are not")
    expect_error(removeSwappedDrops(cells=list(c("A", "C")), umis=list(1L), genes=list(1L), nreads=list(1L), ref.genes=1:10), "list vectors are not")
})
