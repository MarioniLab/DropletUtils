# A big data mock-up:
library(Matrix)
set.seed(100)
ngenes <- 100
ambient.prof <- 1:100
ambient.prof <- ambient.prof/sum(ambient.prof)
total.count <- c(runif(10000, 0, 100), # to get a stable estimate of the ambient profile.
                 runif(1000, 100, 500)) 

ambient.counts <- matrix(rpois(ngenes*length(total.count), lambda=outer(ambient.prof, total.count)), 
                         ncol=length(total.count), nrow=ngenes)
ambient.counts <- as(ambient.counts, "dgCMatrix")

cell.counts <- matrix(rpois(ngenes*100, lambda=ambient.prof*1000), ncol=100, nrow=ngenes)
cell.counts <- as(cell.counts, "dgCMatrix")

new.prof <- runif(ngenes)
new.prof <- new.prof/sum(new.prof)
cell.counts2 <- matrix(rpois(ngenes*100, lambda=new.prof*200), ncol=100, nrow=ngenes)
cell.counts2 <- as(cell.counts2, "dgCMatrix")

my.counts <- cbind(ambient.counts, cell.counts, cell.counts2)

