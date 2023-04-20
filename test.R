## dyn.load("src/di_nuc_regions.so")
source("di_nuc_regions.R")

## LG_chr1.fa is a two line file containing the sequence on the second line.

lp.seq <- readLines("LP_chr1.fa")[2]
nchar(lp.seq)
## [1] 40732676
## 40 million base pairs in a single read.
## this is too big for me to put in github. I will find a better example later on.

## for troubleshooting:
cpg.regions <- .Call("di_nuc_regions", c(substring(lp.seq, 1, 100), "TA"), 14L)

cpg.regions <- .Call("di_nuc_regions", c(substring(lp.seq, 1, 100), "TA"), 6L)

system.time(
    cpg.regions <- diNucRegions(lp.seq, "CG", 14)
)
##  user  system elapsed 
## 0.315   0.000   0.319 

system.time(
    cpg.regions <- .Call("di_nuc_regions", c(lp.seq, "CG"), 14L, to.df=TRUE)
)
##  user  system elapsed 
## 0.276   0.000   0.277 
## timing more dependant on what else is going on on the machine.

dim(cpg.regions)
## [1] 519374      7

## defining 519 thousand CpG regions.. 
