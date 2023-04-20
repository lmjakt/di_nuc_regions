source("di_nuc_regions.R")

chr1 <- read.fa("chr1.fa")

cpg.regions <- diNucRegions(chr1$seq, "CG", 14, to.df=TRUE)
gpc.regions <- diNucRegions(chr1$seq, "GC", 14, to.df=TRUE)

## check that it reports correctly:
with(cpg.regions, substring(chr1$seq, start[1], end[1]))
##         lp_chr1 
## "CGGTAGAAGATCCG" 

with(gpc.regions, substring(chr1$seq, start[1], end[1]))

par(mfrow=c(1,2))
hist(log2(cpg.regions$score), breaks=40)
hist(log2(gpc.regions$score), breaks=40)

sum( cpg.regions$dinuc_N )
## 2812

sum( gpc.regions$dinuc_N )
## 4230

## suggesting that CpG dinucleotides appear to be depleted from
## the genome.

## to count the number of nucleotides in the aboe sequence we can do:
## note that this only works for a single sequence..
tmp <- utf8ToInt( chr1$seq )
nuc.n <- table( tmp )
names(nuc.n) <- strsplit(intToUtf8( as.integer(names(nuc.n)) ), "")[[1]]
nuc.f <- nuc.n / sum(nuc.n)

## note:
## there are some Ns
## A/T and C/G are not completely symmetrical
## sequence is perhaps too short.

with( cpg.regions, sum(dinuc_N) / (sum(nuc.n) * nuc.f['C'] * nuc.f['G']) )
##         C 
## 0.6230237 

with( gpc.regions, sum(dinuc_N) / (sum(nuc.n) * nuc.f['C'] * nuc.f['G']) )
##         C 
## 0.9050103




