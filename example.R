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
## 2912

sum( gpc.regions$dinuc_N )
## 4230

## suggesting that CpG dinucleotides appear to be depleted from
## the genome.

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

## We can count all nucleotides, dinucleotides and trinucleotides in the sequence using a few bitshifts
## and table: (see di_nuc_regions.R)
n.count <- count.nucleotides(chr1$seq)

dn.count <- count.dinucleotides( chr1$seq )
barplot(sort(dn.count[ !grepl("N", names(dn.count)) ]), las=2)

tri.count <- count.trinucleotides( chr1$seq )
barplot(sort(tri.count[ !grepl("N", names(tri.count)) ]), las=2)

exp.df <- exp.dinuc(n.count)

exp.df * sum(dn.count)

## check observed / expected ratios for each CpG dinucleotide
dn.nm <- strsplit( names(dn.count), "" )
names(dn.nm) <- names(dn.count)

obs.exp <- sapply( 1:length(dn.count), function(i){
    dn.count[i] / (sum(dn.count) * exp.df[ dn.nm[[i]][1], dn.nm[[i]][2] ])
})

## make a plot!
pdf("dn_obs_expected.pdf", width=7, height=7)
barplot( sort(obs.exp[ !grepl("N", names(obs.exp)) ]), las=2 )
abline(h=1, col='red', lty=2)
dev.off()
