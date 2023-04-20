## some magic that determines the location of this file
dyn.load( paste(dirname(sys.frame(1)$ofile), "src", "di_nuc_regions.so", sep="/") )

## takes a single sequence to scan
diNucRegions <- function(seq, dinuc, score, to.df=FALSE){
    tmp <- .Call("di_nuc_regions", c(seq, dinuc), as.integer(score));
    ## convert to base 1 (this should perhaps be done in the C-code)
    tmp[ ,c('start', 'end') ] <- tmp[,c('start', 'end') ] + 1
    if(!to.df)
        return(tmp)
    as.data.frame(tmp)
}

count.nucleotides <- function(seq){
    seq.n <- utf8ToInt(seq)
    tbl <- table(seq.n)
    names(tbl) <- strsplit( intToUtf8( as.integer(names(tbl)) ), "" )[[1]]
    tbl
}

## this assumes that the sequences don't actually use
## much unicode but are restricted to 8 bit representations
count.dinucleotides <- function(seq){
    seq.n <- utf8ToInt(seq)
    n <- length(seq.n)
    seq.dn <- bitwOr( bitwShiftL(seq.n[-n], 8), seq.n[-1] )
    tbl <- table(seq.dn)
    nm <- sapply(names(tbl), function(x){
        x <- as.integer(x)
        paste( intToUtf8(bitwShiftR(x, 8)), intToUtf8(bitwAnd(x, 0xFF)), sep="" )
    })
    names(tbl) <- nm
    tbl
}

## we can also count trinucleotides using this
## methods
count.trinucleotides <- function(seq){
    seq.n <- utf8ToInt(seq)
    n <- length(seq.n)
    i <- seq(1, n-2, 3)
    seq.dn <- bitwOr( bitwOr( bitwShiftL(seq.n[i], 16), bitwShiftL(seq.n[i+1], 8)), seq.n[i+2] )
    tbl <- table(seq.dn)
    nm <- sapply(names(tbl), function(x){
        x <- as.integer(x)
        paste(intToUtf8(bitwShiftR(x, 16)),
              intToUtf8(bitwAnd(bitwShiftR(x, 8), 0xFF)),
              intToUtf8(bitwAnd(x, 0xFF)), sep="" )
    })
    names(tbl) <- nm
    tbl
}

## given counts of individual nucleotides provide the
## expected frequency of dinucleotides
exp.dinuc <- function(n.count){
    s <- sum(n.count)
    n.f <- n.count / s
    n.f %*% t(n.f)
}

## an inefficient function to read fasta files:
## This probably has the same name as the similar
## function in many packages. It may be a good idea
## to comment this out.
read.fa <- function(fn){
    lines <- readLines(fn)
    id.l <- grep(">", lines)
    beg.l <- id.l + 1
    end.l <- c(id.l-1, length(lines))[-1]
    seq <- sapply( 1:length(beg.l), function(i){
        paste( lines[beg.l[i]:end.l[i]], collapse="" )
    })
    id <- sub(">([^ ]+).*", "\\1", lines[id.l])
    titles <- sub(">(.+)", "\\1", lines[id.l])
    names(seq) <- id
    list(seq=seq, titles=titles)
}
              
