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
              
