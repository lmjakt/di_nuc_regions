# Find dinucleotide islands in sequence data

This code provides a function that scans for regions that contain
high numbers of a specified nucleotide. It essentially implements
the same algorithm as "cpg_distribute"
(<https://github.com/lmjakt/cpg_distribute>, Tim Cutts and Gos Micklem) that
was written to identify CpG islands except that it allows the user to specify
the dinucleotide of interest.

To use the code it is necessary to first compile the shared object (dll) from
`di_nuc_regions.c`:

```
cd src
R CMD SHLIB di_nuc_regions.c
```

To use the function `dyn.load()` must be called on the resulting `.so`
file. This is done in the `di_nuc_regions.R` code which contains a wrapper
function for the compiled code. Hence you can simply source this file and
you should be good to go.

For an example of usage see the `examples.R` file.

## Functions implemented:

### `diNucRegions(seq, dinuc, score, to.df=FALSE)`

Identifies regions with elevated frequencies of the dinucleotide specified by
`dinuc` in the sequence `seq`, using a running sum which is incremented by
`score` everytime a `dinuc` is observed and otherwise decremented by 1. The
score must be positive and regions are defined as starting from a 0 to
positive transition and ending at a maximum. This is equivalent to finding
local alignments using the *Smith Waterman* algorithm.

The function calls the compiled function `di_nuc_regions` through the `.Call`
interface. The resulting matrix is converted to a dataframe if `to.df` is
`TRUE`.

### `count.nucleotides(seq)`

Enumerates the numbers of different residues in a `seq`. Uses `utf8ToInt()`
and `intToUtf8()` for efficiency. This means that the `seq` must be a character
vector of length 1.

### `count.dinucleotides(seq)`

Counts dinucleotides in `seq` using `utf8ToInt()` in combination with bitwise
shifts to encode dinucleotides as integers. `seq` must be a character vector
of length 1.

### `count.trinucleotides(seq)`

Counts trinucleotides in `seq` using `utf8ToInt()` in combination with bitwise
shifts to encode trinucleotides as integers. `seq` must be a character vector
of length 1.

### `exp.dinuc(n.count)`

Calculates the expected dinucleotide frequencies based on counts of individual
nucleotides as returned by `count.nucleotides()`.

### `read.fa(fn)`

Parses the fasta file given by `fn`. Returns a list containing two character
vectors: `$seq` and `$title`. `$seq` contains the sequences with the
identifiers as names. `$title` contains the full sequence titles minus the
`>`.

Note that this is not an efficient function.
