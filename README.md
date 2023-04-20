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

More documentation may follow.

