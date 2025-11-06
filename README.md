# clc-viz-experiments
Tested on GCC version 15.2.1. Compile with
```console
git submodule update --init ext/mummer
make -j
```

```
Usage: clc-viz [-m global/semiglobal] [-t text.fasta(.gz)] [-q query.fasta(.gz)]
[--random-anchors ANCHORNUM] [-g gap-gap-ld.bmp] [-r INT]
Verification, WIP implementation, and visualization of linearithmic-time
colinear chaining

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
  -t, --text=PATH               Text sequences file
  -q, --query=PATH              Query sequences file
  -m, --mode=MODE               Chaining mode (global/semiglobal)
                                  (default=`global')
  -a, --anchor-type=ANCHOR      (MUM/MEM)  (default=`MUM')
  -l, --anchor-length=LENGTH    Minimum anchor length  (default=`20')
      --random-anchors=ANCHORNUM
                                Number of random anchors to generate
                                  (default=`-1')
  -g, --debug-case-two-output-file=BMPFILE
                                Visualize case 2 in this file (BMP format)
                                  (default=`')
  -r, --random-seed=INT         Seed for the PRNG (-1 is different at every
                                  invocation)  (default=`-1')
```

## External libraries
- [kseq](https://github.com/lh3/seqtk) for FASTA parsing
- [mummer/essaMEM](https://github.com/mummer4/mummer.git) for seed finding
- [grid_to_bmp](https://people.sc.fsu.edu/~jburkardt/cpp_src/grid_to_bmp/grid_to_bmp.html) that we adapted for debugging

## Dependencies
`gengetopt` for development

## TODOs
- MEM tests from at-cg/ChainX
- PacBio HiFi long read test from algbio/ChainX
- implement all-to-all comparison
- investigate edge case for case 2 (--random-anchors 100 -r 49929335 semi-global mode)
