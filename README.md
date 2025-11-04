# clc-viz-experiments
Tested on GCC version 15.2.1.
```
Usage: clc-viz [-t text.fasta(.gz)] [-q query.fasta(.gz)] [--random-anchors ANCHORNUM]
[-g gap-gap-ld.bmp] [-r INT]
Verification, WIP implementation, and visualization of linearithmic-time
colinear chaining

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
  -t, --text=PATH               Text sequences file
  -q, --query=PATH              Query sequences file
  -n, --random-anchors=ANCHORNUM
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
- [grid_to_bmp](https://people.sc.fsu.edu/~jburkardt/cpp_src/grid_to_bmp/grid_to_bmp.html) that we adapted for debugging

## Dependencies
`gengetopt` for development

## TODOs
- reach feature-parity with [ChainX](https://github.com/algbio/ChainX/)
- stream queries from disk
- compare results with ChainX and test performance
- investigate edge case for case 2 (--random-anchors 100 -r 49929335 semi-global mode)
