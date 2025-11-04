# clc-viz-experiments
Tested on GCC version 15.2.1.
```
Usage: clc-viz [-n ANCHORNUM] [-g gap-gap-ld.bmp]
Verification, WIP implementation, and visualization of linearithmic-time
colinear chaining

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
  -n, --anchors=ANCHORNUM       Number of anchors  (default=`5')
  -g, --gap-gap-lower-diagonal-output-file=BMPFILE
                                Visualize the gap-gap, lower diagonal case in
                                  this file (BMP format)
                                  (default=`gap-gap-ld.bmp')
  -r, --random-seed=INT         Seed for the PRNG (-1 is different at every
                                  invocation)  (default=`-1')
```

## External libraries
- [grid_to_bmp](https://people.sc.fsu.edu/~jburkardt/cpp_src/grid_to_bmp/grid_to_bmp.html) that we adapted for debugging

## Dependencies
`gengetopt` for development

## TODOs
- reach feature-parity with [ChainX](https://github.com/algbio/ChainX/)
- compare results with ChainX and test performance
- investigate edge case for case 2 (-n 100 -r 49929335 semi-global mode)
