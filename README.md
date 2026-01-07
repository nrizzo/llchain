# clc-viz-experiments
Tested on GCC version 15.2.1. Compile with
```console
git submodule update --init ext/mummer
make -j
```

```
Usage: clc-viz [-m global/semiglobal] [-a MUM/MEM] [-l minanchorlength] [-t
text.fasta(.gz)] [-q query.fasta(.gz)] [--chainx|--chainx-opt]
[--chainx-original-magic-numbers] [--custom-anchors anchors.mummer]
[--all-to-all] [--random-anchors ANCHORNUM] [-g gap-gap-ld.bmp] [-r INT]
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
      --all-to-all              Pairwise comparisons (queries)  (default=off)
      --chainx                  Chain with at-cg/ChainX original algorithm
                                  (default=off)
      --chainx-opt              Chain with algbio/ChainX optimal algorithm
                                  (default=off)
      --chainx-original-magic-numbers
                                In ChainX mode, use original magic numbers B =
                                  100, alpha = 4 instead of variable B >= 100,
                                  alpha = 4  (default=off)
      --custom-anchors=PATH     Do not index/query but read the anchors from
                                  this file (NB it should respect the same
                                  order as query file)
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
- [algbio/ChainX](https://github.com/algbio/ChainX) adapted as C++ module

## Dependencies
`gengetopt` for development

## Visualizations
- debug mode (`-g file.bmp`) creates an image file showing the randomly generated anchors (black), an optimal chain (red), and the recursions for case 2 projected forward (colored).
- the result of `--all-to-all` mode is a distance matrix in phylip format, that on Linux could be visualized as follows with `gnuplot`:
```
./clc-viz --all-to-all -a MEM -l 10 --query test/queries.fa > test/matrix.phylip
cat test/matrix.phylip | tail -n +2 | cut -d' ' -f2- | gnuplot -p -e "set view map; set size square; set yrange reverse; splot '-' matrix with image"
```

## TODOs
- split chaining time and backtracking time
- add minimizer anchors
- investigate edge case for case 2 (--random-anchors 100 -r 49929335 semi-global mode)
