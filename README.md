# llchain
Tested on GCC version >= 15.2.1. Build with
```console
git submodule update --init ext/mummer
make -j $(nproc)
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
./llchain --all-to-all -a MEM -l 10 --query test/queries.fa > test/matrix.phylip
cat test/matrix.phylip | tail -n +2 | cut -d' ' -f2- | gnuplot -p -e "set view map; set size square; set yrange reverse; splot '-' matrix with image"
```

## TODOs
- check -Wall warnings
- investigate edge case for case 2 (--random-anchors 100 -r 49929335 semi-global mode)
- investigate edge case for case 2 (--random-anchors=100 -r 696403780)
- see if the ChainX-opt "bug fix" about strict ChainX precedence affects performance
- use [CLI11](https://github.com/CLIUtils/CLI11) instead of gengetopt
- avoid using a list in case 2
