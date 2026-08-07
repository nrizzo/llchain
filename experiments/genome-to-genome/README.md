# genome-to-genome experiments
**WIP/unpublished** experiment on computing the anchored edit distance of human vs chimp, plant vs plant genomes. First, make sure you have the [NCBI datasets command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) installed (and it is visible by your PATH environment variable). Then get the datasets (~9GB) with script
```console
./get_datasets.sh
```

To run the experiment, make sure you have [`seqtk`](https://github.com/lh3/seqtk) installed and execute commands
```console
git submodule update --init ../../ext/mummer && make -j $(nproc) -C ../../
./run_experiment.sh
```
