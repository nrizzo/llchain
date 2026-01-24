# genome-to-genome experiments
First, make sure you have the [NCBI datasets command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) installed and visible by your PATH environment variable.
Then get the datasets (~2GB of disk space) with script
```console
./get_datasets.sh
```

Then, make sure you have (`seqtk`)[https://github.com/lh3/seqtk] installed and visible by your `PATH` environment variable. Finally, run the experiment with commands
```console
git submodule update --init ../../ext/mummer && make -j4 -C ../../
./run_experiment.sh
```

Afterwards, the experiment results (in folder `output`) can be recalculated and shown again with command
```console
./show_results.sh
```
