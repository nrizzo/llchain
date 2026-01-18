# ChainX-human experiments
First, make sure you have [`seqtk`](https://github.com/lh3/seqtk) installed and visible by your PATH environment variable.
Then get the T2T-CHM13v2.0 reference and concatenate it (~2GB of disk space) with commands
```console
git submodule update --init ../../ext/ChainX
../../ext/ChainX/data/human/get_dataset.sh
```

Then, run the experiment (~47GB of RAM) with commands
```console
git submodule update --init ../../ext/mummer && make -j4 -C ../../
./run_experiment.sh
```

Afterwards, the experiment results (in folder `output`) can be recalculated and shown again with command
```console
./show_results.sh
```
