# ChainX-human experiments
This experiment replicates the experiment of the [`ChainX-opt` paper](https://doi.org/10.48550/arXiv.2506.11750), that is, chaining the MUMs or MEMs of 100k HiFi reads (HG002) against the T2T-CHM13 reference (semiglobal mode). First, make sure you have [`seqtk`](https://github.com/lh3/seqtk) installed (and it is visible by your PATH environment variable). Get the (concatenated) T2T-CHM13v2.0 reference and HG002 reads (~26GB of disk space) with commands
```console
git submodule update --init ../../ext/ChainX
sed -i 's/-s 1/-s 2/' ../../ext/ChainX/data/human/get_dataset.sh
../../ext/ChainX/data/human/get_dataset.sh
```

Then, make sure you have a working `python3` installation with `numpy`, `matplotlib`, and run the experiment (~47GB of RAM are needed) with commands
```console
git submodule update --init ../../ext/mummer
make -j $(nproc) ../../
./run_experiment.sh
```
Then `output/plot.svg` will be the plot from the paper. Afterwards, the experiment results (in folder `output`) can be recalculated/replotted with command
```console
./show_results.sh
```
