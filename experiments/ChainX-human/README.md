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

## Plot
Input: files
```
./data/chaining_times_human_mem_chainx
./data/chaining_times_human_mem_chainx-opt
./data/chaining_times_human_mem_llchain
./data/seeding_times_human_mem_chainx
./data/seeding_times_human_mem_chainx-opt
./data/seeding_times_human_mem_llchain
```
containing tab separated anchor-count and seeding/chaining time, like:
```
20922	0.0185096
36	1.1607e-05
8980	0.00180272
90932	0.0387818
34319	0.00626129
754	0.000107567
20657	0.00337901
6660	0.00807215
28253	0.00524187
166	3.774e-05
```

Then run `./plot.py` to generate `plot.pdf`.
