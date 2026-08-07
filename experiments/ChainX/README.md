# ChainX experiments
This experiment replicates Tables 1 and 3 from the [`ChainX` paper](https://doi.org/10.1089/cmb.2022.0266) and compares the cost of the output chains between `ChainX`, `ChainX-opt`, and `llchain`. First, compile [`algbio/ChainX`](https://github.com/algbio/ChainX), `llchain`, and obtain the test datasets with commands
```console
git submodule update --init ../../ext/mummer
make -j $(nproc) -C ../../ llchain
git submodule update --init ../../ext/ChainX
make -C ../../ext/ChainX chainX
```

Then, run the experiment with command
```console
./run_experiment.sh
```

Afterwards, the experiment results (in folder `output`) can be checked and shown again with command
```console
./show_results.sh
```
