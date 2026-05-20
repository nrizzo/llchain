# ChainX experiments
First, compile [`algbio/ChainX`](https://github.com/algbio/ChainX), `llchain`, and obtain the test datasets with commands
```console
git submodule update --init ../../ext/mummer
make -j4 -C ../../ llchain
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
