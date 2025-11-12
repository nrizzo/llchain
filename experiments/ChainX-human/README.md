# ChainX-human experiments
Takes ~47GBs of RAM, requires `seqtk` installed.
```console
git submodule update --init ../../ext/ChainX && make -C ../../ext/ChainX chainX
../../ext/ChainX/data/human/get_dataset.sh
./run_experiment.sh
```
