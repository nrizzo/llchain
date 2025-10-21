FLAGS=-std=c++20 -fmodules -O0 -g
CFLAGS=-std=c++20 -O0 -g
.PHONY : clean

clc-viz : src/clc-viz.cpp src/utils.o src/MinSegmentTree.o src/algo.o src/command-line-parsing/cmdline.o
	g++ -I src $(FLAGS) $^ -o clc-viz

src/command-line-parsing/cmdline.o : src/command-line-parsing/cmdline.h src/command-line-parsing/cmdline.c
	g++ $(CFLAGS) -c $^ -o $@
src/command-line-parsing/cmdline.c src/command-line-parsing/cmdline.h : src/command-line-parsing/config.ggo
	gengetopt \
		--input=./src/command-line-parsing/config.ggo \
		--output-dir=./src/command-line-parsing/ \
		--unnamed-opts

src/utils.o : src/utils.cppm ext/grid_to_bmp.hpp
	g++ $(FLAGS) -I ext -c $< -o $@
src/algo.o : src/algo.cppm src/utils.o src/MinSegmentTree.o
	g++ $(FLAGS) -c $< -o $@
src/MinSegmentTree.o : src/MinSegmentTree.cppm src/utils.o
	g++ $(FLAGS) -c $< -o $@

clean :
	-rm -f src/utils.o src/algo.o src/MinSegmentTree.o src/command-line-parsing/cmdline.o
	-rm -Rf gcm.cache
