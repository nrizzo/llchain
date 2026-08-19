VERSION=$(shell git rev-parse --short HEAD)
FLAGS=-std=c++20 -O3 -Wall -fmodules
CFLAGS=-std=c++20 -O3 -Wall
LDFLAGS=-lz
STDMODULES=iostream vector fstream tuple random algorithm limits cassert map set list string numeric sstream cmath array filesystem thread syncstream latch barrier optional atomic
STDMODULES_FAKEFILES=$(foreach s,$(STDMODULES),makefile.cache/$(s))
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
HTSLIB_LIB=$(ROOT_DIR)ext/htslib
HTSLIB_FLAGS=-L $(HTSLIB_LIB) -lhts -Wl,-rpath $(HTSLIB_LIB)
.PHONY : clean

LLCHAIN_OBJS=src/utils.o src/MinSegmentTree.o src/algo.o ext/grid_to_bmp.o ext/mummer_essaMEM_wrapper.o ext/kseq.o ext/chainx.o ext/cli11.o ext/htslib_wrapper.o
llchain : src/llchain.cpp $(LLCHAIN_OBJS) $(STDMODULES_FAKEFILES)
	g++ -I ext/concurrentqueue $(FLAGS) $(HTSLIB_FLAGS) -DVERSION="\"$(VERSION)\"" $< $(LLCHAIN_OBJS) -o llchain $(LDFLAGS) -Wno-global-module

src/utils.o : src/utils.cppm ext/grid_to_bmp.o $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
src/algo.o : src/algo.cppm src/utils.o src/MinSegmentTree.o $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
src/MinSegmentTree.o : src/MinSegmentTree.cppm src/utils.o $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
CLI11_INCL=ext/CLI11/include
ext/cli11.o : ext/CLI11/src/modules/CLI11.cppm $(shell find $(CLI11_INCL) -name *.hpp)
	g++ $(FLAGS) -I $(CLI11_INCL) -c $< -o $@
ext/grid_to_bmp.o : ext/grid_to_bmp.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
ext/kseq.o : ext/kseq.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
MUMMER_CPPS=ext/mummer/src/essaMEM/
MUMMER_INCL=ext/mummer/include/
MUMMER_CPP_FLAGS=-DNDEBUG
ext/mummer_essaMEM_wrapper.o : ext/mummer_essaMEM_wrapper.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) $(MUMMER_CPP_FLAGS) -I $(MUMMER_INCL) -I $(MUMMER_CPPS) -c $< -o $@
HTSLIB_INCL=$(ROOT_DIR)ext/htslib
ext/htslib_wrapper.o : ext/htslib_wrapper.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -I $(HTSLIB_INCL) -c $< -o $@
ext/chainx.o : ext/chainx.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -I $(MUMMER_INCL) -I $(MUMMER_CPPS) -c $< -o $@
$(STDMODULES_FAKEFILES) :
	g++ $(FLAGS) -xc++-system-header $(notdir $@)
	mkdir -p makefile.cache && touch $@

clean :
	-rm -f $(LLCHAIN_OBJS)
	-rm -Rf gcm.cache makefile.cache
