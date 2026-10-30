VERSION=$(shell git rev-parse --short HEAD)
COMPILER=$(shell $(CXX) --version | head -n 1)

ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
HTSLIB_LIB=$(ROOT_DIR)ext/htslib
HTSLIB_CFLAGS=-L $(HTSLIB_LIB) -lhts -Wl,-rpath $(HTSLIB_LIB)

LDFLAGS=-lz
CXXFLAGS=-std=c++20 -O3 -flto -march=native -DNDEBUG -Wall
GCC_CXXFLAGS=-fmodules
CLANG_CXXFLAGS=-stdlib=libc++ -fexperimental-library -fprebuilt-module-path=clang.cache

ifneq '' '$(findstring clang,$(COMPILER))'
	CXXFLAGS += $(CLANG_CXXFLAGS)
	CLANG_CACHE=clang.cache
	CLANG_STDMODULE=clang.cache/std.pcm
	CLANG_STD_MANIFEST=$(shell $(CXX) -stdlib=libc++ -print-file-name=libc++.modules.json)
	CLANG_STDMODULE_CPPM=$(realpath $(dir $(CLANG_STD_MANIFEST))/$(shell cat "$(CLANG_STD_MANIFEST)" | jq -r '.modules[] | select(.["logical-name"] == "std") | .["source-path"]'))
	STDMODULES_FILES=$(CLANG_STDMODULE)
	CXXFLAGS_MODULE=-fmodule-output=clang.cache/$(basename $(notdir $@)).pcm # interpreted when called
else ifneq '' '$(findstring g++,$(COMPILER))'
	CXXFLAGS += $(GCC_CXXFLAGS)
	EXTRA_CXXFLAGS=-flto-partition=none 
	STDMODULES=iostream vector fstream tuple random algorithm limits cassert map set list string numeric sstream cmath array filesystem thread syncstream latch barrier optional atomic iterator
	STDMODULES_FAKEFILES=$(foreach s,$(STDMODULES),makefile.cache/$(s))
	STDMODULES_FILES=$(STDMODULES_FAKEFILES)
else
	$(warning Only GCC and Clang are supported!)
	exit 1
endif

.PHONY : clean

LLCHAIN_OBJS=src/utils.o src/MinSegmentTree.o src/algo.o ext/grid_to_bmp.o ext/mummer_essaMEM_wrapper.o ext/kseq.o ext/chainx.o ext/cli11.o ext/htslib_wrapper.o
llchain : src/llchain.cpp $(LLCHAIN_OBJS) $(STDMODULES_FILES)
	$(CXX) -I ext/concurrentqueue $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(HTSLIB_CFLAGS) -DVERSION="\"$(VERSION)\"" -DCOMPILER="\"$(COMPILER)\"" $< $(LLCHAIN_OBJS) -o llchain $(LDFLAGS)

src/utils.o : src/utils.cppm ext/grid_to_bmp.o $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(CXXFLAGS_MODULE)
src/algo.o : src/algo.cppm src/utils.o src/MinSegmentTree.o $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(CXXFLAGS_MODULE)
src/MinSegmentTree.o : src/MinSegmentTree.cppm src/utils.o $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(CXXFLAGS_MODULE)
CLI11_INCL=ext/CLI11/include
ext/cli11.o : ext/CLI11/src/modules/CLI11.cppm $(shell find $(CLI11_INCL) -name *.hpp)
	$(CXX) $(CXXFLAGS) -I $(CLI11_INCL) -c $< -o $@ $(CXXFLAGS_MODULE)
ext/grid_to_bmp.o : ext/grid_to_bmp.cppm $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(CXXFLAGS_MODULE)
ext/kseq.o : ext/kseq.cppm $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(CXXFLAGS_MODULE)
MUMMER_CPPS=ext/mummer/src/essaMEM/
MUMMER_INCL=ext/mummer/include/
MUMMER_CPP_FLAGS=-DNDEBUG
ext/mummer_essaMEM_wrapper.o : ext/mummer_essaMEM_wrapper.cppm $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) $(MUMMER_CPP_FLAGS) -I $(MUMMER_INCL) -I $(MUMMER_CPPS) -c $< -o $@ $(CXXFLAGS_MODULE)

HTSLIB_INCL=$(ROOT_DIR)ext/htslib
ext/htslib_wrapper.o : ext/htslib_wrapper.cppm $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -I $(HTSLIB_INCL) -c $< -o $@ $(CXXFLAGS_MODULE)
ext/chainx.o : ext/chainx.cppm $(STDMODULES_FILES)
	$(CXX) $(CXXFLAGS) -I $(MUMMER_INCL) -I $(MUMMER_CPPS) -c $< -o $@ $(CXXFLAGS_MODULE)
$(STDMODULES_FAKEFILES) :
	mkdir -p makefile.cache
	$(CXX) $(CXXFLAGS) -xc++-system-header $(notdir $@)
	touch $@
$(CLANG_STDMODULE) :
	mkdir -p clang.cache
	$(CXX) $(CXXFLAGS) -Wno-reserved-module-identifier -c $(CLANG_STDMODULE_CPPM) -o $(basename $@).o -fmodule-output=$@

clean :
	-rm -f $(LLCHAIN_OBJS)
	-rm -Rf gcm.cache clang.cache makefile.cache
