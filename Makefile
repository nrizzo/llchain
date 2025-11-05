FLAGS=-std=c++20 -O0 -g -fmodules
CFLAGS=-std=c++20 -O0 -g
LDFLAGS=-lz
STDMODULES=iostream vector fstream stdio.h tuple random algorithm limits cassert map set list string ctype.h string.h stdlib.h zlib.h
STDMODULES_FAKEFILES=$(foreach s,$(STDMODULES),gcm.cache/$(s))
.PHONY : clean

clc-viz : src/clc-viz.cpp src/utils.o src/MinSegmentTree.o src/algo.o ext/grid_to_bmp.o src/command-line-parsing/cmdline.o ext/kseq.o ext/mummer_essaMEM_wrapper.o $(STDMODULES_FAKEFILES)
	g++ -I src $(FLAGS) $^ -o clc-viz $(LDFLAGS)

src/command-line-parsing/cmdline.o : src/command-line-parsing/cmdline.h src/command-line-parsing/cmdline.c
	g++ $(CFLAGS) -c $^ -o $@
src/command-line-parsing/cmdline.c src/command-line-parsing/cmdline.h : src/command-line-parsing/config.ggo
	gengetopt \
		--input=./src/command-line-parsing/config.ggo \
		--output-dir=./src/command-line-parsing/ \
		--unnamed-opts

src/utils.o : src/utils.cppm ext/grid_to_bmp.o $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
src/algo.o : src/algo.cppm src/utils.o src/MinSegmentTree.o $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
src/MinSegmentTree.o : src/MinSegmentTree.cppm src/utils.o $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
ext/grid_to_bmp.o : ext/grid_to_bmp.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
ext/kseq.o : ext/kseq.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -c $< -o $@
MUMMER_CPPS=ext/mummer/src/essaMEM/
MUMMER_INCL=ext/mummer/include/
ext/mummer_essaMEM_wrapper.o : ext/mummer_essaMEM_wrapper.cppm $(STDMODULES_FAKEFILES)
	g++ $(FLAGS) -I $(MUMMER_INCL) -I $(MUMMER_CPPS) -c $< -o $@
$(STDMODULES_FAKEFILES) :
	g++ $(FLAGS) -xc++-system-header $(notdir $@)
	touch $@

clean :
	-rm -f src/utils.o src/algo.o src/MinSegmentTree.o src/command-line-parsing/cmdline.o ext/grid_to_bmp.o
	-rm -Rf gcm.cache
