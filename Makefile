DEPS=ext/grid_to_bmp.hpp command-line-parsing/cmdline.h command-line-parsing/cmdline.c

clc-viz : clc-viz.cpp $(DEPS)
	g++ -O0 -g \
	clc-viz.cpp command-line-parsing/cmdline.c \
	-o clc-viz

command-line-parsing/cmdline%c command-line-parsing/cmdline%h : command-line-parsing/config.ggo
	gengetopt \
		--input=./command-line-parsing/config.ggo \
		--output-dir=./command-line-parsing/ \
		--unnamed-opts
