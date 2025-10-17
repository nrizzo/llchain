DEPS=src/command-line-parsing/cmdline.h src/command-line-parsing/cmdline.c
LIB_DIRS=ext

clc-viz : src/clc-viz.cpp $(DEPS)
	g++ -std=c++20 -O0 -g \
	-I $(LIB_DIRS) \
	src/clc-viz.cpp src/command-line-parsing/cmdline.c \
	-o clc-viz

src/command-line-parsing/cmdline%c src/command-line-parsing/cmdline%h : src/command-line-parsing/config.ggo
	gengetopt \
		--input=./src/command-line-parsing/config.ggo \
		--output-dir=./src/command-line-parsing/ \
		--unnamed-opts
