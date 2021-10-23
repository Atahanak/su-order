MAIN_SOURCE=./src/main.cpp

BRANCH := $(shell git rev-parse --abbrev-ref HEAD)

TARGET=./bin/${BRANCH}_order.out

MKDIR_C=mkdir -p
DIRS=bin

all: $(TARGET) .intr

GCC=g++
GCC_FLAGS= -std=c++14 -fopenmp -O3

all : $(TARGET)

$(TARGET): dir
	$(GCC) $(GCC_FLAGS) -o $@ $(MAIN_SOURCE)

dir:
	$(MKDIR_C) $(DIRS)

.intr:
	rm -rf $(OBJECTS) $(ARGPARSE_MAIN) $(NONARGPARSE_MAIN)

.PHONY: clean

clean:
	@rm -rfv $(TARGET)

