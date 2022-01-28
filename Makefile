MAIN_SOURCE=./src/main.cpp
TEST_SOURCE=./tests/test.cpp

BRANCH := $(shell git rev-parse --abbrev-ref HEAD)

TARGET=./bin/${BRANCH}_order.out
TEST_TARGET=./bin/${BRANCH}_test.out

MKDIR_C=mkdir -p
DIRS=bin
METIS_INCLUDE=lib/include/
METIS_LIB=lib/lib/

GCC=g++
GCC_FLAGS= -std=c++14 -fopenmp -O3 -I${METIS_INCLUDE} -L${METIS_LIB} -lmetis
GCC_DEBUG_FLAGS= -std=c++14 -fopenmp -g 

all : $(TARGET)
test: $(TEST_TARGET)

$(TARGET): dir
	$(GCC) -o $@ $(MAIN_SOURCE) $(GCC_FLAGS) 
$(TEST_TARGET): dir
	$(GCC) $(GCC_DEBUG_FLAGS) -o $@ $(TEST_SOURCE)

dir:
	$(MKDIR_C) $(DIRS)

.PHONY: clean

clean:
	@rm -rfv $(TARGET) $(TEST_TARGET)

