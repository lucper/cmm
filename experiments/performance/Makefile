CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I../../algorithm
LDFLAGS  := -fopenmp

ifdef DEBUG
  CXXFLAGS += -g
endif

LIB := ../../algorithm/libcmm.a

all: cmm_perf

cmm_perf: cmm_perf.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^

cmm_perf.o: cmm_perf.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $^

$(LIB):
	$(MAKE) -C ../../src

clean:
	rm -f cmm_perf.o cmm_perf
