CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I../../algorithm -I.
LDFLAGS  := -fopenmp

ifdef DEBUG
  CXXFLAGS += -g
endif

LIB := ../../algorithm/librank.a

all: cmm_perf

cmm_perf: cmm_perf.o motifs_search_test.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $^

$(LIB):
	$(MAKE) -C ../../algorithm

clean:
	rm -f *.o cmm_perf
