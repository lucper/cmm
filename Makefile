CXX      := g++
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native
CFLAGS   := -O3 -march=native
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP

ifdef DEBUG
  CXXFLAGS += -g
endif

OBJS := main.o motifs_search.o rank_table.o data_import.o libsais64.o libsais.o

all: cmm

cmm: $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -fopenmp -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d)

# Runs the tests of cmm and cmm_test (see test/).
.PHONY: test
test:
	$(MAKE) -C test

clean:
	rm -f $(OBJS) $(OBJS:.o=.d) cmm
