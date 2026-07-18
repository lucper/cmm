CXX      := g++
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native
CFLAGS   := -O3 -march=native
LDFLAGS  := -fopenmp

ifdef DEBUG
  CXXFLAGS += -g
endif

OBJS := main.o motifs_search.o rank_table.o data_import.o libsais64.o libsais.o

all: cmm

cmm: $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) cmm
