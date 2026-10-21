CXX      := g++
MPICXX   := mpicxx
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native
CFLAGS   := -O3 -march=native
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP
# The SKIP_MPICXX macros leave out the C++ bindings of MPI (removed in MPI 3.0), which we do not use.
MPIFLAGS := -DCMM_MPI -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX

ifdef DEBUG
  CXXFLAGS += -g
endif

CORE_OBJS := motifs_search.o motif_pair_searcher.o motif_nodes.o csr.o progress_bar.o rank_table.o data_import.o libsais64.o libsais.o
OBJS      := main.o options_parsing.o $(CORE_OBJS)

all: cmm

cmm: $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

cmm_mpi: main_mpi.o options_parsing.o $(CORE_OBJS)
	$(MPICXX) $(LDFLAGS) -o $@ $^

main_mpi.o: main.cpp
	$(MPICXX) $(CXXFLAGS) $(DEPFLAGS) $(MPIFLAGS) -fopenmp -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -fopenmp -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d) main_mpi.d

# Runs the tests of cmm and cmm_test (see test/).
.PHONY: test
test:
	$(MAKE) -C test

clean:
	rm -f $(OBJS) $(OBJS:.o=.d) cmm main_mpi.o main_mpi.d cmm_mpi
