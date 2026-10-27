CXX      := g++
MPICXX   := mpicxx
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native
CFLAGS   := -O3 -march=native
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP
INCLUDES := -Isrc -Iexternal
# The SKIP_MPICXX macros leave out the C++ bindings of MPI (removed in MPI 3.0), which we do not use.
MPIFLAGS := -DCMM_MPI -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX

ifdef DEBUG
  CXXFLAGS += -g
else
  CXXFLAGS += -DNDEBUG
endif

# The objects are built in build/, from the sources of src/ and external/.
CORE_OBJS := build/exact_algorithm.o build/motif_pair_searcher.o build/occurrence_matrix.o \
             build/csr.o build/progress_bar.o build/rank_table.o build/esa.o build/data_import.o \
             build/graph.o build/libsais64.o build/libsais.o
OBJS      := build/main.o build/options_parsing.o $(CORE_OBJS)

all: cmm

cmm: $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

cmm_mpi: build/main_mpi.o build/options_parsing.o $(CORE_OBJS)
	$(MPICXX) $(LDFLAGS) -o $@ $^

build/main_mpi.o: src/main.cpp | build
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) $(DEPFLAGS) $(MPIFLAGS) -fopenmp -c -o $@ $<

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPFLAGS) -fopenmp -c -o $@ $<

build/%.o: external/%.c | build
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

build:
	mkdir -p $@

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d) build/main_mpi.d

# Runs the tests of cmm and cmm_test (see test/).
.PHONY: test
test:
	$(MAKE) -C test

clean:
	rm -rf build cmm cmm_mpi
