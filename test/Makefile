CXX      := g++
MPICXX   := mpicxx
CXXFLAGS := -std=c++17 -Wall -Wextra -O3 -march=native -I../src -I../external -I../external/doctest
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP

ALG_DIR   := ..
BUILD_DIR := $(ALG_DIR)/build
CORE_OBJS := $(BUILD_DIR)/exact_algorithm.o $(BUILD_DIR)/motif_pair_searcher.o \
             $(BUILD_DIR)/occurrence_matrix.o $(BUILD_DIR)/csr.o $(BUILD_DIR)/progress_bar.o \
             $(BUILD_DIR)/rank_table.o $(BUILD_DIR)/esa.o $(BUILD_DIR)/data_import.o \
             $(BUILD_DIR)/graph.o $(BUILD_DIR)/libsais64.o $(BUILD_DIR)/libsais.o
OBJS      := main.o test_search.o test_errors.o test_programs.o test_csr.o \
             test_occurrence_matrix.o test_motif_pair_searcher.o test_graph.o
# The tests of the MPI version (mpicxx; see check_mpi).
MPI_CORE_OBJS := $(BUILD_DIR)/exact_algorithm_mpi.o
MPI_OBJS      := main_mpi.o test_exact_algorithm_mpi.o
# The SKIP_MPICXX macros leave out the C++ bindings of MPI (removed in MPI 3.0), which we do not use.
MPIFLAGS      := -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
# The processes of the MPI tests; --oversubscribe allows more than the cores of the machine.
MPI_NP        := 4
MPI_TIMEOUT   := 120

# Runs the tests, of the search and of the programs cmm and cmm_test.
.PHONY: check
check: tests
	$(MAKE) -s -C $(ALG_DIR)/experiments/src
	./tests

tests: $(OBJS) $(CORE_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

# Runs the tests of the MPI version on MPI_NP processes. A process that fails alone (e.g. an
# exception) leaves the others waiting in a collective, so mpirun stops them after MPI_TIMEOUT seconds.
.PHONY: check_mpi
check_mpi: tests_mpi
	mpirun --oversubscribe --timeout $(MPI_TIMEOUT) -np $(MPI_NP) ./tests_mpi

tests_mpi: $(MPI_OBJS) $(MPI_CORE_OBJS) $(CORE_OBJS)
	$(MPICXX) $(LDFLAGS) -o $@ $^

$(MPI_OBJS): %.o: %.cpp
	$(MPICXX) $(CXXFLAGS) $(DEPFLAGS) $(MPIFLAGS) -fopenmp -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -fopenmp -c -o $@ $<

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d) $(MPI_OBJS:.o=.d)

# The objects of build/ and cmm are brought up to date by the root Makefile.
.PHONY: core
core:
	$(MAKE) -s -C $(ALG_DIR) $(patsubst $(ALG_DIR)/%,%,$(CORE_OBJS)) cmm

$(CORE_OBJS): core ;

.PHONY: core_mpi
core_mpi:
	$(MAKE) -s -C $(ALG_DIR) $(patsubst $(ALG_DIR)/%,%,$(MPI_CORE_OBJS))

$(MPI_CORE_OBJS): core_mpi ;

# Writes the expected outputs of the test cases (see generate_expected in test_search.cpp).
.PHONY: expected
expected: tests
	./tests --generate

clean:
	rm -f $(OBJS) $(OBJS:.o=.d) tests $(MPI_OBJS) $(MPI_OBJS:.o=.d) tests_mpi
