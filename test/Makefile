CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O3 -march=native -I../src -I../external
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP

ALG_DIR   := ..
BUILD_DIR := $(ALG_DIR)/build
CORE_OBJS := $(BUILD_DIR)/motifs_search.o $(BUILD_DIR)/motif_pair_searcher.o $(BUILD_DIR)/motif_nodes.o \
             $(BUILD_DIR)/csr.o $(BUILD_DIR)/progress_bar.o $(BUILD_DIR)/rank_table.o $(BUILD_DIR)/data_import.o \
             $(BUILD_DIR)/libsais64.o $(BUILD_DIR)/libsais.o
OBJS      := main.o test_search.o test_errors.o test_programs.o test_csr.o test_motif_nodes.o \
             test_motif_pair_searcher.o

# Runs the tests, of the search and of the programs cmm and cmm_test.
.PHONY: check
check: tests
	$(MAKE) -s -C $(ALG_DIR)/experiments/src
	./tests

tests: $(OBJS) $(CORE_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -fopenmp -c -o $@ $<

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d)

# The objects of build/ and cmm are brought up to date by the root Makefile.
.PHONY: core
core:
	$(MAKE) -s -C $(ALG_DIR) $(patsubst $(ALG_DIR)/%,%,$(CORE_OBJS)) cmm

$(CORE_OBJS): core ;

# Writes the expected outputs of the test cases (see generate_expected in test_search.cpp).
.PHONY: expected
expected: tests
	./tests --generate

clean:
	rm -f $(OBJS) $(OBJS:.o=.d) tests
