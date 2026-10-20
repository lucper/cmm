CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O3 -march=native -I..
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP

ALG_DIR   := ..
CORE_OBJS := $(ALG_DIR)/motifs_search.o $(ALG_DIR)/rank_table.o $(ALG_DIR)/data_import.o \
             $(ALG_DIR)/libsais64.o $(ALG_DIR)/libsais.o
OBJS      := main.o test_search.o test_errors.o test_programs.o

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

# The objects of the root directory and cmm are brought up to date by the root Makefile.
.PHONY: core
core:
	$(MAKE) -s -C $(ALG_DIR) $(notdir $(CORE_OBJS)) cmm

$(CORE_OBJS): core ;

# Writes the expected outputs of the test cases (see generate_expected in test_search.cpp).
.PHONY: expected
expected: tests
	./tests --generate

clean:
	rm -f $(OBJS) $(OBJS:.o=.d) tests
