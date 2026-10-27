CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I../../src -I../../external -I.
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP

ifdef DEBUG
  CXXFLAGS += -g
else
  CXXFLAGS += -DNDEBUG
endif

ALG_DIR   := ../..
BUILD_DIR := $(ALG_DIR)/build
CORE_OBJS := $(BUILD_DIR)/exact_algorithm.o $(BUILD_DIR)/motif_pair_searcher.o \
             $(BUILD_DIR)/occurrence_matrix.o $(BUILD_DIR)/csr.o $(BUILD_DIR)/progress_bar.o \
             $(BUILD_DIR)/rank_table.o $(BUILD_DIR)/esa.o $(BUILD_DIR)/data_import.o \
             $(BUILD_DIR)/graph.o $(BUILD_DIR)/libsais64.o $(BUILD_DIR)/libsais.o

BIN_DIR := ../bin
BINS    := $(BIN_DIR)/cmm_test $(BIN_DIR)/cmm_eval $(BIN_DIR)/cmm_dedup

all: $(BINS)

$(BIN_DIR)/cmm_test: cmm_test.o $(CORE_OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/cmm_eval: cmm_eval.o cmm_acc.o $(BUILD_DIR)/data_import.o $(BUILD_DIR)/graph.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/cmm_dedup: cmm_dedup.o cmm_acc.o $(BUILD_DIR)/data_import.o $(BUILD_DIR)/graph.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

OBJS := cmm_test.o cmm_eval.o cmm_dedup.o cmm_acc.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -fopenmp -c -o $@ $<

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d)

# The objects of build/ are brought up to date by the root Makefile, which knows their sources
# and headers, every time these binaries are built.
.PHONY: core
core:
	$(MAKE) -C $(ALG_DIR) $(patsubst $(ALG_DIR)/%,%,$(CORE_OBJS))

$(CORE_OBJS): core ;

clean:
	rm -f $(OBJS) $(OBJS:.o=.d)
	rm -f $(BINS)
