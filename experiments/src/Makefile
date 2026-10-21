CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I../.. -I.
LDFLAGS  := -fopenmp
DEPFLAGS := -MMD -MP

ifdef DEBUG
  CXXFLAGS += -g
else
  CXXFLAGS += -DNDEBUG
endif

ALG_DIR   := ../..
CORE_OBJS := $(ALG_DIR)/motifs_search.o $(ALG_DIR)/motif_pair_searcher.o $(ALG_DIR)/motif_nodes.o $(ALG_DIR)/csr.o \
             $(ALG_DIR)/progress_bar.o $(ALG_DIR)/rank_table.o $(ALG_DIR)/data_import.o $(ALG_DIR)/libsais64.o \
             $(ALG_DIR)/libsais.o

BIN_DIR := ../bin
BINS    := $(BIN_DIR)/cmm_test $(BIN_DIR)/cmm_eval $(BIN_DIR)/cmm_dedup

all: $(BINS)

$(BIN_DIR)/cmm_test: cmm_test.o $(CORE_OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/cmm_eval: cmm_eval.o cmm_acc.o $(ALG_DIR)/data_import.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/cmm_dedup: cmm_dedup.o cmm_acc.o $(ALG_DIR)/data_import.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

OBJS := cmm_test.o cmm_eval.o cmm_dedup.o cmm_acc.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -fopenmp -c -o $@ $<

# The headers each object depends on, written by the compiler (-MMD) when it is built.
-include $(OBJS:.o=.d)

# The objects of the root directory are brought up to date by its own Makefile, which knows
# their sources and headers, every time these binaries are built.
.PHONY: core
core:
	$(MAKE) -C $(ALG_DIR) $(notdir $(CORE_OBJS))

$(CORE_OBJS): core ;

clean:
	rm -f $(OBJS) $(OBJS:.o=.d)
	rm -f $(BINS)
