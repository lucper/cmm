CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I../.. -I.
LDFLAGS  := -fopenmp

ifdef DEBUG
  CXXFLAGS += -g
else
  CXXFLAGS += -DNDEBUG
endif

ALG_DIR   := ../..
CORE_OBJS := $(ALG_DIR)/rank_table.o $(ALG_DIR)/data_import.o \
             $(ALG_DIR)/libsais64.o $(ALG_DIR)/libsais.o

BIN_DIR := ../bin
BINS    := $(BIN_DIR)/cmm_test $(BIN_DIR)/cmm_eval $(BIN_DIR)/cmm_dedup

all: $(BINS)

$(BIN_DIR)/cmm_test: cmm_test.o motifs_search_test.o $(CORE_OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/cmm_eval: cmm_eval.o cmm_acc.o $(ALG_DIR)/data_import.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/cmm_dedup: cmm_dedup.o cmm_acc.o $(ALG_DIR)/data_import.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $<

## Header dependencies.
cmm_acc.o:            cmm_acc.cpp   cmm_acc.hpp
cmm_eval.o:           cmm_eval.cpp  cmm_acc.hpp
cmm_dedup.o:          cmm_dedup.cpp cmm_acc.hpp
cmm_test.o:           cmm_test.cpp  motifs_search_test.hpp
motifs_search_test.o: motifs_search_test.cpp motifs_search_test.hpp

$(CORE_OBJS):
	$(MAKE) -C $(ALG_DIR)

clean:
	rm -f *.o
	rm -f $(BINS)
