CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I../../algorithm -I.
LDFLAGS  := -fopenmp

ifdef DEBUG
  CXXFLAGS += -g
else
  CXXFLAGS += -DNDEBUG
endif

ALG_DIR  := ../../algorithm
ALG_OBJS := $(ALG_DIR)/data_import.o
OBJS     := cmm_eval.o

all: cmm_eval

cmm_eval: $(OBJS) $(ALG_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $<

$(ALG_OBJS):
	$(MAKE) -C $(ALG_DIR)

clean:
	rm -f $(OBJS) cmm_perf
