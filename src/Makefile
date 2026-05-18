CXX      := g++
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native -I.
CFLAGS   := -O3 -march=native -I.
LDFLAGS  := -fopenmp
LIBSAIS_OBJ := libsais64.o libsais.o

ifdef DEBUG
  CXXFLAGS += -g
endif

CORE_SRCS := motifs_search.cpp rank_table.cpp data_import.cpp
CORE_OBJS := $(CORE_SRCS:.cpp=.o)

all: libcmm.a cmm

libcmm.a: $(CORE_OBJS) $(LIBSAIS_OBJ)
	ar rcs $@ $^

cmm: main.o libcmm.a
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $^

clean:
	rm -f *.o libcmm.a cmm
