CXX      := g++
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -fno-omit-frame-pointer -O3 -march=native
CFLAGS   := -O3 -march=native
LDFLAGS  := -fopenmp
LIBSAIS_OBJ := libsais64.o libsais.o

ifdef DEBUG
  CXXFLAGS += -g
endif

CORE_SRCS := rank_table.cpp data_import.cpp
CORE_OBJS := $(CORE_SRCS:.cpp=.o)

all: librank.a cmm

librank.a: $(CORE_OBJS) $(LIBSAIS_OBJ)
	ar rcs $@ $^

cmm: main.o motifs_search.cpp librank.a
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $^

clean:
	rm -f *.o librank.a cmm
