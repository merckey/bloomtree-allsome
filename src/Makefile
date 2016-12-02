CXXFLAGS+=-Wall -std=c++11 -I$(HOME)/include/roaring
LDFLAGS+=-lroaring -lsdsl -ldivsufsort -ldivsufsort64 -ljellyfish-2.0 -lz -lpthread #-static

# these shouldnt be necessary, but it helps users who have not set $HOME/[include/lib] in compiler flags
CXXFLAGS+=-I$(HOME)/include -I$(HOME)/include/jellyfish-2.2.6/ 
LDFLAGS+=-L$(HOME)/lib -L$(HOME)/lib/jellyfish-2.2.6

CPP_FILES := $(wildcard *.cc)
OBJ_FILES := $(addprefix ./,$(notdir $(CPP_FILES:.cc=.o)))

debug: CXXFLAGS += -DDEBUG -g -O1
all:   CXXFLAGS += -DNDEBUG -O3 

all: bt

bt: $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

debug: main.o Build.o Query.o Kmers.o BloomTree.o BF.o util.o Count.o
	$(CXX) -o debug $^ $(LDFLAGS)

clean:
	rm -f bt *.o
# DO NOT DELETE
