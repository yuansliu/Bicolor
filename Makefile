CC = g++
GATB=../gatb-core-1.2.3-bin-Linux
# CFLAG = -O3 -Wall -I$(GATB)/include/ -L$(GATB)/lib/ -lgatbcore -lhdf5 -ldl -lz -lpthread -std=c++0x
CFLAG = -Wno-deprecated -I$(GATB)/include -L$(GATB)/lib -lgatbcore -lhdf5 -ldl -lz -lpthread -std=c++0x -O3

all:
	$(CC) iterative.cpp $(CFLAG) -o iterative
	$(CC) combine.cpp $(CFLAG) -o combine

iterative: iterative.cpp
	$(CC) iterative.cpp $(CFLAG) -o iterative

combine: combine.cpp
	$(CC) combine.cpp $(CFLAG) -o combine
