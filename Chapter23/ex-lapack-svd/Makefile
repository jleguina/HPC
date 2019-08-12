CC=g++
CXXFLAGS=-Wall -O0 -pedantic
LDLIBS=-llapack -lblas -ljpeg

default: svd

blas_l1: svd.o

.PHONY: clean

clean:
	rm -f svd *.o output.jpg

