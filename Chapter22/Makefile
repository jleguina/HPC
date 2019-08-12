CC=g++
CXXFLAGS=-Wall -O0 -pedantic
LDLIBS=-lblas
TARGET=blas_l1

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean

clean:
	rm -f $(TARGET) *.o
