CC=mpicxx
CXX=mpicxx
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=
TARGET=mpi_bcast

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean run

# An example of using the code
run: $(TARGET)
	mpiexec -np 2 ./$(TARGET)

clean:
	rm -f $(TARGET) *.o
