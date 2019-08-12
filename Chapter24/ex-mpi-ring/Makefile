CC=mpicxx
CXX=mpicxx
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=
TARGET=ring

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean run

# An example of using the code
run: $(TARGET)
	mpiexec -np 8 ./$(TARGET)

clean:
	rm -f $(TARGET) *.o
