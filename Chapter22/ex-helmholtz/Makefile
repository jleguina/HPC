CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-llapack -lblas
TARGET=cg_helmholtz

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean doc

doc:
	doxygen Doxyfile

clean:
	rm -rf $(TARGET) *.o html latex
