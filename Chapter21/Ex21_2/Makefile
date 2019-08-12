CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-lfftw3
TARGET=fftw

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

generate: generate.o

.PHONY: clean signal spectrum

signal: generate
	./generate > signal.txt

spectrum: fftw signal
	./fftw > spectrum.txt

clean:
	rm -f $(TARGET) generate *.o *.txt
