CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-lboost_program_options
TARGET=sort_boost

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean run

# An example of using the code
run: $(TARGET)
	./$(TARGET) --size 15 --min 10 --max 50 --descending

clean:
	rm -f $(TARGET) *.o
