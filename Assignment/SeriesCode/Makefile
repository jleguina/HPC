CXX=g++
CXXFLAGS=-std=c++11 -Wall -O0 
DIFF=0 0 10 10 1 7 5 6 0 0 0 1
ADVX=0 0 10 10 1 7 5 6  1 0 0 0
ADVY=0 0 10 10 1 7 5 6  0 1 0 0
BURG=0 0 10 10 1 7 5 6  1.00 0.50 1.00 0.02
TARGET=compile
TARGET_DIFF=diff
TARGET_AVDX=advx
TARGET_ADVY=advy
TARGET_BURG=burg

default: $(TARGET)
all: $(TARGET_DIFF) $(TARGET_AVDX) $(TARGET_AVDY) $(TARGET_BURG) 

Model.o: Model.cpp Model.h 
	$(CXX) $(CXXFLAGS) -o Model.o -c Model.cpp

burgers.p:  Burgers.cpp Burgers.h
	$(CXX) $(CXXFLAGS) -o Burgers.o -c Burgers.cpp

burguers.o: burguers.cpp
	$(CXX) $(CXXFLAGS) -o burguers.o -c burguers.cpp

$(TARGET): burguers.o Burgers.o Model.o 
	$(CXX) $(CXXFLAGS) -o Assignment burguers.o Burgers.o Model.o 

$(TARGET_DIFF): $(TARGET)
	./Assignment $(DIFF)

$(TARGET_ADVX): $(TARGET)
	./Assignment $(ADVX)

$(TARGET_ADVY): $(TARGET)
	./Assignment $(ADVY)

$(TARGET_BURG): $(TARGET)
	./Assignment $(BURG)

.PHONY: clean

clean:
	rm -f $(TARGET) *.o
