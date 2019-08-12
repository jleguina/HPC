CXX=mpicxx
CXXFLAGS= -Wall -O0 
DIFF=0 0 10 10 1 20001 20001 4000 0 0 0 1 4 4
ADVX=0 0 10 10 1 20001 20001 4000  1 0 0 0 4 4
ADVY=0 0 10 10 1 20001 20001 4000  0 1 0 0 4 4
BURG=0 0 10 10 1 20001 20001 4000  1.00 0.50 1.00 0.02 4 4
TARGET=compile
TARGET_DIFF=diff
TARGET_AVDX=advx
TARGET_ADVY=advy
TARGET_BURG=burg

default: $(TARGET)
all: $(TARGET_DIFF) $(TARGET_AVDX) $(TARGET_AVDY) $(TARGET_BURG) 

Model.o: Model.cpp Model.h 
	$(CXX) $(CXXFLAGS) -o Model.o -c Model.cpp

Burgers.o:  Burgers.cpp Burgers.h
	$(CXX) $(CXXFLAGS) -o Burgers.o -c Burgers.cpp
	
BurgersSingle.o:  BurgersSingle.cpp BurgersSingle.h
	$(CXX) $(CXXFLAGS) -o BurgersSingle.o -c BurgersSingle.cpp

burguers.o: burguers.cpp
	$(CXX) $(CXXFLAGS) -o burguers.o -c burguers.cpp

$(TARGET): burguers.o Burgers.o Model.o 
	$(CXX) $(CXXFLAGS) -o prog burguers.o Model.o Burgers.o BurgersSingle.o

$(TARGET_DIFF): $(TARGET)
	mpiexec -np 16 ./prog $(DIFF)

$(TARGET_ADVX): $(TARGET)
	mpiexec -np 16 ./prog $(ADVX)

$(TARGET_ADVY): $(TARGET)
	mpiexec -np 16 ./prog $(ADVY)

$(TARGET_BURG): $(TARGET)
	mpiexec -np 16 ./prog $(BURG)

.PHONY: clean

clean:
	rm -f $(TARGET) *.o
