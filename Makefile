CC=g++
CFLAGS=-Wall -fopenmp -pthread -g
CFLAGS2=-c
LDFLAGS=
SOURCES=
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=SPCE

all:	o2

o3:
	rm -rf $(EXECUTABLE)	
	$(CC) $(CFLAGS) -O3 *.cpp -o $(EXECUTABLE)

o2:
	rm -rf $(EXECUTABLE)	
	$(CC) $(CFLAGS) -O2 *.cpp -o $(EXECUTABLE)
	
o1:
	rm -rf $(EXECUTABLE)	
	$(CC) $(CFLAGS) -O1 *.cpp -o $(EXECUTABLE)

o0:
	rm -rf $(EXECUTABLE)	
	$(CC) $(CFLAGS) *.cpp -o $(EXECUTABLE)


verbose: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
	rm -rf *.o *~

.cpp.o:
	$(CC) $(CFLAGS2) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *~
