SRC_DIR=src
INCLUDE_DIR=include
OBJ_DIR=obj

CXX=g++
PARALLEL_FLAGS=-fopenmp -pthread
DEBUG_FLAG=-g
CXXFLAGS=-c -Wall -I $(INCLUDE_DIR)
LDFLAGS=
SOURCES= $(shell ls $(SRC_DIR)/*.cpp)
OBJECTS=$(addprefix $(OBJ_DIR)/,$(notdir $(SOURCES:.cpp=.o)))
EXECUTABLE=SPCE

all:	o2

o3: CXXFLAGS += -O3
o3: verbose

o2: CXXFLAGS += -O2
o2: verbose

o1: CXXFLAGS += -O1
o1: verbose

o0: verbose

verbose: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

remove:
	rm -rf $(EXECUTABLE)
	
clean: remove
	cd obj && rm -rf * && cd ..
	