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

o0: CXXFLAGS += -g
o0: verbose

verbose: init $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

mbar:
	$(CXX) -lm $(SRC_DIR)/mbar.c -o mbar

remove:
	rm -rf $(EXECUTABLE) mbar

clean: remove
	rm -rf *~
	rm -rf $(OBJ_DIR)/*
	rm -rf $(SRC_DIR)/*~
	rm -rf $(INCLUDE_DIR)/*~

init:
	mkdir -p $(OBJ_DIR)