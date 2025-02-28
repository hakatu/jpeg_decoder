# Define the compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall

# Include paths (adjust if headers are in a different directory)
INCLUDE_PATH = .
CXXFLAGS += -I$(INCLUDE_PATH)

# Target executable
TARGET = idct_test

# Source file
SRC = main.cpp

# Object file
OBJ = $(SRC:.cpp=.o)

# Default target: build the executable
all: $(TARGET)

# Compile main.cpp into main.o
$(OBJ): $(SRC) idct_fast.h idct_aan.h
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $(OBJ)

# Link the object file to create the executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Clean target: remove object files and executable
clean:
	rm -f $(OBJ) $(TARGET)