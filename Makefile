# CHIC - Makefile (C++ library build preserved)
CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -funroll-loops -DNDEBUG -fPIC
INCLUDES = -I/usr/include/eigen3 -I. -Isrc

# C++ sources
LIB_SOURCE = src/CHIC.cpp src/CHIC_EARTH.cpp
HEADER = src/CHIC.h

# C++ build artifacts
BUILD_DIR = build
LIB_OBJ = $(BUILD_DIR)/CHIC.o $(BUILD_DIR)/CHIC_EARTH.o
LIB_TARGET = $(BUILD_DIR)/libchic.a

.PHONY: all lib clean dirs

all: lib

dirs:
	@mkdir -p $(BUILD_DIR)

# Build static library (C++-only)
lib: dirs $(LIB_TARGET)

$(LIB_TARGET): $(LIB_OBJ)
	@echo "Creating static library: $@"
	ar rcs $@ $^

$(BUILD_DIR)/%.o: src/%.cpp $(HEADER) | dirs
	@echo "Compiling library object: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.a
	@echo "Cleaned C++ build artifacts"
