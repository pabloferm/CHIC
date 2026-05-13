# CHIC - Makefile
CXX      = g++
CXXFLAGS = -std=c++26 -O3 -march=native -ffast-math -funroll-loops \
           -fPIC -DNDEBUG \
           -fomit-frame-pointer \
           -fno-trapping-math \
           -fassociative-math \
           -freciprocal-math \
           -ffinite-math-only
INCLUDES = -I/usr/include/eigen3 -I. -Isrc

# C++ sources
LIB_SOURCE = src/CHIC.cpp
HEADER     = src/CHIC.h

# Build artifacts
BUILD_DIR  = build
LIB_OBJ    = $(BUILD_DIR)/CHIC.o
LIB_SOURCE = src/CHIC.cpp src/CHIC_EARTH.cpp
HEADER     = src/CHIC.h src/CHIC_EARTH.h

# C++ build
BUILD_DIR = build
LIB_OBJ   = $(BUILD_DIR)/CHIC.o $(BUILD_DIR)/CHIC_EARTH.o
LIB_TARGET = $(BUILD_DIR)/libchic.a

# Examples
EXAMPLE_DIR  = examples/cpp
EXAMPLE_BINS = $(EXAMPLE_DIR)/chic_benchmark $(EXAMPLE_DIR)/example

# Prevent make from deleting intermediate files (library, objects)
.SECONDARY:

EXAMPLE_BINS = $(EXAMPLE_DIR)/example $(EXAMPLE_DIR)/test_chicearth $(EXAMPLE_DIR)/speed

.SECONDARY:
.PHONY: all lib examples clean dirs

all: lib

dirs:
	@mkdir -p $(BUILD_DIR)

# Build static library
lib: dirs $(LIB_TARGET)

$(LIB_TARGET): $(LIB_OBJ)
	@echo "Creating static library: $@"
	ar rcs $@ $^

# Compile each library object explicitly
$(BUILD_DIR)/CHIC.o: src/CHIC.cpp $(HEADER) | dirs
	@echo "Compiling: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/CHIC_EARTH.o: src/CHIC_EARTH.cpp $(HEADER) | dirs
	@echo "Compiling: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Build all examples
examples: lib $(EXAMPLE_BINS)

# Compile each example .cpp -> build/.o
$(BUILD_DIR)/%.o: $(EXAMPLE_DIR)/%.cpp $(HEADER) | dirs
	@echo "Compiling: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link each example binary
$(EXAMPLE_DIR)/%: $(BUILD_DIR)/%.o $(LIB_TARGET)
	@echo "Linking: $@"
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.a $(EXAMPLE_BINS)
	@echo "Cleaned compilation artifacts."
