# CHIC - Makefile
CXX      = g++
CXXFLAGS = -std=c++26 -O3 -march=native -ffast-math -funroll-loops \
           -fPIC -DNDEBUG \
           -fomit-frame-pointer \
           -fno-trapping-math \
           -fassociative-math \
           -freciprocal-math \
           -ffinite-math-only

# --- Eigen auto-detection ---
EIGEN_INCLUDE :=
 
# 1. Try pkg-config 
ifneq ($(shell pkg-config --exists eigen3 2>/dev/null && echo yes),)
    EIGEN_INCLUDE := $(shell pkg-config --cflags-only-I eigen3)
    EIGEN_VERSION := $(shell pkg-config --modversion eigen3)
    $(info Eigen $(EIGEN_VERSION) found via pkg-config: $(EIGEN_INCLUDE))
else
    # 2. Common installation paths
    EIGEN_SEARCH_PATHS := \
        /usr/include/eigen3 \
        /usr/local/include/eigen3 \
        /opt/homebrew/include/eigen3 \
        /opt/local/include/eigen3 \
        /usr/include/Eigen \
        /usr/local/include/Eigen
 
    EIGEN_INCLUDE := $(firstword $(foreach p,$(EIGEN_SEARCH_PATHS),\
        $(if $(wildcard $(p)/Eigen/Dense),-I$(p),)))
 
    ifeq ($(EIGEN_INCLUDE),)
        $(error Eigen not found. Install it (e.g. apt install libeigen3-dev) \
or set EIGEN_INCLUDE manually: make EIGEN_INCLUDE=-I/your/path/to/eigen3)
    endif
 
    # Display Eigen version used
    EIGEN_WORLD  := $(shell grep -r 'define EIGEN_WORLD_VERSION' \
        $(patsubst -I%,%,$(EIGEN_INCLUDE))/Eigen/src/Core/util/Macros.h 2>/dev/null \
        | awk '{print $$3}')
    EIGEN_MAJOR  := $(shell grep -r 'define EIGEN_MAJOR_VERSION' \
        $(patsubst -I%,%,$(EIGEN_INCLUDE))/Eigen/src/Core/util/Macros.h 2>/dev/null \
        | awk '{print $$3}')
    EIGEN_MINOR  := $(shell grep -r 'define EIGEN_MINOR_VERSION' \
        $(patsubst -I%,%,$(EIGEN_INCLUDE))/Eigen/src/Core/util/Macros.h 2>/dev/null \
        | awk '{print $$3}')
    EIGEN_VERSION := $(EIGEN_WORLD).$(EIGEN_MAJOR).$(EIGEN_MINOR)
    $(info Eigen $(EIGEN_VERSION) found at: $(EIGEN_INCLUDE))
endif

INCLUDES = $(EIGEN_INCLUDE) -I. -Isrc

# Build artifacts
BUILD_DIR  = build
LIB_OBJ    = $(BUILD_DIR)/CHIC.o
LIB_SOURCE = src/CHIC.cpp src/CHIC_EARTH.cpp src/CHIC_BATCH.cpp
HEADER     = src/CHIC.h src/CHIC_EARTH.h src/CHIC_BATCH.h

# C++ build
BUILD_DIR = build
LIB_OBJ   = $(BUILD_DIR)/CHIC.o $(BUILD_DIR)/CHIC_EARTH.o $(BUILD_DIR)/CHIC_BATCH.o
LIB_TARGET = $(BUILD_DIR)/libchic.a

# Python artifacts
BUILD_PY = bindings/build 
EGG_PY = bindings/chic.egg-info

# Examples
EXAMPLE_DIR  = examples/cpp
EXAMPLE_BINS = $(EXAMPLE_DIR)/chic_benchmark $(EXAMPLE_DIR)/example $(EXAMPLE_DIR)/chicdiff_benchmark $(EXAMPLE_DIR)/chic_profile $(EXAMPLE_DIR)/chicdiff2_benchmark $(EXAMPLE_DIR)/chic_batch_example

# Python
PYTHON = bindings

# Prevent make from deleting intermediate files (library, objects)
.SECONDARY:

EXAMPLE_BINS = $(EXAMPLE_DIR)/example $(EXAMPLE_DIR)/test_chicearth $(EXAMPLE_DIR)/speed $(EXAMPLE_DIR)/chic_batch_example

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

$(BUILD_DIR)/CHIC_BATCH.o: src/CHIC_BATCH.cpp $(HEADER) | dirs
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
	rm -rf $(BUILD_DIR)
	@echo "Cleaned C++ build artifacts"
	rm -rf $(BUILD_PY) $(EGG_PY)
	@echo "Cleaned Python build artifacts"
	rm -f $(EXAMPLE_BINS)
	@echo "Cleaned example binaries"
