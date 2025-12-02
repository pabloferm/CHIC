# CHIC - Makefile (fixed)
CXX = g++
PYTHON_EXEC = python3

PYBIND_INCLUDES := $(shell $(PYTHON_EXEC) -m pybind11 --includes 2>/dev/null)
PYTHON_LDFLAGS := $(shell $(PYTHON_EXEC)-config --ldflags 2>/dev/null)
PYTHON_EXT_SUFFIX := $(shell $(PYTHON_EXEC)-c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))" 2>/dev/null)

# Compiler settings (LTO disabled while debugging)
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -funroll-loops -DNDEBUG -fPIC
INCLUDES = -I/usr/include/eigen3 -I. -Isrc $(PYBIND_INCLUDES)

# Sources
LIB_SOURCE = src/CHIC.cpp
BINDING_SOURCE = bindings/pybind_CHIC.cpp
HEADER = src/CHIC.h

# Build artifacts
BUILD_DIR = build
LIB_OBJ = $(BUILD_DIR)/CHIC.o
BINDING_OBJ = $(BUILD_DIR)/pybind_CHIC.o
LIB_TARGET = $(BUILD_DIR)/libchic.a
PYTHON_TARGET = bindings/pychic$(PYTHON_EXT_SUFFIX)

# Ensure dirs safely via a recipe (do not use $(shell ...) that may clobber)
.PHONY: all lib python clean dirs
all: dirs $(LIB_TARGET)

dirs:
	@mkdir -p $(BUILD_DIR)

lib: $(LIB_TARGET)
python: dirs $(PYTHON_TARGET)

# Build static library
$(LIB_TARGET): $(LIB_OBJ)
	@echo "Creating static library: $@"
	ar rcs $@ $^

$(LIB_OBJ): $(LIB_SOURCE) $(HEADER) | dirs
	@echo "Compiling library object: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Build binding object
$(BINDING_OBJ): $(BINDING_SOURCE) $(HEADER) | dirs
	@echo "Compiling Python bindings object: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link Python extension — link the static archive directly (no -lchic)
$(PYTHON_TARGET): $(LIB_TARGET) $(BINDING_OBJ)
	@echo "Linking Python module: $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -shared -o $@ $(LIB_TARGET) $(BINDING_OBJ) $(PYTHON_LDFLAGS)

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.a bindings/*.so
	@echo "Cleaned build artifacts"

