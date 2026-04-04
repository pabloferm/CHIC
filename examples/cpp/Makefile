CXX      := g++
CXXFLAGS := -std=c++17 -Ofast -march=native -ffast-math -funroll-loops -DNDEBUG -fPIC\
            -I../../src -I. \
            $(shell pkg-config --cflags eigen3 2>/dev/null || echo "-I/usr/include/eigen3")
LDFLAGS  :=

# Add CHIC.cpp and any other translation units from your src directory here
DEPS     := ../../src/CHIC_EARTH.cpp ../../src/CHIC.cpp

all: test_chicearth benchmark_chicearth

# ---- Full test: spot checks + unitarity + gnuplot dump ----
test_chicearth: test_chicearth.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ---- Performance-only benchmark ----
test_perf: benchmark_chicearth.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run: test_chicearth
	./test_chicearth

perf: benchmark_chicearth
	./benchmark_chicearth

plot: run
	gnuplot plot_oscillations.gp

clean:
	rm -f test_chicearth benchmark_chicearth oscillations.dat oscillations.png plot_oscillations.gp
