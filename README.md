<img src="docs/CHIC.png" width="200" />

# CHIC

C++ library for computing three/flavor neutrino oscillation probabilities and their derivatives.

## Structure

- `src/`: C++ source and headers
- `bindings/`: Python package and bindings
- `tests/`: Tests for C++ and Python

## Requirements

- Eigen library: https://libeigen.gitlab.io/?title=Main_Page
- Only if building python bindings:
	- pybind11: https://pybind11.readthedocs.io/en/stable/installing.html
	- numpy: https://numpy.org/install/

## Build instructions

- For the C++ CHIC library only
```
make clean
make 
```
or 
```
make lib
```

- For the python bindings (optional)
```
pip install ./bindings
```
and making it available in system-wide
```
export PYTHONPATH=$PWD/bindings/:$PYTHONPATH
```

## Examples

- For C++ library, you can check the following minimal example:
```
examples/cpp/example.cpp
```
compile and run it as
```
g++ -o example examples/cpp/example.cpp -Lbuild -lchic -I src -I/usr/include/eigen3
./example
```

- For the python module, you can check the following minimal example:
```
python3 examples/python/example.py
```
