![Alt text](docs/CHIC.png)

# CHIC

C++ project with Python bindings (pybind11) and tests for both C++ and Python components.

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

- For the CHIC library
```
make clean
make
```

- For the python bindings (optional)
```
make python
```
