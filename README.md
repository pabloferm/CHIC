<img src="docs/CHIC.png" width="200" />

# CHIC

C++ library for computing three/flavor neutrino oscillation probabilities and their derivatives.

## Structure

- `src/`: C++ source and headers
- `bindings/`: Python package and bindings
- `examples/`: Minimal examples for C++ and Python

## Requirements

- Eigen library: https://libeigen.gitlab.io/?title=Main_Page
(in the following a standard installation path (i.e. /usr/include/eigen3) is assumed)
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
or run a simple benchmark program.compile and run it as
```
examples/cpp/chic_benchmark.cpp
```
To compile them, just do
```
make examples
```
The executables will appear under the `examples/cpp/` folder.

- For the python module, you can check the following minimal example:
```
python3 examples/python/example.py
```

- Additionally, there is a Jupyter notebook:
```
examples/python/example.ipynb
```

## Reference
You can check and cite the associated paper [here](https://arxiv.org/abs/2512.16427).
