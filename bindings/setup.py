import setuptools
from setuptools import setup, Extension
import pybind11
import os

# --- Configuration ---
# 1. Source files to be compiled
sources = ['pybind_CHIC.cpp', '../src/CHIC.cpp']
sources = ["../bindings/pybind_CHIC.cpp", "../src/CHIC.cpp"]

# 2. Compiler flags
compile_args = ['-std=c++17', '-O3', '-march=native', '-ffast-math']

# 3. Include directories
# NOTE: This assumes Eigen is installed in a standard location like /usr/include/eigen3
# or that your system's compiler can find it automatically.
# Adjust EIGEN_INCLUDE_DIR if Eigen is located elsewhere (e.g., in a local vcpkg install).
try:
    # Try finding Eigen using a common system path
    EIGEN_INCLUDE_DIR = '/usr/include/eigen3'
    if not os.path.isdir(EIGEN_INCLUDE_DIR):
        EIGEN_INCLUDE_DIR = '/usr/local/include/eigen3'
except:
    EIGEN_INCLUDE_DIR = '/usr/include/eigen3' # Default fallback
include_dirs = [
    # pybind11 headers are mandatory
    pybind11.get_include(),
    # Eigen headers are mandatory for matrix conversion
    EIGEN_INCLUDE_DIR
]
include_dirs = ["../src", pybind11.get_include(), "/usr/include/eigen3"]
# ---------------------

chicos_module = Extension(
    'pychic',
    sources=sources,
    include_dirs=include_dirs,
    language='c++17',
    extra_compile_args=compile_args,
)

setup(
    name='chic',
    version='1.0.0',
    description='Python bindings for CHIC',
    ext_modules=[chicos_module],
    # Required dependencies for the Python module
    install_requires=[
        'numpy',
        'pybind11>=2.6'
    ],
    zip_safe=False,
)
