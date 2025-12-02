# bindings/setup.py
import os
from setuptools import setup, Extension
import pybind11
import sysconfig

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Force same compiler as your Makefile
os.environ.setdefault("CXX", "g++")
os.environ.setdefault("CC", "gcc")

# Flags exactly matching your Makefile
compile_args = [
    "-std=c++17",
    "-O3",
    "-march=native",
    "-ffast-math",
    "-funroll-loops",
    "-DNDEBUG",
    "-fPIC",
]

# Make sure extension suffix is the same used by your Python
ext_name = "pychic"

sources = [
    os.path.join(ROOT, "bindings", "pybind_CHIC.cpp"),
    os.path.join(ROOT, "src", "CHIC.cpp"),
]

include_dirs = [
    os.path.join(ROOT, "src"), 
    pybind11.get_include(),
    "/usr/include/eigen3",
]

# Force linking against libstdc++ and keep lib paths discovered by python
extra_link_args = [
    "-lstdc++",
]

ext_modules = [
    Extension(
        ext_name,
        sources=sources,
        include_dirs=include_dirs,
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=extra_link_args,
    )
]

setup(
    name="chic",
    version="1.0.0",
    description="Python bindings for CHIC",
    ext_modules=ext_modules,
    install_requires=["numpy", "pybind11>=2.6"],
    zip_safe=False,
)
