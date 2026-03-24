# bindings/setup.py
import os
import sys
from setuptools import setup, Extension

# At this point pybind11 IS available — it's in [build-system].requires
# and setuptools guarantees build deps are installed before setup.py runs.
# The previous failures were because the *old* setup.py was cached / still
# present without a pyproject.toml alongside it. With both files present,
# setuptools installs build-system.requires first, then runs setup.py.
#
# If you still hit issues, the nuclear option is:
#   pip-3.14 install pybind11 first, then pip-3.14 install bindings/

from pybind11.setup_helpers import Pybind11Extension, build_ext

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

os.environ.setdefault("CXX", "g++")
os.environ.setdefault("CC",  "gcc")

ext_modules = [
    Pybind11Extension(
        "pychic",
        sources=[
            os.path.join(ROOT, "bindings", "pybind_CHIC.cpp"),
            os.path.join(ROOT, "src", "CHIC.cpp"),
        ],
        include_dirs=[
            os.path.join(ROOT, "src"),
            "/usr/include/eigen3",
        ],
        extra_compile_args=[
            "-std=c++17",
            "-O3",
            "-march=native",
            "-ffast-math",
            "-funroll-loops",
            "-DNDEBUG",
            "-fPIC",
        ],
        extra_link_args=["-lstdc++"],
        cxx_std=17,
    )
]

setup(
    name="chic",
    version="1.0.1",
    description="Python bindings for CHIC",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy", "pybind11>=2.6"],
    zip_safe=False,
)
