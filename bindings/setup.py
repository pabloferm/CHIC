# bindings/setup.py
import os
import re
import subprocess
import sys
from setuptools import setup, Extension

from pybind11.setup_helpers import Pybind11Extension, build_ext

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

os.environ.setdefault("CXX", "g++")
os.environ.setdefault("CC",  "gcc")


# ---------------------------------------------------------------------------
# Eigen auto-detection (mirrors the Makefile logic so both use the same path)
# ---------------------------------------------------------------------------
def _read_eigen_version(eigen_root: str) -> str:
    """Parse EIGEN_*_VERSION defines from Eigen's Macros.h."""
    macros_h = os.path.join(eigen_root, "Eigen", "src", "Core", "util", "Macros.h")
    if not os.path.isfile(macros_h):
        return "unknown"
    values = {}
    pattern = re.compile(r"#define\s+(EIGEN_WORLD_VERSION|EIGEN_MAJOR_VERSION|EIGEN_MINOR_VERSION)\s+(\d+)")
    with open(macros_h) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                values[m.group(1)] = m.group(2)
    return "{}.{}.{}".format(
        values.get("EIGEN_WORLD_VERSION", "?"),
        values.get("EIGEN_MAJOR_VERSION", "?"),
        values.get("EIGEN_MINOR_VERSION", "?"),
    )


def find_eigen() -> str:
    """
    Return the Eigen include directory (e.g. '/usr/include/eigen3').
    
    Raises SystemExit with an actionable message if Eigen cannot be found.
    """
    # 1. Explicit override
    override = os.environ.get("EIGEN_INCLUDE", "").strip()
    if override:
        path = override.lstrip("-I")          # accept both raw path and -I<path>
        if os.path.isfile(os.path.join(path, "Eigen", "Dense")):
            print(f"Eigen found via EIGEN_INCLUDE override: {path}")
            return path
        sys.exit(f"[setup.py] EIGEN_INCLUDE is set to '{path}' but "
                 "Eigen/Dense was not found there.")

    # 2. pkg-config
    try:
        result = subprocess.run(
            ["pkg-config", "--cflags-only-I", "eigen3"],
            capture_output=True, text=True, check=True,
        )
        # pkg-config returns e.g. "-I/usr/include/eigen3 "
        include_flag = result.stdout.strip().split()[0]   # take first token
        path = include_flag.lstrip("-I")
        version_result = subprocess.run(
            ["pkg-config", "--modversion", "eigen3"],
            capture_output=True, text=True,
        )
        version = version_result.stdout.strip() or _read_eigen_version(path)
        print(f"Eigen {version} found via pkg-config: {path}")
        return path
    except (subprocess.CalledProcessError, FileNotFoundError, IndexError):
        pass

    # 3. Common install paths
    search_paths = [
        "/usr/include/eigen3",
        "/usr/local/include/eigen3",
        "/opt/homebrew/include/eigen3",   # macOS Homebrew (Apple Silicon / Intel)
        "/opt/local/include/eigen3",      # MacPorts
        "/usr/include/Eigen",
        "/usr/local/include/Eigen",
    ]
    for path in search_paths:
        if os.path.isfile(os.path.join(path, "Eigen", "Dense")):
            version = _read_eigen_version(path)
            print(f"Eigen {version} found at: {path}")
            return path



EIGEN_INCLUDE = find_eigen()

# ---------------------------------------------------------------------------

ext_modules = [
    Pybind11Extension(
        "pychic",
        sources=[
            os.path.join(ROOT, "bindings", "pybind_CHIC.cpp"),
            os.path.join(ROOT, "src", "CHIC.cpp"),
        ],
        include_dirs=[
            os.path.join(ROOT, "src"),
            EIGEN_INCLUDE,
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
    version="1.1.0",
    description="Python bindings for CHIC",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy", "pybind11>=2.6"],
    zip_safe=False,
)