"""
setup.py — build script for the CHIC pybind11 extensions.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# ── Resolve project root (one level above bindings/) ─────────────────────────

HERE       = Path(__file__).parent.resolve()   # …/CHIC/bindings
ROOT       = HERE.parent                       # …/CHIC
SRC_DIR    = ROOT / "src"
BUILD_DIR  = ROOT / "build"

# ── Locate Eigen headers ─────────────────────────────────────────────────────

def _find_eigen() -> str:
    """
    Return the include path that exposes <Eigen/...> headers.

    Search order:
      1. pkg-config  (covers dnf/apt/brew system installs)
      2. Common hard-coded system paths
    """
    # 1. pkg-config
    try:
        result = subprocess.run(
            ["pkg-config", "--cflags-only-I", "eigen3"],
            capture_output=True, text=True, check=True,
        )
        # output is like "-I/usr/include/eigen3"
        path = result.stdout.strip().lstrip("-I").split()[0]
        if path and Path(path).is_dir():
            return path
    except (FileNotFoundError, subprocess.CalledProcessError, IndexError):
        pass

    # 2. Hard-coded fallback paths
    candidates = [
        "/usr/include/eigen3",                   # Fedora / Debian / Ubuntu
        "/usr/local/include/eigen3",
        "/opt/homebrew/include/eigen3",          # macOS Homebrew (Apple Silicon)
        "/opt/homebrew/opt/eigen/include/eigen3",
        "/usr/local/opt/eigen/include/eigen3",   # macOS Homebrew (Intel)
    ]
    for p in candidates:
        if Path(p).is_dir():
            return p

    raise RuntimeError(
        "Could not locate Eigen3 headers.\n"
        "  Fedora/RHEL:    sudo dnf install eigen3-devel\n"
        "  Debian/Ubuntu:  sudo apt install libeigen3-dev\n"
        "  macOS:          brew install eigen\n"
    )


# ── Sanity-check that libchic.a exists ───────────────────────────────────────

_LIB = BUILD_DIR / "libchic.a"
if not _LIB.exists():
    raise RuntimeError(
        f"Static library not found: {_LIB}\n"
        "Run make in the project root before installing the bindings."
    )

# ── Compiler / linker configuration ──────────────────────────────────────────

EIGEN_INCLUDE = _find_eigen()

# Mirror the flags used in the project Makefile so the ABI matches.
COMPILE_ARGS = [
    "-O3",
    "-std=c++26",
    "-march=native",
    "-ffast-math",
    "-funroll-loops",
    "-fPIC",
    "-DNDEBUG",
    "-fomit-frame-pointer",
    "-fno-trapping-math",
    "-fassociative-math",
    "-freciprocal-math",
    "-ffinite-math-only",
    "-fvisibility=hidden",
]
if sys.platform == "darwin":
    COMPILE_ARGS += ["-stdlib=libc++"]

# Link against the pre-built static library instead of recompiling the C++ sources.
LINK_ARGS = [str(_LIB)]

INCLUDE_DIRS = [
    EIGEN_INCLUDE,
    str(ROOT),       # for headers included as "CHIC.h" / "CHIC_EARTH.h"
    str(SRC_DIR),    # for headers included as "src/..."
]

# ── Extension definitions ─────────────────────────────────────────────────────

_common = dict(
    sources=[str(HERE / "pybind_CHIC.cpp")],
    include_dirs=INCLUDE_DIRS,
    extra_compile_args=COMPILE_ARGS,
    extra_link_args=LINK_ARGS,
    cxx_std=17,   # pybind11 helpers accept an int/str here; actual std set above
)

ext_modules = [
    Pybind11Extension(name="pychic",       **_common),
    Pybind11Extension(name="pychic_earth", **_common),
]

# ── Setup call ────────────────────────────────────────────────────────────────

setup(
    name="chic",
    version="2.0.0",
    description="Python bindings for CHIC",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy", "pybind11>=2.6"],
    zip_safe=False,
)