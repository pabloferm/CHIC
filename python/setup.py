from setuptools import setup, Extension
import pybind11

setup(
    name="chic",
    packages=["chic"],
    ext_modules=[
        Extension(
            "chic_cpp",
            ["chic/bindings.cpp"],
            include_dirs=[
                pybind11.get_include(),
                "../../cpp/include"
            ],
            language="c++"
        )
    ],
    zip_safe=False,
)