#include <pybind11/pybind11.h>
#include "chic.hpp"

namespace py = pybind11;

PYBIND11_MODULE(chic_cpp, m) {
    py::class_<chic::Example>(m, "Example")
        .def(py::init<>())
        .def("add", &chic::Example::add);
}