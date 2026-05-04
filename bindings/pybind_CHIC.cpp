#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "CHIC.h"  // Include your CHIC headers

namespace py = pybind11;
using Matrix3d = Eigen::Matrix<double, 3, 3>; // Explicitly define Matrix3d

PYBIND11_MODULE(pychic, m) {
    m.doc() = "Python bindings for the CHIC (Caley-Hamilton Invariants and Constants for neutrino oscillations) C++ library.";
    // Bind the CHIC class
    py::class_<CHIC>(m, "CHIC")
        .def(py::init<std::string, double, double, double, double, double, double, double, double>(),
             py::arg("mode") = "neutrino",
             py::arg("theta_12") = 0.5836381018669037,  // = 33.44°
             py::arg("theta_23") = 0.8587019919812102,  // = 49.2°,
             py::arg("theta_13") = 0.14957471689591406, // = 8.57°
             py::arg("delta_cp") = 4.084070449666731,   // = 234°
             py::arg("dm2_21") = 7.42e-5,
             py::arg("dm2_31") = 2.51e-3,
             py::arg("density") = 2.8,
             py::arg("Y_e") = 0.5)
        .def("update_dcp", &CHIC::update_dcp)
        .def("update_th23", &CHIC::update_th23)
        .def("update_th12", &CHIC::update_th12)
        .def("update_th13", &CHIC::update_th13)
        .def("update_dm221", &CHIC::update_dm221)
        .def("update_dm231", &CHIC::update_dm231)
        .def("update_density", &CHIC::update_density)
        .def("compute_oscillations", &CHIC::compute_oscillations)
        .def("get_hamiltonian", &CHIC::get_hamiltonian)
        .def("get_amplitude", &CHIC::get_amplitude)
        .def("get_eigenvalues", &CHIC::get_eigenvalues);

    // Bind the CHICDIFF class
    py::class_<CHICDIFF, CHIC>(m, "CHICDIFF")
        .def(py::init<std::string, double, double, double, double, double, double, double, double>(),
             py::arg("mode") = "neutrino",
             py::arg("theta_12") = 0.5836381018669037,  // = 33.44°
             py::arg("theta_23") = 0.8587019919812102,  // = 49.2°,
             py::arg("theta_13") = 0.14957471689591406, // = 8.57°
             py::arg("delta_cp") = 4.084070449666731,   // = 234°
             py::arg("dm2_21") = 7.42e-5,
             py::arg("dm2_31") = 2.51e-3,
             py::arg("density") = 2.8,
             py::arg("Y_e") = 0.5)
        .def("compute_oscillations_derivatives", &CHICDIFF::compute_oscillations_derivatives)
        .def("get_diff_hamiltonian", &CHICDIFF::get_diff_hamiltonian)
        .def("get_diff_amplitude", &CHICDIFF::get_diff_amplitude);
}

