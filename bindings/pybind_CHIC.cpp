#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "CHIC.h"
#include "CHIC_EARTH.h"

namespace py = pybind11;
using Matrix3d = Eigen::Matrix<double, 3, 3>;

// ======================================================================== //
//  pychic — constant-density neutrino oscillations                         //
// ======================================================================== //

PYBIND11_MODULE(pychic, m) {
    m.doc() = "Python bindings for the CHIC (Cayley-Hamilton Invariants and "
              "Constants) C++ library. Provides fast neutrino oscillation "
              "probabilities and their analytic derivatives for a medium of "
              "constant density, using the Cayley-Hamilton matrix formalism.";

    // ------------------------------------------------------------------ //
    //  CHIC                                                               //
    // ------------------------------------------------------------------ //
    py::class_<CHIC>(m, "CHIC",
        R"doc(
        Neutrino oscillations in matter of constant density.

        Computes the 3×3 oscillation-probability matrix for a neutrino
        travelling through a uniform-density medium using the
        Cayley-Hamilton method.

        Parameters
        ----------
        mode : str
            "neutrino" or "antineutrino" (default "neutrino").
        theta_12, theta_23, theta_13 : float
            Mixing angles in degrees.
        delta_cp : float
            CP-violating phase in degrees.
        dm2_21, dm2_31 : float
            Mass-squared differences in eV².
        density : float
            Matter density in g/cm³ (default 2.8).
        Y_e : float
            Electron fraction (default 0.5).
        )doc")

        .def(py::init<std::string, double, double, double, double,
                      double, double, double, double>(),
             py::arg("mode")     = "neutrino",
             py::arg("theta_12")       = 0.5836381018669037,
             py::arg("theta_23")       = 0.8587019919812102,
             py::arg("theta_13")       = 0.14957471689591406,
             py::arg("delta_cp")       = 4.084070449666731,
             py::arg("dm2_21")   = 7.42e-5,
             py::arg("dm2_31")   = 2.51e-3,
             py::arg("density")  = 2.8,
             py::arg("Y_e")      = 0.5)

        .def("update_dcp", &CHIC::update_dcp,
             py::arg("delta_cp"),
             "Update the CP-violating phase (degrees).")
        .def("update_th23", &CHIC::update_th23,
             py::arg("theta_23"),
             "Update the theta_23 mixing angle (degrees).")
        .def("update_th12", &CHIC::update_th12,
             py::arg("theta_12"),
             "Update the theta_12 mixing angle (degrees).")
        .def("update_th13", &CHIC::update_th13,
             py::arg("theta_13"),
             "Update the theta_13 mixing angle (degrees).")
        .def("update_dm221", &CHIC::update_dm221,
             py::arg("dm2_21"),
             "Update the solar mass-squared difference dm²₂₁ (eV²).")
        .def("update_dm231", &CHIC::update_dm231,
             py::arg("dm2_31"),
             "Update the atmospheric mass-squared difference dm²₃₁ (eV²).")
        .def("update_density", &CHIC::update_density,
             py::arg("density"),
             "Update the matter density (g/cm³).")

        .def("compute_oscillations", &CHIC::compute_oscillations,
             R"doc(
             Compute the 3×3 oscillation probability matrix.

             Parameters
             ----------
             E : float
                 Neutrino energy in GeV.
             L : float
                 Baseline length in km.

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype float64
                 P[α][β] = probability that a flavour-α neutrino is
                 detected as flavour β.
             )doc")

        // ---- accessors ------------------------------------------------ //
        .def("get_hamiltonian", &CHIC::get_hamiltonian,
             R"doc(
             Return the effective Hamiltonian in the flavour basis.

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype complex128
             )doc")
        .def("get_amplitude", &CHIC::get_amplitude,
             R"doc(
             Return the transition amplitude matrix after the last call
             to compute_oscillations().

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype complex128
             )doc")
        .def("get_eigenvalues", &CHIC::get_eigenvalues,
             R"doc(
             Return the three eigenvalues of the effective Hamiltonian.

             Returns
             -------
             numpy.ndarray, shape (3,), dtype float64
             )doc");

    // ------------------------------------------------------------------ //
    //  CHICDIFF  (derives from CHIC)                                      //
    // ------------------------------------------------------------------ //
    py::class_<CHICDIFF, CHIC>(m, "CHICDIFF",
        R"doc(
        Neutrino oscillations in constant-density matter with analytic
        derivatives.

        Extends CHIC to additionally compute the derivative of the
        oscillation probability matrix with respect to any oscillation
        parameter, using the Cayley-Hamilton formalism.

        Constructor arguments are identical to CHIC.
        )doc")

        .def(py::init<std::string, double, double, double, double,
                      double, double, double, double>(),
             py::arg("mode")     = "neutrino",
             py::arg("theta_12")       = 0.5836381018669037,
             py::arg("theta_23")       = 0.8587019919812102,
             py::arg("theta_13")       = 0.14957471689591406,
             py::arg("delta_cp")       = 4.084070449666731,
             py::arg("dm2_21")   = 7.42e-5,
             py::arg("dm2_31")   = 2.51e-3,
             py::arg("density")  = 2.8,
             py::arg("Y_e")      = 0.5)

        .def("compute_oscillations_derivatives",
             &CHICDIFF::compute_oscillations_derivatives,
             R"doc(
             Compute the derivative of the oscillation probability matrix
             with respect to one oscillation parameter.

             Parameters
             ----------
             param : str
                 Parameter name: one of "theta_12", "theta_23", "theta_13",
                 "delta_cp", "dm2_21", "dm2_31".
             E : float
                 Neutrino energy in GeV.
             L : float
                 Baseline length in km.

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype float64
                 dP[α][β] / d(param).
             )doc")
        .def("get_diff_hamiltonian", &CHICDIFF::get_diff_hamiltonian,
             R"doc(
             Return the derivative of the effective Hamiltonian with respect
             to the last requested parameter.

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype complex128
             )doc")
        .def("get_diff_amplitude", &CHICDIFF::get_diff_amplitude,
             R"doc(
             Return the derivative of the transition amplitude matrix with
             respect to the last requested parameter.

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype complex128
             )doc");
}

// ======================================================================== //
//  pychic_earth — neutrino propagation through the Earth                   //
// ======================================================================== //

PYBIND11_MODULE(pychic_earth, m) {
    m.doc() = "Python bindings for CHIC EARTH — neutrino oscillation "
              "propagation through the Earth using the Cayley-Hamilton "
              "method with a layered density profile.";

    // ------------------------------------------------------------------ //
    //  CHICEARTH                                                          //
    // ------------------------------------------------------------------ //
    py::class_<CHICEARTH>(m, "CHICEARTH",
        R"doc(
        Neutrino propagation through the Earth.

        Computes the full oscillation-probability matrix for a neutrino
        travelling through the Earth at a given energy and zenith angle,
        using a layered Earth density model (e.g. PREM).

        Parameters
        ----------
        mode : str
            "neutrino" or "antineutrino" (default "neutrino").
        theta_12, theta_23, theta_13 : float
            Mixing angles in radians.
        delta_cp : float
            CP-violating phase in radians.
        dm2_21, dm2_31 : float
            Mass-squared differences in eV².
        model : str
            Earth density model identifier, e.g. "PREM10" (default).
        detector_depth : float
            Depth of the detector below the surface in km (default 0.0).
        )doc")

        .def(py::init<std::string_view,
                      double, double, double, double,
                      double, double,
                      std::string_view, double>(),
             py::arg("mode")           = "neutrino",
             py::arg("theta_12")       = 0.5836381018669037,
             py::arg("theta_23")       = 0.8587019919812102,
             py::arg("theta_13")       = 0.14957471689591406,
             py::arg("delta_cp")       = 4.084070449666731,
             py::arg("dm2_21")         = 7.42e-5,
             py::arg("dm2_31")         = 2.51e-3,
             py::arg("model")          = "PREM10",
             py::arg("detector_depth") = 0.0)

        .def("update_dcp", &CHICEARTH::update_dcp,
             py::arg("delta_cp"),
             "Update the CP-violating phase (radians).")
        .def("update_th23", &CHICEARTH::update_th23,
             py::arg("theta_23"),
             "Update the theta_23 mixing angle (radians).")
        .def("update_th12", &CHICEARTH::update_th12,
             py::arg("theta_12"),
             "Update the theta_12 mixing angle (radians).")
        .def("update_th13", &CHICEARTH::update_th13,
             py::arg("theta_13"),
             "Update the theta_13 mixing angle (radians).")
        .def("update_dm221", &CHICEARTH::update_dm221,
             py::arg("dm2_21"),
             "Update the solar mass-squared difference dm²₂₁ (eV²).")
        .def("update_dm231", &CHICEARTH::update_dm231,
             py::arg("dm2_31"),
             "Update the atmospheric mass-squared difference dm²₃₁ (eV²).")

        .def("compute_oscillations", &CHICEARTH::compute_oscillations,
             py::arg("E"), py::arg("cos_zenith"), py::arg("h") = 0.0,
             R"doc(
             Compute the 3×3 oscillation probability matrix.

             Parameters
             ----------
             E : float
                 Neutrino energy in GeV.
             cos_zenith : float
                 Cosine of the zenith angle (−1 = upward-going through the
                 full Earth diameter, 0 = horizontal).
             h : float, optional
                 Production height above the surface in km (default 0.0).

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype float64
                 P[α][β] = probability that a flavour-α neutrino arrives
                 as flavour β.
             )doc")

        .def("get_amplitude", &CHICEARTH::get_amplitude,
             R"doc(
             Return the accumulated transition amplitude matrix J.

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype complex128
                 The full propagation amplitude after the last call to
                 compute_oscillations().
             )doc");

    // ------------------------------------------------------------------ //
    //  CHICEARTHDIFF  (derives from CHICEARTH)                           //
    // ------------------------------------------------------------------ //
    py::class_<CHICEARTHDIFF, CHICEARTH>(m, "CHICEARTHDIFF",
        R"doc(
        Neutrino propagation through the Earth with analytic derivatives.

        Extends CHICEARTH to also compute the derivative of the oscillation
        probability matrix with respect to oscillation parameters, using the
        Cayley-Hamilton formalism.

        Constructor arguments are identical to CHICEARTH.
        )doc")

        .def(py::init<std::string_view,
                      double, double, double, double,
                      double, double,
                      std::string_view, double>(),
             py::arg("mode")           = "neutrino",
             py::arg("theta_12")       = 0.5836381018669037,
             py::arg("theta_23")       = 0.8587019919812102,
             py::arg("theta_13")       = 0.14957471689591406,
             py::arg("delta_cp")       = 4.084070449666731,
             py::arg("dm2_21")         = 7.42e-5,
             py::arg("dm2_31")         = 2.51e-3,
             py::arg("model")          = "PREM10",
             py::arg("detector_depth") = 0.0)

        .def("compute_oscillations",
             &CHICEARTHDIFF::compute_oscillations,
             py::arg("E"), py::arg("cos_zenith"), py::arg("h") = 0.0,
             R"doc(
             Compute the 3×3 oscillation probability matrix (and cache
             intermediate quantities needed for subsequent derivative calls).

             Parameters
             ----------
             E : float
                 Neutrino energy in GeV.
             cos_zenith : float
                 Cosine of the zenith angle.
             h : float, optional
                 Production height in km (default 0.0).

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype float64
             )doc")

        .def("compute_oscillations_derivatives",
             &CHICEARTHDIFF::compute_oscillations_derivatives,
             py::arg("param"), py::arg("E"), py::arg("cos_zenith"),
             py::arg("h") = 0.0,
             R"doc(
             Compute the derivative of the oscillation probability matrix
             with respect to one oscillation parameter.

             Parameters
             ----------
             param : str
                 Parameter name: one of "theta_12", "theta_23", "theta_13",
                 "delta_cp", "dm2_21", "dm2_31", "E".
             E : float
                 Neutrino energy in GeV.
             cos_zenith : float
                 Cosine of the zenith angle.
             h : float, optional
                 Production height in km (default 0.0).

             Returns
             -------
             numpy.ndarray, shape (3, 3), dtype float64
                 dP[α][β] / d(param).
             )doc");
}