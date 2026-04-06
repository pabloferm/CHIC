#include "CHIC.h"
#include <complex>

// =================================================================
// == Base class for CHIC neutrino oscillation probabilities only ==
// =================================================================

CHIC::CHIC(std::string_view mode,
           double theta_12,
           double theta_23,
           double theta_13,
           double delta_cp,
           double dm2_21,
           double dm2_31,
           double density,
           double Y_e) {

    // Account for neutrino or antineutrino oscillations
    if (mode == "neutrino" || mode == "nu") {
        flip = 1.0;
    } else if (mode == "antineutrino" || mode == "nub" || mode == "nubar" ||
            mode == "anti-neutrino" || mode == "antinu") {
        flip = -1.0;
    } else {
        throw std::invalid_argument(
            "Not a valid (neutrino or anti-neutrino) mode selected.");
    }

    m2 = zero_cache;

    // Apply values
    update_th12(theta_12);
    update_th23(theta_23);
    update_th13(theta_13);
    update_dcp(delta_cp);
    update_dm221(dm2_21);
    update_dm231(dm2_31);
    update_Ye(Y_e);
    update_density(density);

    // Construct basic matrices
    _set_vacuum();
    _set_matter();
}

// PMNS matrix
void CHIC::pmns_matrix() {
    U(0, 0) = c12 * c13;
    U(0, 1) = s12 * c13;
    U(0, 2) = s13 * e_midelta;

    U(1, 0) = -s12 * c23 - c12 * s23 * s13 * e_idelta;
    U(1, 1) =  c12 * c23 - s12 * s23 * s13 * e_idelta;
    U(1, 2) =  s23 * c13;

    U(2, 0) =  s12 * s23 - c12 * c23 * s13 * e_idelta;
    U(2, 1) = -c12 * s23 - s12 * c23 * s13 * e_idelta;
    U(2, 2) =  c23 * c13;
}

// Vacuum Hamiltonian and related
void CHIC::_set_vacuum() {
    pmns_matrix();

    m2(1, 1) = dm2_21;
    m2(2, 2) = dm2_31;

    // Traceless vacuum Hamiltonian w/o energy-independent part
    trK0 = dm2_21 + dm2_31;
    Hs0  = U * m2 * U.adjoint() - (trK0 * OptConstants::INV_3) * identity_cache;

    // Terms of the shifted Hamiltonian squared
    Hs0_2 = Hs0 * Hs0;
    re_Hs0Vs0 = OptConstants::INV_3 * Hs0;
    re_Hs0Vs0(0,0) *= 4.0; 
    re_Hs0Vs0(1,1) *= -2.0; 
    re_Hs0Vs0(1,2) *= -2.0; 
    re_Hs0Vs0(2,1) *= -2.0; 
    re_Hs0Vs0(2,2) *= -2.0; 

    // Determinant term calculation independent of energy and density
    const double det_term = std::norm(U(1, 1) * U(2, 2) - U(1, 2) * U(2, 1));
    detH0 = dm2_21 * dm2_31 * det_term;

    update_pmns   = false;
    update_matter = true;   // V hasn't been applied yet — force matter rebuild
}

// Matter potential and related
void CHIC::_set_matter() {
    V = flip * OptConstants::WEAK * Y_e * density;
    const double V3 = V * OptConstants::INV_3;

    // Full cross-term for single-density use: re_Hs0Vs = re_Hs0Vs_unit * V
    re_Hs0Vs = re_Hs0Vs0 * V;

    // Vs² is diagonal: diag(4V²/9, V²/9, V²/9) — store the two distinct entries
    Vs2_diag0 = 4.0 * V3 * V3;   // (2V/3)²
    Vs2_diag1 = V3 * V3;         // (V/3)²
    
    update_matter = false;
}

void CHIC::_compute_hamiltonians() {
    // Full traceless Hamiltonian
    Hs = Hs0 * inv_2E;
    Hs(0, 0) += 2.0 * V * OptConstants::INV_3;
    Hs(1, 1) -= V * OptConstants::INV_3;
    Hs(2, 2) -= V * OptConstants::INV_3;
    // Full traceless Hamiltonian squared:
    Hs2.noalias() = Hs0_2 * inv_4E_squared + re_Hs0Vs * inv_2E;
    Hs2(0, 0) += Vs2_diag0;
    Hs2(1, 1) += Vs2_diag1;
    Hs2(2, 2) += Vs2_diag1;

    // Trace of full Hamiltonian (used for DetHs)
    TrH       = inv_2E * trK0 + V;
    TrHs2 = std::real(Hs2.trace()); // !?
    DetHs = detH0 * V * inv_4E_squared
          + TrHs2 * TrH  * OptConstants::INV_6
          - TrH * TrH * TrH * OptConstants::INV_27;

    // Calculation of eigenvalues: trigonometric-quadratic-linear method
    if (TrHs2 < EPSILON) {
        throw std::invalid_argument(
            "Eigenvalues are (almost) degenerate, please check your input parameters.");
    } else {
        lambdas[0] = std::sqrt(OptConstants::TWO_OVER_3 * TrHs2) *
            std::cos(std::acos(OptConstants::SQRT_54 * DetHs / TrHs2 / std::sqrt(TrHs2)) * OptConstants::INV_3);
        lambdas[1] = -(0.5 * lambdas[0] -
            std::sqrt(0.25 * lambdas[0] * lambdas[0] - DetHs / lambdas[0]));
        lambdas[2] = -(lambdas[0] + lambdas[1]);
    }

    // Cyclic product and difference of eigenvalues
    prod_lambdas = DetHs / lambdas.array();
    diff_lambdas = 1.0 / (2.0 * lambdas.array().square() + prod_lambdas.array());
}

void CHIC::compute_hamiltonians(double E) {
    if (update_pmns)   _set_vacuum();
    if (update_matter) _set_matter();
    if (E != E0) { E0 = E; inv_2E = 0.5 / E0; inv_4E_squared = inv_2E * inv_2E; _compute_hamiltonians();}
}

void CHIC::_exponential() {
    // Exponentials
    iL = std::complex<double>(0, -L0 * OptConstants::BASELINE_FACTOR);
    exp_lambdas = diff_lambdas.array() * (iL*lambdas).array().exp();

    // Amplitude coefficients
    cJ[0] = prod_lambdas.dot(exp_lambdas);
    cJ[1] = lambdas.dot(exp_lambdas);
    cJ[2] = exp_lambdas.sum();

    // Amplitude
    J.noalias() = cJ[2] * Hs2 + cJ[1] * Hs;
    J.diagonal().array() += cJ[0];
}

// Compute oscillations for a given energy E and baseline L
Eigen::Matrix3d CHIC::compute_oscillations(double E, double L) {
    _amplitude(E, L);
    return J.cwiseAbs2();
}

void CHIC::_amplitude(double E, double L) {
    if (update_pmns)          { _set_vacuum();           E0 = E; _compute_hamiltonians(); L0 = L; }
    else if (update_matter)   { _set_matter();           E0 = E; _compute_hamiltonians(); L0 = L; }
    else if (E != E0)        { E0 = E; _compute_hamiltonians(); L0 = L; }
    else if (L == L0)        { return; }
    else                     { L0 = L; }

    _exponential();
}


// ====================================================================== \\
// = Derived class for derivative calculations with lazy initialization = \\
// ====================================================================== \\

CHICDIFF::CHICDIFF(std::string_view mode,
                   double theta_12,
                   double theta_23,
                   double theta_13,
                   double delta_cp,
                   double dm2_21,
                   double dm2_31,
                   double density,
                   double Y_e)
    : CHIC(mode, theta_12, theta_23, theta_13, delta_cp,
           dm2_21, dm2_31, density, Y_e) {}


Eigen::Matrix3d CHICDIFF::compute_oscillations_derivatives(std::string_view param, double E, double L) {
    _amplitude(E, L);

    param0 = param;
    _set_dHs();

    _amplitude_derivative();

    return 2.0 * dJ.cwiseProduct(J.conjugate()).real();
}

void CHICDIFF::_set_dHs() {
    dHs = zero_cache;
    if (param0 == "density") {
        dHs_drho();
    } else if (param0 == "dm231") {
        dHs_ddm231();
    } else if (param0 == "dm221") {
        dHs_ddm221();
    } else if (param0 == "dcp") {
        dHs_ddcp();
    } else if (param0 == "th23") {
        dHs_dth23();
    } else if (param0 == "th13") {
        dHs_dth13();
    } else if (param0 == "th12") {
        dHs_dth12();
    } else if (param0 == "E") {
        dHs_dE();
    } else {
        throw std::invalid_argument(std::string("Unknown parameter: ") + std::string(param0));
    }
    if (param0 == "E") {
      dHs = -2 * inv_4E_squared * dHs;
    } else if (param0 != "density") {
      dHs = inv_2E * dHs;
    }
    comm_dHH.noalias()  = dHs * Hs  + Hs  * dHs;
    comm_dHH2.noalias() = dHs * Hs2 + Hs2 * dHs;
}

// Computes amplitude derivative matrix
void CHICDIFF::_amplitude_derivative() {
    // Integral matrix (symmetric)
    I_ijk.diagonal() = iL * diff_lambdas.cwiseProduct(exp_lambdas);
    I_ijk(0, 1) = I_ijk(1, 0) = (diff_lambdas[0]*exp_lambdas[1] - diff_lambdas[1]*exp_lambdas[0])
                / (lambdas[1] - lambdas[0]);
    I_ijk(0, 2) = I_ijk(2, 0) = (diff_lambdas[2]*exp_lambdas[0] - diff_lambdas[0]*exp_lambdas[2])
                / (lambdas[0] - lambdas[2]);
    I_ijk(1, 2) = I_ijk(2, 1) = (diff_lambdas[1]*exp_lambdas[2] - diff_lambdas[2]*exp_lambdas[1])
                / (lambdas[2] - lambdas[1]);

    Eigen::Vector3cd I_prod = I_ijk * prod_lambdas;

    cdJ_diag[0] = prod_lambdas.dot(I_prod);    // dH

    I_prod = I_ijk * lambdas;
    cdJ_diag[1] = lambdas.dot(I_prod);       // H dH H
    cdJ_off[0] = prod_lambdas.dot(I_prod);   // dH H

    I_prod = I_ijk * unit;
    cdJ_diag[2] = unit.dot(I_prod);            // H2 dH H2
    cdJ_off[1] = prod_lambdas.dot(I_prod);     // dH H2
    cdJ_off[2] = lambdas.dot(I_prod);          // H dH H2

    dJ.noalias() = cdJ_diag[0] * dHs;
    dJ.noalias()+= cdJ_diag[1] * Hs * dHs * Hs;
    dJ.noalias()+= cdJ_diag[2] * Hs2 * dHs * Hs2;
    dJ.noalias()+= cdJ_off[0] * comm_dHH;
    dJ.noalias()+= cdJ_off[1] * comm_dHH2;
    dJ.noalias()+= cdJ_off[2] * Hs * comm_dHH * Hs;
}

// === Derivative Hamiltonian builders ===

void CHICDIFF::dHs_drho() {
    dHs(0, 0) = flip * OptConstants::WEAK * Y_e;
    dHs.diagonal().array() -= flip * OptConstants::WEAK * Y_e * OptConstants::INV_3;
}

void CHICDIFF::dHs_ddm221() {
    Eigen::Matrix3d dm = zero_cache;
    dm(1, 1) = 1.0;
    dHs = U * dm * U.adjoint();
    dHs.diagonal().array() -= OptConstants::INV_3;
}

void CHICDIFF::dHs_ddm231() {
    Eigen::Matrix3d dm = zero_cache;
    dm(2, 2) = 1.0;
    dHs = U * dm * U.adjoint();
    dHs.diagonal().array() -= OptConstants::INV_3;
}

void CHICDIFF::dHs_dth23() {
    Eigen::Matrix3cd dU = zero_cache;
    dU(1, 0) = U(2, 0);
    dU(1, 1) =  U(2, 1);
    dU(1, 2) = U(2, 2);
    dU(2, 0) = - U(1, 0);
    dU(2, 1) = - U(1, 1);
    dU(2, 2) = - U(1, 2);
    dHs.noalias() = dU * m2 * U.adjoint() + U * m2 * dU.adjoint();
}

void CHICDIFF::dHs_dth12() {
    Eigen::Matrix3cd dU = zero_cache;
    dU(0, 0) = - U(0, 1);
    dU(0, 1) =  U(0, 0);
    dU(1, 0) = - U(1, 1);
    dU(1, 1) = U(1, 0);
    dU(2, 0) = - U(2, 1);
    dU(2, 1) = U(2, 0);
    dHs.noalias() = dU * m2 * U.adjoint() + U * m2 * dU.adjoint();
}

void CHICDIFF::dHs_dth13() {
    Eigen::Matrix3cd dU = zero_cache;
    dU(0, 0) = - c12 * s13;
    dU(0, 1) = - s12 * s13;
    dU(0, 2) = c13 * e_midelta;
    dU(1, 0) = - c12 * s23 * c13 * e_idelta;
    dU(1, 1) = - s12 * s23 * c13 * e_idelta;
    dU(1, 2) = - s23 * s13;
    dU(2, 0) = - c12 * c23 * c13 * e_idelta;
    dU(2, 1) = - s12 * c23 * c13 * e_idelta;
    dU(2, 2) = - c23 * s13;
    dHs.noalias() = dU * m2 * U.adjoint() + U * m2 * dU.adjoint();
}

void CHICDIFF::dHs_ddcp() {
    Eigen::Matrix3cd dU = zero_cache;
    std::complex<double> mi_flip(0.0, -flip);
    dU(0, 2) = mi_flip * s13 * e_midelta;
    dU(1, 0) = mi_flip * c12 * s23 * s13 * e_idelta;
    dU(1, 1) = mi_flip * s12 * s23 * s13 * e_idelta;
    dU(2, 0) = mi_flip * c12 * c23 * s13 * e_idelta;
    dU(2, 1) = mi_flip * s12 * c23 * s13 * e_idelta;
    dHs.noalias() = dU * m2 * U.adjoint() + U * m2 * dU.adjoint();
}

void CHICDIFF::dHs_dE() {
    dHs = Hs0;
}