#include "CHIC.h"

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

    // Traceless vacuum Hamiltonian w/o energy-independent part
    trK0 = dm2_21 + dm2_31;
    Eigen::Matrix3cd Um2 = U;
    Um2.col(0).setZero();
    Um2.col(1) *= dm2_21;
    Um2.col(2) *= dm2_31;
    Hs0.noalias() = Um2 * U.adjoint();
    Hs0.diagonal().array() -= trK0 * OptConstants::INV_3;

    // Terms of the shifted Hamiltonian squared
    Hs0_2.noalias() = Hs0 * Hs0;
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

    // Vs² is diagonal
    Vs2_diag0 = OptConstants::FOUR_OVER_NINE * V * V;   // (2V/3)²
    Vs2_diag1 = OptConstants::INV_9 * V * V;         // (V/3)²
    
    update_matter = false;
}

void CHIC::_compute_hamiltonians() {
    // Full traceless Hamiltonian
    inv_2E = 0.5 / E0; 
    inv_4E_squared = inv_2E * inv_2E;
    Hs = Hs0 * inv_2E;
    Hs(0, 0) += 2.0 * V * OptConstants::INV_3;
    Hs(1, 1) -= V * OptConstants::INV_3;
    Hs(2, 2) -= V * OptConstants::INV_3;

    // Full traceless Hamiltonian squared:
    Hs2.noalias() = Hs0_2 * inv_4E_squared + re_Hs0Vs0 * V * inv_2E;
    Hs2(0, 0) += Vs2_diag0;
    Hs2(1, 1) += Vs2_diag1;
    Hs2(2, 2) += Vs2_diag1;

    // Trace of full Hamiltonian (used for DetHs)
    TrH       = inv_2E * trK0 + V;
    TrHs2 = Hs2(0,0).real() + Hs2(1,1).real() + Hs2(2,2).real();
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
    bool recompute = (E != E0);
    if (update_pmns)   { _set_vacuum();  recompute = true; }
    if (update_matter) { _set_matter();  recompute = true; }
    if (recompute)     { E0 = E; _compute_hamiltonians();  }
}

void CHIC::_exponential() {
    // Exponentials
    iL = -L0 * OptConstants::BASELINE_FACTOR;

    double s, c;
    for (int k = 0; k < 3; ++k) {
        __builtin_sincos(iL * lambdas[k], &s, &c);
        exp_lambdas[k] = diff_lambdas[k] * std::complex<double>(c, s);
    }

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
    bool recompute = (E != E0);
    if (update_pmns)   { _set_vacuum();  recompute = true; }
    if (update_matter) { _set_matter();  recompute = true; }
    if (recompute)     { E0 = E; _compute_hamiltonians(); L0 = L; }
    else if (L == L0)  { return; }
    else               { L0 = L; }

    _exponential();
}


// ====================================================================== 
// = Derived class for derivative calculations with lazy initialization = 
// ====================================================================== 

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
    
    if (param0 == "L") {
        dHs.setZero();
        dJ = std::complex<double>(0.0, - OptConstants::BASELINE_FACTOR) * (Hs*J);
    }
    else {
        _set_dHs();
        _amplitude_derivative();
        if (param0 == "EL") {
            Eigen::Matrix3cd dJL = std::complex<double>(0.0, - OptConstants::BASELINE_FACTOR) * (Hs*J);
            Eigen::Matrix3d ddprob = dJ.cwiseProduct(dJL.conjugate()).real();
            dJ = std::complex<double>(0.0, - OptConstants::BASELINE_FACTOR) * (-2.0*Hs0/inv_4E_squared * J + Hs*dJ);
            ddprob += dJ.cwiseProduct(J.conjugate()).real();
            return 2.0 * ddprob;
        }
    }
    
    return 2.0 * dJ.cwiseProduct(J.conjugate()).real();
}

void CHICDIFF::compute_hamiltonians_and_derivatives(std::string_view param, double E) {
    compute_hamiltonians(E);
    param0 = param;
    _set_dHs();
}

void CHICDIFF::_set_dHs() {
    dHs.setZero();
    if      (param0 == "density") { dHs_drho();   }
    else if (param0 == "dm231")   { dHs_ddm231(); }
    else if (param0 == "dm221")   { dHs_ddm221(); }
    else if (param0 == "dcp")     { dHs_ddcp();   }
    else if (param0 == "th23")    { dHs_dth23();  }
    else if (param0 == "th13")    { dHs_dth13();  }
    else if (param0 == "th12")    { dHs_dth12();  }
    else if ((param0 == "E") || (param0 == "EL")) { dHs_dE();     }
    else if (param0 == "L")       { dHs.setZero();}
    else {
        throw std::invalid_argument(std::string("Unknown parameter: ") + std::string(param0));
    }

    if ((param0 == "E") || (param0 == "EL")) {
      dHs = -2 * inv_4E_squared * dHs;
    } else if (param0 != "density") {
      dHs = inv_2E * dHs;
    }
    
    Eigen::Matrix3cd dHs_Hs, dHs_Hs2;
    dHs_Hs.noalias()   = dHs * Hs;
    dHs_Hs2.noalias()  = dHs * Hs2;
    
    comm_dHsHs.noalias()  = dHs_Hs + dHs_Hs.adjoint();    // = dHs·Hs + Hs·dHs

    comm_dHsHs2.noalias() = dHs_Hs2 + dHs_Hs2.adjoint();  // = dHs·Hs2 + Hs2·dHs

    Hs_dHs_Hs.noalias()   = Hs  * dHs_Hs;   // cached: Hs·dHs·Hs
    Hs2_dHs_Hs2.noalias() = Hs2 * dHs_Hs2;  // cached: Hs2·dHs·Hs2

    // Hs_comm_dHsHs_Hs.noalias() = Hs * comm_dHsHs * Hs;
    Hs_comm_dHsHs_Hs.noalias()  = Hs_dHs_Hs * Hs;
    Hs_comm_dHsHs_Hs.noalias() += Hs2 * dHs_Hs;
}

// Computes amplitude derivative matrix
void CHICDIFF::_amplitude_derivative() {
    // Integral matrix (symmetric)
    I_ijk.diagonal() = std::complex<double>(0.0, iL) * diff_lambdas.cwiseProduct(exp_lambdas);
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

    I_prod = I_ijk * Eigen::Vector3cd::Ones();
    cdJ_diag[2] = Eigen::Vector3cd::Ones().dot(I_prod);            // H2 dH H2
    cdJ_off[1] = prod_lambdas.dot(I_prod);     // dH H2
    cdJ_off[2] = lambdas.dot(I_prod);          // H dH H2

    dJ.noalias() = cdJ_diag[0] * dHs;
    dJ.noalias()+= cdJ_diag[1] * Hs_dHs_Hs;
    dJ.noalias()+= cdJ_diag[2] * Hs2_dHs_Hs2;
    dJ.noalias()+= cdJ_off[0] * comm_dHsHs;
    dJ.noalias()+= cdJ_off[1] * comm_dHsHs2;
    dJ.noalias()+= cdJ_off[2] * Hs_comm_dHsHs_Hs;
}

// === Derivative Hamiltonian builders ===

void CHICDIFF::dHs_drho() {
    dHs(0, 0) = flip * OptConstants::WEAK * Y_e;
    dHs.diagonal().array() -= flip * OptConstants::WEAK * Y_e * OptConstants::INV_3;
}

void CHICDIFF::dHs_ddm221() {
    dHs.noalias() = U.col(1) * U.col(1).adjoint();
    dHs.diagonal().array() -= OptConstants::INV_3;
}

void CHICDIFF::dHs_ddm231() {
    dHs.noalias() = U.col(2) * U.col(2).adjoint();
    dHs.diagonal().array() -= OptConstants::INV_3;
}

void CHICDIFF::dHs_dth23() {
    dHs.noalias()  = Eigen::Vector3cd(0.0,  U(2,1)*dm2_21, -U(1,1)*dm2_21) * U.col(1).adjoint();
    dHs.noalias() += Eigen::Vector3cd(0.0,  U(2,2)*dm2_31, -U(1,2)*dm2_31) * U.col(2).adjoint();
    const Eigen::Matrix3cd tmp = dHs; dHs.noalias() += tmp.adjoint();
}

void CHICDIFF::dHs_dth12() {
    dHs.noalias()  = (U.col(0) * dm2_21) * U.col(1).adjoint();  // rank-1: 9 cmul
    const Eigen::Matrix3cd tmp = dHs; dHs.noalias() += tmp.adjoint();
}

void CHICDIFF::dHs_dth13() {
    dHs.noalias()  = Eigen::Vector3cd(-s12*s13               * dm2_21,
                                        -s12*s23*c13*e_idelta  * dm2_21,
                                        -s12*c23*c13*e_idelta  * dm2_21) * U.col(1).adjoint();
    dHs.noalias() += Eigen::Vector3cd( c13*e_midelta            * dm2_31,
                                       -s23*s13                 * dm2_31,
                                       -c23*s13                 * dm2_31) * U.col(2).adjoint();
    const Eigen::Matrix3cd tmp = dHs; dHs.noalias() += tmp.adjoint();
}

void CHICDIFF::dHs_ddcp() {

    const std::complex<double> mi_flip(0.0, -flip);
    dHs.noalias()   = Eigen::Vector3cd(0.0,
                                         mi_flip*s12*s23*s13*e_idelta * dm2_21,
                                         mi_flip*s12*c23*s13*e_idelta * dm2_21) * U.col(1).adjoint();
    dHs.row(0)     += (mi_flip * s13 * e_midelta * dm2_31) * U.col(2).adjoint();
    const Eigen::Matrix3cd tmp = dHs; dHs.noalias() += tmp.adjoint();
}

void CHICDIFF::dHs_dE() {
    dHs = Hs0;
}
