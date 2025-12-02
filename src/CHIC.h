#ifndef CHIC_H
#define CHIC_H

#include <Eigen/Dense>
#include "opt_constants.h"


// ================================================================= \\
// == Base class for CHIC neutrino oscillation probabilities only == \\
// ================================================================= \\

class CHIC {
 public:
    // Constructor with default parameters (angles in radians)
    explicit CHIC(std::string_view mode = "neutrino",
        double theta_12 = 0.5836381018669037,  // theta_12_deg = 33.44
        double theta_23 = 0.8587019919812102,  // theta_23_deg = 49.2
        double theta_13 = 0.14957471689591406, // theta_13_deg = 8.57
        double delta_cp = 4.084070449666731,   // delta_cp_deg = 234
        double dm2_21 = 7.42e-5,  // MeV^2/c^4
        double dm2_31 = 2.51e-3,  // MeV^2/c^4
        double density = 2.8,     // g/cm^3
        double Y_e = 0.5   // effective electron fraction
        );

    // Update dcp
    inline void update_dcp(double delta_cp) {
        e_idelta = std::exp(std::complex<double>(0, flip * delta_cp));
        e_midelta = 1.0 / e_idelta;
        need_update = true;
    }
    
    // Update theta_23
    inline void update_th23(double theta_23) {
        c23 = std::cos(theta_23);
        s23 = std::sin(theta_23);
        need_update = true;
    }

    // Update theta_12
    inline void update_th12(double theta_12) {
        c12 = std::cos(theta_12);
        s12 = std::sin(theta_12);
        need_update = true;
    }

    // Update theta_13
    inline void update_th13(double theta_13) {
        c13 = std::cos(theta_13);
        s13 = std::sin(theta_13);
        need_update = true;
    }

    // Update dm221
    inline void update_dm221(double dm_2_21) {
        dm2_21 = dm_2_21;
        need_update = true;
    }

    // Update dm231
    inline void update_dm231(double dm_2_31) {
        dm2_31 = dm_2_31;
        need_update = true;
    }

    // Update density
    inline void update_density(double rho) {
        density = rho;
        need_update = true;
    }

    // Called to compute oscillations
    Eigen::Matrix3d compute_oscillations(double E, double L);

    // Access to results
    Eigen::Matrix3cd get_hamiltonian();            // Reduced Hamiltonian
    Eigen::Matrix3cd get_amplitude();              // Coefficients of the amplitude
    Eigen::Vector3d get_eigenvalues();             // Eigenvalues

 protected:
    static constexpr double EPSILON = 1e-20;

    // Scalars for the Hamiltonian
    double s12{}, c12{}, s13{}, c13{}, s23{}, c23{};
    double dm2_21{}, dm2_31{}, V{}, Y_e{}, density{};
    double E0{}, L0{}, inv_2E{};
    double flip = 1.0;
    std::complex<double> iL{}, e_idelta{}, e_midelta{};
    
    // Pre-computed values for optimization
    double trH_cubed, detHs_squared, trHs2_cubed;  // Cache expensive determinant calculation
    
    // Invariants
    double TrH, TrH2, DetH, TrHs2, DetHs, detH0, trK0;
    double inv_4E_squared{};
    
    // Eigenvalues (lambdas)
    Eigen::Vector3d prod_lambdas, diff_lambdas, lambdas;
    Eigen::Vector3cd exp_lambdas, cJ;
    
    const Eigen::Matrix3d identity_cache = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d zero_cache = Eigen::Matrix3d::Zero();

    // Pre-computed matrices
    Eigen::Matrix3cd U;                        // PMNS matrix
    Eigen::Matrix3cd Hs0, Hs;                  // Reduced Hamiltonian
    Eigen::Matrix3cd Hs2, Hs0_2, re_Hs0Vs;     // Reduced squared Hamiltonian
    Eigen::Matrix3d Vs_2;                      // Matter potential matrix squared
    Eigen::Matrix3d m2 = zero_cache;           // Mass squared difference matrix
    Eigen::Matrix3d v = zero_cache;            // Matter potential matrix
    Eigen::Matrix3d Vs = zero_cache;           // Reduced matter potential matrix
    Eigen::Matrix3cd J;                        // Amplitude matrix
    
    // Flag to denote if matrix update is needed.
    bool need_update = true;
    
    void _compute_hamiltonians();
    void _amplitude();
    void pmns_matrix();
    void _set_matrices();
};




// =================================================== \\
// = Derived class for probabilities and derivatives = \\
// =================================================== \\

class CHICDIFF : public CHIC {
 public:
    // Constructor - inherits from base class
    explicit CHICDIFF(std::string_view mode = "neutrino",
                      double theta_12 = 0.5836381018669037,
                      double theta_23 = 0.8587019919812102,
                      double theta_13 = 0.14957471689591406,
                      double delta_cp = 4.084070449666731,
                      double dm2_21 = 7.42e-5,
                      double dm2_31 = 2.51e-3,
                      double density = 2.8,
                      double Y_e = 0.5);
    
    // Derivative calculation method
    Eigen::Matrix3d compute_oscillations_derivatives(std::string_view param, double E, double L);
    // Access to derivative results
    Eigen::Matrix3cd get_diff_hamiltonian();      // Derivative of the reduced Hamiltonian
    Eigen::Matrix3cd get_diff_amplitude();        // Derivative of the amplitude    

 private:
    std::string param0{"none"};
    const Eigen::Vector3d unit = Eigen::Vector3d::Ones(3);

    Eigen::Matrix3cd I_ijk;                 // Integral matrix
    Eigen::Matrix3cd comm_dHH, comm_dHH2;   // Commutator matrices
    Eigen::Matrix3cd dJ;              // Derivative amplitude matrix
    Eigen::Matrix3cd dHs;             // Derivative of reduced Hamiltonian
    std::complex<double> cdJ_diag[3], cdJ_off[3];  // Coefficients
    
    void _amplitude_derivative();
    void _set_dHs();
    
    // Build derivatives of the reduced Hamiltonian
    void dHs_drho();
    void dHs_ddm221();
    void dHs_ddm231();
    void dHs_dth23();
    void dHs_dth13();
    void dHs_dth12();
    void dHs_ddcp();

};

#endif // CHIC_H
