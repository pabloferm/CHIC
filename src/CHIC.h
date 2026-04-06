#ifndef CHIC_H
#define CHIC_H

#include <Eigen/Dense>
#include <cstddef>
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
        double dm2_21 = 7.42e-5,  // eV^2/c^4
        double dm2_31 = 2.51e-3,  // eV^2/c^4
        double density = 2.8,     // g/cm^3
        double Y_e = 0.5   // effective electron fraction
        );

    // Update dcp
    inline void update_dcp(double delta_cp) {
        e_idelta = std::exp(std::complex<double>(0, flip * delta_cp));
        e_midelta = 1.0 / e_idelta;
        update_pmns = true;
    }
    // Update theta_23
    inline void update_th23(double theta_23) {c23 = std::cos(theta_23); s23 = std::sin(theta_23); update_pmns = true; }
    // Update theta_12
    inline void update_th12(double theta_12) {c12 = std::cos(theta_12); s12 = std::sin(theta_12); update_pmns = true; }
    // Update theta_13
    inline void update_th13(double theta_13) {c13 = std::cos(theta_13); s13 = std::sin(theta_13); update_pmns = true; }

    // Update dm221
    inline void update_dm221(double dm_2_21) {dm2_21 = dm_2_21; update_pmns = true;}
    // Update dm231
    inline void update_dm231(double dm_2_31) {dm2_31 = dm_2_31; update_pmns = true;}
    // Update density
    inline void update_density(double rho) { 
        density = rho; 
        // if (density == 0.0) {
        //     vacuum = true;
        //     update_matter = false;
        //     return
        // }
        update_matter = true;
    }
    // Update Ye
    inline void update_Ye(double Ye) { 
        Y_e = Ye;      
        // if (Y_e == 0.0) {
        //     vacuum = true;
        //     update_matter = false;
        //     return
        // }
        update_matter = true; }

    // Called to compute oscillations
    Eigen::Matrix3d compute_oscillations(double E, double L);
    // Compute hamiltonian and eigenvalues
    void compute_hamiltonians(double E);

    // Access to main results
    const Eigen::Matrix3cd& get_hamiltonian        () const noexcept { return Hs;          }
    const Eigen::Matrix3cd& get_hamiltonian_squared() const noexcept { return Hs2;         }
    const Eigen::Matrix3cd& get_amplitude          () const noexcept { return J;           }
    const Eigen::Vector3d&  get_eigenvalues        () const noexcept { return lambdas;     }
    const Eigen::Vector3d&  get_prod_eigenvalues   () const noexcept { return prod_lambdas;}
    const Eigen::Vector3d&  get_diff_eigenvalues   () const noexcept { return diff_lambdas;}

    // Access to vacuum variables
    const Eigen::Matrix3cd& get_Hs0            () const noexcept { return Hs0;             }
    const Eigen::Matrix3cd& get_Hs0_squared    () const noexcept { return Hs0_2;           }
    const Eigen::Matrix3cd& get_re_Hs0Vs0      () const noexcept { return re_Hs0Vs0;       }
    double                  get_trK0           () const noexcept { return trK0;            }
    double                  get_detH0          () const noexcept { return detH0;           }
    double                  get_flip           () const noexcept { return flip;            }

 protected:
    static constexpr double EPSILON = 1e-20;

    // Scalars for the Hamiltonian
    double s12{}, c12{}, s13{}, c13{}, s23{}, c23{};
    double dm2_21{}, dm2_31{}, V{}, Y_e{}, density{};
    double E0{}, L0{}, inv_2E{}, inv_4E_squared{};
    double flip = 1.0;
    double Vs2_diag0{}, Vs2_diag1{};
    std::complex<double> iL{}, e_idelta{}, e_midelta{};
    
    // Invariants
    double TrH, TrHs2, DetHs, detH0, trK0;

    // Eigenvalues (lambdas)
    Eigen::Vector3d prod_lambdas, diff_lambdas, lambdas;
    Eigen::Vector3cd exp_lambdas, cJ;
    
    const Eigen::Matrix3d identity_cache = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d zero_cache = Eigen::Matrix3d::Zero();
    
    // Pre-computed matrices
    Eigen::Matrix3cd U;                        // PMNS matrix
    Eigen::Matrix3cd Hs0, Hs;                  // Reduced Hamiltonian
    Eigen::Matrix3cd Hs2, Hs0_2;               // Reduced squared Hamiltonian
    Eigen::Matrix3cd re_Hs0Vs, re_Hs0Vs0;      // Vacuum hamiltonian and matter potential product
    Eigen::Matrix3d m2 = zero_cache;           // Mass squared difference matrix
    Eigen::Matrix3d v = zero_cache;            // Matter potential matrix
    Eigen::Matrix3d Vs = zero_cache;           // Reduced matter potential matrix
    Eigen::Matrix3d Vs2 = zero_cache;          // Reduced matter potential squared
    Eigen::Matrix3cd J;                        // Amplitude matrix
    
    // Flag to denote if matrix update is needed.
    bool update_pmns = true;
    bool update_matter = true;
    bool vacuum = false;
    
    void pmns_matrix();
    void _set_vacuum();
    void _set_matter();
    void _exponential();
    void _compute_hamiltonians();
    void _amplitude(double E, double L);

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
    const Eigen::Matrix3cd& get_diff_hamiltonian() const noexcept { return dHs; } // Derivative of the reduced Hamiltonian
    const Eigen::Matrix3cd& get_diff_amplitude  () const noexcept { return dJ;  } // Derivative of the amplitude    

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
    void dHs_dE();

};

#endif // CHIC_H