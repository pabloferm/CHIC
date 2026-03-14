#ifndef CHICEARTH_H
#define CHICEARTH_H

#include <Eigen/Dense>
#include <cstddef>
#include "opt_constants.h"
#include "EARTH.h"
#include "CHIC.h"


// =================================================== \\
// ====== Class for Earth neutrino propagation ======= \\
// =================================================== \\

class CHICEARTH {
 public:
    // Constructor - inherits from base class
    explicit CHICEARTH(std::string_view mode = "neutrino",
                      double theta_12 = 0.5836381018669037,
                      double theta_23 = 0.8587019919812102,
                      double theta_13 = 0.14957471689591406,
                      double delta_cp = 4.084070449666731,
                      double dm2_21 = 7.42e-5,
                      double dm2_31 = 2.51e-3,
                      std::string_view model = "PREM4");

    // Called to compute oscillations
    Eigen::Matrix3d compute_oscillations(double E, double cos_zenith);

    ~CHICEARTH() {
        delete chic;
        delete[] tracks;
    }
    
    // Update dcp
    inline void update_dcp(double delta_cp) {
        chic->update_dcp(delta_cp);
    }
    
    // Update theta_23
    inline void update_th23(double theta_23) {
        chic->update_th23(theta_23);
    }

    // Update theta_12
    inline void update_th12(double theta_12) {
        chic->update_th12(theta_12);
    }

    // Update theta_13
    inline void update_th13(double theta_13) {
        chic->update_th13(theta_13);
    }

    // Update dm221
    inline void update_dm221(double dm_2_21) {
        chic->update_dm221(dm_2_21);
    }

    // Update dm231
    inline void update_dm231(double dm_2_31) {
        chic->update_dm231(dm_2_31);
    }

 private:
    void _compute_hamiltonians(double E);
    void _amplitude();
    void _build_track(double cos_zenith);
    // we need to compute eigenvalues and Hs for a given energy
    // that is several calls to _compute_hamiltonians
    // once we have that for each density layer, cache linked to energy
    // use the same logic as in CHIC for the baseline

    std::vector<Eigen::Vector3d> lambdas;
    std::vector<Eigen::Vector3d> prod_lambdas;
    Eigen::Matrix3cd J;
    double* tracks;
    bool need_update = true;
    const EarthModel& Earth;
    CHIC* chic = nullptr;
    int deepest;
    double cos_zenith0, E0;

};

#endif // CHICEARTH_H
