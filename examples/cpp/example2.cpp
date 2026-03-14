#include "CHICEARTH.h"

int main() {
    double energy = 0.7;  // GeV
    double cos_zenith = 0.5;  // km

    // CHIC 
    CHICEARTH chic;
    
    chic.update_dcp(4.1);  // dcp in rad (234 deg)
    chic.update_dm231(2.5e-3); // eV^2/c^4
    chic.update_dm221(7.5e-5); // eV^2/c^4
    chic.update_th12(0.583638); // theta_12 in rad (33.44 deg)
    chic.update_th13(0.149574); // theta_13 in rad (8.57 deg)
    chic.update_th23(0.858701); // theta_23 in rad (49.2 deg)
    
    Eigen::Matrix3d prob;
    prob = chic.compute_oscillations(
        energy, cos_zenith);


    // // CHICDIFF 
    // CHICDIFF dchic;
    
    // Eigen::Matrix3d dprob;
    // prob = dchic.compute_oscillations(
    //     energy, baseline);
    // dprob = dchic.compute_oscillations_derivatives(
    //     "dcp", energy, baseline);
    // dprob = dchic.compute_oscillations_derivatives(
    //     "density", energy, baseline);
    // dprob = dchic.compute_oscillations_derivatives(
    //     "th23", energy, baseline);

    // return 0;
}

