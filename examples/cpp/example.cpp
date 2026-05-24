#include "CHIC.h"
#include <iostream>

int main() {
    double energy = 0.2;  // GeV
    double baseline = 295.;  // km

    // CHIC 
    CHIC chic;
    
    chic.update_dcp(4.1);  // dcp in rad (234 deg)
    chic.update_dm231(2.5e-3); // eV^2/c^4
    chic.update_dm221(7.5e-5); // eV^2/c^4
    chic.update_th12(0.583638); // theta_12 in rad (33.44 deg)
    chic.update_th13(0.149574); // theta_13 in rad (8.57 deg)
    chic.update_th23(0.858701); // theta_23 in rad (49.2 deg)
    chic.update_density(3.0); // g/cm^3
    
    for (int i=0; i<100; i++) {
    Eigen::Matrix3d prob;
    double en = energy + i*0.8/100.;
    prob = chic.compute_oscillations(
        en, baseline);
    std::cout << en << " " << prob(1,0) << " " << prob(1,1) << " " << prob(1,2) << " " << prob(1,0) + prob(1,1) + prob(1,2)<< std::endl;
    }

    // CHICDIFF 
    CHICDIFF dchic;
    
    Eigen::Matrix3d prob;
    Eigen::Matrix3d dprob;
    prob = dchic.compute_oscillations(
        energy, baseline);
    dprob = dchic.compute_oscillations_derivatives(
        "dcp", energy, baseline);
    dprob = dchic.compute_oscillations_derivatives(
        "density", energy, baseline);
    dprob = dchic.compute_oscillations_derivatives(
        "th23", energy, baseline);

    return 0;
}

