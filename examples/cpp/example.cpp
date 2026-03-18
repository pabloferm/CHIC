#include "CHIC.h"
#include <iostream>

int main() {
    double energy = 0.125603;  // GeV
    double baseline = 250.;  // km

    // CHIC
    CHIC chic;

    chic.update_dcp(0.066);  // dcp in rad (234 deg)
    chic.update_dm231(-2.509e-3); // eV^2/c^4
    chic.update_dm221(7.53e-5); // eV^2/c^4
    chic.update_th12(0.5872523687443223); // theta_12 in rad (33.44 deg)
    chic.update_th13(0.14819001778459273); // theta_13 in rad (8.57 deg)
    chic.update_th23(0.8134128187551903); // theta_23 in rad (49.2 deg)
    chic.update_density(2.6); // g/cm^3

    Eigen::Matrix3d prob;
    prob = chic.compute_oscillations(
        energy, baseline);
    std::cout <<prob<<std::endl;


    // CHICDIFF
    CHICDIFF dchic;

    Eigen::Matrix3d dprob;
    prob = dchic.compute_oscillations(
        energy, baseline);
    dprob = dchic.compute_oscillations_derivatives(
        "dcp", energy, baseline);
    dprob = dchic.compute_oscillations_derivatives(
        "density", energy, baseline);
    dprob = dchic.compute_oscillations_derivatives(
        "th23", energy, baseline);
    std::cout <<prob<<std::endl;

    return 0;
} 
