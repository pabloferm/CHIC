#include "CHICEARTH.h"

// =================================================== \\
// ====== Class for Earth neutrino propagation ======= \\
// =================================================== \\

CHICEARTH::CHICEARTH(std::string_view mode,
           double theta_12,
           double theta_23,
           double theta_13,
           double delta_cp,
           double dm2_21,
           double dm2_31,
           std::string_view model) {

            // Initialize CHIC
            chic = new CHIC(mode, theta_12, theta_23, theta_13, delta_cp, dm2_21, dm2_31);

            // Get density layers
            if (model == "PREM4") {
                Earth = &PREM4;
            } else {
                throw std::invalid_argument("Not a valid Earth model.");
            }
            lambdas.resize(Earth->Nlayers);
            prod_lambdas.resize(Earth->Nlayers);
            tracks = new double[Earth->Nlayers];
            tracks[Earth->Nlayers - 1] = 1;
           }

// Computes the eigenvalues for each layer at a given energy
void CHICEARTH::_compute_hamiltonians(double E) {
    chic->compute_hamiltonians(E0);
    for (int i=0; i<Earth->Nlayers; i++){
        chic->update_density(Earth->density[i]);
        chic->update_Ye(Earth->Ye[i]);
        lambdas[i] = chic->get_eigenvalues();
        prod_lambdas[i] = chic->get_prod_eigenvalues();
        }
}

void CHICEARTH::_build_track(double cos_zenith) {
    if (cos_zenith0 == cos_zenith) {
        return;
    }
    double sin_zenith = std::sin(std::acos(cos_zenith0));
    double seg2;
    for (int i = 0; i < Earth->Nlayers; i++) {
        seg2 = Earth->radii[i] * Earth->radii[i] - sin_zenith * sin_zenith;
        if (seg2 > 0) {
            tracks[i] = std::sqrt(seg2);
            deepest = i;
        } else {
            tracks[i] = 0;
        }
    }
}

void CHICEARTH::_amplitude() {
    chic->_amplitude(2 * R_EARTH * tracks[deepest], lambdas[deepest], prod_lambdas[deepest]);
    J = chic->get_amplitude();
    for (int i = deepest; i>0; i--) {
        chic->_amplitude(R_EARTH * std::abs(tracks[i] - tracks[i-1]), lambdas[i], prod_lambdas[i]);
        J = chic->get_amplitude() * J * chic->get_amplitude(); // Symmetric of layers' path
    }
}

Eigen::Matrix3d CHICEARTH::compute_oscillations(double E, double cos_zenith) {
    _compute_hamiltonians(E);
        
    _build_track(cos_zenith);

    _amplitude();

    return J.cwiseAbs2();
}