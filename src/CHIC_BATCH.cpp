#include "CHIC_BATCH.h"
#include <stdexcept>

CHIC_BATCH::CHIC_BATCH(bool derivatives) 
    : diff(derivatives) {
    start_chic();
}

// Store initial neutrino flux before oscillations

void CHIC_BATCH::set_initial_flux(const std::vector<std::vector<double>>& flux) {
    
    initial_flux = flux;
}

// Calling CHIC(DIFF)
void CHIC_BATCH::start_chic() {
    // Check if CHIC object for this mode already exists
    chic_objects[0] = std::make_unique<CHICDIFF>("neutrino");
    chic_objects[1] = std::make_unique<CHICDIFF>("antineutrino");
}

void CHIC_BATCH::update_parameter(const std::string& param, double value) {
    
    // Update all CHIC objects
    for (int nu=0; nu<2; ++nu) {
        
        if (param == "theta_12" || param == "th12") {
            chic_objects[nu]->update_th12(value);
        } 
        else if (param == "theta_23" || param == "th23") {
            chic_objects[nu]->update_th23(value);
        } 
        else if (param == "theta_13" || param == "th13") {
            chic_objects[nu]->update_th13(value);
        } 
        else if (param == "delta_cp" || param == "dcp") {
            chic_objects[nu]->update_dcp(value);
        } 
        else if (param == "dm2_21" || param == "dm221") {
            chic_objects[nu]->update_dm221(value);
        } 
        else if (param == "dm2_31" || param == "dm231") {
            chic_objects[nu]->update_dm231(value);
        } 
        else if (param == "density" || param == "rho") {
            chic_objects[nu]->update_density(value);
        } 
        else if (param == "y_e" || param == "ye") {
            chic_objects[nu]->update_Ye(value);
        } 
        else {
            throw std::invalid_argument(
                "Unknown parameter: '" + param + "'. Supported parameters: "
                "theta_12, theta_23, theta_13, delta_cp, dm2_21, dm2_31, density, Y_e"
            );
        }
    }
}

std::vector<double> CHIC_BATCH::compute_oscillations(
    const std::vector<int>& flavors,
    const std::vector<double>& energies,
    const std::vector<double>& baselines) {
    
    if (initial_flux.empty() || initial_flux[0].size() != 3) {
    throw std::runtime_error("Initial flux not set. Call set_initial_flux().");
}
    
    weights.resize(flavors.size());
    Eigen::Matrix3d prob_matrix;

    // Iterate over all events
    for (size_t i = 0; i < flavors.size(); ++i) {
        if (flavors[i] < 0) {
            prob_matrix = chic_objects[0]->compute_oscillations(energies[i], baselines[i]);
        }
        else {
            prob_matrix = chic_objects[1]->compute_oscillations(energies[i], baselines[i]);
        }
        weights[i] = 0.0;
        int target_flavor = int(0.5 * std::abs(flavors[i]) - 6.0);
        for (int initial_flavor = 0; initial_flavor < 3; ++initial_flavor) {
            weights[i] += prob_matrix(initial_flavor, target_flavor) * initial_flux[i][initial_flavor];
        }
    }
    return weights;
}

std::vector<double> CHIC_BATCH::compute_oscillations_derivatives(
    const std::string& param,
    const std::vector<int>& flavors,
    const std::vector<double>& energies,
    const std::vector<double>& baselines) {
    
    if (!diff) {
        throw std::runtime_error(
            "Derivatives not enabled. Create CHIC_BATCH with derivatives=true"
        );
    }
    
    if (energies.empty()) throw std::invalid_argument("energies array cannot be empty");
    if (baselines.empty()) throw std::invalid_argument("baselines array cannot be empty");
    
    dweights.resize(flavors.size());
    Eigen::Matrix3d prob_matrix;
    
    // Iterate over all events
    for (size_t i = 0; i < flavors.size(); ++i) {
        if (flavors[i] < 0) {
            prob_matrix = chic_objects[0]->compute_oscillations_derivatives(param, energies[i], baselines[i]);
        }
        else {
            prob_matrix = chic_objects[1]->compute_oscillations_derivatives(param, energies[i], baselines[i]);
        }
        dweights[i] = 0.0;
        int target_flavor = int(0.5 * std::abs(flavors[i]) - 6.0);
        for (int initial_flavor = 0; initial_flavor < 3; ++initial_flavor) {
            dweights[i] += prob_matrix(initial_flavor, target_flavor) * initial_flux[i][initial_flavor];
        }
    }
    return dweights;
}
