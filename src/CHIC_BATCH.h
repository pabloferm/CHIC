#ifndef CHIC_BATCH_H
#define CHIC_BATCH_H

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "CHIC.h"

// CHIC_BATCH: Batch computation for CHIC

class CHIC_BATCH {
public:
    // Constructor: enables/disables derivative calculations
    explicit CHIC_BATCH(bool derivatives = true);

    // Store initial neutrino flux before oscillations
    void set_initial_flux(const std::vector<std::vector<double>>& flux);

    // Batch calculation
    std::vector<double> compute_oscillations(
        const std::vector<int>& flavors,
        const std::vector<double>& energies,
        const std::vector<double>& baselines
    );

    std::vector<double> compute_oscillations_derivatives(
        const std::string& param,
        const std::vector<int>& flavors,
        const std::vector<double>& energies,
        const std::vector<double>& baselines
    );

    void update_parameter(const std::string& param, double value);    void update_parameters(const std::map<std::string, double>& params);

    // Accessor methods
    const std::vector<std::vector<double>>& get_initial_flux() const noexcept { return initial_flux; }
    const std::vector<double>& get_weights() const noexcept { return weights; }
    bool derivatives_enabled() const noexcept { return diff; }

private:
    // Data members
    std::vector<std::vector<double>> initial_flux;                              // Initial flux per flavor
    std::vector<double> weights, dweights;                                   // Optional weighting factors
    bool diff;                                                      // Enable/disable derivatives
    std::unique_ptr<CHICDIFF> chic_objects[2]; // One CHICDIFF per mode

    // Helper methods
    void start_chic();
    
   
};

#endif // CHIC_BATCH_H