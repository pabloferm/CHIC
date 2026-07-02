#include "CHIC_BATCH_GRID_derived.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>

/**
 * Constructor: Initialize CHIC_BATCH_GRID (calls parent constructor)
 */
CHIC_BATCH_GRID::CHIC_BATCH_GRID() 
    : CHIC_BATCH(true) {
}

/**
 * Create and precompute oscillation grid with Hermite derivatives
 */
void CHIC_BATCH_GRID::create_and_compute_grid(
    int n_e, double E_min_val, double E_max_val,
    int n_l, double L_min_val, double L_max_val) {
    
    if (n_e < 2) throw std::invalid_argument("n_energy must be >= 2");
    if (n_l < 2) throw std::invalid_argument("n_baseline must be >= 2");
    if (E_min_val >= E_max_val) throw std::invalid_argument("E_min_val must be < E_max_val");
    if (L_min_val >= L_max_val) throw std::invalid_argument("L_min_val must be < L_max");
    if (E_min_val < 0.0) throw std::invalid_argument("Energy must be positive");
    if (L_min_val < 0.0) throw std::invalid_argument("Baseline must be positive");
        
    // Store parameters
    n_modes = mode_types.size();
    n_flavors = target_flavors.size();
    n_energy = n_e;
    n_baseline = n_l;
    E_min = E_min_val;
    E_max = E_max_val;
    L_min = L_min_val;
    L_max = L_max_val;
    
    // Create grid axes
    energy_grid.resize(n_energy);
    for (size_t i = 0; i < n_energy; ++i) {
        energy_grid[i] = E_min + (i * (E_max - E_min)) / (n_energy - 1);
    }
    
    baseline_grid.resize(n_baseline);
    for (size_t i = 0; i < n_baseline; ++i) {
        baseline_grid[i] = L_min + (i * (L_max - L_min)) / (n_baseline - 1);
    }
        
    // Allocate storage
    size_t total_size = n_energy * n_baseline;
    grid_P.resize(total_size);
    grid_dE.resize(total_size);
    grid_dL.resize(total_size);
    grid_d2EdL.resize(total_size);
    
    std::cout << "  Total grid points: " << total_size << std::endl;
    std::cout << "  Memory per variable: " << (total_size * sizeof(float) / 1e6) << " MB" << std::endl;
    
    // Compute grid values and derivatives
    for (size_t e = 0; e < n_energy; ++e) {
        double E = energy_grid[e];
                
        for (size_t l = 0; l < n_baseline; ++l) {
            double L = baseline_grid[l];
                    
            // Compute probability
            Eigen::Matrix3d prob = chic->compute_oscillations(E, L);
                    
            // Compute derivatives
            Eigen::Matrix3d dE = chic->compute_oscillations_derivatives("E", E, L);
            Eigen::Matrix3d dL = chic->compute_oscillations_derivatives("L", E, L);
                    
            // Approximate mixed partial via numerical differentiation
            Eigen::Matrix3d d2EdL = chic->compute_oscillations_derivatives("EL", E, L)
            // Eigen::Matrix3d d2EdL = chic->compute_oscillations_derivatives("E", E + 0.001, L) - 
            //                                chic->compute_oscillations_derivatives("E", E - 0.001, L);
            // d2EdL /= 0.002;
                    
            size_t idx = grid_index(m, f, e, l);
            grid_P[idx] = (float)prob(0, flavor);
            grid_dE[idx] = (float)dE(0, flavor);
            grid_dL[idx] = (float)dL(0, flavor);
            grid_d2EdL[idx] = (float)d2EdL(0, flavor);
        }
    }
}

/**
 * Find grid cell containing (energy, baseline)
 */
void CHIC_BATCH_GRID::find_grid_cell(double energy, double baseline) const {
    e_idx = -1;
    l_idx = -1;
    
    // Check bounds
    if (energy < E_min || energy > E_max) return;
    if (baseline < L_min || baseline > L_max) return;
    
    // Find energy index
    if (energy == energy_grid[n_energy - 1]) {
        e_idx = n_energy - 2;
    } else {
        for (size_t i = 0; i < n_energy - 1; ++i) {
            if (energy >= energy_grid[i] && energy < energy_grid[i + 1]) {
                e_idx = i;
                break;
            }
        }
    }
    
    // Find baseline index
    if (baseline == baseline_grid[n_baseline - 1]) {
        l_idx = n_baseline - 2;
    } else {
        for (size_t i = 0; i < n_baseline - 1; ++i) {
            if (baseline >= baseline_grid[i] && baseline < baseline_grid[i + 1]) {
                l_idx = i;
                break;
            }
        }
    }
}

/**
 * Hermite basis functions for 1D interpolation
 */
void CHIC_BATCH_GRID::hermite_basis_1d(double t, double& h00, double& h10, 
                                         double& h01, double& h11) {
    double t2 = t * t;
    double t3 = t2 * t;
    
    h00 = 1.0 - 3.0*t2 + 2.0*t3;
    h10 = t - 2.0*t2 + t3;
    h01 = 3.0*t2 - 2.0*t3;
    h11 = -t2 + t3;
}

/**
 * Bicubic Hermite interpolation in 2D
 */
double CHIC_BATCH_GRID::hermite_interpolate_2d(
    const std::string& mode,
    double energy,
    double baseline) const {
    
    if (!grid_ready) {
        throw std::runtime_error("Grid not ready. Call create_and_compute_grid() first.");
    }
    
    // Find grid cell
    int e_idx, l_idx;
    e_idx, l_idx = find_grid_cell(energy, baseline);
    
    if (e_idx < 0 || l_idx < 0) {
        throw std::out_of_range(
            "Energy=" + std::to_string(energy) + " GeV or Baseline=" + 
            std::to_string(baseline) + " km out of grid bounds"
        );
    }
    
    // Get spacing in this cell
    double dE = energy_grid[e_idx + 1] - energy_grid[e_idx];
    double dL = baseline_grid[l_idx + 1] - baseline_grid[l_idx];
    
    // Normalized coordinates [0,1] within cell
    double tE = (energy - energy_grid[e_idx]) / dE;
    double tL = (baseline - baseline_grid[l_idx]) / dL;
    
    // Hermite basis functions
    double hE00, hE10, hE01, hE11;
    double hL00, hL10, hL01, hL11;
    hermite_basis_1d(tE, hE00, hE10, hE01, hE11);
    hermite_basis_1d(tL, hL00, hL10, hL01, hL11);
    
    
    // Get corner values and derivatives
    Eigen::Matrix3d P00 = (*grid_data)[grid_index(e_idx, l_idx)];
    Eigen::Matrix3d P01 = (*grid_data)[grid_index(e_idx, l_idx + 1)];
    Eigen::Matrix3d P10 = (*grid_data)[grid_index(e_idx + 1, l_idx)];
    Eigen::Matrix3d P11 = (*grid_data)[grid_index(e_idx + 1, l_idx + 1)];
    
    Eigen::Matrix3d dE00 = grid_dE[grid_index(e_idx, l_idx)] * dE;
    Eigen::Matrix3d dE01 = grid_dE[grid_index(e_idx, l_idx + 1)] * dE;
    Eigen::Matrix3d dE10 = grid_dE[grid_index(e_idx + 1, l_idx)] * dE;
    Eigen::Matrix3d dE11 = grid_dE[grid_index(e_idx + 1, l_idx + 1)] * dE;
    
    Eigen::Matrix3d dL00 = grid_dL[grid_index(e_idx, l_idx)] * dL;
    Eigen::Matrix3d dL01 = grid_dL[grid_index(e_idx, l_idx + 1)] * dL;
    Eigen::Matrix3d dL10 = grid_dL[grid_index(e_idx + 1, l_idx)] * dL;
    Eigen::Matrix3d dL11 = grid_dL[grid_index(e_idx + 1, l_idx + 1)] * dL;
    
    Eigen::Matrix3d d2EdL00 = grid_d2EdL[grid_index(e_idx, l_idx)] * dE * dL;
    Eigen::Matrix3d d2EdL01 = grid_d2EdL[grid_index(e_idx, l_idx + 1)] * dE * dL;
    Eigen::Matrix3d d2EdL10 = grid_d2EdL[grid_index(e_idx + 1, l_idx)] * dE * dL;
    Eigen::Matrix3d d2EdL11 = grid_d2EdL[grid_index(e_idx + 1, l_idx + 1)] * dE * dL;
    
    // 2D Hermite interpolation
    Eigen::Matrix3d result = 0.0;
    result += hE00 * hL00 * P00;
    result += hE10 * hL00 * dE00;
    result += hE01 * hL00 * P10;
    result += hE11 * hL00 * dE10;
    
    result += hE00 * hL10 * dL00;
    result += hE10 * hL10 * d2EdL00;
    result += hE01 * hL10 * dL10;
    result += hE11 * hL10 * d2EdL10;
    
    result += hE00 * hL01 * P01;
    result += hE10 * hL01 * dE01;
    result += hE01 * hL01 * P11;
    result += hE11 * hL01 * dE11;
    
    result += hE00 * hL11 * dL01;
    result += hE10 * hL11 * d2EdL01;
    result += hE01 * hL11 * dL11;
    result += hE11 * hL11 * d2EdL11;
    
    return result;
}

/**
 * Interpolate single weight using precomputed grid
 */
double CHIC_BATCH_GRID::interpolate_weight(
    int flavor,
    double energy,
    double baseline) const {
    
    if (!grid_ready) {
        throw std::runtime_error("Grid not ready. Call create_and_compute_grid() first.");
    }
    
    if (initial_flux.empty()) {
        throw std::runtime_error("Initial flux not set");
    }
    
    // Get probability via Hermite interpolation
    Eigen::Matrix3d prob = hermite_interpolate_2d(energy, baseline, 0);
    
    // Apply initial flux (use first event's flux)
    double weight = 0.0;
    for (int flavor = 0; flavor < 3; ++flavor) {
        weight += prob * initial_flux[0][flavor];
    }
    
    return weight;
}

/**
 * Batch interpolation for multiple events
 */
std::vector<double> CHIC_BATCH_GRID::interpolate_weights_batch(
    const std::vector<std::tuple<std::string, int, double, double>>& events) const {
    
    if (!grid_ready) {
        throw std::runtime_error("Grid not ready. Call create_and_compute_grid() first.");
    }
    
    std::vector<double> results;
    results.reserve(events.size());
    
    for (size_t i = 0; i < events.size(); ++i) {
        const std::string& mode = std::get<0>(events[i]);
        int flavor = std::get<1>(events[i]);
        double energy = std::get<2>(events[i]);
        double baseline = std::get<3>(events[i]);
        
        double prob = hermite_interpolate_2d(energy, baseline, 0);
        
        // Use per-event flux if available
        double weight = 0.0;
        size_t flux_idx = (i < initial_flux.size()) ? i : 0;
        
        for (int f = 0; f < 3; ++f) {
            weight += prob * initial_flux[flux_idx][f];
        }
        
        results.push_back(weight);
    }
    
    return results;
}