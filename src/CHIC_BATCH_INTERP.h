#ifndef CHIC_BATCH_GRID_H
#define CHIC_BATCH_GRID_H

#include "CHIC_BATCH.h"
#include <vector>
#include <map>
#include <string>


class CHIC_BATCH_GRID : public CHIC_BATCH {
public:
    // Constructor: enables/disables derivative calculations
    explicit CHIC_BATCH_GRID();
    
    // Destructor
    ~CHIC_BATCH_GRID() override = default;

    // Create and precompute oscillation grid with Hermite derivatives
    void create_and_compute_grid(
        int n_energy, double E_min, double E_max,
        int n_baseline, double L_min, double L_max
    );

    // Interpolate oscillation weight at single point using precomputed grid
    double interpolate_weight(
        int flavor,
        double energy,
        double baseline
    );

    // Batch interpolation using precomputed grid (MC-style)
    std::vector<double> compute_oscillations(
        const std::vector<int>& flavors,
        const std::vector<double>& energies,
        const std::vector<double>& baselines);

    // Check if grid is precomputed and ready for interpolation
    bool is_grid_ready() const noexcept { return grid_ready; }


private:
    // ================================================================
    // Grid data storage
    // ================================================================
    bool grid_ready = false;
    
    // Grid axes
    std::vector<double> energy_grid;                // Energy points
    std::vector<double> baseline_grid;              // Baseline points
    
    // Grid dimensions
    size_t n_energy = 0, n_baseline = 0, n_modes = 0, n_flavors = 0;
    double E_min = 0, E_max = 0, L_min = 0, L_max = 0;
    
    // Grid data storage (compact float32 format)
    // Index: [mode_idx][flavor_idx][energy_idx][baseline_idx]
    // Flattened to 1D: idx = m*nf*ne*nl + f*ne*nl + e*nl + l
    std::vector<Eigen::Matrix3d> grid_P;      // Probability values
    std::vector<Eigen::Matrix3d> grid_dE;     // ∂P/∂E (energy derivative)
    std::vector<Eigen::Matrix3d> grid_dL;     // ∂P/∂L (baseline derivative)
    std::vector<Eigen::Matrix3d> grid_d2EdL;  // ∂²P/∂E∂L (mixed partial)
    
    // Mode and flavor lookup
    std::map<std::string, size_t> mode_index;
    std::map<int, size_t> flavor_index;

    // ================================================================
    // Private helper methods
    // ================================================================
    
    /**
     * Get flat index in grid arrays
     * Index: [mode_idx * nf*ne*nl + flavor_idx * ne*nl + energy_idx * nl + baseline_idx]
     */
    inline size_t grid_index(size_t m, size_t f, size_t e, size_t l) const noexcept {
        return m * n_flavors * n_energy * n_baseline + 
               f * n_energy * n_baseline + 
               e * n_baseline + l;
    }
    
    // Bicubic Hermite interpolation in 1D
    static void hermite_basis_1d(double t, double& h00, double& h10, 
                                  double& h01, double& h11);
    
    // Bicubic Hermite interpolation in 2D (E, L)
    double hermite_interpolate_2d(
        int flavor,
        double energy,
        double baseline,
        int grid_var
    );
    
    // Find grid cell containing (energy, baseline)
    void find_grid_cell(double energy, double baseline,
                       int& e_idx, int& l_idx) const;
};

#endif // CHIC_BATCH_GRID_H