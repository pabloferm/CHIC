#include "CHIC_BATCH.h"
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>

// Simple MC event structure
struct MCEvent {
    int event_id;
    int true_flavor;            // PDG
    double true_energy;         // GeV
    double baseline;            // km
    std::vector<double> flux;   // [e, mu, tau] composition
    double oscillation_weight;  // Output: computed weight
};

int main() {
    std::cout << "=== CHIC_BATCH MC Example ===" << std::endl;
    std::cout << "Processing simulated oscillation events\n" << std::endl;

    // ============================================================
    // Create batch processor
    // ============================================================
    CHIC_BATCH batch(true);  // Derivatives disabled for speed
    batch.update_parameter("theta_23", 0.85);
    std::cout << "✓ Created CHIC_BATCH for MC processing\n" << std::endl;

    // ============================================================
    // Generate simulated MC events
    // ============================================================
    std::vector<MCEvent> mc_events;
    std::mt19937 rng(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<double> energy_dist(0.5, 5.0);
    
    // Simulate 1000 events
    const int n_events = 1000000;
    
    std::cout << "Generating " << n_events << " simulated events..." << std::endl;
    for (int i = 0; i < n_events; ++i) {
        MCEvent evt;
        evt.event_id = i;
        
        // Randomly assign flavor (equal distribution)
        evt.true_flavor = 2 * (rng() % 3 + 6);
        std::cout << evt.true_flavor << std::endl;
        
        // Random energy
        evt.true_energy = energy_dist(rng);
        
        // Baseline: T2K (295 km) or NOvA (810 km)
        evt.baseline = (rng() % 2 == 0) ? 295.0 : 810.0;
        
        // Random initial flux composition
        double f_e = rng() % 100 / 100.0;
        double f_mu = (rng() % (int)(100 - f_e * 100)) / 100.0;
        double f_tau = 1.0 - f_e - f_mu;
        evt.flux = {f_e, f_mu, f_tau};
        
        mc_events.push_back(evt);
    }
    std::cout << "✓ Generated " << mc_events.size() << " events\n" << std::endl;

    // ============================================================
    // Extract arrays for batch processing
    // ============================================================
    std::vector<int> flavors;
    std::vector<double> energies;
    std::vector<double> baselines;
    std::vector<std::vector<double>> per_event_fluxes;
    
    for (const auto& evt : mc_events) {
        flavors.push_back(evt.true_flavor);
        energies.push_back(evt.true_energy);
        baselines.push_back(evt.baseline);
        per_event_fluxes.push_back(evt.flux);
    }
    
    // Set per-event flux
    batch.set_initial_flux(per_event_fluxes);
    std::cout << "✓ Loaded per-event flux for all events\n" << std::endl;

    // ============================================================
    // Compute oscillation weights (all events at once!)
    // ============================================================
    std::cout << "Computing oscillation weights for all events..." << std::endl;
    std::vector<double> weights = batch.compute_oscillations(
        flavors,
        energies,
        baselines
    );
    std::cout << "✓ Computed " << weights.size() << " oscillation weights\n" << std::endl;

    // ============================================================
    // Store results back into MC events
    // ============================================================
    for (size_t i = 0; i < mc_events.size(); ++i) {
        mc_events[i].oscillation_weight = weights[i];
    }

    // ============================================================
    // Analyze results
    // ============================================================
    std::cout << "\n=== Event Statistics ===" << std::endl;
    
    double sum_weights = 0.0;
    double max_weight = 0.0;
    double min_weight = 1.0;
    int nu_count = 0, nubar_count = 0;
    
    for (const auto& evt : mc_events) {
        if (evt.true_flavor < 0) nubar_count += 1;
        else nu_count += 1;
        sum_weights += evt.oscillation_weight;
        max_weight = std::max(max_weight, evt.oscillation_weight);
        min_weight = std::min(min_weight, evt.oscillation_weight);
    }
    
    std::cout << "Neutrino events:     " << nu_count << std::endl;
    std::cout << "Antineutrino events: " << nubar_count << std::endl;
    std::cout << "Sum of weights:      " << std::fixed << std::setprecision(6) 
              << sum_weights << std::endl;
    std::cout << "Mean weight:         " << sum_weights / mc_events.size() << std::endl;
    std::cout << "Min weight:          " << min_weight << std::endl;
    std::cout << "Max weight:          " << max_weight << std::endl;

    // ============================================================
    // Show sample events
    // ============================================================
    std::cout << "\n=== Sample Events ===" << std::endl;
    std::cout << std::left << std::setw(8) << "ID"
              << std::setw(10) << "Flavor"
              << std::setw(12) << "Energy(GeV)"
              << std::setw(12) << "Baseline(km)"
              << std::setw(12) << "Weight" << std::endl;
    std::cout << std::string(66, '-') << std::endl;
    
    for (int i = 0; i < std::min(10, (int)mc_events.size()); ++i) {
        const auto& evt = mc_events[i];
        std::string flavor_name = (std::abs(evt.true_flavor) == 12) ? "e" : 
                                  (std::abs(evt.true_flavor) == 14) ? "mu" : "tau";
        
        std::cout << std::left << std::setw(8) << evt.event_id
                  << std::setw(10) << flavor_name
                  << std::fixed << std::setprecision(3)
                  << std::setw(12) << evt.true_energy
                  << std::setw(12) << evt.baseline
                  << std::setw(12) << evt.oscillation_weight << std::endl;
    }
    std::cout << "... and " << (mc_events.size() - 10) << " more events\n" << std::endl;

    // ============================================================
    // MC Parameter scan
    // ============================================================
    std::cout << "=== Parameter Scan ===" << std::endl;
    std::cout << "Scanning theta_23 effect on oscillation weights\n" << std::endl;
    
    std::vector<double> theta_23_values = {0.7, 0.75, 0.775, 0.8, 0.85};
    
    for (double theta_23 : theta_23_values) {
        batch.update_parameter("theta_23", theta_23);
        
        auto scan_weights = batch.compute_oscillations(
            flavors, energies, baselines
        );
        
        double scan_sum = 0.0;
        for (double w : scan_weights) sum_weights += w;
        
        std::cout << "θ₂₃ = " << std::fixed << std::setprecision(3) 
                  << theta_23 << " → Sum = " << scan_sum << std::endl;
    }

    return 0;
}
