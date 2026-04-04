#include "CHICEARTH.h"
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

// ------------------------------------------------------------------ //
//  Helpers
// ------------------------------------------------------------------ //

static void print_result(const char* label, int N, double total_ms) {
  const double per_us = total_ms * 1e3 / N;
  const double per_ns = total_ms * 1e6 / N;
  std::cout << std::left  << std::setw(56) << label
            << std::right
            << std::setw(8)  << N           << " calls  |"
            << std::setw(10) << std::fixed  << std::setprecision(2) << total_ms << " ms  |"
            << std::setw(8)  << std::setprecision(2) << per_us    << " us/call  |"
            << std::setw(8)  << std::setprecision(0) << per_ns    << " ns/call\n";
}

// cos_zenith threshold below which all 4 PREM layers are crossed.
// The innermost shell has radius 1221.5 km, Earth radius 6371 km.
// cos_zenith_max = -sqrt(1 - (R_inner/R_earth)^2) ~ -0.9816
static constexpr double CZ_CORE = -0.9816;

using hrclock = std::chrono::high_resolution_clock;

// ------------------------------------------------------------------ //
//  A) Fixed energy, varying cos_zenith — only core-crossing tracks
// ------------------------------------------------------------------ //
static void bench_fixed_energy(CHICEARTH& earth,
                                const std::vector<double>& czs_core) {
  const double E_values[] = {1.0, 10.0, 100.0};
  double sink = 0.0;

  for (double E : E_values) {
    earth.compute_oscillations(E, czs_core[0]); // warm-up / cache H

    auto t0 = hrclock::now();
    for (double cz : czs_core)
      sink += earth.compute_oscillations(E, cz).sum();
    auto t1 = hrclock::now();

    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    char label[80];
    std::snprintf(label, sizeof(label),
                  "A) Fixed E=%.0f GeV, core-crossing cos_z (cached H)", E);
    print_result(label, (int)czs_core.size(), ms);
  }
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  B) Varying energy AND cos_zenith — core-crossing only, full cost
// ------------------------------------------------------------------ //
static void bench_full_grid(CHICEARTH& earth,
                             const std::vector<double>& energies,
                             const std::vector<double>& czs_core) {
  double sink = 0.0;
  earth.compute_oscillations(energies[0], czs_core[0]); // warm-up

  auto t0 = hrclock::now();
  for (double E : energies)
    for (double cz : czs_core)
      sink += earth.compute_oscillations(E, cz).sum();
  auto t1 = hrclock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  print_result("B) Varying E + core-crossing cos_z (full cost)",
               (int)(energies.size() * czs_core.size()), ms);
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  C) Fixed cos_zenith, varying energy — track cached, H recomputed
// ------------------------------------------------------------------ //
static void bench_fixed_cz(CHICEARTH& earth,
                            const std::vector<double>& energies,
                            double cz) {
  double sink = 0.0;
  earth.compute_oscillations(energies[0], cz); // warm-up

  auto t0 = hrclock::now();
  for (double E : energies)
    sink += earth.compute_oscillations(E, cz).sum();
  auto t1 = hrclock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  char label[80];
  std::snprintf(label, sizeof(label),
                "C) Fixed cos_z=%.2f, varying E (track cached, H recomputed)", cz);
  print_result(label, (int)energies.size(), ms);
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  D) Scan over dm231 — fixed E and cos_zenith, varying dm231
// ------------------------------------------------------------------ //
static void bench_dm231(CHICEARTH& earth,
                         const std::vector<double>& dm231_values,
                         double E, double cz) {
  double sink = 0.0;
  earth.update_dm231(dm231_values[0]);
  earth.compute_oscillations(E, cz); // warm-up

  auto t0 = hrclock::now();
  for (double dm : dm231_values) {
    earth.update_dm231(dm);
    sink += earth.compute_oscillations(E, cz).sum();
  }
  auto t1 = hrclock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  print_result("D) Varying dm231, fixed E=10 GeV cos_z=-1",
               (int)dm231_values.size(), ms);
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  D) Scan over theta_23 — fixed E and cos_zenith, varying theta_23
// ------------------------------------------------------------------ //
static void bench_th23(CHICEARTH& earth,
                        const std::vector<double>& th23_values,
                        double E, double cz) {
  double sink = 0.0;
  earth.update_th23(th23_values[0]);
  earth.compute_oscillations(E, cz); // warm-up

  auto t0 = hrclock::now();
  for (double th : th23_values) {
    earth.update_th23(th);
    sink += earth.compute_oscillations(E, cz).sum();
  }
  auto t1 = hrclock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  print_result("E) Varying theta_23, fixed E=10 GeV cos_z=-1",
               (int)th23_values.size(), ms);
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  E) Joint scan over dm231 x theta_23 — core-crossing trajectory
// ------------------------------------------------------------------ //
static void bench_dm231_th23(CHICEARTH& earth,
                               const std::vector<double>& dm231_values,
                               const std::vector<double>& th23_values,
                               double E, double cz) {
  double sink = 0.0;
  earth.update_dm231(dm231_values[0]);
  earth.update_th23(th23_values[0]);
  earth.compute_oscillations(E, cz); // warm-up

  auto t0 = hrclock::now();
  for (double dm : dm231_values) {
    earth.update_dm231(dm);
    for (double th : th23_values) {
      earth.update_th23(th);
      sink += earth.compute_oscillations(E, cz).sum();
    }
  }
  auto t1 = hrclock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  print_result("F) dm231 x theta_23 grid, fixed E=10 GeV cos_z=-1",
               (int)(dm231_values.size() * th23_values.size()), ms);
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  main
// ------------------------------------------------------------------ //
int main() {
  // Default NuFIT 5.2 best-fit parameters (NO)
  CHICEARTH earth("neutrino");

  const int NE   = 2000;
  const int Ncz  = 2000;
  const int Nparam = 200;

  const double log_E_min = std::log10(0.1);
  const double log_E_max = std::log10(100.0);
  const double cz_min    = -1.0;

  // Energies: log-spaced
  std::vector<double> energies(NE);
  for (int ie = 0; ie < NE; ++ie)
    energies[ie] = std::pow(10.0, log_E_min +
                            (log_E_max - log_E_min) * ie / (NE - 1));

  // cos_zenith: only core-crossing (below CZ_CORE threshold)
  std::vector<double> czs_core;
  for (int icz = 0; icz < Ncz; ++icz) {
    const double cz = cz_min + (CZ_CORE - cz_min) * icz / (Ncz - 1);
    czs_core.push_back(cz);
  }

  // dm231: range 2.3e-3 to 2.7e-3 eV^2 (3σ interval around best-fit)
  std::vector<double> dm231_values(Nparam);
  for (int i = 0; i < Nparam; ++i)
    dm231_values[i] = 2.3e-3 + (2.7e-3 - 2.3e-3) * i / (Nparam - 1);

  // theta_23: range 0.6 to 1.0 rad (covers full octant range)
  std::vector<double> th23_values(Nparam);
  for (int i = 0; i < Nparam; ++i)
    th23_values[i] = 0.6 + (1.0 - 0.6) * i / (Nparam - 1);

  const double E_fixed  = 10.0;  // GeV
  const double cz_fixed = -1.0;  // straight down, core-crossing

  const std::string sep(104, '-');
  std::cout << "\n" << sep << "\n"
            << " CHICEARTH performance benchmark  "
            << "(core-crossing tracks only, cos_zenith < " << CZ_CORE << ")\n"
            << sep << "\n";

  bench_fixed_energy(earth, czs_core);
  std::cout << sep << "\n";
  bench_full_grid(earth, energies, czs_core);
  std::cout << sep << "\n";
  bench_fixed_cz(earth, energies, -1.0);
  bench_fixed_cz(earth, energies, -0.99);  // just inside core
  bench_fixed_cz(earth, energies, CZ_CORE + 0.001); // just outside core
  std::cout << sep << "\n";
  bench_dm231(earth, dm231_values, E_fixed, cz_fixed);
  bench_th23(earth, th23_values, E_fixed, cz_fixed);
  std::cout << sep << "\n";
  bench_dm231_th23(earth, dm231_values, th23_values, E_fixed, cz_fixed);
  std::cout << sep << "\n\n";

  return 0;
}