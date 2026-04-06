#include "CHIC_EARTH.h"
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

// ------------------------------------------------------------------ //
//  Helpers
// ------------------------------------------------------------------ //

bool check_unitarity(const Eigen::Matrix3d& P, double tol = 1e-6) {
  bool ok = true;
  for (int i = 0; i < 3; ++i) {
    const double row_sum = P.row(i).sum();
    const double col_sum = P.col(i).sum();
    if (std::abs(row_sum - 1.0) > tol) {
      std::cerr << "  [FAIL] row " << i << " sum = " << row_sum << "\n";
      ok = false;
    }
    if (std::abs(col_sum - 1.0) > tol) {
      std::cerr << "  [FAIL] col " << i << " sum = " << col_sum << "\n";
      ok = false;
    }
  }
  return ok;
}

static const char* flavour[3] = {"e", "mu", "tau"};

void print_matrix(const Eigen::Matrix3d& P) {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "         -> e        -> mu       -> tau\n";
  for (int a = 0; a < 3; ++a) {
    std::cout << "  nu_" << flavour[a] << " ";
    for (int b = 0; b < 3; ++b)
      std::cout << std::setw(12) << P(a, b);
    std::cout << "\n";
  }
}

// ------------------------------------------------------------------ //
//  1. Spot-check
// ------------------------------------------------------------------ //
void test_spot_checks(CHICEARTH& earth) {
  std::cout << "\n========== Spot checks ==========\n";

  const double cases[][2] = {
      {1.0,  -1.0},
      {1.0,  -0.5},
      {10.0, -1.0},
      {10.0, -0.5},
      {10.0, -1.0},
  };

  bool all_ok = true;
  for (auto& c : cases) {
    const double E          = c[0];
    const double cos_zenith = c[1];
    Eigen::Matrix3d P = earth.compute_oscillations(E, cos_zenith);
    std::cout << "\nE = " << E << " GeV,  cos_zenith = " << cos_zenith << "\n";
    print_matrix(P);
    const bool ok = check_unitarity(P);
    std::cout << "  Unitarity: " << (ok ? "PASS" : "FAIL") << "\n";
    all_ok &= ok;
  }
  std::cout << "\nSpot-check unitarity overall: " << (all_ok ? "PASS" : "FAIL") << "\n";
}

// ------------------------------------------------------------------ //
//  2. Unitarity sweep
// ------------------------------------------------------------------ //
void test_unitarity_grid(CHICEARTH& earth, int NE = 40, int Ncz = 40,
                         double tol = 1e-6) {
  std::cout << "\n========== Unitarity grid sweep ("
            << NE << " x " << Ncz << " points) ==========\n";

  const double log_E_min = std::log10(0.1);
  const double log_E_max = std::log10(100.0);
  const double cz_min    = -1.0;
  const double cz_max    = -0.01;

  int failures = 0;
  for (int ie = 0; ie < NE; ++ie) {
    const double E = std::pow(10.0, log_E_min +
                              (log_E_max - log_E_min) * ie / (NE - 1));
    for (int icz = 0; icz < Ncz; ++icz) {
      const double cz = cz_min + (cz_max - cz_min) * icz / (Ncz - 1);
      Eigen::Matrix3d P = earth.compute_oscillations(E, cz);
      if (!check_unitarity(P, tol)) {
        ++failures;
        std::cerr << "  Unitarity FAIL at E=" << E
                  << " GeV, cos_zenith=" << cz << "\n";
      }
    }
  }
  std::cout << "Failures: " << failures << " / " << NE * Ncz << "  -> "
            << (failures == 0 ? "PASS" : "FAIL") << "\n";
}

// ------------------------------------------------------------------ //
//  3. Performance benchmark
//
//  Two scenarios timed separately:
//
//  A) Fixed energy, varying cos_zenith:
//     Hamiltonians are cached — measures only _build_track + _amplitude.
//     Typical use case: scanning zenith at a fixed energy bin.
//
//  B) Varying energy AND cos_zenith:
//     Full cost per call including _compute_hamiltonians.
//     Worst-case / uncached scenario.
// ------------------------------------------------------------------ //
void test_performance(CHICEARTH& earth, int NE = 200, int Ncz = 200) {
  using clock = std::chrono::high_resolution_clock;

  const double log_E_min = std::log10(0.1);
  const double log_E_max = std::log10(100.0);
  const double cz_min    = -1.0;
  const double cz_max    = -0.01;

  // Pre-build grids — keep log/pow out of the timed region
  std::vector<double> energies(NE), czs(Ncz);
  for (int ie  = 0; ie  < NE;  ++ie)
    energies[ie] = std::pow(10.0, log_E_min + (log_E_max - log_E_min) * ie / (NE - 1));
  for (int icz = 0; icz < Ncz; ++icz)
    czs[icz] = cz_min + (cz_max - cz_min) * icz / (Ncz - 1);

  // Accumulator — prevents the compiler from dead-store-eliminating the calls
  double sink = 0.0;

  // ---- A) Fixed energy, varying cos_zenith ----
  {
    const double E_fixed = energies[NE / 2];
    earth.compute_oscillations(E_fixed, czs[0]); // warm-up: cache Hamiltonians

    auto t0 = clock::now();
    for (int icz = 0; icz < Ncz; ++icz)
      sink += earth.compute_oscillations(E_fixed, czs[icz]).sum();
    auto t1 = clock::now();

    const double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double per_us   = std::chrono::duration<double, std::micro>(t1 - t0).count() / Ncz;
    std::cout << "\n========== Performance A: fixed E=" << E_fixed
              << " GeV, varying cos_zenith ==========\n"
              << "  Trajectories  : " << Ncz << "\n"
              << "  Total         : " << std::fixed << std::setprecision(3)
              << total_ms << " ms\n"
              << "  Per trajectory: " << std::setprecision(3) << per_us << " us\n";
  }

  // ---- B) Varying energy AND cos_zenith ----
  {
    earth.compute_oscillations(energies[0], czs[0]); // warm-up

    auto t0 = clock::now();
    for (int ie = 0; ie < NE; ++ie)
      for (int icz = 0; icz < Ncz; ++icz)
        sink += earth.compute_oscillations(energies[ie], czs[icz]).sum();
    auto t1 = clock::now();

    const int    Ntraj    = NE * Ncz;
    const double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double per_us   = std::chrono::duration<double, std::micro>(t1 - t0).count() / Ntraj;
    std::cout << "\n========== Performance B: varying E and cos_zenith ==========\n"
              << "  Trajectories  : " << Ntraj << "\n"
              << "  Total         : " << std::fixed << std::setprecision(3)
              << total_ms << " ms\n"
              << "  Per trajectory: " << std::setprecision(3) << per_us << " us\n";
  }

  // Prevent sink from being optimised away
  if (sink == 0.0) std::cout << "";
}

// ------------------------------------------------------------------ //
//  4. Dump grid + gnuplot script
// ------------------------------------------------------------------ //
void dump_and_plot(CHICEARTH& earth,
                   int NE = 800, int Ncz = 50,
                   const std::string& datafile  = "oscillations.dat",
                   const std::string& scriptfile = "plot_oscillations.gp") {
  std::cout << "\n========== Generating gnuplot data: " << datafile << " ==========\n";

  const double log_E_min = std::log10(0.1);
  const double log_E_max = std::log10(100.0);
  const double cz_min    = -1.0;
  const double cz_max    =  0.0;

  std::ofstream dat(datafile);
  if (!dat) { std::cerr << "Cannot open " << datafile << "\n"; return; }

  dat << "# E[GeV]  cos_zenith  "
         "Pee  Pemu  Petau  "
         "Pmue  Pmumu  Pmutau  "
         "Ptaue  Ptaumu  Ptautau\n";

  for (int icz = 0; icz < Ncz; ++icz) {
    const double cz = cz_min + (cz_max - cz_min) * icz / (Ncz - 1);
    if (icz > 0) dat << "\n";
    for (int ie = 0; ie < NE; ++ie) {
      const double E = std::pow(10.0, log_E_min +
                                (log_E_max - log_E_min) * ie / (NE - 1));
      Eigen::Matrix3d P = earth.compute_oscillations(E, cz);
      dat << std::scientific << std::setprecision(6) << E << "  " << cz;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          dat << "  " << P(a, b);
      dat << "\n";
    }
  }
  dat.close();
  std::cout << "Data written to " << datafile << "\n";

  std::ofstream gp(scriptfile);
  if (!gp) { std::cerr << "Cannot open " << scriptfile << "\n"; return; }

  const char* titles[3][3] = {
    {"P_{ee}",    "P_{e#mu}",    "P_{e#tau}"},
    {"P_{#mue}",  "P_{#mu#mu}",  "P_{#mu#tau}"},
    {"P_{#taue}", "P_{#tau#mu}", "P_{#tau#tau}"},
  };

  gp << "set terminal pngcairo size 1200,1100 enhanced font 'Arial,11'\n"
     << "set output 'oscillations.png'\n"
     << "set multiplot layout 3,3 title 'Neutrino oscillation probabilities (Earth)'\n"
     << "set xlabel 'E [GeV]'\n"
     << "set ylabel 'cos(zenith)'\n"
     << "set cbrange [0:1]\n"
     << "set palette rgbformulae 33,13,10\n"
     << "set logscale x\n"
     << "set xrange [0.1:100]\n"
     << "set yrange [-1:0]\n";

  int col = 3;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b, ++col)
      gp << "set title '" << titles[a][b] << "'\n"
         << "plot '" << datafile << "' using 1:2:" << col << " with image notitle\n";

  gp << "unset multiplot\n";
  gp.close();
  std::cout << "Gnuplot script written to " << scriptfile << "\n"
            << "Run with:  gnuplot " << scriptfile << "\n";
}

// ------------------------------------------------------------------ //
//  main
// ------------------------------------------------------------------ //
int main() {
  CHICEARTH earth("neutrino");

  test_spot_checks(earth);
  test_unitarity_grid(earth);
  test_performance(earth);
  dump_and_plot(earth);

  return 0;
}