/**
 * chic_profile.cpp  —  Fine-grained timing harness for CHIC / CHICDIFF
 *
 * Strategy
 * --------
 * We time every logical stage individually using std::chrono::steady_clock.
 * Each stage is repeated N_INNER times inside a hot loop (to amortise clock
 * overhead) and that loop is itself repeated N_OUTER times so we can report
 * mean ± stddev and spot jitter.  All results are printed as ns/call.
 *
 * Stages timed
 * ────────────
 *  Base class (CHIC)
 *   1. pmns_matrix()
 *   2. _set_vacuum()          (includes pmns + Hs0, Hs0_2, re_Hs0Vs0, detH0)
 *   3. _set_matter()
 *   4. _compute_hamiltonians()
 *   5. _exponential()
 *   6. compute_oscillations()  end-to-end (E + L sweep)
 *
 *  Derived class (CHICDIFF) — per-parameter derivative cost
 *   For each param in {density, dm221, dm231, th12, th13, th23, dcp, E}:
 *   7a. dHs builder only
 *   7b. _set_dHs() (builder + scaling + comm_dHH / comm_dHH2)
 *   7c. _amplitude_derivative()
 *   7d. compute_oscillations_derivatives()  end-to-end
 *
 * Compile (example):
 *   g++ -O2 -march=native -std=c++17 chic_profile.cpp CHIC.cpp \
 *       -I<eigen_path> -o chic_profile
 *   ./chic_profile
 *
 * If you have access to Linux perf / VTune, run with:
 *   perf stat -e cycles,instructions,cache-misses ./chic_profile
 */

#include "CHIC.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// ── tuneable knobs ──────────────────────────────────────────────────────────
static constexpr int    N_INNER   = 1'000;   // repetitions inside one timing block
static constexpr int    N_OUTER   = 50;      // independent timing blocks → stddev
static constexpr double E_GEV     = 1.0;     // neutrino energy [GeV]
static constexpr double L_KM      = 1300.0;  // baseline [km]
static constexpr double DENSITY   = 2.848;   // matter density [g/cm³]
// ────────────────────────────────────────────────────────────────────────────

// PDG 2023 central values
static constexpr double TH12  = 0.5836;
static constexpr double TH23  = 0.8587;
static constexpr double TH13  = 0.1503;
static constexpr double DCP   = 1.36 * M_PI;
static constexpr double DM221 = 7.53e-5;
static constexpr double DM231 = 2.453e-3;
static constexpr double YE    = 0.5;

// ── timing helpers ──────────────────────────────────────────────────────────
using Clock = std::chrono::steady_clock;
using NS    = std::chrono::duration<double, std::nano>;

struct Stats { double mean_ns, std_ns, min_ns; };

// Time a callable F() repeated N_INNER × N_OUTER times.
// Returns per-call statistics in nanoseconds.
template <typename F>
Stats measure(F&& f) {
    std::vector<double> samples(N_OUTER);
    for (int o = 0; o < N_OUTER; ++o) {
        auto t0 = Clock::now();
        for (int i = 0; i < N_INNER; ++i) f();
        auto t1 = Clock::now();
        samples[o] = NS(t1 - t0).count() / N_INNER;
    }
    double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / N_OUTER;
    double var  = 0.0;
    for (double s : samples) var += (s - mean) * (s - mean);
    var /= N_OUTER;
    double mn = *std::min_element(samples.begin(), samples.end());
    return {mean, std::sqrt(var), mn};
}

// Pretty printer
void print_row(const std::string& label, const Stats& s) {
    std::cout << std::left  << std::setw(42) << label
              << std::right << std::setw(10) << std::fixed << std::setprecision(2) << s.mean_ns
              << std::setw(10) << s.std_ns
              << std::setw(10) << s.min_ns
              << "\n";
}

void print_header(const std::string& section) {
    std::cout << "\n"
              << "══════════════════════════════════════════════════════════════════\n"
              << "  " << section << "\n"
              << "══════════════════════════════════════════════════════════════════\n"
              << std::left  << std::setw(42) << "Stage"
              << std::right << std::setw(10) << "mean ns"
              << std::setw(10) << "std ns"
              << std::setw(10) << "min ns"
              << "\n"
              << std::string(72, '─') << "\n";
}

// ── subclass that exposes protected internals for profiling ─────────────────
// We inherit CHICDIFF and add thin public wrappers around every private stage.
struct CHICProbe : public CHICDIFF {
    using CHICDIFF::CHICDIFF;

    // Base stages
    void pub_pmns()           { pmns_matrix(); }
    void pub_set_vacuum()     { _set_vacuum(); }
    void pub_set_matter()     { _set_matter(); }
    void pub_compute_H()      { _compute_hamiltonians(); }
    void pub_exponential()    { _exponential(); }

    // Force a state that makes _compute_hamiltonians and _exponential runnable
    void prime(double E, double L) {
        // Run the full pipeline once so all caches are valid
        _amplitude(E, L);
    }

    // Individual dHs builders (called after dHs.setZero() so state is clean)
    void pub_dHs_drho()   { dHs.setZero(); dHs_drho();   }
    void pub_dHs_ddm221() { dHs.setZero(); dHs_ddm221(); }
    void pub_dHs_ddm231() { dHs.setZero(); dHs_ddm231(); }
    void pub_dHs_dth12()  { dHs.setZero(); dHs_dth12();  }
    void pub_dHs_dth13()  { dHs.setZero(); dHs_dth13();  }
    void pub_dHs_dth23()  { dHs.setZero(); dHs_dth23();  }
    void pub_dHs_ddcp()   { dHs.setZero(); dHs_ddcp();   }
    void pub_dHs_dE()     { dHs.setZero(); dHs_dE();     }

    // Full _set_dHs (builder + scaling + comms)
    void pub_set_dHs(std::string_view p) { param0 = p; _set_dHs(); }

    // _amplitude_derivative alone (requires _set_dHs to have been called)
    void pub_amplitude_derivative() { _amplitude_derivative(); }
};

// ── main ────────────────────────────────────────────────────────────────────
int main() {
    std::cout << "CHIC performance profile\n"
              << "  N_INNER=" << N_INNER << "  N_OUTER=" << N_OUTER
              << "  E=" << E_GEV << " GeV  L=" << L_KM << " km"
              << "  density=" << DENSITY << " g/cm³\n";

    CHICProbe p("nu", TH12, TH23, TH13, DCP, DM221, DM231, DENSITY, YE);
    p.prime(E_GEV, L_KM);

    // ── 1. Base-class pipeline stages ───────────────────────────────────────
    print_header("Base class  CHIC — pipeline stages");

    print_row("pmns_matrix()",
        measure([&]{ p.pub_pmns(); }));

    print_row("_set_vacuum()  (pmns + Hs0 + Hs0_2 + re_Hs0Vs0 + detH0)",
        measure([&]{ p.pub_set_vacuum(); }));

    print_row("_set_matter()",
        measure([&]{ p.pub_set_matter(); }));

    // _compute_hamiltonians needs E0 set; it's already set by prime()
    print_row("_compute_hamiltonians()  (Hs, Hs2, eigenvalues, lambdas)",
        measure([&]{ p.pub_compute_H(); }));

    print_row("_exponential()  (exp_lambdas, cJ, J)",
        measure([&]{ p.pub_exponential(); }));

    // End-to-end oscillation probability (fixed E, fixed L — exercises caching)
    {
        double E = E_GEV, L = L_KM;
        print_row("compute_oscillations()  [cached E+L — only _exponential]",
            measure([&]{ volatile auto r = p.compute_oscillations(E, L); (void)r; }));
    }
    {
        // Vary E slightly each call to defeat the E-cache and force full recompute
        double E = E_GEV;
        int    call = 0;
        print_row("compute_oscillations()  [new E each call — full pipeline]",
            measure([&]{
                E = E_GEV + 1e-9 * (call++ & 0xFFF);   // tiny variation, same physics
                volatile auto r = p.compute_oscillations(E, L_KM);
                (void)r;
            }));
    }

    // ── 2. dHs builders (isolated) ──────────────────────────────────────────
    print_header("CHICDIFF — dHs builders  (isolated, after setZero)");

    print_row("dHs_drho()",   measure([&]{ p.pub_dHs_drho();   }));
    print_row("dHs_ddm221()", measure([&]{ p.pub_dHs_ddm221(); }));
    print_row("dHs_ddm231()", measure([&]{ p.pub_dHs_ddm231(); }));
    print_row("dHs_dth12()",  measure([&]{ p.pub_dHs_dth12();  }));
    print_row("dHs_dth13()",  measure([&]{ p.pub_dHs_dth13();  }));
    print_row("dHs_dth23()",  measure([&]{ p.pub_dHs_dth23();  }));
    print_row("dHs_ddcp()",   measure([&]{ p.pub_dHs_ddcp();   }));
    print_row("dHs_dE()",     measure([&]{ p.pub_dHs_dE();     }));

    // ── 3. _set_dHs (builder + scaling + comm_dHH + comm_dHH2) ─────────────
    print_header("CHICDIFF — _set_dHs()  (builder + scale + 2 × matmul pairs)");

    for (const char* param : {"density","dm221","dm231","th12","th13","th23","dcp","E"}) {
        std::string label = std::string("_set_dHs(\"") + param + "\")";
        print_row(label, measure([&]{ p.pub_set_dHs(param); }));
    }

    // ── 4. _amplitude_derivative  (after a representative _set_dHs) ─────────
    print_header("CHICDIFF — _amplitude_derivative()  (I_ijk + cdJ + dJ assembly)");

    // Pre-set dHs for each parameter then time _amplitude_derivative alone
    for (const char* param : {"density","dm221","dm231","th12","th13","th23","dcp","E"}) {
        // Warm up dHs / comm matrices for this param
        p.pub_set_dHs(param);
        std::string label = std::string("_amplitude_derivative()  [dHs=") + param + "]";
        print_row(label,
            measure([&]{
                // Re-set dHs so the function sees consistent inputs each iteration
                p.pub_set_dHs(param);
                p.pub_amplitude_derivative();
            }));
    }

    // ── 5. End-to-end derivative (full pipeline per param) ──────────────────
    print_header("CHICDIFF — compute_oscillations_derivatives()  end-to-end");

    {
        double E = E_GEV, L = L_KM;
        for (const char* param : {"density","dm221","dm231","th12","th13","th23","dcp","E"}) {
            std::string label = std::string("param=\"") + param + "\"  [cached E+L]";
            print_row(label,
                measure([&]{
                    volatile auto r = p.compute_oscillations_derivatives(param, E, L);
                    (void)r;
                }));
        }
    }
    {
        // Force full recompute by varying E
        int call = 0;
        double E = E_GEV;
        for (const char* param : {"density","dm221","dm231","th12","th13","th23","dcp","E"}) {
            call = 0;
            std::string label = std::string("param=\"") + param + "\"  [new E each call]";
            print_row(label,
                measure([&]{
                    E = E_GEV + 1e-9 * (call++ & 0xFFF);
                    volatile auto r = p.compute_oscillations_derivatives(param, E, L_KM);
                    (void)r;
                }));
        }
    }

    // ── 6. Cost breakdown summary ────────────────────────────────────────────
    std::cout << "\n"
              << "══════════════════════════════════════════════════════════════════\n"
              << "  Cost breakdown guide\n"
              << "══════════════════════════════════════════════════════════════════\n"
              << "  _set_dHs − dHs_builder   ≈ cost of 2× matmul pairs (comm_dHH)\n"
              << "  _amplitude_derivative     ≈ I_ijk assembly + 6 matmul terms\n"
              << "  _compute_hamiltonians     ≈ Hs, Hs2 build + trig eigenvalues\n"
              << "  _exponential              ≈ 3× std::polar + cJ dots + J build\n"
              << "\n"
              << "  Bottleneck signals:\n"
              << "    • If _compute_hamiltonians >> _exponential: trig (acos/sqrt)\n"
              << "      is dominating — consider a lookup or approximation.\n"
              << "    • If _amplitude_derivative >> _set_dHs: the 6-term dJ\n"
              << "      assembly (matmul chains) is the bottleneck.\n"
              << "    • If comm_dHH cost (_set_dHs − builder) is large relative\n"
              << "      to _amplitude_derivative: consider the commented-out\n"
              << "      adjoint shortcut (saves one matmul per pair).\n"
              << "    • If std ns is large (> 20% of mean): cache pressure or\n"
              << "      OS scheduling — pin the process and disable turbo.\n"
              << "\n";

    return 0;
}
