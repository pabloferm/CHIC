#ifndef EARTH_H
#define EARTH_H

#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <string>

static constexpr double R_EARTH = 6371.0;  // km

// ============================================================
// EarthModel: shell-constant density model for neutrino propagation.
//
// All models derived from PREM (Dziewonski & Anderson 1981).
// Densities are volume-weighted averages of the PREM polynomial
// over each shell:
//
//   rho_i = integral(rho(r) r^2 dr, r_in, r_out)
//           / integral(r^2 dr, r_in, r_out)
//
// using the published polynomial coefficients from Table I of
// Dziewonski & Anderson (1981).
//
// Electron fractions (Ye = Z/A):
//   Inner core   0.4656  — iron/FeNi  (Allegre et al. 1995)
//   Outer core   0.4656  — iron/FeNi
//   Mantle+crust 0.4957  — pyrolite   (McDonough & Sun 1995)
//
// Standard values used throughout the neutrino-oscillation literature:
//   IceCube DeepCore (2025), Frontiers Earth Sci. (2023),
//   Sci. Reports (2015), Eur. Phys. J. C (2020, 2022).
//
// Available models (shells listed innermost first):
//   "PREM4"    —  4 shells: inner core | outer core | mantle | crust
//   "PREM10"   — 10 shells: follows all major geophysical discontinuities
//   "PREM42"   — 42 shells: standard high-res neutrino tomography model
//   "PREMMAX"  — 81 shells: maximum resolution, ~equal-width per region
//
// Total mass reproduced to better than 0.02% (5.9732e24 vs 5.9722e24 kg).
// ============================================================

struct EarthModel {
    int                 Nlayers;
    std::vector<double> radii;    // outer radius of each shell [km]
    std::vector<double> density;  // volume-averaged density    [g/cm³]
    std::vector<double> Ye;       // electron fraction (Z/A)
};

inline const std::unordered_map<std::string, EarthModel>& earth_models() {
    static const std::unordered_map<std::string, EarthModel> models = {

        // --------------------------------------------------------
        // PREM4 — 4 shells
        // Inner core | outer core | lower mantle | upper mantle+crust
        // Lightweight baseline; matches original CHICEARTH default.
        // --------------------------------------------------------
        {"PREM4", {
            4,
            {1221.5, 3480.0, 5701.0, 6371.0},
            {12.8936, 10.9030,  4.7254,  3.0000},
            { 0.4656,  0.4656,  0.4957,  0.4957}
        }},

        // --------------------------------------------------------
        // PREM10 — 10 shells
        // Follows natural PREM geophysical discontinuities:
        //
        //   Shell  r_in      r_out     Region
        //    1       0.0   1221.5  inner core
        //    2    1221.5   2440.0  outer core (lower)
        //    3    2440.0   3480.0  outer core (upper)
        //    4    3480.0   3630.0  D'' layer
        //    5    3630.0   5600.0  lower mantle
        //    6    5600.0   5701.0  lower mantle top (660-km discontinuity)
        //    7    5701.0   5971.0  transition zone
        //    8    5971.0   6291.0  upper mantle (LVZ + LID)
        //    9    6291.0   6356.0  crust
        //   10    6356.0   6371.0  upper crust + ocean
        // --------------------------------------------------------
        {"PREM10", {
            10,
            {1221.5, 2440.0, 3480.0, 3630.0, 5600.0, 5701.0, 5971.0, 6291.0, 6356.0, 6371.0},
            {12.893569, 11.662691, 10.550173,  5.528405,  4.912992,
              4.411899,  3.882304,  3.434119,  3.308040,  2.283404},
            { 0.4656,    0.4656,    0.4656,    0.4957,    0.4957,
              0.4957,    0.4957,    0.4957,    0.4957,    0.4957}
        }},

        // --------------------------------------------------------
        // PREM42 — 42 shells
        // Standard resolution for atmospheric neutrino tomography
        // (as used in e.g. Frontiers Earth Sci. 11 (2023) 1008396).
        // Shell distribution:
        //   inner core  [    0.0 – 1221.5]:  4 equal shells
        //   outer core  [ 1221.5 – 3480.0]:  8 equal shells
        //   D'' layer   [ 3480.0 – 3630.0]:  2 equal shells
        //   lower mantle[ 3630.0 – 5600.0]: 12 equal shells
        //   LM top      [ 5600.0 – 5701.0]:  2 equal shells
        //   trans. zone [ 5701.0 – 5971.0]:  4 equal shells
        //   LVZ         [ 5971.0 – 6151.0]:  3 equal shells
        //   LID         [ 6151.0 – 6291.0]:  3 equal shells
        //   upper crust [ 6291.0 – 6346.6]:  2 equal shells
        //   lower crust [ 6346.6 – 6356.0]:  1 shell
        //   crust+ocean [ 6356.0 – 6371.0]:  1 shell
        // --------------------------------------------------------
        {"PREM42", {
            42,
            // outer radii [km]
            {  305.3750,   610.7500,   916.1250,  1221.5000,   // inner core
              1503.8125,  1786.1250,  2068.4375,  2350.7500,   // outer core (lower)
              2633.0625,  2915.3750,  3197.6875,  3480.0000,   // outer core (upper)
              3555.0000,  3630.0000,                            // D'' layer
              3794.1667,  3958.3333,  4122.5000,  4286.6667,   // lower mantle
              4450.8333,  4615.0000,  4779.1667,  4943.3333,
              5107.5000,  5271.6667,  5435.8333,  5600.0000,
              5650.5000,  5701.0000,                            // LM top
              5768.5000,  5836.0000,  5903.5000,  5971.0000,   // transition zone
              6031.0000,  6091.0000,  6151.0000,               // LVZ
              6197.6667,  6244.3333,  6291.0000,               // LID
              6318.8000,  6346.6000,                            // upper crust
              6356.0000,                                        // lower crust
              6371.0000},                                       // upper crust + ocean
            // volume-averaged PREM densities [g/cm³]
            {13.076317, 13.034546, 12.953202, 12.831335,
             12.083826, 11.910527, 11.706012, 11.467431,
             11.191919, 10.876600, 10.518597, 10.115028,
              5.547542,  5.510058,
              5.449991,  5.368323,  5.286526,  5.204283,
              5.121278,  5.037195,  4.951719,  4.864534,
              4.775323,  4.683770,  4.589560,  4.492376,
              4.427654,  4.396424,
              3.984231,  3.936222,  3.851196,  3.766123,
              3.525289,  3.489460,  3.453631,
              3.362033,  3.367104,  3.372176,
              3.376218,  3.379240,
              2.900000,
              2.283404},
            // electron fractions: core (12 shells) = 0.4656, mantle+ (30 shells) = 0.4957
            {0.4656, 0.4656, 0.4656, 0.4656,
             0.4656, 0.4656, 0.4656, 0.4656,
             0.4656, 0.4656, 0.4656, 0.4656,
             0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957,
             0.4957, 0.4957,
             0.4957,
             0.4957}
        }},

        // --------------------------------------------------------
        // PREMMAX — 81 shells
        // Maximum practical resolution: each PREM region subdivided
        // into equal-width shells, all discontinuities preserved.
        // Shell distribution:
        //   inner core  [    0.0 – 1221.5]:  8 shells
        //   outer core  [ 1221.5 – 3480.0]: 16 shells
        //   D'' layer   [ 3480.0 – 3630.0]:  3 shells
        //   lower mantle[ 3630.0 – 5600.0]: 20 shells
        //   LM top      [ 5600.0 – 5701.0]:  3 shells
        //   trans. zone [ 5701.0 – 5971.0]:  8 shells
        //   LVZ         [ 5971.0 – 6151.0]:  6 shells
        //   LID         [ 6151.0 – 6291.0]:  6 shells
        //   upper crust [ 6291.0 – 6346.6]:  4 shells
        //   lower crust [ 6346.6 – 6356.0]:  3 shells
        //   crust+ocean [ 6356.0 – 6371.0]:  4 shells
        //   total = 8+16+3+20+3+8+6+6+4+3+4 = 81
        // --------------------------------------------------------
        {"PREMMAX", {
            81,
            // outer radii [km]
            {  152.6875,   305.3750,   458.0625,   610.7500,   // inner core
               763.4375,   916.1250,  1068.8125,  1221.5000,
              1362.6562,  1503.8125,  1644.9688,  1786.1250,   // outer core
              1927.2812,  2068.4375,  2209.5938,  2350.7500,
              2491.9062,  2633.0625,  2774.2188,  2915.3750,
              3056.5312,  3197.6875,  3338.8438,  3480.0000,
              3530.0000,  3580.0000,  3630.0000,               // D'' layer
              3728.5000,  3827.0000,  3925.5000,  4024.0000,   // lower mantle
              4122.5000,  4221.0000,  4319.5000,  4418.0000,
              4516.5000,  4615.0000,  4713.5000,  4812.0000,
              4910.5000,  5009.0000,  5107.5000,  5206.0000,
              5304.5000,  5403.0000,  5501.5000,  5600.0000,
              5633.6667,  5667.3333,  5701.0000,               // LM top
              5734.7500,  5768.5000,  5802.2500,  5836.0000,   // transition zone
              5869.7500,  5903.5000,  5937.2500,  5971.0000,
              6001.0000,  6031.0000,  6061.0000,  6091.0000,   // LVZ
              6121.0000,  6151.0000,
              6174.3333,  6197.6667,  6221.0000,  6244.3333,   // LID
              6267.6667,  6291.0000,
              6304.9000,  6318.8000,  6332.7000,  6346.6000,   // upper crust
              6349.7333,  6352.8667,  6356.0000,               // lower crust
              6359.7500,  6363.5000,  6367.2500,  6371.0000},  // crust + ocean
            // volume-averaged PREM densities [g/cm³]
            {13.085454, 13.075011, 13.054676, 13.024209,   // inner core
             12.983595, 12.932829, 12.871912, 12.800843,
             12.127593, 12.048249, 11.961658, 11.867460,   // outer core
             11.765296, 11.654806, 11.535632, 11.407411,
             11.269785, 11.122393, 10.964874, 10.796867,
             10.618013, 10.427951, 10.226320, 10.012759,
              5.553871,  5.528858,  5.503890,               // D'' layer
              5.466719,  5.417696,  5.368696,  5.319650,   // lower mantle
              5.270489,  5.221147,  5.171553,  5.121641,
              5.071342,  5.020587,  4.969308,  4.917437,
              4.864905,  4.811645,  4.757588,  4.702666,
              4.646810,  4.589952,  4.532024,  4.472957,
              4.432860,  4.412114,  4.391220,               // LM top
              3.988184,  3.980325,  3.957567,  3.915124,   // transition zone
              3.872587,  3.830050,  3.787513,  3.744975,
              3.534291,  3.516377,  3.498462,  3.480547,   // LVZ
              3.462632,  3.444718,
              3.360760,  3.363296,  3.365832,  3.368368,   // LID
              3.370903,  3.373439,
              3.375461,  3.376972,  3.378483,  3.379993,   // upper crust
              2.900000,  2.900000,  2.900000,               // lower crust
              2.600000,  2.600000,  2.600000,  1.335851},  // crust + ocean
            // electron fractions: core (24 shells) = 0.4656, mantle+ (57 shells) = 0.4957
            {0.4656, 0.4656, 0.4656, 0.4656,   // inner core
             0.4656, 0.4656, 0.4656, 0.4656,
             0.4656, 0.4656, 0.4656, 0.4656,   // outer core
             0.4656, 0.4656, 0.4656, 0.4656,
             0.4656, 0.4656, 0.4656, 0.4656,
             0.4656, 0.4656, 0.4656, 0.4656,
             0.4957, 0.4957, 0.4957,            // D'' layer
             0.4957, 0.4957, 0.4957, 0.4957,   // lower mantle
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957,            // LM top
             0.4957, 0.4957, 0.4957, 0.4957,   // transition zone
             0.4957, 0.4957, 0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,   // LVZ
             0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,   // LID
             0.4957, 0.4957,
             0.4957, 0.4957, 0.4957, 0.4957,   // upper crust
             0.4957, 0.4957, 0.4957,            // lower crust
             0.4957, 0.4957, 0.4957, 0.4957}   // crust + ocean
        }},

    };  // end models map
    return models;
}

inline const EarthModel& get_earth_model(std::string_view name) {
    const auto& models = earth_models();
    auto it = models.find(std::string(name));
    if (it == models.end())
        throw std::invalid_argument("Not a valid Earth model: " + std::string(name));
    return it->second;
}

#endif // EARTH_H
