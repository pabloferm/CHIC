#pragma once
#include <vector>
#include <string>

#ifndef EARTH_H
#define EARTH_H

static constexpr double R_EARTH = 6371.0;  // km

struct EarthModel {
    int Nlayers;
    std::vector<double> radii;
    std::vector<double> density;
    std::vector<double> Ye;
};

EarthModel make_PREM4() {
    return {
        4,
        {1221.5, 3480.0, 5701.0, 6371.0},   // radii
        {13.0,   11.0,   5.0,    3.0},      // density
        {0.467,  0.467,  0.494,  0.494}     // Ye
    };
}



// template <size_t _Nlayers>
// struct EarthModel {
//     const int Nlayers = int(Nlayers);
//     double density[_Nlayers];
//     double Ye[_Nlayers];
//     double radii[_Nlayers];

//     // Constructor to initialize the arrays with values
//     EarthModel(const double (&d)[_Nlayers], const double (&y)[_Nlayers], const double (&r)[_Nlayers])
//         : density(), Ye(), radii() {
//         for (int i = 0; i < Nlayers; ++i) {
//             density[i] = d[i];
//             Ye[i] = y[i];
//             radii[i] = r[i];
//         }
//     }
// };

// //=================//
// //= PREM 4 layers =//
// //=================//
// constexpr double density4[4] = {2.8, 3.3, 10.0, 12.8};  // density in g/cm^3
// constexpr double Ye4[4] = {0.5, 0.5, 0.5, 0.5};        // effective electron density
// constexpr double radii4[4] = {1.0, 0.9945, 0.5384, 0.1931};  // outer radii fraction
// EarthModel<4> PREM4(density4, Ye4, radii4);

// #endif // EARTH_H
