#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "CHIC.h"

int main()
{
    using namespace std::chrono;

    // ------------------------------------- //
    // Set the vacuum oscillation parameters //
    // ------------------------------------- //
    CHIC chic;
    chic.update_dcp(-2.2);
    chic.update_dm231(2.5e-3);
    chic.update_dm221(7.5e-5);
    chic.update_th12(0.590500015160318);   // sin^2(th12)=0.31
    chic.update_th13(0.1418970546041639);   // sin^2(th13)=0.02
    chic.update_th23(0.8354818739782283);   // sin^2(th23)=0.55
    chic.update_density(3.0);     // g/cm^3

    // ------------------------------- //
    // Set the sweep parameters        //
    // ------------------------------- //
    const int    numPoints = 500;
    const double Emin      = 0.05;    // GeV
    const double Emax      = 50.0;    // GeV
    const double Lmin      = 0.1;     // km
    const double Lmax      = 10000.;  // km
    const double L_fixed   = 295.;    // km  (used in energy-only loop)
    const double E_fixed   = 2.5;     // GeV (used in baseline-only loop)
    const int    warmup    = 200;     // calls to skip before timing

    long long total_duration_ns;
    long long count;

    // ================================================ //
    //  1. Double loop: E x L                           //
    // ================================================ //
    total_duration_ns = 0;
    count             = 0;

    std::cout << "# Double loop (E x L)\n";
    std::cout << "E(GeV),L(km),"
              << "P(nue->nue),P(numu->nue),P(nutau->nue),"
              << "P(nue->numu),P(numu->numu),P(nutau->numu),"
              << "P(nue->nutau),P(numu->nutau),P(nutau->nutau)\n";

    for (int i = 0; i < numPoints; i++)
    {
        double t = static_cast<double>(i) / (numPoints - 1);
        double E = Emin * std::pow(Emax / Emin, t);

        for (int j = 0; j < numPoints; j++)
        {
            double s = static_cast<double>(j) / (numPoints - 1);
            double L = Lmin * std::pow(Lmax / Lmin, s);

            auto start = high_resolution_clock::now();
            Eigen::Matrix3d prob = chic.compute_oscillations(E, L);
            auto end   = high_resolution_clock::now();

            if (count >= warmup)
                total_duration_ns += duration_cast<nanoseconds>(end - start).count();
            ++count;

            std::cout << std::scientific << std::setprecision(15)
                      << E << "," << L << ","
                      << prob(0,0) << "," << prob(1,0) << "," << prob(2,0) << ","
                      << prob(0,1) << "," << prob(1,1) << "," << prob(2,1) << ","
                      << prob(0,2) << "," << prob(1,2) << "," << prob(2,2) << "\n";
        }
    }
    {
        long long effective = count - warmup;
        std::cout << "# Double loop  | calls: " << count
                  << "  | avg: " << static_cast<double>(total_duration_ns) / effective
                  << " ns\n\n";
    }

    // ================================================ //
    //  2. Energy loop only (fixed L = L_fixed)         //
    // ================================================ //
    total_duration_ns = 0;
    count             = 0;

    std::cout << "# Energy loop only (L = " << L_fixed << " km)\n";
    std::cout << "E(GeV),L(km),"
              << "P(nue->nue),P(numu->nue),P(nutau->nue),"
              << "P(nue->numu),P(numu->numu),P(nutau->numu),"
              << "P(nue->nutau),P(numu->nutau),P(nutau->nutau)\n";

    for (int i = 0; i < numPoints; i++)
    {
        double t = static_cast<double>(i) / (numPoints - 1);
        double E = Emin * std::pow(Emax / Emin, t);

        auto start = high_resolution_clock::now();
        Eigen::Matrix3d prob = chic.compute_oscillations(E, L_fixed);
        auto end   = high_resolution_clock::now();

        if (count >= warmup)
            total_duration_ns += duration_cast<nanoseconds>(end - start).count();
        ++count;

        std::cout << std::scientific << std::setprecision(15)
                  << E << "," << L_fixed << ","
                  << prob(0,0) << "," << prob(1,0) << "," << prob(2,0) << ","
                  << prob(0,1) << "," << prob(1,1) << "," << prob(2,1) << ","
                  << prob(0,2) << "," << prob(1,2) << "," << prob(2,2) << "\n";
    }
    {
        long long effective = count - warmup;
        std::cout << "# Energy loop  | calls: " << count
                  << "  | avg: " << static_cast<double>(total_duration_ns) / effective
                  << " ns\n\n";
    }

    // ================================================ //
    //  3. Baseline loop only (fixed E = E_fixed)       //
    // ================================================ //
    total_duration_ns = 0;
    count             = 0;

    std::cout << "# Baseline loop only (E = " << E_fixed << " GeV)\n";
    std::cout << "E(GeV),L(km),"
              << "P(nue->nue),P(numu->nue),P(nutau->nue),"
              << "P(nue->numu),P(numu->numu),P(nutau->numu),"
              << "P(nue->nutau),P(numu->nutau),P(nutau->nutau)\n";

    for (int j = 0; j < numPoints; j++)
    {
        double s = static_cast<double>(j) / (numPoints - 1);
        double L = Lmin * std::pow(Lmax / Lmin, s);

        auto start = high_resolution_clock::now();
        Eigen::Matrix3d prob = chic.compute_oscillations(E_fixed, L);
        auto end   = high_resolution_clock::now();

        if (count >= warmup)
            total_duration_ns += duration_cast<nanoseconds>(end - start).count();
        ++count;

        std::cout << std::scientific << std::setprecision(15)
                  << E_fixed << "," << L << ","
                  << prob(0,0) << "," << prob(1,0) << "," << prob(2,0) << ","
                  << prob(0,1) << "," << prob(1,1) << "," << prob(2,1) << ","
                  << prob(0,2) << "," << prob(1,2) << "," << prob(2,2) << "\n";
    }
    {
        long long effective = count - warmup;
        std::cout << "# Baseline loop | calls: " << count
                  << "  | avg: " << static_cast<double>(total_duration_ns) / effective
                  << " ns\n\n";
    }

    return 0;
}
