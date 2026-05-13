#include "CHIC_EARTH.h"

#include <Eigen/Dense>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#define sq(x) ((x)*(x))

std::string optimization_level = "3";

using clock_type = std::chrono::high_resolution_clock;

// ============================================================
// Atmospheric benchmark
// ============================================================

void Atmospheric_Speed(
    CHICEARTH& earth,
    bool quick_osc_param,
    int n_1,
    int n_2,
    FILE* data)
{
    const double E  = 3.0;
    const double cz = -1.0;

    Eigen::Matrix3d probs;

    double speed, speed_sum, speedsq_sum;

    speed_sum   = 0;
    speedsq_sum = 0;

    // warmup
    probs = earth.compute_oscillations(E, cz);

    for (int i = 0; i < n_1; i++)
    {
        auto start_time = clock_type::now();

        for (int j = 0; j < n_2; j++)
        {
            if (quick_osc_param)
            {
                double theta23 =
                        0.72 + i * 0.1 / n_2;

                earth.update_th23(theta23);
            }
            else
            {
                double dm31 =
                    2.5e-3 + i * 0.1e-3 / n_2;

                earth.update_dm231(dm31);
            }

            probs = earth.compute_oscillations(E, cz);
        }

        auto end_time = clock_type::now();

        speed =
            std::chrono::duration_cast<
                std::chrono::duration<double>
            >(end_time - start_time).count()
            / n_2 * 1e9;

        speed_sum   += speed;
        speedsq_sum += sq(speed);
    }

    if (quick_osc_param)
    {
        fprintf(data, "quick ");
        printf("quick ");
    }
    else
    {
        fprintf(data, "full ");
        printf("full ");
    }

    printf(
        "time = %g +- %g ns\n",
        speed_sum / n_1,
        sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1))
    );

    fprintf(
        data,
        "%g %g\n",
        speed_sum / n_1,
        sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1))
    );
}

void Atmospheric_Speed()
{
    std::string fname =
        "data/speed/Atmospherics_"
        + optimization_level
        + ".txt";

    FILE* data = fopen(fname.c_str(), "w");

    if (!data)
    {
        perror("fopen failed");
        return;
    }

    CHICEARTH earth("neutrino");

    for (int i = 0; i < 2; i++)
    {
        Atmospheric_Speed(
            earth,
            i == 0,
            1e4,
            1e2,
            data
        );
    }

    fclose(data);
}

// ============================================================
// E vs cosz benchmark
// ============================================================

void E_vs_cosz_Speed(
    CHICEARTH& earth,
    int nE,
    int ncosz,
    FILE* data)
{
    std::vector<double> Es;
    std::vector<double> coszs;

    double Emin     = 2;
    double Emax     = 40;

    double coszmin  = -1;
    double coszmax  = -0.9999999;

    double Estep    = (Emax - Emin) / nE;
    double coszstep = (coszmax - coszmin) / ncosz;

    Es.reserve(nE);
    coszs.reserve(ncosz);

    for (int i = 0; i < nE; i++)
        Es.emplace_back(Emin + i * Estep);

    for (int i = 0; i < ncosz; i++)
        coszs.emplace_back(coszmin + i * coszstep);

    // warmup
    earth.compute_oscillations(Es[0], coszs[0]);

    int n_1 = 5e1;
    int n_2 = 5e1;

    double speed;
    double speed_sum;
    double speedsq_sum;
    double speed_mean;
    double speed_std;

    speed_sum   = 0;
    speedsq_sum = 0;

    Eigen::Matrix3d probs;

    for (int i = 0; i < n_1; i++)
    {
        auto start_time = clock_type::now();

        for (int j = 0; j < n_2; j++)
        {
            earth.update_dm231(
                2.5e-3 + i * 0.1e-3 / n_2
            );

            for (int ie = 0; ie < nE; ie++)
            {
                for (int icz = 0; icz < ncosz; icz++)
                {
                    probs =
                        earth.compute_oscillations(
                            Es[ie],
                            coszs[icz]
                        );
                }
            }
        }

        auto end_time = clock_type::now();

        speed =
            std::chrono::duration_cast<
                std::chrono::duration<double>
            >(end_time - start_time).count()
            / n_2 * 1e9;

        speed_sum   += speed;
        speedsq_sum += sq(speed);
    }

    speed_mean =
        speed_sum / n_1 / nE / ncosz;

    speed_std =
        sqrt(
            speedsq_sum / n_1
            - sq(speed_sum / n_1)
        ) / nE / ncosz;

    printf(
        "nE = %d, ncosz = %d, time = %g +- %g ns\n",
        nE,
        ncosz,
        speed_mean,
        speed_std
    );

    fprintf(
        data,
        "%d %d %g %g\n",
        nE,
        ncosz,
        speed_mean,
        speed_std
    );
}

void E_vs_cosz_Speed()
{
    std::string fname =
        "data/speed/E_vs_cosz_"
        + optimization_level
        + ".txt";

    FILE* data = fopen(fname.c_str(), "w");

    if (!data)
    {
        perror("fopen failed");
        return;
    }

    CHICEARTH earth("neutrino");

    for (int j = 1; j <= 3; j++)
    {
        E_vs_cosz_Speed(
            earth,
            1,
            pow(10, j),
            data
        );

        E_vs_cosz_Speed(
            earth,
            pow(10, j),
            1,
            data
        );
    }

    fclose(data);
}

// ============================================================
// main
// ============================================================

int main()
{
    system("mkdir -p data/speed");

    Atmospheric_Speed();

    E_vs_cosz_Speed();

    return 0;
}
