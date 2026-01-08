#ifndef OPT_CONSTANTS_H
#define OPT_CONSTANTS_H

#include <cmath>

// Common constants used by both CHIC implementations
namespace OptConstants {
static constexpr double PI = std::acos(-1);
static constexpr double WEAK = 7.632470714045e-5;
static constexpr double BASELINE_FACTOR = 5.067730716156394;
static constexpr double INV_3 = 1.0 / 3.0;
static constexpr double INV_9 = 1.0 / 9.0;
static constexpr double TWO_THIRDS = 2.0 / 3.0;
static constexpr double INV_27 = 1.0 / 27.0;
static constexpr double SQRT_54 = std::sqrt(54.0);
static constexpr double INV_6 = 1.0 / 6.0;
static constexpr double TWO_OVER_3 = 2.0 / 3.0;
} // namespace OptConstants

#endif // OPT_CONSTANTS_H
