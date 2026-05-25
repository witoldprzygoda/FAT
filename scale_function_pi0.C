// scale_function_pi0.C — finalna gładka funkcja skalująca R(OA) z fit_pi0.C (pi0 Dalitz).
//
// R(x) = 1 + A * t^alpha * (1-t)^beta,   t = (x - xL)/(xR - xL)
// Parametry z fitu: A = 7.50596, alpha = 1.45695, beta = 2.37547
// xL = 0.50 deg, xR = 6.00 deg
//
// Użycie:
//   root -l
//   .L scale_function_pi0.C
//   scale_function_pi0(2.5);

#include <cmath>

double scale_function_pi0(double x)
{
    const double xL = 0.50;
    const double xR = 6.00;

    const double A     = 7.50596;
    const double alpha = 1.45695;
    const double beta  = 2.37547;

    if (x <= xL) return 1.0;
    if (x >= xR) return 1.0;

    const double t = (x - xL) / (xR - xL);

    return 1.0 + A * std::pow(t, alpha) * std::pow(1.0 - t, beta);
}
