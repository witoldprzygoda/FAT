// scale_function.C — finalna gładka funkcja skalująca R(OA) z fit.C (eta Dalitz).
//
// R(x) = 1 + N * t^a * (1-t)^b,   t = (x - xL)/(xR - xL)
// Parametry z fitu: N = 10.3171, a = 2.12720, b = 2.41167
// xL = 0.25 deg, xR = 5.00 deg
//
// Użycie:
//   root -l
//   .L scale_function.C
//   scale_function(2.5);   // -> ~1.45

#include <cmath>

double scale_function(double x)
{
    const double xL = 0.25;
    const double xR = 5.00;

    const double N = 10.3171;
    const double a = 2.12720;
    const double b = 2.41167;

    if (x <= xL) return 1.0;
    if (x >= xR) return 1.0;

    const double t = (x - xL) / (xR - xL);

    return 1.0 + N * std::pow(t, a) * std::pow(1.0 - t, b);
}
