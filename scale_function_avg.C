// scale_function_avg.C — średnia i błąd z dwóch funkcji skalujących (pi0 + eta Dalitz).
//
// avg(x)  = (R_pi0(x) + R_eta(x)) / 2
// err(x)  = (max - min) / 2 = |R_pi0(x) - R_eta(x)| / 2
//
// Użycie:
//   root -l
//   .L scale_function_avg.C
//   double v, e;
//   scale_function_avg(2.5, v, e);
//   // albo:
//   scale_function_avg(2.5);   // wypisze v i e na ekran

#include <cmath>
#include <iostream>

// --- pi0 Dalitz: xL=0.50, xR=6.00,  A=7.506, alpha=1.457, beta=2.375 ---
double scale_function_pi0(double x)
{
    const double xL = 0.50, xR = 6.00;
    const double A = 7.50596, alpha = 1.45695, beta = 2.37547;
    if (x <= xL || x >= xR) return 1.0;
    const double t = (x - xL) / (xR - xL);
    return 1.0 + A * std::pow(t, alpha) * std::pow(1.0 - t, beta);
}

// --- eta Dalitz: xL=0.25, xR=5.00,  N=10.317, a=2.127, b=2.412 ---
double scale_function_eta(double x)
{
    const double xL = 0.25, xR = 5.00;
    const double N = 10.3171, a = 2.12720, b = 2.41167;
    if (x <= xL || x >= xR) return 1.0;
    const double t = (x - xL) / (xR - xL);
    return 1.0 + N * std::pow(t, a) * std::pow(1.0 - t, b);
}

// Wersja zwracająca przez referencję.
void scale_function_avg(double x, double &val, double &err)
{
    const double r_pi0 = scale_function_pi0(x);
    const double r_eta = scale_function_eta(x);

    val = 0.5 * (r_pi0 + r_eta);
    err = 0.5 * std::fabs(r_pi0 - r_eta);
}

// Wersja interaktywna — drukuje wynik.
void scale_function_avg(double x)
{
    double v, e;
    scale_function_avg(x, v, e);
    const double r_pi0 = scale_function_pi0(x);
    const double r_eta = scale_function_eta(x);

    std::cout << "x = " << x << " deg\n"
              << "  R_pi0 = " << r_pi0 << "\n"
              << "  R_eta = " << r_eta << "\n"
              << "  avg   = " << v << "\n"
              << "  err   = " << e << "\n"
              << "  R     = " << v << " +/- " << e << "\n";
}
