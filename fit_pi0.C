// ============================================================================
// fit_pi0.C
//
// Wersja wygładzona:
//
// 1. Yield EXP/SIM:
//      Y(x) = exp[Chebyshev(log x)]
//    ale z mniejszą liczbą parametrów:
//      YIELD_DEG = 10
//
// 2. Yield fit/drawing zaczyna się od OA = 0.50,
//    żeby uniknąć artefaktu brzegowego przy bardzo niskich OA.
//
// 3. Common anchors:
//      - pierwsze trzy punkty: OA = 0.10, 0.30, 0.50
//      - wszystkie punkty OA >= 6.0
//
// 4. Ratio:
//      prosty dodatni bump:
//
//      R(x) = 1,                                  x <= 0.50
//      R(x) = 1 + A*t^alpha*(1-t)^beta,           0.50 < x < 6.00
//      R(x) = 1,                                  x >= 6.00
//
//      t = (x - 0.50)/(6.00 - 0.50)
//
//    To jest dużo mniej elastyczne niż poprzedni Chebyshev dla ratio,
//    więc nie powinno robić sztucznie połamanej krzywej.
//
// Uruchomienie:
//
//   root -l fit_pi0.C
//
// albo:
//
//   root -l
//   .L fit_pi0.C
//   fit_pi0()
//
// ============================================================================

#include <iostream>
#include <cmath>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

// ============================================================================
// Configuration
// ============================================================================

static const int NPTS = 75;

static const int YIELD_DEG  = 10;
static const int YIELD_NPAR = YIELD_DEG + 1;

static const double X_MIN_DATA  = 0.10;
static const double X_MIN_FIT   = 0.50;
static const double X_MAX       = 14.90;

static const double X_RATIO_LEFT = 0.50;
static const double X_RATIO_FLAT = 6.00;

// ============================================================================
// Input data
// ============================================================================

double OA[NPTS] = {
     0.10,  0.30,  0.50,  0.70,  0.90,
     1.10,  1.30,  1.50,  1.70,  1.90,
     2.10,  2.30,  2.50,  2.70,  2.90,
     3.10,  3.30,  3.50,  3.70,  3.90,
     4.10,  4.30,  4.50,  4.70,  4.90,
     5.10,  5.30,  5.50,  5.70,  5.90,
     6.10,  6.30,  6.50,  6.70,  6.90,
     7.10,  7.30,  7.50,  7.70,  7.90,
     8.10,  8.30,  8.50,  8.70,  8.90,
     9.10,  9.30,  9.50,  9.70,  9.90,
    10.10, 10.30, 10.50, 10.70, 10.90,
    11.10, 11.30, 11.50, 11.70, 11.90,
    12.10, 12.30, 12.50, 12.70, 12.90,
    13.10, 13.30, 13.50, 13.70, 13.90,
    14.10, 14.30, 14.50, 14.70, 14.90
};

double EX[NPTS] = {
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10,
    0.10, 0.10, 0.10, 0.10, 0.10
};

double Y_EXP[NPTS] = {
      4442.858,  14721.985,  27340.256,  43082.391,  61373.211,
     76537.125,  90357.109,  99895.836, 116665.508, 117769.375,
    116828.281, 112120.281, 102283.625,  95549.562,  86950.664,
     78048.031,  68912.695,  61566.848,  54915.379,  47510.820,
     42399.441,  37664.594,  34015.105,  30561.301,  27725.137,
     25139.730,  23332.826,  21286.613,  19872.275,  18262.000,
     17453.473,  16333.997,  15769.065,  14949.738,  14204.017,
     13519.814,  13023.897,  12483.342,  11894.700,  11438.716,
     10165.751,   9708.549,  10060.918,   9168.984,   8739.830,
      8444.884,   8502.866,   7950.846,   7593.554,   7554.796,
      7398.008,   7118.540,   7127.146,   6851.275,   6697.279,
      6728.697,   6304.032,   5838.902,   6004.288,   5920.101,
      5637.324,   5636.839,   5480.261,   5216.864,   5288.264,
      5010.884,   4883.057,   4690.852,   4696.927,   4598.343,
      4648.832,   4417.613,   4455.843,   3876.733,   3959.274
};

double EY_EXP[NPTS] = {
     77.434, 140.730, 188.780, 237.013, 281.764,
    315.935, 343.210, 362.192, 389.778, 393.685,
    392.219, 385.031, 373.723, 360.993, 345.375,
    327.758, 309.086, 291.506, 275.303, 257.194,
    242.911, 229.111, 216.878, 205.587, 195.809,
    186.896, 179.767, 171.724, 166.000, 159.484,
    155.322, 150.250, 147.099, 143.360, 140.132,
    136.415, 133.993, 131.042, 127.973, 125.678,
    118.727, 116.233, 118.794, 112.801, 109.900,
    108.697, 108.333, 105.991, 104.082, 103.359,
    102.162, 101.143, 100.025,  98.874,  98.229,
     96.974,  95.828,  89.961,  93.461,  92.644,
     91.198,  91.209,  89.883,  88.448,  89.079,
     86.319,  86.209,  85.691,  84.273,  84.682,
     83.923,  82.849,  83.162,  76.870,  80.218
};

double Y_SIM[NPTS] = {
      2915.547,  10815.353,  23035.730,  37878.557,  52193.160,
     64313.478,  73029.315,  75185.509,  76745.422,  75852.717,
     73125.684,  69366.526,  64936.993,  60127.063,  55369.449,
     50811.749,  46675.071,  42725.114,  39153.470,  35809.016,
     32954.444,  29197.349,  27080.080,  25141.965,  23547.194,
     21936.649,  20748.626,  19377.831,  18355.819,  17431.173,
     16460.739,  15620.450,  14869.760,  14271.526,  13592.401,
     13013.727,  12500.667,  12030.392,  11581.697,  11179.031,
     10809.220,  10388.612,  10098.497,   9771.726,   9435.577,
      9122.468,   8863.351,   8565.816,   7924.511,   8051.886,
      7440.551,   7218.630,   7086.588,   6872.613,   6961.963,
      6501.651,   6347.500,   6393.941,   6258.123,   5841.801,
      5711.923,   5536.440,   5108.441,   4947.845,   4928.737,
      4978.642,   4638.572,   4509.764,   4385.223,   4312.650,
      4195.839,   4114.923,   4011.447,   3891.846,   3809.436
};

double EY_SIM[NPTS] = {
    15.864, 30.538, 44.495, 57.053, 67.026,
    74.477, 79.444, 80.621, 81.581, 81.257,
    79.921, 77.977, 75.566, 72.847, 69.995,
    67.129, 64.154, 61.415, 58.850, 56.376,
    54.138, 50.914, 49.057, 47.312, 45.746,
    44.298, 42.971, 41.657, 40.478, 39.440,
    38.346, 37.451, 36.555, 35.860, 35.002,
    34.290, 33.630, 32.982, 32.379, 31.811,
    31.235, 30.690, 30.244, 29.725, 29.217,
    28.777, 28.365, 27.895, 26.772, 27.049,
    25.957, 25.572, 25.305, 24.925, 25.161,
    24.230, 23.936, 24.127, 23.837, 22.949,
    22.693, 22.363, 21.449, 21.123, 21.016,
    21.178, 20.405, 20.119, 19.855, 19.663,
    19.394, 19.200, 18.955, 18.654, 18.454
};

// ============================================================================
// Parametrization arrays
// ============================================================================

double Y_EXP_PAR[NPTS];
double EY_EXP_PAR[NPTS];

double Y_SIM_PAR[NPTS];
double EY_SIM_PAR[NPTS];

// ============================================================================
// Utility functions
// ============================================================================

void WeightedAverage2(
    double y1,
    double e1,
    double y2,
    double e2,
    double &yavg,
    double &eyavg
)
{
    const double w1 = 1.0 / (e1 * e1);
    const double w2 = 1.0 / (e2 * e2);

    yavg  = (y1 * w1 + y2 * w2) / (w1 + w2);
    eyavg = 1.0 / std::sqrt(w1 + w2);
}

double MapLogXToUnit(double x)
{
    const double lx    = std::log(x);
    const double lxMin = std::log(X_MIN_FIT);
    const double lxMax = std::log(X_MAX);

    return 2.0 * (lx - lxMin) / (lxMax - lxMin) - 1.0;
}

double ChebyshevEval(double u, const double *par, int npar)
{
    if (npar <= 0) return 0.0;
    if (npar == 1) return par[0];

    double t0 = 1.0;
    double t1 = u;

    double result = par[0] * t0 + par[1] * t1;

    for (int i = 2; i < npar; ++i)
    {
        const double t2 = 2.0 * u * t1 - t0;
        result += par[i] * t2;

        t0 = t1;
        t1 = t2;
    }

    return result;
}

// ============================================================================
// Build parametrization points
//
// Common weighted averages:
//   - first three points,
//   - OA >= 6.0.
//
// Independent region:
//   - 0.50 < OA < 6.0.
// ============================================================================

void BuildParametrizationArrays()
{
    std::cout << "\n============================================================\n";
    std::cout << "Preparing pi0 parametrization points\n";
    std::cout << "Common anchors: first three bins and OA >= 6 deg\n";
    std::cout << "Independent EXP/SIM region: 0.50 < OA < 6 deg\n";
    std::cout << "Yield model: exp(Chebyshev polynomial in log(OA)), degree "
              << YIELD_DEG << "\n";
    std::cout << "============================================================\n\n";

    for (int i = 0; i < NPTS; ++i)
    {
        const bool useCommonPoint = (i < 3) || (OA[i] >= 6.0);

        if (useCommonPoint)
        {
            double yavg  = 0.0;
            double eyavg = 0.0;

            WeightedAverage2(
                Y_EXP[i], EY_EXP[i],
                Y_SIM[i], EY_SIM[i],
                yavg, eyavg
            );

            Y_EXP_PAR[i]  = yavg;
            EY_EXP_PAR[i] = eyavg;

            Y_SIM_PAR[i]  = yavg;
            EY_SIM_PAR[i] = eyavg;
        }
        else
        {
            Y_EXP_PAR[i]  = Y_EXP[i];
            EY_EXP_PAR[i] = EY_EXP[i];

            Y_SIM_PAR[i]  = Y_SIM[i];
            EY_SIM_PAR[i] = EY_SIM[i];
        }
    }
}

// ============================================================================
// Yield model in log space
// ============================================================================

double LogYieldCheb(double *xx, double *par)
{
    const double x = xx[0];

    if (x <= 0.0) return 0.0;

    const double u = MapLogXToUnit(x);
    return ChebyshevEval(u, par, YIELD_NPAR);
}

// ============================================================================
// Yield model in normal space
// ============================================================================

double YieldCheb(double *xx, double *par)
{
    const double x = xx[0];

    // Nie rysujemy / nie ekstrapolujemy Chebysheva poniżej 0.50,
    // bo tam był artefakt brzegowy.
    if (x < X_MIN_FIT) return 0.0;

    const double logy = LogYieldCheb(xx, par);

    if (logy > 30.0) return std::exp(30.0);
    if (logy < -30.0) return std::exp(-30.0);

    return std::exp(logy);
}

// ============================================================================
// Smooth positive ratio bump
//
// R(x) = 1                                         for x <= xL
// R(x) = 1 + A * t^alpha * (1-t)^beta             for xL < x < xR
// R(x) = 1                                         for x >= xR
//
// t = (x-xL)/(xR-xL)
//
// par[0] = A
// par[1] = alpha
// par[2] = beta
// ============================================================================

double RatioBetaBump(double *xx, double *par)
{
    const double x = xx[0];

    const double A     = par[0];
    const double alpha = par[1];
    const double beta  = par[2];

    const double xL = X_RATIO_LEFT;
    const double xR = X_RATIO_FLAT;

    if (x <= xL) return 1.0;
    if (x >= xR) return 1.0;

    const double t = (x - xL) / (xR - xL);

    if (t <= 0.0) return 1.0;
    if (t >= 1.0) return 1.0;

    return 1.0 + A * std::pow(t, alpha) * std::pow(1.0 - t, beta);
}

// ============================================================================
// Main macro
// ============================================================================

void fit_pi0()
{
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    BuildParametrizationArrays();

    // ========================================================================
    // Original graphs
    // ========================================================================

    TGraphErrors *gExpOriginal =
        new TGraphErrors(NPTS, OA, Y_EXP, EX, EY_EXP);

    gExpOriginal->SetName("gExpOriginal");
    gExpOriginal->SetTitle(
        "#pi^{0} yield vs OA;OA(e^{+}e^{-}) [deg];yield(#mu #pm 1#sigma)"
    );

    gExpOriginal->SetMarkerStyle(20);
    gExpOriginal->SetMarkerSize(0.75);
    gExpOriginal->SetMarkerColor(kBlack);
    gExpOriginal->SetLineColor(kBlack);

    TGraphErrors *gSimOriginal =
        new TGraphErrors(NPTS, OA, Y_SIM, EX, EY_SIM);

    gSimOriginal->SetName("gSimOriginal");
    gSimOriginal->SetMarkerStyle(24);
    gSimOriginal->SetMarkerSize(0.75);
    gSimOriginal->SetMarkerColor(kRed + 1);
    gSimOriginal->SetLineColor(kRed + 1);

    // ========================================================================
    // Graphs used for yield parametrization
    // ========================================================================

    TGraphErrors *gExpPar =
        new TGraphErrors(NPTS, OA, Y_EXP_PAR, EX, EY_EXP_PAR);

    gExpPar->SetName("gExpPar");
    gExpPar->SetMarkerStyle(20);
    gExpPar->SetMarkerSize(0.45);
    gExpPar->SetMarkerColor(kBlue + 1);
    gExpPar->SetLineColor(kBlue + 1);

    TGraphErrors *gSimPar =
        new TGraphErrors(NPTS, OA, Y_SIM_PAR, EX, EY_SIM_PAR);

    gSimPar->SetName("gSimPar");
    gSimPar->SetMarkerStyle(24);
    gSimPar->SetMarkerSize(0.45);
    gSimPar->SetMarkerColor(kMagenta + 1);
    gSimPar->SetLineColor(kMagenta + 1);

    // ========================================================================
    // Build log-yield graphs for fitting
    // Only OA >= 0.50 enters Chebyshev fit.
    // ========================================================================

    double X_LOG[NPTS];
    double EX_LOG[NPTS];

    double LOG_EXP[NPTS];
    double ELOG_EXP[NPTS];

    double LOG_SIM[NPTS];
    double ELOG_SIM[NPTS];

    int nLog = 0;

    for (int i = 0; i < NPTS; ++i)
    {
        if (OA[i] >= X_MIN_FIT)
        {
            X_LOG[nLog]  = OA[i];
            EX_LOG[nLog] = 0.0;

            LOG_EXP[nLog] = std::log(Y_EXP_PAR[i]);
            LOG_SIM[nLog] = std::log(Y_SIM_PAR[i]);

            ELOG_EXP[nLog] = EY_EXP_PAR[i] / Y_EXP_PAR[i];
            ELOG_SIM[nLog] = EY_SIM_PAR[i] / Y_SIM_PAR[i];

            // Nieco większa waga głównego piku, ale łagodniej niż wcześniej.
            if (OA[i] < 6.0)
            {
                ELOG_EXP[nLog] *= 0.45;
                ELOG_SIM[nLog] *= 0.45;
            }

            ++nLog;
        }
    }

    TGraphErrors *gExpLog =
        new TGraphErrors(nLog, X_LOG, LOG_EXP, EX_LOG, ELOG_EXP);

    TGraphErrors *gSimLog =
        new TGraphErrors(nLog, X_LOG, LOG_SIM, EX_LOG, ELOG_SIM);

    // ========================================================================
    // Fit log-yield parametrizations
    // ========================================================================

    TF1 *fLogExp = new TF1(
        "fLogExp",
        LogYieldCheb,
        X_MIN_FIT,
        X_MAX,
        YIELD_NPAR
    );

    TF1 *fLogSim = new TF1(
        "fLogSim",
        LogYieldCheb,
        X_MIN_FIT,
        X_MAX,
        YIELD_NPAR
    );

    for (int ip = 0; ip < YIELD_NPAR; ++ip)
    {
        fLogExp->SetParameter(ip, 0.0);
        fLogSim->SetParameter(ip, 0.0);
    }

    fLogExp->SetParameter(0, std::log(30000.0));
    fLogSim->SetParameter(0, std::log(25000.0));

    std::cout << "\n\n========== FIT LOG-YIELD EXP ==========\n";
    gExpLog->Fit(fLogExp, "RME");

    std::cout << "\n\n========== FIT LOG-YIELD SIM ==========\n";
    gSimLog->Fit(fLogSim, "RME");

    TF1 *fYieldExp = new TF1(
        "fYieldExp",
        YieldCheb,
        X_MIN_FIT,
        X_MAX,
        YIELD_NPAR
    );

    TF1 *fYieldSim = new TF1(
        "fYieldSim",
        YieldCheb,
        X_MIN_FIT,
        X_MAX,
        YIELD_NPAR
    );

    for (int ip = 0; ip < YIELD_NPAR; ++ip)
    {
        fYieldExp->SetParameter(ip, fLogExp->GetParameter(ip));
        fYieldSim->SetParameter(ip, fLogSim->GetParameter(ip));
    }

    fYieldExp->SetLineColor(kBlue + 1);
    fYieldExp->SetLineWidth(3);

    fYieldSim->SetLineColor(kRed + 1);
    fYieldSim->SetLineWidth(3);

    // ========================================================================
    // Canvas 1: yields and analytical parametrizations
    // ========================================================================

    TCanvas *c1 =
        new TCanvas("c1", "pi0 EXP and SIM analytical parametrizations", 1200, 800);

    c1->SetGrid();

    gExpOriginal->Draw("AP");
    gExpOriginal->GetXaxis()->SetLimits(0.0, 15.2);
    gExpOriginal->GetYaxis()->SetRangeUser(0.0, 130000.0);

    gSimOriginal->Draw("P SAME");

    gExpPar->Draw("P SAME");
    gSimPar->Draw("P SAME");

    fYieldExp->Draw("SAME");
    fYieldSim->Draw("SAME");

    TLine *lineCommonLow = new TLine(X_RATIO_LEFT, 0.0, X_RATIO_LEFT, 130000.0);
    lineCommonLow->SetLineStyle(3);
    lineCommonLow->SetLineColor(kGray + 1);
    lineCommonLow->Draw("SAME");

    TLine *lineCommonHigh = new TLine(X_RATIO_FLAT, 0.0, X_RATIO_FLAT, 130000.0);
    lineCommonHigh->SetLineStyle(3);
    lineCommonHigh->SetLineColor(kGray + 1);
    lineCommonHigh->Draw("SAME");

    TLegend *leg1 = new TLegend(0.47, 0.57, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    leg1->AddEntry(gExpOriginal, "EXP original", "lep");
    leg1->AddEntry(gSimOriginal, "SIM original (#times 6.37149)", "lep");
    leg1->AddEntry(gExpPar, "EXP parametrization points", "p");
    leg1->AddEntry(gSimPar, "SIM parametrization points", "p");
    leg1->AddEntry(fYieldExp, "EXP exp(Chebyshev(log OA)), deg=10", "l");
    leg1->AddEntry(fYieldSim, "SIM exp(Chebyshev(log OA)), deg=10", "l");
    leg1->AddEntry(lineCommonHigh, "common anchors: OA #leq 0.5 and OA #geq 6", "l");

    leg1->Draw();

    c1->Update();
    c1->SaveAs("pi0_yield_analytic_logcheb_deg10.png");

    // ========================================================================
    // Original point ratio
    // ========================================================================

    double R_ORIG[NPTS];
    double ER_ORIG[NPTS];

    for (int i = 0; i < NPTS; ++i)
    {
        R_ORIG[i] = Y_EXP[i] / Y_SIM[i];

        const double relExp = EY_EXP[i] / Y_EXP[i];
        const double relSim = EY_SIM[i] / Y_SIM[i];

        ER_ORIG[i] = R_ORIG[i] * std::sqrt(relExp * relExp + relSim * relSim);
    }

    TGraphErrors *gRatioOriginal =
        new TGraphErrors(NPTS, OA, R_ORIG, EX, ER_ORIG);

    gRatioOriginal->SetName("gRatioOriginal");
    gRatioOriginal->SetTitle(
        "#pi^{0} ratio EXP/SIM;OA(e^{+}e^{-}) [deg];EXP / SIM"
    );

    gRatioOriginal->SetMarkerStyle(20);
    gRatioOriginal->SetMarkerSize(0.75);
    gRatioOriginal->SetMarkerColor(kBlack);
    gRatioOriginal->SetLineColor(kBlack);

    // ========================================================================
    // Ratio graph for direct fit
    //
    // We include:
    //   - first three points as artificial anchors at R=1,
    //   - original point ratios for 0.50 < OA < 6.00.
    // ========================================================================

    double X_RFIT[NPTS];
    double EX_RFIT[NPTS];

    double R_FIT[NPTS];
    double ER_FIT[NPTS];

    int nRFit = 0;

    for (int i = 0; i < NPTS; ++i)
    {
        if (OA[i] <= X_RATIO_FLAT)
        {
            X_RFIT[nRFit]  = OA[i];
            EX_RFIT[nRFit] = 0.0;

            if (OA[i] <= X_RATIO_LEFT)
            {
                R_FIT[nRFit]  = 1.0;
                ER_FIT[nRFit] = 0.015;
            }
            else
            {
                R_FIT[nRFit]  = R_ORIG[i];
                ER_FIT[nRFit] = ER_ORIG[i];

                // Łagodniejsze ważenie ratio, żeby nie gonić każdego punktu.
                ER_FIT[nRFit] *= 1.8;
            }

            ++nRFit;
        }
    }

    TGraphErrors *gRatioFit =
        new TGraphErrors(nRFit, X_RFIT, R_FIT, EX_RFIT, ER_FIT);

    gRatioFit->SetName("gRatioFit");
    gRatioFit->SetMarkerStyle(24);
    gRatioFit->SetMarkerSize(0.8);
    gRatioFit->SetMarkerColor(kRed + 1);
    gRatioFit->SetLineColor(kRed + 1);

    // ========================================================================
    // Fit simple smooth positive ratio
    // ========================================================================

    TF1 *fRatio = new TF1(
        "fRatio",
        RatioBetaBump,
        0.0,
        15.0,
        3
    );

    fRatio->SetParNames("A", "alpha", "beta");

    fRatio->SetParameters(
        2.8,   // A
        1.0,   // alpha
        1.6    // beta
    );

    fRatio->SetParLimits(0, 0.0, 20.0);   // A >= 0
    fRatio->SetParLimits(1, 0.2, 8.0);    // alpha
    fRatio->SetParLimits(2, 0.2, 8.0);    // beta

    fRatio->SetLineColor(kBlue + 1);
    fRatio->SetLineWidth(3);

    std::cout << "\n\n========== FIT SIMPLE POSITIVE RATIO ==========\n";
    gRatioFit->Fit(fRatio, "RME");

    // ========================================================================
    // Canvas 2: original ratio + smooth ratio fit
    // ========================================================================

    TCanvas *c2 =
        new TCanvas("c2", "pi0 smooth constrained ratio EXP/SIM", 1200, 720);

    c2->SetGrid();

    gRatioOriginal->Draw("AP");
    gRatioOriginal->GetXaxis()->SetLimits(0.0, 15.2);
    gRatioOriginal->GetYaxis()->SetRangeUser(0.8, 1.9);

    gRatioFit->Draw("P SAME");
    fRatio->Draw("SAME");

    TLine *line1 = new TLine(0.0, 1.0, 15.2, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray + 2);
    line1->SetLineWidth(2);
    line1->Draw("SAME");

    TLine *lLow = new TLine(X_RATIO_LEFT, 0.8, X_RATIO_LEFT, 1.9);
    lLow->SetLineStyle(3);
    lLow->SetLineColor(kGray + 1);
    lLow->Draw("SAME");

    TLine *lHigh = new TLine(X_RATIO_FLAT, 0.8, X_RATIO_FLAT, 1.9);
    lHigh->SetLineStyle(3);
    lHigh->SetLineColor(kGray + 1);
    lHigh->Draw("SAME");

    TLegend *leg2 = new TLegend(0.45, 0.63, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg2->AddEntry(gRatioOriginal, "original point ratio EXP/SIM", "lep");
    leg2->AddEntry(gRatioFit, "points used in ratio fit", "lep");
    leg2->AddEntry(fRatio, "smooth positive beta-bump ratio", "l");
    leg2->AddEntry(line1, "ratio = 1", "l");
    leg2->AddEntry(lHigh, "R=1 for OA #leq 0.5 and OA #geq 6 deg", "l");

    leg2->Draw();

    c2->Update();
    c2->SaveAs("pi0_ratio_smooth_beta_bump.png");

    // ========================================================================
    // Canvas 3: final scale function only
    // ========================================================================

    TCanvas *c3 =
        new TCanvas("c3", "pi0 final smooth scale function", 1200, 620);

    c3->SetGrid();

    TF1 *fScale = (TF1*)fRatio->Clone("fScale");

    fScale->SetTitle(
        "#pi^{0} final smooth scale function;OA(e^{+}e^{-}) [deg];scale factor R(OA)"
    );

    fScale->SetMinimum(0.95);
    fScale->SetMaximum(1.9);
    fScale->Draw();

    TLine *line1b = new TLine(0.0, 1.0, 15.2, 1.0);
    line1b->SetLineStyle(2);
    line1b->SetLineColor(kGray + 2);
    line1b->SetLineWidth(2);
    line1b->Draw("SAME");

    TLine *lLow2 = new TLine(X_RATIO_LEFT, 0.95, X_RATIO_LEFT, 1.9);
    lLow2->SetLineStyle(3);
    lLow2->SetLineColor(kGray + 1);
    lLow2->Draw("SAME");

    TLine *lHigh2 = new TLine(X_RATIO_FLAT, 0.95, X_RATIO_FLAT, 1.9);
    lHigh2->SetLineStyle(3);
    lHigh2->SetLineColor(kGray + 1);
    lHigh2->Draw("SAME");

    c3->Update();
    c3->SaveAs("pi0_final_scale_function_beta_bump.png");

    // ========================================================================
    // Summary
    // ========================================================================

    std::cout << "\n============================================================\n";
    std::cout << "PI0 summary\n";
    std::cout << "============================================================\n";

    std::cout << "Yield parametrization: exp(Chebyshev polynomial in log(OA))\n";
    std::cout << "Yield fit starts at OA = " << X_MIN_FIT << "\n";
    std::cout << "YIELD_DEG = " << YIELD_DEG << "\n";
    std::cout << "Ratio parametrization: beta-bump, 3 parameters\n";
    std::cout << "Common anchors: first three bins and OA >= 6 deg\n\n";

    std::cout << "EXP log-yield chi2/ndf = "
              << fLogExp->GetChisquare()
              << " / " << fLogExp->GetNDF()
              << " = " << fLogExp->GetChisquare() / fLogExp->GetNDF()
              << "\n";

    std::cout << "SIM log-yield chi2/ndf = "
              << fLogSim->GetChisquare()
              << " / " << fLogSim->GetNDF()
              << " = " << fLogSim->GetChisquare() / fLogSim->GetNDF()
              << "\n";

    std::cout << "Ratio beta-bump chi2/ndf = "
              << fRatio->GetChisquare()
              << " / " << fRatio->GetNDF()
              << " = " << fRatio->GetChisquare() / fRatio->GetNDF()
              << "\n";

    std::cout << "\nFinal ratio parameters:\n";
    std::cout << "A     = " << fRatio->GetParameter(0)
              << " +/- " << fRatio->GetParError(0) << "\n";
    std::cout << "alpha = " << fRatio->GetParameter(1)
              << " +/- " << fRatio->GetParError(1) << "\n";
    std::cout << "beta  = " << fRatio->GetParameter(2)
              << " +/- " << fRatio->GetParError(2) << "\n";

    std::cout << "\nSaved files:\n";
    std::cout << "  pi0_yield_analytic_logcheb_deg10.png\n";
    std::cout << "  pi0_ratio_smooth_beta_bump.png\n";
    std::cout << "  pi0_final_scale_function_beta_bump.png\n";
    std::cout << "============================================================\n\n";
}
