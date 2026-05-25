// ============================================================================
// fit.C
//
// Makro ROOT do:
//   1. fitu EXP i SIM z uwspólnionymi punktami kotwiczącymi,
//   2. narysowania f_EXP i f_SIM,
//   3. wyznaczenia finalnej funkcji skalującej EXP/SIM przez bezpośredni fit ratio.
//
// Najważniejsza część:
//   Finalne ratio NIE jest brane jako f_EXP/f_SIM,
//   tylko jest fitowane bezpośrednio funkcją dodatniego garbu:
//
//      R(x) = 1,                                      x <= xL
//      R(x) = 1 + N * t^a * (1-t)^b,                  xL < x < xR
//      R(x) = 1,                                      x >= xR
//
//   gdzie:
//      xL = 0.25 deg
//      xR = 5.00 deg
//      t  = (x - xL)/(xR - xL)
//
// Własności finalnego ratio:
//   - R(0.25) = 1,
//   - R(x >= 5 deg) = 1,
//   - R(x) >= 1 wszędzie,
//   - gładki bump tylko w zakresie 0.25-5 deg.
//
// Uruchomienie:
//   root -l fit.C
//
// albo:
//   root -l
//   .L fit.C
//   fit()
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
// Dane wejściowe
// ============================================================================

static const int NPTS = 30;

double OA[NPTS] = {
    0.25,  0.75,  1.25,  1.75,  2.25,
    2.75,  3.25,  3.75,  4.25,  4.75,
    5.25,  5.75,  6.25,  6.75,  7.25,
    7.75,  8.25,  8.75,  9.25,  9.75,
    10.25, 10.75, 11.25, 11.75, 12.25,
    12.75, 13.25, 13.75, 14.25, 14.75
};

double EX[NPTS] = {
    0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25
};

double Y_EXP[NPTS] = {
    1444.467,  5177.358, 10008.392, 11708.789, 11122.247,
    9232.160,  7268.656,  4278.196,  3491.841,  3275.125,
    2166.496,  1873.992,  2129.470,  1741.683,  1235.898,
    1317.266,  1414.660,  1254.628,  1236.897,   879.671,
     886.802,  1175.600,   902.032,   684.849,   864.040,
     676.138,   839.332,   654.945,   996.267,   973.130
};

double EY_EXP[NPTS] = {
     61.131, 124.089, 168.831, 186.743, 181.433,
    166.499, 145.887, 119.113, 107.443,  98.975,
     80.864,  73.560,  74.793,  65.901,  58.532,
     60.117,  59.666,  55.937,  54.277,  49.900,
     46.840,  52.010,  49.071,  44.463,  50.150,
     42.450,  44.306,  40.398,  45.760,  45.508
};

double Y_SIM[NPTS] = {
    1077.131, 4818.426, 7794.102, 8835.195, 7793.852,
    6455.475, 5149.369, 4211.167, 3370.398, 2948.859,
    2505.472, 2225.088, 1954.286, 1770.050, 1556.445,
    1528.479, 1371.111, 1217.892, 1175.670, 1112.282,
    1062.485, 1012.026,  978.710,  888.434,  818.776,
     792.706,  806.860,  712.259,  693.663,  715.831
};

double EY_SIM[NPTS] = {
    12.031, 24.686, 31.855, 34.136, 31.993,
    29.119, 26.078, 23.379, 20.667, 19.044,
    17.536, 16.638, 15.602, 14.520, 13.720,
    13.408, 12.596, 11.989, 11.630, 11.151,
    10.801, 10.495, 10.253,  9.683,  9.172,
     9.068,  8.955,  8.382,  8.248,  8.199
};

// ============================================================================
// Tablice używane do fitu po zastosowaniu wspólnych punktów kotwiczących
// ============================================================================

double Y_EXP_FIT[NPTS];
double EY_EXP_FIT[NPTS];

double Y_SIM_FIT[NPTS];
double EY_SIM_FIT[NPTS];

// ============================================================================
// Średnia ważona dwóch punktów
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

// ============================================================================
// Przygotowanie punktów do fitu
//
// Zasada:
//   - i == 0        : wspólny punkt EXP/SIM,
//   - OA[i] >= 5.0 : wspólny punkt EXP/SIM,
//   - reszta       : oryginalny EXP albo SIM.
//
// Czyli różnice między EXP i SIM zostają w obszarze:
//
//   0.5 < OA < 5.0 deg
//
// ============================================================================

void BuildCommonAnchorFitArrays()
{
    std::cout << "\n============================================================\n";
    std::cout << "Preparing fit points with common anchors\n";
    std::cout << "Common anchors: first bin and OA >= 5 deg\n";
    std::cout << "Independent EXP/SIM region: 0.5 < OA < 5 deg\n";
    std::cout << "============================================================\n\n";

    for (int i = 0; i < NPTS; ++i)
    {
        const bool useCommonPoint = (i == 0) || (OA[i] >= 5.0);

        if (useCommonPoint)
        {
            double yavg  = 0.0;
            double eyavg = 0.0;

            WeightedAverage2(
                Y_EXP[i], EY_EXP[i],
                Y_SIM[i], EY_SIM[i],
                yavg, eyavg
            );

            Y_EXP_FIT[i]  = yavg;
            EY_EXP_FIT[i] = eyavg;

            Y_SIM_FIT[i]  = yavg;
            EY_SIM_FIT[i] = eyavg;

            std::cout << "OA = " << OA[i]
                      << "  common = " << yavg
                      << " +/- " << eyavg
                      << "   from EXP " << Y_EXP[i] << " +/- " << EY_EXP[i]
                      << " and SIM " << Y_SIM[i] << " +/- " << EY_SIM[i]
                      << "\n";
        }
        else
        {
            Y_EXP_FIT[i]  = Y_EXP[i];
            EY_EXP_FIT[i] = EY_EXP[i];

            Y_SIM_FIT[i]  = Y_SIM[i];
            EY_SIM_FIT[i] = EY_SIM[i];
        }
    }

    std::cout << "\n";
}

// ============================================================================
// Model fitu yield(OA)
//
// par[0] = A      - amplituda części peakowej
// par[1] = p      - potęga x^p
// par[2] = k      - stała spadku exp(-k*x)
//
// par[3] = D      - amplituda ogona potęgowego
// par[4] = s      - przesunięcie w ogonie: (x+s)^q
// par[5] = q      - potęga ogona
//
// par[6] = B      - wspólny poziom bazowy / asymptota
//
// par[7] = x0     - środek przejścia peak -> tail
// par[8] = Delta  - szerokość przejścia
//
// f_left  = B + A*x^p*exp(-k*x)
// f_right = B + D/(x+s)^q
//
// w(x) = 1/(1 + exp((x-x0)/Delta))
//
// f(x) = w*f_left + (1-w)*f_right
//
// ============================================================================

double SmoothPiecewiseModel(double *xx, double *par)
{
    const double x = xx[0];

    const double A     = par[0];
    const double p     = par[1];
    const double k     = par[2];

    const double D     = par[3];
    const double s     = par[4];
    const double q     = par[5];

    const double B     = par[6];

    const double x0    = par[7];
    const double Delta = par[8];

    if (x <= 0.0) return 0.0;
    if (Delta <= 0.0) return 0.0;
    if (x + s <= 0.0) return 0.0;

    const double arg = (x - x0) / Delta;
    const double w = 1.0 / (1.0 + std::exp(arg));

    const double left  = B + A * std::pow(x, p) * std::exp(-k * x);
    const double right = B + D / std::pow(x + s, q);

    return w * left + (1.0 - w) * right;
}

// ============================================================================
// Bezpośredni model ratio
//
// par[0] = N  - amplituda garbu, N >= 0
// par[1] = a  - kształt lewego zbocza
// par[2] = b  - kształt prawego zbocza
//
// R(x) = 1,                                      x <= xL
// R(x) = 1 + N * t^a * (1-t)^b,                  xL < x < xR
// R(x) = 1,                                      x >= xR
//
// xL = 0.25
// xR = 5.00
// t = (x-xL)/(xR-xL)
//
// Dzięki temu:
//   R(x) >= 1 wszędzie,
//   R(0.25) = 1,
//   R(x >= 5) = 1.
//
// ============================================================================

double PositiveRatioBump(double *xx, double *par)
{
    const double x = xx[0];

    const double N = par[0];
    const double a = par[1];
    const double b = par[2];

    const double xL = 0.25;
    const double xR = 5.00;

    if (x <= xL) return 1.0;
    if (x >= xR) return 1.0;

    const double t = (x - xL) / (xR - xL);

    if (t <= 0.0) return 1.0;
    if (t >= 1.0) return 1.0;

    return 1.0 + N * std::pow(t, a) * std::pow(1.0 - t, b);
}

// ============================================================================
// Konfiguracja funkcji fitującej yield
// ============================================================================

void ConfigureFitFunction(TF1 *f, bool isExp)
{
    f->SetParNames(
        "A",
        "p",
        "k",
        "D",
        "s",
        "q",
        "B",
        "x0",
        "Delta"
    );

    if (isExp)
    {
        f->SetParameters(
            3.6e4,  // A
            3.0,    // p
            1.6,    // k
            9.0e4,  // D
            0.05,   // s
            2.35,   // q
            600.0,  // B
            4.0,    // x0
            0.45    // Delta
        );
    }
    else
    {
        f->SetParameters(
            2.5e4,  // A
            2.8,    // p
            1.55,   // k
            8.0e4,  // D
            0.05,   // s
            2.35,   // q
            600.0,  // B
            4.0,    // x0
            0.45    // Delta
        );
    }

    f->SetParLimits(0, 0.0, 1.0e7);    // A
    f->SetParLimits(1, 0.1, 8.0);      // p
    f->SetParLimits(2, 0.01, 10.0);    // k

    f->SetParLimits(3, 0.0, 1.0e9);    // D
    f->SetParLimits(4, 0.001, 10.0);   // s
    f->SetParLimits(5, 0.1, 10.0);     // q

    f->SetParLimits(6, 0.0, 5000.0);   // B

    // Wspólna geometria przejścia dla EXP i SIM.
    // To stabilizuje fity yieldów.
    f->FixParameter(7, 4.0);   // x0
    f->FixParameter(8, 0.45);  // Delta
}

// ============================================================================
// Funkcja pomocnicza: chi2 ręcznie
// ============================================================================

double ManualChi2(
    TF1 *f,
    const double *x,
    const double *y,
    const double *ey,
    int n
)
{
    double chi2 = 0.0;

    for (int i = 0; i < n; ++i)
    {
        const double yf = f->Eval(x[i]);
        const double pull = (y[i] - yf) / ey[i];
        chi2 += pull * pull;
    }

    return chi2;
}

// ============================================================================
// Główne makro
// ============================================================================

void fit()
{
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    BuildCommonAnchorFitArrays();

    const double fitMin = 0.25;
    const double fitMax = 14.75;

    // ========================================================================
    // Grafy oryginalne
    // ========================================================================

    TGraphErrors *gExpOriginal =
        new TGraphErrors(NPTS, OA, Y_EXP, EX, EY_EXP);

    gExpOriginal->SetName("gExpOriginal");
    gExpOriginal->SetTitle(
        "EXP/SIM yield vs OA;OA(e^{+}e^{-}) [deg];yield(#mu #pm 1#sigma)"
    );

    gExpOriginal->SetMarkerStyle(20);
    gExpOriginal->SetMarkerSize(0.9);
    gExpOriginal->SetMarkerColor(kBlack);
    gExpOriginal->SetLineColor(kBlack);

    TGraphErrors *gSimOriginal =
        new TGraphErrors(NPTS, OA, Y_SIM, EX, EY_SIM);

    gSimOriginal->SetName("gSimOriginal");
    gSimOriginal->SetMarkerStyle(24);
    gSimOriginal->SetMarkerSize(0.9);
    gSimOriginal->SetMarkerColor(kRed + 1);
    gSimOriginal->SetLineColor(kRed + 1);

    // ========================================================================
    // Grafy używane do fitu yieldów
    // ========================================================================

    TGraphErrors *gExpFit =
        new TGraphErrors(NPTS, OA, Y_EXP_FIT, EX, EY_EXP_FIT);

    gExpFit->SetName("gExpFit");
    gExpFit->SetMarkerStyle(20);
    gExpFit->SetMarkerSize(0.75);
    gExpFit->SetMarkerColor(kBlue + 1);
    gExpFit->SetLineColor(kBlue + 1);

    TGraphErrors *gSimFit =
        new TGraphErrors(NPTS, OA, Y_SIM_FIT, EX, EY_SIM_FIT);

    gSimFit->SetName("gSimFit");
    gSimFit->SetMarkerStyle(24);
    gSimFit->SetMarkerSize(0.75);
    gSimFit->SetMarkerColor(kMagenta + 1);
    gSimFit->SetLineColor(kMagenta + 1);

    // ========================================================================
    // Fit EXP
    // ========================================================================

    TF1 *fExp = new TF1(
        "fExp",
        SmoothPiecewiseModel,
        fitMin,
        fitMax,
        9
    );

    ConfigureFitFunction(fExp, true);

    fExp->SetLineColor(kBlue + 1);
    fExp->SetLineWidth(3);

    std::cout << "\n\n========== FIT EXP, common anchors ==========\n";
    gExpFit->Fit(fExp, "RME");

    // ========================================================================
    // Fit SIM
    // ========================================================================

    TF1 *fSim = new TF1(
        "fSim",
        SmoothPiecewiseModel,
        fitMin,
        fitMax,
        9
    );

    ConfigureFitFunction(fSim, false);

    fSim->SetLineColor(kRed + 1);
    fSim->SetLineWidth(3);

    std::cout << "\n\n========== FIT SIM, common anchors ==========\n";
    gSimFit->Fit(fSim, "RME");

    // ========================================================================
    // Canvas 1: oryginalne punkty + punkty fitowe + oba fity
    // ========================================================================

    TCanvas *c1 =
        new TCanvas("c1", "EXP and SIM fits with common anchors", 1150, 800);

    c1->SetGrid();

    gExpOriginal->Draw("AP");
    gExpOriginal->GetXaxis()->SetLimits(0.0, 15.2);
    gExpOriginal->GetYaxis()->SetRangeUser(0.0, 13000.0);

    gSimOriginal->Draw("P SAME");

    // Punkty, które faktycznie poszły do fitu.
    gExpFit->Draw("P SAME");
    gSimFit->Draw("P SAME");

    fExp->Draw("SAME");
    fSim->Draw("SAME");

    TLegend *leg1 = new TLegend(0.50, 0.58, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    leg1->AddEntry(gExpOriginal, "EXP original", "lep");
    leg1->AddEntry(gSimOriginal, "SIM original", "lep");
    leg1->AddEntry(gExpFit, "EXP points used in fit", "p");
    leg1->AddEntry(gSimFit, "SIM points used in fit", "p");
    leg1->AddEntry(fExp, "fit EXP with common anchors", "l");
    leg1->AddEntry(fSim, "fit SIM with common anchors", "l");

    leg1->Draw();

    c1->Update();
    c1->SaveAs("fit_exp_sim_commonAnchors.png");

    // ========================================================================
    // Punktowe ratio z oryginalnych danych
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
        "Ratio EXP/SIM;OA(e^{+}e^{-}) [deg];EXP / SIM"
    );

    gRatioOriginal->SetMarkerStyle(20);
    gRatioOriginal->SetMarkerSize(0.9);
    gRatioOriginal->SetMarkerColor(kBlack);
    gRatioOriginal->SetLineColor(kBlack);

    // ========================================================================
    // Graf ratio używany do bezpośredniego fitu
    //
    // Bierzemy punkty do x <= 5.0.
    // Punkt x=0.25 zastępujemy wartością 1, bo to jest wymuszony anchor.
    //
    // Powyżej 5.0 finalna funkcja ratio ma być z definicji 1,
    // więc tych punktów NIE używamy do fitu garbu.
    // ========================================================================

    double xRatioFit[NPTS];
    double yRatioFit[NPTS];
    double exRatioFit[NPTS];
    double eyRatioFit[NPTS];

    int nRatioFit = 0;

    for (int i = 0; i < NPTS; ++i)
    {
        if (OA[i] <= 5.0)
        {
            xRatioFit[nRatioFit]  = OA[i];
            exRatioFit[nRatioFit] = EX[i];

            if (i == 0)
            {
                // Lewy anchor: ratio = 1 dokładnie.
                // Dajemy mały błąd, żeby fit mocno trzymał punkt startowy.
                yRatioFit[nRatioFit]  = 1.0;
                eyRatioFit[nRatioFit] = 0.01;
            }
            else
            {
                yRatioFit[nRatioFit]  = R_ORIG[i];
                eyRatioFit[nRatioFit] = ER_ORIG[i];
            }

            ++nRatioFit;
        }
    }

    TGraphErrors *gRatioFit =
        new TGraphErrors(nRatioFit, xRatioFit, yRatioFit, exRatioFit, eyRatioFit);

    gRatioFit->SetName("gRatioFit");
    gRatioFit->SetTitle(
        "Ratio points used in direct constrained fit;OA(e^{+}e^{-}) [deg];EXP / SIM"
    );

    gRatioFit->SetMarkerStyle(24);
    gRatioFit->SetMarkerSize(0.9);
    gRatioFit->SetMarkerColor(kRed + 1);
    gRatioFit->SetLineColor(kRed + 1);

    // ========================================================================
    // Bezpośredni fit ratio
    // ========================================================================

    TF1 *fRatioDirect = new TF1(
        "fRatioDirect",
        PositiveRatioBump,
        0.25,
        15.0,
        3
    );

    fRatioDirect->SetParNames("N", "a", "b");

    // Starty dobrane pod bump ok. 1-4 deg.
    fRatioDirect->SetParameters(
        6.0,   // N
        1.5,   // a
        1.5    // b
    );

    fRatioDirect->SetParLimits(0, 0.0, 100.0);  // N >= 0
    fRatioDirect->SetParLimits(1, 0.1, 20.0);   // a > 0
    fRatioDirect->SetParLimits(2, 0.1, 20.0);   // b > 0

    fRatioDirect->SetLineColor(kBlue + 1);
    fRatioDirect->SetLineWidth(3);

    std::cout << "\n\n========== DIRECT CONSTRAINED RATIO FIT ==========\n";
    gRatioFit->Fit(fRatioDirect, "RME");

    // ========================================================================
    // Canvas 2: ratio oryginalne + bezpośredni fit ratio
    // ========================================================================

    TCanvas *c2 =
        new TCanvas("c2", "Direct constrained ratio EXP/SIM", 1150, 700);

    c2->SetGrid();

    gRatioOriginal->Draw("AP");
    gRatioOriginal->GetXaxis()->SetLimits(0.0, 15.2);
    gRatioOriginal->GetYaxis()->SetRangeUser(0.8, 1.8);

    gRatioFit->Draw("P SAME");
    fRatioDirect->Draw("SAME");

    TLine *line1 = new TLine(0.0, 1.0, 15.2, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray + 2);
    line1->SetLineWidth(2);
    line1->Draw("SAME");

    // Granice regionu, gdzie finalna funkcja ratio może się różnić od 1.
    TLine *lLow = new TLine(0.25, 0.8, 0.25, 1.8);
    TLine *lHigh = new TLine(5.0, 0.8, 5.0, 1.8);

    lLow->SetLineStyle(3);
    lHigh->SetLineStyle(3);

    lLow->SetLineColor(kGray + 1);
    lHigh->SetLineColor(kGray + 1);

    lLow->Draw("SAME");
    lHigh->Draw("SAME");

    TLegend *leg2 = new TLegend(0.46, 0.64, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg2->AddEntry(gRatioOriginal, "original point ratio EXP/SIM", "lep");
    leg2->AddEntry(gRatioFit, "points used in direct ratio fit", "lep");
    leg2->AddEntry(fRatioDirect, "constrained ratio fit: R #geq 1, flat above 5 deg", "l");
    leg2->AddEntry(line1, "ratio = 1", "l");
    leg2->AddEntry(lHigh, "R=1 for OA #geq 5 deg", "l");

    leg2->Draw();

    c2->Update();
    c2->SaveAs("ratio_direct_positive_flatAbove5.png");

    // ========================================================================
    // Canvas 3: sama finalna funkcja skalująca
    // ========================================================================

    TCanvas *c3 =
        new TCanvas("c3", "Final smooth scale function", 1150, 600);

    c3->SetGrid();

    TF1 *fScale = (TF1*)fRatioDirect->Clone("fScale");

    fScale->SetTitle(
        "Final smooth scale function;OA(e^{+}e^{-}) [deg];scale factor R(OA)"
    );

    fScale->SetMinimum(0.95);
    fScale->SetMaximum(1.8);
    fScale->Draw();

    TLine *line1b = new TLine(0.0, 1.0, 15.2, 1.0);
    line1b->SetLineStyle(2);
    line1b->SetLineColor(kGray + 2);
    line1b->SetLineWidth(2);
    line1b->Draw("SAME");

    TLine *lLow2 = new TLine(0.25, 0.95, 0.25, 1.8);
    TLine *lHigh2 = new TLine(5.0, 0.95, 5.0, 1.8);

    lLow2->SetLineStyle(3);
    lHigh2->SetLineStyle(3);

    lLow2->SetLineColor(kGray + 1);
    lHigh2->SetLineColor(kGray + 1);

    lLow2->Draw("SAME");
    lHigh2->Draw("SAME");

    c3->Update();
    c3->SaveAs("final_scale_function_positive_flatAbove5.png");

    // ========================================================================
    // Canvas 4: pulle dla fitów yieldów
    // ========================================================================

    double PULL_EXP[NPTS];
    double PULL_SIM[NPTS];
    double PULL_EY[NPTS];

    for (int i = 0; i < NPTS; ++i)
    {
        PULL_EXP[i] = (Y_EXP_FIT[i] - fExp->Eval(OA[i])) / EY_EXP_FIT[i];
        PULL_SIM[i] = (Y_SIM_FIT[i] - fSim->Eval(OA[i])) / EY_SIM_FIT[i];
        PULL_EY[i]  = 0.0;
    }

    TGraphErrors *gPullExp =
        new TGraphErrors(NPTS, OA, PULL_EXP, EX, PULL_EY);

    gPullExp->SetName("gPullExp");
    gPullExp->SetTitle(
        "Pulls for yield fit points;OA(e^{+}e^{-}) [deg];(data-fit)/#sigma"
    );

    gPullExp->SetMarkerStyle(20);
    gPullExp->SetMarkerSize(0.85);
    gPullExp->SetMarkerColor(kBlack);
    gPullExp->SetLineColor(kBlack);

    TGraphErrors *gPullSim =
        new TGraphErrors(NPTS, OA, PULL_SIM, EX, PULL_EY);

    gPullSim->SetName("gPullSim");
    gPullSim->SetMarkerStyle(24);
    gPullSim->SetMarkerSize(0.85);
    gPullSim->SetMarkerColor(kRed + 1);
    gPullSim->SetLineColor(kRed + 1);

    TCanvas *c4 =
        new TCanvas("c4", "Pulls EXP and SIM yield fits", 1150, 560);

    c4->SetGrid();

    gPullExp->Draw("AP");
    gPullExp->GetXaxis()->SetLimits(0.0, 15.2);
    gPullExp->GetYaxis()->SetRangeUser(-8.0, 8.0);

    gPullSim->Draw("P SAME");

    TF1 *zero = new TF1("zero", "0", 0.0, 15.2);
    zero->SetLineColor(kBlue + 1);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("SAME");

    TLegend *leg4 = new TLegend(0.62, 0.74, 0.88, 0.88);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);

    leg4->AddEntry(gPullExp, "EXP fit-point pulls", "lep");
    leg4->AddEntry(gPullSim, "SIM fit-point pulls", "lep");

    leg4->Draw();

    c4->Update();
    c4->SaveAs("pulls_commonAnchors.png");

    // ========================================================================
    // Podsumowanie w terminalu
    // ========================================================================

    const double chi2ExpManual =
        ManualChi2(fExp, OA, Y_EXP_FIT, EY_EXP_FIT, NPTS);

    const double chi2SimManual =
        ManualChi2(fSim, OA, Y_SIM_FIT, EY_SIM_FIT, NPTS);

    const int nFreeExp = fExp->GetNumberFreeParameters();
    const int nFreeSim = fSim->GetNumberFreeParameters();

    const int ndfExp = NPTS - nFreeExp;
    const int ndfSim = NPTS - nFreeSim;

    std::cout << "\n============================================================\n";
    std::cout << "Fit summary\n";
    std::cout << "============================================================\n";

    std::cout << "EXP ROOT chi2/ndf   = "
              << fExp->GetChisquare()
              << " / " << fExp->GetNDF()
              << " = " << fExp->GetChisquare() / fExp->GetNDF()
              << "\n";

    std::cout << "EXP manual chi2/ndf = "
              << chi2ExpManual
              << " / " << ndfExp
              << " = " << chi2ExpManual / ndfExp
              << "\n\n";

    std::cout << "SIM ROOT chi2/ndf   = "
              << fSim->GetChisquare()
              << " / " << fSim->GetNDF()
              << " = " << fSim->GetChisquare() / fSim->GetNDF()
              << "\n";

    std::cout << "SIM manual chi2/ndf = "
              << chi2SimManual
              << " / " << ndfSim
              << " = " << chi2SimManual / ndfSim
              << "\n\n";

    std::cout << "Direct ratio fit chi2/ndf = "
              << fRatioDirect->GetChisquare()
              << " / " << fRatioDirect->GetNDF()
              << " = " << fRatioDirect->GetChisquare() / fRatioDirect->GetNDF()
              << "\n";

    std::cout << "\nFinal ratio parameters:\n";
    std::cout << "N = " << fRatioDirect->GetParameter(0)
              << " +/- " << fRatioDirect->GetParError(0) << "\n";
    std::cout << "a = " << fRatioDirect->GetParameter(1)
              << " +/- " << fRatioDirect->GetParError(1) << "\n";
    std::cout << "b = " << fRatioDirect->GetParameter(2)
              << " +/- " << fRatioDirect->GetParError(2) << "\n";

    std::cout << "\nSaved files:\n";
    std::cout << "  fit_exp_sim_commonAnchors.png\n";
    std::cout << "  ratio_direct_positive_flatAbove5.png\n";
    std::cout << "  final_scale_function_positive_flatAbove5.png\n";
    std::cout << "  pulls_commonAnchors.png\n";
    std::cout << "============================================================\n\n";
}
