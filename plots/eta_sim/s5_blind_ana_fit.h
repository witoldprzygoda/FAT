// ============================================================================
// s5_blind_ana_fit.h  --  fit utilities of the blind analysis (label blind_ana)
// ============================================================================
// Header-only, used by s5_blind_ana_dev.C, s5_blind_ana_eta.C, s5_blind_ana_x17.C.
//
//   * Lin    : linear model y_i = sum_k a_k F_ik, fitted by iterated generalised
//              least squares with variance V_i = max(mu_i, 1) + sum_k a_k^2 T_ik
//              (mu_i = current model, T_ik = template variance of column k:
//              Barlow-Beeston-lite). With T = 0 the converged solution is the
//              binned Poisson maximum likelihood (identity link, IRLS).
//              Optional Gaussian constraints on single parameters.
//              Returns parameters, covariance, Baker-Cousins likelihood-ratio
//              chi2 and the weights w of any linear functional g^T a
//              (w = V^-1 F C g), used for the event-group (sandwich) error.
//   * Gaussian fraction, erfc step, bin-integrated shapes
//   * plotting style helpers (large fonts, EPS -> ps2pdf -dEPSCrop + PNG)
// Double precision everywhere.
// ============================================================================
#ifndef S5_BLIND_ANA_FIT_H
#define S5_BLIND_ANA_FIT_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"

namespace bla {

// ---------------------------------------------------------------------------
// small dense linear algebra
// ---------------------------------------------------------------------------
typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;

inline Mat Zero(int n, int m) { return Mat(n, Vec(m, 0.)); }

// Cholesky inversion of a symmetric positive definite matrix (in place)
inline bool CholInv(Mat& a) {
    const int n = (int)a.size();
    Mat L = Zero(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = a[i][j];
            for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
            if (i == j) { if (!(s > 0)) return false; L[i][i] = std::sqrt(s); }
            else L[i][j] = s / L[j][j];
        }
    Mat Li = Zero(n, n);
    for (int i = 0; i < n; ++i) {
        Li[i][i] = 1.0 / L[i][i];
        for (int j = 0; j < i; ++j) {
            double s = 0;
            for (int k = j; k < i; ++k) s -= L[i][k] * Li[k][j];
            Li[i][j] = s / L[i][i];
        }
    }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = std::max(i, j); k < n; ++k) s += Li[k][i] * Li[k][j];
            a[i][j] = s;
        }
    return true;
}

// ---------------------------------------------------------------------------
// linear fit with model variance (IRLS = Poisson ML) and template variance
// ---------------------------------------------------------------------------
struct Constraint { int k; double val, sig; };

struct LinProblem {
    int n = 0, np = 0;
    Vec y;                 // data (counts) per point
    Mat F;                 // [point][param] model columns
    Mat T;                 // [point][param] template variance per unit a_k^2 (may be empty)
    std::vector<Constraint> cons;
    double floorVar = 1.0; // minimal Poisson variance
    bool poissonVar = true;// V = mu (true) or V = given var0 (false)
    Vec var0;              // fixed variance (if poissonVar == false)
    int maxIt = 30;
};

struct LinResult {
    bool ok = false;
    int np = 0, n = 0;
    Vec a;                 // parameters
    Mat C;                 // covariance (from the final variance)
    Vec V;                 // final variance per point
    Vec mu;                // final model per point
    double chi2 = 0;       // Baker-Cousins likelihood-ratio chi2 (Poisson part)
    double chi2P = 0;      // Pearson chi2 with the final variance
    int ndf = 0;
    double prob = 0;
    int nIt = 0;
    // weights of the functional g^T a: w_i = (V^-1 F C g)_i
    Vec Weights(const Vec& g) const {
        Vec Cg(np, 0.), w(n, 0.);
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) Cg[k] += C[k][l] * g[l];
        for (int i = 0; i < n; ++i) {
            double s = 0;
            for (int k = 0; k < np; ++k) s += Fc[i][k] * Cg[k];
            w[i] = s / V[i];
        }
        return w;
    }
    double Var(const Vec& g) const {
        double s = 0;
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) s += g[k] * C[k][l] * g[l];
        return s;
    }
    Mat Fc;                // copy of F (for the weights)
};

inline LinResult LinFit(const LinProblem& P, const Vec* start = nullptr) {
    LinResult r; r.np = P.np; r.n = P.n;
    const int n = P.n, np = P.np;
    if (n <= 0 || np <= 0) return r;
    r.Fc = P.F;
    Vec a(np, 1.0);
    if (start) a = *start;
    Vec V(n, 1.0);
    auto modelAt = [&](int i, const Vec& p) { double s = 0; for (int k = 0; k < np; ++k) s += p[k] * P.F[i][k]; return s; };
    // first variance: observed counts (or fixed)
    for (int i = 0; i < n; ++i) {
        double v = P.poissonVar ? std::max(P.y[i], P.floorVar) : P.var0[i];
        if (!P.T.empty()) for (int k = 0; k < np; ++k) v += a[k] * a[k] * P.T[i][k];
        V[i] = v;
    }
    Mat A; Vec b;
    bool conv = false;
    for (int it = 0; it < P.maxIt; ++it) {
        A = Zero(np, np); b.assign(np, 0.);
        for (int i = 0; i < n; ++i) {
            const double wi = 1.0 / V[i];
            for (int k = 0; k < np; ++k) {
                const double fk = P.F[i][k] * wi;
                if (fk == 0) continue;
                b[k] += fk * P.y[i];
                for (int l = 0; l <= k; ++l) A[k][l] += fk * P.F[i][l];
            }
        }
        for (const auto& c : P.cons) { A[c.k][c.k] += 1.0 / (c.sig * c.sig); b[c.k] += c.val / (c.sig * c.sig); }
        for (int k = 0; k < np; ++k) for (int l = k + 1; l < np; ++l) A[k][l] = A[l][k];
        // scale for conditioning
        Vec sc(np, 1.0);
        for (int k = 0; k < np; ++k) sc[k] = A[k][k] > 0 ? 1.0 / std::sqrt(A[k][k]) : 1.0;
        Mat As = A;
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) As[k][l] *= sc[k] * sc[l];
        if (!CholInv(As)) return r;
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) As[k][l] *= sc[k] * sc[l];
        Vec an(np, 0.);
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) an[k] += As[k][l] * b[l];
        double dmax = 0;
        for (int k = 0; k < np; ++k) dmax = std::max(dmax, std::fabs(an[k] - a[k]) / std::sqrt(std::max(As[k][k], 1e-300)));
        a = an; r.C = As; r.nIt = it + 1;
        // new variance
        for (int i = 0; i < n; ++i) {
            double v = P.poissonVar ? std::max(modelAt(i, a), P.floorVar) : P.var0[i];
            if (!P.T.empty()) for (int k = 0; k < np; ++k) v += a[k] * a[k] * P.T[i][k];
            V[i] = v;
        }
        if (it > 0 && dmax < 1e-6) { conv = true; break; }
    }
    (void)conv;
    r.a = a; r.V = V; r.mu.assign(n, 0.);
    r.chi2 = 0; r.chi2P = 0;
    for (int i = 0; i < n; ++i) {
        const double m = modelAt(i, a);
        r.mu[i] = m;
        const double mm = std::max(m, 1e-6);
        r.chi2 += 2.0 * (mm - P.y[i] + (P.y[i] > 0 ? P.y[i] * std::log(P.y[i] / mm) : 0.0));
        r.chi2P += (P.y[i] - m) * (P.y[i] - m) / V[i];
    }
    r.ndf = n - np + (int)P.cons.size();
    r.prob = r.ndf > 0 ? TMath::Prob(r.chi2P, r.ndf) : 0;
    r.ok = true;
    return r;
}

// ---------------------------------------------------------------------------
// shapes
// ---------------------------------------------------------------------------
inline double GaussFrac(double x1, double x2, double m, double s) {
    const double q = 1.0 / (std::sqrt(2.0) * s);
    return 0.5 * (std::erf((x2 - m) * q) - std::erf((x1 - m) * q));
}
// bifurcated Gaussian, unit area, integral over [x1, x2]
inline double BifurFrac(double x1, double x2, double m, double sL, double sR) {
    const double nL = sL / (sL + sR), nR = sR / (sL + sR);   // area left / right of m
    double s = 0;
    // left part: 2 nL * N(m, sL) restricted to x < m
    {
        const double a = std::min(x1, m), b = std::min(x2, m);
        if (b > a) s += 2.0 * nL * GaussFrac(a, b, m, sL);
    }
    {
        const double a = std::max(x1, m), b = std::max(x2, m);
        if (b > a) s += 2.0 * nR * GaussFrac(a, b, m, sR);
    }
    return s;
}
// step 0.5 erfc((x - e)/(sqrt2 s)) averaged over a bin [x1, x2] (Simpson, 8 points)
inline double StepBin(double x1, double x2, double e, double s) {
    const int n = 8; double sum = 0;
    for (int i = 0; i < n; ++i) {
        const double x = x1 + (i + 0.5) * (x2 - x1) / n;
        sum += 0.5 * std::erfc((x - e) / (std::sqrt(2.0) * s));
    }
    return sum / n;
}

// content of a histogram between x1 and x2 with fractional edge bins (uniform inside a bin)
inline double FracIntegral(const TH1D* h, double x1, double x2, double* err = nullptr) {
    double s = 0, e2 = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b), hi = h->GetXaxis()->GetBinUpEdge(b);
        const double a = std::max(lo, x1), c = std::min(hi, x2);
        if (c <= a) continue;
        const double f = (c - a) / (hi - lo);
        s += f * h->GetBinContent(b);
        e2 += f * f * h->GetBinError(b) * h->GetBinError(b);
    }
    if (err) *err = std::sqrt(e2);
    return s;
}

// ---------------------------------------------------------------------------
// figure style: large fonts, vector PDF with tight page + PNG
// ---------------------------------------------------------------------------
inline void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.050, "XYZ");
    gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetTitleOffset(1.10, "X");
    gStyle->SetTitleOffset(1.25, "Y");
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendTextSize(0.045);
    gStyle->SetEndErrorSize(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);
}
// save as EPS -> PDF (ps2pdf -dEPSCrop) and PNG; remove the EPS
inline void SaveFig(TCanvas* c, const std::string& base) {
    c->SaveAs((base + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", base.c_str(), base.c_str(), base.c_str()));
    c->SaveAs((base + ".png").c_str());
}

} // namespace bla

#endif
