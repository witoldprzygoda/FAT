// ============================================================================
// x17c_model.h  --  FROZEN X17 model of X-c (SPEC report_all WP9; prefix x17c; 2026-11-14)
// ============================================================================
// This file is x17_model.h of the x17ab package (namespace x17m, X-a/X-b, validated on toys) copied verbatim into
// namespace x17c, with FOUR additions (marked "x17c:") needed by the frozen procedure (x17_frozen.txt):
//   1. Config::xNorm  0 = normalisation (1): X17 rows = r N5 R_X        (r = BR_X BR_ee / BR_dir, N5 fitted)
//                     1 = normalisations (2)/(3): X17 rows = r xA N5ref R_X with a FIXED eta scale xA (a_eta of WP5,
//                         or the peak-based a_P) and N5ref = direct rows at SMASH (a_eta = 1); the uncertainty of xA
//                         goes into the X17 yield nuisance sX. (The a_eta parameter of the eta-related background terms
//                         is NOT used for the X17 normalisation, so the side-band fit cannot move it.)
//      The X17 yield nuisance is LOG-NORMAL in x17c: factor (1 + sX)^thX (x17_model.h: linear 1 + sX thX).
//   2. line-shape nuisances: X_j -> X_j + thW dXw_j + thS dXs_j, thW, thS ~ N(0, 1) (Config::wSig, sSig > 0 switch them
//      on; Inputs::dXw, dXs = half differences of the line at width x (1 +- 0.10) and M(ee) scale x (1 +- 0.02)).
//   3. Config::qScale: the test statistic q~ is divided by it (event-group variance factor phi, or the Pearson
//      chi2/ndf fall-back); ToyData then draws over-dispersed pseudo-data qScale x Poisson(mu / qScale).
//   4. Glob::w0, s0: centres of the line-shape constraints (global observables of the toys).
// Everything else (model, fitter, q~, asymptotic CLs, FillMaskedRows, builders, LineFor) is unchanged.
// --- original header of x17_model.h follows ---
// x17_model.h  --  X17 search: 2D model in M4 x M(ee), Poisson likelihood, 90 % CLs (asymptotic and toys)
// SPEC report_all, WP9 X-a (model) and X-b (validation); reused by X-c (side-band fit, Asimov, freeze) and X-d.
// Written by the x17ab package (2026-11-14). Namespace x17m (here: x17c). Header only; ACLiC-compatible.
//
// THE MODEL (bin i in M4, bin j in M(ee); all terms in expected RAW rows of the fitted sample):
//
//   mu_ij = k_e(j) CB_ij                                             combinatorial background, measured (like-sign)
//         + n_j  B_c(j)(i) [1 + b_c t_i (+ c_c t_i^2)]               other background: free yield per M(ee) bin,
//                                                                    M4 shape per class c, tilt per class
//         + gP_c(j)(i) [ N5 D_j + N4 C_j (1 + sC thC T_j) ]          eta peak: direct decays (D) and conversions (C)
//         + r N5 R_X (1 + sX thX) gX(i) X_j                          X17 line, POI r
//         + a_eta [ E23_ij + g Geta_ij ]                             eta-related non-peak terms (C2 + C3, G_eta)
//
//   t_i = (M4_i - 0.55)/0.1;  D_j, C_j, X_j = fractions of ALL direct / conversion / X17 rows of the M4 range that fall in
//   M(ee) bin j (the M(ee) flows count in the total, so N5 and N4 are totals over the full M(ee) range);
//   R_X = eps_X / eps_5 (efficiency of eta -> pi+pi-X17(ee) over that of eta -> pi+pi-e+e-, same selection and range);
//   POI  r = BR(eta -> pi+pi- X17) BR(X17 -> e+e-) / BR(eta -> pi+pi- e+e-)   => BR_X BR_ee = r x 2.68e-4.
//
//   Constraints (Gaussian penalty terms in -2 ln L): k_E2 = 1 +- 0.05 (k_E1 free by default); a_eta = a0 +- sa;
//   g = 1 +- 0.30; thC ~ N(0, 1) (relative uncertainty sC of the conversion tail T_j = M(ee) >= 12 MeV);
//   thX ~ N(0, 1) (relative uncertainty sX of the X17 yield: eps_X/eps_5, N5 bias, BR_dir, ...).
//   Background-shape modes (Config::bkgMode): 0 B (1 + b t), b free per class (the literal spec form); 1 B fixed;
//   2 flat (1 + b t + c t^2) free per class ("pol2" variant); 3 B (1 + b t + c t^2);
//   4 B (1 + (b + b' u_j) t + (c + c' u_j) t^2), u_j = (M(ee)_j - 30 MeV)/30 MeV, all four free per class (NOMINAL of X-a/X-b:
//     the only form that passed the toy tests with a background shape that changes with M(ee); see work_all/x17ab_results.md).
//
// LIKELIHOOD: -2 ln L = 2 sum_used [mu - n ln mu] + sum_constraints ((p - p0)/sigma)^2, minimised by Fisher scoring
//   (Newton with the expected Hessian 2 sum J J^T / mu) with step halving; every parameter that has no information
//   (zero Fisher diagonal, e.g. n_j of a fully masked row) is frozen at its start value.
// TEST STATISTIC: q~_r (Cowan, Cranmer, Gross, Vitells, EPJC 71 (2011) 1554), r_hat >= 0 in the reference fit.
// LIMITS: CLs = CL_{s+b}/CL_b = 0.10.
//   asymptotic: sigma_A(r) from the background-only Asimov data set (the model at the r = 0 fit of the observed data);
//   toys: Poisson pseudo-data from the conditional fits at r and at 0; the centres of the Gaussian constraints are
//         drawn around the conditional estimates (global observables), q~_r computed for every toy.
//
// INPUT BUILDERS: BuildFromHists (hists_v1.root, data scale, layout X2D, data CB + simulation templates) and
//   BuildToyGrid (blind-test half-A templates of the simulation study, toy grid in simulated M4). Line of an arbitrary
//   mass m_X: LineFromFile (x17_lines.root of x17ab_line.C).
// ============================================================================
#ifndef ALLREPORT_X17C_MODEL_H
#define ALLREPORT_X17C_MODEL_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TTree.h"

namespace x17c {

const double kBRdir = 2.68e-4;               // PDG BR(eta -> pi+pi-e+e-)
const double kAlpha = 0.10;                  // CLs level (90 % CL)
const double kTailMee = 0.012;               // GeV/c2: the conversion-tail nuisance acts at M(ee) >= 12 MeV/c2
const char* const kLinesFile = "/hdd3/przygoda/HADES/FAT_out/all_report/x17ab/x17_lines.root";

// ============================================================================ grid
struct Grid {
    int nx = 0, ny = 0;
    std::vector<double> xe, ye;   // bin edges, GeV/c2
    double xc(int i) const { return 0.5 * (xe[i] + xe[i + 1]); }
    double yc(int j) const { return 0.5 * (ye[j] + ye[j + 1]); }
    int idx(int i, int j) const { return j * nx + i; }
    int n() const { return nx * ny; }
};

// ============================================================================ model inputs
struct Inputs {
    Grid G;
    std::vector<double> CB;      // [nx*ny] measured CB (k = 1), raw rows
    std::vector<double> E23;     // [nx*ny] C2 + C3 at a_eta = 1
    std::vector<double> Geta;    // [nx*ny] G_eta at a_eta = 1
    std::vector<int> kcl;        // [ny] CB class: 0 E1 (M(ee) < 10 MeV), 1 E2
    std::vector<int> bcl;        // [ny] background- and peak-shape class
    int nbcl = 2;
    std::vector<std::vector<double>> B;    // [nbcl][nx] background M4 shape (unit sum)
    std::vector<std::vector<double>> gP;   // [nbcl][nx] peak M4 shape (unit sum)
    std::vector<double> gX;                // [nx] X17 M4 shape (unit sum)
    std::vector<double> D, C, X;           // [ny] fractions of the totals (flows included in the totals)
    std::vector<double> Ctail;             // [ny] 1 where the conversion-tail nuisance acts
    double RX = 1.;                        // eps_X / eps_5
    std::vector<double> tx;                // [nx] tilt variable
    double N5start = 1000., N4start = 1000.;
    std::string label;
    double N5ref = 0.;                     // x17c: direct rows at a_eta = 1 (SMASH) for normalisation (2)
    std::vector<double> dXw, dXs;          // x17c: [ny] line-shape derivatives (width, M(ee) scale), per 1 sigma
};

// ============================================================================ configuration
struct Config {
    double kSig[2] = {-1., 0.05};   // per CB class: < 0 free, 0 fixed at k0, > 0 Gaussian constraint
    double k0[2] = {1., 1.};
    double a0 = 1., aSig = 0.10;    // a_eta (aSig < 0 free, 0 fixed)
    double gSig = 0.30;             // G_eta factor around 1 (0 fixed)
    int bkgMode = 0;                // 0 B(1+bt); 1 B fixed; 2 flat(1+bt+ct^2); 3 B(1+bt+ct^2); 4 B(1+(b+b'u)t+(c+c'u)t^2), u = (M(ee)-30 MeV)/30 MeV
    double cTailSig = 0.;           // relative uncertainty of the conversion tail (<= 0: no nuisance)
    double xSig = 0.;               // relative uncertainty of the X17 yield (<= 0: no nuisance)
    std::vector<char> use;          // [nx*ny] bins in the likelihood (empty = all)
    bool allowNegR = true;          // r_hat may be negative in the free fit
    int xNorm = 0;                  // x17c: 0 X17 rows = r N5 R_X (normalisation 1); 1 = r xA N5ref R_X (normalisations 2, 3)
    double xA = 1.;                 // x17c: fixed eta scale of normalisations (2)/(3)
    double wSig = 0., sSig = 0.;    // x17c: > 0 switches the line-width / M(ee)-scale nuisance on (size is in Inputs::dXw, dXs)
    double qScale = 1.;             // x17c: q~ divided by qScale (variance factor); toys over-dispersed by the same factor
    double gC = 1.;                 // x17c: centre of the g constraint (WP5 value); GlobOf copies it
};
// centres of the Gaussian constraints (global observables)
struct Glob {
    double k0[2] = {1., 1.};
    double a0 = 1., g0 = 1., c0 = 0., x0 = 0.;
    double w0 = 0., s0 = 0.;        // x17c: line-shape constraints
};
inline Glob GlobOf(const Config& c) { Glob g; g.k0[0] = c.k0[0]; g.k0[1] = c.k0[1]; g.a0 = c.a0; g.g0 = c.gC; return g; }

// ============================================================================ parameters
// layout: n_j (ny) | N5 | N4 | k_E1 | k_E2 | b_c (nbcl) | c_c (nbcl) | b'_c (nbcl) | c'_c (nbcl) | a | g | thC | thX | thW | thS | r
struct Layout {
    int ny = 0, nb = 0;
    int iN5 = 0, iN4 = 0, iK = 0, iB = 0, iC = 0, iBy = 0, iCy = 0, iA = 0, iG = 0, iTC = 0, iTX = 0, iTW = 0, iTS = 0, iR = 0, np = 0;
    void Set(int ny_, int nb_) {
        ny = ny_; nb = nb_;
        iN5 = ny; iN4 = ny + 1; iK = ny + 2; iB = ny + 4; iC = iB + nb; iBy = iC + nb; iCy = iBy + nb; iA = iCy + nb; iG = iA + 1; iTC = iA + 2; iTX = iA + 3;
        iTW = iA + 4; iTS = iA + 5; iR = iA + 6;
        np = iR + 1;
    }
};
struct FitResult {
    bool ok = false;
    int iter = 0;
    double F = std::numeric_limits<double>::infinity();   // -2 ln L (up to a data constant)
    std::vector<double> p;
    std::vector<double> err;    // sqrt of the diagonal of the inverse Fisher matrix (free parameters)
    double rErr = 0.;
    std::string msg;
};

class Model {
public:
    const Inputs* in = nullptr;
    Config cfg;
    Layout L;
    explicit Model(const Inputs& i, const Config& c = Config()) : in(&i), cfg(c) { L.Set(i.G.ny, i.nbcl); }

    bool Used(int b) const { return cfg.use.empty() || cfg.use[b]; }
    // which parameters float (r handled by the caller through rFree)
    std::vector<char> FreeMask(bool rFree) const {
        std::vector<char> f(L.np, 0);
        for (int j = 0; j < L.ny; ++j) f[j] = 1;
        f[L.iN5] = 1; f[L.iN4] = 1;
        for (int e = 0; e < 2; ++e) f[L.iK + e] = cfg.kSig[e] != 0.;
        const bool tilt = cfg.bkgMode != 1, quad = cfg.bkgMode >= 2, lin = cfg.bkgMode == 4;
        for (int c = 0; c < L.nb; ++c) { f[L.iB + c] = tilt; f[L.iC + c] = quad; f[L.iBy + c] = lin; f[L.iCy + c] = lin; }
        f[L.iA] = cfg.aSig != 0.;
        f[L.iG] = cfg.gSig != 0.;
        f[L.iTC] = cfg.cTailSig > 0.;
        f[L.iTX] = cfg.xSig > 0.;
        f[L.iTW] = cfg.wSig > 0. && !in->dXw.empty();
        f[L.iTS] = cfg.sSig > 0. && !in->dXs.empty();
        f[L.iR] = rFree;
        return f;
    }
    // x17c: effective line fraction of M(ee) bin j and its derivatives
    inline double Xeff(const std::vector<double>& p, int j) const {
        const Inputs& I = *in;
        double x = I.X[j];
        if (cfg.wSig > 0 && !I.dXw.empty()) x += p[L.iTW] * I.dXw[j];
        if (cfg.sSig > 0 && !I.dXs.empty()) x += p[L.iTS] * I.dXs[j];
        return x;
    }
    // x17c: the yield factor that multiplies r (N5 for normalisation 1, a_eta N5ref for normalisation 2)
    inline double Nnorm(const std::vector<double>& p) const { return cfg.xNorm == 1 ? cfg.xA * in->N5ref : p[L.iN5]; }
    std::vector<double> Start(const std::vector<double>& n, double r, const Glob& g) const {
        const Inputs& I = *in;
        std::vector<double> p(L.np, 0.);
        p[L.iN5] = I.N5start; p[L.iN4] = I.N4start;
        p[L.iK] = g.k0[0]; p[L.iK + 1] = g.k0[1];
        p[L.iA] = g.a0; p[L.iG] = g.g0; p[L.iTC] = 0; p[L.iTX] = 0; p[L.iTW] = 0; p[L.iTS] = 0; p[L.iR] = r;
        for (int j = 0; j < L.ny; ++j) {
            double s = 0, o = 0;
            const int c = I.bcl[j];
            for (int i = 0; i < I.G.nx; ++i) {
                const int b = I.G.idx(i, j);
                if (!Used(b)) continue;
                s += n[b];
                o += p[L.iK + I.kcl[j]] * I.CB[b] + g.a0 * (I.E23[b] + I.Geta[b]) +
                     I.gP[c][i] * (I.N5start * I.D[j] + I.N4start * I.C[j]);
            }
            double sh = 0;
            for (int i = 0; i < I.G.nx; ++i) if (Used(I.G.idx(i, j))) sh += (cfg.bkgMode == 2 ? 1. / I.G.nx : I.B[c][i]);
            p[j] = sh > 0 ? std::max(1e-3 * std::max(s, 1.), (s - o)) / sh : 1.;
        }
        return p;
    }
    // expected content of bin (i, j) and (optionally) its derivatives (dense vector of size np, only touched entries)
    inline double Mu(const std::vector<double>& p, int i, int j, double* d = nullptr) const {
        const Inputs& I = *in;
        const int b = I.G.idx(i, j), c = I.bcl[j], ke = I.kcl[j];
        const double t = I.tx[i];
        const double shp = cfg.bkgMode == 2 ? 1. / I.G.nx : I.B[c][i];
        double poly = 1.;
        const double u = cfg.bkgMode == 4 ? (I.G.yc(j) - 0.030) / 0.030 : 0.;
        if (cfg.bkgMode != 1) poly += (p[L.iB + c] + p[L.iBy + c] * u) * t;
        if (cfg.bkgMode >= 2) poly += (p[L.iC + c] + p[L.iCy + c] * u) * t * t;
        const double ctf = 1. + (cfg.cTailSig > 0 ? cfg.cTailSig * p[L.iTC] * I.Ctail[j] : 0.);
        const double lxs = cfg.xSig > 0 ? std::log(1. + cfg.xSig) : 0.;   // x17c: log-normal X17 yield nuisance
        const double xf = cfg.xSig > 0 ? std::exp(lxs * p[L.iTX]) : 1.;
        const double gp = I.gP[c][i];
        const double xl = I.gX[i] * Xeff(p, j) * I.RX;   // x17c: line-shape nuisances inside Xeff
        const double nn = Nnorm(p);                        // x17c: N5 (normalisation 1) or xA N5ref (normalisations 2, 3)
        const double mu = p[L.iK + ke] * I.CB[b] + p[j] * shp * poly + gp * (p[L.iN5] * I.D[j] + p[L.iN4] * I.C[j] * ctf) +
                          p[L.iR] * nn * xl * xf + p[L.iA] * (I.E23[b] + p[L.iG] * I.Geta[b]);
        if (d) {
            d[j] = shp * poly;
            d[L.iN5] = gp * I.D[j] + (cfg.xNorm == 1 ? 0. : p[L.iR] * xl * xf);
            d[L.iN4] = gp * I.C[j] * ctf;
            d[L.iK + ke] = I.CB[b];
            d[L.iB + c] = p[j] * shp * t;
            d[L.iC + c] = p[j] * shp * t * t;
            d[L.iBy + c] = p[j] * shp * t * u;
            d[L.iCy + c] = p[j] * shp * t * t * u;
            d[L.iA] = I.E23[b] + p[L.iG] * I.Geta[b];
            d[L.iG] = p[L.iA] * I.Geta[b];
            d[L.iTC] = cfg.cTailSig > 0 ? gp * p[L.iN4] * I.C[j] * cfg.cTailSig * I.Ctail[j] : 0.;
            d[L.iTX] = cfg.xSig > 0 ? p[L.iR] * nn * xl * xf * lxs : 0.;
            d[L.iTW] = (cfg.wSig > 0 && !I.dXw.empty()) ? p[L.iR] * nn * I.gX[i] * I.dXw[j] * I.RX * xf : 0.;
            d[L.iTS] = (cfg.sSig > 0 && !I.dXs.empty()) ? p[L.iR] * nn * I.gX[i] * I.dXs[j] * I.RX * xf : 0.;
            d[L.iR] = nn * xl * xf;
        }
        return mu;
    }
    // indices touched by bin (i, j)
    inline int Touched(int i, int j, int* idx) const {
        const Inputs& I = *in;
        int n = 0;
        idx[n++] = j; idx[n++] = L.iN5; idx[n++] = L.iN4; idx[n++] = L.iK + I.kcl[j];
        idx[n++] = L.iB + I.bcl[j]; idx[n++] = L.iC + I.bcl[j]; idx[n++] = L.iBy + I.bcl[j]; idx[n++] = L.iCy + I.bcl[j];
        idx[n++] = L.iA; idx[n++] = L.iG; idx[n++] = L.iTC; idx[n++] = L.iTX; idx[n++] = L.iTW; idx[n++] = L.iTS; idx[n++] = L.iR;
        (void)i;
        return n;
    }
    // constraint terms
    double Penalty(const std::vector<double>& p, const Glob& g, std::vector<double>* grad = nullptr, std::vector<double>* hdiag = nullptr) const {
        double F = 0;
        auto add = [&](int k, double c0, double s) {
            if (!(s > 0)) return;
            const double z = (p[k] - c0) / s;
            F += z * z;
            if (grad) (*grad)[k] += 2 * z / s;
            if (hdiag) (*hdiag)[k] += 2 / (s * s);
        };
        for (int e = 0; e < 2; ++e) add(L.iK + e, g.k0[e], cfg.kSig[e]);
        add(L.iA, g.a0, cfg.aSig);
        add(L.iG, g.g0, cfg.gSig);
        if (cfg.cTailSig > 0) add(L.iTC, g.c0, 1.);
        if (cfg.xSig > 0) add(L.iTX, g.x0, 1.);
        if (cfg.wSig > 0 && !in->dXw.empty()) add(L.iTW, g.w0, 1.);   // x17c
        if (cfg.sSig > 0 && !in->dXs.empty()) add(L.iTS, g.s0, 1.);   // x17c
        return F;
    }
    // -2 ln L (data constant omitted: sum n ln n - n added for readability of the deviance)
    double NLL(const std::vector<double>& p, const std::vector<double>& n, const Glob& g) const {
        const Inputs& I = *in;
        double F = 0;
        for (int j = 0; j < I.G.ny; ++j)
            for (int i = 0; i < I.G.nx; ++i) {
                const int b = I.G.idx(i, j);
                if (!Used(b)) continue;
                const double mu = Mu(p, i, j);
                const double nn = n[b];
                if (!(mu > 0)) { if (nn > 0 || mu < 0) return std::numeric_limits<double>::infinity(); continue; }
                F += 2 * (mu - nn + (nn > 0 ? nn * std::log(nn / mu) : 0.));
            }
        return F + Penalty(p, g);
    }
    // expected data (model prediction) for every bin
    std::vector<double> Expected(const std::vector<double>& p) const {
        const Inputs& I = *in;
        std::vector<double> m(I.G.n(), 0.);
        for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) m[I.G.idx(i, j)] = std::max(0., Mu(p, i, j));
        return m;
    }

    // ------------------------------------------------------------------ the fit
    FitResult Fit(const std::vector<double>& n, double r, bool rFree, const Glob& g, const std::vector<double>* start = nullptr,
                  int maxIter = 300) const {
        const Inputs& I = *in;
        FitResult R;
        std::vector<double> p = start ? *start : Start(n, r, g);
        if (!rFree) p[L.iR] = r;
        std::vector<char> fr = FreeMask(rFree);
        const int np = L.np;
        std::vector<double> d(np, 0.), grad(np), H((size_t)np * np), step(np);
        std::vector<double> hdiag(np);
        int idx[24];
        double F = NLL(p, n, g);
        if (!std::isfinite(F)) {   // restart from the default start
            p = Start(n, r, g);
            if (!rFree) p[L.iR] = r;
            F = NLL(p, n, g);
            if (!std::isfinite(F)) { R.msg = "start not finite"; R.p = p; return R; }
        }
        double lambda = 0.;   // Levenberg-Marquardt damping (0 = pure Fisher scoring)
        int it = 0;
        bool conv = false;
        std::vector<char> act = fr;
        for (it = 0; it < maxIter; ++it) {
            std::fill(grad.begin(), grad.end(), 0.);
            std::fill(H.begin(), H.end(), 0.);
            std::fill(hdiag.begin(), hdiag.end(), 0.);
            for (int j = 0; j < I.G.ny; ++j)
                for (int i = 0; i < I.G.nx; ++i) {
                    const int b = I.G.idx(i, j);
                    if (!Used(b)) continue;
                    const int nt = Touched(i, j, idx);
                    for (int q = 0; q < nt; ++q) d[idx[q]] = 0.;
                    const double mu = Mu(p, i, j, d.data());
                    const double mm = std::max(mu, 1e-9);
                    const double gfac = 2. * (1. - n[b] / mm);
                    const double w = 2. / mm;
                    for (int q = 0; q < nt; ++q) {
                        const int a = idx[q];
                        if (!act[a] || d[a] == 0.) continue;
                        grad[a] += gfac * d[a];
                        for (int s = 0; s < nt; ++s) {
                            const int c = idx[s];
                            if (!act[c] || d[c] == 0.) continue;
                            H[(size_t)a * np + c] += w * d[a] * d[c];
                        }
                    }
                }
            Penalty(p, g, &grad, &hdiag);
            for (int a = 0; a < np; ++a) H[(size_t)a * np + a] += hdiag[a];
            // freeze parameters without information
            for (int a = 0; a < np; ++a) if (act[a] && !(H[(size_t)a * np + a] > 0)) act[a] = 0;
            // solve (H + lambda diag H) step = -grad on the active set (Cholesky)
            std::vector<int> ai;
            for (int a = 0; a < np; ++a) if (act[a]) ai.push_back(a);
            const int m = (int)ai.size();
            std::vector<double> A((size_t)m * m), bvec(m);
            for (int u = 0; u < m; ++u) {
                bvec[u] = -grad[ai[u]];
                for (int v = 0; v < m; ++v) A[(size_t)u * m + v] = H[(size_t)ai[u] * np + ai[v]];
                A[(size_t)u * m + u] *= (1. + lambda);
            }
            if (!Chol(A, bvec, m)) { lambda = lambda > 0 ? lambda * 10 : 1e-3; if (lambda > 1e6) { R.msg = "singular"; break; } continue; }
            std::fill(step.begin(), step.end(), 0.);
            for (int u = 0; u < m; ++u) step[ai[u]] = bvec[u];
            // line search
            double t = 1.;
            std::vector<double> pn(np);
            double Fn = F;
            bool okStep = false;
            for (int ls = 0; ls < 40; ++ls) {
                for (int a = 0; a < np; ++a) pn[a] = p[a] + t * step[a];
                Fn = NLL(pn, n, g);
                if (std::isfinite(Fn) && Fn <= F + 1e-10) { okStep = true; break; }
                t *= 0.5;
            }
            if (!okStep) {
                if (lambda < 1e6) { lambda = lambda > 0 ? lambda * 10 : 1e-3; continue; }
                R.msg = "line search failed";
                break;
            }
            const double dF = F - Fn;
            p = pn; F = Fn;
            lambda = lambda > 1e-8 ? lambda * 0.1 : 0.;
            // expected decrease (Newton decrement)
            double dec = 0;
            for (int a = 0; a < np; ++a) dec += -grad[a] * step[a];
            if (dF < 1e-7 && dec < 1e-6) { conv = true; break; }
        }
        R.iter = it;
        R.ok = conv;
        if (!conv && R.msg.empty()) R.msg = "max iterations";
        R.F = F;
        R.p = p;
        // errors from the inverse Fisher matrix at the optimum (active parameters)
        R.err.assign(np, 0.);
        {
            std::fill(H.begin(), H.end(), 0.);
            std::fill(hdiag.begin(), hdiag.end(), 0.);
            for (int j = 0; j < I.G.ny; ++j)
                for (int i = 0; i < I.G.nx; ++i) {
                    const int b = I.G.idx(i, j);
                    if (!Used(b)) continue;
                    const int nt = Touched(i, j, idx);
                    for (int q = 0; q < nt; ++q) d[idx[q]] = 0.;
                    const double mu = std::max(Mu(p, i, j, d.data()), 1e-9);
                    for (int q = 0; q < nt; ++q) {
                        const int a = idx[q];
                        if (!act[a] || d[a] == 0.) continue;
                        for (int s = 0; s < nt; ++s) {
                            const int c = idx[s];
                            if (!act[c] || d[c] == 0.) continue;
                            H[(size_t)a * np + c] += 2. / mu * d[a] * d[c];
                        }
                    }
                }
            Penalty(p, g, nullptr, &hdiag);
            for (int a = 0; a < np; ++a) H[(size_t)a * np + a] += hdiag[a];
            std::vector<int> ai;
            for (int a = 0; a < np; ++a) if (act[a]) ai.push_back(a);
            const int m = (int)ai.size();
            std::vector<double> A((size_t)m * m);
            for (int u = 0; u < m; ++u) for (int v = 0; v < m; ++v) A[(size_t)u * m + v] = 0.5 * H[(size_t)ai[u] * np + ai[v]];
            std::vector<double> inv;
            if (Inverse(A, m, inv)) for (int u = 0; u < m; ++u) R.err[ai[u]] = std::sqrt(std::max(0., inv[(size_t)u * m + u]));
            R.rErr = R.err[L.iR];
        }
        return R;
    }

    // ------------------------------------------------------------------ linear algebra
    static bool Chol(std::vector<double>& A, std::vector<double>& b, int n) {   // solves A x = b in place (b <- x)
        for (int j = 0; j < n; ++j) {
            double s = A[(size_t)j * n + j];
            for (int k = 0; k < j; ++k) s -= A[(size_t)j * n + k] * A[(size_t)j * n + k];
            if (!(s > 0)) return false;
            const double l = std::sqrt(s);
            A[(size_t)j * n + j] = l;
            for (int i = j + 1; i < n; ++i) {
                double t = A[(size_t)i * n + j];
                for (int k = 0; k < j; ++k) t -= A[(size_t)i * n + k] * A[(size_t)j * n + k];
                A[(size_t)i * n + j] = t / l;
            }
        }
        for (int i = 0; i < n; ++i) { double t = b[i]; for (int k = 0; k < i; ++k) t -= A[(size_t)i * n + k] * b[k]; b[i] = t / A[(size_t)i * n + i]; }
        for (int i = n - 1; i >= 0; --i) { double t = b[i]; for (int k = i + 1; k < n; ++k) t -= A[(size_t)k * n + i] * b[k]; b[i] = t / A[(size_t)i * n + i]; }
        return true;
    }
    static bool Inverse(const std::vector<double>& A0, int n, std::vector<double>& inv) {
        inv.assign((size_t)n * n, 0.);
        for (int c = 0; c < n; ++c) {
            std::vector<double> A = A0, b(n, 0.);
            b[c] = 1.;
            if (!Chol(A, b, n)) return false;
            for (int r = 0; r < n; ++r) inv[(size_t)r * n + c] = b[r];
        }
        return true;
    }
};

// ============================================================================ test statistic and limits
struct Profile {
    FitResult hat;        // free fit (r may be negative)
    FitResult zero;       // r = 0
    double rhat = 0., rhatErr = 0.;
    double Fref = 0.;     // reference: F at max(r_hat, 0)
};
inline Profile ProfileFits(const Model& M, const std::vector<double>& n, const Glob& g) {
    Profile P;
    P.zero = M.Fit(n, 0., false, g);
    P.hat = M.Fit(n, 0., true, g, &P.zero.p);
    if (!P.hat.ok && P.zero.ok) {   // retry from the default start
        FitResult h2 = M.Fit(n, 0., true, g);
        if (h2.ok && h2.F < P.hat.F) P.hat = h2;
    }
    P.rhat = P.hat.p.empty() ? 0. : P.hat.p[M.L.iR];
    P.rhatErr = P.hat.rErr;
    P.Fref = (P.rhat >= 0. && P.hat.F <= P.zero.F) ? P.hat.F : P.zero.F;
    if (P.rhat < 0.) P.Fref = P.zero.F;
    return P;
}
// q~_mu for one data set (the profile P of that data set)
inline double QTilde(const Model& M, const std::vector<double>& n, const Glob& g, const Profile& P, double mu, FitResult* fmu = nullptr) {
    if (P.rhat > mu) return 0.;
    const std::vector<double>* st = (P.rhat >= 0 && !P.hat.p.empty()) ? &P.hat.p : &P.zero.p;
    FitResult f = M.Fit(n, mu, false, g, st);
    if (!f.ok) { FitResult f2 = M.Fit(n, mu, false, g); if (f2.F < f.F) f = f2; }
    if (fmu) *fmu = f;
    return std::max(0., f.F - P.Fref) / (M.cfg.qScale > 0 ? M.cfg.qScale : 1.);   // x17c: variance factor
}
inline double Phi(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.)); }
inline double PhiInv(double p) { return TMath::NormQuantile(p); }
// asymptotic CLs for q~ (Cowan et al. 2011, eqs. 65-66), qA = q_mu of the background-only Asimov set
inline double CLsAsym(double q, double qA, double* clsb = nullptr, double* clb = nullptr) {
    if (!(qA > 0)) { if (clsb) *clsb = 1; if (clb) *clb = 1; return 1.; }
    const double sq = std::sqrt(std::max(q, 0.)), sqA = std::sqrt(qA);
    double ps, pb;
    if (q <= qA) { ps = 1 - Phi(sq); pb = Phi(sqA - sq); }
    else { ps = 1 - Phi((q + qA) / (2 * sqA)); pb = 1 - Phi((q - qA) / (2 * sqA)); }
    if (clsb) *clsb = ps;
    if (clb) *clb = pb;
    return pb > 0 ? ps / pb : 1.;
}
struct Limit {
    bool ok = false;
    double rhat = 0, rhatErr = 0;
    double ul = 0;                 // observed (asymptotic CLs)
    double exp[5] = {0, 0, 0, 0, 0};   // expected -2, -1, median, +1, +2 sigma (Asimov)
    double sigmaA = 0;             // sigma_A at the observed limit
    double ulPL = 0;               // one-sided profile-likelihood limit without CLs (q~ = 2.706)
    Profile prof;
    std::vector<double> asimov;
    Glob glob;
};
// Asimov background-only data set: the model at the r = 0 fit of the data
inline std::vector<double> AsimovB(const Model& M, const Profile& P) { return M.Expected(P.zero.p); }
// root finder on a monotone function f(mu) - target in [lo, hi]
inline double Solve(const std::function<double(double)>& f, double target, double lo, double hi, int nit = 60, double tol = 1e-4) {
    double flo = f(lo) - target, fhi = f(hi) - target;
    int grow = 0;
    while (flo * fhi > 0 && grow < 30) { hi *= 2; fhi = f(hi) - target; ++grow; }
    if (flo * fhi > 0) return std::numeric_limits<double>::quiet_NaN();
    for (int it = 0; it < nit; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double fm = f(mid) - target;
        if (fm * flo > 0) { lo = mid; flo = fm; } else { hi = mid; fhi = fm; }
        if (hi - lo < tol * std::max(1e-12, std::fabs(mid))) break;
    }
    return 0.5 * (lo + hi);
}
// observed and expected asymptotic CLs limits on r; scale = a typical size of r (start of the bracket)
inline Limit LimitAsym(const Model& M, const std::vector<double>& n, const Glob& g, double scale, bool doExp = true) {
    Limit Lr;
    Lr.glob = g;
    Lr.prof = ProfileFits(M, n, g);
    const Profile& P = Lr.prof;
    Lr.rhat = P.rhat; Lr.rhatErr = P.rhatErr;
    Lr.asimov = AsimovB(M, P);
    Profile PA;   // Asimov profile: r_hat = 0 exactly; reference = r = 0 fit
    PA.zero = M.Fit(Lr.asimov, 0., false, g, &P.zero.p);
    PA.hat = PA.zero; PA.rhat = 0; PA.Fref = PA.zero.F;
    auto qA = [&](double mu) { return QTilde(M, Lr.asimov, g, PA, mu); };
    auto qObs = [&](double mu) { return QTilde(M, n, g, P, mu); };
    auto cls = [&](double mu) {
        if (mu <= 0) return 1.;
        return CLsAsym(qObs(mu), qA(mu));
    };
    const double lo = std::max(1e-9, 1e-4 * scale);
    double hi = std::max(scale, P.rhat > 0 ? P.rhat + 3 * std::max(P.rhatErr, 1e-6) : scale);
    // CLs decreases with mu: solve CLs(mu) = alpha
    Lr.ul = Solve([&](double mu) { return -cls(mu); }, -kAlpha, lo, hi);
    {
        const double qa = qA(Lr.ul);
        Lr.sigmaA = qa > 0 ? Lr.ul / std::sqrt(qa) : 0.;
    }
    Lr.ulPL = Solve([&](double mu) { return qObs(mu); }, 2.70554, std::max(lo, P.rhat > 0 ? P.rhat : lo), hi);
    if (doExp) {
        for (int k = 0; k < 5; ++k) {
            const double N = k - 2.;
            const double z = PhiInv(1 - kAlpha * Phi(N)) + N;
            // mu_N = sigma(mu_N) z, sigma(mu) = mu / sqrt(qA(mu)) -> solve sqrt(qA(mu)) = z
            if (z <= 0) { Lr.exp[k] = 0; continue; }
            Lr.exp[k] = Solve([&](double mu) { return std::sqrt(qA(mu)); }, z, lo, hi);
        }
    }
    Lr.ok = P.zero.ok && std::isfinite(Lr.ul);
    return Lr;
}

// ------------------------------------------------------------------ rows without data (blinding, X-c)
// After a fit with Config::use masking whole M(ee) rows (e.g. 12-28 MeV/c2 before x17_frozen.txt), the n_j of those rows
// carry no information and stay at their start values. Before an Asimov set is built from such a fit (X-c), they must be
// predicted: n_j = sim_j x exp(pol(y_j)), with sim_j a smooth reference per row (e.g. the simulated G_D + G_C row sums, or
// the like-sign rows) and the polynomial (order 'order' in y = M(ee)) fitted by least squares to ln(n_j/sim_j) of the
// rows that have data (weights: the relative Hessian errors of n_j if given, else equal). Returns the number of rows filled.
inline int FillMaskedRows(const Model& M, std::vector<double>& p, const std::vector<double>& simRow, int order = 2,
                          const std::vector<double>* nErr = nullptr) {
    const Inputs& I = *M.in;
    std::vector<char> rowUsed(I.G.ny, 0);
    for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) if (M.Used(I.G.idx(i, j))) rowUsed[j] = 1;
    const int nc = order + 1;
    std::vector<double> A((size_t)nc * nc, 0.), b(nc, 0.);
    int nfit = 0;
    for (int j = 0; j < I.G.ny; ++j) {
        if (!rowUsed[j] || !(p[j] > 0) || !(simRow[j] > 0)) continue;
        const double y = (I.G.yc(j) - 0.030) / 0.030, z = std::log(p[j] / simRow[j]);
        double w = 1.;
        if (nErr && (*nErr)[j] > 0) { const double e = (*nErr)[j] / p[j]; w = 1. / (e * e); }
        std::vector<double> xp(nc, 1.);
        for (int k = 1; k < nc; ++k) xp[k] = xp[k - 1] * y;
        for (int u = 0; u < nc; ++u) { b[u] += w * xp[u] * z; for (int v = 0; v < nc; ++v) A[(size_t)u * nc + v] += w * xp[u] * xp[v]; }
        ++nfit;
    }
    if (nfit < nc || !Model::Chol(A, b, nc)) return 0;
    int nf = 0;
    for (int j = 0; j < I.G.ny; ++j) {
        if (rowUsed[j] || !(simRow[j] > 0)) continue;
        const double y = (I.G.yc(j) - 0.030) / 0.030;
        double z = 0, xp = 1;
        for (int k = 0; k < nc; ++k) { z += b[k] * xp; xp *= y; }
        p[j] = simRow[j] * std::exp(z);
        ++nf;
    }
    return nf;
}

// ------------------------------------------------------------------ toy-based CLs
struct ToyCLsPoint { double mu = 0, qobs = 0, clsb = 0, clb = 0, cls = 0, clsErr = 0; int nt = 0; };
// Poisson pseudo-data and randomised global observables around the conditional estimate p (constrained parameters)
inline void ToyData(const Model& M, const std::vector<double>& p, std::mt19937_64& rng, std::vector<double>& n, Glob& g, const Glob& gobs) {
    const std::vector<double> mu = M.Expected(p);
    n.resize(mu.size());
    const double S = M.cfg.qScale > 1. ? M.cfg.qScale : 1.;   // x17c: over-dispersed toys S x Poisson(mu / S)
    for (size_t b = 0; b < mu.size(); ++b) {
        if (mu[b] <= 0) { n[b] = 0; continue; }
        std::poisson_distribution<long long> P(mu[b] / S);
        n[b] = S * (double)P(rng);
    }
    std::normal_distribution<double> N01(0., 1.);
    const Config& c = M.cfg;
    const Layout& L = M.L;
    g = gobs;
    for (int e = 0; e < 2; ++e) if (c.kSig[e] > 0) g.k0[e] = p[L.iK + e] + c.kSig[e] * N01(rng);
    if (c.aSig > 0) g.a0 = p[L.iA] + c.aSig * N01(rng);
    if (c.gSig > 0) g.g0 = p[L.iG] + c.gSig * N01(rng);
    if (c.cTailSig > 0) g.c0 = p[L.iTC] + N01(rng);
    if (c.xSig > 0) g.x0 = p[L.iTX] + N01(rng);
    if (c.wSig > 0 && !M.in->dXw.empty()) g.w0 = p[L.iTW] + N01(rng);   // x17c
    if (c.sSig > 0 && !M.in->dXs.empty()) g.s0 = p[L.iTS] + N01(rng);   // x17c
}
inline ToyCLsPoint ToyCLs(const Model& M, const std::vector<double>& nobs, const Glob& g, const Profile& P, double mu, int ntoys,
                          unsigned long long seed) {
    ToyCLsPoint T;
    T.mu = mu;
    FitResult fmu;
    T.qobs = QTilde(M, nobs, g, P, mu, &fmu);
    if (fmu.p.empty()) fmu = M.Fit(nobs, mu, false, g);
    const std::vector<double> pmu = fmu.p, p0 = P.zero.p;
    std::mt19937_64 rng(seed);
    int ns = 0, nb = 0;
    std::vector<double> n;
    Glob gt;
    for (int t = 0; t < ntoys; ++t) {
        for (int h = 0; h < 2; ++h) {
            ToyData(M, h == 0 ? pmu : p0, rng, n, gt, g);
            const Profile Pt = ProfileFits(M, n, gt);
            const double q = QTilde(M, n, gt, Pt, mu);
            if (q >= T.qobs - 1e-9) (h == 0 ? ns : nb)++;
        }
    }
    T.nt = ntoys;
    T.clsb = (double)ns / ntoys;
    T.clb = (double)nb / ntoys;
    T.cls = T.clb > 0 ? T.clsb / T.clb : 1.;
    // binomial errors propagated
    const double es = std::sqrt(std::max(T.clsb * (1 - T.clsb), 1. / ntoys) / ntoys), eb = std::sqrt(std::max(T.clb * (1 - T.clb), 1. / ntoys) / ntoys);
    T.clsErr = T.clb > 0 ? T.cls * std::sqrt(std::pow(es / std::max(T.clsb, 1e-9), 2) + std::pow(eb / T.clb, 2)) : 0.;
    return T;
}
// toy UL: CLs scanned on a grid of mu; weighted least-squares fit of ln CLs = a + b mu + c mu^2 to the points with
// 0.005 < CLs < 0.6 (weights from the binomial errors), solved for CLs = 0.10; fallback: first crossing, log-linear
inline double ToyUL(const std::vector<ToyCLsPoint>& pts) {
    double S[5] = {0, 0, 0, 0, 0}, T[3] = {0, 0, 0};
    int n = 0;
    for (const auto& p : pts) {
        if (!(p.cls > 0.005 && p.cls < 0.6)) continue;
        const double e = p.clsErr > 0 ? p.clsErr / p.cls : 0.1;
        const double w = 1. / (e * e), y = std::log(p.cls), x = p.mu;
        double xp = 1;
        for (int k = 0; k < 5; ++k) { S[k] += w * xp; xp *= x; }
        T[0] += w * y; T[1] += w * y * x; T[2] += w * y * x * x;
        ++n;
    }
    if (n >= 4) {
        std::vector<double> A = {S[0], S[1], S[2], S[1], S[2], S[3], S[2], S[3], S[4]}, b = {T[0], T[1], T[2]};
        if (Model::Chol(A, b, 3)) {
            const double a = b[0] - std::log(kAlpha), bb = b[1], c = b[2];
            double mu = std::numeric_limits<double>::quiet_NaN();
            if (std::fabs(c) < 1e-12 * std::fabs(bb)) mu = -a / bb;
            else {
                const double D = bb * bb - 4 * c * a;
                if (D >= 0) {
                    const double r1 = (-bb - std::sqrt(D)) / (2 * c), r2 = (-bb + std::sqrt(D)) / (2 * c);
                    const double lo = pts.front().mu, hi = pts.back().mu;
                    for (double r : {r1, r2}) if (r >= lo && r <= hi) { mu = r; break; }
                }
            }
            if (std::isfinite(mu)) return mu;
        }
    }
    for (size_t k = 1; k < pts.size(); ++k) {
        const double a = pts[k - 1].cls, b = pts[k].cls;
        if (a >= kAlpha && b < kAlpha) {
            const double la = std::log(std::max(a, 1e-6)), lb = std::log(std::max(b, 1e-6)), lt = std::log(kAlpha);
            return pts[k - 1].mu + (pts[k].mu - pts[k - 1].mu) * (la - lt) / (la - lb);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// ============================================================================ helpers for the builders
// exp(pol3) smoothing of a 1D shape (content with errors); returns the smoothed vector (same normalisation) on the bins
inline std::vector<double> SmoothExpPol3(const std::vector<double>& y, const std::vector<double>& ey, const std::vector<double>& x, int order = 3) {
    const int n = (int)y.size();
    TH1D h("x17c_sm", "", n, 0, n);
    h.SetDirectory(nullptr);
    double s = 0;
    for (int i = 0; i < n; ++i) { h.SetBinContent(i + 1, y[i]); h.SetBinError(i + 1, ey[i] > 0 ? ey[i] : 1e-9 + 0.1 * std::fabs(y[i])); s += y[i]; }
    std::vector<double> out(n, 0.);
    if (!(s > 0)) return out;
    // fit in the variable x (tilt variable) through an explicit TF1 on the bin index -> x map
    const double x0 = x.front(), dx = (x.back() - x.front()) / std::max(1, n - 1);
    TF1 f("x17c_ep3", Form("exp([0]+[1]*(%g+%g*(x-0.5))+[2]*pow(%g+%g*(x-0.5),2)+[3]*pow(%g+%g*(x-0.5),3))", x0, dx, x0, dx, x0, dx), 0, n);
    f.SetParameters(std::log(s / n), 0, 0, 0);
    if (order < 3) f.FixParameter(3, 0.);
    if (order < 2) f.FixParameter(2, 0.);
    h.Fit(&f, "QN0");
    double t = 0;
    for (int i = 0; i < n; ++i) { out[i] = f.Integral(i, i + 1); t += out[i]; }
    if (t > 0) for (auto& v : out) v *= s / t;
    return out;
}
inline void Normalise(std::vector<double>& v) { double s = 0; for (double x : v) s += x; if (s > 0) for (auto& x : v) x /= s; }
// factorise and smooth a 2D template T (x fastest): T_ij -> (row sum_j) x s_c(i), s_c = exp(pol3)-smoothed M4 shape of the
// class c(j) (sum of the rows of the class, errors from err2). Removes the single-cell spikes of heavy simulated rows
// (one std-set row = 27 expected data rows) while keeping the M(ee) distribution and the per-class M4 shape.
inline void SmoothFactorise(std::vector<double>& T, const std::vector<double>& err2, const Inputs& I, int order = 3) {
    const int nx = I.G.nx, ny = I.G.ny;
    std::vector<std::vector<double>> s(I.nbcl, std::vector<double>(nx, 0.)), e(I.nbcl, std::vector<double>(nx, 0.));
    std::vector<double> row(ny, 0.);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int b = I.G.idx(i, j);
            s[I.bcl[j]][i] += T[b];
            e[I.bcl[j]][i] += err2.empty() ? std::fabs(T[b]) : err2[b];
            row[j] += T[b];
        }
    for (int c = 0; c < I.nbcl; ++c) {
        double tot = 0; for (double v : s[c]) tot += v;
        if (!(tot > 0)) continue;
        for (auto& v : e[c]) v = std::sqrt(v);
        std::vector<double> sm = SmoothExpPol3(s[c], e[c], I.tx, order);
        double t2 = 0; for (double v : sm) t2 += v;
        if (!(t2 > 0)) sm = s[c], t2 = tot;
        for (auto& v : sm) v /= t2;
        s[c] = sm;
    }
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) T[I.G.idx(i, j)] = row[j] * s[I.bcl[j]][i];
}

// a line (fine histogram in MeV, any binning) mapped to the M(ee) grid [GeV] with x' = shift + scale x and a resolution
// change around the line mean: x' = m0 + (1 + res) (x - m0); returns the fraction of the histogram TOTAL (incl. flows)
// in every grid bin; tot = normalisation (sum of all fine bins incl. flows) if > 0
inline std::vector<double> LineToGrid(const TH1D* h, const Grid& G, double scale = 1., double res = 0., double tot = -1.) {
    std::vector<double> X(G.ny, 0.);
    if (!h) return X;
    double T = tot > 0 ? tot : h->Integral(0, h->GetNbinsX() + 1);
    if (!(T > 0)) return X;
    const double m0 = h->GetMean();
    for (int k = 1; k <= h->GetNbinsX(); ++k) {
        const double c = h->GetBinContent(k);
        if (c == 0) continue;
        double a = h->GetXaxis()->GetBinLowEdge(k), b = h->GetXaxis()->GetBinUpEdge(k);
        auto tr = [&](double x) { return scale * (m0 + (1 + res) * (x - m0)) * 1e-3; };   // MeV -> GeV
        double ta = tr(a), tb = tr(b);
        if (tb < ta) std::swap(ta, tb);
        const double w = tb - ta;
        for (int j = 0; j < G.ny; ++j) {
            const double lo = std::max(ta, G.ye[j]), hi = std::min(tb, G.ye[j + 1]);
            if (hi > lo) X[j] += c * (hi - lo) / w;
        }
    }
    for (auto& x : X) x /= T;
    return X;
}

// ============================================================================ builder: hists_v1.root (data scale, X2D)
// CB options: 0 = WP0 spread rule (CB_X2D, nominal); 1 = factorised over the X2D M4 range: CB_fine(i) x f(j) with f = the
// like-sign share (N++ + N--) of M(ee) bin j summed over the X2D M4 range; 2 = factorised geometric mean: the row total
// 2 sqrt(N++_j N--_j) (summed over the X2D M4 range) spread in M4 like CB_fine(i)
struct HistOpt {
    std::string sel = "S0", var = "nom";
    int cbMode = 0;
    int cbRep = -1;          // >= 0: use CB replica (CBrep_X2D) number cbRep instead of the nominal CB (mode 0 only)
    int nbcl = 2;            // 2: E1/E2; 3: E1, 10-30, >= 30 MeV
    bool smoothB = true;     // exp(pol3) smoothing of the background shape per class
    bool smoothEta = true;   // factorise and smooth C2+C3 and G_eta (SmoothFactorise)
    int smoothOrder = 2;     // order of the exp(pol) smoothing of the background shape B (2 = exp(pol2), NOMINAL)
    double kappaB = 1.;      // G_D + kappa_B G_C
    double bTilt = 0.;       // fixed tilt multiplied into B (from WP5), b per 0.1 GeV/c2
    double aEta = 1.;        // scale of the eta templates used for the N5/N4 start values
    std::string simSuffix = "";   // "", "_hA", "_hB" (halves: x 2 applied)
    int simRep = -1;         // >= 0: use simulation replica number simRep (TH3D _rep) for C4 (conversion) only
};
inline TH1* GetH(TFile* f, const std::string& path) {
    TH1* h = f ? (TH1*)f->Get(path.c_str()) : nullptr;
    if (!h) { printf("x17c: %s not found\n", path.c_str()); return nullptr; }
    return h;
}
// fill the per-class vectors of an Inputs from 2D arrays (x fastest) of D, C, X totals etc. (shared by the builders)
inline void ClassesByMee(Inputs& I, int nbcl) {
    I.nbcl = nbcl;
    I.kcl.assign(I.G.ny, 0); I.bcl.assign(I.G.ny, 0);
    for (int j = 0; j < I.G.ny; ++j) {
        const double y = I.G.yc(j);
        I.kcl[j] = y < 0.010 ? 0 : 1;
        if (nbcl == 2) I.bcl[j] = y < 0.010 ? 0 : 1;
        else I.bcl[j] = y < 0.010 ? 0 : (y < 0.030 ? 1 : 2);
    }
    I.Ctail.assign(I.G.ny, 0.);
    for (int j = 0; j < I.G.ny; ++j) I.Ctail[j] = I.G.ye[j] >= kTailMee - 1e-9 ? 1. : 0.;
}
// sum over the rows j of class c of a 2D array (nx x (ny+2) with flows at j = -1 -> index 0 ...)
inline bool BuildFromHists(TFile* f, Inputs& I, const HistOpt& o, double nbToRaw) {
    const std::string dd = "data/" + o.sel + "/all/", sd = "sim/" + o.sel + "__" + o.var + "/";
    TH2* hcb = (TH2*)GetH(f, dd + "CB_X2D");
    TH2* hpp = (TH2*)GetH(f, dd + "Npp_X2D");
    TH2* hmm = (TH2*)GetH(f, dd + "Nmm_X2D");
    const char* comps[8] = {"C2", "C3", "C4", "C5", "C6", "Geta", "GD", "Gconv"};
    TH2* hs[8];
    for (int c = 0; c < 8; ++c) { hs[c] = (TH2*)GetH(f, sd + comps[c] + "_X2D" + o.simSuffix); if (!hs[c]) return false; }
    if (!hcb || !hpp || !hmm) return false;
    const double sfac = nbToRaw * (o.simSuffix.empty() ? 1. : 2.);
    Grid& G = I.G;
    G.nx = hcb->GetNbinsX(); G.ny = hcb->GetNbinsY();
    G.xe.clear(); G.ye.clear();
    for (int i = 1; i <= G.nx + 1; ++i) G.xe.push_back(hcb->GetXaxis()->GetBinLowEdge(i));
    for (int j = 1; j <= G.ny + 1; ++j) G.ye.push_back(hcb->GetYaxis()->GetBinLowEdge(j));
    const int nx = G.nx, ny = G.ny;
    I.tx.resize(nx);
    for (int i = 0; i < nx; ++i) I.tx[i] = (G.xc(i) - 0.55) / 0.1;
    ClassesByMee(I, o.nbcl);
    // CB
    I.CB.assign(G.n(), 0.);
    TH3* hrep = nullptr;
    if (o.cbRep >= 0) hrep = (TH3*)GetH(f, dd + "CBrep_X2D");
    if (o.cbMode == 0) {
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
            I.CB[G.idx(i, j)] = hrep ? hrep->GetBinContent(i + 1, j + 1, o.cbRep + 1) : hcb->GetBinContent(i + 1, j + 1);
    } else {
        std::vector<double> fine(nx, 0.), ls(ny + 2, 0.), pp(ny + 2, 0.), mm(ny + 2, 0.);
        for (int i = 0; i < nx; ++i) for (int j = 0; j <= ny + 1; ++j) fine[i] += hcb->GetBinContent(i + 1, j);
        double lt = 0;
        for (int j = 0; j <= ny + 1; ++j) for (int i = 0; i < nx; ++i) {
            pp[j] += hpp->GetBinContent(i + 1, j); mm[j] += hmm->GetBinContent(i + 1, j);
        }
        for (int j = 0; j <= ny + 1; ++j) { ls[j] = pp[j] + mm[j]; lt += ls[j]; }
        double ftot = 0; for (double x : fine) ftot += x;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            double v;
            if (o.cbMode == 1) v = lt > 0 ? fine[i] * ls[j + 1] / lt : 0.;
            else v = ftot > 0 ? 2. * std::sqrt(pp[j + 1] * mm[j + 1]) * fine[i] / ftot : 0.;
            I.CB[G.idx(i, j)] = v;
        }
    }
    // eta-related terms
    I.E23.assign(G.n(), 0.); I.Geta.assign(G.n(), 0.);
    std::vector<double> e23(G.n(), 0.), egt(G.n(), 0.);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const int b = G.idx(i, j);
        I.E23[b] = sfac * (hs[0]->GetBinContent(i + 1, j + 1) + hs[1]->GetBinContent(i + 1, j + 1));
        I.Geta[b] = sfac * hs[5]->GetBinContent(i + 1, j + 1);
        e23[b] = sfac * sfac * (std::pow(hs[0]->GetBinError(i + 1, j + 1), 2) + std::pow(hs[1]->GetBinError(i + 1, j + 1), 2));
        egt[b] = sfac * sfac * std::pow(hs[5]->GetBinError(i + 1, j + 1), 2);
    }
    if (o.smoothEta) { SmoothFactorise(I.E23, e23, I); SmoothFactorise(I.Geta, egt, I); }
    // conversion template (optionally a replica)
    TH2* hC4 = hs[2];
    TH2D* c4rep = nullptr;
    if (o.simRep >= 0) {
        TH3* h3 = (TH3*)GetH(f, sd + "C4_X2D_rep");
        if (h3) {
            c4rep = new TH2D("x17c_c4rep", "", nx, G.xe.data(), ny, G.ye.data());
            c4rep->SetDirectory(nullptr);
            for (int i = 0; i <= nx + 1; ++i) for (int j = 0; j <= ny + 1; ++j) c4rep->SetBinContent(i, j, h3->GetBinContent(i, j, o.simRep + 1));
            hC4 = c4rep;
        }
    }
    // M(ee) fractions (totals include the M(ee) flows; M4 inside the range)
    auto frac = [&](TH2* h, std::vector<double>& v) {
        double tot = 0;
        v.assign(ny, 0.);
        for (int j = 0; j <= ny + 1; ++j) for (int i = 1; i <= nx; ++i) tot += h->GetBinContent(i, j);
        for (int j = 0; j < ny; ++j) for (int i = 1; i <= nx; ++i) v[j] += h->GetBinContent(i, j + 1);
        if (tot > 0) for (auto& x : v) x /= tot;
        return tot;
    };
    const double t5 = frac(hs[3], I.D), t4 = frac(hC4, I.C), t6 = frac(hs[4], I.X);
    I.RX = (t5 > 0) ? t6 * kBRdir / t5 : 1.;
    I.N5start = o.aEta * t5 * sfac; I.N4start = o.aEta * t4 * sfac;
    // peak M4 shapes per class from C5 (+ C4 in the same proportion as their sums in the class)
    I.gP.assign(I.nbcl, std::vector<double>(nx, 0.));
    I.B.assign(I.nbcl, std::vector<double>(nx, 0.));
    std::vector<std::vector<double>> be(I.nbcl, std::vector<double>(nx, 0.));
    for (int j = 0; j < ny; ++j) {
        const int c = I.bcl[j];
        for (int i = 0; i < nx; ++i) {
            I.gP[c][i] += hs[3]->GetBinContent(i + 1, j + 1);
            I.B[c][i] += hs[6]->GetBinContent(i + 1, j + 1) + o.kappaB * hs[7]->GetBinContent(i + 1, j + 1);
            const double e1 = hs[6]->GetBinError(i + 1, j + 1), e2 = o.kappaB * hs[7]->GetBinError(i + 1, j + 1);
            be[c][i] += e1 * e1 + e2 * e2;
        }
    }
    for (int c = 0; c < I.nbcl; ++c) {
        for (auto& e : be[c]) e = std::sqrt(e);
        if (o.smoothB) I.B[c] = SmoothExpPol3(I.B[c], be[c], I.tx, o.smoothOrder);
        if (o.bTilt != 0.) for (int i = 0; i < nx; ++i) I.B[c][i] *= (1 + o.bTilt * I.tx[i]);
        Normalise(I.B[c]);
        Normalise(I.gP[c]);
    }
    // X17 M4 shape (C6, all M(ee))
    I.gX.assign(nx, 0.);
    for (int i = 0; i < nx; ++i) for (int j = 0; j <= ny + 1; ++j) I.gX[i] += hs[4]->GetBinContent(i + 1, j);
    Normalise(I.gX);
    delete c4rep;
    I.label = o.sel + "__" + o.var;
    return true;
}

// ============================================================================ builder: toy grid (simulation blind test)
// templates.root of the simulation study (half A x 2, rows at L_sim, M4 0.4830-0.6080 in 2.5 MeV x M(ee) 0-100 MeV in
// 0.5 MeV, pre + both pions in the vertex fit, simulated M4 without mass transfer). The toy grid: M4 bins i0..i1 (1-based
// of the template axis), M(ee) bins merged by 'merge' from y0 (1-based) for ny bins.
struct ToyGridOpt {
    int i0 = 17, i1 = 50;      // M4 0.5230-0.6080 (simulation scale)
    int y0 = 5, merge = 4, ny = 30;   // M(ee) 2-62 MeV/c2 in 2 MeV/c2
    double scale = 6080.5 / 221.05;   // half A x 2 at L_sim -> toys at L_data (6080.5) of the simulation study
    int nbcl = 2;
    bool smoothB = true;
    bool smoothEta = true;     // factorise and smooth C2+C3 and G_eta
    int smoothOrder = 2;       // order of the exp(pol) smoothing of B (NOMINAL 2)
    int lineMode = 0;          // 0: X from the half-A C6 template; 1: from the slice line of x17_lines.root at 17.0
};
const char* const kToyTemplates = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/templates/templates.root";
const char* const kToyFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/pseudo/x17_toys.root";
const char* const kToyKey = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/sealed/x17_key.root";
// rebin a toy-format TH2 (50 x 200) onto the toy grid: out[j*nx+i]
inline std::vector<double> ToyRebin(const TH2* h, const ToyGridOpt& o, std::vector<double>* err2 = nullptr) {
    const int nx = o.i1 - o.i0 + 1;
    std::vector<double> v((size_t)nx * o.ny, 0.);
    if (err2) err2->assign(v.size(), 0.);
    for (int j = 0; j < o.ny; ++j)
        for (int k = 0; k < o.merge; ++k) {
            const int by = o.y0 + j * o.merge + k;
            for (int i = 0; i < nx; ++i) {
                v[(size_t)j * nx + i] += h->GetBinContent(o.i0 + i, by);
                if (err2) (*err2)[(size_t)j * nx + i] += std::pow(h->GetBinError(o.i0 + i, by), 2);
            }
        }
    return v;
}
inline bool BuildToyGrid(Inputs& I, const ToyGridOpt& o) {
    TFile* f = TFile::Open(kToyTemplates, "READ");
    if (!f || f->IsZombie()) return false;
    auto T = [&](const char* c) { return (TH2*)GetH(f, std::string("m4mee_X17_") + c); };
    TH2 *d5 = T("direct5"), *c4 = T("conv4"), *b23 = T("bump23"), *ge = T("Geta"), *gp = T("Gpi0D"), *gc = T("Gconv"), *gv = T("Gvec"),
        *go = T("Goth"), *x = T("x17perBR");
    if (!d5 || !c4 || !b23 || !ge || !gp || !gc || !gv || !go || !x) return false;
    Grid& G = I.G;
    G.nx = o.i1 - o.i0 + 1; G.ny = o.ny;
    G.xe.clear(); G.ye.clear();
    for (int i = 0; i <= G.nx; ++i) G.xe.push_back(d5->GetXaxis()->GetBinLowEdge(o.i0 + i));
    for (int j = 0; j <= G.ny; ++j) G.ye.push_back(1e-3 * d5->GetYaxis()->GetBinLowEdge(o.y0 + j * o.merge));
    const int nx = G.nx, ny = G.ny;
    I.tx.resize(nx);
    for (int i = 0; i < nx; ++i) I.tx[i] = (G.xc(i) - 0.55) / 0.1;
    ClassesByMee(I, o.nbcl);
    auto RB = [&](TH2* h, std::vector<double>* e2 = nullptr) {
        std::vector<double> v = ToyRebin(h, o, e2);
        for (auto& a : v) a *= o.scale;
        if (e2) for (auto& a : *e2) a *= o.scale * o.scale;
        return v;
    };
    std::vector<double> e23, egt;
    I.E23 = RB(b23, &e23);
    I.Geta = RB(ge, &egt);
    if (o.smoothEta) { SmoothFactorise(I.E23, e23, I); SmoothFactorise(I.Geta, egt, I); }
    I.CB.assign(G.n(), 0.);
    // totals over the M4 range incl. all M(ee) (flows of the 0-100 MeV axis included)
    auto total = [&](TH2* h) { double t = 0; for (int j = 0; j <= h->GetNbinsY() + 1; ++j) for (int i = o.i0; i <= o.i1; ++i) t += h->GetBinContent(i, j); return t; };
    auto fracs = [&](TH2* h, std::vector<double>& v) {
        const double t = total(h);
        std::vector<double> a = ToyRebin(h, o);
        v.assign(ny, 0.);
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) v[j] += a[(size_t)j * nx + i];
        if (t > 0) for (auto& q : v) q /= t;
        return t;
    };
    const double t5 = fracs(d5, I.D), t4 = fracs(c4, I.C), t6 = fracs(x, I.X);
    I.RX = t5 > 0 ? t6 * kBRdir / t5 : 1.;
    I.N5start = t5 * o.scale; I.N4start = t4 * o.scale;
    std::vector<double> vd5 = ToyRebin(d5, o), vx = ToyRebin(x, o);
    std::vector<double> e2p, e2c, e2v, e2o;
    std::vector<double> vp = ToyRebin(gp, o, &e2p), vc = ToyRebin(gc, o, &e2c), vv = ToyRebin(gv, o, &e2v), vo = ToyRebin(go, o, &e2o);
    I.gP.assign(I.nbcl, std::vector<double>(nx, 0.));
    I.B.assign(I.nbcl, std::vector<double>(nx, 0.));
    std::vector<std::vector<double>> be(I.nbcl, std::vector<double>(nx, 0.));
    I.gX.assign(nx, 0.);
    for (int j = 0; j < ny; ++j) {
        const int c = I.bcl[j];
        for (int i = 0; i < nx; ++i) {
            const size_t b = (size_t)j * nx + i;
            I.gP[c][i] += vd5[b];
            I.B[c][i] += vp[b] + vc[b] + vv[b] + vo[b];
            be[c][i] += e2p[b] + e2c[b] + e2v[b] + e2o[b];
            I.gX[i] += vx[b];
        }
    }
    for (int c = 0; c < I.nbcl; ++c) {
        for (auto& e : be[c]) e = std::sqrt(e);
        if (o.smoothB) I.B[c] = SmoothExpPol3(I.B[c], be[c], I.tx, o.smoothOrder);
        Normalise(I.B[c]);
        Normalise(I.gP[c]);
    }
    Normalise(I.gX);
    f->Close();
    delete f;
    I.label = "toy grid (half-A templates)";
    return true;
}

// ============================================================================ X17 line for an arbitrary mass
// reads x17_lines.root (x17ab_line.C). NOMINAL line at mass m: the simulated X17 line (C6, 17.0 MeV/c2, 0.69 M rows)
// morphed: x' = mean(m) + (rms(m)/rms_C6) (x - mean_C6), with the truncated mean and RMS of the reweighted C5 slices
// smoothed by pol2 in m (8-45 MeV/c2); at 17.0 the C6 line itself. R_X(m) = R_X(17, C6) x eps(m)/eps(17), eps = pol2 fit
// of the slice efficiencies (variant 0 nom, 1 DL, 2 C2, 3 NP). Variant rawSlice: the slice histogram itself (shifted to m,
// linear interpolation between the 0.5 MeV/c2 grid masses) and the raw slice R_X.
struct LineSet {
    bool ok = false;
    std::vector<double> m;            // MeV
    std::vector<TH1D*> h;             // fine slice lines (nominal variant)
    std::vector<double> RX[4];        // raw slice R_X
    TH1D* c6 = nullptr;
    double RX17c6 = 0, mean6 = 17, rms6 = 4.5;
    TF1 *fMean = nullptr, *fRms = nullptr, *fEps[4] = {nullptr, nullptr, nullptr, nullptr};
};
inline LineSet& Lines() { static LineSet s; return s; }
inline bool LoadLines(const char* file = kLinesFile, const char* var = "nom") {
    LineSet& S = Lines();
    if (S.ok) return true;
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) return false;
    TTree* t = (TTree*)f->Get("lines");
    TH1D* hc = (TH1D*)f->Get("consts");
    TH1D* hm = (TH1D*)f->Get("c6moments");
    if (!t || !hc || !hm) return false;
    double m, RX[4];
    t->SetBranchAddress("m", &m);
    t->SetBranchAddress("RX", RX);
    for (Long64_t k = 0; k < t->GetEntries(); ++k) {
        t->GetEntry(k);
        S.m.push_back(m);
        for (int v = 0; v < 4; ++v) S.RX[v].push_back(RX[v]);
        TH1D* h = (TH1D*)f->Get(Form("ee_m%03d_%s", (int)std::lround(10 * m), var));
        if (h) { h = (TH1D*)h->Clone(Form("x17c_line_%d", (int)k)); h->SetDirectory(nullptr); }
        S.h.push_back(h);
    }
    S.c6 = (TH1D*)f->Get("C6_ee");
    if (S.c6) { S.c6 = (TH1D*)S.c6->Clone("x17c_c6line"); S.c6->SetDirectory(nullptr); }
    S.RX17c6 = hc->GetBinContent(3);
    S.mean6 = hm->GetBinContent(1); S.rms6 = hm->GetBinContent(2);
    S.fMean = (TF1*)f->Get("fMean"); S.fRms = (TF1*)f->Get("fRms");
    const char* vn[4] = {"nom", "DL", "C2", "NP"};
    for (int v = 0; v < 4; ++v) S.fEps[v] = (TF1*)f->Get(Form("fEps_%s", vn[v]));
    f->Close(); delete f;
    S.ok = S.c6 && S.fMean && S.fRms && S.fEps[0];
    return S.ok;
}
// fractions of a fine line (MeV) in the M(ee) bins of G (GeV) after x' = scale (c1 + w (x - c0))
inline std::vector<double> LineMap(const TH1D* h, const Grid& G, double c0, double c1, double w, double scale, double tot = -1.) {
    std::vector<double> X(G.ny, 0.);
    const double T = tot > 0 ? tot : h->Integral(0, h->GetNbinsX() + 1);
    if (!(T > 0)) return X;
    for (int k = 1; k <= h->GetNbinsX(); ++k) {
        const double c = h->GetBinContent(k);
        if (c == 0) continue;
        double ta = 1e-3 * scale * (c1 + w * (h->GetXaxis()->GetBinLowEdge(k) - c0));
        double tb = 1e-3 * scale * (c1 + w * (h->GetXaxis()->GetBinUpEdge(k) - c0));
        if (tb < ta) std::swap(ta, tb);
        for (int j = 0; j < G.ny; ++j) {
            const double lo = std::max(ta, G.ye[j]), hi = std::min(tb, G.ye[j + 1]);
            if (hi > lo) X[j] += c * (hi - lo) / (tb - ta);
        }
    }
    for (auto& x : X) x /= T;
    return X;
}
// X_j and R_X for mass mX [MeV/c2]; rxVar 0 nom, 1 DL, 2 C2, 3 NP; scale = M(ee) scale (1 +- 0.02), res = relative change of
// the line width (+-0.10); rawSlice = use the slice histogram instead of the morphed C6 line
inline bool LineFor(double mX, const Grid& G, std::vector<double>& X, double& RX, int rxVar = 0, double scale = 1., double res = 0.,
                    bool rawSlice = false) {
    LineSet& S = Lines();
    if (!S.ok && !LoadLines()) return false;
    const double e17 = S.fEps[0]->Eval(17.);   // anchored through the NOMINAL slice efficiency at 17
    if (!rawSlice) {
        const double c1 = (std::fabs(mX - 17.) < 1e-9) ? S.mean6 : mX + S.fMean->Eval(mX) - (17. + S.fMean->Eval(17.)) + S.mean6;
        const double wr = (std::fabs(mX - 17.) < 1e-9) ? 1. : S.fRms->Eval(mX) / S.fRms->Eval(17.);
        X = LineMap(S.c6, G, S.mean6, c1, wr * (1 + res), scale);
        RX = S.RX17c6 * (S.fEps[rxVar] ? S.fEps[rxVar]->Eval(mX) / e17 : 1.);
        return true;
    }
    const double dm = S.m[1] - S.m[0];
    const double u = (mX - S.m[0]) / dm;
    const int k = (int)std::floor(u + 1e-9);
    if (k < 0 || k + 1 >= (int)S.m.size()) return false;
    const double w = u - k;
    auto one = [&](int q) {
        const double m0 = S.h[q]->GetMean();
        return LineMap(S.h[q], G, m0, m0 + (mX - S.m[q]), 1 + res, scale);
    };
    std::vector<double> Xa = one(k), Xb = w > 1e-9 ? one(k + 1) : Xa;
    X.assign(G.ny, 0.);
    for (int j = 0; j < G.ny; ++j) X[j] = (1 - w) * Xa[j] + w * Xb[j];
    RX = (1 - w) * S.RX[rxVar][k] + w * S.RX[rxVar][k + 1];
    return true;
}

}  // namespace x17c

#endif  // ALLREPORT_X17C_MODEL_H
