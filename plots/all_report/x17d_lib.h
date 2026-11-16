// ============================================================================
// x17d_lib.h  --  WP9 X-d (prefix x17d): X17 UNBLINDING helpers. Applies the FROZEN procedure of
//                 /hdd3/przygoda/HADES/FAT_out/all_report/x17/x17_frozen.txt (x17c_lib.h / x17c_model.h, unchanged and
//                 only included) to the full data X2D. Written 2026-11-14.
// ============================================================================
// What is new here (nothing of the frozen model, fitter, CLs or configuration is changed):
//   LoadFull   the full data X2D: side-band cells from hists_v1.root (frozen loader x17cl::LoadSideBands, including its
//              mask check) + the masked rows 12-28 MeV/c2 from the sealed file (x17_frozen.txt Sect. 10, step 1); the
//              side-band cells of the sealed file are compared with hists_v1.root (must be identical).
//   LocalP0    local p0 of the frozen model: q0 = (F(r = 0) - F(r_hat)) / qScale for r_hat > 0 (x17_frozen.txt Sect. 6).
//   TStat      two-sided profile-likelihood ratio with the physical boundary r >= 0 (Feldman-Cousins ordering):
//              t_r = (F(r) - F(max(r_hat, 0))) / qScale (the q~ of the frozen procedure without the one-sided zeroing).
//   FCToys     Neyman construction with t_r: toys from the conditional fit of the DATA at r (profile construction, as the
//              frozen ToyCLs), global observables drawn (x17c::ToyData, over-dispersed by qScale).
//   Comps      the model split into components (CB, other background, direct D, conversion C, eta-related, X17) for the
//              figures.
// ============================================================================
#ifndef ALLREPORT_X17D_LIB_H
#define ALLREPORT_X17D_LIB_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "x17c_lib.h"

namespace x17d {
using namespace x17cl;

const char* const kOutD = "/hdd3/przygoda/HADES/FAT_out/all_report/x17d";
const char* const kFigD = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTabD = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
// frozen numbers of x17_frozen.txt that are not in x17c_config.txt (Sect. 6 and 3)
const double kBpredX17c = 6311.7;     // background-only expectation of the frozen model in X17c (side-band fit, rows predicted)
const double kSigBX17c = 1195.5;      // frozen sigma_B(X17c)
const double kS1X17c = 2.0954e7;      // C6 per unit BR in X17c, raw rows at a_eta = 1
const double kN5ref = 7458.8;         // direct rows at SMASH (a_eta = 1) in the X2D M4 range
const double kFrozenExpMed17 = 2.759e-05;   // frozen Asimov median at 17 MeV/c2, normalisation (1) (reproduction check)

// ---------------------------------------------------------------- data: side bands (hists_v1) + masked rows (sealed)
struct Data {
    bool ok = false;
    std::vector<double> n;      // [nx*ny] all 1050 cells
    std::vector<char> useSB;    // side-band cells (as in X-c)
    std::vector<double> nSBonly;   // side-band rows only, masked rows 0 (as X-c used them)
    double nSB = 0, nMasked = 0, maxDiffSB = 0;
};
inline bool LoadFull(TFile* f, Data& D, const char* sel = "S0") {
    if (!LoadSideBands(f, D.n, D.useSB, sel)) return false;   // frozen loader (checks that the mask is intact)
    D.nSBonly = D.n;
    TFile* fs = TFile::Open(alldefs::kSealedFile, "READ");
    if (!fs || fs->IsZombie()) { printf("x17d: cannot open the sealed file\n"); return false; }
    TH2* h = (TH2*)fs->Get(Form("data/%s/all/N_X2D", sel));
    TH2* hm = (TH2*)f->Get(Form("data/%s/all/N_X2D", sel));
    if (!h || !hm) { printf("x17d: N_X2D missing in the sealed file\n"); return false; }
    const int nx = hm->GetNbinsX(), ny = hm->GetNbinsY();
    if (h->GetNbinsX() != nx || h->GetNbinsY() != ny) { printf("x17d: sealed X2D binning differs\n"); return false; }
    for (int i = 1; i <= nx + 1; ++i) if (std::fabs(h->GetXaxis()->GetBinLowEdge(i) - hm->GetXaxis()->GetBinLowEdge(i)) > 1e-12) { printf("x17d: x edges differ\n"); return false; }
    for (int j = 1; j <= ny + 1; ++j) if (std::fabs(h->GetYaxis()->GetBinLowEdge(j) - hm->GetYaxis()->GetBinLowEdge(j)) > 1e-12) { printf("x17d: y edges differ\n"); return false; }
    D.nSB = 0; D.nMasked = 0; D.maxDiffSB = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const size_t b = (size_t)j * nx + i;
            const double v = h->GetBinContent(i + 1, j + 1);
            if (D.useSB[b]) { D.maxDiffSB = std::max(D.maxDiffSB, std::fabs(v - D.n[b])); D.nSB += D.n[b]; }
            else { D.n[b] = v; D.nMasked += v; }
        }
    fs->Close();
    delete fs;
    D.ok = D.maxDiffSB == 0.;
    if (!D.ok) printf("x17d: WARNING sealed side-band cells differ from hists_v1.root (max %.3g)\n", D.maxDiffSB);
    return true;
}

// ---------------------------------------------------------------- context: frozen inputs + full data
struct Ctx {
    TFile* f = nullptr;
    Wp5 w;
    Frozen z;
    x17c::Inputs I;
    Data D;
};
inline bool Init(Ctx& c) {
    c.f = TFile::Open(alldefs::kHistFile, "READ");
    if (!c.f || c.f->IsZombie()) return false;
    c.w = ReadWp5();
    if (!c.w.ok || !ReadConfig(c.z)) return false;
    if (!Build(c.f, c.I, c.w, BuildOpt(), c.z.a0)) return false;   // frozen nominal build (x17c_freeze.C step 1)
    if (!SetLine(c.I, 17.0)) return false;
    return LoadFull(c.f, c.D);
}

// ---------------------------------------------------------------- local p0 (asymptotic)
struct P0 { bool ok = false; double rhat = 0, rerr = 0, q0 = 0, Z = 0, p0 = 0.5, F0 = 0, Fhat = 0; };
inline P0 LocalP0(const x17c::Model& M, const std::vector<double>& n, const x17c::Glob& g) {
    P0 r;
    const x17c::Profile P = x17c::ProfileFits(M, n, g);
    r.ok = P.zero.ok && P.hat.ok;
    r.rhat = P.rhat; r.rerr = P.rhatErr; r.F0 = P.zero.F; r.Fhat = P.hat.F;
    const double qs = M.cfg.qScale > 0 ? M.cfg.qScale : 1.;
    r.q0 = (P.rhat > 0) ? std::max(0., P.zero.F - P.hat.F) / qs : 0.;
    r.Z = std::sqrt(r.q0);
    r.p0 = 1. - x17c::Phi(r.Z);
    return r;
}

// ---------------------------------------------------------------- Feldman-Cousins: two-sided t_r with r >= 0
inline double TStat(const x17c::Model& M, const std::vector<double>& n, const x17c::Glob& g, const x17c::Profile& P, double r,
                    x17c::FitResult* fr = nullptr) {
    const std::vector<double>* st = (P.rhat >= 0 && !P.hat.p.empty()) ? &P.hat.p : &P.zero.p;
    x17c::FitResult f = M.Fit(n, r, false, g, st);
    if (!f.ok) { x17c::FitResult f2 = M.Fit(n, r, false, g); if (f2.F < f.F) f = f2; }
    if (fr) *fr = f;
    const double qs = M.cfg.qScale > 0 ? M.cfg.qScale : 1.;
    return std::max(0., f.F - P.Fref) / qs;
}
struct FCPoint { double r = 0, tobs = 0, t68 = 0, t90 = 0, t95 = 0, pval = 0; int nt = 0, nok = 0; };
inline double Quantile(std::vector<double> v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double x = q * (v.size() - 1);
    const size_t i = (size_t)std::floor(x);
    const double w = x - i;
    return i + 1 < v.size() ? (1 - w) * v[i] + w * v[i + 1] : v.back();
}
inline FCPoint FCToys(const x17c::Model& M, const std::vector<double>& nobs, const x17c::Glob& g, const x17c::Profile& P, double r, int ntoys,
                      unsigned long long seed) {
    FCPoint T;
    T.r = r;
    x17c::FitResult fr;
    T.tobs = TStat(M, nobs, g, P, r, &fr);
    const std::vector<double> pr = fr.p;   // nuisances profiled on the data at r (profile construction)
    std::mt19937_64 rng(seed);
    std::vector<double> n, ts;
    x17c::Glob gt;
    int nge = 0;
    for (int t = 0; t < ntoys; ++t) {
        x17c::ToyData(M, pr, rng, n, gt, g);
        const x17c::Profile Pt = x17c::ProfileFits(M, n, gt);
        if (!Pt.zero.ok) continue;
        const double tt = TStat(M, n, gt, Pt, r);
        ts.push_back(tt);
        if (tt >= T.tobs - 1e-9) ++nge;
    }
    T.nt = ntoys; T.nok = (int)ts.size();
    T.t68 = Quantile(ts, 0.68); T.t90 = Quantile(ts, 0.90); T.t95 = Quantile(ts, 0.95);
    T.pval = ts.empty() ? 0. : (double)nge / ts.size();
    return T;
}

// ---------------------------------------------------------------- model components of a cell
// out: 0 CB, 1 other background, 2 direct D, 3 conversion C, 4 eta-related (C2 + C3, G_eta), 5 X17
inline void Comps(const x17c::Model& M, const std::vector<double>& p, int i, int j, double out[6]) {
    const x17c::Inputs& I = *M.in;
    const x17c::Layout& L = M.L;
    const int b = I.G.idx(i, j), c = I.bcl[j];
    const double mu = M.Mu(p, i, j);
    std::vector<double> q = p;
    for (int jj = 0; jj < L.ny; ++jj) q[jj] = 0.;
    out[1] = mu - M.Mu(q, i, j);
    q = p; q[L.iR] = 0.;
    out[5] = mu - M.Mu(q, i, j);
    out[0] = p[L.iK + I.kcl[j]] * I.CB[b];
    out[2] = I.gP[c][i] * p[L.iN5] * I.D[j];
    const double ctf = 1. + (M.cfg.cTailSig > 0 ? M.cfg.cTailSig * p[L.iTC] * I.Ctail[j] : 0.);
    out[3] = I.gP[c][i] * p[L.iN4] * I.C[j] * ctf;
    out[4] = p[L.iA] * (I.E23[b] + p[L.iG] * I.Geta[b]);
}
// Pearson chi2 and deviance of a fit over a set of cells (all cells if use is empty)
struct Gof { double chi2 = 0, dev = 0; int ncell = 0, nfree = 0; std::vector<double> chi2row, resPeak, varPeak; };
inline Gof GoodnessOfFit(const x17c::Model& M, const x17c::FitResult& fit, const std::vector<double>& n, const std::vector<char>* use = nullptr,
                         double pkLo = 0.5325, double pkHi = 0.5575) {
    Gof G;
    const x17c::Inputs& I = *M.in;
    G.chi2row.assign(I.G.ny, 0.); G.resPeak.assign(I.G.ny, 0.); G.varPeak.assign(I.G.ny, 0.);
    for (int j = 0; j < I.G.ny; ++j)
        for (int i = 0; i < I.G.nx; ++i) {
            const int b = I.G.idx(i, j);
            if (use && !(*use)[b]) continue;
            const double mu = M.Mu(fit.p, i, j), nn = n[b];
            if (!(mu > 0)) continue;
            const double c = (nn - mu) * (nn - mu) / mu;
            G.chi2 += c; G.chi2row[j] += c;
            G.dev += 2 * (mu - nn + (nn > 0 ? nn * std::log(nn / mu) : 0.));
            ++G.ncell;
            if (I.G.xc(i) > pkLo && I.G.xc(i) < pkHi) { G.resPeak[j] += nn - mu; G.varPeak[j] += mu; }
        }
    for (size_t k = 0; k < fit.err.size(); ++k) if (fit.err[k] > 0) ++G.nfree;
    return G;
}

}  // namespace x17d

#endif  // ALLREPORT_X17D_LIB_H
