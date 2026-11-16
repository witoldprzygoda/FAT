// ============================================================================
// flex_lib.h  --  flexible SMASH description of the data (package flex of report_all)
// ============================================================================
// Model, bin i of fit histogram h = (M(pi+pi-) slice j, M(ee) class e), RAW data rows:
//
//   mu_i = k_e CB_i + sum_t coef_t(theta) f_t,i T_t,i
//
//   T_t   = SMASH components (sigma x PDG BR x V3, s = 0.99615): C2, C3, C4, C5, G_eta (raw simulated rows, Barlow-Beeston-lite
//           variance), G_D, G_C, their omega/rho0 split, T3 (smoothed exp(pol3) per histogram as in Part III, no variance)
//   coef  = theta[p1] (x theta[p2]): one free scale per component in the flexible model (products only for the Part III
//           configuration, which this library reproduces as a check: a_eta, a kappa, a rho_d, a g, kappa_B)
//   f     = 1 for the eta classes; for a "background" term c_j (1 + b_j x), x = (M4 - 0.55)/0.1 at the bin centre, with
//           optional per-slice normalisations c_j (Part III) and a common (b) or per-slice (b_j) tilt
//   k_e   = CB scale per class: k_E2 = 1 +- 0.05 (Part III D5), k_E1 free only when E1 is fitted
//
// Fitter: the WP4/WP5 estimator in general form (wp5_lib.h): Poisson / IRLS fixed point by Fisher scoring with the variance
// V = max(mu, 1) + k^2 var(CB) + sum coef^2 f^2 var(T_raw) (BB-lite); errors: sandwich with the event-group covariance of the
// counts (scaled to the model) + the CB covariance from the 200 like-sign replicas + one part per Gaussian constraint; MC =
// spread of 100 refits with the WP0 bootstrap replicas of the simulation (smoothing redone).
// Inputs: /hdd3/.../all_report/flex/flex_data.root, flex_sim.root (layout X of flex_fill.C: the F10 bins 0.40-0.70 + 10 MeV
// bins 0.70-1.00; classes E1 < 10, E2a 10-20, E2b >= 20 MeV). Summed over E2a + E2b, these reproduce hists_v1.root F10 E2
// exactly (flex_fill_checks.txt); flex_fit.C checks this again at run time.
// ============================================================================
#ifndef FLEX_LIB_H
#define FLEX_LIB_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TRandom3.h"

#include "s5_blind_ana_fit.h"
#include "hists_io.h"
#include "wp4_lib.h"

namespace flex {
using namespace alldefs;
typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;

const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/flex";
const char* const kDataFile = "/hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_data.root";
const char* const kSimFile = "/hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_sim.root";
const char* const kFigDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTabDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const int kNX = 100, kNCls = 3, kNHX = 15, kNRepS = 100, kNRepC = 200;
inline const char* ClsTag(int c) { static const char* n[kNCls] = {"E1", "E2a", "E2b"}; return n[c]; }
inline std::string XName(int j, int c) { return Form("X_S%d%s", j + 1, ClsTag(c)); }
inline std::vector<double> XEdges() {
    std::vector<double> e = F10Edges();
    for (int i = 1; i <= 30; ++i) e.push_back(std::round((0.700 + 0.010 * i) * 1e6) / 1e6);
    return e;
}

inline Mat Zero(int n, int m) { return Mat(n, Vec(m, 0.)); }
inline double Quad(const Mat& C, const Vec& g) {
    double s = 0;
    for (size_t k = 0; k < g.size(); ++k) { if (g[k] == 0) continue; for (size_t l = 0; l < g.size(); ++l) s += g[k] * C[k][l] * g[l]; }
    return s;
}
inline Mat Add(const Mat& a, const Mat& b) { Mat c = a; for (size_t i = 0; i < a.size(); ++i) for (size_t j = 0; j < a[i].size(); ++j) c[i][j] += b[i][j]; return c; }

// ===========================================================================
// 1. components
// ===========================================================================
struct CompDef {
    std::string name, text;
    std::vector<int> parts;        // alldefs::Comp summed
    bool smooth = false;           // smoothed exp(pol3) per histogram (non-eta) or raw with BB-lite (eta classes)
    std::vector<int> parent;       // parent for sparse smoothing (N_eff < 20): its shape in the same histogram
};
enum CompId { cC2 = 0, cC3, cC4, cC5, cGeta, cGD, cGC, cGDom, cGDnom, cGCom, cGCnom, cT3, cGpi0D, cGoth, kNCompF };
inline std::vector<CompDef> AllComps() {
    std::vector<CompDef> v(kNCompF);
    v[cC2] = {"C2", "bump, eta -> pi+pi-pi0, pi0 Dalitz", {kC2}, false, {}};
    v[cC3] = {"C3", "bump, eta -> pi+pi-pi0, photon conversion", {kC3}, false, {}};
    v[cC4] = {"C4", "peak, eta -> pi+pi-gamma, photon conversion", {kC4}, false, {}};
    v[cC5] = {"C5", "peak, direct eta -> pi+pi-e+e-", {kC5}, false, {}};
    v[cGeta] = {"Geta", "other eta (eta in the event, pair not of C2-C5)", {kGeta}, false, {}};
    v[cGD] = {"GD", "non-eta, Dalitz-type pair (pi0 Dalitz + vector + other)", {kGD}, true, {}};
    v[cGC] = {"GC", "non-eta, conversion pair", {kGconv}, true, {}};
    v[cGDom] = {"GD_om", "G_D with both pions from one omega or rho0", {kGpi0Dom, kGvecom, kGothom}, true, {kGD}};
    v[cGDnom] = {"GD_nom", "G_D, other pions", {kGpi0Dnom, kGvecnom, kGothnom}, true, {kGD}};
    v[cGCom] = {"GC_om", "G_C with both pions from one omega or rho0", {kGconvom}, true, {kGconv}};
    v[cGCnom] = {"GC_nom", "G_C, other pions", {kGconvnom}, true, {kGconv}};
    v[cT3] = {"T3", "fake lepton (not in the SMASH sum of Part III, D6)", {kT3}, true, {}};
    v[cGpi0D] = {"Gpi0D", "non-eta pi0 Dalitz pair", {kGpi0D}, true, {kGD}};
    v[cGoth] = {"Goth", "non-eta other / vector pair", {kGoth, kGvec}, true, {kGD}};
    return v;
}

// ===========================================================================
// 2. layout of a fit: fit classes x slices x M4 range
// ===========================================================================
struct FitClass { std::string tag; std::vector<int> xcls; int kidx; };   // kidx 0: k_E1, 1: k_E2
inline FitClass ClsE1() { return {"E1", {0}, 0}; }
inline FitClass ClsE2() { return {"E2", {1, 2}, 1}; }
inline FitClass ClsE2a() { return {"E2a", {1}, 1}; }
inline FitClass ClsE2b() { return {"E2b", {2}, 1}; }
struct Layout {
    std::vector<FitClass> fc;
    double mlo = 0.40, mhi = 0.70;
    int nh() const { return (int)fc.size() * kNSlice; }
    int J(int h) const { return h % kNSlice; }
    int F(int h) const { return h / kNSlice; }
    std::string HTag(int h) const { return Form("S%d%s", J(h) + 1, fc[F(h)].tag.c_str()); }
};

// ===========================================================================
// 3. templates (expected raw data rows) per component x fit histogram x X bin (1..100)
// ===========================================================================
struct Templ {
    int nc = 0, nh = 0;
    std::vector<Vec> c, v, sm;   // [comp * nh + h][x], x = 0..kNX+1
    Vec& C(int comp, int h) { return c[comp * nh + h]; }
    const Vec& C(int comp, int h) const { return c[comp * nh + h]; }
    const Vec& V(int comp, int h) const { return v[comp * nh + h]; }
    const Vec& S(int comp, int h) const { return sm[comp * nh + h]; }
};
// smoothing of one histogram in the pieces 0.40-0.70 (F10 edges) and 0.70-1.00 (10 MeV), wp4::SmoothFit (exp(pol3) in
// x = (M4 - 0.55)/0.1, one common scale, parent shape for N_eff < 20); returns bin integrals on the X axis
struct SmoothInfo { double neff[2] = {0, 0}; int order[2] = {-1, -1}; bool scaled[2] = {false, false}; double chi2[2] = {0, 0}; int nb[2] = {0, 0}; };
inline Vec SmoothX(const Vec& c, const Vec& v, const Vec* pc, const Vec* pv, bool full, SmoothInfo* info = nullptr) {
    const std::vector<double> E = XEdges();
    Vec out(kNX + 2, 0.);
    for (int piece = 0; piece < (full ? 2 : 1); ++piece) {
        const int x0 = piece == 0 ? 1 : 71, x1 = piece == 0 ? 70 : 100, n = x1 - x0 + 1;
        std::vector<double> e(E.begin() + (x0 - 1), E.begin() + x1 + 1);
        Vec cc(n + 1, 0.), vv(n + 1, 0.), pcc(n + 1, 0.), pvv(n + 1, 0.);
        for (int i = 1; i <= n; ++i) { cc[i] = c[x0 + i - 1]; vv[i] = v[x0 + i - 1]; if (pc) { pcc[i] = (*pc)[x0 + i - 1]; pvv[i] = (*pv)[x0 + i - 1]; } }
        wp4::Smooth par;
        if (pc) par = wp4::SmoothFit(pcc, pvv, e, nullptr, -1);
        const wp4::Smooth s = wp4::SmoothFit(cc, vv, e, pc ? &par : nullptr, -1);
        // diagnostic chi2 of the smoothing under its own scaled-Poisson model: var(c_i) = sbar f_i, sbar = sum v / sum c
        double ch = 0, C = 0, Vs = 0;
        for (int i = 1; i <= n; ++i) { C += cc[i]; Vs += vv[i]; }
        const double sbar = C > 0 ? Vs / C : 1.;
        for (int i = 1; i <= n; ++i) {
            out[x0 + i - 1] = s.Int(e[i - 1], e[i]);
            const double f = out[x0 + i - 1];
            if (f > 0) ch += std::pow(cc[i] - f, 2) / (sbar * std::max(f, sbar));
        }
        if (info) { info->neff[piece] = s.neff; info->order[piece] = s.order; info->scaled[piece] = s.scaled; info->chi2[piece] = ch; info->nb[piece] = n; }
    }
    return out;
}

struct SimSrc {   // all component histograms of the X layout (nb per bin) and replicas, loaded once
    TFile* f = nullptr;
    std::map<std::string, Vec> c, v;
    std::map<std::string, std::vector<Vec>> rep;   // [r][x]
    bool Open() { f = TFile::Open(kSimFile); return f && !f->IsZombie(); }
    void Need(int comp, int j, int k, bool withRep) {
        const std::string n = std::string(CompTag(comp)) + "_" + XName(j, k);
        if (!c.count(n)) {
            TH1* h = (TH1*)f->Get(n.c_str());
            Vec a(kNX + 2, 0.), b(kNX + 2, 0.);
            for (int x = 1; x <= kNX; ++x) { a[x] = h->GetBinContent(x) * kNbToRaw; b[x] = std::pow(h->GetBinError(x) * kNbToRaw, 2); }
            c[n] = a; v[n] = b;
        }
        if (withRep && !rep.count(n)) {
            TH2* h = (TH2*)f->Get((n + "_rep").c_str());
            std::vector<Vec> r(kNRepS, Vec(kNX + 2, 0.));
            for (int q = 0; q < kNRepS; ++q) for (int x = 1; x <= kNX; ++x) r[q][x] = h->GetBinContent(x, q + 1) * kNbToRaw;
            rep[n] = r;
        }
    }
};
// sum of the parts over the classes of the fit class; replica r < 0 = full sample (the variance is always the full one)
// (the variance of the merged class E2 = E2a + E2b is the event-group variance of the merged cell, histogram <comp>_X_S<j>E2;
//  the variance of a sum of components is the sum of their variances, as the hists_v1 GD histogram is not used for parts)
inline void SumParts(SimSrc& S, const std::vector<int>& parts, int j, const FitClass& fc, int r, Vec& c, Vec& v) {
    c.assign(kNX + 2, 0.); v.assign(kNX + 2, 0.);
    const bool merged = fc.xcls.size() == 2 && fc.xcls[0] == 1 && fc.xcls[1] == 2;
    for (int p : parts) {
        for (int k : fc.xcls) {
            S.Need(p, j, k, r >= 0);
            const std::string n = std::string(CompTag(p)) + "_" + XName(j, k);
            const Vec& a = r < 0 ? S.c[n] : S.rep[n][r];
            const Vec& b = S.v[n];
            for (int x = 1; x <= kNX; ++x) { c[x] += a[x]; if (!merged) v[x] += b[x]; }
        }
        if (merged) {
            const std::string n = std::string(CompTag(p)) + Form("_X_S%dE2", j + 1);
            if (!S.v.count(n)) {
                TH1* h = (TH1*)S.f->Get(n.c_str());
                Vec b(kNX + 2, 0.);
                for (int x = 1; x <= kNX; ++x) b[x] = std::pow(h->GetBinError(x) * kNbToRaw, 2);
                S.v[n] = b;
            }
            const Vec& b = S.v[n];
            for (int x = 1; x <= kNX; ++x) v[x] += b[x];
        }
    }
}
// templates of the components "use" (indices into AllComps) for a layout; r < 0 full sample, else replica r
inline Templ MakeTempl(SimSrc& S, const Layout& L, const std::vector<int>& use, int r, std::vector<SmoothInfo>* info = nullptr) {
    const std::vector<CompDef> D = AllComps();
    Templ T; T.nc = kNCompF; T.nh = L.nh();
    T.c.assign(T.nc * T.nh, Vec()); T.v.assign(T.nc * T.nh, Vec()); T.sm.assign(T.nc * T.nh, Vec());
    const bool full = L.mhi > 0.70 + 1e-9;
    if (info) info->assign(T.nc * T.nh, SmoothInfo());
    for (int comp : use)
        for (int h = 0; h < T.nh; ++h) {
            const int j = L.J(h);
            const FitClass& fc = L.fc[L.F(h)];
            Vec c, v;
            SumParts(S, D[comp].parts, j, fc, r, c, v);
            T.c[comp * T.nh + h] = c; T.v[comp * T.nh + h] = v;
            if (D[comp].smooth) {
                if (!D[comp].parent.empty()) {
                    Vec pc, pv;
                    SumParts(S, D[comp].parent, j, fc, r, pc, pv);
                    T.sm[comp * T.nh + h] = SmoothX(c, v, &pc, &pv, full, info ? &(*info)[comp * T.nh + h] : nullptr);
                } else T.sm[comp * T.nh + h] = SmoothX(c, v, nullptr, nullptr, full, info ? &(*info)[comp * T.nh + h] : nullptr);
            }
        }
    return T;
}

// ===========================================================================
// 4. data of a layout: bins, N, CB, covariances
// ===========================================================================
struct Bins {
    int n = 0;
    std::vector<int> h, j, f, e, x;   // fit histogram, slice, fit class, k index, X bin
    Vec y, cb, vcb, xc, lo, hi;
    std::vector<double> eg, ccb;      // dense n x n: event-group covariance of N, CB covariance (replicas)
    std::vector<Vec> cbrep;           // [r][i]
};
inline bool MakeBins(const Layout& L, Bins& B) {
    B = Bins();
    TFile* f = TFile::Open(kDataFile);
    if (!f || f->IsZombie()) return false;
    const std::vector<double> E = XEdges();
    TMatrixDSym* M = (TMatrixDSym*)f->Get("egcov_X");
    std::map<std::string, TH1*> hN, hC;
    std::map<std::string, TH2*> hR;
    std::vector<std::vector<int>> cellsOf;   // X cells of each fit bin
    for (int h = 0; h < L.nh(); ++h) {
        const int j = L.J(h);
        const FitClass& fc = L.fc[L.F(h)];
        for (int x = 1; x <= kNX; ++x) {
            if (E[x - 1] < L.mlo - 1e-9 || E[x] > L.mhi + 1e-9) continue;
            double y = 0, cb = 0;
            std::vector<int> cl;
            for (int k : fc.xcls) {
                const std::string n = XName(j, k);
                if (!hN.count(n)) { hN[n] = (TH1*)f->Get(("N_" + n).c_str()); hC[n] = (TH1*)f->Get(("CB_" + n).c_str()); hR[n] = (TH2*)f->Get(("CBrep_" + n).c_str()); }
                y += hN[n]->GetBinContent(x); cb += hC[n]->GetBinContent(x);
                cl.push_back((kNCls * j + k) * kNX + x - 1);
            }
            B.h.push_back(h); B.j.push_back(j); B.f.push_back(L.F(h)); B.e.push_back(fc.kidx); B.x.push_back(x);
            B.y.push_back(y); B.cb.push_back(cb); B.lo.push_back(E[x - 1]); B.hi.push_back(E[x]); B.xc.push_back(TiltX(0.5 * (E[x - 1] + E[x])));
            cellsOf.push_back(cl);
            ++B.n;
        }
    }
    const int n = B.n;
    B.eg.assign((size_t)n * n, 0.);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            double s = 0;
            for (int p : cellsOf[a]) for (int q : cellsOf[b]) s += (*M)(p, q);
            B.eg[(size_t)a * n + b] = s;
        }
    B.cbrep.assign(kNRepC, Vec(n, 0.));
    for (int i = 0; i < n; ++i)
        for (int k : L.fc[B.f[i]].xcls) {
            TH2* r = hR[XName(B.j[i], k)];
            for (int q = 0; q < kNRepC; ++q) B.cbrep[q][i] += r->GetBinContent(B.x[i], q + 1);
        }
    Vec m(n, 0.);
    for (int q = 0; q < kNRepC; ++q) for (int i = 0; i < n; ++i) m[i] += B.cbrep[q][i] / kNRepC;
    B.ccb.assign((size_t)n * n, 0.);
    for (int a = 0; a < n; ++a)
        for (int b = a; b < n; ++b) {
            double s = 0;
            for (int q = 0; q < kNRepC; ++q) s += (B.cbrep[q][a] - m[a]) * (B.cbrep[q][b] - m[b]);
            s /= kNRepC;   // 1/n as wp5_lib.h
            B.ccb[(size_t)a * n + b] = B.ccb[(size_t)b * n + a] = s;
        }
    B.vcb.assign(n, 0.);
    for (int i = 0; i < n; ++i) B.vcb[i] = B.ccb[(size_t)i * n + i];
    f->Close();
    return true;
}

// ===========================================================================
// 5. model
// ===========================================================================
struct Term { int comp, p1, p2; bool bg; int emask; };   // emask: bit k = applies to k index k (1 = E1, 2 = E2, 3 = both)
struct Model {
    std::string tag, text;
    std::vector<std::string> names;
    Vec start;
    std::vector<int> fixed;
    std::vector<bla::Constraint> cons;
    std::vector<Term> terms;
    int iK[2] = {-1, -1}; double kFix[2] = {1, 1};
    int iC[kNSlice] = {-1, -1, -1, -1, -1};
    int iB = -1;
    int iBj[kNSlice] = {-1, -1, -1, -1, -1};
    std::vector<int> scaleOf;   // parameter index of the free scale of each component (-1 if none), size kNCompF
    int np() const { return (int)names.size(); }
    int Idx(const std::string& n) const { for (int k = 0; k < np(); ++k) if (names[k] == n) return k; return -1; }
    int Add(const std::string& n, double st) { names.push_back(n); start.push_back(st); fixed.push_back(0); return np() - 1; }
    std::vector<int> Comps() const { std::vector<int> v; for (const auto& t : terms) if (std::find(v.begin(), v.end(), t.comp) == v.end()) v.push_back(t.comp); return v; }
};
struct ModelOpt {
    std::vector<int> comps;       // components, each with its own free scale
    int tilt = 0;                 // 0 none, 1 common b, 2 per-slice b_j (background terms)
    bool perSliceNorm = false;    // c_j on the background terms (then the scale of the first background component is fixed to 1)
    bool kE1free = true;          // when E1 is fitted: k_E1 free (CB-like term), else fixed to 1
    double kE2sig = 0.05;
    std::vector<int> e1only;      // components that apply to E1 only (e.g. T3 as the fake-lepton term)
    std::map<int, double> fixScale;   // components with the scale held at a value
    std::map<int, std::pair<double, double>> consScale;   // Gaussian constraint on a scale
};
inline bool IsBg(int comp) { const std::vector<CompDef> D = AllComps(); return D[comp].smooth; }
inline Model MakeFlex(const ModelOpt& o, bool hasE1, bool hasE2, const std::string& tag, const std::string& text) {
    Model M; M.tag = tag; M.text = text;
    M.scaleOf.assign(kNCompF, -1);
    const std::vector<CompDef> D = AllComps();
    bool firstBg = true;
    for (int c : o.comps) {
        const int p = M.Add("s_" + D[c].name, 1.0);
        M.scaleOf[c] = p;
        int em = 3;
        if (std::find(o.e1only.begin(), o.e1only.end(), c) != o.e1only.end()) em = 1;
        M.terms.push_back({c, p, -1, D[c].smooth, em});
        if (o.fixScale.count(c)) { M.start[p] = o.fixScale.at(c); M.fixed[p] = 1; }
        if (o.consScale.count(c)) M.cons.push_back({p, o.consScale.at(c).first, o.consScale.at(c).second});
        if (o.perSliceNorm && D[c].smooth && firstBg) { M.fixed[p] = 1; M.start[p] = 1.0; firstBg = false; }
    }
    if (o.perSliceNorm) for (int j = 0; j < kNSlice; ++j) M.iC[j] = M.Add(Form("c_%d", j + 1), 1.0);
    if (o.tilt == 1) M.iB = M.Add("b", 0.0);
    if (o.tilt == 2) for (int j = 0; j < kNSlice; ++j) M.iBj[j] = M.Add(Form("b_%d", j + 1), 0.0);
    if (hasE1) { M.iK[0] = M.Add("k_E1", 1.0); if (!o.kE1free) M.fixed[M.iK[0]] = 1; }
    if (hasE2) { M.iK[1] = M.Add("k_E2", 1.0); M.cons.push_back({M.iK[1], 1.0, o.kE2sig}); }
    return M;
}
// Part III F-main (check): a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta] + c_j (1 + b x)[G_D + kappa_B G_C]
inline Model MakePartIII(bool hasE1, bool hasE2) {
    Model M; M.tag = "P3"; M.text = "Part III F-main (kappa = 1 +- 0.2, g = 1 +- 0.3, k_E2 = 1 +- 0.05)";
    M.scaleOf.assign(kNCompF, -1);
    const int ia = M.Add("a_eta", 0.4), ik = M.Add("kappa", 1.0), ir = M.Add("rho_d", 1.0), ig = M.Add("g", 1.0);
    M.cons.push_back({ik, 1.0, 0.2}); M.cons.push_back({ig, 1.0, 0.3});
    M.terms = {{cC2, ia, -1, false, 3}, {cC3, ia, ik, false, 3}, {cC4, ia, ik, false, 3}, {cC5, ia, ir, false, 3}, {cGeta, ia, ig, false, 3}};
    for (int j = 0; j < kNSlice; ++j) M.iC[j] = M.Add(Form("c_%d", j + 1), 1.0);
    M.iB = M.Add("b", -0.07);
    const int ikb = M.Add("kappa_B", 2.0);
    M.terms.push_back({cGD, -1, -1, true, 3});
    M.terms.push_back({cGC, ikb, -1, true, 3});
    if (hasE1) M.iK[0] = M.Add("k_E1", 1.0);
    if (hasE2) { M.iK[1] = M.Add("k_E2", 1.0); M.cons.push_back({M.iK[1], 1.0, 0.05}); }
    return M;
}

// per-bin template values of the model's components (smoothed or raw)
struct BinT { std::vector<Vec> t, v; };   // [comp][i]
inline BinT BinTempl(const Bins& B, const Templ& T, const Model& M) {
    BinT o; o.t.assign(kNCompF, Vec()); o.v.assign(kNCompF, Vec());
    const std::vector<CompDef> D = AllComps();
    for (int c : M.Comps()) {
        o.t[c].assign(B.n, 0.); o.v[c].assign(B.n, 0.);
        for (int i = 0; i < B.n; ++i) {
            if (D[c].smooth) o.t[c][i] = T.S(c, B.h[i])[B.x[i]];
            else { o.t[c][i] = T.C(c, B.h[i])[B.x[i]]; o.v[c][i] = T.V(c, B.h[i])[B.x[i]]; }
        }
    }
    return o;
}
// evaluation: mu, Jacobian, V (with BB-lite), Vd (data part)
inline void Eval(const Bins& B, const BinT& BT, const Model& M, const Vec& th, Vec& mu, Mat* J, Vec* V, Vec* Vd,
                 std::vector<Vec>* parts = nullptr) {
    const int n = B.n, np = M.np();
    mu.assign(n, 0.);
    if (J) J->assign(n, Vec(np, 0.));
    if (V) V->assign(n, 0.);
    if (Vd) Vd->assign(n, 0.);
    if (parts) parts->assign(kNCompF + 1, Vec(n, 0.));   // last = k CB
    for (int i = 0; i < n; ++i) {
        const int e = B.e[i], j = B.j[i];
        const double k = M.iK[e] >= 0 ? th[M.iK[e]] : M.kFix[e];
        const double c = M.iC[j] >= 0 ? th[M.iC[j]] : 1.;
        const double b = M.iBj[j] >= 0 ? th[M.iBj[j]] : (M.iB >= 0 ? th[M.iB] : 0.);
        const double x = B.xc[i], tl = 1. + b * x;
        double m = k * B.cb[i], vs = 0;
        if (parts) (*parts)[kNCompF][i] = k * B.cb[i];
        for (const auto& T : M.terms) {
            if (!((T.emask >> e) & 1)) continue;
            const double a1 = T.p1 >= 0 ? th[T.p1] : 1., a2 = T.p2 >= 0 ? th[T.p2] : 1.;
            const double co = a1 * a2;
            const double f = T.bg ? c * tl : 1.;
            const double tv = BT.t[T.comp][i];
            m += co * f * tv;
            vs += co * co * f * f * BT.v[T.comp][i];
            if (parts) (*parts)[T.comp][i] += co * f * tv;
            if (J) {
                if (T.p1 >= 0) (*J)[i][T.p1] += a2 * f * tv;
                if (T.p2 >= 0) (*J)[i][T.p2] += a1 * f * tv;
                if (T.bg) {
                    if (M.iC[j] >= 0) (*J)[i][M.iC[j]] += co * tl * tv;
                    if (M.iBj[j] >= 0) (*J)[i][M.iBj[j]] += co * c * x * tv;
                    else if (M.iB >= 0) (*J)[i][M.iB] += co * c * x * tv;
                }
            }
        }
        if (J && M.iK[e] >= 0) (*J)[i][M.iK[e]] = B.cb[i];
        mu[i] = m;
        const double vb = std::max(m, 1.0) + k * k * B.vcb[i];
        if (V) (*V)[i] = vb + vs;
        if (Vd) (*Vd)[i] = vb;
    }
}

// ===========================================================================
// 6. fitter (wp5_lib.h Fit, generalised)
// ===========================================================================
struct FitRes {
    bool ok = false;
    int nit = 0, n = 0, nfree = 0, ndf = 0;
    Vec th;
    Mat C, Cn, Ccb, Cstat;
    std::vector<Mat> Ccons;
    Vec mu, V, Vd;
    double chi2P = 0, chi2BC = 0;
};
// sandwich: N part uses the event-group covariance eg (scaled to the model: eg_ab sqrt(mu_a mu_b / (y_a y_b)));
// for toys eg = nullptr: Poisson (diag mu)
inline FitRes FitOnce(const Bins& B, const BinT& BT, const Model& M, const Vec* start, bool sandwich,
                      const std::vector<double>* eg, int maxIt, double damp);
// Fisher scoring; if the plain iteration does not converge (oscillating IRLS fixed point), it is repeated with damped steps
// (x 0.5, then x 0.25, more iterations); the fixed point is the same, only the path changes
inline FitRes Fit(const Bins& B, const BinT& BT, const Model& M, const Vec* start = nullptr, bool sandwich = true,
                  const std::vector<double>* eg = nullptr, int maxIt = 400) {
    FitRes R = FitOnce(B, BT, M, start, sandwich, eg, maxIt, 1.0);
    if (!R.ok) R = FitOnce(B, BT, M, start, sandwich, eg, 3 * maxIt, 0.5);
    if (!R.ok) R = FitOnce(B, BT, M, start, sandwich, eg, 8 * maxIt, 0.25);
    return R;
}
inline FitRes FitOnce(const Bins& B, const BinT& BT, const Model& M, const Vec* start, bool sandwich,
                      const std::vector<double>* eg, int maxIt, double damp) {
    FitRes R;
    const int n = B.n, np = M.np();
    R.n = n;
    std::vector<int> fr;
    for (int k = 0; k < np; ++k) if (!M.fixed[k]) fr.push_back(k);
    const int nf = (int)fr.size();
    R.nfree = nf;
    Vec th = start ? *start : M.start;
    for (int k = 0; k < np; ++k) if (M.fixed[k]) th[k] = M.start[k];
    Vec mu, V;
    Mat J;
    Vec Pd(nf, 0.), P0(nf, 0.);
    for (const auto& c : M.cons) for (int a = 0; a < nf; ++a) if (fr[a] == c.k) { Pd[a] += 1. / (c.sig * c.sig); P0[a] = c.val; }
    bool conv = false;
    for (int it = 0; it < maxIt; ++it) {
        Eval(B, BT, M, th, mu, &J, &V, nullptr);
        Mat A = Zero(nf, nf); Vec b(nf, 0.);
        for (int i = 0; i < n; ++i) {
            const double w = 1. / V[i], r = B.y[i] - mu[i];
            for (int a = 0; a < nf; ++a) {
                const double ja = J[i][fr[a]] * w;
                if (ja == 0) continue;
                b[a] += ja * r;
                for (int c = 0; c <= a; ++c) A[a][c] += ja * J[i][fr[c]];
            }
        }
        for (int a = 0; a < nf; ++a) { A[a][a] += Pd[a]; b[a] -= Pd[a] * (th[fr[a]] - P0[a]); }
        for (int a = 0; a < nf; ++a) for (int c = a + 1; c < nf; ++c) A[a][c] = A[c][a];
        Vec sc(nf, 1.);
        for (int a = 0; a < nf; ++a) sc[a] = A[a][a] > 0 ? 1. / std::sqrt(A[a][a]) : 1.;
        Mat As = A;
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) As[a][c] *= sc[a] * sc[c];
        if (!bla::CholInv(As)) { R.nit = it; return R; }
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) As[a][c] *= sc[a] * sc[c];
        Vec d(nf, 0.);
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) d[a] += As[a][c] * b[c];
        double dmax = 0, lim = 1.;
        for (int a = 0; a < nf; ++a) {
            const double z = std::fabs(d[a]) / std::sqrt(std::max(As[a][a], 1e-300));
            dmax = std::max(dmax, z);
            if (z > 3.) lim = std::min(lim, 3. / z);
        }
        for (int a = 0; a < nf; ++a) th[fr[a]] += damp * lim * d[a];
        R.nit = it + 1;
        if (it > 0 && dmax < 1e-5) { conv = true; break; }
    }
    if (!conv) return R;
    Vec Vd;
    Eval(B, BT, M, th, mu, &J, &V, &Vd);
    Mat A = Zero(nf, nf);
    for (int i = 0; i < n; ++i) for (int a = 0; a < nf; ++a) { const double ja = J[i][fr[a]] / V[i]; if (ja == 0) continue; for (int c = 0; c < nf; ++c) A[a][c] += ja * J[i][fr[c]]; }
    for (int a = 0; a < nf; ++a) A[a][a] += Pd[a];
    Vec sc(nf, 1.);
    for (int a = 0; a < nf; ++a) sc[a] = A[a][a] > 0 ? 1. / std::sqrt(A[a][a]) : 1.;
    for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) A[a][c] *= sc[a] * sc[c];
    if (!bla::CholInv(A)) return R;
    for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) A[a][c] *= sc[a] * sc[c];
    auto embed = [&](const Mat& m) { Mat o = Zero(np, np); for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) o[fr[a]][fr[c]] = m[a][c]; return o; };
    R.C = embed(A);
    auto sand = [&](const Mat& Mm) {
        Mat T = Zero(nf, nf), O = Zero(nf, nf);
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { double s = 0; for (int k = 0; k < nf; ++k) s += A[a][k] * Mm[k][c]; T[a][c] = s; }
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { double s = 0; for (int k = 0; k < nf; ++k) s += T[a][k] * A[k][c]; O[a][c] = s; }
        return embed(O);
    };
    if (sandwich) {
        Mat U = Zero(n, nf);
        Vec sf(n, 0.), kk(n, 0.);
        for (int i = 0; i < n; ++i) {
            for (int a = 0; a < nf; ++a) U[i][a] = J[i][fr[a]] / V[i];
            sf[i] = B.y[i] > 0 ? std::sqrt(std::max(mu[i], 1e-9) / B.y[i]) : 0.;
            kk[i] = M.iK[B.e[i]] >= 0 ? th[M.iK[B.e[i]]] : M.kFix[B.e[i]];
        }
        Mat Mn = Zero(nf, nf), Mc = Zero(nf, nf);
        Vec tn(nf), tc(nf);
        for (int i = 0; i < n; ++i) {
            std::fill(tn.begin(), tn.end(), 0.); std::fill(tc.begin(), tc.end(), 0.);
            const size_t ri = (size_t)i * n;
            for (int l = 0; l < n; ++l) {
                double sn;
                if (eg) {
                    sn = (*eg)[ri + l];
                    if (B.y[i] > 0 && B.y[l] > 0) sn *= sf[i] * sf[l];
                    else sn = (i == l) ? std::max(mu[i], 0.) : 0.;
                } else sn = (i == l) ? std::max(mu[i], 0.) : 0.;
                const double scb = B.ccb[ri + l] * kk[i] * kk[l];
                if (sn == 0 && scb == 0) continue;
                for (int a = 0; a < nf; ++a) { tn[a] += sn * U[l][a]; tc[a] += scb * U[l][a]; }
            }
            for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { Mn[a][c] += U[i][a] * tn[c]; Mc[a][c] += U[i][a] * tc[c]; }
        }
        R.Cn = sand(Mn); R.Ccb = sand(Mc);
        R.Cstat = Add(R.Cn, R.Ccb);
        for (const auto& c : M.cons) {
            Mat Mp = Zero(nf, nf);
            for (int a = 0; a < nf; ++a) if (fr[a] == c.k) Mp[a][a] = 1. / (c.sig * c.sig);
            R.Ccons.push_back(sand(Mp));
            R.Cstat = Add(R.Cstat, R.Ccons.back());
        }
    } else {
        R.Cn = R.C; R.Ccb = Zero(np, np); R.Cstat = R.C;
    }
    R.th = th; R.mu = mu; R.V = V; R.Vd = Vd;
    for (int i = 0; i < n; ++i) {
        const double m = std::max(mu[i], 1e-9), y = B.y[i];
        R.chi2BC += 2. * (m - y + (y > 0 ? y * std::log(y / m) : 0.));
        R.chi2P += (y - mu[i]) * (y - mu[i]) / V[i];
    }
    for (const auto& c : M.cons) { const double z = (th[c.k] - c.val) / c.sig; R.chi2BC += z * z; R.chi2P += z * z; }
    R.ndf = n - nf + (int)M.cons.size();
    R.ok = true;
    return R;
}

// ===========================================================================
// 7. replicas (MC error), hybrid GLS chi2
// ===========================================================================
struct RepSet { std::vector<BinT> bt; std::vector<Templ> T; };
inline RepSet MakeReplicaSet(SimSrc& S, const Layout& L, const Bins& B, const Model& M, int nrep = kNRepS) {
    RepSet R;
    const std::vector<int> use = M.Comps();
    for (int r = 0; r < nrep; ++r) {
        Templ T = MakeTempl(S, L, use, r);
        R.bt.push_back(BinTempl(B, T, M));
        R.T.push_back(T);
    }
    return R;
}
struct RepFits { std::vector<FitRes> fit; Mat cov; Vec mean; int nok = 0; };
inline RepFits RunReplicas(const Bins& B, const Model& M, const RepSet& R, const Vec& th0) {
    RepFits o;
    const int np = M.np();
    o.fit.resize(R.bt.size());
    for (size_t r = 0; r < R.bt.size(); ++r) o.fit[r] = Fit(B, R.bt[r], M, &th0, false);
    o.mean.assign(np, 0.); o.cov = Zero(np, np);
    for (auto& f : o.fit) if (f.ok) { ++o.nok; for (int k = 0; k < np; ++k) o.mean[k] += f.th[k]; }
    for (int k = 0; k < np; ++k) o.mean[k] /= std::max(1, o.nok);
    for (auto& f : o.fit) if (f.ok) for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) o.cov[k][l] += (f.th[k] - o.mean[k]) * (f.th[l] - o.mean[l]);
    for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) o.cov[k][l] /= std::max(1, o.nok - 1);
    return o;
}
inline bool InvSym(std::vector<double>& Mv, int n) {
    TMatrixDSym m(n);
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) m(a, b) = Mv[(size_t)a * n + b];
    TDecompChol ch(m);
    if (!ch.Decompose()) return false;
    TMatrixDSym inv(n);
    ch.Invert(inv);
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) Mv[(size_t)a * n + b] = inv(a, b);
    return true;
}
// hybrid GLS (v5 verification): data covariance (event groups scaled to the model + CB replicas) + exact diagonal BB-lite
// variance of the raw templates + replica covariance of the smoothed (background) part of the prediction
struct GLSRes { double chi2 = 0, chi2data = 0; int ndf = 0; bool ok = false; };
inline GLSRes HybridGLS(const Bins& B, const BinT& BT, const Model& M, const FitRes& F, const RepSet& R, bool bgRepCov = true) {
    GLSRes g;
    const int n = B.n;
    const std::vector<CompDef> D = AllComps();
    std::vector<double> W((size_t)n * n, 0.);
    Vec sf(n), kk(n);
    for (int i = 0; i < n; ++i) { sf[i] = B.y[i] > 0 ? std::sqrt(std::max(F.mu[i], 1e-9) / B.y[i]) : 0.; kk[i] = M.iK[B.e[i]] >= 0 ? F.th[M.iK[B.e[i]]] : M.kFix[B.e[i]]; }
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            double s = B.eg[(size_t)a * n + b];
            if (B.y[a] > 0 && B.y[b] > 0) s *= sf[a] * sf[b]; else s = (a == b) ? std::max(F.mu[a], 1.) : 0.;
            s += kk[a] * kk[b] * B.ccb[(size_t)a * n + b];
            W[(size_t)a * n + b] = s;
        }
    std::vector<double> Wd = W;
    // BB-lite diagonal (raw templates): V - Vd of the fit
    for (int i = 0; i < n; ++i) W[(size_t)i * n + i] += F.V[i] - F.Vd[i];
    // replica covariance of the smoothed part
    std::vector<Vec> bg;
    for (const auto& bt : R.bt) {
        Vec mu; std::vector<Vec> parts;
        Eval(B, bt, M, F.th, mu, nullptr, nullptr, nullptr, &parts);
        Vec s(n, 0.);
        if (bgRepCov) for (int c = 0; c < kNCompF; ++c) if (D[c].smooth) for (int i = 0; i < n; ++i) s[i] += parts[c][i];
        bg.push_back(s);
    }
    const int nr = (int)bg.size();
    Vec m(n, 0.);
    for (auto& v : bg) for (int i = 0; i < n; ++i) m[i] += v[i] / nr;
    for (auto& v : bg) for (int a = 0; a < n; ++a) { const double da = v[a] - m[a]; if (da == 0) continue; for (int b = 0; b < n; ++b) W[(size_t)a * n + b] += da * (v[b] - m[b]) / (nr - 1); }
    if (!InvSym(W, n) || !InvSym(Wd, n)) return g;
    Vec r(n);
    for (int i = 0; i < n; ++i) r[i] = B.y[i] - F.mu[i];
    for (int a = 0; a < n; ++a) { double s = 0, sd = 0; for (int b = 0; b < n; ++b) { s += W[(size_t)a * n + b] * r[b]; sd += Wd[(size_t)a * n + b] * r[b]; } g.chi2 += r[a] * s; g.chi2data += r[a] * sd; }
    for (const auto& c : M.cons) { const double z = (F.th[c.k] - c.val) / c.sig; g.chi2 += z * z; g.chi2data += z * z; }
    g.ndf = F.ndf;
    g.ok = true;
    return g;
}

// errors of a fit: stat (sandwich), MC (replicas), total
struct Summary {
    Vec v, eStat, eMC, eTot, shift, eN, eCB;
    Mat covTot, corTot, covStat;
};
inline Summary Summarise(const Model& M, const FitRes& F, const RepFits* R) {
    Summary s;
    const int np = M.np();
    s.v = F.th; s.eStat.assign(np, 0.); s.eMC.assign(np, 0.); s.eTot.assign(np, 0.); s.shift.assign(np, 0.); s.eN.assign(np, 0.); s.eCB.assign(np, 0.);
    s.covStat = F.Cstat;
    s.covTot = F.Cstat;
    if (R) s.covTot = Add(F.Cstat, R->cov);
    s.corTot = Zero(np, np);
    for (int k = 0; k < np; ++k) {
        s.eStat[k] = std::sqrt(std::max(0., F.Cstat[k][k]));
        s.eN[k] = std::sqrt(std::max(0., F.Cn[k][k]));
        s.eCB[k] = std::sqrt(std::max(0., F.Ccb[k][k]));
        if (R) { s.eMC[k] = std::sqrt(std::max(0., R->cov[k][k])); s.shift[k] = R->mean[k] - F.th[k]; }
        s.eTot[k] = std::sqrt(std::max(0., s.covTot[k][k]));
    }
    for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) s.corTot[k][l] = (s.eTot[k] > 0 && s.eTot[l] > 0) ? s.covTot[k][l] / (s.eTot[k] * s.eTot[l]) : 0.;
    return s;
}

}  // namespace flex

#endif
