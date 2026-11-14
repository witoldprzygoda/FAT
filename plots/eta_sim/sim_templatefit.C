// ============================================================================
// sim_templatefit.C  --  PRELIMINARY template fit of the data M(pi+pi-e+e-)
// ============================================================================
// Global fit of the data spectrum N - CB in M4 = M(pi+pi-e+e-) with simulated
// templates (one scale factor per component), as in the X17 analysis of this
// channel. PRELIMINARY: the eta -> pi+pi-e+e- template (eta set) is missing,
// and L_data = 6080.5 nb^-1 is unconfirmed.
//
// Data   : N - CB per bin, CB = 2 sqrt(N++ N--) in 10 MeV bins interpolated to
//          2.5 MeV (ScanEngine); covariance = event-group covariance of N
//          (rows of one event in the same or in different bins) + covariance of
//          the interpolated CB (like-sign event groups, linear propagation).
// Model  : templates in data rows (sum w * L_data / L_sim, L_data = 6080.5):
//          alpha = 1 means agreement with SMASH.
//            bump  = signal classes 2 + 3 (pi0 Dalitz + conversion; their
//                    relative share fixed by the simulation)
//            peak  = signal class 4 (eta -> pi+pi-gamma, converted photon)
//            bkg   = template 2 (physical background); split variant:
//                    omega / rho0 pion pairs (sim_pcls 3, 4) and the rest
//            G     = Gaussian (unit area) at the data eta position 0.5439 with
//                    the width of the peak-conversion template, free content:
//                    stands for the missing eta -> pi+pi-e+e-
//          M1 3 templates; M2 M1 + G; M3 M1 + G with a_peak = a_bump (one eta
//          factor); M4 split background; M5 M4 + G, a_peak = a_bump;
//          M6 M1 with a linear tilt of the background, a_bkg + b (M - 0.6)/0.1;
//          M7 M6 + G, a_peak = a_bump.
// chi2   : (y - F p)^T V^-1 (y - F p), V = V_data + sum_kl c_k c_l C_kl
//          (C_kl = event-group covariance of the templates k, l; c = their
//          fitted coefficients, bin dependent with the tilt); V re-evaluated at
//          the current p until p is stable (iterated generalised least
//          squares; no Poisson likelihood on weighted samples).
// Implied contents (data rows):
//          S = eta peak content by the data-study fit P_main (Gaussian 0.5439 /
//              5.4 MeV + pol2 in 0.5225-0.61) applied to the fitted model;
//          E = edge-step sum of the data study (M(pi+pi-) slices 0.30-0.34,
//              0.34-0.37, 0.37-0.39, 0.39-0.40; scan::FitEdgeStep) applied to
//              the fitted model; both against the same measures on the data.
//          Errors: (a) from the fitted parameters (covariance), (b) from the
//          template statistics at fixed parameters (linear propagation of the
//          template covariance through the numerical gradient of S and E).
//
// Outputs: figs/fit_*.png (.pdf for the main ones), tables/fit_*.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_templatefit", kTRUE)' 'plots/eta_sim/sim_templatefit.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"

#include "SimEngine.h"

namespace stf {
using namespace simeng;

// ---------------------------------------------------------------- constants
const int kNB = 320;                        // fine grid: 2.5 MeV, 0.30-1.10
const double kLo = 0.30, kHi = 1.10;
const double kBW = (kHi - kLo) / kNB;
const double kU = kLdata / kLsim;           // sum w -> expected data rows
const double kTilt0 = 0.6, kTiltW = 0.1;    // tilt variable t = (M - 0.6) / 0.1
inline double TiltT(double x) { return (x - kTilt0) / kTiltW; }

enum { cBump = 0, cPeak, cVM, cRest, kNC };
const char* const kCTag[kNC] = {"bump", "peak", "bkgVM", "bkgRest"};
inline bool IsBkg(int k) { return k == cVM || k == cRest; }
inline bool IsVM(const Row& r) { return r.pcls == 3 || r.pcls == 4; }
std::vector<Comp> FitComps() {
    std::vector<Comp> c(kNC);
    c[cBump] = {"bump (pi0 Dalitz + conversion)", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); }};
    c[cPeak] = {"peak (conversion)", [](const Row& r) { return r.tmpl == 1 && r.cls == 4; }};
    c[cVM]   = {"phys. bkg, omega/rho0 pions", [](const Row& r) { return r.tmpl == 2 && IsVM(r); }};
    c[cRest] = {"phys. bkg, other pions", [](const Row& r) { return r.tmpl == 2 && !IsVM(r); }};
    return c;
}

const int kNP = 4;
const int kPresetId[kNP] = {kPre, kPreOA2, kSlides, kPeakNom};
const char* const kPTag[kNP] = {"pre", "preOA2", "slides", "peaknom"};
const char* const kPName[kNP] = {"pre", "pre+OA>2", "slides", "peak nominal"};

// template options: 0 mode A, w; 1 mode B, w; 2 mode A, w_db; 3 mode A, w, simulated M4 scaled
const int kNO = 4;
// review finding 8: the scale variant uses the K0S M(pi+pi-) scale (data 0.9970 / sim 1.0018, tables/spectra_eta.txt),
// not the offset of the conversion-only simulated peak against the all-decay data peak; its transfer to M4 is UNVERIFIED
double gSimScale = 0.9970 / 1.0018;
inline SimOpt OptOf(int o) { SimOpt s; if (o == 1) s.mode = 1; if (o == 2) s.db = true; return s; }
inline double ScaleOf(int o) { return o == 3 ? gSimScale : 1.0; }
const char* const kOTag[kNO] = {"A_w", "B_w", "A_wdb", "A_w_sc"};

const int kNR = 3;
const double kRLo[kNR] = {0.43, 0.45, 0.40}, kRHi[kNR] = {0.75, 0.70, 0.85};
const int kNBW = 2;
const int kReb[kNBW] = {1, 2};               // 2.5 and 5 MeV

// M(pi+pi-) slices of the edge-step sum E (as plots/scan_bump.C)
const int kNSl = 4;
const double kSlLo[kNSl] = {0.30, 0.34, 0.37, 0.39}, kSlHi[kNSl] = {0.34, 0.37, 0.39, 0.40};
inline scan::StepCfg SliceCfg(int i) {
    scan::StepCfg c; c.mpLo = kSlLo[i]; c.mpHi = kSlHi[i];
    if (i == 3) c.halfRange = 0.035;
    return c;
}
const int kSlOff = 20;                       // fine bin of 0.35 (start of the ScanEngine grid)

struct Model {
    std::string tag, label;
    std::vector<std::vector<int>> cols;      // templates multiplied by each alpha
    bool tilt, gauss;
    std::vector<std::string> pname;          // alphas, [b_tilt], [A_G]
};
std::vector<Model> Models() {
    std::vector<Model> m(7);
    m[0] = {"M1", "3 templates", {{cBump}, {cPeak}, {cVM, cRest}}, false, false, {"a_bump", "a_peak", "a_bkg"}};
    m[1] = {"M2", "3 templates + G", {{cBump}, {cPeak}, {cVM, cRest}}, false, true, {"a_bump", "a_peak", "a_bkg", "A_G"}};
    m[2] = {"M3", "3 templates + G, a_peak = a_bump", {{cBump, cPeak}, {cVM, cRest}}, false, true, {"a_eta", "a_bkg", "A_G"}};
    m[3] = {"M4", "split background", {{cBump}, {cPeak}, {cVM}, {cRest}}, false, false, {"a_bump", "a_peak", "a_VM", "a_rest"}};
    m[4] = {"M5", "split background + G, a_peak = a_bump", {{cBump, cPeak}, {cVM}, {cRest}}, false, true, {"a_eta", "a_VM", "a_rest", "A_G"}};
    m[5] = {"M6", "3 templates, background tilt", {{cBump}, {cPeak}, {cVM, cRest}}, true, false, {"a_bump", "a_peak", "a_bkg", "b_tilt"}};
    m[6] = {"M7", "background tilt + G, a_peak = a_bump", {{cBump, cPeak}, {cVM, cRest}}, true, true, {"a_eta", "a_bkg", "b_tilt", "A_G"}};
    return m;
}

// coefficients of the templates for a parameter vector
struct Coef {
    double c[kNC] = {0, 0, 0, 0};
    double beta = 0, AG = 0;
    double at(int k, double x) const { return c[k] + (IsBkg(k) ? beta * TiltT(x) : 0.); }
};
inline Coef CoefOf(const Model& m, const std::vector<double>& p) {
    Coef q;
    const int nt = (int)m.cols.size();
    for (int j = 0; j < nt; ++j) for (int k : m.cols[j]) q.c[k] = p[j];
    int j = nt;
    if (m.tilt) q.beta = p[j++];
    if (m.gauss) q.AG = p[j++];
    return q;
}

// ---------------------------------------------------------------- templates
struct TplSet {
    std::vector<double> s[kNC];               // data rows per 2.5 MeV bin (alpha = 1)
    std::vector<double> C[kNC][kNC];          // event-group covariance, data rows^2, kNB x kNB
    double Sall[kNC], Vall[kNC], Nrow[kNC];   // all selected rows (M4 0.30-1.10 loaded)
    std::vector<double> sl[kNSl][kNC];        // per M(pi+pi-) slice, data rows per 2.5 MeV bin
    std::vector<double> slv[kNSl][kNC];       // their diagonal group variances
    double slPeakFrac[kNSl];                  // share of the peak template in each slice
};

void FillTpl(const SimData& d, const Sel& sel, const SimOpt& opt, double scale, TplSet& t) {
    const std::vector<Comp> comps = FitComps();
    for (int c = 0; c < kNC; ++c) {
        t.s[c].assign(kNB, 0.);
        for (int c2 = 0; c2 < kNC; ++c2) t.C[c][c2].assign(kNB * kNB, 0.);
        t.Sall[c] = t.Vall[c] = t.Nrow[c] = 0;
    }
    std::vector<std::pair<int, double>> e[kNC], ag[kNC];
    auto flush = [&]() {
        bool any = false;
        for (int c = 0; c < kNC; ++c) {
            ag[c].clear();
            if (e[c].empty()) continue;
            any = true;
            std::sort(e[c].begin(), e[c].end());
            double tot = 0;
            size_t i = 0;
            while (i < e[c].size()) {
                const int b = e[c][i].first; double acc = 0;
                while (i < e[c].size() && e[c][i].first == b) { acc += e[c][i].second; ++i; }
                tot += acc;
                if (b >= 0) ag[c].push_back(std::make_pair(b, acc));
            }
            t.Vall[c] += tot * tot;
            e[c].clear();
        }
        if (!any) return;
        for (int c1 = 0; c1 < kNC; ++c1)
            for (const auto& a : ag[c1])
                for (int c2 = 0; c2 < kNC; ++c2)
                    for (const auto& b : ag[c2]) t.C[c1][c2][a.first * kNB + b.first] += a.second * b.second;
    };
    for (int set = 1; set <= 2; ++set) {
        if (!(opt.sets & set)) continue;
        const std::vector<Row>& v = d.tree(set, opt.hyp);
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const Row& r = v[i];
            if (r.grp != cur) { flush(); cur = r.grp; }
            if (!Owned(r, opt) || !sel(r)) continue;
            const double w = opt.db ? r.w_db : r.w;
            const double m = r.m4 * scale;
            const int ib = (m >= kLo && m < kHi) ? std::min(kNB - 1, (int)((m - kLo) / kBW)) : -1;
            for (int c = 0; c < kNC; ++c) {
                if (!comps[c].f(r)) continue;
                if (ib >= 0) t.s[c][ib] += w;
                t.Sall[c] += w; t.Nrow[c] += 1;
                e[c].push_back(std::make_pair(ib, w));
            }
        }
        flush();
    }
    for (int c = 0; c < kNC; ++c) {
        for (auto& x : t.s[c]) x *= kU;
        for (int c2 = 0; c2 < kNC; ++c2) for (auto& x : t.C[c][c2]) x *= kU * kU;
        t.Sall[c] *= kU; t.Vall[c] *= kU * kU;
    }
    // slices (sums and diagonal variances)
    Binning b; b.nbin = kNB; b.lo = kLo; b.hi = kHi;
    if (scale != 1.0) b.x = [scale](const Row& r) { return (double)r.m4 * scale; };
    double tot = 0; for (double x : t.s[cPeak]) tot += x;
    for (int k = 0; k < kNSl; ++k) {
        const float a = kSlLo[k], bb = kSlHi[k];
        const Sel ss = SelAnd(sel, [a, bb](const Row& r) { return r.mpipi >= a && r.mpipi < bb; });
        const std::vector<WSpec> w = FillMany(d, ss, comps, opt, b);
        double sp = 0;
        for (int c = 0; c < kNC; ++c) {
            t.sl[k][c] = w[c].s; t.slv[k][c] = w[c].v;
            for (auto& x : t.sl[k][c]) x *= kU;
            for (auto& x : t.slv[k][c]) x *= kU * kU;
        }
        for (double x : t.sl[k][cPeak]) sp += x;
        t.slPeakFrac[k] = tot > 0 ? sp / tot : 0;
    }
}

// ---------------------------------------------------------------- data
struct DSet {
    scan::Spectra sp;               // 0.30-1.10, 2.5 MeV, CB in 10 MeV bins
    std::vector<double> cov;        // covariance of N - CB per 2.5 MeV bin, kNB x kNB
    double diagRatio = 0;           // sum of group variances of N / sum of N (0.40-0.85)
    scan::Spectra slm[kNSl];        // slices on the ScanEngine grid (0.35-0.95), like-sign lists dropped
    scan::PeakFit pf; scan::Spread spr;
    double E = 0, eE = 0, h[kNSl], eh[kNSl];
};

void FillData(const scan::Data& dd, const scan::Sel& sel, DSet& D) {
    scan::FillOpt fo; fo.nbin = kNB; fo.lo = kLo; fo.hi = kHi; fo.reb = 4;
    D.sp = scan::Fill(dd, sel, fo);
    const scan::Spectra& sp = D.sp;
    D.cov.assign(kNB * kNB, 0.);
    auto groups = [&](int sample, int div, int nb, std::vector<double>& G) {
        std::vector<int> b;
        auto flush = [&]() {
            if (b.empty()) return;
            std::sort(b.begin(), b.end());
            std::vector<std::pair<int, double>> a;
            for (int x : b) { if (!a.empty() && a.back().first == x) a.back().second += 1; else a.push_back(std::make_pair(x, 1.)); }
            for (const auto& p : a) for (const auto& q : a) G[p.first * nb + q.first] += p.second * q.second;
            b.clear();
        };
        int cur = -1;
        for (const scan::Row& r : dd.s[sample]) {
            if (!r.cor_ok && !dd.fallback) continue;
            if (r.grp != cur) { flush(); cur = r.grp; }
            const double m = r.m4;
            if (!(m == m)) continue;
            const int ib = sp.bin(m);
            if (ib < 0 || ib >= sp.nbin) continue;
            if (!sel(r)) continue;
            b.push_back(ib / div);
        }
        flush();
    };
    groups(0, 1, kNB, D.cov);                       // N: event groups
    {
        double sv = 0, sn = 0;
        for (int i = 0; i < kNB; ++i) if (sp.x(i) > 0.40 && sp.x(i) < 0.85) { sv += D.cov[i * kNB + i]; sn += sp.n[0][i]; }
        D.diagRatio = sn > 0 ? sv / sn : 0;
    }
    // CB: like-sign event groups per wide bin, propagated to the interpolated fine bins
    const int nw = sp.nw, reb = sp.reb;
    for (int q = 0; q < 2; ++q) {
        std::vector<double> G(nw * nw, 0.);
        groups(1 + q, reb, nw, G);
        std::vector<double> Dq(nw, 1.);
        for (int j = 0; j < nw; ++j) {
            const double own = sp.w[1 + q][j], oth = sp.w[2 - q][j];
            Dq[j] = (own > 0 && oth > 0) ? std::sqrt(oth / own) : 1.0;
        }
        std::vector<int> j0(kNB), j1(kNB);
        std::vector<double> a0(kNB), a1(kNB);
        for (int i = 0; i < kNB; ++i) {
            double f; sp.Interp(i, j0[i], j1[i], f);
            a0[i] = (1 - f) / reb * Dq[j0[i]]; a1[i] = f / reb * Dq[j1[i]];
        }
        for (int i = 0; i < kNB; ++i)
            for (int k = 0; k < kNB; ++k) {
                const int ji[2] = {j0[i], j1[i]}, jk[2] = {j0[k], j1[k]};
                const double ai[2] = {a0[i], a1[i]}, ak[2] = {a0[k], a1[k]};
                double s = 0;
                for (int p = 0; p < 2; ++p) for (int qq = 0; qq < 2; ++qq) s += ai[p] * ak[qq] * G[ji[p] * nw + jk[qq]];
                D.cov[i * kNB + k] += s;
            }
    }
    D.pf = scan::FitPeakMain(sp);
    D.spr = scan::PeakModelSpread(sp);
    std::vector<scan::Sel> ss;
    for (int k = 0; k < kNSl; ++k) {
        const float a = kSlLo[k], bb = kSlHi[k];
        ss.push_back([sel, a, bb](const scan::Row& r) { return r.mpipi >= a && r.mpipi < bb && sel(r); });
    }
    const std::vector<scan::Spectra> sv = scan::FillMany(dd, ss);   // default grid, as the data study
    double v = 0;
    D.E = 0;
    for (int k = 0; k < kNSl; ++k) {
        const scan::StepFit f = scan::FitEdgeStep(sv[k], SliceCfg(k));
        D.h[k] = f.h; D.eh[k] = f.eh;
        D.E += f.h; v += f.eh * f.eh;
        D.slm[k] = sv[k]; D.slm[k].ls[0].clear(); D.slm[k].ls[1].clear(); D.slm[k].win.clear();
    }
    D.eE = std::sqrt(v);
}

// ---------------------------------------------------------------- fit
struct FitIn {
    int reb = 1, npts = 0;
    double lo = 0, hi = 0;
    std::vector<int> pt;                     // coarse bin index of each point
    std::vector<double> x, y, Vd;            // centre, N - CB, data covariance
    std::vector<double> T[kNC];
    std::vector<double> CT[kNC][kNC];
    std::vector<double> G;                   // unit Gaussian per point
};

FitIn MakeIn(const DSet& D, const TplSet& t, int reb, double lo, double hi, double gm, double gs) {
    FitIn in; in.reb = reb; in.lo = lo; in.hi = hi;
    const int nc = kNB / reb;
    const double bw = kBW * reb;
    for (int I = 0; I < nc; ++I) {
        const double xc = kLo + (I + 0.5) * bw;
        if (xc >= lo && xc < hi) in.pt.push_back(I);
    }
    const int n = in.npts = (int)in.pt.size();
    in.x.resize(n); in.y.assign(n, 0.); in.Vd.assign(n * n, 0.); in.G.resize(n);
    for (int c = 0; c < kNC; ++c) {
        in.T[c].assign(n, 0.);
        for (int c2 = 0; c2 < kNC; ++c2) in.CT[c][c2].assign(n * n, 0.);
    }
    for (int q = 0; q < n; ++q) {
        const int I = in.pt[q];
        const double x1 = kLo + I * bw;
        in.x[q] = x1 + 0.5 * bw;
        in.G[q] = scan::GaussFrac(x1, x1 + bw, gm, gs);
        for (int a = 0; a < reb; ++a) {
            const int i = I * reb + a;
            in.y[q] += D.sp.sig[i];
            for (int c = 0; c < kNC; ++c) in.T[c][q] += t.s[c][i];
        }
        for (int p = 0; p < n; ++p) {
            const int J = in.pt[p];
            for (int a = 0; a < reb; ++a)
                for (int b = 0; b < reb; ++b) {
                    const int i = I * reb + a, k = J * reb + b;
                    in.Vd[q * n + p] += D.cov[i * kNB + k];
                    for (int c = 0; c < kNC; ++c)
                        for (int c2 = 0; c2 < kNC; ++c2) in.CT[c][c2][q * n + p] += t.C[c][c2][i * kNB + k];
                }
        }
    }
    return in;
}

struct FitOut {
    bool ok = false, conv = false;
    int np = 0, npts = 0, ndf = 0, iter = 0;
    double chi2 = 0, prob = 0, tvShare = 0;
    std::vector<double> p, e, cov;          // parameters, errors, covariance (np x np)
    Coef q;                                 // template coefficients
    std::vector<double> fit, sig, pull;     // per point: model, sqrt(V_qq), pull
    // implied contents: value, error from the parameters, error from the template statistics
    double S = 0, eS = 0, tS = 0, Speak = 0, E = 0, eE = 0, tE = 0, Ebump = 0;
    double bumpAll = 0, eBumpAll = 0, peakAll = 0, ePeakAll = 0;
};

FitOut DoFit(const FitIn& in, const Model& m) {
    FitOut o;
    const int n = in.npts, nt = (int)m.cols.size(), np = nt + (m.tilt ? 1 : 0) + (m.gauss ? 1 : 0);
    o.np = np; o.npts = n; o.ndf = n - np;
    std::vector<double> F(n * np, 0.);
    for (int q = 0; q < n; ++q) {
        for (int j = 0; j < nt; ++j) for (int k : m.cols[j]) F[q * np + j] += in.T[k][q];
        int j = nt;
        if (m.tilt) F[q * np + j++] = (in.T[cVM][q] + in.T[cRest][q]) * TiltT(in.x[q]);
        if (m.gauss) F[q * np + j++] = in.G[q];
    }
    std::vector<double> p(np, 1.);
    for (int j = nt; j < np; ++j) p[j] = 0.;
    std::vector<double> V(n * n), cq(kNC * n);
    auto buildV = [&](const std::vector<double>& pp) {
        const Coef c = CoefOf(m, pp);
        for (int k = 0; k < kNC; ++k) for (int q = 0; q < n; ++q) cq[k * n + q] = c.at(k, in.x[q]);
        for (int q = 0; q < n; ++q)
            for (int r = 0; r < n; ++r) {
                const int z = q * n + r;
                double v = in.Vd[z];
                for (int k = 0; k < kNC; ++k) {
                    const double a = cq[k * n + q];
                    if (a == 0) continue;
                    for (int l = 0; l < kNC; ++l) { const double b = cq[l * n + r]; if (b != 0) v += a * b * in.CT[k][l][z]; }
                }
                V[z] = v;
            }
    };
    scan::LinRes r;
    for (int it = 0; it < 200; ++it) {
        buildV(p);
        std::vector<double> Vc = V;
        r = scan::LinFitCov(n, np, in.y.data(), Vc, F.data());
        o.iter = it + 1;
        if (!r.ok) return o;
        double dmax = 0;
        for (int j = 0; j < np; ++j) dmax = std::max(dmax, std::fabs(r.par[j] - p[j]) / std::sqrt(r.cov[j][j]));
        for (int j = 0; j < np; ++j) p[j] = r.par[j];
        if (dmax < 1e-6) { o.conv = true; break; }
    }
    o.ok = true;
    o.p = p; o.e.resize(np); o.cov.resize(np * np);
    for (int j = 0; j < np; ++j) { o.e[j] = std::sqrt(r.cov[j][j]); for (int l = 0; l < np; ++l) o.cov[j * np + l] = r.cov[j][l]; }
    o.q = CoefOf(m, p);
    buildV(p);
    {
        std::vector<double> Vc = V;
        o.chi2 = scan::LinFitCov(n, np, in.y.data(), Vc, F.data()).chi2;
    }
    o.prob = TMath::Prob(o.chi2, o.ndf);
    o.fit.resize(n); o.sig.resize(n); o.pull.resize(n);
    double sh = 0;
    for (int q = 0; q < n; ++q) {
        double f = 0;
        for (int j = 0; j < np; ++j) f += p[j] * F[q * np + j];
        o.fit[q] = f; o.sig[q] = std::sqrt(V[q * n + q]);
        o.pull[q] = (in.y[q] - f) / o.sig[q];
        sh += (V[q * n + q] - in.Vd[q * n + q]) / V[q * n + q];
    }
    o.tvShare = sh / n;
    return o;
}

// model on the fine grid (data rows per 2.5 MeV); slice >= 0: the M(pi+pi-) slice templates
std::vector<double> ModelFine(const TplSet& t, const Coef& c, double gm, double gs, int slice = -1) {
    std::vector<double> y(kNB, 0.);
    for (int i = 0; i < kNB; ++i) {
        const double x1 = kLo + i * kBW, xc = x1 + 0.5 * kBW;
        for (int k = 0; k < kNC; ++k) y[i] += c.at(k, xc) * (slice < 0 ? t.s[k][i] : t.sl[slice][k][i]);
        if (c.AG != 0) y[i] += c.AG * (slice < 0 ? 1. : t.slPeakFrac[slice]) * scan::GaussFrac(x1, x1 + kBW, gm, gs);
    }
    return y;
}
// P_main of the data study applied to a model (weights as for the data: model + data CB, CB covariance)
double SOf(const std::vector<double>& y, const DSet& D) {
    const scan::PeakFit f = scan::FitPeakArrays(y.data(), D.sp.cb.data(), &D.sp, kNB, kLo, kBW, scan::PeakCfg(), false);
    return f.ok ? f.S : std::nan("");
}
// step height of slice k for a fine-grid model of that slice
double HOf(const std::vector<double>& y, const DSet& D, int k) {
    scan::Spectra ms = D.slm[k];
    for (int i = 0; i < ms.nbin; ++i) ms.sig[i] = y[i + kSlOff];
    const scan::StepFit f = scan::FitEdgeStep(ms, SliceCfg(k));
    return f.ok ? f.h : std::nan("");
}
double EOf(const TplSet& t, const Coef& c, double gm, double gs, const DSet& D) {
    double E = 0;
    for (int k = 0; k < kNSl; ++k) E += HOf(ModelFine(t, c, gm, gs, k), D, k);
    return E;
}

void Implied(FitOut& o, const Model& m, const TplSet& t, const DSet& D, double gm, double gs) {
    if (!o.ok) return;
    auto evalS = [&](const std::vector<double>& p) { return SOf(ModelFine(t, CoefOf(m, p), gm, gs), D); };
    auto evalE = [&](const std::vector<double>& p) { return EOf(t, CoefOf(m, p), gm, gs, D); };
    o.S = evalS(o.p); o.E = evalE(o.p);
    const int np = o.np;
    std::vector<double> gS(np), gE(np);
    for (int j = 0; j < np; ++j) {
        const double d = 0.2 * o.e[j] + 1e-9;
        std::vector<double> a = o.p, b = o.p;
        a[j] += d; b[j] -= d;
        gS[j] = (evalS(a) - evalS(b)) / (2 * d);
        gE[j] = (evalE(a) - evalE(b)) / (2 * d);
    }
    double vS = 0, vE = 0;
    for (int j = 0; j < np; ++j) for (int l = 0; l < np; ++l) { vS += gS[j] * gS[l] * o.cov[j * np + l]; vE += gE[j] * gE[l] * o.cov[j * np + l]; }
    o.eS = std::sqrt(std::max(0., vS)); o.eE = std::sqrt(std::max(0., vE));
    // template statistics at fixed parameters: gradient of S (E) with respect to the model in each bin,
    // variance = g^T Cov(model) g, Cov(model)_ij = sum_kl c_k(x_i) c_l(x_j) C_kl,ij (slices: diagonal only)
    {
        const std::vector<double> y0 = ModelFine(t, o.q, gm, gs);
        std::vector<int> idx; std::vector<double> g;
        for (int i = 0; i < kNB; ++i) {
            const double xc = kLo + (i + 0.5) * kBW;
            if (xc < 0.515 || xc > 0.615) continue;
            std::vector<double> y = y0; const double dl = std::max(1.0, 0.01 * std::fabs(y0[i]));
            y[i] += dl;
            idx.push_back(i); g.push_back((SOf(y, D) - o.S) / dl);
        }
        double v = 0;
        for (size_t a = 0; a < idx.size(); ++a)
            for (size_t b = 0; b < idx.size(); ++b) {
                const int i = idx[a], j = idx[b];
                const double xi = kLo + (i + 0.5) * kBW, xj = kLo + (j + 0.5) * kBW;
                double c = 0;
                for (int k = 0; k < kNC; ++k) for (int l = 0; l < kNC; ++l) c += o.q.at(k, xi) * o.q.at(l, xj) * t.C[k][l][i * kNB + j];
                v += g[a] * g[b] * c;
            }
        o.tS = std::sqrt(std::max(0., v));
        double vE2 = 0;
        for (int k = 0; k < kNSl; ++k) {
            const std::vector<double> ys = ModelFine(t, o.q, gm, gs, k);
            const double h0 = HOf(ys, D, k);
            for (int i = kSlOff; i < kSlOff + 240; ++i) {
                const double xc = kLo + (i + 0.5) * kBW;
                if (xc < 0.42 || xc > 0.62) continue;
                double var = 0;
                for (int c = 0; c < kNC; ++c) var += std::pow(o.q.at(c, xc), 2) * t.slv[k][c][i];
                if (var <= 0) continue;
                std::vector<double> y = ys; const double dl = std::max(1.0, 0.01 * std::fabs(ys[i]));
                y[i] += dl;
                const double gi = (HOf(y, D, k) - h0) / dl;
                vE2 += gi * gi * var;
            }
        }
        o.tE = std::sqrt(vE2);
    }
    // parts: S of the peak template + Gaussian, E of the bump template
    {
        Coef c0 = o.q; c0.c[cPeak] = 0; c0.AG = 0;
        o.Speak = o.S - SOf(ModelFine(t, c0, gm, gs), D);
        Coef c1 = o.q; c1.c[cBump] = 0;
        o.Ebump = o.E - EOf(t, c1, gm, gs, D);
    }
    double ba = 0, pa = 0;
    for (int i = 0; i < kNB; ++i) { ba += t.s[cBump][i]; pa += t.s[cPeak][i]; }
    auto colOf = [&](int comp) { for (size_t j = 0; j < m.cols.size(); ++j) for (int k : m.cols[j]) if (k == comp) return (int)j; return -1; };
    const int jb = colOf(cBump), jp = colOf(cPeak);
    o.bumpAll = o.q.c[cBump] * ba; o.peakAll = o.q.c[cPeak] * pa;
    o.eBumpAll = jb >= 0 ? o.e[jb] * ba : 0; o.ePeakAll = jp >= 0 ? o.e[jp] * pa : 0;
}

// ---------------------------------------------------------------- plotting helpers
void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy, double lm, double tm, double bm = 0.13,
            double rm = 0.02) {
    TPad* p = new TPad(name, "", x1, y1, x2, y2);
    p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
    p->SetLogy(logy); p->Draw(); p->cd();
    return p;
}
void Axes(TH1D* h, const char* xt, const char* yt, double ts = 0.055, double yoff = 1.2) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 0.88); h->GetYaxis()->SetLabelSize(ts * 0.88);
    h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(yoff);
    h->GetYaxis()->SetMaxDigits(4);
}
void Line(TH1D* h, int col, int lw, int ls, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerSize(0); h->SetMarkerColor(col);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); } else h->SetFillStyle(0);
}
void Mk(TH1D* h, int col, int st, double sz = 1.0) {
    h->SetMarkerStyle(st); h->SetMarkerSize(sz); h->SetMarkerColor(col); h->SetLineColor(col); h->SetLineWidth(1);
}
const char* kXt = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const char* kYt = "d#sigma/dM [nb/(GeV/c^{2})]";

// fine-grid array (data rows per 2.5 MeV) -> histogram in nb/(GeV/c2), rebinned by reb; cov: kNB x kNB (optional)
TH1D* HFine(const char* name, const std::vector<double>& y, const std::vector<double>* cov, int reb) {
    const int nc = kNB / reb;
    TH1D* h = new TH1D(name, "", nc, kLo, kHi);
    h->SetDirectory(nullptr); h->Sumw2();
    const double f = 1.0 / (kLdata * kBW * reb);
    for (int I = 0; I < nc; ++I) {
        double s = 0, v = 0;
        for (int a = 0; a < reb; ++a) {
            s += y[I * reb + a];
            if (cov) for (int b = 0; b < reb; ++b) v += (*cov)[(I * reb + a) * kNB + I * reb + b];
        }
        h->SetBinContent(I + 1, s * f); h->SetBinError(I + 1, std::sqrt(std::max(0., v)) * f);
    }
    return h;
}
std::string ParStr(const Model& m, const FitOut& r, bool latex) {
    std::string s;
    for (int j = 0; j < r.np; ++j) {
        std::string n = m.pname[j];
        if (latex) {
            if (n == "A_G") n = "A_{G}";
            else if (n == "b_tilt") n = "#beta";
            else n = "#alpha_{" + n.substr(2) + "}";
            if (n == "#alpha_{eta}") n = "#alpha_{#eta}";
        }
        s += Form(latex ? "%s %s  " : "%s=%s ", n.c_str(), m.pname[j] == "A_G" ? Form(latex ? "%.0f#pm%.0f" : "%.0f(%.0f)", r.p[j], r.e[j])
                                                                         : Form(latex ? "%.2f#pm%.2f" : "%.3f(%.3f)", r.p[j], r.e[j]));
    }
    return s;
}

}  // namespace stf

// ============================================================================
void sim_templatefit() {
    using namespace simeng;
    using namespace stf;
    TStopwatch swAll;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
    const std::string tab = std::string(kTabDir) + "/";
    const std::string fig = std::string(kFigDir) + "/";
    const std::vector<Model> models = Models();
    const int nM = (int)models.size();
    FILE* fchk = fopen((tab + "fit_checks.txt").c_str(), "w");
    auto both = [&](FILE* f, const char* s) { printf("%s", s); fflush(stdout); if (f) { fprintf(f, "%s", s); fflush(f); } };

    // (1) event groups of the simulation: one isBest == 1 row per group; grp changes iff (event, vz) changes
    {
        TFile* f = TFile::Open(kCompact, "READ");
        for (int k : {kStdUS, kLepUS}) {
            TTree* t = (TTree*)f->Get(TreeName(k));
            Int_t grp, isBest, event; Float_t vz;
            t->SetBranchStatus("*", 0);
            for (const char* b : {"grp", "isBest", "event", "vz"}) t->SetBranchStatus(b, 1);
            t->SetBranchAddress("grp", &grp); t->SetBranchAddress("isBest", &isBest);
            t->SetBranchAddress("event", &event); t->SetBranchAddress("vz", &vz);
            const Long64_t n = t->GetEntries();
            long ngrp = 0, bad = 0, keyMismatch = 0; int cur = -1, nb = 0, pe = -1; float pv = 0;
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                const bool newKey = (i == 0) || event != pe || vz != pv;
                const bool newGrp = grp != cur;
                if (newKey != newGrp) keyMismatch++;
                if (newGrp) { if (i > 0 && nb != 1) bad++; ngrp++; nb = 0; cur = grp; }
                if (isBest == 1) nb++;
                pe = event; pv = vz;
            }
            if (nb != 1) bad++;
            both(fchk, Form("sim groups %s: rows %lld, groups %ld, groups without exactly one isBest==1 row: %ld, "
                            "rows where a new (event, vz) key and a new grp disagree: %ld\n", TreeName(k), n, ngrp, bad, keyMismatch));
            t->ResetBranchAddresses();
        }
        f->Close(); delete f;
    }

    // ---------------- load ----------------
    SimData d;
    if (!Load(d, kCompact, (1 << kStdUS) | (1 << kLepUS), kLo, kHi)) return;
    scan::Data dd;
    if (!scan::Load(dd, scan::kCompact, kLo, kHi)) return;

    // (2) reproduction of the data-study numbers (no vz cut, isBest >= 0, cor_ok == 1)
    {
        const scan::Spectra sp = scan::Fill(dd, scan::SelNone());
        const scan::PeakFit pf = scan::FitPeakMain(sp);
        const scan::Spread spr = scan::PeakModelSpread(sp);
        double E = 0, v = 0;
        std::vector<scan::Sel> ss;
        for (int k = 0; k < kNSl; ++k) {
            const float a = kSlLo[k], b = kSlHi[k];
            ss.push_back([a, b](const scan::Row& r) { return r.mpipi >= a && r.mpipi < b; });
        }
        const std::vector<scan::Spectra> sv = scan::FillMany(dd, ss);
        for (int k = 0; k < kNSl; ++k) { const scan::StepFit f = scan::FitEdgeStep(sv[k], SliceCfg(k)); E += f.h; v += f.eh * f.eh; }
        both(fchk, Form("data study reproduction (no vz cut, data loaded in 0.30-1.10): S = %.0f +- %.0f +- %.0f (reported 4657 +- 390 +- 1233), "
                        "E = %.0f +- %.0f (reported 1241 +- 76)\n", pf.S, pf.eS, spr.half, E, std::sqrt(v)));
    }

    // ---------------- Gaussian shape and position of the simulated conversion peak (pre, mode A, w) ----------------
    double gm = scan::kPeakMean, gs = 0, gmSim = 0, egm = 0, egs = 0;
    {
        TplSet t0;
        FillTpl(d, SelPreset(kPre), OptOf(0), 1.0, t0);
        TH1D* h = new TH1D("hPeakTpl", "", kNB, kLo, kHi);
        h->SetDirectory(nullptr);
        for (int i = 0; i < kNB; ++i) { h->SetBinContent(i + 1, t0.s[cPeak][i]); h->SetBinError(i + 1, std::sqrt(t0.C[cPeak][cPeak][i * kNB + i])); }
        TF1* g = new TF1("gPeak", "gaus", 0.53, 0.56);
        g->SetParameters(h->GetMaximum(), 0.5445, 0.006);
        h->Fit(g, "QNR0");
        for (int it = 0; it < 3; ++it) {
            g->SetRange(g->GetParameter(1) - 2 * g->GetParameter(2), g->GetParameter(1) + 2 * g->GetParameter(2));
            h->Fit(g, "QNR0");
        }
        gmSim = g->GetParameter(1); gs = std::fabs(g->GetParameter(2)); egm = g->GetParError(1); egs = g->GetParError(2);
        double tot = 0, core = 0;
        for (int i = 0; i < kNB; ++i) { tot += t0.s[cPeak][i]; const double x = kLo + (i + 0.5) * kBW; if (std::fabs(x - gmSim) < 2 * gs) core += t0.s[cPeak][i]; }
        both(fchk, Form("Gaussian of the peak-conversion template (pre, mode A, w; +-2 sigma iterated fit): mean %.5f +- %.5f, sigma %.5f +- %.5f GeV/c2, "
                        "chi2/ndf %.1f/%d; share within +-2 sigma %.3f; data peak 0.5439, 5.4 MeV.\n  -> G of the fit at %.4f (data), sigma %.5f; "
                        "option A_w_sc scales the simulated M4 by %.5f (K0S M(pipi) scale data/sim; not this peak offset, review finding 8)\n", gmSim, egm, gs, egs, g->GetChisquare(), g->GetNDF(), core / tot, gm, gs, gSimScale));
    }

    // ---------------- templates and data per preset ----------------
    std::vector<TplSet> T(kNP * kNO);
    std::vector<DSet> D(kNP);
    TStopwatch sw;
    for (int ip = 0; ip < kNP; ++ip) {
        FillData(dd, DataSelPreset(kPresetId[ip]), D[ip]);
        for (int io = 0; io < kNO; ++io) FillTpl(d, SelPreset(kPresetId[ip]), OptOf(io), ScaleOf(io), T[ip * kNO + io]);
        printf("preset %s filled (%.1f s)\n", kPName[ip], sw.RealTime()); sw.Continue(); fflush(stdout);
    }
    // (3) template sums and variances against the engine (diagonal of the group covariance = WSpec::v)
    {
        Binning b; b.nbin = kNB; b.lo = kLo; b.hi = kHi;
        double worst = 0;
        for (int io = 0; io < 3; ++io) {
            const std::vector<WSpec> w = FillMany(d, SelPreset(kPre), FitComps(), OptOf(io), b);
            const TplSet& t = T[0 * kNO + io];
            for (int c = 0; c < kNC; ++c)
                for (int i = 0; i < kNB; ++i) {
                    const double a = t.s[c][i], bb = w[c].s[i] * kU, va = t.C[c][c][i * kNB + i], vb = w[c].v[i] * kU * kU;
                    if (bb > 0) worst = std::max(worst, std::fabs(a / bb - 1));
                    if (vb > 0) worst = std::max(worst, std::fabs(va / vb - 1));
                }
        }
        both(fchk, Form("template sums and diagonal group variances against simeng::FillMany (pre, options A_w, B_w, A_wdb, 4 templates, 320 bins): "
                        "largest relative difference %.2e\n", worst));
        for (int ip = 0; ip < kNP; ++ip)
            both(fchk, Form("data %s: sum of group variances of N / sum of N in 0.40-0.85 = %.4f (Poisson would be 1)\n", kPName[ip], D[ip].diagRatio));
        for (int ip = 0; ip < kNP; ++ip) {
            const std::vector<WSum> s = SumMany(d, SelPreset(kPresetId[ip]), {CompBumpDalitz(), CompBumpConv()}, SimOpt());
            both(fchk, Form("bump template %s (A, w): pi0 Dalitz share %.3f, conversion share %.3f (fixed in the fit)\n", kPName[ip],
                            s[0].S / (s[0].S + s[1].S), s[1].S / (s[0].S + s[1].S)));
        }
        for (int ip = 0; ip < kNP; ++ip)
            both(fchk, Form("data %s: S(P_main) = %.0f +- %.0f +- %.0f (model spread), E = %.0f +- %.0f; h per slice %.0f %.0f %.0f %.0f%s\n", kPName[ip],
                            D[ip].pf.S, D[ip].pf.eS, D[ip].spr.half, D[ip].E, D[ip].eE, D[ip].h[0], D[ip].h[1], D[ip].h[2], D[ip].h[3],
                            kPresetId[ip] == kSlides ? " (E not valid: eta-rest cut, data study)" : ""));
    }

    // ---------------- template table ----------------
    {
        FILE* f = fopen((tab + "fit_templates.txt").c_str(), "w");
        for (int pass = 0; pass < 2; ++pass) {
            FILE* o = pass ? f : stdout;
            if (!o) continue;
            fprintf(o, "# templates of the fit, e+e-, std + lep owned rows, in expected data rows (sum w * L_data / L_sim, L_data = %.1f nb^-1 UNCONFIRMED, L_sim = %.2f)\n", kLdata, kLsim);
            fprintf(o, "# N_eff = (sum w)^2 / V, V from event groups (rows of one event in different bins included in the range sums)\n");
            fprintf(o, "# eta -> pi+pi-e+e- (direct decay) is NOT in any template (eta set missing). A_w_sc: simulated M4 x %.5f (K0S M(pipi) scale data/sim)\n", gSimScale);
            fprintf(o, "%-13s %-7s %-8s %10s %12s %10s | %12s %10s %8s\n", "preset", "option", "template", "rows", "rows(data)", "N_eff",
                    "0.43-0.75", "+-", "N_eff");
            for (int ip = 0; ip < kNP; ++ip)
                for (int io = 0; io < kNO; ++io) {
                    const TplSet& t = T[ip * kNO + io];
                    std::vector<int> idx;
                    for (int i = 0; i < kNB; ++i) { const double x = kLo + (i + 0.5) * kBW; if (x >= 0.43 && x < 0.75) idx.push_back(i); }
                    for (int c = 0; c < kNC; ++c) {
                        double s = 0, v = 0;
                        for (int i : idx) { s += t.s[c][i]; for (int k : idx) v += t.C[c][c][i * kNB + k]; }
                        fprintf(o, "%-13s %-7s %-8s %10.0f %12.1f %10.1f | %12.1f %10.1f %8.1f\n", kPName[ip], kOTag[io], kCTag[c], t.Nrow[c], t.Sall[c],
                                t.Vall[c] > 0 ? t.Sall[c] * t.Sall[c] / t.Vall[c] : 0., s, std::sqrt(v), v > 0 ? s * s / v : 0.);
                    }
                }
            for (int ip = 0; ip < kNP; ++ip) {
                const TplSet& t = T[ip * kNO];
                fprintf(o, "# %s (A, w): share of the peak template in the M(pipi) slices 0.30-0.34/0.34-0.37/0.37-0.39/0.39-0.40: %.3f %.3f %.3f %.3f\n",
                        kPName[ip], t.slPeakFrac[0], t.slPeakFrac[1], t.slPeakFrac[2], t.slPeakFrac[3]);
            }
        }
        if (f) fclose(f);
    }

    // ---------------- decomposition at alpha = 1 (pure SMASH, mode A, w) ----------------
    // which simulated components carry S (P_main) and E (edge steps): contribution of a component =
    // measure(sum of all) - measure(sum without it); also rows in the bump and peak windows
    {
        std::vector<Comp> dc = {CompBumpDalitz(), CompBumpConv(), CompPeakConv()};
        for (int l : {1, 3, 4, 5, 6}) dc.push_back({Form("T2 %s", LclsName(l)), [l](const Row& r) { return r.tmpl == 2 && r.lcls == l; }});
        for (int p : {0, 1, 2, 3, 4, 5, 6, 7, 8})
            dc.push_back({Form("T2 %s", PclsName(p)), [p](const Row& r) { return r.tmpl == 2 && r.pcls == p; }});
        dc.push_back({"T2 std-set rows", [](const Row& r) { return r.tmpl == 2 && r.set == 1; }});
        dc.push_back({"T2 lep-set rows", [](const Row& r) { return r.tmpl == 2 && r.set == 2; }});
        const int nd = (int)dc.size();
        const int nLcls = 3 + 5;   // first block: signals + T2 by lepton class (a partition of T1 + T2 except other lcls)
        FILE* f = fopen((tab + "fit_decomposition.txt").c_str(), "w");
        Binning b; b.nbin = kNB; b.lo = kLo; b.hi = kHi;
        for (int ip = 0; ip < kNP; ++ip) {
            const Sel sel = SelPreset(kPresetId[ip]);
            const std::vector<WSpec> w = FillMany(d, sel, dc, SimOpt(), b);
            std::vector<std::vector<WSpec>> ws(kNSl);
            for (int k = 0; k < kNSl; ++k) {
                const float a = kSlLo[k], bb = kSlHi[k];
                ws[k] = FillMany(d, SelAnd(sel, [a, bb](const Row& r) { return r.mpipi >= a && r.mpipi < bb; }), dc, SimOpt(), b);
            }
            // total = T1 + T2 (from the fit templates, alpha = 1)
            const TplSet& t = T[ip * kNO];
            Coef one; for (int k = 0; k < kNC; ++k) one.c[k] = 1.;
            const std::vector<double> yAll = ModelFine(t, one, gm, gs);
            std::vector<std::vector<double>> ysl(kNSl);
            for (int k = 0; k < kNSl; ++k) ysl[k] = ModelFine(t, one, gm, gs, k);
            const double S0 = SOf(yAll, D[ip]);
            double E0 = 0; for (int k = 0; k < kNSl; ++k) E0 += HOf(ysl[k], D[ip], k);
            double bw0 = 0, pw0 = 0;
            for (int i = 0; i < kNB; ++i) {
                const double x = kLo + (i + 0.5) * kBW;
                if (x >= scan::kBumpLo && x < scan::kBumpHi) bw0 += yAll[i];
                if (x >= scan::kPeakLo && x < scan::kPeakHi) pw0 += yAll[i];
            }
            for (int pass = 0; pass < 2; ++pass) {
                FILE* o = pass ? f : stdout;
                if (!o) continue;
                fprintf(o, "%s## %s, mode A, w, alpha = 1 (SMASH, L_data = %.1f unconfirmed): S(P_main) of T1+T2 = %.0f (data %.0f +- %.0f), E = %.0f (data %.0f +- %.0f);"
                           " rows in bump window %.0f, peak window %.0f (data N - CB: %.0f, %.0f)\n", ip ? "\n" : "", kPName[ip], kLdata, S0, D[ip].pf.S, D[ip].pf.eS, E0,
                        D[ip].E, D[ip].eE, bw0, pw0, D[ip].sp.WindowApprox(scan::kBumpLo, scan::kBumpHi).sig(), D[ip].sp.WindowApprox(scan::kPeakLo, scan::kPeakHi).sig());
                fprintf(o, "%-40s %12s %10s %12s %12s %10s %10s\n", "component (removed from T1+T2)", "rows 0.3-1.1", "N_eff", "bump window", "peak window", "dS", "dE");
            }
            for (int j = 0; j < nd; ++j) {
                std::vector<double> y = yAll;
                double bwj = 0, pwj = 0;
                for (int i = 0; i < kNB; ++i) {
                    const double v = w[j].s[i] * kU, x = kLo + (i + 0.5) * kBW;
                    y[i] -= v;
                    if (x >= scan::kBumpLo && x < scan::kBumpHi) bwj += v;
                    if (x >= scan::kPeakLo && x < scan::kPeakHi) pwj += v;
                }
                const double dS = S0 - SOf(y, D[ip]);
                double dE = E0;
                for (int k = 0; k < kNSl; ++k) {
                    std::vector<double> ys = ysl[k];
                    for (int i = 0; i < kNB; ++i) ys[i] -= ws[k][j].s[i] * kU;
                    dE -= HOf(ys, D[ip], k);
                }
                for (int pass = 0; pass < 2; ++pass) {
                    FILE* o = pass ? f : stdout;
                    if (!o) continue;
                    if (j == nLcls) fprintf(o, "-- T2 by pion-pair class --\n");
                    if (j == nd - 2) fprintf(o, "-- T2 by set --\n");
                    fprintf(o, "%-40s %12.0f %10.0f %12.0f %12.0f %10.0f %10.0f\n", dc[j].name.c_str(), w[j].S * kU, w[j].NeffIn(), bwj, pwj, dS, dE);
                }
            }
        }
        if (f) fclose(f);
    }

    // ---------------- all fits ----------------
    auto IDX =[&](int ip, int io, int ib, int ir, int im) { return (((ip * kNO + io) * kNBW + ib) * kNR + ir) * nM + im; };
    std::vector<FitOut> R(kNP * kNO * kNBW * kNR * nM);
    std::vector<FitIn> nominalIn(kNP);
    sw.Start();
    for (int ip = 0; ip < kNP; ++ip) {
        for (int io = 0; io < kNO; ++io)
            for (int ib = 0; ib < kNBW; ++ib)
                for (int ir = 0; ir < kNR; ++ir) {
                    const FitIn in = MakeIn(D[ip], T[ip * kNO + io], kReb[ib], kRLo[ir], kRHi[ir], gm, gs);
                    if (io == 0 && ib == 0 && ir == 0) nominalIn[ip] = in;
                    for (int im = 0; im < nM; ++im) {
                        FitOut& o = R[IDX(ip, io, ib, ir, im)];
                        o = DoFit(in, models[im]);
                        Implied(o, models[im], T[ip * kNO + io], D[ip], gm, gs);
                    }
                }
        printf("fits of %s done (%.1f s)\n", kPName[ip], sw.RealTime()); sw.Continue(); fflush(stdout);
    }

    // ---------------- results table ----------------
    {
        FILE* f = fopen((tab + "fit_results.txt").c_str(), "w");
        for (int pass = 0; pass < 2; ++pass) {
            FILE* o = pass ? f : stdout;
            if (!o) continue;
            fprintf(o, "# PRELIMINARY template fit of data N - CB in M4 (e+e-); templates in data rows with L_data = %.1f nb^-1 (UNCONFIRMED): alpha = 1 = SMASH\n", kLdata);
            fprintf(o, "# eta -> pi+pi-e+e- (direct decay) is MISSING in the simulation; a_peak and A_G absorb it\n");
            fprintf(o, "# M1 a_bump*bump + a_peak*peak + a_bkg*bkg; M2 M1 + A_G*G; M3 a_eta*(bump+peak) + a_bkg*bkg + A_G*G; M4 a_bump, a_peak, a_VM*bkg(omega/rho0 pions)\n");
            fprintf(o, "#  + a_rest*bkg(other pions); M5 a_eta*(bump+peak) + a_VM + a_rest + A_G*G; M6 M1 with bkg coefficient a_bkg + b_tilt*(M-0.6)/0.1; M7 M6 with a_peak = a_bump + A_G*G\n");
            fprintf(o, "# G: mean %.4f (data peak), sigma %.5f GeV/c2 (peak-conversion template, its mean %.5f); A_G = whole Gaussian in data rows\n", gm, gs, gmSim);
            fprintf(o, "# options: A_w mode A, w (PDG); B_w mode B; A_wdb w_db; A_w_sc mode A, w, simulated M4 x %.5f (K0S M(pipi) scale data/sim)\n", gSimScale);
            fprintf(o, "# chi2: data covariance (event groups, CB propagation) + template covariance (event groups) with the fitted coefficients; tv = mean template share of V\n");
            fprintf(o, "# S = P_main (data-study peak fit) on the fitted model; Speak = part from peak template + G; E = edge-step sum on the model; Ebump = part from the bump\n");
            fprintf(o, "#   errors: first = from the parameters, second = template statistics at fixed parameters\n");
            fprintf(o, "%-12s %-6s %-3s %-9s %-3s %-62s %10s %5s %8s %4s %22s %7s %22s %7s\n", "preset", "opt", "bin", "range", "mod",
                    "parameters", "chi2/ndf", "", "prob", "tv", "S +- par +- tmpl", "Speak", "E +- par +- tmpl", "Ebump");
            for (int ip = 0; ip < kNP; ++ip)
                for (int io = 0; io < kNO; ++io)
                    for (int ib = 0; ib < kNBW; ++ib)
                        for (int ir = 0; ir < kNR; ++ir)
                            for (int im = 0; im < nM; ++im) {
                                const FitOut& r = R[IDX(ip, io, ib, ir, im)];
                                fprintf(o, "%-12s %-6s %-3s %.2f-%.2f %-3s %-62s %6.1f/%-3d %5.2f %8.1e %4.2f %7.0f +-%5.0f +-%5.0f %7.0f %7.0f +-%5.0f +-%5.0f %7.0f%s\n",
                                        kPName[ip], kOTag[io], kReb[ib] == 1 ? "2.5" : "5", kRLo[ir], kRHi[ir], models[im].tag.c_str(),
                                        ParStr(models[im], r, false).c_str(), r.chi2, r.ndf, r.ndf > 0 ? r.chi2 / r.ndf : 0., r.prob, r.tvShare, r.S, r.eS,
                                        r.tS, r.Speak, r.E, r.eE, r.tE, r.Ebump, r.conv ? "" : " NOT CONVERGED");
                            }
        }
        if (f) fclose(f);
    }
    // ---------------- correlations of the nominal fits ----------------
    {
        FILE* f = fopen((tab + "fit_correlations.txt").c_str(), "w");
        if (f) {
            fprintf(f, "# correlation coefficients of the fitted parameters, nominal configuration (mode A, w, 2.5 MeV, 0.43-0.75)\n");
            for (int ip = 0; ip < kNP; ++ip)
                for (int im = 0; im < nM; ++im) {
                    const FitOut& r = R[IDX(ip, 0, 0, 0, im)];
                    fprintf(f, "\n%s %s:\n", kPName[ip], models[im].tag.c_str());
                    for (int j = 0; j < r.np; ++j) {
                        fprintf(f, "  %-7s", models[im].pname[j].c_str());
                        for (int l = 0; l < r.np; ++l) fprintf(f, " %6.3f", r.cov[j * r.np + l] / (r.e[j] * r.e[l]));
                        fprintf(f, "\n");
                    }
                }
            fclose(f);
        }
    }
    // ---------------- contents table ----------------
    {
        FILE* f = fopen((tab + "fit_contents.txt").c_str(), "w");
        for (int pass = 0; pass < 2; ++pass) {
            FILE* o = pass ? f : stdout;
            if (!o) continue;
            fprintf(o, "# eta peak and bump contents implied by the fit (data rows); PRELIMINARY: no eta -> pi+pi-e+e- template; L_data = %.1f unconfirmed\n", kLdata);
            fprintf(o, "# data-study values (no vz cut): peak S = 4657 +- 390 +- 1233 (no cut), E = 1241 +- 76 +- 319 (step sum, not a yield)\n");
            fprintf(o, "# per preset: data S (P_main) +- stat +- model spread; data E +- stat (slides: E not valid, eta-rest cut)\n");
            fprintf(o, "# S_fit, E_fit: the same measures on the fitted model, +- parameters +- template statistics; ratios with the quadratic sum of all errors\n");
            fprintf(o, "# bump = a_bump x bump template (M4 0.30-1.10); peak = a_peak x peak template (M4 0.30-1.10); A_G = Gaussian content\n");
            for (int ip = 0; ip < kNP; ++ip) {
                fprintf(o, "\n## %s: data S = %.0f +- %.0f +- %.0f; data E = %.0f +- %.0f%s\n", kPName[ip], D[ip].pf.S, D[ip].pf.eS, D[ip].spr.half, D[ip].E,
                        D[ip].eE, kPresetId[ip] == kSlides ? " (not valid)" : "");
                fprintf(o, "%-6s %-3s %-9s %-3s %22s %14s %22s %14s %15s %15s %12s\n", "opt", "bin", "range", "mod", "S_fit", "S_fit/S_data", "E_fit",
                        "E_fit/E_data", "bump rows", "peak rows", "A_G");
                for (int io = 0; io < kNO; ++io)
                    for (int ib = 0; ib < kNBW; ++ib)
                        for (int ir = 0; ir < kNR; ++ir)
                            for (int im = 0; im < nM; ++im) {
                                const FitOut& r = R[IDX(ip, io, ib, ir, im)];
                                const double rs = r.S / D[ip].pf.S, rse = rs * std::sqrt((r.eS * r.eS + r.tS * r.tS) / (r.S * r.S) + std::pow(D[ip].pf.eS / D[ip].pf.S, 2));
                                const double re = r.E / D[ip].E, ree = re * std::sqrt((r.eE * r.eE + r.tE * r.tE) / (r.E * r.E) + std::pow(D[ip].eE / D[ip].E, 2));
                                fprintf(o, "%-6s %-3s %.2f-%.2f %-3s %6.0f +-%5.0f +-%5.0f %6.2f +- %-5.2f %6.0f +-%5.0f +-%5.0f %6.2f +- %-5.2f %7.0f+-%-6.0f %7.0f+-%-6.0f %6.0f+-%-5.0f\n",
                                        kOTag[io], kReb[ib] == 1 ? "2.5" : "5", kRLo[ir], kRHi[ir], models[im].tag.c_str(), r.S, r.eS, r.tS, rs, rse, r.E,
                                        r.eE, r.tE, re, ree, r.bumpAll, r.eBumpAll, r.peakAll, r.ePeakAll, r.q.AG, models[im].gauss ? r.e[r.np - 1] : 0.);
                            }
            }
        }
        if (f) fclose(f);
    }

    // =====================================================================
    // figures
    // =====================================================================
    const std::string dataLine = Form("data: N - CB, CB = 2#sqrt{N_{++}N_{--}} (10 MeV bins, no k factor); rows / L_{data} / bin width, %s", kLdataNote);
    const std::string simLine = Form("templates: SMASH std + lep owned rows, mode A, w (PDG); #alpha = 1 means SMASH (#sigma = #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1}); no efficiency correction", kLsim);
    const char* direct = "PRELIMINARY: #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) is NOT in the simulation; L_{data} unconfirmed";
    Binning bDisp; bDisp.nbin = kNB; bDisp.lo = kLo; bDisp.hi = kHi;
    double yHiLog = 0;
    for (int ip = 0; ip < kNP; ++ip)
        for (int i = 0; i + 1 < kNB; i += 2) {
            const double x = kLo + (i + 1) * kBW;
            if (x > 0.38 && x < 0.90) yHiLog = std::max(yHiLog, (D[ip].sp.n[0][i] + D[ip].sp.n[0][i + 1]) / (kLdata * 2 * kBW));
        }
    const double yLoLog = yHiLog * 2e-5, yHiLogP = yHiLog * 400;

    // ---------------- per preset: nominal fit M1 ----------------
    for (int ip = 0; ip < kNP; ++ip) {
        const int io = 0, ib = 0, ir = 0, im = 0;
        const FitOut& r = R[IDX(ip, io, ib, ir, im)];
        const TplSet& t = T[ip * kNO + io];
        const DSet& Ds = D[ip];
        const std::vector<WSpec> disp = FillMany(d, SelPreset(kPresetId[ip]), {CompFake(), CompDiff()}, SimOpt(), bDisp);
        auto comp = [&](std::vector<int> ks) {
            std::vector<double> y(kNB, 0.);
            for (int i = 0; i < kNB; ++i) { const double xc = kLo + (i + 0.5) * kBW; for (int k : ks) y[i] += r.q.at(k, xc) * t.s[k][i]; }
            return y;
        };
        const std::vector<double> yB = comp({cVM, cRest}), yBump = comp({cBump}), yPeak = comp({cPeak});
        std::vector<double> yTot(kNB);
        for (int i = 0; i < kNB; ++i) yTot[i] = yB[i] + yBump[i] + yPeak[i];
        std::vector<double> covN(kNB * kNB, 0.), covC(kNB * kNB, 0.);
        for (int i = 0; i < kNB; ++i) { covN[i * kNB + i] = Ds.sp.n[0][i]; covC[i * kNB + i] = Ds.sp.cbVar[i]; }
        std::vector<double> n0(Ds.sp.n[0].begin(), Ds.sp.n[0].end()), cb(Ds.sp.cb.begin(), Ds.sp.cb.end()), sg(Ds.sp.sig.begin(), Ds.sp.sig.end());
        TH1D* hN = HFine(Form("hN%d", ip), n0, &covN, 2);
        TH1D* hC = HFine(Form("hC%d", ip), cb, &covC, 2);
        TH1D* hS = HFine(Form("hS%d", ip), sg, &Ds.cov, 2);
        TH1D* hBk = HFine(Form("hBk%d", ip), yB, nullptr, 2);
        TH1D* hBu = HFine(Form("hBu%d", ip), yBump, nullptr, 2);
        TH1D* hPk = HFine(Form("hPk%d", ip), yPeak, nullptr, 2);
        TH1D* hTo = HFine(Form("hTo%d", ip), yTot, nullptr, 2);
        std::vector<double> yF(kNB), yD(kNB);
        for (int i = 0; i < kNB; ++i) { yF[i] = disp[0].s[i] * kU; yD[i] = disp[1].s[i] * kU; }
        TH1D* hF = HFine(Form("hF%d", ip), yF, nullptr, 2);
        TH1D* hDf = HFine(Form("hDf%d", ip), yD, nullptr, 2);
        Mk(hN, 1, 20, 0.8); Mk(hC, 2, 21, 0.8); Mk(hS, 4, 20, 0.8);
        Line(hBk, kColPhysBkg, 1, 1, kColPhysFill); Line(hTo, kColSim, 4, 1); Line(hBu, kColBumpDalitz, 3, 1); Line(hPk, kColPeakConv, 3, 1);
        Line(hF, kColFake, 2, 2); Line(hDf, kColDiff, 2, 2);

        TCanvas* c = new TCanvas(Form("cFit%d", ip), "", 2000, 1800);
        c->cd(); MkPad(Form("pa%d", ip), 0.0, 0.50, 1.0, 1.0, true, 0.075, 0.27);
        TH1D* fr = (TH1D*)hN->Clone(Form("fr%d", ip)); fr->Reset(); Axes(fr, kXt, kYt, 0.055, 0.62);
        fr->GetXaxis()->SetRangeUser(0.38, 0.90);
        fr->SetMinimum(yLoLog); fr->SetMaximum(yHiLogP); fr->Draw("axis");
        hBk->Draw("hist same"); hTo->Draw("hist same"); hBu->Draw("hist same"); hPk->Draw("hist same");
        hF->Draw("hist same"); hDf->Draw("hist same");
        hN->Draw("E1 P same"); hC->Draw("E1 P same"); hS->Draw("E1 P same");
        for (double xl : {kRLo[ir], kRHi[ir]}) { TLine* l = new TLine(xl, yLoLog, xl, yHiLog * 3); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
        fr->Draw("axis same");
        TLegend* lg = new TLegend(0.08, 0.52, 0.97, 0.725);
        lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
        lg->AddEntry(hN, "data: all e^{+}e^{-} combinations N", "p");
        lg->AddEntry(hTo, Form("fit total (#chi^{2}/ndf = %.0f/%d)", r.chi2, r.ndf), "l");
        lg->AddEntry(hBu, Form("bump: #alpha_{bump} = %.3f #pm %.3f", r.p[0], r.e[0]), "l");
        lg->AddEntry(hC, "data: CB = 2#sqrt{N_{++}N_{--}}", "p");
        lg->AddEntry(hBk, Form("phys. bkg: #alpha_{bkg} = %.3f #pm %.3f", r.p[2], r.e[2]), "f");
        lg->AddEntry(hPk, Form("peak conv.: #alpha_{peak} = %.2f #pm %.2f", r.p[1], r.e[1]), "l");
        lg->AddEntry(hS, "data: N - CB (fitted)", "p");
        lg->AddEntry(hF, "fake lepton, #alpha = 1, not in fit", "l");
        lg->AddEntry(hDf, "different sources, #alpha = 1, not in fit", "l");
        lg->Draw();
        Txt(0.075, 0.945, Form("(a) PRELIMINARY fit of N - CB, \"%s\", e^{+}e^{-}; model M1 (3 templates), fit %.2f-%.2f GeV/c^{2} (dashed), 2.5 MeV bins (shown: 5 MeV)",
                               kPName[ip], kRLo[ir], kRHi[ir]), 0.036);
        Txt(0.075, 0.903, PresetLabel(kPresetId[ip]), 0.032);
        Txt(0.075, 0.866, dataLine.c_str(), 0.030, 4);
        Txt(0.075, 0.829, simLine.c_str(), 0.030);
        Txt(0.075, 0.792, direct, 0.032, 2);
        // (b) linear zoom 0.43-0.62, 2.5 MeV, stacked fitted components
        c->cd(); MkPad(Form("pb%d", ip), 0.0, 0.0, 0.5, 0.50, false, 0.15, 0.17);
        TH1D* zS = HFine(Form("zS%d", ip), sg, &Ds.cov, 1);
        Mk(zS, 4, 20, 0.9);
        TH1D* zB = HFine(Form("zB%d", ip), yB, nullptr, 1);
        TH1D* zB1 = HFine(Form("zB1%d", ip), yBump, nullptr, 1); zB1->Add(zB);
        TH1D* zB2 = HFine(Form("zB2%d", ip), yPeak, nullptr, 1); zB2->Add(zB1);
        Line(zB, kColPhysBkg, 1, 1, kColPhysFill); Line(zB1, kColBumpDalitz, 1, 1, kColBumpDalitz); Line(zB2, kColPeakConv, 1, 1, kColPeakConv);
        TH1D* zT = (TH1D*)zB2->Clone(Form("zT%d", ip)); Line(zT, kColSim, 3, 1);
        double ym = 0;
        for (int i = 1; i <= zS->GetNbinsX(); ++i) { const double x = zS->GetBinCenter(i); if (x > 0.43 && x < 0.62) ym = std::max(ym, std::max(zS->GetBinContent(i) + zS->GetBinError(i), zT->GetBinContent(i))); }
        TH1D* fz = (TH1D*)zS->Clone(Form("fz%d", ip)); fz->Reset(); Axes(fz, kXt, kYt, 0.050, 1.45);
        fz->GetXaxis()->SetRangeUser(0.43, 0.62); fz->SetMinimum(0); fz->SetMaximum(ym * 1.30); fz->Draw("axis");
        zB2->Draw("hist same"); zB1->Draw("hist same"); zB->Draw("hist same"); zT->Draw("hist same"); zS->Draw("E1 P same");
        fz->Draw("axis same");
        TLegend* lz = new TLegend(0.17, 0.62, 0.62, 0.82); lz->SetBorderSize(0); lz->SetFillStyle(0); lz->SetTextSize(0.036);
        lz->AddEntry(zS, "data N - CB", "p"); lz->AddEntry(zB2, "+ #alpha_{peak} #times peak conversion", "f");
        lz->AddEntry(zB1, "+ #alpha_{bump} #times bump (#pi^{0} Dalitz + conv.)", "f"); lz->AddEntry(zB, "#alpha_{bkg} #times physical background", "f");
        lz->Draw();
        Txt(0.15, 0.93, "(b) linear zoom, 2.5 MeV bins, fitted parts stacked", 0.042);
        Txt(0.15, 0.885, "L_{data} = 6080.5 nb^{-1}, unconfirmed; no direct #eta decay in sim", 0.034, 2);
        // (c) pulls
        c->cd(); MkPad(Form("pc%d", ip), 0.5, 0.25, 1.0, 0.50, false, 0.13, 0.20, 0.03);
        const FitIn& in = nominalIn[ip];
        TH1D* hp = new TH1D(Form("hp%d", ip), "", kNB, kLo, kHi); hp->SetDirectory(nullptr);
        for (int q = 0; q < in.npts; ++q) hp->SetBinContent(in.pt[q] + 1, r.pull[q]);
        hp->SetFillColor(kAzure - 4); hp->SetLineColor(kAzure + 2);
        TH1D* fp = (TH1D*)hp->Clone(Form("fp%d", ip)); fp->Reset(); Axes(fp, "", "(data - fit) / #sigma", 0.085, 0.65);
        fp->GetXaxis()->SetLabelSize(0); fp->GetXaxis()->SetRangeUser(0.40, 0.85);
        fp->SetMinimum(-6); fp->SetMaximum(6); fp->Draw("axis"); hp->Draw("hist same");
        for (double yl : {-2., 0., 2.}) { TLine* l = new TLine(0.40, yl, 0.85, yl); l->SetLineStyle(yl == 0 ? 1 : 2); l->SetLineColor(kGray + 2); l->Draw(); }
        fp->Draw("axis same");
        Txt(0.13, 0.86, "(c) pulls; #sigma^{2} = data + #alpha^{2} #times template variance", 0.075);
        // (d) ratio data / fit
        c->cd(); MkPad(Form("pd%d", ip), 0.5, 0.0, 1.0, 0.25, false, 0.13, 0.15, 0.30);
        TH1D* hr = (TH1D*)hS->Clone(Form("hr%d", ip)); hr->Divide(hTo);
        for (int I = 1; I <= hr->GetNbinsX(); ++I) hr->SetBinError(I, hTo->GetBinContent(I) > 0 ? hS->GetBinError(I) / hTo->GetBinContent(I) : 0);
        Mk(hr, 4, 20, 0.8);
        TH1D* fq = (TH1D*)hr->Clone(Form("fq%d", ip)); fq->Reset(); Axes(fq, kXt, "data / fit", 0.10, 0.55);
        fq->GetXaxis()->SetRangeUser(0.40, 0.85); fq->SetMinimum(0.6); fq->SetMaximum(1.4); fq->GetYaxis()->SetNdivisions(504); fq->Draw("axis");
        TLine* one = new TLine(0.40, 1, 0.85, 1); one->SetLineStyle(2); one->Draw();
        for (double xl : {kRLo[ir], kRHi[ir]}) { TLine* l = new TLine(xl, 0.6, xl, 1.4); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
        hr->Draw("E1 P same");
        Txt(0.13, 0.89, "(d) data (N - CB) / fit, 5 MeV bins, data errors only", 0.085);
        c->SaveAs((fig + "fit_spectra_" + kPTag[ip] + ".png").c_str());
        c->SaveAs((fig + "fit_spectra_" + kPTag[ip] + ".pdf").c_str());
    }

    // ---------------- per preset: the models compared ----------------
    {
        const int show[6] = {0, 2, 3, 4, 5, 6};
        const int mcol[7] = {kColSim, kMagenta + 1, kRed + 1, kOrange + 7, kBlue + 1, kViolet + 1, kGray + 3};
        const int msty[7] = {1, 2, 1, 7, 9, 1, 2};
        for (int ip = 0; ip < kNP; ++ip) {
            const int io = 0, ib = 0, ir = 0;
            const TplSet& t = T[ip * kNO + io];
            const DSet& Ds = D[ip];
            const FitIn& in = nominalIn[ip];
            TCanvas* c = new TCanvas(Form("cModels%d", ip), "", 2000, 1800);
            c->cd(); MkPad(Form("mL%d", ip), 0.0, 0.0, 0.5, 1.0, false, 0.15, 0.24, 0.08, 0.05);
            std::vector<double> sg(Ds.sp.sig.begin(), Ds.sp.sig.end());
            TH1D* zS = HFine(Form("mzS%d", ip), sg, &Ds.cov, 2);
            Mk(zS, 4, 20, 1.0);
            double ym = 0;
            for (int i = 1; i <= zS->GetNbinsX(); ++i) { const double x = zS->GetBinCenter(i); if (x > 0.43 && x < 0.75) ym = std::max(ym, zS->GetBinContent(i) + zS->GetBinError(i)); }
            TH1D* fz = (TH1D*)zS->Clone(Form("mfz%d", ip)); fz->Reset(); Axes(fz, kXt, kYt, 0.040, 1.6);
            fz->GetXaxis()->SetRangeUser(0.43, 0.75); fz->SetMinimum(0); fz->SetMaximum(ym * 1.45); fz->Draw("axis");
            TLegend* lg = new TLegend(0.17, 0.60, 0.97, 0.75); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.024);
            lg->AddEntry(zS, "data N - CB", "p");
            for (int im = 0; im < nM; ++im) {
                const FitOut& r = R[IDX(ip, io, ib, ir, im)];
                TH1D* h = HFine(Form("mod%d_%d", ip, im), ModelFine(t, r.q, gm, gs), nullptr, 2);
                Line(h, mcol[im], 2, msty[im]); h->Draw("hist same");
                lg->AddEntry(h, Form("%s %s: #chi^{2}/ndf %.0f/%d", models[im].tag.c_str(), models[im].label.c_str(), r.chi2, r.ndf), "l");
            }
            zS->Draw("E1 P same"); fz->Draw("axis same"); lg->Draw();
            Txt(0.15, 0.965, Form("PRELIMINARY: fit models, \"%s\", e^{+}e^{-}, fit 0.43-0.75", kPName[ip]), 0.030);
            Txt(0.15, 0.94, "fit in 2.5 MeV bins; spectra shown in 5 MeV bins", 0.022);
            Txt(0.15, 0.915, "templates std + lep, mode A, w (PDG); #alpha = 1 means SMASH", 0.022);
            Txt(0.15, 0.89, "data: L_{data} = 6080.5 nb^{-1}, unconfirmed", 0.022, 4);
            Txt(0.15, 0.865, Form("G: Gaussian at %.4f GeV/c^{2}, #sigma = %.1f MeV, free content", gm, gs * 1000), 0.022);
            Txt(0.15, 0.84, "M3, M5, M7: #alpha_{peak} = #alpha_{bump}; M4, M5: #omega/#rho^{0} pions vs rest", 0.022);
            Txt(0.15, 0.815, "M6, M7: background #times (#alpha_{bkg} + #beta (M - 0.6)/0.1)", 0.022);
            Txt(0.15, 0.79, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation (G stands for it)", 0.022, 2);
            for (int s = 0; s < 6; ++s) {
                const int im = show[s];
                const FitOut& r = R[IDX(ip, io, ib, ir, im)];
                const double h = 0.94 / 6, y2 = 0.98 - s * h, y1 = y2 - h;
                const bool last = s == 5;
                c->cd(); MkPad(Form("mR%d_%d", ip, s), 0.5, last ? 0.0 : y1, 1.0, y2, false, 0.10, 0.02, last ? 0.27 : 0.02);
                TH1D* hp = new TH1D(Form("mp%d_%d", ip, s), "", kNB, kLo, kHi); hp->SetDirectory(nullptr);
                for (int q = 0; q < in.npts; ++q) hp->SetBinContent(in.pt[q] + 1, r.pull[q]);
                hp->SetFillColor(mcol[im]); hp->SetLineColor(mcol[im]);
                TH1D* fp = (TH1D*)hp->Clone(Form("mfp%d_%d", ip, s)); fp->Reset();
                Axes(fp, last ? kXt : "", "pull", last ? 0.085 : 0.11, 0.40);
                fp->GetXaxis()->SetRangeUser(0.43, 0.75); fp->SetMinimum(-6); fp->SetMaximum(6); fp->GetYaxis()->SetNdivisions(505);
                if (!last) fp->GetXaxis()->SetLabelSize(0);
                fp->Draw("axis"); hp->Draw("hist same");
                for (double yl : {-2., 2.}) { TLine* l = new TLine(0.43, yl, 0.75, yl); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
                Txt(0.11, last ? 0.88 : 0.84, Form("%s  #chi^{2}/ndf = %.0f/%d:  %s", models[im].tag.c_str(), r.chi2, r.ndf, ParStr(models[im], r, true).c_str()),
                    last ? 0.052 : 0.066);
            }
            c->SaveAs((fig + "fit_models_" + kPTag[ip] + ".png").c_str());
            if (ip == 0) c->SaveAs((fig + "fit_models_" + kPTag[ip] + ".pdf").c_str());
        }
    }

    // ---------------- templates and their precision ("pre") ----------------
    {
        const TplSet& t = T[0];
        const DSet& Ds = D[0];
        TCanvas* c = new TCanvas("cTpl", "", 2000, 1600);
        c->cd(); MkPad("tA", 0.0, 0.0, 0.5, 1.0, true, 0.15, 0.22, 0.10);
        std::vector<double> sg(Ds.sp.sig.begin(), Ds.sp.sig.end());
        TH1D* hS = HFine("tS", sg, &Ds.cov, 2); Mk(hS, 4, 20, 0.8);
        const int col[kNC] = {kColBumpDalitz, kColPeakConv, kViolet + 1, kColPhysBkg};
        const char* nm[kNC] = {"bump (#pi^{0} Dalitz + conversion)", "peak (conversion)", "phys. bkg, #omega/#rho^{0} pions", "phys. bkg, other pions"};
        TH1D* fr = (TH1D*)hS->Clone("tfr"); fr->Reset(); Axes(fr, kXt, kYt, 0.040, 1.6);
        fr->GetXaxis()->SetRangeUser(0.38, 0.90); fr->SetMinimum(yLoLog); fr->SetMaximum(yHiLogP); fr->Draw("axis");
        TLegend* lg = new TLegend(0.17, 0.62, 0.97, 0.77); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.026);
        lg->AddEntry(hS, "data N - CB", "p");
        std::vector<double> sum(kNB, 0.);
        for (int k = 0; k < kNC; ++k) {
            for (int i = 0; i < kNB; ++i) sum[i] += t.s[k][i];
            TH1D* h = HFine(Form("tT%d", k), t.s[k], &t.C[k][k], 2); Line(h, col[k], 3, 1); h->Draw("hist same");
            lg->AddEntry(h, Form("%s: %.0f rows, N_{eff} %.0f", nm[k], t.Sall[k], t.Vall[k] > 0 ? t.Sall[k] * t.Sall[k] / t.Vall[k] : 0.), "l");
        }
        TH1D* hSum = HFine("tSum", sum, nullptr, 2); Line(hSum, kColSim, 3, 2); hSum->Draw("hist same");
        lg->AddEntry(hSum, "sum of the four (= signal + physical background)", "l");
        hS->Draw("E1 P same"); fr->Draw("axis same"); lg->Draw();
        Txt(0.15, 0.95, "(a) templates at #alpha = 1, \"pre\", e^{+}e^{-}, 5 MeV bins", 0.032);
        Txt(0.15, 0.915, "std + lep owned rows, mode A, w (PDG); rows = expected data rows", 0.024);
        Txt(0.15, 0.885, "data and templates: L_{data} = 6080.5 nb^{-1}, unconfirmed", 0.024, 4);
        Txt(0.15, 0.855, "edge of the #omega/#rho^{0}-pion template at 0.745 (#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}, #pi^{0} Dalitz)", 0.024);
        Txt(0.15, 0.825, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation", 0.024, 2);
        // (b) relative precision per 2.5 MeV bin
        c->cd(); MkPad("tB", 0.5, 0.0, 1.0, 1.0, true, 0.15, 0.22, 0.10);
        TH1D* rD = new TH1D("rD", "", kNB, kLo, kHi); rD->SetDirectory(nullptr);
        TH1D* rT = new TH1D("rT", "", kNB, kLo, kHi); rT->SetDirectory(nullptr);
        TH1D* rR = new TH1D("rR", "", kNB, kLo, kHi); rR->SetDirectory(nullptr);
        TH1D* rB = new TH1D("rB", "", kNB, kLo, kHi); rB->SetDirectory(nullptr);
        for (int i = 0; i < kNB; ++i) {
            double vs = 0;
            for (int k = 0; k < kNC; ++k) for (int l = 0; l < kNC; ++l) vs += t.C[k][l][i * kNB + i];
            if (Ds.sp.sig[i] > 0) rD->SetBinContent(i + 1, std::sqrt(Ds.cov[i * kNB + i]) / Ds.sp.sig[i]);
            if (sum[i] > 0) rT->SetBinContent(i + 1, std::sqrt(vs) / sum[i]);
            if (t.s[cRest][i] > 0) rR->SetBinContent(i + 1, std::sqrt(t.C[cRest][cRest][i * kNB + i]) / t.s[cRest][i]);
            if (t.s[cBump][i] > 0 && kLo + (i + 0.5) * kBW < 0.53) rB->SetBinContent(i + 1, std::sqrt(t.C[cBump][cBump][i * kNB + i]) / t.s[cBump][i]);
        }
        Mk(rD, 4, 20, 0.9); rD->SetLineWidth(0);
        Line(rT, kColSim, 3, 1); Line(rR, kColPhysBkg, 2, 2); Line(rB, kColBumpDalitz, 2, 1);
        TH1D* fb = (TH1D*)rD->Clone("tfb"); fb->Reset(); Axes(fb, kXt, "relative error per 2.5 MeV bin", 0.040, 1.6);
        fb->GetXaxis()->SetRangeUser(0.40, 0.85); fb->SetMinimum(0.003); fb->SetMaximum(1.0); fb->Draw("axis");
        rT->Draw("hist same"); rR->Draw("hist same"); rB->Draw("hist same"); rD->Draw("P same");
        TLegend* lb = new TLegend(0.17, 0.66, 0.97, 0.77); lb->SetBorderSize(0); lb->SetFillStyle(0); lb->SetTextSize(0.026);
        lb->AddEntry(rD, "data N - CB (event groups + CB propagation)", "p");
        lb->AddEntry(rT, "sum of the templates (event-group covariance)", "l");
        lb->AddEntry(rR, "phys. bkg, other pions", "l");
        lb->AddEntry(rB, "bump template (below 0.53)", "l");
        lb->Draw();
        Txt(0.15, 0.95, "(b) statistical precision per bin, \"pre\"", 0.032);
        Txt(0.15, 0.915, "the template errors are 3-6 times larger than the data errors:", 0.024);
        Txt(0.15, 0.885, "the template variance is about 95 % of the variance in the fit", 0.024);
        c->SaveAs((fig + "fit_templates_pre.png").c_str());
        c->SaveAs((fig + "fit_templates_pre.pdf").c_str());
    }

    // ---------------- stability summary (forest plots) ----------------
    struct VarDef { const char* label; int io, ib, ir, im; };
    const std::vector<VarDef> vars = {
        {"M1 nominal (0.43-0.75, 2.5 MeV, A, w)", 0, 0, 0, 0}, {"M1, range 0.45-0.70", 0, 0, 1, 0}, {"M1, range 0.40-0.85", 0, 0, 2, 0},
        {"M1, 5 MeV bins", 0, 1, 0, 0}, {"M1, mode B", 1, 0, 0, 0}, {"M1, w_{db}", 2, 0, 0, 0}, {"M1, sim M4 #times 0.995 (K0S)", 3, 0, 0, 0},
        {"M2: + G", 0, 0, 0, 1}, {"M3: + G, #alpha_{peak} = #alpha_{bump}", 0, 0, 0, 2},
        {"M4: split background", 0, 0, 0, 3}, {"M4, range 0.45-0.70", 0, 0, 1, 3}, {"M4, range 0.40-0.85", 0, 0, 2, 3}, {"M4, mode B", 1, 0, 0, 3},
        {"M5: split + G, #alpha_{peak} = #alpha_{bump}", 0, 0, 0, 4}, {"M5, range 0.40-0.85", 0, 0, 2, 4},
        {"M6: background tilt", 0, 0, 0, 5}, {"M6, range 0.45-0.70", 0, 0, 1, 5}, {"M6, range 0.40-0.85", 0, 0, 2, 5}, {"M6, mode B", 1, 0, 0, 5},
        {"M7: tilt + G, #alpha_{peak} = #alpha_{bump}", 0, 0, 0, 6}, {"M7, range 0.40-0.85", 0, 0, 2, 6}, {"M7, sim M4 #times 0.995 (K0S)", 3, 0, 0, 6}};
    const int nv = (int)vars.size();
    const int pcol[kNP] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
    const int pmk[kNP] = {20, 21, 22, 33};
    auto getPar = [&](const FitOut& r, const Model& m, const std::string& name, double& v, double& e) {
        for (int j = 0; j < r.np; ++j) if (m.pname[j] == name) { v = r.p[j]; e = r.e[j]; return true; }
        return false;
    };
    auto forest = [&](const char* cname, const char* file, int npan, const char* const* ptitle,
                      std::function<bool(int, int, int, double&, double&)> get, const double* xlo, const double* xhi, const double* xref,
                      const char* head1, const char* head2) {
        TCanvas* c = new TCanvas(cname, "", 2600, 1900);
        const double lab = 0.17, w = (1.0 - lab) / npan, tm = 0.10, bm = 0.05;
        c->cd(); MkPad(Form("%s_lab", cname), 0.0, 0.0, lab, 1.0, false, 0.0, tm, bm, 0.0);
        for (int iv = 0; iv < nv; ++iv) {
            const double y = bm + (1.0 - tm - bm) * (nv - iv - 0.5) / nv;
            Txt(0.97, y, vars[iv].label, 0.060, 1, 32);
        }
        for (int k = 0; k < npan; ++k) {
            c->cd(); MkPad(Form("%s_%d", cname, k), lab + k * w, 0.0, lab + (k + 1) * w, 1.0, false, 0.05, tm, bm, 0.05);
            TH1D* fr = new TH1D(Form("%s_fr%d", cname, k), "", 10, xlo[k], xhi[k]); fr->SetDirectory(nullptr);
            fr->SetMinimum(0); fr->SetMaximum(nv); fr->GetYaxis()->SetLabelSize(0); fr->GetYaxis()->SetTickLength(0);
            fr->GetXaxis()->SetLabelSize(24.0 / (2600 * w)); fr->GetXaxis()->SetNdivisions(504); fr->GetXaxis()->SetLabelOffset(-0.01); fr->GetXaxis()->SetTickLength(0.01); fr->Draw("axis");
            for (int iv = 0; iv < nv; ++iv) { TLine* l = new TLine(xlo[k], iv, xhi[k], iv); l->SetLineColor(kGray); l->SetLineStyle(3); l->Draw(); }
            if (xref && xref[k] > xlo[k] && xref[k] < xhi[k]) { TLine* l = new TLine(xref[k], 0, xref[k], nv); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
            for (int ip = 0; ip < kNP; ++ip)
                for (int iv = 0; iv < nv; ++iv) {
                    double v, e;
                    if (!get(k, ip, iv, v, e)) continue;
                    const double y = nv - iv - 0.5 + (1.5 - ip) * 0.2;
                    const double vv = std::max(xlo[k], std::min(xhi[k], v));
                    TLine* l = new TLine(std::max(xlo[k], v - e), y, std::min(xhi[k], v + e), y); l->SetLineColor(pcol[ip]); l->SetLineWidth(2); l->Draw();
                    TMarker* mk = new TMarker(vv, y, (v < xlo[k] || v > xhi[k]) ? 24 + ip % 3 : pmk[ip]);
                    mk->SetMarkerColor(pcol[ip]); mk->SetMarkerSize(1.5); mk->Draw();
                }
            Txt(0.06, 0.92, ptitle[k], 28.0 / (2600 * w));
        }
        c->cd();
        Txt(0.005, 0.975, head1, 0.0145);
        Txt(0.005, 0.950, head2, 0.0145);
        TLegend* lg = new TLegend(0.005, 0.90, 0.17, 0.94); lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.014);
        for (int ip = 0; ip < kNP; ++ip) { TMarker* mk = new TMarker(0, 0, pmk[ip]); mk->SetMarkerColor(pcol[ip]); mk->SetMarkerSize(1.5); lg->AddEntry(mk, kPName[ip], "p"); }
        lg->Draw();
        c->SaveAs((fig + file + ".png").c_str());
        c->SaveAs((fig + file + ".pdf").c_str());
    };
    {
        const char* titles[7] = {"#alpha_{bump} or #alpha_{#eta}", "#alpha_{peak}", "#alpha_{bkg} or #alpha_{rest}", "#alpha_{VM}", "tilt #beta", "A_{G} [rows]", "#chi^{2}/ndf"};
        const double xlo[7] = {0.0, -1.0, 0.8, -0.5, -0.3, -6000, 0.5}, xhi[7] = {1.2, 4.0, 1.8, 2.0, 0.3, 16000, 6.0};
        const double xref[7] = {1, 1, 1, 1, 0, 0, 1};
        auto get = [&](int k, int ip, int iv, double& v, double& e) -> bool {
            const VarDef& V = vars[iv];
            const FitOut& r = R[IDX(ip, V.io, V.ib, V.ir, V.im)];
            const Model& m = models[V.im];
            if (!r.ok) return false;
            switch (k) {
                case 0: return getPar(r, m, "a_bump", v, e) || getPar(r, m, "a_eta", v, e);
                case 1: return getPar(r, m, "a_peak", v, e);
                case 2: return getPar(r, m, "a_bkg", v, e) || getPar(r, m, "a_rest", v, e);
                case 3: return getPar(r, m, "a_VM", v, e);
                case 4: return getPar(r, m, "b_tilt", v, e);
                case 5: return getPar(r, m, "A_G", v, e);
                case 6: v = r.chi2 / r.ndf; e = std::sqrt(2.0 / r.ndf); return true;
            }
            return false;
        };
        forest("cAlpha", "fit_alpha_summary", 7, titles, get, xlo, xhi, xref,
               "PRELIMINARY template fit of data N - CB (e^{+}e^{-}): fitted parameters per variant (rows) and selection (colours); #alpha = 1 means SMASH with L_{data} = 6080.5 nb^{-1} (unconfirmed: every #alpha scales with 1/L_{data})",
               "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation (G stands for it); open markers: value outside the axis range; errors: fit (data + template variance)");
    }
    {
        const char* titles[4] = {"S_{fit} / S_{data} (P_{main})", "E_{fit} / E_{data} (edge steps)", "#alpha_{bump} #times bump [rows]", "peak tmpl. + G [rows]"};
        const double xlo[4] = {-0.5, 0.0, 0, 0}, xhi[4] = {3.0, 2.0, 100000, 30000};
        const double xref[4] = {1, 1, 0, 0};
        auto get = [&](int k, int ip, int iv, double& v, double& e) -> bool {
            const VarDef& V = vars[iv];
            const FitOut& r = R[IDX(ip, V.io, V.ib, V.ir, V.im)];
            if (!r.ok) return false;
            switch (k) {
                case 0: v = r.S / D[ip].pf.S;
                        e = std::fabs(v) * std::sqrt((r.eS * r.eS + r.tS * r.tS) / (r.S * r.S) + std::pow(D[ip].pf.eS / D[ip].pf.S, 2)); return true;
                case 1: if (kPresetId[ip] == kSlides) return false;
                        v = r.E / D[ip].E; e = std::fabs(v) * std::sqrt((r.eE * r.eE + r.tE * r.tE) / (r.E * r.E) + std::pow(D[ip].eE / D[ip].E, 2)); return true;
                case 2: v = r.bumpAll; e = r.eBumpAll; return true;
                case 3: v = r.peakAll + r.q.AG; e = std::sqrt(r.ePeakAll * r.ePeakAll + (models[V.im].gauss ? std::pow(r.e[r.np - 1], 2) : 0.)); return true;
            }
            return false;
        };
        forest("cCont", "fit_contents_summary", 4, titles, get, xlo, xhi, xref,
               "PRELIMINARY: eta peak (P_{main} of the data study) and bump (edge-step sum E) measured on the fitted model, against the same measures on the data (ratio 1 = agreement)",
               "errors: fit parameters + template statistics + data statistics; slides: E not valid (#eta-rest cut); L_{data} = 6080.5 nb^{-1} unconfirmed; no direct #eta decay in sim");
    }
    if (fchk) fclose(fchk);
    printf("sim_templatefit: total %.1f s\n", swAll.RealTime());
}
