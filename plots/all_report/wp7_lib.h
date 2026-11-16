// ============================================================================
// wp7_lib.h  --  WP7 systematics of SPEC report_all (prefix wp7): generic layouts, data from rows, all refits
// ============================================================================
// Written by WP7, 2026-11-14. Everything here is a GENERALISATION of the validated WP4/WP5/WP6 code (no new method):
//   * data: N, CB (spread rule, 2 sqrt(N++ N--) in 10 MeV cells, interpolated, spread by f(cell | 50 MeV band)),
//     200 CB replicas (wp0 generator and seed) and the event-group covariance are rebuilt from the rows
//     (wp7_rows.root) for ANY layout, selection, subset and CB variant (wp6_ana.C LoadDataW / CBWin generalised);
//   * simulation: F10/F5g/win/tot of hists_v1.root (every sim directory) or the G3 layout of hists_v1.root
//     (S0__nom, pre__nom) summed into any layout (range, binning, slice edges, M(ee) split, K0S veto);
//   * smoothing: wp4::SmoothFit (F-main/F-gold, WP5 rule) and the WP3 smoothing of wp6_ana.C (stage A, w3::Smooth,
//     copied below and generalised to exp(pol n), n = 2..4);
//   * fit: copy of wp5::Eval / wp5::Fit (Fisher scoring, BB-lite, sandwich with the event-group and CB-replica
//     covariances) with optional model terms: exp(b x) tilt, omega/rho0 split, common c (1 + c1 y), kappa_B = kappa,
//     BB-lite on raw background templates;
//   * estimators: B4 and the peak side bands exactly as wp6_ana.C (MakeB4, MakeB4All, MakePSB mode 0), errors N + CB +
//     stage A + F-main as wp6::Evaluate.
// The nominal configuration reproduces WP5 (F-main, F-gold) and WP6 (B4, side bands) to rounding (checked in wp7_syst.C).
// ============================================================================
#ifndef WP7_LIB_H
#define WP7_LIB_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TMath.h"

#include "wp5_lib.h"

namespace wp7 {
using namespace alldefs;
typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;
using wp5::Spec;
using wp5::Opt;
using wp5::FitRes;
using wp5::Term;
using wp5::Zero;
using wp5::Quad;

const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/wp7";
const char* const kRowsFile = "/hdd3/przygoda/HADES/FAT_out/all_report/wp7/wp7_rows.root";
const uint64_t kSeedCB = 0xCB2026ULL;   // = wp0_fill.C kSeedCB, wp6_ana.C kSeedCB6: the CB replicas of hists_v1
const int kNRepCB = 200;
const double kYf = kPT3;                // yields = raw rows x 64/63
const int kNSg = 5;                     // C2 C3 C4 C5 G_eta (index = wp5::SigT)
const int kSigComp[kNSg] = {kC2, kC3, kC4, kC5, kGeta};
const double kSigEta = 1.0665;          // mb, SMASH sigma(pp -> eta X)
const double kNormEta = 1.0110 * 0.97;  // C_sec x (1 - f_eta), spec AMENDMENT 1 (0.98 +0.03/-0.04)

// ---------------------------------------------------------------------------------------------------------------
// replica generator: copy of wp0_fill.C (same seed -> same CB replicas as hists_v1 CBrep)
// ---------------------------------------------------------------------------------------------------------------
inline uint64_t Mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline int Pois1(uint64_t h) {
    const double u = (double)(h >> 11) * (1.0 / 9007199254740992.0);
    double p = 0.36787944117144233, c = p;
    int k = 0;
    while (u >= c && k < 30) { ++k; p /= k; c += p; }
    return k;
}
inline void GroupPois(uint64_t seed, int tree, int grp, int n, double* k) {
    const uint64_t base = Mix(seed ^ Mix(((uint64_t)(uint32_t)tree << 32) | (uint64_t)(uint32_t)grp));
    for (int r = 0; r < n; ++r) k[r] = Pois1(Mix(base + 0x632BE59BD9B4E019ULL * (uint64_t)(r + 1)));
}

// ---------------------------------------------------------------------------------------------------------------
// 1. rows
// ---------------------------------------------------------------------------------------------------------------
struct Rows {
    size_t n = 0;
    std::vector<int8_t> k, isBest, cor, vtx, per, fold;
    std::vector<int> grp;
    std::vector<float> m4, mpipi, mee, mm4, vz, dz;
};
inline bool LoadRows(Rows& R) {
    TFile* f = TFile::Open(kRowsFile);
    if (!f) return false;
    TTree* t = (TTree*)f->Get("rows");
    Char_t k, isBest, cor, vtx, per, fold; Int_t grp; Float_t m4, mpipi, mee, mm4, vz, dz;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"k", "isBest", "cor", "vtx", "per", "fold", "grp", "m4", "mpipi", "mee", "mm4", "vz", "dz"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("k", &k); t->SetBranchAddress("isBest", &isBest); t->SetBranchAddress("cor", &cor); t->SetBranchAddress("vtx", &vtx);
    t->SetBranchAddress("per", &per); t->SetBranchAddress("fold", &fold); t->SetBranchAddress("grp", &grp);
    t->SetBranchAddress("m4", &m4); t->SetBranchAddress("mpipi", &mpipi); t->SetBranchAddress("mee", &mee);
    t->SetBranchAddress("mm4", &mm4); t->SetBranchAddress("vz", &vz); t->SetBranchAddress("dz", &dz);
    const Long64_t n = t->GetEntries();
    R.n = n;
    R.k.resize(n); R.isBest.resize(n); R.cor.resize(n); R.vtx.resize(n); R.per.resize(n); R.fold.resize(n); R.grp.resize(n);
    R.m4.resize(n); R.mpipi.resize(n); R.mee.resize(n); R.mm4.resize(n); R.vz.resize(n); R.dz.resize(n);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        R.k[i] = k; R.isBest[i] = isBest; R.cor[i] = cor; R.vtx[i] = vtx; R.per[i] = per; R.fold[i] = fold; R.grp[i] = grp;
        R.m4[i] = m4; R.mpipi[i] = mpipi; R.mee[i] = mee; R.mm4[i] = mm4; R.vz[i] = vz; R.dz[i] = dz;
    }
    f->Close();
    return true;
}
// selection of a row (= alldefs::PassSel on the stored columns; corfb: cor_ok == 0 rows with the _fb kinematics kept)
struct SelCfg {
    int sel = kSelS0;
    bool best = false;     // isBest == 1
    bool corfb = false;    // cor_ok fallback
    int sub = kSubAll;     // kSubAll, kSubFold0, kSubFold1, kSubPerI + p
    std::string Tag() const {
        std::string t = SelTag(sel);
        if (best) t += "_best";
        if (corfb) t += "_corfb";
        if (sub != kSubAll) t += std::string("_") + SubTag(sub);
        return t;
    }
};
inline bool PassRow(const Rows& R, size_t i, const SelCfg& s) {
    if (!s.corfb && !R.cor[i]) return false;
    if (s.best && R.isBest[i] != 1) return false;
    const float vz = R.vz[i];
    const bool vtx = R.vtx[i] != 0;
    bool ok = false;
    switch (s.sel) {
        case kSelS0:    ok = vz > -160.f && vtx; break;
        case kSelPre:   ok = vz > -160.f; break;
        case kSelDz40:  ok = vz > -160.f && R.dz[i] < 40.f; break;
        case kSelPeak:  ok = vz > -160.f && R.dz[i] < 40.f && R.mm4[i] < 2.85f; break;
        case kSelVz150: ok = vz > -150.f && vtx; break;
        case kSelVz170: ok = vz > -170.f && vtx; break;
        default: ok = false;
    }
    if (!ok) return false;
    if (s.sub == kSubFold0 && R.fold[i] != 0) return false;
    if (s.sub == kSubFold1 && R.fold[i] != 1) return false;
    if (s.sub >= kSubPerI && R.per[i] != s.sub - kSubPerI) return false;
    return true;
}

// ---------------------------------------------------------------------------------------------------------------
// 2. layouts and windows
// ---------------------------------------------------------------------------------------------------------------
struct Layout {
    std::string tag = "std";
    Vec e;                                   // M4 bin edges of the 10 main (slice x class) and 5 F5g histograms
    double se[4] = {0.34, 0.37, 0.39, 0.43}; // slice edges (data scale)
    double eSplit = 0.010;                   // E1 < eSplit <= E2
    bool veto = false;                       // K0S veto: M(pipi) 0.485-0.507 removed from every histogram
    double simEdgeRel = 0;                   // simulation only: slice edges x (1 + simEdgeRel) (G3 source, linear sharing)
    int nb() const { return (int)e.size() - 1; }
};
inline Vec EdgesF10like(double hi, double w1 = 0.005, double w2 = 0.0025, double w3 = 0.005) {
    Vec e;
    for (double x = 0.400; x < 0.520 - 1e-9; x += w1) e.push_back(x);
    for (double x = 0.520; x < 0.570 - 1e-9; x += w2) e.push_back(x);
    for (double x = 0.570; x < hi - 1e-9; x += w3) e.push_back(x);
    e.push_back(hi);
    for (auto& x : e) x = std::round(x * 1e6) / 1e6;
    return e;
}
inline Layout StdLayout() { Layout L; L.e = F10Edges(); return L; }
inline int LSlice(const Layout& L, double m) {
    if (!(m == m)) return -1;
    if (L.veto && m >= 0.485 && m < 0.507) return -1;
    int j = 0;
    while (j < 4 && m >= L.se[j]) ++j;
    return j;
}
inline int LClass(const Layout& L, double q) { if (!(q == q)) return -1; return q < L.eSplit ? 0 : 1; }

struct WinDef { const char* tag; double lo, hi; bool aligned; };
enum { kWB = 0, kWP, kWSl, kWSh, kWU575, kWU565, kWU560, kWU650, kWU750, kWCRL, kWK400, kWR0, kWU850, kWtot, kNW };
const WinDef kWins[kNW] = {
    {"B", 0.4550, 0.5225, true},      {"P", 0.53364, 0.55416, false},  {"SBlo", 0.5250, 0.5336, false}, {"SBhi", 0.5542, 0.5628, false},
    {"U575", 0.575, 0.700, true},     {"U565", 0.565, 0.700, true},    {"U560", 0.560, 0.700, true},    {"U575_650", 0.575, 0.650, true},
    {"U575_750", 0.575, 0.750, true}, {"CRL", 0.5225, 0.5300, true},   {"K400", 0.400, 0.5225, true},   {"R0", 0.400, 0.455, true},
    {"U575_850", 0.575, 0.850, true}, {"tot", -1e9, 1e9, false}};
inline double FracFine(int i, double a, double b) {   // copy of wp6_ana.C
    if (a < -1e8) return 1.0;
    const double x1 = kFineLo + i * kFineW, x2 = x1 + kFineW;
    const double o = std::min(x2, b) - std::max(x1, a);
    return o > 0 ? o / kFineW : 0.;
}

// ---------------------------------------------------------------------------------------------------------------
// 3. data in a layout (raw rows): items = 10 main histograms x nb, 5 F5g x nb, 10 cells x kNW windows
// ---------------------------------------------------------------------------------------------------------------
struct CBOpt {
    int rule = 0;      // 0: spread rule (nominal); 1: like-sign arithmetic sum N++ + N-- (same spread); 2: per-cell geometric mean
    int reb = kCBreb;  // fine bins per CB cell (4 = 10 MeV, 8 = 20 MeV)
    double k = 1.0;    // overall factor on the CB
    std::string Tag() const { return Form("cb%d_r%d_k%.3f", rule, reb, k); }
};
struct DataL {
    Layout L; SelCfg S; CBOpt C;
    int nb = 0, nI = 0;
    Vec N, CB;
    std::vector<double> covN, covCB;
    double mpp[kNSlice] = {0, 0, 0, 0, 0};   // mean M(pipi) of the e+e- rows in 0.40-0.70 per slice
    long nUS = 0, nLS = 0;
    double Ntot = 0;                           // e+e- rows in the layout cells at any M4
    int IMain(int h, int x) const { return h * nb + x; }
    int IF5(int j, int x) const { return 10 * nb + j * nb + x; }
    int IWin(int h, int w) const { return 15 * nb + h * kNW + w; }
};
struct ItemFine { std::vector<std::pair<int, double>> f; int cell; };   // fine bins and overlap fraction; CB cell
inline void ItemDefs(const DataL& D, std::vector<ItemFine>& it) {
    it.assign(D.nI, ItemFine());
    const Vec& e = D.L.e;
    auto fill = [&](int idx, int cell, double lo, double hi) {
        it[idx].cell = cell;
        for (int i = 0; i < kNfine; ++i) { const double u = FracFine(i, lo, hi); if (u > 0) it[idx].f.push_back({i, u}); }
    };
    for (int h = 0; h < 10; ++h) for (int x = 0; x < D.nb; ++x) fill(D.IMain(h, x), h, e[x], e[x + 1]);
    for (int j = 0; j < 5; ++j) for (int x = 0; x < D.nb; ++x) fill(D.IF5(j, x), 10 + j, e[x], e[x + 1]);
    for (int h = 0; h < 10; ++h) for (int w = 0; w < kNW; ++w) fill(D.IWin(h, w), h, kWins[w].lo, kWins[w].hi);
}
// CB of every item from like-sign counts: tp, tm [kNfine]; cp, cm [15 cells x kNfine]
inline void CBItems(const DataL& D, const std::vector<ItemFine>& it, const double* tp, const double* tm, const double* cp, const double* cm, Vec& out) {
    const CBOpt& C = D.C;
    out.assign(D.nI, 0.);
    if (C.rule == 2) {
        // per-cell geometric mean in 10 MeV cells of the cell (WP3 MakeCBgm), uniform inside the cell
        const int nJ = kNfine / kCBreb;
        std::vector<double> g((size_t)15 * nJ, 0.);
        for (int c = 0; c < 15; ++c)
            for (int J = 0; J < nJ; ++J) {
                double a = 0, b = 0;
                for (int i = J * kCBreb; i < (J + 1) * kCBreb; ++i) { a += cp[c * kNfine + i]; b += cm[c * kNfine + i]; }
                g[(size_t)c * nJ + J] = 2. * std::sqrt(a * b) / kCBreb;   // per fine bin
            }
        for (int q = 0; q < D.nI; ++q) {
            double s = 0;
            for (const auto& p : it[q].f) s += g[(size_t)it[q].cell * nJ + p.first / kCBreb] * p.second;
            out[q] = C.k * s;
        }
        return;
    }
    const int reb = C.reb, nw = kNfine / reb;
    std::vector<double> cbW(nw, 0.), cbF(kNfine, 0.);
    for (int j = 0; j < nw; ++j) {
        double a = 0, b = 0;
        for (int i = j * reb; i < (j + 1) * reb; ++i) { a += tp[i]; b += tm[i]; }
        cbW[j] = C.rule == 1 ? a + b : 2. * std::sqrt(a * b);
    }
    for (int i = 0; i < kNfine; ++i) {   // linear interpolation between the cell centres (wp6_ana.C CBWin)
        const double x = kFineLo + (i + 0.5) * kFineW;
        const double t = (x - (kFineLo + 0.5 * reb * kFineW)) / (reb * kFineW);
        double v;
        if (t <= 0) v = cbW[0];
        else if (t >= nw - 1) v = cbW[nw - 1];
        else { const int j0 = (int)t; const double f = t - j0; v = (1 - f) * cbW[j0] + f * cbW[j0 + 1]; }
        cbF[i] = v / reb;
    }
    std::vector<double> La(kNBand, 0.), Lc((size_t)15 * kNBand, 0.);
    for (int i = 0; i < kNfine; ++i) La[BandOfFine(i)] += tp[i] + tm[i];
    for (int c = 0; c < 15; ++c) for (int i = 0; i < kNfine; ++i) Lc[(size_t)c * kNBand + BandOfFine(i)] += cp[c * kNfine + i] + cm[c * kNfine + i];
    for (int q = 0; q < D.nI; ++q) {
        double s = 0;
        const int c = it[q].cell;
        for (const auto& p : it[q].f) {
            const int b = BandOfFine(p.first);
            if (La[b] > 0) s += cbF[p.first] * p.second * Lc[(size_t)c * kNBand + b] / La[b];
        }
        out[q] = C.k * s;
    }
}
inline bool BuildData(const Rows& R, DataL& D, bool withCov = true) {
    const Layout& L = D.L;
    D.nb = L.nb();
    D.nI = 15 * D.nb + 10 * kNW;
    const int nI = D.nI;
    D.N.assign(nI, 0.);
    if (withCov) D.covN.assign((size_t)nI * nI, 0.);
    std::vector<int> fx(kNfine, -1);
    for (int i = 0; i < kNfine; ++i) {
        const double c = kFineLo + (i + 0.5) * kFineW;
        if (c >= L.e.front() && c < L.e.back()) fx[i] = (int)(std::upper_bound(L.e.begin(), L.e.end(), c) - L.e.begin()) - 1;
    }
    std::vector<double> tp(kNfine, 0.), tm(kNfine, 0.), cp((size_t)15 * kNfine, 0.), cm((size_t)15 * kNfine, 0.);
    std::vector<double> rtp, rtm, rcp, rcm;
    if (withCov) { rtp.assign((size_t)kNRepCB * kNfine, 0.); rtm = rtp; rcp.assign((size_t)kNRepCB * 15 * kNfine, 0.); rcm = rcp; }
    std::vector<int> items;
    int curK = -1, curG = -1;
    bool haveK = false;
    std::vector<double> kr(kNRepCB, 1.);
    double smp[kNSlice] = {0, 0, 0, 0, 0}, cnp[kNSlice] = {0, 0, 0, 0, 0};
    auto flush = [&]() {
        if (withCov) for (int a : items) { double* row = &D.covN[(size_t)a * nI]; for (int b : items) row[b] += 1.; }
        items.clear();
    };
    D.nUS = D.nLS = 0; D.Ntot = 0;
    for (size_t i = 0; i < R.n; ++i) {
        if (R.k[i] != curK || R.grp[i] != curG) {
            if (curK == 0) flush();
            curK = R.k[i]; curG = R.grp[i]; haveK = false;
        }
        if (!PassRow(R, i, D.S)) continue;
        const double m4 = R.m4[i];
        const int j = LSlice(L, R.mpipi[i]), ec = LClass(L, R.mee[i]);
        const int h = (j >= 0 && ec >= 0) ? 2 * j + ec : -1;
        const int j5 = (j >= 0 && R.mee[i] >= kF5gMee) ? j : -1;
        const int fi = FineBin(m4);
        if (curK == 0) {
            ++D.nUS;
            const int x = fi >= 0 ? fx[fi] : -1;
            if (h >= 0 && x >= 0) { const int a = D.IMain(h, x); items.push_back(a); D.N[a] += 1; }
            if (j5 >= 0 && x >= 0) { const int a = D.IF5(j5, x); items.push_back(a); D.N[a] += 1; }
            if (h >= 0) {
                D.Ntot += 1;
                for (int w = 0; w < kNW; ++w)
                    if (m4 >= kWins[w].lo && m4 < kWins[w].hi) { const int a = D.IWin(h, w); items.push_back(a); D.N[a] += 1; }
            }
            if (j >= 0 && m4 >= 0.40 && m4 < 0.70) { smp[j] += R.mpipi[i]; cnp[j] += 1; }
        } else {
            ++D.nLS;
            if (fi < 0) continue;
            if (withCov && !haveK) { GroupPois(kSeedCB, curK, curG, kNRepCB, kr.data()); haveK = true; }
            const bool pp = curK == 1;
            (pp ? tp : tm)[fi] += 1;
            if (h >= 0) (pp ? cp : cm)[(size_t)h * kNfine + fi] += 1;
            if (j5 >= 0) (pp ? cp : cm)[(size_t)(10 + j5) * kNfine + fi] += 1;
            if (withCov) {
                std::vector<double>& A = pp ? rtp : rtm;
                std::vector<double>& B = pp ? rcp : rcm;
                for (int r = 0; r < kNRepCB; ++r) {
                    const double w = kr[r];
                    if (w == 0) continue;
                    A[(size_t)r * kNfine + fi] += w;
                    if (h >= 0) B[((size_t)r * 15 + h) * kNfine + fi] += w;
                    if (j5 >= 0) B[((size_t)r * 15 + 10 + j5) * kNfine + fi] += w;
                }
            }
        }
    }
    if (curK == 0) flush();
    for (int j = 0; j < kNSlice; ++j) D.mpp[j] = cnp[j] > 0 ? smp[j] / cnp[j] : 0.;
    std::vector<ItemFine> it;
    ItemDefs(D, it);
    CBItems(D, it, tp.data(), tm.data(), cp.data(), cm.data(), D.CB);
    if (withCov) {
        std::vector<Vec> rep(kNRepCB);
        Vec mean(nI, 0.);
        for (int r = 0; r < kNRepCB; ++r) {
            CBItems(D, it, &rtp[(size_t)r * kNfine], &rtm[(size_t)r * kNfine], &rcp[(size_t)r * 15 * kNfine], &rcm[(size_t)r * 15 * kNfine], rep[r]);
            for (int a = 0; a < nI; ++a) mean[a] += rep[r][a] / kNRepCB;
        }
        D.covCB.assign((size_t)nI * nI, 0.);
        Vec d(nI);
        for (int r = 0; r < kNRepCB; ++r) {
            for (int a = 0; a < nI; ++a) d[a] = rep[r][a] - mean[a];
            for (int a = 0; a < nI; ++a) {
                if (d[a] == 0) continue;
                double* row = &D.covCB[(size_t)a * nI];
                const double da = d[a] / (kNRepCB - 1);
                for (int b = 0; b < nI; ++b) row[b] += da * d[b];
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------------------
// 4. simulation templates in a layout (expected raw data rows = nb x kNbToRaw x scale)
// ---------------------------------------------------------------------------------------------------------------
struct Tmpl {
    int nb = 0;
    Vec s[kNSg], vs[kNSg];                              // main [h nb + x]
    Vec gd, vgd, gc, vgc, gdo, vgdo, gco, vgco;         // G_D, G_C and their omega/rho0 parts
    Vec s5[kNSg], vs5[kNSg], gd5, vgd5, gc5, vgc5;      // F5g [j nb + x]
    Vec ws[kNSg];                                       // aligned windows per cell [h kNW + w]; NaN if not available
    double wa[kNSg][kNW];                               // all cells, all M(ee); NaN if not available
};
inline void AllocT(Tmpl& T, int nb) {
    T.nb = nb;
    const size_t n = (size_t)10 * nb, n5 = (size_t)5 * nb;
    for (int t = 0; t < kNSg; ++t) { T.s[t].assign(n, 0.); T.vs[t].assign(n, 0.); T.s5[t].assign(n5, 0.); T.vs5[t].assign(n5, 0.); T.ws[t].assign(10 * kNW, NAN); }
    T.gd.assign(n, 0.); T.vgd = T.gd; T.gc = T.gd; T.vgc = T.gd; T.gdo = T.gd; T.vgdo = T.gd; T.gco = T.gd; T.vgco = T.gd;
    T.gd5.assign(n5, 0.); T.vgd5 = T.gd5; T.gc5 = T.gd5; T.vgc5 = T.gd5;
    for (int t = 0; t < kNSg; ++t) for (int w = 0; w < kNW; ++w) T.wa[t][w] = NAN;
}
inline void ScaleT(Tmpl& T, double f) {
    auto sc = [&](Vec& v, double g) { for (auto& x : v) x *= g; };
    for (int t = 0; t < kNSg; ++t) { sc(T.s[t], f); sc(T.vs[t], f * f); sc(T.s5[t], f); sc(T.vs5[t], f * f); sc(T.ws[t], f); for (int w = 0; w < kNW; ++w) T.wa[t][w] *= f; }
    sc(T.gd, f); sc(T.gc, f); sc(T.gdo, f); sc(T.gco, f); sc(T.gd5, f); sc(T.gc5, f);
    sc(T.vgd, f * f); sc(T.vgc, f * f); sc(T.vgdo, f * f); sc(T.vgco, f * f); sc(T.vgd5, f * f); sc(T.vgc5, f * f);
}
// window sums of aligned windows from the main bins (window inside the layout range, edges on bin edges)
inline void WinFromBins(const Layout& L, Tmpl& T) {
    const Vec& e = L.e;
    const int nb = L.nb();
    for (int w = 0; w < kNW; ++w) {
        if (!kWins[w].aligned) continue;
        const double lo = kWins[w].lo, hi = kWins[w].hi;
        if (lo < e.front() - 1e-9 || hi > e.back() + 1e-9) continue;
        bool okl = false, okh = false;
        for (double x : e) { if (std::fabs(x - lo) < 1e-7) okl = true; if (std::fabs(x - hi) < 1e-7) okh = true; }
        if (!okl || !okh) continue;
        for (int t = 0; t < kNSg; ++t) {
            double tot = 0;
            for (int h = 0; h < 10; ++h) {
                double s = 0;
                for (int x = 0; x < nb; ++x) if (e[x] >= lo - 1e-9 && e[x + 1] <= hi + 1e-9) s += T.s[t][h * nb + x];
                T.ws[t][h * kNW + w] = s; tot += s;
            }
            T.wa[t][w] = tot;
        }
    }
}
// F10 / F5g / win / tot of sim/<sel>__<var> (standard layout only)
inline bool LoadTmplF10(TFile* fh, const char* sel, const char* var, Tmpl& T, FILE* log = nullptr) {
    const Layout L = StdLayout();
    const int nb = L.nb();
    AllocT(T, nb);
    const double k = kNbToRaw;
    auto get = [&](int comp, const std::string& lay) -> TH1* { return allio::Sim(fh, sel, var, comp, lay); };
    auto rd = [&](TH1* x, Vec& c, Vec& v, size_t off, double sign = 1.) {
        for (int b = 1; b <= nb; ++b) { c[off + b - 1] += sign * x->GetBinContent(b) * k; v[off + b - 1] += std::pow(x->GetBinError(b) * k, 2); }
    };
    for (int h = 0; h < 10; ++h) {
        const std::string lay = F10Name(h / 2, h % 2);
        for (int t = 0; t < kNSg; ++t) { TH1* x = get(kSigComp[t], lay); if (!x) return false; rd(x, T.s[t], T.vs[t], (size_t)h * nb); delete x; }
        TH1* a = get(kGD, lay); TH1* b = get(kGconv, lay); if (!a || !b) return false;
        rd(a, T.gd, T.vgd, (size_t)h * nb); rd(b, T.gc, T.vgc, (size_t)h * nb); delete a; delete b;
        for (int c : {kGpi0Dom, kGvecom, kGothom}) { TH1* x = get(c, lay); if (!x) return false; rd(x, T.gdo, T.vgdo, (size_t)h * nb); delete x; }
        TH1* o = get(kGconvom, lay); if (!o) return false; rd(o, T.gco, T.vgco, (size_t)h * nb); delete o;
    }
    for (int j = 0; j < 5; ++j) {
        const std::string lay = F5gName(j);
        for (int t = 0; t < kNSg; ++t) { TH1* x = get(kSigComp[t], lay); if (!x) return false; rd(x, T.s5[t], T.vs5[t], (size_t)j * nb); delete x; }
        TH1* a = get(kGD, lay); TH1* b = get(kGconv, lay); if (!a || !b) return false;
        rd(a, T.gd5, T.vgd5, (size_t)j * nb); rd(b, T.gc5, T.vgc5, (size_t)j * nb); delete a; delete b;
    }
    WinFromBins(L, T);
    for (int t = 0; t < kNSg; ++t) {
        TH1* w = get(kSigComp[t], "win"); TH1* tt = get(kSigComp[t], "tot"); if (!w || !tt) return false;
        T.wa[t][kWB] = w->GetBinContent(kWB + 1) * k;       // = sum of the F10 cells (checked by WP6)
        T.wa[t][kWP] = w->GetBinContent(alldefs::kWP + 1) * k;
        T.wa[t][kWSl] = w->GetBinContent(kWSBlo + 1) * k;
        T.wa[t][kWSh] = w->GetBinContent(kWSBhi + 1) * k;
        T.wa[t][kWtot] = tt->GetBinContent(1) * k;
        delete w; delete tt;
    }
    if (log) fprintf(log, "  templates F10 sim/%s__%s loaded\n", sel, var);
    return true;
}
// G3 of sim/<sel>__nom summed into the layout L (edges on the 2.5 MeV grid of 0.40-0.85; slice edges on the G3 grid or
// shared linearly; M(ee) split on a G3 edge). Variances: sums of the G3 cell variances (approximation, flagged).
// Non-aligned all-cell windows (P, side bands, tot) from win/tot of the same directory (x the K0S-veto ratio).
inline bool LoadTmplG3(TFile* fh, const char* sel, const Layout& L, Tmpl& T, FILE* log = nullptr) {
    const int nb = L.nb();
    AllocT(T, nb);
    const double k = kNbToRaw;
    const Vec ye = G3PiPiEdges(), ze = G3MeeEdges();
    const int nx = kG3nx, ny = (int)ye.size() - 1, nz = (int)ze.size() - 1;
    // x-bin -> layout bin
    std::vector<int> xb(nx + 1, -1);
    for (int x = 1; x <= nx; ++x) {
        const double c = kG3xlo + (x - 0.5) * 0.0025;
        if (c >= L.e.front() && c < L.e.back()) xb[x] = (int)(std::upper_bound(L.e.begin(), L.e.end(), c) - L.e.begin()) - 1;
    }
    // y-bin -> slice weights (simulation edges x (1 + simEdgeRel)); veto bin 0.485-0.507
    std::vector<std::array<double, 5>> wy(ny + 1);
    double sedge[4]; for (int q = 0; q < 4; ++q) sedge[q] = L.se[q] * (1. + L.simEdgeRel);
    for (int y = 1; y <= ny; ++y) {
        const double a = ye[y - 1], b = ye[y];
        for (int j = 0; j < 5; ++j) wy[y][j] = 0.;
        if (L.veto && std::fabs(a - 0.485) < 1e-9) continue;
        for (int j = 0; j < 5; ++j) {
            const double lo = j == 0 ? -1e9 : sedge[j - 1], hi = j == 4 ? 1e9 : sedge[j];
            const double o = std::min(b, hi) - std::max(a, lo);
            if (o > 0) wy[y][j] = o / (b - a);
        }
    }
    std::vector<int> zc(nz + 1, -1); std::vector<bool> z5(nz + 1, false);
    for (int z = 1; z <= nz; ++z) { const double c = 0.5 * (ze[z - 1] + std::min(ze[z], 0.1)); zc[z] = c < L.eSplit ? 0 : 1; z5[z] = ze[z - 1] >= kF5gMee - 1e-9; }
    auto add3 = [&](int comp, Vec* c, Vec* v, Vec* c5, Vec* v5, double sign) -> bool {
        TH3* g = (TH3*)allio::Sim(fh, sel, "nom", comp, "G3");
        if (!g) return false;
        for (int x = 1; x <= nx; ++x) {
            if (xb[x] < 0) continue;
            for (int y = 1; y <= ny; ++y)
                for (int z = 1; z <= nz; ++z) {
                    const double val = g->GetBinContent(x, y, z) * k, var = std::pow(g->GetBinError(x, y, z) * k, 2);
                    if (val == 0 && var == 0) continue;
                    for (int j = 0; j < 5; ++j) {
                        const double wgt = wy[y][j];
                        if (wgt <= 0) continue;
                        const int h = 2 * j + zc[z];
                        if (c) { (*c)[(size_t)h * nb + xb[x]] += sign * wgt * val; (*v)[(size_t)h * nb + xb[x]] += wgt * wgt * var; }
                        if (c5 && z5[z]) { (*c5)[(size_t)j * nb + xb[x]] += sign * wgt * val; (*v5)[(size_t)j * nb + xb[x]] += wgt * wgt * var; }
                    }
                }
        }
        delete g;
        return true;
    };
    for (int t = 0; t < kNSg; ++t) if (!add3(kSigComp[t], &T.s[t], &T.vs[t], &T.s5[t], &T.vs5[t], 1.)) return false;
    if (!add3(kGD, &T.gd, &T.vgd, &T.gd5, &T.vgd5, 1.)) return false;
    if (!add3(kGconv, &T.gc, &T.vgc, &T.gc5, &T.vgc5, 1.)) return false;
    for (int c : {kGpi0Dom, kGvecom, kGothom}) if (!add3(c, &T.gdo, &T.vgdo, nullptr, nullptr, 1.)) return false;
    if (!add3(kGconvom, &T.gco, &T.vgco, nullptr, nullptr, 1.)) return false;
    WinFromBins(L, T);
    // aligned windows outside the layout range: directly from G3 (no layout binning needed)
    for (int t = 0; t < kNSg; ++t) {
        TH3* g = (TH3*)allio::Sim(fh, sel, "nom", kSigComp[t], "G3");
        for (int w = 0; w < kNW; ++w) {
            if (!kWins[w].aligned || T.wa[t][w] == T.wa[t][w]) continue;
            if (kWins[w].lo < kG3xlo - 1e-9 || kWins[w].hi > kG3xhi + 1e-9) continue;
            double tot = 0;
            for (int h = 0; h < 10; ++h) T.ws[t][h * kNW + w] = 0.;
            for (int x = 1; x <= nx; ++x) {
                const double c = kG3xlo + (x - 0.5) * 0.0025;
                if (!(c >= kWins[w].lo && c < kWins[w].hi)) continue;
                for (int y = 1; y <= ny; ++y) for (int z = 1; z <= nz; ++z) for (int j = 0; j < 5; ++j) {
                    if (wy[y][j] <= 0) continue;
                    const double v = wy[y][j] * g->GetBinContent(x, y, z) * k;
                    T.ws[t][(2 * j + zc[z]) * kNW + w] += v; tot += v;
                }
            }
            T.wa[t][w] = tot;
        }
        // veto ratio for the non-aligned all-cell windows (fractional 2.5 MeV bins)
        double rat[kNW];
        for (int w = 0; w < kNW; ++w) rat[w] = 1.;
        if (L.veto) {
            for (int w : {(int)kWP, (int)kWSl, (int)kWSh, (int)kWtot}) {
                double all = 0, kept = 0;
                for (int x = 1; x <= nx; ++x) {
                    const double a = kG3xlo + (x - 1) * 0.0025, b = a + 0.0025;
                    const double o = w == kWtot ? 1. : std::max(0., std::min(b, kWins[w].hi) - std::max(a, kWins[w].lo)) / 0.0025;
                    if (o <= 0) continue;
                    for (int y = 1; y <= ny; ++y) for (int z = 1; z <= nz; ++z) {
                        const double v = g->GetBinContent(x, y, z) * o;
                        all += v;
                        if (!(std::fabs(ye[y - 1] - 0.485) < 1e-9)) kept += v;
                    }
                }
                rat[w] = all > 0 ? kept / all : 1.;
            }
        }
        delete g;
        TH1* wh = allio::Sim(fh, sel, "nom", kSigComp[t], "win"); TH1* tt = allio::Sim(fh, sel, "nom", kSigComp[t], "tot");
        if (!wh || !tt) return false;
        T.wa[t][kWP] = wh->GetBinContent(alldefs::kWP + 1) * k * rat[kWP];
        T.wa[t][kWSl] = wh->GetBinContent(kWSBlo + 1) * k * rat[kWSl];
        T.wa[t][kWSh] = wh->GetBinContent(kWSBhi + 1) * k * rat[kWSh];
        T.wa[t][kWtot] = tt->GetBinContent(1) * k * rat[kWtot];
        delete wh; delete tt;
    }
    if (log) fprintf(log, "  templates G3 sim/%s__nom in layout %s\n", sel, L.tag.c_str());
    return true;
}
// the 100 WP0 bootstrap replicas of sim/<sel>__nom (standard layout); variances = those of the full sample (wp5::LoadReplicas)
inline bool LoadTmplF10Rep(TFile* fh, const char* sel, const Tmpl& Tf, int nrep, std::vector<Tmpl>& R) {
    const Layout L = StdLayout();
    const int nb = L.nb();
    const double k = kNbToRaw;
    R.assign(nrep, Tmpl());
    for (auto& T : R) AllocT(T, nb);
    for (auto& T : R) {
        for (int t = 0; t < kNSg; ++t) { T.vs[t] = Tf.vs[t]; T.vs5[t] = Tf.vs5[t]; }
        T.vgd = Tf.vgd; T.vgc = Tf.vgc; T.vgdo = Tf.vgdo; T.vgco = Tf.vgco; T.vgd5 = Tf.vgd5; T.vgc5 = Tf.vgc5;
    }
    auto rd2 = [&](int comp, const std::string& lay, std::function<Vec&(Tmpl&)> dst, size_t off, bool opt) -> bool {
        TH2* x = (TH2*)fh->Get((std::string("sim/") + sel + "__nom/" + CompTag(comp) + "_" + lay + "_rep").c_str());
        if (!x) { if (!opt) printf("wp7: replica %s %s missing\n", CompTag(comp), lay.c_str()); return opt; }
        for (int r = 0; r < nrep; ++r) { Vec& v = dst(R[r]); for (int b = 1; b <= nb; ++b) v[off + b - 1] += x->GetBinContent(b, r + 1) * k; }
        return true;
    };
    for (int h = 0; h < 10; ++h) {
        const std::string lay = F10Name(h / 2, h % 2);
        const size_t off = (size_t)h * nb;
        for (int t = 0; t < kNSg; ++t) if (!rd2(kSigComp[t], lay, [t](Tmpl& T) -> Vec& { return T.s[t]; }, off, false)) return false;
        if (!rd2(kGD, lay, [](Tmpl& T) -> Vec& { return T.gd; }, off, false)) return false;
        if (!rd2(kGconv, lay, [](Tmpl& T) -> Vec& { return T.gc; }, off, false)) return false;
        for (int c : {kGpi0Dom, kGvecom, kGothom}) rd2(c, lay, [](Tmpl& T) -> Vec& { return T.gdo; }, off, true);
        rd2(kGconvom, lay, [](Tmpl& T) -> Vec& { return T.gco; }, off, true);
    }
    for (int j = 0; j < 5; ++j) {
        const std::string lay = F5gName(j);
        const size_t off = (size_t)j * nb;
        for (int t = 0; t < kNSg; ++t) if (!rd2(kSigComp[t], lay, [t](Tmpl& T) -> Vec& { return T.s5[t]; }, off, false)) return false;
        if (!rd2(kGD, lay, [](Tmpl& T) -> Vec& { return T.gd5; }, off, false)) return false;
        if (!rd2(kGconv, lay, [](Tmpl& T) -> Vec& { return T.gc5; }, off, false)) return false;
    }
    for (auto& T : R) WinFromBins(L, T);
    for (int t = 0; t < kNSg; ++t) {
        TH2* w = (TH2*)fh->Get((std::string("sim/") + sel + "__nom/" + CompTag(kSigComp[t]) + "_win_rep").c_str());
        TH2* tt = (TH2*)fh->Get((std::string("sim/") + sel + "__nom/" + CompTag(kSigComp[t]) + "_tot_rep").c_str());
        if (!w || !tt) return false;
        for (int r = 0; r < nrep; ++r) {
            R[r].wa[t][kWB] = w->GetBinContent(kWB + 1, r + 1) * k;
            R[r].wa[t][kWP] = w->GetBinContent(alldefs::kWP + 1, r + 1) * k;
            R[r].wa[t][kWSl] = w->GetBinContent(kWSBlo + 1, r + 1) * k;
            R[r].wa[t][kWSh] = w->GetBinContent(kWSBhi + 1, r + 1) * k;
            R[r].wa[t][kWtot] = tt->GetBinContent(1, r + 1) * k;
        }
    }
    return true;
}
// S9 variant: factorised signal templates T(j, e, x) = total x f(e) x g(j) x h(x | j), per class; F5g: f(>= 20 MeV) per slice
inline void FactoriseSignal(Tmpl& T, const Layout& L) {
    const int nb = T.nb;
    for (int t = 0; t < kNSg; ++t) {
        double tot = 0, fe[2] = {0, 0}, gj[5] = {0, 0, 0, 0, 0};
        Vec hx((size_t)5 * nb, 0.), vx((size_t)5 * nb, 0.);
        for (int h = 0; h < 10; ++h) for (int x = 0; x < nb; ++x) {
            const double c = T.s[t][h * nb + x];
            tot += c; fe[h % 2] += c; gj[h / 2] += c; hx[(h / 2) * nb + x] += c; vx[(h / 2) * nb + x] += T.vs[t][h * nb + x];
        }
        if (!(tot > 0)) continue;
        for (int h = 0; h < 10; ++h) {
            const int j = h / 2;
            const double f = fe[h % 2] / tot;
            for (int x = 0; x < nb; ++x) {
                T.s[t][h * nb + x] = f * hx[j * nb + x];          // = tot f(e) g(j) h(x | j)
                T.vs[t][h * nb + x] = f * f * vx[j * nb + x];
            }
        }
        // F5g: the >= 20 MeV share per slice from the F5g sums, shape h(x | j)
        for (int j = 0; j < 5; ++j) {
            double s5 = 0, sj = 0;
            for (int x = 0; x < nb; ++x) { s5 += T.s5[t][j * nb + x]; sj += hx[j * nb + x]; }
            const double f = sj > 0 ? s5 / sj : 0.;
            for (int x = 0; x < nb; ++x) { T.s5[t][j * nb + x] = f * hx[j * nb + x]; T.vs5[t][j * nb + x] = f * f * vx[j * nb + x]; }
        }
    }
    WinFromBins(L, T);
}

// ---------------------------------------------------------------------------------------------------------------
// 5. smoothing of the non-eta templates
// ---------------------------------------------------------------------------------------------------------------
const double kGL4x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
const double kGL4w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
inline double Tsm(double m) { return (m - 0.55) / 0.15; }
struct SD {          // smoothed density
    int kind = -1;   // -1 zero; 0 wp4::Smooth; 1 exp(pol n) in t = (M4 - 0.55)/0.15 from 'start' (WP3 form); 2 raw histogram (piecewise flat)
    wp4::Smooth s4;
    double start = 1e9; int ord = 3; double p[6] = {0, 0, 0, 0, 0, 0}; int dom = -1;
    const Vec* raw = nullptr; const Vec* e = nullptr; size_t off = 0;
    double Dens(double m) const {
        if (kind == 0) return s4.Dens(m);
        if (kind == 1) {
            if (m < start) return 0.;
            const double t = Tsm(m);
            double a = 0; for (int k = ord; k >= 0; --k) a = a * t + p[k];
            return std::exp(std::min(700., a));
        }
        if (kind == 2) {
            if (!(m >= e->front() && m < e->back())) return 0.;
            const int x = (int)(std::upper_bound(e->begin(), e->end(), m) - e->begin()) - 1;
            return (*raw)[off + x] / ((*e)[x + 1] - (*e)[x]);
        }
        return 0.;
    }
};
// WP3 smoothing (wp6_ana.C w3::Smooth, copied and generalised to order n <= 4): Poisson likelihood of the expected rows,
// fitted from the first bin with content (dom) or from 'dom' given; replicas start from the nominal parameters.
struct SmR { bool ok = false; int dom = 0; int ord = 3; double p[6] = {0, 0, 0, 0, 0, 0}; Vec f; double chi2 = 0; int ndf = 0; };
inline void BinMomN(double lo, double hi, const double* p, int ord, double* mom, int nm) {
    for (int k = 0; k < nm; ++k) mom[k] = 0;
    const double h = 0.5 * (hi - lo), c = 0.5 * (hi + lo);
    for (int g = 0; g < 4; ++g) {
        const double m = c + h * kGL4x[g], t = Tsm(m);
        double a = 0; for (int k = ord; k >= 0; --k) a = a * t + p[k];
        const double e = std::exp(std::min(700., a)) * h * kGL4w[g];
        double tk = 1;
        for (int k = 0; k < nm; ++k) { mom[k] += tk * e; tk *= t; }
    }
}
inline SmR SmoothW3(const Vec& y, const Vec& e2, const Vec& lo, const Vec& hi, int ord = 3, int dom = -1, const double* pstart = nullptr) {
    SmR r; r.ord = ord;
    const int n = (int)y.size(), np = ord + 1, nm = 2 * ord + 1;
    if (dom < 0) { dom = n; for (int i = 0; i < n; ++i) if (y[i] > 0) { dom = i; break; } }
    r.dom = dom;
    r.f.assign(n, 0.);
    double tot = 0; for (int i = dom; i < n; ++i) tot += y[i];
    if (dom >= n - (ord + 1) || !(tot > 0)) return r;
    double p[6] = {0, 0, 0, 0, 0, 0};
    if (pstart) for (int k = 0; k < np; ++k) p[k] = pstart[k];
    else {
        bla::Mat A = bla::Zero(np, np); Vec b(np, 0.);
        for (int i = dom; i < n; ++i) {
            if (!(y[i] > 0)) continue;
            const double t = Tsm(0.5 * (lo[i] + hi[i])), w = y[i], z = std::log(y[i] / (hi[i] - lo[i]));
            double tk[6]; tk[0] = 1; for (int k = 1; k < np; ++k) tk[k] = tk[k - 1] * t;
            for (int k = 0; k < np; ++k) { b[k] += w * tk[k] * z; for (int l = 0; l < np; ++l) A[k][l] += w * tk[k] * tk[l]; }
        }
        for (int k = 0; k < np; ++k) A[k][k] += 1e-9;
        if (bla::CholInv(A)) for (int k = 0; k < np; ++k) { p[k] = 0; for (int l = 0; l < np; ++l) p[k] += A[k][l] * b[l]; }
        else { p[0] = std::log(tot / (hi[n - 1] - lo[dom])); }
    }
    auto F = [&](const double* q) {
        double s = 0, mom[11];
        for (int i = dom; i < n; ++i) {
            BinMomN(lo[i], hi[i], q, ord, mom, 1);
            const double f = std::max(mom[0], 1e-300);
            s += f - (y[i] > 0 ? y[i] * std::log(f) : 0.);
        }
        return s;
    };
    double Fc = F(p);
    for (int it = 0; it < 100; ++it) {
        bla::Mat H = bla::Zero(np, np); Vec g(np, 0.);
        double mom[11];
        for (int i = dom; i < n; ++i) {
            BinMomN(lo[i], hi[i], p, ord, mom, nm);
            const double f = std::max(mom[0], 1e-300), a = 1.0 - y[i] / f, c2 = y[i] / (f * f);
            for (int k = 0; k < np; ++k) {
                g[k] += a * mom[k];
                for (int l = 0; l < np; ++l) H[k][l] += a * mom[k + l] + c2 * mom[k] * mom[l];
            }
        }
        bla::Mat Hi = H;
        double lam = 0;
        bool inv = false;
        for (int tr = 0; tr < 20 && !inv; ++tr) {
            Hi = H;
            for (int k = 0; k < np; ++k) Hi[k][k] += lam * std::max(1e-12, std::fabs(H[k][k]));
            inv = bla::CholInv(Hi);
            lam = lam == 0 ? 1e-6 : lam * 10;
        }
        if (!inv) break;
        double d[6];
        for (int k = 0; k < np; ++k) { d[k] = 0; for (int l = 0; l < np; ++l) d[k] -= Hi[k][l] * g[l]; }
        double st = 1.0, Fn = 0, pn[6] = {0, 0, 0, 0, 0, 0};
        for (int h = 0; h < 40; ++h) {
            for (int k = 0; k < np; ++k) pn[k] = p[k] + st * d[k];
            Fn = F(pn);
            if (Fn <= Fc + 1e-12 * std::fabs(Fc)) break;
            st *= 0.5;
        }
        const double dF = Fc - Fn;
        if (Fn <= Fc) { for (int k = 0; k < np; ++k) p[k] = pn[k]; Fc = Fn; }
        if (dF >= 0 && dF < 1e-10 * std::max(1.0, std::fabs(Fc))) break;
        if (st < 1e-10) break;
    }
    for (int k = 0; k < np; ++k) r.p[k] = p[k];
    double sy = 0, se = 0;
    for (int i = dom; i < n; ++i) { double mom[11]; BinMomN(lo[i], hi[i], p, ord, mom, 1); r.f[i] = mom[0]; sy += y[i]; se += e2[i]; }
    const double wbar = sy > 0 ? se / sy : 1.0;
    r.chi2 = 0; r.ndf = n - dom - np;
    for (int i = dom; i < n; ++i) if (r.f[i] > 0) r.chi2 += (y[i] - r.f[i]) * (y[i] - r.f[i]) / (wbar * r.f[i]);
    r.ok = true;
    return r;
}
// smoothing sets: 15 histograms (10 main, 5 F5g) x 4 templates (G_D, G_C, G_D om part, G_C om part)
struct SmSet {
    SD d[15][4];
    Vec I[4];        // smoothed bin integrals [h nb + x], h = 0..14
    double chi2[15][4];
};
// mode 0 = WP4/WP5 rule (wp4::SmoothFit from 0.400, N_eff order rule, F5g parent = main E2 of the slice; ord 4: own fit,
// same domain and rules); mode 1 = WP3 rule (exp(pol ord) from the first bin with content; F5g: same rule);
// raw = no smoothing (S4 "raw BB-lite": the bin contents themselves)
inline void SmoothT(const Tmpl& T, const Layout& L, int mode, int ord, bool raw, SmSet& S, bool omega, const SmSet* ref = nullptr) {
    const int nb = T.nb;
    const Vec& e = L.e;
    Vec lo(nb), hi(nb);
    for (int x = 0; x < nb; ++x) { lo[x] = e[x]; hi[x] = e[x + 1]; }
    for (int q = 0; q < 4; ++q) S.I[q].assign((size_t)15 * nb, 0.);
    for (int h = 0; h < 15; ++h)
        for (int q = 0; q < 4; ++q) {
            S.chi2[h][q] = 0;
            S.d[h][q] = SD();
            if (q >= 2 && (h >= 10 || !omega)) continue;
            const Vec* src; const Vec* var; size_t off;
            if (h < 10) { src = q == 0 ? &T.gd : q == 1 ? &T.gc : q == 2 ? &T.gdo : &T.gco; var = q == 0 ? &T.vgd : q == 1 ? &T.vgc : q == 2 ? &T.vgdo : &T.vgco; off = (size_t)h * nb; }
            else { src = q == 0 ? &T.gd5 : &T.gc5; var = q == 0 ? &T.vgd5 : &T.vgc5; off = (size_t)(h - 10) * nb; }
            SD& d = S.d[h][q];
            if (raw) {
                d.kind = 2; d.raw = src; d.e = &L.e; d.off = off;
                for (int x = 0; x < nb; ++x) S.I[q][(size_t)h * nb + x] = (*src)[off + x];
                continue;
            }
            if (mode == 0) {
                Vec c(nb + 1, 0.), v(nb + 1, 0.);
                for (int x = 0; x < nb; ++x) { c[x + 1] = (*src)[off + x]; v[x + 1] = (*var)[off + x]; }
                const wp4::Smooth* par = (h >= 10 && S.d[2 * (h - 10) + 1][q].kind == 0) ? &S.d[2 * (h - 10) + 1][q].s4 : nullptr;
                if (ord <= 3) {
                    d.kind = 0;
                    d.s4 = wp4::SmoothFit(c, v, e, par, ord == 3 ? -1 : ord);
                    for (int x = 0; x < nb; ++x) S.I[q][(size_t)h * nb + x] = d.s4.Int(e[x], e[x + 1]);
                } else {
                    // exp(pol 4) from 0.400 with the WP4 N_eff rule (>= 50: 4, 20-50: 2, < 20 in F5g: parent shape)
                    double C = 0, V = 0; for (int x = 1; x <= nb; ++x) { C += c[x]; V += v[x]; }
                    const double neff = V > 0 ? C * C / V : 0;
                    if (!(C > 0)) continue;
                    if (neff < 50) {
                        d.kind = 0;
                        d.s4 = wp4::SmoothFit(c, v, e, par, -1);
                        for (int x = 0; x < nb; ++x) S.I[q][(size_t)h * nb + x] = d.s4.Int(e[x], e[x + 1]);
                    } else {
                        Vec y(nb), vv(nb); for (int x = 0; x < nb; ++x) { y[x] = c[x + 1]; vv[x] = v[x + 1]; }
                        SmR r = SmoothW3(y, vv, lo, hi, 4, 0);
                        d.kind = 1; d.start = lo[0]; d.ord = 4; for (int k = 0; k < 6; ++k) d.p[k] = r.p[k];
                        for (int x = 0; x < nb; ++x) S.I[q][(size_t)h * nb + x] = r.f[x];
                    }
                }
            } else {
                Vec y(nb), vv(nb); for (int x = 0; x < nb; ++x) { y[x] = (*src)[off + x]; vv[x] = (*var)[off + x]; }
                const SD* rd = ref ? &ref->d[h][q] : nullptr;
                SmR r = (rd && rd->kind == 1) ? SmoothW3(y, vv, lo, hi, ord, rd->dom, rd->p) : SmoothW3(y, vv, lo, hi, ord);
                if (!r.ok) continue;
                d.kind = 1; d.ord = ord; d.dom = r.dom; d.start = r.dom < nb ? lo[r.dom] : 1e9; for (int k = 0; k < 6; ++k) d.p[k] = r.p[k];
                for (int x = 0; x < nb; ++x) S.I[q][(size_t)h * nb + x] = r.f[x];
                S.chi2[h][q] = r.ndf > 0 ? r.chi2 / r.ndf : 0.;
            }
        }
}
// integral over [lo, hi) of density x tilt(m), with the tilt form of the model (0: 1 + b x + q x^2, 1: exp(b x))
inline double IntTiltSD(const SD& d, double lo, double hi, double b, double q, int form) {
    if (d.kind < 0 || !(hi > lo)) return 0.;
    if (d.kind == 1) lo = std::max(lo, d.start);
    if (!(hi > lo)) return 0.;
    const int n = std::max(1, (int)std::ceil((hi - lo) / 0.0025 - 1e-9));
    const double w = (hi - lo) / n;
    double tot = 0;
    for (int k = 0; k < n; ++k) {
        const double c = lo + (k + 0.5) * w;
        for (int g = 0; g < 8; ++g) {
            const double m = c + 0.5 * w * wp4::kGLx[g];
            const double x = TiltX(m);
            const double tl = form == 1 ? std::exp(b * x) : 1. + b * x + q * x * x;
            tot += d.Dens(m) * tl * wp4::kGLw[g] * 0.5 * w;
        }
    }
    return tot;
}

// ---------------------------------------------------------------------------------------------------------------
// 6. fit bins, model, fitter (copy of wp5::Eval / wp5::Fit with extra model terms and generic covariance items)
// ---------------------------------------------------------------------------------------------------------------
struct GBins {
    int n = 0;
    Vec y, cb, vcb, x, gd, gc, gdo, gco, vbd, vbc, fix, lo, hi, yj;
    std::vector<int> j, e, item;
    Vec s[kNSg], vs[kNSg];
};
struct Model {
    int tilt = 0;            // 0: 1 + b x [+ q x^2]; 1: exp(b x)
    int iOm = -1;            // omega/rho0 scale c_om (background = [G_D + (c_om - 1) G_D,om] + kB [G_C + (c_om - 1) G_C,om])
    bool kBtie = false;      // kappa_B = kappa (parameter iKap)
    int iKap = -1;
    int iC0 = -1, iC1 = -1;  // c_j = c0 (1 + c1 y_j)
    bool bbBkg = false;      // BB-lite variance of the (raw) background templates in V
};
// bins of the main layout: classes useE (bit 1 E1, bit 2 E2), M4 range [mlo, mhi); keep(j, e, lo, hi) optional extra filter;
// or the F5g layout (f5 = true)
inline void MakeGBins(const DataL& D, const Tmpl& T, const SmSet& S, bool f5, int useE, double mlo, double mhi, GBins& B,
                      const std::function<bool(int, int, double, double)>& keep = nullptr,
                      const std::function<double(int, int, int)>& fixfn = nullptr) {
    B = GBins();
    const Vec& e = D.L.e;
    const int nb = D.nb;
    const int nh = f5 ? 5 : 10;
    for (int q = 0; q < nh; ++q) {
        const int jj = f5 ? q : q / 2, ec = f5 ? 1 : q % 2;
        if (!f5 && !((useE >> ec) & 1)) continue;
        for (int x = 0; x < nb; ++x) {
            const bool inR = e[x] >= mlo - 1e-9 && e[x + 1] <= mhi + 1e-9;
            const bool kp = keep ? keep(jj, ec, e[x], e[x + 1]) : false;
            if (!inR && !kp) continue;
            if (keep && !kp && !inR) continue;
            const int it = f5 ? D.IF5(q, x) : D.IMain(q, x);
            const int hs = f5 ? 10 + q : q;   // smoothing index
            const size_t ti = (size_t)q * nb + x;
            // CB variance in the weights: replica variance with 1/n as wp5::LoadData (the covariance itself is stored with 1/(n-1))
            B.y.push_back(D.N[it]); B.cb.push_back(D.CB[it]); B.vcb.push_back(D.covCB.empty() ? 0. : D.covCB[(size_t)it * D.nI + it] * (kNRepCB - 1.) / kNRepCB);
            B.x.push_back(TiltX(0.5 * (e[x] + e[x + 1]))); B.lo.push_back(e[x]); B.hi.push_back(e[x + 1]);
            B.gd.push_back(S.I[0][(size_t)hs * nb + x]); B.gc.push_back(S.I[1][(size_t)hs * nb + x]);
            B.gdo.push_back(f5 ? 0. : S.I[2][(size_t)hs * nb + x]); B.gco.push_back(f5 ? 0. : S.I[3][(size_t)hs * nb + x]);
            B.vbd.push_back(f5 ? T.vgd5[ti] : T.vgd[ti]); B.vbc.push_back(f5 ? T.vgc5[ti] : T.vgc[ti]);
            B.fix.push_back(fixfn ? fixfn(f5 ? 1 : 0, q, x) : 0.);
            B.yj.push_back((D.mpp[jj] - 0.38) / 0.05);
            B.j.push_back(jj); B.e.push_back(ec); B.item.push_back(it);
            for (int t = 0; t < kNSg; ++t) { B.s[t].push_back(f5 ? T.s5[t][ti] : T.s[t][ti]); B.vs[t].push_back(f5 ? T.vs5[t][ti] : T.vs[t][ti]); }
            ++B.n;
        }
    }
}
inline void EvalX(const GBins& B, const Spec& S, const Model& M, const Vec& th, Vec& mu, Mat* J, Vec* V, Vec* Vd) {
    const int n = B.n, np = S.np();
    mu.assign(n, 0.);
    if (J) J->assign(n, Vec(np, 0.));
    if (V) V->assign(n, 0.);
    if (Vd) Vd->assign(n, 0.);
    const double kB = M.kBtie ? th[M.iKap] : (S.iKB >= 0 ? th[S.iKB] : S.kBfix);
    const double com = M.iOm >= 0 ? th[M.iOm] : 1.;
    for (int i = 0; i < n; ++i) {
        const int e = B.e[i], j = B.j[i];
        const double k = S.iK[e] >= 0 ? th[S.iK[e]] : S.kFix[e];
        double c;
        if (M.iC0 >= 0) c = th[M.iC0] * (1. + th[M.iC1] * B.yj[i]);
        else c = S.iC[j] >= 0 ? th[S.iC[j]] : 1.;
        const double b = S.iBj[j] >= 0 ? th[S.iBj[j]] : (S.iB >= 0 ? th[S.iB] : S.bFix);
        const double q = S.iQ >= 0 ? th[S.iQ] : 0.;
        const double bD = B.gd[i] + (com - 1.) * B.gdo[i], bC = B.gc[i] + (com - 1.) * B.gco[i];
        const double bg = bD + kB * bC;
        const double x = B.x[i];
        const double tl = M.tilt == 1 ? std::exp(b * x) : 1. + b * x + q * x * x;
        double m = k * B.cb[i] + c * tl * bg + B.fix[i], vs = 0;
        for (const auto& T : S.terms) {
            if (T.slice >= 0 && T.slice != j) continue;
            const double c2 = T.p2 >= 0 ? th[T.p2] : 1.;
            const double co = th[T.p1] * c2;
            m += co * B.s[T.t][i];
            vs += co * co * B.vs[T.t][i];
            if (J) { (*J)[i][T.p1] += c2 * B.s[T.t][i]; if (T.p2 >= 0) (*J)[i][T.p2] += th[T.p1] * B.s[T.t][i]; }
        }
        mu[i] = m;
        if (J) {
            if (S.iK[e] >= 0) (*J)[i][S.iK[e]] = B.cb[i];
            if (M.iC0 >= 0) { (*J)[i][M.iC0] += (1. + th[M.iC1] * B.yj[i]) * tl * bg; (*J)[i][M.iC1] += th[M.iC0] * B.yj[i] * tl * bg; }
            else if (S.iC[j] >= 0) (*J)[i][S.iC[j]] = tl * bg;
            const double dtdb = M.tilt == 1 ? x * tl : x;
            if (S.iBj[j] >= 0) (*J)[i][S.iBj[j]] = c * dtdb * bg;
            else if (S.iB >= 0) (*J)[i][S.iB] = c * dtdb * bg;
            if (S.iQ >= 0) (*J)[i][S.iQ] = c * x * x * bg;
            if (M.kBtie) (*J)[i][M.iKap] += c * tl * bC;
            else if (S.iKB >= 0) (*J)[i][S.iKB] = c * tl * bC;
            if (M.iOm >= 0) (*J)[i][M.iOm] = c * tl * (B.gdo[i] + kB * B.gco[i]);
        }
        double vb = std::max(m, 1.0) + k * k * B.vcb[i];
        if (M.bbBkg) vb += c * c * tl * tl * (B.vbd[i] + kB * kB * B.vbc[i]);
        if (V) (*V)[i] = vb + vs;
        if (Vd) (*Vd)[i] = std::max(m, 1.0) + k * k * B.vcb[i];
    }
}
inline FitRes FitX(const GBins& B, const Spec& S, const Model& M, const DataL& D, const Vec* start = nullptr, bool sandwich = true, int maxIt = 300) {
    FitRes R;
    const int n = B.n, np = S.np();
    R.n = n;
    std::vector<int> fr;
    for (int k = 0; k < np; ++k) if (!S.fixed[k]) fr.push_back(k);
    const int nf = (int)fr.size();
    R.nfree = nf;
    Vec th = start ? *start : S.start;
    for (int k = 0; k < np; ++k) if (S.fixed[k]) th[k] = S.start[k];
    Vec mu, V;
    Mat J;
    Vec Pd(nf, 0.), P0(nf, 0.);
    for (const auto& c : S.cons) for (int a = 0; a < nf; ++a) if (fr[a] == c.k) { Pd[a] += 1. / (c.sig * c.sig); P0[a] = c.val; }
    bool conv = false;
    for (int it = 0; it < maxIt; ++it) {
        EvalX(B, S, M, th, mu, &J, &V, nullptr);
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
        for (int a = 0; a < nf; ++a) th[fr[a]] += lim * d[a];
        R.nit = it + 1;
        if (it > 0 && dmax < 1e-5) { conv = true; break; }
    }
    if (!conv) return R;
    Vec Vd;
    EvalX(B, S, M, th, mu, &J, &V, &Vd);
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
    if (sandwich && !D.covN.empty()) {
        const int nI = D.nI;
        Mat U = Zero(n, nf);
        Vec sf(n, 0.), kk(n, 0.);
        for (int i = 0; i < n; ++i) {
            for (int a = 0; a < nf; ++a) U[i][a] = J[i][fr[a]] / V[i];
            sf[i] = B.y[i] > 0 ? std::sqrt(std::max(mu[i], 1e-9) / B.y[i]) : 0.;
            const int e = B.e[i];
            kk[i] = S.iK[e] >= 0 ? th[S.iK[e]] : S.kFix[e];
        }
        Mat Mn = Zero(nf, nf), Mc = Zero(nf, nf);
        Vec tn(nf), tc(nf);
        for (int i = 0; i < n; ++i) {
            std::fill(tn.begin(), tn.end(), 0.); std::fill(tc.begin(), tc.end(), 0.);
            const size_t ri = (size_t)B.item[i] * nI;
            for (int l = 0; l < n; ++l) {
                double sn = D.covN[ri + B.item[l]];
                if (B.y[i] > 0 && B.y[l] > 0) sn *= sf[i] * sf[l];
                else sn = (i == l) ? std::max(mu[i], 0.) : 0.;
                const double scb = D.covCB[ri + B.item[l]] * kk[i] * kk[l];
                if (sn == 0 && scb == 0) continue;
                for (int a = 0; a < nf; ++a) { tn[a] += sn * U[l][a]; tc[a] += scb * U[l][a]; }
            }
            for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { Mn[a][c] += U[i][a] * tn[c]; Mc[a][c] += U[i][a] * tc[c]; }
        }
        R.Cn = sand(Mn); R.Ccb = sand(Mc);
        R.Cstat = wp5::Add(R.Cn, R.Ccb);
        R.Ccons.clear();
        for (const auto& c : S.cons) {
            Mat Mp = Zero(nf, nf);
            for (int a = 0; a < nf; ++a) if (fr[a] == c.k) Mp[a][a] = 1. / (c.sig * c.sig);
            R.Ccons.push_back(sand(Mp));
            R.Cstat = wp5::Add(R.Cstat, R.Ccons.back());
        }
    } else {
        R.Cn = R.C; R.Ccb = Zero(np, np); R.Cstat = R.C;
        for (size_t q = 0; q < S.cons.size(); ++q) R.Ccons.push_back(Zero(np, np));
    }
    R.th = th; R.mu = mu; R.V = V; R.Vd = Vd;
    R.chi2BC = R.chi2P = R.chi2E = 0;
    for (int i = 0; i < n; ++i) {
        const double m = std::max(mu[i], 1e-9), y = B.y[i];
        R.chi2BC += 2. * (m - y + (y > 0 ? y * std::log(y / m) : 0.));
        R.chi2P += (y - mu[i]) * (y - mu[i]) / V[i];
    }
    for (const auto& c : S.cons) { const double z = (th[c.k] - c.val) / c.sig; R.chi2BC += z * z; R.chi2P += z * z; }
    R.ndf = n - nf + (int)S.cons.size();
    R.ok = true;
    return R;
}

// ---------------------------------------------------------------------------------------------------------------
// 7. specs with the model options
// ---------------------------------------------------------------------------------------------------------------
struct BkgOpt {
    bool quad = false, tiltj = false, expb = false, omega = false, commonC = false, kBtie = false;
    bool fixKB = false; double kBval = 1.;
    double kE2fix = -1;     // >= 0: k_E2 fixed to this value; < 0: constrained 1 +- 0.05
    bool useE1 = false;     // k_E1 free (spec configuration)
};
inline void AddBkgX(Spec& S, Model& M, const BkgOpt& o) {
    if (o.commonC) { M.iC0 = S.Add("c0", 1.0); M.iC1 = S.Add("c1", 0.0); }
    else for (int j = 0; j < kNSlice; ++j) S.iC[j] = S.Add(Form("c_%d", j + 1), 1.0);
    if (o.tiltj) { for (int j = 0; j < kNSlice; ++j) S.iBj[j] = S.Add(Form("b_%d", j + 1), -0.07); }
    else S.iB = S.Add("b", -0.07);
    if (o.quad) S.iQ = S.Add("q", 0.0);
    M.tilt = o.expb ? 1 : 0;
    S.iKB = S.Add("kappa_B", o.fixKB ? o.kBval : 2.0);
    if (o.fixKB || o.kBtie) S.fixed[S.iKB] = 1;
    if (o.omega) M.iOm = S.Add("c_om", 1.0);
    if (o.useE1) S.iK[0] = S.Add("k_E1", 1.0);
    S.iK[1] = S.Add("k_E2", o.kE2fix >= 0 ? o.kE2fix : 1.0);
    if (o.kE2fix >= 0) S.fixed[S.iK[1]] = 1; else S.cons.push_back({S.iK[1], 1.0, 0.05});
}
struct EtaOpt { double kappaSig = 0.2; bool kappaFree = false; double gSig = 0.30; bool gFree = false; };
inline Spec SpecMainX(const EtaOpt& eo, const BkgOpt& bo, Model& M) {
    Spec S; M = Model();
    const int ia = S.Add("a_eta", 0.4), ik = S.Add("kappa", 1.0), ir = S.Add("rho_d", 1.0), ig = S.Add("g", 1.0);
    if (!eo.kappaFree) S.cons.push_back({ik, 1.0, eo.kappaSig});
    if (!eo.gFree) S.cons.push_back({ig, 1.0, eo.gSig});
    S.terms.push_back({wp5::kS2, ia, -1, -1});
    S.terms.push_back({wp5::kS3, ia, ik, -1});
    S.terms.push_back({wp5::kS4, ia, ik, -1});
    S.terms.push_back({wp5::kS5, ia, ir, -1});
    S.terms.push_back({wp5::kSeta, ia, ig, -1});
    AddBkgX(S, M, bo);
    if (bo.kBtie) { M.kBtie = true; M.iKap = ik; }
    return S;
}
inline Spec SpecGoldX(const EtaOpt& eo, const BkgOpt& bo, double kB, Model& M) {
    Spec S; M = Model();
    const int i2 = S.Add("alpha_2", 0.4), i5 = S.Add("alpha_5", 0.4), ig = S.Add("g", 1.0);
    if (!eo.gFree) S.cons.push_back({ig, 1.0, eo.gSig});
    S.terms = {{wp5::kS2, i2, -1, -1}, {wp5::kS5, i5, -1, -1}, {wp5::kSeta, i2, ig, -1}};
    BkgOpt b = bo; b.fixKB = true; b.kBval = kB; b.kBtie = false; b.omega = false; b.useE1 = false;
    AddBkgX(S, M, b);
    return S;
}
inline Spec SpecStageAX(const BkgOpt& bo, Model& M) {
    Spec S; M = Model();
    BkgOpt b = bo; b.kBtie = false; b.useE1 = false;
    AddBkgX(S, M, b);
    return S;
}

// ---------------------------------------------------------------------------------------------------------------
// 8. stage-A background in windows and the estimators (wp6_ana.C, generalised)
// ---------------------------------------------------------------------------------------------------------------
// B_A of cell h (main histogram) over [lo, hi): c_j int tilt(m) [G_D + (c_om - 1) G_D,om + kB (G_C + (c_om - 1) G_C,om)]
inline double BAint(const SmSet& S, const Spec& Sp, const Model& M, const Vec& th, const DataL& D, int h, double lo, double hi) {
    const int j = h / 2;
    const double kB = M.kBtie ? th[M.iKap] : (Sp.iKB >= 0 ? th[Sp.iKB] : Sp.kBfix);
    const double b = Sp.iBj[j] >= 0 ? th[Sp.iBj[j]] : (Sp.iB >= 0 ? th[Sp.iB] : Sp.bFix);
    const double q = Sp.iQ >= 0 ? th[Sp.iQ] : 0.;
    const double com = M.iOm >= 0 ? th[M.iOm] : 1.;
    const double c = M.iC0 >= 0 ? th[M.iC0] * (1. + th[M.iC1] * (D.mpp[j] - 0.38) / 0.05) : th[Sp.iC[j]];
    double v = IntTiltSD(S.d[h][0], lo, hi, b, q, M.tilt) + kB * IntTiltSD(S.d[h][1], lo, hi, b, q, M.tilt);
    if (M.iOm >= 0) v += (com - 1.) * (IntTiltSD(S.d[h][2], lo, hi, b, q, M.tilt) + kB * IntTiltSD(S.d[h][3], lo, hi, b, q, M.tilt));
    return c * v;
}
inline Vec BATable(const SmSet& S, const Spec& Sp, const Model& M, const Vec& th, const DataL& D) {
    Vec t((size_t)10 * kNW, 0.);
    for (int h = 0; h < 10; ++h)
        for (int w = 0; w < kNW; ++w) {
            if (kWins[w].lo < -1e8) continue;
            t[h * kNW + w] = BAint(S, Sp, M, th, D, h, kWins[w].lo, kWins[w].hi);
        }
    return t;
}
// estimator input state
struct St { const Vec* N; const Vec* CB; const Vec* BA; const Tmpl* T; Vec f; const DataL* D; };
inline std::vector<int> Cells(int e, int slice = -1) {
    std::vector<int> c;
    for (int j = 0; j < kNSlice; ++j) {
        if (slice >= 0 && j != slice) continue;
        if (e == 0 || e == 2) c.push_back(2 * j);
        if (e == 1 || e == 2) c.push_back(2 * j + 1);
    }
    return c;
}
inline double SN(const St& s, const std::vector<int>& c, int w) { double v = 0; for (int h : c) v += (*s.N)[s.D->IWin(h, w)]; return v; }
inline double SC(const St& s, const std::vector<int>& c, int w) { double v = 0; for (int h : c) v += (*s.CB)[s.D->IWin(h, w)]; return v; }
inline double SB(const St& s, const std::vector<int>& c, int w) { double v = 0; for (int h : c) v += (*s.BA)[h * kNW + w]; return v; }
inline double SS(const St& s, int t, const std::vector<int>& c, int w) {
    if (c.size() == 10) return s.T->wa[t][w];
    double v = 0; for (int h : c) v += s.T->ws[t][h * kNW + w];
    return v;
}
inline double EtaAll(const St& s, const std::vector<int>& c, int w) {
    const Vec& f = s.f;
    return f[0] * SS(s, 0, c, w) + f[0] * f[1] * (SS(s, 1, c, w) + SS(s, 2, c, w)) + f[0] * f[2] * SS(s, 3, c, w) + f[0] * f[3] * SS(s, 4, c, w);
}
inline double EtaOth(const St& s, const std::vector<int>& c, int w) {
    const Vec& f = s.f;
    return f[0] * f[1] * SS(s, 2, c, w) + f[0] * f[2] * SS(s, 3, c, w) + f[0] * f[3] * SS(s, 4, c, w);
}
inline double Tail(const St& s, const std::vector<int>& c, int w) {
    const Vec& f = s.f;
    return f[0] * SS(s, 0, c, w) + f[0] * f[1] * SS(s, 1, c, w) + f[0] * f[3] * SS(s, 4, c, w);
}
inline double U45(const St& s, const std::vector<int>& c, int w) { return SS(s, 2, c, w) + SS(s, 3, c, w); }
typedef std::function<double(const St&)> Est;
inline Est MakeB4(const std::vector<int>& cells, int wU) {
    return [=](const St& s) {
        const double sB = SN(s, cells, kWB) - SC(s, cells, kWB);
        const double sU = SN(s, cells, wU) - SC(s, cells, wU);
        const double R = SB(s, cells, kWB) / SB(s, cells, wU);
        return kYf * (sB - (sU - EtaAll(s, cells, wU)) * R - EtaOth(s, cells, kWB));
    };
}
inline double FracE2(const St& s) {
    const std::vector<int> e2 = Cells(1), all = Cells(2);
    return (SS(s, 0, e2, kWB) + s.f[1] * SS(s, 1, e2, kWB)) / (SS(s, 0, all, kWB) + s.f[1] * SS(s, 1, all, kWB));
}
inline Est MakeB4All(int wU) { return [=](const St& s) { return MakeB4(Cells(1), wU)(s) / FracE2(s); }; }
inline Est MakePSB() {
    return [=](const St& s) {
        const std::vector<int> c = Cells(2);
        const double sP = SN(s, c, kWP) - SC(s, c, kWP);
        const double sS = SN(s, c, kWSl) - SC(s, c, kWSl) + SN(s, c, kWSh) - SC(s, c, kWSh);
        const double bP = SB(s, c, kWP), bS = SB(s, c, kWSl) + SB(s, c, kWSh);
        const double tP = Tail(s, c, kWP), tS = Tail(s, c, kWSl) + Tail(s, c, kWSh);
        const double f = (U45(s, c, kWSl) + U45(s, c, kWSh)) / U45(s, c, kWP);
        const double R = (bP + tP) / (bS + tS);
        return kYf * (sP - R * sS) / (1. - R * f);
    };
}
inline Est MakeAP() { return [=](const St& s) { return MakePSB()(s) / (kYf * U45(s, Cells(2), kWP)); }; }
inline Est MakeAB(int wU) { return [=](const St& s) { return MakeB4(Cells(1), wU)(s) / (kYf * (SS(s, 0, Cells(1), kWB) + SS(s, 1, Cells(1), kWB))); }; }
inline Est MakeDR(int wU) {
    return [=](const St& s) {
        const double yb = MakeB4All(wU)(s), yp = MakePSB()(s);
        const std::vector<int> a = Cells(2);
        const double sb = SS(s, 0, a, kWB) + SS(s, 1, a, kWB), sp = U45(s, a, kWP);
        return (yp / yb) / (sp / sb);
    };
}

}  // namespace wp7
#endif
