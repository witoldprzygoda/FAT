// ============================================================================
// wp0_fill.C  --  WP0 histogram factory of SPEC report_all (prefix wp0)
// ============================================================================
// One pass over the data compact rows (scan_v2: us, pp, mm) and one pass over the
// simulated e+e- rows (std_us, lep_us, eta_us, owned), filling every histogram of
// WP0 into /hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root (directory and
// naming scheme: README.md of this directory and the TNamed "README" in the file).
//
//   data : N, N++, N--, CB (= 2 sqrt(N++ N--) in 10 MeV cells, interpolated, spread
//          over the cells by f(cell | 50 MeV band) from N++ + N--), CB20 (20 MeV
//          cells), 200 Poisson event-group replicas of the CB; subsets all / folds /
//          periods; selections S0, pre, dz40, PEAK, vz150, vz170, none, S0_best,
//          pre_best, S0_corfb, pre_corfb
//   sim  : components C2-C6, G_eta, G_pi0D, G_conv, G_vec, G_oth (+ omega/rho0 split),
//          G_D, T12, T3, T4 in nb per bin (sum w / L_sim) with event-group errors;
//          halves A/B (even/odd event group) and 100 Poisson event-group replicas
//          for S0__nom and pre__nom; variants of weights, ownership, ghosts, mass
//          transfer and selection
//   blinding (D19): the data N of the 2 MeV M(ee) bins overlapping 13-27 MeV in X2D
//          (all M4) and in B_mee_R2..R4 is replaced by -1 in hists_v1.root; the
//          unmasked cells are in sealed/hists_v1_x17_sealed.root (open only after
//          x17_frozen.txt exists)
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp0", kTRUE)' \
//        'plots/all_report/wp0_fill.C+O(-1)'
// Argument: maximum number of rows per tree (-1 = all; a test run writes
// hists_v1_test.root instead).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TNamed.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "defs.h"
#include "sim_eta_physics.h"

namespace wp0 {
using namespace alldefs;

// ---------------------------------------------------------------------------
// counter-based random numbers: Poisson(1) weight of an event group in a replica
// ---------------------------------------------------------------------------
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
const uint64_t kSeedSim = 0x5157A11ULL, kSeedCB = 0xCB2026ULL;

// ---------------------------------------------------------------------------
// binning helpers
// ---------------------------------------------------------------------------
inline std::vector<double> Uniform(int n, double lo, double hi) {
    std::vector<double> e(n + 1);
    for (int i = 0; i <= n; ++i) e[i] = lo + (hi - lo) * i / n;
    return e;
}
// 0 underflow, 1..n, n + 1 overflow (NaN -> overflow)
inline int FindBin(const std::vector<double>& e, double x) {
    const int n = (int)e.size() - 1;
    if (!(x == x)) return n + 1;
    if (x < e[0]) return 0;
    if (x >= e[n]) return n + 1;
    return (int)(std::upper_bound(e.begin(), e.end(), x) - e.begin());
}
inline int UBin(int n, double lo, double hi, double x) {
    if (!(x == x)) return n + 1;
    if (x < lo) return 0;
    if (x >= hi) return n + 1;
    return std::min(n, (int)((x - lo) / (hi - lo) * n) + 1);
}

struct HShape {
    std::string name;
    int dim = 1;
    std::vector<double> ex, ey, ez;
    int nx() const { return (int)ex.size() - 1; }
    int ny() const { return dim >= 2 ? (int)ey.size() - 1 : 0; }
    int nz() const { return dim >= 3 ? (int)ez.size() - 1 : 0; }
    int ncell() const { return (nx() + 2) * (dim >= 2 ? ny() + 2 : 1) * (dim >= 3 ? nz() + 2 : 1); }
};
inline HShape Shape1(const std::string& n, const std::vector<double>& ex) { HShape s; s.name = n; s.dim = 1; s.ex = ex; return s; }
inline HShape Shape2(const std::string& n, const std::vector<double>& ex, const std::vector<double>& ey) {
    HShape s; s.name = n; s.dim = 2; s.ex = ex; s.ey = ey; return s;
}
inline HShape Shape3(const std::string& n, const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez) {
    HShape s; s.name = n; s.dim = 3; s.ex = ex; s.ey = ey; s.ez = ez; return s;
}

// global edge tables
std::vector<double> gF10e, gX2Dex, gX2Dey, gG3ex, gG3ey, gG3ez, gFinee, gBM4e, gMapx, gMapy, gWine, gTote;
std::vector<double> gBVe[kNBVar];
std::vector<int> gF10x, gX2Dx, gG3x;   // fine bin -> axis bin
void InitEdges() {
    gF10e = F10Edges();
    gX2Dex = Uniform(kX2Dnx, kX2Dxlo, kX2Dxhi);
    gX2Dey = Uniform(kX2Dny, kX2Dylo, kX2Dyhi);
    gG3ex = Uniform(kG3nx, kG3xlo, kG3xhi);
    gG3ey = G3PiPiEdges();
    gG3ez = G3MeeEdges();
    gFinee = Uniform(kNfine, kFineLo, kFineHi);
    gBM4e = Uniform(kNbefM4, kFineLo, kFineHi);
    gMapx = Uniform(kMapNx, kMapXlo, kMapXhi);
    gMapy = Uniform(kMapNy, kMapYlo, kMapYhi);
    gWine = Uniform(kNWin, 0., (double)kNWin);
    gTote = Uniform(1, 0., 1.);
    for (int v = 0; v < kNBVar; ++v) gBVe[v] = Uniform(kBVnb[v], kBVlo[v], kBVhi[v]);
    gF10x.resize(kNfine); gX2Dx.resize(kNfine); gG3x.resize(kNfine);
    for (int i = 0; i < kNfine; ++i) {
        const double c = kFineLo + (i + 0.5) * kFineW;
        gF10x[i] = FindBin(gF10e, c);
        gX2Dx[i] = FindBin(gX2Dex, c);
        gG3x[i] = FindBin(gG3ex, c);
    }
}

// ---------------------------------------------------------------------------
// layouts (sets of histograms of one block)
// ---------------------------------------------------------------------------
enum LType { kLstd = 0, kLbef, kLfull, kLrepSim, kLrepData, kNLType };
struct Layout {
    std::vector<HShape> h;
    std::vector<int> off;
    int stride = 0;
    int iF10[10], iF5g[5], iX2D = -1, iFine = -1, iBM4 = -1, iWin = -1, iTot = -1, iG3 = -1, iMap = -1;
    int iBV[kNBVar][kNReg];
    bool hasBef = false;
    int Add(const HShape& s) { h.push_back(s); off.push_back(stride); stride += s.ncell(); return (int)h.size() - 1; }
};
Layout MakeLayout(int type) {
    Layout L;
    for (int k = 0; k < 10; ++k) L.iF10[k] = -1;
    for (int k = 0; k < 5; ++k) L.iF5g[k] = -1;
    for (int v = 0; v < kNBVar; ++v) for (int k = 0; k < kNReg; ++k) L.iBV[v][k] = -1;
    if (type == kLrepData) {
        L.iFine = L.Add(Shape1("fine_M4", gFinee));
        L.iBM4 = L.Add(Shape1("B_M4", gBM4e));
    }
    for (int s = 0; s < kNSlice; ++s)
        for (int e = 0; e < 2; ++e) L.iF10[2 * s + e] = L.Add(Shape1(F10Name(s, e), gF10e));
    for (int s = 0; s < kNSlice; ++s) L.iF5g[s] = L.Add(Shape1(F5gName(s), gF10e));
    L.iX2D = L.Add(Shape2("X2D", gX2Dex, gX2Dey));
    if (type == kLstd || type == kLbef || type == kLfull) {
        L.iFine = L.Add(Shape1("fine_M4", gFinee));
        L.iBM4 = L.Add(Shape1("B_M4", gBM4e));
    }
    L.iWin = L.Add(Shape1("win", gWine));
    if (type != kLrepData) L.iTot = L.Add(Shape1("tot", gTote));
    if (type == kLbef || type == kLfull) {
        L.hasBef = true;
        L.iMap = L.Add(Shape2("B_map", gMapx, gMapy));
        for (int v = 0; v < kNBVar; ++v)
            for (int k = 0; k < kNReg; ++k) L.iBV[v][k] = L.Add(Shape1(BVarName(v, k), gBVe[v]));
    }
    if (type == kLfull) L.iG3 = L.Add(Shape3("G3", gG3ex, gG3ey, gG3ez));
    return L;
}
Layout gL[kNLType];

// ---------------------------------------------------------------------------
// cells of one row (for one value of M4)
// ---------------------------------------------------------------------------
struct Cells {
    int fi, slice, ecl, f10h, f10c, f5gh, x2d, x2y, fine, bm4, g3, g3y, g3z, map, my, reg, nwin;
    int bv[kNBVar];
    int win[kNWin];
};
inline void MakeCells(double m4, double mpipi, double mee, const double* bvals, Cells& c) {
    const int fi = FineBin(m4);
    c.fi = fi;
    const bool under = (m4 < kFineLo);   // NaN -> overflow
    c.slice = Slice(mpipi);
    c.ecl = Eclass(mee);
    c.f10c = fi >= 0 ? gF10x[fi] : (under ? 0 : kNF10 + 1);
    c.f10h = (c.slice >= 0 && c.ecl >= 0) ? 2 * c.slice + c.ecl : -1;
    c.f5gh = (c.slice >= 0 && mee >= kF5gMee) ? c.slice : -1;
    const int x2x = fi >= 0 ? gX2Dx[fi] : (under ? 0 : kX2Dnx + 1);
    c.x2y = FindBin(gX2Dey, mee);
    c.x2d = x2x + (kX2Dnx + 2) * c.x2y;
    c.fine = fi >= 0 ? fi + 1 : (under ? 0 : kNfine + 1);
    c.bm4 = fi >= 0 ? fi / 2 + 1 : (under ? 0 : kNbefM4 + 1);
    const int g3x = fi >= 0 ? gG3x[fi] : (under ? 0 : kG3nx + 1);
    c.g3y = FindBin(gG3ey, mpipi);
    c.g3z = FindBin(gG3ez, mee);
    c.g3 = g3x + (kG3nx + 2) * (c.g3y + ((int)gG3ey.size() + 1) * c.g3z);
    const int mx = fi >= 0 ? fi / 4 + 1 : (under ? 0 : kMapNx + 1);
    c.my = UBin(kMapNy, kMapYlo, kMapYhi, mpipi);
    c.map = mx + (kMapNx + 2) * c.my;
    c.reg = Region(m4);
    for (int v = 0; v < kNBVar; ++v) c.bv[v] = UBin(kBVnb[v], kBVlo[v], kBVhi[v], bvals[v]);
    c.nwin = 0;
    for (int w = 0; w < kNWin; ++w) if (InWindow(w, m4, mee, mpipi)) c.win[c.nwin++] = w;
}
template <class F> inline void ForCells(const Layout& L, const Cells& c, F f) {
    if (c.f10h >= 0 && L.iF10[0] >= 0) f(L.off[L.iF10[c.f10h]] + c.f10c);
    if (c.f5gh >= 0 && L.iF5g[0] >= 0) f(L.off[L.iF5g[c.f5gh]] + c.f10c);
    if (L.iX2D >= 0) f(L.off[L.iX2D] + c.x2d);
    if (L.iFine >= 0) f(L.off[L.iFine] + c.fine);
    if (L.iBM4 >= 0) f(L.off[L.iBM4] + c.bm4);
    if (L.iWin >= 0) for (int k = 0; k < c.nwin; ++k) f(L.off[L.iWin] + c.win[k] + 1);
    if (L.iTot >= 0) f(L.off[L.iTot] + 1);
    if (L.iG3 >= 0) f(L.off[L.iG3] + c.g3);
    if (L.iMap >= 0) f(L.off[L.iMap] + c.map);
    if (L.hasBef && c.reg >= 0)
        for (int v = 0; v < kNBVar; ++v) f(L.off[L.iBV[v][c.reg]] + c.bv[v]);
}

// ---------------------------------------------------------------------------
// accumulators
// ---------------------------------------------------------------------------
// sums with event-group variance (and optional half-A sums): v[i] = sum over groups of (sum w in cell i)^2
struct Acc {
    std::vector<double> s, v, sA, vA, acc;
    std::vector<int> stamp, touched;
    int gid = 0;
    bool a = false, halves = false;
    void Init(size_t n, bool h) {
        halves = h;
        s.assign(n, 0.); v.assign(n, 0.); acc.assign(n, 0.); stamp.assign(n, 0);
        if (h) { sA.assign(n, 0.); vA.assign(n, 0.); }
        touched.reserve(1 << 16);
    }
    inline void Begin(bool halfA) { Flush(); ++gid; a = halfA; }
    inline void Add(int i, double w) {
        s[i] += w;
        if (halves && a) sA[i] += w;
        if (stamp[i] != gid) { stamp[i] = gid; acc[i] = 0.; touched.push_back(i); }
        acc[i] += w;
    }
    void Flush() {
        for (int i : touched) {
            const double q = acc[i] * acc[i];
            v[i] += q;
            if (halves && a) vA[i] += q;
        }
        touched.clear();
    }
};
struct RepAcc {
    int nrep = 0;
    std::vector<double> r;
    void Init(size_t ncell, int n) { nrep = n; r.assign(ncell * n, 0.); }
    inline void Add(size_t i, double w, const double* k) {
        double* p = &r[i * nrep];
        for (int q = 0; q < nrep; ++q) p[q] += w * k[q];
    }
};

// ---------------------------------------------------------------------------
// histogram writing
// ---------------------------------------------------------------------------
// content c[cell] x scale, error sqrt(v[cell]) x scale (v may be null: error 0)
TH1* MakeHist(const HShape& s, const std::string& name, const double* c, const double* v, double scale, const char* title = "") {
    TH1* h = nullptr;
    if (s.dim == 1) h = new TH1D(name.c_str(), title, s.nx(), s.ex.data());
    else if (s.dim == 2) h = new TH2D(name.c_str(), title, s.nx(), s.ex.data(), s.ny(), s.ey.data());
    else h = new TH3D(name.c_str(), title, s.nx(), s.ex.data(), s.ny(), s.ey.data(), s.nz(), s.ez.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    const int n = s.ncell();
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        h->SetBinContent(i, c[i] * scale);
        h->SetBinError(i, v ? std::sqrt(std::max(0., v[i])) * scale : 0.);
        sum += c[i];
    }
    h->SetEntries(sum);
    if (s.name == "win")
        for (int w = 0; w < kNWin; ++w) h->GetXaxis()->SetBinLabel(w + 1, kWin[w].tag);
    return h;
}
// replicas: one more axis (replica index 0..nrep-1); rep[cell * nrep + r]
TH1* MakeRepHist(const HShape& s, const std::string& name, const double* rep, int nrep, double scale) {
    std::vector<double> re = Uniform(nrep, -0.5, nrep - 0.5);
    TH1* h = nullptr;
    if (s.dim == 1) h = new TH2D(name.c_str(), "replicas (y = replica index)", s.nx(), s.ex.data(), nrep, re.data());
    else h = new TH3D(name.c_str(), "replicas (z = replica index)", s.nx(), s.ex.data(), s.ny(), s.ey.data(), nrep, re.data());
    h->SetDirectory(nullptr);
    const int n = s.ncell(), nxy = n;
    for (int i = 0; i < n; ++i)
        for (int r = 0; r < nrep; ++r) h->SetBinContent(i + nxy * (r + 1), rep[(size_t)i * nrep + r] * scale);
    if (s.name == "win")
        for (int w = 0; w < kNWin; ++w) h->GetXaxis()->SetBinLabel(w + 1, kWin[w].tag);
    return h;
}
void WriteObj(TDirectory* d, TObject* o) { d->WriteTObject(o, o->GetName(), "Overwrite"); delete o; }

// ===========================================================================
// CB construction
// ===========================================================================
// 1D CB on the fine grid from the like-sign fine counts: wide cells of reb fine bins,
// 2 sqrt(N++ N--) per wide cell, linear interpolation to the fine-bin centres
// (ScanEngine.h Spectra::Interp). Keeps the interpolation weights for the exact
// linear error propagation (variance of a wide cell = N++ + N--).
struct CB1D {
    int reb = 4, nw = 80;
    std::vector<double> cbW, varW, cb;
    std::vector<int> j0, j1;
    std::vector<double> a0, a1;
    void Make(const double* pp, const double* mm, int r) {   // pp, mm: fine counts [kNfine]
        reb = r; nw = kNfine / reb;
        std::vector<double> wp(nw, 0.), wm(nw, 0.);
        for (int i = 0; i < kNfine; ++i) { wp[i / reb] += pp[i]; wm[i / reb] += mm[i]; }
        cbW.assign(nw, 0.); varW.assign(nw, 0.);
        for (int j = 0; j < nw; ++j) { cbW[j] = 2. * std::sqrt(wp[j] * wm[j]); varW[j] = wp[j] + wm[j]; }
        cb.assign(kNfine, 0.); j0.assign(kNfine, 0); j1.assign(kNfine, 0); a0.assign(kNfine, 0.); a1.assign(kNfine, 0.);
        for (int i = 0; i < kNfine; ++i) {
            const double t = (i + 0.5) / reb - 0.5;
            int k0, k1; double f;
            if (t <= 0) { k0 = k1 = 0; f = 0; }
            else if (t >= nw - 1) { k0 = k1 = nw - 1; f = 0; }
            else { k0 = (int)t; k1 = k0 + 1; f = t - k0; }
            j0[i] = k0; j1[i] = k1; a0[i] = (1 - f) / reb; a1[i] = f / reb;
            cb[i] = a0[i] * cbW[k0] + a1[i] * cbW[k1];
        }
    }
    // sum and linear variance of sum_i u_i cb_i
    std::pair<double, double> Sum(const std::vector<std::pair<int, double>>& u) const {
        std::vector<double> c(nw, 0.);
        double s = 0;
        for (const auto& p : u) {
            const int i = p.first;
            s += p.second * cb[i];
            c[j0[i]] += p.second * a0[i];
            c[j1[i]] += p.second * a1[i];
        }
        double v = 0;
        for (int j = 0; j < nw; ++j) v += c[j] * c[j] * varW[j];
        return std::make_pair(s, v);
    }
};
// cell families of the CB spreading (cells independent of M4)
enum Fam { kFamF10 = 0, kFamF5g, kFamX2D, kFamG3, kFamMap, kFamX17, kFamCRK, kNFam };
const int kFamN[kNFam] = {10, 6, kX2Dny + 2, 29 * 7, kMapNy + 2, 2, 2};
const bool kFamRep[kNFam] = {true, true, true, false, false, true, true};
inline int FamCell(int fam, const Cells& c, double mee) {
    switch (fam) {
        case kFamF10: return c.f10h;
        case kFamF5g: return c.slice < 0 ? -1 : (mee >= kF5gMee ? c.slice : 5);
        case kFamX2D: return c.x2y;
        case kFamG3: return c.g3y + 29 * c.g3z;
        case kFamMap: return c.my;
        case kFamX17: return (mee >= kWin[kWX17c].meeLo && mee < kWin[kWX17c].meeHi) ? 0 : 1;
        case kFamCRK: return c.slice < 0 ? -1 : (c.slice == kWin[kWCRK].slice ? 0 : 1);
    }
    return -1;
}
// f(cell | band) from the counts cnt[band * n + cell]; empty bands take the pooled fraction (counted in nEmpty)
std::vector<double> Fractions(const double* cnt, int n, int* nEmpty = nullptr) {
    std::vector<double> f(kNBand * n, 0.), pool(n, 0.);
    double ptot = 0;
    for (int b = 0; b < kNBand; ++b) for (int c = 0; c < n; ++c) { pool[c] += cnt[b * n + c]; ptot += cnt[b * n + c]; }
    for (int b = 0; b < kNBand; ++b) {
        double t = 0;
        for (int c = 0; c < n; ++c) t += cnt[b * n + c];
        for (int c = 0; c < n; ++c) f[b * n + c] = t > 0 ? cnt[b * n + c] / t : (ptot > 0 ? pool[c] / ptot : (c == 0 ? 1. : 0.));
        if (!(t > 0) && nEmpty) ++*nEmpty;
    }
    return f;
}
// fraction of fine bin i inside [a, b) (edges within 1e-6 of a fine-bin edge are snapped to it)
inline double FracFine(int i, double a, double b) {
    double fa = (a - kFineLo) / kFineW, fb = (b - kFineLo) / kFineW;
    if (std::fabs(fa - std::round(fa)) < 1e-6) fa = std::round(fa);
    if (std::fabs(fb - std::round(fb)) < 1e-6) fb = std::round(fb);
    const double o = std::min((double)(i + 1), fb) - std::max((double)i, fa);
    return o > 0 ? std::min(1., o) : 0.;
}

// ===========================================================================
// data
// ===========================================================================
struct DSel { const char* tag; int sel; bool best; bool fb; const char* text; };
const int kNDSel = 11;
const DSel kDSel[kNDSel] = {
    {"S0", kSelS0, false, false, "S0 = pre + both pions in the vertex fit (NOMINAL)"},
    {"pre", kSelPre, false, false, "pre"},
    {"dz40", kSelDz40, false, false, "pre + dz_pi < 40 mm (DIRECT)"},
    {"PEAK", kSelPeak, false, false, "pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2"},
    {"vz150", kSelVz150, false, false, "S0 with vz > -150 mm"},
    {"vz170", kSelVz170, false, false, "S0 with vz > -170 mm"},
    {"none", kSelNone, false, false, "none = isBest >= 0, vz > -500, PT3, cor_ok (no vz > -160, no pion-vertex cut)"},
    {"S0_best", kSelS0, true, false, "S0 and isBest == 1"},
    {"pre_best", kSelPre, true, false, "pre and isBest == 1"},
    {"S0_corfb", kSelS0, false, true, "S0, cor_ok fallback (cor_ok == 0 rows with the reconstructed pi-)"},
    {"pre_corfb", kSelPre, false, true, "pre, cor_ok fallback"}};
const char* const kWhat[3] = {"N", "Npp", "Nmm"};

struct DataOut {
    // blocks: (dsel, sub, sample) -> base offset into the Acc arrays
    int base[kNDSel][kNSub][3];
    Acc A;
    // like-sign band-cell counts: [dsel][sub][fam] -> vector kNBand * kFamN
    std::vector<double> lsbc[kNDSel][kNSub][kNFam];
    // replicas (sub all): fine like-sign counts [dsel][pp/mm] (kNfine * kNRepCB), band-cell [dsel][fam]
    std::vector<double> lsFineRep[kNDSel][2];
    std::vector<double> lsbcRep[kNDSel][kNFam];
    long long nRows[3] = {0, 0, 0}, nPass[kNDSel][3];
    long long nNaNcell[kNDSel] = {0};
};

struct DVals { double m4, mee, mpipi, oa, mm4, pt4, y4; };
inline DVals ValsOf(const scan::Row& r, bool fb) {
    DVals d;
    if (fb && !r.cor_ok) { d.m4 = r.m4_fb; d.mee = r.mee_fb; d.mpipi = r.mpipi_fb; d.oa = r.oa_fb; d.mm4 = r.mm4_fb; d.pt4 = r.pt4_fb; d.y4 = r.y4_fb; }
    else { d.m4 = r.m4; d.mee = r.mee; d.mpipi = r.mpipi; d.oa = r.oa; d.mm4 = r.mm4; d.pt4 = r.pt4; d.y4 = r.y4; }
    return d;
}

void FillData(const scan::Data& D, DataOut& O, long long maxRows) {
    TStopwatch sw;
    size_t tot = 0;
    for (int ds = 0; ds < kNDSel; ++ds)
        for (int sb = 0; sb < kNSub; ++sb)
            for (int k = 0; k < 3; ++k) {
                O.base[ds][sb][k] = (int)tot;
                tot += gL[sb == kSubAll ? kLfull : kLstd].stride;
            }
    printf("wp0: data cells %zu (%.1f MB per array)\n", tot, tot * 8e-6);
    O.A.Init(tot, false);
    for (int ds = 0; ds < kNDSel; ++ds) {
        for (int sb = 0; sb < kNSub; ++sb)
            for (int f = 0; f < kNFam; ++f) O.lsbc[ds][sb][f].assign(kNBand * kFamN[f], 0.);
        for (int q = 0; q < 2; ++q) O.lsFineRep[ds][q].assign((size_t)kNfine * kNRepCB, 0.);
        for (int f = 0; f < kNFam; ++f) if (kFamRep[f]) O.lsbcRep[ds][f].assign((size_t)kNBand * kFamN[f] * kNRepCB, 0.);
        for (int k = 0; k < 3; ++k) O.nPass[ds][k] = 0;
    }
    std::vector<double> kr(kNRepCB, 0.);
    for (int k = 0; k < 3; ++k) {
        const std::vector<scan::Row>& v = D.s[k];
        const size_t n = maxRows > 0 ? std::min((size_t)maxRows, v.size()) : v.size();
        int cur = -1;
        bool haveK = false;
        for (size_t i = 0; i < n; ++i) {
            const scan::Row& r = v[i];
            if (r.grp != cur) { O.A.Begin(false); cur = r.grp; haveK = false; }
            O.nRows[k]++;
            bool pass[kNDSel];
            bool any = false, anyFb = false, anyNo = false;
            const DVals v0 = ValsOf(r, false), v1 = ValsOf(r, true);
            for (int ds = 0; ds < kNDSel; ++ds) {
                const DSel& d = kDSel[ds];
                const DVals& vv = d.fb ? v1 : v0;
                const bool corOk = d.fb ? true : (r.cor_ok == 1);
                pass[ds] = PassSel(d.sel, r, corOk, r.dz_pi, (float)vv.mm4) && (!d.best || r.isBest == 1) && (vv.m4 == vv.m4);
                if (pass[ds]) { any = true; if (d.fb) anyFb = true; else anyNo = true; }
            }
            if (!any) continue;
            Cells c0, c1;
            if (anyNo) { const double b0[kNBVar] = {v0.mpipi, v0.mee, v0.oa, r.vz, v0.mm4, v0.pt4, v0.y4}; MakeCells(v0.m4, v0.mpipi, v0.mee, b0, c0); }
            if (anyFb) { const double b1[kNBVar] = {v1.mpipi, v1.mee, v1.oa, r.vz, v1.mm4, v1.pt4, v1.y4}; MakeCells(v1.m4, v1.mpipi, v1.mee, b1, c1); }
            const int per = PeriodOfDay(r.file_id);
            const int subs[3] = {kSubAll, kSubFold0 + (r.fold & 1), per >= 0 ? kSubPerI + per : -1};
            for (int ds = 0; ds < kNDSel; ++ds) {
                if (!pass[ds]) continue;
                O.nPass[ds][k]++;
                const Cells& c = kDSel[ds].fb ? c1 : c0;
                const double mee = kDSel[ds].fb ? v1.mee : v0.mee;
                if (c.f10h < 0) O.nNaNcell[ds]++;
                for (int q = 0; q < 3; ++q) {
                    const int sb = subs[q];
                    if (sb < 0) continue;
                    const Layout& L = gL[sb == kSubAll ? kLfull : kLstd];
                    const int base = O.base[ds][sb][k];
                    ForCells(L, c, [&](int idx) { O.A.Add(base + idx, 1.); });
                    if (k == 0 || c.fi < 0) continue;
                    const int band = BandOfFine(c.fi);
                    for (int f = 0; f < kNFam; ++f) {
                        const int cell = FamCell(f, c, mee);
                        if (cell >= 0) O.lsbc[ds][sb][f][band * kFamN[f] + cell] += 1.;
                    }
                }
                // replicas of the like-sign rows (sub all)
                if (k > 0 && c.fi >= 0) {
                    if (!haveK) { GroupPois(kSeedCB, k, r.grp, kNRepCB, kr.data()); haveK = true; }
                    double* p = &O.lsFineRep[ds][k - 1][(size_t)c.fi * kNRepCB];
                    for (int q = 0; q < kNRepCB; ++q) p[q] += kr[q];
                    const int band = BandOfFine(c.fi);
                    for (int f = 0; f < kNFam; ++f) {
                        if (!kFamRep[f]) continue;
                        const int cell = FamCell(f, c, mee);
                        if (cell < 0) continue;
                        double* pf = &O.lsbcRep[ds][f][((size_t)band * kFamN[f] + cell) * kNRepCB];
                        for (int q = 0; q < kNRepCB; ++q) pf[q] += kr[q];
                    }
                }
            }
        }
        O.A.Flush();
        printf("wp0: data sample %s: %lld rows in %.1f s\n", scan::SampleName(k), O.nRows[k], sw.RealTime());
        sw.Continue();
    }
}

// CB of every histogram of a data block (dsel, sub): out arrays cb[], var[] in the layout of the block
// (rep: also 200 replicas in the data replica layout). Returns the CB lost in empty bands (0.40-0.70).
struct CBOut { std::vector<double> c, v, c20, v20, rep; int nEmpty = 0; };
void MakeCBBlock(const DataOut& O, int ds, int sb, bool withRep, CBOut& out) {
    const Layout& L = gL[sb == kSubAll ? kLfull : kLstd];
    const double* S = O.A.s.data();
    const int bpp = O.base[ds][sb][1], bmm = O.base[ds][sb][2];
    const double* ppFine = S + bpp + L.off[L.iFine] + 1;
    const double* mmFine = S + bmm + L.off[L.iFine] + 1;
    std::vector<std::vector<double>> F(kNFam);
    for (int f = 0; f < kNFam; ++f) F[f] = Fractions(O.lsbc[ds][sb][f].data(), kFamN[f], &out.nEmpty);
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1 && sb != kSubAll) break;   // CB20 only for the full sample
        CB1D cb;
        cb.Make(ppFine, mmFine, pass == 0 ? kCBreb : kCBreb20);
        std::vector<double>& C = pass == 0 ? out.c : out.c20;
        std::vector<double>& V = pass == 0 ? out.v : out.v20;
        C.assign(L.stride, 0.); V.assign(L.stride, 0.);
        // M4-axis bins of width wf fine bins starting at fine index i0: for each axis bin x (1-based) the fine set
        auto axisSum = [&](int i0, int nb, int wf, std::vector<double>& s, std::vector<double>& var) {
            s.assign(nb + 2, 0.); var.assign(nb + 2, 0.);
            for (int x = 1; x <= nb; ++x) {
                std::vector<std::pair<int, double>> u;
                for (int q = 0; q < wf; ++q) u.push_back(std::make_pair(i0 + (x - 1) * wf + q, 1.));
                const auto p = cb.Sum(u);
                s[x] = p.first; var[x] = p.second;
            }
        };
        std::vector<double> s, var;
        // fine and B_M4
        axisSum(0, kNfine, 1, s, var);
        for (int x = 1; x <= kNfine; ++x) { C[L.off[L.iFine] + x] = s[x]; V[L.off[L.iFine] + x] = var[x]; }
        axisSum(0, kNbefM4, 2, s, var);
        for (int x = 1; x <= kNbefM4; ++x) { C[L.off[L.iBM4] + x] = s[x]; V[L.off[L.iBM4] + x] = var[x]; }
        // F10 / F5g: axis bins from the F10 edges (fine index ranges)
        {
            std::vector<double> fs(kNF10 + 2, 0.), fv(kNF10 + 2, 0.);
            std::vector<int> fband(kNF10 + 2, -1);
            for (int x = 1; x <= kNF10; ++x) {
                std::vector<std::pair<int, double>> u;
                for (int i = 0; i < kNfine; ++i) if (gF10x[i] == x) { u.push_back(std::make_pair(i, 1.)); fband[x] = BandOfFine(i); }
                const auto p = cb.Sum(u);
                fs[x] = p.first; fv[x] = p.second;
            }
            for (int h = 0; h < 10; ++h)
                for (int x = 1; x <= kNF10; ++x) {
                    const double f = F[kFamF10][fband[x] * kFamN[kFamF10] + h];
                    C[L.off[L.iF10[h]] + x] = f * fs[x]; V[L.off[L.iF10[h]] + x] = f * f * fv[x];
                }
            for (int h = 0; h < 5; ++h)
                for (int x = 1; x <= kNF10; ++x) {
                    const double f = F[kFamF5g][fband[x] * kFamN[kFamF5g] + h];
                    C[L.off[L.iF5g[h]] + x] = f * fs[x]; V[L.off[L.iF5g[h]] + x] = f * f * fv[x];
                }
        }
        // X2D: x = 2.5 MeV bins from fine index 89
        {
            const int i0 = (int)std::lround((kX2Dxlo - kFineLo) / kFineW);
            axisSum(i0, kX2Dnx, 1, s, var);
            for (int x = 1; x <= kX2Dnx; ++x) {
                const int band = BandOfFine(i0 + x - 1);
                for (int y = 0; y < kX2Dny + 2; ++y) {
                    const double f = F[kFamX2D][band * kFamN[kFamX2D] + y];
                    const int idx = L.off[L.iX2D] + x + (kX2Dnx + 2) * y;
                    C[idx] = f * s[x]; V[idx] = f * f * var[x];
                }
            }
        }
        // windows: fractional fine bins; X17c and CRK with the band fraction of their extra condition
        for (int w = 0; w < kNWin; ++w) {
            std::vector<std::pair<int, double>> u;
            for (int i = 0; i < kNfine; ++i) {
                double fr = FracFine(i, kWin[w].lo, kWin[w].hi);
                if (fr <= 0) continue;
                const int band = BandOfFine(i);
                if (w == kWX17c || w == kWX17P) fr *= F[kFamX17][band * 2 + 0];
                if (w == kWCRK) fr *= F[kFamCRK][band * 2 + 0];
                u.push_back(std::make_pair(i, fr));
            }
            const auto p = cb.Sum(u);
            C[L.off[L.iWin] + w + 1] = p.first; V[L.off[L.iWin] + w + 1] = p.second;
        }
        // tot: CB in 0.30-1.10
        {
            std::vector<std::pair<int, double>> u;
            for (int i = 0; i < kNfine; ++i) u.push_back(std::make_pair(i, 1.));
            const auto p = cb.Sum(u);
            C[L.off[L.iTot] + 1] = p.first; V[L.off[L.iTot] + 1] = p.second;
        }
        if (sb != kSubAll) continue;
        // G3
        {
            const int i0 = (int)std::lround((kG3xlo - kFineLo) / kFineW);
            axisSum(i0, kG3nx, 1, s, var);
            for (int x = 1; x <= kG3nx; ++x) {
                const int band = BandOfFine(i0 + x - 1);
                for (int yz = 0; yz < kFamN[kFamG3]; ++yz) {
                    const double f = F[kFamG3][band * kFamN[kFamG3] + yz];
                    if (f == 0) continue;
                    const int idx = L.off[L.iG3] + x + (kG3nx + 2) * yz;
                    C[idx] = f * s[x]; V[idx] = f * f * var[x];
                }
            }
        }
        // B_map: 10 MeV x bins
        axisSum(0, kMapNx, 4, s, var);
        for (int x = 1; x <= kMapNx; ++x) {
            const int band = BandOfFine((x - 1) * 4);
            for (int y = 0; y < kMapNy + 2; ++y) {
                const double f = F[kFamMap][band * kFamN[kFamMap] + y];
                const int idx = L.off[L.iMap] + x + (kMapNx + 2) * y;
                C[idx] = f * s[x]; V[idx] = f * f * var[x];
            }
        }
        // BEFORE variables: CB(region) x f(bin | region), f from N++ + N-- of the region histogram
        for (int k = 0; k < kNReg; ++k) {
            std::vector<std::pair<int, double>> u;
            for (int i = 0; i < kNfine; ++i) {
                const double fr = FracFine(i, kWin[kWR0 + k].lo, kWin[kWR0 + k].hi);
                if (fr > 0) u.push_back(std::make_pair(i, fr));
            }
            const auto p = cb.Sum(u);
            for (int vv = 0; vv < kNBVar; ++vv) {
                const int o = L.off[L.iBV[vv][k]], nc = L.h[L.iBV[vv][k]].ncell();
                double t = 0;
                for (int b = 0; b < nc; ++b) t += S[bpp + o + b] + S[bmm + o + b];
                for (int b = 0; b < nc; ++b) {
                    const double f = t > 0 ? (S[bpp + o + b] + S[bmm + o + b]) / t : (b == 0 ? 1. : 0.);
                    C[o + b] = f * p.first; V[o + b] = f * f * p.second;
                }
            }
        }
    }
    if (!withRep) return;
    // ---- replicas (data replica layout: fine, B_M4, F10, F5g, X2D, win)
    const Layout& R = gL[kLrepData];
    out.rep.assign((size_t)R.stride * kNRepCB, 0.);
    std::vector<double> pp(kNfine), mm(kNfine);
    std::vector<double> cnt;
    for (int r = 0; r < kNRepCB; ++r) {
        for (int i = 0; i < kNfine; ++i) {
            pp[i] = O.lsFineRep[ds][0][(size_t)i * kNRepCB + r];
            mm[i] = O.lsFineRep[ds][1][(size_t)i * kNRepCB + r];
        }
        CB1D cb;
        cb.Make(pp.data(), mm.data(), kCBreb);
        std::vector<std::vector<double>> Fr(kNFam);
        for (int f = 0; f < kNFam; ++f) {
            if (!kFamRep[f]) continue;
            cnt.assign(kNBand * kFamN[f], 0.);
            for (int b = 0; b < kNBand * kFamN[f]; ++b) cnt[b] = O.lsbcRep[ds][f][(size_t)b * kNRepCB + r];
            Fr[f] = Fractions(cnt.data(), kFamN[f]);
        }
        auto put = [&](int hidx, int cell, double val) { out.rep[((size_t)R.off[hidx] + cell) * kNRepCB + r] = val; };
        for (int i = 0; i < kNfine; ++i) put(R.iFine, i + 1, cb.cb[i]);
        for (int x = 1; x <= kNbefM4; ++x) put(R.iBM4, x, cb.cb[2 * (x - 1)] + cb.cb[2 * (x - 1) + 1]);
        std::vector<double> fs(kNF10 + 2, 0.);
        std::vector<int> fband(kNF10 + 2, 0);
        for (int i = 0; i < kNfine; ++i) { const int x = gF10x[i]; if (x >= 1 && x <= kNF10) { fs[x] += cb.cb[i]; fband[x] = BandOfFine(i); } }
        for (int h = 0; h < 10; ++h)
            for (int x = 1; x <= kNF10; ++x) put(R.iF10[h], x, Fr[kFamF10][fband[x] * kFamN[kFamF10] + h] * fs[x]);
        for (int h = 0; h < 5; ++h)
            for (int x = 1; x <= kNF10; ++x) put(R.iF5g[h], x, Fr[kFamF5g][fband[x] * kFamN[kFamF5g] + h] * fs[x]);
        const int i0 = (int)std::lround((kX2Dxlo - kFineLo) / kFineW);
        for (int x = 1; x <= kX2Dnx; ++x) {
            const int i = i0 + x - 1, band = BandOfFine(i);
            for (int y = 0; y < kX2Dny + 2; ++y) put(R.iX2D, x + (kX2Dnx + 2) * y, Fr[kFamX2D][band * kFamN[kFamX2D] + y] * cb.cb[i]);
        }
        for (int w = 0; w < kNWin; ++w) {
            double sw = 0;
            for (int i = 0; i < kNfine; ++i) {
                double fr = FracFine(i, kWin[w].lo, kWin[w].hi);
                if (fr <= 0) continue;
                const int band = BandOfFine(i);
                if (w == kWX17c || w == kWX17P) fr *= Fr[kFamX17][band * 2 + 0];
                if (w == kWCRK) fr *= Fr[kFamCRK][band * 2 + 0];
                sw += fr * cb.cb[i];
            }
            put(R.iWin, w + 1, sw);
        }
    }
}

// ===========================================================================
// simulation
// ===========================================================================
struct SimDir { std::string tag; int sel, var, ltype; bool rep, halves; size_t base, repBase; };
const int kNWv = 9;    // weight variants: 0 nom, 1 V2, 2 1aY, 3 PS, 4 convOff, 5 GS, 6 off, 7 modeB, 8 wdb
const int kNSv = 6;    // M4 variants: 0 nom, 1 sK0, 2 sUp, 3 sDn, 4 r089, 5 r110
inline int VarW(int v) {
    switch (v) {
        case kVV2: return 1; case kVbump1aY: return 2; case kVbumpPS: return 3; case kVconvOff: return 4; case kVGS: return 5;
        case kVwtsOff: return 6; case kVmodeB: return 7; case kVwdb: return 8;
    }
    return 0;
}
inline int VarS(int v) {
    switch (v) { case kVsK0: return 1; case kVsUp: return 2; case kVsDn: return 3; case kVr089: return 4; case kVr110: return 5; }
    return 0;
}
const int kPhysOfW[6] = {kPhysNom, kPhysV2, kPhysBump1aY, kPhysBumpPS, kPhysConvOff, kPhysGS};
const double kSofS[kNSv] = {kS, kSK0, kS + kSerr, kS - kSerr, kS, kS};
const double kRofS[kNSv] = {kR, kR, kR, kR, kRlo, kRhi};

// MM4 and y4 of a simulated row recomputed for a beam kinetic energy T [MeV] (corrected momenta, PDG masses)
inline void Kin4(const simeng::Row& r, double T, double& mm4, double& y4) {
    const double mp = simeng::kMpMeV;
    const double E = T + mp, pz = std::sqrt(E * E - mp * mp);
    TLorentzVector ini(0, 0, pz, E + mp);
    const TLorentzVector f = simeng::P4(r.ep_p, r.ep_theta, r.ep_phi, simeng::kMeMeV) + simeng::P4(r.em_p, r.em_theta, r.em_phi, simeng::kMeMeV) +
                             simeng::P4(r.pip_p, r.pip_theta, r.pip_phi, simeng::kMpiMeV) +
                             simeng::P4(r.pim_p, r.pim_theta, r.pim_phi, simeng::kMpiMeV);
    const TLorentzVector miss = ini - f;
    const double m2 = miss.M2();
    mm4 = (m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2)) * 1e-3;
    y4 = f.Rapidity() - ini.Rapidity();
}

struct SimOut {
    std::vector<SimDir> dirs;
    Acc A;
    RepAcc R;
    // checks (pre, nominal weights, unscaled M4)
    double c23SimBump = 0, c23SimBumpV = 0, c45SimPeak = 0, mm4diffMax = 0, y4diffMax = 0;
    long long nRows = 0, nPassAny = 0, nNoComp = 0, nTmpl[5] = {0, 0, 0, 0, 0}, nKinChk = 0;
};

void BuildSimDirs(SimOut& O) {
    const int selMain[2] = {kSelS0, kSelPre};
    for (int s = 0; s < 2; ++s)
        for (int v = 0; v < kNSimVar; ++v) {
            SimDir d;
            d.tag = std::string(SelTag(selMain[s])) + "__" + SimVarTag(v);
            d.sel = selMain[s]; d.var = v;
            d.ltype = v == kVnom ? kLfull : (v == kVghostA || v == kVwtsOff) ? kLbef : kLstd;
            d.rep = d.halves = (v == kVnom);
            O.dirs.push_back(d);
        }
    const int selVar[4] = {kSelDz40, kSelPeak, kSelVz150, kSelVz170};
    for (int s = 0; s < 4; ++s) {
        SimDir d;
        d.tag = std::string(SelTag(selVar[s])) + "__nom";
        d.sel = selVar[s]; d.var = kVnom; d.ltype = kLstd; d.rep = d.halves = false;
        O.dirs.push_back(d);
    }
    size_t tot = 0, rtot = 0;
    for (auto& d : O.dirs) {
        d.base = tot; tot += (size_t)gL[d.ltype].stride * kNComp;
        d.repBase = rtot; if (d.rep) rtot += (size_t)gL[kLrepSim].stride * kNComp;
    }
    printf("wp0: sim %zu directories, %zu cells (%.1f MB per array), replica cells %zu (%.1f MB)\n", O.dirs.size(), tot,
           tot * 8e-6, rtot, rtot * kNRepSim * 8e-6);
    O.A.Init(tot, true);
    O.R.Init(rtot, kNRepSim);
}

void FillSim(const simeng::SimData& D, SimOut& O, long long maxRows) {
    TStopwatch sw;
    std::vector<double> kr(kNRepSim, 0.);
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    const int nd = (int)O.dirs.size();
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& v = D.t[trees[t]];
        const size_t n = maxRows > 0 ? std::min((size_t)maxRows, v.size()) : v.size();
        int cur = -1;
        bool haveK = false;
        for (size_t i = 0; i < n; ++i) {
            const simeng::Row& r = v[i];
            if (r.grp != cur) { O.A.Begin(r.grp % 2 == 0); cur = r.grp; haveK = false; }
            O.nRows++;
            if (r.tmpl >= 0 && r.tmpl <= 4) O.nTmpl[r.tmpl]++;
            // selections (sim rows: cor_ok == 1 always; MM4 at 4500 MeV for PEAK)
            const float dzpi = DzPi(r.pip_z, r.pim_z, r.vz);
            double mm4s = r.mm4, y4s = r.y4;
            bool pass[kNSel];
            bool any = false;
            const bool base = PassBase(r, r.cor_ok == 1);
            if (base) {
                Kin4(r, kTsim, mm4s, y4s);
                if (O.nKinChk < 200000) {   // validation of the recomputation at 4530 MeV against the stored values
                    double m30, y30;
                    Kin4(r, kTdata, m30, y30);
                    O.mm4diffMax = std::max(O.mm4diffMax, std::fabs(m30 - r.mm4));
                    O.y4diffMax = std::max(O.y4diffMax, std::fabs(y30 - r.y4));
                    O.nKinChk++;
                }
            }
            for (int s = 0; s < kNSel; ++s) { pass[s] = base && s != kSelNone && PassSel(s, r, r.cor_ok == 1, dzpi, (float)mm4s); any |= pass[s]; }
            if (!any) continue;
            O.nPassAny++;
            int cn[6], cg[6];
            const int ncn = CompsOf(r, false, cn), ncg = CompsOf(r, true, cg);
            if (ncn == 0 && ncg == 0) { O.nNoComp++; continue; }
            // weights
            const bool oA = r.ownA == 1, oB = r.ownB == 1;
            double wv[kNWv];
            for (int q = 0; q < 6; ++q) wv[q] = oA ? r.w * simphys::RowFactor(r, kPhysOfW[q]) : 0.;
            wv[6] = oA ? r.w : 0.;
            const double f108 = (oA || oB) ? simphys::RowFactor(r, kPhysNom) : 0.;
            wv[7] = oB ? r.w * f108 : 0.;
            wv[8] = oA ? r.w_db * f108 : 0.;
            const double c6A = oA ? 2. : 0., c6B = oB ? 2. : 0.;
            // acceptance-1 checks (pre, nominal, s = 1)
            if (pass[kSelPre]) {
                for (int q = 0; q < ncn; ++q) {
                    if ((cn[q] == kC2 || cn[q] == kC3) && r.m4 >= 0.4576f && r.m4 < 0.5255f) { O.c23SimBump += wv[0]; }
                    if ((cn[q] == kC4 || cn[q] == kC5) && r.m4 >= 0.5333f && r.m4 < 0.5587f) O.c45SimPeak += wv[0];
                }
            }
            // cells per M4 variant (lazy)
            const double mpp = kSPiPi * r.mpipi;
            const double bvals[kNBVar] = {mpp, r.mee, r.oa, r.vz, mm4s, r.pt4, y4s};
            Cells cs[kNSv];
            bool have[kNSv] = {false, false, false, false, false, false};
            for (int d = 0; d < nd; ++d) {
                const SimDir& dir = O.dirs[d];
                if (!pass[dir.sel]) continue;
                const int si = VarS(dir.var), wi = VarW(dir.var);
                if (!have[si]) { MakeCells(M4sim(r.m4, r.m4_true, kSofS[si], kRofS[si]), mpp, r.mee, bvals, cs[si]); have[si] = true; }
                const Cells& c = cs[si];
                const bool gh = dir.var == kVghostA;
                const int nc = gh ? ncg : ncn;
                const int* cl = gh ? cg : cn;
                const Layout& L = gL[dir.ltype];
                for (int q = 0; q < nc; ++q) {
                    const int comp = cl[q];
                    const double w = comp == kC6 ? (wi == 7 ? c6B : c6A) : wv[wi];
                    if (w == 0.) continue;
                    const size_t b = dir.base + (size_t)comp * L.stride;
                    ForCells(L, c, [&](int idx) { O.A.Add((int)(b + idx), w); });
                    if (dir.rep) {
                        if (!haveK) { GroupPois(kSeedSim, trees[t], r.grp, kNRepSim, kr.data()); haveK = true; }
                        const Layout& Rl = gL[kLrepSim];
                        const size_t rb = dir.repBase + (size_t)comp * Rl.stride;
                        ForCells(Rl, c, [&](int idx) { O.R.Add(rb + idx, w, kr.data()); });
                    }
                }
            }
        }
        O.A.Flush();
        printf("wp0: sim tree %s: %zu rows, %.1f s\n", simeng::TreeName(trees[t]), n, sw.RealTime());
        sw.Continue();
    }
}

// ===========================================================================
// acceptance helpers
// ===========================================================================
struct Chk { std::string what; double val, ref, tol; bool rel; bool ok; };
std::vector<Chk> gChk;
void AddChk(const std::string& w, double val, double ref, double tol, bool rel) {
    Chk c{w, val, ref, tol, rel, false};
    const double d = rel ? (ref != 0 ? std::fabs(val - ref) / std::fabs(ref) : std::fabs(val)) : std::fabs(val - ref);
    c.ok = d <= tol;
    gChk.push_back(c);
    printf("wp0 check %-80s %18.10g ref %18.10g %s\n", w.c_str(), val, ref, c.ok ? "OK" : "FAIL");
}
// sum of the in-range cells of a 1D histogram (bins lo..hi, 1-based) of a block
inline double Sum1(const double* a, size_t off, int b1, int b2) { double s = 0; for (int b = b1; b <= b2; ++b) s += a[off + b]; return s; }

// cell-sum closure of one block (contents a[], layout L): returns the largest relative deviation
double Closure(const double* a, const Layout& L, size_t base, std::string& worst) {
    double dmax = 0;
    auto test = [&](const std::string& n, double x, double ref) {
        const double d = std::fabs(x - ref) / std::max(1.0, std::fabs(ref));
        if (d > dmax) { dmax = d; worst = n; }
    };
    const size_t of = base + L.off[L.iFine];
    // F10: sum over the 10 histograms and 70 bins = fine 0.40-0.70 (fine bins 41..160)
    double s10 = 0;
    for (int h = 0; h < 10; ++h) s10 += Sum1(a, base + L.off[L.iF10[h]], 1, kNF10);
    test("F10", s10, Sum1(a, of, 41, 160));
    // F10 per M4 bin and F5g per bin <= F10 (F5g is a subset): check F5g x-sum against G3 later
    // X2D (with the M(ee) flows): sum over all y of x bins 1..35 = fine 0.5225-0.6100 (fine bins 90..124)
    double sx = 0;
    for (int y = 0; y < kX2Dny + 2; ++y) for (int x = 1; x <= kX2Dnx; ++x) sx += a[base + L.off[L.iX2D] + x + (kX2Dnx + 2) * y];
    test("X2D", sx, Sum1(a, of, 90, 124));
    // B_M4 = fine
    test("B_M4", Sum1(a, base + L.off[L.iBM4], 1, kNbefM4), Sum1(a, of, 1, kNfine));
    // aligned windows: B, Fit, Bef, CRU, CRL, Jlo, Jhi
    const int aw[7] = {kWB, kWFit, kWBef, kWCRU, kWCRL, kWJlo, kWJhi};
    for (int k = 0; k < 7; ++k) {
        const int w = aw[k];
        const int b1 = (int)std::lround((kWin[w].lo - kFineLo) / kFineW) + 1, b2 = (int)std::lround((kWin[w].hi - kFineLo) / kFineW);
        test(std::string("win ") + kWin[w].tag, a[base + L.off[L.iWin] + w + 1], Sum1(a, of, b1, b2));
    }
    // CRK window = S5 F10 bins 0.400-0.5225
    {
        double s = 0;
        for (int e = 0; e < 2; ++e) for (int x = 1; x <= kNF10; ++x) if (gF10e[x] <= 0.5225 + 1e-9) s += a[base + L.off[L.iF10[8 + e]] + x];
        test("win CRK vs F10", a[base + L.off[L.iWin] + kWCRK + 1], s);
    }
    if (L.iG3 >= 0) {
        // G3 (all y, z) per x bin summed over 0.40-0.85 = fine 41..220
        double sg = 0, s5g = 0, s10g = 0;
        const int ny = (int)gG3ey.size() + 1, nz = (int)gG3ez.size() + 1;
        for (int z = 0; z < nz; ++z) for (int y = 0; y < ny; ++y) for (int x = 1; x <= kG3nx; ++x) {
            const double c = a[base + L.off[L.iG3] + x + (kG3nx + 2) * (y + ny * z)];
            sg += c;
            if (z == 5 && x <= 120) s5g += c;    // M(ee) >= 20 MeV, M4 0.40-0.70
            if (x <= 120) s10g += c;
        }
        test("G3", sg, Sum1(a, of, 41, 220));
        double s5 = 0;
        for (int h = 0; h < 5; ++h) s5 += Sum1(a, base + L.off[L.iF5g[h]], 1, kNF10);
        test("F5g vs G3", s5, s5g);
        test("F10 vs G3", s10, s10g);
    }
    if (L.hasBef) {
        // map (with M(pipi) flows) = B_M4 0.30-1.10
        double sm = 0;
        for (int y = 0; y < kMapNy + 2; ++y) for (int x = 1; x <= kMapNx; ++x) sm += a[base + L.off[L.iMap] + x + (kMapNx + 2) * y];
        test("B_map", sm, Sum1(a, of, 1, kNfine));
        // BEFORE variables in each region (with flows) = window R_k
        for (int k = 0; k < kNReg; ++k)
            for (int v = 0; v < kNBVar; ++v) {
                const HShape& s = L.h[L.iBV[v][k]];
                test(BVarName(v, k), Sum1(a, base + L.off[L.iBV[v][k]], 0, s.nx() + 1), a[base + L.off[L.iWin] + kWR0 + k + 1]);
            }
    }
    return dmax;
}

}  // namespace wp0

// ===========================================================================
// main
// ===========================================================================
void wp0_fill(long long maxRows = -1) {
    using namespace wp0;
    TStopwatch total;
    InitEdges();
    for (int t = 0; t < kNLType; ++t) gL[t] = MakeLayout(t);
    printf("wp0: layout strides std %d, bef %d, full %d, repSim %d, repData %d\n", gL[kLstd].stride, gL[kLbef].stride,
           gL[kLfull].stride, gL[kLrepSim].stride, gL[kLrepData].stride);
    // edge function checks (Sect. 2)
    AddChk("E(0.30, 0) = 0.5186", Edge(0.30, 0.), 0.5186, 0.00006, false);
    AddChk("E(0.39, 0) = 0.5008", Edge(0.39, 0.), 0.5008, 0.00006, false);
    AddChk("E(0.39, 0.060) = 0.5096", Edge(0.39, 0.060), 0.5096, 0.00006, false);
    AddChk("E(0.4116, 0) = 0.4739", Edge(0.4116, 0.), 0.4739, 0.0003, false);
    AddChk("E(0.30, 0.100) = 0.533", Edge(0.30, 0.100), 0.533, 0.0006, false);

    const bool test = maxRows > 0;
    const std::string outName = test ? std::string("/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1_test.root") : std::string(kHistFile);
    const std::string sealName = test ? std::string("/hdd3/przygoda/HADES/FAT_out/all_report/sealed/hists_v1_x17_sealed_test.root")
                                      : std::string(kSealedFile);
    gSystem->mkdir("/hdd3/przygoda/HADES/FAT_out/all_report/sealed", kTRUE);
    TFile* fo = TFile::Open(outName.c_str(), "RECREATE");
    TFile* fs = TFile::Open(sealName.c_str(), "RECREATE");
    if (!fo || fo->IsZombie() || !fs || fs->IsZombie()) { printf("wp0: cannot open the output files\n"); return; }
    fo->SetCompressionSettings(505);
    FILE* tab = fopen(test ? "/hdd3/przygoda/HADES/FAT_out/all_report/wp0/tables/wp0_numbers_test.txt"
                           : "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp0_numbers.txt", "w");
    fprintf(tab, "# wp0_numbers.txt -- plots/all_report/wp0_fill.C (WP0 histogram factory), %s\n", test ? "TEST RUN" : "full run");
    fprintf(tab, "# output %s\n", outName.c_str());

    // ------------------------------------------------------------------ data
    {
        scan::Data D;
        if (!scan::Load(D, scan::kCompact)) return;
        DataOut* O = new DataOut();
        FillData(D, *O, maxRows);
        D.s[0].clear(); D.s[0].shrink_to_fit(); D.s[1].clear(); D.s[1].shrink_to_fit(); D.s[2].clear(); D.s[2].shrink_to_fit();
        TDirectory* dd = fo->mkdir("data");
        TDirectory* sd = fs->mkdir("data");
        // masks (blinding, D19)
        std::vector<double> maskX(gL[kLstd].h[gL[kLstd].iX2D].ncell(), 0.), maskB(kBVnb[kBmee] + 2, 0.);
        std::vector<int> maskYbins, maskBbins;
        for (int y = 1; y <= kX2Dny; ++y) if (gX2Dey[y] > 0.013 + 1e-12 && gX2Dey[y - 1] < 0.027 - 1e-12) maskYbins.push_back(y);
        for (int b = 1; b <= kBVnb[kBmee]; ++b) if (gBVe[kBmee][b] > 0.013 + 1e-12 && gBVe[kBmee][b - 1] < 0.027 - 1e-12) maskBbins.push_back(b);
        for (int y : maskYbins) for (int x = 0; x <= kX2Dnx + 1; ++x) maskX[x + (kX2Dnx + 2) * y] = 1.;
        for (int b : maskBbins) maskB[b] = 1.;
        WriteObj(dd, MakeHist(gL[kLstd].h[gL[kLstd].iX2D], "mask_X2D", maskX.data(), nullptr, 1.,
                              "1 = data N masked (-1) until x17_frozen.txt exists (D19)"));
        WriteObj(dd, MakeHist(gL[kLfull].h[gL[kLfull].iBV[kBmee][3]], "mask_B_mee", maskB.data(), nullptr, 1.,
                              "1 = data N of B_mee_R2, R3, R4 masked (-1) (D19)"));
        fprintf(tab, "\n[data]\n# rows read: us %lld, pp %lld, mm %lld\n", O->nRows[0], O->nRows[1], O->nRows[2]);
        double worstAll = 0;
        std::string worstName;
        for (int ds = 0; ds < kNDSel; ++ds) {
            TDirectory* d1 = dd->mkdir(kDSel[ds].tag);
            TDirectory* s1 = sd->mkdir(kDSel[ds].tag);
            d1->WriteTObject(new TNamed("info", Form("data, selection %s: %s. Contents: raw rows (no 64/63, no L); N = e+e- rows, Npp, Nmm = "
                                                     "like-sign rows, errors = event-group (sum over groups of k^2); CB = 2 sqrt(N++ N--) "
                                                     "(10 MeV cells, interpolated, spread by f(cell | 50 MeV band)), k = 1; CB20 = 20 MeV cells; "
                                                     "CBrep_* = 200 Poisson event-group replicas of the like-sign rows (y/z axis = replica).",
                                                     kDSel[ds].tag, kDSel[ds].text)), "info", "Overwrite");
            for (int sb = 0; sb < kNSub; ++sb) {
                TDirectory* d2 = d1->mkdir(SubTag(sb));
                const Layout& L = gL[sb == kSubAll ? kLfull : kLstd];
                CBOut cb;
                MakeCBBlock(*O, ds, sb, sb == kSubAll, cb);
                // closure checks (before masking)
                for (int k = 0; k < 3; ++k) {
                    std::string wn;
                    const double dm = Closure(O->A.s.data(), L, O->base[ds][sb][k], wn);
                    if (dm > worstAll) { worstAll = dm; worstName = Form("data %s %s %s %s", kDSel[ds].tag, SubTag(sb), kWhat[k], wn.c_str()); }
                }
                {
                    std::string wn;
                    const double dm = Closure(cb.c.data(), L, 0, wn);
                    if (dm > worstAll) { worstAll = dm; worstName = Form("data %s %s CB %s", kDSel[ds].tag, SubTag(sb), wn.c_str()); }
                }
                for (int hi = 0; hi < (int)L.h.size(); ++hi) {
                    const HShape& s = L.h[hi];
                    for (int k = 0; k < 3; ++k) {
                        const size_t o = O->base[ds][sb][k] + L.off[hi];
                        std::vector<double> c(O->A.s.begin() + o, O->A.s.begin() + o + s.ncell());
                        std::vector<double> v(O->A.v.begin() + o, O->A.v.begin() + o + s.ncell());
                        const bool isX = (hi == L.iX2D), isB = (L.hasBef && (hi == L.iBV[kBmee][2] || hi == L.iBV[kBmee][3] || hi == L.iBV[kBmee][4]));
                        if (k == 0 && (isX || isB)) {
                            // sealed copy, then mask
                            TDirectory* s2 = s1->GetDirectory(SubTag(sb)) ? s1->GetDirectory(SubTag(sb)) : s1->mkdir(SubTag(sb));
                            WriteObj(s2, MakeHist(s, std::string("N_") + s.name, c.data(), v.data(), 1.));
                            const std::vector<double>& m = isX ? maskX : maskB;
                            for (int q = 0; q < s.ncell(); ++q) if (m[q] > 0) { c[q] = -1.; v[q] = 0.; }
                        }
                        WriteObj(d2, MakeHist(s, std::string(kWhat[k]) + "_" + s.name, c.data(), v.data(), 1.));
                    }
                    WriteObj(d2, MakeHist(s, std::string("CB_") + s.name, cb.c.data() + L.off[hi], cb.v.data() + L.off[hi], 1.));
                    if (sb == kSubAll) WriteObj(d2, MakeHist(s, std::string("CB20_") + s.name, cb.c20.data() + L.off[hi], cb.v20.data() + L.off[hi], 1.));
                }
                if (sb == kSubAll) {
                    const Layout& R = gL[kLrepData];
                    for (int hi = 0; hi < (int)R.h.size(); ++hi)
                        WriteObj(d2, MakeRepHist(R.h[hi], std::string("CBrep_") + R.h[hi].name, cb.rep.data() + (size_t)R.off[hi] * kNRepCB, kNRepCB, 1.));
                    // replica RMS vs linear propagation (windows)
                    fprintf(tab, "# %s CB windows: value, linear error, replica mean, replica rms\n", kDSel[ds].tag);
                    for (int w = 0; w < kNWin; ++w) {
                        double m = 0, m2 = 0;
                        for (int r = 0; r < kNRepCB; ++r) { const double x = cb.rep[((size_t)R.off[R.iWin] + w + 1) * kNRepCB + r]; m += x; m2 += x * x; }
                        m /= kNRepCB; m2 = std::sqrt(std::max(0., m2 / kNRepCB - m * m));
                        const size_t on = O->base[ds][sb][0] + L.off[L.iWin] + w + 1;
                        fprintf(tab, "%-10s %-5s N %10.0f (+- %7.1f)  CB %10.2f lin +- %7.2f  rep mean %10.2f rms %7.2f  N-CB %10.1f\n", kDSel[ds].tag,
                                kWin[w].tag, O->A.s[on], std::sqrt(O->A.v[on]), cb.c[L.off[L.iWin] + w + 1], std::sqrt(cb.v[L.off[L.iWin] + w + 1]), m, m2,
                                O->A.s[on] - cb.c[L.off[L.iWin] + w + 1]);
                    }
                    fprintf(tab, "# %s: rows passing us %lld, pp %lld, mm %lld; rows without an F10 cell %lld; empty CB bands (all families) %d\n",
                            kDSel[ds].tag, O->nPass[ds][0], O->nPass[ds][1], O->nPass[ds][2], O->nNaNcell[ds], cb.nEmpty);
                    // acceptance 2 and 3
                    const size_t wP = O->base[ds][sb][0] + L.off[L.iWin] + kWP + 1, wB = O->base[ds][sb][0] + L.off[L.iWin] + kWB + 1;
                    if (ds == 0) {
                        AddChk("A2: data S0 peak-window N = 55947", O->A.s[wP], 55947, 0.5, false);
                        AddChk("A2: data S0 peak-window CB/N = 0.059", cb.c[L.off[L.iWin] + kWP + 1] / O->A.s[wP], 0.059, 0.0005, false);
                    }
                    if (std::string(kDSel[ds].tag) == "none")
                        AddChk("A3: data none bump-window W = N - CB = 148091", O->A.s[wB] - cb.c[L.off[L.iWin] + kWB + 1], 148091, 0.5, false);
                }
            }
            printf("wp0: data %s written (%.1f s)\n", kDSel[ds].tag, total.RealTime());
            total.Continue();
        }
        AddChk("A4: data cell sums = 1D (max relative deviation; " + worstName + ")", worstAll, 0, 1e-9, false);
        delete O;
    }

    // ------------------------------------------------------------------ simulation
    {
        simphys::Init();
        simphys::Install();
        simphys::LoadAngles();
        simeng::SimData* D = new simeng::SimData();
        if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
        SimOut* O = new SimOut();
        BuildSimDirs(*O);
        FillSim(*D, *O, maxRows);
        {
            ProcInfo_t pi;
            gSystem->GetProcInfo(&pi);
            printf("wp0: memory after the simulation pass: resident %.2f GB, virtual %.2f GB\n", pi.fMemResident / 1048576., pi.fMemVirtual / 1048576.);
        }
        delete D;
        printf("wp0: sim filled; rows %lld, passing a selection %lld, without a component %lld; tmpl 0-4: %lld %lld %lld %lld %lld; "
               "MM4/y4 recomputation at 4530 vs stored: max |dMM4| %.2e GeV/c2, max |dy| %.2e (%lld rows); angles hit %lld miss %lld\n",
               O->nRows, O->nPassAny, O->nNoComp, O->nTmpl[0], O->nTmpl[1], O->nTmpl[2], O->nTmpl[3], O->nTmpl[4], O->mm4diffMax,
               O->y4diffMax, O->nKinChk, simphys::AngHit(), simphys::AngMiss());
        fprintf(tab, "\n[simulation]\n# rows %lld, passing a selection %lld, no component %lld; tmpl 0-4: %lld %lld %lld %lld %lld\n", O->nRows,
                O->nPassAny, O->nNoComp, O->nTmpl[0], O->nTmpl[1], O->nTmpl[2], O->nTmpl[3], O->nTmpl[4]);
        fprintf(tab, "# MM4, y4 recomputed at 4530 MeV vs stored: max |dMM4| = %.3e GeV/c2, max |dy4| = %.3e (%lld rows)\n", O->mm4diffMax,
                O->y4diffMax, O->nKinChk);
        TDirectory* sdir = fo->mkdir("sim");
        const double sc = 1. / kLsim;
        double worstAll = 0;
        std::string worstName;
        for (const SimDir& dir : O->dirs) {
            TDirectory* d1 = sdir->mkdir(dir.tag.c_str());
            d1->WriteTObject(new TNamed("info", Form("simulation, selection %s, variant %s: %s. Units: nb per bin = sum w / L_sim "
                                                     "(L_sim = 221.05 nb^-1), errors = event groups; expected raw data rows = nb x %.4f "
                                                     "(L = %.1f nb^-1, x 63/64). M4 = M4' (mass transfer), M(pipi) x %.4f. %s",
                                                     SelTag(dir.sel), SimVarTag(dir.var), SimVarText(dir.var), kNbToRaw, kL, kSPiPi,
                                                     dir.rep ? "Halves _hA/_hB (even/odd event group, NOT rescaled) and 100 replicas _rep." : "")),
                             "info", "Overwrite");
            const Layout& L = gL[dir.ltype];
            for (int comp = 0; comp < kNComp; ++comp) {
                const size_t b = dir.base + (size_t)comp * L.stride;
                std::string wn;
                const double dm = Closure(O->A.s.data(), L, b, wn);
                if (dm > worstAll) { worstAll = dm; worstName = "sim " + dir.tag + " " + CompTag(comp) + " " + wn; }
                for (int hi = 0; hi < (int)L.h.size(); ++hi) {
                    const HShape& s = L.h[hi];
                    const size_t o = b + L.off[hi];
                    const std::string nm = std::string(CompTag(comp)) + "_" + s.name;
                    WriteObj(d1, MakeHist(s, nm, O->A.s.data() + o, O->A.v.data() + o, sc));
                    if (dir.halves && hi != L.iG3 && hi != L.iMap && !(L.hasBef && hi >= L.iBV[0][0] && hi <= L.iBV[kNBVar - 1][kNReg - 1])) {
                        std::vector<double> sB(s.ncell()), vB(s.ncell());
                        for (int q = 0; q < s.ncell(); ++q) { sB[q] = O->A.s[o + q] - O->A.sA[o + q]; vB[q] = O->A.v[o + q] - O->A.vA[o + q]; }
                        WriteObj(d1, MakeHist(s, nm + "_hA", O->A.sA.data() + o, O->A.vA.data() + o, sc, "half A (even event groups)"));
                        WriteObj(d1, MakeHist(s, nm + "_hB", sB.data(), vB.data(), sc, "half B (odd event groups)"));
                    }
                }
                if (dir.rep) {
                    const Layout& R = gL[kLrepSim];
                    const size_t rb = dir.repBase + (size_t)comp * R.stride;
                    for (int hi = 0; hi < (int)R.h.size(); ++hi)
                        WriteObj(d1, MakeRepHist(R.h[hi], std::string(CompTag(comp)) + "_" + R.h[hi].name + "_rep",
                                                 O->R.r.data() + (rb + R.off[hi]) * kNRepSim, kNRepSim, sc));
                }
            }
            // window and total table
            fprintf(tab, "# sim %s [nb]: comp, tot (all M4), windows B P SBlo SBhi X17c CRU CRK\n", dir.tag.c_str());
            for (int comp = 0; comp < kNComp; ++comp) {
                const size_t b = dir.base + (size_t)comp * L.stride;
                const size_t ot = b + L.off[L.iTot] + 1, ow = b + L.off[L.iWin];
                fprintf(tab, "%-14s %-10s %12.5g +- %9.3g |", dir.tag.c_str(), CompTag(comp), O->A.s[ot] * sc, std::sqrt(O->A.v[ot]) * sc);
                const int ws[7] = {kWB, kWP, kWSBlo, kWSBhi, kWX17c, kWCRU, kWCRK};
                for (int q = 0; q < 7; ++q) fprintf(tab, " %10.4g", O->A.s[ow + ws[q] + 1] * sc);
                fprintf(tab, "\n");
            }
            if (dir.tag == "pre__nom") {
                auto tot = [&](int comp) { return O->A.s[dir.base + (size_t)comp * L.stride + L.off[L.iTot] + 1] * sc; };
                AddChk("A1: sim pre T1+T2 all M4 = 701.083 nb", tot(kT12), 701.083, 0.0005, true);
                AddChk("A1: sim pre C2+C3 in sim bump window 0.4576-0.5255 (s = 1) = 9.689 nb", O->c23SimBump * sc, 9.689, 0.0005, true);
                AddChk("extra: sim pre C4+C5 in sim peak window (s = 1) = 2.262 nb", O->c45SimPeak * sc, 2.262, 0.001, true);
                AddChk("extra: sim pre C2 all M4 = 7.046 nb", tot(kC2), 7.046, 0.001, true);
                AddChk("extra: sim pre C3 all M4 = 5.19 nb", tot(kC3), 5.19, 0.002, true);
                AddChk("extra: sim pre C4 all M4 = 1.242 nb", tot(kC4), 1.242, 0.001, true);
                AddChk("extra: sim pre C5 all M4 = 1.417 nb", tot(kC5), 1.417, 0.001, true);
                AddChk("extra: sim pre G_eta all M4 = 24.27 nb", tot(kGeta), 24.27, 0.001, true);
                AddChk("extra: sim pre G_pi0D all M4 = 394.5 nb", tot(kGpi0D), 394.5, 0.001, true);
                AddChk("extra: sim pre G_conv all M4 = 265.3 nb", tot(kGconv), 265.3, 0.001, true);
                AddChk("extra: sim pre G_vec all M4 = 0.1066 nb", tot(kGvec), 0.1066, 0.001, true);
                AddChk("extra: sim pre G_oth all M4 = 1.967 nb", tot(kGoth), 1.967, 0.001, true);
                AddChk("extra: sim pre T3 all M4 = 125.0 nb", tot(kT3), 125.0, 0.001, true);
                AddChk("extra: sim pre T4 all M4 = 12.51 nb", tot(kT4), 12.51, 0.001, true);
                AddChk("extra: sim pre X17 x 1e-5 all M4 = 0.07127 nb", tot(kC6) * 1e-5, 0.07127, 0.001, true);
                AddChk("extra: sim pre G_pi0D = om + nom", tot(kGpi0Dom) + tot(kGpi0Dnom), tot(kGpi0D), 1e-9, true);
                AddChk("extra: sim pre G_D = G_pi0D + G_vec + G_oth", tot(kGD), tot(kGpi0D) + tot(kGvec) + tot(kGoth), 1e-9, true);
                AddChk("extra: sim pre T12 = C2..C5 + G_eta + G_D + G_conv", tot(kT12),
                       tot(kC2) + tot(kC3) + tot(kC4) + tot(kC5) + tot(kGeta) + tot(kGD) + tot(kGconv), 1e-9, true);
            }
            printf("wp0: sim %s written (%.1f s)\n", dir.tag.c_str(), total.RealTime());
            total.Continue();
        }
        AddChk("A4: sim cell sums = 1D (max relative deviation; " + worstName + ")", worstAll, 0, 1e-9, false);
        delete O;
    }

    // ------------------------------------------------------------------ README and checks
    fo->cd();
    fo->WriteTObject(new TNamed("README", "hists_v1.root of SPEC report_all (WP0, plots/all_report/wp0_fill.C). Naming scheme: see "
                                          "plots/all_report/README.md. data/<sel>/<sub>/<what>_<layout>, what = N, Npp, Nmm, CB, CB20, "
                                          "CBrep (raw rows); sim/<sel>__<variant>/<comp>_<layout>[_hA|_hB|_rep] (nb per bin). Layouts: "
                                          "F10_S<j>E<e>, F5g_S<j>, X2D, G3, fine_M4, B_M4, B_map, B_<var>_R<k>, win, tot. Data N of X2D "
                                          "(M(ee) 12-28 MeV) and B_mee_R2-R4 (12-28 MeV) masked (-1) until x17_frozen.txt (D19)."),
                     "README", "Overwrite");
    fprintf(tab, "\n[acceptance and checks]\n");
    int nfail = 0;
    for (const auto& c : gChk) {
        fprintf(tab, "%-4s %-96s value %16.8g ref %16.8g tol %g%s\n", c.ok ? "PASS" : "FAIL", c.what.c_str(), c.val, c.ref, c.tol, c.rel ? " (rel)" : "");
        if (!c.ok) ++nfail;
    }
    fprintf(tab, "# %d checks, %d failed; run time %.0f s\n", (int)gChk.size(), nfail, total.RealTime());
    fclose(tab);
    fo->Close();
    fs->Close();
    printf("wp0: done, %d checks, %d failed, %.0f s; output %s\n", (int)gChk.size(), nfail, total.RealTime(), outName.c_str());
}
