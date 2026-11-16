// ============================================================================
// flex_fill.C  --  package flex of report_all: histograms of the "X" layout for the flexible SMASH description
// ============================================================================
// Why: the flexible fit (flex_fit.C) runs nominally on F10 E2 of hists_v1.root (M4 0.40-0.70). Two of its
// variants need cells that hists_v1.root does not have:
//   (1) the M4 range extended to 1.00 GeV/c2 (F10 bins 0.40-0.70, then 10 MeV bins 0.70-1.00: 100 bins);
//   (2) the M(ee) class E2 split into E2a = 10-20 MeV and E2b >= 20 MeV (with E1 < 10 MeV: three classes),
//       with ONE event-group covariance over all cells (wp5_egcov.root has F10 and F5g as separate matrices).
// Layout X: histogram h = 3 j + c (slice j = 0..4: S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5 >= 0.43;
// class c = 0 E1 < 10 MeV, 1 E2a 10-20 MeV, 2 E2b >= 20 MeV), M4 axis: the 70 F10 bins + 30 bins of 10 MeV.
//
// Conventions: those of wp0_fill.C (Part III), followed line by line (code copied from wp0_fill.C / e1_fill.C):
//   data : rows of scan_v2 (us, pp, mm), selection S0 = alldefs::PassSel(kSelS0); M4 bin = axis bin of the 2.5 MeV fine
//          bin (fine grid 0.30-1.10), exactly as wp0 (gF10x); N with event-group variances and the full event-group
//          covariance (rows grouped by adjacency of grp, as wp5_rows.C); CB = 2 sqrt(N++ N--) in 10 MeV cells,
//          interpolated to the fine bins, spread over the 15 (slice x class) cells by f(cell | 50 MeV band) from N++ + N--;
//          200 Poisson(1) event-group replicas of the like-sign rows (wp0 generator and seed kSeedCB).
//          Units: RAW ROWS (no 64/63, no L).
//   sim  : SMASH std + lep + eta e+e- rows, owned, mode A, w x simphys::RowFactor(r, 108); M4' = M4sim(m4, m4_true)
//          (s = 0.99615, r = 1), M(pipi)' = 0.9952 M(pipi); components alldefs::CompsOf; nb per bin (sum w / L_sim) with
//          event-group errors; 100 Poisson(1) event-group replicas (wp0 generator and seed kSeedSim, so the replicas of
//          the 0.40-0.70 cells are those of hists_v1.root).
// Checks (written to flex_fill_checks.txt): the X cells summed to E1 / E2 over 0.40-0.70 reproduce hists_v1.root F10
// (data N, N errors, CB, CB replicas; sim contents, errors, replicas) and F5g (E2b).
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_data.root, flex_sim.root
// Run (from /home/przygoda/HADES/FAT/FAT):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_flex", kTRUE)'
//   root -l -b -q -e "$B" 'plots/all_report/flex_fill.C+O(0)'    # compile
//   root -l -b -q -e "$B" 'plots/all_report/flex_fill.C+O(1)'    # data (about 2 min, 4 GB)
//   root -l -b -q -e "$B" 'plots/all_report/flex_fill.C+O(2)'    # simulation (about 3 min, 9 GB)
//   root -l -b -q -e "$B" 'plots/all_report/flex_fill.C+O(4)'    # checks against hists_v1.root
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixDSym.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "defs.h"
#include "hists_io.h"
#include "sim_eta_physics.h"

namespace flexf {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/flex";

// ---------------------------------------------------------------- random numbers (copy of wp0_fill.C)
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

// ---------------------------------------------------------------- layout X
const int kNX = 100;       // M4 bins per histogram
const int kNCls = 3;       // E1, E2a, E2b
const int kNH = kNSlice * kNCls;
const int kNCell = kNH * kNX;
inline const char* ClsTag(int c) { static const char* n[kNCls] = {"E1", "E2a", "E2b"}; return n[c]; }
inline std::string XName(int j, int c) { return Form("X_S%d%s", j + 1, ClsTag(c)); }
std::vector<double> gXe;   // 101 edges
std::vector<int> gXx;      // fine bin -> axis bin 1..100 (0: below, kNX + 1: above)
void InitLayout() {
    gXe = F10Edges();
    for (int i = 1; i <= 30; ++i) gXe.push_back(std::round((0.700 + 0.010 * i) * 1e6) / 1e6);
    gXx.assign(kNfine, 0);
    for (int i = 0; i < kNfine; ++i) {
        const double c = kFineLo + (i + 0.5) * kFineW;
        if (c < gXe.front()) gXx[i] = 0;
        else if (c >= gXe.back()) gXx[i] = kNX + 1;
        else gXx[i] = (int)(std::upper_bound(gXe.begin(), gXe.end(), c) - gXe.begin());
    }
}
inline int Cls(double mee) {
    if (!(mee == mee)) return -1;
    return mee < kEsplit ? 0 : (mee < kF5gMee ? 1 : 2);
}
// cell index (0..kNCell-1) of a row, -1 if outside; fineOut = fine bin (for the CB)
inline int CellOf(double m4, double mpipi, double mee, int& fineOut, int& hOut) {
    fineOut = FineBin(m4);
    hOut = -1;
    const int j = Slice(mpipi), c = Cls(mee);
    if (j < 0 || c < 0) return -1;
    hOut = kNCls * j + c;
    if (fineOut < 0) return -1;
    const int x = gXx[fineOut];
    if (x < 1 || x > kNX) return -1;
    return hOut * kNX + (x - 1);
}

// ---------------------------------------------------------------- accumulator with event-group variance (wp0 Acc)
struct Acc {
    std::vector<double> s, v, acc;
    std::vector<int> stamp, touched;
    int gid = 0;
    void Init(size_t n) { s.assign(n, 0.); v.assign(n, 0.); acc.assign(n, 0.); stamp.assign(n, 0); touched.reserve(1 << 16); }
    inline void Begin() { Flush(); ++gid; }
    inline void Add(size_t i, double w) {
        s[i] += w;
        if (stamp[i] != gid) { stamp[i] = gid; acc[i] = 0.; touched.push_back((int)i); }
        acc[i] += w;
    }
    void Flush() { for (int i : touched) v[i] += acc[i] * acc[i]; touched.clear(); }
};

// ---------------------------------------------------------------- CB (copy of wp0 CB1D)
struct CB1D {
    int reb = 4, nw = 80;
    std::vector<double> cbW, varW, cb;
    std::vector<int> j0, j1;
    std::vector<double> a0, a1;
    void Make(const double* pp, const double* mm, int r) {
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
std::vector<double> Fractions(const double* cnt, int n) {
    std::vector<double> f(kNBand * n, 0.), pool(n, 0.);
    double ptot = 0;
    for (int b = 0; b < kNBand; ++b) for (int c = 0; c < n; ++c) { pool[c] += cnt[b * n + c]; ptot += cnt[b * n + c]; }
    for (int b = 0; b < kNBand; ++b) {
        double t = 0;
        for (int c = 0; c < n; ++c) t += cnt[b * n + c];
        for (int c = 0; c < n; ++c) f[b * n + c] = t > 0 ? cnt[b * n + c] / t : (ptot > 0 ? pool[c] / ptot : (c == 0 ? 1. : 0.));
    }
    return f;
}
// CB of every X cell from a CB1D and the fractions f[band * kNH + h]
void CellCB(const CB1D& cb, const std::vector<double>& f, std::vector<double>& out, std::vector<double>* var) {
    out.assign(kNCell, 0.);
    if (var) var->assign(kNCell, 0.);
    for (int h = 0; h < kNH; ++h)
        for (int x = 1; x <= kNX; ++x) {
            std::vector<std::pair<int, double>> u;
            for (int i = 0; i < kNfine; ++i) {
                if (gXx[i] != x) continue;
                const double fr = f[BandOfFine(i) * kNH + h];
                if (fr > 0) u.push_back(std::make_pair(i, fr));
            }
            const auto p = cb.Sum(u);
            out[h * kNX + x - 1] = p.first;
            if (var) (*var)[h * kNX + x - 1] = p.second;
        }
}

TH1D* H1(const std::string& n, const double* c, const double* v, double sc) {
    TH1D* h = new TH1D(n.c_str(), "", kNX, gXe.data());
    h->SetDirectory(nullptr); h->Sumw2();
    for (int x = 1; x <= kNX; ++x) { h->SetBinContent(x, c[x - 1] * sc); h->SetBinError(x, v ? std::sqrt(std::max(0., v[x - 1])) * sc : 0.); }
    return h;
}
TH2D* HRep(const std::string& n, const double* rep, int nrep, double sc) {   // rep[(x - 1) * nrep + r]
    TH2D* h = new TH2D(n.c_str(), "replicas (y = replica index)", kNX, gXe.data(), nrep, -0.5, nrep - 0.5);
    h->SetDirectory(nullptr);
    for (int x = 1; x <= kNX; ++x) for (int r = 0; r < nrep; ++r) h->SetBinContent(x, r + 1, rep[(size_t)(x - 1) * nrep + r] * sc);
    return h;
}

// ===========================================================================
// data
// ===========================================================================
void RunData() {
    TStopwatch sw;
    scan::Data D;
    if (!scan::Load(D, scan::kCompact)) return;
    const int NR = kNRepCB;
    Acc A; A.Init(kNCell);
    std::vector<double> lsAcc[2];                    // like-sign rows per cell (information)
    std::vector<double> ls[2], lsRep[2];             // fine like-sign counts and replicas [i * NR + r]
    std::vector<double> fam(kNBand * kNH, 0.), famRep((size_t)kNBand * kNH * NR, 0.);
    for (int q = 0; q < 2; ++q) { ls[q].assign(kNfine, 0.); lsRep[q].assign((size_t)kNfine * NR, 0.); lsAcc[q].assign(kNCell, 0.); }
    std::vector<double> eg((size_t)kNCell * kNCell, 0.);
    std::vector<int> cells;
    long long nPass[3] = {0, 0, 0}, nCell[3] = {0, 0, 0};
    std::vector<double> kr(NR, 0.);
    for (int k = 0; k < 3; ++k) {
        const std::vector<scan::Row>& v = D.s[k];
        int cur = -1;
        bool haveK = false;
        auto flush = [&]() {
            for (size_t a = 0; a < cells.size(); ++a) for (size_t b = 0; b < cells.size(); ++b) eg[(size_t)cells[a] * kNCell + cells[b]] += 1.;
            cells.clear();
        };
        for (size_t i = 0; i < v.size(); ++i) {
            const scan::Row& r = v[i];
            if (r.grp != cur) { if (k == 0) { A.Begin(); flush(); } cur = r.grp; haveK = false; }
            if (!(r.m4 == r.m4)) continue;
            if (!PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4)) continue;
            nPass[k]++;
            int fi, h;
            const int c = CellOf(r.m4, r.mpipi, r.mee, fi, h);
            if (k == 0) {
                if (c >= 0) { A.Add(c, 1.); cells.push_back(c); nCell[0]++; }
                continue;
            }
            if (c >= 0) { lsAcc[k - 1][c] += 1.; nCell[k]++; }
            if (fi < 0) continue;
            ls[k - 1][fi] += 1.;
            const int band = BandOfFine(fi);
            if (h >= 0) fam[band * kNH + h] += 1.;
            if (!haveK) { GroupPois(kSeedCB, k, r.grp, NR, kr.data()); haveK = true; }
            double* p = &lsRep[k - 1][(size_t)fi * NR];
            for (int q = 0; q < NR; ++q) p[q] += kr[q];
            if (h >= 0) { double* pf = &famRep[((size_t)band * kNH + h) * NR]; for (int q = 0; q < NR; ++q) pf[q] += kr[q]; }
        }
        if (k == 0) { A.Flush(); flush(); }
        printf("flex: data sample %s: %lld rows pass S0, %lld in X cells (%.1f s)\n", scan::SampleName(k), nPass[k], nCell[k], sw.RealTime());
        sw.Continue();
    }
    D.s[0].clear(); D.s[1].clear(); D.s[2].clear();
    // CB
    CB1D cb; cb.Make(ls[0].data(), ls[1].data(), kCBreb);
    std::vector<double> C, V;
    CellCB(cb, Fractions(fam.data(), kNH), C, &V);
    std::vector<double> rep((size_t)kNCell * NR, 0.), pp(kNfine), mm(kNfine), cnt(kNBand * kNH), Cr;
    for (int r = 0; r < NR; ++r) {
        for (int i = 0; i < kNfine; ++i) { pp[i] = lsRep[0][(size_t)i * NR + r]; mm[i] = lsRep[1][(size_t)i * NR + r]; }
        CB1D cbr; cbr.Make(pp.data(), mm.data(), kCBreb);
        for (int b = 0; b < kNBand * kNH; ++b) cnt[b] = famRep[(size_t)b * NR + r];
        CellCB(cbr, Fractions(cnt.data(), kNH), Cr, nullptr);
        for (int q = 0; q < kNCell; ++q) rep[(size_t)q * NR + r] = Cr[q];
    }
    TFile* fo = TFile::Open(Form("%s/flex_data.root", kOutDir), "RECREATE");
    for (int j = 0; j < kNSlice; ++j)
        for (int c = 0; c < kNCls; ++c) {
            const int h = kNCls * j + c, o = h * kNX;
            fo->WriteTObject(H1("N_" + XName(j, c), A.s.data() + o, A.v.data() + o, 1.));
            fo->WriteTObject(H1("Npp_" + XName(j, c), lsAcc[0].data() + o, lsAcc[0].data() + o, 1.));
            fo->WriteTObject(H1("Nmm_" + XName(j, c), lsAcc[1].data() + o, lsAcc[1].data() + o, 1.));
            fo->WriteTObject(H1("CB_" + XName(j, c), C.data() + o, V.data() + o, 1.));
            fo->WriteTObject(HRep("CBrep_" + XName(j, c), rep.data() + (size_t)o * NR, NR, 1.));
        }
    TMatrixDSym M(kNCell);
    for (int a = 0; a < kNCell; ++a) for (int b = 0; b < kNCell; ++b) M(a, b) = eg[(size_t)a * kNCell + b];
    fo->WriteTObject(&M, "egcov_X");
    fo->Close();
    printf("flex: data written (%.1f s)\n", sw.RealTime());
}

// ===========================================================================
// simulation
// ===========================================================================
void RunSim() {
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    const int NR = kNRepSim;
    Acc A; A.Init((size_t)kNComp * kNCell);
    Acc A2; A2.Init((size_t)kNComp * kNSlice * kNX);   // merged class E2 = E2a + E2b (its own event-group variance)
    std::vector<double> R((size_t)kNComp * kNCell * NR, 0.);
    std::vector<double> kr(NR, 0.);
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    long long nPass = 0;
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& v = D->t[trees[t]];
        int cur = -1;
        bool haveK = false;
        for (size_t i = 0; i < v.size(); ++i) {
            const simeng::Row& r = v[i];
            if (r.grp != cur) { A.Begin(); A2.Begin(); cur = r.grp; haveK = false; }
            if (!PassSel(kSelS0, r, r.cor_ok == 1, 0.f, 0.f)) continue;
            if (r.ownA != 1) continue;   // mode A
            int cn[8];
            const int ncn = CompsOf(r, false, cn);
            if (ncn == 0) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom);
            int fi, h;
            const int c = CellOf(M4sim(r.m4, r.m4_true), kSPiPi * r.mpipi, r.mee, fi, h);
            if (c < 0) continue;
            nPass++;
            for (int q = 0; q < ncn; ++q) {
                const double wc = cn[q] == kC6 ? 2. : w;
                if (wc == 0.) continue;
                A.Add((size_t)cn[q] * kNCell + c, wc);
                if (h % kNCls != 0) A2.Add((size_t)cn[q] * kNSlice * kNX + (h / kNCls) * kNX + (c % kNX), wc);
                if (!haveK) { GroupPois(kSeedSim, trees[t], r.grp, NR, kr.data()); haveK = true; }
                double* p = &R[((size_t)cn[q] * kNCell + c) * NR];
                for (int k = 0; k < NR; ++k) p[k] += wc * kr[k];
            }
        }
        A.Flush(); A2.Flush();
        printf("flex: sim tree %s: %zu rows (%.1f s)\n", simeng::TreeName(trees[t]), v.size(), sw.RealTime());
        sw.Continue();
    }
    printf("flex: angles hit %lld miss %lld; rows in X cells %lld\n", simphys::AngHit(), simphys::AngMiss(), nPass);
    delete D;
    TFile* fo = TFile::Open(Form("%s/flex_sim.root", kOutDir), "RECREATE");
    const double sc = 1. / kLsim;
    for (int comp = 0; comp < kNComp; ++comp)
        for (int j = 0; j < kNSlice; ++j)
            for (int c = 0; c < kNCls; ++c) {
                const size_t o = (size_t)comp * kNCell + (kNCls * j + c) * kNX;
                const std::string n = std::string(CompTag(comp)) + "_" + XName(j, c);
                fo->WriteTObject(H1(n, A.s.data() + o, A.v.data() + o, sc));
                fo->WriteTObject(HRep(n + "_rep", R.data() + o * NR, NR, sc));
            }
    // merged E2 (contents = E2a + E2b; error = event-group error of the merged cell, as hists_v1 F10 E2)
    for (int comp = 0; comp < kNComp; ++comp)
        for (int j = 0; j < kNSlice; ++j) {
            const size_t o = (size_t)comp * kNSlice * kNX + j * kNX;
            fo->WriteTObject(H1(std::string(CompTag(comp)) + Form("_X_S%dE2", j + 1), A2.s.data() + o, A2.v.data() + o, sc));
        }
    fo->Close();
    printf("flex: sim written (%.1f s)\n", sw.RealTime());
}

// ===========================================================================
// checks against hists_v1.root
// ===========================================================================
void Check() {
    TFile* fh = TFile::Open(kHistFile);
    TFile* fd = TFile::Open(Form("%s/flex_data.root", kOutDir));
    TFile* fs = TFile::Open(Form("%s/flex_sim.root", kOutDir));
    FILE* out = fopen(Form("%s/flex_fill_checks.txt", kOutDir), "w");
    auto P = [&](const char* fmt, double a, double b, const char* what) {
        fprintf(out, fmt, what, a, b); printf(fmt, what, a, b);
    };
    double worst[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // data: F10 E1 = X E1, F10 E2 = X E2a + X E2b, F5g = X E2b, bins 1..70
    for (int j = 0; j < kNSlice; ++j) {
        for (int e = 0; e < 3; ++e) {   // 0 F10 E1, 1 F10 E2, 2 F5g
            const std::string ref = e < 2 ? F10Name(j, e) : F5gName(j);
            TH1* n = allio::Data(fh, "S0", "all", "N", ref);
            TH1* cb = allio::Data(fh, "S0", "all", "CB", ref);
            TH2* cr = (TH2*)fh->Get(Form("data/S0/all/CBrep_%s", ref.c_str()));
            std::vector<int> cl = e == 0 ? std::vector<int>{0} : e == 1 ? std::vector<int>{1, 2} : std::vector<int>{2};
            for (int x = 1; x <= kNF10; ++x) {
                double sn = 0, scb = 0, sv = 0;
                std::vector<double> srep(kNRepCB, 0.);
                for (int c : cl) {
                    TH1* a = (TH1*)fd->Get(("N_" + XName(j, c)).c_str());
                    TH1* b = (TH1*)fd->Get(("CB_" + XName(j, c)).c_str());
                    TH2* r = (TH2*)fd->Get(("CBrep_" + XName(j, c)).c_str());
                    sn += a->GetBinContent(x); scb += b->GetBinContent(x);
                    if (cl.size() == 1) sv = a->GetBinError(x) * a->GetBinError(x);
                    for (int k = 0; k < kNRepCB; ++k) srep[k] += r->GetBinContent(x, k + 1);
                }
                worst[0] = std::max(worst[0], std::fabs(sn - n->GetBinContent(x)));
                worst[1] = std::max(worst[1], std::fabs(scb - cb->GetBinContent(x)) / std::max(1., cb->GetBinContent(x)));
                if (cl.size() == 1) worst[2] = std::max(worst[2], std::fabs(sv - n->GetBinError(x) * n->GetBinError(x)));
                if (cr && e < 2) for (int k = 0; k < kNRepCB; ++k) worst[3] = std::max(worst[3], std::fabs(srep[k] - cr->GetBinContent(x, k + 1)) / std::max(1., cr->GetBinContent(x, k + 1)));
            }
            delete n; delete cb;
        }
    }
    // event-group covariance diagonal vs N errors
    TMatrixDSym* M = (TMatrixDSym*)fd->Get("egcov_X");
    for (int h = 0; h < kNH; ++h) {
        TH1* a = (TH1*)fd->Get(("N_" + XName(h / kNCls, h % kNCls)).c_str());
        for (int x = 1; x <= kNX; ++x) worst[4] = std::max(worst[4], std::fabs((*M)(h * kNX + x - 1, h * kNX + x - 1) - a->GetBinError(x) * a->GetBinError(x)));
    }
    // event-group covariance of F10 (wp5_egcov.root) = X summed
    TFile* f5 = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/wp5_egcov.root");
    TMatrixDSym* M10 = (TMatrixDSym*)f5->Get("egcov_S0_all_F10");
    {
        auto xcells = [&](int h10, int x) { std::vector<int> v; const int j = h10 / 2, e = h10 % 2; if (e == 0) v.push_back((kNCls * j) * kNX + x); else { v.push_back((kNCls * j + 1) * kNX + x); v.push_back((kNCls * j + 2) * kNX + x); } return v; };
        for (int a = 0; a < 700; ++a)
            for (int b = 0; b < 700; ++b) {
                const auto va = xcells(a / 70, a % 70), vb = xcells(b / 70, b % 70);
                double s = 0;
                for (int p : va) for (int q : vb) s += (*M)(p, q);
                worst[5] = std::max(worst[5], std::fabs(s - (*M10)(a, b)));
            }
    }
    // simulation
    for (int comp = 0; comp < kNComp; ++comp)
        for (int j = 0; j < kNSlice; ++j)
            for (int e = 0; e < 2; ++e) {
                TH1* ref = allio::Sim(fh, "S0", "nom", comp, F10Name(j, e));
                TH2* rr = (TH2*)fh->Get(Form("sim/S0__nom/%s_%s_rep", CompTag(comp), F10Name(j, e).c_str()));
                std::vector<int> cl = e == 0 ? std::vector<int>{0} : std::vector<int>{1, 2};
                for (int x = 1; x <= kNF10; ++x) {
                    double s = 0, v = 0;
                    std::vector<double> sr(kNRepSim, 0.);
                    for (int c : cl) {
                        TH1* a = (TH1*)fs->Get((std::string(CompTag(comp)) + "_" + XName(j, c)).c_str());
                        TH2* r = (TH2*)fs->Get((std::string(CompTag(comp)) + "_" + XName(j, c) + "_rep").c_str());
                        s += a->GetBinContent(x);
                        if (cl.size() == 1) v = a->GetBinError(x) * a->GetBinError(x);
                        for (int k = 0; k < kNRepSim; ++k) sr[k] += r->GetBinContent(x, k + 1);
                    }
                    if (cl.size() == 2) {
                        TH1* m2 = (TH1*)fs->Get((std::string(CompTag(comp)) + Form("_X_S%dE2", j + 1)).c_str());
                        v = m2->GetBinError(x) * m2->GetBinError(x);
                        worst[6] = std::max(worst[6], std::fabs(m2->GetBinContent(x) - s) / std::max(1e-6, ref->GetBinContent(x)));
                    }
                    const double sc = std::max(1e-6, ref->GetBinContent(x));
                    worst[6] = std::max(worst[6], std::fabs(s - ref->GetBinContent(x)) / sc);
                    worst[6] = std::max(worst[6], std::fabs(v - ref->GetBinError(x) * ref->GetBinError(x)) / (sc * sc));
                    if (rr) for (int k = 0; k < kNRepSim; ++k) worst[7] = std::max(worst[7], std::fabs(sr[k] - rr->GetBinContent(x, k + 1)) / std::max(1e-6, rr->GetBinContent(x, k + 1)));
                }
                delete ref;
            }
    fprintf(out, "# flex_fill.C checks against hists_v1.root (X layout, 0.40-0.70, S0)\n");
    P("%-60s %.3g (ref %.3g)\n", worst[0], 0, "data N: max |X - hists_v1| (F10 E1, E2, F5g) [rows]");
    P("%-60s %.3g (ref %.3g)\n", worst[1], 0, "data CB: max relative difference");
    P("%-60s %.3g (ref %.3g)\n", worst[2], 0, "data N error^2: max |X - hists_v1| (E1, F5g)");
    P("%-60s %.3g (ref %.3g)\n", worst[3], 0, "data CB replicas: max relative difference (F10)");
    P("%-60s %.3g (ref %.3g)\n", worst[4], 0, "egcov_X diagonal vs N error^2");
    P("%-60s %.3g (ref %.3g)\n", worst[5], 0, "egcov_X summed to F10 vs wp5_egcov egcov_S0_all_F10");
    P("%-60s %.3g (ref %.3g)\n", worst[6], 0, "sim contents and errors^2 (all comps, F10): max relative");
    P("%-60s %.3g (ref %.3g)\n", worst[7], 0, "sim replicas (all comps, F10): max relative");
    fclose(out);
}

}  // namespace flexf

void flex_fill(int mode = 0) {
    flexf::InitLayout();
    gSystem->mkdir(flexf::kOutDir, kTRUE);
    printf("flex: X layout %d bins x %d histograms, edges %.4f ... %.4f\n", flexf::kNX, flexf::kNH, flexf::gXe.front(), flexf::gXe.back());
    if (mode & 1) flexf::RunData();
    if (mode & 2) flexf::RunSim();
    if (mode & 4) flexf::Check();
}
