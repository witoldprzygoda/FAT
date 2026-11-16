// ============================================================================
// e1_fill.C  --  package e1 of report_all: data and SMASH after the "typical" eta cuts (E1-E4) and the
//                X17-oriented OA < 6 deg cuts (X1-X4), histogram factory (prefix e1)
// ============================================================================
// Why: the analyst asked to see the standard comparison data / SMASH AFTER the typical cuts of the earlier
// analyses (slides / old code; their purpose is to enhance a possible eta -> pi+pi- X17 channel), and with
// an extra OA(e+e-) < 6 deg cut (X17 -> e+e- pairs are expected at small opening angles).
//
// Conventions: those of wp0_fill.C (Part III), followed line by line:
//   data : rows of scan_v2 (us, pp, mm); pre = alldefs::PassSel(kSelPre), S0 = alldefs::PassSel(kSelS0);
//          CB = 2 sqrt(N++ N--) in 10 MeV M4 cells (fine grid 0.30-1.10 in 2.5 MeV), linear interpolation to the
//          fine bins, windows with fractional fine bins (wp0 CB1D / FracFine); CB of a cell of another variable
//          (M(ee), OA) in the peak band = sum over the fine M4 bins of CB(fine) x f(cell | 50 MeV band), f = share
//          of the like-sign rows (N++ + N--) of the band in the cell (wp0 X2D rule); k = 1; errors: N event-group
//          (sum over groups of k^2), CB linear and RMS of 200 Poisson(1) event-group replicas of the like-sign rows
//          (wp0 generator and seed). Units: RAW ROWS (no 64/63, no L).
//   sim  : SMASH std + lep + eta e+e- rows, owned, mode A, w x simphys::RowFactor(r, 108) (V3 + KLOE + conv F_V
//          + bkg rules); M4' = alldefs::M4sim(m4, m4_true) (s = 0.99615, r = 1); M(pipi)' = 0.9952 M(pipi) (data
//          scale, used in the M(pipi) cuts); M(ee), OA, oa_lab, oa_eta unchanged; components alldefs::CompsOf
//          (C2..C5, C6 = X17 per unit BR with weight 2 per owned row, G_eta, G_pi0D, G_conv, G_vec, G_oth, G_D,
//          T12, T3, T4); units nb = sum w / L_sim (221.05 nb^-1), event-group errors.
//   old  : reproduction of the earlier study (sim_rw_compare.C, tables/rw_compare_enriched.txt): data
//          scan::FillMany with its FillOpt (0.35-1.00, 260 bins, reb 4), (N - CB) / 6080.5; sim T1 + T2
//          (simeng::CompSim) at the UNSCALED M4 and the raw M(pipi), w x factor 108, mode A.
// Selections: base (pre | S0) x cut set:
//   ""    none                 E1  oa_lab < 50, oa_eta > 140, M(pipi) < 0.46, OA > 2 (slides)
//   E2    E1 without OA > 2    E3  oa_lab < 50, oa_eta > 140, M(pipi) < 0.42 (old code)     E3oa E3 + OA > 2
//   E4    oa_lab < 50, oa_eta > 140 (eta band)                 E4oa E4 + OA > 2
//   X1    E2 + OA < 6          X2  E2 + 2 < OA < 6 (= E1 kinematic cuts + 2 < OA < 6)
//   O6    OA < 6 (pre_O6 = X4, S0_O6 = X3)                     X2old E3 + 2 < OA < 6 (old code with its OA flags)
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/e1/e1_<data|sim>.json (+ e1_<data|sim>.root)
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_e1", kTRUE)' 'plots/all_report/e1_fill.C+O(0)'  (compile)
//   ... 'plots/all_report/e1_fill.C+O(1)'   data (about 40 s with the load)
//   ... 'plots/all_report/e1_fill.C+O(2)'   simulation (about 2 min, about 8 GB); the two can run in parallel
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "defs.h"
#include "sim_eta_physics.h"

namespace e1 {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/e1";

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
const uint64_t kSeedCB = 0xCB2026ULL;

// ---------------------------------------------------------------- selections
const int kNCut = 11;
const char* const kCutTag[kNCut] = {"", "E1", "E2", "E3", "E3oa", "E4", "E4oa", "X1", "X2", "O6", "X2old"};
const char* const kCutText[kNCut] = {
    "no further cut",
    "E1 (slides): oa_lab < 50, oa_eta > 140, M(pipi) < 0.46, OA > 2",
    "E2: E1 without OA > 2",
    "E3 (old code): oa_lab < 50, oa_eta > 140, M(pipi) < 0.42",
    "E3oa: E3 + OA > 2",
    "E4 (eta band): oa_lab < 50, oa_eta > 140",
    "E4oa: E4 + OA > 2",
    "X1: E2 + OA < 6",
    "X2: E2 + 2 < OA < 6",
    "O6: OA < 6 (pre: X4, S0: X3)",
    "X2old: E3 + 2 < OA < 6"};
inline bool PassCut(int c, float mpp, float oa, float oalab, float oaeta) {
    if (c == 0) return true;
    if (c == 9) return oa < 6.f;
    if (!(oalab < 50.f && oaeta > 140.f)) return false;
    switch (c) {
        case 1: return mpp < 0.46f && oa > 2.f;
        case 2: return mpp < 0.46f;
        case 3: return mpp < 0.42f;
        case 4: return mpp < 0.42f && oa > 2.f;
        case 5: return true;
        case 6: return oa > 2.f;
        case 7: return mpp < 0.46f && oa < 6.f;
        case 8: return mpp < 0.46f && oa > 2.f && oa < 6.f;
        case 10: return mpp < 0.42f && oa > 2.f && oa < 6.f;
    }
    return false;
}
const int kNSelE = 2 * kNCut;   // 0..10 pre base, 11..21 S0 base
inline std::string SelName(int s) {
    const int b = s / kNCut, c = s % kNCut;
    std::string n = b == 0 ? "pre" : "S0";
    if (c) n += std::string("_") + kCutTag[c];
    return n;
}
// old study: pre, E1, E2, E3, E3oa, E4, E4oa (raw M(pipi))
const int kNOld = 7;
const int kOldCut[kNOld] = {0, 1, 2, 3, 4, 5, 6};
const char* const kOldTag[kNOld] = {"pre", "E1", "E2", "E3", "E3oa", "E4", "E4oa"};
// old regions (seen::kReg + bands of sim_rw_compare.C)
const int kNOldReg = 10;
const char* const kOldRegTag[kNOldReg] = {"below", "bump", "gap", "peak", "post", "r56", "r70", "all", "band049", "band056"};
const double kOldRegLo[kNOldReg] = {0.35, 0.4550, 0.5225, scan::kPeakLo, scan::kPeakHi, 0.56, 0.70, 0.35, 0.4875, 0.56};
const double kOldRegHi[kNOldReg] = {0.455, 0.5225, scan::kPeakLo, scan::kPeakHi, 0.56, 0.70, 1.00, 1.00, 0.5225, 0.595};

// ---------------------------------------------------------------- windows (by the row's exact M4; CB with fractional fine bins)
struct Win { const char* tag; double lo, hi; const char* text; };
const int kNW = 20;
const Win kW[kNW] = {
    {"below", 0.35, 0.455, "below the bump 0.35-0.455"},
    {"R0", 0.0, 0.455, "M4 < 0.455 (Table B1 R0)"},
    {"B", 0.4550, 0.5225, "bump window B"},
    {"R2", 0.5225, 0.53364, "0.5225-0.5336"},
    {"P", 0.53364, 0.55416, "peak window P"},
    {"post", 0.55416, 0.56, "0.5542-0.56"},
    {"R4", 0.55416, 0.565, "0.5542-0.565"},
    {"r56", 0.56, 0.70, "0.56-0.70"},
    {"R5", 0.565, 0.70, "0.565-0.70 (CR-U)"},
    {"CR", 0.575, 0.70, "0.575-0.70 (background normalisation of the scaled display)"},
    {"r70", 0.70, 1.00, "0.70-1.00"},
    {"R6", 0.70, 0.90, "0.70-0.90"},
    {"R7", 0.90, 1.10, "0.90-1.10"},
    {"Fit", 0.40, 0.70, "0.40-0.70 (display fit range)"},
    {"all", 0.35, 1.00, "0.35-1.00"},
    {"PB", 0.5325, 0.5575, "eta peak band (x17d Fig. X1)"},
    {"X2D", 0.5225, 0.6100, "X2D M4 range (X17 limit range)"},
    {"Jlo", 0.4875, 0.5225, "jump denominator band"},
    {"Jhi", 0.560, 0.595, "jump numerator band"},
    {"Bef", 0.30, 1.10, "0.30-1.10"}};
const int kWPB = 15, kWX2D = 16;
// peak band on the fine grid: fine bins 93..102 (0.5325-0.5575)
const int kPBi0 = 93, kPBi1 = 102;
// X17 M(ee) range 15.5-24.5 MeV; X2D M(ee) range 2-62 MeV
const double kXmeeLo = 0.0155, kXmeeHi = 0.0245, kX2DmeeLo = 0.002, kX2DmeeHi = 0.062;

// ---------------------------------------------------------------- cell layout of one (selection, component) block
// fine 0.30-1.10 2.5 MeV (322 with flows), b5 0.30-1.10 5 MeV (162), win (kNW), x17 window (1), x2dm (X2D M4 and M(ee) 2-62, 1),
// mee in the peak band (32: 2 MeV 0-62, overflow), oa in the peak band (41: 1 deg 0-40, overflow), oaX = oa in the X17 window (41),
// pee = pair momentum in the peak band (sim only; 101: 20 MeV/c 0-2000, overflow), meeW = M(ee) 0-62 MeV in 0.40-0.70 (32)
const int kNfineC = kNfine + 2, kNb5 = 160, kNb5C = kNb5 + 2, kNmee = 32, kNoa = 41, kNpee = 101;
enum Part { pFine = 0, pB5, pWin, pX17, pX2Dm, pMee, pOa, pOaX, pPee, pMeeW, kNPart };
const int kPartN[kNPart] = {kNfineC, kNb5C, kNW, 1, 1, kNmee, kNoa, kNoa, kNpee, kNmee};
const char* const kPartTag[kNPart] = {"fine", "b5", "win", "x17", "x2dm", "mee", "oa", "oaX", "pee", "meeW"};
int gOff[kNPart], gStride = 0;
void InitLayout() {
    gStride = 0;
    for (int p = 0; p < kNPart; ++p) { gOff[p] = gStride; gStride += kPartN[p]; }
}
inline int MeeBin(double mee) { if (!(mee == mee) || mee < 0) return -1; const int j = (int)(mee / 0.002); return std::min(j, kNmee - 1); }
inline int OaBin(double oa) { if (!(oa == oa) || oa < 0) return -1; const int j = (int)oa; return std::min(j, kNoa - 1); }
inline int PeeBin(double p) { if (!(p == p) || p < 0) return -1; const int j = (int)(p / 20.); return std::min(j, kNpee - 1); }
inline int B5Bin(int fi, double m4) { return fi >= 0 ? fi / 2 + 1 : (m4 < kFineLo ? 0 : kNb5 + 1); }
inline int FineCell(int fi, double m4) { return fi >= 0 ? fi + 1 : (m4 < kFineLo ? 0 : kNfine + 1); }

// list of (part offset + index) cells of a row (pee only when p >= 0)
inline int CellsOf(double m4, double mee, double oa, double pee, int* out) {
    int n = 0;
    const int fi = FineBin(m4);
    out[n++] = gOff[pFine] + FineCell(fi, m4);
    out[n++] = gOff[pB5] + B5Bin(fi, m4);
    for (int w = 0; w < kNW; ++w) if (m4 >= kW[w].lo && m4 < kW[w].hi) out[n++] = gOff[pWin] + w;
    const bool pb = m4 >= kW[kWPB].lo && m4 < kW[kWPB].hi;
    const bool xm = mee >= kXmeeLo && mee < kXmeeHi;
    if (pb && xm) out[n++] = gOff[pX17];
    if (m4 >= kW[kWX2D].lo && m4 < kW[kWX2D].hi && mee >= kX2DmeeLo && mee < kX2DmeeHi) out[n++] = gOff[pX2Dm];
    if (pb) {
        const int jm = MeeBin(mee), jo = OaBin(oa);
        if (jm >= 0) out[n++] = gOff[pMee] + jm;
        if (jo >= 0) out[n++] = gOff[pOa] + jo;
        if (jo >= 0 && xm) out[n++] = gOff[pOaX] + jo;
        if (pee >= 0) { const int jp = PeeBin(pee); if (jp >= 0) out[n++] = gOff[pPee] + jp; }
    }
    if (m4 >= 0.40 && m4 < 0.70) { const int jm = MeeBin(mee); if (jm >= 0) out[n++] = gOff[pMeeW] + jm; }
    return n;
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
    void Flush() {
        for (int i : touched) v[i] += acc[i] * acc[i];
        touched.clear();
    }
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
inline double FracFine(int i, double a, double b) {
    double fa = (a - kFineLo) / kFineW, fb = (b - kFineLo) / kFineW;
    if (std::fabs(fa - std::round(fa)) < 1e-6) fa = std::round(fa);
    if (std::fabs(fb - std::round(fb)) < 1e-6) fb = std::round(fb);
    const double o = std::min((double)(i + 1), fb) - std::max((double)i, fa);
    return o > 0 ? std::min(1., o) : 0.;
}
// families of the CB spreading in the peak band: mee (32 cells), oa (41), x17 (2: in / out of 15.5-24.5 MeV), oaX (42: OA bin
// in the X17 M(ee) range, 41 = outside), x2dm (2: M(ee) 2-62 / outside)
enum Fam { fMee = 0, fOa, fX17, fOaX, fX2Dm, kNFam };
const int kFamN[kNFam] = {kNmee, kNoa, 2, kNoa + 1, 2};
inline int FamCell(int f, double mee, double oa) {
    const bool xm = mee >= kXmeeLo && mee < kXmeeHi;
    switch (f) {
        case fMee: return MeeBin(mee);
        case fOa: return OaBin(oa);
        case fX17: return xm ? 0 : 1;
        case fOaX: { const int j = OaBin(oa); return j < 0 ? -1 : (xm ? j : kNoa); }
        case fX2Dm: return (mee >= kX2DmeeLo && mee < kX2DmeeHi) ? 0 : 1;
    }
    return -1;
}
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

// ---------------------------------------------------------------- json writer
struct JW {
    FILE* f = nullptr;
    bool first = true;
    void Open(const std::string& n) { f = fopen(n.c_str(), "w"); }
    void Arr(const char* key, const double* a, int n, bool comma = true) {
        fprintf(f, "%s\"%s\": [", comma ? ", " : "", key);
        for (int i = 0; i < n; ++i) fprintf(f, i ? ",%.10g" : "%.10g", std::isfinite(a[i]) ? a[i] : 0.);
        fprintf(f, "]");
    }
};

// ===========================================================================
// data
// ===========================================================================
void RunData() {
    TStopwatch sw;
    scan::Data D;
    if (!scan::Load(D, scan::kCompact)) return;
    const int NR = kNRepCB;
    // N (us) with group variance; like-sign counts (fine and family band-cells) and replicas
    Acc A; A.Init((size_t)kNSelE * gStride);
    std::vector<double> ls[kNSelE][2];            // like-sign fine counts (pp, mm), plain
    std::vector<double> lsAll[kNSelE][2];         // like-sign rows per block cell (all parts; plain counts) for information
    std::vector<double> fam[kNSelE][kNFam];       // band-cell like-sign counts
    std::vector<double> lsRep[kNSelE][2];         // replicas of the fine like-sign counts [i * NR + r]
    std::vector<double> famRep[kNSelE][kNFam];    // replicas of the band-cell counts
    for (int s = 0; s < kNSelE; ++s) {
        for (int q = 0; q < 2; ++q) { ls[s][q].assign(kNfine, 0.); lsAll[s][q].assign(gStride, 0.); lsRep[s][q].assign((size_t)kNfine * NR, 0.); }
        for (int f = 0; f < kNFam; ++f) { fam[s][f].assign(kNBand * kFamN[f], 0.); famRep[s][f].assign((size_t)kNBand * kFamN[f] * NR, 0.); }
    }
    long long nPass[kNSelE][3] = {};
    std::vector<double> kr(NR, 0.);
    int cells[64];
    for (int k = 0; k < 3; ++k) {
        const std::vector<scan::Row>& v = D.s[k];
        int cur = -1;
        bool haveK = false;
        for (size_t i = 0; i < v.size(); ++i) {
            const scan::Row& r = v[i];
            if (r.grp != cur) { if (k == 0) A.Begin(); cur = r.grp; haveK = false; }
            if (!(r.m4 == r.m4)) continue;
            const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, r.dz_pi, r.mm4);
            if (!pre) continue;
            const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4);
            bool pc[kNCut];
            bool any = false;
            for (int c = 0; c < kNCut; ++c) { pc[c] = PassCut(c, r.mpipi, r.oa, r.oa_lab, r.oa_eta); any |= pc[c]; }
            if (!any) continue;
            const int nc = CellsOf(r.m4, r.mee, r.oa, -1., cells);
            const int fi = FineBin(r.m4);
            for (int s = 0; s < kNSelE; ++s) {
                const int b = s / kNCut, c = s % kNCut;
                if (!pc[c] || (b == 1 && !s0)) continue;
                nPass[s][k]++;
                if (k == 0) {
                    const size_t base = (size_t)s * gStride;
                    for (int q = 0; q < nc; ++q) A.Add(base + cells[q], 1.);
                    continue;
                }
                for (int q = 0; q < nc; ++q) lsAll[s][k - 1][cells[q]] += 1.;
                if (fi < 0) continue;
                ls[s][k - 1][fi] += 1.;
                const int band = BandOfFine(fi);
                for (int f = 0; f < kNFam; ++f) {
                    const int cell = FamCell(f, r.mee, r.oa);
                    if (cell >= 0) fam[s][f][band * kFamN[f] + cell] += 1.;
                }
                if (!haveK) { GroupPois(kSeedCB, k, r.grp, NR, kr.data()); haveK = true; }
                double* p = &lsRep[s][k - 1][(size_t)fi * NR];
                for (int q = 0; q < NR; ++q) p[q] += kr[q];
                for (int f = 0; f < kNFam; ++f) {
                    const int cell = FamCell(f, r.mee, r.oa);
                    if (cell < 0) continue;
                    double* pf = &famRep[s][f][((size_t)band * kFamN[f] + cell) * NR];
                    for (int q = 0; q < NR; ++q) pf[q] += kr[q];
                }
            }
        }
        if (k == 0) A.Flush();
        printf("e1: data sample %s done (%.1f s)\n", scan::SampleName(k), sw.RealTime());
        sw.Continue();
    }

    // ---- old-study reproduction: scan::FillMany with the FillOpt of sim_rw_compare.C (EnrichedSpectra)
    std::vector<scan::Sel> osel;
    for (int j = 0; j < kNOld; ++j) {
        const int c = kOldCut[j];
        osel.push_back([c](const scan::Row& r) { return simeng::PassPre(r) && PassCut(c, r.mpipi, r.oa, r.oa_lab, r.oa_eta); });
    }
    scan::FillOpt fo; fo.nbin = 260; fo.lo = 0.35; fo.hi = 1.00; fo.reb = 4;
    for (int k = 0; k < kNOldReg; ++k) fo.windows.push_back(std::make_pair(kOldRegLo[k], kOldRegHi[k]));
    const std::vector<scan::Spectra> osp = scan::FillMany(D, osel, fo);
    D.s[0].clear(); D.s[0].shrink_to_fit(); D.s[1].clear(); D.s[1].shrink_to_fit(); D.s[2].clear(); D.s[2].shrink_to_fit();

    // ---- CB per selection: nominal, linear variance, replica mean / rms
    JW j; j.Open(std::string(kOutDir) + "/e1_data.json");
    fprintf(j.f, "{\"meta\": {\"what\": \"e1 data, raw rows; CB k = 1\", \"L\": %.4f, \"PT3\": %.10g, \"nrep\": %d, \"stride\": %d", kL, kPT3, NR, gStride);
    fprintf(j.f, ", \"parts\": {");
    for (int p = 0; p < kNPart; ++p) fprintf(j.f, "%s\"%s\": [%d, %d]", p ? ", " : "", kPartTag[p], gOff[p], kPartN[p]);
    fprintf(j.f, "}, \"windows\": [");
    for (int w = 0; w < kNW; ++w) fprintf(j.f, "%s[\"%s\", %.6f, %.6f, \"%s\"]", w ? ", " : "", kW[w].tag, kW[w].lo, kW[w].hi, kW[w].text);
    fprintf(j.f, "], \"cuts\": [");
    for (int c = 0; c < kNCut; ++c) fprintf(j.f, "%s[\"%s\", \"%s\"]", c ? ", " : "", kCutTag[c], kCutText[c]);
    fprintf(j.f, "]}, \"sel\": {");
    TFile* fr = TFile::Open((std::string(kOutDir) + "/e1_data.root").c_str(), "RECREATE");
    for (int s = 0; s < kNSelE; ++s) {
        // nominal CB
        CB1D cb; cb.Make(ls[s][0].data(), ls[s][1].data(), kCBreb);
        std::vector<std::vector<double>> F(kNFam);
        for (int f = 0; f < kNFam; ++f) F[f] = Fractions(fam[s][f].data(), kFamN[f]);
        // CB of every block cell: a vector of (fine index, weight) per cell, for the families with f
        std::vector<double> C(gStride, 0.), V(gStride, 0.);
        auto cellCB = [&](const CB1D& cbx, const std::vector<std::vector<double>>& Fx, std::vector<double>& out, std::vector<double>* var) {
            out.assign(gStride, 0.);
            if (var) var->assign(gStride, 0.);
            auto put = [&](int cell, const std::vector<std::pair<int, double>>& u) {
                const auto p = cbx.Sum(u);
                out[cell] = p.first;
                if (var) (*var)[cell] = p.second;
            };
            for (int i = 0; i < kNfine; ++i) put(gOff[pFine] + i + 1, {{i, 1.}});
            for (int x = 1; x <= kNb5; ++x) put(gOff[pB5] + x, {{2 * (x - 1), 1.}, {2 * (x - 1) + 1, 1.}});
            for (int w = 0; w < kNW; ++w) {
                std::vector<std::pair<int, double>> u;
                for (int i = 0; i < kNfine; ++i) { const double fr2 = FracFine(i, kW[w].lo, kW[w].hi); if (fr2 > 0) u.push_back({i, fr2}); }
                put(gOff[pWin] + w, u);
            }
            auto famPut = [&](int part, int fam, int ncell, int i0, int i1, int famCellOffset) {
                for (int cc = 0; cc < ncell; ++cc) {
                    std::vector<std::pair<int, double>> u;
                    for (int i = i0; i <= i1; ++i) {
                        const int band = BandOfFine(i);
                        const double fr2 = Fx[fam][band * kFamN[fam] + cc + famCellOffset];
                        if (fr2 > 0) u.push_back({i, fr2});
                    }
                    put(gOff[part] + cc, u);
                }
            };
            famPut(pX17, fX17, 1, kPBi0, kPBi1, 0);
            famPut(pMee, fMee, kNmee, kPBi0, kPBi1, 0);
            famPut(pOa, fOa, kNoa, kPBi0, kPBi1, 0);
            famPut(pOaX, fOaX, kNoa, kPBi0, kPBi1, 0);
            // x2dm: X2D M4 range (fine 89..123) x M(ee) 2-62
            famPut(pX2Dm, fX2Dm, 1, 89, 123, 0);
            // meeW: M(ee) in 0.40-0.70 (fine 40..159)
            famPut(pMeeW, fMee, kNmee, 40, 159, 0);
        };
        cellCB(cb, F, C, &V);
        // replicas
        std::vector<double> m1(gStride, 0.), m2(gStride, 0.), pp(kNfine), mm(kNfine), cnt, Cr;
        for (int r = 0; r < NR; ++r) {
            for (int i = 0; i < kNfine; ++i) { pp[i] = lsRep[s][0][(size_t)i * NR + r]; mm[i] = lsRep[s][1][(size_t)i * NR + r]; }
            CB1D cbr; cbr.Make(pp.data(), mm.data(), kCBreb);
            std::vector<std::vector<double>> Fr(kNFam);
            for (int f = 0; f < kNFam; ++f) {
                cnt.assign(kNBand * kFamN[f], 0.);
                for (int b = 0; b < kNBand * kFamN[f]; ++b) cnt[b] = famRep[s][f][(size_t)b * NR + r];
                Fr[f] = Fractions(cnt.data(), kFamN[f]);
            }
            cellCB(cbr, Fr, Cr, nullptr);
            for (int q = 0; q < gStride; ++q) { m1[q] += Cr[q]; m2[q] += Cr[q] * Cr[q]; }
        }
        std::vector<double> rms(gStride, 0.), mean(gStride, 0.);
        for (int q = 0; q < gStride; ++q) { mean[q] = m1[q] / NR; rms[q] = std::sqrt(std::max(0., m2[q] / NR - mean[q] * mean[q])); }
        std::vector<double> lin(gStride, 0.);
        for (int q = 0; q < gStride; ++q) lin[q] = std::sqrt(std::max(0., V[q]));
        const size_t base = (size_t)s * gStride;
        std::vector<double> N(A.s.begin() + base, A.s.begin() + base + gStride), NV(A.v.begin() + base, A.v.begin() + base + gStride);
        fprintf(j.f, "%s\"%s\": {\"npass\": [%lld, %lld, %lld]", s ? ", " : "", SelName(s).c_str(), nPass[s][0], nPass[s][1], nPass[s][2]);
        j.Arr("N", N.data(), gStride); j.Arr("NV", NV.data(), gStride);
        j.Arr("Npp", lsAll[s][0].data(), gStride); j.Arr("Nmm", lsAll[s][1].data(), gStride);
        j.Arr("CB", C.data(), gStride); j.Arr("CBlin", lin.data(), gStride); j.Arr("CBrepmean", mean.data(), gStride); j.Arr("CBrep", rms.data(), gStride);
        fprintf(j.f, "}");
        // ROOT copies (fine, b5)
        TDirectory* d = fr->mkdir(SelName(s).c_str());
        for (int p : {pFine, pB5}) {
            const int nb = p == pFine ? kNfine : kNb5;
            for (int what = 0; what < 2; ++what) {
                TH1D* h = new TH1D(Form("%s_%s", what ? "CB" : "N", kPartTag[p]), "", nb, kFineLo, kFineHi);
                h->Sumw2();
                for (int b = 0; b <= nb + 1; ++b) {
                    const int q = gOff[p] + b;
                    h->SetBinContent(b, what ? C[q] : N[q]);
                    h->SetBinError(b, what ? rms[q] : std::sqrt(NV[q]));
                }
                d->WriteTObject(h); delete h;
            }
        }
        printf("e1: data %-10s us %9lld pp %7lld mm %7lld  P: N %8.0f CB %8.1f (lin %.1f, rep %.1f)  X17: N %6.0f CB %6.1f (rep %.1f)\n",
               SelName(s).c_str(), nPass[s][0], nPass[s][1], nPass[s][2], N[gOff[pWin] + 4], C[gOff[pWin] + 4], lin[gOff[pWin] + 4],
               rms[gOff[pWin] + 4], N[gOff[pX17]], C[gOff[pX17]], rms[gOff[pX17]]);
    }
    fprintf(j.f, "}, \"old\": {");
    for (int o = 0; o < kNOld; ++o) {
        fprintf(j.f, "%s\"%s\": {", o ? ", " : "", kOldTag[o]);
        for (int k = 0; k < kNOldReg; ++k) {
            const scan::WinSum& w = osp[o].win[2 + k];
            fprintf(j.f, "%s\"%s\": [%.10g, %.10g, %.10g, %.10g]", k ? ", " : "", kOldRegTag[k], w.n[0], w.cb, w.sig(), w.sigErr());
        }
        fprintf(j.f, "}");
    }
    fprintf(j.f, "}}\n");
    fclose(j.f);
    fr->Close(); delete fr;
    printf("e1: data written (%.1f s)\n", sw.RealTime());
}

// ===========================================================================
// simulation
// ===========================================================================
const int kNC = 14;
const int kCompOf[kNC] = {kC2, kC3, kC4, kC5, kC6, kGeta, kGpi0D, kGconv, kGvec, kGoth, kGD, kT12, kT3, kT4};
inline int LocalComp(int c) { for (int q = 0; q < kNC; ++q) if (kCompOf[q] == c) return q; return -1; }

void RunSim() {
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    Acc A; A.Init((size_t)kNSelE * kNC * gStride);
    // old-study reproduction (T1 + T2 = CompSim, raw M4, raw M(pipi)): [old sel][region]
    Acc O; O.Init((size_t)kNOld * kNOldReg);
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    int cells[64];
    long long nPass[kNSelE] = {};
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& v = D->t[trees[t]];
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const simeng::Row& r = v[i];
            if (r.grp != cur) { A.Begin(); O.Begin(); cur = r.grp; }
            const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, 0.f);
            if (!pre) continue;
            const bool oA = r.ownA == 1;
            if (!oA) continue;   // mode A: rows not owned in A carry weight 0 in every component
            const double f108 = simphys::RowFactor(r, kPhysNom);
            const double w = r.w * f108;
            // old-study sums
            {
                const bool cs = (r.tmpl == 1 && !simeng::X17Row(r)) || r.tmpl == 2;
                if (cs && w != 0.)
                    for (int o = 0; o < kNOld; ++o) {
                        if (!PassCut(kOldCut[o], r.mpipi, r.oa, r.oa_lab, r.oa_eta)) continue;
                        for (int k = 0; k < kNOldReg; ++k)
                            if (r.m4 >= kOldRegLo[k] && r.m4 < kOldRegHi[k]) O.Add((size_t)o * kNOldReg + k, w);
                    }
            }
            const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, 0.f);
            const float mpp = (float)(kSPiPi * r.mpipi);
            bool pc[kNCut];
            bool any = false;
            for (int c = 0; c < kNCut; ++c) { pc[c] = PassCut(c, mpp, r.oa, r.oa_lab, r.oa_eta); any |= pc[c]; }
            if (!any) continue;
            int cn[8];
            const int ncn = CompsOf(r, false, cn);
            if (ncn == 0) continue;
            const double m4p = M4sim(r.m4, r.m4_true);
            const int nc = CellsOf(m4p, r.mee, r.oa, r.p_ee, cells);
            for (int s = 0; s < kNSelE; ++s) {
                const int b = s / kNCut, c = s % kNCut;
                if (!pc[c] || (b == 1 && !s0)) continue;
                nPass[s]++;
                for (int q = 0; q < ncn; ++q) {
                    const int lc = LocalComp(cn[q]);
                    if (lc < 0) continue;
                    const double wc = cn[q] == kC6 ? 2. : w;
                    if (wc == 0.) continue;
                    const size_t base = ((size_t)s * kNC + lc) * gStride;
                    for (int k = 0; k < nc; ++k) A.Add(base + cells[k], wc);
                }
            }
        }
        A.Flush(); O.Flush();
        printf("e1: sim tree %s: %zu rows (%.1f s)\n", simeng::TreeName(trees[t]), v.size(), sw.RealTime());
        sw.Continue();
    }
    printf("e1: angles hit %lld miss %lld\n", simphys::AngHit(), simphys::AngMiss());
    delete D;
    JW j; j.Open(std::string(kOutDir) + "/e1_sim.json");
    fprintf(j.f, "{\"meta\": {\"what\": \"e1 simulation, nb per cell = sum w / L_sim; C6 per unit BR\", \"Lsim\": %.4f, \"stride\": %d, \"comps\": [", kLsim, gStride);
    for (int q = 0; q < kNC; ++q) fprintf(j.f, "%s\"%s\"", q ? ", " : "", CompTag(kCompOf[q]));
    fprintf(j.f, "]}, \"sel\": {");
    TFile* fr = TFile::Open((std::string(kOutDir) + "/e1_sim.root").c_str(), "RECREATE");
    const double sc = 1. / kLsim;
    for (int s = 0; s < kNSelE; ++s) {
        fprintf(j.f, "%s\"%s\": {\"npass\": %lld", s ? ", " : "", SelName(s).c_str(), nPass[s]);
        TDirectory* d = fr->mkdir(SelName(s).c_str());
        for (int q = 0; q < kNC; ++q) {
            const size_t base = ((size_t)s * kNC + q) * gStride;
            std::vector<double> c(gStride), e(gStride);
            for (int k = 0; k < gStride; ++k) { c[k] = A.s[base + k] * sc; e[k] = A.v[base + k] * sc * sc; }
            fprintf(j.f, ", \"%s\": {", CompTag(kCompOf[q]));
            j.Arr("s", c.data(), gStride, false); j.Arr("v", e.data(), gStride);
            fprintf(j.f, "}");
            TH1D* h = new TH1D(Form("%s_b5", CompTag(kCompOf[q])), "nb per bin", kNb5, kFineLo, kFineHi);
            h->Sumw2();
            for (int b = 0; b <= kNb5 + 1; ++b) { h->SetBinContent(b, c[gOff[pB5] + b]); h->SetBinError(b, std::sqrt(e[gOff[pB5] + b])); }
            d->WriteTObject(h); delete h;
        }
        fprintf(j.f, "}");
        const size_t bT = ((size_t)s * kNC + LocalComp(kT12)) * gStride, bX = ((size_t)s * kNC + LocalComp(kC6)) * gStride;
        printf("e1: sim %-10s rows %9lld  T12 P %.4f nb  C6 X17win %.2f nb  C6 X2Dm %.2f nb\n", SelName(s).c_str(), nPass[s],
               A.s[bT + gOff[pWin] + 4] * sc, A.s[bX + gOff[pX17]] * sc, A.s[bX + gOff[pX2Dm]] * sc);
    }
    fprintf(j.f, "}, \"old\": {");
    for (int o = 0; o < kNOld; ++o) {
        fprintf(j.f, "%s\"%s\": {", o ? ", " : "", kOldTag[o]);
        for (int k = 0; k < kNOldReg; ++k)
            fprintf(j.f, "%s\"%s\": [%.10g, %.10g]", k ? ", " : "", kOldRegTag[k], O.s[(size_t)o * kNOldReg + k] * sc,
                    std::sqrt(O.v[(size_t)o * kNOldReg + k]) * sc);
        fprintf(j.f, "}");
    }
    fprintf(j.f, "}}\n");
    fclose(j.f);
    fr->Close(); delete fr;
    printf("e1: sim written (%.1f s)\n", sw.RealTime());
}

}  // namespace e1

void e1_fill(int mode = 0) {
    e1::InitLayout();
    gSystem->mkdir(e1::kOutDir, kTRUE);
    printf("e1: stride %d cells\n", e1::gStride);
    if (mode & 1) e1::RunData();
    if (mode & 2) e1::RunSim();
}
