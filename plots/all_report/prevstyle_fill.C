// ============================================================================
// prevstyle_fill.C  --  package prevstyle of report_all: data and SMASH after the "typical" eta cuts (E1) in the
//                       style of the earlier analysis (HADES collaboration-meeting slides), histogram factory
// ============================================================================
// Written 2026-11-16 (prefix prevstyle). New file; nothing else is changed. Conventions of e1_fill.C (package e1),
// followed line by line, with two selections instead of 22:
//   sel 0 "best": pre + E1 + isBest == 1 (best combination only, as in the earlier analysis; data AND simulation)
//   sel 1 "all" : pre + E1, isBest >= 0 (our nominal; = e1 "pre_E1")
//   E1 = oa_lab < 50 deg, oa_eta > 140 deg, M(pipi) < 0.46 GeV/c2, OA(e+e-) > 2 deg
//   pre = isBest >= 0, vz > -500, PT3, vz > -160 mm, cor_ok == 1 (alldefs::PassSel(kSelPre))
// data : rows of scan_v2 (us, pp, mm); N with event-group variance; CB = 2 sqrt(N++ N--) of the like-sign rows with the
//        same cuts (for sel 0 the like-sign rows with isBest == 1, i.e. the naive like-sign formula of the earlier analysis),
//        10 MeV M4 cells interpolated to the 2.5 MeV fine bins (wp0 CB1D), k = 1; windows with fractional fine bins;
//        CB of an M(ee) cell in an M4 window = sum over its fine M4 bins of frac x CB(fine) x f(cell | 50 MeV band) (wp0
//        X2D rule); CB errors = RMS of 200 Poisson(1) event-group replicas of the like-sign rows (wp0 generator and seed).
//        Units: RAW COUNTS (no 64/63, no L).
// sim  : SMASH std + lep + eta e+e- rows, owned, mode A, w x simphys::RowFactor(r, 108); C2 rows also x their QED gamma*
//        weight q_A of package c2recon (c2recon_weights.txt, joined by (tree, index) of the same loader call; "C2iso" keeps
//        the isotropic C2 for comparison); M4' = 0.99615 M4; M(pipi)' = 0.9952 M(pipi) in the M(pipi) cut; components C2,
//        C2iso, C3, C4, C5, C6 (X17 per unit BR, weight 2 per owned row), G_eta, G_D = G_pi0D + G_vec + G_oth, G_C = G_conv,
//        T12 (C2 with q_A); units nb = sum w / 221.05 nb^-1, event-group variance.
// cells of one block: fine (322: flows + 320 x 2.5 MeV, 0.30-1.10), b5 (162: 5 MeV), windows (kNW), and M(ee) in 4 M4 windows
//        (SBlo 0.518-0.530, Wpk 0.530-0.557, SBhi 0.557-0.569, PB 0.5325-0.5575; 2 MeV bins 0-62 MeV + overflow = 32 each),
//        and the two X17 M(ee) ranges (15.5-24.5 MeV, ours; 12-22 MeV, the earlier analysis') in the same 4 M4 windows (8 cells;
//        CB with the like-sign share of the range in the 50 MeV band).
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_<data|sim>.json (+ .log)
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/aclic", kTRUE)' 'plots/all_report/prevstyle_fill.C+O(0)'  (compile)
//   ... 'plots/all_report/prevstyle_fill.C+O(1)'   data (about 40 s)
//   ... 'plots/all_report/prevstyle_fill.C+O(2)'   simulation (about 2 min, about 8 GB); the two can run in parallel
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TStopwatch.h"
#include "TSystem.h"

#include "defs.h"
#include "sim_eta_physics.h"

namespace prevst {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/prevstyle";
const char* const kC2Dir = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon";
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) { fprintf(gLog, __VA_ARGS__); fflush(gLog); } fflush(stdout); } while (0)

// ---------------------------------------------------------------- random numbers (copy of wp0_fill.C / e1_fill.C)
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
const int kNSelP = 2;
const char* const kSelName[kNSelP] = {"best", "all"};
const char* const kSelText[kNSelP] = {"pre + E1 + isBest == 1 (best combination only)", "pre + E1, isBest >= 0 (all stored combinations)"};
inline bool PassE1(float mpp, float oa, float oalab, float oaeta) {
    return oalab < 50.f && oaeta > 140.f && mpp < 0.46f && oa > 2.f;
}

// ---------------------------------------------------------------- windows (exact M4 for N and sim; CB with fractional fine bins)
struct Win { const char* tag; double lo, hi; const char* text; };
const int kNW = 13;
const Win kW[kNW] = {
    {"all", 0.35, 0.90, "0.35-0.90 (integral normalisation of the earlier style)"},
    {"fit1", 0.40, 0.70, "0.40-0.70 (chi2 range of the M4 panels)"},
    {"gfit", 0.43, 0.75, "0.43-0.75 (three-scale global fit range)"},
    {"SBlo", 0.518, 0.530, "lower side band 518-530 MeV/c2"},
    {"Wpk", 0.530, 0.557, "eta peak window 530-557 MeV/c2"},
    {"SBhi", 0.557, 0.569, "upper side band 557-569 MeV/c2"},
    {"PB", 0.5325, 0.5575, "e1 eta peak band (check)"},
    {"P", 0.53364, 0.55416, "peak window P (check)"},
    {"CR", 0.575, 0.70, "0.575-0.70 (c_b normalisation)"},
    {"B", 0.4550, 0.5225, "bump window B"},
    {"r56", 0.56, 0.70, "0.56-0.70"},
    {"r70", 0.70, 0.90, "0.70-0.90"},
    {"below", 0.35, 0.455, "below the bump 0.35-0.455"}};
const int kNMW = 4;
const int kMeeWin[kNMW] = {3, 4, 5, 6};   // SBlo, Wpk, SBhi, PB

// ---------------------------------------------------------------- cell layout of one block
const int kNfineC = kNfine + 2, kNb5 = 160, kNb5C = kNb5 + 2, kNmee = 32;
enum Part { pFine = 0, pB5, pWin, pMee, pXr, kNPart };
const int kNXr = 2;
const double kXrLo[kNXr] = {0.0155, 0.012}, kXrHi[kNXr] = {0.0245, 0.022};
const int kPartN[kNPart] = {kNfineC, kNb5C, kNW, kNMW * kNmee, kNMW * kNXr};
const char* const kPartTag[kNPart] = {"fine", "b5", "win", "mee", "xr"};
inline bool InXr(int x, double mee) { return mee >= kXrLo[x] && mee < kXrHi[x]; }
int gOff[kNPart], gStride = 0;
void InitLayout() {
    gStride = 0;
    for (int p = 0; p < kNPart; ++p) { gOff[p] = gStride; gStride += kPartN[p]; }
}
inline int MeeBin(double mee) { if (!(mee == mee) || mee < 0) return -1; const int j = (int)(mee / 0.002); return std::min(j, kNmee - 1); }
inline int B5Bin(int fi, double m4) { return fi >= 0 ? fi / 2 + 1 : (m4 < kFineLo ? 0 : kNb5 + 1); }
inline int FineCell(int fi, double m4) { return fi >= 0 ? fi + 1 : (m4 < kFineLo ? 0 : kNfine + 1); }
inline int CellsOf(double m4, double mee, int* out) {
    int n = 0;
    const int fi = FineBin(m4);
    out[n++] = gOff[pFine] + FineCell(fi, m4);
    out[n++] = gOff[pB5] + B5Bin(fi, m4);
    for (int w = 0; w < kNW; ++w) if (m4 >= kW[w].lo && m4 < kW[w].hi) out[n++] = gOff[pWin] + w;
    const int jm = MeeBin(mee);
    if (jm >= 0)
        for (int q = 0; q < kNMW; ++q) {
            const Win& x = kW[kMeeWin[q]];
            if (m4 >= x.lo && m4 < x.hi) {
                out[n++] = gOff[pMee] + q * kNmee + jm;
                for (int xr = 0; xr < kNXr; ++xr) if (InXr(xr, mee)) out[n++] = gOff[pXr] + q * kNXr + xr;
            }
        }
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

// ---------------------------------------------------------------- CB (copy of wp0 CB1D / e1_fill.C)
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
// like-sign share of an M(ee) cell in its 50 MeV band (wp0 X2D rule; empty band -> pooled share)
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
    void Open(const std::string& n) { f = fopen(n.c_str(), "w"); }
    void Arr(const char* key, const double* a, int n, bool comma = true) {
        fprintf(f, "%s\"%s\": [", comma ? ", " : "", key);
        for (int i = 0; i < n; ++i) fprintf(f, i ? ",%.10g" : "%.10g", std::isfinite(a[i]) ? a[i] : 0.);
        fprintf(f, "]");
    }
};
void WriteMeta(JW& j, const char* what) {
    fprintf(j.f, "{\"meta\": {\"what\": \"%s\", \"L\": %.4f, \"PT3\": %.10g, \"Lsim\": %.4f, \"nrep\": %d, \"stride\": %d", what, kL, kPT3,
            kLsim, kNRepCB, gStride);
    fprintf(j.f, ", \"parts\": {");
    for (int p = 0; p < kNPart; ++p) fprintf(j.f, "%s\"%s\": [%d, %d]", p ? ", " : "", kPartTag[p], gOff[p], kPartN[p]);
    fprintf(j.f, "}, \"windows\": [");
    for (int w = 0; w < kNW; ++w) fprintf(j.f, "%s[\"%s\", %.6f, %.6f, \"%s\"]", w ? ", " : "", kW[w].tag, kW[w].lo, kW[w].hi, kW[w].text);
    fprintf(j.f, "], \"meewin\": [");
    for (int q = 0; q < kNMW; ++q) fprintf(j.f, "%s\"%s\"", q ? ", " : "", kW[kMeeWin[q]].tag);
    fprintf(j.f, "], \"sels\": [");
    for (int s = 0; s < kNSelP; ++s) fprintf(j.f, "%s[\"%s\", \"%s\"]", s ? ", " : "", kSelName[s], kSelText[s]);
    fprintf(j.f, "]}");
}

// ===========================================================================
// data
// ===========================================================================
void RunData() {
    TStopwatch sw;
    scan::Data D;
    if (!scan::Load(D, scan::kCompact)) return;
    const int NR = kNRepCB;
    Acc A; A.Init((size_t)kNSelP * gStride);
    std::vector<double> ls[kNSelP][2], lsAll[kNSelP][2], fam[kNSelP], lsRep[kNSelP][2], famRep[kNSelP];
    std::vector<double> famX[kNSelP][kNXr], famXRep[kNSelP][kNXr];   // X17 ranges: 2 cells (in, out) per band
    for (int s = 0; s < kNSelP; ++s) {
        for (int q = 0; q < 2; ++q) { ls[s][q].assign(kNfine, 0.); lsAll[s][q].assign(gStride, 0.); lsRep[s][q].assign((size_t)kNfine * NR, 0.); }
        fam[s].assign(kNBand * kNmee, 0.); famRep[s].assign((size_t)kNBand * kNmee * NR, 0.);
        for (int x = 0; x < kNXr; ++x) { famX[s][x].assign(kNBand * 2, 0.); famXRep[s][x].assign((size_t)kNBand * 2 * NR, 0.); }
    }
    long long nPass[kNSelP][3] = {};
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
            if (!PassSel(kSelPre, r, r.cor_ok == 1, r.dz_pi, r.mm4)) continue;
            if (!PassE1(r.mpipi, r.oa, r.oa_lab, r.oa_eta)) continue;
            const int nc = CellsOf(r.m4, r.mee, cells);
            const int fi = FineBin(r.m4);
            for (int s = 0; s < kNSelP; ++s) {
                if (s == 0 && r.isBest != 1) continue;
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
                const int jm = MeeBin(r.mee);
                if (jm >= 0) fam[s][band * kNmee + jm] += 1.;
                if (jm >= 0) for (int x = 0; x < kNXr; ++x) famX[s][x][band * 2 + (InXr(x, r.mee) ? 0 : 1)] += 1.;
                if (!haveK) { GroupPois(kSeedCB, k, r.grp, NR, kr.data()); haveK = true; }
                double* p = &lsRep[s][k - 1][(size_t)fi * NR];
                for (int q = 0; q < NR; ++q) p[q] += kr[q];
                if (jm >= 0) {
                    double* pf = &famRep[s][((size_t)band * kNmee + jm) * NR];
                    for (int q = 0; q < NR; ++q) pf[q] += kr[q];
                    for (int x = 0; x < kNXr; ++x) {
                        double* px = &famXRep[s][x][((size_t)band * 2 + (InXr(x, r.mee) ? 0 : 1)) * NR];
                        for (int q = 0; q < NR; ++q) px[q] += kr[q];
                    }
                }
            }
        }
        if (k == 0) A.Flush();
        PL("prevstyle: data sample %s done (%.1f s)\n", scan::SampleName(k), sw.RealTime());
        sw.Continue();
    }
    D.s[0].clear(); D.s[0].shrink_to_fit(); D.s[1].clear(); D.s[1].shrink_to_fit(); D.s[2].clear(); D.s[2].shrink_to_fit();

    JW j; j.Open(std::string(kOutDir) + "/prevstyle_data.json");
    WriteMeta(j, "prevstyle data, raw counts; CB k = 1");
    fprintf(j.f, ", \"sel\": {");
    auto cellCB = [&](const CB1D& cbx, const std::vector<double>& F, const std::vector<double>* FX, std::vector<double>& out,
                      std::vector<double>* var) {
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
        for (int q = 0; q < kNMW; ++q) {
            const Win& x = kW[kMeeWin[q]];
            for (int c = 0; c < kNmee; ++c) {
                std::vector<std::pair<int, double>> u;
                for (int i = 0; i < kNfine; ++i) {
                    const double fr2 = FracFine(i, x.lo, x.hi);
                    if (!(fr2 > 0)) continue;
                    const double f = F[BandOfFine(i) * kNmee + c];
                    if (f > 0) u.push_back({i, fr2 * f});
                }
                put(gOff[pMee] + q * kNmee + c, u);
            }
            for (int xr = 0; xr < kNXr; ++xr) {
                std::vector<std::pair<int, double>> u;
                for (int i = 0; i < kNfine; ++i) {
                    const double fr2 = FracFine(i, x.lo, x.hi);
                    if (!(fr2 > 0)) continue;
                    const double f = FX[xr][BandOfFine(i) * 2 + 0];
                    if (f > 0) u.push_back({i, fr2 * f});
                }
                put(gOff[pXr] + q * kNXr + xr, u);
            }
        }
    };
    for (int s = 0; s < kNSelP; ++s) {
        CB1D cb; cb.Make(ls[s][0].data(), ls[s][1].data(), kCBreb);
        const std::vector<double> F = Fractions(fam[s].data(), kNmee);
        std::vector<double> FX[kNXr];
        for (int x = 0; x < kNXr; ++x) FX[x] = Fractions(famX[s][x].data(), 2);
        std::vector<double> C, V;
        cellCB(cb, F, FX, C, &V);
        std::vector<double> m1(gStride, 0.), m2(gStride, 0.), pp(kNfine), mm(kNfine), cnt(kNBand * kNmee), Cr;
        for (int r = 0; r < NR; ++r) {
            for (int i = 0; i < kNfine; ++i) { pp[i] = lsRep[s][0][(size_t)i * NR + r]; mm[i] = lsRep[s][1][(size_t)i * NR + r]; }
            CB1D cbr; cbr.Make(pp.data(), mm.data(), kCBreb);
            for (int b = 0; b < kNBand * kNmee; ++b) cnt[b] = famRep[s][(size_t)b * NR + r];
            const std::vector<double> Fr = Fractions(cnt.data(), kNmee);
            std::vector<double> FXr[kNXr], cx(kNBand * 2);
            for (int x = 0; x < kNXr; ++x) {
                for (int b = 0; b < kNBand * 2; ++b) cx[b] = famXRep[s][x][(size_t)b * NR + r];
                FXr[x] = Fractions(cx.data(), 2);
            }
            cellCB(cbr, Fr, FXr, Cr, nullptr);
            for (int q = 0; q < gStride; ++q) { m1[q] += Cr[q]; m2[q] += Cr[q] * Cr[q]; }
        }
        std::vector<double> rms(gStride, 0.), mean(gStride, 0.), lin(gStride, 0.);
        for (int q = 0; q < gStride; ++q) {
            mean[q] = m1[q] / NR; rms[q] = std::sqrt(std::max(0., m2[q] / NR - mean[q] * mean[q]));
            lin[q] = std::sqrt(std::max(0., V[q]));
        }
        const size_t base = (size_t)s * gStride;
        std::vector<double> N(A.s.begin() + base, A.s.begin() + base + gStride), NV(A.v.begin() + base, A.v.begin() + base + gStride);
        fprintf(j.f, "%s\"%s\": {\"npass\": [%lld, %lld, %lld]", s ? ", " : "", kSelName[s], nPass[s][0], nPass[s][1], nPass[s][2]);
        j.Arr("N", N.data(), gStride); j.Arr("NV", NV.data(), gStride);
        j.Arr("Npp", lsAll[s][0].data(), gStride); j.Arr("Nmm", lsAll[s][1].data(), gStride);
        j.Arr("CB", C.data(), gStride); j.Arr("CBlin", lin.data(), gStride); j.Arr("CBrepmean", mean.data(), gStride); j.Arr("CBrep", rms.data(), gStride);
        fprintf(j.f, "}");
        const int wP = gOff[pWin] + 7, wW = gOff[pWin] + 4;
        PL("prevstyle: data %-4s us %9lld pp %7lld mm %7lld | P: N %8.0f CB %8.1f (rep %.1f) | Wpk: N %8.0f CB %8.1f (rep %.1f)\n", kSelName[s],
           nPass[s][0], nPass[s][1], nPass[s][2], N[wP], C[wP], rms[wP], N[wW], C[wW], rms[wW]);
    }
    fprintf(j.f, "}}\n");
    fclose(j.f);
    PL("prevstyle: data written (%.1f s)\n", sw.RealTime());
}

// ===========================================================================
// simulation
// ===========================================================================
enum LComp { lC2 = 0, lC2iso, lC3, lC4, lC5, lC6, lGeta, lGD, lGC, lT12, kNL };
const char* const kLTag[kNL] = {"C2", "C2iso", "C3", "C4", "C5", "C6", "Geta", "GD", "GC", "T12"};

// QED gamma* weights of the C2 rows (c2recon): per tree, per row index; NaN = no weight
bool ReadQED(const simeng::SimData& D, std::vector<double> q[], std::vector<float> m4ref[], long& nRows) {
    for (int t = 0; t < simeng::kNTreeAll; ++t) { q[t].assign(D.t[t].size(), std::numeric_limits<double>::quiet_NaN()); m4ref[t].assign(D.t[t].size(), -1.f); }
    std::ifstream in(Form("%s/c2recon_rows.txt", kC2Dir));
    std::ifstream iw(Form("%s/c2recon_weights.txt", kC2Dir));
    if (!in || !iw) { PL("prevstyle: cannot open the c2recon files\n"); return false; }
    std::string line;
    nRows = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        int t; long ix; double x[11];
        s >> t >> ix;
        for (int k = 0; k < 11; ++k) s >> x[k];   // set pat_entry grp s0 best cor_ok W108 W0 W116 W100 m4
        if (t < 0 || t >= simeng::kNTreeAll || ix < 0 || ix >= (long)D.t[t].size()) { PL("prevstyle: bad c2recon row %d %ld\n", t, ix); return false; }
        m4ref[t][ix] = (float)x[10];
        ++nRows;
    }
    long nw = 0;
    while (std::getline(iw, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        int t; long ix; double qa;
        s >> t >> ix >> qa;
        if (t < 0 || t >= simeng::kNTreeAll || ix < 0 || ix >= (long)D.t[t].size()) { PL("prevstyle: bad c2recon weight %d %ld\n", t, ix); return false; }
        q[t][ix] = qa;
        ++nw;
    }
    PL("prevstyle: c2recon rows %ld, weights %ld\n", nRows, nw);
    return nw == nRows;
}

void RunSim() {
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    std::vector<double> qed[simeng::kNTreeAll];
    std::vector<float> m4ref[simeng::kNTreeAll];
    long nC2ref = 0;
    if (!ReadQED(*D, qed, m4ref, nC2ref)) return;
    Acc A; A.Init((size_t)kNSelP * kNL * gStride);
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    int cells[64];
    long long nPass[kNSelP] = {};
    long nC2pre = 0, nC2miss = 0, nC2m4bad = 0;
    double sumQ = 0, sumW = 0;
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& v = D->t[trees[t]];
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const simeng::Row& r = v[i];
            if (r.grp != cur) { A.Begin(); cur = r.grp; }
            if (!PassSel(kSelPre, r, r.cor_ok == 1, 0.f, 0.f)) continue;
            if (r.ownA != 1) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom);
            int cn[8];
            const int ncn = CompsOf(r, false, cn);
            if (ncn == 0) continue;
            bool isC2 = false;
            for (int q = 0; q < ncn; ++q) if (cn[q] == kC2) isC2 = true;
            double qa = 1.;
            if (isC2) {
                ++nC2pre;
                qa = qed[trees[t]][i];
                if (!(qa == qa)) { ++nC2miss; qa = 1.; }
                else if (std::fabs(m4ref[trees[t]][i] - r.m4) > 1e-6f) ++nC2m4bad;
                sumQ += w * qa; sumW += w;
            }
            const float mpp = (float)(kSPiPi * r.mpipi);
            if (!PassE1(mpp, r.oa, r.oa_lab, r.oa_eta)) continue;
            const double m4p = M4sim(r.m4, r.m4_true);
            const int nc = CellsOf(m4p, r.mee, cells);
            for (int s = 0; s < kNSelP; ++s) {
                if (s == 0 && r.isBest != 1) continue;
                nPass[s]++;
                for (int q = 0; q < ncn; ++q) {
                    int lc[2] = {-1, -1};
                    double wc[2] = {w, w};
                    switch (cn[q]) {
                        case kC2: lc[0] = lC2; wc[0] = w * qa; lc[1] = lC2iso; break;
                        case kC3: lc[0] = lC3; break;
                        case kC4: lc[0] = lC4; break;
                        case kC5: lc[0] = lC5; break;
                        case kC6: lc[0] = lC6; wc[0] = 2.; break;
                        case kGeta: lc[0] = lGeta; break;
                        case kGD: lc[0] = lGD; break;
                        case kGconv: lc[0] = lGC; break;
                        case kT12: lc[0] = lT12; wc[0] = isC2 ? w * qa : w; break;
                        default: break;
                    }
                    for (int u = 0; u < 2; ++u) {
                        if (lc[u] < 0 || wc[u] == 0.) continue;
                        const size_t base = ((size_t)s * kNL + lc[u]) * gStride;
                        for (int k = 0; k < nc; ++k) A.Add(base + cells[k], wc[u]);
                    }
                }
            }
        }
        A.Flush();
        PL("prevstyle: sim tree %s: %zu rows (%.1f s)\n", simeng::TreeName(trees[t]), v.size(), sw.RealTime());
        sw.Continue();
    }
    PL("prevstyle: C2 rows (pre, ownA): %ld (c2recon list %ld), without a QED weight %ld, m4 mismatch %ld; weighted mean q_A %.5f\n",
       nC2pre, nC2ref, nC2miss, nC2m4bad, sumW > 0 ? sumQ / sumW : 0.);
    if (nC2miss || nC2m4bad || nC2pre != nC2ref) PL("prevstyle: WARNING: the QED join is not exact\n");
    delete D;
    JW j; j.Open(std::string(kOutDir) + "/prevstyle_sim.json");
    WriteMeta(j, "prevstyle simulation, nb per cell = sum w / L_sim; C6 per unit BR; C2 and T12 with the c2recon QED weights");
    fprintf(j.f, ", \"comps\": [");
    for (int q = 0; q < kNL; ++q) fprintf(j.f, "%s\"%s\"", q ? ", " : "", kLTag[q]);
    fprintf(j.f, "], \"qed\": {\"nC2pre\": %ld, \"nC2ref\": %ld, \"miss\": %ld, \"m4bad\": %ld, \"meanq\": %.8g}", nC2pre, nC2ref, nC2miss,
            nC2m4bad, sumW > 0 ? sumQ / sumW : 0.);
    fprintf(j.f, ", \"sel\": {");
    const double sc = 1. / kLsim;
    for (int s = 0; s < kNSelP; ++s) {
        fprintf(j.f, "%s\"%s\": {\"npass\": %lld", s ? ", " : "", kSelName[s], nPass[s]);
        for (int q = 0; q < kNL; ++q) {
            const size_t base = ((size_t)s * kNL + q) * gStride;
            std::vector<double> c(gStride), e(gStride);
            for (int k = 0; k < gStride; ++k) { c[k] = A.s[base + k] * sc; e[k] = A.v[base + k] * sc * sc; }
            fprintf(j.f, ", \"%s\": {", kLTag[q]);
            j.Arr("s", c.data(), gStride, false); j.Arr("v", e.data(), gStride);
            fprintf(j.f, "}");
        }
        fprintf(j.f, "}");
        const size_t bT = ((size_t)s * kNL + lT12) * gStride;
        PL("prevstyle: sim %-4s rows %9lld  T12 P %.4f nb  Wpk %.4f nb\n", kSelName[s], nPass[s], A.s[bT + gOff[pWin] + 7] * sc,
           A.s[bT + gOff[pWin] + 4] * sc);
    }
    fprintf(j.f, "}}\n");
    fclose(j.f);
    PL("prevstyle: sim written (%.1f s)\n", sw.RealTime());
}

}  // namespace prevst

void prevstyle_fill(int mode = 0) {
    prevst::InitLayout();
    gSystem->mkdir(prevst::kOutDir, kTRUE);
    if (mode) prevst::gLog = fopen(Form("%s/prevstyle_fill_%s.log", prevst::kOutDir, mode == 1 ? "data" : (mode == 2 ? "sim" : "both")), "w");
    printf("prevstyle: stride %d cells\n", prevst::gStride);
    if (mode & 1) prevst::RunData();
    if (mode & 2) prevst::RunSim();
    if (prevst::gLog) fclose(prevst::gLog);
}
