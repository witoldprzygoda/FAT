// ============================================================================
// s5_opt_eta.C  --  optimisation of the eta-selection cuts with the simulation truth, and its bias
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SIMULATION ONLY: no measured data are read
// (SimEngine.h includes ../ScanEngine.h for its types, but no data function is called here).
//
// Simulation: SMASH std + lep + eta owned rows, partition mode A, weight w (PDG) x physics weights of the eta
// decays, main variant (simphys::kPhysMain = 108 = V3 + KLOE bump + conversion F_V + background rules).
// Simulated sum = template 1 (signal classes 2-5) + template 2 (physical background); templates 3 (fake
// lepton) and 4 (different sources) enter B only in the stated variant "+T3+T4". X17 rows (class 6) never
// enter anything. sigma = sum w / L_sim (221.05 nb^-1); expected counts at L_data = 6080.5 nb^-1 (UNCONFIRMED,
// a scale only: the location of the optimum does not depend on it).
//
// Base selection "pre" (SimEngine kPre): isBest >= 0, vz > -500, PT3 bit, vz > -160 mm, cor_ok == 1.
//
// Targets (M4 windows of this study, [lo, hi)):
//   PEAK   = classes 4 + 5 in the peak window 0.5333-0.5587 GeV/c2; B = other T1 + T2 rows of the window
//   DIRECT = class 5 only in the peak window (the conversion class 4 counts as background)
//   BUMP   = classes 2 + 3 in the bump window 0.4576-0.5255 GeV/c2; B = other T1 + T2 rows of the window
// FOM F = S / sqrt(S + B) in expected counts at L_data; purity S / (S + B); efficiency S / S_pre.
//
// Method (exact exhaustive grid search):
//   Nine cut variables with grids that contain the named presets: OA > a, alpha_lab < b, M(pipi) > c,
//   M(pipi) < d, alpha_eta > e, both pions in the vertex fit (off/on), |z_pi - z_vertex| < 40 mm for both pions
//   (off/on), MM4 < f, pT4 > g. Every row of a window gets, per variable, the index of the grid cell its value
//   falls into (bin edges = the cut values, compared in float exactly as the SimEngine presets); the weights
//   are filled into a 9-dimensional histogram (13.8 M cells x 7 quantities x 2 halves, double), which is
//   then summed cumulatively along every axis (suffix sums for lower cuts and flags, prefix sums for upper
//   cuts). After that each cell holds the exact S and B of one cut combination, for both halves. Every cell is
//   evaluated (no approximation at the grid values). Checked against direct row loops and against the engine.
// Plateau: per variable the profile (best F over all other variables) and the slice through the chosen point;
//   the plateau is the set of values with F >= 0.95 F_best.
// Freeze rule ("simple and robust"): among the cells with F >= 0.95 F_best on the full sample AND on each half
//   (relative to that half's own best), take the smallest number of active cuts; among those the cell whose
//   worst grid neighbour (each active threshold moved by one grid step) keeps the highest F; ties by F. A
//   neighbour minimum >= 0.90 F_best is required (else the next number of cuts is tried).
// Purity alternative (rule P; reported, NOT frozen): on the same plateau, at most 3 active cuts and the same
//   neighbour condition, the cell of highest purity S / (S + B).
// Split-half check: half A = event groups with an even ordinal (grp, per tree = per set; the 'fold' column of
//   the compact files), half B = odd. Each half is scaled by 2. Optimise on A (best cell and freeze rule on A
//   alone), evaluate on B, and vice versa: bias = F(optimising half) - F(other half) at the chosen cuts. The
//   same is repeated for K random splits of the event groups (hash of set, grp, seed) to give the mean bias and
//   its spread.
// Errors: event-group variance (SimEngine rule; rows of one event group correlated), delta method for F,
//   group-correlated formulas for purity and efficiency. These are the errors from the simulation statistics.
//
// Output: figs/s5_opt_eta_{profile_peak,profile_direct,profile_bump,presets,split,m4,composition,maps}.png/.pdf,
//         tables/s5_opt_eta_results.txt, tables/s5_opt_eta_frozen.txt
// Run from /home/przygoda/HADES/FAT/FAT (single thread; about 10 GB during the load, then about 2.5 GB; about 30 s + 6 s per split):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_opt_eta", kTRUE)' 'plots/eta_sim/s5_opt_eta.C+(100)'
// argument: number K of random splits (default 100)
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"
#include "THStack.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace s5o {
using namespace simeng;

// ============================================================================
// constants
// ============================================================================
const double kPkLo = 0.5333, kPkHi = 0.5587;   // peak window [GeV/c2]
const double kBuLo = 0.4576, kBuHi = 0.5255;   // bump window [GeV/c2]
const double kLd = 6080.5;                     // nb^-1, UNCONFIRMED (scale only)
const double kK = kLd / kLsim;                 // sum w -> expected counts at L_data
const char* const kTag = "s5_opt_eta";
const char* const kLdNote = "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)";
const double kPlat = 0.95;                     // plateau: F >= 0.95 F_best
const double kRobMin = 0.90;                   // freeze rule: worst neighbour >= 0.90 F_best
const int NT = 3;                              // targets
const char* TName(int t) { static const char* n[NT] = {"PEAK", "DIRECT", "BUMP"}; return n[t]; }
const char* TDesc(int t) {
    static const char* n[NT] = {"classes 4+5 in the peak window 0.5333-0.5587", "class 5 in the peak window 0.5333-0.5587 (class 4 = background)",
                                "classes 2+3 in the bump window 0.4576-0.5255"};
    return n[t];
}
const char* TLab(int t) {
    static const char* n[NT] = {"PEAK: classes 4+5, peak window", "DIRECT: class 5, peak window", "BUMP: classes 2+3, bump window"};
    return n[t];
}

FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(gOut ? gOut : stdout, fmt, ap);
    va_end(ap);
}

// ============================================================================
// groups of rows (fixed for all agents)
// ============================================================================
enum Grp { gS2 = 0, gS3, gS4, gS5, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, NG };
const char* GName(int g) {
    static const char* n[NG] = {"Gsig2 bump pi0 Dalitz", "Gsig3 bump conversion", "Gsig4 peak conversion", "Gsig5 peak direct",
                                "G_eta (P2 or L2)", "G_pi0D (L1)", "G_conv (L3)", "G_vec (L4-L6)", "G_oth", "G_fake (T3)",
                                "G_diff (T4)"};
    return (g >= 0 && g < NG) ? n[g] : "?";
}
const char* GShort(int g) {
    static const char* n[NG] = {"Gsig2", "Gsig3", "Gsig4", "Gsig5", "G_eta", "G_pi0D", "G_conv", "G_vec", "G_oth", "G_fake", "G_diff"};
    return (g >= 0 && g < NG) ? n[g] : "?";
}
inline int GroupOf(const Row& r) {
    if (r.tmpl == 1) return gS2 + (r.cls - 2);   // classes 2-5 (X17 rows are dropped before)
    if (r.tmpl == 2) {
        if (r.pcls == 2 || r.lcls == 2) return gEta;
        if (r.lcls == 1) return gPi0D;
        if (r.lcls == 3) return gConv;
        if (r.lcls >= 4 && r.lcls <= 6) return gVec;
        return gOth;
    }
    if (r.tmpl == 3) return gFake;
    return gDiff;
}

// ============================================================================
// cut dimensions
// ============================================================================
enum DimType { kLower = 0, kUpper = 1, kFlag = 2 };
struct Dim {
    std::string key, title, unit;
    int type;
    std::vector<float> g;   // lower: g[0] = none (-inf), g[1..] ascending; upper: g[0..n-2] ascending, g[n-1] = none; flag: 0, 1
    int n() const { return (int)g.size(); }
    int none() const { return type == kUpper ? n() - 1 : 0; }
};
const int ND = 9;
enum { dOA = 0, dLab, dMlo, dMhi, dEta, dVtx, dDz, dMM, dPt };
std::vector<Dim> gD;
void MakeDims() {
    const float NI = -std::numeric_limits<float>::infinity(), PI = std::numeric_limits<float>::infinity();
    gD.clear();
    gD.push_back({"oa_min", "OA >", "deg", kLower, {NI, 1.f, 1.5f, 2.f, 2.5f, 3.f, 4.f, 5.f, 6.f, 8.f}});
    gD.push_back({"alpha_lab_max", "#alpha_{lab} <", "deg", kUpper, {25.f, 30.f, 35.f, 40.f, 45.f, 50.f, 55.f, 60.f, 70.f, PI}});
    gD.push_back({"mpipi_min", "M_{#pi#pi} >", "GeV/c^{2}", kLower, {NI, 0.28f, 0.30f, 0.32f, 0.34f, 0.36f, 0.38f, 0.40f, 0.42f}});
    gD.push_back({"mpipi_max", "M_{#pi#pi} <", "GeV/c^{2}", kUpper, {0.38f, 0.40f, 0.41f, 0.42f, 0.44f, 0.46f, 0.48f, 0.50f, 0.53f, PI}});
    gD.push_back({"alpha_eta_min", "#alpha_{#eta} >", "deg", kLower, {NI, 100.f, 120.f, 130.f, 140.f, 150.f, 160.f, 170.f}});
    gD.push_back({"pions_in_vertex_fit", "both #pi in vertex fit", "", kFlag, {0.f, 1.f}});
    gD.push_back({"dz_pi_lt_40", "|z_{#pi} - z_{vtx}| < 40 mm", "", kFlag, {0.f, 1.f}});
    gD.push_back({"mm4_max", "MM_{4} <", "GeV/c^{2}", kUpper, {2.50f, 2.60f, 2.70f, 2.75f, 2.80f, 2.85f, 2.90f, PI}});
    gD.push_back({"pt4_min", "p_{T,4} >", "MeV/c", kLower, {NI, 100.f, 200.f, 300.f, 400.f, 500.f}});
}
std::string ValStr(int d, int j) {
    const Dim& D = gD[d];
    if (D.type == kFlag) return j ? "on" : "off";
    if (j == D.none()) return "none";
    return Form("%g", D.g[j]);
}
// index of a value in the grid of dimension d (exact float semantics of the cuts)
inline int IndexOf(int d, float x) {
    const Dim& D = gD[d];
    const int n = D.n();
    if (D.type == kFlag) return x > 0.5f ? 1 : 0;
    if (D.type == kLower) {   // largest j with x > g[j]; NaN -> 0
        int idx = 0;
        for (int j = 1; j < n; ++j) if (x > D.g[j]) idx = j;
        return idx;
    }
    for (int j = 0; j < n - 1; ++j) if (x < D.g[j]) return j;   // upper: smallest j with x < g[j]; NaN or above -> none
    return n - 1;
}
inline bool PassDim(int d, int ri, int cj) { return gD[d].type == kUpper ? ri <= cj : ri >= cj; }

struct Cut { int j[ND]; };
Cut NoCut() { Cut c; for (int d = 0; d < ND; ++d) c.j[d] = gD[d].none(); return c; }
int Find(int d, float v) {   // grid index of a cut value
    for (int j = 0; j < gD[d].n(); ++j) if (gD[d].g[j] == v) return j;
    printf("s5o::Find: value %g not in the grid of %s\n", v, gD[d].key.c_str());
    return gD[d].none();
}
std::string CutStr(const Cut& c) {
    std::string s;
    for (int d = 0; d < ND; ++d) {
        if (c.j[d] == gD[d].none()) continue;
        std::string x;
        switch (d) {
            case dOA: x = "OA>" + ValStr(d, c.j[d]); break;
            case dLab: x = "lab<" + ValStr(d, c.j[d]); break;
            case dMlo: x = "Mpp>" + ValStr(d, c.j[d]); break;
            case dMhi: x = "Mpp<" + ValStr(d, c.j[d]); break;
            case dEta: x = "eta>" + ValStr(d, c.j[d]); break;
            case dVtx: x = "pi-vtx"; break;
            case dDz: x = "dz<40"; break;
            case dMM: x = "MM4<" + ValStr(d, c.j[d]); break;
            case dPt: x = "pT4>" + ValStr(d, c.j[d]); break;
        }
        s += (s.empty() ? "" : ", ") + x;
    }
    return s.empty() ? "pre (no cut)" : s;
}
int Complexity(const Cut& c) { int k = 0; for (int d = 0; d < ND; ++d) k += c.j[d] != gD[d].none(); return k; }

// named presets as grid cells
const int NPS = 7;
const char* PsName(int p) {
    static const char* n[NPS] = {"pre", "pre+OA>2", "pre+pions in vertex fit", "slides", "old code", "peak nominal", "bump nominal"};
    return n[p];
}
Cut PresetCut(int p) {
    Cut c = NoCut();
    switch (p) {
        case kPre: break;
        case kPreOA2: c.j[dOA] = Find(dOA, 2.f); break;
        case kPrePiVtx: c.j[dVtx] = 1; break;
        case kSlides: c.j[dLab] = Find(dLab, 50.f); c.j[dMhi] = Find(dMhi, 0.46f); c.j[dEta] = Find(dEta, 140.f); c.j[dOA] = Find(dOA, 2.f); break;
        case kOldCode: c.j[dLab] = Find(dLab, 50.f); c.j[dMhi] = Find(dMhi, 0.42f); c.j[dEta] = Find(dEta, 140.f); break;
        case kPeakNom: c.j[dLab] = Find(dLab, 55.f); c.j[dMlo] = Find(dMlo, 0.30f); break;
        case kBumpNom: c.j[dLab] = Find(dLab, 40.f); c.j[dMlo] = Find(dMlo, 0.30f); c.j[dMhi] = Find(dMhi, 0.41f); break;
    }
    return c;
}

// ============================================================================
// compact rows (all "pre" rows owned in mode A with w > 0, X17 dropped)
// ============================================================================
struct CRow {
    float m4;
    double w;           // w (PDG) x physics factor (main)
    int grp;
    unsigned char set, fold, g, win, q;
    unsigned char ix[ND];
};
std::vector<CRow> gR;
// quantity index: peak window 0 c4, 1 c5, 2 other T1+T2, 3 T3+T4; bump window 4 c2+c3, 5 other T1+T2, 6 T3+T4; 255 none
inline int QOf(int win, int g) {
    if (win == 1) { if (g == gS4) return 0; if (g == gS5) return 1; if (g <= gOth) return 2; return 3; }
    if (win == 2) { if (g == gS2 || g == gS3) return 4; if (g <= gOth) return 5; return 6; }
    return 255;
}
inline bool PassCut(const CRow& r, const Cut& c) {
    for (int d = 0; d < ND; ++d) if (!PassDim(d, r.ix[d], c.j[d])) return false;
    return true;
}
// signal / background membership of a quantity for a target (v = 1: B includes T3 + T4)
inline bool IsSig(int t, int q) { return t == 0 ? (q == 0 || q == 1) : t == 1 ? q == 1 : q == 4; }
inline bool IsBkg(int t, int q, int v) {
    if (t == 0) return q == 2 || (v && q == 3);
    if (t == 1) return q == 0 || q == 2 || (v && q == 3);
    return q == 5 || (v && q == 6);
}
inline int WinOf(int t) { return t < 2 ? 1 : 2; }

// ============================================================================
// the 9-D table
// ============================================================================
const int NQ = 7, NQ2 = 14;
struct Table {
    int n[ND];
    long st[ND];
    long ncell = 0;
    std::vector<double> T;
    void Init() {
        ncell = 1;
        for (int d = ND - 1; d >= 0; --d) { n[d] = gD[d].n(); st[d] = ncell; ncell *= n[d]; }
        T.assign((size_t)ncell * NQ2, 0.);
    }
    long CellOf(const unsigned char* ix) const { long c = 0; for (int d = 0; d < ND; ++d) c += ix[d] * st[d]; return c; }
    long CellOf(const Cut& k) const { long c = 0; for (int d = 0; d < ND; ++d) c += k.j[d] * st[d]; return c; }
    Cut Decode(long c) const { Cut k; for (int d = 0; d < ND; ++d) k.j[d] = (int)((c / st[d]) % n[d]); return k; }
    // fill the window rows (wi: indices into gR) with the halves h (0 = A, 1 = B)
    void Fill(const std::vector<size_t>& wi, const std::vector<long>& wc, const std::vector<unsigned char>& half) {
        std::fill(T.begin(), T.end(), 0.);
        for (size_t k = 0; k < wi.size(); ++k) {
            const CRow& r = gR[wi[k]];
            T[(size_t)wc[k] * NQ2 + half[k] * NQ + r.q] += r.w;
        }
    }
    void Cumulate() {
        for (int d = 0; d < ND; ++d) {
            const long s = st[d];
            const int nd = n[d];
            if (gD[d].type == kUpper) {
                for (long c = 0; c < ncell; ++c) {
                    if ((c / s) % nd == 0) continue;
                    double* a = &T[(size_t)c * NQ2];
                    const double* b = &T[(size_t)(c - s) * NQ2];
                    for (int q = 0; q < NQ2; ++q) a[q] += b[q];
                }
            } else {
                for (long c = ncell - 1; c >= 0; --c) {
                    if ((c / s) % nd == nd - 1) continue;
                    double* a = &T[(size_t)c * NQ2];
                    const double* b = &T[(size_t)(c + s) * NQ2];
                    for (int q = 0; q < NQ2; ++q) a[q] += b[q];
                }
            }
        }
    }
    // sum w (unscaled) of quantity q in sample s (0 full, 1 A, 2 B)
    inline double Q(long c, int s, int q) const {
        const double* a = &T[(size_t)c * NQ2];
        return s == 0 ? a[q] + a[NQ + q] : s == 1 ? a[q] : a[NQ + q];
    }
    // S, B in expected counts at L_data (halves x 2)
    inline void SB(long c, int t, int s, int v, double& S, double& B) const {
        const double* a = &T[(size_t)c * NQ2];
        double q[NQ];
        for (int i = 0; i < NQ; ++i) q[i] = s == 0 ? a[i] + a[NQ + i] : s == 1 ? 2 * a[i] : 2 * a[NQ + i];
        if (t == 0) { S = q[0] + q[1]; B = q[2] + (v ? q[3] : 0.); }
        else if (t == 1) { S = q[1]; B = q[0] + q[2] + (v ? q[3] : 0.); }
        else { S = q[4]; B = q[5] + (v ? q[6] : 0.); }
        S *= kK; B *= kK;
    }
    inline double F(long c, int t, int s, int v = 0) const {
        double S, B;
        SB(c, t, s, v, S, B);
        return (S > 0 && S + B > 0) ? S / std::sqrt(S + B) : 0.;
    }
};
Table gT;

// ============================================================================
// evaluation of one cut with event-group errors (direct row loop)
// ============================================================================
struct Ev {
    double S = 0, B = 0, F = 0, eF = 0, pur = 0, epur = 0, eff = 0, eeff = 0;  // S, B: expected counts at L_data
    double sS = 0, sB = 0, VS = 0, VB = 0, CSB = 0, sSpre = 0;                // sums of w (unscaled, sample weight included)
    double nS = 0, nB = 0;
    double sSs[4] = {0, 0, 0, 0}, vSs[4] = {0, 0, 0, 0}, nSs[4] = {0, 0, 0, 0};
    double sBs[4] = {0, 0, 0, 0}, vBs[4] = {0, 0, 0, 0}, nBs[4] = {0, 0, 0, 0};
    double sG[NG], vG[NG], nG[NG];
    Ev() { for (int g = 0; g < NG; ++g) sG[g] = vG[g] = nG[g] = 0; }
    double NeffS() const { return VS > 0 ? sS * sS / VS : 0; }
    double NeffB() const { return VB > 0 ? sB * sB / VB : 0; }
};
// half: per gR row (0 A, 1 B) or nullptr; h = -1 full sample, 0 or 1 = that half (scaled x 2)
Ev Evaluate(const Cut& c, int t, int v, const std::vector<unsigned char>* half = nullptr, int h = -1) {
    Ev e;
    const int win = WinOf(t);
    const double hs = h >= 0 ? 2. : 1.;
    double sP = 0, sF = 0, bP = 0, gs[NG];
    for (int g = 0; g < NG; ++g) gs[g] = 0;
    int curSet = -1, curGrp = -1, gset = 0;
    double Vpur = 0, Veff = 0;
    std::vector<double> gsP, gsF, gbP;   // per group (for the purity / efficiency variance after the means are known)
    auto flush = [&]() {
        if (sP != 0 || sF != 0 || bP != 0) {
            e.VS += sP * sP; e.VB += bP * bP; e.CSB += sP * bP;
            e.vSs[gset] += sP * sP; e.vBs[gset] += bP * bP;
            gsP.push_back(sP); gsF.push_back(sF); gbP.push_back(bP);
        }
        for (int g = 0; g < NG; ++g) if (gs[g] != 0) { e.vG[g] += gs[g] * gs[g]; gs[g] = 0; }
        sP = sF = bP = 0;
    };
    for (size_t i = 0; i < gR.size(); ++i) {
        const CRow& r = gR[i];
        if (r.set != curSet || r.grp != curGrp) { flush(); curSet = r.set; curGrp = r.grp; gset = r.set; }
        if (r.win != win) continue;
        if (h >= 0 && (*half)[i] != h) continue;
        const double w = r.w * hs;
        const bool ps = PassCut(r, c);
        const bool sig = IsSig(t, r.q), bkg = IsBkg(t, r.q, v);
        if (sig) {
            e.sSpre += w;
            if (ps) { sP += w; e.sS += w; e.nS += 1; e.sSs[r.set] += w; e.nSs[r.set] += 1; }
            else sF += w;
        }
        if (bkg && ps) { bP += w; e.sB += w; e.nB += 1; e.sBs[r.set] += w; e.nBs[r.set] += 1; }
        if (ps) { e.sG[r.g] += w; e.nG[r.g] += 1; gs[r.g] += w; }
    }
    flush();
    e.S = kK * e.sS; e.B = kK * e.sB;
    const double SB = e.S + e.B;
    if (e.S > 0 && SB > 0) {
        e.F = e.S / std::sqrt(SB);
        const double dS = (e.S + 2 * e.B) / (2 * std::pow(SB, 1.5)), dB = -e.S / (2 * std::pow(SB, 1.5));
        const double v2 = dS * dS * e.VS + dB * dB * e.VB + 2 * dS * dB * e.CSB;
        e.eF = kK * std::sqrt(std::max(0., v2));
    }
    const double tot = e.sS + e.sB;
    if (tot > 0) e.pur = e.sS / tot;
    if (e.sSpre > 0) e.eff = e.sS / e.sSpre;
    for (size_t k = 0; k < gsP.size(); ++k) {
        const double a = (1 - e.pur) * gsP[k] - e.pur * gbP[k];
        Vpur += a * a;
        const double b = (1 - e.eff) * gsP[k] - e.eff * gsF[k];
        Veff += b * b;
    }
    e.epur = tot > 0 ? std::sqrt(Vpur) / tot : 0;
    e.eeff = e.sSpre > 0 ? std::sqrt(Veff) / e.sSpre : 0;
    return e;
}

// difference F(a) - F(b) of two cuts on the full sample with its event-group error (rows shared by both cuts are
// correlated; delta method per event group). Added 2026-11-13 (review finding M5).
struct DF { double d = 0, e = 0; };
DF DiffF(const Cut& ca, const Cut& cb, int t) {
    const Ev ea = Evaluate(ca, t, 0), eb = Evaluate(cb, t, 0);
    auto der = [](const Ev& e, double& dS, double& dB) {
        const double SB = e.S + e.B;
        dS = (e.S + 2 * e.B) / (2 * std::pow(SB, 1.5)); dB = -e.S / (2 * std::pow(SB, 1.5));
    };
    double dSa, dBa, dSb, dBb;
    der(ea, dSa, dBa); der(eb, dSb, dBb);
    const int win = WinOf(t);
    double sa = 0, ba = 0, sb = 0, bb = 0, V = 0;
    int curSet = -1, curGrp = -1;
    auto flush = [&]() { const double x = dSa * sa + dBa * ba - dSb * sb - dBb * bb; V += x * x; sa = ba = sb = bb = 0; };
    for (size_t i = 0; i < gR.size(); ++i) {
        const CRow& r = gR[i];
        if (r.set != curSet || r.grp != curGrp) { flush(); curSet = r.set; curGrp = r.grp; }
        if (r.win != win) continue;
        const bool sig = IsSig(t, r.q), bkg = IsBkg(t, r.q, 0);
        if (!sig && !bkg) continue;
        if (PassCut(r, ca)) { if (sig) sa += r.w; else ba += r.w; }
        if (PassCut(r, cb)) { if (sig) sb += r.w; else bb += r.w; }
    }
    flush();
    DF out; out.d = ea.F - eb.F; out.e = kK * std::sqrt(V);
    return out;
}

// ============================================================================
// search over the table
// ============================================================================
struct Best { long cell = 0; double F = -1; };
struct Search {
    Best best[NT][3][2];                  // target, sample, variant
    std::vector<double> prof[NT][3][ND];  // profile: max F per grid value (variant 0)
    std::vector<long> profC[NT][ND];      // cell of the full-sample profile maximum (for its error band; 2026-11-13)
    std::vector<double> profV[NT][ND];    // profile of the variant (full sample)
    std::vector<double> prof2[NT][2];     // 2-D profiles (full): (lab, eta) and (Mpp lo, Mpp hi)
};
const int kPair[2][2] = {{dLab, dEta}, {dMlo, dMhi}};
// full = true: all samples, variants, profiles; false: halves only, variant 0, no profiles (random splits)
void RunSearch(const Table& tb, Search& R, bool full) {
    for (int t = 0; t < NT; ++t) {
        for (int s = 0; s < 3; ++s) {
            for (int v = 0; v < 2; ++v) R.best[t][s][v] = Best();
            for (int d = 0; d < ND; ++d) R.prof[t][s][d].assign(gD[d].n(), 0.);
        }
        for (int d = 0; d < ND; ++d) R.profC[t][d].assign(gD[d].n(), 0);
        for (int d = 0; d < ND; ++d) R.profV[t][d].assign(gD[d].n(), 0.);
        for (int p = 0; p < 2; ++p) R.prof2[t][p].assign(gD[kPair[p][0]].n() * gD[kPair[p][1]].n(), 0.);
    }
    int ix[ND] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (long c = 0; c < tb.ncell; ++c) {
        for (int t = 0; t < NT; ++t) {
            for (int s = full ? 0 : 1; s < 3; ++s) {
                const double f = tb.F(c, t, s, 0);
                if (f > R.best[t][s][0].F) { R.best[t][s][0].F = f; R.best[t][s][0].cell = c; }
                if (full) for (int d = 0; d < ND; ++d) { double& p = R.prof[t][s][d][ix[d]]; if (f > p) { p = f; if (s == 0) R.profC[t][d][ix[d]] = c; } }
                if (full && s == 0)
                    for (int p = 0; p < 2; ++p) {
                        double& q = R.prof2[t][p][ix[kPair[p][0]] * gD[kPair[p][1]].n() + ix[kPair[p][1]]];
                        if (f > q) q = f;
                    }
            }
            if (full) {
                const double f = tb.F(c, t, 0, 1);
                if (f > R.best[t][0][1].F) { R.best[t][0][1].F = f; R.best[t][0][1].cell = c; }
                for (int d = 0; d < ND; ++d) { double& p = R.profV[t][d][ix[d]]; if (f > p) p = f; }
            }
        }
        for (int d = ND - 1; d >= 0; --d) { if (++ix[d] < gD[d].n()) break; ix[d] = 0; }
    }
}
// worst grid neighbour of the active threshold cuts (flags have no neighbour), relative to fbest
double Robustness(const Table& tb, long c, int t, int s, double fbest) {
    const Cut k = tb.Decode(c);
    double m = tb.F(c, t, s) / fbest;
    for (int d = 0; d < ND; ++d) {
        if (gD[d].type == kFlag || k.j[d] == gD[d].none()) continue;
        if (k.j[d] > 0) m = std::min(m, tb.F(c - tb.st[d], t, s) / fbest);
        if (k.j[d] < gD[d].n() - 1) m = std::min(m, tb.F(c + tb.st[d], t, s) / fbest);
    }
    return m;
}
struct Pick { long cell = 0; double F = 0, rob = 0; int cx = 0; bool halvesOK = false, robOK = false; long nPlat = 0, nPlatK = 0; };
// freeze rule on sample s; halves = true: the cell must also be on the plateau of each half
Pick RuleR(const Table& tb, int t, int s, const double fb[3], bool halves) {
    Pick out;
    for (int pass = 0; pass < 2; ++pass) {   // pass 1 (only if pass 0 found nothing): without the half condition
        const bool useH = halves && pass == 0;
        int kmin = 99;
        long nplat = 0;
        int ix[ND] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<int> cxOf;   // not stored: recompute
        for (long c = 0; c < tb.ncell; ++c) {
            if (tb.F(c, t, s) >= kPlat * fb[s] && (!useH || (tb.F(c, t, 1) >= kPlat * fb[1] && tb.F(c, t, 2) >= kPlat * fb[2]))) {
                ++nplat;
                int k = 0;
                for (int d = 0; d < ND; ++d) k += ix[d] != gD[d].none();
                kmin = std::min(kmin, k);
            }
            for (int d = ND - 1; d >= 0; --d) { if (++ix[d] < gD[d].n()) break; ix[d] = 0; }
        }
        if (!nplat) continue;
        for (int K = kmin; K <= ND; ++K) {
            Pick best;
            best.rob = -1;
            long nk = 0;
            int jx[ND] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (long c = 0; c < tb.ncell; ++c) {
                int k = 0;
                for (int d = 0; d < ND; ++d) k += jx[d] != gD[d].none();
                if (k == K) {
                    const double f = tb.F(c, t, s);
                    if (f >= kPlat * fb[s] && (!useH || (tb.F(c, t, 1) >= kPlat * fb[1] && tb.F(c, t, 2) >= kPlat * fb[2]))) {
                        ++nk;
                        const double rb = Robustness(tb, c, t, s, fb[s]);
                        if (rb > best.rob + 1e-12 || (std::fabs(rb - best.rob) <= 1e-12 && f > best.F)) {
                            best.rob = rb; best.F = f; best.cell = c; best.cx = K;
                        }
                    }
                }
                for (int d = ND - 1; d >= 0; --d) { if (++jx[d] < gD[d].n()) break; jx[d] = 0; }
            }
            if (nk && best.rob >= kRobMin) {
                best.halvesOK = useH; best.robOK = true; best.nPlat = nplat; best.nPlatK = nk;
                return best;
            }
            if (nk && best.rob > out.rob) { out = best; out.halvesOK = useH; out.nPlat = nplat; out.nPlatK = nk; }
        }
        if (out.rob > 0) return out;   // plateau found but no robust cell: best robustness overall
    }
    return out;
}
// purity alternative (rule P, reported, not frozen): on the same plateau (and half condition), at most maxCuts active
// cuts and a worst neighbour >= kRobMin, the cell of highest purity S / (S + B)
inline double Purity(const Table& tb, long c, int t, int s) {
    double S, B;
    tb.SB(c, t, s, 0, S, B);
    return S + B > 0 ? S / (S + B) : 0.;
}
Pick RuleP(const Table& tb, int t, int s, const double fb[3], bool halves, int maxCuts) {
    Pick out;
    double bestP = -1;
    int ix[ND] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (long c = 0; c < tb.ncell; ++c) {
        int k = 0;
        for (int d = 0; d < ND; ++d) k += ix[d] != gD[d].none();
        if (k <= maxCuts) {
            const double f = tb.F(c, t, s);
            if (f >= kPlat * fb[s] && (!halves || (tb.F(c, t, 1) >= kPlat * fb[1] && tb.F(c, t, 2) >= kPlat * fb[2]))) {
                const double p = Purity(tb, c, t, s);
                if (p > bestP) {
                    const double rb = Robustness(tb, c, t, s, fb[s]);
                    if (rb >= kRobMin) { bestP = p; out.cell = c; out.F = f; out.rob = rb; out.cx = k; out.halvesOK = halves; out.robOK = true; }
                }
            }
        }
        for (int d = ND - 1; d >= 0; --d) { if (++ix[d] < gD[d].n()) break; ix[d] = 0; }
    }
    if (bestP < 0) { out.cell = 0; out.F = 0; }
    return out;
}
const int kMaxCutsP = 3;
const int NRU = 3;   // 0 best cell, 1 freeze rule R, 2 purity rule P
const char* RuName(int r) { static const char* n[NRU] = {"best cell", "freeze rule", "purity rule"}; return n[r]; }

// ============================================================================
// random splits
// ============================================================================
inline uint64_t SplitMix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline unsigned char HalfOf(const CRow& r, int seed) {
    if (seed == 0) return r.fold;   // deterministic: parity of the event-group ordinal per set
    return (unsigned char)(SplitMix(((uint64_t)r.set << 40) ^ (uint64_t)(uint32_t)r.grp ^ ((uint64_t)seed << 48)) & 1ULL);
}
struct SplitRes {
    // [rule][target][direction 0 = opt A eval B, 1 = opt B eval A]: F on the optimising half, F on the other half,
    // purity on both, and the chosen cell
    double fOpt[NRU][NT][2], fOth[NRU][NT][2], pOpt[NRU][NT][2], pOth[NRU][NT][2];
    double fFull[NRU][NT][2];   // F of the full sample (A + B; independent of the split) at the half pick
    long cell[NRU][NT][2];
};
SplitRes DoSplit(int seed, const std::vector<size_t>& wi, const std::vector<long>& wc, Search* keep = nullptr) {
    std::vector<unsigned char> half(wi.size());
    for (size_t k = 0; k < wi.size(); ++k) half[k] = HalfOf(gR[wi[k]], seed);
    gT.Fill(wi, wc, half);
    gT.Cumulate();
    Search S;
    RunSearch(gT, S, keep != nullptr);
    SplitRes o;
    for (int t = 0; t < NT; ++t) {
        double fb[3] = {keep ? S.best[t][0][0].F : 1., S.best[t][1][0].F, S.best[t][2][0].F};
        for (int dir = 0; dir < 2; ++dir) {
            const int so = 1 + dir, sx = 2 - dir;
            long cs[NRU];
            cs[0] = S.best[t][so][0].cell;
            cs[1] = RuleR(gT, t, so, fb, false).cell;
            cs[2] = RuleP(gT, t, so, fb, false, kMaxCutsP).cell;
            for (int r = 0; r < NRU; ++r) {
                o.cell[r][t][dir] = cs[r];
                o.fOpt[r][t][dir] = gT.F(cs[r], t, so);
                o.fOth[r][t][dir] = gT.F(cs[r], t, sx);
                o.fFull[r][t][dir] = gT.F(cs[r], t, 0);
                o.pOpt[r][t][dir] = Purity(gT, cs[r], t, so);
                o.pOth[r][t][dir] = Purity(gT, cs[r], t, sx);
            }
        }
    }
    if (keep) *keep = S;
    return o;
}

// ============================================================================
// plotting helpers
// ============================================================================
void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetHistLineWidth(3);
}
void Save(TCanvas* c, const char* name) {
    const TString base = Form("%s/%s", kFigDir, name);
    c->SaveAs(base + ".png");
    c->SaveAs(base + ".eps");
    gSystem->Exec(Form("cd %s && ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", kFigDir, name, name, name));
}
void Txt(double x, double y, double size, const char* s, int align = 11, int col = 1, bool ndc = true) {
    TLatex l;
    if (ndc) l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(size);
    l.SetTextAlign(align);
    l.SetTextColor(col);
    l.DrawLatex(x, y, s);
}
void AxisStyle(TH1* h, double ls, double ts, double xo = 1.0, double yo = 1.0) {
    for (TAxis* a : {h->GetXaxis(), h->GetYaxis()}) {
        a->SetLabelSize(ls); a->SetTitleSize(ts); a->SetLabelFont(42); a->SetTitleFont(42);
    }
    h->GetXaxis()->SetTitleOffset(xo);
    h->GetYaxis()->SetTitleOffset(yo);
}

// ============================================================================
// main
// ============================================================================
void Run(int K) {
    TStopwatch sw;
    Style();
    MakeDims();
    // ---------------------------------------------------------------- weights and rows
    if (!simphys::Init()) { printf("simphys::Init failed\n"); return; }
    simphys::Install();
    if (!simphys::LoadAngles()) { printf("simphys::LoadAngles failed\n"); return; }
    SimOpt o;
    o.phys = simphys::kPhysMain;   // mode A, w (PDG), sets 7, e+e-
    SimData* d = new SimData;
    if (!Load(*d, kCompact, kMaskUS, -1e9, 1e9, 1)) return;
    long nZero = 0, nX17 = 0, nPre = 0, nZeroSet[4] = {0, 0, 0, 0}, nZeroTmpl[5] = {0, 0, 0, 0, 0}, nKeptSet[4] = {0, 0, 0, 0};
    for (int set = 1; set <= 3; ++set) {
        for (const Row& r : d->tree(set, kUS)) {
            if (!PassPre(r) || r.ownA != 1) continue;
            ++nPre;
            if (r.tmpl == 1 && r.cls == 6) { ++nX17; continue; }
            const double w = Weight(r, o);
            if (!(w > 0)) { ++nZero; ++nZeroSet[set]; if (r.tmpl >= 1 && r.tmpl <= 4) ++nZeroTmpl[r.tmpl]; continue; }
            ++nKeptSet[set];
            CRow c;
            c.m4 = r.m4; c.w = w; c.grp = r.grp; c.set = (unsigned char)set; c.fold = (unsigned char)(r.fold & 1);
            c.g = (unsigned char)GroupOf(r);
            c.win = (r.m4 >= kPkLo && r.m4 < kPkHi) ? 1 : (r.m4 >= kBuLo && r.m4 < kBuHi) ? 2 : 0;
            c.q = (unsigned char)QOf(c.win, c.g);
            const float dz = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
            const float vals[ND] = {r.oa, r.oa_lab, r.mpipi, r.mpipi, r.oa_eta,
                                    (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1) ? 1.f : 0.f,
                                    dz < 40.f ? 1.f : 0.f, r.mm4, r.pt4};
            for (int k = 0; k < ND; ++k) c.ix[k] = (unsigned char)IndexOf(k, vals[k]);
            gR.push_back(c);
        }
    }
    printf("compact rows: %zu (pre and owned A: %ld; X17 dropped %ld; w <= 0 dropped %ld); angles hit %lld miss %lld\n", gR.size(), nPre,
           nX17, nZero, simphys::AngHit(), simphys::AngMiss());
    const std::string dropNote = Form("kept rows per set std %ld lep %ld eta %ld; w <= 0 rows per set std %ld lep %ld eta %ld, per template T1 %ld T2 %ld T3 %ld T4 %ld",
                                      nKeptSet[1], nKeptSet[2], nKeptSet[3], nZeroSet[1], nZeroSet[2], nZeroSet[3], nZeroTmpl[1], nZeroTmpl[2],
                                      nZeroTmpl[3], nZeroTmpl[4]);
    printf("%s\n", dropNote.c_str());

    // ---------------------------------------------------------------- engine cross-checks (before freeing the rows)
    std::vector<std::string> checks;
    int nCheckBad = 0;
    auto chk = [&](const char* what, double a, double b, double tol) {
        const double rel = std::fabs(a - b) / std::max(1e-300, std::fabs(b));
        const bool ok = rel <= tol;
        if (!ok) ++nCheckBad;
        checks.push_back(Form("%-70s mine %.10g engine %.10g rel %.2e %s", what, a, b, rel, ok ? "OK" : "FAIL"));
    };
    // engine sums at pre and at the presets in both windows; compared after the table is built
    struct EngRef { double S[NT], B[NT], V[NT]; };
    std::vector<EngRef> eng(NPS);
    {
        Comp cPk{"peak 4+5", [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5); }};
        Comp cC5{"c5", [](const Row& r) { return r.tmpl == 1 && r.cls == 5; }};
        Comp cBu{"bump 2+3", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); }};
        Comp cSim = CompSim();
        for (int p = 0; p < NPS; ++p) {
            Sel sPk = [p](const Row& r) { return PassPreset(p, r) && r.m4 >= kPkLo && r.m4 < kPkHi; };
            Sel sBu = [p](const Row& r) { return PassPreset(p, r) && r.m4 >= kBuLo && r.m4 < kBuHi; };
            const std::vector<WSum> a = SumMany(*d, sPk, {cPk, cC5, cSim}, o);
            const std::vector<WSum> b = SumMany(*d, sBu, {cBu, cSim}, o);
            eng[p].S[0] = a[0].S; eng[p].B[0] = a[2].S - a[0].S; eng[p].V[0] = a[0].V;
            eng[p].S[1] = a[1].S; eng[p].B[1] = a[2].S - a[1].S; eng[p].V[1] = a[1].V;
            eng[p].S[2] = b[0].S; eng[p].B[2] = b[1].S - b[0].S; eng[p].V[2] = b[0].V;
        }
    }
    delete d;
    d = nullptr;
    printf("rows ready after %.0f s\n", sw.RealTime()); sw.Continue();

    // ---------------------------------------------------------------- window rows and cells
    gT.Init();
    std::vector<size_t> wi;
    std::vector<long> wc;
    for (size_t i = 0; i < gR.size(); ++i)
        if (gR[i].win && gR[i].q != 255) { wi.push_back(i); wc.push_back(gT.CellOf(gR[i].ix)); }
    printf("table: %ld cells x %d doubles (%.2f GB); window rows %zu\n", gT.ncell, NQ2, gT.ncell * NQ2 * 8e-9, wi.size());
    std::vector<unsigned char> halfAll(gR.size());
    for (size_t i = 0; i < gR.size(); ++i) halfAll[i] = gR[i].fold;

    // ---------------------------------------------------------------- deterministic split + full search
    Search S0;
    const SplitRes R0 = DoSplit(0, wi, wc, &S0);
    printf("deterministic split and full search done after %.0f s\n", sw.RealTime()); sw.Continue();

    // checks: table against direct loops and the engine
    for (int p = 0; p < NPS; ++p) {
        const long c = gT.CellOf(PresetCut(p));
        for (int t = 0; t < NT; ++t) {
            double S, B;
            gT.SB(c, t, 0, 0, S, B);
            chk(Form("%s %s: S table vs engine", PsName(p), TName(t)), S / kK, eng[p].S[t], 1e-9);
            chk(Form("%s %s: B table vs engine", PsName(p), TName(t)), B / kK, eng[p].B[t], 1e-9);
            const Ev e = Evaluate(PresetCut(p), t, 0);
            chk(Form("%s %s: S direct loop vs engine", PsName(p), TName(t)), e.sS, eng[p].S[t], 1e-9);
            chk(Form("%s %s: V(S) group variance vs engine", PsName(p), TName(t)), e.VS, eng[p].V[t], 1e-9);
        }
    }
    {
        uint64_t seed = 12345;
        for (int k = 0; k < 30; ++k) {
            seed = SplitMix(seed);
            const long c = (long)(seed % (uint64_t)gT.ncell);
            const Cut cut = gT.Decode(c);
            for (int t = 0; t < NT; ++t)
                for (int h = 0; h < 2; ++h) {
                    const Ev e = Evaluate(cut, t, 1, &halfAll, h);
                    double S, B;
                    gT.SB(c, t, 1 + h, 1, S, B);
                    if (e.S > 0 || S > 0) chk(Form("random cell %ld %s half %c: S table vs loop", c, TName(t), 'A' + h), S, e.S, 1e-9);
                    if (e.B > 0 || B > 0) chk(Form("random cell %ld %s half %c: B(+T3T4) table vs loop", c, TName(t), 'A' + h), B, e.B, 1e-9);
                }
        }
    }
    printf("checks: %zu, failed %d\n", checks.size(), nCheckBad);

    // ---------------------------------------------------------------- freeze rule on the full sample
    Pick frz[NT], pur[NT];
    Ev evBest[NT], evFrz[NT], evPur[NT], evPre[NT], evPs[NT][NPS], evBestV[NT], evFrzV[NT], evPsV[NT][NPS];
    for (int t = 0; t < NT; ++t) {
        const double fb[3] = {S0.best[t][0][0].F, S0.best[t][1][0].F, S0.best[t][2][0].F};
        frz[t] = RuleR(gT, t, 0, fb, true);
        pur[t] = RuleP(gT, t, 0, fb, true, kMaxCutsP);
        evBest[t] = Evaluate(gT.Decode(S0.best[t][0][0].cell), t, 0);
        evFrz[t] = Evaluate(gT.Decode(frz[t].cell), t, 0);
        evPur[t] = Evaluate(gT.Decode(pur[t].cell), t, 0);
        evPre[t] = Evaluate(NoCut(), t, 0);
        evBestV[t] = Evaluate(gT.Decode(S0.best[t][0][1].cell), t, 1);
        evFrzV[t] = Evaluate(gT.Decode(frz[t].cell), t, 1);
        for (int p = 0; p < NPS; ++p) { evPs[t][p] = Evaluate(PresetCut(p), t, 0); evPsV[t][p] = Evaluate(PresetCut(p), t, 1); }
        printf("%s: best F %.3f (%s); frozen F %.3f (%s), %d cuts, robustness %.3f\n", TName(t), S0.best[t][0][0].F,
               CutStr(gT.Decode(S0.best[t][0][0].cell)).c_str(), frz[t].F, CutStr(gT.Decode(frz[t].cell)).c_str(), frz[t].cx, frz[t].rob);
    }
    // statistical error of the full-sample profile points (evaluated at the cell of each maximum) and the significance of
    // the differences frozen - pre, best - pre, pions-in-vertex - pre (2026-11-13, review finding M5)
    std::vector<double> profE[NT][ND];
    DF dFrzPre[NT], dBestPre[NT], dVtxPre[NT], dFrzVtx[NT];
    for (int t = 0; t < NT; ++t) {
        for (int dd = 0; dd < ND; ++dd) {
            profE[t][dd].assign(gD[dd].n(), 0.);
            for (int j = 0; j < gD[dd].n(); ++j) profE[t][dd][j] = Evaluate(gT.Decode(S0.profC[t][dd][j]), t, 0).eF;
        }
        dFrzPre[t] = DiffF(gT.Decode(frz[t].cell), NoCut(), t);
        dBestPre[t] = DiffF(gT.Decode(S0.best[t][0][0].cell), NoCut(), t);
        dVtxPre[t] = DiffF(PresetCut(kPrePiVtx), NoCut(), t);
        dFrzVtx[t] = DiffF(gT.Decode(frz[t].cell), PresetCut(kPrePiVtx), t);
        printf("%s: F(frozen) - F(pre) = %+.3f +- %.3f; F(best) - F(pre) = %+.3f +- %.3f; F(pions in vertex) - F(pre) = %+.3f +- %.3f; "
               "F(frozen) - F(pions in vertex) = %+.3f +- %.3f\n", TName(t), dFrzPre[t].d, dFrzPre[t].e, dBestPre[t].d, dBestPre[t].e,
               dVtxPre[t].d, dVtxPre[t].e, dFrzVtx[t].d, dFrzVtx[t].e);
    }
    // slices through the frozen point and the best point
    std::vector<double> slF[NT][ND], slB[NT][ND];
    for (int t = 0; t < NT; ++t)
        for (int dd = 0; dd < ND; ++dd) {
            slF[t][dd].assign(gD[dd].n(), 0.); slB[t][dd].assign(gD[dd].n(), 0.);
            Cut a = gT.Decode(frz[t].cell), b = gT.Decode(S0.best[t][0][0].cell);
            for (int j = 0; j < gD[dd].n(); ++j) {
                a.j[dd] = j; b.j[dd] = j;
                slF[t][dd][j] = gT.F(gT.CellOf(a), t, 0);
                slB[t][dd][j] = gT.F(gT.CellOf(b), t, 0);
            }
        }
    // deterministic split: evaluate with errors at the chosen cells
    Ev evSplit[NT][NRU][2][2];   // target, rule, direction, optimising / other half
    for (int t = 0; t < NT; ++t)
        for (int ru = 0; ru < NRU; ++ru)
            for (int dir = 0; dir < 2; ++dir) {
                const int ho = dir, hx = 1 - dir;
                evSplit[t][ru][dir][0] = Evaluate(gT.Decode(R0.cell[ru][t][dir]), t, 0, &halfAll, ho);
                evSplit[t][ru][dir][1] = Evaluate(gT.Decode(R0.cell[ru][t][dir]), t, 0, &halfAll, hx);
            }
    // values needed later from the deterministic table (it is overwritten by the random splits)
    double fHalfAtFrz[NT][2], fHalfBest[NT][2], fFullAtHalfPick[NT][NRU][2], fVarAtFrzRel[NT], fHalfAtPur[NT][2], fHalfAtFrzPre[NT][2];
    for (int t = 0; t < NT; ++t) {
        for (int h = 0; h < 2; ++h) {
            fHalfAtFrz[t][h] = gT.F(frz[t].cell, t, 1 + h);
            fHalfAtPur[t][h] = gT.F(pur[t].cell, t, 1 + h);
            fHalfAtFrzPre[t][h] = gT.F(gT.CellOf(NoCut()), t, 1 + h);
            fHalfBest[t][h] = S0.best[t][1 + h][0].F;
            for (int ru = 0; ru < NRU; ++ru) fFullAtHalfPick[t][ru][h] = gT.F(R0.cell[ru][t][h], t, 0);
        }
        fVarAtFrzRel[t] = gT.F(frz[t].cell, t, 0, 1) / S0.best[t][0][1].F;
    }
    // presets relative
    double fPs[NT][NPS], fPsV[NT][NPS];
    for (int t = 0; t < NT; ++t)
        for (int p = 0; p < NPS; ++p) { fPs[t][p] = gT.F(gT.CellOf(PresetCut(p)), t, 0); fPsV[t][p] = gT.F(gT.CellOf(PresetCut(p)), t, 0, 1); }

    // ---------------------------------------------------------------- spectra at pre and the frozen cuts (before the random splits)
    const int nb = 130;
    const double mlo = 0.35, mhi = 1.00, bw = (mhi - mlo) / nb;
    std::vector<double> spec[4][NG], specSum[4], specSumV[4];
    std::vector<Cut> specCut = {NoCut(), gT.Decode(frz[0].cell), gT.Decode(frz[1].cell), gT.Decode(frz[2].cell)};
    for (int k = 0; k < 4; ++k) {
        for (int g = 0; g < NG; ++g) spec[k][g].assign(nb, 0.);
        specSum[k].assign(nb, 0.); specSumV[k].assign(nb, 0.);
        int cs = -1, cg = -1;
        std::vector<std::pair<int, double>> buf;
        auto flush = [&]() {
            std::sort(buf.begin(), buf.end());
            size_t i = 0;
            while (i < buf.size()) { const int b = buf[i].first; double a = 0; while (i < buf.size() && buf[i].first == b) { a += buf[i].second; ++i; } specSumV[k][b] += a * a; }
            buf.clear();
        };
        for (const CRow& r : gR) {
            if (r.set != cs || r.grp != cg) { flush(); cs = r.set; cg = r.grp; }
            if (!(r.m4 >= mlo && r.m4 < mhi)) continue;
            if (!PassCut(r, specCut[k])) continue;
            const int b = std::min(nb - 1, (int)((r.m4 - mlo) / bw));
            spec[k][r.g][b] += r.w;
            if (r.g <= gOth) { specSum[k][b] += r.w; buf.push_back(std::make_pair(b, r.w)); }
        }
        flush();
    }

    // ---------------------------------------------------------------- random splits
    std::vector<SplitRes> RS;
    TStopwatch swr;
    for (int k = 1; k <= K; ++k) {
        RS.push_back(DoSplit(k, wi, wc));
        if (k == 1 || k % 5 == 0)
            printf("random split %d / %d done (%.0f s per split)\n", k, K, swr.RealTime() / k), swr.Continue();
    }
    // bias statistics: per split the mean of the two directions (the half-to-half offset at fixed cuts cancels to first
    // order); absolute (F units) and relative to the other half; purity bias; how often the half pick equals the full pick
    double mB[NT][NRU], eB[NT][NRU], rmsB[NT][NRU], mR[NT][NRU], eR[NT][NRU], rmsR[NT][NRU], mP[NT][NRU], eP[NT][NRU];
    double mLoss[NT][NRU], eLoss[NT][NRU];   // F_full(half pick) / F_full(best) - 1 (indicative; the half is part of the full sample)
    int nSame[NT][NRU];
    std::vector<double> relAvg[NT][NRU];     // per split, two-direction mean of the relative bias (histograms)
    auto ms = [](const std::vector<double>& v, double& m, double& e, double& r) {
        m = 0; for (double x : v) m += x; m /= std::max<size_t>(1, v.size());
        double s2 = 0; for (double x : v) s2 += (x - m) * (x - m);
        r = v.size() > 1 ? std::sqrt(s2 / (v.size() - 1)) : 0; e = v.size() ? r / std::sqrt((double)v.size()) : 0;
    };
    for (int t = 0; t < NT; ++t)
        for (int ru = 0; ru < NRU; ++ru) {
            std::vector<double> a, rel, pb, ls;
            nSame[t][ru] = 0;
            const long ref = ru == 0 ? S0.best[t][0][0].cell : ru == 1 ? frz[t].cell : pur[t].cell;
            for (const SplitRes& x : RS) {
                const double* opt = x.fOpt[ru][t];
                const double* oth = x.fOth[ru][t];
                a.push_back(0.5 * ((opt[0] - oth[0]) + (opt[1] - oth[1])));
                rel.push_back(0.5 * ((opt[0] / oth[0] - 1) + (opt[1] / oth[1] - 1)));
                pb.push_back(0.5 * ((x.pOpt[ru][t][0] - x.pOth[ru][t][0]) + (x.pOpt[ru][t][1] - x.pOth[ru][t][1])));
                relAvg[t][ru].push_back(rel.back());
                ls.push_back(0.5 * (x.fFull[ru][t][0] + x.fFull[ru][t][1]) / S0.best[t][0][0].F - 1);
                for (int dir = 0; dir < 2; ++dir) nSame[t][ru] += x.cell[ru][t][dir] == ref;
            }
            double m, e, r;
            ms(a, m, e, r);
            mB[t][ru] = m; eB[t][ru] = e; rmsB[t][ru] = r;
            ms(rel, m, e, r);
            mR[t][ru] = m; eR[t][ru] = e; rmsR[t][ru] = r;
            ms(pb, m, e, r);
            mP[t][ru] = m; eP[t][ru] = e;
            ms(ls, m, e, r);
            mLoss[t][ru] = m; eLoss[t][ru] = e;
        }
    printf("random splits done after %.0f s\n", sw.RealTime()); sw.Continue();

    // ============================================================================ tables
    gOut = fopen(Form("%s/%s_results.txt", kTabDir, kTag), "w");
    P("# %s_results.txt -- plots/eta_sim/%s.C (2026-11-12): optimisation of the eta-selection cuts with the simulation truth\n", kTag, kTag);
    P("# SIMULATION ONLY (no measured data read; no combinatorial background exists in simulation and none is added).\n");
    P("# Simulation: SMASH std + lep + eta owned rows (mode A), w (PDG) x physics weights of the eta decays, main variant\n");
    P("#   simphys::kPhysMain = %d (V3 direct decay + KLOE-2 bump + F_V conversion + background rules); %s\n", simphys::kPhysMain,
      OptLabel(o).c_str());
    P("# sigma = sum w / L_sim, L_sim = %.2f nb^-1; S, B = expected counts at L_data = %.1f nb^-1 (UNCONFIRMED, a scale only)\n", kLsim, kLd);
    P("# base 'pre' = isBest >= 0, vz > -500, PT3 bit, vz > -160 mm, cor_ok == 1; all cut variables from corrected momenta\n");
    P("# windows (M4, [lo, hi)): PEAK %.4f-%.4f, BUMP %.4f-%.4f GeV/c2\n", kPkLo, kPkHi, kBuLo, kBuHi);
    P("# targets: PEAK = classes 4+5 in the peak window; DIRECT = class 5 in the peak window (class 4 is background);\n");
    P("#          BUMP = classes 2+3 in the bump window; B = all other T1 + T2 rows of the window (variant '+T3T4': + templates 3 + 4)\n");
    P("# FOM F = S / sqrt(S + B); purity S / (S + B); efficiency S / S_pre. Errors: simulation statistics, event-group variance\n");
    P("#   (rows of one event group correlated), delta method for F. X17 rows (class 6) never enter.\n");
    P("# rows: %zu compact pre rows owned in mode A (X17 dropped %ld, w <= 0 dropped %ld); window rows %zu; true-angle lookups hit %lld miss %lld\n",
      gR.size(), nX17, nZero, wi.size(), simphys::AngHit(), simphys::AngMiss());
    P("# %s\n", dropNote.c_str());
    P("#   (w <= 0: eta-set rows whose only forced decay is X17 carry W = 2 BR(X17) = 0, see SimEngine.h)\n");
    P("\n## 0. Method\n");
    P("# Exhaustive grid search, exact at the grid values: every window row gets per variable the index of the grid interval of its\n");
    P("# value (edges = cut values, float comparison as the SimEngine presets: x > a for lower cuts, x < b for upper cuts, NaN fails\n");
    P("# every real cut); weights are filled into a 9-D histogram (7 quantities x 2 halves, double) and summed cumulatively along each\n");
    P("# axis (suffix sums for lower cuts and flags, prefix sums for upper cuts), so each cell holds S and B of one cut combination.\n");
    P("# Cells: %ld. Plateau = F >= %.2f F_best. Profile = best F over all other variables at a fixed value; slice = F with the other\n", gT.ncell, kPlat);
    P("# variables at the chosen point.\n");
    P("# Freeze rule: among cells with F >= %.2f F_best on the full sample AND on each half (relative to that half's best), the fewest\n", kPlat);
    P("# active cuts; among those the highest worst-neighbour F (each active threshold moved one grid step; flags have no neighbour),\n");
    P("# required >= %.2f F_best (else one more cut is allowed); ties by F.\n", kRobMin);
    P("# Split: half A = event groups with an even ordinal grp (per tree = per set; 'fold' column), half B = odd; each half x 2.\n");
    P("# Bias = F(optimising half) - F(other half) at the cells chosen on the optimising half (best cell, and the freeze rule applied\n");
    P("# to that half alone). Random splits: half = bit 0 of splitmix64((set << 40) ^ grp ^ (seed << 48)), seeds 1..%d.\n", K);
    P("# grids (none = no cut; flags off/on):\n");
    for (int dd = 0; dd < ND; ++dd) {
        P("#   %-20s %-6s:", gD[dd].key.c_str(), gD[dd].type == kLower ? "x > " : gD[dd].type == kUpper ? "x < " : "flag");
        for (int j = 0; j < gD[dd].n(); ++j) P(" %s", ValStr(dd, j).c_str());
        P("  [%s]\n", gD[dd].unit.c_str());
    }

    auto evLine = [&](const char* lab, const Ev& e, double fbest) {
        P("%-34s F = %8.3f +- %6.3f (%.3f of best) | S = %9.1f B = %10.1f | purity %.4f +- %.4f | eff %.4f +- %.4f | N_eff S %8.1f B %8.1f\n", lab,
          e.F, e.eF, fbest > 0 ? e.F / fbest : 0., e.S, e.B, e.pur, e.epur, e.eff, e.eeff, e.NeffS(), e.NeffB());
    };
    auto setLine = [&](const Ev& e) {
        P("%-34s", "   per set (rows / N_eff):");
        for (int s = 1; s <= 3; ++s) {
            const double ns = e.vSs[s] > 0 ? e.sSs[s] * e.sSs[s] / e.vSs[s] : 0, nbk = e.vBs[s] > 0 ? e.sBs[s] * e.sBs[s] / e.vBs[s] : 0;
            P(" %s: S %.0f / %.1f, B %.0f / %.1f, B share %.3f;", SetName(s), e.nSs[s], ns, e.nBs[s], nbk, e.sB > 0 ? e.sBs[s] / e.sB : 0.);
        }
        P("\n");
    };
    auto compLines = [&](const Ev& e) {
        P("   composition of the window rows passing (sigma [nb] +- err; rows; N_eff):\n");
        for (int g = 0; g < NG; ++g)
            P("     %-24s %10.5f +- %8.5f  rows %8.0f  N_eff %10.1f\n", GName(g), e.sG[g] / kLsim, std::sqrt(e.vG[g]) / kLsim, e.nG[g],
              e.vG[g] > 0 ? e.sG[g] * e.sG[g] / e.vG[g] : 0.);
    };

    P("\n## 1. Base 'pre' per target\n");
    for (int t = 0; t < NT; ++t) {
        P("\n# %s: %s\n", TName(t), TDesc(t));
        evLine("pre", evPre[t], S0.best[t][0][0].F);
        setLine(evPre[t]);
        compLines(evPre[t]);
    }

    P("\n## 2. Optimum per target (full sample, halves, variant +T3T4)\n");
    for (int t = 0; t < NT; ++t) {
        P("\n# %s: %s\n", TName(t), TDesc(t));
        P("best full sample : %s\n", CutStr(gT.Decode(S0.best[t][0][0].cell)).c_str());
        evLine("   best (full)", evBest[t], S0.best[t][0][0].F);
        setLine(evBest[t]);
        compLines(evBest[t]);
        P("best on half A   : %s  (F_A = %.3f, F_B there = %.3f)\n", CutStr(gT.Decode(R0.cell[0][t][0])).c_str(), R0.fOpt[0][t][0], R0.fOth[0][t][0]);
        P("best on half B   : %s  (F_B = %.3f, F_A there = %.3f)\n", CutStr(gT.Decode(R0.cell[0][t][1])).c_str(), R0.fOpt[0][t][1], R0.fOth[0][t][1]);
        P("variant +T3T4    : best %s, F = %.3f; frozen cut there %.3f of the variant best\n", CutStr(gT.Decode(S0.best[t][0][1].cell)).c_str(),
          S0.best[t][0][1].F, fVarAtFrzRel[t]);
        evLine("   best (variant +T3T4)", evBestV[t], S0.best[t][0][1].F);
        P("frozen (rule)    : %s\n", CutStr(gT.Decode(frz[t].cell)).c_str());
        P("   rule: %d active cuts, worst neighbour %.3f of best, plateau cells %ld (with this number of cuts %ld), on both half plateaus: %s\n",
          frz[t].cx, frz[t].rob, frz[t].nPlat, frz[t].nPlatK, frz[t].halvesOK ? "yes" : "NO (condition dropped)");
        evLine("   frozen (full)", evFrz[t], S0.best[t][0][0].F);
        setLine(evFrz[t]);
        compLines(evFrz[t]);
        evLine("   frozen (variant +T3T4)", evFrzV[t], S0.best[t][0][1].F);
        P("   frozen on the halves: F_A %.3f (%.3f of the A best), F_B %.3f (%.3f of the B best)\n", fHalfAtFrz[t][0], fHalfAtFrz[t][0] / fHalfBest[t][0],
          fHalfAtFrz[t][1], fHalfAtFrz[t][1] / fHalfBest[t][1]);
        P("purity rule (alternative, not frozen): %s\n", CutStr(gT.Decode(pur[t].cell)).c_str());
        P("   rule P: highest purity on the same plateau (full + both halves) with <= %d cuts and worst neighbour >= %.2f: %d cuts, worst neighbour %.3f\n",
          kMaxCutsP, kRobMin, pur[t].cx, pur[t].rob);
        evLine("   purity rule (full)", evPur[t], S0.best[t][0][0].F);
        setLine(evPur[t]);
        compLines(evPur[t]);
        P("   purity rule on the halves: F_A %.3f (%.3f of the A best), F_B %.3f (%.3f of the B best)\n", fHalfAtPur[t][0], fHalfAtPur[t][0] / fHalfBest[t][0],
          fHalfAtPur[t][1], fHalfAtPur[t][1] / fHalfBest[t][1]);
    }

    P("\n## 3. Plateau per variable (F / F_best of the full sample; '*' = on the plateau F >= %.2f)\n", kPlat);
    P("# profile = best over all other variables; slice(frozen) and slice(best) = the other variables at the frozen / best point\n");
    P("# profA, profB = profiles of the halves relative to their own best; profV = profile of the variant +T3T4 (to its best)\n");
    for (int t = 0; t < NT; ++t) {
        P("\n# %s (F_best = %.3f)\n", TName(t), S0.best[t][0][0].F);
        const double fb = S0.best[t][0][0].F;
        for (int dd = 0; dd < ND; ++dd) {
            P("%-20s %-9s", gD[dd].key.c_str(), "value");
            for (int j = 0; j < gD[dd].n(); ++j) P(" %7s", ValStr(dd, j).c_str());
            P("\n");
            auto row = [&](const char* lab, const std::vector<double>& v, double norm, bool mark) {
                P("%-20s %-9s", "", lab);
                for (double x : v) P(" %6.3f%s", x / norm, mark && x >= kPlat * norm ? "*" : " ");
                P("\n");
            };
            row("profile", S0.prof[t][0][dd], fb, true);
            row("profA", S0.prof[t][1][dd], S0.best[t][1][0].F, true);
            row("profB", S0.prof[t][2][dd], S0.best[t][2][0].F, true);
            row("profV", S0.profV[t][dd], S0.best[t][0][1].F, true);
            row("sl.frozen", slF[t][dd], fb, true);
            row("sl.best", slB[t][dd], fb, true);
            // ranges
            std::string pr, sr;
            for (int j = 0; j < gD[dd].n(); ++j) {
                if (S0.prof[t][0][dd][j] >= kPlat * fb) pr += (pr.empty() ? "" : " ") + ValStr(dd, j);
                if (slF[t][dd][j] >= kPlat * fb) sr += (sr.empty() ? "" : " ") + ValStr(dd, j);
            }
            P("%-20s plateau (profile): {%s}; plateau (slice at frozen): {%s}; frozen value %s\n", "", pr.c_str(), sr.c_str(),
              ValStr(dd, gT.Decode(frz[t].cell).j[dd]).c_str());
        }
    }

    P("\n## 4. Named presets (full sample)\n");
    for (int t = 0; t < NT; ++t) {
        P("\n# %s (F_best = %.3f; variant +T3T4 best %.3f)\n", TName(t), S0.best[t][0][0].F, S0.best[t][0][1].F);
        for (int p = 0; p < NPS; ++p) {
            evLine(Form("%s", PsName(p)), evPs[t][p], S0.best[t][0][0].F);
            P("%-34s variant +T3T4: F = %.3f +- %.3f (%.3f of the variant best), purity %.4f\n", "", evPsV[t][p].F, evPsV[t][p].eF,
              evPsV[t][p].F / S0.best[t][0][1].F, evPsV[t][p].pur);
        }
        evLine("best", evBest[t], S0.best[t][0][0].F);
        evLine("frozen", evFrz[t], S0.best[t][0][0].F);
    }

    P("\n## 5. Split-half check\n");
    P("# deterministic split (fold = grp parity per set). F with errors (simulation statistics of that half). Per direction the\n");
    P("# difference F_opt - F_other contains the half-to-half offset at fixed cuts (e.g. 'pre' F_A %.3f / F_B %.3f for PEAK);\n",
      fHalfAtFrzPre[0][0], fHalfAtFrzPre[0][1]);
    P("# the mean of the two directions cancels it to first order and is the bias estimate ('2-dir mean').\n");
    for (int t = 0; t < NT; ++t) {
        P("\n# %s (pre on the halves: F_A %.3f, F_B %.3f)\n", TName(t), fHalfAtFrzPre[t][0], fHalfAtFrzPre[t][1]);
        for (int ru = 0; ru < NRU; ++ru) {
            double sumAbs = 0, sumRel = 0;
            for (int dir = 0; dir < 2; ++dir) {
                const Ev& a = evSplit[t][ru][dir][0];
                const Ev& b = evSplit[t][ru][dir][1];
                sumAbs += a.F - b.F; sumRel += a.F / b.F - 1;
                P("%-11s optimise %c evaluate %c: %-62s F_opt %8.3f +- %6.3f, F_other %8.3f +- %6.3f, diff %+7.3f (%+6.2f %%), F_full %8.3f\n",
                  RuName(ru), 'A' + dir, 'B' - dir, CutStr(gT.Decode(R0.cell[ru][t][dir])).c_str(), a.F, a.eF, b.F, b.eF, a.F - b.F,
                  100 * (a.F / b.F - 1), fFullAtHalfPick[t][ru][dir]);
            }
            P("%-11s 2-dir mean bias %+.3f (%+.2f %%)\n", RuName(ru), 0.5 * sumAbs, 50 * sumRel);
        }
    }
    P("\n# %d random splits (two directions each): bias = F(optimising half) - F(other half) at the cells chosen on the optimising half\n", K);
    P("# per split the mean of both directions; mean over splits +- error of the mean (rms / sqrt(K)); rms over splits.\n");
    P("# purity bias = purity(optimising half) - purity(other half), same averaging. same pick = how often (of 2K) the half pick equals the\n");
    P("# full-sample pick of the same rule. loss = F_full(half pick) / F_full(best) - 1 (indicative only: the half is part of the full sample).\n");
    for (int t = 0; t < NT; ++t)
        for (int ru = 0; ru < NRU; ++ru)
            P("%-7s %-11s: bias %+.3f +- %.3f (rms %.3f) in F units; relative %+.2f +- %.2f %% (rms %.2f %%); purity bias %+.4f +- %.4f; same pick %d/%d; loss %+.2f +- %.2f %%\n",
              TName(t), RuName(ru), mB[t][ru], eB[t][ru], rmsB[t][ru], 100 * mR[t][ru], 100 * eR[t][ru], 100 * rmsR[t][ru], mP[t][ru], eP[t][ru],
              nSame[t][ru], 2 * K, 100 * mLoss[t][ru], 100 * eLoss[t][ru]);
    P("\n# per random split: target, rule, F_opt(A), F_oth(B), F_opt(B), F_oth(A), cell chosen on A, cell chosen on B\n");
    for (size_t k = 0; k < RS.size(); ++k)
        for (int t = 0; t < NT; ++t)
            for (int ru = 0; ru < NRU; ++ru)
                P("split %3zu %-6s %-11s %8.3f %8.3f %8.3f %8.3f  A: %-58s B: %s\n", k + 1, TName(t), RuName(ru), RS[k].fOpt[ru][t][0],
                  RS[k].fOth[ru][t][0], RS[k].fOpt[ru][t][1], RS[k].fOth[ru][t][1], CutStr(gT.Decode(RS[k].cell[ru][t][0])).c_str(),
                  CutStr(gT.Decode(RS[k].cell[ru][t][1])).c_str());
    P("\n## 5b. Differences of F with the event-group error, correlated rows included (full sample; L_data UNCONFIRMED)\n");
    P("# the same rows enter both selections, so the error of a difference is much smaller than the errors of the two F values\n");
    P("%-8s %22s %22s %22s %22s\n", "target", "frozen - pre", "best - pre", "pions-in-vtx - pre", "frozen - pions-in-vtx");
    for (int t = 0; t < NT; ++t)
        P("%-8s %12.3f +- %6.3f %12.3f +- %6.3f %12.3f +- %6.3f %12.3f +- %6.3f\n", TName(t), dFrzPre[t].d, dFrzPre[t].e, dBestPre[t].d,
          dBestPre[t].e, dVtxPre[t].d, dVtxPre[t].e, dFrzVtx[t].d, dFrzVtx[t].e);
    P("# note: 'best - pre' is biased upwards by the choice of the best cell (split-half bias, section 5)\n");
    P("\n## 6. Checks (%zu, failed %d)\n", checks.size(), nCheckBad);
    for (const std::string& s : checks) P("%s\n", s.c_str());
    P("\n# run time %.0f s\n", sw.RealTime()); sw.Continue();
    fclose(gOut);
    gOut = nullptr;

    // ---------------------------------------------------------------- frozen file
    gOut = fopen(Form("%s/%s_frozen.txt", kTabDir, kTag), "w");
    P("# %s_frozen.txt -- frozen eta-selection cuts from the simulation-truth optimisation (plots/eta_sim/%s.C, 2026-11-12)\n", kTag, kTag);
    P("# SIMULATION ONLY (no CB, physics weights V3 main); FOM = S/sqrt(S+B) at L_data = 6080.5 nb^-1 (UNCONFIRMED, scale only)\n");
    P("# base: pre = isBest >= 0, vz > -500, PT3 bit, vz > -160 mm, corrected momenta, cor_ok == 1\n");
    P("# cuts: *_min means value > cut, *_max means value < cut; none = not applied; flags 0 = off, 1 = on\n");
    P("#   pions_in_vertex_fit = 1: pip_isUsedVertex == 1 && pim_isUsedVertex == 1\n");
    P("#   dz_pi_max = 40: max(|pip_z - vz|, |pim_z - vz|) < 40 mm; units: deg, GeV/c2, MeV/c\n");
    P("# rule: fewest active cuts on the 5%% plateau of the full sample and of both halves, highest worst-neighbour FOM\n");
    P("# NOTE (2026-11-13, review M4): dz_pi needs the pion track z (pip_z, pim_z). The data compact rows (ScanEngine.h, scan::Row) do not\n");
    P("#   have them, so PEAK and DIRECT cannot yet be applied to data. Data-ready alternative on the plateau: pre + both pions in the\n");
    P("#   vertex fit (F/F_best = %.3f PEAK, %.3f DIRECT, %.3f BUMP).\n", evPs[0][kPrePiVtx].F / S0.best[0][0][0].F,
      evPs[1][kPrePiVtx].F / S0.best[1][0][0].F, evPs[2][kPrePiVtx].F / S0.best[2][0][0].F);
    P("# NOTE (review M5): differences with the correlated event-group error, F(frozen) - F(pre): PEAK %+.2f +- %.2f, DIRECT %+.2f +- %.2f;\n",
      dFrzPre[0].d, dFrzPre[0].e, dFrzPre[1].d, dFrzPre[1].e);
    P("#   F(frozen) - F(pre + pions in vertex fit): PEAK %+.2f +- %.2f, DIRECT %+.2f +- %.2f, BUMP %+.2f +- %.2f. The numbers below are\n",
      dFrzVtx[0].d, dFrzVtx[0].e, dFrzVtx[1].d, dFrzVtx[1].e, dFrzVtx[2].d, dFrzVtx[2].e);
    P("#   full-sample (half A + half B) values: do not give them to the analyst of a blind test (review H1).\n\n");
    P("[windows]\npeak_m4_min = %.4f\npeak_m4_max = %.4f\nbump_m4_min = %.4f\nbump_m4_max = %.4f\n", kPkLo, kPkHi, kBuLo, kBuHi);
    for (int t = 0; t < NT; ++t) {
        const Cut c = gT.Decode(frz[t].cell);
        P("\n[%s]\n", TName(t));
        P("target = %s\n", TDesc(t));
        P("oa_min = %s\n", ValStr(dOA, c.j[dOA]).c_str());
        P("alpha_lab_max = %s\n", ValStr(dLab, c.j[dLab]).c_str());
        P("mpipi_min = %s\n", ValStr(dMlo, c.j[dMlo]).c_str());
        P("mpipi_max = %s\n", ValStr(dMhi, c.j[dMhi]).c_str());
        P("alpha_eta_min = %s\n", ValStr(dEta, c.j[dEta]).c_str());
        P("pions_in_vertex_fit = %d\n", c.j[dVtx]);
        P("dz_pi_max = %s\n", c.j[dDz] ? "40" : "none");
        P("mm4_max = %s\n", ValStr(dMM, c.j[dMM]).c_str());
        P("pt4_min = %s\n", ValStr(dPt, c.j[dPt]).c_str());
        P("fom = %.3f +- %.3f\nfom_over_best = %.3f\npurity = %.4f +- %.4f\nefficiency = %.4f +- %.4f\n", evFrz[t].F, evFrz[t].eF,
          evFrz[t].F / S0.best[t][0][0].F, evFrz[t].pur, evFrz[t].epur, evFrz[t].eff, evFrz[t].eeff);
        P("best_point = %s (fom %.3f +- %.3f, purity %.4f)\n", CutStr(gT.Decode(S0.best[t][0][0].cell)).c_str(), S0.best[t][0][0].F,
          evBest[t].eF, evBest[t].pur);
        P("optimisation_bias_frozen_rule = %+.2f +- %.2f %% (best cell %+.2f +- %.2f %%; %d random split halves)\n", 100 * mR[t][1], 100 * eR[t][1],
          100 * mR[t][0], 100 * eR[t][0], K);
        P("# alternative, not frozen (highest purity on the same plateau, <= %d cuts): %s; fom %.3f (%.3f of best), purity %.4f, efficiency %.4f\n",
          kMaxCutsP, CutStr(gT.Decode(pur[t].cell)).c_str(), evPur[t].F, evPur[t].F / S0.best[t][0][0].F, evPur[t].pur, evPur[t].eff);
        P("justification = %d active cut%s, the fewest on the 5%% plateau (full sample %.3f, half A %.3f, half B %.3f of their best); "
          "every one-step change of a threshold keeps >= %.3f of the best FOM\n",
          frz[t].cx, frz[t].cx == 1 ? "" : "s", evFrz[t].F / S0.best[t][0][0].F, fHalfAtFrz[t][0] / fHalfBest[t][0], fHalfAtFrz[t][1] / fHalfBest[t][1],
          frz[t].rob);
    }
    fclose(gOut);
    gOut = nullptr;

    // ============================================================================ figures
    const int cFull = kBlack, cA = kAzure + 2, cB = kOrange + 7, cSl = kGreen + 2, cFrz = kMagenta + 1;
    // ---------- 1. profiles per target
    for (int t = 0; t < NT; ++t) {
        TCanvas* c = new TCanvas(Form("cprof%d", t), "", 2100, 2000);
        TPad* top = new TPad("top", "", 0, 0.86, 1, 1);
        top->Draw();
        top->cd();
        const double ts = 0.105;   // 29 px
        Txt(0.01, 0.84, ts, Form("%s: F = S/#sqrt{S+B} relative to the best (F_{best} = %.1f, %s)", TLab(t), S0.best[t][0][0].F, kLdNote));
        Txt(0.01, 0.66, ts, Form("best: %s;   frozen: %s", CutStr(gT.Decode(S0.best[t][0][0].cell)).c_str(), CutStr(gT.Decode(frz[t].cell)).c_str()));
        Txt(0.01, 0.48, ts, "SMASH std + lep + eta (mode A), w (PDG) #times physics weights V3; base 'pre'; B = other T1+T2 rows of the window");
        {
            TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.40);
            lg->SetTextSize(ts); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetNColumns(4);
            TGraph* d1 = new TGraph(); d1->SetLineColor(cFull); d1->SetMarkerColor(cFull); d1->SetMarkerStyle(20); d1->SetLineWidth(3);
            TGraph* d2 = new TGraph(); d2->SetLineColor(cA); d2->SetMarkerColor(cA); d2->SetMarkerStyle(24); d2->SetLineStyle(2); d2->SetLineWidth(3);
            TGraph* d3 = new TGraph(); d3->SetLineColor(cB); d3->SetMarkerColor(cB); d3->SetMarkerStyle(25); d3->SetLineStyle(2); d3->SetLineWidth(3);
            TGraph* d4 = new TGraph(); d4->SetLineColor(cSl); d4->SetMarkerColor(cSl); d4->SetMarkerStyle(21); d4->SetLineWidth(3);
            TMarker* d5 = new TMarker(0, 0, 23); d5->SetMarkerColor(cFrz); d5->SetMarkerSize(2.5);
            TBox* d6 = new TBox(); d6->SetFillColor(kGray);
            lg->AddEntry(d1, "profile, full sample", "lp");
            lg->AddEntry(d4, "slice through frozen", "lp");
            lg->AddEntry(d2, "profile, half A", "lp");
            lg->AddEntry(d3, "profile, half B", "lp");
            lg->AddEntry(d5, "frozen value", "p");
            lg->AddEntry(d6, "plateau #geq 0.95", "f");
            TBox* d7 = new TBox(); d7->SetFillColor(kGray + 3); d7->SetFillStyle(3345); d7->SetLineColor(kGray + 3);
            lg->AddEntry(d7, "#pm 1#sigma (sim. stat.) of the full profile", "f");
            lg->Draw();
        }
        c->cd();
        TPad* bot = new TPad("bot", "", 0, 0, 1, 0.86);
        bot->Draw();
        bot->cd();
        bot->Divide(3, 3, 0.002, 0.002);
        const double fb = S0.best[t][0][0].F;
        const Cut kf = gT.Decode(frz[t].cell);
        for (int dd = 0; dd < ND; ++dd) {
            bot->cd(dd + 1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.27); gPad->SetTopMargin(0.05);
            const int n = gD[dd].n();
            double ymin = 1.0;
            for (int j = 0; j < n; ++j) {
                ymin = std::min(ymin, S0.prof[t][0][dd][j] / fb);
                ymin = std::min(ymin, slF[t][dd][j] / fb);
            }
            ymin = std::max(0., std::floor((ymin - 0.05) * 10.) / 10.);
            TH1D* fr = new TH1D(Form("fr%d_%d", t, dd), "", n, 0, n);
            fr->SetDirectory(nullptr);
            for (int j = 0; j < n; ++j) fr->GetXaxis()->SetBinLabel(j + 1, ValStr(dd, j).c_str());
            fr->SetMinimum(ymin); fr->SetMaximum(1.12);
            fr->GetXaxis()->SetTitle(Form("%s%s", gD[dd].title.c_str(), gD[dd].unit.empty() ? "" : Form(" [%s]", gD[dd].unit.c_str())));
            fr->GetYaxis()->SetTitle("F / F_{best}");
            AxisStyle(fr, 0.065, 0.065, 1.75, 1.05);
            fr->GetXaxis()->SetLabelSize(0.07);
            fr->GetXaxis()->LabelsOption("v");
            fr->GetYaxis()->SetNdivisions(505);
            fr->Draw("AXIS");
            TBox* bx = new TBox(0, kPlat, n, 1.0);
            bx->SetFillColor(kGray);
            bx->Draw();
            TLine* l1 = new TLine(0, 1, n, 1); l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
            auto gr = [&](const std::vector<double>& v, double norm, int col, int sty, int mk, double off) {
                TGraph* g = new TGraph(n);
                for (int j = 0; j < n; ++j) g->SetPoint(j, j + 0.5 + off, v[j] / norm);
                g->SetLineColor(col); g->SetMarkerColor(col); g->SetLineStyle(sty); g->SetLineWidth(3); g->SetMarkerStyle(mk); g->SetMarkerSize(1.8);
                g->Draw("LP");
                return g;
            };
            gr(S0.prof[t][1][dd], S0.best[t][1][0].F, cA, 2, 24, -0.08);
            gr(S0.prof[t][2][dd], S0.best[t][2][0].F, cB, 2, 25, 0.08);
            {
                TGraphErrors* ge = new TGraphErrors(n);
                for (int j = 0; j < n; ++j) { ge->SetPoint(j, j + 0.5, S0.prof[t][0][dd][j] / fb); ge->SetPointError(j, 0.38, profE[t][dd][j] / fb); }
                ge->SetFillColor(kGray + 3); ge->SetFillStyle(3345); ge->SetLineColor(kGray + 3);
                ge->Draw("2 SAME");
            }
            gr(slF[t][dd], fb, cSl, 1, 21, 0.);
            gr(S0.prof[t][0][dd], fb, cFull, 1, 20, 0.);
            TMarker* m = new TMarker(kf.j[dd] + 0.5, 1.07, 23);
            m->SetMarkerColor(cFrz); m->SetMarkerSize(3.0); m->Draw();
            fr->Draw("AXIS SAME");
        }
        Save(c, Form("%s_profile_%s", kTag, t == 0 ? "peak" : t == 1 ? "direct" : "bump"));
    }

    // ---------- 2. presets, best and frozen: FOM ratio, purity, efficiency
    {
        TCanvas* c = new TCanvas("cps", "", 2100, 1900);
        TPad* top = new TPad("top", "", 0, 0.88, 1, 1);
        top->Draw(); top->cd();
        const double ts = 0.125;   // 28.5 px
        Txt(0.01, 0.78, ts, Form("Named selections against the optimum; simulation only, physics weights V3; F = S/#sqrt{S+B} at %s", kLdNote));
        Txt(0.01, 0.56, ts, "filled: B = other T1+T2 rows of the window; open grey: variant B + T3 + T4; error bars: simulation statistics (event groups)");
        for (int t = 0; t < NT; ++t)
            Txt((t + 0.5) / 3., 0.12, 0.15, Form("%s (F_{best} = %.1f)", TName(t), S0.best[t][0][0].F), 21);
        c->cd();
        TPad* bot = new TPad("bot", "", 0, 0, 1, 0.88);
        bot->Draw(); bot->cd();
        bot->Divide(3, 3, 0.002, 0.002);
        const int nsel = NPS + 2;
        const char* slab[NPS + 2] = {"pre", "+OA>2", "+#pi vtx", "slides", "old code", "peak nom.", "bump nom.", "best", "frozen"};
        for (int row = 0; row < 3; ++row)
            for (int t = 0; t < NT; ++t) {
                bot->cd(row * 3 + t + 1);
                gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.31); gPad->SetTopMargin(0.04);
                TH1D* fr = new TH1D(Form("fps%d_%d", row, t), "", nsel, 0, nsel);
                fr->SetDirectory(nullptr);
                for (int k = 0; k < nsel; ++k) fr->GetXaxis()->SetBinLabel(k + 1, slab[k]);
                AxisStyle(fr, 0.065, 0.068, 1.0, 1.15);
                fr->GetXaxis()->SetLabelSize(0.075);
                fr->GetXaxis()->LabelsOption("v");
                fr->GetYaxis()->SetTitle(row == 0 ? "F / F_{best}" : row == 1 ? "purity S/(S+B)" : "efficiency S/S_{pre}");
                std::vector<const Ev*> ev, evv;
                for (int p = 0; p < NPS; ++p) { ev.push_back(&evPs[t][p]); evv.push_back(&evPsV[t][p]); }
                ev.push_back(&evBest[t]); ev.push_back(&evFrz[t]);
                evv.push_back(&evBestV[t]); evv.push_back(&evFrzV[t]);
                double ymax = 0, ymin = 1e30;
                TGraphErrors* g = new TGraphErrors(nsel);
                TGraphErrors* gv = new TGraphErrors(nsel);
                for (int k = 0; k < nsel; ++k) {
                    double y, e, yv, ev2;
                    if (row == 0) { y = ev[k]->F / S0.best[t][0][0].F; e = ev[k]->eF / S0.best[t][0][0].F; yv = evv[k]->F / S0.best[t][0][1].F; ev2 = evv[k]->eF / S0.best[t][0][1].F; }
                    else if (row == 1) { y = ev[k]->pur; e = ev[k]->epur; yv = evv[k]->pur; ev2 = evv[k]->epur; }
                    else { y = ev[k]->eff; e = ev[k]->eeff; yv = evv[k]->eff; ev2 = evv[k]->eeff; }
                    g->SetPoint(k, k + 0.4, y); g->SetPointError(k, 0, e);
                    gv->SetPoint(k, k + 0.62, yv); gv->SetPointError(k, 0, ev2);
                    ymax = std::max(ymax, std::max(y + e, yv + ev2));
                    ymin = std::min(ymin, std::min(y - e, yv - ev2));
                }
                fr->SetMinimum(row == 0 ? std::max(0., std::floor((ymin - 0.05) * 10.) / 10.) : 0.);
                fr->SetMaximum(row == 0 ? 1.08 : ymax * 1.22);
                fr->Draw("AXIS");
                if (row == 0) {
                    TBox* bx = new TBox(0, kPlat, nsel, 1.0); bx->SetFillColor(kGray); bx->Draw();
                }
                g->SetMarkerStyle(20); g->SetMarkerSize(2.2); g->SetLineWidth(3);
                gv->SetMarkerStyle(24); gv->SetMarkerSize(2.2); gv->SetLineWidth(2); gv->SetMarkerColor(kGray + 2); gv->SetLineColor(kGray + 2);
                gv->Draw("P");
                g->Draw("P");
                TMarker* mf = new TMarker(nsel - 1 + 0.4, g->GetY()[nsel - 1], 20); mf->SetMarkerColor(cFrz); mf->SetMarkerSize(2.4); mf->Draw();
                fr->Draw("AXIS SAME");
            }
        Save(c, Form("%s_presets", kTag));
    }

    // ---------- 3. split-half check
    {
        TCanvas* c = new TCanvas("csp", "", 2100, 1750);
        TPad* top = new TPad("top", "", 0, 0.86, 1, 1);
        top->Draw(); top->cd();
        const double ts = 0.105;
        Txt(0.01, 0.82, ts, Form("Split-half check of the optimisation bias; each half #times 2; F = S/#sqrt{S+B} at %s", kLdNote));
        Txt(0.01, 0.62, ts, Form("top: deterministic split (event-group parity), errors = simulation statistics; bottom: %d random splits, mean of both directions", K));
        const int rcol[NRU] = {kAzure + 2, kOrange + 7, kGreen + 2};
        TLegend* lg = new TLegend(0.0, 0.25, 1.0, 0.47);
        lg->SetNColumns(4); lg->SetTextSize(ts); lg->SetBorderSize(0); lg->SetFillStyle(0);
        TGraphErrors* lgo = new TGraphErrors(); lgo->SetMarkerStyle(20); lgo->SetMarkerSize(2.2); lgo->SetMarkerColor(kBlack); lgo->SetLineWidth(3);
        TGraphErrors* lgx = new TGraphErrors(); lgx->SetMarkerStyle(24); lgx->SetMarkerSize(2.2); lgx->SetMarkerColor(kGray + 2); lgx->SetLineColor(kGray + 2); lgx->SetLineWidth(3);
        TGraph* lgf = new TGraph(); lgf->SetMarkerStyle(29); lgf->SetMarkerSize(2.6);
        TLine* lgb = new TLine(); lgb->SetLineStyle(2); lgb->SetLineWidth(2);
        lg->AddEntry(lgo, "F on the optimising half", "pe");
        lg->AddEntry(lgx, "F on the other half", "pe");
        lg->AddEntry(lgf, "F of the full sample, same cuts", "p");
        lg->AddEntry(lgb, "best F, full sample", "l");
        lg->Draw();
        for (int t = 0; t < NT; ++t) Txt((t + 0.5) / 3., 0.04, 0.13, TName(t), 21);
        c->cd();
        TPad* bot = new TPad("bot", "", 0, 0, 1, 0.86);
        bot->Draw(); bot->cd();
        bot->Divide(3, 2, 0.002, 0.002);
        for (int t = 0; t < NT; ++t) {
            bot->cd(t + 1);
            gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.15); gPad->SetTopMargin(0.12);
            const int nx = 2 * NRU;
            TH1D* fr = new TH1D(Form("fsp%d", t), "", nx, 0, nx);
            fr->SetDirectory(nullptr);
            for (int k = 0; k < nx; ++k) fr->GetXaxis()->SetBinLabel(k + 1, k % 2 ? "B#rightarrowA" : "A#rightarrowB");
            AxisStyle(fr, 0.06, 0.065, 1.0, 1.45);
            fr->GetXaxis()->SetLabelSize(0.07);
            fr->GetYaxis()->SetTitle("F = S/#sqrt{S+B}");
            TGraphErrors* go = new TGraphErrors(nx);
            TGraphErrors* gx = new TGraphErrors(nx);
            TGraph* gf = new TGraph(nx);
            double lo = 1e30, hi = 0;
            for (int k = 0; k < nx; ++k) {
                const int ru = k / 2, dir = k % 2;
                const Ev& a = evSplit[t][ru][dir][0];
                const Ev& b = evSplit[t][ru][dir][1];
                go->SetPoint(k, k + 0.3, a.F); go->SetPointError(k, 0, a.eF);
                gx->SetPoint(k, k + 0.7, b.F); gx->SetPointError(k, 0, b.eF);
                gf->SetPoint(k, k + 0.5, fFullAtHalfPick[t][ru][dir]);
                lo = std::min(lo, std::min(a.F - a.eF, b.F - b.eF)); hi = std::max(hi, std::max(a.F + a.eF, b.F + b.eF));
            }
            lo = std::min(lo, S0.best[t][0][0].F); hi = std::max(hi, S0.best[t][0][0].F);
            const double pad = 0.12 * (hi - lo);
            fr->SetMinimum(lo - pad); fr->SetMaximum(hi + 2.2 * pad);
            fr->Draw("AXIS");
            TLine* lb = new TLine(0, S0.best[t][0][0].F, nx, S0.best[t][0][0].F); lb->SetLineStyle(2); lb->SetLineWidth(2); lb->Draw();
            for (int r = 1; r < NRU; ++r) { TLine* sep = new TLine(2 * r, lo - pad, 2 * r, hi + 2.2 * pad); sep->SetLineColor(kGray + 1); sep->Draw(); }
            const char* rshort[NRU] = {"best cell", "freeze", "purity"};
            for (int r = 0; r < NRU; ++r) Txt(2 * r + 1, hi + 1.25 * pad, 0.055, rshort[r], 21, rcol[r], false);
            go->SetMarkerStyle(20); go->SetMarkerSize(2.2); go->SetLineWidth(3);
            gx->SetMarkerStyle(24); gx->SetMarkerSize(2.2); gx->SetMarkerColor(kGray + 2); gx->SetLineColor(kGray + 2); gx->SetLineWidth(3);
            gf->SetMarkerStyle(29); gf->SetMarkerSize(2.6); gf->SetMarkerColor(kBlack);
            gx->Draw("P"); go->Draw("P"); gf->Draw("P");
            fr->Draw("AXIS SAME");
            bot->cd(t + 4);
            gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.06);
            double amin = 0, amax = 0;
            for (int ru = 0; ru < NRU; ++ru) for (double x : relAvg[t][ru]) { amin = std::min(amin, 100 * x); amax = std::max(amax, 100 * x); }
            const double span = std::max(1., amax - amin);
            std::vector<TH1D*> hh(NRU);
            double hmax = 0;
            for (int ru = 0; ru < NRU; ++ru) {
                hh[ru] = new TH1D(Form("hb%d_%d", ru, t), "", 30, amin - 0.08 * span, amax + 0.08 * span);
                hh[ru]->SetDirectory(nullptr);
                for (double x : relAvg[t][ru]) hh[ru]->Fill(100 * x);
                hh[ru]->SetLineColor(rcol[ru]); hh[ru]->SetLineWidth(4 - ru % 2); hh[ru]->SetFillStyle(0);
                hmax = std::max(hmax, hh[ru]->GetMaximum());
            }
            hh[0]->SetFillColor(kAzure - 9); hh[0]->SetFillStyle(1001);
            AxisStyle(hh[0], 0.06, 0.065, 1.15, 1.45);
            hh[0]->GetXaxis()->SetTitle("bias (F_{opt} - F_{other}) / F_{other} [%]");
            hh[0]->GetYaxis()->SetTitle("random splits");
            hh[0]->SetMaximum(1.75 * hmax);
            hh[0]->Draw("HIST");
            for (int ru = 1; ru < NRU; ++ru) hh[ru]->Draw("HIST SAME");
            TLine* z = new TLine(0, 0, 0, 1.75 * hmax); z->SetLineStyle(2); z->Draw();
            for (int ru = 0; ru < NRU; ++ru)
                Txt(0.30, 0.87 - 0.075 * ru, 0.055, Form("%s: %+.2f #pm %.2f %%", RuName(ru), 100 * mR[t][ru], 100 * eR[t][ru]), 11, rcol[ru]);
        }
        Save(c, Form("%s_split", kTag));
    }

    // ---------- 4. M4 spectra at pre and the frozen cuts
    {
        TCanvas* c = new TCanvas("cm4", "", 2100, 1900);
        TPad* top = new TPad("top", "", 0, 0.88, 1, 1);
        top->Draw(); top->cd();
        Txt(0.01, 0.76, 0.125, "Simulated d#sigma/dM_{4} per group at 'pre' and at the frozen cuts (SMASH std + lep + eta, mode A, w #times physics weights V3; no CB)");
        Txt(0.01, 0.56, 0.125, "warm: signal classes 2-5; cold / grey: physical background groups; dashed: T3 fake lepton, T4 different sources (not in the sum)");
        TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.48);
        lg->SetNColumns(6); lg->SetTextSize(0.125); lg->SetBorderSize(0); lg->SetFillStyle(0);
        c->cd();
        TPad* bot = new TPad("bot", "", 0, 0, 1, 0.88);
        bot->Draw(); bot->cd();
        bot->Divide(2, 2, 0.002, 0.002);
        const int gcol[NG] = {kOrange - 3, kOrange + 7, kRed, kRed + 3, kViolet + 1, kAzure + 2, kCyan + 2, kBlue + 2, kGray + 2, kTeal + 3, kGray + 1};
        const int gsty[NG] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};
        const char* ttl[4] = {"pre", "PEAK frozen", "DIRECT frozen", "BUMP frozen"};
        for (int k = 0; k < 4; ++k) {
            bot->cd(k + 1);
            gPad->SetLogy(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.15);
            TH1D* hs = new TH1D(Form("hsum%d", k), "", nb, mlo, mhi);
            hs->SetDirectory(nullptr);
            double ymax = 0, ymin = 1e30;
            for (int b = 0; b < nb; ++b) {
                hs->SetBinContent(b + 1, specSum[k][b] / (kLsim * bw));
                hs->SetBinError(b + 1, std::sqrt(specSumV[k][b]) / (kLsim * bw));
                ymax = std::max(ymax, hs->GetBinContent(b + 1));
            }
            std::vector<TH1D*> hg(NG);
            for (int g = 0; g < NG; ++g) {
                hg[g] = new TH1D(Form("hg%d_%d", k, g), "", nb, mlo, mhi);
                hg[g]->SetDirectory(nullptr);
                for (int b = 0; b < nb; ++b) {
                    const double y = spec[k][g][b] / (kLsim * bw);
                    hg[g]->SetBinContent(b + 1, y);
                    if (y > 0) ymin = std::min(ymin, y);
                }
                hg[g]->SetLineColor(gcol[g]); hg[g]->SetLineStyle(gsty[g]); hg[g]->SetLineWidth(3);
            }
            hs->SetMaximum(ymax * 8); hs->SetMinimum(std::max(ymax * 2e-5, ymin * 0.5));
            AxisStyle(hs, 0.055, 0.06, 1.05, 1.15);
            hs->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            hs->GetYaxis()->SetTitle("d#sigma/dM_{4} [nb/(GeV/c^{2})]");
            hs->SetLineColor(kBlack); hs->SetLineWidth(4); hs->SetFillColor(kGray); hs->SetMarkerSize(0);
            hs->Draw("E2");
            TH1D* hl = (TH1D*)hs->Clone(Form("hsl%d", k)); hl->SetFillStyle(0); hl->Draw("HIST SAME");
            for (int g = 0; g < NG; ++g) hg[g]->Draw("HIST SAME");
            for (double x : {kBuLo, kBuHi, kPkLo, kPkHi}) {
                TLine* l = new TLine(x, hs->GetMinimum(), x, ymax * 3); l->SetLineStyle(3); l->SetLineWidth(2); l->SetLineColor(kGray + 3); l->Draw();
            }
            Txt(kBuLo - 0.004, ymax * 2.2, 0.045, "bump #rightarrow", 32, kGray + 3, false);
            Txt(kPkHi + 0.004, ymax * 2.2, 0.045, "#leftarrow peak", 12, kGray + 3, false);
            Txt(0.15, 0.935, 0.05, ttl[k]);
            Txt(0.15, 0.875, 0.045, CutStr(specCut[k]).c_str());
            hs->Draw("AXIS SAME");
            if (k == 0) {
                lg->AddEntry(hl, "T1+T2 sum", "l");
                for (int g = 0; g < NG; ++g) lg->AddEntry(hg[g], GShort(g), "l");
            }
        }
        top->cd();
        lg->Draw();
        Save(c, Form("%s_m4", kTag));
    }

    // ---------- 4b. window composition by group (2026-11-13: T3 and T4 drawn separately, not stacked; plain names)
    {
        TCanvas* c = new TCanvas("ccomp", "", 2100, 1150);
        TPad* top = new TPad("top", "", 0, 0.76, 1, 1);
        top->Draw(); top->cd();
        const int gcol[NG] = {kOrange - 3, kOrange + 7, kRed, kRed + 3, kViolet + 1, kAzure + 2, kCyan + 2, kBlue + 2, kGray + 2, kBlack, kBlack};
        const char* gnice[NG] = {"bump, #pi^{0} Dalitz (cl. 2)", "bump, conversion (cl. 3)", "peak, conversion (cl. 4)", "peak, direct (cl. 5)",
                                 "#eta produced, not signal", "#pi^{0} Dalitz pair", "conversion pair", "#rho^{0}, #omega, #phi pairs", "other T2",
                                 "T3 fake lepton (not in sum)", "T4 two sources (not in sum)"};
        TLegend* lg = new TLegend(0.0, 0.02, 1.0, 0.98);
        lg->SetNColumns(4); lg->SetTextSize(0.125); lg->SetBorderSize(0); lg->SetFillStyle(0);
        TLine* lT3 = new TLine(); lT3->SetLineStyle(7); lT3->SetLineWidth(5); lT3->SetLineColor(kBlack);
        TLine* lT4 = new TLine(); lT4->SetLineStyle(10); lT4->SetLineWidth(5); lT4->SetLineColor(kGray + 2);
        c->cd();
        TPad* bot = new TPad("bot", "", 0, 0, 1, 0.76);
        bot->Draw(); bot->cd();
        bot->Divide(3, 1, 0.002, 0.002);
        const char* sl[5] = {"pre", "best", "frozen", "purity", "slides"};
        for (int t = 0; t < NT; ++t) {
            bot->cd(t + 1);
            gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.10); gPad->SetTopMargin(0.10);
            const Ev* ev[5] = {&evPre[t], &evBest[t], &evFrz[t], &evPur[t], &evPs[t][kSlides]};
            THStack* st = new THStack(Form("st%d", t), "");
            for (int g = 0; g < gFake; ++g) {
                TH1D* h = new TH1D(Form("hc%d_%d", t, g), "", 5, 0, 5);
                h->SetDirectory(nullptr);
                for (int k = 0; k < 5; ++k) h->SetBinContent(k + 1, ev[k]->sG[g] / kLsim);
                h->SetFillColor(gcol[g]); h->SetLineColor(gcol[g]);
                st->Add(h);
                if (t == 0) lg->AddEntry(h, gnice[g], "f");
            }
            if (t == 0) { lg->AddEntry(lT3, gnice[gFake], "l"); lg->AddEntry(lT4, gnice[gDiff], "l"); }
            double mx = 0;
            for (int k = 0; k < 5; ++k) { double s2 = 0; for (int g = 0; g < gFake; ++g) s2 += ev[k]->sG[g] / kLsim; mx = std::max(mx, s2); }
            TH1D* fr = new TH1D(Form("fc%d", t), "", 5, 0, 5);
            fr->SetDirectory(nullptr);
            for (int k = 0; k < 5; ++k) fr->GetXaxis()->SetBinLabel(k + 1, sl[k]);
            fr->SetMaximum(mx * 1.16); fr->SetMinimum(0);
            AxisStyle(fr, 0.06, 0.065, 1.0, 1.35);
            fr->GetXaxis()->SetLabelSize(0.075);
            fr->GetXaxis()->LabelsOption("h");
            fr->GetYaxis()->SetTitle("#sigma in the window [nb]");
            fr->Draw("AXIS");
            st->Draw("HIST SAME");
            for (int k = 0; k < 5; ++k) {
                // white halo under the T3 / T4 lines, so that they stay visible on dark bars
                for (int q = 0; q < 2; ++q) {
                    const double y = ev[k]->sG[q == 0 ? gFake : gDiff] / kLsim;
                    TLine* h = new TLine(k + 0.12, y, k + 0.88, y); h->SetLineWidth(11); h->SetLineColor(kWhite); h->Draw();
                    TLine* a = new TLine(k + 0.12, y, k + 0.88, y);
                    a->SetLineStyle(q == 0 ? 7 : 10); a->SetLineWidth(5); a->SetLineColor(q == 0 ? kBlack : kGray + 2); a->Draw();
                }
                Txt(k + 0.5, mx * 1.03, 0.058, Form("%.2f", ev[k]->pur), 21, 1, false);
            }
            fr->Draw("AXIS SAME");
            Txt(0.20, 0.925, 0.065, Form("%s (%s window)", TName(t), t < 2 ? "peak" : "bump"));
        }
        top->cd();
        lg->Draw();
        Txt(0.755, 0.12, 0.12, "numbers: purity S/(S+B)");
        Save(c, Form("%s_composition", kTag));
    }

    // ---------- 5. 2-D profile maps
    {
        TCanvas* c = new TCanvas("cmap", "", 2100, 1500);
        TPad* top = new TPad("top", "", 0, 0.92, 1, 1);
        top->Draw(); top->cd();
        Txt(0.01, 0.55, 0.24, Form("2-D profiles of F/F_{best} (best over the other seven variables); red contour: 0.95; %s", kLdNote));
        Txt(0.01, 0.12, 0.24, "black star: best cell; magenta triangle: frozen cut; simulation only, physics weights V3");
        c->cd();
        TPad* bot = new TPad("bot", "", 0, 0, 1, 0.92);
        bot->Draw(); bot->cd();
        bot->Divide(3, 2, 0.002, 0.002);
        gStyle->SetPalette(kBird);
        for (int p = 0; p < 2; ++p)
            for (int t = 0; t < NT; ++t) {
                bot->cd(p * 3 + t + 1);
                gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.15); gPad->SetBottomMargin(0.22); gPad->SetTopMargin(0.08);
                const int dx = kPair[p][0], dy = kPair[p][1];
                const int nx = gD[dx].n(), ny = gD[dy].n();
                TH2D* h = new TH2D(Form("hm%d_%d", p, t), "", nx, 0, nx, ny, 0, ny);
                h->SetDirectory(nullptr);
                for (int i = 0; i < nx; ++i) h->GetXaxis()->SetBinLabel(i + 1, ValStr(dx, i).c_str());
                for (int j = 0; j < ny; ++j) h->GetYaxis()->SetBinLabel(j + 1, ValStr(dy, j).c_str());
                for (int i = 0; i < nx; ++i)
                    for (int j = 0; j < ny; ++j) h->SetBinContent(i + 1, j + 1, S0.prof2[t][p][i * ny + j] / S0.best[t][0][0].F);
                h->SetMinimum(0.4); h->SetMaximum(1.0);
                AxisStyle(h, 0.06, 0.062, 1.55, 1.2);
                h->GetZaxis()->SetLabelSize(0.05);
                h->GetXaxis()->SetLabelSize(0.065); h->GetYaxis()->SetLabelSize(0.065);
                h->GetXaxis()->LabelsOption("v");
                h->GetXaxis()->SetTitle(Form("%s [%s]", gD[dx].title.c_str(), gD[dx].unit.c_str()));
                h->GetYaxis()->SetTitle(Form("%s [%s]", gD[dy].title.c_str(), gD[dy].unit.c_str()));
                h->Draw("COLZ");
                TH2D* hc = (TH2D*)h->Clone(Form("hc%d_%d", p, t));
                double lev[1] = {kPlat};
                hc->SetContour(1, lev); hc->SetLineColor(kRed); hc->SetLineWidth(4);
                hc->Draw("CONT3 SAME");
                const Cut kb = gT.Decode(S0.best[t][0][0].cell), kf = gT.Decode(frz[t].cell);
                TMarker* mb = new TMarker(kb.j[dx] + 0.5, kb.j[dy] + 0.5, 29); mb->SetMarkerSize(3.2); mb->SetMarkerColor(kBlack); mb->Draw();
                TMarker* mf = new TMarker(kf.j[dx] + 0.5, kf.j[dy] + 0.5, 23); mf->SetMarkerSize(3.0); mf->SetMarkerColor(cFrz); mf->Draw();
                Txt(0.17, 0.94, 0.06, TName(t));
            }
        Save(c, Form("%s_maps", kTag));
    }
    printf("done after %.0f s\n", sw.RealTime());
}
}  // namespace s5o

void s5_opt_eta(int K = 100) { s5o::Run(K); }
