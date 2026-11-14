// ============================================================================
// s5_spec.C  --  the whole simulated M(pi+pi-e+e-) spectrum and how the eta-selection cuts shape it
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SIMULATION ONLY: no measured data are read and no
// combinatorial background is added (the simulation has none).
//
// Rows: std + lep + eta sets, owned in partition mode A (SimOpt::sets = 7), weight w (PDG) times the physics weights
// of the eta decays, main variant (simphys::Init(); simphys::Install(); SimOpt::phys = simphys::kPhysMain = 108 = V3:
// direct decay QED x form factors x decay angles, bump KLOE-2 Dalitz plot, conversion F_V, background rules).
// sigma = sum w / L_sim, L_sim = 221.05 nb^-1. Counts at data luminosity (S/sqrt(S+B)) use L_data = 6080.5 nb^-1 as a
// scale only: L_data = 6080.5 nb^-1 (UNCONFIRMED).
// Base selection "pre" = SimEngine kPre (isBest >= 0, vz > -500 mm, PT3 bit, vz > -160 mm, corrected momenta, cor_ok).
// Simulated sum = template 1 (signal classes 2-5) + template 2 (physical background). Template 3 (fake lepton) and
// template 4 (leptons from two different sources) are drawn dashed and never summed; X17 (class 6, eta set) carries
// W = 2 BR per row with BR = 1e-5 (hypothesis), dashed, never summed.
//
// Groups (common definition of step s5; exactly these):
//   G_sig2..G_sig5  signal classes 2 (bump, pi0 Dalitz), 3 (bump, pi0 photon converted), 4 (peak, eta -> pi+pi-gamma
//                   converted), 5 (peak, direct eta -> pi+pi-e+e-, eta set)
//   G_eta   template 2 with pion-pair class 2 (both pions from one eta) or lepton-pair class 2 (eta Dalitz pair)
//   G_pi0D  remaining template 2 with a pi0 Dalitz pair (L1)
//   G_conv  remaining template 2 with a conversion pair (L3)
//   G_vec   remaining template 2 with rho0/omega/phi -> e+e-, omega -> pi0 e+e- or a second-decay pair (L4-L6)
//   G_oth   all other template 2 rows (by the class rules these are the eta-set pairs of a forced decay, L10/L11,
//           with pions that are not the pions of the same decay)
//   G_fake  template 3;  G_diff  template 4
//   "multi-pion background" in this macro = G_pi0D + G_conv + G_vec
// M4 windows (GeV/c2, fixed for all agents): below 0.4576 | bump 0.4576-0.5255 | gap 0.5255-0.5333 | peak
//   0.5333-0.5587 | 0.5587-0.70 | 0.70-1.00 | above 1.00 (the gap appears only in the tables, so the windows add up).
// Cuts (corrected momenta): OA > 2 deg; alpha_lab (oa_lab) < 50 deg; M(pipi) < 0.46 GeV/c2 (variant 0.42);
//   alpha_eta (oa_eta) > 140 deg. Sequence pre -> +OA -> +alpha_lab -> +M(pipi) -> +alpha_eta (= "slides"), the 0.42
//   variant (final = "old code" + OA > 2), and the order "most effective first" (cuts ranked by the background
//   B_peak = (T1 + T2) - (classes 4 + 5) each removes alone from "pre" in the peak window).
//
// Statistics: one pass over the rows accumulates, for every selection (the 16 subsets of the four cuts, the 8
// subsets with M(pipi) < 0.42, peak nominal, bump nominal, pre + pions in the vertex fit) x M4 window x group, the
// sums of w and their second moments over event groups (SimEngine grouping: adjacent rows with the same grp in one
// tree; the sets are independent). Every error (cross sections, kept fractions, removals, purities, S/sqrt(S+B)) is
// the event-group error of SimEngine, ratios by the delta method. Spectra: per bin sum w and the sum over groups of
// the squared group sums (as SimEngine WSpec). Cross-checks against SimEngine FillMany / Efficiency /
// SumMany in tables/s5_spec_checks.txt. Double precision everywhere.
//
// Output (prefix s5_spec):
//   figs/s5_spec_full_pre       full spectrum after pre: log 0.28-2.50 stacked, linear zoom 0.40-0.70
//   figs/s5_spec_steps_log      the cut sequence (6 panels: 5 steps + final with 0.42), log 0.28-2.50, same axes
//   figs/s5_spec_steps_zoom     the same, linear zoom 0.40-0.70
//   figs/s5_spec_final          final selection only, large, with the window composition
//   figs/s5_spec_ribbon_nominal fraction of each group kept after each step, per M4 window (nominal order)
//   figs/s5_spec_ribbon_reverse the same in the reverse order (the ranking "most effective first" in the peak window equals the
//                               nominal order; the reverse order = most effective first over all M4); a ribbon_ranked /
//                               ribbon_rankall figure is drawn only when that order differs from these two
//   figs/s5_spec_attribution    removal attributed to each cut (bump window, peak window, all M4): the orders, alone, last
//   figs/s5_spec_purity         purity of classes 2-5 and 4+5, multi-pion share, S/sqrt(S+B) at L_data (UNCONFIRMED)
//   figs/s5_spec_presets        final-selection spectra of slides, old code, peak nominal, bump nominal
//   tables/s5_spec_windows.txt  sigma per group and window per selection (+ N rows, N_eff), composition by set
//   tables/s5_spec_flow.txt     ranking, cut flows (sigma, kept fraction) per group and window, attribution
//   tables/s5_spec_purity.txt   purities, multi-pion share, dominant group, S/sqrt(S+B) per step
//   tables/s5_spec_checks.txt   cross-checks against SimEngine and reference numbers
//   /hdd3/.../scan/s5_spec_hists.root  dsigma/dM4 histograms (TH1D) of the drawn selections
// Run from /home/przygoda/HADES/FAT/FAT (about 5 min, about 9 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_spec", kTRUE)' 'plots/eta_sim/s5_spec.C+'
// ============================================================================
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace s5spec {
using namespace simeng;

// ---------------------------------------------------------------- constants
const double kBuLo = 0.4576, kBuHi = 0.5255, kPkLo = 0.5333, kPkHi = 0.5587;
const double kBRx17 = 1e-5;
const double kM4Lim = 1.5907;   // sqrt(s) - 2 m_p for T = 4.53 GeV: largest M4 of p p -> p p pi+ pi- e+ e- (X)
const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan";
const char* const kLdTex = "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)";

std::string Fmt(const char* fmt, ...) {
    char buf[8192];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    return std::string(buf);
}
FILE* gT = nullptr;
void T(const char* fmt, ...) {
    if (!gT) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(gT, fmt, ap);
    va_end(ap);
}

// ---------------------------------------------------------------- M4 windows
enum { wBelow = 0, wBump, wGap, wPeak, wMid, wHigh, wTop, wAll, kNW };
const char* const kWinTxt[kNW] = {"<0.4576", "bump", "gap", "peak", "0.5587-0.70", "0.70-1.00", ">1.00", "all M4"};
const char* const kWinTex[kNW] = {"M_{4} < 0.4576 GeV/c^{2}", "bump window 0.4576-0.5255", "gap 0.5255-0.5333",
                                  "peak window 0.5333-0.5587", "M_{4} 0.5587-0.70 GeV/c^{2}", "M_{4} 0.70-1.00 GeV/c^{2}",
                                  "M_{4} > 1.00 GeV/c^{2}", "all M_{4}"};
inline int WinOf(float m4) {
    const double m = m4;
    if (!(m == m)) return -1;
    if (m < kBuLo) return wBelow;
    if (m < kBuHi) return wBump;
    if (m < kPkLo) return wGap;
    if (m < kPkHi) return wPeak;
    if (m < 0.70) return wMid;
    if (m < 1.00) return wHigh;
    return wTop;
}
const int kWinFig[6] = {wBelow, wBump, wPeak, wMid, wHigh, wTop};

// ---------------------------------------------------------------- groups
enum { lS2 = 0, lS3, lS4, lS5, lEta, lPi0D, lConv, lVec, lOth, lFake, lDiff, lX17, kNB };
enum { cT12 = kNB, cSig, cPeak, cBump, cBkg, cMulti, cBPk, cBBu, cAll, cBPkV, cBBuV, kNL };
inline unsigned Mem(int L) {
    if (L < kNB) return 1u << L;
    const unsigned t12 = 0x1FFu, pk = (1u << lS4) | (1u << lS5), bu = (1u << lS2) | (1u << lS3), all = 0x7FFu;
    switch (L) {
        case cT12:   return t12;
        case cSig:   return 0xFu;
        case cPeak:  return pk;
        case cBump:  return bu;
        case cBkg:   return 0x1F0u;
        case cMulti: return (1u << lPi0D) | (1u << lConv) | (1u << lVec);
        case cBPk:   return t12 & ~pk;
        case cBBu:   return t12 & ~bu;
        case cAll:   return all;
        case cBPkV:  return all & ~pk;
        case cBBuV:  return all & ~bu;
    }
    return 0u;
}
const char* const kLabTxt[kNL] = {
    "G_sig2 bump pi0 Dalitz", "G_sig3 bump conversion", "G_sig4 peak conversion", "G_sig5 peak direct",
    "G_eta", "G_pi0D", "G_conv", "G_vec", "G_oth", "G_fake (T3)", "G_diff (T4)", "X17 BR=1e-5 (not in sum)",
    "SUM T1+T2", "eta signal 2-5", "peak 4+5", "bump 2+3", "T2 phys. bkg", "multi-pion pi0D+conv+vec",
    "B_peak=T12-(4+5)", "B_bump=T12-(2+3)", "T1+T2+T3+T4", "B'_peak=T1234-(4+5)", "B'_bump=T1234-(2+3)"};
const char* const kLabLeg[kNL] = {
    "cl. 2 bump: #pi^{0} Dalitz", "cl. 3 bump: #pi^{0} #gamma conv.", "cl. 4 peak: #pi^{+}#pi^{-}#gamma conv.",
    "cl. 5 peak: direct #pi^{+}#pi^{-}e^{+}e^{-}", "G_{#eta}: #eta produced, not signal", "G_{#pi^{0}D}: #pi^{0} Dalitz pair",
    "G_{conv}: conversion pair", "G_{vec}: #rho^{0}, #omega, #phi pairs", "G_{oth}: other T2 (#eta-set pairs)",
    "G_{fake} = T3 (not in sum)", "G_{diff} = T4 (not in sum)", "X17, BR = 10^{-5} (not in sum)",
    "sum T1 + T2", "classes 2-5", "peak 4+5", "bump 2+3", "T2", "multi-pion", "B_{peak}", "B_{bump}", "T1-T4", "B'_{peak}", "B'_{bump}"};
const char* const kLabShort[kNL] = {"cl. 2", "cl. 3", "cl. 4", "cl. 5", "G_{#eta}", "G_{#pi^{0}D}", "G_{conv}", "G_{vec}",
                                    "G_{oth}", "G_{fake}", "G_{diff}", "X17 (10^{-5})", "T1+T2", "2-5", "4+5", "2+3", "T2",
                                    "multi-#pi", "B_{peak}", "B_{bump}", "T1-T4", "B'_{peak}", "B'_{bump}"};

inline int LabelOf(const Row& r) {
    if (r.tmpl == 1) {
        switch (r.cls) {
            case 2: return lS2;
            case 3: return lS3;
            case 4: return lS4;
            case 5: return lS5;
            case 6: return lX17;
        }
        return -1;
    }
    if (r.tmpl == 2) {
        if (r.pcls == 2 || r.lcls == 2) return lEta;
        if (r.lcls == 1) return lPi0D;
        if (r.lcls == 3) return lConv;
        if (r.lcls >= 4 && r.lcls <= 6) return lVec;
        return lOth;
    }
    if (r.tmpl == 3) return lFake;
    if (r.tmpl == 4) return lDiff;
    return -1;
}

// ---------------------------------------------------------------- selections
// id 0-15: subsets of the cuts (bit 0 OA > 2, bit 1 alpha_lab < 50, bit 2 M(pipi) < 0.46, bit 3 alpha_eta > 140);
// id 16-23: the subsets containing bit 2 with M(pipi) < 0.42; 24 peak nominal; 25 bump nominal; 26 pions in vertex fit
const int sPeakNom = 24, sBumpNom = 25, sPiVtx = 26, kNS = 27;
inline int Id42(int mask) { return 16 + ((mask & 3) | ((mask & 8) >> 1)); }
inline int Mask42(int id) { const int i = id - 16; return (i & 3) | 4 | ((i & 4) << 1); }
inline int SelOfMask(int mask, bool v42) { return (v42 && (mask & 4)) ? Id42(mask) : mask; }
const char* const kCutTxt[4] = {"OA>2", "alpha_lab<50", "Mpipi<0.46", "alpha_eta>140"};
const char* const kCutTex[4] = {"OA > 2#circ", "#alpha_{lab} < 50#circ", "M_{#pi#pi} < 0.46 GeV/c^{2}", "#alpha_{#eta} > 140#circ"};
const char* const kCutShort[4] = {"+OA", "+#alpha_{lab}", "+M_{#pi#pi}", "+#alpha_{#eta}"};
std::string SelTxt(int s) {
    if (s == sPeakNom) return "peak nominal (pre + alpha_lab<55 + Mpipi>0.30)";
    if (s == sBumpNom) return "bump nominal (pre + alpha_lab<40 + 0.30<Mpipi<0.41)";
    if (s == sPiVtx) return "pre + both pions in vertex fit";
    const bool v42 = s >= 16;
    const int m = v42 ? Mask42(s) : s;
    std::string n = "pre";
    for (int c = 0; c < 4; ++c)
        if (m >> c & 1) n += std::string(" + ") + (c == 2 && v42 ? "Mpipi<0.42" : kCutTxt[c]);
    if (s == 15) n += "  [= slides]";
    if (s == Id42(14)) n += "  [= old code]";
    if (s == 1) n += "  [= pre+OA>2]";
    return n;
}
inline unsigned SelBits(const Row& r) {
    const int c0 = r.oa > 2.f, c1 = r.oa_lab < 50.f, c2 = r.mpipi < 0.46f, c2b = r.mpipi < 0.42f, c3 = r.oa_eta > 140.f;
    const int rb = c0 | (c1 << 1) | (c2 << 2) | (c3 << 3);
    const int rb42 = c0 | (c1 << 1) | (c2b << 2) | (c3 << 3);
    unsigned m = 0;
    for (int k = 0; k < 16; ++k)
        if ((k & ~rb) == 0) m |= 1u << k;
    for (int k = 0; k < 16; ++k)
        if ((k & 4) && (k & ~rb42) == 0) m |= 1u << Id42(k);
    if (r.oa_lab < 55.f && r.mpipi > 0.30f) m |= 1u << sPeakNom;
    if (r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f) m |= 1u << sBumpNom;
    if (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1) m |= 1u << sPiVtx;
    return m;
}

// ---------------------------------------------------------------- compact rows and accumulator
struct CR {
    int g;            // global event-group index
    float m4;
    double w;         // w (PDG) x physics factor (V3); X17: 2 BR
    unsigned sel;     // bit s: passes selection s
    signed char lab, win, set;
};

const int kNH = kNB + 1;   // histogram labels: base groups + T1+T2 sum
const int kHT12 = kNB;
struct HB {
    int n;
    double lo, hi;
    double bw() const { return (hi - lo) / n; }
};
const HB kHB[2] = {{222, 0.28, 2.50}, {120, 0.40, 0.70}};   // 10 MeV full range, 2.5 MeV zoom
const int kMaxBin = 222;
inline int HBin(int b, double m) {
    const HB& h = kHB[b];
    return (m >= h.lo && m < h.hi) ? std::min(h.n - 1, (int)((m - h.lo) / h.bw())) : -1;
}

struct Acc {
    std::vector<double> S, N, Mlab, Msel, hs, hv, hn, hacc;
    std::vector<int> hlast;
    Acc()
        : S(kNS * kNW * kNL, 0.), N(kNS * kNW * kNL, 0.), Mlab((size_t)kNS * kNW * kNB * kNB, 0.),
          Msel((size_t)kNL * kNW * kNS * kNS, 0.), hs((size_t)kNS * 2 * kNH * kMaxBin, 0.), hv(hs.size(), 0.),
          hn(hs.size(), 0.), hacc(hs.size(), 0.), hlast(hs.size(), -1) {}
    static size_t iS(int s, int k, int L) { return ((size_t)s * kNW + k) * kNL + L; }
    static size_t iML(int s, int k, int a, int b) { return (((size_t)s * kNW + k) * kNB + a) * kNB + b; }
    static size_t iMS(int L, int k, int s1, int s2) { return (((size_t)L * kNW + k) * kNS + s1) * kNS + s2; }
    static size_t iH(int s, int b, int h, int bin) { return (((size_t)s * 2 + b) * kNH + h) * kMaxBin + bin; }
    double Sum(int s, int k, int L) const { return S[iS(s, k, L)]; }
    double Rows(int s, int k, int L) const { return N[iS(s, k, L)]; }
    double Var(int s, int k, int L) const { return Msel[iMS(L, k, s, s)]; }
    double Xs(int s, int k, int L) const { return Sum(s, k, L) / kLsim; }
    double XsErr(int s, int k, int L) const { return std::sqrt(std::max(0., Var(s, k, L))) / kLsim; }
    double Neff(int s, int k, int L) const {
        const double v = Var(s, k, L);
        return v > 0 ? Sum(s, k, L) * Sum(s, k, L) / v : 0.;
    }
};

void Accumulate(const std::vector<CR>& rows, Acc& A) {
    std::vector<std::vector<int>> labsOf(kNB), memOf(kNL);
    for (int l = 0; l < kNB; ++l)
        for (int L = 0; L < kNL; ++L)
            if (Mem(L) >> l & 1u) labsOf[l].push_back(L);
    for (int L = 0; L < kNL; ++L)
        for (int l = 0; l < kNB; ++l)
            if (Mem(L) >> l & 1u) memOf[L].push_back(l);
    const unsigned t12 = Mem(cT12);
    std::vector<double> G((size_t)kNS * kNW * kNB, 0.), C((size_t)kNS * kNW * kNL, 0.);
    std::vector<char> touched(kNS * kNW, 0);
    std::vector<int> tl;
    tl.reserve(512);
    int nsk[kNW];
    int sk[kNW][kNS];
    for (int k = 0; k < kNW; ++k) nsk[k] = 0;
    auto flush = [&]() {
        if (tl.empty()) return;
        for (int id : tl) {
            const int s = id / kNW, k = id % kNW;
            const double* v = &G[(size_t)id * kNB];
            for (int a = 0; a < kNB; ++a) {
                if (v[a] == 0) continue;
                for (int b = 0; b < kNB; ++b)
                    if (v[b] != 0) A.Mlab[Acc::iML(s, k, a, b)] += v[a] * v[b];
            }
            double* c = &C[(size_t)id * kNL];
            for (int L = 0; L < kNL; ++L) {
                double x = 0;
                for (int l : memOf[L]) x += v[l];
                c[L] = x;
            }
            sk[k][nsk[k]++] = s;
        }
        for (int k = 0; k < kNW; ++k) {
            if (!nsk[k]) continue;
            for (int L = 0; L < kNL; ++L)
                for (int i = 0; i < nsk[k]; ++i) {
                    const int s1 = sk[k][i];
                    const double x1 = C[(size_t)(s1 * kNW + k) * kNL + L];
                    if (x1 == 0) continue;
                    for (int j = 0; j < nsk[k]; ++j) {
                        const int s2 = sk[k][j];
                        const double x2 = C[(size_t)(s2 * kNW + k) * kNL + L];
                        if (x2 != 0) A.Msel[Acc::iMS(L, k, s1, s2)] += x1 * x2;
                    }
                }
            nsk[k] = 0;
        }
        for (int id : tl) {
            std::fill(G.begin() + (size_t)id * kNB, G.begin() + (size_t)(id + 1) * kNB, 0.);
            touched[id] = 0;
        }
        tl.clear();
    };
    auto hadd = [&](size_t cell, double w, int g) {
        if (A.hlast[cell] != g) {
            A.hv[cell] += A.hacc[cell] * A.hacc[cell];
            A.hacc[cell] = 0;
            A.hlast[cell] = g;
        }
        A.hacc[cell] += w;
        A.hs[cell] += w;
        A.hn[cell] += 1;
    };
    int cur = INT_MIN;
    for (const CR& r : rows) {
        if (r.g != cur) {
            flush();
            cur = r.g;
        }
        const int l = r.lab;
        const int ks[2] = {r.win, wAll};
        int bins[2];
        for (int b = 0; b < 2; ++b) bins[b] = HBin(b, (double)r.m4);
        for (int s = 0; s < kNS; ++s) {
            if (!(r.sel >> s & 1u)) continue;
            for (int q = 0; q < 2; ++q) {
                const int k = ks[q];
                if (k < 0) continue;
                const int id = s * kNW + k;
                if (!touched[id]) {
                    touched[id] = 1;
                    tl.push_back(id);
                }
                G[(size_t)id * kNB + l] += r.w;
                for (int L : labsOf[l]) {
                    A.S[Acc::iS(s, k, L)] += r.w;
                    A.N[Acc::iS(s, k, L)] += 1;
                }
            }
            for (int b = 0; b < 2; ++b) {
                if (bins[b] < 0) continue;
                hadd(Acc::iH(s, b, l, bins[b]), r.w, r.g);
                if (t12 >> l & 1u) hadd(Acc::iH(s, b, kHT12, bins[b]), r.w, r.g);
            }
        }
    }
    flush();
    for (size_t c = 0; c < A.hv.size(); ++c) {
        A.hv[c] += A.hacc[c] * A.hacc[c];
        A.hacc[c] = 0;
    }
}

// ---------------------------------------------------------------- derived quantities with event-group errors
struct VE {
    double v = 0, e = 0;
};
// variance of sum_i c_i S(s_i) for one group label L and window k
double VarSelLin(const Acc& A, int L, int k, std::vector<std::pair<int, double>> c) {
    std::sort(c.begin(), c.end());
    std::vector<std::pair<int, double>> m;
    for (const auto& x : c) {
        if (!m.empty() && m.back().first == x.first) m.back().second += x.second;
        else m.push_back(x);
    }
    double v = 0;
    for (const auto& a : m)
        for (const auto& b : m) v += a.second * b.second * A.Msel[Acc::iMS(L, k, a.first, b.first)];
    return v;
}
// fraction of group L in window k kept at selection s, relative to s0
VE Kept(const Acc& A, int L, int k, int s, int s0 = 0) {
    VE r;
    const double S0 = A.Sum(s0, k, L);
    if (!(S0 > 0)) return r;
    r.v = A.Sum(s, k, L) / S0;
    if (s == s0) return r;
    r.e = std::sqrt(std::max(0., VarSelLin(A, L, k, {{s, 1. / S0}, {s0, -r.v / S0}})));
    return r;
}
// (S(a) - S(b)) / S(s0)
VE Removed(const Acc& A, int L, int k, int a, int b, int s0 = 0) {
    VE r;
    const double S0 = A.Sum(s0, k, L);
    if (!(S0 > 0)) return r;
    r.v = (A.Sum(a, k, L) - A.Sum(b, k, L)) / S0;
    r.e = std::sqrt(std::max(0., VarSelLin(A, L, k, {{a, 1. / S0}, {b, -1. / S0}, {s0, -r.v / S0}})));
    return r;
}
double VarLabLin(const Acc& A, int s, int k, const double* c) {
    double v = 0;
    for (int a = 0; a < kNB; ++a) {
        if (c[a] == 0) continue;
        for (int b = 0; b < kNB; ++b)
            if (c[b] != 0) v += c[a] * c[b] * A.Mlab[Acc::iML(s, k, a, b)];
    }
    return v;
}
double SumMask(const Acc& A, int s, int k, unsigned m) {
    double x = 0;
    for (int l = 0; l < kNB; ++l)
        if (m >> l & 1u) x += A.Sum(s, k, l);
    return x;
}
// ratio sum(num) / sum(den), num within den
VE Ratio(const Acc& A, int s, int k, unsigned num, unsigned den) {
    VE r;
    const double Sn = SumMask(A, s, k, num), Sd = SumMask(A, s, k, den);
    if (!(Sd > 0)) return r;
    r.v = Sn / Sd;
    double c[kNB];
    for (int l = 0; l < kNB; ++l) c[l] = (((num >> l) & 1u) - r.v * ((den >> l) & 1u)) / Sd;
    r.e = std::sqrt(std::max(0., VarLabLin(A, s, k, c)));
    return r;
}
// S / sqrt(S + B) at data luminosity: counts = sigma x L_data (L_data UNCONFIRMED)
struct Sig {
    double S = 0, B = 0, F = 0, eF = 0;   // S, B in expected data rows
};
Sig Signif(const Acc& A, int s, int k, unsigned sig, unsigned bkg) {
    Sig o;
    const double r = kLdata / kLsim;
    const double S = SumMask(A, s, k, sig), Tt = SumMask(A, s, k, sig | bkg);
    o.S = S * r;
    o.B = (Tt - S) * r;
    if (!(Tt > 0)) return o;
    const double F = S / std::sqrt(Tt);
    double c[kNB];
    for (int l = 0; l < kNB; ++l)
        c[l] = ((sig >> l) & 1u) / std::sqrt(Tt) - S * (((sig | bkg) >> l) & 1u) / (2. * Tt * std::sqrt(Tt));
    o.F = F * std::sqrt(r);
    o.eF = std::sqrt(std::max(0., VarLabLin(A, s, k, c))) * std::sqrt(r);
    return o;
}

// ---------------------------------------------------------------- sequences
struct Seq {
    std::string key, title;
    std::vector<int> cuts;
    bool v42 = false;
    std::vector<int> sel;
    std::vector<std::string> stepTex;
};
Seq MakeSeq(const std::string& key, const std::string& title, const std::vector<int>& cuts, bool v42) {
    Seq q;
    q.key = key;
    q.title = title;
    q.cuts = cuts;
    q.v42 = v42;
    q.sel.push_back(0);
    q.stepTex.push_back("pre");
    int m = 0;
    for (int c : cuts) {
        m |= 1 << c;
        q.sel.push_back(SelOfMask(m, v42));
        q.stepTex.push_back(c == 2 && v42 ? "+M_{#pi#pi}^{0.42}" : kCutShort[c]);
    }
    return q;
}

// ---------------------------------------------------------------- plotting helpers
int HexCol(const char* h) { return TColor::GetColor(h); }
int FillCol(int L) {
    static const char* f[kNB] = {"#FD8D3C", "#FECC5C", "#BD0026", "#F03B20", "#9E9AC8", "#C6DBEF",
                                 "#6BAED6", "#2171B5", "#D9D9D9", "#525252", "#238B45", "#E7298A"};
    if (L < kNB) return HexCol(f[L]);
    return kBlack;
}
int LineCol(int L) {
    static const char* f[kNB] = {"#F16913", "#D4A017", "#99000D", "#EF3B2C", "#6A51A3", "#4EB3D3",
                                 "#2171B5", "#08306B", "#969696", "#525252", "#238B45", "#E7298A"};
    if (L < kNB) return HexCol(f[L]);
    if (L == cMulti) return HexCol("#08519C");
    return kBlack;
}
const char* const kWinCol[kNW] = {"#969696", "#F16913", "#BDBDBD", "#99000D", "#4EB3D3", "#2171B5", "#08306B", "#000000"};
const int kStack[9] = {lPi0D, lConv, lVec, lOth, lEta, lS2, lS3, lS4, lS5};
int gId = 0;
std::string Nm(const char* p) { return Fmt("%s_%d", p, gId++); }

void Txt(double x, double y, const std::string& s, double size, int col = kBlack, int align = 11) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size);
    t.SetTextColor(col);
    t.SetTextAlign(align);
    t.DrawLatex(x, y, s.c_str());
}
TPad* MkPad(double x1, double y1, double x2, double y2, bool logy, double l = 0.13, double r = 0.03, double t = 0.10,
            double b = 0.15) {
    TPad* p = new TPad(Nm("pad").c_str(), "", x1, y1, x2, y2);
    p->SetLeftMargin(l);
    p->SetRightMargin(r);
    p->SetTopMargin(t);
    p->SetBottomMargin(b);
    p->SetLogy(logy);
    p->SetTicks(1, 1);
    p->Draw();
    p->cd();
    return p;
}
void Axes(TH1D* h, const char* xt, const char* yt, double ts, double yoff = 1.15) {
    h->GetXaxis()->SetTitle(xt);
    h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(ts);
    h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 0.92);
    h->GetYaxis()->SetLabelSize(ts * 0.92);
    h->GetXaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetTitleOffset(yoff);
    h->GetXaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetMaxDigits(4);
}
const char* const kXM4 = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const char* const kYds = "d#sigma/dM_{4} [nb/(GeV/c^{2})]";
TH1D* HSpec(const Acc& A, int s, int b, int h) {
    const HB& hb = kHB[b];
    TH1D* x = new TH1D(Nm("h").c_str(), "", hb.n, hb.lo, hb.hi);
    x->SetDirectory(nullptr);
    x->Sumw2();
    const double f = 1. / (kLsim * hb.bw());
    for (int i = 0; i < hb.n; ++i) {
        const size_t c = Acc::iH(s, b, h, i);
        x->SetBinContent(i + 1, A.hs[c] * f);
        x->SetBinError(i + 1, std::sqrt(A.hv[c]) * f);
    }
    return x;
}
double HMax(const Acc& A, int s, int b, double xlo = -1e9, double xhi = 1e9) {
    TH1D* h = HSpec(A, s, b, kHT12);
    double m = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        if (h->GetBinCenter(i) > xlo && h->GetBinCenter(i) < xhi) m = std::max(m, h->GetBinContent(i) + h->GetBinError(i));
    delete h;
    return m;
}
double HMinPos(const Acc& A, int s, int b) {
    TH1D* h = HSpec(A, s, b, kHT12);
    double m = 1e300;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        if (h->GetBinContent(i) > 0) m = std::min(m, h->GetBinContent(i));
    delete h;
    return m;
}
void VLine(double x, double y1, double y2, int col, int ls, int lw = 2) {
    TLine* l = new TLine(x, y1, x, y2);
    l->SetLineColor(col);
    l->SetLineStyle(ls);
    l->SetLineWidth(lw);
    l->Draw();
}
// stacked spectrum of selection s in binning b (0 full, 1 zoom)
void DrawSpectrum(const Acc& A, int s, int b, bool logy, double ylo, double yhi, double ts, double yoff = 1.15, int ndivFull = 512) {
    const HB& hb = kHB[b];
    TH1D* fr = new TH1D(Nm("fr").c_str(), "", hb.n, hb.lo, hb.hi);
    fr->SetDirectory(nullptr);
    Axes(fr, kXM4, kYds, ts, yoff);
    fr->GetXaxis()->SetNdivisions(b == 0 ? ndivFull : 506);
    fr->SetMinimum(ylo);
    fr->SetMaximum(yhi);
    fr->Draw("axis");
    TH1D* cum[9];
    for (int k = 0; k < 9; ++k) {
        cum[k] = HSpec(A, s, b, kStack[k]);
        if (k) cum[k]->Add(cum[k - 1]);
    }
    for (int k = 8; k >= 0; --k) {
        cum[k]->SetFillColor(FillCol(kStack[k]));
        cum[k]->SetFillStyle(1001);
        cum[k]->SetLineColor(FillCol(kStack[k]));
        cum[k]->SetLineWidth(1);
        cum[k]->Draw("hist same");
    }
    TH1D* hs = HSpec(A, s, b, kHT12);
    hs->SetLineColor(kBlack);
    hs->SetLineWidth(2);
    hs->SetMarkerSize(0);
    hs->SetFillStyle(0);
    hs->Draw("hist same");
    TH1D* he = (TH1D*)hs->Clone(Nm("he").c_str());
    he->SetLineWidth(1);
    he->Draw("E1 same");
    const int dl[3] = {lFake, lDiff, lX17};
    const int st[3] = {9, 10, 9};
    for (int q = 0; q < 3; ++q) {
        TH1D* h = HSpec(A, s, b, dl[q]);
        h->SetLineColor(LineCol(dl[q]));
        h->SetLineStyle(st[q]);
        h->SetLineWidth(dl[q] == lX17 ? 4 : 3);
        h->SetFillStyle(0);
        h->Draw("hist same");
    }
    const int wc = HexCol("#636363");
    if (b == 1) {
        for (double x : {kBuLo, kBuHi, kPkLo, kPkHi}) VLine(x, ylo, yhi, wc, 3);
    } else {
        for (double x : {kBuLo, kBuHi, kPkLo, kPkHi, 0.70, 1.00}) VLine(x, ylo, logy ? ylo * std::pow(yhi / ylo, 0.80) : 0.8 * yhi, wc, 3, 1);
        VLine(kM4Lim, ylo, logy ? ylo * std::pow(yhi / ylo, 0.62) : 0.6 * yhi, HexCol("#08306B"), 2, 2);
    }
    gPad->RedrawAxis();
}
TH1D* Dummy(int lc, int lw, int ls, int fc) {
    TH1D* h = new TH1D(Nm("dum").c_str(), "", 1, 0, 1);
    h->SetDirectory(nullptr);
    h->SetLineColor(lc);
    h->SetLineWidth(lw);
    h->SetLineStyle(ls);
    if (fc >= 0) {
        h->SetFillColor(fc);
        h->SetFillStyle(1001);
    } else
        h->SetFillStyle(0);
    return h;
}
void LegendStack(double x1, double y1, double x2, double y2, double ts, int ncol) {
    TLegend* lg = new TLegend(x1, y1, x2, y2);
    lg->SetNColumns(ncol);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->SetTextSize(ts);
    lg->SetTextFont(42);
    for (int k = 8; k >= 0; --k) lg->AddEntry(Dummy(FillCol(kStack[k]), 1, 1, FillCol(kStack[k])), kLabLeg[kStack[k]], "f");
    lg->AddEntry(Dummy(kBlack, 3, 1, -1), "sum T1 + T2 (#pm event-group error)", "l");
    lg->AddEntry(Dummy(LineCol(lFake), 3, 9, -1), kLabLeg[lFake], "l");
    lg->AddEntry(Dummy(LineCol(lDiff), 3, 10, -1), kLabLeg[lDiff], "l");
    lg->AddEntry(Dummy(LineCol(lX17), 4, 9, -1), kLabLeg[lX17], "l");
    lg->AddEntry(Dummy(HexCol("#636363"), 2, 3, -1), "window edges (bump, peak, 0.70, 1.00)", "l");
    lg->Draw();
}
// 2026-11-13 (review: header text 3-4 pt): one short note line, drawn larger; the second line is empty
std::string NoteSim(int part) {
    if (part == 0) return Fmt("SMASH std+lep+eta, mode A, w(PDG) #times V3; d#sigma/dM = #Sigmaw/L_{sim}/bin, L_{sim} = %.2f nb^{-1}; no CB", kLsim);
    return "";
}
void SaveFig(TCanvas* c, const std::string& name) {
    const std::string base = std::string(kFigDir) + "/" + name;
    c->SaveAs((base + ".png").c_str());
    c->SaveAs((base + ".eps").c_str());
    gSystem->Exec(Fmt("ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", base.c_str(), base.c_str(), base.c_str()).c_str());
}
TPad* HeaderPad(TCanvas* c, double y0) {
    c->cd();
    TPad* p = new TPad(Nm("hdr").c_str(), "", 0, y0, 1, 1);
    p->SetMargin(0, 0, 0, 0);
    p->Draw();
    p->cd();
    return p;
}
std::string Pm(double v, double e, int prec = 3) {
    if (!(e > 0)) return Fmt("%.*g", prec, v);
    return Fmt("%.*g #pm %.2g", prec, v, e);
}

}  // namespace s5spec

// ============================================================================
void s5_spec() {
    using namespace simeng;
    using namespace s5spec;
    TStopwatch swAll;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");

    // ------------------------------------------------------------ physics weights and rows
    if (!simphys::Init()) { printf("s5_spec: simphys::Init failed\n"); return; }
    simphys::Install();
    simphys::LoadAngles();
    SimData d;
    if (!Load(d, kCompact, kMaskUS)) return;
    SimOpt opt;   // sets 7, mode A, w (PDG), e+e-
    opt.phys = simphys::kPhysMain;
    printf("s5_spec: options %s (%s)\n", OptLabel(opt).c_str(), simphys::PhysLabel(opt.phys).c_str());

    std::vector<CR> rows;
    rows.reserve(11000000);
    long long nBad = 0, nRowsSet[4] = {0, 0, 0, 0};
    double wSet[4][kNB];
    for (int i = 0; i < 4; ++i)
        for (int l = 0; l < kNB; ++l) wSet[i][l] = 0;
    double wSetFin[4][kNB];
    for (int i = 0; i < 4; ++i)
        for (int l = 0; l < kNB; ++l) wSetFin[i][l] = 0;
    // rows above the kinematic limit M4 = sqrt(s) - 2 m_p (pre and final): sum w per pion-pair class and per set,
    // and how many of them have a true M4 below the limit (m4_true < 0: a misidentified track)
    double wLim[2][10], wLimSet[2][4], wLimTrueBelow[2], wLimTrueNone[2], wLimAll[2];
    for (int i = 0; i < 2; ++i) {
        for (int p = 0; p < 10; ++p) wLim[i][p] = 0;
        for (int q = 0; q < 4; ++q) wLimSet[i][q] = 0;
        wLimTrueBelow[i] = wLimTrueNone[i] = wLimAll[i] = 0;
    }
    {
        int g = -1;
        for (int set = 1; set <= 3; ++set) {
            const std::vector<Row>& v = d.tree(set, kUS);
            bool first = true;
            int cur = 0;
            for (const Row& r : v) {
                if (first || r.grp != cur) {
                    ++g;
                    cur = r.grp;
                    first = false;
                }
                if (!Owned(r, opt) || !PassPre(r)) continue;
                const int lab = LabelOf(r);
                if (lab < 0) {
                    ++nBad;
                    continue;
                }
                CR c;
                c.g = g;
                c.m4 = r.m4;
                c.w = lab == lX17 ? 2. * kBRx17 : r.w * simphys::RowFactor(r, opt.phys);
                c.sel = SelBits(r);
                c.lab = (signed char)lab;
                c.win = (signed char)WinOf(r.m4);
                c.set = (signed char)set;
                rows.push_back(c);
                ++nRowsSet[set];
                wSet[set][lab] += c.w;
                if (c.sel >> 15 & 1u) wSetFin[set][lab] += c.w;
                if (r.m4 > kM4Lim && (Mem(cT12) >> lab & 1u))
                    for (int i = 0; i < 2; ++i) {
                        if (i == 1 && !(c.sel >> 15 & 1u)) continue;
                        if (r.pcls >= 0 && r.pcls < 10) wLim[i][r.pcls] += c.w;
                        wLimSet[i][set] += c.w;
                        wLimAll[i] += c.w;
                        if (r.m4_true < 0) wLimTrueNone[i] += c.w;
                        else if (r.m4_true < kM4Lim) wLimTrueBelow[i] += c.w;
                    }
            }
        }
    }
    printf("s5_spec: %zu rows after pre (std %lld, lep %lld, eta %lld), %lld without a group label\n", rows.size(), nRowsSet[1],
           nRowsSet[2], nRowsSet[3], nBad);

    Acc A;
    {
        TStopwatch sw;
        Accumulate(rows, A);
        printf("s5_spec: accumulation %.1f s\n", sw.RealTime());
    }

    // ------------------------------------------------------------ checks against SimEngine
    gT = fopen(Fmt("%s/s5_spec_checks.txt", kTabDir).c_str(), "w");
    T("# s5_spec_checks.txt -- plots/eta_sim/s5_spec.C (simulation only)\n");
    T("# Cross-checks of the one-pass accumulator of s5_spec.C against SimEngine.h (FillMany, Efficiency, EfficiencyVs, SumMany)\n");
    T("# and against reference numbers of the earlier steps. Options: %s; %s.\n", OptLabel(opt).c_str(), simphys::PhysLabel(opt.phys).c_str());
    T("# rows after pre (owned mode A, e+e-): std %lld, lep %lld, eta %lld, total %zu; rows without a group label: %lld (expected 0)\n",
      nRowsSet[1], nRowsSet[2], nRowsSet[3], rows.size(), nBad);
    T("# angle lookups of the V3 weight: hit %lld, miss %lld\n\n", simphys::AngHit(), simphys::AngMiss());
    int nFail = 0;
    auto Check = [&](const std::string& what, double mine, double ref, double tolRel) {
        const double rel = ref != 0 ? std::fabs(mine - ref) / std::fabs(ref) : std::fabs(mine);   // ref 0: mine is a relative difference
        const bool ok = rel <= tolRel;
        if (!ok) ++nFail;
        T("%-4s %-92s mine %-14.8g ref %-14.8g rel %.2e (tol %.0e)\n", ok ? "OK" : "FAIL", what.c_str(), mine, ref, rel, tolRel);
    };
    {
        TStopwatch sw;
        Binning bW;
        bW.nbin = 7;
        bW.lo = 0;
        bW.hi = 7;
        bW.x = [](const Row& r) {
            const int k = WinOf(r.m4);
            return k < 0 ? -1. : k + 0.5;
        };
        std::vector<Comp> comps;
        for (int L = 0; L < kNL; ++L) {
            const unsigned m = Mem(L);
            Comp c(kLabTxt[L], [m](const Row& r) {
                const int l = LabelOf(r);
                return l >= 0 && ((m >> l) & 1u);
            });
            if (L == lX17) c.wf = [](const Row&) { return 2. * kBRx17; };
            comps.push_back(c);
        }
        const int pp[7][2] = {{kPre, 0}, {kPreOA2, 1}, {kPrePiVtx, sPiVtx}, {kSlides, 15}, {kOldCode, Id42(14)}, {kPeakNom, sPeakNom}, {kBumpNom, sBumpNom}};
        T("## 1. window sums and event-group variances per group: SimEngine FillMany (window index as binning) against the accumulator\n");
        for (int q = 0; q < 7; ++q) {
            const std::vector<WSpec> sp = FillMany(d, SelPreset(pp[q][0]), comps, opt, bW);
            double mS = 0, mV = 0;
            for (int L = 0; L < kNL; ++L) {
                for (int k = 0; k < 7; ++k) {
                    const double a = A.Sum(pp[q][1], k, L), b = sp[L].s[k];
                    const double va = A.Var(pp[q][1], k, L), vb = sp[L].v[k];
                    mS = std::max(mS, std::fabs(a - b) / std::max(std::fabs(b), 1e-300) * (b != 0 || a != 0));
                    mV = std::max(mV, std::fabs(va - vb) / std::max(std::fabs(vb), 1e-300) * (vb != 0 || va != 0));
                }
                mS = std::max(mS, std::fabs(A.Sum(pp[q][1], wAll, L) - sp[L].Sall) / std::max(std::fabs(sp[L].Sall), 1e-300) * (sp[L].Sall != 0));
                mV = std::max(mV, std::fabs(A.Var(pp[q][1], wAll, L) - sp[L].Vall) / std::max(std::fabs(sp[L].Vall), 1e-300) * (sp[L].Vall != 0));
            }
            Check(Fmt("%s: max rel. diff of sum w over 23 groups x 8 windows", PresetName(pp[q][0])), mS, 0, 1e-9);
            Check(Fmt("%s: max rel. diff of the group variance", PresetName(pp[q][0])), mV, 0, 1e-9);
            if (q == 0 || q == 3) {
                Check(Fmt("%s: sigma(T1+T2) all M4 [nb]", PresetName(pp[q][0])), A.Xs(pp[q][1], wAll, cT12), sp[cT12].Sall / kLsim, 1e-12);
                Check(Fmt("%s: sigma(4+5) peak window [nb]", PresetName(pp[q][0])), A.Xs(pp[q][1], wPeak, cPeak), sp[cPeak].s[wPeak] / kLsim, 1e-12);
            }
        }
        // selection identity
        long long mis[7] = {0, 0, 0, 0, 0, 0, 0};
        for (int set = 1; set <= 3; ++set)
            for (const Row& r : d.tree(set, kUS)) {
                if (!Owned(r, opt) || !PassPre(r)) continue;
                const unsigned b = SelBits(r);
                for (int q = 0; q < 7; ++q)
                    if (((b >> pp[q][1]) & 1u) != (unsigned)PassPreset(pp[q][0], r)) ++mis[q];
            }
        T("\n## 2. selection bits of the accumulator against SimEngine PassPreset (rows passing pre)\n");
        for (int q = 0; q < 7; ++q) Check(Fmt("%s = selection id %d (%s): mismatching rows", PresetName(pp[q][0]), pp[q][1], SelTxt(pp[q][1]).c_str()), (double)mis[q], 0, 0);
        // efficiency (kept fraction) and purity errors
        T("\n## 3. kept fractions and purities with event-group errors against SimEngine Efficiency / EfficiencyVs\n");
        const std::vector<EffRes> ev = EfficiencyVs(d, SelPreset(kPre), SelPreset(kSlides), comps[cPeak], opt, bW);
        for (int k : {wBump, wPeak, wMid}) {
            const VE m = Kept(A, cPeak, k, 15, 0);
            Check(Fmt("peak 4+5 kept pre -> slides, window %s: value", kWinTxt[k]), m.v, ev[k].eps, 1e-10);
            Check(Fmt("peak 4+5 kept pre -> slides, window %s: error", kWinTxt[k]), m.e, ev[k].err, 1e-8);
        }
        {
            const VE m = Kept(A, cPeak, wAll, 15, 0);
            Check("peak 4+5 kept pre -> slides, all M4: value", m.v, ev.back().eps, 1e-10);
            Check("peak 4+5 kept pre -> slides, all M4: error", m.e, ev.back().err, 1e-8);
        }
        const EffRes et = Efficiency(d, SelPreset(kPre), SelPreset(kSlides), CompSim(), opt);
        const VE mt = Kept(A, cT12, wAll, 15, 0);
        Check("T1+T2 kept pre -> slides, all M4: value", mt.v, et.eps, 1e-10);
        Check("T1+T2 kept pre -> slides, all M4: error", mt.e, et.err, 1e-8);
        const Sel slPk = [](const Row& r) { return PassPreset(kSlides, r) && WinOf(r.m4) == wPeak; };
        const Sel isPk = [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5); };
        const EffRes ep = Efficiency(d, slPk, isPk, CompSim(), opt);
        const VE mp = Ratio(A, 15, wPeak, Mem(cPeak), Mem(cT12));
        Check("purity (4+5)/(T1+T2), slides, peak window: value", mp.v, ep.eps, 1e-10);
        Check("purity (4+5)/(T1+T2), slides, peak window: error", mp.e, ep.err, 1e-8);
        const Sel prBu = [](const Row& r) { return PassPreset(kPre, r) && WinOf(r.m4) == wBump; };
        const Sel isBu = [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); };
        const EffRes eb = Efficiency(d, prBu, isBu, CompSim(), opt);
        const VE mb = Ratio(A, 0, wBump, Mem(cBump), Mem(cT12));
        Check("purity (2+3)/(T1+T2), pre, bump window: value", mb.v, eb.eps, 1e-10);
        Check("purity (2+3)/(T1+T2), pre, bump window: error", mb.e, eb.err, 1e-8);
        // spectra
        T("\n## 4. spectra: SimEngine FillMany (CompSim, CompX17) in the binnings of the figures against the accumulator\n");
        for (int b = 0; b < 2; ++b)
            for (int s : {0, 15}) {
                Binning bb;
                bb.nbin = kHB[b].n;
                bb.lo = kHB[b].lo;
                bb.hi = kHB[b].hi;
                const std::vector<WSpec> sp = FillMany(d, SelPreset(s == 0 ? kPre : kSlides), {CompSim(), CompX17(kBRx17)}, opt, bb);
                double mS = 0, mV = 0, mX = 0;
                for (int i = 0; i < bb.nbin; ++i) {
                    const size_t c = Acc::iH(s, b, kHT12, i), cx = Acc::iH(s, b, lX17, i);
                    if (sp[0].s[i] != 0) mS = std::max(mS, std::fabs(A.hs[c] - sp[0].s[i]) / std::fabs(sp[0].s[i]));
                    if (sp[0].v[i] != 0) mV = std::max(mV, std::fabs(A.hv[c] - sp[0].v[i]) / std::fabs(sp[0].v[i]));
                    if (sp[1].s[i] != 0) mX = std::max(mX, std::fabs(A.hs[cx] - sp[1].s[i]) / std::fabs(sp[1].s[i]));
                }
                Check(Fmt("%s, %d bins %.2f-%.2f: max rel. diff of T1+T2 per bin", s == 0 ? "pre" : "slides", bb.nbin, bb.lo, bb.hi), mS, 0, 1e-9);
                Check(Fmt("%s, %d bins %.2f-%.2f: max rel. diff of the T1+T2 variance per bin", s == 0 ? "pre" : "slides", bb.nbin, bb.lo, bb.hi), mV, 0, 1e-9);
                Check(Fmt("%s, %d bins %.2f-%.2f: max rel. diff of X17 per bin", s == 0 ? "pre" : "slides", bb.nbin, bb.lo, bb.hi), mX, 0, 1e-9);
            }
        // references of the earlier steps
        T("\n## 5. reference numbers of the earlier steps (pre, all M4)\n");
        SimOpt o0 = opt;
        o0.phys = 0;
        const std::vector<WSum> r0 = SumMany(d, SelPreset(kPre), {CompSim(), CompSigCls(2), CompSigCls(3), CompSigCls(4), CompSigCls(5)}, o0);
        Check("no physics weights: sigma(T1+T2) [nb] (README key numbers 701.8 +- 1.5)", r0[0].Xsec(), 701.8, 0.0005);
        // classes 2-4: rules v2 (2026-11-13) values of report.tex tab:cls (rules v1: 8.05, 5.93, 1.18 and 7.16, 5.29, 1.20)
        Check("no physics weights: sigma class 2 [nb] (report.tex 7.92)", r0[1].Xsec(), 7.92, 0.001);
        Check("no physics weights: sigma class 3 [nb] (report.tex 5.82)", r0[2].Xsec(), 5.82, 0.001);
        Check("no physics weights: sigma class 4 [nb] (report.tex 1.23)", r0[3].Xsec(), 1.23, 0.005);
        Check("no physics weights: sigma class 5 [nb] (report.tex 1.385)", r0[4].Xsec(), 1.385, 0.0005);
        Check("V3: sigma class 2 [nb] (report.tex 7.05)", A.Xs(0, wAll, lS2), 7.05, 0.001);
        Check("V3: sigma class 3 [nb] (report.tex 5.19)", A.Xs(0, wAll, lS3), 5.19, 0.001);
        Check("V3: sigma class 4 [nb] (report.tex 1.24)", A.Xs(0, wAll, lS4), 1.24, 0.005);
        Check("V3: sigma class 5 [nb] (report.tex 1.417)", A.Xs(0, wAll, lS5), 1.417, 0.0005);
        Check("V3 / no weights: T1+T2 all M4 (README x 0.999)", A.Xs(0, wAll, cT12) / r0[0].Xsec(), 0.999, 0.0006);
        Check("X17 BR 1e-5, pre, all M4 [nb] (README 0.0713)", A.Xs(0, wAll, lX17), 0.0713, 0.001);
        // group consistency
        T("\n## 6. internal consistency\n");
        double dW = 0;
        for (int s = 0; s < kNS; ++s)
            for (int L = 0; L < kNL; ++L) {
                double sw = 0;
                for (int k = 0; k < 7; ++k) sw += A.Sum(s, k, L);
                dW = std::max(dW, std::fabs(sw - A.Sum(s, wAll, L)) / std::max(1e-300, std::fabs(A.Sum(s, wAll, L))) * (A.Sum(s, wAll, L) != 0));
            }
        Check("windows (incl. gap) add up to all M4, every selection and group (max rel. diff, summation order)", dW, 0, 1e-9);
        double dG = 0;
        for (int s = 0; s < kNS; ++s)
            for (int k = 0; k < kNW; ++k) {
                double sw = 0;
                for (int l = 0; l <= lOth; ++l) sw += A.Sum(s, k, l);
                dG = std::max(dG, std::fabs(sw - A.Sum(s, k, cT12)) / std::max(1e-300, std::fabs(A.Sum(s, k, cT12))) * (A.Sum(s, k, cT12) != 0));
            }
        Check("classes 2-5 + G_eta + G_pi0D + G_conv + G_vec + G_oth = T1+T2 (max rel. diff, summation order)", dG, 0, 1e-9);
        T("\n# %d checks outside the tolerance\n", nFail);
        printf("s5_spec: checks done in %.0f s, %d outside the tolerance\n", sw.RealTime(), nFail);
    }
    if (gT) fclose(gT);
    gT = nullptr;
    // free the trees
    for (int k = 0; k < kNTreeAll; ++k) {
        d.t[k].clear();
        d.t[k].shrink_to_fit();
    }

    // ------------------------------------------------------------ sequences and ranking
    const Seq qNom = MakeSeq("nominal", "nominal order: OA, #alpha_{lab}, M_{#pi#pi} < 0.46, #alpha_{#eta}", {0, 1, 2, 3}, false);
    const Seq q42 = MakeSeq("v042", "variant M_{#pi#pi} < 0.42: OA, #alpha_{lab}, M_{#pi#pi} < 0.42, #alpha_{#eta}", {0, 1, 2, 3}, true);
    VE remAlone[4];
    for (int c = 0; c < 4; ++c) remAlone[c] = Removed(A, cBPk, wPeak, 0, 1 << c, 0);
    std::vector<int> rk = {0, 1, 2, 3};
    std::stable_sort(rk.begin(), rk.end(), [&](int a, int b) { return remAlone[a].v > remAlone[b].v; });
    std::vector<int> greedy;
    {
        int m = 0;
        std::vector<int> left = {0, 1, 2, 3};
        while (!left.empty()) {
            int best = -1;
            double bv = -1;
            for (int c : left) {
                const double r = A.Sum(m, wPeak, cBPk) - A.Sum(m | (1 << c), wPeak, cBPk);
                if (r > bv) {
                    bv = r;
                    best = c;
                }
            }
            greedy.push_back(best);
            m |= 1 << best;
            left.erase(std::find(left.begin(), left.end(), best));
        }
    }
    // the same ranking with the background of the whole spectrum: T2 (= T1+T2 - classes 2-5) over all M4
    VE remAllM4[4];
    for (int c = 0; c < 4; ++c) remAllM4[c] = Removed(A, cBkg, wAll, 0, 1 << c, 0);
    std::vector<int> rkA = {0, 1, 2, 3};
    std::stable_sort(rkA.begin(), rkA.end(), [&](int a, int b) { return remAllM4[a].v > remAllM4[b].v; });
    const std::vector<int> rev(rk.rbegin(), rk.rend());
    const Seq qRk = MakeSeq("ranked", "ranked by B_{peak} removed in the peak window", rk, false);
    const Seq qRev = MakeSeq("reverse", "reverse of the ranked order (least effective first)", rev, false);
    const Seq qRkA = MakeSeq("rankall", "ranked by the T2 background removed over all M_{4}", rkA, false);
    const bool rkIsNom = rk == qNom.cuts;
    const bool rkADistinct = rkA != rk && rkA != rev && rkA != qNom.cuts;
    std::vector<const Seq*> seqs = {&qNom, &q42};
    if (!rkIsNom) seqs.push_back(&qRk);
    seqs.push_back(&qRev);
    if (rkADistinct) seqs.push_back(&qRkA);
    // orders compared in the attribution (first = nominal; the ranked order is shown separately only if it differs)
    std::vector<std::pair<std::string, const Seq*>> ords = {{rkIsNom ? "nominal = ranked" : "nominal", &qNom}};
    if (!rkIsNom) ords.push_back({"ranked (peak B)", &qRk});
    ords.push_back({rkA == rev ? "reverse = all-M4 ranked" : "reverse", &qRev});
    if (rkADistinct) ords.push_back({"ranked (all-M4 B)", &qRkA});
    printf("s5_spec: ranked order (peak-window B_peak):");
    for (int c : rk) printf(" %s (%.3f)", kCutTxt[c], remAlone[c].v);
    printf("%s; greedy:", rkIsNom ? " = nominal" : "");
    for (int c : greedy) printf(" %s", kCutTxt[c]);
    printf("; all-M4 T2 ranking:");
    for (int c : rkA) printf(" %s (%.3f)", kCutTxt[c], remAllM4[c].v);
    printf("\n");

    // ------------------------------------------------------------ table: windows
    std::vector<int> tabSel = {0, 1, 3, 7, 15, Id42(7), Id42(15)};
    for (const Seq* q : seqs)
        for (int s : q->sel)
            if (std::find(tabSel.begin(), tabSel.end(), s) == tabSel.end()) tabSel.push_back(s);
    for (int s : {Id42(14), sPeakNom, sBumpNom, sPiVtx})
        if (std::find(tabSel.begin(), tabSel.end(), s) == tabSel.end()) tabSel.push_back(s);
    gT = fopen(Fmt("%s/s5_spec_windows.txt", kTabDir).c_str(), "w");
    T("# s5_spec_windows.txt -- plots/eta_sim/s5_spec.C (SIMULATION ONLY, no CB)\n");
    T("# sigma [nb] = sum w / L_sim (L_sim = %.2f nb^-1) per group and M4 window; std + lep + eta owned rows, mode A, e+e-,\n", kLsim);
    T("# w (PDG) x physics weights V3 (%s). Errors: event groups. [N_eff] = (sum w)^2 / group variance.\n", simphys::PhysLabel(opt.phys).c_str());
    T("# Windows [GeV/c2]: <0.4576 | bump 0.4576-0.5255 | gap 0.5255-0.5333 | peak 0.5333-0.5587 | 0.5587-0.70 | 0.70-1.00 | >1.00 | all.\n");
    T("# Sum T1+T2 = classes 2-5 + G_eta + G_pi0D + G_conv + G_vec + G_oth; G_fake (T3), G_diff (T4), X17 (BR 1e-5) are not in the sum.\n");
    T("# Composite rows below the groups: classes 2-5, peak 4+5, bump 2+3, T2, multi-pion = G_pi0D+G_conv+G_vec, B_peak = T1+T2-(4+5),\n");
    T("# B_bump = T1+T2-(2+3), T1+T2+T3+T4 and the B' with T3+T4.\n");
    for (int s : tabSel) {
        T("\n## selection id %d: %s\n", s, SelTxt(s).c_str());
        T("%-26s", "group");
        for (int k = 0; k < kNW; ++k) T(" %30s", kWinTxt[k]);
        T(" %10s\n", "rows(all)");
        for (int L = 0; L < kNL; ++L) {
            T("%-26s", kLabTxt[L]);
            for (int k = 0; k < kNW; ++k)
                T(" %10.4g +- %-8.2g [%7.0f]", A.Xs(s, k, L), A.XsErr(s, k, L), A.Neff(s, k, L));
            T(" %10.0f\n", A.Rows(s, wAll, L));
        }
    }
    T("\n## composition by set: sigma [nb] of each group, all M4 (no errors; sums of w / L_sim)\n");
    T("%-26s %12s %12s %12s | %12s %12s %12s\n", "group", "pre: std", "lep", "eta", "slides: std", "lep", "eta");
    for (int l = 0; l < kNB; ++l)
        T("%-26s %12.5g %12.5g %12.5g | %12.5g %12.5g %12.5g\n", kLabTxt[l], wSet[1][l] / kLsim, wSet[2][l] / kLsim, wSet[3][l] / kLsim,
          wSetFin[1][l] / kLsim, wSetFin[2][l] / kLsim, wSetFin[3][l] / kLsim);
    T("# G_oth: by the class rules only eta-set rows (lepton-pair class L10/L11 of a forced decay, pions not of the same decay).\n");
    T("\n## T1+T2 above the kinematic limit M4 > sqrt(s) - 2 m_p = %.4f GeV/c2 (p p -> p p pi+ pi- e+ e-): sigma [nb] (sums of w / L_sim)\n", kM4Lim);
    for (int i = 0; i < 2; ++i) {
        T("%-8s all %.4g nb (of T1+T2 %.4g nb); true M4 below the limit %.4g nb, a misidentified track (no true M4) %.4g nb; by set std %.4g, lep %.4g, eta %.4g\n",
          i == 0 ? "pre" : "slides", wLimAll[i] / kLsim, A.Xs(i == 0 ? 0 : 15, wAll, cT12), wLimTrueBelow[i] / kLsim, wLimTrueNone[i] / kLsim,
          wLimSet[i][1] / kLsim, wLimSet[i][2] / kLsim, wLimSet[i][3] / kLsim);
        T("         by pion-pair class:");
        for (int p = 0; p < 10; ++p) T(" P%d %.4g;", p, wLim[i][p] / kLsim);
        T("\n");
    }
    fclose(gT);
    gT = nullptr;

    // ------------------------------------------------------------ table: flow
    const int flowLab[] = {lS2, lS3, lS4, lS5, lEta, lPi0D, lConv, lVec, lOth, lFake, lDiff, lX17, cT12, cSig, cPeak, cBump, cMulti, cBPk, cBBu};
    gT = fopen(Fmt("%s/s5_spec_flow.txt", kTabDir).c_str(), "w");
    T("# s5_spec_flow.txt -- plots/eta_sim/s5_spec.C (SIMULATION ONLY, no CB)\n");
    T("# Cut flows per group and M4 window; same rows, weights and errors as s5_spec_windows.txt (event groups, delta method).\n");
    T("# Cuts: OA > 2 deg; alpha_lab < 50 deg; M(pipi) < 0.46 GeV/c2 (variant 0.42); alpha_eta > 140 deg. kept = sigma(step) / sigma(pre).\n");
    T("\n## 1. ranking 'most effective first': fraction of B_peak = (T1+T2) - (classes 4+5) in the peak window removed by each cut alone (from pre)\n");
    for (int c = 0; c < 4; ++c) {
        const VE ks = Kept(A, cPeak, wPeak, 1 << c, 0);
        T("  %-14s removes B_peak %.4f +- %.4f   (keeps peak 4+5: %.4f +- %.4f; B_peak N_eff at pre %.0f)\n", kCutTxt[c], remAlone[c].v, remAlone[c].e,
          ks.v, ks.e, A.Neff(0, wPeak, cBPk));
    }
    T("  ranked order:");
    for (int c : rk) T(" %s", kCutTxt[c]);
    T("%s\n  greedy order (largest additional B_peak removal at each step):", rkIsNom ? "  = the nominal order" : "");
    for (int c : greedy) T(" %s", kCutTxt[c]);
    T("\n  => the ranking in the peak window reproduces the nominal order; to show how the attribution depends on the order, the\n");
    T("     reverse order (least effective first) is added%s.\n", rkADistinct ? ", and the order ranked by the T2 background removed over all M4" : "");
    T("\n## 1b. the same, background removed by each cut alone from pre, per window (B = T1+T2 minus the window's signal in the bump and peak\n");
    T("##     windows, T2 elsewhere) and over all M4 (T2)\n");
    T("%-14s", "cut");
    for (int k = 0; k < kNW; ++k) T(" %18s", kWinTxt[k]);
    T("\n");
    for (int c = 0; c < 4; ++c) {
        T("%-14s", kCutTxt[c]);
        for (int k = 0; k < kNW; ++k) {
            const int L = k == wBump ? cBBu : k == wPeak ? cBPk : cBkg;
            const VE r = Removed(A, L, k, 0, 1 << c, 0);
            T(" %8.4f +- %-7.4f", r.v, r.e);
        }
        T("\n");
    }
    T("  all-M4 ranking (T2):");
    for (int c : rkA) T(" %s", kCutTxt[c]);
    T("%s\n", rkADistinct ? "" : "  (equal to the nominal, ranked or reverse order: not shown separately)");
    for (const Seq* q : seqs) {
        T("\n## 2. sequence '%s': steps", q->key.c_str());
        for (size_t j = 0; j < q->sel.size(); ++j) T(" [%zu] %s;", j, SelTxt(q->sel[j]).c_str());
        T("\n");
        for (int k = 0; k < kNW; ++k) {
            T("\n### window %s: sigma [nb] +- err (top) and kept fraction +- err (bottom) per step\n", kWinTxt[k]);
            T("%-26s", "group");
            for (size_t j = 0; j < q->sel.size(); ++j) T(" %22s", Fmt("step %zu", j).c_str());
            T(" %9s\n", "Neff(pre)");
            for (int L : flowLab) {
                if (A.Sum(0, k, L) == 0) continue;
                T("%-26s", kLabTxt[L]);
                for (int s : q->sel) T(" %10.4g +- %-9.2g", A.Xs(s, k, L), A.XsErr(s, k, L));
                T(" %9.0f\n", A.Neff(0, k, L));
                T("%-26s", "   kept");
                for (int s : q->sel) {
                    const VE e = Kept(A, L, k, s, 0);
                    T(" %10.4f +- %-9.4f", e.v, e.e);
                }
                T(" %9s\n", Fmt("Neff(fin) %.0f", A.Neff(q->sel.back(), k, L)).c_str());
            }
        }
    }
    T("\n## 3. attribution: fraction of the group (sigma at pre) removed by each cut, in the window\n");
    T("#  columns: the cut at its step in each order (nominal = ranked by B_peak; reverse%s); 0.42 = nominal order with M(pipi) < 0.42;\n",
      rkADistinct ? "; ranked by the all-M4 T2" : "");
    T("#  alone = the cut alone on pre; last = unique removal (all four cuts against the three others). The removals of one order add up\n");
    T("#  to 1 - kept(final); alone and last bracket the attribution.\n");
    auto stepOf = [](const Seq& q, int c) {
        for (size_t j = 0; j < q.cuts.size(); ++j)
            if (q.cuts[j] == c) return (int)j + 1;
        return -1;
    };
    for (int k : {wBelow, wBump, wPeak, wMid, wHigh, wTop, wAll}) {
        T("\n### window %s\n", kWinTxt[k]);
        T("%-26s %-14s", "group", "cut");
        for (const auto& o : ords) T(" %18s", o.first.c_str());
        T(" %18s %18s %18s\n", "0.42 order", "alone", "last (unique)");
        for (int L : flowLab) {
            if (A.Sum(0, k, L) == 0 || L == lX17) continue;
            for (int c = 0; c < 4; ++c) {
                T("%-26s %-14s", c == 0 ? kLabTxt[L] : "", kCutTxt[c]);
                for (const auto& o : ords) {
                    const int j = stepOf(*o.second, c);
                    const VE r = Removed(A, L, k, o.second->sel[j - 1], o.second->sel[j], 0);
                    T(" %8.4f +- %-7.4f", r.v, r.e);
                }
                const int jn = stepOf(qNom, c);
                const VE r4 = Removed(A, L, k, q42.sel[jn - 1], q42.sel[jn], 0);
                const VE ra = Removed(A, L, k, 0, SelOfMask(1 << c, false), 0);
                const VE rl = Removed(A, L, k, SelOfMask(15 & ~(1 << c), false), 15, 0);
                T(" %8.4f +- %-7.4f %8.4f +- %-7.4f %8.4f +- %-7.4f\n", r4.v, r4.e, ra.v, ra.e, rl.v, rl.e);
            }
        }
    }
    fclose(gT);
    gT = nullptr;

    // ------------------------------------------------------------ table: purity and S/sqrt(S+B)
    gT = fopen(Fmt("%s/s5_spec_purity.txt", kTabDir).c_str(), "w");
    T("# s5_spec_purity.txt -- plots/eta_sim/s5_spec.C (SIMULATION ONLY, no CB)\n");
    T("# Purity = sigma(group) / sigma(T1+T2) in the window (errors: event groups, delta method). eta = classes 2-5, peak = 4+5,\n");
    T("# bump = 2+3, multi-pion = G_pi0D + G_conv + G_vec, G_eta share. 'dominant' = largest single group of T1+T2 in the window.\n");
    T("# S/sqrt(S+B): peak = classes 4+5 in the peak window, bump = classes 2+3 in the bump window; B = everything else of T1+T2 in that\n");
    T("# window (B' adds T3 + T4); counts = sigma x L_data, L_data = 6080.5 nb^-1 (UNCONFIRMED, a scale only; F scales with sqrt(L_data)).\n");
    T("# Errors of F: simulation statistics only (event groups, delta method).\n");
    for (const Seq* q : seqs) {
        T("\n## sequence '%s'\n", q->key.c_str());
        for (size_t j = 0; j < q->sel.size(); ++j) {
            const int s = q->sel[j];
            T("\n### step %zu: %s\n", j, SelTxt(s).c_str());
            T("%-12s %18s %18s %18s %18s %18s  %-24s %10s\n", "window", "P(2-5)", "P(4+5)", "P(2+3)", "multi-pion", "G_eta", "dominant (share)", "Neff T12");
            for (int k = 0; k < kNW; ++k) {
                if (A.Sum(s, k, cT12) <= 0) continue;
                const VE pe = Ratio(A, s, k, Mem(cSig), Mem(cT12)), pp = Ratio(A, s, k, Mem(cPeak), Mem(cT12)),
                         pb = Ratio(A, s, k, Mem(cBump), Mem(cT12)), pm = Ratio(A, s, k, Mem(cMulti), Mem(cT12)),
                         pg = Ratio(A, s, k, 1u << lEta, Mem(cT12));
                int dom = 0;
                for (int l = 0; l <= lOth; ++l)
                    if (A.Sum(s, k, l) > A.Sum(s, k, dom)) dom = l;
                T("%-12s %8.4f +- %-7.4f %8.4f +- %-7.4f %8.4f +- %-7.4f %8.4f +- %-7.4f %8.4f +- %-7.4f  %-14s (%.3f) %10.0f\n", kWinTxt[k], pe.v, pe.e,
                  pp.v, pp.e, pb.v, pb.e, pm.v, pm.e, pg.v, pg.e, kLabTxt[dom], A.Sum(s, k, dom) / A.Sum(s, k, cT12), A.Neff(s, k, cT12));
            }
            const Sig fp = Signif(A, s, wPeak, Mem(cPeak), Mem(cBPk)), fpv = Signif(A, s, wPeak, Mem(cPeak), Mem(cBPkV));
            const Sig fb = Signif(A, s, wBump, Mem(cBump), Mem(cBBu)), fbv = Signif(A, s, wBump, Mem(cBump), Mem(cBBuV));
            const Sig fp0 = Signif(A, 0, wPeak, Mem(cPeak), Mem(cBPk)), fb0 = Signif(A, 0, wBump, Mem(cBump), Mem(cBBu));
            T("  peak: S %.1f  B %.1f  S/sqrt(S+B) %.3f +- %.3f  (F/F_pre %.3f) | with T3+T4: B' %.1f  F' %.3f +- %.3f   [N_eff S %.0f, B %.0f]\n", fp.S,
              fp.B, fp.F, fp.eF, fp.F / fp0.F, fpv.B, fpv.F, fpv.eF, A.Neff(s, wPeak, cPeak), A.Neff(s, wPeak, cBPk));
            T("  bump: S %.1f  B %.1f  S/sqrt(S+B) %.3f +- %.3f  (F/F_pre %.3f) | with T3+T4: B' %.1f  F' %.3f +- %.3f   [N_eff S %.0f, B %.0f]\n", fb.S,
              fb.B, fb.F, fb.eF, fb.F / fb0.F, fbv.B, fbv.F, fbv.eF, A.Neff(s, wBump, cBump), A.Neff(s, wBump, cBBu));
        }
    }
    T("\n## named presets (final-selection comparison)\n");
    for (int s : {15, Id42(14), sPeakNom, sBumpNom, sPiVtx}) {
        const Sig fp = Signif(A, s, wPeak, Mem(cPeak), Mem(cBPk)), fb = Signif(A, s, wBump, Mem(cBump), Mem(cBBu));
        const Sig fpv = Signif(A, s, wPeak, Mem(cPeak), Mem(cBPkV)), fbv = Signif(A, s, wBump, Mem(cBump), Mem(cBBuV));
        const VE pp = Ratio(A, s, wPeak, Mem(cPeak), Mem(cT12)), pb = Ratio(A, s, wBump, Mem(cBump), Mem(cT12));
        T("%-58s peak: P %.4f +- %.4f, F %.2f +- %.2f (F' %.2f) | bump: P %.4f +- %.4f, F %.2f +- %.2f (F' %.2f)\n", SelTxt(s).c_str(), pp.v, pp.e, fp.F,
          fp.eF, fpv.F, pb.v, pb.e, fb.F, fb.eF, fbv.F);
    }
    fclose(gT);
    gT = nullptr;

    // ------------------------------------------------------------ histograms to a ROOT file
    {
        TFile f(Fmt("%s/s5_spec_hists.root", kOutDir).c_str(), "RECREATE");
        std::vector<int> hs = tabSel;
        for (int s : hs)
            for (int b = 0; b < 2; ++b)
                for (int h = 0; h < kNH; ++h) {
                    TH1D* x = HSpec(A, s, b, h);
                    static const char* hn[kNH] = {"sig2", "sig3", "sig4", "sig5", "Geta", "Gpi0D", "Gconv", "Gvec", "Goth", "Gfake", "Gdiff", "X17_br1e-5", "sumT12"};
                    x->SetName(Fmt("h_sel%02d_%s_%s", s, b == 0 ? "full" : "zoom", hn[h]).c_str());
                    x->SetTitle(Fmt("%s; M4 [GeV/c2]; dsigma/dM4 [nb/(GeV/c2)]", SelTxt(s).c_str()).c_str());
                    x->Write();
                    delete x;
                }
        f.Close();
    }

    // ============================================================ figures
    const double kYlo = 1e-3;   // lower edge of the log panels [nb/(GeV/c2)] (2026-11-13: 1e-3 and top x 1e3, review: 4 empty decades)
    // ---- 1. full spectrum after pre
    {
        TCanvas* c = new TCanvas("cFull", "", 2600, 2500);
        HeaderPad(c, 0.80);
        Txt(0.012, 0.87, "Simulated d#sigma/dM_{4} after \"pre\" (isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta), stacked by group", 0.095);
        Txt(0.012, 0.765, NoteSim(0), 0.088);
        Txt(0.012, 0.675, NoteSim(1), 0.060);
        Txt(0.012, 0.615, "Stack: T1 + T2 (signals warm, backgrounds cold/grey); dashed, not stacked: T3, T4, X17", 0.078);
        LegendStack(0.005, 0.0, 0.995, 0.585, 0.074, 4);
        c->cd();
        const double ymx = HMax(A, 0, 0);
        MkPad(0.0, 0.40, 1.0, 0.80, true, 0.09, 0.02, 0.09, 0.14);
        DrawSpectrum(A, 0, 0, true, kYlo, ymx * 1e3, 0.060, 0.72);
        Txt(0.10, 0.93, "(a) pre, log scale, 10 MeV/c^{2} bins, M_{4} 0.28-2.50 GeV/c^{2}", 0.060);
        Txt(0.36, 0.82, Fmt("#sigma_{T1+T2} = %.1f #pm %.1f nb  (T3: %.1f nb, T4: %.1f nb, dashed)", A.Xs(0, wAll, cT12), A.XsErr(0, wAll, cT12),
                            A.Xs(0, wAll, lFake), A.Xs(0, wAll, lDiff)), 0.046);
        Txt(0.36, 0.755, Fmt("classes 2-5: %.2f nb; multi-pion (G_{#pi^{0}D}+G_{conv}+G_{vec}): %.1f nb", A.Xs(0, wAll, cSig), A.Xs(0, wAll, cMulti)), 0.046);
        Txt(0.36, 0.69, Fmt("dashed blue: M_{4} = #sqrt{s} - 2m_{p} = %.3f GeV/c^{2}; above it %.2f nb (mostly misidentified tracks)", kM4Lim, wLimAll[0] / kLsim), 0.046,
            HexCol("#08306B"));
        c->cd();
        MkPad(0.0, 0.0, 1.0, 0.40, false, 0.09, 0.02, 0.09, 0.14);
        const double zm = HMax(A, 0, 1);
        DrawSpectrum(A, 0, 1, false, 0, zm * 1.32, 0.060, 0.72);
        Txt(0.10, 0.93, "(b) pre, linear zoom 0.40-0.70 GeV/c^{2}, 2.5 MeV/c^{2} bins; dotted: bump and peak windows", 0.060);
        Txt(0.12, 0.80, Fmt("bump window: #sigma_{2+3} = %s nb, P_{2+3} = %.3f", Pm(A.Xs(0, wBump, cBump), A.XsErr(0, wBump, cBump)).c_str(),
                            Ratio(A, 0, wBump, Mem(cBump), Mem(cT12)).v), 0.050);
        Txt(0.12, 0.73, Fmt("peak window: #sigma_{4+5} = %s nb, P_{4+5} = %.3f", Pm(A.Xs(0, wPeak, cPeak), A.XsErr(0, wPeak, cPeak)).c_str(),
                            Ratio(A, 0, wPeak, Mem(cPeak), Mem(cT12)).v), 0.050);
        Txt(0.12, 0.66, Fmt("X17 (BR 10^{-5}, dashed, at the bottom): %.3f nb in the peak window", A.Xs(0, wPeak, lX17)), 0.050, LineCol(lX17));
        SaveFig(c, "s5_spec_full_pre");
    }

    // ---- 2./3. cut sequence: log and zoom, six panels
    const int stepSel[6] = {qNom.sel[0], qNom.sel[1], qNom.sel[2], qNom.sel[3], qNom.sel[4], q42.sel[4]};
    const char* stepTitle[6] = {"(a) pre", "(b) pre + OA > 2#circ", "(c) + #alpha_{lab} < 50#circ", "(d) + M_{#pi#pi} < 0.46 GeV/c^{2}",
                                "(e) + #alpha_{#eta} > 140#circ = final (= slides)", "(f) final with M_{#pi#pi} < 0.42 (= old code + OA > 2#circ)"};
    for (int fig = 0; fig < 2; ++fig) {
        const int b = fig == 0 ? 0 : 1;
        TCanvas* c = new TCanvas(Fmt("cSteps%d", fig).c_str(), "", 2800, 2900);
        HeaderPad(c, 0.865);
        Txt(0.010, 0.86, fig == 0 ? "Cut sequence: simulated d#sigma/dM_{4} after each step, log scale, M_{4} 0.28-2.50 GeV/c^{2} (same axes in all panels)"
                                  : "Cut sequence: simulated d#sigma/dM_{4} after each step, linear zoom 0.40-0.70 GeV/c^{2} (y range of each panel differs)",
            0.105);
        Txt(0.010, 0.73, NoteSim(0), 0.088);
        Txt(0.010, 0.63, NoteSim(1), 0.064);
        LegendStack(0.005, 0.0, 0.995, 0.66, 0.082, 4);
        double yhi = 0;
        for (int i = 0; i < 6; ++i) yhi = std::max(yhi, HMax(A, stepSel[i], b));
        for (int i = 0; i < 6; ++i) {
            const int s = stepSel[i];
            c->cd();
            const int row = i / 2, col = i % 2;
            const double y2 = 0.865 - row * (0.865 / 3.), y1 = y2 - 0.865 / 3.;
            MkPad(col * 0.5, y1, col * 0.5 + 0.5, y2, fig == 0, 0.14, 0.025, 0.10, 0.15);
            if (fig == 0) DrawSpectrum(A, s, 0, true, kYlo, yhi * 1e3, 0.058, 1.12);
            else DrawSpectrum(A, s, 1, false, 0, HMax(A, s, 1) * 1.55, 0.058, 1.12);
            Txt(0.15, 0.925, stepTitle[i], 0.058);
            const VE pk = Ratio(A, s, wPeak, Mem(cPeak), Mem(cT12)), pb = Ratio(A, s, wBump, Mem(cBump), Mem(cT12));
            const double tx = fig == 0 ? 0.36 : 0.17;
            Txt(tx, 0.82, Fmt("#sigma_{T1+T2} = %s nb", Pm(A.Xs(s, wAll, cT12), A.XsErr(s, wAll, cT12), 4).c_str()), 0.048);
            Txt(tx, 0.755, Fmt("peak: #sigma_{4+5} = %s nb, P = %.3f", Pm(A.Xs(s, wPeak, cPeak), A.XsErr(s, wPeak, cPeak)).c_str(), pk.v), 0.048);
            Txt(tx, 0.69, Fmt("bump: #sigma_{2+3} = %s nb, P = %.3f", Pm(A.Xs(s, wBump, cBump), A.XsErr(s, wBump, cBump)).c_str(), pb.v), 0.048);
        }
        SaveFig(c, fig == 0 ? "s5_spec_steps_log" : "s5_spec_steps_zoom");
    }

    // ---- 4. final selection, large
    {
        const int s = 15;
        TCanvas* c = new TCanvas("cFinal", "", 3000, 2500);
        HeaderPad(c, 0.82);
        Txt(0.010, 0.86, "Final selection: pre + OA > 2#circ + #alpha_{lab} < 50#circ + M_{#pi#pi} < 0.46 GeV/c^{2} + #alpha_{#eta} > 140#circ (= \"slides\")", 0.100);
        Txt(0.010, 0.74, NoteSim(0), 0.088);
        Txt(0.010, 0.645, NoteSim(1), 0.060);
        LegendStack(0.005, 0.0, 0.995, 0.66, 0.076, 4);
        c->cd();
        const double ymx = HMax(A, s, 0);
        MkPad(0.0, 0.41, 0.66, 0.82, true, 0.12, 0.02, 0.09, 0.14);
        DrawSpectrum(A, s, 0, true, kYlo, ymx * 3e3, 0.058, 0.95);
        Txt(0.13, 0.93, "(a) final, log scale, 10 MeV/c^{2} bins", 0.058);
        Txt(0.42, 0.82, Fmt("#sigma_{T1+T2} = %s nb", Pm(A.Xs(s, wAll, cT12), A.XsErr(s, wAll, cT12), 4).c_str()), 0.052);
        Txt(0.42, 0.75, Fmt("M_{4} < 0.70: %.2f nb, > 0.70: %.2f nb", A.Xs(s, wBelow, cT12) + A.Xs(s, wBump, cT12) + A.Xs(s, wGap, cT12) + A.Xs(s, wPeak, cT12) + A.Xs(s, wMid, cT12),
                            A.Xs(s, wHigh, cT12) + A.Xs(s, wTop, cT12)), 0.052);
        c->cd();
        const double zm = HMax(A, s, 1);
        MkPad(0.0, 0.0, 0.66, 0.41, false, 0.12, 0.02, 0.09, 0.14);
        DrawSpectrum(A, s, 1, false, 0, zm * 1.30, 0.058, 0.95);
        Txt(0.13, 0.93, "(b) final, linear zoom, 2.5 MeV/c^{2} bins; dotted: bump and peak windows", 0.058);
        for (int t = 0; t < 2; ++t) {
            const int k = t == 0 ? wBump : wPeak;
            c->cd();
            MkPad(0.66, t == 0 ? 0.41 : 0.0, 1.0, t == 0 ? 0.82 : 0.41, false, 0.02, 0.02, 0.02, 0.02);
            Txt(0.04, 0.92, Fmt("(%c) %s", t == 0 ? 'c' : 'd', t == 0 ? "bump window 0.4576-0.5255" : "peak window 0.5333-0.5587"), 0.058);
            Txt(0.04, 0.845, "group", 0.048);
            Txt(0.66, 0.845, "#sigma [nb]", 0.048, kBlack, 31);
            Txt(0.96, 0.845, "share", 0.048, kBlack, 31);
            double y = 0.78;
            const int ord[13] = {lS2, lS3, lS4, lS5, lEta, lPi0D, lConv, lVec, lOth, cT12, lFake, lDiff, lX17};
            for (int L : ord) {
                const int colr = L < kNB ? LineCol(L) : kBlack;
                Txt(0.04, y, kLabShort[L], 0.046, colr);
                Txt(0.66, y, Pm(A.Xs(s, k, L), A.XsErr(s, k, L), 3), 0.046, colr, 31);
                if (L <= lOth || L == cT12) Txt(0.96, y, Fmt("%.3f", A.Sum(s, k, L) / A.Sum(s, k, cT12)), 0.046, colr, 31);
                y -= L == cT12 ? 0.062 : 0.056;
            }
            const Sig f = t == 0 ? Signif(A, s, wBump, Mem(cBump), Mem(cBBu)) : Signif(A, s, wPeak, Mem(cPeak), Mem(cBPk));
            Txt(0.04, 0.035, Fmt("S/#sqrt{S+B} = %.1f #pm %.1f at L_{data} (UNCONF.)", f.F, f.eF), 0.044);
        }
        SaveFig(c, "s5_spec_final");
    }

    // ---- 5./6. ribbons
    const int ribLab[13] = {lS2, lS3, lS4, lS5, lEta, lPi0D, lConv, lVec, lOth, lFake, lDiff, cT12, cMulti};
    std::vector<const Seq*> ribSeq = {&qNom};
    if (!rkIsNom) ribSeq.push_back(&qRk);
    ribSeq.push_back(&qRev);
    if (rkADistinct) ribSeq.push_back(&qRkA);
    for (const Seq* q : ribSeq) {
        TCanvas* c = new TCanvas(Nm("cRib").c_str(), "", 3000, 2300);
        HeaderPad(c, 0.83);
        std::string what = q->key == "nominal" ? (rkIsNom ? "nominal order (= ranked by the B_{peak} each cut removes alone in the peak window)" : "nominal order")
                         : q->key == "reverse" ? (rkA == rev ? "reverse order: least effective first in the peak window = most effective first over all M_{4} (T2 removed alone)"
                                                             : "reverse order (least effective first in the peak window)")
                         : q->key == "ranked" ? "order ranked by B_{peak} removed in the peak window" : "order ranked by the T2 background removed over all M_{4}";
        Txt(0.010, 0.88, "Ribbon cut flow: fraction of each group kept after each step (relative to pre), per M_{4} window", 0.100);
        Txt(0.010, 0.765, what, 0.080);
        std::string steps = "steps:";
        for (size_t j = 1; j < q->sel.size(); ++j) steps += std::string(" ") + q->stepTex[j] + " = " + kCutTex[q->cuts[j - 1]] + (j + 1 < q->sel.size() ? ";" : "");
        Txt(0.010, 0.655, steps, 0.068);
        Txt(0.010, 0.56, "Errors: event groups (cumulative steps, correlated); groups with N_{eff}(pre) < 3 or zero kept not drawn (log y). Simulation only.", 0.064);
        TLegend* lg = new TLegend(0.005, 0.0, 0.995, 0.50);
        lg->SetNColumns(5);
        lg->SetBorderSize(0);
        lg->SetFillStyle(0);
        lg->SetTextSize(0.072);
        lg->SetTextFont(42);
        for (int L : ribLab) {
            TH1D* dm = Dummy(LineCol(L), L == cT12 ? 5 : 3, (L == lFake ? 9 : L == lDiff ? 10 : L == cMulti ? 7 : 1), -1);
            dm->SetMarkerColor(LineCol(L));
            dm->SetMarkerStyle(L == cT12 ? 29 : L == cMulti ? 33 : L <= lS5 ? 20 : (L == lFake || L == lDiff) ? 24 : 21);
            dm->SetMarkerSize(2.4);
            lg->AddEntry(dm, L == cT12 ? "sum T1 + T2" : L == cMulti ? "multi-pion (#pi^{0}D + conv + vec)" : kLabLeg[L], "lp");
        }
        lg->Draw();
        const int n = (int)q->sel.size();
        for (int i = 0; i < 6; ++i) {
            const int k = kWinFig[i];
            c->cd();
            const int row = i / 3, col = i % 3;
            const double y2 = 0.83 - row * 0.415, y1 = y2 - 0.415;
            MkPad(col / 3., y1, (col + 1) / 3., y2, true, 0.19, 0.03, 0.155, 0.12);
            double ymin = 1;
            for (int L : ribLab)
                for (int s : q->sel) {
                    const VE e = Kept(A, L, k, s, 0);
                    if (e.v > 0 && A.Neff(0, k, L) >= 3) ymin = std::min(ymin, e.v - e.e > 0 ? e.v - e.e : e.v * 0.5);
                }
            TH1D* fr = new TH1D(Nm("frR").c_str(), "", n, 0, n);
            fr->SetDirectory(nullptr);
            Axes(fr, "", "kept fraction (rel. to pre)", 0.054, 1.55);
            for (int j = 0; j < n; ++j) fr->GetXaxis()->SetBinLabel(j + 1, q->stepTex[j].c_str());
            fr->GetXaxis()->SetLabelSize(0.075);
            fr->GetXaxis()->SetLabelOffset(0.012);
            fr->GetXaxis()->LabelsOption("h");
            double lo = std::max(1e-5, ymin * 0.4);
            {
                const double nice[] = {0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01};
                double l2 = std::pow(10., std::floor(std::log10(lo)));
                for (double x : nice)
                    if (x <= lo) { l2 = std::max(l2, x); break; }
                lo = l2;
            }
            fr->SetMinimum(lo);
            fr->SetMaximum(lo >= 0.1 ? 1.6 : 4.0);
            if (lo >= 0.1) fr->GetYaxis()->SetMoreLogLabels();
            fr->GetYaxis()->SetNoExponent(lo > 1e-3);
            fr->Draw("axis");
            TLine* l1 = new TLine(0, 1, n, 1);
            l1->SetLineStyle(3);
            l1->Draw();
            int idx = 0;
            for (int L : ribLab) {
                ++idx;
                if (!(A.Sum(0, k, L) > 0) || A.Neff(0, k, L) < 3) continue;
                TGraphErrors* g = new TGraphErrors();
                for (int j = 0; j < n; ++j) {
                    const VE e = Kept(A, L, k, q->sel[j], 0);
                    if (!(e.v > 0)) continue;
                    const int p = g->GetN();
                    g->SetPoint(p, j + 0.5 + (idx - 7) * 0.022, e.v);
                    g->SetPointError(p, 0, e.e);
                }
                g->SetLineColor(LineCol(L));
                g->SetMarkerColor(LineCol(L));
                g->SetLineWidth(L == cT12 ? 5 : 3);
                g->SetLineStyle(L == lFake ? 9 : L == lDiff ? 10 : L == cMulti ? 7 : 1);
                g->SetMarkerStyle(L == cT12 ? 29 : L == cMulti ? 33 : L <= lS5 ? 20 : (L == lFake || L == lDiff) ? 24 : 21);
                g->SetMarkerSize(L == cT12 ? 2.8 : 2.0);
                g->Draw("LP same");
            }
            Txt(0.19, 0.93, Fmt("(%c) %s", 'a' + i, kWinTex[k]), 0.058);
            Txt(0.19, 0.865, Fmt("#sigma_{T1+T2}: pre %.4g nb, final %.4g nb", A.Xs(0, k, cT12), A.Xs(q->sel.back(), k, cT12)), 0.046);
            gPad->RedrawAxis();
        }
        SaveFig(c, Fmt("s5_spec_ribbon_%s", q->key.c_str()));
    }

    // ---- 7. attribution
    {
        TCanvas* c = new TCanvas("cAttr", "", 3000, 3100);
        HeaderPad(c, 0.865);
        Txt(0.010, 0.88, "Which cut removes what: fraction of the group's #sigma at pre removed by each cut, in the bump window, the peak window and over all M_{4}", 0.090);
        Txt(0.010, 0.77, "Cuts: OA > 2#circ, #alpha_{lab} < 50#circ, M_{#pi#pi} < 0.46 GeV/c^{2}, #alpha_{#eta} > 140#circ. Bars: the cut at its step in each order; the cut alone on pre;", 0.070);
        std::string ordS = "last = unique removal (all four against the other three). Orders:";
        for (const auto& o : ords) {
            ordS += std::string(" ") + o.first + " (";
            for (size_t j = 0; j < o.second->cuts.size(); ++j) ordS += std::string(j ? ", " : "") + kCutShort[o.second->cuts[j]];
            ordS += ")";
        }
        Txt(0.010, 0.68, ordS, 0.070);
        Txt(0.010, 0.59, "B_{bump} = (T1+T2) - (2+3), B_{peak} = (T1+T2) - (4+5), T2 = physical background, multi-pion = G_{#pi^{0}D}+G_{conv}+G_{vec}. The bars of one order add up to 1 - kept(final). Errors: event groups (orange).", 0.066);
        const char* vcol[5] = {"#252525", "#737373", "#BDBDBD", "#3182BD", "#9ECAE1"};
        std::vector<std::string> vname;
        for (const auto& o : ords) vname.push_back(o.first + " order");
        vname.push_back("alone (on pre)");
        vname.push_back("last (unique)");
        const int nv = (int)vname.size();
        std::vector<int> vc;
        for (int v = 0; v < nv; ++v) vc.push_back(HexCol(v < (int)ords.size() ? vcol[v] : vcol[3 + (v - (int)ords.size())]));
        TLegend* lg = new TLegend(0.01, 0.02, 0.99, 0.50);
        lg->SetNColumns(nv);
        lg->SetBorderSize(0);
        lg->SetFillStyle(0);
        lg->SetTextSize(0.080);
        lg->SetTextFont(42);
        for (int v = 0; v < nv; ++v) lg->AddEntry(Dummy(vc[v], 1, 1, vc[v]), vname[v].c_str(), "f");
        lg->Draw();
        const int pl[3][3] = {{cBump, cBBu, cMulti}, {cPeak, cBPk, cMulti}, {cSig, cBkg, cMulti}};
        const char* pn[3][3] = {{"bump 2+3", "B_{bump}", "multi-pion"}, {"peak 4+5", "B_{peak}", "multi-pion"}, {"classes 2-5", "T2", "multi-pion"}};
        const char* wn[3] = {"bump window", "peak window", "all M_{4}"};
        const char* cl[4] = {"OA", "#alpha_{lab}", "M_{#pi#pi}", "#alpha_{#eta}"};
        for (int r = 0; r < 3; ++r)
            for (int q = 0; q < 3; ++q) {
                const int k = r == 0 ? wBump : r == 1 ? wPeak : wAll, L = pl[r][q];
                c->cd();
                const double y2 = 0.865 - r * (0.865 / 3.), y1 = y2 - 0.865 / 3.;
                MkPad(q / 3., y1, (q + 1) / 3., y2, false, 0.17, 0.03, 0.11, 0.12);
                std::vector<std::vector<VE>> x(4, std::vector<VE>(nv));
                double vmax = 0;
                for (int ct = 0; ct < 4; ++ct) {
                    for (size_t o = 0; o < ords.size(); ++o) {
                        const Seq& sq = *ords[o].second;
                        int j = 0;
                        for (int jj = 0; jj < 4; ++jj)
                            if (sq.cuts[jj] == ct) j = jj + 1;
                        x[ct][o] = Removed(A, L, k, sq.sel[j - 1], sq.sel[j], 0);
                    }
                    x[ct][nv - 2] = Removed(A, L, k, 0, 1 << ct, 0);
                    x[ct][nv - 1] = Removed(A, L, k, 15 & ~(1 << ct), 15, 0);
                    for (int v = 0; v < nv; ++v) vmax = std::max(vmax, x[ct][v].v + x[ct][v].e);
                }
                TH1D* fr = new TH1D(Nm("frA").c_str(), "", 4, 0, 4);
                fr->SetDirectory(nullptr);
                Axes(fr, "", "removed fraction", 0.056, 1.35);
                for (int ct = 0; ct < 4; ++ct) fr->GetXaxis()->SetBinLabel(ct + 1, cl[ct]);
                fr->GetXaxis()->SetLabelSize(0.080);
                fr->GetXaxis()->LabelsOption("h");
                fr->SetMinimum(0);
                fr->SetMaximum(std::max(0.05, vmax * 1.32));
                fr->Draw("axis");
                TGraphErrors* ge = new TGraphErrors();
                const double bw = 0.8 / nv;
                for (int ct = 0; ct < 4; ++ct)
                    for (int v = 0; v < nv; ++v) {
                        const double x0 = ct + 0.1 + bw * v, x1 = x0 + 0.92 * bw;
                        TBox* bx = new TBox(x0, 0, x1, std::max(0., x[ct][v].v));
                        bx->SetFillColor(vc[v]);
                        bx->SetLineColor(vc[v]);
                        bx->Draw();
                        const int p = ge->GetN();
                        ge->SetPoint(p, 0.5 * (x0 + x1), x[ct][v].v);
                        ge->SetPointError(p, 0, x[ct][v].e);
                    }
                ge->SetLineColor(HexCol("#E6550D"));
                ge->SetLineWidth(3);
                ge->SetMarkerSize(0);
                ge->Draw("P same");
                gPad->RedrawAxis();
                Txt(0.18, 0.92, Fmt("(%c) %s, %s", 'a' + r * 3 + q, pn[r][q], wn[r]), 0.058);
                const VE kf = Kept(A, L, k, 15, 0);
                Txt(0.36, 0.81, Fmt("kept (all four): %.3f #pm %.3f", kf.v, kf.e), 0.050);
            }
        SaveFig(c, "s5_spec_attribution");
    }

    // ---- 8. purity and S/sqrt(S+B)
    {
        TCanvas* c = new TCanvas("cPur", "", 2800, 2500);
        HeaderPad(c, 0.86);
        Txt(0.010, 0.84, "The resulting shape: purity per M_{4} window after each step, multi-pion share, and S/#sqrt{S+B} of the peak and the bump", 0.100);
        Txt(0.010, 0.68, "Purity = #sigma(group)/#sigma(T1+T2) in the window. Filled markers, solid: nominal order (M_{#pi#pi} < 0.46); open markers, dashed: M_{#pi#pi} < 0.42 at the last two steps.", 0.066);
        Txt(0.010, 0.56, Fmt("(d): counts = #sigma #times L_{data}, %s, a scale only (F #propto #sqrt{L_{data}}); B = rest of T1+T2 in the window; dashed: B + T3 + T4. Errors: simulation statistics.", kLdTex), 0.066);
        TLegend* lg = new TLegend(0.01, 0.0, 0.99, 0.48);
        lg->SetNColumns(6);
        lg->SetBorderSize(0);
        lg->SetFillStyle(0);
        lg->SetTextSize(0.078);
        lg->SetTextFont(42);
        for (int k : kWinFig) {
            TH1D* dm = Dummy(HexCol(kWinCol[k]), 4, 1, -1);
            dm->SetMarkerColor(HexCol(kWinCol[k]));
            dm->SetMarkerStyle(20);
            dm->SetMarkerSize(2.4);
            lg->AddEntry(dm, k == wBump ? "bump window" : k == wPeak ? "peak window" : kWinTxt[k], "lp");
        }
        lg->Draw();
        const int n = (int)qNom.sel.size();
        auto frame = [&](const char* yt, double lo, double hi, bool logy) {
            TH1D* fr = new TH1D(Nm("frP").c_str(), "", n, 0, n);
            fr->SetDirectory(nullptr);
            Axes(fr, "", yt, 0.055, 1.30);
            for (int j = 0; j < n; ++j) fr->GetXaxis()->SetBinLabel(j + 1, qNom.stepTex[j].c_str());
            fr->GetXaxis()->SetLabelSize(0.072);
            fr->GetXaxis()->LabelsOption("h");
            fr->SetMinimum(lo);
            fr->SetMaximum(hi);
            gPad->SetLogy(logy);
            fr->Draw("axis");
        };
        auto graph = [&](std::function<VE(int)> f, const Seq& q, int j0, int col, int mst, int ls, double off) {
            TGraphErrors* g = new TGraphErrors();
            for (int j = j0; j < (int)q.sel.size(); ++j) {
                const VE e = f(q.sel[j]);
                if (!(e.v > 0)) continue;
                const int p = g->GetN();
                g->SetPoint(p, j + 0.5 + off, e.v);
                g->SetPointError(p, 0, e.e);
            }
            g->SetLineColor(col);
            g->SetMarkerColor(col);
            g->SetLineWidth(3);
            g->SetLineStyle(ls);
            g->SetMarkerStyle(mst);
            g->SetMarkerSize(2.0);
            g->Draw("LP same");
        };
        const char* ptitle[3] = {"(a) purity of the #eta contributions, classes 2-5", "(b) purity of the peak, classes 4+5",
                                 "(c) multi-pion share, G_{#pi^{0}D}+G_{conv}+G_{vec}"};
        const unsigned pnum[3] = {Mem(cSig), Mem(cPeak), Mem(cMulti)};
        for (int p = 0; p < 3; ++p) {
            c->cd();
            const int row = p / 2, col = p % 2;
            const double y2 = 0.86 - row * 0.43, y1 = y2 - 0.43;
            MkPad(col * 0.5, y1, col * 0.5 + 0.5, y2, p < 2, 0.15, 0.03, 0.10, 0.12);
            if (p < 2) frame(p == 0 ? "#sigma(2-5) / #sigma(T1+T2)" : "#sigma(4+5) / #sigma(T1+T2)", 1e-4, 3.0, true);
            else frame("#sigma(multi-pion) / #sigma(T1+T2)", 0, 1.18, false);
            if (p == 2) {
                TLine* l = new TLine(0, 0.5, n, 0.5);
                l->SetLineStyle(2);
                l->SetLineWidth(2);
                l->Draw();
            }
            int idx = 0;
            for (int k : kWinFig) {
                ++idx;
                const unsigned num = pnum[p];
                auto f = [&, k, num](int s) { return Ratio(A, s, k, num, Mem(cT12)); };
                graph(f, qNom, 0, HexCol(kWinCol[k]), 20, 1, (idx - 3.5) * 0.02);
                graph(f, q42, 3, HexCol(kWinCol[k]), 24, 2, (idx - 3.5) * 0.02 + 0.06);
            }
            Txt(0.17, 0.92, ptitle[p], 0.056);
            if (p < 2) Txt(0.18, 0.835, "windows above 0.70 GeV/c^{2}: < 10^{-6}, not drawn (see table)", 0.044);
            gPad->RedrawAxis();
        }
        c->cd();
        MkPad(0.5, 0.0, 1.0, 0.43, false, 0.15, 0.03, 0.10, 0.12);
        double fmax = 0;
        for (const Seq* q : {&qNom, &q42})
            for (int s : q->sel)
                for (int w = 0; w < 2; ++w) {
                    const Sig f = w == 0 ? Signif(A, s, wPeak, Mem(cPeak), Mem(cBPk)) : Signif(A, s, wBump, Mem(cBump), Mem(cBBu));
                    fmax = std::max(fmax, f.F + f.eF);
                }
        frame("S/#sqrt{S+B} at L_{data}", 0, fmax * 1.6, false);
        for (int w = 0; w < 2; ++w)
            for (int v = 0; v < 2; ++v) {
                const int k = w == 0 ? wPeak : wBump;
                const unsigned sg = w == 0 ? Mem(cPeak) : Mem(cBump);
                const unsigned bg = w == 0 ? (v == 0 ? Mem(cBPk) : Mem(cBPkV)) : (v == 0 ? Mem(cBBu) : Mem(cBBuV));
                auto f = [&, k, sg, bg](int s) {
                    const Sig x = Signif(A, s, k, sg, bg);
                    VE e;
                    e.v = x.F;
                    e.e = x.eF;
                    return e;
                };
                const int col = HexCol(w == 0 ? kWinCol[wPeak] : kWinCol[wBump]);
                graph(f, qNom, 0, col, v == 0 ? 20 : 22, v == 0 ? 1 : 2, (w - 0.5) * 0.04);
                graph(f, q42, 3, col, v == 0 ? 24 : 26, v == 0 ? 1 : 2, (w - 0.5) * 0.04 + 0.06);
            }
        Txt(0.17, 0.92, "(d) S/#sqrt{S+B} at L_{data} (UNCONFIRMED)", 0.056);
        Txt(0.40, 0.84, "dark red: peak 4+5, peak window", 0.046, HexCol(kWinCol[wPeak]));
        Txt(0.40, 0.785, "orange: bump 2+3, bump window", 0.046, HexCol(kWinCol[wBump]));
        Txt(0.40, 0.73, "dashed: B + T3 + T4; open: M_{#pi#pi} < 0.42", 0.046);
        gPad->RedrawAxis();
        SaveFig(c, "s5_spec_purity");
    }

    // ---- 9. presets
    {
        const int ps[4] = {15, Id42(14), sPeakNom, sBumpNom};
        const char* pt[4] = {"OA > 2#circ, #alpha_{lab} < 50#circ, M_{#pi#pi} < 0.46, #alpha_{#eta} > 140#circ",
                             "#alpha_{lab} < 50#circ, M_{#pi#pi} < 0.42, #alpha_{#eta} > 140#circ", "#alpha_{lab} < 55#circ, M_{#pi#pi} > 0.30",
                             "#alpha_{lab} < 40#circ, 0.30 < M_{#pi#pi} < 0.41"};
        TCanvas* c = new TCanvas("cPre", "", 4000, 2500);
        HeaderPad(c, 0.83);
        Txt(0.008, 0.87, "Named selections (all = pre + the listed cuts, M_{#pi#pi} in GeV/c^{2}): simulated d#sigma/dM_{4}, log 0.28-2.50 (top, same axes) and linear zoom 0.40-0.70 (bottom, own y range)", 0.085);
        Txt(0.008, 0.755, NoteSim(0), 0.088);
        Txt(0.008, 0.66, NoteSim(1), 0.066);
        LegendStack(0.005, 0.0, 0.995, 0.66, 0.080, 4);
        double yhi = 0;
        for (int i = 0; i < 4; ++i) yhi = std::max(yhi, HMax(A, ps[i], 0));
        for (int i = 0; i < 4; ++i) {
            const int s = ps[i];
            for (int r = 0; r < 2; ++r) {
                c->cd();
                const double y2 = r == 0 ? 0.83 : 0.415, y1 = r == 0 ? 0.415 : 0.0;
                MkPad(i * 0.25, y1, (i + 1) * 0.25, y2, r == 0, 0.17, 0.045, 0.12, 0.14);
                if (r == 0) DrawSpectrum(A, s, 0, true, kYlo, yhi * 1e3, 0.056, 1.40, 505);
                else DrawSpectrum(A, s, 1, false, 0, HMax(A, s, 1) * 1.55, 0.056, 1.40);
                Txt(0.17, 0.955, Fmt("(%c) %s", 'a' + i + 4 * r, i == 0 ? "slides" : i == 1 ? "old code" : i == 2 ? "peak nominal" : "bump nominal"), 0.058);
                Txt(0.17, 0.90, pt[i], 0.046);
                if (r == 0) {
                    Txt(0.20, 0.80, Fmt("#sigma_{T1+T2} = %.4g nb", A.Xs(s, wAll, cT12)), 0.050);
                    Txt(0.20, 0.735, Fmt("multi-pion share, M_{4} > 0.70: %.2f", (A.Sum(s, wHigh, cMulti) + A.Sum(s, wTop, cMulti)) /
                                                                                  std::max(1e-300, A.Sum(s, wHigh, cT12) + A.Sum(s, wTop, cT12))), 0.050);
                } else {
                    const VE pk = Ratio(A, s, wPeak, Mem(cPeak), Mem(cT12)), pb = Ratio(A, s, wBump, Mem(cBump), Mem(cT12));
                    Txt(0.20, 0.80, Fmt("peak: #sigma_{4+5} %.3g nb, P %.3f", A.Xs(s, wPeak, cPeak), pk.v), 0.050);
                    Txt(0.20, 0.735, Fmt("bump: #sigma_{2+3} %.3g nb, P %.3f", A.Xs(s, wBump, cBump), pb.v), 0.050);
                }
            }
        }
        SaveFig(c, "s5_spec_presets");
    }
    printf("s5_spec: done in %.0f s\n", swAll.RealTime());
}
