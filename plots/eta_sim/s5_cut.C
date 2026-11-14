// ============================================================================
// s5_cut.C  --  cut flow on pure signal and on the simulated backgrounds
//               (report section 5, steps 3 and 4: what each cut does to each
//               signal class and to each background group)
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SIMULATION ONLY: no measured data are read.
//
// Rows: std + lep + eta sets, owned in partition mode A (SimOpt::sets = 7), weight w (PDG) times the physics
// weights of the eta decays, main variant V3 (simphys::kPhysMain = 108; sim_eta_physics.h).
// sigma = sum w / L_sim, L_sim = 221.05 nb^-1. No combinatorial background exists in simulation and none is added.
// Simulated sum = T1 (signal classes 2-5) + T2 (physical background). T3 (fake lepton) and T4 (leptons from two
// different sources) are shown separately (dashed) and enter a background only in the stated variant T2+T3+T4.
// X17 (class 6, eta set) is a hypothesis: W = 2 BR per row with BR = 1e-5 (illustration), never in a sum.
//
// Groups (the same for every agent of the study):
//   S2, S3, S4, S5 = signal classes 2 (bump, pi0 Dalitz), 3 (bump, conversion of a pi0 photon), 4 (peak, conversion
//   of the eta -> pi+pi-gamma photon), 5 (peak, direct eta -> pi+pi-e+e-); X17 = class 6;
//   G_eta  = T2 rows with both pions from one eta (pion-pair class P2) or an eta Dalitz pair (L2);
//   G_pi0D = other T2 rows with a pi0 Dalitz pair (L1); G_conv = other T2 rows with a conversion pair (L3);
//   G_vec  = other T2 rows with rho0/omega/phi -> e+e-, omega -> pi0 e+e- or a second-decay pair (L4-L6);
//   G_oth  = all other T2 rows; G_fake = T3; G_diff = T4.
//
// Base "pre" = SimEngine preset kPre (isBest >= 0, vz > -500 mm, PT3 bit, vz > -160 mm, cor_ok). Cuts, from the
// corrected momenta, in the order of the sequential flow:
//   OA(e+e-) > 2 deg, alpha_lab < 50 deg, M(pi+pi-) < 0.46 GeV/c2, alpha_eta > 140 deg   (these four = "slides"),
//   both pions in the vertex fit, |z_pi - v_z| < 40 mm for both pions, MM4 < 2.85 GeV/c2, pT4 > 400 MeV/c
//   (the last two: the values of the earlier eta studies, sim_eta_enriched.C and sim_eta_summary.C).
// M4 windows (fixed for all agents): peak 0.5333-0.5587, bump 0.4576-0.5255 GeV/c2.
//
// Parts (bits of the argument):
//   1  N-1 distributions of the eight cut variables in the bump and in the peak window ("slides" without the cut,
//      or "slides" for the four other variables): normalised shapes and absolute dsigma/dx;
//      figures figs/s5_cut_nm1_<var>, table tables/s5_cut_nm1.txt
//   2  efficiencies per class and group: each cut alone from pre, sequential (step and cumulative), N-1 (slides) and
//      N-1 (all eight cuts), with event-group errors; tables/s5_cut_eff.txt, figures s5_cut_eff_{bump,peak}, s5_cut_flow
//   4  composition per M4 region after pre and after slides (cross section and share per group, lepton-pair and
//      pion-pair classes); tables/s5_cut_comp.txt, figures s5_cut_comp, s5_cut_m4
//   8  ranking: Q = eps_S / sqrt(eps_B) and eps_S per cut for the peak, the bump and X17; what each cut removes;
//      expected counts at L_data = 6080.5 nb^-1 (UNCONFIRMED, scale only); tables/s5_cut_rank.txt, figure s5_cut_rank
// Always: checks against the engine (SumMany, Efficiency) and the earlier tables -> tables/s5_cut_checks.txt.
// The optional second argument names a Markdown file for the key tables of the report.
//
// Errors: event groups (rows of one event correlated; pass and fail rows correlated in an efficiency); std, lep and
// eta sets independent. All sums in double precision; histograms TH1D filled from the double sums.
// Figures: vector PDF (EPS -> ps2pdf -dEPSCrop) and PNG.
//
// Run from /home/przygoda/HADES/FAT/FAT (about 3 min, about 9 GB of memory, one process):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_cut", kTRUE)' \
//        'plots/eta_sim/s5_cut.C+(15)'
// ============================================================================
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace s5c {
using namespace simeng;

// ============================================================================
// constants: windows, regions
// ============================================================================
const double kBuLo = 0.4576, kBuHi = 0.5255, kPkLo = 0.5333, kPkHi = 0.5587;
const double kBR17 = 1e-5;          // illustrative BR(eta -> pi+pi- X17, X17 -> e+e-)
const double kLd = 6080.5;          // nb^-1, L_data, UNCONFIRMED, used as a scale only

// windows of the efficiency accumulators: the two analysis windows, all M4, and the other M4 regions
enum { wBump = 0, wPeak = 1, wAll = 2, wLow = 3, wGap = 4, wR5 = 5, wR6 = 6, wR7 = 7, kNW = 8 };
const double kWLo[kNW] = {kBuLo, kPkLo, -1e30, -1e30, kBuHi, kPkHi, 0.70, 1.00};
const double kWHi[kNW] = {kBuHi, kPkHi, 1e30, kBuLo, kPkLo, 0.70, 1.00, 1e30};
const char* const kWTxt[kNW] = {"bump window 0.4576-0.5255", "peak window 0.5333-0.5587", "all M4", "M4 < 0.4576",
                                "gap 0.5255-0.5333", "M4 0.5587-0.70", "M4 0.70-1.00", "M4 > 1.00"};
const char* const kWTag[kNW] = {"bump", "peak", "all", "low", "gap", "r5", "r6", "r7"};
// the seven M4 regions in mass order (window index)
const int kRegW[7] = {wLow, wBump, wGap, wPeak, wR5, wR6, wR7};

const int kNR = 7;                  // M4 regions; index kNR = all M4
const double kRLo[kNR] = {-1e30, kBuLo, kBuHi, kPkLo, kPkHi, 0.70, 1.00};
const double kRHi[kNR] = {kBuLo, kBuHi, kPkLo, kPkHi, 0.70, 1.00, 1e30};
const char* const kRName[kNR + 1] = {"< 0.4576", "bump 0.4576-0.5255", "gap 0.5255-0.5333", "peak 0.5333-0.5587",
                                     "0.5587-0.70", "0.70-1.00", "> 1.00", "all M4"};
const char* const kRShort[kNR + 1] = {"< 0.458", "bump", "gap", "peak", "0.559-0.70", "0.70-1.00", "> 1.00", "all M4"};

// ============================================================================
// components: 12 groups + aggregates
// ============================================================================
enum { gS2 = 0, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, kNG,
       aSB = kNG, aSP, aSig, aT2, aSim, aT234, aBX, kNC };
const char* const kCName[kNC] = {
    "S2 bump, pi0 Dalitz (cls 2)", "S3 bump, conversion (cls 3)", "S4 peak, pipigamma conv. (cls 4)",
    "S5 peak, direct pipiee (cls 5)", "X17 (cls 6), BR = 1e-5", "G_eta (P2 or L2)", "G_pi0D (L1)", "G_conv (L3)",
    "G_vec (L4-L6)", "G_oth (other T2)", "G_fake (T3)", "G_diff (T4)", "bump signal S2+S3", "peak signal S4+S5",
    "signal T1 = S2..S5", "phys. background T2", "simulated sum T1+T2", "bkg variant T2+T3+T4", "X17 bkg T2+S4+S5"};
const char* const kCShort[kNC] = {"S2", "S3", "S4", "S5", "X17", "G_eta", "G_pi0D", "G_conv", "G_vec", "G_oth",
                                  "G_fake", "G_diff", "S2+S3", "S4+S5", "T1", "T2", "T1+T2", "T2+T3+T4", "T2+S4+S5"};
const char* const kCLeg[kNC] = {
    "S2 bump: #pi^{0} Dalitz", "S3 bump: #pi^{0}, #gamma conversion", "S4 peak: #pi^{+}#pi^{-}#gamma, conversion",
    "S5 peak: #pi^{+}#pi^{-}e^{+}e^{-} direct", "X17 (BR 10^{-5}, not summed)", "G_{#eta}: #eta made, not signal",
    "G_{#pi^{0}D}: #pi^{0} Dalitz pair (L1)", "G_{conv}: conversion pair (L3)", "G_{vec}: #rho#omega#phi #rightarrow ee (L4-6)",
    "G_{oth}: #eta-set ee (L10), other #pi", "G_{fake}: fake lepton (T3)", "G_{diff}: two sources (T4)",
    "S2 + S3", "S4 + S5", "T1", "T2 total (phys. bkg)", "T1 + T2 (simulated sum)", "T2 + T3 + T4", "T2 + S4 + S5"};
// warm colours for signals, cold / grey for backgrounds
const char* const kHex[kNC] = {"#F57C00", "#8D4E1A", "#D62728", "#C2185B", "#7B0000", "#5E35B1", "#1565C0",
                               "#00ACC1", "#00796B", "#78909C", "#9E9E9E", "#424242", "#F57C00", "#D62728",
                               "#D62728", "#000000", "#000000", "#455A64", "#455A64"};
int Col(int c) { return TColor::GetColor(kHex[c]); }
int LSty(int c) { return c == gX17 ? 2 : c == gFake ? 7 : c == gDiff ? 3 : 1; }
int LWid(int c) { return (c == aT2 || c == aSim) ? 5 : 4; }
int MSty(int c) {
    switch (c) {
        case gS2: return 20; case gS3: return 21; case gS4: return 22; case gS5: return 23; case gX17: return 29;
        case gEta: return 24; case gPi0D: return 25; case gConv: return 26; case gVec: return 32; case gOth: return 28;
        case gFake: return 27; case gDiff: return 30; case aT2: return 34; case aSim: return 34;
    }
    return 20;
}
bool IsSignal(int c) { return c <= gX17 || c == aSB || c == aSP || c == aSig; }

// ============================================================================
// cuts
// ============================================================================
enum { cOA = 0, cLab, cMpp, cEta, cVtx, cDz, cMm4, cPt4, kNCut };
const char* const kCutName[kNCut] = {"OA > 2", "alpha_lab < 50", "M(pipi) < 0.46", "alpha_eta > 140", "pi in vtx fit",
                                     "abs(z_pi-vz) < 40", "MM4 < 2.85", "pT4 > 400"};
const char* const kCutTex[kNCut] = {"OA > 2#circ", "#alpha_{lab} < 50#circ", "M_{#pi#pi} < 0.46", "#alpha_{#eta} > 140#circ",
                                    "both #pi in vertex fit", "|z_{#pi} - v_{z}| < 40 mm", "MM_{4} < 2.85", "p_{T,4} > 400 MeV/c"};
const char* const kCutShort[kNCut] = {"OA", "#alpha_{lab}", "M_{#pi#pi}", "#alpha_{#eta}", "#pi vtx", "|#Deltaz_{#pi}|",
                                      "MM_{4}", "p_{T,4}"};
const unsigned kSlM = 0x0F, kAllC = 0xFF;
inline unsigned Bit(int c) { return 1u << c; }

// ============================================================================
// compact rows
// ============================================================================
struct MRow {
    double w;
    float m4, oa, lab, mpp, eta, mm4, pt4, dz;
    int gid;
    unsigned mask;
    unsigned short cb;
    unsigned char g, nv, set;
    signed char lcls, pcls, tmpl;
};
std::vector<MRow> gR;

unsigned MaskOf(int g) {
    unsigned m = 1u << g;
    if (g == gS2 || g == gS3) m |= (1u << aSB) | (1u << aSig) | (1u << aSim);
    else if (g == gS4 || g == gS5) m |= (1u << aSP) | (1u << aSig) | (1u << aSim) | (1u << aBX);
    else if (g >= gEta && g <= gOth) m |= (1u << aT2) | (1u << aSim) | (1u << aT234) | (1u << aBX);
    else if (g == gFake || g == gDiff) m |= (1u << aT234);
    return m;
}
int GroupOf(const Row& r) {
    if (r.tmpl == 1) {
        if (r.cls >= 2 && r.cls <= 5) return gS2 + (r.cls - 2);
        if (r.cls == 6) return gX17;
        return -1;
    }
    if (r.tmpl == 2) {
        if (r.pcls == 2 || r.lcls == 2) return gEta;
        if (r.lcls == 1) return gPi0D;
        if (r.lcls == 3) return gConv;
        if (r.lcls >= 4 && r.lcls <= 6) return gVec;
        return gOth;
    }
    if (r.tmpl == 3) return gFake;
    if (r.tmpl == 4) return gDiff;
    return -1;
}
unsigned short CutBits(const MRow& r) {
    unsigned short b = 0;
    if (r.oa > 2.f) b |= 1;
    if (r.lab < 50.f) b |= 2;
    if (r.mpp < 0.46f) b |= 4;
    if (r.eta > 140.f) b |= 8;
    if (r.nv == 2) b |= 16;
    if (r.dz < 40.f) b |= 32;
    if (r.mm4 < 2.85f) b |= 64;
    if (r.pt4 > 400.f) b |= 128;
    return b;
}

// ============================================================================
// output helpers
// ============================================================================
FILE* gOut = nullptr;
FILE* gMd = nullptr;
void P(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt); vfprintf(gOut ? gOut : stdout, fmt, ap); va_end(ap);
}
void M(const char* fmt, ...) {
    if (!gMd) return;
    va_list ap; va_start(ap, fmt); vfprintf(gMd, fmt, ap); va_end(ap);
}
std::string Fs(const char* fmt, ...) {
    char buf[8192];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    return buf;
}
FILE* OpenTab(const char* name) {
    FILE* f = fopen(Form("%s/%s", kTabDir, name), "w");
    if (!f) printf("cannot open %s/%s\n", kTabDir, name);
    return f;
}
void CloseTab() { if (gOut) fclose(gOut); gOut = nullptr; }
// value +- error with a precision that follows the error
std::string VE(double v, double e, int minDig = 3) {
    if (!(v == v)) return "-";
    int d = minDig;
    if (e > 0) { const int k = (int)std::ceil(-std::log10(e)) + 1; d = std::max(minDig, std::min(6, k)); }
    return Fs("%.*f +- %.*f", d, v, d, e);
}
std::string VEg(double v, double e) { return Fs("%.4g +- %.2g", v, e); }

// ============================================================================
// weighted sums with event-group variance (pass and fail rows of one group correlated)
// ============================================================================
struct Acc {
    double Sa = 0, Sp = 0, aa = 0, ab = 0, bb = 0, Na = 0, Np = 0;
    double eps() const { return Sa > 0 ? Sp / Sa : 0.; }
    double err() const {
        if (Sa <= 0) return 0.;
        const double e = Sp / Sa, v = (1 - e) * (1 - e) * aa - 2 * e * (1 - e) * ab + e * e * bb;
        return v > 0 ? std::sqrt(v) / Sa : 0.;
    }
    double Va() const { return aa + 2 * ab + bb; }
    double NeffA() const { const double v = Va(); return v > 0 ? Sa * Sa / v : 0.; }
    double NeffP() const { return aa > 0 ? Sp * Sp / aa : 0.; }
    double xsA() const { return Sa / kLsim; }
    double xsAe() const { return std::sqrt(std::max(0., Va())) / kLsim; }
    double xsP() const { return Sp / kLsim; }
    double xsPe() const { return std::sqrt(std::max(0., aa)) / kLsim; }
};

// efficiency configurations: type x cut, plus "slides from pre"
enum { tAlone = 0, tStep, tCum, tNm1S, tNm1A, kNTy };
const char* const kTyName[kNTy] = {"each cut alone, from pre", "sequential step (cut given all previous cuts)",
                                   "sequential, cumulative from pre", "N-1 (slides): cut given the other slides cuts; "
                                   "last four: on top of slides", "N-1 (all eight): cut given the other seven"};
const int kNCfg = kNTy * kNCut + 1;
const int kCfgSlides = kNTy * kNCut;
unsigned gBase[kNCfg], gPass[kNCfg];
inline int Cfg(int t, int c) { return t * kNCut + c; }
void InitCfg() {
    for (int c = 0; c < kNCut; ++c) {
        const unsigned b = Bit(c), prev = b - 1;
        gBase[Cfg(tAlone, c)] = 0;                               gPass[Cfg(tAlone, c)] = b;
        gBase[Cfg(tStep, c)] = prev;                             gPass[Cfg(tStep, c)] = b;
        gBase[Cfg(tCum, c)] = 0;                                 gPass[Cfg(tCum, c)] = prev | b;
        gBase[Cfg(tNm1S, c)] = c < 4 ? (kSlM & ~b) : kSlM; gPass[Cfg(tNm1S, c)] = b;
        gBase[Cfg(tNm1A, c)] = kAllC & ~b;                       gPass[Cfg(tNm1A, c)] = b;
    }
    gBase[kCfgSlides] = 0; gPass[kCfgSlides] = kSlM;
}
std::vector<Acc> gE;   // [window][cfg][component]
inline Acc& E(int w, int k, int c) { return gE[((size_t)w * kNCfg + k) * kNC + c]; }

void EvalEff() {
    TStopwatch sw;
    gE.assign((size_t)kNW * kNCfg * kNC, Acc());
    const size_t n = gE.size();
    std::vector<double> ga(n, 0.), gb(n, 0.);
    std::vector<char> tc(n, 0);
    std::vector<int> touched;
    touched.reserve(4096);
    int cur = INT_MIN;
    for (size_t i = 0; i <= gR.size(); ++i) {
        const bool end = i == gR.size();
        if (end || gR[i].gid != cur) {
            for (int idx : touched) {
                Acc& a = gE[idx];
                a.aa += ga[idx] * ga[idx]; a.ab += ga[idx] * gb[idx]; a.bb += gb[idx] * gb[idx];
                ga[idx] = gb[idx] = 0; tc[idx] = 0;
            }
            touched.clear();
            if (end) break;
            cur = gR[i].gid;
        }
        const MRow& r = gR[i];
        int wins[2], nw = 0;
        for (int k = 0; k < 7; ++k) if (r.m4 >= kWLo[kRegW[k]] && r.m4 < kWHi[kRegW[k]]) { wins[nw++] = kRegW[k]; break; }
        wins[nw++] = wAll;
        for (int iw = 0; iw < nw; ++iw) {
            const int w = wins[iw];
            for (int k = 0; k < kNCfg; ++k) {
                if ((r.cb & gBase[k]) != gBase[k]) continue;
                const bool p = (r.cb & gPass[k]) == gPass[k];
                unsigned m = r.mask;
                while (m) {
                    const int c = __builtin_ctz(m);
                    m &= m - 1;
                    const int idx = (int)(((size_t)w * kNCfg + k) * kNC + c);
                    Acc& a = gE[idx];
                    a.Sa += r.w; a.Na += 1;
                    if (p) { a.Sp += r.w; a.Np += 1; ga[idx] += r.w; } else gb[idx] += r.w;
                    if (!tc[idx]) { tc[idx] = 1; touched.push_back(idx); }
                }
            }
        }
    }
    printf("s5c::EvalEff: %.1f s\n", sw.RealTime());
}

// ============================================================================
// weighted histograms with event-group variance
// ============================================================================
struct HS {
    std::vector<double> s, v;   // bins 0..nb-1, nb = underflow, nb+1 = overflow
    double S = 0, V = 0, N = 0, Sx = 0;
};
struct Job {
    int var = 0;
    unsigned base = 0;
    double mlo = -1e30, mhi = 1e30;
    int nb = 0;
    double lo = 0, hi = 1;
    std::vector<HS> h;
    std::vector<std::vector<std::pair<int, double>>> buf;
    void Init(int n, double a, double b) {
        nb = n; lo = a; hi = b;
        h.assign(kNC, HS());
        for (HS& x : h) { x.s.assign(nb + 2, 0.); x.v.assign(nb + 2, 0.); }
        buf.assign(kNC, std::vector<std::pair<int, double>>());
    }
    int Bin(double x) const {
        if (!(x == x)) return nb + 1;
        if (x < lo) return nb;
        if (x >= hi) return nb + 1;
        return std::min(nb - 1, (int)((x - lo) / (hi - lo) * nb));
    }
    double bw() const { return (hi - lo) / nb; }
    void Flush() {
        for (int c = 0; c < kNC; ++c) {
            std::vector<std::pair<int, double>>& e = buf[c];
            if (e.empty()) continue;
            std::sort(e.begin(), e.end());
            double all = 0;
            size_t i = 0;
            while (i < e.size()) {
                const int b = e[i].first;
                double acc = 0;
                while (i < e.size() && e[i].first == b) { acc += e[i].second; ++i; }
                h[c].v[b] += acc * acc;
                all += acc;
            }
            h[c].V += all * all;
            e.clear();
        }
    }
    double Neff(int c) const { return h[c].V > 0 ? h[c].S * h[c].S / h[c].V : 0.; }
    // weighted quantile from the histogram (linear inside the bin)
    double Quantile(int c, double q) const {
        const HS& x = h[c];
        if (x.S <= 0) return std::numeric_limits<double>::quiet_NaN();
        double t = q * x.S - x.s[nb];
        if (t < 0) return lo;
        for (int i = 0; i < nb; ++i) {
            if (t <= x.s[i] && x.s[i] > 0) return lo + (i + t / x.s[i]) * bw();
            t -= x.s[i];
        }
        return hi;
    }
};

const int kNVar = 9;   // eight cut variables + M4
struct VarDef {
    const char* tag; const char* title; const char* yShape; const char* yAbs;
    int nb; double lo, hi, cut; int keep;   // keep: +1 keeps x > cut, -1 keeps x < cut
    std::vector<double> alt;
};
std::vector<VarDef> gVar;
void InitVars() {
    gVar.clear();
    gVar.push_back({"oa", "OA(e^{+}e^{-}) [deg]", "1/#sigma d#sigma/dOA [1/deg]", "d#sigma/dOA [nb/deg]", 40, 0, 20, 2, +1, {}});
    gVar.push_back({"alab", "#alpha_{lab} [deg]", "1/#sigma d#sigma/d#alpha_{lab} [1/deg]", "d#sigma/d#alpha_{lab} [nb/deg]", 36, 0, 180, 50, -1, {40, 55}});
    gVar.push_back({"mpipi", "M(#pi^{+}#pi^{-}) [GeV/c^{2}]", "1/#sigma d#sigma/dM_{#pi#pi} [(GeV/c^{2})^{-1}]",
                    "d#sigma/dM_{#pi#pi} [nb/(GeV/c^{2})]", 30, 0.27, 0.57, 0.46, -1, {0.30, 0.42}});
    gVar.push_back({"aeta", "#alpha_{#eta} [deg]", "1/#sigma d#sigma/d#alpha_{#eta} [1/deg]", "d#sigma/d#alpha_{#eta} [nb/deg]", 30, 90, 180, 140, +1, {}});
    gVar.push_back({"vtx", "pions used in the vertex fit", "fraction of #sigma", "#sigma [nb]", 3, -0.5, 2.5, 1.5, +1, {}});
    gVar.push_back({"dzpi", "max |z_{#pi} - v_{z}| [mm]", "1/#sigma d#sigma/d|#Deltaz| [1/mm]", "d#sigma/d|#Deltaz| [nb/mm]", 40, 0, 120, 40, -1, {}});
    gVar.push_back({"mm4", "MM_{4} [GeV/c^{2}]", "1/#sigma d#sigma/dMM_{4} [(GeV/c^{2})^{-1}]", "d#sigma/dMM_{4} [nb/(GeV/c^{2})]", 32, 1.8, 3.4, 2.85, -1, {}});
    gVar.push_back({"pt4", "p_{T,4} [MeV/c]", "1/#sigma d#sigma/dp_{T,4} [(MeV/c)^{-1}]", "d#sigma/dp_{T,4} [nb/(MeV/c)]", 40, 0, 1000, 400, +1, {}});
    gVar.push_back({"m4", "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "", "d#sigma/dM_{4} [nb/(GeV/c^{2})]", 160, 0.30, 1.10, -1, 0, {}});
}
inline double VarVal(const MRow& r, int v) {
    switch (v) {
        case 0: return r.oa; case 1: return r.lab; case 2: return r.mpp; case 3: return r.eta; case 4: return r.nv;
        case 5: return r.dz; case 6: return r.mm4; case 7: return r.pt4;
    }
    return r.m4;
}

void RunJobs(std::vector<Job>& jobs) {
    TStopwatch sw;
    std::vector<char> tj(jobs.size(), 0);
    std::vector<int> touched;
    int cur = INT_MIN;
    for (size_t i = 0; i <= gR.size(); ++i) {
        const bool end = i == gR.size();
        if (end || gR[i].gid != cur) {
            for (int j : touched) { jobs[j].Flush(); tj[j] = 0; }
            touched.clear();
            if (end) break;
            cur = gR[i].gid;
        }
        const MRow& r = gR[i];
        for (size_t j = 0; j < jobs.size(); ++j) {
            Job& jb = jobs[j];
            if ((r.cb & jb.base) != jb.base) continue;
            if (!(r.m4 >= jb.mlo && r.m4 < jb.mhi)) continue;
            const double x = VarVal(r, jb.var);
            const int b = jb.Bin(x);
            unsigned m = r.mask;
            while (m) {
                const int c = __builtin_ctz(m);
                m &= m - 1;
                HS& h = jb.h[c];
                h.s[b] += r.w; h.S += r.w; h.N += 1;
                if (x == x) h.Sx += r.w * x;
                jb.buf[c].push_back(std::make_pair(b, r.w));
            }
            if (!tj[j]) { tj[j] = 1; touched.push_back((int)j); }
        }
    }
    printf("s5c::RunJobs: %zu jobs in %.1f s\n", jobs.size(), sw.RealTime());
}

// ============================================================================
// loading: compact rows with the weights of the main variant
// ============================================================================
double gOthL[12], gOthN[12];   // G_oth: sum w and rows per lepton-pair class (check of the group definition)
long long gNoGroup = 0, gZeroW = 0;
bool Build(SimData& d) {
    TStopwatch sw;
    SimOpt o; o.phys = simphys::kPhysMain;
    for (int k = 0; k < 12; ++k) gOthL[k] = gOthN[k] = 0;
    gR.clear();
    gR.reserve(8000000);
    int gid = 0;
    for (int set = 1; set <= 3; ++set) {
        const std::vector<Row>& v = d.tree(set, kUS);
        int cur = INT_MIN;
        for (const Row& r : v) {
            if (r.grp != cur) { ++gid; cur = r.grp; }
            if (!PassPre(r)) continue;
            if (!Owned(r, o)) continue;
            const int g = GroupOf(r);
            if (g < 0) { ++gNoGroup; continue; }
            const double w = g == gX17 ? 2. * kBR17 : Weight(r, o);
            if (!(w > 0)) { ++gZeroW; continue; }
            MRow m;
            m.w = w;
            m.m4 = r.m4; m.oa = r.oa; m.lab = r.oa_lab; m.mpp = r.mpipi; m.eta = r.oa_eta; m.mm4 = r.mm4; m.pt4 = r.pt4;
            m.dz = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
            m.nv = (unsigned char)((r.pip_isUsedVertex == 1) + (r.pim_isUsedVertex == 1));
            m.gid = gid;
            m.g = (unsigned char)g;
            m.mask = MaskOf(g);
            m.set = (unsigned char)set;
            m.lcls = (signed char)r.lcls; m.pcls = (signed char)r.pcls; m.tmpl = (signed char)r.tmpl;
            m.cb = CutBits(m);
            gR.push_back(m);
            if (g == gOth && r.lcls >= 0 && r.lcls < 12) { gOthL[r.lcls] += w; gOthN[r.lcls] += 1; }
        }
    }
    printf("s5c::Build: %zu rows (pre, owned mode A, w > 0), %lld without group, %lld with w = 0, %.1f s\n", gR.size(),
           gNoGroup, gZeroW, sw.RealTime());
    return !gR.empty();
}

// ============================================================================
// checks against the engine and the earlier tables
// ============================================================================
int gNchk = 0, gNbad = 0;
// column 'sigma main' [nb] of tables/rw_weights_xsec.txt (sim_rw_build.C): section "## <sec>", line "<comp> <win> | off +- e
// main +- e ..."; NaN if not found (the check then differs). Replaces the rules-v1 numbers written here before 2026-11-13.
double XsecMain(const std::string& sec, const std::string& comp, const std::string& win) {
    std::ifstream in(Form("%s/rw_weights_xsec.txt", kTabDir));
    std::string l;
    bool inSec = false;
    while (std::getline(in, l)) {
        if (l.compare(0, 3, "## ") == 0) { inSec = l.substr(3) == sec; continue; }
        if (!inSec || l.compare(0, comp.size(), comp) != 0) continue;
        const size_t w = l.find_first_not_of(' ', comp.size()), bar = l.find('|');
        if (w == std::string::npos || bar == std::string::npos || l.compare(w, win.size(), win) != 0) continue;
        double off = 0, eoff = 0, main = 0;
        if (sscanf(l.c_str() + bar + 1, "%lf +- %lf %lf", &off, &eoff, &main) == 3) return main;
    }
    printf("s5_cut: no '%s %s' in section '%s' of tables/rw_weights_xsec.txt\n", comp.c_str(), win.c_str(), sec.c_str());
    return std::numeric_limits<double>::quiet_NaN();
}
void Check(const char* what, double mine, double ref, double tol) {
    ++gNchk;
    const double d = mine - ref;
    const bool ok = std::fabs(d) <= tol;
    if (!ok) ++gNbad;
    P("  %-74s this %12.6g  ref %12.6g  diff %+9.2g  %s\n", what, mine, ref, d, ok ? "OK" : "DIFFERS");
}
void Checks(SimData& d) {
    gOut = OpenTab("s5_cut_checks.txt");
    P("# s5_cut_checks: checks of s5_cut.C (compact-row accumulators against the engine and the earlier tables)\n");
    P("# weights: std + lep + eta owned rows (mode A), w (PDG) x phys %d (%s); X17: W = 2 x %.0e per class-6 row\n\n",
      simphys::kPhysMain, simphys::PhysLabel(simphys::kPhysMain).c_str(), kBR17);
    SimOpt o; o.phys = simphys::kPhysMain;
    // 1. totals against simeng::SumMany (all M4)
    P("## 1. sigma [nb], all M4: this macro against simeng::SumMany (same rows, same weights)\n");
    std::vector<Comp> comps = {CompSigCls(2), CompSigCls(3), CompSigCls(4), CompSigCls(5), CompX17(kBR17),
                               CompPhysBkg(), CompFake(), CompDiff(), CompSim()};
    const int mine[9] = {gS2, gS3, gS4, gS5, gX17, aT2, gFake, gDiff, aSim};
    for (int s : {(int)simeng::kPre, (int)simeng::kSlides}) {
        const std::vector<WSum> sm = SumMany(d, SelPreset(s), comps, o);
        for (int i = 0; i < 9; ++i) {
            const Acc& a = s == kPre ? E(wAll, Cfg(tAlone, 0), mine[i]) : E(wAll, kCfgSlides, mine[i]);
            const double my = s == kPre ? a.xsA() : a.xsP(), mye = s == kPre ? a.xsAe() : a.xsPe();
            Check(Form("%s: %s sigma", PresetName(s), kCName[mine[i]]), my, sm[i].Xsec(), 1e-9 * std::max(1., sm[i].Xsec()));
            Check(Form("%s: %s error", PresetName(s), kCName[mine[i]]), mye, sm[i].XsecErr(), 1e-9 * std::max(1., sm[i].Xsec()));
        }
    }
    // 2. efficiencies against simeng::Efficiency
    P("\n## 2. efficiencies: this macro against simeng::Efficiency (event-group errors)\n");
    {
        Sel basePk = [](const Row& r) { return PassPre(r) && r.m4 >= kPkLo && r.m4 < kPkHi; };
        Sel baseBu = [](const Row& r) { return PassPre(r) && r.m4 >= kBuLo && r.m4 < kBuHi; };
        Sel baseX = [](const Row& r) {
            return PassPre(r) && r.m4 >= kPkLo && r.m4 < kPkHi && r.oa_lab < 50.f && r.mpipi < 0.46f && r.oa_eta > 140.f;
        };
        EffRes e1 = Efficiency(d, basePk, [](const Row& r) { return r.oa > 2.f; }, CompSigCls(4), o);
        const Acc& a1 = E(wPeak, Cfg(tAlone, cOA), gS4);
        Check("S4 in the peak window, OA > 2 alone: eps", a1.eps(), e1.eps, 1e-12);
        Check("S4 in the peak window, OA > 2 alone: error", a1.err(), e1.err, 1e-12);
        EffRes e2 = Efficiency(d, baseBu, [](const Row& r) { return r.oa_lab < 50.f; }, CompPhysBkg(), o);
        const Acc& a2 = E(wBump, Cfg(tAlone, cLab), aT2);
        Check("T2 in the bump window, alpha_lab < 50 alone: eps", a2.eps(), e2.eps, 1e-12);
        Check("T2 in the bump window, alpha_lab < 50 alone: error", a2.err(), e2.err, 1e-12);
        EffRes e3 = Efficiency(d, baseX, [](const Row& r) { return r.oa > 2.f; }, CompX17(kBR17), o);
        const Acc& a3 = E(wPeak, Cfg(tNm1S, cOA), gX17);
        Check("X17 in the peak window, OA > 2 N-1 (slides): eps", a3.eps(), e3.eps, 1e-12);
        Check("X17 in the peak window, OA > 2 N-1 (slides): error", a3.err(), e3.err, 1e-12);
        EffRes e4 = Efficiency(d, basePk, [](const Row& r) { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; },
                               CompSigCls(5), o);
        const Acc& a4 = E(wPeak, Cfg(tAlone, cVtx), gS5);
        Check("S5 in the peak window, both pions in the vertex fit alone: eps", a4.eps(), e4.eps, 1e-12);
        Check("S5 in the peak window, both pions in the vertex fit alone: error", a4.err(), e4.err, 1e-12);
    }
    // 3. earlier tables (printing precision)
    P("\n## 3. against the earlier tables (tolerance = 0.6 x last printed digit)\n");
    P("# rw_weights_xsec.txt, column 'sigma main' (phys 108), sections pre and E1 (slides), all M4; eta_engine_reference.txt X17\n");
    // read from the table (rules v1 values written here before: pre 7.1591, 5.2907, 1.1952, 1.4168, 701.1347;
    // slides 3.9915, 2.0896, 0.5703, 1.0294, 68.0034)
    const char* xsComp[5] = {"bump pi0 Dalitz (2)", "bump conversion (3)", "peak conversion (4)", "peak direct (5)", "T1 + T2 (simulated sum)"};
    double refPre[5], refSl[5];
    for (int i = 0; i < 5; ++i) { refPre[i] = XsecMain("pre", xsComp[i], "all M4"); refSl[i] = XsecMain("E1 (slides)", xsComp[i], "all M4"); }
    const int cc[5] = {gS2, gS3, gS4, gS5, aSim};
    for (int i = 0; i < 5; ++i) Check(Form("pre, all M4: %s", kCName[cc[i]]), E(wAll, Cfg(tAlone, 0), cc[i]).xsA(), refPre[i], 6e-5);
    for (int i = 0; i < 5; ++i) Check(Form("slides, all M4: %s", kCName[cc[i]]), E(wAll, kCfgSlides, cc[i]).xsP(), refSl[i], 6e-5);
    Check("pre, all M4: X17 BR 1e-5", E(wAll, Cfg(tAlone, 0), gX17).xsA(), 0.07127, 6e-6);
    Check("slides, all M4: X17 BR 1e-5", E(wAll, kCfgSlides, gX17).xsP(), 0.0589, 6e-5);
    // 4. content of G_oth
    P("\n## 4. content of G_oth (T2 rows not in G_eta, L1, L3, L4-L6), pre, all M4: lepton-pair class\n");
    double tot = 0;
    for (int k = 0; k < 12; ++k) tot += gOthL[k];
    for (int k = 0; k < 12; ++k)
        if (gOthN[k] > 0) P("  %-28s rows %9.0f  sigma %.5f nb  share %.4f\n", LclsName(k), gOthN[k], gOthL[k] / kLsim, gOthL[k] / tot);
    Check("G_oth: share outside L10/L11", 1. - (gOthL[10] + gOthL[11]) / tot, 0., 1e-12);
    P("\n# rows dropped: %lld without a group (no template), %lld with W = 0\n", gNoGroup, gZeroW);
    P("# %d checks, %d outside the tolerance\n", gNchk, gNbad);
    CloseTab();
    printf("s5c::Checks: %d checks, %d outside the tolerance\n", gNchk, gNbad);
}

// ============================================================================
// plotting helpers
// ============================================================================
void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLegendFont(42); gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetEndErrorSize(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetCanvasColor(0); gStyle->SetPadColor(0); gStyle->SetFrameFillColor(0);
    gStyle->SetHatchesLineWidth(2);
    gStyle->SetLineScalePS(1.0);   // vector output: line widths as on screen
    TGaxis::SetMaxDigits(3);
}
void AxisStyle(TH1* h, double s = 0.055, double yoff = 1.3, double xoff = 1.1) {
    h->GetXaxis()->SetLabelSize(s); h->GetYaxis()->SetLabelSize(s);
    h->GetXaxis()->SetTitleSize(s * 1.05); h->GetYaxis()->SetTitleSize(s * 1.05);
    h->GetXaxis()->SetTitleOffset(xoff); h->GetYaxis()->SetTitleOffset(yoff);
    h->GetXaxis()->SetNdivisions(506); h->GetYaxis()->SetNdivisions(505);
    h->GetXaxis()->SetLabelOffset(0.012);
}
TLatex* Txt(double x, double y, const char* s, double size, int align = 11, int col = 1, int font = 42) {
    TLatex* t = new TLatex(x, y, s);
    t->SetNDC(); t->SetTextSize(size); t->SetTextAlign(align); t->SetTextColor(col); t->SetTextFont(font);
    t->Draw();
    return t;
}
TLatex* TxtPx(double x, double y, const char* s, double px, int align = 11, int col = 1) { return Txt(x, y, s, px, align, col, 43); }
void SaveFig(TCanvas* c, const char* name) {
    const TString base = Form("%s/%s", kFigDir, name);
    c->SaveAs(base + ".eps");
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf", base.Data(), base.Data()));
    gSystem->Unlink(base + ".eps");
    // PNG rendered from the PDF (ROOT's own PNG backend misplaces sub- and superscripts after #circ and ^{2})
    gSystem->Exec(Form("pdftoppm -png -singlefile -scale-to-x %d -scale-to-y -1 %s.pdf %s", c->GetWw(), base.Data(), base.Data()));
}
int& Uid() { static int n = 0; return n; }
TH1D* NewH(const char* pfx, int nb, double lo, double hi) {
    TH1D* h = new TH1D(Form("%s_%d", pfx, ++Uid()), "", nb, lo, hi);
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}
void LineStyle(TH1* h, int c) {
    h->SetLineColor(Col(c)); h->SetLineStyle(LSty(c)); h->SetLineWidth(LWid(c));
    h->SetMarkerColor(Col(c)); h->SetMarkerStyle(MSty(c)); h->SetMarkerSize(0);
    h->SetFillStyle(0);
}
// histogram of component c of a job: mode 0 = shape (1/sigma dsigma/dx, normalised to all rows of the component,
// under- and overflow included), 1 = dsigma/dx [nb / unit of x]; for the vertex flag (bw = 1) per bin
TH1D* MakeH(const Job& j, int c, int mode) {
    TH1D* h = NewH("h", j.nb, j.lo, j.hi);
    const HS& x = j.h[c];
    const double f = mode == 0 ? (x.S > 0 ? 1. / (x.S * j.bw()) : 0.) : 1. / (kLsim * j.bw());
    for (int i = 0; i < j.nb; ++i) { h->SetBinContent(i + 1, x.s[i] * f); h->SetBinError(i + 1, std::sqrt(x.v[i]) * f); }
    LineStyle(h, c);
    return h;
}
// frame y limits for a set of histograms
void YLim(const std::vector<TH1D*>& hh, bool logy, double& ymin, double& ymax) {
    ymax = 0; double mpos = 1e300;
    for (TH1D* h : hh)
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            const double y = h->GetBinContent(i);
            ymax = std::max(ymax, y);
            if (y > 0) mpos = std::min(mpos, y);
        }
    if (ymax <= 0) ymax = 1;
    if (logy) { ymin = std::max(mpos * 0.5, ymax * 3e-5); ymax *= 8; }
    else { ymin = 0; ymax *= 1.22; }
}
double YAt(double ymin, double ymax, bool logy, double f) {
    return logy ? std::exp(std::log(ymin) + f * (std::log(ymax) - std::log(ymin))) : ymin + f * (ymax - ymin);
}

// ============================================================================
// part 1: N-1 distributions
// ============================================================================
// components drawn per window
const std::vector<int> kDrawBump = {gS2, gS3, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2};
const std::vector<int> kDrawPeak = {gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2};
const std::vector<int> kLegAll = {gS2, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2};
const double kNeffMin = 10.;      // efficiency points and absolute curves
const double kNeffShape = 30.;    // normalised shapes

void DrawNm1(int v, const Job& jb, const Job& jp) {
    const VarDef& V = gVar[v];
    TCanvas* c = new TCanvas(Form("c_nm1_%s", V.tag), "", 2400, 1400);
    const double xL = 0.74;
    for (int w = 0; w < 2; ++w)
        for (int kind = 0; kind < 2; ++kind) {
            c->cd();
            const double x0 = kind == 0 ? 0. : xL / 2, x1 = kind == 0 ? xL / 2 : xL;
            const double y0 = w == 0 ? 0.5 : 0., y1 = w == 0 ? 1. : 0.5;
            TPad* p = new TPad(Form("p_nm1_%s_%d_%d", V.tag, w, kind), "", x0, y0, x1, y1);
            p->SetMargin(0.17, 0.05, 0.16, 0.10);
            p->SetLogy(kind == 1);
            p->Draw(); p->cd();
            const Job& j = w == 0 ? jb : jp;
            const std::vector<int>& cs = w == 0 ? kDrawBump : kDrawPeak;
            std::vector<TH1D*> hh; std::vector<int> hc;
            for (int cc : cs) {
                if (j.h[cc].S <= 0 || j.Neff(cc) < (kind == 0 ? kNeffShape : kNeffMin)) continue;
                hh.push_back(MakeH(j, cc, kind)); hc.push_back(cc);
            }
            double ymin, ymax;
            YLim(hh, kind == 1, ymin, ymax);
            TH1D* fr = NewH("fr", j.nb, j.lo, j.hi);
            if (v == 4) for (int i = 1; i <= 3; ++i) fr->GetXaxis()->SetBinLabel(i, Form("%d", i - 1));
            fr->SetMinimum(ymin); fr->SetMaximum(ymax);
            fr->GetXaxis()->SetTitle(V.title);
            fr->GetYaxis()->SetTitle(kind == 0 ? V.yShape : V.yAbs);
            AxisStyle(fr, 0.055, 1.45);
            if (v == 4) fr->GetXaxis()->SetLabelSize(0.075);
            fr->Draw("AXIS");
            // alternative cut values (thin, grey) and the cut in use
            for (double a : V.alt) {
                TLine* l = new TLine(a, ymin, a, YAt(ymin, ymax, kind == 1, 0.93));
                l->SetLineColor(kGray + 1); l->SetLineStyle(3); l->SetLineWidth(3); l->Draw();
                TLatex* t = new TLatex(a, YAt(ymin, ymax, kind == 1, 0.945), Form("%g", a));
                t->SetTextSize(0.045); t->SetTextAlign(21); t->SetTextColor(kGray + 2); t->Draw();
            }
            for (size_t i = 0; i < hh.size(); ++i) hh[i]->Draw("hist same");
            for (size_t i = 0; i < hh.size(); ++i)
                if (IsSignal(hc[i])) { hh[i]->SetMarkerSize(0); hh[i]->Draw("E0 X0 same"); }
            const double yTop = YAt(ymin, ymax, kind == 1, 0.88);
            TLine* lc = new TLine(V.cut, ymin, V.cut, yTop);
            lc->SetLineColor(kBlack); lc->SetLineStyle(2); lc->SetLineWidth(5); lc->Draw();
            const double dx = 0.09 * (j.hi - j.lo) * V.keep;
            TArrow* ar = new TArrow(V.cut, yTop, V.cut + dx, yTop, 0.022, "|>");
            ar->SetLineWidth(4); ar->SetLineColor(kBlack); ar->SetFillColor(kBlack); ar->Draw();
            fr->Draw("AXIS same");
            Txt(0.17, 0.925, Form("%s, %s", w == 0 ? "bump window" : "peak window",
                                  kind == 0 ? "normalised shapes" : "absolute"), 0.058);
        }
    // legend pad
    c->cd();
    TPad* pl = new TPad(Form("pl_nm1_%s", V.tag), "", xL, 0, 1, 1);
    pl->SetMargin(0, 0, 0, 0);
    pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.04, 0.53, 0.99, 0.99);
    lg->SetTextFont(43); lg->SetTextSize(30); lg->SetFillStyle(0); lg->SetBorderSize(0);
    for (int cc : kLegAll) {
        TH1D* h = NewH("leg", 1, 0, 1);
        LineStyle(h, cc);
        lg->AddEntry(h, kCLeg[cc], "l");
    }
    lg->Draw();
    double y = 0.50;
    const double dy = 0.029;
    TxtPx(0.04, y, Form("N-1 efficiency of %s:", kCutTex[v]), 30); y -= dy * 1.15;
    const int cut = v;
    auto line = [&](int w, const std::vector<int>& cs) {
        std::string s = w == 0 ? "bump: " : "peak: ";
        for (size_t i = 0; i < cs.size(); ++i) {
            const Acc& a = E(w, Cfg(tNm1S, cut), cs[i]);
            s += Fs("%s%s %.2f", i ? ", " : "", kCShort[cs[i]], a.eps());
        }
        TxtPx(0.04, y, s.c_str(), 27); y -= dy;
    };
    line(wBump, {gS2, gS3, aT2});
    line(wPeak, {gS4, gS5, gX17, aT2});
    y -= dy * 0.4;
    TxtPx(0.04, y, v < 4 ? "N-1: 'slides' without this cut" : "N-1: 'slides' (this cut not part of it)", 27); y -= dy;
    TxtPx(0.04, y, "slides = pre + #alpha_{lab} < 50#circ,", 27); y -= dy;
    TxtPx(0.04, y, "   M_{#pi#pi} < 0.46, #alpha_{#eta} > 140#circ, OA > 2#circ", 27); y -= dy;
    TxtPx(0.04, y, "dashed black: cut in use, arrow: kept side", 27); y -= dy;
    if (!V.alt.empty()) { TxtPx(0.04, y, "dotted grey: other values used", 27); y -= dy; }
    TxtPx(0.04, y, "SIMULATION ONLY (SMASH), no CB", 27); y -= dy;
    TxtPx(0.04, y, "std + lep + eta, mode A, w(PDG) #times phys V3", 27); y -= dy;
    TxtPx(0.04, y, "#sigma = #Sigma w / L_{sim}, L_{sim} = 221.05 nb^{-1}", 27); y -= dy;
    TxtPx(0.04, y, "shapes: each normalised to its own #sigma", 27); y -= dy;
    TxtPx(0.04, y, "N_{eff} #geq 30 (shapes), #geq 10 (absolute)", 27); y -= dy;
    TxtPx(0.04, y, "error bars (event groups): signals only", 27); y -= dy;
    TxtPx(0.04, y, "T3, T4 (dashed): not in any sum", 27);
    SaveFig(c, Form("s5_cut_nm1_%s", V.tag));
}

void WriteNm1Table(const std::vector<Job>& jobs) {
    gOut = OpenTab("s5_cut_nm1.txt");
    P("# s5_cut_nm1: N-1 distributions of the cut variables (figures s5_cut_nm1_<var>)\n");
    P("# SIMULATION ONLY. std + lep + eta owned rows (mode A), w (PDG) x phys V3 (108); sigma = sum w / 221.05 nb^-1; no CB\n");
    P("# N-1 selection: 'slides' (pre + OA > 2, alpha_lab < 50, M(pipi) < 0.46, alpha_eta > 140) without the cut of the\n");
    P("#   variable (OA, alpha_lab, M(pipi), alpha_eta), or 'slides' itself (vertex flag, |z_pi - vz|, MM4, pT4)\n");
    P("# columns: sigma in the N-1 selection and window [nb] (event-group error), N_eff, weighted mean, median, 10 %% and 90 %%\n");
    P("#   quantiles (from the histogram), fraction kept by the cut (= N-1 efficiency, event-group error), fraction outside\n");
    P("#   the histogram range; X17 for BR = 1e-5 (W = 2e-5 per row)\n");
    const std::vector<int> cs = {gS2, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aSB, aSP, aT2, aSim, aT234, aBX};
    for (const Job& j : jobs) {
        if (j.var >= 8) continue;
        const int w = j.mlo == kBuLo ? wBump : wPeak;
        P("\n## %s, %s (cut %s), histogram %d bins %g-%g\n", gVar[j.var].tag, kWTxt[w], kCutName[j.var], j.nb, j.lo, j.hi);
        P("%-26s %24s %10s %9s %9s %9s %9s %22s %8s\n", "component", "sigma [nb]", "N_eff", "mean", "median", "q10", "q90",
          "kept by cut", "outside");
        for (int c : cs) {
            const HS& h = j.h[c];
            const Acc& a = E(w, Cfg(tNm1S, j.var), c);
            if (h.S <= 0) { P("%-26s %24s\n", kCShort[c], "-"); continue; }
            P("%-26s %24s %10.1f %9.4g %9.4g %9.4g %9.4g %22s %8.4f\n", kCName[c] + 0, VEg(h.S / kLsim, std::sqrt(h.V) / kLsim).c_str(),
              j.Neff(c), h.Sx / h.S, j.Quantile(c, 0.5), j.Quantile(c, 0.1), j.Quantile(c, 0.9),
              VE(a.eps(), a.err()).c_str(), (h.s[j.nb] + h.s[j.nb + 1]) / h.S);
        }
    }
    CloseTab();
}

// ============================================================================
// part 2: efficiency tables and figures
// ============================================================================
const std::vector<int> kTabC = {gS2, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff,
                                aSB, aSP, aSig, aT2, aSim, aT234, aBX};
void WriteEffTable() {
    gOut = OpenTab("s5_cut_eff.txt");
    P("# s5_cut_eff: efficiency of every cut per signal class and background group, in the bump and the peak window\n");
    P("# SIMULATION ONLY. std + lep + eta owned rows (mode A), w (PDG) x phys V3 (108); sigma = sum w / 221.05 nb^-1; no CB\n");
    P("# base 'pre' (isBest >= 0, vz > -500, PT3 bit, vz > -160, cor_ok); window on the reconstructed M4\n");
    P("# cuts: 1 OA > 2 deg, 2 alpha_lab < 50 deg, 3 M(pipi) < 0.46 GeV/c2, 4 alpha_eta > 140 deg (1-4 = slides),\n");
    P("#   5 both pions in the vertex fit, 6 |z_pi - vz| < 40 mm (both pions), 7 MM4 < 2.85 GeV/c2, 8 pT4 > 400 MeV/c\n");
    P("# errors: event groups, pass and fail rows of one group correlated; X17: BR = 1e-5 (W = 2e-5 per row, not in sums)\n");
    P("# groups: G_eta = T2 with P2 or L2; G_pi0D = other T2 with L1; G_conv = L3; G_vec = L4-L6; G_oth = other T2\n");
    P("#   (here: only eta-set direct-decay e+e- pairs L10 with other pions, see s5_cut_checks.txt); G_fake = T3; G_diff = T4\n");
    for (int w = 0; w < 3; ++w) {
        P("\n\n############################## %s ##############################\n", kWTxt[w]);
        // cut flow
        P("\n## sigma [nb] after each step of the sequential flow (cumulative), and N_eff\n");
        P("%-26s %20s", "component", "pre");
        for (int c = 0; c < kNCut; ++c) P(" %20s", Fs("+%s", kCutName[c]).c_str());
        P(" %9s %9s %9s\n", "Neff pre", "Neff sl.", "Neff all");
        for (int cc : kTabC) {
            const Acc& a0 = E(w, Cfg(tAlone, 0), cc);
            P("%-26s %20s", kCShort[cc], VEg(a0.xsA(), a0.xsAe()).c_str());
            for (int c = 0; c < kNCut; ++c) { const Acc& a = E(w, Cfg(tCum, c), cc); P(" %20s", VEg(a.xsP(), a.xsPe()).c_str()); }
            P(" %9.1f %9.1f %9.1f\n", a0.NeffA(), E(w, kCfgSlides, cc).NeffP(), E(w, Cfg(tCum, kNCut - 1), cc).NeffP());
        }
        for (int t : {(int)tAlone, (int)tStep, (int)tCum, (int)tNm1S, (int)tNm1A}) {
            P("\n## efficiency: %s\n", kTyName[t]);
            P("%-26s", "component");
            for (int c = 0; c < kNCut; ++c) P(" %17s", kCutName[c]);
            P(" %10s\n", "Neff base");
            for (int cc : kTabC) {
                P("%-26s", kCShort[cc]);
                for (int c = 0; c < kNCut; ++c) {
                    const Acc& a = E(w, Cfg(t, c), cc);
                    P(" %17s", a.Sa > 0 ? VE(a.eps(), a.err()).c_str() : "-");
                }
                // N_eff of the base of the last cut of this type (smallest base)
                P(" %10.1f\n", E(w, Cfg(t, t == tAlone || t == tCum ? 0 : kNCut - 1), cc).NeffA());
            }
        }
        P("\n## slides from pre (all four slides cuts together)\n");
        for (int cc : kTabC) { const Acc& a = E(w, kCfgSlides, cc); P("%-26s %17s  N_eff(pass) %.1f\n", kCShort[cc], VE(a.eps(), a.err()).c_str(), a.NeffP()); }
    }
    // where in M4 each cut acts: efficiency per M4 region
    P("\n\n############################## efficiency per M4 region ##############################\n");
    P("# regions of M4 [GeV/c2]: < 0.4576, bump 0.4576-0.5255, gap 0.5255-0.5333, peak 0.5333-0.5587, 0.5587-0.70, 0.70-1.00, > 1.00\n");
    for (int cc : std::vector<int>{aT2, gPi0D, gConv, gEta, gFake, gDiff, aSim}) {
        for (int t : {(int)tAlone, (int)tNm1S}) {
            P("\n## %s: %s\n", kCName[cc], kTyName[t]);
            P("%-18s", "cut");
            for (int k = 0; k < 7; ++k) P(" %17s", kRShort[k]);
            P("\n");
            for (int c = 0; c <= kNCut; ++c) {
                if (c == kNCut && t != tAlone) continue;
                P("%-18s", c == kNCut ? "slides (4 cuts)" : kCutName[c]);
                for (int k = 0; k < 7; ++k) {
                    const Acc& a = c == kNCut ? E(kRegW[k], kCfgSlides, cc) : E(kRegW[k], Cfg(t, c), cc);
                    P(" %17s", a.Sa > 0 ? VE(a.eps(), a.err()).c_str() : "-");
                }
                P("\n");
            }
            if (t == tAlone) {
                P("%-18s", "sigma pre [nb]");
                for (int k = 0; k < 7; ++k) P(" %17s", VEg(E(kRegW[k], Cfg(tAlone, 0), cc).xsA(), E(kRegW[k], Cfg(tAlone, 0), cc).xsAe()).c_str());
                P("\n");
            }
        }
    }
    CloseTab();
}

void DrawEffWindow(int w) {
    TCanvas* c = new TCanvas(Form("c_eff_%d", w), "", 2400, 1800);
    const double xL = 0.76;
    const int types[3] = {tAlone, tNm1S, tStep};
    const char* tt[3] = {"(a) each cut alone, from pre", "(b) N-1: cut given the other 'slides' cuts (last four: on top of slides)",
                         "(c) sequential step: cut given all cuts to its left"};
    const std::vector<int> cs = w == wBump ? std::vector<int>{gS2, gS3, gS5, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2}
                                           : std::vector<int>{gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2};
    for (int k = 0; k < 3; ++k) {
        c->cd();
        TPad* p = new TPad(Form("p_eff_%d_%d", w, k), "", 0, 1. - (k + 1) / 3., xL, 1. - k / 3.);
        p->SetMargin(0.085, 0.015, 0.15, 0.11);
        p->SetGridy(1);
        p->Draw(); p->cd();
        TH1D* fr = NewH("fre", kNCut, 0, kNCut);
        for (int i = 0; i < kNCut; ++i) fr->GetXaxis()->SetBinLabel(i + 1, kCutShort[i]);
        fr->SetMinimum(0); fr->SetMaximum(1.08);
        fr->GetYaxis()->SetTitle("efficiency");
        AxisStyle(fr, 0.065, 0.55);
        fr->GetXaxis()->SetLabelSize(0.085);
        fr->Draw("AXIS");
        const int n = (int)cs.size();
        for (int i = 0; i < n; ++i) {
            const int cc = cs[i];
            TGraphErrors* g = new TGraphErrors();
            for (int cu = 0; cu < kNCut; ++cu) {
                const Acc& a = E(w, Cfg(types[k], cu), cc);
                if (a.Sa <= 0 || a.NeffA() < kNeffMin) continue;
                const int ip = g->GetN();
                g->SetPoint(ip, cu + 0.5 + (i - (n - 1) / 2.) * 0.066, a.eps());
                g->SetPointError(ip, 0, a.err());
            }
            g->SetMarkerStyle(MSty(cc)); g->SetMarkerColor(Col(cc)); g->SetLineColor(Col(cc));
            g->SetMarkerSize(2.4); g->SetLineWidth(3);
            g->Draw("P same");
        }
        for (int cu = 1; cu < kNCut; ++cu) {
            TLine* l = new TLine(cu, 0, cu, 1.08);
            l->SetLineColor(kGray); l->SetLineStyle(1); l->Draw();
        }
        fr->Draw("AXIS same");
        Txt(0.085, 0.915, tt[k], 0.07);
    }
    c->cd();
    TPad* pl = new TPad(Form("pl_eff_%d", w), "", xL, 0, 1, 1);
    pl->SetMargin(0, 0, 0, 0);
    pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.04, 0.52, 0.99, 0.99);
    lg->SetTextFont(43); lg->SetTextSize(32); lg->SetFillStyle(0); lg->SetBorderSize(0);
    for (int cc : cs) {
        TGraphErrors* g = new TGraphErrors();
        g->SetMarkerStyle(MSty(cc)); g->SetMarkerColor(Col(cc)); g->SetLineColor(Col(cc)); g->SetMarkerSize(2.4);
        lg->AddEntry(g, kCLeg[cc], "p");
    }
    lg->Draw();
    double y = 0.49;
    const double dy = 0.021;
    TxtPx(0.04, y, w == wBump ? "BUMP window" : "PEAK window", 34); y -= dy * 1.2;
    TxtPx(0.04, y, w == wBump ? "M_{4} = 0.4576-0.5255 GeV/c^{2}" : "M_{4} = 0.5333-0.5587 GeV/c^{2}", 30); y -= dy * 1.3;
    TxtPx(0.04, y, "SIMULATION ONLY (SMASH), no CB", 30); y -= dy;
    TxtPx(0.04, y, "std + lep + eta, mode A,", 30); y -= dy;
    TxtPx(0.04, y, "w(PDG) #times phys V3", 30); y -= dy;
    TxtPx(0.04, y, "errors: event groups", 30); y -= dy;
    TxtPx(0.04, y, "points shown if N_{eff}(base) #geq 10", 30); y -= dy;
    TxtPx(0.04, y, "X17: BR = 10^{-5}, never summed", 30); y -= dy;
    TxtPx(0.04, y, "T3, T4: not in any sum", 30); y -= dy * 1.3;
    TxtPx(0.04, y, "filled: signals; open: backgrounds", 30); y -= dy * 1.3;
    TxtPx(0.04, y, "cuts: OA > 2#circ", 30); y -= dy;
    TxtPx(0.04, y, "#alpha_{lab} < 50#circ", 30); y -= dy;
    TxtPx(0.04, y, "M_{#pi#pi} < 0.46 GeV/c^{2}", 30); y -= dy;
    TxtPx(0.04, y, "#alpha_{#eta} > 140#circ", 30); y -= dy;
    TxtPx(0.04, y, "#pi vtx: both #pi in the vertex fit", 30); y -= dy;
    TxtPx(0.04, y, "|#Deltaz_{#pi}| = |z_{#pi} - v_{z}| < 40 mm", 30); y -= dy;
    TxtPx(0.04, y, "MM_{4} < 2.85 GeV/c^{2}", 30); y -= dy;
    TxtPx(0.04, y, "p_{T,4} > 400 MeV/c", 30);
    SaveFig(c, Form("s5_cut_eff_%s", kWTag[w]));
}

void DrawFlow() {
    TCanvas* c = new TCanvas("c_flow", "", 2400, 1600);
    const double xL = 0.76;
    for (int w = 0; w < 2; ++w) {
        c->cd();
        TPad* p = new TPad(Form("p_flow_%d", w), "", 0, w == 0 ? 0.5 : 0., xL, w == 0 ? 1. : 0.5);
        p->SetMargin(0.10, 0.015, 0.14, 0.10);
        p->SetLogy(1); p->SetGridy(1);
        p->Draw(); p->cd();
        const std::vector<int> cs = w == wBump ? std::vector<int>{gS2, gS3, gS5, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2}
                                               : std::vector<int>{gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2};
        double ymax = 0, ymin = 1e300;
        std::vector<TGraphErrors*> gg;
        for (int cc : cs) {
            TGraphErrors* g = new TGraphErrors();
            for (int s = 0; s <= kNCut; ++s) {
                const Acc& a = s == 0 ? E(w, Cfg(tAlone, 0), cc) : E(w, Cfg(tCum, s - 1), cc);
                const double v = s == 0 ? a.xsA() : a.xsP(), e = s == 0 ? a.xsAe() : a.xsPe();
                if (!(v > 0)) continue;
                const int ip = g->GetN();
                g->SetPoint(ip, s + 0.5, v); g->SetPointError(ip, 0, e);
                ymax = std::max(ymax, v); ymin = std::min(ymin, v);
            }
            g->SetMarkerStyle(MSty(cc)); g->SetMarkerColor(Col(cc)); g->SetLineColor(Col(cc));
            g->SetLineStyle(LSty(cc)); g->SetLineWidth(cc == aT2 ? 5 : 3); g->SetMarkerSize(2.4);
            gg.push_back(g);
        }
        TH1D* fr = NewH("frf", kNCut + 1, 0, kNCut + 1);
        fr->GetXaxis()->SetBinLabel(1, "pre");
        for (int i = 0; i < kNCut; ++i) fr->GetXaxis()->SetBinLabel(i + 2, Form("+%s", kCutShort[i]));
        fr->SetMinimum(ymin * 0.3); fr->SetMaximum(ymax * 4);
        fr->GetYaxis()->SetTitle(Form("#sigma in the %s window [nb]", w == 0 ? "bump" : "peak"));
        AxisStyle(fr, 0.06, 0.75);
        fr->GetXaxis()->SetLabelSize(0.075);
        fr->Draw("AXIS");
        for (TGraphErrors* g : gg) g->Draw("LP same");
        fr->Draw("AXIS same");
        Txt(0.10, 0.92, Form("sequential cut flow, %s window M_{4} = %s GeV/c^{2}", w == 0 ? "bump" : "peak",
                             w == 0 ? "0.4576-0.5255" : "0.5333-0.5587"), 0.065);
    }
    c->cd();
    TPad* pl = new TPad("pl_flow", "", xL, 0, 1, 1);
    pl->SetMargin(0, 0, 0, 0); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.04, 0.40, 0.99, 0.99);
    lg->SetTextFont(43); lg->SetTextSize(31); lg->SetFillStyle(0); lg->SetBorderSize(0);
    for (int cc : std::vector<int>{gS2, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2}) {
        TGraphErrors* g = new TGraphErrors();
        g->SetMarkerStyle(MSty(cc)); g->SetMarkerColor(Col(cc)); g->SetLineColor(Col(cc)); g->SetLineStyle(LSty(cc));
        g->SetLineWidth(3); g->SetMarkerSize(2.4);
        lg->AddEntry(g, kCLeg[cc], "lp");
    }
    lg->Draw();
    double y = 0.36;
    const double dy = 0.026;
    TxtPx(0.04, y, "order: pre #rightarrow OA > 2#circ", 28); y -= dy;
    TxtPx(0.04, y, "#rightarrow #alpha_{lab} < 50#circ", 28); y -= dy;
    TxtPx(0.04, y, "#rightarrow M_{#pi#pi} < 0.46 #rightarrow #alpha_{#eta} > 140#circ", 28); y -= dy;
    TxtPx(0.04, y, "#rightarrow both #pi in vertex fit", 28); y -= dy;
    TxtPx(0.04, y, "#rightarrow |z_{#pi} - v_{z}| < 40 mm", 28); y -= dy;
    TxtPx(0.04, y, "#rightarrow MM_{4} < 2.85 GeV/c^{2}", 28); y -= dy;
    TxtPx(0.04, y, "#rightarrow p_{T,4} > 400 MeV/c", 28); y -= dy * 1.3;
    TxtPx(0.04, y, "SIMULATION ONLY (SMASH), no CB", 28); y -= dy;
    TxtPx(0.04, y, "std + lep + eta, mode A, w(PDG) #times phys V3", 28); y -= dy;
    TxtPx(0.04, y, "#sigma = #Sigma w / 221.05 nb^{-1}; errors: event groups", 28); y -= dy;
    TxtPx(0.04, y, "X17: BR = 10^{-5}, never summed", 28); y -= dy;
    TxtPx(0.04, y, "T3, T4 (dashed): not in any sum", 28);
    SaveFig(c, "s5_cut_flow");
}

// where in M4 each cut acts: efficiency of T2 per M4 region, each cut alone and N-1 (slides)
const char* const kCutHex[kNCut + 1] = {"#1565C0", "#D84315", "#2E7D32", "#6A1B9A", "#00838F", "#AD1457", "#795548", "#F9A825", "#000000"};
void DrawRegions() {
    TCanvas* c = new TCanvas("c_reg", "", 2400, 1500);
    const double xL = 0.76;
    const int ms[kNCut + 1] = {20, 21, 22, 23, 33, 34, 29, 47, 43};
    for (int k = 0; k < 2; ++k) {
        const int t = k == 0 ? tAlone : tNm1S;
        c->cd();
        TPad* p = new TPad(Form("p_reg_%d", k), "", 0, k == 0 ? 0.5 : 0., xL, k == 0 ? 1. : 0.5);
        p->SetMargin(0.09, 0.015, 0.15, 0.11);
        p->SetLogy(1); p->SetGridy(1);
        p->Draw(); p->cd();
        TH1D* fr = NewH("frg", 7, 0, 7);
        for (int i = 0; i < 7; ++i) fr->GetXaxis()->SetBinLabel(i + 1, kRShort[i]);
        fr->SetMinimum(k == 0 ? 5e-5 : 1e-3); fr->SetMaximum(1.6);
        fr->GetYaxis()->SetTitle("efficiency of T2");
        AxisStyle(fr, 0.065, 0.62);
        fr->GetXaxis()->SetLabelSize(0.075);
        fr->Draw("AXIS");
        for (int w : {wBump, wPeak}) {
            const int ib = w == wBump ? 1 : 3;
            TBox* b = new TBox(ib, fr->GetMinimum(), ib + 1, 1.6);
            b->SetFillColor(TColor::GetColor(w == wBump ? "#FFE0B2" : "#FFCDD2")); b->SetLineWidth(0); b->Draw();
        }
        for (int cu = 0; cu <= kNCut; ++cu) {
            if (cu == kNCut && k == 1) continue;
            TGraphErrors* g = new TGraphErrors();
            for (int i = 0; i < 7; ++i) {
                const Acc& a = cu == kNCut ? E(kRegW[i], kCfgSlides, aT2) : E(kRegW[i], Cfg(t, cu), aT2);
                if (a.Sa <= 0 || !(a.eps() > 0)) continue;
                const int ip = g->GetN();
                g->SetPoint(ip, i + 0.5 + (cu - 4) * 0.045, a.eps()); g->SetPointError(ip, 0, a.err());
            }
            const int col = TColor::GetColor(kCutHex[cu]);
            g->SetMarkerStyle(ms[cu]); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(2.4);
            g->SetLineWidth(cu == kNCut ? 5 : 3); g->SetLineStyle(cu == kNCut ? 2 : 1);
            g->Draw("LP same");
        }
        fr->Draw("AXIS same");
        Txt(0.09, 0.915, k == 0 ? "(a) T2 (physical background): each cut alone, from pre, per M_{4} region"
                                : "(b) T2: N-1 (slides; the last four on top of slides), per M_{4} region", 0.068);
    }
    c->cd();
    TPad* pl = new TPad("pl_reg", "", xL, 0, 1, 1);
    pl->SetMargin(0, 0, 0, 0); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.04, 0.50, 0.99, 0.99);
    lg->SetTextFont(43); lg->SetTextSize(32); lg->SetFillStyle(0); lg->SetBorderSize(0);
    for (int cu = 0; cu <= kNCut; ++cu) {
        TGraphErrors* g = new TGraphErrors();
        const int col = TColor::GetColor(kCutHex[cu]);
        g->SetMarkerStyle(ms[cu]); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(2.4);
        g->SetLineWidth(cu == kNCut ? 5 : 3); g->SetLineStyle(cu == kNCut ? 2 : 1);
        lg->AddEntry(g, cu == kNCut ? "slides (4 cuts), panel (a)" : kCutTex[cu], "lp");
    }
    lg->Draw();
    double y = 0.46;
    const double dy = 0.028;
    TxtPx(0.04, y, "regions of M_{4} [GeV/c^{2}]:", 30); y -= dy;
    TxtPx(0.04, y, "bump 0.4576-0.5255 (orange band)", 30); y -= dy;
    TxtPx(0.04, y, "gap 0.5255-0.5333", 30); y -= dy;
    TxtPx(0.04, y, "peak 0.5333-0.5587 (red band)", 30); y -= dy * 1.3;
    TxtPx(0.04, y, "SIMULATION ONLY (SMASH), no CB", 30); y -= dy;
    TxtPx(0.04, y, "std + lep + eta, mode A,", 30); y -= dy;
    TxtPx(0.04, y, "w(PDG) #times phys V3", 30); y -= dy;
    TxtPx(0.04, y, "errors: event groups", 30);
    SaveFig(c, "s5_cut_regions");
}

// ============================================================================
// part 4: composition per M4 region
// ============================================================================
const int kNI = kNC + 72;   // components + for T2, T3, T4: lepton-pair class (12 slots) + pion-pair class (12 slots)
inline int ItemL(int t, int l) { return kNC + (t - 2) * 24 + l; }
inline int ItemP(int t, int p) { return kNC + (t - 2) * 24 + 12 + p; }
const int kNDen = 4;        // denominators of the shares: T1+T2, T2, T3, T4
const int kDenIdx[kNDen] = {aSim, aT2, gFake, gDiff};
struct RA { double S = 0, aa = 0, N = 0; double ad[kNDen] = {0, 0, 0, 0}; };
std::vector<RA> gC;             // [sel][region][item]
std::vector<double> gCellS, gCellV, gCellN;   // [sel][region][lcls 12][pcls 10], T2 only
inline RA& C(int s, int rg, int it) { return gC[((size_t)s * (kNR + 1) + rg) * kNI + it]; }
inline size_t CellIdx(int s, int rg, int l, int p) { return (((size_t)s * (kNR + 1) + rg) * 12 + l) * 10 + p; }
// share of item it in denominator k with the delta-method event-group error
std::pair<double, double> Share(int s, int rg, int it, int k) {
    const RA& a = C(s, rg, it);
    const RA& d = C(s, rg, kDenIdx[k]);
    if (d.S <= 0) return std::make_pair(0., 0.);
    const double R = a.S / d.S;
    const double v = a.aa - 2 * R * a.ad[k] + R * R * d.aa;
    return std::make_pair(R, v > 0 ? std::sqrt(v) / d.S : 0.);
}
void RunComp() {
    TStopwatch sw;
    gC.assign((size_t)2 * (kNR + 1) * kNI, RA());
    gCellS.assign((size_t)2 * (kNR + 1) * 120, 0.); gCellV = gCellS; gCellN = gCellS;
    std::vector<double> ga((size_t)2 * (kNR + 1) * kNI, 0.), gc((size_t)2 * (kNR + 1) * 120, 0.);
    std::vector<int> touchedSR;
    std::vector<char> tsr(2 * (kNR + 1), 0);
    int cur = INT_MIN;
    for (size_t i = 0; i <= gR.size(); ++i) {
        const bool end = i == gR.size();
        if (end || gR[i].gid != cur) {
            for (int sr : touchedSR) {
                double* a = &ga[(size_t)sr * kNI];
                double den[kNDen];
                for (int k = 0; k < kNDen; ++k) den[k] = a[kDenIdx[k]];
                for (int it = 0; it < kNI; ++it) {
                    if (a[it] == 0) continue;
                    RA& r = gC[(size_t)sr * kNI + it];
                    r.aa += a[it] * a[it];
                    for (int k = 0; k < kNDen; ++k) r.ad[k] += a[it] * den[k];
                    a[it] = 0;
                }
                double* cc = &gc[(size_t)sr * 120];
                for (int q = 0; q < 120; ++q) if (cc[q] != 0) { gCellV[(size_t)sr * 120 + q] += cc[q] * cc[q]; cc[q] = 0; }
                tsr[sr] = 0;
            }
            touchedSR.clear();
            if (end) break;
            cur = gR[i].gid;
        }
        const MRow& r = gR[i];
        int rg = -1;
        for (int k = 0; k < kNR; ++k) if (r.m4 >= kRLo[k] && r.m4 < kRHi[k]) rg = k;
        for (int s = 0; s < 2; ++s) {
            if (s == 1 && (r.cb & kSlM) != kSlM) continue;
            for (int rr : {rg, kNR}) {
                if (rr < 0) continue;
                const int sr = s * (kNR + 1) + rr;
                double* a = &ga[(size_t)sr * kNI];
                auto add = [&](int it) { RA& x = gC[(size_t)sr * kNI + it]; x.S += r.w; x.N += 1; a[it] += r.w; };
                unsigned m = r.mask;
                while (m) { const int c = __builtin_ctz(m); m &= m - 1; add(c); }
                if (r.tmpl >= 2 && r.tmpl <= 4 && r.lcls >= 0 && r.lcls < 12 && r.pcls >= 0 && r.pcls < 10) {
                    add(ItemL(r.tmpl, r.lcls));
                    add(ItemP(r.tmpl, r.pcls));
                    if (r.tmpl == 2) {
                        const size_t q = (size_t)sr * 120 + r.lcls * 10 + r.pcls;
                        gCellS[q] += r.w; gCellN[q] += 1; gc[q] += r.w;
                    }
                }
                if (!tsr[sr]) { tsr[sr] = 1; touchedSR.push_back(sr); }
            }
        }
    }
    printf("s5c::RunComp: %.1f s\n", sw.RealTime());
}

void WriteCompTable() {
    gOut = OpenTab("s5_cut_comp.txt");
    P("# s5_cut_comp: composition of the simulated pi+pi-e+e- rows per M4 region after 'pre' and after 'slides'\n");
    P("# SIMULATION ONLY. std + lep + eta owned rows (mode A), w (PDG) x phys V3 (108); sigma = sum w / 221.05 nb^-1; no CB\n");
    P("# share of T1+T2 (the simulated sum) and share of T2 (physical background) with delta-method event-group errors;\n");
    P("#   T3 and T4 are not in T1+T2: their 'share of T1+T2' is a ratio to it; X17 (BR 1e-5) only as a ratio\n");
    P("# regions of M4 [GeV/c2]: < 0.4576, bump 0.4576-0.5255, gap 0.5255-0.5333, peak 0.5333-0.5587, 0.5587-0.70, 0.70-1.00, > 1.00\n");
    const std::vector<int> cs = {gS2, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aSig, aT2, aSim, aT234};
    for (int s = 0; s < 2; ++s) {
        P("\n\n############################## selection %s ##############################\n", s == 0 ? "pre" : "slides");
        for (int rg = 0; rg <= kNR; ++rg) {
            P("\n## %s, M4 %s\n", s == 0 ? "pre" : "slides", kRName[rg]);
            P("%-26s %26s %24s %24s %10s %10s\n", "component", "sigma [nb]", "share of T1+T2", "share of T2", "N_eff", "rows");
            for (int c : cs) {
                const RA& a = C(s, rg, c);
                const double neff = a.aa > 0 ? a.S * a.S / a.aa : 0.;
                const std::pair<double, double> sh = Share(s, rg, c, 0), s2 = Share(s, rg, c, 1);
                const bool inT2 = c >= gEta && c <= gOth;
                P("%-26s %26s %24s %24s %10.1f %10.0f\n", kCName[c], VEg(a.S / kLsim, std::sqrt(a.aa) / kLsim).c_str(),
                  VE(sh.first, sh.second, 4).c_str(), (inT2 || c == aT2) ? VE(s2.first, s2.second, 4).c_str() : "", neff, a.N);
            }
            for (int t = 2; t <= 4; ++t) {
                const int k = t == 2 ? 1 : t == 3 ? 2 : 3;
                P("# %s: lepton-pair class (sigma [nb], share of %s)\n", t == 2 ? "T2" : t == 3 ? "T3" : "T4", t == 2 ? "T2" : t == 3 ? "T3" : "T4");
                for (int l = 0; l < 12; ++l) {
                    const RA& a = C(s, rg, ItemL(t, l));
                    if (a.S <= 0) continue;
                    const std::pair<double, double> sh = Share(s, rg, ItemL(t, l), k);
                    P("    %-26s %26s %24s  N_eff %.1f\n", LclsName(l), VEg(a.S / kLsim, std::sqrt(a.aa) / kLsim).c_str(),
                      VE(sh.first, sh.second, 4).c_str(), a.aa > 0 ? a.S * a.S / a.aa : 0.);
                }
                P("# %s: pion-pair class\n", t == 2 ? "T2" : t == 3 ? "T3" : "T4");
                for (int p = 0; p < 10; ++p) {
                    const RA& a = C(s, rg, ItemP(t, p));
                    if (a.S <= 0) continue;
                    const std::pair<double, double> sh = Share(s, rg, ItemP(t, p), k);
                    P("    %-26s %26s %24s  N_eff %.1f\n", PclsName(p), VEg(a.S / kLsim, std::sqrt(a.aa) / kLsim).c_str(),
                      VE(sh.first, sh.second, 4).c_str(), a.aa > 0 ? a.S * a.S / a.aa : 0.);
                }
            }
            // top cells of T2
            std::vector<std::pair<double, int>> cells;
            const double t2 = C(s, rg, aT2).S;
            for (int q = 0; q < 120; ++q) {
                const double v = gCellS[CellIdx(s, rg, 0, 0) + q];
                if (v > 0) cells.push_back(std::make_pair(-v, q));
            }
            std::sort(cells.begin(), cells.end());
            P("# T2: largest (lepton-pair, pion-pair) cells: sigma [nb], share of T2, N_eff, rows\n");
            for (size_t i = 0; i < cells.size() && i < 6; ++i) {
                const int q = cells[i].second;
                const double S = gCellS[CellIdx(s, rg, 0, 0) + q], V = gCellV[CellIdx(s, rg, 0, 0) + q], N = gCellN[CellIdx(s, rg, 0, 0) + q];
                P("    %-24s x %-22s %22s  %.4f  N_eff %8.1f  rows %8.0f\n", LclsName(q / 10), PclsName(q % 10),
                  VEg(S / kLsim, std::sqrt(V) / kLsim).c_str(), t2 > 0 ? S / t2 : 0., V > 0 ? S * S / V : 0., N);
            }
        }
    }
    // compact summary
    P("\n\n############################## summary ##############################\n");
    P("# per selection and region: sigma(T1+T2), sigma(T2) [nb]; shares of T2 per group; dominant lepton-pair and pion-pair class of T2\n");
    for (int s = 0; s < 2; ++s) {
        P("\n## %s\n", s == 0 ? "pre" : "slides");
        P("%-20s %20s %20s %7s %7s %7s %7s %7s %7s %7s  %-28s %-26s\n", "region", "T1+T2 [nb]", "T2 [nb]", "G_eta", "G_pi0D",
          "G_conv", "G_vec", "G_oth", "T3/T12", "T4/T12", "dominant ee (share of T2)", "dominant pions (share of T2)");
        for (int rg = 0; rg <= kNR; ++rg) {
            const RA& sim = C(s, rg, aSim);
            const RA& t2 = C(s, rg, aT2);
            int bl = 0, bp = 0;
            for (int l = 0; l < 12; ++l) if (C(s, rg, ItemL(2, l)).S > C(s, rg, ItemL(2, bl)).S) bl = l;
            for (int p = 0; p < 10; ++p) if (C(s, rg, ItemP(2, p)).S > C(s, rg, ItemP(2, bp)).S) bp = p;
            P("%-20s %20s %20s", kRName[rg], VEg(sim.S / kLsim, std::sqrt(sim.aa) / kLsim).c_str(),
              VEg(t2.S / kLsim, std::sqrt(t2.aa) / kLsim).c_str());
            for (int g = gEta; g <= gOth; ++g) P(" %7.4f", Share(s, rg, g, 1).first);
            P(" %7.4f %7.4f  %-28s %-26s\n", Share(s, rg, gFake, 0).first, Share(s, rg, gDiff, 0).first,
              Fs("%s (%.3f)", LclsName(bl), Share(s, rg, ItemL(2, bl), 1).first).c_str(),
              Fs("%s (%.3f)", PclsName(bp), Share(s, rg, ItemP(2, bp), 1).first).c_str());
        }
    }
    CloseTab();
}

// palettes (cold) for the class bars
const char* const kHexL[12] = {"#BDBDBD", "#1565C0", "#5E35B1", "#00ACC1", "#00796B", "#26A69A", "#80CBC4", "#616161",
                               "#9E9E9E", "#E0E0E0", "#90A4AE", "#37474F"};
const char* const kHexP[10] = {"#E0E0E0", "#BDBDBD", "#5E35B1", "#283593", "#1E88E5", "#4FC3F7", "#00838F", "#78909C",
                               "#B0BEC5", "#004D40"};
void DrawComp() {
    TCanvas* c = new TCanvas("c_comp", "", 2400, 1800);
    const double yB = 0.25;
    const char* colT[3] = {"share of T1 + T2 per group", "T2: lepton-pair class", "T2: pion-pair class"};
    std::vector<int> grp1 = {gS2, gS3, gS4, gS5, gEta, gPi0D, gConv, gVec, gOth};
    for (int s = 0; s < 2; ++s)
        for (int col = 0; col < 3; ++col) {
            c->cd();
            TPad* p = new TPad(Form("p_comp_%d_%d", s, col), "", col / 3., s == 0 ? yB + (1 - yB) / 2 : yB, (col + 1) / 3.,
                               s == 0 ? 1. : yB + (1 - yB) / 2);
            p->SetMargin(0.26, 0.03, 0.13, 0.11);
            p->Draw(); p->cd();
            TH2D* fr = new TH2D(Form("frc_%d_%d", s, col), "", 10, 0, 1, kNR + 1, 0, kNR + 1);
            fr->SetDirectory(nullptr);
            for (int rg = 0; rg <= kNR; ++rg) fr->GetYaxis()->SetBinLabel(kNR + 1 - rg, kRShort[rg]);
            fr->GetXaxis()->SetTitle(col == 0 ? "share of T1 + T2" : "share of T2");
            AxisStyle(fr, 0.062, 1.0, 0.95);
            fr->GetYaxis()->SetLabelSize(0.068);
            fr->Draw("AXIS");
            for (int rg = 0; rg <= kNR; ++rg) {
                const double yc = kNR + 0.5 - rg, h = 0.36;
                double x = 0;
                auto seg = [&](double f, int colr) {
                    if (f <= 0) return;
                    TBox* b = new TBox(x, yc - h, std::min(1., x + f), yc + h);
                    b->SetFillColor(colr); b->SetLineColor(kWhite); b->SetLineWidth(1); b->Draw("l");
                    b = new TBox(x, yc - h, std::min(1., x + f), yc + h);
                    b->SetFillStyle(0); b->SetLineColor(kWhite); b->SetLineWidth(1); b->Draw("l");
                    if (f > 0.13) {
                        TLatex* t = new TLatex(x + f / 2, yc, Form("%.0f", 100 * f));
                        t->SetTextAlign(22); t->SetTextSize(0.05); t->SetTextColor(kWhite); t->Draw();
                    }
                    x += f;
                };
                if (col == 0) for (int g : grp1) seg(Share(s, rg, g, 0).first, Col(g));
                else if (col == 1) for (int l = 0; l < 12; ++l) seg(Share(s, rg, ItemL(2, l), 1).first, TColor::GetColor(kHexL[l]));
                else for (int q = 0; q < 10; ++q) seg(Share(s, rg, ItemP(2, q), 1).first, TColor::GetColor(kHexP[q]));
            }
            fr->Draw("AXIS same");
            Txt(0.26, 0.92, Form("%s: %s", s == 0 ? "pre" : "slides", colT[col]), 0.062);
        }
    // legends
    c->cd();
    TPad* pl = new TPad("pl_comp", "", 0, 0, 1, yB);
    pl->SetMargin(0, 0, 0, 0); pl->Draw(); pl->cd();
    TLegend* l1 = new TLegend(0.01, 0.13, 0.33, 0.99);
    l1->SetTextFont(43); l1->SetTextSize(30); l1->SetNColumns(2); l1->SetFillStyle(0); l1->SetBorderSize(0);
    const char* shortLeg[9] = {"S2 bump #pi^{0} Dalitz", "S3 bump conversion", "S4 peak conversion", "S5 peak direct",
                               "G_{#eta}", "G_{#pi^{0}D} (L1)", "G_{conv} (L3)", "G_{vec} (L4-6)", "G_{oth} (L10)"};
    for (size_t i = 0; i < grp1.size(); ++i) {
        TBox* b = new TBox(0, 0, 1, 1); b->SetFillColor(Col(grp1[i])); b->SetLineColor(Col(grp1[i]));
        l1->AddEntry(b, shortLeg[i], "f");
    }
    l1->Draw();
    TLegend* l2 = new TLegend(0.34, 0.13, 0.66, 0.99);
    l2->SetTextFont(43); l2->SetTextSize(30); l2->SetNColumns(2); l2->SetFillStyle(0); l2->SetBorderSize(0);
    for (int l : {1, 2, 3, 4, 5, 6, 10, 11}) {
        TBox* b = new TBox(0, 0, 1, 1); b->SetFillColor(TColor::GetColor(kHexL[l])); b->SetLineColor(TColor::GetColor(kHexL[l]));
        l2->AddEntry(b, LclsName(l), "f");
    }
    l2->Draw();
    TLegend* l3 = new TLegend(0.67, 0.13, 0.99, 0.99);
    l3->SetTextFont(43); l3->SetTextSize(30); l3->SetNColumns(2); l3->SetFillStyle(0); l3->SetBorderSize(0);
    for (int q = 0; q < 10; ++q) {
        TBox* b = new TBox(0, 0, 1, 1); b->SetFillColor(TColor::GetColor(kHexP[q])); b->SetLineColor(TColor::GetColor(kHexP[q]));
        l3->AddEntry(b, PclsName(q), "f");
    }
    l3->Draw();
    c->cd();
    TxtPx(0.005, 0.006, "SIMULATION ONLY (SMASH std + lep + eta, mode A, w(PDG) #times phys V3), no CB; numbers in %; "
          "regions in M_{4} [GeV/c^{2}]: bump 0.4576-0.5255, gap 0.5255-0.5333, peak 0.5333-0.5587", 30);
    SaveFig(c, "s5_cut_comp");
}

void DrawM4(const Job& jPre, const Job& jSl) {
    TCanvas* c = new TCanvas("c_m4", "", 2400, 1600);
    const double xL = 0.76;
    const std::vector<int> cs = {gS2, gS3, gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aSim};
    for (int s = 0; s < 2; ++s) {
        c->cd();
        TPad* p = new TPad(Form("p_m4_%d", s), "", 0, s == 0 ? 0.5 : 0., xL, s == 0 ? 1. : 0.5);
        p->SetMargin(0.10, 0.015, 0.15, 0.10);
        p->SetLogy(1);
        p->Draw(); p->cd();
        const Job& j = s == 0 ? jPre : jSl;
        std::vector<TH1D*> hh;
        for (int cc : cs) if (j.h[cc].S > 0) hh.push_back(MakeH(j, cc, 1));
        double ymin, ymax;
        YLim(hh, true, ymin, ymax);
        ymin = std::max(ymin, ymax * 1e-6);
        TH1D* fr = NewH("frm", j.nb, j.lo, j.hi);
        fr->SetMinimum(ymin); fr->SetMaximum(ymax);
        fr->GetXaxis()->SetTitle(gVar[8].title); fr->GetYaxis()->SetTitle(gVar[8].yAbs);
        AxisStyle(fr, 0.06, 0.78);
        fr->Draw("AXIS");
        for (int w = 0; w < 2; ++w) {
            TBox* b = new TBox(kWLo[w], ymin, kWHi[w], ymax);
            b->SetFillColor(TColor::GetColor(w == 0 ? "#FFE0B2" : "#FFCDD2")); b->SetLineWidth(0); b->Draw();
        }
        for (TH1D* h : hh) h->Draw("hist same");
        fr->Draw("AXIS same");
        Txt(0.10, 0.92, s == 0 ? "pre: isBest #geq 0, PT3, v_{z} > -160 mm" : "slides: pre + OA > 2 deg, #alpha_{lab} < 50 deg, "
            "M_{#pi#pi} < 0.46, #alpha_{#eta} > 140 deg", 0.065);
    }
    c->cd();
    TPad* pl = new TPad("pl_m4", "", xL, 0, 1, 1);
    pl->SetMargin(0, 0, 0, 0); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.04, 0.40, 0.99, 0.99);
    lg->SetTextFont(43); lg->SetTextSize(31); lg->SetFillStyle(0); lg->SetBorderSize(0);
    for (int cc : cs) { TH1D* h = NewH("legm", 1, 0, 1); LineStyle(h, cc); lg->AddEntry(h, kCLeg[cc], "l"); }
    lg->Draw();
    double y = 0.36;
    const double dy = 0.026;
    TxtPx(0.04, y, "bands: bump (orange) and peak (red) windows", 28); y -= dy;
    TxtPx(0.04, y, "SIMULATION ONLY (SMASH), no CB", 28); y -= dy;
    TxtPx(0.04, y, "std + lep + eta, mode A, w(PDG) #times phys V3", 28); y -= dy;
    TxtPx(0.04, y, "#sigma = #Sigma w / 221.05 nb^{-1}; 5 MeV/c^{2} bins", 28); y -= dy;
    TxtPx(0.04, y, "X17: BR = 10^{-5}, never summed", 28); y -= dy;
    TxtPx(0.04, y, "T3, T4 (dashed): not in T1 + T2", 28);
    SaveFig(c, "s5_cut_m4");
}

// ============================================================================
// part 8: ranking
// ============================================================================
struct Tgt { const char* name; const char* tex; int w, s, b; };
const std::vector<Tgt> kTgt = {
    {"peak: S4+S5 vs T2", "peak: S = S4 + S5, B = T2", wPeak, aSP, aT2},
    {"bump: S2+S3 vs T2", "bump: S = S2 + S3, B = T2", wBump, aSB, aT2},
    {"X17: X17 vs T2+S4+S5", "X17: S = X17, B = T2 + S4 + S5", wPeak, gX17, aBX},
    {"peak, variant B = T2+T3+T4", "", wPeak, aSP, aT234},
    {"bump, variant B = T2+T3+T4", "", wBump, aSB, aT234},
    {"S4 (conversion peak) vs T2", "", wPeak, gS4, aT2},
    {"S5 (direct peak) vs T2", "", wPeak, gS5, aT2},
    {"S2 (bump, Dalitz) vs T2", "", wBump, gS2, aT2},
    {"S3 (bump, conversion) vs T2", "", wBump, gS3, aT2}};
struct QV { double es, ees, eb, eeb, q, eq; };
QV Qof(int w, int k, int s, int b) {
    const Acc& as = E(w, k, s);
    const Acc& ab = E(w, k, b);
    QV o;
    o.es = as.eps(); o.ees = as.err(); o.eb = ab.eps(); o.eeb = ab.err();
    o.q = o.eb > 0 ? o.es / std::sqrt(o.eb) : 0.;
    o.eq = (o.es > 0 && o.eb > 0) ? o.q * std::hypot(o.ees / o.es, 0.5 * o.eeb / o.eb) : 0.;
    return o;
}
const char* Verdict(const QV& v) {
    const double tol = std::max(0.02, 2 * v.eq);
    if (v.q > 1 + tol) return v.es < 0.9 ? "informative, costs signal" : "informative";
    if (v.q < 1 - tol) return "harmful (Q < 1)";
    return v.es < 0.9 ? "costs signal, no gain" : "neutral";
}
void WriteRankTable() {
    gOut = OpenTab("s5_cut_rank.txt");
    P("# s5_cut_rank: information carried by each cut; Q = eps_S / sqrt(eps_B) (factor on S/sqrt(B) when B dominates),\n");
    P("#   purity gain eps_S / eps_B; errors: event groups of S and of B, combined as independent (rows of one event can\n");
    P("#   be in S and in B: the small correlation is neglected)\n");
    P("# SIMULATION ONLY, no CB (the data have a CB that is not in B): Q, purity and significances here are upper bounds\n");
    P("# verdict: informative if Q > 1 + max(0.02, 2 dQ); harmful if Q < 1 - max(0.02, 2 dQ); 'costs signal' if eps_S < 0.9\n");
    P("# variants: alone = from pre; N-1 (slides) = given the other slides cuts (vertex, |z|, MM4, pT4: on top of slides);\n");
    P("#   N-1 (all) = given the other seven cuts\n");
    for (const Tgt& t : kTgt) {
        P("\n## %s (%s)\n", t.name, kWTxt[t.w]);
        P("%-16s | %-60s | %-60s | %-60s\n", "cut", "alone from pre: eps_S, eps_B, Q", "N-1 (slides): eps_S, eps_B, Q",
          "N-1 (all eight): eps_S, eps_B, Q");
        for (int c = 0; c < kNCut; ++c) {
            P("%-16s", kCutName[c]);
            for (int ty : {(int)tAlone, (int)tNm1S, (int)tNm1A}) {
                const QV v = Qof(t.w, Cfg(ty, c), t.s, t.b);
                P(" | %6.3f(%3.0f) %6.3f(%3.0f) Q %6.3f(%3.0f) %-22s", v.es, 1000 * v.ees, v.eb, 1000 * v.eeb, v.q, 1000 * v.eq, Verdict(v));
            }
            P("\n");
        }
        // full selections
        const QV vs = Qof(t.w, kCfgSlides, t.s, t.b), va = Qof(t.w, Cfg(tCum, kNCut - 1), t.s, t.b);
        P("%-16s | eps_S %.3f(%.0f) eps_B %.4f(%.0f) Q %.3f(%.0f) purity gain %.2f\n", "slides from pre", vs.es, 1000 * vs.ees, vs.eb,
          1e4 * vs.eeb, vs.q, 1000 * vs.eq, vs.eb > 0 ? vs.es / vs.eb : 0.);
        P("%-16s | eps_S %.3f(%.0f) eps_B %.4f(%.0f) Q %.3f(%.0f) purity gain %.2f\n", "all 8 from pre", va.es, 1000 * va.ees, va.eb,
          1e4 * va.eeb, va.q, 1000 * va.eq, va.eb > 0 ? va.es / va.eb : 0.);
    }
    // what each cut removes (N-1 slides): removed sigma per background group
    P("\n\n## what each cut removes, N-1 (slides): sigma removed [nb] per group and its share of the removed T2\n");
    for (int w = 0; w < 2; ++w) {
        P("\n# %s\n", kWTxt[w]);
        P("%-16s", "cut");
        for (int g : std::vector<int>{gS2, gS3, gS4, gS5, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2}) P(" %15s", kCShort[g]);
        P("\n");
        for (int c = 0; c < kNCut; ++c) {
            P("%-16s", kCutName[c]);
            const Acc& t2 = E(w, Cfg(tNm1S, c), aT2);
            const double rt2 = t2.Sa - t2.Sp;
            for (int g : std::vector<int>{gS2, gS3, gS4, gS5, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2}) {
                const Acc& a = E(w, Cfg(tNm1S, c), g);
                const double r = (a.Sa - a.Sp) / kLsim;
                if (g >= gEta && g <= gOth && rt2 > 0) P(" %8.4f(%4.2f)", r, (a.Sa - a.Sp) / rt2);
                else P(" %15.4f", r);
            }
            P("\n");
        }
    }
    // expected counts at L_data (UNCONFIRMED)
    P("\n\n## expected rows at L_data = %.1f nb^-1 (UNCONFIRMED; scale only), no CB, no efficiency correction\n", kLd);
    P("%-20s %-8s %14s %14s %14s %9s %12s %14s %12s\n", "selection", "window", "S", "B (T2)", "B' (T2+T3+T4)", "S/B",
      "S/sqrt(S+B)", "S/sqrt(S+B')", "N_eff(S)");
    for (int w = 0; w < 2; ++w) {
        const int s = w == wPeak ? aSP : aSB;
        for (int sel = 0; sel < 3; ++sel) {
            const int k = sel == 0 ? Cfg(tAlone, 0) : sel == 1 ? kCfgSlides : Cfg(tCum, kNCut - 1);
            const Acc& as = E(w, k, s);
            const Acc& ab = E(w, k, aT2);
            const Acc& ab2 = E(w, k, aT234);
            const double S = (sel == 0 ? as.Sa : as.Sp) * kLd / kLsim, B = (sel == 0 ? ab.Sa : ab.Sp) * kLd / kLsim,
                         B2 = (sel == 0 ? ab2.Sa : ab2.Sp) * kLd / kLsim;
            P("%-20s %-8s %14.0f %14.0f %14.0f %9.4f %12.1f %14.1f %12.1f\n", sel == 0 ? "pre" : sel == 1 ? "slides" : "all 8 cuts",
              kWTag[w], S, B, B2, B > 0 ? S / B : 0., S / std::sqrt(S + B), S / std::sqrt(S + B2), sel == 0 ? as.NeffA() : as.NeffP());
        }
    }
    P("\n# X17 (BR = 1e-5) in the peak window against B_X = T2 + S4 + S5 (the eta peak is background for X17)\n");
    for (int sel = 0; sel < 3; ++sel) {
        const int k = sel == 0 ? Cfg(tAlone, 0) : sel == 1 ? kCfgSlides : Cfg(tCum, kNCut - 1);
        const Acc& as = E(wPeak, k, gX17);
        const Acc& ab = E(wPeak, k, aBX);
        const double S = (sel == 0 ? as.Sa : as.Sp) * kLd / kLsim, B = (sel == 0 ? ab.Sa : ab.Sp) * kLd / kLsim;
        P("%-20s X17 %10.1f  B_X %12.0f  S/B %.5f  S/sqrt(B) %.3f\n", sel == 0 ? "pre" : sel == 1 ? "slides" : "all 8 cuts", S, B,
          S / B, S / std::sqrt(B));
    }
    CloseTab();
}

void DrawRank() {
    TCanvas* c = new TCanvas("c_rank", "", 2400, 1900);
    const double yT = 0.93;
    const int ty[3] = {tAlone, tNm1S, tNm1A};
    const int ms[3] = {24, 21, 22};
    const char* tyLeg[3] = {"alone, from pre", "N-1 (slides; last four on top of slides)", "N-1 (all eight cuts)"};
    const int tcol[3] = {TColor::GetColor("#616161"), TColor::GetColor("#1565C0"), TColor::GetColor("#D62728")};
    for (int r = 0; r < 3; ++r)
        for (int col = 0; col < 2; ++col) {
            const Tgt& t = kTgt[r];
            c->cd();
            TPad* p = new TPad(Form("p_rank_%d_%d", r, col), "", col / 2., yT - (r + 1) * yT / 3, (col + 1) / 2., yT - r * yT / 3);
            p->SetMargin(0.12, 0.02, 0.17, 0.12);
            p->SetGridy(1);
            p->Draw(); p->cd();
            TH1D* fr = NewH("frr", kNCut, 0, kNCut);
            for (int i = 0; i < kNCut; ++i) fr->GetXaxis()->SetBinLabel(i + 1, kCutShort[i]);
            double lo = col == 0 ? 0.85 : 0., hi = col == 0 ? 1.15 : 1.08;
            std::vector<TGraphErrors*> gg;
            for (int k = 0; k < 3; ++k) {
                TGraphErrors* g = new TGraphErrors();
                for (int cu = 0; cu < kNCut; ++cu) {
                    const QV v = Qof(t.w, Cfg(ty[k], cu), t.s, t.b);
                    const double y = col == 0 ? v.q : v.es, e = col == 0 ? v.eq : v.ees;
                    if (!(y > 0)) continue;
                    if (col == 0) { lo = std::min(lo, y - e - 0.03); hi = std::max(hi, y + e + 0.06); }
                    const int ip = g->GetN();
                    g->SetPoint(ip, cu + 0.5 + (k - 1) * 0.22, y); g->SetPointError(ip, 0, e);
                }
                g->SetMarkerStyle(ms[k]); g->SetMarkerColor(tcol[k]); g->SetLineColor(tcol[k]);
                g->SetMarkerSize(2.6); g->SetLineWidth(3);
                gg.push_back(g);
            }
            fr->SetMinimum(lo); fr->SetMaximum(hi);
            fr->GetYaxis()->SetTitle(col == 0 ? "Q = #varepsilon_{S}/#sqrt{#varepsilon_{B}}" : "#varepsilon_{S}");
            AxisStyle(fr, 0.07, 0.75);
            fr->GetXaxis()->SetLabelSize(0.085);
            fr->Draw("AXIS");
            TLine* l1 = new TLine(0, 1, kNCut, 1);
            l1->SetLineColor(kBlack); l1->SetLineStyle(2); l1->SetLineWidth(3); l1->Draw();
            for (TGraphErrors* g : gg) g->Draw("P same");
            fr->Draw("AXIS same");
            Txt(0.12, 0.905, Form("%s (%s)", t.tex, t.w == wPeak ? "peak window" : "bump window"), 0.072);
            if (col == 0) Txt(0.97, 0.80, "Q > 1: the cut carries information", 0.06, 31, kGray + 2);
            else Txt(0.97, 0.22, "#varepsilon_{S} < 1: signal lost", 0.06, 31, kGray + 2);
        }
    c->cd();
    TPad* pl = new TPad("pl_rank", "", 0, yT, 1, 1);
    pl->SetMargin(0, 0, 0, 0); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.01, 0.05, 0.69, 0.95);
    lg->SetTextFont(43); lg->SetTextSize(32); lg->SetNColumns(3); lg->SetFillStyle(0); lg->SetBorderSize(0);
    for (int k = 0; k < 3; ++k) {
        TGraphErrors* g = new TGraphErrors();
        g->SetMarkerStyle(ms[k]); g->SetMarkerColor(tcol[k]); g->SetLineColor(tcol[k]); g->SetMarkerSize(2.6);
        lg->AddEntry(g, tyLeg[k], "p");
    }
    lg->Draw();
    TxtPx(0.70, 0.55, "SIMULATION ONLY, no CB; std + lep + eta, mode A,", 28);
    TxtPx(0.70, 0.12, "w(PDG) #times phys V3; X17: BR = 10^{-5}; errors: event groups", 28);
    SaveFig(c, "s5_cut_rank");
}

// ============================================================================
// Markdown key tables for the report (optional)
// ============================================================================
void WriteMd() {
    if (!gMd) return;
    M("<!-- written by plots/eta_sim/s5_cut.C -->\n\n");
    for (int w = 0; w < 2; ++w) {
        const std::vector<int> cs = w == wBump ? std::vector<int>{gS2, gS3, gS5, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2}
                                               : std::vector<int>{gS4, gS5, gX17, gEta, gPi0D, gConv, gVec, gOth, gFake, gDiff, aT2, aBX};
        M("### %s: sigma [nb] along the sequential flow\n\n| component | pre |", kWTxt[w]);
        for (int c = 0; c < kNCut; ++c) M(" +%s |", kCutName[c]);
        M(" N_eff pre / end |\n|---|---|");
        for (int c = 0; c < kNCut; ++c) M("---|");
        M("---|\n");
        for (int cc : cs) {
            const Acc& a0 = E(w, Cfg(tAlone, 0), cc);
            M("| %s | %s |", kCShort[cc], VEg(a0.xsA(), a0.xsAe()).c_str());
            for (int c = 0; c < kNCut; ++c) { const Acc& a = E(w, Cfg(tCum, c), cc); M(" %s |", VEg(a.xsP(), a.xsPe()).c_str()); }
            M(" %.0f / %.0f |\n", a0.NeffA(), E(w, Cfg(tCum, kNCut - 1), cc).NeffP());
        }
        for (int t : {(int)tAlone, (int)tNm1S, (int)tStep}) {
            M("\n### %s: efficiency, %s\n\n| component |", kWTxt[w], kTyName[t]);
            for (int c = 0; c < kNCut; ++c) M(" %s |", kCutName[c]);
            M("\n|---|");
            for (int c = 0; c < kNCut; ++c) M("---|");
            M("\n");
            for (int cc : cs) {
                M("| %s |", kCShort[cc]);
                for (int c = 0; c < kNCut; ++c) { const Acc& a = E(w, Cfg(t, c), cc); M(" %s |", a.Sa > 0 ? Fs("%.3f(%.0f)", a.eps(), 1000 * a.err()).c_str() : "-"); }
                M("\n");
            }
        }
        M("\n");
    }
    // ranking
    for (int it = 0; it < 3; ++it) {
        const Tgt& t = kTgt[it];
        M("\n### ranking, %s (%s)\n\n| cut | alone: eps_S / eps_B / Q | N-1 (slides): eps_S / eps_B / Q | N-1 (all): eps_S / eps_B / Q | verdict (N-1 slides) |\n|---|---|---|---|---|\n",
          t.name, kWTxt[t.w]);
        for (int c = 0; c < kNCut; ++c) {
            M("| %s |", kCutName[c]);
            QV vs{};
            for (int ty : {(int)tAlone, (int)tNm1S, (int)tNm1A}) {
                const QV v = Qof(t.w, Cfg(ty, c), t.s, t.b);
                if (ty == tNm1S) vs = v;
                M(" %.3f(%.0f) / %.3f(%.0f) / %.3f(%.0f) |", v.es, 1000 * v.ees, v.eb, 1000 * v.eeb, v.q, 1000 * v.eq);
            }
            M(" %s |\n", Verdict(vs));
        }
    }
    // where in M4 each cut acts
    for (int t : {(int)tAlone, (int)tNm1S}) {
        M("\n### T2 efficiency per M4 region: %s\n\n| cut |", kTyName[t]);
        for (int k = 0; k < 7; ++k) M(" %s |", kRName[k]);
        M("\n|---|");
        for (int k = 0; k < 7; ++k) M("---|");
        M("\n");
        for (int c = 0; c <= kNCut; ++c) {
            if (c == kNCut && t != tAlone) continue;
            M("| %s |", c == kNCut ? "slides (4 cuts)" : kCutName[c]);
            for (int k = 0; k < 7; ++k) {
                const Acc& a = c == kNCut ? E(kRegW[k], kCfgSlides, aT2) : E(kRegW[k], Cfg(t, c), aT2);
                M(" %s |", a.eps() < 0.01 ? Fs("%.2g(%.1g)", a.eps(), a.err()).c_str() : Fs("%.3f(%.0f)", a.eps(), 1000 * a.err()).c_str());
            }
            M("\n");
        }
    }
    // composition summary
    for (int s = 0; s < 2; ++s) {
        M("\n### composition, %s: sigma [nb] and shares of T2 per group; T3, T4 relative to T1+T2\n\n", s == 0 ? "pre" : "slides");
        M("| M4 region | T1+T2 [nb] | T2 [nb] | signal share of T1+T2 | G_eta | G_pi0D | G_conv | G_vec | G_oth | T3/(T1+T2) | T4/(T1+T2) | dominant ee | dominant pions |\n");
        M("|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
        for (int rg = 0; rg <= kNR; ++rg) {
            const RA& sim = C(s, rg, aSim);
            const RA& t2 = C(s, rg, aT2);
            int bl = 0, bp = 0;
            for (int l = 0; l < 12; ++l) if (C(s, rg, ItemL(2, l)).S > C(s, rg, ItemL(2, bl)).S) bl = l;
            for (int p = 0; p < 10; ++p) if (C(s, rg, ItemP(2, p)).S > C(s, rg, ItemP(2, bp)).S) bp = p;
            const std::pair<double, double> sg = Share(s, rg, aSig, 0);
            M("| %s | %s | %s | %s |", kRName[rg], VEg(sim.S / kLsim, std::sqrt(sim.aa) / kLsim).c_str(),
              VEg(t2.S / kLsim, std::sqrt(t2.aa) / kLsim).c_str(), Fs("%.4f(%.0f)", sg.first, 1e4 * sg.second).c_str());
            for (int g = gEta; g <= gOth; ++g) { const std::pair<double, double> x = Share(s, rg, g, 1); M(" %.3f(%.0f) |", x.first, 1000 * x.second); }
            M(" %.3f | %.3f | %s %.2f | %s %.2f |\n", Share(s, rg, gFake, 0).first, Share(s, rg, gDiff, 0).first, LclsName(bl),
              Share(s, rg, ItemL(2, bl), 1).first, PclsName(bp), Share(s, rg, ItemP(2, bp), 1).first);
        }
    }
}

}  // namespace s5c

// ============================================================================
void s5_cut(int parts = 15, const char* mdOut = "") {
    using namespace s5c;
    TStopwatch sw;
    Style();
    InitCfg();
    InitVars();
    if (!simphys::Init()) { printf("simphys::Init failed\n"); return; }
    simphys::Install();
    if (!simphys::LoadAngles()) { printf("simphys::LoadAngles failed\n"); return; }
    {
        SimData* d = new SimData;
        if (!Load(*d, kCompact, kMaskUS, -1e9, 1e9, 1)) return;
        simphys::AngHit() = 0; simphys::AngMiss() = 0;
        if (!Build(*d)) return;
        printf("V3 angles: found %lld, missing %lld\n", simphys::AngHit(), simphys::AngMiss());
        EvalEff();
        Checks(*d);
        delete d;
    }
    if (mdOut && mdOut[0]) gMd = fopen(mdOut, "w");
    std::vector<Job> jobs;
    if (parts & 1 || parts & 4) {
        for (int v = 0; v < 8; ++v)
            for (int w = 0; w < 2; ++w) {
                Job j;
                j.var = v; j.base = v < 4 ? (kSlM & ~Bit(v)) : kSlM; j.mlo = kWLo[w]; j.mhi = kWHi[w];
                j.Init(gVar[v].nb, gVar[v].lo, gVar[v].hi);
                jobs.push_back(j);
            }
        for (int s = 0; s < 2; ++s) {
            Job j;
            j.var = 8; j.base = s == 0 ? 0 : kSlM;
            j.Init(gVar[8].nb, gVar[8].lo, gVar[8].hi);
            jobs.push_back(j);
        }
        RunJobs(jobs);
    }
    if (parts & 1) {
        WriteNm1Table(jobs);
        for (int v = 0; v < 8; ++v) DrawNm1(v, jobs[2 * v], jobs[2 * v + 1]);
    }
    if (parts & 2) {
        WriteEffTable();
        DrawEffWindow(wBump);
        DrawEffWindow(wPeak);
        DrawFlow();
        DrawRegions();
    }
    if (parts & 4) {
        RunComp();
        WriteCompTable();
        DrawComp();
        DrawM4(jobs[16], jobs[17]);
    }
    if (parts & 8) {
        WriteRankTable();
        DrawRank();
    }
    if (gMd) {
        if (!(parts & 4)) RunComp();
        WriteMd();
        fclose(gMd); gMd = nullptr;
    }
    printf("s5_cut: done in %.1f s (checks: %d, outside tolerance: %d)\n", sw.RealTime(), gNchk, gNbad);
}
