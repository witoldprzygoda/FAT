// ============================================================================
// s5_opt_x17.C  --  selection for an X17 estimate in eta -> pi+pi- X17, X17 -> e+e- (simulation only, truth known)
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-. SMASH simulation, std + lep + eta sets, owned rows in partition mode A,
// weight w (PDG) times the physics weights of the eta decays (sim_eta_physics.h, main variant V3 = simphys::kPhysMain).
// NO measured data are read (no ScanEngine data loading); no combinatorial background exists in simulation.
//
// Signal  : X17 rows (four-particle class 6, eta set), W = 2 BR per row (BR = BR(eta -> pi+pi- X17) x BR(X17 -> e+e-)).
// Background: every other row of templates 1 + 2 (with the physics weights), split into the agreed groups
//   Gsig2/3/4/5 (signal classes 2-5), G_eta, G_pi0D, G_conv, G_vec, G_oth (template 2), and as a variant G_fake (T3)
//   and G_diff (T4).
// Figure of merit (v2, 2026-11-13, after the review findings H2 and M2): expected 90 % CL upper limit on BR from
//   counting, UL = sqrt(UL_stat^2 + (1.645 sigma_B,MC / s1)^2): UL_stat = Poisson CLs with continuous n = B (Asimov;
//   1.64 sqrt(B) for large B) and sigma_B,MC = the simulation (event-group) error of B, because the background of a
//   counting analysis is taken from this simulation; s1 = selected X17 rows per unit BR at L_data = 6080.5 nb^-1
//   (UNCONFIRMED). gFomMC = 0 restores v1 (B known exactly). The M4 window is FIXED to the study peak window
//   0.5333-0.5587 (v1 optimised it, a rule deviation); its edges are only scanned at the optimum (extra variables).
//   Also: median expected CLs (integer n), Feldman-Cousins, Z for BR = 1e-5 (all with B known).
// Optimisation: coordinate descent over the cut grids, choosing in every one-dimensional step the loosest value
//   whose figure of merit is within 2 % of the step minimum (plateau rule). Split-half check: event groups with even
//   grp (half A) and odd grp (half B) of each e+e- tree (grp = index of the event group in the tree, file order);
//   halves are evaluated with weights x 2.
// Shape analysis: binned Poisson likelihood of M(ee) in the M4 window (background shape from simulation, X17 as a
//   Gaussian line or as the simulated template), Asimov B-only data, CLs median expected limit from the asymptotic
//   formula q(mu_up) = 2.706, background normalisation known or profiled.
//
// Parts (bit mask, default 0 = all): 1 line shape, 2 background distributions, 4 optimisation + split half,
//   8 final numbers and shape analysis, 16 figures. Every part needs the loading step (about 2 min, 4 GB).
// Outputs: plots/eta_sim/figs/s5_opt_x17_*.{png,pdf}, plots/eta_sim/tables/s5_opt_x17_*.txt,
//   plots/eta_sim/tables/s5_opt_x17_frozen.txt.
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_opt_x17/aclic", kTRUE)' 'plots/eta_sim/s5_opt_x17.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Math/ProbFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TBox.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace s5x {
using namespace simeng;

// ---------------------------------------------------------------- output
FILE* gF = nullptr;
void P(const char* fmt, ...) {
    char buf[32768];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    fputs(buf, stdout);
    if (gF) fputs(buf, gF);
}
void OpenTab(const char* name) {
    const std::string p = std::string(kTabDir) + "/" + name;
    gF = fopen(p.c_str(), "w");
    if (!gF) printf("cannot write %s\n", p.c_str());
}
void CloseTab() { if (gF) fclose(gF); gF = nullptr; }
int gCnt = 0;
std::string Nm(const char* s) { return Form("%s_%d", s, gCnt++); }
// vector PDF with a tight page (EPS -> ps2pdf -dEPSCrop) and a PNG rasterised from that PDF (ROOT's own PNG renderer
// misplaces TLatex sub- and superscripts in NDC text)
void Save(TCanvas* c, const char* name) {
    const std::string base = std::string(kFigDir) + "/" + name;
    c->SaveAs((base + ".eps").c_str());
    gSystem->Exec(Form("cd %s && ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps && pdftoppm -png -singlefile -scale-to-x %d -scale-to-y -1 %s.pdf %s",
                       kFigDir, name, name, name, c->GetWw(), name, name));
}

// ---------------------------------------------------------------- constants
const double kLd = kLdata;                // 6080.5 nb^-1, UNCONFIRMED
const double kToRows = kLd / kLsim;       // sum w -> expected data rows
const char* const kLdNote = "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)";
const double kPkLo = 0.5333, kPkHi = 0.5587;   // M4 peak window of this study (simulation-centred)
const double kBuLo = 0.4576, kBuHi = 0.5255;   // M4 bump window
const double kMX = 17.0;                        // X17 mass [MeV/c2]
const double kNo = -999.;                       // "no cut"
const double kBRill = 1e-5;                     // illustrative BR for figures and Z
const double kTol = 0.02;                       // plateau rule
const double kM4Lo = 0.40, kM4Hi = 0.70;        // M4 range loaded for the lep set
inline bool On(double v) { return v > -998.; }

// ---------------------------------------------------------------- categories
enum Cat { cX17 = 0, cS2, cS3, cS4, cS5, cEta, cPi0D, cConv, cVec, cOth, cFake, cDiff, cX17o, kNCat };
const char* kCatName[kNCat] = {"X17 (class 6)", "Gsig2 bump pi0 Dalitz", "Gsig3 bump conversion", "Gsig4 peak conversion",
                               "Gsig5 direct eta->pipiee", "G_eta", "G_pi0D", "G_conv", "G_vec", "G_oth",
                               "G_fake (T3)", "G_diff (T4)", "X17 other rows"};
const char* kCatTex[kNCat] = {"X17 (class 6)", "Gsig2: bump, #pi^{0} Dalitz", "Gsig3: bump, conversion",
                              "Gsig4: peak, conversion", "Gsig5: direct #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}",
                              "G_{#eta}: #eta produced, no signal", "G_{#pi^{0}D}: #pi^{0} Dalitz pair",
                              "G_{conv}: conversion pair", "G_{vec}: #rho^{0}/#omega/#phi pairs", "G_{oth}: other",
                              "G_{fake} (T3, variant)", "G_{diff} (T4, variant)", "X17 other rows"};
// colours (v2): eta signals warm (orange, gold, dark red, red), X17 magenta; backgrounds cold or grey
const int kCatCol[kNCat] = {kMagenta + 1, kOrange + 1, kOrange - 2, kRed + 3, kRed, kViolet + 1, kGray + 1,
                            kTeal - 6, kGreen + 3, kGray, kAzure - 9, kGray + 3, kPink + 1};
inline bool IsBkg(int c) { return c >= cS2 && c <= cOth; }
inline bool IsBkg34(int c) { return c >= cS2 && c <= cDiff; }
int CatOf(const Row& r) {
    if (r.tmpl == 1 && r.cls == 6) return cX17;
    if (r.set == 3 && r.w == 0.) return cX17o;     // eta-set row with an X17 forced-decay track (BR(X17) = 0 in production)
    if (r.tmpl == 1) {
        if (r.cls >= 2 && r.cls <= 5) return cS2 + (r.cls - 2);
        return -1;
    }
    if (r.tmpl == 2) {
        if (r.pcls == 2 || r.lcls == 2) return cEta;
        if (r.lcls == 1) return cPi0D;
        if (r.lcls == 3) return cConv;
        if (r.lcls >= 4 && r.lcls <= 6) return cVec;
        return cOth;
    }
    if (r.tmpl == 3) return cFake;
    if (r.tmpl == 4) return cDiff;
    return -1;
}

// ---------------------------------------------------------------- compact rows
struct CR {
    float m4, mee, oa, lab, eta, mpp, mm4, pt4, dz, pee, rvee, docaee;   // mee [MeV/c2], pee [MeV/c]; rv_ee, doca_ee [mm]
    unsigned char vtx, cat, half, set;
    int gid;
    double w;    // background: w (PDG) x physics factor; X17 rows: 2 (per unit BR)
    double w0;   // the same without the physics weights (variant)
};
bool gW0 = false;   // evaluate with w0 (physics weights OFF, variant)
std::vector<CR> gAll;      // pre, owned (mode A), M4 in [0.40, 0.70)
std::vector<CR> gScan;     // gAll with M4 in [0.50, 0.60) and M(ee) < 60 MeV (optimisation)
std::vector<CR> gX17all;   // X17 rows, pre, all M4 (eta set)

// ---------------------------------------------------------------- cuts
struct Cuts {
    double m4lo = kPkLo, m4hi = kPkHi;
    double meelo = kNo, meehi = kNo;     // MeV/c2
    double oalo = kNo, oahi = kNo;       // deg
    double labhi = kNo;                  // deg
    double mpplo = kNo, mpphi = kNo;     // GeV/c2
    double etalo = kNo;                  // deg
    double vtx = kNo;                    // 1 = both pions in the vertex fit
    double dzhi = kNo;                   // mm
    double mm4hi = kNo;                  // GeV/c2
    double pt4lo = kNo;                  // MeV/c
    double rveehi = kNo;                 // mm, e+e- crossing point, transverse distance from the vertex
    double docaeehi = kNo;               // mm, distance of closest approach of the two lepton tracks
};
inline bool PassC(const CR& r, const Cuts& c, bool useMee = true, bool useM4 = true) {
    if (useM4 && !(r.m4 >= c.m4lo && r.m4 < c.m4hi)) return false;
    if (useMee && On(c.meelo) && !(r.mee >= c.meelo)) return false;
    if (useMee && On(c.meehi) && !(r.mee < c.meehi)) return false;
    if (On(c.oalo) && !(r.oa > c.oalo)) return false;
    if (On(c.oahi) && !(r.oa < c.oahi)) return false;
    if (On(c.labhi) && !(r.lab < c.labhi)) return false;
    if (On(c.mpplo) && !(r.mpp > c.mpplo)) return false;
    if (On(c.mpphi) && !(r.mpp < c.mpphi)) return false;
    if (On(c.etalo) && !(r.eta > c.etalo)) return false;
    if (On(c.vtx) && !r.vtx) return false;
    if (On(c.dzhi) && !(r.dz < c.dzhi)) return false;
    if (On(c.mm4hi) && !(r.mm4 < c.mm4hi)) return false;
    if (On(c.pt4lo) && !(r.pt4 > c.pt4lo)) return false;
    if (On(c.rveehi) && !(r.rvee >= 0 && r.rvee < c.rveehi)) return false;
    if (On(c.docaeehi) && !(r.docaee < c.docaeehi)) return false;
    return true;
}
// the same selection on engine rows (pre included)
inline Sel SelOf(const Cuts& c, bool useMee = true, bool useM4 = true) {
    return [c, useMee, useM4](const Row& r) {
        if (!PassPre(r)) return false;
        CR x;
        x.m4 = r.m4; x.mee = r.mee * 1000.f; x.oa = r.oa; x.lab = r.oa_lab; x.eta = r.oa_eta; x.mpp = r.mpipi;
        x.mm4 = r.mm4; x.pt4 = r.pt4; x.dz = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
        x.vtx = (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1);
        x.rvee = r.rv_ee; x.docaee = r.doca_ee;
        return PassC(x, c, useMee, useM4);
    };
}
std::string CutStr(const Cuts& c) {
    std::string s = Form("M4 %.4f-%.4f", c.m4lo, c.m4hi);
    if (On(c.meelo) || On(c.meehi)) s += Form(", M(ee) %s-%s MeV", On(c.meelo) ? Form("%.1f", c.meelo) : "0", On(c.meehi) ? Form("%.1f", c.meehi) : "inf");
    if (On(c.oalo)) s += Form(", OA > %.1f", c.oalo);
    if (On(c.oahi)) s += Form(", OA < %.1f", c.oahi);
    if (On(c.labhi)) s += Form(", a_lab < %.0f", c.labhi);
    if (On(c.mpplo)) s += Form(", M(pipi) > %.3f", c.mpplo);
    if (On(c.mpphi)) s += Form(", M(pipi) < %.3f", c.mpphi);
    if (On(c.etalo)) s += Form(", a_eta > %.0f", c.etalo);
    if (On(c.vtx)) s += ", pions in vertex";
    if (On(c.dzhi)) s += Form(", dz_pi < %.0f", c.dzhi);
    if (On(c.mm4hi)) s += Form(", MM4 < %.2f", c.mm4hi);
    if (On(c.pt4lo)) s += Form(", pT4 > %.0f", c.pt4lo);
    if (On(c.rveehi)) s += Form(", rv_ee < %.0f", c.rveehi);
    if (On(c.docaeehi)) s += Form(", doca_ee < %.1f", c.docaeehi);
    return s;
}
Cuts SlidesCuts() {   // preset kSlides in the peak window
    Cuts c; c.labhi = 50; c.mpphi = 0.46; c.etalo = 140; c.oalo = 2; return c;
}

// ---------------------------------------------------------------- evaluation
struct Ev {
    double s = 0, vs = 0, ns = 0;                 // X17 sum (per unit BR), group variance, rows
    double b[kNCat] = {0}, vb[kNCat] = {0}, nb[kNCat] = {0};
    double B = 0, VB = 0, B34 = 0, VB34 = 0;      // T1 (no X17) + T2; + T3 + T4
    double NeffB() const { return VB > 0 ? B * B / VB : 0.; }
    double NeffB34() const { return VB34 > 0 ? B34 * B34 / VB34 : 0.; }
};
// half: 0 all, 1 half A (even grp), 2 half B (odd grp); halves are scaled by 2
Ev Eval(const std::vector<CR>& rows, const Cuts& c, int half = 0, bool useMee = true, bool useM4 = true) {
    Ev e;
    const double f = half ? 2. : 1.;
    double a[kNCat];
    std::fill(a, a + kNCat, 0.);
    bool any = false;
    int cur = -1;
    auto flush = [&]() {
        if (!any) return;
        double sb = 0, sb34 = 0;
        for (int k = 0; k < kNCat; ++k) {
            if (a[k] == 0.) continue;
            e.vb[k] += a[k] * a[k];
            if (IsBkg(k)) sb += a[k];
            if (IsBkg34(k)) sb34 += a[k];
        }
        e.vs += a[cX17] * a[cX17];
        e.VB += sb * sb; e.VB34 += sb34 * sb34;
        std::fill(a, a + kNCat, 0.);
        any = false;
    };
    for (const CR& r : rows) {
        if (r.gid != cur) { flush(); cur = r.gid; }
        if (half && r.half != half - 1) continue;
        if (!PassC(r, c, useMee, useM4)) continue;
        const double w = (gW0 ? r.w0 : r.w) * f;
        a[r.cat] += w; any = true;
        e.b[r.cat] += w; e.nb[r.cat] += 1;
    }
    flush();
    e.s = e.b[cX17]; e.ns = e.nb[cX17];
    e.vb[cX17] = e.vs;
    for (int k = 0; k < kNCat; ++k) { if (IsBkg(k)) e.B += e.b[k]; if (IsBkg34(k)) e.B34 += e.b[k]; }
    return e;
}

// ---------------------------------------------------------------- statistics
// P(N <= n | mu) for real n >= 0: regularised upper incomplete gamma Q(n + 1, mu)
inline double PCdf(double n, double mu) { return ROOT::Math::inc_gamma_c(n + 1., mu); }
// CLs upper limit on s for observed n (real) and known b at confidence cl
double ClsUL(double n, double b, double cl = 0.90) {
    if (b < 0) b = 0;
    const double den = PCdf(n, b);
    if (!(den > 0)) return 0.;
    auto f = [&](double s) { return PCdf(n, s + b) / den - (1. - cl); };
    double lo = 0, hi = std::max(5., 3. * std::sqrt(n + b + 1.) + std::max(0., n - b) + 5.);
    int it = 0;
    while (f(hi) > 0 && it < 200) { hi *= 2; ++it; }
    for (int i = 0; i < 200; ++i) { const double m = 0.5 * (lo + hi); if (f(m) > 0) lo = m; else hi = m; }
    return 0.5 * (lo + hi);
}
// asymptotic (Cowan et al.) median expected CLs limit for one bin: 2 [s - b ln(1 + s/b)] = 2.706
double AsyUL(double b) {
    if (b <= 0) return 0.;
    const double q = 2.705543;
    auto f = [&](double s) { return 2. * (s - b * std::log1p(s / b)) - q; };
    double lo = 0, hi = 10. + 5. * std::sqrt(b);
    while (f(hi) < 0) hi *= 2;
    for (int i = 0; i < 200; ++i) { const double m = 0.5 * (lo + hi); if (f(m) < 0) lo = m; else hi = m; }
    return 0.5 * (lo + hi);
}
// Poisson quantile: smallest integer n with P(N <= n | b) >= p
int PQuant(double b, double p) {
    int n = std::max(0, (int)std::floor(b - 10 * std::sqrt(b + 1.)));
    while (PCdf(n, b) < p) ++n;
    while (n > 0 && PCdf(n - 1, b) >= p) --n;
    return n;
}
// Feldman-Cousins (known b): lower edge of the acceptance region in n at signal mu
inline double LogP(int n, double m) { return m > 0 ? n * std::log(m) - m - std::lgamma(n + 1.) : (n == 0 ? 0. : -1e300); }
int FCLow(double mu, double b, double cl) {
    const double m = mu + b;
    const int nlo = std::max(0, (int)std::floor(m - 12 * std::sqrt(m + 1.) - 10));
    const int nhi = (int)std::ceil(m + 12 * std::sqrt(m + 1.) + 10);
    std::vector<std::pair<double, int>> R;
    R.reserve(nhi - nlo + 1);
    for (int n = nlo; n <= nhi; ++n) {
        const double mub = std::max(0., n - b);
        R.push_back(std::make_pair(LogP(n, m) - LogP(n, mub + b), n));
    }
    std::sort(R.begin(), R.end(), [](const std::pair<double, int>& x, const std::pair<double, int>& y) { return x.first > y.first; });
    double acc = 0; int nmin = 1 << 30;
    for (const auto& x : R) { acc += std::exp(LogP(x.second, m)); nmin = std::min(nmin, x.second); if (acc >= cl) break; }
    return nmin;
}
double FCUpper(int n0, double b, double cl = 0.90) {
    double lo = 0, hi = std::max(5., n0 - b + 10 * std::sqrt(n0 + b + 1.) + 10);
    while (FCLow(hi, b, cl) <= n0) hi *= 2;
    for (int i = 0; i < 60; ++i) { const double m = 0.5 * (lo + hi); if (FCLow(m, b, cl) > n0) hi = m; else lo = m; }
    return 0.5 * (lo + hi);
}
// Asimov discovery significance for s on known b
inline double ZA(double s, double b) { return b > 0 ? std::sqrt(std::max(0., 2 * ((s + b) * std::log1p(s / b) - s))) : 0.; }

// figure of merit: expected 90 % CL UL on BR (Asimov CLs with continuous n = B)
// v2: ul = UL including the MC error of B (gFomMC = 1, default) or UL_stat (B known exactly, gFomMC = 0)
int gFomMC = 1;
struct Fom { double ul = 1e9, s1 = 0, B = 0, sup = 0, ulStat = 1e9, ulTot = 1e9, sBmc = 0; };
Fom FomOf(const Ev& e, bool with34 = false) {
    Fom f;
    f.s1 = e.s * kToRows;                       // expected rows per unit BR
    f.B = (with34 ? e.B34 : e.B) * kToRows;     // expected background rows
    f.sBmc = std::sqrt(with34 ? e.VB34 : e.VB) * kToRows;   // simulation (event-group) error of B in rows
    if (f.s1 <= 0) return f;
    f.sup = ClsUL(f.B, f.B);
    f.ulStat = f.sup / f.s1;
    f.ulTot = std::sqrt(f.ulStat * f.ulStat + std::pow(1.6449 * f.sBmc / f.s1, 2));
    f.ul = gFomMC ? f.ulTot : f.ulStat;
    return f;
}
double FomCuts(const std::vector<CR>& rows, const Cuts& c, int half = 0) { return FomOf(Eval(rows, c, half)).ul; }

// ---------------------------------------------------------------- optimisation variables
struct OVar {
    const char* key;
    const char* label;      // axis label
    double Cuts::*field;
    int dir;                // +1: larger value looser (upper cut); -1: smaller value looser (lower cut)
    std::vector<double> vals;
    double scale;           // display scale
};
std::vector<double> Rng(double a, double b, double d) {
    std::vector<double> v;
    for (double x = a; x <= b + 1e-9; x += d) v.push_back(std::round(x * 1e6) / 1e6);
    return v;
}
std::vector<OVar> gVars;
void InitVars() {
    gVars.clear();
    // v2: the M4 window is fixed to the study peak window (kPkLo, kPkHi); its edges are extra variables (scan only)
    gVars.push_back({"meelo", "M(e^{+}e^{-}) lower edge [MeV/c^{2}]", &Cuts::meelo, -1, Rng(6, 17, 0.5), 1});
    gVars.push_back({"meehi", "M(e^{+}e^{-}) upper edge [MeV/c^{2}]", &Cuts::meehi, +1, Rng(17, 34, 0.5), 1});
    {
        std::vector<double> v = {kNo}; for (double x : Rng(0.5, 8, 0.5)) v.push_back(x);
        gVars.push_back({"oalo", "OA lower cut [deg]", &Cuts::oalo, -1, v, 1});
    }
    {
        std::vector<double> v = Rng(5, 24, 1); v.push_back(kNo);
        gVars.push_back({"oahi", "OA upper cut [deg]", &Cuts::oahi, +1, v, 1});
    }
    {
        std::vector<double> v = Rng(20, 90, 5); v.push_back(kNo);
        gVars.push_back({"labhi", "#alpha_{lab} upper cut [deg]", &Cuts::labhi, +1, v, 1});
    }
    {
        std::vector<double> v = {kNo}; for (double x : Rng(0.28, 0.40, 0.01)) v.push_back(x);
        gVars.push_back({"mpplo", "M(#pi^{+}#pi^{-}) lower cut [GeV/c^{2}]", &Cuts::mpplo, -1, v, 1});
    }
    {
        std::vector<double> v = Rng(0.40, 0.53, 0.01); v.push_back(kNo);
        gVars.push_back({"mpphi", "M(#pi^{+}#pi^{-}) upper cut [GeV/c^{2}]", &Cuts::mpphi, +1, v, 1});
    }
    {
        std::vector<double> v = {kNo}; for (double x : Rng(90, 175, 5)) v.push_back(x);
        gVars.push_back({"etalo", "#alpha_{#eta} lower cut [deg]", &Cuts::etalo, -1, v, 1});
    }
    gVars.push_back({"vtx", "both pions in the vertex fit (0 no, 1 yes)", &Cuts::vtx, -1, {kNo, 1.}, 1});
}
// extra variables, scanned at the optimum only (not part of the descent)
std::vector<OVar> gExtra;
void InitExtra() {
    gExtra.clear();
    { std::vector<double> v = Rng(0.5180, 0.5455, 0.0025); v.push_back(kPkLo); std::sort(v.begin(), v.end());
      gExtra.push_back({"m4lo", "M_{4} window, lower edge [GeV/c^{2}]", &Cuts::m4lo, -1, v, 1}); }
    { std::vector<double> v = Rng(0.5480, 0.5755, 0.0025); v.push_back(kPkHi); std::sort(v.begin(), v.end());
      gExtra.push_back({"m4hi", "M_{4} window, upper edge [GeV/c^{2}]", &Cuts::m4hi, +1, v, 1}); }
    { std::vector<double> v = {10, 20, 30, 40, 60, 100}; v.push_back(kNo); gExtra.push_back({"dzhi", "|z_{#pi} - v_{z}| upper cut (both #pi) [mm]", &Cuts::dzhi, +1, v, 1}); }
    { std::vector<double> v = Rng(2.0, 3.0, 0.05); v.push_back(kNo); gExtra.push_back({"mm4hi", "MM_{4} upper cut [GeV/c^{2}]", &Cuts::mm4hi, +1, v, 1}); }
    { std::vector<double> v = {kNo}; for (double x : Rng(50, 500, 50)) v.push_back(x); gExtra.push_back({"pt4lo", "p_{T,4} lower cut [MeV/c]", &Cuts::pt4lo, -1, v, 1}); }
    { std::vector<double> v = {2, 4, 6, 8, 10, 15, 20, 30, 50}; v.push_back(kNo); gExtra.push_back({"rveehi", "e^{+}e^{-} crossing point, r from vertex, upper cut [mm]", &Cuts::rveehi, +1, v, 1}); }
    { std::vector<double> v = {1, 2, 3, 5, 8, 12, 20}; v.push_back(kNo); gExtra.push_back({"docaeehi", "doca(e^{+},e^{-}) upper cut [mm]", &Cuts::docaeehi, +1, v, 1}); }
}
inline double Loose(const OVar& v, double x) { return !On(x) ? 1e30 : (v.dir > 0 ? x : -x); }

struct StepLog { int pass; std::string var; double from, to, fom; };
// coordinate descent with the plateau rule; returns the cuts, appends the path to lg
Cuts Descent(const std::vector<CR>& rows, Cuts c, int half, std::vector<StepLog>* lg, int maxPass = 8) {
    for (int pass = 0; pass < maxPass; ++pass) {
        bool changed = false;
        for (const OVar& v : gVars) {
            std::vector<double> f(v.vals.size());
            double fmin = 1e30;
            for (size_t i = 0; i < v.vals.size(); ++i) {
                Cuts t = c; t.*(v.field) = v.vals[i];
                if (t.m4hi <= t.m4lo || (On(t.meelo) && On(t.meehi) && t.meehi <= t.meelo) ||
                    (On(t.oalo) && On(t.oahi) && t.oahi <= t.oalo) || (On(t.mpplo) && On(t.mpphi) && t.mpphi <= t.mpplo)) { f[i] = 1e30; continue; }
                f[i] = FomCuts(rows, t, half);
                fmin = std::min(fmin, f[i]);
            }
            int best = -1; double bl = -1e31;
            for (size_t i = 0; i < v.vals.size(); ++i)
                if (f[i] <= fmin * (1 + kTol) && Loose(v, v.vals[i]) > bl) { bl = Loose(v, v.vals[i]); best = (int)i; }
            if (best < 0) continue;
            const double old = c.*(v.field);
            if (std::fabs(old - v.vals[best]) > 1e-9) {
                c.*(v.field) = v.vals[best];
                changed = true;
                if (lg) lg->push_back({pass, v.key, old, v.vals[best], f[best]});
            }
        }
        if (!changed) break;
    }
    return c;
}

// ---------------------------------------------------------------- shape analysis (binned likelihood in M(ee))
struct ShapeIn {
    std::vector<double> b0, b1, s;   // rows per bin: b0 = Gsig5 (direct decay), b1 = other background; s = X17 per unit BR
    double lo = 0, hi = 0; int nb = 0;
};
// NLL = sum (lam - n ln lam); n ln lam = 0 for n = 0
double NLL(const ShapeIn& in, const std::vector<double>& n, double mu, double t0, double t1) {
    double L = 0;
    for (int i = 0; i < in.nb; ++i) {
        const double lam = t0 * in.b0[i] + t1 * in.b1[i] + mu * in.s[i];
        if (lam <= 0) { if (n[i] > 0) L += 1e10; continue; }
        L += lam - (n[i] > 0 ? n[i] * std::log(lam) : 0.);
    }
    return L;
}
// profile the background normalisations: mode 0 none, 1 one common factor, 2 two factors (direct decay, rest)
double ProfNLL(const ShapeIn& in, const std::vector<double>& n, double mu, int mode, double* t0o = nullptr, double* t1o = nullptr) {
    double t0 = 1, t1 = 1;
    if (mode == 1) {
        for (int it = 0; it < 100; ++it) {
            double g = 0, h = 0;
            for (int i = 0; i < in.nb; ++i) {
                const double B = in.b0[i] + in.b1[i];
                const double lam = t0 * B + mu * in.s[i];
                if (lam <= 0) continue;
                g += B * (1 - n[i] / lam); h += n[i] * B * B / (lam * lam);
            }
            if (h <= 0) break;
            const double d = g / h;
            t0 = std::max(0., t0 - d);
            if (std::fabs(d) < 1e-13) break;
        }
        t1 = t0;
    } else if (mode == 2) {
        for (int it = 0; it < 200; ++it) {
            double g0 = 0, g1 = 0, h00 = 0, h01 = 0, h11 = 0;
            for (int i = 0; i < in.nb; ++i) {
                const double lam = t0 * in.b0[i] + t1 * in.b1[i] + mu * in.s[i];
                if (lam <= 0) continue;
                g0 += in.b0[i] * (1 - n[i] / lam); g1 += in.b1[i] * (1 - n[i] / lam);
                const double q = n[i] / (lam * lam);
                h00 += q * in.b0[i] * in.b0[i]; h01 += q * in.b0[i] * in.b1[i]; h11 += q * in.b1[i] * in.b1[i];
            }
            const double det = h00 * h11 - h01 * h01;
            if (!(det > 0)) break;
            const double d0 = (h11 * g0 - h01 * g1) / det, d1 = (h00 * g1 - h01 * g0) / det;
            t0 = std::max(0., t0 - d0); t1 = std::max(0., t1 - d1);
            if (std::fabs(d0) < 1e-13 && std::fabs(d1) < 1e-13) break;
        }
    }
    if (t0o) *t0o = t0;
    if (t1o) *t1o = t1;
    return NLL(in, n, mu, t0, t1);
}
// q(mu) on the Asimov B-only data (mu-hat = 0, theta-hat = 1)
double QAsimov(const ShapeIn& in, double mu, int mode) {
    std::vector<double> n(in.nb);
    for (int i = 0; i < in.nb; ++i) n[i] = in.b0[i] + in.b1[i];
    return 2 * (ProfNLL(in, n, mu, mode) - NLL(in, n, 0, 1, 1));
}
double ShapeUL(const ShapeIn& in, int mode) {
    const double q = 2.705543;
    double lo = 0, hi = 1e-6;
    while (QAsimov(in, hi, mode) < q && hi < 1) hi *= 2;
    for (int i = 0; i < 100; ++i) { const double m = 0.5 * (lo + hi); if (QAsimov(in, m, mode) < q) lo = m; else hi = m; }
    return 0.5 * (lo + hi);
}
double ShapeZ(const ShapeIn& in, double mu1, int mode) {
    std::vector<double> n(in.nb);
    for (int i = 0; i < in.nb; ++i) n[i] = in.b0[i] + in.b1[i] + mu1 * in.s[i];
    const double q0 = 2 * (ProfNLL(in, n, 0, mode) - NLL(in, n, mu1, 1, 1));
    return q0 > 0 ? std::sqrt(q0) : 0.;
}
// histograms (rows at L_data) of the background and of X17 per unit BR in M(ee) for the cuts without the M(ee) window
void ShapeHists(const Cuts& c, double lo, double hi, int nb, ShapeIn& in, std::vector<double>* sTpl, std::vector<std::vector<double>>* grp) {
    in.lo = lo; in.hi = hi; in.nb = nb;
    in.b0.assign(nb, 0.); in.b1.assign(nb, 0.); in.s.assign(nb, 0.);
    if (sTpl) sTpl->assign(nb, 0.);
    if (grp) grp->assign(kNCat, std::vector<double>(nb, 0.));
    const double bw = (hi - lo) / nb;
    for (const CR& r : gAll) {
        if (!PassC(r, c, false, true)) continue;
        if (!(r.mee >= lo && r.mee < hi)) continue;
        const int i = std::min(nb - 1, (int)((r.mee - lo) / bw));
        const double w = r.w * kToRows;
        if (grp) (*grp)[r.cat][i] += w;
        if (r.cat == cX17) { if (sTpl) (*sTpl)[i] += w; continue; }
        if (r.cat == cS5) in.b0[i] += w;
        else if (IsBkg(r.cat)) in.b1[i] += w;
    }
}
void FillGauss(ShapeIn& in, double total, double mean, double sigma) {
    const double bw = (in.hi - in.lo) / in.nb;
    const double norm = ROOT::Math::gaussian_cdf(in.hi, sigma, mean) - ROOT::Math::gaussian_cdf(in.lo, sigma, mean);
    for (int i = 0; i < in.nb; ++i) {
        const double a = in.lo + i * bw, b = a + bw;
        in.s[i] = total * (ROOT::Math::gaussian_cdf(b, sigma, mean) - ROOT::Math::gaussian_cdf(a, sigma, mean)) / norm;
    }
}

// ---------------------------------------------------------------- line-shape fits
struct LineFit {
    double mean = 0, meanE = 0, sig = 0, sigE = 0;      // Gaussian core (iterative +-1.5 sigma)
    double q16 = 0, q50 = 0, q84 = 0, q025 = 0, q975 = 0;
    double f1 = 0, f2 = 0, f3 = 0, flo3 = 0, fhi3 = 0;  // fractions within +-1,2,3 sigma_core; tails beyond 3 sigma
    double m1 = 0, s1 = 0, m2 = 0, s2 = 0, fr2 = 0, chi2 = 0; int ndf = 0;   // double Gaussian
    double neff = 0, tot = 0;
};
// weighted quantile of (x, w) pairs
double WQuant(std::vector<std::pair<double, double>>& v, double p) {
    std::sort(v.begin(), v.end());
    double tot = 0; for (auto& x : v) tot += x.second;
    double acc = 0;
    for (auto& x : v) { acc += x.second; if (acc >= p * tot) return x.first; }
    return v.empty() ? 0. : v.back().first;
}
LineFit FitLine(const std::vector<std::pair<double, double>>& xw0, TH1D* h) {
    LineFit L;
    std::vector<std::pair<double, double>> xw = xw0;
    double tot = 0, t2 = 0; for (auto& x : xw) { tot += x.second; t2 += x.second * x.second; }
    L.tot = tot; L.neff = t2 > 0 ? tot * tot / t2 : 0;
    L.q16 = WQuant(xw, 0.1587); L.q50 = WQuant(xw, 0.5); L.q84 = WQuant(xw, 0.8413);
    L.q025 = WQuant(xw, 0.025); L.q975 = WQuant(xw, 0.975);
    double m = L.q50, s = 0.5 * (L.q84 - L.q16);
    TF1 g(Nm("g").c_str(), "gaus", 0, 100);
    for (int it = 0; it < 6; ++it) {
        g.SetParameters(h->GetMaximum(), m, s);
        h->Fit(&g, "QN0", "", m - 1.5 * s, m + 1.5 * s);
        m = g.GetParameter(1); s = std::fabs(g.GetParameter(2));
    }
    L.mean = m; L.meanE = g.GetParError(1); L.sig = s; L.sigE = g.GetParError(2);
    double a1 = 0, a2 = 0, a3 = 0, lo3 = 0, hi3 = 0;
    for (auto& x : xw) {
        const double d = (x.first - m) / s;
        if (std::fabs(d) < 1) a1 += x.second;
        if (std::fabs(d) < 2) a2 += x.second;
        if (std::fabs(d) < 3) a3 += x.second;
        if (d <= -3) lo3 += x.second;
        if (d >= 3) hi3 += x.second;
    }
    L.f1 = a1 / tot; L.f2 = a2 / tot; L.f3 = a3 / tot; L.flo3 = lo3 / tot; L.fhi3 = hi3 / tot;
    TF1 dg(Nm("dg").c_str(), "[0]*((1-[5])*TMath::Gaus(x,[1],[2],1)+[5]*TMath::Gaus(x,[3],[4],1))", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    dg.SetParameters(h->Integral("width"), m, s, m, 2.5 * s, 0.3);
    dg.SetParLimits(5, 0., 1.); dg.SetParLimits(2, 0.1, 20); dg.SetParLimits(4, 0.1, 60);
    h->Fit(&dg, "QN0", "", std::max(h->GetXaxis()->GetXmin(), m - 8 * s), std::min(h->GetXaxis()->GetXmax(), m + 8 * s));
    L.m1 = dg.GetParameter(1); L.s1 = std::fabs(dg.GetParameter(2)); L.m2 = dg.GetParameter(3); L.s2 = std::fabs(dg.GetParameter(4));
    L.fr2 = dg.GetParameter(5); L.chi2 = dg.GetChisquare(); L.ndf = dg.GetNDF();
    if (L.s2 < L.s1) { std::swap(L.m1, L.m2); std::swap(L.s1, L.s2); L.fr2 = 1 - L.fr2; }
    return L;
}

// ---------------------------------------------------------------- style
void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetLabelSize(0.05, "xyz"); gStyle->SetTitleSize(0.055, "xyz");
    gStyle->SetTitleOffset(1.15, "y"); gStyle->SetTitleOffset(1.05, "x");
    gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.09);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetNumberContours(60); gStyle->SetPalette(kBird);
}
void Txt(double x, double y, const char* s, double size = 0.048, int col = kBlack, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.SetTextFont(42);
    t.DrawLatex(x, y, s);
}
void Frame(TH1* h, const char* xt, const char* yt) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetTitleOffset(1.05); h->GetYaxis()->SetTitleOffset(1.25);
}

// ---------------------------------------------------------------- global state of the run
SimData* gD = nullptr;
SimOpt gOpt;
Cuts gOpt0, gOptA, gOptB;    // full-sample optimum, half A optimum, half B optimum
std::vector<StepLog> gLog, gLogA, gLogB;
LineFit gLinePk, gLineAll, gLineOpt;

// group components for the engine
std::vector<Comp> GroupComps(double brX17) {
    std::vector<Comp> v;
    for (int k = 0; k < kNCat; ++k) {
        Comp c{kCatName[k], [k](const Row& r) { return CatOf(r) == k; }};
        if (k == cX17 || k == cX17o) c.wf = [brX17](const Row&) { return 2. * brX17; };
        v.push_back(c);
    }
    return v;
}

// ============================================================================ loading
bool LoadAll() {
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    gOpt = SimOpt(); gOpt.sets = kSetsAll; gOpt.phys = simphys::kPhysMain; gOpt.mode = 0;
    gD = new SimData();
    SimData dl;
    if (!Load(*gD, kCompact, (1 << kStdUS) | (1 << kEtaUS))) return false;      // std and eta e+e- trees, all M4
    if (!Load(dl, kCompact, (1 << kLepUS), kM4Lo, kM4Hi)) return false;         // lep e+e- tree, M4 in [0.40, 0.70)
    gD->t[kLepUS] = std::move(dl.t[kLepUS]);
    gD->loaded |= (1 << kLepUS);
    printf("loaded in %.0f s\n", sw.RealTime());
    // compact rows
    long long nX17o = 0;
    for (int set = 1; set <= 3; ++set) {
        for (const Row& r : gD->tree(set, kUS)) {
            if (!Owned(r, gOpt)) continue;
            if (!PassPre(r)) continue;
            const int cat = CatOf(r);
            if (cat < 0) continue;
            CR c;
            c.m4 = r.m4; c.mee = r.mee * 1000.f; c.oa = r.oa; c.lab = r.oa_lab; c.eta = r.oa_eta; c.mpp = r.mpipi;
            c.mm4 = r.mm4; c.pt4 = r.pt4; c.dz = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz)); c.pee = r.p_ee;
            c.rvee = r.rv_ee; c.docaee = r.doca_ee;
            c.vtx = (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1);
            c.cat = (unsigned char)cat; c.half = (unsigned char)(r.grp & 1); c.set = (unsigned char)set;
            c.gid = r.grp + (set - 1) * (1 << 28);
            c.w = (cat == cX17 || cat == cX17o) ? 2. : r.w * simphys::RowFactor(r, simphys::kPhysMain);
            c.w0 = (cat == cX17 || cat == cX17o) ? 2. : r.w;
            if (cat == cX17o) ++nX17o;
            if (cat == cX17) gX17all.push_back(c);
            if (!(r.m4 >= kM4Lo && r.m4 < kM4Hi)) continue;
            if (!(c.w > 0) && !(c.w0 > 0)) continue;
            gAll.push_back(c);
            if (r.m4 >= 0.50f && r.m4 < 0.60f && c.mee < 60.f) gScan.push_back(c);
        }
    }
    printf("compact rows: all %zu, scan %zu, X17 all M4 %zu (X17-other rows %lld); %.0f s\n", gAll.size(), gScan.size(),
           gX17all.size(), nX17o, sw.RealTime());
    return true;
}

// ============================================================================ part 1: X17 line shape
void Part1() {
    OpenTab("s5_opt_x17_lineshape.txt");
    P("# s5_opt_x17_lineshape: reconstructed M(ee) line of X17 (eta -> pi+pi- X17, X17 -> e+e-, m = 17 MeV/c2; class 6, eta set)\n");
    P("# selection pre (isBest >= 0, vz > -500, PT3 bit, vz > -160, cor_ok); M(ee) from corrected momenta; per unit BR weights (W = 2)\n");
    P("# core = Gaussian fitted iteratively in +-1.5 sigma; q = weighted quantiles; f(k) = fraction within +-k sigma_core;\n");
    P("# double Gaussian over +-8 sigma (narrow 1 / wide 2, fr2 = wide fraction). Errors: fit errors (statistics of the class-6 rows).\n\n");
    // true M(ee)
    {
        double s = 0, s2 = 0, n = 0, mn = 1e9, mx = -1e9;
        for (const Row& r : gD->tree(3, kUS)) {
            if (!(r.tmpl == 1 && r.cls == 6) || !PassPre(r) || !r.ownA) continue;
            const double m = r.mee_true * 1000.;
            s += m; s2 += m * m; n += 1; mn = std::min(mn, m); mx = std::max(mx, m);
        }
        P("true M(ee) of the accepted X17 rows: mean %.4f MeV, rms %.4f MeV, min %.4f, max %.4f (rows %.0f)\n\n", s / n,
          std::sqrt(std::max(0., s2 / n - (s / n) * (s / n))), mn, mx, n);
    }
    struct Sl { const char* name; std::function<bool(const CR&)> f; };
    std::vector<Sl> sls = {
        {"pre, all M4", [](const CR&) { return true; }},
        {"pre, M4 peak window", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi; }},
        {"pre, M4 bump window", [](const CR& r) { return r.m4 >= kBuLo && r.m4 < kBuHi; }},
        {"peak window, OA < 3", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.oa < 3; }},
        {"peak window, OA 3-5", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.oa >= 3 && r.oa < 5; }},
        {"peak window, OA 5-8", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.oa >= 5 && r.oa < 8; }},
        {"peak window, OA > 8", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.oa >= 8; }},
        {"peak window, p_ee < 300", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.pee < 300; }},
        {"peak window, p_ee 300-500", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.pee >= 300 && r.pee < 500; }},
        {"peak window, p_ee > 500", [](const CR& r) { return r.m4 >= kPkLo && r.m4 < kPkHi && r.pee >= 500; }},
    };
    P("%-28s %9s %9s %16s %16s %8s %8s %8s %6s %6s %6s %7s %7s | %7s %6s %7s %6s %5s %9s\n", "selection", "sig/BR nb", "N_eff",
      "mean [MeV]", "sigma_core", "q16", "q50", "q84", "f1", "f2", "f3", "<-3s", ">+3s", "m1", "s1", "m2", "s2", "fr2", "chi2/ndf");
    for (size_t k = 0; k < sls.size(); ++k) {
        std::vector<std::pair<double, double>> xw;
        TH1D h(Nm("hls").c_str(), "", 240, 0, 60);
        h.Sumw2();
        for (const CR& r : gX17all) if (sls[k].f(r)) { xw.push_back(std::make_pair((double)r.mee, r.w)); h.Fill(r.mee, r.w); }
        if (xw.size() < 50) continue;
        const LineFit L = FitLine(xw, &h);
        if (k == 0) gLineAll = L;
        if (k == 1) gLinePk = L;
        P("%-28s %9.5f %9.0f %7.3f+-%.3f %7.3f+-%.3f %8.3f %8.3f %8.3f %6.3f %6.3f %6.3f %7.4f %7.4f | %7.3f %6.3f %7.3f %6.3f %5.3f %5.0f/%d\n",
          sls[k].name, L.tot / kLsim, L.neff, L.mean, L.meanE, L.sig, L.sigE, L.q16, L.q50, L.q84, L.f1, L.f2, L.f3, L.flo3, L.fhi3,
          L.m1, L.s1, L.m2, L.s2, L.fr2, L.chi2, L.ndf);
    }
    P("\n# sig/BR = X17 cross section per unit BR [nb] (sum W / L_sim with W = 2): multiply by BR for the cross section.\n");
    // M4 of X17
    {
        std::vector<std::pair<double, double>> xw;
        TH1D h(Nm("hm4").c_str(), "", 300, 0.40, 0.70);
        for (const CR& r : gX17all) if (r.m4 >= 0.40 && r.m4 < 0.70) { xw.push_back(std::make_pair((double)r.m4 * 1000., r.w)); h.Fill(r.m4 * 1000., r.w); }
        TH1D h2(Nm("hm4b").c_str(), "", 300, 400, 700);
        for (auto& x : xw) h2.Fill(x.first, x.second);
        const LineFit L = FitLine(xw, &h2);
        double tot = 0, inPk = 0, inBu = 0;
        for (const CR& r : gX17all) { tot += r.w; if (r.m4 >= kPkLo && r.m4 < kPkHi) inPk += r.w; if (r.m4 >= kBuLo && r.m4 < kBuHi) inBu += r.w; }
        P("\nX17 reconstructed M4 (pre, 0.40-0.70): core mean %.2f +- %.2f MeV/c2, sigma_core %.2f +- %.2f MeV/c2; q16/q50/q84 %.1f/%.1f/%.1f;"
          " fraction of all-M4 X17 in the peak window %.4f, in the bump window %.4f\n", L.mean, L.meanE, L.sig, L.sigE, L.q16, L.q50, L.q84,
          inPk / tot, inBu / tot);
    }
    CloseTab();
}

// ============================================================================ part 2: backgrounds near 17 MeV
void Part2() {
    OpenTab("s5_opt_x17_bkg.txt");
    P("# s5_opt_x17_bkg: background groups near M(ee) = 17 MeV/c2 in the M4 peak window %.4f-%.4f (pre)\n", kPkLo, kPkHi);
    P("# std + lep + eta owned rows, mode A, w (PDG) x physics weights V3 (kPhysMain); sigma = sum w / 221.05 nb^-1; errors from event groups\n");
    P("# X17 columns: per unit BR (W = 2); multiply by BR. 'X17 other rows' = eta-set rows with an X17 track that are not class 6\n");
    P("#   (W = 0 in production; shown with W = 2 BR, natural factors neglected; NOT part of the agreed signal)\n\n");
    struct W { const char* name; double lo, hi; };
    std::vector<W> ws = {{"M(ee) all", -1, 1e9}, {"M(ee) 0-10", 0, 10}, {"M(ee) 10-12", 10, 12}, {"M(ee) 12-22", 12, 22}, {"M(ee) 14-20", 14, 20},
                         {"M(ee) 22-30", 22, 30}, {"M(ee) 30-60", 30, 60}};
    for (int pass = 0; pass < 2; ++pass) {
        Cuts c; if (pass == 1) c = gOpt0;
        P("## %s\n", pass == 0 ? "selection pre + M4 peak window" : Form("selection optimum without its M(ee) window: %s", CutStr(c).c_str()));
        P("%-26s", "group [nb] (N_eff)");
        for (const W& w : ws) P(" %24s", w.name);
        P("\n");
        std::vector<Ev> evs;
        for (const W& w : ws) {
            Cuts t = c; t.meelo = w.lo < 0 ? kNo : w.lo; t.meehi = w.hi > 1e8 ? kNo : w.hi;
            evs.push_back(Eval(gAll, t));
        }
        for (int k = 0; k < kNCat; ++k) {
            P("%-26s", kCatName[k]);
            for (const Ev& e : evs) {
                const double v = e.b[k] / kLsim, er = std::sqrt(e.vb[k]) / kLsim, ne = e.vb[k] > 0 ? e.b[k] * e.b[k] / e.vb[k] : 0;
                P(" %9.4g+-%-7.2g(%5.0f)", v, er, ne);
            }
            P("\n");
        }
        P("%-26s", "T1(no X17)+T2 sum");
        for (const Ev& e : evs) P(" %9.4g+-%-7.2g(%5.0f)", e.B / kLsim, std::sqrt(e.VB) / kLsim, e.NeffB());
        P("\n%-26s", "+ T3 + T4 (variant)");
        for (const Ev& e : evs) P(" %9.4g+-%-7.2g(%5.0f)", e.B34 / kLsim, std::sqrt(e.VB34) / kLsim, e.NeffB34());
        P("\n%-26s", "S/B at BR = 1e-5");
        for (const Ev& e : evs) P(" %24.4f", e.B > 0 ? e.s * kBRill / e.B : 0.);
        P("\n%-26s", "Gsig5 share of B");
        for (const Ev& e : evs) P(" %24.3f", e.B > 0 ? e.b[cS5] / e.B : 0.);
        P("\n\n");
    }
    // OA distribution summaries (quantiles) near 17 MeV, pre + peak window + M(ee) 12-22
    P("## OA quantiles (deg) of each group in the peak window with M(ee) 12-22 MeV/c2 (pre): q10 / q50 / q90 and the fraction with OA > 2, > 4, > 8\n");
    for (int k = 0; k < kNCat; ++k) {
        std::vector<std::pair<double, double>> xw;
        for (const CR& r : gAll) if (r.cat == k && r.m4 >= kPkLo && r.m4 < kPkHi && r.mee >= 12 && r.mee < 22) xw.push_back(std::make_pair((double)r.oa, r.w));
        if (xw.size() < 3) { P("%-26s (fewer than 3 rows)\n", kCatName[k]); continue; }
        double t = 0, a2 = 0, a4 = 0, a8 = 0;
        for (auto& x : xw) { t += x.second; if (x.first > 2) a2 += x.second; if (x.first > 4) a4 += x.second; if (x.first > 8) a8 += x.second; }
        const double q10 = WQuant(xw, 0.1), q50 = WQuant(xw, 0.5), q90 = WQuant(xw, 0.9);
        P("%-26s rows %7zu  q10 %6.2f  q50 %6.2f  q90 %6.2f   >2 %.3f  >4 %.3f  >8 %.3f\n", kCatName[k], xw.size(), q10, q50, q90, a2 / t, a4 / t, a8 / t);
    }
    CloseTab();
}

// ============================================================================ part 4: optimisation and split half
std::vector<double> gScanFull[32], gScanA[32], gScanB[32];   // 1D scans at the optimum (descent vars then extras)
void Scan1D(const Cuts& c0) {
    int iv = 0;
    for (int set = 0; set < 2; ++set) {
        const std::vector<OVar>& vs = set == 0 ? gVars : gExtra;
        for (const OVar& v : vs) {
            gScanFull[iv].clear(); gScanA[iv].clear(); gScanB[iv].clear();
            for (double x : v.vals) {
                Cuts t = c0; t.*(v.field) = x;
                bool bad = t.m4hi <= t.m4lo || (On(t.meelo) && On(t.meehi) && t.meehi <= t.meelo) || (On(t.oalo) && On(t.oahi) && t.oahi <= t.oalo) ||
                           (On(t.mpplo) && On(t.mpphi) && t.mpphi <= t.mpplo);
                gScanFull[iv].push_back(bad ? NAN : FomCuts(gScan, t, 0));
                gScanA[iv].push_back(bad ? NAN : FomCuts(gScan, t, 1));
                gScanB[iv].push_back(bad ? NAN : FomCuts(gScan, t, 2));
            }
            ++iv;
        }
    }
}
void PrintCuts(const char* name, const Cuts& c) {
    P("%-14s m4lo = %.4f, m4hi = %.4f, meelo = %s, meehi = %s, oalo = %s, oahi = %s, labhi = %s, mpplo = %s, mpphi = %s, etalo = %s, vtx = %s\n",
      name, c.m4lo, c.m4hi, On(c.meelo) ? Form("%.1f", c.meelo) : "none", On(c.meehi) ? Form("%.1f", c.meehi) : "none",
      On(c.oalo) ? Form("%.1f", c.oalo) : "none", On(c.oahi) ? Form("%.0f", c.oahi) : "none", On(c.labhi) ? Form("%.0f", c.labhi) : "none",
      On(c.mpplo) ? Form("%.2f", c.mpplo) : "none", On(c.mpphi) ? Form("%.2f", c.mpphi) : "none", On(c.etalo) ? Form("%.0f", c.etalo) : "none",
      On(c.vtx) ? "yes" : "no");
}
void Part4() {
    InitVars(); InitExtra();
    Cuts start; start.meelo = 12; start.meehi = 22;
    TStopwatch sw;
    gOpt0 = Descent(gScan, start, 0, &gLog);
    gOptA = Descent(gScan, start, 1, &gLogA);
    gOptB = Descent(gScan, start, 2, &gLogB);
    printf("descents done in %.0f s\n", sw.RealTime());
    // second start (looser M(ee) window, slides cuts) to test the uniqueness of the optimum
    Cuts start2 = SlidesCuts(); start2.meelo = 10; start2.meehi = 26;
    std::vector<StepLog> lg2;
    const Cuts opt2 = Descent(gScan, start2, 0, &lg2);
    Scan1D(gOpt0);
    OpenTab("s5_opt_x17_scan.txt");
    P("# s5_opt_x17_scan: optimisation of the X17 selection (figure of merit = expected 90%% CL UL on BR(eta->pipiX17) x BR(X17->ee))\n");
    P("# FOM (v2) = sqrt(UL_stat^2 + (1.645 sigma_B,MC / s1)^2), UL_stat = s_up(B) / s1: s_up = CLs upper limit for n = B (Asimov,\n");
    P("#   continuous Poisson), sigma_B,MC = simulation (event-group) error of B, s1 = selected X17 rows per unit BR (FOM mode %d);\n", gFomMC);
    P("#   the M4 window is fixed to %.4f-%.4f (study peak window); its edges are scanned at the optimum only\n", kPkLo, kPkHi);
    P("#   rows at L_data = 6080.5 nb^-1 (UNCONFIRMED); B = T1 (without X17) + T2 with physics weights V3; pre always applied\n");
    P("# descent: coordinate descent over the variables below, plateau rule = loosest value within %.0f %% of the step minimum\n", kTol * 100);
    P("# split half: half A = event groups with even grp, half B = odd grp (each e+e- tree); half results use weights x 2\n");
    P("# start: pre + M4 %.4f-%.4f + M(ee) 12-22 MeV; second start: slides cuts + M(ee) 10-26 MeV\n\n", kPkLo, kPkHi);
    auto plog = [&](const char* nm, const std::vector<StepLog>& L) {
        P("## path of the descent (%s): pass, variable, from -> to, FOM after the step\n", nm);
        for (const StepLog& s : L) P("  %d  %-6s %9s -> %-9s  %.4g\n", s.pass, s.var.c_str(), On(s.from) ? Form("%.4g", s.from) : "none",
                                      On(s.to) ? Form("%.4g", s.to) : "none", s.fom);
    };
    plog("full sample", gLog); plog("half A", gLogA); plog("half B", gLogB); plog("full sample, second start", lg2);
    P("\n## optima\n");
    PrintCuts("full", gOpt0); PrintCuts("half A", gOptA); PrintCuts("half B", gOptB); PrintCuts("full, start 2", opt2);
    P("\n## figure of merit (UL on BR, x 1e-6) of each optimum evaluated on the full sample and on each half\n");
    P("%-16s %10s %10s %10s\n", "optimum", "full", "half A", "half B");
    for (int k = 0; k < 4; ++k) {
        const Cuts& c = k == 0 ? gOpt0 : k == 1 ? gOptA : k == 2 ? gOptB : opt2;
        P("%-16s %10.4f %10.4f %10.4f\n", k == 0 ? "full" : k == 1 ? "half A" : k == 2 ? "half B" : "full, start 2", FomCuts(gScan, c, 0) * 1e6,
          FomCuts(gScan, c, 1) * 1e6, FomCuts(gScan, c, 2) * 1e6);
    }
    {
        const double aa = FomCuts(gScan, gOptA, 1), ab = FomCuts(gScan, gOptA, 2), bb = FomCuts(gScan, gOptB, 2), ba = FomCuts(gScan, gOptB, 1);
        P("\nsplit-half bias: optimum A on A %.4f, on B %.4f (ratio %.3f); optimum B on B %.4f, on A %.4f (ratio %.3f) [x 1e-6]\n", aa * 1e6, ab * 1e6,
          ab / aa, bb * 1e6, ba * 1e6, ba / bb);
        P("the full-sample optimum on the halves: A %.4f, B %.4f [x 1e-6]; spread of the half results measures the MC statistical noise of the FOM\n",
          FomCuts(gScan, gOpt0, 1) * 1e6, FomCuts(gScan, gOpt0, 2) * 1e6);
    }
    // 1D scans at the optimum
    P("\n## one-dimensional scans at the full-sample optimum (UL x 1e-6; full / half A / half B); '*' = chosen value\n");
    int iv = 0;
    for (int set = 0; set < 2; ++set) {
        const std::vector<OVar>& vs = set == 0 ? gVars : gExtra;
        for (const OVar& v : vs) {
            P("# %s%s\n", v.key, set ? " (extra variable, not in the descent)" : "");
            for (size_t i = 0; i < v.vals.size(); ++i) {
                const bool ch = std::fabs(gOpt0.*(v.field) - v.vals[i]) < 1e-9;
                P("   %8s %c  %9.4f %9.4f %9.4f\n", On(v.vals[i]) ? Form("%.4g", v.vals[i]) : "none", ch ? '*' : ' ', gScanFull[iv][i] * 1e6,
                  gScanA[iv][i] * 1e6, gScanB[iv][i] * 1e6);
            }
            ++iv;
        }
    }
    CloseTab();
}

// ============================================================================ part 8: final numbers, shape analysis, frozen file
struct Final {
    Ev e; Fom f, f34;
    double ulAsy = 0, ulMedCls = 0, ulLoCls = 0, ulHiCls = 0, ulFC = 0, z = 0, z34 = 0;
    int nmed = 0, n16 = 0, n84 = 0;
    double ulBp = 0, ulBm = 0;     // B +- MC error
};
Final Summarise(const Cuts& c, int half = 0) {
    Final F;
    F.e = Eval(gAll, c, half);
    F.f = FomOf(F.e); F.f34 = FomOf(F.e, true);
    const double B = F.f.B, s1 = F.f.s1;
    if (s1 <= 0) return F;
    F.ulAsy = AsyUL(B) / s1;
    F.nmed = PQuant(B, 0.5); F.n16 = PQuant(B, 0.1587); F.n84 = PQuant(B, 0.8413);
    F.ulMedCls = ClsUL(F.nmed, B) / s1; F.ulLoCls = ClsUL(F.n16, B) / s1; F.ulHiCls = ClsUL(F.n84, B) / s1;
    F.ulFC = FCUpper(F.nmed, B) / s1;
    F.z = ZA(kBRill * s1, B); F.z34 = ZA(kBRill * s1, F.f34.B);
    const double sB = std::sqrt(F.e.VB) * kToRows;
    F.ulBp = ClsUL(B + sB, B + sB) / s1; F.ulBm = ClsUL(std::max(0., B - sB), std::max(0., B - sB)) / s1;
    return F;
}
struct ShapeRes { double ul[3] = {0, 0, 0}, ulT[3] = {0, 0, 0}, z[3] = {0, 0, 0}, zT[3] = {0, 0, 0}; double sTot = 0, mean = 0, sig = 0; };
ShapeRes gShOpt, gShSl, gShOptPk;
ShapeIn gInOpt, gInSl;
std::vector<double> gTplOpt, gTplSl;
std::vector<std::vector<double>> gGrpOpt, gGrpSl;
const double kShLo = 2, kShHi = 62; const int kShN = 60;
ShapeRes DoShape(const Cuts& c, ShapeIn& in, std::vector<double>& tpl, std::vector<std::vector<double>>& grp, const char* name) {
    ShapeRes R;
    ShapeHists(c, kShLo, kShHi, kShN, in, &tpl, &grp);
    // line shape at this selection (X17 rows of gAll passing the cuts without the M(ee) window)
    std::vector<std::pair<double, double>> xw;
    TH1D h(Nm("hsh").c_str(), "", 240, 0, 60);
    for (const CR& r : gAll) if (r.cat == cX17 && PassC(r, c, false, true)) { xw.push_back(std::make_pair((double)r.mee, r.w)); h.Fill(r.mee, r.w); }
    const LineFit L = FitLine(xw, &h);
    double tot = 0; for (double x : tpl) tot += x;
    R.sTot = tot; R.mean = L.mean; R.sig = L.sig;
    for (int mode = 0; mode < 3; ++mode) {
        FillGauss(in, tot, L.mean, L.sig);
        R.ul[mode] = ShapeUL(in, mode); R.z[mode] = ShapeZ(in, kBRill, mode);
        in.s = tpl;
        R.ulT[mode] = ShapeUL(in, mode); R.zT[mode] = ShapeZ(in, kBRill, mode);
    }
    FillGauss(in, tot, L.mean, L.sig);   // leave the Gaussian in place for the figures
    P("%-34s X17 in %g-%g MeV: %.2f rows per unit BR (x 1e-5: %.2f); line mean %.2f, sigma %.2f MeV\n", name, kShLo, kShHi, tot, tot * 1e-5, L.mean, L.sig);
    const char* mn[3] = {"B known", "B norm. profiled", "Gsig5 and rest profiled"};
    for (int mode = 0; mode < 3; ++mode)
        P("    %-26s UL(Gauss) %.4f  UL(template) %.4f [x 1e-6];  Z(1e-5) Gauss %.2f  template %.2f\n", mn[mode], R.ul[mode] * 1e6, R.ulT[mode] * 1e6,
          R.z[mode], R.zT[mode]);
    return R;
}
void PrintFinal(const char* name, const Cuts& c, const Final& F) {
    P("## %s: %s\n", name, CutStr(c).c_str());
    P("X17 per unit BR: sigma %.5f +- %.5f nb (N_eff %.0f, rows %.0f) = %.1f expected rows per unit BR; at BR = 1e-5: %.2f rows\n", F.e.s / kLsim,
      std::sqrt(F.e.vs) / kLsim, F.e.vs > 0 ? F.e.s * F.e.s / F.e.vs : 0, F.e.ns, F.f.s1, F.f.s1 * 1e-5);
    P("%-26s %12s %10s %8s %9s %12s\n", "background group", "sigma [nb]", "error", "N_eff", "rows", "data rows");
    for (int k = 1; k < kNCat; ++k) {
        const double ne = F.e.vb[k] > 0 ? F.e.b[k] * F.e.b[k] / F.e.vb[k] : 0;
        P("%-26s %12.5g %10.3g %8.1f %9.0f %12.2f%s\n", kCatName[k], F.e.b[k] / kLsim, std::sqrt(F.e.vb[k]) / kLsim, ne, F.e.nb[k], F.e.b[k] * kToRows,
          k == cX17o ? "  (per unit BR, x BR; not in B)" : (k >= cFake ? "  (variant only)" : ""));
    }
    P("%-26s %12.5g %10.3g %8.1f %9s %12.2f +- %.2f\n", "B = T1(no X17) + T2", F.e.B / kLsim, std::sqrt(F.e.VB) / kLsim, F.e.NeffB(), "", F.f.B,
      std::sqrt(F.e.VB) * kToRows);
    P("%-26s %12.5g %10.3g %8.1f %9s %12.2f +- %.2f\n", "B + T3 + T4 (variant)", F.e.B34 / kLsim, std::sqrt(F.e.VB34) / kLsim, F.e.NeffB34(), "",
      F.f34.B, std::sqrt(F.e.VB34) * kToRows);
    P("expected 90%% CL UL on BR, B known exactly (UL_stat) [x 1e-6]: Asimov CLs %.4f (1.64 sqrt(B) approx %.4f; asymptotic %.4f); median CLs (n = %d) %.4f, "
      "band n = %d..%d: %.4f..%.4f; FC (n = %d) %.4f\n", F.f.ulStat * 1e6, 1.6449 * std::sqrt(F.f.B) / F.f.s1 * 1e6, F.ulAsy * 1e6, F.nmed,
      F.ulMedCls * 1e6, F.n16, F.n84, F.ulLoCls * 1e6, F.ulHiCls * 1e6, F.nmed, F.ulFC * 1e6);
    P("  MC statistics of B: UL with B +- 1 sigma_MC: %.4f / %.4f; variant B + T3 + T4: UL %.4f, Z(1e-5) %.2f\n", F.ulBp * 1e6, F.ulBm * 1e6,
      F.f34.ulStat * 1e6, F.z34);
    P("  X17 other rows (not class 6, W = 2 BR approx.): %.2f rows per unit BR = %.3f of the class-6 yield (would lower the UL by about this fraction)\n",
      F.e.b[cX17o] * kToRows, F.e.s > 0 ? F.e.b[cX17o] / F.e.s : 0);
    P("discovery significance for BR = 1e-5 (Asimov, B known): Z = %.2f (s = %.1f, B = %.1f, s/B = %.3f)\n", F.z, F.f.s1 * 1e-5, F.f.B,
      F.f.B > 0 ? F.f.s1 * 1e-5 / F.f.B : 0);
    P("v2 limits [x 1e-6]: UL_stat (B known exactly) %.4f; MC error of B %.1f rows -> 1.645 sigma_B,MC / s1 = %.4f; UL including it %.4f "
      "(ratio %.2f); Z(1e-5) with the MC error of B: %.2f\n", F.f.ulStat * 1e6, F.f.sBmc, 1.6449 * F.f.sBmc / F.f.s1 * 1e6, F.f.ulTot * 1e6,
      F.f.ulTot / F.f.ulStat, F.f.s1 * 1e-5 / std::sqrt(F.f.B + F.f.sBmc * F.f.sBmc));
    P("share of the MC variance of B by group (diagonal terms; rows of one event in two groups give the rest): ");
    for (int k = cS2; k <= cOth; ++k) if (F.e.vb[k] > 0) P("%s %.3f; ", kCatName[k], F.e.vb[k] / F.e.VB);
    P("\n\n");
}
Final gFinOpt, gFinSl, gFinSlW, gFinPre, gFinV1;
Cuts gSlW;
void Part8() {
    OpenTab("s5_opt_x17_result.txt");
    P("# s5_opt_x17_result: X17 selection, final numbers (simulation only; L_data = 6080.5 nb^-1 UNCONFIRMED for all rows and limits)\n");
    P("# std + lep + eta owned rows, mode A, w (PDG) x physics weights V3; signal = class 6 with W = 2 BR; B = T1 (no X17) + T2\n");
    P("# errors: event-group variance; N_eff = S^2 / V. Limits: 90 %% CL, B known exactly (no systematic uncertainty)\n\n");
    // engine cross-check of the compact evaluation (optimum)
    {
        std::vector<Comp> comps = {CompX17(1.0), CompSim(), CompFake(), CompDiff()};
        const std::vector<WSum> ws = SumMany(*gD, SelOf(gOpt0), comps, gOpt);
        const Ev e = Eval(gAll, gOpt0);
        P("engine cross-check (optimum): X17 per unit BR %.6f vs %.6f; T1+T2 %.6f vs %.6f; T3 %.6f vs %.6f; T4 %.6f vs %.6f (sum w)\n",
          ws[0].S, e.s, ws[1].S, e.B, ws[2].S, e.b[cFake], ws[3].S, e.b[cDiff]);
        P("  group variance: X17 %.6g vs %.6g; T1+T2 %.6g vs %.6g\n", ws[0].V, e.vs, ws[1].V, e.VB);
        const std::vector<WSum> wp = SumMany(*gD, SelOf(Cuts(), false, true), {CompX17(1e-5)}, gOpt);
        P("  reference: X17 at BR 1e-5 in pre + peak window %.4f-%.4f: %.5f nb (eta_engine_reference uses 0.5336-0.5542: 0.05575)\n\n", kPkLo, kPkHi, wp[0].S / kLsim);
    }
    gFinPre = Summarise([] { Cuts c; c.meelo = 12; c.meehi = 22; return c; }());
    gFinOpt = Summarise(gOpt0);
    gFinSl = Summarise(SlidesCuts());
    gSlW = SlidesCuts(); gSlW.meelo = gOpt0.meelo; gSlW.meehi = gOpt0.meehi;
    gFinSlW = Summarise(gSlW);
    PrintFinal("pre + peak window + M(ee) 12-22 (start)", [] { Cuts c; c.meelo = 12; c.meehi = 22; return c; }(), gFinPre);
    PrintFinal("OPTIMUM (frozen, v2)", gOpt0, gFinOpt);
    {
        // v1 selection (frozen 2026-11-12, used by the blind test): M4 0.5330-0.5555, M(ee) 13.0-26.5, pions in the vertex fit
        Cuts v1; v1.m4lo = 0.5330; v1.m4hi = 0.5555; v1.meelo = 13.0; v1.meehi = 26.5; v1.vtx = 1.;
        gFinV1 = Summarise(v1);
        PrintFinal("v1 selection (frozen 2026-11-12, used by the blind test)", v1, gFinV1);
        Cuts v1p = v1; v1p.m4lo = kPkLo; v1p.m4hi = kPkHi;
        const Final F = Summarise(v1p);
        P("v1 cuts in the fixed peak window %.4f-%.4f: UL_stat %.4f, UL with MC error %.4f [x 1e-6]\n\n", kPkLo, kPkHi, F.f.ulStat * 1e6, F.f.ulTot * 1e6);
    }
    PrintFinal("slides + peak window (no M(ee) window)", SlidesCuts(), gFinSl);
    PrintFinal("slides + peak window + optimum M(ee) window", gSlW, gFinSlW);
    // optimum on the halves
    for (int h = 1; h <= 2; ++h) {
        const Final F = Summarise(gOpt0, h);
        P("optimum evaluated on half %c (weights x 2): UL %.4f e-6 (UL_stat %.4f), B %.1f +- %.1f rows, s1 %.1f, Z(1e-5) %.2f\n", h == 1 ? 'A' : 'B',
          F.f.ul * 1e6, F.f.ulStat * 1e6, F.f.B, std::sqrt(F.e.VB) * kToRows, F.f.s1, F.z);
    }
    for (int k = 0; k < 2; ++k) {
        const Cuts& c = k == 0 ? gOptA : gOptB;
        const Final Fs = Summarise(c, k + 1), Fo = Summarise(c, 2 - k);
        P("optimum of half %c: on its own half UL %.4f e-6 (B %.1f), on the other half UL %.4f e-6 (B %.1f), full %.4f e-6\n", k == 0 ? 'A' : 'B',
          Fs.f.ul * 1e6, Fs.f.B, Fo.f.ul * 1e6, Fo.f.B, Summarise(c).f.ul * 1e6);
    }
    // cut flow
    P("\n## cut flow towards the optimum (UL x 1e-6, s1 = X17 rows per unit BR, B rows; L_data UNCONFIRMED)\n");
    {
        std::vector<std::pair<std::string, Cuts>> fl;
        Cuts c; c.m4lo = gOpt0.m4lo; c.m4hi = gOpt0.m4hi;
        fl.push_back(std::make_pair(std::string("pre + M4 window"), c));
        c.meelo = gOpt0.meelo; c.meehi = gOpt0.meehi; fl.push_back(std::make_pair(std::string("+ M(ee) window"), c));
        c.oalo = gOpt0.oalo; c.oahi = gOpt0.oahi; fl.push_back(std::make_pair(std::string("+ OA cuts"), c));
        c.labhi = gOpt0.labhi; fl.push_back(std::make_pair(std::string("+ alpha_lab"), c));
        c.mpplo = gOpt0.mpplo; c.mpphi = gOpt0.mpphi; fl.push_back(std::make_pair(std::string("+ M(pipi)"), c));
        c.etalo = gOpt0.etalo; fl.push_back(std::make_pair(std::string("+ alpha_eta"), c));
        c.vtx = gOpt0.vtx; fl.push_back(std::make_pair(std::string("+ pions in vertex"), c));
        P("%-20s %10s %10s %10s %9s %9s %8s %7s %7s %7s\n", "step", "UL", "UL_stat", "s1", "B", "dB(MC)", "N_eff B", "eps_S", "eps_B", "Z(1e-5)");
        const Ev e0 = Eval(gAll, fl[0].second);
        for (auto& x : fl) {
            const Ev e = Eval(gAll, x.second);
            const Fom f = FomOf(e);
            P("%-20s %10.4f %10.4f %10.1f %9.1f %9.1f %8.0f %7.3f %7.4f %7.2f\n", x.first.c_str(), f.ul * 1e6, f.ulStat * 1e6, f.s1, f.B, std::sqrt(e.VB) * kToRows, e.NeffB(),
              e.s / e0.s, e.B / e0.B, ZA(f.s1 * 1e-5, f.B));
        }
    }
    // extra variables at the optimum
    P("\n## extra variables at the optimum (not adopted unless stated): best value and UL\n");
    {
        int iv = (int)gVars.size();
        for (const OVar& v : gExtra) {
            double best = 1e30; double bx = kNo;
            for (size_t i = 0; i < v.vals.size(); ++i) if (gScanFull[iv][i] < best) { best = gScanFull[iv][i]; bx = v.vals[i]; }
            double none = gScanFull[iv][0];
            for (size_t i = 0; i < v.vals.size(); ++i) if (!On(v.vals[i])) none = gScanFull[iv][i];
            P("%-8s no cut %.4f e-6, best %s -> %.4f e-6 (gain %.1f %%)\n", v.key, none * 1e6, On(bx) ? Form("%.4g", bx) : "none", best * 1e6,
              100 * (1 - best / none));
            ++iv;
        }
    }
    // variants of the counting result at the optimum
    P("\n## variants at the optimum (counting, Asimov CLs; UL x 1e-6; L_data UNCONFIRMED)\n");
    {
        const Final& F = gFinOpt;
        const double s1 = F.f.s1, B = F.f.B;
        P("baseline: UL (with MC error of B) %.4f, UL_stat %.4f, B = %.1f +- %.1f (MC) rows -> UL_stat +- %.4f from a 1-sigma shift of B; Z(1e-5) %.2f\n",
          F.f.ul * 1e6, F.f.ulStat * 1e6, B, std::sqrt(F.e.VB) * kToRows, 0.5 * (F.ulBp - F.ulBm) * 1e6, F.z);
        P("# the variants below are UL_stat (B known exactly), as in v1\n");
        gW0 = true; const Ev e0 = Eval(gAll, gOpt0); gW0 = false;
        const Fom f0 = FomOf(e0);
        P("physics weights OFF (same cuts): B %.1f rows (Gsig5 %.1f instead of %.1f; bump/conv. classes unweighted), UL_stat %.4f, UL with MC error %.4f, Z(1e-5) %.2f\n", f0.B,
          e0.b[cS5] * kToRows, F.e.b[cS5] * kToRows, f0.ulStat * 1e6, f0.ulTot * 1e6, ZA(1e-5 * f0.s1, f0.B));
        {
            Cuts start; start.meelo = 12; start.meehi = 22;
            gW0 = true; const Cuts c0 = Descent(gScan, start, 0, nullptr); const double u0 = FomCuts(gScan, c0, 0); gW0 = false;
            P("optimum re-found with the physics weights OFF (FOM mode %d): %s; FOM (weights off) %.4f, the same cuts with weights ON %.4f\n", gFomMC, CutStr(c0).c_str(),
              u0 * 1e6, FomCuts(gScan, c0, 0) * 1e6);
        }
        for (double sc : {0.5, 0.75, 1.25, 1.5, 2.0}) {
            const double Bs = B * sc;
            P("B x %.2f: UL %.4f (x %.3f), Z(1e-5) %.2f\n", sc, ClsUL(Bs, Bs) / s1 * 1e6, ClsUL(Bs, Bs) / s1 / F.f.ulStat, ZA(1e-5 * s1, Bs));
        }
        {
            double bs = 0; for (int k = cS2; k <= cS5; ++k) bs += F.e.b[k];
            for (double sc : {0.7, 1.3}) {
                const double Bs = (bs + sc * (F.e.B - bs)) * kToRows;
                P("template-2 groups (G_eta, G_pi0D, G_conv, G_vec, G_oth) x %.1f, Gsig2-5 unchanged: B %.1f, UL %.4f\n", sc, Bs, ClsUL(Bs, Bs) / s1 * 1e6);
            }
        }
        const double s1o = s1 + F.e.b[cX17o] * kToRows;
        P("signal including the X17 other rows (W = 2 BR, approximation): s1 %.4g instead of %.4g, UL %.4f, Z(1e-5) %.2f\n", s1o, s1,
          ClsUL(B, B) / s1o * 1e6, ZA(1e-5 * s1o, B));
        P("background composition (share of B): ");
        for (int k = cS2; k <= cOth; ++k) P("%s %.3f; ", kCatName[k], F.e.b[k] / F.e.B);
        P("\n");
    }
    // shape analysis
    P("\n## shape analysis: binned Poisson likelihood of M(ee) (%g-%g MeV, %d bins) in the M4 window, all other cuts applied\n", kShLo, kShHi, kShN);
    P("# Asimov B-only data; median expected CLs limit from q(mu_up) = 2.706 (asymptotic); X17 = Gaussian line (core mean, sigma; carrying the\n");
    P("# whole selected X17 yield in the fit range) or the simulated X17 template; background shape = simulation (Gsig5 and the rest separately\n");
    P("# in the two-factor mode). Counting at the same selection with the optimum M(ee) window is given above.\n");
    gShOpt = DoShape([] { Cuts c = gOpt0; c.meelo = kNo; c.meehi = kNo; return c; }(), gInOpt, gTplOpt, gGrpOpt, "optimum (no M(ee) window)");
    {
        ShapeIn tmp; std::vector<double> t; std::vector<std::vector<double>> g;
        Cuts c = gOpt0; c.meelo = kNo; c.meehi = kNo; c.m4lo = kPkLo; c.m4hi = kPkHi;
        gShOptPk = DoShape(c, tmp, t, g, "optimum cuts, M4 peak window");
    }
    gShSl = DoShape(SlidesCuts(), gInSl, gTplSl, gGrpSl, "slides (peak window)");
    {
        P("\n# robustness of the shape result at the optimum cuts (Gaussian line, mean %.2f, sigma %.2f MeV): binning and fit range\n", gShOpt.mean, gShOpt.sig);
        P("%-22s %10s %12s %12s %12s %8s\n", "range, bin", "X17 rows", "B known", "B profiled", "2 factors", "Z(1e-5)");
        Cuts c = gOpt0; c.meelo = kNo; c.meehi = kNo;
        struct BV { double lo, hi, bw; };
        for (const BV& b : std::vector<BV>{{2, 62, 0.5}, {2, 62, 1}, {2, 62, 2}, {2, 62, 4}, {8, 40, 1}, {8, 40, 2}, {10, 30, 1}}) {
            ShapeIn in; std::vector<double> tpl;
            ShapeHists(c, b.lo, b.hi, (int)std::lround((b.hi - b.lo) / b.bw), in, &tpl, nullptr);
            double tot = 0; for (double x : tpl) tot += x;
            FillGauss(in, tot, gShOpt.mean, gShOpt.sig);
            P("%5.0f-%-4.0f %4.1f MeV    %10.4g %12.4f %12.4f %12.4f %8.2f\n", b.lo, b.hi, b.bw, tot * 1e-5, ShapeUL(in, 0) * 1e6, ShapeUL(in, 1) * 1e6,
              ShapeUL(in, 2) * 1e6, ShapeZ(in, kBRill, 1));
        }
        P("# X17 rows column: expected X17 rows at BR = 1e-5 inside the fit range; Z with the background normalisation profiled\n");
    }
    {
        // counting cross-check through the shape machinery: one bin = the optimum M(ee) window
        ShapeIn one; one.nb = 1; one.b0 = {gFinOpt.e.b[cS5] * kToRows}; one.b1 = {(gFinOpt.e.B - gFinOpt.e.b[cS5]) * kToRows}; one.s = {gFinOpt.f.s1};
        P("check: one-bin likelihood at the optimum: UL %.4f e-6 (asymptotic counting %.4f e-6)\n", ShapeUL(one, 0) * 1e6, gFinOpt.ulAsy * 1e6);
    }
    CloseTab();
    // frozen file
    OpenTab("s5_opt_x17_frozen.txt");
    const Cuts& c = gOpt0;
    P("# s5_opt_x17_frozen: chosen X17 selection, version 2 (frozen %s), from s5_opt_x17.C\n", "2026-11-13");
    P("# v2 (review findings H2, M2): the figure of merit includes the simulation (MC) error of the background, and the M4 window is fixed to the\n");
    P("#   study peak window 0.5333-0.5587. The blind test of 2026-11-13 used v1 (tables/s5_opt_x17_frozen_v1.txt: M4 0.5330-0.5555, M(ee)\n");
    P("#   13.0-26.5 MeV/c2, pions in the vertex fit; UL_stat 5.48e-6). UL_stat assumes the background known exactly; it is a statistical limit only.\n");
    P("# simulation only; std + lep + eta owned rows, mode A, w (PDG) x physics weights V3; signal = class 6 (W = 2 BR); B = T1 (no X17) + T2\n");
    P("# all cut variables from corrected momenta; 'none' = no cut; M in GeV/c2 except M(ee) in MeV/c2; angles in deg\n");
    P("base = pre (isBest >= 0, vertex z > -500 mm, PT3 bit, vertex z > -160 mm, cor_ok == 1)\n");
    P("m4_lo = %.4f\nm4_hi = %.4f\n", c.m4lo, c.m4hi);
    P("mee_lo_MeV = %s\nmee_hi_MeV = %s\n", On(c.meelo) ? Form("%.1f", c.meelo) : "none", On(c.meehi) ? Form("%.1f", c.meehi) : "none");
    P("oa_lo_deg = %s\noa_hi_deg = %s\n", On(c.oalo) ? Form("%.1f", c.oalo) : "none", On(c.oahi) ? Form("%.0f", c.oahi) : "none");
    P("alpha_lab_hi_deg = %s\n", On(c.labhi) ? Form("%.0f", c.labhi) : "none");
    P("mpipi_lo = %s\nmpipi_hi = %s\n", On(c.mpplo) ? Form("%.2f", c.mpplo) : "none", On(c.mpphi) ? Form("%.2f", c.mpphi) : "none");
    P("alpha_eta_lo_deg = %s\n", On(c.etalo) ? Form("%.0f", c.etalo) : "none");
    P("pions_in_vertex_fit = %s\n", On(c.vtx) ? "yes" : "no");
    P("dz_pi_hi_mm = none\nmm4_hi = none\npt4_lo = none\n");
    P("x17_sigma_per_unit_BR_nb = %.5f +- %.5f\n", gFinOpt.e.s / kLsim, std::sqrt(gFinOpt.e.vs) / kLsim);
    P("x17_rows_per_unit_BR = %.1f   # L_data = 6080.5 nb^-1 (UNCONFIRMED)\n", gFinOpt.f.s1);
    P("background_sigma_nb = %.5f +- %.5f   # N_eff %.0f\n", gFinOpt.e.B / kLsim, std::sqrt(gFinOpt.e.VB) / kLsim, gFinOpt.e.NeffB());
    P("background_rows = %.1f +- %.1f   # L_data = 6080.5 nb^-1 (UNCONFIRMED)\n", gFinOpt.f.B, std::sqrt(gFinOpt.e.VB) * kToRows);
    P("expected_UL90_BR_with_MC_error_of_B = %.3e   # PRIMARY (v2): counting, sqrt(UL_stat^2 + (1.645 sigma_B,MC / s1)^2), sigma_B,MC = %.1f rows\n",
      gFinOpt.f.ulTot, gFinOpt.f.sBmc);
    P("expected_UL90_BR_asimov_CLs = %.3e +- %.1e   # UL_stat: counting, B known exactly; +- = shift for a 1-sigma change of B\n", gFinOpt.f.ulStat,
      0.5 * (gFinOpt.ulBp - gFinOpt.ulBm));
    P("expected_UL90_BR_median_CLs = %.3e   # band %.3e .. %.3e\n", gFinOpt.ulMedCls, gFinOpt.ulLoCls, gFinOpt.ulHiCls);
    P("expected_UL90_BR_FC = %.3e\n", gFinOpt.ulFC);
    P("expected_UL90_BR_shape_Bknown = %.3e   # binned M(ee) likelihood, Gaussian line\n", gShOpt.ul[0]);
    P("expected_UL90_BR_shape_Bprofiled = %.3e\n", gShOpt.ul[1]);
    P("Z_BR_1e-5 = %.2f   # Asimov counting, B known\n", gFinOpt.z);
    P("expected_UL90_BR_slides_counting = %.3e   # slides + peak window + the M(ee) window above, with the MC error of B (UL_stat %.3e)\n",
      gFinSlW.f.ulTot, gFinSlW.f.ulStat);
    P("expected_UL90_BR_v1_selection_with_MC_error_of_B = %.3e   # the v1 selection evaluated with the v2 figure of merit (UL_stat %.3e)\n",
      gFinV1.f.ulTot, gFinV1.f.ulStat);
    CloseTab();
}

// ============================================================================ part 16: figures
TH1D* HistOf(const std::vector<CR>& rows, std::function<bool(const CR&)> f, std::function<double(const CR&)> x, int cat, int nb, double lo, double hi,
             double scale, const char* nm) {
    // weighted histogram of one category with event-group errors (rows of a group adjacent)
    TH1D* h = new TH1D(Nm(nm).c_str(), "", nb, lo, hi);
    h->SetDirectory(nullptr); h->Sumw2();
    std::vector<double> s(nb + 2, 0.), v(nb + 2, 0.);
    std::vector<std::pair<int, double>> g;
    int cur = -1;
    auto flush = [&]() {
        std::sort(g.begin(), g.end());
        size_t i = 0;
        while (i < g.size()) { const int b = g[i].first; double a = 0; while (i < g.size() && g[i].first == b) { a += g[i].second; ++i; } v[b] += a * a; }
        g.clear();
    };
    const double bw = (hi - lo) / nb;
    for (const CR& r : rows) {
        if (r.gid != cur) { flush(); cur = r.gid; }
        if (cat >= 0 && r.cat != cat) continue;
        if (cat == -2 && !IsBkg(r.cat)) continue;
        if (!f(r)) continue;
        const double xv = x(r);
        if (!(xv >= lo && xv < hi)) continue;
        const int b = std::min(nb - 1, (int)((xv - lo) / bw)) + 1;
        s[b] += r.w; g.push_back(std::make_pair(b, r.w));
    }
    flush();
    for (int b = 1; b <= nb; ++b) { h->SetBinContent(b, s[b] * scale); h->SetBinError(b, std::sqrt(v[b]) * scale); }
    return h;
}
void DrawGroups(TPad* p, std::function<bool(const CR&)> f, std::function<double(const CR&)> x, int nb, double lo, double hi, double scale,
                const char* xt, const char* yt, bool logy, bool legend, const char* title, double xline1 = -1, double xline2 = -1,
                double brShow = kBRill, double ymin = -1) {
    p->cd(); p->SetLogy(logy);
    std::vector<TH1D*> hs(kNCat, nullptr);
    double mx = 0, mn = 1e30;
    for (int k = 0; k < kNCat; ++k) {
        hs[k] = HistOf(gAll, f, x, k, nb, lo, hi, scale * ((k == cX17 || k == cX17o) ? brShow : 1.), "hg");
        for (int b = 1; b <= nb; ++b) if (hs[k]->GetBinContent(b) > 0) { mx = std::max(mx, hs[k]->GetBinContent(b)); mn = std::min(mn, hs[k]->GetBinContent(b)); }
    }
    TH1D* tot = HistOf(gAll, f, x, -2, nb, lo, hi, scale, "htot");
    mx = std::max(mx, tot->GetMaximum());
    TH1D* fr = new TH1D(Nm("fr").c_str(), "", nb, lo, hi);
    fr->SetDirectory(nullptr);
    Frame(fr, xt, yt);
    if (logy) { fr->SetMaximum(mx * 30); fr->SetMinimum(ymin > 0 ? ymin : std::max(mn * 0.5, mx * 1e-5)); }
    else { fr->SetMaximum(mx * 1.35); fr->SetMinimum(0); }
    fr->Draw("AXIS");
    tot->SetLineColor(kBlack); tot->SetLineWidth(3); tot->SetMarkerSize(0);
    tot->Draw("HIST SAME");
    TH1D* tc = (TH1D*)tot->Clone(Nm("tc").c_str()); tc->SetFillColor(kGray); tc->SetFillStyle(3001); tc->SetMarkerSize(0); tc->Draw("E2 SAME");
    for (int k = 1; k < kNCat; ++k) {
        if (k == cX17o) continue;
        hs[k]->SetLineColor(kCatCol[k]); hs[k]->SetLineWidth(2); hs[k]->SetMarkerSize(0);
        if (k >= cFake) hs[k]->SetLineStyle(2);
        if (hs[k]->GetMaximum() > 0) hs[k]->Draw("HIST SAME");
    }
    hs[cX17]->SetLineColor(kCatCol[cX17]); hs[cX17]->SetLineWidth(4); hs[cX17]->SetLineStyle(7); hs[cX17]->Draw("HIST SAME");
    if (xline1 > lo) { TLine* l = new TLine(xline1, fr->GetMinimum(), xline1, logy ? mx * 3 : mx * 1.1); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw(); }
    if (xline2 > lo) { TLine* l = new TLine(xline2, fr->GetMinimum(), xline2, logy ? mx * 3 : mx * 1.1); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw(); }
    Txt(0.17, 0.93, title, 0.05);
    (void)legend;
}
// legend of DrawGroups in its own pad (text 0.045 of the pad height or more)
void GroupLegend(TPad* p, double ts, double brShow = kBRill, const char* extra = nullptr) {
    p->cd();
    TLegend* lg = new TLegend(0.02, 0.06, 0.99, 0.97);
    lg->SetTextSize(ts); lg->SetFillStyle(0); lg->SetTextFont(42);
    TH1D* t = new TH1D(Nm("lt").c_str(), "", 1, 0, 1); t->SetLineColor(kBlack); t->SetLineWidth(3);
    t->SetFillColor(kGray); t->SetFillStyle(3001);
    lg->AddEntry(t, "background sum T1 (no X17) + T2", "lf");
    TH1D* x = new TH1D(Nm("lx").c_str(), "", 1, 0, 1); x->SetLineColor(kCatCol[cX17]); x->SetLineWidth(4); x->SetLineStyle(7);
    lg->AddEntry(x, Form("X17, BR = %.0e (not in the sum)", brShow), "l");
    for (int k = 1; k < kNCat; ++k) {
        if (k == cX17o) continue;
        TH1D* h = new TH1D(Nm("lh").c_str(), "", 1, 0, 1); h->SetLineColor(kCatCol[k]); h->SetLineWidth(3);
        if (k >= cFake) h->SetLineStyle(2);
        lg->AddEntry(h, kCatTex[k], "l");
    }
    if (extra) lg->AddEntry((TObject*)nullptr, extra, "");
    lg->Draw();
}
// width of the X17 line in slices of a variable of the compact row (peak window, pre)
TGraphErrors* SigmaVs(std::function<double(const CR&)> x, const std::vector<double>& e, int col, int sty) {
    TGraphErrors* g = new TGraphErrors();
    for (size_t i = 0; i + 1 < e.size(); ++i) {
        std::vector<std::pair<double, double>> xw;
        TH1D hh(Nm("hh").c_str(), "", 240, 0, 60);
        double sx = 0, sw = 0;
        for (const CR& r : gX17all) {
            if (!(r.m4 >= kPkLo && r.m4 < kPkHi)) continue;
            const double v = x(r);
            if (!(v >= e[i] && v < e[i + 1])) continue;
            xw.push_back(std::make_pair((double)r.mee, r.w)); hh.Fill(r.mee, r.w); sx += v * r.w; sw += r.w;
        }
        if (xw.size() < 50) continue;
        const LineFit L = FitLine(xw, &hh);
        const int n = g->GetN();
        g->SetPoint(n, sx / sw, L.sig); g->SetPointError(n, 0, L.sigE);
    }
    g->SetMarkerStyle(sty); g->SetMarkerSize(1.8); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
    return g;
}
void Part16() {
    Style();
    // ---------------- F1: X17 line shape (3 x 2)
    {
        TCanvas* c = new TCanvas(Nm("c1").c_str(), "", 2400, 1400);
        c->Divide(3, 2);
        TH1D* h = new TH1D(Nm("hl").c_str(), "", 120, 0, 60); h->SetDirectory(nullptr); h->Sumw2();
        TH1D* hp = new TH1D(Nm("hlp").c_str(), "", 120, 0, 60); hp->SetDirectory(nullptr); hp->Sumw2();
        for (const CR& r : gX17all) { h->Fill(r.mee, r.w / kLsim / 0.5); if (r.m4 >= kPkLo && r.m4 < kPkHi) hp->Fill(r.mee, r.w / kLsim / 0.5); }
        for (int pad = 1; pad <= 2; ++pad) {
            c->cd(pad); gPad->SetLogy(pad == 2);
            TH1D* hd = (TH1D*)hp->Clone(Nm("hd").c_str()); hd->SetDirectory(nullptr);
            Frame(hd, "M(e^{+}e^{-}) [MeV/c^{2}]", "d#sigma/dM(ee) per unit BR [nb/(MeV/c^{2})]");
            hd->SetLineColor(kRed + 1); hd->SetLineWidth(3); hd->SetMarkerColor(kRed + 1); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.9);
            h->SetLineColor(kOrange + 7); h->SetLineWidth(2); h->SetLineStyle(2);
            if (pad == 2) { hd->SetMinimum(h->GetMaximum() * 2e-4); hd->SetMaximum(h->GetMaximum() * 4); } else { hd->SetMinimum(0); hd->SetMaximum(h->GetMaximum() * 1.3); }
            hd->Draw("E1"); h->Draw("HIST SAME"); hd->Draw("E1 SAME");
            TF1* g = new TF1(Nm("fg").c_str(), "gaus", 0, 60);
            g->SetParameters(hp->GetMaximum(), gLinePk.mean, gLinePk.sig);
            hp->Fit(g, "QN0", "", gLinePk.mean - 1.5 * gLinePk.sig, gLinePk.mean + 1.5 * gLinePk.sig);
            g->SetLineColor(kBlack); g->SetLineWidth(2); g->SetRange(0, 60); g->SetNpx(500); g->Draw("SAME");
            TF1* dg = new TF1(Nm("fdg").c_str(), "[0]*((1-[5])*TMath::Gaus(x,[1],[2],1)+[5]*TMath::Gaus(x,[3],[4],1))", 0, 60);
            dg->SetParameters(hp->Integral("width"), gLinePk.m1, gLinePk.s1, gLinePk.m2, gLinePk.s2, gLinePk.fr2);
            dg->SetParLimits(5, 0., 1.);
            hp->Fit(dg, "QN0", "", std::max(0., gLinePk.mean - 8 * gLinePk.sig), std::min(60., gLinePk.mean + 8 * gLinePk.sig));
            dg->SetLineColor(kBlue + 1); dg->SetLineStyle(2); dg->SetLineWidth(3); dg->SetNpx(500); dg->Draw("SAME");
            if (pad == 1) {
                TLegend* lg = new TLegend(0.50, 0.58, 0.99, 0.90); lg->SetTextSize(0.046); lg->SetFillStyle(0);
                lg->AddEntry(hd, "X17, peak window", "lep");
                lg->AddEntry(h, "X17, all M_{4}", "l");
                lg->AddEntry(g, Form("core: #sigma = %.2f", gLinePk.sig), "l");
                lg->AddEntry(dg, Form("2 Gauss.: #sigma %.1f, %.1f", gLinePk.s1, gLinePk.s2), "l");
                lg->Draw();
                Txt(0.17, 0.93, "(a) X17 line in reconstructed M(e^{+}e^{-}), pre", 0.05);
            } else {
                Txt(0.17, 0.93, "(b) the same, log scale: the tails", 0.05);
                Txt(0.56, 0.83, Form("mean %.2f MeV/c^{2}", gLinePk.mean), 0.046);
                Txt(0.56, 0.76, Form("#pm2#sigma: %.3f", gLinePk.f2), 0.046);
                Txt(0.56, 0.69, Form("< -3#sigma: %.3f", gLinePk.flo3), 0.046);
                Txt(0.56, 0.62, Form("> +3#sigma: %.3f", gLinePk.fhi3), 0.046);
                Txt(0.56, 0.55, Form("wide part: %.0f %%", 100 * gLinePk.fr2), 0.046);
            }
        }
        // (c) sigma vs OA, (d) sigma vs p_ee
        for (int k = 0; k < 2; ++k) {
            c->cd(3 + k);
            TGraphErrors* g = k == 0 ? SigmaVs([](const CR& r) { return (double)r.oa; }, {0, 2, 3, 4, 5, 6, 8, 11, 20}, kRed + 1, 20)
                                     : SigmaVs([](const CR& r) { return (double)r.pee; }, {0, 200, 300, 400, 500, 650, 850, 1500}, kRed + 1, 21);
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", 10, 0, k == 0 ? 14 : 1100); fr->SetDirectory(nullptr);
            Frame(fr, k == 0 ? "OA(e^{+}e^{-}) [deg]" : "p(e^{+}e^{-}) [MeV/c]", "#sigma_{core} of M(e^{+}e^{-}) [MeV/c^{2}]");
            fr->SetMinimum(0); fr->SetMaximum(7); fr->Draw("AXIS");
            g->Draw("P SAME");
            TLine* l = new TLine(0, gLinePk.sig, k == 0 ? 14 : 1100, gLinePk.sig); l->SetLineStyle(2); l->Draw();
            Txt(0.17, 0.93, k == 0 ? "(c) width of the X17 line against OA" : "(d) width of the X17 line against p(ee)", 0.05);
            Txt(0.20, 0.20, "M_{4} peak window; dashed: all", 0.046);
        }
        // (e) M4 of X17
        c->cd(5);
        {
            TH1D* hm = new TH1D(Nm("hm").c_str(), "", 120, 0.40, 0.70); hm->SetDirectory(nullptr); hm->Sumw2();
            for (const CR& r : gX17all) hm->Fill(r.m4, r.w / kLsim / 0.0025);
            Frame(hm, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "d#sigma/dM_{4} per unit BR [nb/(GeV/c^{2})]");
            hm->SetLineColor(kRed + 1); hm->SetLineWidth(3); hm->SetMarkerSize(0);
            gPad->SetLogy(1); hm->SetMinimum(hm->GetMaximum() * 1e-4); hm->SetMaximum(hm->GetMaximum() * 5);
            hm->Draw("HIST");
            for (double x : {kPkLo, kPkHi}) { TLine* l = new TLine(x, hm->GetMinimum(), x, hm->GetMaximum()); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
            for (double x : {gOpt0.m4lo, gOpt0.m4hi}) { TLine* l = new TLine(x, hm->GetMinimum(), x, hm->GetMaximum()); l->SetLineStyle(1); l->SetLineColor(kGreen + 3); l->SetLineWidth(2); l->Draw(); }
            Txt(0.17, 0.93, "(e) X17 in M_{4}, pre", 0.05);
            Txt(0.62, 0.83, "dashed: peak", 0.046);
            Txt(0.62, 0.77, "window", 0.046);
            Txt(0.62, 0.69, "green: optimum", 0.046, kGreen + 3);
        }
        // (f) OA shapes of X17 and of the direct decay in the peak window, M(ee) 12-22
        c->cd(6);
        {
            TH1D* a = new TH1D(Nm("oa1").c_str(), "", 40, 0, 20); a->SetDirectory(nullptr);
            TH1D* b = new TH1D(Nm("oa2").c_str(), "", 40, 0, 20); b->SetDirectory(nullptr);
            TH1D* d = new TH1D(Nm("oa3").c_str(), "", 40, 0, 20); d->SetDirectory(nullptr);
            for (const CR& r : gAll) {
                if (!(r.m4 >= kPkLo && r.m4 < kPkHi && r.mee >= 12 && r.mee < 22)) continue;
                if (r.cat == cX17) a->Fill(r.oa, r.w);
                else if (r.cat == cS5) b->Fill(r.oa, r.w);
                else if (IsBkg(r.cat)) d->Fill(r.oa, r.w);
            }
            for (TH1D* x : {a, b, d}) if (x->Integral() > 0) x->Scale(1. / x->Integral() / 0.5);
            Frame(a, "OA(e^{+}e^{-}) [deg]", "normalised [1/deg]");
            a->SetLineColor(kRed + 1); a->SetLineWidth(3);
            b->SetLineColor(kBlue + 1); b->SetLineWidth(3); b->SetLineStyle(2);
            d->SetLineColor(kGray + 2); d->SetLineWidth(3); d->SetLineStyle(3);
            a->SetMaximum(1.3 * std::max(a->GetMaximum(), std::max(b->GetMaximum(), d->GetMaximum()))); a->SetMinimum(0);
            a->Draw("HIST"); b->Draw("HIST SAME"); d->Draw("HIST SAME");
            TLegend* lg = new TLegend(0.45, 0.66, 0.97, 0.90); lg->SetTextSize(0.046); lg->SetFillStyle(0);
            lg->AddEntry(a, "X17", "l"); lg->AddEntry(b, "Gsig5 direct decay", "l"); lg->AddEntry(d, "all background", "l"); lg->Draw();
            Txt(0.17, 0.93, "(f) OA, peak window, M(ee) 12-22 MeV/c^{2}", 0.05);
        }
        c->cd(0);
        Save(c, "s5_opt_x17_lineshape");
    }
    // ---------------- F2: backgrounds in M(ee) and OA near 17 MeV (2 x 2 plots + legend column)
    {
        TCanvas* c = new TCanvas(Nm("c2").c_str(), "", 2500, 1400);
        TPad* pp[4];
        for (int i = 0; i < 4; ++i) {
            const double x0 = (i % 2) * 0.38, y0 = i < 2 ? 0.5 : 0.0;
            pp[i] = new TPad(Nm("p").c_str(), "", x0, y0, x0 + 0.38, y0 + 0.5);
            c->cd(); pp[i]->Draw();
        }
        TPad* pl = new TPad(Nm("pl").c_str(), "", 0.76, 0.0, 1.0, 1.0); c->cd(); pl->Draw();
        Cuts cpre; Cuts copt = gOpt0;
        DrawGroups(pp[0], [&](const CR& r) { return PassC(r, cpre, false, true); }, [](const CR& r) { return (double)r.mee; }, 30, 0, 60, 1. / kLsim / 2.,
                   "M(e^{+}e^{-}) [MeV/c^{2}]", "d#sigma/dM(ee) [nb/(MeV/c^{2})]", true, false, "(a) pre + M_{4} peak window", gOpt0.meelo, gOpt0.meehi, kBRill, 3e-4);
        DrawGroups(pp[1], [&](const CR& r) { return PassC(r, copt, false, true); }, [](const CR& r) { return (double)r.mee; }, 30, 0, 60, 1. / kLsim / 2.,
                   "M(e^{+}e^{-}) [MeV/c^{2}]", "d#sigma/dM(ee) [nb/(MeV/c^{2})]", true, false, "(b) optimum without its M(ee) window", gOpt0.meelo, gOpt0.meehi, kBRill, 3e-4);
        Cuts cpw; cpw.meelo = 12; cpw.meehi = 22;
        Cuts cow = gOpt0; cow.oalo = kNo; cow.oahi = kNo;
        DrawGroups(pp[2], [&](const CR& r) { return PassC(r, cpw, true, true); }, [](const CR& r) { return (double)r.oa; }, 25, 0, 25, 1. / kLsim,
                   "OA(e^{+}e^{-}) [deg]", "d#sigma/dOA [nb/deg]", true, false, "(c) pre + peak window + M(ee) 12-22", gOpt0.oalo, On(gOpt0.oahi) ? gOpt0.oahi : -1, kBRill, 3e-4);
        DrawGroups(pp[3], [&](const CR& r) { return PassC(r, cow, true, true); }, [](const CR& r) { return (double)r.oa; }, 25, 0, 25, 1. / kLsim,
                   "OA(e^{+}e^{-}) [deg]", "d#sigma/dOA [nb/deg]", true, false, "(d) optimum (no OA cut is chosen)", gOpt0.oalo, On(gOpt0.oahi) ? gOpt0.oahi : -1, kBRill, 3e-4);
        GroupLegend(pl, 0.036, kBRill, "dotted lines: optimum M(ee) window");
        c->cd(0);
        Save(c, "s5_opt_x17_bkg_mee_oa");
    }
    // ---------------- F3: 1D scans
    {
        TCanvas* c = new TCanvas(Nm("c3").c_str(), "", 2400, 2100);
        c->Divide(4, 5);
        int iv = 0;
        for (int set = 0; set < 2; ++set) {
            const std::vector<OVar>& vs = set == 0 ? gVars : gExtra;
            for (const OVar& v : vs) {
                if (iv >= 19) break;
                TPad* p = (TPad*)c->cd(iv + 1);
                p->SetLeftMargin(0.17); p->SetBottomMargin(0.19); p->SetTopMargin(0.12);
                bool hasNone = false;
                double lo = 1e30, hi = -1e30;
                for (double x : v.vals) { if (On(x)) { lo = std::min(lo, x); hi = std::max(hi, x); } else hasNone = true; }
                const bool isVtx = std::string(v.key) == "vtx";
                double step = (v.vals.size() > 2) ? (hi - lo) / std::max(1, (int)v.vals.size() - 2) : 1;
                if (isVtx) { lo = 0; hi = 1; step = 1; }
                std::vector<double> xs;
                for (double x : v.vals) xs.push_back(On(x) ? x : (isVtx ? 0. : (v.dir > 0 ? hi + step : lo - step)));
                TGraph* g[3];
                double ymn = 1e30, ymx = 0;
                for (int k = 0; k < 3; ++k) {
                    g[k] = new TGraph();
                    const std::vector<double>& y = k == 0 ? gScanFull[iv] : k == 1 ? gScanA[iv] : gScanB[iv];
                    for (size_t i = 0; i < xs.size(); ++i)
                        if (y[i] == y[i] && y[i] < 1) { g[k]->SetPoint(g[k]->GetN(), xs[i], y[i] * 1e6); ymn = std::min(ymn, y[i] * 1e6); ymx = std::max(ymx, y[i] * 1e6); }
                    g[k]->Sort();
                }
                const double xlo = *std::min_element(xs.begin(), xs.end()), xhi = *std::max_element(xs.begin(), xs.end());
                TH1D* fr = new TH1D(Nm("fr").c_str(), "", 10, xlo - 0.05 * (xhi - xlo + 1e-9), xhi + 0.05 * (xhi - xlo + 1e-9));
                fr->SetDirectory(nullptr);
                Frame(fr, v.label, "UL(BR) [10^{-6}]");
                fr->GetXaxis()->SetLabelSize(0.065); fr->GetYaxis()->SetLabelSize(0.065);
                fr->GetXaxis()->SetTitleSize(0.062); fr->GetYaxis()->SetTitleSize(0.07); fr->GetYaxis()->SetTitleOffset(1.1);
                fr->GetXaxis()->SetNdivisions(505);
                const double ycap = std::min(ymx, ymn * 2.5);
                const double ylo = std::max(0., ymn - 0.1 * (ycap - ymn) - 0.02 * ymn), yhi = ycap + 0.15 * (ycap - ymn) + 0.02 * ymn;
                fr->SetMinimum(ylo); fr->SetMaximum(yhi);
                fr->Draw("AXIS");
                const int col[3] = {kBlack, kBlue + 1, kGreen + 3};
                for (int k = 2; k >= 0; --k) {
                    g[k]->SetLineColor(col[k]); g[k]->SetMarkerColor(col[k]); g[k]->SetLineWidth(k == 0 ? 3 : 2); g[k]->SetMarkerStyle(k == 0 ? 20 : 24);
                    g[k]->SetMarkerSize(0.9); if (k) g[k]->SetLineStyle(2);
                    g[k]->Draw("LP SAME");
                }
                const double xc = On(gOpt0.*(v.field)) ? gOpt0.*(v.field) : (isVtx ? 0. : (v.dir > 0 ? hi + step : lo - step));
                TLine* l = new TLine(xc, ylo, xc, yhi); l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->Draw();
                Txt(0.19, 0.91, set ? Form("%s (extra, not in the descent)", v.key) : v.key, 0.07);
                if (hasNone && !isVtx) Txt(0.93, 0.80, v.dir > 0 ? "right end: no cut" : "left end: no cut", 0.06, kGray + 2, 31);
                ++iv;
            }
        }
        TPad* p = (TPad*)c->cd(20);
        p->cd();
        TLegend* lg = new TLegend(0.03, 0.30, 0.98, 0.97); lg->SetTextSize(0.075); lg->SetFillStyle(0);
        TGraph* d0 = new TGraph(); d0->SetLineWidth(3); d0->SetMarkerStyle(20);
        TGraph* d1 = new TGraph(); d1->SetLineColor(kBlue + 1); d1->SetLineStyle(2); d1->SetMarkerStyle(24); d1->SetMarkerColor(kBlue + 1);
        TGraph* d2 = new TGraph(); d2->SetLineColor(kGreen + 3); d2->SetLineStyle(2); d2->SetMarkerStyle(24); d2->SetMarkerColor(kGreen + 3);
        TLine* dl = new TLine(); dl->SetLineColor(kRed + 1); dl->SetLineWidth(2);
        lg->AddEntry(d0, "full sample", "lp"); lg->AddEntry(d1, "half A (even grp), w x 2", "lp"); lg->AddEntry(d2, "half B (odd grp), w x 2", "lp");
        lg->AddEntry(dl, "chosen value", "l");
        lg->Draw();
        Txt(0.04, 0.17, gFomMC ? "exp. 90% CL UL on BR, incl. MC error of B" : "expected 90% CL UL on BR, B known", 0.068);
        Txt(0.04, 0.06, kLdNote, 0.072);
        c->cd(0);
        Save(c, "s5_opt_x17_scan1d");
    }
    // ---------------- F4: 2D maps
    {
        TCanvas* c = new TCanvas(Nm("c4").c_str(), "", 1900, 1580);
        TPad* pm[4];
        for (int i = 0; i < 4; ++i) {
            const double x0 = (i % 2) * 0.5, y0 = i < 2 ? 0.53 : 0.06;
            pm[i] = new TPad(Nm("pm").c_str(), "", x0, y0, x0 + 0.5, y0 + 0.47);
            c->cd(); pm[i]->Draw();
        }
        TPad* ps = new TPad(Nm("ps").c_str(), "", 0, 0, 1, 0.05); c->cd(); ps->Draw();
        ps->cd(); { TLatex t; t.SetNDC(); t.SetTextSize(0.42); t.SetTextFont(42); t.SetTextAlign(22);
                    t.DrawLatex(0.5, 0.5, Form("UL on BR relative to the optimum; other cuts at the optimum; %s", kLdNote)); }
        struct M2 { const char* t; double Cuts::*fx; double Cuts::*fy; std::vector<double> xv, yv; const char* xt; const char* yt; };
        std::vector<M2> ms = {
            {"(a) M(ee) window", &Cuts::meelo, &Cuts::meehi, Rng(6, 17, 0.5), Rng(17, 34, 0.5), "M(ee) lower edge [MeV/c^{2}]", "M(ee) upper edge [MeV/c^{2}]"},
            {"(b) OA cuts", &Cuts::oalo, &Cuts::oahi, Rng(0, 8, 0.5), Rng(5, 26, 1), "OA lower cut [deg] (0: none)", "OA upper cut [deg] (26: none)"},
            {"(c) M_{4} window", &Cuts::m4lo, &Cuts::m4hi, Rng(0.518, 0.5455, 0.0025), Rng(0.548, 0.5755, 0.0025), "M_{4} lower edge [GeV/c^{2}]", "M_{4} upper edge [GeV/c^{2}]"},
            {"(d) #alpha_{lab}, #alpha_{#eta}", &Cuts::labhi, &Cuts::etalo, Rng(20, 95, 5), Rng(85, 175, 5), "#alpha_{lab} upper cut [deg] (95: none)", "#alpha_{#eta} lower cut [deg] (85: none)"}};
        const double best = FomCuts(gScan, gOpt0, 0) * 1e6;
        for (size_t m = 0; m < ms.size(); ++m) {
            TPad* p = pm[m]; p->cd();
            p->SetRightMargin(0.18); p->SetLeftMargin(0.19); p->SetTopMargin(0.13); p->SetBottomMargin(0.16);
            const M2& M = ms[m];
            const int nx = (int)M.xv.size(), ny = (int)M.yv.size();
            const double dx = M.xv[1] - M.xv[0], dy = M.yv[1] - M.yv[0];
            TH2D* h = new TH2D(Nm("h2").c_str(), "", nx, M.xv[0] - dx / 2, M.xv[nx - 1] + dx / 2, ny, M.yv[0] - dy / 2, M.yv[ny - 1] + dy / 2);
            h->SetDirectory(nullptr);
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j) {
                    Cuts t = gOpt0;
                    double xv = M.xv[i], yv = M.yv[j];
                    if (m == 1) { if (xv == 0) xv = kNo; if (yv >= 26) yv = kNo; }
                    if (m == 3) { if (xv >= 95) xv = kNo; if (yv <= 85) yv = kNo; }
                    if ((m == 0 && yv <= xv) || (m == 2 && yv <= xv) || (m == 1 && On(xv) && On(yv) && yv <= xv)) continue;
                    t.*(M.fx) = xv; t.*(M.fy) = yv;
                    const double f = FomCuts(gScan, t, 0) * 1e6;
                    if (f < 1e3) h->SetBinContent(i + 1, j + 1, f / best);
                }
            Frame(h, M.xt, M.yt);
            h->GetYaxis()->SetTitleOffset(1.6); h->GetXaxis()->SetNdivisions(506);
            h->GetZaxis()->SetTitle("UL / UL(optimum)"); h->GetZaxis()->SetTitleSize(0.05); h->GetZaxis()->SetLabelSize(0.045);
            h->GetZaxis()->SetTitleOffset(1.25);
            h->SetMinimum(0.95); h->SetMaximum(1.6);
            h->Draw("COLZ");
            TH2D* hc = (TH2D*)h->Clone(Nm("hc").c_str());
            double lv[3] = {1.02, 1.05, 1.10};
            hc->SetContour(3, lv); hc->SetLineColor(kBlack); hc->SetLineWidth(2); hc->Draw("CONT3 SAME");
            double ox = gOpt0.*(M.fx), oy = gOpt0.*(M.fy);
            if (m == 1) { if (!On(ox)) ox = 0; if (!On(oy)) oy = 26; }
            if (m == 3) { if (!On(ox)) ox = 95; if (!On(oy)) oy = 85; }
            TGraph* gp = new TGraph(1); gp->SetPoint(0, ox, oy); gp->SetMarkerStyle(29); gp->SetMarkerSize(3.5); gp->SetMarkerColor(kRed + 1); gp->Draw("P SAME");
            Txt(0.19, 0.935, M.t, 0.055);
            Txt(0.19, 0.885, "contours 1.02, 1.05, 1.10; star: optimum", 0.046);
        }
        c->cd(0);
        Save(c, "s5_opt_x17_maps");
    }
    // ---------------- F5: M4 at the optimum (no M4 window) and slides
    {
        TCanvas* c = new TCanvas(Nm("c5").c_str(), "", 2500, 900);
        TPad* p1 = new TPad(Nm("p").c_str(), "", 0.0, 0.0, 0.38, 1.0); c->cd(); p1->Draw();
        TPad* p2 = new TPad(Nm("p").c_str(), "", 0.38, 0.0, 0.76, 1.0); c->cd(); p2->Draw();
        TPad* pl = new TPad(Nm("pl").c_str(), "", 0.76, 0.0, 1.0, 1.0); c->cd(); pl->Draw();
        Cuts co = gOpt0; Cuts cs = gSlW;
        DrawGroups(p1, [&](const CR& r) { return PassC(r, co, true, false); }, [](const CR& r) { return (double)r.m4; }, 60, 0.40, 0.70,
                   1. / kLsim / 0.005, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "d#sigma/dM_{4} [nb/(GeV/c^{2})]", true, false,
                   "(a) optimum without its M_{4} window", gOpt0.m4lo, gOpt0.m4hi, kBRill, 2e-2);
        DrawGroups(p2, [&](const CR& r) { return PassC(r, cs, true, false); }, [](const CR& r) { return (double)r.m4; }, 60, 0.40, 0.70,
                   1. / kLsim / 0.005, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "d#sigma/dM_{4} [nb/(GeV/c^{2})]", true, false,
                   Form("(b) slides + M(ee) %.1f-%.1f MeV/c^{2}", gOpt0.meelo, gOpt0.meehi), kPkLo, kPkHi, kBRill, 2e-2);
        GroupLegend(pl, 0.042, kBRill, "dotted lines: M_{4} windows");
        c->cd(0);
        Save(c, "s5_opt_x17_m4");
    }
    // ---------------- F6: shape analysis (row 1: two fits + legend; row 2: q curves, limit summary)
    {
        TCanvas* c = new TCanvas(Nm("c6").c_str(), "", 2500, 1500);
        TPad* pa = new TPad(Nm("p").c_str(), "", 0.00, 0.5, 0.38, 1.0); c->cd(); pa->Draw();
        TPad* pb = new TPad(Nm("p").c_str(), "", 0.38, 0.5, 0.76, 1.0); c->cd(); pb->Draw();
        TPad* pl = new TPad(Nm("p").c_str(), "", 0.76, 0.5, 1.00, 1.0); c->cd(); pl->Draw();
        TPad* pc = new TPad(Nm("p").c_str(), "", 0.00, 0.0, 0.45, 0.5); c->cd(); pc->Draw();
        TPad* pd = new TPad(Nm("p").c_str(), "", 0.45, 0.0, 1.00, 0.5); c->cd(); pd->Draw();
        const int order[9] = {cOth, cVec, cConv, cEta, cPi0D, cS2, cS3, cS4, cS5};
        TH1D* hxL = nullptr; TH1D* hxtL = nullptr; std::vector<TH1D*> hL;
        for (int k = 0; k < 2; ++k) {
            ShapeIn& in = k == 0 ? gInOpt : gInSl;
            const ShapeRes& R = k == 0 ? gShOpt : gShSl;
            std::vector<std::vector<double>>& grp = k == 0 ? gGrpOpt : gGrpSl;
            TPad* p = k == 0 ? pa : pb;
            p->cd(); p->SetLogy(1);
            THStack* st = new THStack(Nm("st").c_str(), "");
            double ymax = 0;
            std::vector<double> tot(in.nb, 0.);
            for (int oi = 0; oi < 9; ++oi) {
                const int cat = order[oi];
                TH1D* h = new TH1D(Nm("hs").c_str(), "", in.nb, in.lo, in.hi); h->SetDirectory(nullptr);
                for (int i = 0; i < in.nb; ++i) { h->SetBinContent(i + 1, grp[cat][i]); tot[i] += grp[cat][i]; }
                h->SetFillColor(kCatCol[cat]); h->SetLineColor(kCatCol[cat]);
                st->Add(h);
                if (k == 0) hL.push_back(h);
            }
            for (double x : tot) ymax = std::max(ymax, x);
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", in.nb, in.lo, in.hi); fr->SetDirectory(nullptr);
            Frame(fr, "M(e^{+}e^{-}) [MeV/c^{2}]", "expected rows / (1 MeV/c^{2})");
            fr->SetMaximum(ymax * 30); fr->SetMinimum(std::max(0.05, ymax * 2e-4));
            fr->Draw("AXIS");
            st->Draw("HIST SAME");
            TH1D* hx = new TH1D(Nm("hx").c_str(), "", in.nb, in.lo, in.hi); hx->SetDirectory(nullptr);
            TH1D* hxt = new TH1D(Nm("hxt").c_str(), "", in.nb, in.lo, in.hi); hxt->SetDirectory(nullptr);
            const std::vector<double>& tpl = k == 0 ? gTplOpt : gTplSl;
            for (int i = 0; i < in.nb; ++i) { hx->SetBinContent(i + 1, in.s[i] * R.ul[1]); hxt->SetBinContent(i + 1, tpl[i] * kBRill); }
            hx->SetLineColor(kRed + 1); hx->SetLineWidth(4);
            hxt->SetLineColor(kOrange + 7); hxt->SetLineWidth(3); hxt->SetLineStyle(7);
            hx->Draw("HIST SAME"); hxt->Draw("HIST SAME");
            if (k == 0) { hxL = hx; hxtL = hxt; }
            gPad->RedrawAxis();
            Txt(0.17, 0.93, k == 0 ? "(a) optimum without its M(ee) window" : "(b) slides, M_{4} peak window", 0.05);
            Txt(0.17, 0.85, Form("UL (shape, B norm. free) = %.2f #times 10^{-6}", R.ul[1] * 1e6), 0.046, kRed + 1);
            Txt(0.17, 0.78, kLdNote, 0.046);
        }
        pl->cd();
        {
            TLegend* lg = new TLegend(0.02, 0.04, 0.99, 0.97); lg->SetTextSize(0.05); lg->SetFillStyle(0);
            lg->AddEntry(hxL, "X17 Gaussian at the expected UL", "l");
            lg->AddEntry(hxtL, "X17 simulated, BR = 10^{-5}", "l");
            for (int oi = 8; oi >= 0; --oi) lg->AddEntry(hL[oi], kCatTex[order[oi]], "f");
            lg->Draw();
        }
        // (c) q(BR) curves
        pc->cd();
        {
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", 10, 0, 12); fr->SetDirectory(nullptr);
            Frame(fr, "BR(#eta#rightarrow#pi^{+}#pi^{-}X17)#timesBR(X17#rightarrowe^{+}e^{-}) [10^{-6}]", "q(BR), Asimov B-only data");
            fr->SetMinimum(0); fr->SetMaximum(8); fr->Draw("AXIS");
            TLine* l = new TLine(0, 2.7055, 12, 2.7055); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
            TLegend* lg = new TLegend(0.16, 0.52, 0.58, 0.90); lg->SetTextSize(0.046); lg->SetFillStyle(0);
            struct Cv { ShapeIn* in; int mode; int col; int sty; const char* nm; };
            static ShapeIn one, oneS;
            one.nb = 1; one.b0 = {gFinOpt.e.b[cS5] * kToRows}; one.b1 = {(gFinOpt.e.B - gFinOpt.e.b[cS5]) * kToRows}; one.s = {gFinOpt.f.s1};
            oneS.nb = 1; oneS.b0 = {gFinSlW.e.b[cS5] * kToRows}; oneS.b1 = {(gFinSlW.e.B - gFinSlW.e.b[cS5]) * kToRows}; oneS.s = {gFinSlW.f.s1};
            std::vector<Cv> cv = {{&one, 0, kBlack, 1, "opt., counting"}, {&gInOpt, 0, kRed + 1, 1, "opt., shape, B known"},
                                  {&gInOpt, 1, kRed + 1, 2, "opt., shape, B free"}, {&oneS, 0, kGray + 2, 1, "slides + M(ee) win., count."},
                                  {&gInSl, 0, kBlue + 1, 1, "slides, shape, B known"}, {&gInSl, 1, kBlue + 1, 2, "slides, shape, B free"}};
            for (const Cv& x : cv) {
                TGraph* g = new TGraph();
                for (int i = 0; i <= 120; ++i) { const double br = i * 0.1e-6; g->SetPoint(i, br * 1e6, QAsimov(*x.in, br, x.mode)); }
                g->SetLineColor(x.col); g->SetLineStyle(x.sty); g->SetLineWidth(3); g->Draw("L SAME");
                lg->AddEntry(g, x.nm, "l");
            }
            lg->Draw();
            Txt(0.17, 0.93, "(c) q(BR); the expected CLs limit is at q = 2.71", 0.05);
        }
        // (d) summary of limits
        pd->cd();
        {
            struct It { const char* nm; double v; int col; };
            std::vector<It> it = {
                {"opt., count., + MC error of B", gFinOpt.f.ulTot, kBlack},
                {"opt., count., B known", gFinOpt.f.ulStat, kBlack},
                {"v1, count., + MC error of B", gFinV1.f.ulTot, kGray + 2},
                {"v1, count., B known", gFinV1.f.ulStat, kGray + 2},
                {"opt., count., median CLs", gFinOpt.ulMedCls, kBlack},
                {"opt., count., Feldman-Cousins", gFinOpt.ulFC, kBlack},
                {"opt., count., B + T3 + T4", gFinOpt.f34.ulStat, kGray + 2},
                {"opt., shape (Gauss), B known", gShOpt.ul[0], kRed + 1},
                {"opt., shape (Gauss), B free", gShOpt.ul[1], kRed + 1},
                {"opt., shape (templ.), B free", gShOpt.ulT[1], kRed + 1},
                {"slides + M(ee) win., count.", gFinSlW.f.ulStat, kBlue + 1},
                {"slides, shape (Gauss), B known", gShSl.ul[0], kBlue + 1},
                {"slides, shape (Gauss), B free", gShSl.ul[1], kBlue + 1},
                {"pre + M(ee) 12-22, count.", gFinPre.f.ulStat, kGray + 1}};
            const int n = (int)it.size();
            double vmax = 0; for (auto& x : it) vmax = std::max(vmax, x.v);
            pd->SetLeftMargin(0.45); pd->SetRightMargin(0.03); pd->SetTopMargin(0.16);
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", 10, 0, 1.2 * vmax * 1e6); fr->SetDirectory(nullptr);
            Frame(fr, "expected 90% CL UL on BR [10^{-6}]", "");
            fr->SetMinimum(0.3); fr->SetMaximum(n + 0.7);
            fr->GetYaxis()->SetLabelSize(0); fr->GetYaxis()->SetTickLength(0);
            fr->Draw("AXIS");
            for (int i = 0; i < n; ++i) {
                const double y = n - i;
                TBox* b = new TBox(0, y - 0.33, it[i].v * 1e6, y + 0.33); b->SetFillColor(it[i].col); b->SetFillStyle(1001); b->Draw();
                TLatex t; t.SetTextAlign(32); t.SetTextSize(0.05); t.SetTextFont(42);
                t.DrawLatex(-0.02 * fr->GetXaxis()->GetXmax(), y, it[i].nm);
                TLatex v; v.SetTextAlign(12); v.SetTextSize(0.05); v.SetTextFont(42);
                v.DrawLatex(it[i].v * 1e6 + 0.01 * fr->GetXaxis()->GetXmax(), y, Form("%.2f", it[i].v * 1e6));
            }
            Txt(0.46, 0.93, "(d) expected limits (90% CL); B known unless stated", 0.05);
            Txt(0.46, 0.87, kLdNote, 0.045);
        }
        c->cd(0);
        Save(c, "s5_opt_x17_shape");
    }
}

}  // namespace s5x

void s5_opt_x17(int parts = 0) {
    using namespace s5x;
    if (parts == 0) parts = 31;
    TStopwatch sw;
    if (!LoadAll()) return;
    if (parts & 1) Part1();
    if (parts & (2 | 4 | 8 | 16)) Part4();       // the optimum is needed by every later part
    if (parts & 2) Part2();
    if (parts & (8 | 16)) Part8();
    if (parts & 16) Part16();
    printf("s5_opt_x17 done in %.0f s\n", sw.RealTime());
}
