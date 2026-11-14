// ============================================================================
// sim_cuts.C  --  cut efficiencies and cut optimisation from the simulation truth
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SMASH std + lep
// owned rows (partition mode A, weight w = PDG BR) unless a variant is named.
// No combinatorial background in simulation: signal = template 1 split by the
// four-particle class (bump pi0 Dalitz, bump conversion, peak conversion);
// physical background = template 2. Templates 3 and 4 only for information.
// eta -> pi+pi-e+e- (direct) is NOT in the simulation (eta set missing).
//
// 1. Efficiency of every cut alone (base "pre" = isBest >= 0, PT3, vz > -160,
//    z > -500, cor_ok; the vz family starts from "pre0" = pre without
//    vz > -160) for the signal classes and for the physical background in the
//    two windows, weighted, with event-group errors (pass and fail rows of one
//    group correlated). Compared with the data-derived efficiencies of the
//    data study (tables of plots/scan_peak.C, scan_bump.C, scan_extra.C).
// 2. Optimisation from truth: S and B in the peak window and in the bump
//    window, FOM = S / sqrt(S + B) (expected counts at L_data) and S/B, per cut
//    value and in the maps (oa_lab, M(pipi) lower) and (OA lower, OA upper);
//    repeated with B multiplied by R = data (N - CB) / simulation (T1 + T2) in
//    0.60 - 0.80 GeV/c2, R evaluated for every selection.
// 3. Ribbon cut flows of the data study (same order, nominal values and bands).
//
// Windows: bump 0.4550 - 0.5225 (data); peak: the simulated peak conversion is
// fitted (Gaussian, chi2) and the simulation window is mean +- 1.9 sigma (the
// rule of the data window); the data window 0.53364 - 0.55416 is a variant.
//
// Absolute scale: sigma = sum w / L_sim (L_sim = 221.05 nb^-1); data
// (N - CB) / L_data, L_data = 6080.5 nb^-1 (UNCONFIRMED).
//
// Output: figs/cuts_*.png (.pdf for the main ones), tables/cuts_*.txt
// Run from /home/przygoda/HADES/FAT/FAT (about 2 min, about 7 GB of memory):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_cuts", kTRUE)' 'plots/eta_sim/sim_cuts.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"

namespace scut {
using namespace simeng;

// ============================================================================
// constants, windows
// ============================================================================
const double kNoCut = -999.;
const double kPkLoD = scan::kPeakLo, kPkHiD = scan::kPeakHi;   // data peak window 0.53364 - 0.55416
const double kBuLo = scan::kBumpLo, kBuHi = scan::kBumpHi;      // bump window 0.4550 - 0.5225
const double kRLo = 0.60, kRHi = 0.80;                          // window of R = data / simulation
double gPkLo = kPkLoD, gPkHi = kPkHiD;                          // simulation peak window (from the fit)
double gPkMean = 0, gPkMeanE = 0, gPkSig = 0, gPkSigE = 0, gPkChi2 = 0;
int gPkNdf = 0;
const char* const kDataScan = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan";
const double kTabLsim = kLsim;

FILE* gOut = nullptr;   // current table
void P(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(gOut ? gOut : stdout, fmt, ap);
    va_end(ap);
}

// ============================================================================
// targets: (class, window) combinations accumulated in one pass
// ============================================================================
enum Tgt { tBD = 0, tBC, tBDC, tPC, tTag, tSbu, tSpk, tSpkD, tB2bu, tB2pk, tB2pkD, tSimBu, tSimPk, tSimPkD, tSimR,
           tT3bu, tT3pk, tT4bu, tT4pk, tPCinBu, tBUinPk, NT };
const char* TgtName(int t) {
    static const char* n[NT] = {"bump pi0 Dalitz (all M4)", "bump conversion (all M4)", "bump D+C (all M4)", "peak conversion (all M4)",
                                "bump D+C, photon-tagged (all M4)", "S_bump: bump D+C in bump window",
                                "S_peak: peak conv. in sim peak window", "peak conv. in data peak window",
                                "B_bump: phys. bkg in bump window", "B_peak: phys. bkg in sim peak window",
                                "phys. bkg in data peak window", "T1+T2 in bump window", "T1+T2 in sim peak window",
                                "T1+T2 in data peak window", "T1+T2 in 0.60-0.80", "fake lepton in bump window",
                                "fake lepton in sim peak window", "diff. sources in bump window",
                                "diff. sources in sim peak window", "peak conv. in bump window",
                                "bump D+C in sim peak window"};
    return (t >= 0 && t < NT) ? n[t] : "?";
}

// compact simulation row (only rows that belong to at least one target)
struct MRow {
    float m4, oa, oa_lab, mpipi, oa_eta, mm4, pt4, vz, dzpi;
    signed char nvtxpi, vtx4;
    unsigned char own, tmpl, cls, tag;
    unsigned int mask;
    int gid;
    float dal;          // eta -> pi+pi-pi0 Dalitz-plot factor (bump classes), 1 otherwise
    double w, wdb;
};
// compact data row (window of R only)
struct DRow {
    float m4, oa, oa_lab, mpipi, oa_eta, mm4, pt4, vz, dzpi;
    signed char nvtxpi, vtx4;
    unsigned char s;    // 0 us, 1 pp, 2 mm
};
std::vector<MRow> gRows;
std::vector<DRow> gData;

inline unsigned TgtMask(int tmpl, int cls, bool tag, double m4) {
    unsigned m = 0;
    const bool sig = tmpl == 1, bu = sig && (cls == 2 || cls == 3), pc = sig && cls == 4, sim = tmpl == 1 || tmpl == 2;
    const bool inBu = m4 >= kBuLo && m4 < kBuHi, inPk = m4 >= gPkLo && m4 < gPkHi;
    const bool inPkD = m4 >= kPkLoD && m4 < kPkHiD, inR = m4 >= kRLo && m4 < kRHi;
    if (sig && cls == 2) m |= 1u << tBD;
    if (sig && cls == 3) m |= 1u << tBC;
    if (bu) m |= 1u << tBDC;
    if (pc) m |= 1u << tPC;
    if (bu && tag) m |= 1u << tTag;
    if (inBu) {
        if (bu) m |= 1u << tSbu;
        if (tmpl == 2) m |= 1u << tB2bu;
        if (sim) m |= 1u << tSimBu;
        if (tmpl == 3) m |= 1u << tT3bu;
        if (tmpl == 4) m |= 1u << tT4bu;
        if (pc) m |= 1u << tPCinBu;
    }
    if (inPk) {
        if (pc) m |= 1u << tSpk;
        if (tmpl == 2) m |= 1u << tB2pk;
        if (sim) m |= 1u << tSimPk;
        if (tmpl == 3) m |= 1u << tT3pk;
        if (tmpl == 4) m |= 1u << tT4pk;
        if (bu) m |= 1u << tBUinPk;
    }
    if (inPkD) {
        if (pc) m |= 1u << tSpkD;
        if (tmpl == 2) m |= 1u << tB2pkD;
        if (sim) m |= 1u << tSimPkD;
    }
    if (inR && sim) m |= 1u << tSimR;
    return m;
}

// eta -> pi+ pi- pi0 Dalitz plot: |A|^2 = 1 + a y + b y^2 + d x^2 + f y^3 (x terms dropped: x is not in the
// compact file), y = 3 T_pi0 / Q - 1. Parameters: KLOE-2 2016 (values quoted from memory, UNVERIFIED);
// normalised to the phase-space average (TGenPhaseSpace) so that the total eta -> pi+pi-pi0 rate is unchanged.
const double kDalA = -1.095, kDalB = 0.145, kDalF = 0.141;
const double kMetaG = 0.547862, kMpi0G = 0.1349768, kMpiG = 0.13957;
const double kQeta = kMetaG - 2 * kMpiG - kMpi0G;
double gDalNorm = 1.0;
inline double DalY(double mpp) {
    const double Tpi0 = (kMetaG * kMetaG + kMpi0G * kMpi0G - mpp * mpp) / (2 * kMetaG) - kMpi0G;
    return 3 * Tpi0 / kQeta - 1;
}
inline double DalW(double mpp) { const double y = DalY(mpp); return 1 + kDalA * y + kDalB * y * y + kDalF * y * y * y; }

// ============================================================================
// cuts (one code for simulation and data rows)
// ============================================================================
struct Cuts {
    double vzLo = -160;       // "pre": vz > -160; -1e9 = binding z > -500 only ("pre0")
    double oaLo = kNoCut, oaHi = kNoCut, labHi = kNoCut, mppLo = kNoCut, mppHi = kNoCut, etaLo = kNoCut, mm4Hi = kNoCut,
           pt4Lo = kNoCut, dzHi = kNoCut;
    int vtx = 0;              // 0 none, 1 >= 1 pion, 2 both pions, 3 all four tracks in the vertex fit
};
inline bool On(double v) { return v != kNoCut; }
template <class R> inline bool Pass(const R& r, const Cuts& c) {
    if (c.vzLo > -1e8 && !(r.vz > c.vzLo)) return false;
    if (On(c.oaLo) && c.oaLo > 0 && !(r.oa > c.oaLo)) return false;
    if (On(c.oaHi) && !(r.oa < c.oaHi)) return false;
    if (On(c.labHi) && !(r.oa_lab < c.labHi)) return false;
    if (On(c.mppLo) && !(r.mpipi > c.mppLo)) return false;
    if (On(c.mppHi) && !(r.mpipi < c.mppHi)) return false;
    if (On(c.etaLo) && !(r.oa_eta > c.etaLo)) return false;
    if (On(c.mm4Hi) && !(r.mm4 < c.mm4Hi)) return false;
    if (On(c.pt4Lo) && c.pt4Lo > 0 && !(r.pt4 > c.pt4Lo)) return false;
    if (On(c.dzHi) && !(r.dzpi < c.dzHi)) return false;
    if (c.vtx == 1 && r.nvtxpi < 1) return false;
    if (c.vtx == 2 && r.nvtxpi < 2) return false;
    if (c.vtx == 3 && !r.vtx4) return false;
    return true;
}
Cuts Pre() { return Cuts(); }
Cuts Pre0() { Cuts c; c.vzLo = -1e9; return c; }
Cuts PeakStart() { Cuts c = Pre0(); c.labHi = 55; c.mppLo = 0.30; return c; }           // scan_extra.C (txt), peak
Cuts BumpStart() { Cuts c = Pre0(); c.labHi = 40; c.mppLo = 0.30; c.mppHi = 0.41; return c; }  // scan_extra.C (txt), bump
Cuts PeakStartMD() { Cuts c = Pre0(); c.labHi = 50; c.mppHi = 0.46; return c; }          // scan_extra.md, peak
Cuts BumpStartMD() { Cuts c = Pre0(); c.labHi = 50; c.mppHi = 0.42; return c; }          // scan_extra.md, bump

// ============================================================================
// weighted sums with event-group variance (pass / fail correlated)
// ============================================================================
struct Acc {
    double Sa = 0, Sp = 0, aa = 0, ab = 0, bb = 0, Na = 0, Np = 0;
    double eps() const { return Sa > 0 ? Sp / Sa : 0.; }
    double epsErr() const {
        if (Sa <= 0) return 0.;
        const double e = Sp / Sa, v = (1 - e) * (1 - e) * aa - 2 * e * (1 - e) * ab + e * e * bb;
        return v > 0 ? std::sqrt(v) / Sa : 0.;
    }
    double Va() const { return aa + 2 * ab + bb; }
    double NeffP() const { return aa > 0 ? Sp * Sp / aa : 0.; }
    double NeffA() const { const double v = Va(); return v > 0 ? Sa * Sa / v : 0.; }
    double xs() const { return Sp / kLsim; }                 // nb, rows passing the selection
    double xsErr() const { return std::sqrt(aa) / kLsim; }
};
struct EOpt { int mode = 0; bool db = false; bool dal = false; };
std::vector<Acc> Eval(const Cuts& base, const Cuts& sel, const EOpt& o = EOpt()) {
    std::vector<Acc> acc(NT);
    double ga[NT], gb[NT];
    for (int t = 0; t < NT; ++t) ga[t] = gb[t] = 0;
    unsigned touched = 0;
    int cur = std::numeric_limits<int>::min();
    const unsigned char ob = o.mode ? 2 : 1;
    for (size_t i = 0; i <= gRows.size(); ++i) {
        const bool end = i == gRows.size();
        if (end || gRows[i].gid != cur) {
            unsigned m = touched;
            while (m) {
                const int t = __builtin_ctz(m);
                m &= m - 1;
                acc[t].aa += ga[t] * ga[t]; acc[t].ab += ga[t] * gb[t]; acc[t].bb += gb[t] * gb[t];
                ga[t] = gb[t] = 0;
            }
            touched = 0;
            if (end) break;
            cur = gRows[i].gid;
        }
        const MRow& r = gRows[i];
        if (!(r.own & ob)) continue;
        if (!Pass(r, base)) continue;
        const bool p = Pass(r, sel);
        const double w = (o.db ? r.wdb : r.w) * (o.dal ? (double)r.dal : 1.0);
        unsigned m = r.mask;
        while (m) {
            const int t = __builtin_ctz(m);
            m &= m - 1;
            Acc& a = acc[t];
            a.Sa += w; a.Na += 1;
            if (p) { a.Sp += w; a.Np += 1; ga[t] += w; } else gb[t] += w;
            touched |= 1u << t;
        }
    }
    return acc;
}

// data (N - CB) in 0.60 - 0.80 [nb], CB = 2 sqrt(N++ N--) in 10 MeV bins
struct DW { double sig = 0, err = 0, N = 0, CB = 0; bool ok = false; };
DW DataWin(const Cuts& c) {
    DW o;
    if (On(c.dzHi)) return o;      // pion z not in the data compact file
    const int nb = 20;
    double n[3][20];
    for (int k = 0; k < 3; ++k) for (int b = 0; b < nb; ++b) n[k][b] = 0;
    for (const DRow& r : gData) {
        if (!Pass(r, c)) continue;
        const int b = (int)((r.m4 - kRLo) / 0.01);
        if (b < 0 || b >= nb) continue;
        n[r.s][b] += 1;
    }
    double s = 0, v = 0;
    for (int b = 0; b < nb; ++b) {
        const double cb = 2 * std::sqrt(n[1][b] * n[2][b]);
        s += n[0][b] - cb; v += n[0][b] + n[1][b] + n[2][b];
        o.N += n[0][b]; o.CB += cb;
    }
    o.sig = s / kLdata; o.err = std::sqrt(v) / kLdata; o.ok = true;
    return o;
}

// one evaluated selection
struct Res {
    std::vector<Acc> a;
    DW dw, dw0;                 // data at the selection and at its base
    double R = 1, Rerr = 0, R0 = 1, R0err = 0;
    bool Rok = false;
};
Res Run(const Cuts& base, const Cuts& sel, const EOpt& o = EOpt(), bool data = true) {
    Res r;
    r.a = Eval(base, sel, o);
    if (!data) return r;
    r.dw = DataWin(sel); r.dw0 = DataWin(base);
    const Acc& s = r.a[tSimR];
    if (r.dw0.ok && s.Sa > 0) {
        const double sim0 = s.Sa / kLsim;
        r.R0 = r.dw0.sig / sim0;
        r.R0err = r.R0 * std::sqrt(std::pow(r.dw0.err / r.dw0.sig, 2) + s.Va() / (s.Sa * s.Sa));
    }
    if (r.dw.ok && s.Sp > 0) {
        const double sim = s.Sp / kLsim;
        r.R = r.dw.sig / sim;
        r.Rerr = r.R * std::sqrt(std::pow(r.dw.err / r.dw.sig, 2) + s.aa / (s.Sp * s.Sp));
        r.Rok = true;
    } else { r.R = r.R0; r.Rerr = r.R0err; r.Rok = false; }
    return r;
}

// figure of merit of a window: 0 = peak (sim window), 1 = bump, 2 = peak (data window)
struct FR { double S = 0, Se = 0, B = 0, Be = 0, SB = 0, F = 0, rel = 0, relErr = 0, sbRel = 0; };
inline int TS(int win) { return win == 0 ? tSpk : win == 1 ? tSbu : tSpkD; }
inline int TB(int win) { return win == 0 ? tB2pk : win == 1 ? tB2bu : tB2pkD; }
// R, R0: factors on B at the selection and at the base (1 = simulation only)
FR Fom(const std::vector<Acc>& a, int win, double R = 1., double R0 = 1.) {
    const Acc& s = a[TS(win)];
    const Acc& b = a[TB(win)];
    FR f;
    f.S = s.Sp / kLsim; f.Se = std::sqrt(s.aa) / kLsim;
    f.B = R * b.Sp / kLsim; f.Be = R * std::sqrt(b.aa) / kLsim;
    const double S0 = s.Sa / kLsim, B0 = R0 * b.Sa / kLsim;
    f.SB = f.B > 0 ? f.S / f.B : 0.;
    f.sbRel = (B0 > 0 && S0 > 0 && f.B > 0) ? f.SB / (S0 / B0) : 0.;
    f.F = (f.S + f.B) > 0 ? std::sqrt(kLdata) * f.S / std::sqrt(f.S + f.B) : 0.;
    const double F0 = (S0 + B0) > 0 ? std::sqrt(kLdata) * S0 / std::sqrt(S0 + B0) : 0.;
    f.rel = F0 > 0 ? f.F / F0 : 0.;
    const double es = s.eps(), eb = b.eps();
    if (f.S + f.B > 0) {
        const double kS = 1 - f.S / (2 * (f.S + f.B)), kB = f.B / (2 * (f.S + f.B));
        const double rs = es > 0 ? s.epsErr() / es : 0., rb = eb > 0 ? b.epsErr() / eb : 0.;
        f.relErr = f.rel * std::sqrt(kS * kS * rs * rs + kB * kB * rb * rb);
    }
    return f;
}

// ============================================================================
// cut families (each cut alone)
// ============================================================================
enum Kind { kVz, kOaLo, kOaHi, kLab, kMppLo, kMppHi, kEta, kMm4, kPt4, kDz, kVtx };
void SetCut(Cuts& c, int k, double v) {
    switch (k) {
        case kVz: c.vzLo = (v == kNoCut) ? -1e9 : v; break;
        case kOaLo: c.oaLo = v; break;
        case kOaHi: c.oaHi = v; break;
        case kLab: c.labHi = v; break;
        case kMppLo: c.mppLo = v; break;
        case kMppHi: c.mppHi = v; break;
        case kEta: c.etaLo = v; break;
        case kMm4: c.mm4Hi = v; break;
        case kPt4: c.pt4Lo = v; break;
        case kDz: c.dzHi = v; break;
        case kVtx: c.vtx = (int)v; break;
    }
}
struct Fam {
    int kind;
    std::string key, title;
    std::vector<double> v;
    bool pre0;          // base = pre0 (vz family)
    bool starts;        // also evaluated at the scan_extra start selections
    std::vector<Res> r, rPk, rBu, rD, rBuD;   // rD, rBuD: eta -> 3 pi Dalitz-plot weight on the bump classes
};
std::vector<double> Steps(double a, double b, double d) {
    std::vector<double> v;
    const int n = (int)std::lround((b - a) / d);
    for (int i = 0; i <= n; ++i) v.push_back(std::round((a + i * d) * 1000.) / 1000.);
    return v;
}
std::vector<Fam> MakeFamilies() {
    std::vector<Fam> f;
    auto add = [&](int k, const char* key, const char* tit, std::vector<double> v, bool p0, bool st) {
        Fam x; x.kind = k; x.key = key; x.title = tit; x.v = v; x.pre0 = p0; x.starts = st; f.push_back(x);
    };
    add(kVz, "vz", "v_{z} > [mm]", {kNoCut, -200, -190, -180, -170, -165, -160, -155, -150, -140}, true, true);
    add(kOaLo, "oa_lo", "OA > [deg]", Steps(0, 4, 0.5), false, false);
    add(kOaHi, "oa_hi", "OA < [deg]", {4, 6, 8, 10, 15, 20, 30, kNoCut}, false, false);
    std::vector<double> lab = Steps(25, 70, 5); lab.push_back(kNoCut);
    add(kLab, "oa_lab", "#alpha_{lab} < [deg]", lab, false, false);
    std::vector<double> mlo = {kNoCut}; for (double x : Steps(0.28, 0.38, 0.02)) mlo.push_back(x);
    add(kMppLo, "mpipi_lo", "M_{#pi#pi} > [GeV/c^{2}]", mlo, false, false);
    std::vector<double> mhi = Steps(0.38, 0.50, 0.01); mhi.push_back(kNoCut);
    add(kMppHi, "mpipi_hi", "M_{#pi#pi} < [GeV/c^{2}]", mhi, false, false);
    add(kEta, "eta_rest", "#alpha_{#eta} > [deg]", {kNoCut, 100, 120, 140, 160}, false, false);
    std::vector<double> mm = Steps(2.60, 2.95, 0.05); mm.push_back(kNoCut);
    add(kMm4, "mm4_hi", "MM_{4} < [GeV/c^{2}]", mm, false, false);
    add(kPt4, "pt4_lo", "p_{T,4} > [MeV/c]", Steps(0, 450, 50), false, true);
    add(kDz, "dz_pi", "|z_{#pi} - v_{z}| < [mm], both #pi", {kNoCut, 150, 100, 60, 40, 30, 20}, false, true);
    add(kVtx, "vtx", "tracks in the vertex fit", {0, 1, 2, 3}, false, true);
    return f;
}
std::string ValLab(int k, double v) {
    if (k == kVtx) {
        static const char* n[4] = {"none", "#geq1 #pi", "both #pi", "all 4"};
        const int i = (int)v;
        return (i >= 0 && i < 4) ? n[i] : "?";
    }
    if (v == kNoCut || ((k == kOaLo || k == kPt4) && v == 0)) return "none";
    if (k == kMppLo || k == kMppHi || k == kMm4) return Form("%.2f", v);
    return Form("%g", v);
}
std::string ValTxt(int k, double v) {   // plain text for tables
    if (k == kVtx) {
        static const char* n[4] = {"none", ">=1 pion", "both pions", "all four"};
        const int i = (int)v;
        return (i >= 0 && i < 4) ? n[i] : "?";
    }
    return ValLab(k, v);
}

// ============================================================================
// data tables of the data study (parsed, not recomputed)
// ============================================================================
struct DPt { double v, y, e; };
std::map<std::string, std::vector<DPt>> gD;   // key "family/quantity"
void AddD(const std::string& fam, const std::string& q, double v, double y, double e) { gD[fam + "/" + q].push_back({v, y, e}); }
const std::vector<DPt>* GetD(const std::string& fam, const std::string& q) {
    auto it = gD.find(fam + "/" + q);
    return it == gD.end() ? nullptr : &it->second;
}
bool FindD(const std::string& fam, const std::string& q, double v, DPt& out) {
    const std::vector<DPt>* d = GetD(fam, q);
    if (!d) return false;
    for (const DPt& p : *d)
        if ((p.v == kNoCut && v == kNoCut) || (p.v != kNoCut && v != kNoCut && std::fabs(p.v - v) < 1e-6)) { out = p; return true; }
    return false;
}
std::vector<std::string> Tok(const std::string& s) {
    std::istringstream is(s);
    std::vector<std::string> v;
    std::string t;
    while (is >> t) if (t != "|") v.push_back(t);
    return v;
}
std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}
std::vector<std::string> Split(const std::string& s, char c) {
    std::vector<std::string> v;
    std::string cur;
    for (char x : s) { if (x == c) { v.push_back(cur); cur.clear(); } else cur += x; }
    v.push_back(cur);
    return v;
}
bool IsNumTok(const std::string& s) {
    if (s == "none") return true;
    if (s.empty()) return false;
    const char c = s[0];
    return (c >= '0' && c <= '9') || c == '-' || c == '.';
}
double NumTok(const std::string& s) {
    if (s == "none") return kNoCut;
    std::string u;
    for (char c : s) if (c != '*') u += c;
    return atof(u.c_str());
}
std::vector<double> Nums(const std::string& s) {
    std::vector<double> v;
    for (const std::string& t : Tok(s)) if (IsNumTok(t) && t.find('/') == std::string::npos) v.push_back(NumTok(t));
    return v;
}
int gNDataLines = 0;
void ParsePeakSingle() {
    std::ifstream in(Form("%s/peak_scan_single.txt", kDataScan));
    std::string l;
    const std::map<std::string, std::string> fm = {{"oa_lo", "oa_lo"}, {"oa_hi", "oa_hi"}, {"oa_lab", "oa_lab"},
                                                   {"mpipi_lo", "mpipi_lo"}, {"mpipi_hi", "mpipi_hi"}, {"mm4_hi", "mm4_hi"}};
    struct L { std::string fam; double v, N, cbn, F, sb, eps, err; int flag; };
    std::vector<L> ls;
    while (std::getline(in, l)) {
        const std::vector<std::string> t = Tok(l);
        if (t.size() < 16 || !fm.count(t[0]) || !IsNumTok(t[1])) continue;
        L x; x.fam = t[0]; x.v = NumTok(t[1]); x.flag = atoi(t[2].c_str()); x.N = atof(t[3].c_str());
        x.cbn = atof(t[4].c_str()); x.F = atof(t[7].c_str()); x.sb = atof(t[8].c_str()); x.eps = atof(t[13].c_str()); x.err = atof(t[14].c_str());
        ls.push_back(x);
    }
    double N0 = 0, c0 = 0, F0 = 0;
    for (const L& x : ls) if (x.fam == "oa_lo" && x.v == 0) { N0 = x.N; c0 = x.cbn; F0 = x.F; }
    for (const L& x : ls) {
        const double v = (x.fam == "oa_lo" && x.v == 0) ? 0 : x.v;
        AddD(x.fam, "pkS", v, x.eps, x.err);
        AddD(x.fam, "pkNCB", v, x.N * (1 - x.cbn) / (N0 * (1 - c0)), 0);
        AddD(x.fam, "pkF", v, F0 > 0 ? x.F / F0 : 0, 0);
        AddD(x.fam, "pkFlag", v, x.flag, 0);
        AddD(x.fam, "pkSB", v, x.sb, 0);
        ++gNDataLines;
    }
}
void ParseBumpSingle() {
    std::ifstream in(Form("%s/bump_scan1d.txt", kDataScan));
    std::string l, fam;
    const std::map<std::string, std::string> fm = {{"OA lower", "oa_lo"}, {"OA upper", "oa_hi"}, {"oa_lab upper", "oa_lab"},
                                                   {"M(pipi) upper", "mpipi_hi"}, {"M(pipi) lower", "mpipi_lo"},
                                                   {"mm4 upper", "mm4_hi"}, {"eta-rest lower", "eta_rest"}};
    while (std::getline(in, l)) {
        if (l.find("Plateaus") == 0) break;
        bool hdr = false;
        for (auto& kv : fm)
            if (l.find(kv.first + "  (base") == 0) { fam = kv.second; hdr = true; }
        if (hdr || fam.empty()) continue;
        const std::vector<std::string> t = Tok(l);
        if (t.size() < 26 || !IsNumTok(t[0])) continue;
        double v = NumTok(t[0]);
        if (fam == "oa_lo" && v == kNoCut) v = 0;
        const double effE = atof(t[4].c_str()), errE = atof(t[5].c_str()), effT = atof(t[14].c_str()),
                     errT = atof(t[15].c_str()), effW = atof(t[22].c_str()), effN = atof(t[23].c_str());
        AddD(fam, "buT", v, effT, errT);
        const bool eValid = !(fam == "eta_rest" && v != kNoCut && v >= 120);   // edge fits not valid (bump_scan1d.txt)
        if (eValid) AddD(fam, "buE", v, effE, errE);
        AddD(fam, "buW", v, effW, 0);
        AddD(fam, "buGT", v, effN > 0 ? effT / std::sqrt(effN) : 0, effN > 0 ? errT / std::sqrt(effN) : 0);
        if (eValid) AddD(fam, "buGE", v, effN > 0 ? effE / std::sqrt(effN) : 0, effN > 0 ? errE / std::sqrt(effN) : 0);
        ++gNDataLines;
    }
}
// extra tables: cut text -> fields (numbers of each '|' field)
std::map<std::string, std::vector<std::vector<double>>> ParseExtra(const char* file) {
    std::map<std::string, std::vector<std::vector<double>>> m;
    std::ifstream in(Form("%s/%s", kDataScan, file));
    std::string l;
    while (std::getline(in, l)) {
        if (l.size() < 45 || l[0] == '#' || l[1] != ' ' || l[2] != ' ') continue;
        const size_t p = l.find('|', 30);
        if (p == std::string::npos) continue;
        const std::string cut = Trim(l.substr(3, p - 3));
        std::vector<std::vector<double>> f;
        for (const std::string& s : Split(l.substr(p + 1), '|')) f.push_back(Nums(s));
        m[cut] = f;
    }
    return m;
}
void ParseExtraAll() {
    auto pk = ParseExtra("extra_peak_candidates.txt");
    auto bu = ParseExtra("extra_bump_candidates.txt");
    double Fst = 0;
    if (pk.count("no cut") && pk["no cut"].size() > 3 && !pk["no cut"][3].empty()) Fst = pk["no cut"][3][0];
    struct It { std::string fam, cut; double v; };
    std::vector<It> its;
    its.push_back({"vz", "no cut", kNoCut});
    for (double v : {-200., -190., -180., -170., -165., -160., -155., -150., -140.}) its.push_back({"vz", Form("vz > %g", v), v});
    its.push_back({"pt4_lo", "no cut", 0});
    for (double v : Steps(50, 450, 50)) its.push_back({"pt4_lo", Form("pt4 > %g", v), v});
    its.push_back({"dz_pi", "no cut", kNoCut});
    for (double v : {150., 100., 60., 40., 30., 20.}) its.push_back({"dz_pi", Form("|z - vz| (both pi) < %g", v), v});
    its.push_back({"vtx", "no cut", 0});
    its.push_back({"vtx", "isUsedVertex at least one pion", 1});
    its.push_back({"vtx", "isUsedVertex both pions", 2});
    its.push_back({"vtx", "isUsedVertex all four", 3});
    for (const It& it : its) {
        if (pk.count(it.cut)) {
            const auto& f = pk[it.cut];   // f[0] epsS err, f[3] F dF, f[10] epsTag err Gt
            if (f.size() > 10 && f[0].size() >= 2) AddD(it.fam, "xpkS", it.v, f[0][0], f[0][1]);
            if (f.size() > 3 && !f[3].empty() && Fst > 0) AddD(it.fam, "xpkF", it.v, f[3][0] / Fst, 0);
            ++gNDataLines;
        }
        if (bu.count(it.cut)) {
            const auto& f = bu[it.cut];   // f[0] epsTag err, f[1] epsN sqrtN epsNCB, f[2] Gb err pull
            if (f.size() > 2 && f[0].size() >= 2) AddD(it.fam, "xbuT", it.v, f[0][0], f[0][1]);
            if (f.size() > 2 && f[1].size() >= 3) AddD(it.fam, "xbuNCB", it.v, f[1][2], 0);
            if (f.size() > 2 && f[2].size() >= 2) AddD(it.fam, "xbuG", it.v, f[2][0], f[2][1]);
            ++gNDataLines;
        }
    }
}
// two-dimensional blocks: title line, column-label line, rows "label v v v"
struct Blk { std::vector<double> cols, rows; std::vector<std::vector<double>> v; bool ok = false; };
Blk ParseBlock(const char* file, const char* section, const char* title) {
    Blk b;
    std::ifstream in(Form("%s/%s", kDataScan, file));
    std::string l;
    bool sec = (section == nullptr || !*section);
    while (std::getline(in, l)) {
        if (!sec) { if (l.find(section) != std::string::npos) sec = true; continue; }
        if (l.find(title) == std::string::npos) continue;
        std::getline(in, l);
        for (const std::string& t : Tok(l)) b.cols.push_back(NumTok(t));
        while (std::getline(in, l)) {
            const std::vector<std::string> t = Tok(l);
            if (t.size() != b.cols.size() + 1 || !IsNumTok(t[0])) break;
            b.rows.push_back(NumTok(t[0]));
            std::vector<double> r;
            for (size_t i = 1; i < t.size(); ++i) r.push_back(atof(t[i].c_str()));
            b.v.push_back(r);
        }
        b.ok = !b.rows.empty();
        break;
    }
    return b;
}
// ribbon tables of the data: "central" lines in order
std::vector<std::vector<std::vector<double>>> ParseRibbon(const char* file) {
    std::vector<std::vector<std::vector<double>>> out;
    std::ifstream in(Form("%s/%s", kDataScan, file));
    std::string l;
    while (std::getline(in, l)) {
        if (l.empty() || l[0] == '#') continue;
        const std::vector<std::string> t = Tok(l);
        if (t.empty() || t[0] != "central") continue;
        const std::vector<std::string> f = Split(l, '|');
        std::vector<std::vector<double>> r;
        for (size_t i = 1; i < f.size(); ++i) r.push_back(Nums(f[i]));
        out.push_back(r);
    }
    return out;
}

// ============================================================================
// plotting helpers
// ============================================================================
const int kColData = 601;       // kBlue+1: data signal
const int kColTag = 803;        // kOrange+3: photon-tagged bump (sim)
const int kColB2pk = 923;       // kGray+3
const int kColB2bu = 921;       // kGray+1
struct Ser {
    std::vector<double> x, y, e;
    int col = 1, lsty = 1, lw = 3, mk = 20;
    double ms = 1.3;
    bool line = true, band = true;
    std::string leg;
};
Ser SimSer(const Fam& f, std::function<std::pair<double, double>(const Res&)> q, int col, int lsty, const char* leg,
           int which = 0) {
    Ser s; s.col = col; s.lsty = lsty; s.leg = leg;
    const std::vector<Res>& rr = which == 0 ? f.r : which == 1 ? f.rPk : which == 2 ? f.rBu : f.rD;
    for (size_t i = 0; i < rr.size(); ++i) {
        const std::pair<double, double> v = q(rr[i]);
        s.x.push_back(i + 0.5); s.y.push_back(v.first); s.e.push_back(v.second);
    }
    return s;
}
Ser DatSer(const Fam& f, const char* q, int col, int mk, const char* leg, double ms = 1.5) {
    Ser s; s.col = col; s.mk = mk; s.leg = leg; s.line = false; s.band = false; s.ms = ms;
    for (size_t i = 0; i < f.v.size(); ++i) {
        DPt p;
        if (FindD(f.key, q, f.v[i], p) && p.y > 0) { s.x.push_back(i + 0.5); s.y.push_back(p.y); s.e.push_back(p.e); }
    }
    return s;
}
void DrawSer(const Ser& s, double dx = 0) {
    const int n = (int)s.x.size();
    if (!n) return;
    std::vector<double> x(s.x), ex(n, 0.);
    for (double& v : x) v += dx;
    TGraphErrors* g = new TGraphErrors(n, &x[0], &s.y[0], &ex[0], &s.e[0]);
    if (s.line) {
        if (s.band) {
            TGraphErrors* gb = (TGraphErrors*)g->Clone();
            gb->SetFillColorAlpha(s.col, 0.25); gb->SetLineWidth(0);
            gb->Draw("3");
        }
        g->SetLineColor(s.col); g->SetLineStyle(s.lsty); g->SetLineWidth(s.lw);
        g->SetMarkerColor(s.col); g->SetMarkerStyle(s.mk); g->SetMarkerSize(0.8);
        g->Draw("LP");
    } else {
        g->SetLineColor(s.col); g->SetMarkerColor(s.col); g->SetMarkerStyle(s.mk); g->SetMarkerSize(s.ms);
        g->SetLineWidth(2);
        g->Draw("PZ");
    }
}
TH1D* Frame(const char* name, const Fam& f, double ylo, double yhi, const char* ytit) {
    const int n = (int)f.v.size();
    TH1D* h = new TH1D(name, "", n, 0, n);
    h->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) h->GetXaxis()->SetBinLabel(i + 1, ValLab(f.kind, f.v[i]).c_str());
    h->GetXaxis()->SetLabelSize(n > 12 ? 0.050 : n > 9 ? 0.058 : 0.072);
    h->GetXaxis()->LabelsOption("h");
    h->GetXaxis()->SetTitle(f.title.c_str());
    h->GetXaxis()->SetTitleSize(0.065); h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitle(ytit);
    h->GetYaxis()->SetTitleSize(0.065); h->GetYaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetLabelSize(0.06);
    h->SetMinimum(ylo); h->SetMaximum(yhi);
    h->SetStats(0);
    return h;
}
TPad* GridPad(TCanvas* c, const char* name, int i, int nx, int ny, double ytop) {
    c->cd();
    const int ix = i % nx, iy = i / nx;
    const double w = 1.0 / nx, h = ytop / ny;
    TPad* p = new TPad(name, "", ix * w, std::max(0.0, ytop - (iy + 1) * h), (ix + 1) * w, ytop - iy * h);
    p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.09); p->SetBottomMargin(0.17);
    p->SetGridy(1);
    p->Draw();
    p->cd();
    return p;
}
void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
void Header(TCanvas* c, double ytop, const std::vector<std::string>& lines) {
    c->cd();
    TPad* p = new TPad("hdr", "", 0, ytop, 1, 1);
    p->SetFillStyle(0); p->SetBorderSize(0);
    p->Draw(); p->cd();
    const int n = (int)lines.size();
    const double px = c->GetWh() * (1 - ytop);    // header height in pixels
    for (int i = 0; i < n; ++i) {
        TLatex t; t.SetNDC(); t.SetTextFont(i == 0 ? 63 : 43); t.SetTextAlign(13);
        t.SetTextSize(i == 0 ? 29 : 23);
        t.DrawLatex(0.008, 1.0 - (6 + i * std::min(38.0, (px - 8) / n)) / px, lines[i].c_str());
    }
}
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs(Form("%s/%s.png", kFigDir, name));
    if (pdf) c->SaveAs(Form("%s/%s.pdf", kFigDir, name));
}
std::pair<double, double> Range(const std::vector<Ser>& ss, double lo0, double hi0, double pad = 0.08) {
    double lo = 1e30, hi = -1e30;
    for (const Ser& s : ss)
        for (size_t i = 0; i < s.y.size(); ++i) {
            if (!(s.y[i] > 0)) continue;
            lo = std::min(lo, s.y[i] - s.e[i]); hi = std::max(hi, s.y[i] + s.e[i]);
        }
    if (lo > hi) return std::make_pair(lo0, hi0);
    const double d = hi - lo;
    return std::make_pair(std::min(lo0, lo - pad * d), std::max(hi0, hi + pad * d));
}

}  // namespace scut

// ============================================================================
// main
// ============================================================================
void sim_cuts() {
    using namespace simeng;
    using namespace scut;
    TStopwatch swAll;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gSystem->mkdir(kFigDir, kTRUE);
    gSystem->mkdir(kTabDir, kTRUE);

    // ------------------------------------------------------------------ simulation
    SimData* d = new SimData;
    if (!Load(*d, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;

    // simulated peak position: peak conversion, "pre", mode A, w; chi2 fit with group-variance errors
    TH1D* hPk = nullptr;
    {
        Binning b; b.nbin = 40; b.lo = 0.50; b.hi = 0.60;
        WSpec sp = Fill(*d, SelPreset(kPre), CompPeakConv(), SimOpt(), b);
        hPk = sp.Hist("hPkFit", 0);
        TF1 g("gpk", "gaus", 0.52, 0.57);
        double m = 0.5439, s = 0.0054;
        for (int it = 0; it < 6; ++it) {
            g.SetRange(m - 2 * s, m + 2 * s);
            g.SetParameters(hPk->GetMaximum(), m, s);
            hPk->Fit(&g, "QR0");
            m = g.GetParameter(1); s = std::fabs(g.GetParameter(2));
        }
        gPkMean = m; gPkMeanE = g.GetParError(1); gPkSig = s; gPkSigE = g.GetParError(2);
        gPkChi2 = g.GetChisquare(); gPkNdf = g.GetNDF();
        gPkLo = gPkMean - 1.9 * gPkSig; gPkHi = gPkMean + 1.9 * gPkSig;
        printf("sim peak (conversion): mean %.5f +- %.5f, sigma %.5f +- %.5f, chi2/ndf %.1f/%d; window %.5f - %.5f\n",
               gPkMean, gPkMeanE, gPkSig, gPkSigE, gPkChi2, gPkNdf, gPkLo, gPkHi);
    }
    // spectra for the window figure (engine)
    std::vector<WSpec> wsWin;
    {
        Binning b; b.nbin = 88; b.lo = 0.40; b.hi = 0.62;
        wsWin = FillMany(*d, SelPreset(kPre), {CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompPhysBkg(), CompFake(), CompDiff()},
                         SimOpt(), b);
    }
    // engine cross-check values (compared with this macro's own sums below)
    const Comp cB2bu{"T2 bump window", [](const Row& r) { return r.tmpl == 2 && r.m4 >= kBuLo && r.m4 < kBuHi; }};
    const Comp cSimR{"T1+T2 in R window", [](const Row& r) { return (r.tmpl == 1 || r.tmpl == 2) && r.m4 >= kRLo && r.m4 < kRHi; }};
    const EffRes eng1 = Efficiency(*d, SelPreset(kPre), SelAnd(SelPreset(kPre), [](const Row& r) { return r.oa > 2.f; }), CompPeakConv());
    const EffRes eng2 = Efficiency(*d, SelPreset(kPre), SelAnd(SelPreset(kPre), [](const Row& r) { return r.oa_lab < 50.f; }), cB2bu);
    const EffRes eng3 = Efficiency(*d, SelPreset(kPre),
                                   SelAnd(SelPreset(kPre), [](const Row& r) { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; }),
                                   CompBumpDalitz());
    const std::vector<WSum> eng4 = SumMany(*d, SelPreset(kPre), {cSimR, CompPeakConv()});
    const EffRes eng5 = Efficiency(*d, SelPreset(kPre), SelAnd(SelPreset(kPre), [](const Row& r) { return r.oa > 2.f; }), cB2bu,
                                   [] { SimOpt o; o.mode = 1; o.db = true; return o; }());

    // eta -> pi+pi-pi0 Dalitz-plot weight: normalisation = phase-space average of W(y)
    {
        TGenPhaseSpace ps;
        TLorentzVector eta(0, 0, 0, kMetaG);
        double mm[3] = {kMpiG, kMpiG, kMpi0G};
        ps.SetDecay(eta, 3, mm);
        double sw = 0, sww = 0;
        for (int i = 0; i < 400000; ++i) {
            const double wt = ps.Generate();
            const TLorentzVector pp = *ps.GetDecay(0) + *ps.GetDecay(1);
            sw += wt; sww += wt * DalW(pp.M());
        }
        gDalNorm = sww / sw;
        printf("Dalitz-plot normalisation <W>_PS = %.4f\n", gDalNorm);
    }
    // compact rows
    {
        size_t nIn = 0;
        for (int set = 1; set <= 2; ++set) {
            for (const Row& r : d->tree(set, kUS)) {
                ++nIn;
                if (!(r.isBest >= 0 && r.vz > -500.f && r.cor_ok == 1 && PT3(r))) continue;   // pre0
                const bool tag = r.g_idx > 0 && r.m_eeg >= 0.119f && r.m_eeg < 0.173f;
                const unsigned mask = TgtMask(r.tmpl, r.cls, tag, r.m4);
                if (!mask) continue;
                MRow m;
                m.m4 = r.m4; m.oa = r.oa; m.oa_lab = r.oa_lab; m.mpipi = r.mpipi; m.oa_eta = r.oa_eta; m.mm4 = r.mm4;
                m.pt4 = r.pt4; m.vz = r.vz;
                m.dzpi = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
                m.nvtxpi = (r.pip_isUsedVertex == 1) + (r.pim_isUsedVertex == 1);
                m.vtx4 = (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1 && r.ep_isUsedVertex == 1 && r.em_isUsedVertex == 1);
                m.own = (r.ownA == 1 ? 1 : 0) | (r.ownB == 1 ? 2 : 0);
                m.tmpl = r.tmpl; m.cls = r.cls; m.tag = tag;
                m.mask = mask;
                m.gid = r.grp + (set == 2 ? (1 << 28) : 0);
                m.dal = (r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && r.mpipi_true > 0) ? (float)(DalW(r.mpipi_true) / gDalNorm) : 1.f;
                m.w = r.w; m.wdb = r.w_db;
                gRows.push_back(m);
            }
        }
        printf("compact simulation rows: %zu of %zu\n", gRows.size(), nIn);
    }
    delete d;
    d = nullptr;

    // ------------------------------------------------------------------ data (R window only)
    {
        scan::Data dd;
        if (!scan::Load(dd, scan::kCompact, kRLo - 0.001, kRHi + 0.001)) return;
        for (int k = 0; k < 3; ++k)
            for (const scan::Row& r : dd.s[k]) {
                if (!(r.isBest >= 0 && r.vz > -500.f && r.cor_ok == 1)) continue;
                if (!(r.m4 >= kRLo && r.m4 < kRHi)) continue;
                DRow x;
                x.m4 = r.m4; x.oa = r.oa; x.oa_lab = r.oa_lab; x.mpipi = r.mpipi; x.oa_eta = r.oa_eta; x.mm4 = r.mm4;
                x.pt4 = r.pt4; x.vz = r.vz; x.dzpi = std::numeric_limits<float>::quiet_NaN();
                x.nvtxpi = (r.pip_isUsedVertex == 1) + (r.pim_isUsedVertex == 1);
                x.vtx4 = (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1 && r.ep_isUsedVertex == 1 && r.em_isUsedVertex == 1);
                x.s = (unsigned char)k;
                gData.push_back(x);
            }
        printf("data rows in %.2f-%.2f: %zu\n", kRLo, kRHi, gData.size());
    }
    ParsePeakSingle();
    ParseBumpSingle();
    ParseExtraAll();
    printf("data table lines parsed: %d\n", gNDataLines);

    // ------------------------------------------------------------------ reference "pre" and cross-checks
    const Res rPre = Run(Pre(), Pre());
    const Res rPre0 = Run(Pre0(), Pre0());
    {
        Cuts s1 = Pre(); s1.oaLo = 2;
        Cuts s2 = Pre(); s2.labHi = 50;
        Cuts s3 = Pre(); s3.vtx = 2;
        const std::vector<Acc> a1 = Eval(Pre(), s1), a2 = Eval(Pre(), s2), a3 = Eval(Pre(), s3);
        EOpt ob; ob.mode = 1; ob.db = true;
        const std::vector<Acc> a5 = Eval(Pre(), s1, ob);
        gOut = fopen(Form("%s/cuts_summary.txt", kTabDir), "w");
        P("# cuts_summary.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# simulation: std + lep owned rows, mode A, w (PDG) unless stated; sigma = sum w / L_sim, L_sim = %.2f nb^-1\n", kLsim);
        P("# data: (N - CB) / L_data, L_data = %.1f nb^-1 (UNCONFIRMED), CB = 2 sqrt(N++ N--) in 10 MeV bins, no k factor\n", kLdata);
        P("# eta -> pi+pi-e+e- (direct) is NOT in the simulation (eta set missing)\n#\n");
        P("SIMULATED PEAK (peak conversion, pre, chi2 Gaussian fit, 2.5 MeV bins, iterated range mean +- 2 sigma)\n");
        P("  mean   %.5f +- %.5f GeV/c2   (data 0.5439; difference %+.2f +- %.2f MeV)\n", gPkMean, gPkMeanE,
          1e3 * (gPkMean - 0.5439), 1e3 * gPkMeanE);
        P("  sigma  %.5f +- %.5f GeV/c2   (data 0.0054)\n", gPkSig, gPkSigE);
        P("  chi2/ndf %.1f/%d\n", gPkChi2, gPkNdf);
        P("WINDOWS\n  sim peak window  %.5f - %.5f (mean +- 1.9 sigma, the rule of the data window)\n", gPkLo, gPkHi);
        P("  data peak window %.5f - %.5f (variant)\n  bump window      %.4f - %.4f (data)\n  R window         %.2f - %.2f\n", kPkLoD,
          kPkHiD, kBuLo, kBuHi, kRLo, kRHi);
        P("\nCROSS-CHECK against the engine (simeng::Efficiency / SumMany on the full rows)\n");
        P("  peak conv, OA > 2:        engine %.5f +- %.5f   this macro %.5f +- %.5f\n", eng1.eps, eng1.err, a1[tPC].eps(), a1[tPC].epsErr());
        P("  T2 bump window, lab < 50: engine %.5f +- %.5f   this macro %.5f +- %.5f\n", eng2.eps, eng2.err, a2[tB2bu].eps(), a2[tB2bu].epsErr());
        P("  bump Dalitz, pions vtx:   engine %.5f +- %.5f   this macro %.5f +- %.5f\n", eng3.eps, eng3.err, a3[tBD].eps(), a3[tBD].epsErr());
        P("  T1+T2 in 0.60-0.80, sum w: engine %.4f (V %.4f)   this macro %.4f (V %.4f)\n", eng4[0].S, eng4[0].V, rPre.a[tSimR].Sp, rPre.a[tSimR].aa);
        P("  peak conv, sum w:          engine %.4f (V %.4f)   this macro %.4f (V %.4f)\n", eng4[1].S, eng4[1].V, rPre.a[tPC].Sp, rPre.a[tPC].aa);
        P("  T2 bump, OA > 2, mode B w_db: engine %.5f +- %.5f   this macro %.5f +- %.5f\n", eng5.eps, eng5.err, a5[tB2bu].eps(), a5[tB2bu].epsErr());
        P("\nREFERENCE \"pre\" (mode A, w): sigma [nb] +- stat, rows, N_eff\n");
        for (int t = 0; t < NT; ++t)
            P("  %-42s %9.4f +- %7.4f  rows %8.0f  N_eff %9.1f\n", TgtName(t), rPre.a[t].xs(), rPre.a[t].xsErr(), rPre.a[t].Np, rPre.a[t].NeffP());
        P("\nR = data (N - CB) / simulation (T1 + T2) in %.2f - %.2f GeV/c2 at L_data = %.1f nb^-1 (UNCONFIRMED)\n", kRLo, kRHi, kLdata);
        P("  pre : data %.3f +- %.3f nb (N %.0f, CB %.0f), sim %.3f +- %.3f nb, R = %.4f +- %.4f\n", rPre.dw.sig, rPre.dw.err, rPre.dw.N,
          rPre.dw.CB, rPre.a[tSimR].xs(), rPre.a[tSimR].xsErr(), rPre.R, rPre.Rerr);
        P("  pre0: data %.3f +- %.3f nb (N %.0f, CB %.0f), sim %.3f +- %.3f nb, R = %.4f +- %.4f\n", rPre0.dw.sig, rPre0.dw.err, rPre0.dw.N,
          rPre0.dw.CB, rPre0.a[tSimR].xs(), rPre0.a[tSimR].xsErr(), rPre0.R, rPre0.Rerr);
        for (int w = 0; w < 3; ++w) {
            const FR f = Fom(rPre.a, w), fr = Fom(rPre.a, w, rPre.R, rPre.R);
            P("  %s: S = %.4f +- %.4f nb, B = %.4f +- %.4f nb, S/B = %.4f, F = %.2f (counts at L_data); with R: S/B = %.4f, F = %.2f\n",
              w == 0 ? "peak (sim window) " : w == 1 ? "bump              " : "peak (data window)", f.S, f.Se, f.B, f.Be, f.SB, f.F, fr.SB, fr.F);
        }
        {
            EOpt od; od.dal = true;
            const std::vector<Acc> ad = Eval(Pre(), Pre(), od);
            P("\nETA -> PI+ PI- PI0 DALITZ-PLOT WEIGHT (variant for the bump classes; the simulated decay looks like phase space)\n");
            P("  W = 1 + a y + b y^2 + f y^3, a = %.3f, b = %.3f, f = %.3f (KLOE-2 2016, from memory, UNVERIFIED), x terms dropped,\n", kDalA, kDalB, kDalF);
            P("  y = 3 T_pi0/Q - 1 from mpipi_true; normalised to the phase-space average <W>_PS = %.4f (TGenPhaseSpace, 4e5 decays)\n", gDalNorm);
            P("  pre: bump D+C (all M4) %.4f -> %.4f nb (x %.4f); S_bump %.4f -> %.4f nb; tagged %.4f -> %.4f nb\n", rPre.a[tBDC].xs(), ad[tBDC].xs(),
              ad[tBDC].Sp / rPre.a[tBDC].Sp, rPre.a[tSbu].xs(), ad[tSbu].xs(), rPre.a[tTag].xs(), ad[tTag].xs());
            P("  bump pi0 Dalitz %.4f -> %.4f nb, bump conversion %.4f -> %.4f nb; F(bump) %.2f -> %.2f\n", rPre.a[tBD].xs(), ad[tBD].xs(), rPre.a[tBC].xs(),
              ad[tBC].xs(), Fom(rPre.a, 1).F, Fom(ad, 1).F);
        }
        fclose(gOut); gOut = nullptr;
    }

    // ------------------------------------------------------------------ 1. each cut alone
    std::vector<Fam> fams = MakeFamilies();
    for (Fam& f : fams) {
        for (double v : f.v) {
            Cuts base = f.pre0 ? Pre0() : Pre();
            Cuts sel = base; SetCut(sel, f.kind, v);
            f.r.push_back(Run(base, sel));
            EOpt od; od.dal = true;
            f.rD.push_back(Run(base, sel, od, false));
            if (f.starts) {
                Cuts bp = PeakStart(), sp = bp; SetCut(sp, f.kind, v);
                Cuts bb = BumpStart(), sb = bb; SetCut(sb, f.kind, v);
                f.rPk.push_back(Run(bp, sp, EOpt(), false));
                f.rBu.push_back(Run(bb, sb, EOpt(), false));
                f.rBuD.push_back(Run(bb, sb, od, false));
            }
        }
    }
    printf("1D scans done (%.0f s)\n", swAll.RealTime()); swAll.Continue();

    // table of efficiencies
    {
        const int cols[] = {tBD, tBC, tPC, tTag, tSbu, tSpk, tB2pk, tB2bu, tSimPk, tSimBu, tT3pk, tT3bu};
        const char* cn[] = {"bumpDal", "bumpConv", "peakConv", "bumpTag", "S_bump", "S_peak", "B_peak", "B_bump", "T12_pk", "T12_bu", "fake_pk", "fake_bu"};
        gOut = fopen(Form("%s/cuts_eff_single.txt", kTabDir), "w");
        P("# cuts_eff_single.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# Efficiency of each cut ALONE relative to \"pre\" (isBest >= 0, PT3, vz > -160, z > -500, cor_ok), vz family relative to pre0 (no vz > -160).\n");
        P("# Simulation std + lep, owned rows mode A, weight w (PDG). eps = sum w (pass) / sum w (base); error from event groups with the pass/fail\n");
        P("# correlation of a group. Classes: bumpDal = bump pi0 Dalitz, bumpConv = bump conversion, peakConv = peak conversion (all M4);\n");
        P("# bumpTag = bump D+C with a photon of 0.119 <= M(ee gamma) < 0.173 (emulates the data ECAL tag, UNVERIFIED);\n");
        P("# S_bump = bump D+C in %.4f-%.4f; S_peak = peak conv in the sim peak window %.4f-%.4f; B_* = physical background (T2) in the window;\n", kBuLo, kBuHi, gPkLo, gPkHi);
        P("# T12_* = T1 + T2 (simulated spectrum) in the window; fake_* = template 3 (not part of the simulated spectrum, information only).\n");
        P("# Base N_eff (sum w)^2 / V per class is given in the line 'N_eff(base)'. Last two columns: bump D+C (all M4) and bumpTag with the\n");
        P("# eta -> pi+pi-pi0 Dalitz-plot weight (variant; the simulated decay looks like phase space; see cuts_summary.txt).\n#\n");
        for (const Fam& f : fams) {
            P("\n%-10s %-10s", f.key.c_str(), "value");
            for (const char* c : cn) P(" | %-15s", c);
            P(" | %-15s | %-15s\n", "bumpDC_dalw", "bumpTag_dalw");
            for (size_t i = 0; i < f.v.size(); ++i) {
                P("%-10s %-10s", f.key.c_str(), ValTxt(f.kind, f.v[i]).c_str());
                for (int t : cols) P(" | %6.4f +- %5.4f", f.r[i].a[t].eps(), f.r[i].a[t].epsErr());
                P(" | %6.4f +- %5.4f | %6.4f +- %5.4f\n", f.rD[i].a[tBDC].eps(), f.rD[i].a[tBDC].epsErr(), f.rD[i].a[tTag].eps(), f.rD[i].a[tTag].epsErr());
            }
            P("%-10s %-10s", f.key.c_str(), "N_eff(base)");
            for (int t : cols) P(" | %15.1f", f.r[0].a[t].NeffA());
            P("\n");
            if (f.starts) {
                P("%-10s at the scan_extra start selections (data base): peak start = pre0 + oa_lab < 55 + M(pipi) > 0.30 (peak conv, S_peak, B_peak);\n", "");
                P("%-10s bump start = pre0 + oa_lab < 40 + 0.30 < M(pipi) < 0.41 (bump classes, tagged, B_bump)\n", "");
                for (size_t i = 0; i < f.v.size(); ++i) {
                    P("%-10s %-10s", (f.key + "@st").c_str(), ValTxt(f.kind, f.v[i]).c_str());
                    for (int t : cols) {
                        const bool pk = (t == tPC || t == tSpk || t == tB2pk || t == tSimPk || t == tT3pk);
                        const Acc& a = pk ? f.rPk[i].a[t] : f.rBu[i].a[t];
                        P(" | %6.4f +- %5.4f", a.eps(), a.epsErr());
                    }
                    P("\n");
                }
            }
        }
        fclose(gOut); gOut = nullptr;
    }

    // comparison with the data (all matched points)
    {
        gOut = fopen(Form("%s/cuts_compare_data.txt", kTabDir), "w");
        P("# cuts_compare_data.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# Simulated efficiency (truth classes, mode A, w) against the data-derived efficiency of the data study, same cut, same base.\n");
        P("# Data: peak eps_S = S_pass/(S_pass+S_fail) from fits of the data peak (it CONTAINS eta -> pi+pi-e+e- direct, missing in the simulation);\n");
        P("#       bump eff_T = ECAL-tagged eta (m5c, pi0 window 0.119-0.173), eff_E = edge-step sum; tables in %s.\n", kDataScan);
        P("# Sim classes: peak -> peak conversion (all M4); tagged -> bump D+C with photon 0.119 <= M(ee gamma) < 0.173; bump all -> bump D+C (all M4);\n");
        P("#              edge E -> S_bump (bump D+C in the bump window).  pull = (sim - data) / sqrt(err_sim^2 + err_data^2)\n");
        P("# Last column: sim tagged with the eta -> pi+pi-pi0 Dalitz-plot weight (variant, see cuts_summary.txt).\n");
        P("# Bases: 1D scans of scan_peak.C / scan_bump.C = no cut (sim: pre, i.e. with vz > -160, which keeps 99.6 %% of the signal);\n");
        P("#        scan_extra.C (txt) = start selections (sim: same start, pre0 base)\n#\n");
        P("%-9s %-11s %-8s | %-16s %-16s %6s | %-16s %-16s %6s | %-16s %-16s %6s | %-16s | %-16s\n", "family", "value", "base", "data peak", "sim peak", "pull",
          "data tagged", "sim tagged", "pull", "data edge E", "sim S_bump", "pull", "sim bump all", "sim tagged Dal.w");
        auto pull = [](double a, double ea, double b, double eb) { const double e = std::sqrt(ea * ea + eb * eb); return e > 0 ? (a - b) / e : 0.; };
        for (const Fam& f : fams) {
            for (size_t i = 0; i < f.v.size(); ++i) {
                DPt ps, pt, pe;
                const bool hs = FindD(f.key, "pkS", f.v[i], ps), ht = FindD(f.key, "buT", f.v[i], pt), he = FindD(f.key, "buE", f.v[i], pe);
                if (!hs && !ht && !he) continue;
                const std::vector<Acc>& a = f.r[i].a;
                Acc bu; bu.Sa = a[tBD].Sa + a[tBC].Sa; bu.Sp = a[tBD].Sp + a[tBC].Sp;   // approximate error below
                const double buErr = std::sqrt(std::pow(a[tBD].epsErr() * a[tBD].Sa, 2) + std::pow(a[tBC].epsErr() * a[tBC].Sa, 2)) / bu.Sa;
                P("%-9s %-11s %-8s |", f.key.c_str(), ValTxt(f.kind, f.v[i]).c_str(), "no cut");
                if (hs) P(" %6.3f +- %5.3f   %6.3f +- %5.3f  %+5.1f |", ps.y, ps.e, a[tPC].eps(), a[tPC].epsErr(), pull(a[tPC].eps(), a[tPC].epsErr(), ps.y, ps.e));
                else P(" %-16s %-16s %6s |", "-", "-", "");
                if (ht) P(" %6.3f +- %5.3f   %6.3f +- %5.3f  %+5.1f |", pt.y, pt.e, a[tTag].eps(), a[tTag].epsErr(), pull(a[tTag].eps(), a[tTag].epsErr(), pt.y, pt.e));
                else P(" %-16s %-16s %6s |", "-", "-", "");
                if (he) P(" %6.3f +- %5.3f   %6.3f +- %5.3f  %+5.1f |", pe.y, pe.e, a[tSbu].eps(), a[tSbu].epsErr(), pull(a[tSbu].eps(), a[tSbu].epsErr(), pe.y, pe.e));
                else P(" %-16s %-16s %6s |", "-", "-", "");
                P(" %6.3f +- %5.3f | %6.3f +- %5.3f\n", bu.Sp / bu.Sa, buErr, f.rD[i].a[tTag].eps(), f.rD[i].a[tTag].epsErr());
            }
            if (!f.starts) continue;
            for (size_t i = 0; i < f.v.size(); ++i) {
                DPt ps, pt;
                const bool hs = FindD(f.key, "xpkS", f.v[i], ps), ht = FindD(f.key, "xbuT", f.v[i], pt);
                if (!hs && !ht) continue;
                const Acc& ap = f.rPk[i].a[tPC];
                const Acc& at = f.rBu[i].a[tTag];
                const std::vector<Acc>& ab = f.rBu[i].a;
                const double buEps = (ab[tBD].Sp + ab[tBC].Sp) / (ab[tBD].Sa + ab[tBC].Sa);
                const double buErr = std::sqrt(std::pow(ab[tBD].epsErr() * ab[tBD].Sa, 2) + std::pow(ab[tBC].epsErr() * ab[tBC].Sa, 2)) /
                                     (ab[tBD].Sa + ab[tBC].Sa);
                P("%-9s %-11s %-8s |", f.key.c_str(), ValTxt(f.kind, f.v[i]).c_str(), "start");
                if (hs) P(" %6.3f +- %5.3f   %6.3f +- %5.3f  %+5.1f |", ps.y, ps.e, ap.eps(), ap.epsErr(), pull(ap.eps(), ap.epsErr(), ps.y, ps.e));
                else P(" %-16s %-16s %6s |", "-", "-", "");
                if (ht) P(" %6.3f +- %5.3f   %6.3f +- %5.3f  %+5.1f |", pt.y, pt.e, at.eps(), at.epsErr(), pull(at.eps(), at.epsErr(), pt.y, pt.e));
                else P(" %-16s %-16s %6s |", "-", "-", "");
                P(" %-16s %-16s %6s |", "-", "-", "");
                P(" %6.3f +- %5.3f | %6.3f +- %5.3f\n", buEps, buErr, f.rBuD[i].a[tTag].eps(), f.rBuD[i].a[tTag].epsErr());
            }
        }
        // the numbers quoted in scan_extra.md (start: peak oa_lab < 50, M(pipi) < 0.46; bump oa_lab < 50, M(pipi) < 0.42)
        P("\n# scan_extra.md start selections (peak: pre0 + oa_lab < 50 + M(pipi) < 0.46; bump: pre0 + oa_lab < 50 + M(pipi) < 0.42)\n");
        struct Q { const char* name; int kind; double v; double dpk, dpkE, dtg, dtgE; };
        const Q qs[] = {{"both pions in vertex fit", kVtx, 2, 0.912, 0.045, 0.945, 0.014},
                        {"all four in vertex fit", kVtx, 3, 0.771, 0.049, 0.758, 0.016},
                        {"vz > -160", kVz, -160, 1.008, 0.029, 0.999, 0.007},
                        {"|z_pi - vz| < 40 (peak) / 60 (tagged)", kDz, 40, 0.974, 0.038, 1.001, 0.009}};
        for (const Q& q : qs) {
            Cuts bp = PeakStartMD(), sp = bp; SetCut(sp, q.kind, q.v);
            Cuts bb = BumpStartMD(), sb = bb; SetCut(sb, q.kind, q.kind == kDz ? 60 : q.v);
            const std::vector<Acc> ap = Eval(bp, sp), ab = Eval(bb, sb);
            P("%-40s | data peak %5.3f +- %5.3f  sim peak conv %5.3f +- %5.3f  pull %+5.1f | data tagged %5.3f +- %5.3f  sim tagged %5.3f +- %5.3f  pull %+5.1f | sim bump all %5.3f\n",
              q.name, q.dpk, q.dpkE, ap[tPC].eps(), ap[tPC].epsErr(), pull(ap[tPC].eps(), ap[tPC].epsErr(), q.dpk, q.dpkE), q.dtg, q.dtgE,
              ab[tTag].eps(), ab[tTag].epsErr(), pull(ab[tTag].eps(), ab[tTag].epsErr(), q.dtg, q.dtgE),
              (ab[tBD].Sp + ab[tBC].Sp) / (ab[tBD].Sa + ab[tBC].Sa));
        }
        fclose(gOut); gOut = nullptr;
    }

    // ------------------------------------------------------------------ figures: efficiencies
    const std::vector<std::string> noteSim = {
        Form("SMASH std + lep, owned rows (mode A), weight w (PDG); event-group errors; %s", "no efficiency correction"),
        Form("%s; %s", kLdataNote, NoteNoDirect())};
    {
        TCanvas* c = new TCanvas("cEffS", "", 2400, 1750);
        const double yt = 0.885;
        Header(c, yt, {"Simulation truth: efficiency of each cut alone for the #eta signal classes, compared with the data study",
                       "Base: pre (isBest #geq 0, PT3, v_{z} > -160 mm, cor_ok; v_{z} panel: without v_{z} > -160). Lines: simulation (band = stat. error). "
                       "Markers: data (blue) from scan_peak / scan_bump (base: no cut) and scan_extra (open: its start selections).",
                       "Data peak #varepsilon_{S}: pass/fail fits of the data peak, which contains #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct), "
                       "absent from the simulation. Dashed: simulation at the scan_extra start selections (like-for-like with open markers).",
                       "Long dash: bump with the #eta#rightarrow#pi^{+}#pi^{-}#pi^{0} Dalitz-plot weight (the simulated decay looks like phase space). "
                       "SMASH std + lep, owned rows (mode A), w (PDG); event-group errors.",
                       noteSim[1]});
        std::vector<TPad*> pads;
        for (size_t k = 0; k < fams.size(); ++k) {
            const Fam& f = fams[k];
            TPad* p = GridPad(c, Form("pS%zu", k), (int)k, 4, 3, yt);
            auto ef = [](int t) { return [t](const Res& r) { return std::make_pair(r.a[t].eps(), r.a[t].epsErr()); }; };
            std::vector<Ser> ss;
            ss.push_back(SimSer(f, ef(tBD), kColBumpDalitz, 1, "sim: bump, #pi^{0} Dalitz"));
            ss.push_back(SimSer(f, ef(tBC), kColBumpConv, 1, "sim: bump, conversion"));
            ss.push_back(SimSer(f, ef(tPC), kColPeakConv, 1, "sim: peak, conversion"));
            ss.push_back(SimSer(f, ef(tTag), kColTag, 3, "sim: bump D+C, photon-tagged"));
            { Ser a = SimSer(f, ef(tBDC), kColBumpDalitz, 9, "", 3); a.band = false; a.lw = 3; ss.push_back(a); }
            { Ser a = SimSer(f, ef(tTag), kColTag, 7, "", 3); a.band = false; a.lw = 3; ss.push_back(a); }
            if (f.starts) {
                Ser a = SimSer(f, ef(tPC), kColPeakConv, 2, "", 1); a.band = false; a.lw = 2; ss.push_back(a);
                Ser b = SimSer(f, ef(tTag), kColTag, 2, "", 2); b.band = false; b.lw = 2; ss.push_back(b);
            }
            ss.push_back(DatSer(f, "pkS", kColData, 20, "data: peak #varepsilon_{S}"));
            ss.push_back(DatSer(f, "buT", kColData, 21, "data: tagged #eta (bump)"));
            ss.push_back(DatSer(f, "buE", kColData, 26, "data: bump edge step E"));
            ss.push_back(DatSer(f, "xpkS", kColData, 24, "data: peak, scan_extra start"));
            ss.push_back(DatSer(f, "xbuT", kColData, 25, "data: tagged #eta, scan_extra start"));
            const std::pair<double, double> yr = Range(ss, 0.0, 1.05, 0.05);
            TH1D* fr = Frame(Form("frS%zu", k), f, std::max(0.0, yr.first), std::min(1.25, yr.second), "efficiency");
            fr->Draw("AXIS");
            TLine* l1 = new TLine(0, 1, f.v.size(), 1); l1->SetLineStyle(3); l1->Draw();
            for (const Ser& s : ss) DrawSer(s, s.line ? 0 : 0.0);
            pads.push_back(p);
        }
        // legend panel
        TPad* pl = GridPad(c, "pSleg", 11, 4, 3, yt);
        pl->SetGridy(0);
        TLegend* lg = new TLegend(0.02, 0.02, 0.98, 0.98);
        lg->SetTextSize(0.052); lg->SetBorderSize(0); lg->SetFillStyle(0);
        auto addL = [&](int col, int ls, const char* t) { TGraph* g = new TGraph(); g->SetLineColor(col); g->SetLineStyle(ls); g->SetLineWidth(3); lg->AddEntry(g, t, "l"); };
        auto addM = [&](int col, int mk, const char* t) { TGraph* g = new TGraph(); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.5); lg->AddEntry(g, t, "p"); };
        addL(kColBumpDalitz, 1, "sim: bump, #pi^{0} Dalitz (all M_{4})");
        addL(kColBumpConv, 1, "sim: bump, conversion (all M_{4})");
        addL(kColPeakConv, 1, "sim: peak, conversion (all M_{4})");
        addL(kColTag, 3, "sim: bump D+C, photon in #pi^{0} window");
        addL(kColPeakConv, 2, "sim at scan_extra start (peak / tagged)");
        addL(kColBumpDalitz, 9, "sim: bump D+C, #eta#rightarrow3#pi Dalitz-plot weight");
        addL(kColTag, 7, "sim: tagged, Dalitz-plot weight");
        addM(kColData, 20, "data: peak #varepsilon_{S} (fit, no-cut base)");
        addM(kColData, 21, "data: ECAL-tagged #eta, #varepsilon_{T}");
        addM(kColData, 26, "data: bump edge step, #varepsilon_{E}");
        addM(kColData, 24, "data: peak #varepsilon_{S}, scan_extra start");
        addM(kColData, 25, "data: tagged #eta, scan_extra start");
        lg->Draw();
        Save(c, "cuts_eff_signal", true);
    }
    {
        TCanvas* c = new TCanvas("cEffB", "", 2400, 1700);
        const double yt = 0.90;
        Header(c, yt, {"Simulation truth: efficiency of each cut alone for the physical background in the two windows, compared with data N - CB",
                       Form("Base: pre (v_{z} panel: pre without v_{z} > -160). Peak window (sim) %.4f - %.4f, bump window %.4f - %.4f GeV/c^{2}. "
                            "Physical background = template 2 (correlated true e^{+}e^{-} pair), no CB in simulation.", gPkLo, gPkHi, kBuLo, kBuHi),
                       "Markers: data N - CB kept in the window (black: peak window 0.5336 - 0.5542, open: bump window), CB = 2#sqrt{N_{++}N_{--}}; "
                       "it contains the signal and fake-lepton pairs. Teal: simulated spectrum T1 + T2 in the same windows (like-for-like).",
                       noteSim[0]});
        for (size_t k = 0; k < fams.size(); ++k) {
            const Fam& f = fams[k];
            GridPad(c, Form("pB%zu", k), (int)k, 4, 3, yt);
            auto ef = [](int t) { return [t](const Res& r) { return std::make_pair(r.a[t].eps(), r.a[t].epsErr()); }; };
            std::vector<Ser> ss;
            ss.push_back(SimSer(f, ef(tB2pk), kColB2pk, 1, ""));
            ss.push_back(SimSer(f, ef(tB2bu), kColB2bu, 2, ""));
            Ser a = SimSer(f, ef(tSimPkD), kColSim, 1, ""); a.band = false; a.lw = 2; ss.push_back(a);
            Ser b = SimSer(f, ef(tSimBu), kColSim, 2, ""); b.band = false; b.lw = 2; ss.push_back(b);
            Ser e = SimSer(f, ef(tT3pk), kColFake, 3, ""); e.band = false; e.lw = 2; ss.push_back(e);
            ss.push_back(DatSer(f, "pkNCB", 1, 20, ""));
            ss.push_back(DatSer(f, "buW", 1, 24, ""));
            ss.push_back(DatSer(f, "xbuNCB", 1, 25, ""));
            const std::pair<double, double> yr = Range(ss, 0.0, 1.05, 0.05);
            TH1D* fr = Frame(Form("frB%zu", k), f, std::max(0.0, yr.first), std::min(1.25, yr.second), "efficiency");
            fr->Draw("AXIS");
            TLine* l1 = new TLine(0, 1, f.v.size(), 1); l1->SetLineStyle(3); l1->Draw();
            for (const Ser& s : ss) DrawSer(s);
        }
        TPad* pl = GridPad(c, "pBleg", 11, 4, 3, yt);
        pl->SetGridy(0);
        TLegend* lg = new TLegend(0.02, 0.02, 0.98, 0.98);
        lg->SetTextSize(0.052); lg->SetBorderSize(0); lg->SetFillStyle(0);
        auto addL = [&](int col, int ls, const char* t) { TGraph* g = new TGraph(); g->SetLineColor(col); g->SetLineStyle(ls); g->SetLineWidth(3); lg->AddEntry(g, t, "l"); };
        auto addM = [&](int col, int mk, const char* t) { TGraph* g = new TGraph(); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.5); lg->AddEntry(g, t, "p"); };
        addL(kColB2pk, 1, "sim: phys. bkg, peak window (sim)");
        addL(kColB2bu, 2, "sim: phys. bkg, bump window");
        addL(kColSim, 1, "sim: T1+T2, data peak window");
        addL(kColSim, 2, "sim: T1+T2, bump window");
        addL(kColFake, 3, "sim: fake lepton, peak window (info)");
        addM(1, 20, "data: N - CB kept, peak window");
        addM(1, 24, "data: N - CB kept, bump window");
        addM(1, 25, "data: N - CB, bump window, scan_extra start");
        lg->Draw();
        Save(c, "cuts_eff_background", true);
    }

    // ------------------------------------------------------------------ 2. FOM per cut value
    // per family and window: best value, region within 5 % of the best, region F >= F_best - 1 (counts at L_data)
    struct Opt1 { int best = -1, bestR = -1; std::vector<int> reg5, reg1, reg5R; };
    std::map<std::string, Opt1> opt1[2];
    {
        gOut = fopen(Form("%s/cuts_fom_1d.txt", kTabDir), "w");
        P("# cuts_fom_1d.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# Each cut alone on \"pre\" (vz family: on pre0). S = signal, B = physical background (T2) in the window, sigma [nb] (sum w / L_sim).\n");
        P("# F = S / sqrt(S + B) in expected counts at L_data = %.1f nb^-1 (UNCONFIRMED; the optimum position does not depend on it).\n", kLdata);
        P("# F/F0 relative to the base; error from the efficiency errors of S and B (event groups), S-B correlation neglected.\n");
        P("# R = data (N - CB) / sim (T1 + T2) in %.2f-%.2f at the same selection (n/a: no pion z in the data file -> R of the base).\n", kRLo, kRHi);
        P("# F'/F'0: the same with B -> R * B (R of the selection; R error not propagated).\n");
        P("# Peak: S = peak conversion in the sim peak window %.4f-%.4f (direct decay missing!). Bump: S = bump D+C in %.4f-%.4f.\n", gPkLo, gPkHi, kBuLo, kBuHi);
        P("# Optimum: best F/F0; 'within 5%%' = values with F >= 0.95 F_best; 'F_best-1' = values with F >= F_best - 1 (counts at L_data).\n");
        for (int w = 0; w < 2; ++w) {
            P("\n==================== %s window ====================\n", w == 0 ? "PEAK" : "BUMP");
            for (const Fam& f : fams) {
                P("\n%-9s %-10s | %9s %8s %9s %8s %8s %7s | %7s %6s | %6s %6s | %7s %6s %8s\n", f.key.c_str(), "value", "S [nb]", "err", "B [nb]",
                  "err", "S/B", "F", "F/F0", "err", "R", "err", "F'/F'0", "err", "S/B'");
                Opt1 o;
                std::vector<FR> fs, frs;
                for (size_t i = 0; i < f.v.size(); ++i) {
                    const Res& r = f.r[i];
                    fs.push_back(Fom(r.a, w));
                    frs.push_back(Fom(r.a, w, r.R, r.R0));
                }
                double best = -1, bestR = -1;
                for (size_t i = 0; i < f.v.size(); ++i) {
                    if (fs[i].F > best) { best = fs[i].F; o.best = (int)i; }
                    if (frs[i].F > bestR) { bestR = frs[i].F; o.bestR = (int)i; }
                }
                for (size_t i = 0; i < f.v.size(); ++i) {
                    if (fs[i].F >= 0.95 * best) o.reg5.push_back((int)i);
                    if (fs[i].F >= best - 1) o.reg1.push_back((int)i);
                    if (frs[i].F >= 0.95 * bestR) o.reg5R.push_back((int)i);
                }
                for (size_t i = 0; i < f.v.size(); ++i) {
                    const Res& r = f.r[i];
                    P("%-9s %-10s | %9.4f %8.4f %9.4f %8.4f %8.4f %7.2f | %7.4f %6.4f | %6.3f %6.3f%s | %7.4f %6.4f %8.4f %s%s\n", f.key.c_str(),
                      ValTxt(f.kind, f.v[i]).c_str(), fs[i].S, fs[i].Se, fs[i].B, fs[i].Be, fs[i].SB, fs[i].F, fs[i].rel, fs[i].relErr, r.R, r.Rerr,
                      r.Rok ? " " : "*", frs[i].rel, frs[i].relErr, frs[i].SB, (int)i == o.best ? "<best" : "", (int)i == o.bestR ? " <best'" : "");
                }
                auto lst = [&](const std::vector<int>& v) { std::string s; for (int i : v) s += ValTxt(f.kind, f.v[i]) + " "; return s; };
                P("   optimum: best %s (F/F0 = %.4f +- %.4f); within 5%%: %s; F >= F_best-1: %s\n", ValTxt(f.kind, f.v[o.best]).c_str(), fs[o.best].rel,
                  fs[o.best].relErr, lst(o.reg5).c_str(), lst(o.reg1).c_str());
                P("   with B*R: best %s (F'/F'0 = %.4f); within 5%%: %s\n", ValTxt(f.kind, f.v[o.bestR]).c_str(), frs[o.bestR].rel, lst(o.reg5R).c_str());
                opt1[w][f.key] = o;
            }
        }
        P("\n==================== BUMP window, eta -> 3 pi Dalitz-plot weight on the bump classes (variant) ====================\n");
        for (const Fam& f : fams) {
            std::vector<FR> fs;
            double best = -1; int ib = 0;
            for (size_t i = 0; i < f.v.size(); ++i) { fs.push_back(Fom(f.rD[i].a, 1)); if (fs.back().F > best) { best = fs.back().F; ib = (int)i; } }
            P("%-9s F/F0:", f.key.c_str());
            for (size_t i = 0; i < f.v.size(); ++i) P(" %s:%.4f+-%.4f", ValTxt(f.kind, f.v[i]).c_str(), fs[i].rel, fs[i].relErr);
            std::string reg;
            for (size_t i = 0; i < f.v.size(); ++i) if (fs[i].F >= 0.95 * best) reg += ValTxt(f.kind, f.v[i]) + " ";
            P("\n   best %s (S = %.4f nb, S/B = %.4f); within 5%%: %s\n", ValTxt(f.kind, f.v[ib]).c_str(), fs[ib].S, fs[ib].SB, reg.c_str());
        }
        fclose(gOut); gOut = nullptr;
    }
    for (int w = 0; w < 2; ++w) {
        TCanvas* c = new TCanvas(Form("cF%d", w), "", 2400, 1700);
        const double yt = 0.90;
        const int col = w == 0 ? kColPeakConv : kColBumpDalitz;
        Header(c, yt, {w == 0 ? "Simulation truth: figure of merit of each cut alone, PEAK window (S = peak conversion; #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} direct missing)"
                              : "Simulation truth: figure of merit of each cut alone, BUMP window (S = bump #pi^{0} Dalitz + conversion)",
                       Form("F = S/#sqrt{S+B} relative to pre (B = physical background in the window: %s). Solid: simulation; dashed: B #times R, "
                            "R = data (N-CB)/sim (T1+T2) in 0.60 - 0.80 at the same selection; long dash (bump only): S with the #eta#rightarrow3#pi Dalitz-plot weight.",
                            (w == 0 ? TString::Format("%.4f - %.4f", gPkLo, gPkHi) : TString::Format("%.4f - %.4f", kBuLo, kBuHi)).Data()),
                       w == 0 ? "Blue: data F_{stat}/F_{stat,0} (scan_peak, no-cut base; open: scan_extra start). Shaded: simulation optimum region "
                                "(F within 5 % of the best). Data plateau of the peak study: #alpha_{lab} 35 - none, M_{#pi#pi} > none - 0.34, OA > 0 - 1.5, OA < 8 - none."
                              : "Blue: data gain #varepsilon_{T}/#sqrt{#varepsilon_{N}} (tagged #eta, full) and #varepsilon_{E}/#sqrt{#varepsilon_{N}} (edge step, open triangles); "
                                "open squares: scan_extra G_{b}. Shaded: simulation optimum (within 5 % of the best).",
                       std::string(noteSim[1]) + "; F in counts at L_{data} (its value does not move the optimum)"});
        for (size_t k = 0; k < fams.size(); ++k) {
            const Fam& f = fams[k];
            GridPad(c, Form("pF%d%zu", w, k), (int)k, 4, 3, yt);
            Ser s1, s2, s3;
            s1.col = col; s1.leg = "sim"; s2.col = col; s2.lsty = 2; s2.band = false; s2.lw = 2; s3.col = kColTag; s3.lsty = 9; s3.band = false; s3.lw = 3;
            for (size_t i = 0; i < f.v.size(); ++i) {
                const FR a = Fom(f.r[i].a, w), b = Fom(f.r[i].a, w, f.r[i].R, f.r[i].R0), d3 = Fom(f.rD[i].a, w);
                s1.x.push_back(i + 0.5); s1.y.push_back(a.rel); s1.e.push_back(a.relErr);
                s2.x.push_back(i + 0.5); s2.y.push_back(b.rel); s2.e.push_back(0);
                s3.x.push_back(i + 0.5); s3.y.push_back(d3.rel); s3.e.push_back(0);
            }
            std::vector<Ser> ss = {s1, s2};
            if (w == 1) ss.push_back(s3);
            if (w == 0) { ss.push_back(DatSer(f, "pkF", kColData, 20, "")); ss.push_back(DatSer(f, "xpkF", kColData, 24, "")); }
            else { ss.push_back(DatSer(f, "buGT", kColData, 21, "")); ss.push_back(DatSer(f, "buGE", kColData, 26, "")); ss.push_back(DatSer(f, "xbuG", kColData, 25, "")); }
            std::pair<double, double> yr = Range(ss, 0.9, 1.1, 0.06);
            yr.first = std::max(yr.first, 0.3); yr.second = std::min(yr.second, 1.6);
            TH1D* fr = Frame(Form("frF%d%zu", w, k), f, yr.first, yr.second, "F / F_{pre}");
            fr->Draw("AXIS");
            const Opt1& o = opt1[w][f.key];
            for (int i : o.reg5) { TBox* bx = new TBox(i + 0.02, yr.first, i + 0.98, yr.second); bx->SetFillColorAlpha(col, 0.12); bx->Draw(); }
            TLine* l1 = new TLine(0, 1, f.v.size(), 1); l1->SetLineStyle(3); l1->Draw();
            for (const Ser& s : ss) DrawSer(s);
            fr->Draw("AXIS SAME");
        }
        TPad* pl = GridPad(c, Form("pFleg%d", w), 11, 4, 3, yt);
        pl->SetGridy(0);
        TLegend* lg = new TLegend(0.02, 0.02, 0.98, 0.98);
        lg->SetTextSize(0.052); lg->SetBorderSize(0); lg->SetFillStyle(0);
        auto addL = [&](int cl, int ls, const char* t) { TGraph* g = new TGraph(); g->SetLineColor(cl); g->SetLineStyle(ls); g->SetLineWidth(3); lg->AddEntry(g, t, "l"); };
        auto addM = [&](int cl, int mk, const char* t) { TGraph* g = new TGraph(); g->SetMarkerColor(cl); g->SetMarkerStyle(mk); g->SetMarkerSize(1.5); lg->AddEntry(g, t, "p"); };
        addL(col, 1, "sim: F/F_{pre}, B = phys. bkg");
        addL(col, 2, "sim: F'/F'_{pre}, B #times R");
        if (w == 1) addL(kColTag, 9, "sim: S with #eta#rightarrow3#pi Dalitz-plot weight");
        if (w == 0) { addM(kColData, 20, "data: F_{stat}/F_{stat}(no cut)"); addM(kColData, 24, "data: F_{stat}, scan_extra start"); }
        else { addM(kColData, 21, "data: #varepsilon_{T}/#sqrt{#varepsilon_{N}}"); addM(kColData, 26, "data: #varepsilon_{E}/#sqrt{#varepsilon_{N}}"); addM(kColData, 25, "data: G_{b}, scan_extra start"); }
        TBox* bx = new TBox(); bx->SetFillColorAlpha(col, 0.12); lg->AddEntry(bx, "sim optimum: F #geq 0.95 F_{best}", "f");
        lg->Draw();
        Save(c, w == 0 ? "cuts_fom_peak" : "cuts_fom_bump", true);
    }
    // S/B against each cut
    {
        TCanvas* c = new TCanvas("cSB", "", 2400, 1700);
        const double yt = 0.90;
        Header(c, yt, {"Simulation truth: S/B of each cut alone in the two windows",
                       Form("B = physical background (template 2; no CB in the simulation); dashed: B #times R (R = data (N-CB)/sim (T1+T2) in 0.60 - 0.80). "
                            "Peak window (sim) %.4f - %.4f, bump window %.4f - %.4f GeV/c^{2}.", gPkLo, gPkHi, kBuLo, kBuHi),
                       "Blue: data S/B of the peak study (S in the window over polynomial + CB: its B includes the CB and the fake-lepton pairs, so it is lower by construction). "
                       "Long dash: bump with the #eta#rightarrow3#pi Dalitz-plot weight.",
                       noteSim[1]});
        for (size_t k = 0; k < fams.size(); ++k) {
            const Fam& f = fams[k];
            TPad* p = GridPad(c, Form("pSB%zu", k), (int)k, 4, 3, yt);
            p->SetLogy(1);
            std::vector<Ser> ss;
            for (int w = 0; w < 2; ++w) {
                Ser a, b, d3;
                a.col = b.col = (w == 0 ? kColPeakConv : kColBumpDalitz); b.lsty = 2; b.band = false; b.lw = 2;
                d3.col = kColTag; d3.lsty = 9; d3.band = false; d3.lw = 3;
                for (size_t i = 0; i < f.v.size(); ++i) {
                    const FR x = Fom(f.r[i].a, w), y = Fom(f.r[i].a, w, f.r[i].R, f.r[i].R0), z = Fom(f.rD[i].a, w);
                    const double e = (x.S > 0 && x.B > 0) ? x.SB * std::sqrt(std::pow(x.Se / x.S, 2) + std::pow(x.Be / x.B, 2)) : 0;
                    a.x.push_back(i + 0.5); a.y.push_back(x.SB); a.e.push_back(e);
                    b.x.push_back(i + 0.5); b.y.push_back(y.SB); b.e.push_back(0);
                    d3.x.push_back(i + 0.5); d3.y.push_back(z.SB); d3.e.push_back(0);
                }
                ss.push_back(a); ss.push_back(b);
                if (w == 1) ss.push_back(d3);
            }
            ss.push_back(DatSer(f, "pkSB", kColData, 20, ""));
            TH1D* fr = Frame(Form("frSB%zu", k), f, 0.03, 3.0, "S/B");
            fr->Draw("AXIS");
            for (const Ser& s : ss) DrawSer(s);
        }
        TPad* pl = GridPad(c, "pSBleg", 11, 4, 3, yt);
        pl->SetGridy(0);
        TLegend* lg = new TLegend(0.02, 0.02, 0.98, 0.98);
        lg->SetTextSize(0.052); lg->SetBorderSize(0); lg->SetFillStyle(0);
        auto addL = [&](int cl, int ls, const char* t) { TGraph* g = new TGraph(); g->SetLineColor(cl); g->SetLineStyle(ls); g->SetLineWidth(3); lg->AddEntry(g, t, "l"); };
        TGraph* gm = new TGraph(); gm->SetMarkerColor(kColData); gm->SetMarkerStyle(20); gm->SetMarkerSize(1.5);
        addL(kColPeakConv, 1, "peak window: S = peak conversion");
        addL(kColPeakConv, 2, "peak window, B #times R");
        addL(kColBumpDalitz, 1, "bump window: S = bump D+C");
        addL(kColBumpDalitz, 2, "bump window, B #times R");
        addL(kColTag, 9, "bump, #eta#rightarrow3#pi Dalitz-plot weight");
        lg->AddEntry(gm, "data peak S/B (B incl. CB, fakes)", "p");
        lg->Draw();
        Save(c, "cuts_sb", true);
    }
    // R against each cut
    {
        TCanvas* c = new TCanvas("cR", "", 2400, 1700);
        const double yt = 0.92;
        Header(c, yt, {Form("R = data (N - CB) / simulation (T1 + T2) in %.2f - %.2f GeV/c^{2}, for each cut alone (base pre; v_{z} panel: pre0)", kRLo, kRHi),
                       Form("CB = 2#sqrt{N_{++}N_{--}} in 10 MeV bins, no k factor; simulation std + lep, mode A, w; %s. No point: pion z not in the data file.", kLdataNote)});
        for (size_t k = 0; k < fams.size(); ++k) {
            const Fam& f = fams[k];
            GridPad(c, Form("pR%zu", k), (int)k, 4, 3, yt);
            Ser s; s.line = false; s.band = false; s.col = 1; s.mk = 20; s.ms = 1.3;
            for (size_t i = 0; i < f.v.size(); ++i)
                if (f.r[i].Rok) { s.x.push_back(i + 0.5); s.y.push_back(f.r[i].R); s.e.push_back(f.r[i].Rerr); }
            std::pair<double, double> yr = Range({s}, rPre.R - 0.05, rPre.R + 0.05, 0.1);
            TH1D* fr = Frame(Form("frR%zu", k), f, yr.first, yr.second, "R (0.60 - 0.80)");
            fr->Draw("AXIS");
            TLine* l1 = new TLine(0, rPre.R, f.v.size(), rPre.R); l1->SetLineStyle(2); l1->SetLineColor(kColSim); l1->Draw();
            DrawSer(s);
        }
        GridPad(c, "pRleg", 11, 4, 3, yt);
        Txt(0.05, 0.75, Form("R(pre) = %.3f #pm %.3f", rPre.R, rPre.Rerr), 0.08);
        Txt(0.05, 0.62, Form("R(pre0) = %.3f #pm %.3f", rPre0.R, rPre0.Rerr), 0.08);
        Txt(0.05, 0.49, "dashed: R(pre)", 0.07, kColSim);
        Txt(0.05, 0.36, "error: data (Poisson, CB) and sim (groups)", 0.06);
        Save(c, "cuts_R", false);
    }
    printf("1D figures done (%.0f s)\n", swAll.RealTime()); swAll.Continue();

    // ------------------------------------------------------------------ 2b. maps
    struct Map2 {
        std::string key, name, xt, yt;
        int kx, ky;
        std::vector<double> xv, yv;
        std::vector<std::vector<Res>> r, rD;   // [ix][iy]; rD: Dalitz-plot weight on the bump classes
    };
    std::vector<Map2> maps(2);
    maps[0].key = "lab_mpplo"; maps[0].name = "(#alpha_{lab} <, M_{#pi#pi} >)"; maps[0].xt = "#alpha_{lab} < [deg]"; maps[0].yt = "M_{#pi#pi} > [GeV/c^{2}]";
    maps[0].kx = kLab; maps[0].ky = kMppLo;
    maps[0].xv = fams[3].v; maps[0].yv = fams[4].v;
    maps[1].key = "oalo_oahi"; maps[1].name = "(OA >, OA <)"; maps[1].xt = "OA > [deg]"; maps[1].yt = "OA < [deg]";
    maps[1].kx = kOaLo; maps[1].ky = kOaHi;
    maps[1].xv = fams[1].v; maps[1].yv = fams[2].v;
    for (Map2& m : maps) {
        m.r.assign(m.xv.size(), std::vector<Res>(m.yv.size()));
        m.rD.assign(m.xv.size(), std::vector<Res>(m.yv.size()));
        for (size_t ix = 0; ix < m.xv.size(); ++ix)
            for (size_t iy = 0; iy < m.yv.size(); ++iy) {
                Cuts s = Pre(); SetCut(s, m.kx, m.xv[ix]); SetCut(s, m.ky, m.yv[iy]);
                m.r[ix][iy] = Run(Pre(), s);
                EOpt od; od.dal = true;
                m.rD[ix][iy] = Run(Pre(), s, od, false);
            }
    }
    printf("maps done (%.0f s)\n", swAll.RealTime()); swAll.Continue();
    // data plateaus on the maps: peak F_stat >= best - 1 (scan_peak maps), bump gain_T within 5 % (scan_bump OA map)
    const Blk dPk0 = ParseBlock("peak_maps.txt", "", "F_stat: columns oa_lab <, rows M(pipi) >");
    const Blk dPk1 = ParseBlock("peak_maps.txt", "", "F_stat: columns OA >, rows OA <");
    const Blk dBu1 = ParseBlock("bump_scan2d.txt", "2D: OA lower (rows) against OA upper (columns)", "gain_T = eff_T / sqrt(eff_N)");
    auto dataPlateau = [&](int w, int im, size_t ix, size_t iy, const Map2& m) -> int {
        // 1 inside, 0 outside, -1 no data map
        const Blk* b = nullptr;
        bool rowsAreY = true;
        double thr = 0;
        if (w == 0 && im == 0) b = &dPk0;
        if (w == 0 && im == 1) b = &dPk1;
        if (w == 1 && im == 1) { b = &dBu1; rowsAreY = false; }
        if (!b || !b->ok) return -1;
        double best = -1e30;
        for (auto& r : b->v) for (double x : r) best = std::max(best, x);
        thr = (w == 0) ? best - 1 : 0.95 * best;
        const double xv = m.xv[ix], yv = m.yv[iy];
        auto eq = [](double a, double bb, bool zeroIsNone) {
            if (zeroIsNone) { if (a == kNoCut) a = 0; if (bb == kNoCut) bb = 0; }
            return (a == kNoCut && bb == kNoCut) || (a != kNoCut && bb != kNoCut && std::fabs(a - bb) < 1e-6);
        };
        const bool xz = (m.kx == kOaLo);
        for (size_t i = 0; i < b->rows.size(); ++i)
            for (size_t j = 0; j < b->cols.size(); ++j) {
                const double rv = b->rows[i], cv = b->cols[j];
                const bool match = rowsAreY ? (eq(cv, xv, xz) && eq(rv, yv, false)) : (eq(rv, xv, xz) && eq(cv, yv, false));
                if (match) return b->v[i][j] >= thr ? 1 : 0;
            }
        return -1;
    };
    {
        gOut = fopen(Form("%s/cuts_maps.txt", kTabDir), "w");
        P("# cuts_maps.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# Two-dimensional scans on \"pre\" (other cuts open). Rows = y cut, columns = x cut. F/F0 relative to pre (F = S/sqrt(S+B), B = phys. bkg);\n");
        P("# F'/F'0 with B -> R * B (R = data (N-CB)/sim (T1+T2) in %.2f-%.2f at the same selection, L_data %.1f UNCONFIRMED); S/B in %%.\n", kRLo, kRHi, kLdata);
        P("# Marks: B = best cell, 5 = within 5%% of the best, 1 = F >= F_best - 1 (counts at L_data), D = inside the data plateau\n");
        P("# (peak: F_stat >= best - 1 in peak_maps.txt [the OA map of the data has oa_lab < 50 applied]; bump OA map: gain_T >= 0.95 best).\n");
        for (int w = 0; w < 2; ++w)
            for (int im = 0; im < 2; ++im) {
                const Map2& m = maps[im];
                const char* wn = w == 0 ? "PEAK" : "BUMP";
                double best = -1, bestR = -1; size_t bx = 0, by = 0, bxR = 0, byR = 0;
                for (size_t ix = 0; ix < m.xv.size(); ++ix)
                    for (size_t iy = 0; iy < m.yv.size(); ++iy) {
                        const FR f = Fom(m.r[ix][iy].a, w), g = Fom(m.r[ix][iy].a, w, m.r[ix][iy].R, m.r[ix][iy].R0);
                        if (f.F > best) { best = f.F; bx = ix; by = iy; }
                        if (g.F > bestR) { bestR = g.F; bxR = ix; byR = iy; }
                    }
                for (int q = 0; q < 4; ++q) {
                    const char* qn[4] = {"F/F0", "F'/F'0 (B*R)", "S/B [%]", "R"};
                    P("\n%s window, map %s, %s\n%10s", wn, m.key.c_str(), qn[q], "y \\ x");
                    for (double x : m.xv) P(" %9s", ValTxt(m.kx, x).c_str());
                    P("\n");
                    for (size_t iy = 0; iy < m.yv.size(); ++iy) {
                        P("%10s", ValTxt(m.ky, m.yv[iy]).c_str());
                        for (size_t ix = 0; ix < m.xv.size(); ++ix) {
                            const Res& r = m.r[ix][iy];
                            const FR f = Fom(r.a, w), g = Fom(r.a, w, r.R, r.R0);
                            double val = q == 0 ? f.rel : q == 1 ? g.rel : q == 2 ? 100 * f.SB : r.R;
                            std::string mk;
                            if (q == 0) {
                                if (ix == bx && iy == by) mk += "B";
                                if (f.F >= 0.95 * best) mk += "5";
                                if (f.F >= best - 1) mk += "1";
                                if (dataPlateau(w, im, ix, iy, m) == 1) mk += "D";
                            }
                            if (q == 1) {
                                if (ix == bxR && iy == byR) mk += "B";
                                if (g.F >= 0.95 * bestR) mk += "5";
                            }
                            if (q == 3 && !r.Rok) mk = "*";
                            P(" %5.3f%-4s", val, mk.c_str());
                        }
                        P("\n");
                    }
                }
                const FR fb = Fom(m.r[bx][by].a, w);
                P("best (sim): x %s, y %s, F/F0 = %.4f +- %.4f, S/B = %.4f\n", ValTxt(m.kx, m.xv[bx]).c_str(), ValTxt(m.ky, m.yv[by]).c_str(), fb.rel, fb.relErr, fb.SB);
                const FR gb = Fom(m.r[bxR][byR].a, w, m.r[bxR][byR].R, m.r[bxR][byR].R0);
                P("best (B*R): x %s, y %s, F'/F'0 = %.4f\n", ValTxt(m.kx, m.xv[bxR]).c_str(), ValTxt(m.ky, m.yv[byR]).c_str(), gb.rel);
                if (w == 1) {
                    double bD = -1; size_t dx = 0, dy = 0;
                    for (size_t ix = 0; ix < m.xv.size(); ++ix)
                        for (size_t iy = 0; iy < m.yv.size(); ++iy) { const double F = Fom(m.rD[ix][iy].a, 1).F; if (F > bD) { bD = F; dx = ix; dy = iy; } }
                    P("\nBUMP window, map %s, F/F0 with the eta -> 3 pi Dalitz-plot weight (5 = within 5%% of its best)\n%10s", m.key.c_str(), "y \\ x");
                    for (double x : m.xv) P(" %9s", ValTxt(m.kx, x).c_str());
                    P("\n");
                    for (size_t iy = 0; iy < m.yv.size(); ++iy) {
                        P("%10s", ValTxt(m.ky, m.yv[iy]).c_str());
                        for (size_t ix = 0; ix < m.xv.size(); ++ix) {
                            const FR f = Fom(m.rD[ix][iy].a, 1);
                            P(" %5.3f%-4s", f.rel, (std::string(ix == dx && iy == dy ? "B" : "") + (f.F >= 0.95 * bD ? "5" : "")).c_str());
                        }
                        P("\n");
                    }
                    P("best (Dalitz weight): x %s, y %s, F/F0 = %.4f\n", ValTxt(m.kx, m.xv[dx]).c_str(), ValTxt(m.ky, m.yv[dy]).c_str(), Fom(m.rD[dx][dy].a, 1).rel);
                }
            }
        fclose(gOut); gOut = nullptr;
    }
    for (int w = 0; w < 2; ++w) {
        TCanvas* c = new TCanvas(Form("cM%d", w), "", 2400, 1500);
        const double yt = 0.89;
        Header(c, yt, {w == 0 ? "Simulation truth: two-dimensional cut maps, PEAK window (S = peak conversion; direct decay missing in the simulation)"
                              : "Simulation truth: two-dimensional cut maps, BUMP window (S = bump #pi^{0} Dalitz + conversion)",
                       "Base pre, other cuts open. Left: F/F_{pre}, F = S/#sqrt{S+B}, B = physical background; middle: F'/F'_{pre} with B #times R "
                       "(R = data (N-CB)/sim (T1+T2) in 0.60 - 0.80 at each cell); right: S/B [%].",
                       "Boxes: black solid = sim within 5 % of the best F; black dashed = F #geq F_{best} - 1 (counts at L_{data}); red star = best cell; "
                       "brown (bump, left) = within 5 % with the #eta#rightarrow3#pi Dalitz-plot weight.",
                       "Blue boxes: data plateau (peak: F_{stat} #geq best - 1 in scan_peak; its OA map has #alpha_{lab} < 50 applied; bump OA map: tagged-#eta gain within 5 %).",
                       noteSim[1]});
        for (int im = 0; im < 2; ++im) {
            const Map2& m = maps[im];
            double best = -1, bestR = -1;
            for (size_t ix = 0; ix < m.xv.size(); ++ix)
                for (size_t iy = 0; iy < m.yv.size(); ++iy) {
                    best = std::max(best, Fom(m.r[ix][iy].a, w).F);
                    bestR = std::max(bestR, Fom(m.r[ix][iy].a, w, m.r[ix][iy].R, m.r[ix][iy].R0).F);
                }
            for (int q = 0; q < 3; ++q) {
                c->cd();
                const double h = yt / 2;
                TPad* p = new TPad(Form("pM%d%d%d", w, im, q), "", q / 3.0, yt - (im + 1) * h, (q + 1) / 3.0, yt - im * h);
                p->SetLeftMargin(0.14); p->SetRightMargin(0.13); p->SetTopMargin(0.08); p->SetBottomMargin(0.14);
                p->Draw(); p->cd();
                const int nx = (int)m.xv.size(), ny = (int)m.yv.size();
                TH2D* h2 = new TH2D(Form("hM%d%d%d", w, im, q), "", nx, 0, nx, ny, 0, ny);
                h2->SetDirectory(nullptr);
                for (int ix = 0; ix < nx; ++ix) h2->GetXaxis()->SetBinLabel(ix + 1, ValLab(m.kx, m.xv[ix]).c_str());
                for (int iy = 0; iy < ny; ++iy) h2->GetYaxis()->SetBinLabel(iy + 1, ValLab(m.ky, m.yv[iy]).c_str());
                double zlo = 1e30, zhi = -1e30;
                for (int ix = 0; ix < nx; ++ix)
                    for (int iy = 0; iy < ny; ++iy) {
                        const Res& r = m.r[ix][iy];
                        const FR f = Fom(r.a, w), g = Fom(r.a, w, r.R, r.R0);
                        const double val = q == 0 ? f.rel : q == 1 ? g.rel : 100 * f.SB;
                        h2->SetBinContent(ix + 1, iy + 1, val);
                        zlo = std::min(zlo, val); zhi = std::max(zhi, val);
                    }
                h2->SetMinimum(zlo); h2->SetMaximum(zhi);
                h2->GetXaxis()->SetTitle(m.xt.c_str()); h2->GetYaxis()->SetTitle(m.yt.c_str());
                h2->GetXaxis()->SetTitleSize(0.055); h2->GetYaxis()->SetTitleSize(0.055);
                h2->GetXaxis()->SetLabelSize(0.055); h2->GetYaxis()->SetLabelSize(0.055); h2->GetZaxis()->SetLabelSize(0.045);
                h2->GetYaxis()->SetTitleOffset(1.15);
                h2->SetMarkerSize(nx > 9 ? 1.05 : 1.3);
                gStyle->SetPaintTextFormat(q == 2 ? "4.1f" : "4.2f");
                h2->Draw("COLZ");
                for (int ix = 0; ix < nx; ++ix)
                    for (int iy = 0; iy < ny; ++iy) {
                        TLatex tl; tl.SetTextAlign(22); tl.SetTextFont(42); tl.SetTextSize(nx > 9 ? 0.029 : 0.036);
                        tl.DrawLatex(ix + 0.5, iy + 0.5, Form(q == 2 ? "%.1f" : "%.2f", h2->GetBinContent(ix + 1, iy + 1)));
                    }
                double bestD = -1;
                for (int ix = 0; ix < nx; ++ix) for (int iy = 0; iy < ny; ++iy) bestD = std::max(bestD, Fom(m.rD[ix][iy].a, 1).F);
                for (int ix = 0; ix < nx; ++ix)
                    for (int iy = 0; iy < ny; ++iy) {
                        const Res& r = m.r[ix][iy];
                        const double F = q == 1 ? Fom(r.a, w, r.R, r.R0).F : Fom(r.a, w).F;
                        const double B = q == 1 ? bestR : best;
                        if (F >= 0.95 * B) { TBox* b = new TBox(ix + 0.04, iy + 0.04, ix + 0.96, iy + 0.96); b->SetFillStyle(0); b->SetLineColor(1); b->SetLineWidth(3); b->Draw(); }
                        if (q != 1 && F >= B - 1) { TBox* b = new TBox(ix + 0.08, iy + 0.08, ix + 0.92, iy + 0.92); b->SetFillStyle(0); b->SetLineColor(1); b->SetLineStyle(2); b->SetLineWidth(2); b->Draw(); }
                        if (F >= B - 1e-12) { TMarker* mk = new TMarker(ix + 0.5, iy + 0.80, 29); mk->SetMarkerSize(2.0); mk->SetMarkerColor(kRed + 1); mk->Draw(); }
                        if (dataPlateau(w, im, ix, iy, m) == 1) { TBox* b = new TBox(ix + 0.12, iy + 0.12, ix + 0.88, iy + 0.88); b->SetFillStyle(0); b->SetLineColor(kColData); b->SetLineWidth(3); b->Draw(); }
                        if (w == 1 && q == 0 && Fom(m.rD[ix][iy].a, 1).F >= 0.95 * bestD) { TBox* b = new TBox(ix + 0.16, iy + 0.16, ix + 0.84, iy + 0.84); b->SetFillStyle(0); b->SetLineColor(kColTag); b->SetLineStyle(1); b->SetLineWidth(3); b->Draw(); }
                    }
                const char* qn[3] = {"F/F_{pre} (sim)", "F'/F'_{pre} (B #times R)", "S/B [%] (sim)"};
                Txt(0.14, 0.935, Form("%s  %s", m.name.c_str(), qn[q]), 0.055);
            }
        }
        Save(c, w == 0 ? "cuts_map_peak" : "cuts_map_bump", true);
    }

    // ------------------------------------------------------------------ 3. ribbons
    struct Step { std::string lab; bool applied; int kind; double nom; std::vector<double> alt; };
    struct Flow { std::string key; int win; std::vector<Step> st; };
    std::vector<Flow> flows(2);
    flows[0].key = "peak"; flows[0].win = 0;
    flows[0].st = {{"preselection", true, -1, 0, {}},
                   {"#alpha_{lab} < 55#circ", true, kLab, 55, {50, 60}},
                   {"M_{#pi#pi} < 0.46 (n.a.)", false, kMppHi, 0.46, {0.45, 0.46, 0.47}},
                   {"M_{#pi#pi} > 0.30", true, kMppLo, 0.30, {0.28, 0.32}},
                   {"OA > 2#circ (n.a.)", false, kOaLo, 2, {1.5, 2, 2.5}},
                   {"OA < 20#circ (n.a.)", false, kOaHi, 20, {18, 20, 22}},
                   {"v_{z} > -160", true, kVz, -160, {-170, -150}},
                   {"both #pi in vertex fit", true, kVtx, 2, {1, 3}}};
    flows[1].key = "bump"; flows[1].win = 1;
    flows[1].st = {{"preselection", true, -1, 0, {}},
                   {"#alpha_{lab} < 40#circ", true, kLab, 40, {35, 45}},
                   {"M_{#pi#pi} < 0.41", true, kMppHi, 0.41, {0.40, 0.42}},
                   {"M_{#pi#pi} > 0.30", true, kMppLo, 0.30, {0.28, 0.32}},
                   {"OA > 2#circ (n.a.)", false, kOaLo, 2, {1.5, 2, 2.5}},
                   {"OA < 20#circ (n.a.)", false, kOaHi, 20, {18, 20, 22}},
                   {"v_{z} > -160", true, kVz, -160, {-170, -150}},
                   {"both #pi in vertex fit", true, kVtx, 2, {1, 3}}};
    const auto dRibPk = ParseRibbon("ribbons_peak.txt");
    const auto dRibBu = ParseRibbon("ribbons_bump.txt");
    std::vector<Cuts> finalSel(2);
    {
        gOut = fopen(Form("%s/cuts_ribbons.txt", kTabDir), "w");
        P("# cuts_ribbons.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# Ribbon cut flows of the data study (scan_ribbons.C): same order, nominal values, bands; base pre0 = isBest >= 0, PT3, z > -500, cor_ok.\n");
        P("# central = cuts 1..j at nominal; band = cut j alone moved (a priori: oa_lab +-5 deg, M(pipi) upper +-10 MeV, lower +-20 MeV, OA lower +-0.5,\n");
        P("# OA upper +-2 deg, vz +-10 mm, vertex fit >=1 pion / all four); 'if-appl.' = a cut that is NOT applied, put on top of the cuts so far.\n");
        P("# sigma [nb] = sum w / L_sim (mode A, w); S_pk = peak conversion in the sim peak window %.4f-%.4f, B_pk = phys. bkg there;\n", gPkLo, gPkHi);
        P("# S_bu = bump D+C in %.4f-%.4f, B_bu = phys. bkg there. F/F0 relative to pre0, F = S/sqrt(S+B); F' with B*R. N_eff of S and B in brackets.\n", kBuLo, kBuHi);
        P("# Data columns (from ribbons_*.txt, central lines): peak S/S0 and (N-CB)/(N-CB)0 in the peak window, F_stat/F0; bump eff_T, eff_E, W/W0, gain_T.\n");
    }
    for (int fl = 0; fl < 2; ++fl) {
        Flow& F = flows[fl];
        const int w = F.win, ow = 1 - w;
        const int n = (int)F.st.size();
        std::vector<Cuts> cen(n);
        std::vector<Res> rc(n), rcD(n);
        std::vector<std::vector<Res>> ra(n);
        Cuts cur = Pre0();
        for (int j = 0; j < n; ++j) {
            const Step& s = F.st[j];
            Cuts prev = cur;
            if (j > 0 && s.applied) SetCut(cur, s.kind, s.nom);
            cen[j] = cur;
            rc[j] = Run(Pre0(), cur);
            { EOpt od; od.dal = true; rcD[j] = Run(Pre0(), cur, od, false); }
            for (double v : s.alt) { Cuts a = prev; SetCut(a, s.kind, v); ra[j].push_back(Run(Pre0(), a)); }
        }
        finalSel[fl] = cur;
        const auto& dr = fl == 0 ? dRibPk : dRibBu;
        // table
        P("\n==================== %s flow (window: %s) ====================\n", F.key.c_str(), w == 0 ? "peak" : "bump");
        P("%-5s %-28s | %8s %7s %9s %7s | %8s %7s %9s %7s | %6s %6s | %6s %6s %6s | %6s | %s\n", "role", "cut", "S_pk", "err", "B_pk", "err", "S_bu", "err", "B_bu",
          "err", "S/S0", "B/B0", "F/F0", "err", "F'/F'0", "R", "data");
        auto line = [&](const char* role, const std::string& cut, const Res& r, int j, bool central) {
            const FR a = Fom(r.a, w), b = Fom(r.a, w, r.R, r.R0);
            const Acc& s0 = r.a[TS(w)];
            const Acc& b0 = r.a[TB(w)];
            P("%-5s %-28s | %8.4f %7.4f %9.4f %7.4f | %8.4f %7.4f %9.4f %7.4f | %6.4f %6.4f | %6.4f %6.4f %6.4f | %6.3f |", role, cut.c_str(),
              r.a[tSpk].xs(), r.a[tSpk].xsErr(), r.a[tB2pk].xs(), r.a[tB2pk].xsErr(), r.a[tSbu].xs(), r.a[tSbu].xsErr(), r.a[tB2bu].xs(), r.a[tB2bu].xsErr(),
              s0.eps(), b0.eps(), a.rel, a.relErr, b.rel, r.R);
            if (central && j < (int)dr.size()) {
                const auto& x = dr[j];
                if (fl == 0 && x.size() > 2 && dr[0].size() > 2 && x[0].size() >= 4 && x[1].size() >= 1 && x[2].size() >= 2)
                    P(" S/S0 %.3f, (N-CB)/0 %.3f, F/F0 %.3f", x[1][0] / dr[0][1][0], x[0][3] / dr[0][0][3], x[2][1] / dr[0][2][1]);
                if (fl == 1 && x.size() > 6 && x[0].size() >= 5 && x[2].size() >= 2 && x[4].size() >= 2 && x[6].size() >= 2)
                    P(" eff_T %.3f +- %.3f, eff_E %.3f +- %.3f, W/W0 %.3f, gain_T %.3f", x[4][0], x[4][1], x[2][0], x[2][1], x[0][4], x[6][1]);
            }
            P("  [N_eff S %.0f, B %.0f]\n", s0.NeffP(), b0.NeffP());
        };
        for (int j = 0; j < n; ++j) {
            line("cen", F.st[j].lab, rc[j], j, true);
            if (w == 1) P("      Dalitz-plot weight: S_bu %.4f nb, S/S0 %.4f +- %.4f, tagged kept %.4f +- %.4f, F/F0 %.4f\n", rcD[j].a[tSbu].xs(), rcD[j].a[tSbu].eps(),
                          rcD[j].a[tSbu].epsErr(), rcD[j].a[tTag].eps(), rcD[j].a[tTag].epsErr(), Fom(rcD[j].a, 1).rel);
            for (size_t a = 0; a < ra[j].size(); ++a)
                line(F.st[j].applied ? " band" : " ifap", Form("%s = %s", F.st[j].lab.c_str(), ValTxt(F.st[j].kind, F.st[j].alt[a]).c_str()), ra[j][a], j, false);
        }
        // figure
        TCanvas* c = new TCanvas(Form("cRib%d", fl), "", 2000, 1900);
        const double yt = 0.90;
        const int col = w == 0 ? kColPeakConv : kColBumpDalitz;
        Header(c, yt, {Form("Simulation truth: ribbon cut flow of the data study, %s selection (window %s)", w == 0 ? "PEAK" : "BUMP",
                            (w == 0 ? TString::Format("%.4f - %.4f, sim peak", gPkLo, gPkHi) : TString::Format("%.4f - %.4f", kBuLo, kBuHi)).Data()),
                       "Line: cuts 1..j at nominal values; box at step j: cut j alone moved by its a-priori step (filled: applied; dashed: not applied, range if applied). "
                       "Base pre0 = isBest #geq 0, PT3, z > -500, cor_ok (v_{z} > -160 is a step).",
                       Form("(a) #sigma in the windows: S %s, B = physical background (filled: this window, open: the other window). (b) kept fraction; blue/black: data "
                            "(scan_ribbons). (c) F/F_{0}, F = S/#sqrt{S+B}; dashed: B #times R.", w == 0 ? "= peak conversion" : "= bump D+C"),
                       noteSim[1]});
        const double padH[3] = {0.40, 0.25, 0.25};
        double ytop = yt;
        std::vector<TPad*> pp(3);
        for (int q = 0; q < 3; ++q) {
            c->cd();
            const double hh = padH[q] / 0.90 * yt;
            const double ybot = (q == 2) ? 0.0 : ytop - hh;
            pp[q] = new TPad(Form("pRib%d%d", fl, q), "", 0, ybot, 1, ytop);
            pp[q]->SetLeftMargin(0.10); pp[q]->SetRightMargin(0.02);
            pp[q]->SetTopMargin(0.04); pp[q]->SetBottomMargin(q == 2 ? 0.30 : 0.03);
            pp[q]->SetGridy(1);
            pp[q]->Draw();
            ytop = ybot;
        }
        auto frameR = [&](const char* nm, double lo, double hi, const char* yt2, bool labels, double ls) {
            TH1D* h = new TH1D(nm, "", n, -0.5, n - 0.5);
            h->SetDirectory(nullptr);
            for (int j = 0; j < n; ++j) h->GetXaxis()->SetBinLabel(j + 1, labels ? F.st[j].lab.c_str() : "");
            h->GetXaxis()->SetLabelSize(labels ? 0.085 : 0);
            h->GetYaxis()->SetTitle(yt2); h->GetYaxis()->SetTitleSize(0.075 * ls); h->GetYaxis()->SetLabelSize(0.065 * ls);
            h->GetYaxis()->SetTitleOffset(0.6 / ls);
            h->GetYaxis()->SetNdivisions(ls > 1.2 ? 505 : 510);
            h->SetMinimum(lo); h->SetMaximum(hi);
            h->Draw("AXIS");
            return h;
        };
        // quantity along the flow: central, band (min/max over the step's lines), stat error
        auto ser = [&](std::function<double(const Res&)> q, std::function<double(const Res&)> qe, std::vector<double>& c0, std::vector<double>& e0,
                       std::vector<double>& lo, std::vector<double>& hi) {
            c0.assign(n, 0); e0.assign(n, 0); lo.assign(n, 0); hi.assign(n, 0);
            for (int j = 0; j < n; ++j) {
                c0[j] = q(rc[j]); e0[j] = qe(rc[j]);
                lo[j] = hi[j] = c0[j];
                for (const Res& r : ra[j]) { lo[j] = std::min(lo[j], q(r)); hi[j] = std::max(hi[j], q(r)); }
            }
        };
        auto drawRib = [&](const std::vector<double>& c0, const std::vector<double>& e0, const std::vector<double>& lo, const std::vector<double>& hi,
                           int cl, int mk, int ls, double dx, double ylo, double yhi, bool logy) {
            for (int j = 1; j < n; ++j) {
                if (hi[j] - lo[j] <= 0) continue;
                const double a = std::max(lo[j], logy ? ylo * 1.0001 : ylo), b = std::min(hi[j], yhi);
                TBox* bx = new TBox(j + dx - 0.13, a, j + dx + 0.13, b);
                if (F.st[j].applied) { bx->SetFillColorAlpha(cl, 0.30); bx->SetLineColor(cl); }
                else { bx->SetFillStyle(0); bx->SetLineColor(cl); bx->SetLineStyle(2); bx->SetLineWidth(2); }
                bx->Draw(F.st[j].applied ? "" : "l");
            }
            std::vector<double> x(n), ex(n, 0.);
            for (int j = 0; j < n; ++j) x[j] = j + dx;
            TGraphErrors* g = new TGraphErrors(n, &x[0], &c0[0], &ex[0], &e0[0]);
            g->SetLineColor(cl); g->SetMarkerColor(cl); g->SetMarkerStyle(mk); g->SetMarkerSize(1.6); g->SetLineWidth(3); g->SetLineStyle(ls);
            g->Draw("LP");
        };
        std::vector<double> c0, e0, lo, hi;
        // (a) absolute
        pp[0]->cd();
        pp[0]->SetLogy(1);
        {
            double ymin = 1e30, ymax = -1e30;
            for (int j = 0; j < n; ++j)
                for (int t : {tSpk, tB2pk, tSbu, tB2bu}) {
                    ymin = std::min(ymin, rc[j].a[t].xs()); ymax = std::max(ymax, rc[j].a[t].xs());
                    for (const Res& r : ra[j]) { ymin = std::min(ymin, r.a[t].xs()); ymax = std::max(ymax, r.a[t].xs()); }
                }
            const double ylo = ymin / 2.5, yhi = ymax * 15;
            frameR(Form("frRa%d", fl), ylo, yhi, "#sigma in window [nb]", false, 1.0);
            const int tS = TS(w), tB = TB(w), oS = TS(ow), oB = TB(ow);
            ser([&](const Res& r) { return r.a[tS].xs(); }, [&](const Res& r) { return r.a[tS].xsErr(); }, c0, e0, lo, hi);
            drawRib(c0, e0, lo, hi, col, 20, 1, -0.06, ylo, yhi, true);
            ser([&](const Res& r) { return r.a[tB].xs(); }, [&](const Res& r) { return r.a[tB].xsErr(); }, c0, e0, lo, hi);
            drawRib(c0, e0, lo, hi, kColB2pk, 21, 1, 0.06, ylo, yhi, true);
            const int ocol = w == 0 ? kColBumpDalitz : kColPeakConv;
            ser([&](const Res& r) { return r.a[oS].xs(); }, [&](const Res& r) { return r.a[oS].xsErr(); }, c0, e0, lo, hi);
            drawRib(c0, e0, lo, hi, ocol, 24, 3, -0.18, ylo, yhi, true);
            ser([&](const Res& r) { return r.a[oB].xs(); }, [&](const Res& r) { return r.a[oB].xsErr(); }, c0, e0, lo, hi);
            drawRib(c0, e0, lo, hi, kColB2bu, 25, 3, 0.18, ylo, yhi, true);
            TLegend* lg = new TLegend(0.12, 0.76, 0.98, 0.96);
            lg->SetNColumns(2); lg->SetTextSize(0.045); lg->SetBorderSize(0); lg->SetFillStyle(0);
            auto addP = [&](int cl, int mk, int ls, const char* t) { TGraph* g = new TGraph(); g->SetMarkerColor(cl); g->SetMarkerStyle(mk); g->SetLineColor(cl); g->SetLineStyle(ls); g->SetLineWidth(3); g->SetMarkerSize(1.5); lg->AddEntry(g, t, "lp"); };
            addP(col, 20, 1, w == 0 ? "S: peak conv., peak window" : "S: bump D+C, bump window");
            addP(kColB2pk, 21, 1, w == 0 ? "B: phys. bkg, peak window" : "B: phys. bkg, bump window");
            addP(ocol, 24, 3, w == 0 ? "S: bump D+C, bump window" : "S: peak conv., peak window");
            addP(kColB2bu, 25, 3, w == 0 ? "B: phys. bkg, bump window" : "B: phys. bkg, peak window");
            lg->Draw();
        }
        // (b) kept fraction relative to pre0
        pp[1]->cd();
        {
            frameR(Form("frRb%d", fl), 0.0, 1.15, "kept fraction", false, 1.6);
            const int tS = TS(w), tB = TB(w);
            ser([&](const Res& r) { return r.a[tS].eps(); }, [&](const Res& r) { return r.a[tS].epsErr(); }, c0, e0, lo, hi);
            drawRib(c0, e0, lo, hi, col, 20, 1, -0.06, 0.0, 1.15, false);
            ser([&](const Res& r) { return r.a[tB].eps(); }, [&](const Res& r) { return r.a[tB].epsErr(); }, c0, e0, lo, hi);
            drawRib(c0, e0, lo, hi, kColB2pk, 21, 1, 0.06, 0.0, 1.15, false);
            if (w == 1) {
                ser([&](const Res& r) { return r.a[tTag].eps(); }, [&](const Res& r) { return r.a[tTag].epsErr(); }, c0, e0, lo, hi);
                drawRib(c0, e0, lo, hi, kColTag, 22, 3, -0.18, 0.0, 1.15, false);
                std::vector<double> cd(n), ed(n, 0.);
                for (int j = 0; j < n; ++j) cd[j] = rcD[j].a[tTag].eps();
                drawRib(cd, ed, cd, cd, kColTag, 23, 9, -0.24, 0.0, 1.15, false);
            }
            // data
            TGraphErrors* gS = new TGraphErrors();
            TGraphErrors* gN = new TGraphErrors();
            TGraphErrors* gE = new TGraphErrors();
            for (int j = 0; j < n && j < (int)dr.size(); ++j) {
                const auto& x = dr[j];
                if (fl == 0 && x.size() > 1 && x[0].size() >= 4 && x[1].size() >= 1) {
                    gS->SetPoint(gS->GetN(), j + 0.25, x[1][0] / dr[0][1][0]);
                    gN->SetPoint(gN->GetN(), j + 0.32, x[0][3] / dr[0][0][3]);
                }
                if (fl == 1 && x.size() > 4 && x[4].size() >= 2 && x[0].size() >= 5 && x[2].size() >= 2) {
                    const int k = gS->GetN(); gS->SetPoint(k, j + 0.25, x[4][0]); gS->SetPointError(k, 0, x[4][1]);
                    gN->SetPoint(gN->GetN(), j + 0.32, x[0][4]);
                    const int m = gE->GetN(); gE->SetPoint(m, j + 0.18, x[2][0]); gE->SetPointError(m, 0, x[2][1]);
                }
            }
            gS->SetMarkerStyle(fl == 0 ? 20 : 21); gS->SetMarkerColor(kColData); gS->SetLineColor(kColData); gS->SetMarkerSize(1.6); gS->Draw("PZ");
            gN->SetMarkerStyle(24); gN->SetMarkerColor(1); gN->SetMarkerSize(1.6); gN->Draw("P");
            if (fl == 1) { gE->SetMarkerStyle(26); gE->SetMarkerColor(kColData); gE->SetLineColor(kColData); gE->SetMarkerSize(1.6); gE->Draw("PZ"); }
            TLegend* lg = new TLegend(0.12, 0.05, 0.98, 0.30);
            lg->SetNColumns(3); lg->SetTextSize(0.060); lg->SetBorderSize(0); lg->SetFillStyle(0);
            auto addP = [&](int cl, int mk, int ls, const char* t, const char* o) { TGraph* g = new TGraph(); g->SetMarkerColor(cl); g->SetMarkerStyle(mk); g->SetLineColor(cl); g->SetLineStyle(ls); g->SetLineWidth(3); g->SetMarkerSize(1.5); lg->AddEntry(g, t, o); };
            addP(col, 20, 1, "sim S", "lp");
            addP(kColB2pk, 21, 1, "sim B (phys. bkg)", "lp");
            if (w == 1) addP(kColTag, 22, 3, "sim bump D+C, tagged", "lp");
            if (w == 1) addP(kColTag, 23, 9, "sim tagged, #eta#rightarrow3#pi Dalitz-plot weight", "lp");
            addP(kColData, fl == 0 ? 20 : 21, 1, fl == 0 ? "data S/S_{0} (fit)" : "data tagged #eta #varepsilon_{T}", "p");
            if (fl == 1) addP(kColData, 26, 1, "data edge step #varepsilon_{E}", "p");
            addP(1, 24, 1, "data (N-CB)/(N-CB)_{0} in window", "p");
            lg->Draw();
        }
        // (c) FOM
        pp[2]->cd();
        {
            std::vector<double> fc(n), fe(n), flo(n), fhi(n), gc(n), glo(n), ghi(n);
            for (int j = 0; j < n; ++j) {
                const FR a = Fom(rc[j].a, w), b = Fom(rc[j].a, w, rc[j].R, rc[j].R0);
                fc[j] = a.rel; fe[j] = a.relErr; gc[j] = b.rel;
                flo[j] = fhi[j] = a.rel; glo[j] = ghi[j] = b.rel;
                for (const Res& r : ra[j]) {
                    const double x = Fom(r.a, w).rel, y = Fom(r.a, w, r.R, r.R0).rel;
                    flo[j] = std::min(flo[j], x); fhi[j] = std::max(fhi[j], x); glo[j] = std::min(glo[j], y); ghi[j] = std::max(ghi[j], y);
                }
            }
            std::vector<double> dF;
            for (int j = 0; j < n && j < (int)dr.size(); ++j) {
                const auto& x = dr[j];
                if (fl == 0 && x.size() > 2 && x[2].size() >= 2) dF.push_back(x[2][1] / dr[0][2][1]);
                if (fl == 1 && x.size() > 6 && x[6].size() >= 2) dF.push_back(x[6][1]);
            }
            double lo2 = 1e30, hi2 = -1e30;
            for (int j = 0; j < n; ++j) { lo2 = std::min({lo2, flo[j] - fe[j], glo[j]}); hi2 = std::max({hi2, fhi[j] + fe[j], ghi[j]}); }
            if (w == 1) for (int j = 0; j < n; ++j) { lo2 = std::min(lo2, Fom(rcD[j].a, 1).rel); hi2 = std::max(hi2, Fom(rcD[j].a, 1).rel); }
            for (double x : dF) { lo2 = std::min(lo2, x); hi2 = std::max(hi2, x); }
            const double pad = 0.1 * (hi2 - lo2);
            frameR(Form("frRc%d", fl), lo2 - pad, hi2 + pad, "F / F_{0}", true, 1.6);
            drawRib(fc, fe, flo, fhi, col, 20, 1, -0.06, lo2 - pad, hi2 + pad, false);
            std::vector<double> ze(n, 0.);
            drawRib(gc, ze, glo, ghi, col, 24, 2, 0.06, lo2 - pad, hi2 + pad, false);
            if (w == 1) {
                std::vector<double> cd(n);
                for (int j = 0; j < n; ++j) cd[j] = Fom(rcD[j].a, 1).rel;
                drawRib(cd, ze, cd, cd, kColTag, 23, 9, -0.18, lo2 - pad, hi2 + pad, false);
            }
            TGraph* gd = new TGraph();
            for (size_t j = 0; j < dF.size(); ++j) gd->SetPoint(gd->GetN(), j + 0.25, dF[j]);
            gd->SetMarkerStyle(fl == 0 ? 20 : 21); gd->SetMarkerColor(kColData); gd->SetMarkerSize(1.6); gd->Draw("P");
            TLegend* lg = new TLegend(0.12, 0.70, 0.98, 0.95);
            lg->SetNColumns(3); lg->SetTextSize(0.060); lg->SetBorderSize(0); lg->SetFillStyle(0);
            auto addP = [&](int cl, int mk, int ls, const char* t, const char* o) { TGraph* g = new TGraph(); g->SetMarkerColor(cl); g->SetMarkerStyle(mk); g->SetLineColor(cl); g->SetLineStyle(ls); g->SetLineWidth(3); g->SetMarkerSize(1.5); lg->AddEntry(g, t, o); };
            addP(col, 20, 1, "sim F/F_{0}", "lp");
            addP(col, 24, 2, "sim, B #times R", "lp");
            if (w == 1) addP(kColTag, 23, 9, "sim, Dalitz-plot weight", "lp");
            addP(kColData, fl == 0 ? 20 : 21, 1, fl == 0 ? "data F_{stat}/F_{stat,0}" : "data #varepsilon_{T}/#sqrt{#varepsilon_{N}}", "p");
            lg->Draw();
        }
        Save(c, Form("cuts_ribbon_%s", F.key.c_str()), true);
    }
    fclose(gOut); gOut = nullptr;
    printf("ribbons done (%.0f s)\n", swAll.RealTime()); swAll.Continue();

    // ------------------------------------------------------------------ variants
    {
        gOut = fopen(Form("%s/cuts_variants.txt", kTabDir), "w");
        P("# cuts_variants.txt -- plots/eta_sim/sim_cuts.C\n");
        P("# sigma [nb] +- stat (N_eff) of the window templates in four variants: partition mode A/B, weight w (PDG) / w_db (converter database).\n");
        P("# peak flow final = pre0 + oa_lab < 55 + M(pipi) > 0.30 + vz > -160 + both pions in vertex fit; bump flow final = pre0 + oa_lab < 40 +\n");
        P("# 0.30 < M(pipi) < 0.41 + vz > -160 + both pions in vertex fit. S_pk window = sim peak window %.4f-%.4f.\n", gPkLo, gPkHi);
        struct SV { const char* name; Cuts c; };
        Cuts oa2 = Pre(); oa2.oaLo = 2;
        const SV sv[] = {{"pre", Pre()}, {"pre + OA > 2", oa2}, {"peak flow final", finalSel[0]}, {"bump flow final", finalSel[1]}};
        const int tt[] = {tSpk, tB2pk, tSbu, tB2bu, tTag, tT3pk, tT3bu, tT4pk, tT4bu, tSimR};
        for (const SV& s : sv) {
            P("\n%s\n%-40s", s.name, "template");
            const char* vn[4] = {"A, w", "A, w_db", "B, w", "B, w_db"};
            for (const char* v : vn) P(" | %-28s", v);
            P("\n");
            std::vector<std::vector<Acc>> va;
            for (int v = 0; v < 4; ++v) { EOpt o; o.mode = v / 2; o.db = v % 2; va.push_back(Eval(Pre0(), s.c, o)); }
            for (int t : tt) {
                P("%-40s", TgtName(t));
                for (int v = 0; v < 4; ++v) P(" | %8.4f +- %6.4f (%8.0f)", va[v][t].xs(), va[v][t].xsErr(), va[v][t].NeffP());
                P("\n");
            }
            P("%-40s", "F = S/sqrt(S+B) peak / bump [counts]");
            for (int v = 0; v < 4; ++v) P(" | %8.2f / %8.2f          ", Fom(va[v], 0).F, Fom(va[v], 1).F);
            P("\n");
        }
        fclose(gOut); gOut = nullptr;
    }

    // ------------------------------------------------------------------ window figure
    {
        TCanvas* c = new TCanvas("cWin", "", 2000, 900);
        const double yt = 0.86;
        Header(c, yt, {"Simulated M(#pi^{+}#pi^{-}e^{+}e^{-}) around the #eta: classes and windows (selection pre, mode A, w)",
                       Form("Sim peak (conversion): fit %.4f #pm %.4f GeV/c^{2}, #sigma = %.1f #pm %.1f MeV (data: 0.5439, 5.4 MeV). "
                            "Windows: sim peak %.4f-%.4f (green), data peak %.4f-%.4f (blue), bump %.4f-%.4f (orange).",
                            gPkMean, gPkMeanE, 1e3 * gPkSig, 1e3 * gPkSigE, gPkLo, gPkHi, kPkLoD, kPkHiD, kBuLo, kBuHi),
                       std::string(NoteSimNorm()) + "; " + NoteNoDirect()});
        for (int q = 0; q < 2; ++q) {
            c->cd();
            TPad* p = new TPad(Form("pW%d", q), "", q * 0.5, 0, (q + 1) * 0.5, yt);
            p->SetLeftMargin(0.14); p->SetRightMargin(0.03); p->SetBottomMargin(0.13); p->SetTopMargin(0.04);
            p->SetLogy(q == 0);
            p->Draw(); p->cd();
            if (q == 0) {
                std::vector<TH1D*> hs;
                const int cols[6] = {kColBumpDalitz, kColBumpConv, kColPeakConv, kColPhysBkg, kColFake, kColDiff};
                const int ls[6] = {1, 1, 1, 1, 2, 2};
                for (int k = 0; k < 6; ++k) {
                    TH1D* h = wsWin[k].Hist(Form("hW%d", k), 2);
                    h->SetLineColor(cols[k]); h->SetLineWidth(3); h->SetLineStyle(ls[k]); h->SetMarkerSize(0);
                    hs.push_back(h);
                }
                hs[3]->SetFillColor(kColPhysFill);
                hs[3]->SetMinimum(1.0); hs[3]->SetMaximum(5e4);
                hs[3]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]");
                hs[3]->GetYaxis()->SetTitle("d#sigma/dM [nb/(GeV/c^{2})]");
                hs[3]->GetXaxis()->SetTitleSize(0.05); hs[3]->GetYaxis()->SetTitleSize(0.05);
                hs[3]->GetXaxis()->SetLabelSize(0.045); hs[3]->GetYaxis()->SetLabelSize(0.045);
                hs[3]->Draw("HIST");
                auto band = [&](double a, double b, int cl, double y1, double y2) { TBox* bx = new TBox(a, y1, b, y2); bx->SetFillColorAlpha(cl, 0.15); bx->Draw(); };
                band(kBuLo, kBuHi, kColBumpDalitz, 1.0, 5e4);
                band(gPkLo, gPkHi, kColPeakConv, 1.0, 5e4);
                band(kPkLoD, kPkHiD, kColData, 1.0, 5e4);
                for (int k = 0; k < 6; ++k) if (k != 3) hs[k]->Draw("HIST SAME");
                hs[3]->Draw("AXIS SAME");
                TLegend* lg = new TLegend(0.16, 0.70, 0.97, 0.96);
                lg->SetNColumns(2); lg->SetTextSize(0.032); lg->SetBorderSize(0); lg->SetFillStyle(0);
                const char* nm[6] = {"bump, #pi^{0} Dalitz", "bump, conversion", "peak, conversion", "physical background (T2)", "fake lepton (info)", "different sources (info)"};
                for (int k = 0; k < 6; ++k) lg->AddEntry(hs[k], nm[k], k == 3 ? "f" : "l");
                lg->Draw();
            } else {
                TH1D* h = (TH1D*)hPk->Clone("hPkDraw");
                h->Scale(1.0 / (kLsim * h->GetBinWidth(1)));
                h->SetLineColor(kColPeakConv); h->SetMarkerColor(kColPeakConv); h->SetMarkerStyle(20); h->SetLineWidth(2);
                h->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]");
                h->GetYaxis()->SetTitle("d#sigma/dM [nb/(GeV/c^{2})]");
                h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
                h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
                h->SetMinimum(0); h->SetMaximum(h->GetMaximum() * 1.35);
                h->Draw("E1");
                const double ymax = h->GetMaximum();
                auto band = [&](double a, double b, int cl) { TBox* bx = new TBox(a, 0, b, ymax); bx->SetFillColorAlpha(cl, 0.12); bx->Draw(); };
                band(gPkLo, gPkHi, kColPeakConv);
                band(kPkLoD, kPkHiD, kColData);
                TF1* g = new TF1("gDraw", "gaus", gPkMean - 2 * gPkSig, gPkMean + 2 * gPkSig);
                const double area = rPre.a[tPC].xs() * 0;   // not used
                (void)area;
                // fitted Gaussian in the drawn units
                TF1 gf("gpk2", "gaus", 0.52, 0.57);
                TH1D* hf = (TH1D*)h->Clone("hPkFit2");
                gf.SetParameters(hf->GetMaximum(), gPkMean, gPkSig);
                gf.SetRange(gPkMean - 2 * gPkSig, gPkMean + 2 * gPkSig);
                hf->Fit(&gf, "QR0");
                g->SetParameters(gf.GetParameter(0), gf.GetParameter(1), gf.GetParameter(2));
                g->SetLineColor(kColPeakConv); g->SetLineWidth(3); g->Draw("SAME");
                // data peak shape with the same height, for the position
                TF1* gd = new TF1("gData", "gaus", 0.52, 0.57);
                gd->SetParameters(gf.GetParameter(0), 0.5439, 0.0054);
                gd->SetLineColor(kColData); gd->SetLineStyle(2); gd->SetLineWidth(3); gd->Draw("SAME");
                h->Draw("E1 SAME");
                TLegend* lg = new TLegend(0.16, 0.72, 0.97, 0.96);
                lg->SetTextSize(0.032); lg->SetBorderSize(0); lg->SetFillStyle(0);
                lg->AddEntry(h, "sim: peak conversion (pre), 2.5 MeV bins", "lep");
                lg->AddEntry(g, Form("Gaussian fit: %.4f, #sigma %.1f MeV", gPkMean, 1e3 * gPkSig), "l");
                lg->AddEntry(gd, "data peak: 0.5439, #sigma 5.4 MeV (same height)", "l");
                lg->Draw();
            }
        }
        Save(c, "cuts_windows", false);
    }
    printf("sim_cuts done in %.0f s\n", swAll.RealTime());
}
