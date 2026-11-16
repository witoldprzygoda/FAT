// opt_common_plain.h -- COPY of plots/exp_report/opt_common_cnt.h made by the short-report figure agent (integ6, 2026-11-16):
// internal catalogue labels in the figure text replaced by plain names (glossary_v2.md); outputs redirected:
// figures of the scope -> docs/latex/<dir>/<base>_plain.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/run/opt
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/mk_plain.py; the physics code is unchanged.
// opt_common_cnt.h -- COPY of plots/exp_report/opt_common.h made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/opt
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// opt_common.h  --  shared code of the chapter "Optimisation of the selection on
//                   data" (prefix opt) of the experimental-data report
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Header-only; included by plots/exp_report/opt_named.C, opt_scans.C,
// opt_sculpt.C and opt_bias.C. Builds on plots/ScanEngine.h (frozen peak
// shape, P_main, model variations, edge step).
//
// Contents
//   A. paths, windows of the simulation study
//   B. Cuts: one selection as a set of cut values, Pass(), Describe(), the named
//      selections (pre, peak nominal, bump nominal, vertex, frozen PEAK, DIRECT,
//      BUMP, X17, slides, old code)
//   C. eta peak: P_main + model spread + CB sensitivity (+-3 %)
//   D. bump: edge-step visibility E (+ model spread) and the ECAL-tagged eta T;
//      copies of the procedures of plots/scan_bump.C (SliceCfg, M5c, Tagged,
//      TagCfg, TagOpt, EdgeSpread) so that the numbers are the same
//   E. drawing helpers (style, vector PDF via EPS + ps2pdf -dEPSCrop, PNG)
//
// Conventions (fixed for the report)
//   rows = stored particle combinations; cor_ok == 1 (no fallback); masses from
//   the energy-loss corrected momenta (COR); beam T = 4530 MeV; FAT production
//   already applies vertex z > -500 mm and PT3 (trigbit == 8192);
//   CB = 2 sqrt(N++ N--) per 10 MeV bin, interpolated to 2.5 MeV, k = 1;
//   P_main = chi2 fit of N - CB, frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2,
//   0.5225-0.6100; S = whole Gaussian; model spread = half range of S over the
//   48 variations with p > 0.05; F = S / stat (compares selections, NOT a
//   significance); E = sum of the edge-step heights in four M(pi+pi-) slices
//   (rows per 2.5 MeV; a visibility measure, NOT a yield).
// ============================================================================
#ifndef OPT_COMMON_PLAIN_H
#define OPT_COMMON_PLAIN_H

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../ScanEngine.h"

namespace opt {

using scan::PeakCfg;
using scan::PeakFit;
using scan::Row;
using scan::Sel;
using scan::Spectra;
using scan::StepCfg;
using scan::StepFit;

// ============================================================================
// A. paths and windows
// ============================================================================
const char* const kFigOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp";
const char* const kTabOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/run/opt/tab";
const char* const kWork   = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/run/opt/work";
const char* const kPlainScratch = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/run/opt/figs";
const char* const kPlainScope = " opt_bias=opt_bias_plain ";
const char* const kTabIn = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp";   // read-only inputs
const char* const kScanV1 = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan";   // first cut-scan study

// windows of the simulation study (simulation-centred; frozen there)
const double kSimPeakLo = 0.5333, kSimPeakHi = 0.5587;
const double kSimBumpLo = 0.4576, kSimBumpHi = 0.5255;
// X17 mass band of the frozen simulation selection (v2), MeV/c2
const double kX17Lo = 15.5, kX17Hi = 24.5;

// luminosity (START-free, elastic scattering) and the PT3 trigger-word factor
const double kLumi     = 6079.2;      // nb^-1, runs present in the data
const double kLumiAll  = 6080.5;      // nb^-1, all 17 615 runs
const double kPT3fac   = 64.0 / 63.0; // trigger word 8192 holds 63/64 of the PT3 events

// ============================================================================
// B. selections
// ============================================================================
const float kNoCut = -1e30f;
struct Cuts {
    std::string name;     // short name used in tables
    float vzMin  = -160.f;   // vz > vzMin [mm]; kNoCut = no cut (FAT: vz > -500 always)
    int   best   = 0;        // 1: isBest == 1 only (comparison variant)
    float oaLo   = kNoCut;    // OA(e+e-) > oaLo [deg]
    float oaHi   = kNoCut;    // OA(e+e-) < oaHi [deg]
    float oaLab  = kNoCut;    // oa_lab((pi+pi-),(e+e-)) < oaLab [deg]
    float mpLo   = kNoCut;    // M(pi+pi-) > mpLo [GeV/c2]
    float mpHi   = kNoCut;    // M(pi+pi-) < mpHi [GeV/c2]
    float mm4Hi  = kNoCut;    // MM(pi+pi-e+e-) < mm4Hi [GeV/c2]
    float etaLo  = kNoCut;    // eta-rest angle > etaLo [deg]
    float dzPi   = kNoCut;    // dz_pi = max(|z_pi+ - vz|, |z_pi- - vz|) < dzPi [mm]
    int   vtx    = 0;        // 0 none, 1 at least one pion, 2 both pions, 4 all four tracks in the vertex fit
    float meeLo  = kNoCut;    // M(e+e-) >= meeLo [MeV/c2]
    float meeHi  = kNoCut;    // M(e+e-) <  meeHi [MeV/c2]
    int   fold   = -1;       // -1 all event groups, 0 / 1 one half (parity of the event-group index)
};
inline bool Pass(const Cuts& c, const Row& r) {
    if (!(r.cor_ok == 1 && r.isBest >= 0)) return false;
    if (c.fold >= 0 && r.fold != c.fold) return false;
    if (c.vzMin > kNoCut && !(r.vz > c.vzMin)) return false;
    if (c.best && r.isBest != 1) return false;
    if (c.oaLo > kNoCut && !(r.oa > c.oaLo)) return false;
    if (c.oaHi > kNoCut && !(r.oa < c.oaHi)) return false;
    if (c.oaLab > kNoCut && !(r.oa_lab < c.oaLab)) return false;
    if (c.mpLo > kNoCut && !(r.mpipi > c.mpLo)) return false;
    if (c.mpHi > kNoCut && !(r.mpipi < c.mpHi)) return false;
    if (c.mm4Hi > kNoCut && !(r.mm4 < c.mm4Hi)) return false;
    if (c.etaLo > kNoCut && !(r.oa_eta > c.etaLo)) return false;
    if (c.dzPi > kNoCut && !(r.dz_pi < c.dzPi)) return false;
    if (c.vtx == 1 && !(r.pip_isUsedVertex == 1 || r.pim_isUsedVertex == 1)) return false;
    if (c.vtx == 2 && !(r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1)) return false;
    if (c.vtx == 4 && !(r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1 &&
                        r.ep_isUsedVertex == 1 && r.em_isUsedVertex == 1)) return false;
    if (c.meeLo > kNoCut && !(r.mee * 1000.f >= c.meeLo)) return false;
    if (c.meeHi > kNoCut && !(r.mee * 1000.f < c.meeHi)) return false;
    return true;
}
inline Sel MkSel(const Cuts& c) { return [c](const Row& r) { return Pass(c, r); }; }

inline std::string Describe(const Cuts& c) {
    std::string s;
    auto add = [&](const std::string& t) { if (!s.empty()) s += ", "; s += t; };
    if (c.vzMin > kNoCut) add(Form("vz>%.0f", c.vzMin)); else add("no vz cut");
    if (c.oaLab > kNoCut) add(Form("oa_lab<%g", c.oaLab));
    if (c.mpLo > kNoCut) add(Form("Mpipi>%.2f", c.mpLo));
    if (c.mpHi > kNoCut) add(Form("Mpipi<%.2f", c.mpHi));
    if (c.oaLo > kNoCut) add(Form("OA>%g", c.oaLo));
    if (c.oaHi > kNoCut) add(Form("OA<%g", c.oaHi));
    if (c.mm4Hi > kNoCut) add(Form("MM4<%.2f", c.mm4Hi));
    if (c.etaLo > kNoCut) add(Form("eta_rest>%g", c.etaLo));
    if (c.dzPi > kNoCut) add(Form("dz_pi<%g", c.dzPi));
    if (c.vtx == 1) add(">=1 pion in vertex fit");
    if (c.vtx == 2) add("both pions in vertex fit");
    if (c.vtx == 4) add("all 4 tracks in vertex fit");
    if (c.meeLo > kNoCut || c.meeHi > kNoCut) add(Form("%g<=Mee<%g MeV", c.meeLo > kNoCut ? c.meeLo : 0.f, c.meeHi > kNoCut ? c.meeHi : 1e9f));
    if (c.best) add("isBest==1");
    if (c.fold >= 0) add(Form("fold %d", c.fold));
    return s;
}

// ---- the named selections -----------------------------------------------------
// pre = cor_ok && isBest >= 0 && vz > -160 (FAT: vz > -500, PT3)
inline Cuts Pre()       { Cuts c; c.name = "pre"; return c; }
inline Cuts Pre0()      { Cuts c; c.name = "pre without vz>-160"; c.vzMin = kNoCut; return c; }
inline Cuts PeakNom()   { Cuts c; c.name = "peak nominal"; c.oaLab = 55; c.mpLo = 0.30f; return c; }
inline Cuts PeakNomV()  { Cuts c = PeakNom(); c.name = "peak nominal + vertex"; c.vtx = 2; return c; }
inline Cuts BumpNom()   { Cuts c; c.name = "bump nominal"; c.oaLab = 40; c.mpLo = 0.30f; c.mpHi = 0.41f; return c; }
inline Cuts BumpNomV()  { Cuts c = BumpNom(); c.name = "bump nominal + vertex"; c.vtx = 2; return c; }
inline Cuts PreV()      { Cuts c; c.name = "pre + vertex"; c.vtx = 2; return c; }
// frozen selections of the simulation study (FAT_simwt plots/eta_sim/tables/s5_opt_eta_frozen.txt,
// s5_opt_x17_frozen.txt), chosen there on the simulation alone (these data not used in the choice)
inline Cuts Direct()    { Cuts c; c.name = "DIRECT (frozen)"; c.dzPi = 40; return c; }
inline Cuts PeakF()     { Cuts c; c.name = "PEAK (frozen)"; c.dzPi = 40; c.mm4Hi = 2.85f; return c; }
inline Cuts BumpF()     { Cuts c; c.name = "BUMP (frozen) = pre"; return c; }
inline Cuts X17()       { Cuts c; c.name = "X17 (frozen)"; c.vtx = 2; c.meeLo = kX17Lo; c.meeHi = kX17Hi; return c; }
// combinations and old references
inline Cuts PeakNomDz() { Cuts c = PeakNom(); c.name = "peak nominal + dz_pi<40"; c.dzPi = 40; return c; }
inline Cuts BumpNomDz() { Cuts c = BumpNom(); c.name = "bump nominal + dz_pi<40"; c.dzPi = 40; return c; }
inline Cuts Slides()    { Cuts c; c.name = "slides"; c.oaLab = 50; c.mpHi = 0.46f; c.etaLo = 140; c.oaLo = 2; return c; }
inline Cuts OldCode()   { Cuts c; c.name = "old code"; c.oaLab = 50; c.mpHi = 0.42f; c.etaLo = 140; return c; }
inline Cuts WithBest(Cuts c) { c.best = 1; return c; }
inline Cuts WithOA2(Cuts c)  { c.oaLo = (c.oaLo > kNoCut) ? std::max(c.oaLo, 2.f) : 2.f; return c; }
inline Cuts WithFold(Cuts c, int f) { c.fold = f; return c; }

// the simulation peak window moved to the data mass scale: x (data peak 0.5439 / simulation peak 0.5460)
const double kSimPeakMean = 0.5460;
const double kShPeakLo = kSimPeakLo * scan::kPeakMean / kSimPeakMean, kShPeakHi = kSimPeakHi * scan::kPeakMean / kSimPeakMean;
// fill options with the extra windows: 2 = simulation peak window, 3 = simulation bump window,
// 4 = simulation peak window on the data mass scale
inline scan::FillOpt OptWin() {
    scan::FillOpt o;
    o.windows.push_back(std::make_pair(kSimPeakLo, kSimPeakHi));
    o.windows.push_back(std::make_pair(kSimBumpLo, kSimBumpHi));
    o.windows.push_back(std::make_pair(kShPeakLo, kShPeakHi));
    return o;
}

// ============================================================================
// C. eta peak
// ============================================================================
// scale the CB of a filled spectrum by k (sensitivity of S to the CB normalisation)
inline void ScaleCB(Spectra& sp, double k) {
    for (auto& v : sp.cbW) v *= k;
    for (auto& v : sp.cbWVar) v *= k * k;
    for (int i = 0; i < sp.nbin; ++i) {
        sp.cb[i] *= k; sp.cbVar[i] *= k * k;
        sp.sig[i] = sp.n[0][i] - sp.cb[i];
        sp.sigVar[i] = sp.n[0][i] + sp.cbVar[i];
    }
    for (auto& w : sp.win) { w.cb *= k; w.cbVar *= k * k; }
}

struct PeakRes {
    PeakFit f;                  // P_main nominal
    scan::Spread spr;           // model spread (48 variations)
    double Sm3 = 0, Sp3 = 0;    // S with CB x 0.97 and x 1.03
    double Ftot = 0;            // S / sqrt(stat^2 + spread^2)
    double cbs = 0;             // CB / N in the peak window
    scan::WinSum wSim;          // rows in the simulation peak window (if filled with OptWin)
    scan::WinSum wSh;           // rows in the simulation peak window moved to the data mass scale
    bool hasSim = false;
};
inline PeakRes EvalPeak(const Spectra& sp, bool withSpread = true, bool withCB = true) {
    PeakRes r;
    r.f = scan::FitPeakMain(sp);
    if (withSpread) r.spr = scan::PeakModelSpread(sp);
    if (withCB) {
        Spectra a = sp; ScaleCB(a, 0.97); r.Sm3 = scan::FitPeakMain(a).S;
        Spectra b = sp; ScaleCB(b, 1.03); r.Sp3 = scan::FitPeakMain(b).S;
    }
    const double t = std::sqrt(r.f.eS * r.f.eS + r.spr.half * r.spr.half);
    r.Ftot = t > 0 ? r.f.S / t : 0;
    r.cbs = r.f.N > 0 ? r.f.CB / r.f.N : 0;
    if (sp.win.size() > 2) { r.wSim = sp.win[2]; r.hasSim = true; }
    if (sp.win.size() > 4) r.wSh = sp.win[4];
    return r;
}

// efficiency of a cut from pass / fail fits: eps = p / (p + f) (disjoint samples)
inline void Ratio(double p, double ep, double f, double ef, double& eff, double& err) {
    const double t = p + f;
    if (!(std::fabs(t) > 0)) { eff = 0; err = 0; return; }
    eff = p / t;
    err = std::sqrt(f * f * ep * ep + p * p * ef * ef) / (t * t);
}
// efficiency of a subsample p of b (rows of p are part of b)
inline void SubRatio(double p, double ep, double b, double eb, double& eff, double& err) {
    if (!(b > 0)) { eff = 0; err = 0; return; }
    eff = p / b;
    err = std::sqrt(std::max(0.0, (1.0 - 2.0 * eff) * ep * ep + eff * eff * eb * eb)) / b;
}

// ============================================================================
// D. bump: edge step E and tagged eta T (procedures of plots/scan_bump.C)
// ============================================================================
const int kNSl = 6;
const int kNSig = 4;
const double kSlLo[kNSl] = {0.30, 0.34, 0.37, 0.39, 0.40, 0.42};
const double kSlHi[kNSl] = {0.34, 0.37, 0.39, 0.40, 0.41, 0.46};
const double kCtrlEdge = 0.5150;
inline StepCfg SliceCfg(int i, double lo, double hi) {
    StepCfg c;
    c.mpLo = lo; c.mpHi = hi;
    if (i == 3 || i == 4) c.halfRange = 0.035;
    if (i == 5) { c.useEdge = true; c.edge = kCtrlEdge; }
    return c;
}
// slice bounds of a selection (intersection with its M(pipi) cuts)
inline void PassSlice(const Cuts& w, int i, double& lo, double& hi) {
    lo = kSlLo[i]; hi = kSlHi[i];
    if (w.mpLo > kNoCut) lo = std::max(lo, (double)w.mpLo);
    if (w.mpHi > kNoCut) hi = std::min(hi, (double)w.mpHi);
}
// tagged eta: frozen pi0 window in M(e+e-gamma) and eta shape in m5c (scan_bump.C)
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730, kMeeMax = 0.100;
const double kTagMean = 0.5463, kTagSigma = 0.0062, kTagLo = 0.49, kTagHi = 0.62;
const int kTagOrder = 2, kTagNbin = 160, kTagReb = 8;
const double kTagHLo = 0.40, kTagHHi = 0.80;
inline double M5c(const Row& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
inline bool HasGam(const Row& r) { return r.g_idx > 0 && r.mee < kMeeMax; }
inline bool Tagged(const Row& r) { return HasGam(r) && r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi; }
inline PeakCfg TagCfg() {
    PeakCfg c;
    c.mean = kTagMean; c.sigma = kTagSigma; c.lo = kTagLo; c.hi = kTagHi; c.order = kTagOrder;
    c.winLo = c.mean - scan::kWinNsig * c.sigma; c.winHi = c.mean + scan::kWinNsig * c.sigma;
    return c;
}
inline scan::FillOpt TagOpt() {
    scan::FillOpt o;
    o.nbin = kTagNbin; o.lo = kTagHLo; o.hi = kTagHHi; o.reb = kTagReb;
    o.x = [](const Row& r) { return M5c(r); };
    const PeakCfg c = TagCfg();
    o.windows.push_back(std::make_pair(c.winLo, c.winHi));
    return o;
}
// model spread of T (scan_bump.C section C): orders 2, 3 x ranges 0.49-0.62,
// 0.48-0.66, 0.50-0.60 x width -+0.3 MeV; fits with p > 0.05; nominal always counts
inline double TagSpread(const Spectra& s, const PeakFit& nom) {
    double mn = nom.S, mx = nom.S;
    const double rl[3] = {0.49, 0.48, 0.50}, rh[3] = {0.62, 0.66, 0.60};
    for (int ord = 2; ord <= 3; ++ord) for (int ir = 0; ir < 3; ++ir) for (int iw = -1; iw <= 1; ++iw) {
        PeakCfg c = TagCfg(); c.order = ord; c.lo = rl[ir]; c.hi = rh[ir]; c.sigma = kTagSigma + 0.0003 * iw;
        const PeakFit f = scan::FitPeakMain(s, c);
        if (!f.ok || !(f.prob > 0.05)) continue;
        mn = std::min(mn, f.S); mx = std::max(mx, f.S);
    }
    return 0.5 * (mx - mn);
}
// model spread of E (scan_bump.C EdgeSpread): order 2, 3 x half range 35, 45, 55 MeV x
// smearing 4.4, 5.4, 6.4 MeV; per slice half of (max - min) of h over the fits with
// p > 0.05 (nominal always counts); slices added linearly
inline double EdgeSpread(const std::vector<Spectra>& sl, const Cuts& w, int off = 0) {
    double sum = 0;
    for (int i = 0; i < kNSig; ++i) {
        double lo, hi; PassSlice(w, i, lo, hi);
        if (hi - lo < 0.0045 || sl[off + i].nsel[0] < 50) continue;
        const StepFit nom = scan::FitEdgeStep(sl[off + i], SliceCfg(i, lo, hi));
        if (!nom.ok) continue;
        double mn = nom.h, mx = nom.h;
        for (int ord = 2; ord <= 3; ++ord)
            for (int ih = 0; ih < 3; ++ih)
                for (int is = -1; is <= 1; ++is) {
                    StepCfg c = SliceCfg(i, lo, hi);
                    c.order = ord; c.halfRange = 0.035 + 0.010 * ih; c.s = scan::kEdgeSigma + 0.001 * is;
                    const StepFit f = scan::FitEdgeStep(sl[off + i], c);
                    if (!f.ok || !(TMath::Prob(f.chi2, f.ndf) > 0.05)) continue;
                    mn = std::min(mn, f.h); mx = std::max(mx, f.h);
                }
        sum += 0.5 * (mx - mn);
    }
    return sum;
}

struct BumpRes {
    bool ok = false;
    double E = 0, eE = 0, spreadE = 0;   // sum of step heights (rows / 2.5 MeV)
    double Em3 = 0, Ep3 = 0;             // E with CB x 0.97, x 1.03
    double h[kNSl], eh[kNSl], dchi2[kNSl]; bool okSl[kNSl];
    double chi2 = 0; int ndf = 0;
    scan::WinSum win, winSim;            // data bump window 0.4550-0.5225; simulation window 0.4576-0.5255
    double Q = 0;                        // E / sqrt(N in the bump window)
    PeakFit T; double spreadT = 0;       // tagged eta (if computed)
    bool hasT = false;
};

// E (+ spread) for a list of selections; slices filled in one pass.
// fullSel[j] = spectra of the whole selection (for the bump window), may be filled here.
inline std::vector<BumpRes> EvalBump(const scan::Data& d, const std::vector<Cuts>& cs, bool withSpread,
                                     bool withT, const scan::Data* dT = nullptr, bool withCB = false) {
    const int nc = (int)cs.size();
    std::vector<BumpRes> out(nc);
    std::vector<Sel> sels;
    for (int j = 0; j < nc; ++j) {
        const Cuts w = cs[j];
        for (int i = 0; i < kNSl; ++i) {
            const float lo = kSlLo[i], hi = kSlHi[i];
            sels.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && Pass(w, r); });
        }
        sels.push_back(MkSel(w));
    }
    const std::vector<Spectra> sp = scan::FillMany(d, sels, OptWin());
    const int per = kNSl + 1;
    for (int j = 0; j < nc; ++j) {
        BumpRes& R = out[j];
        const Cuts& w = cs[j];
        double v = 0;
        for (int i = 0; i < kNSl; ++i) {
            R.okSl[i] = false; R.h[i] = R.eh[i] = R.dchi2[i] = 0;
            double lo, hi; PassSlice(w, i, lo, hi);
            if (hi - lo < 0.0045 || sp[j * per + i].nsel[0] < 50) continue;
            const StepFit f = scan::FitEdgeStep(sp[j * per + i], SliceCfg(i, lo, hi));
            if (!f.ok) continue;
            R.okSl[i] = true; R.h[i] = f.h; R.eh[i] = f.eh; R.dchi2[i] = f.dchi2;
            if (i >= kNSig) continue;
            R.E += f.h; v += f.eh * f.eh; R.chi2 += f.chi2; R.ndf += f.ndf;
        }
        R.eE = std::sqrt(v);
        R.ok = true;
        const Spectra& full = sp[j * per + kNSl];
        R.win = full.win[1]; R.winSim = full.win[3];
        if (R.win.n[0] > 0) R.Q = R.E / std::sqrt(R.win.n[0]);
        if (withSpread) {
            std::vector<Spectra> sl(sp.begin() + j * per, sp.begin() + j * per + kNSig);
            R.spreadE = EdgeSpread(sl, w);
        }
        if (withCB) {
            for (int q = 0; q < 2; ++q) {
                const double k = q == 0 ? 0.97 : 1.03;
                double e = 0;
                for (int i = 0; i < kNSig; ++i) {
                    double lo, hi; PassSlice(w, i, lo, hi);
                    if (hi - lo < 0.0045 || sp[j * per + i].nsel[0] < 50) continue;
                    Spectra s = sp[j * per + i]; ScaleCB(s, k);
                    const StepFit f = scan::FitEdgeStep(s, SliceCfg(i, lo, hi));
                    if (f.ok) e += f.h;
                }
                (q == 0 ? R.Em3 : R.Ep3) = e;
            }
        }
    }
    if (withT && dT) {
        std::vector<Sel> ts;
        for (int j = 0; j < nc; ++j) { const Cuts w = cs[j]; ts.push_back([=](const Row& r) { return Tagged(r) && Pass(w, r); }); }
        const std::vector<Spectra> st = scan::FillMany(*dT, ts, TagOpt());
        for (int j = 0; j < nc; ++j) {
            out[j].T = scan::FitPeakMain(st[j], TagCfg());
            out[j].spreadT = withSpread ? TagSpread(st[j], out[j].T) : 0;
            out[j].hasT = true;
        }
    }
    return out;
}

// copy of the rows that pass a filter
inline void Filter(const scan::Data& in, scan::Data& out, const std::function<bool(const Row&)>& keep) {
    for (int k = 0; k < 3; ++k) {
        out.s[k].clear();
        for (const auto& r : in.s[k]) if (keep(r)) out.s[k].push_back(r);
    }
    out.fallback = in.fallback;
}

// ============================================================================
// E. drawing
// ============================================================================
inline void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gStyle->SetLabelSize(0.050, "xyz"); gStyle->SetTitleSize(0.056, "xyz");
    gStyle->SetTitleOffset(1.05, "x"); gStyle->SetTitleOffset(1.25, "y");
    gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.08); gStyle->SetPadBottomMargin(0.14);
    gStyle->SetEndErrorSize(0);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetHatchesLineWidth(1);
}
// vector PDF cropped to the canvas (EPS -> ps2pdf -dEPSCrop) and a PNG preview
inline void Save(TCanvas* c, const char* name) {
    const std::string n = name;
    const std::string sc = kPlainScope;                                              // " name=outbase ... "
    const size_t p = sc.find(" " + n + "=");
    const std::string eb = std::string(kPlainScratch) + "/" + n;                     // EPS always in the scratch folder
    std::string b = eb;
    if (p != std::string::npos) { const size_t q = p + n.size() + 2; b = std::string(kFigOut) + "/" + sc.substr(q, sc.find(' ', q) - q); }
    c->SaveAs((eb + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf", eb.c_str(), b.c_str()));
    c->SaveAs((b + ".png").c_str());
}
inline TLatex* Txt(double x, double y, const char* t, double size = 0.050, int col = kBlack, int align = 11) {
    TLatex* l = new TLatex(x, y, t); l->SetNDC(); l->SetTextFont(42); l->SetTextSize(size);
    l->SetTextColor(col); l->SetTextAlign(align); l->Draw(); return l;
}
inline TLine* HLine(double x1, double x2, double y, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid); l->Draw(); return l;
}
inline TLine* VLine(double x, double y1, double y2, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x, y1, x, y2); l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid); l->Draw(); return l;
}
inline TGraphErrors* GE(int col, int mk, double ms = 1.2) {
    TGraphErrors* g = new TGraphErrors(); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk);
    g->SetMarkerSize(ms); g->SetLineWidth(2); return g;
}
inline void AddPt(TGraphErrors* g, double x, double y, double ey, double ex = 0) {
    const int n = g->GetN(); g->SetPoint(n, x, y); g->SetPointError(n, ex, ey);
}
inline TLegend* Leg(double x1, double y1, double x2, double y2, double size = 0.045) {
    TLegend* l = new TLegend(x1, y1, x2, y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(42); l->SetTextSize(size);
    return l;
}
// empty frame
inline TH1D* Frame(const char* name, double x1, double x2, double y1, double y2, const char* xt, const char* yt) {
    TH1D* h = new TH1D(name, "", 100, x1, x2); h->SetDirectory(nullptr);
    h->SetMinimum(y1); h->SetMaximum(y2);
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->CenterTitle(false);
    h->Draw("AXIS");
    return h;
}

// ---- printing to stdout and a table file at the same time ------------------------
struct Out {
    FILE* f = nullptr;
    void Open(const std::string& path) { Close(); f = fopen(path.c_str(), "w"); }
    void Close() { if (f) fclose(f); f = nullptr; }
    void P(const char* fmt, ...) {
        va_list a, b; va_start(a, fmt); va_copy(b, a);
        vprintf(fmt, a); if (f) vfprintf(f, fmt, b);
        va_end(a); va_end(b);
    }
};

} // namespace opt

#endif // OPT_COMMON_H
