// ============================================================================
// scan_bump.C  --  the bump below the eta mass in M(pi+pi-e+e-) and the cuts
// ============================================================================
// HADES p+p, T = 4.53 GeV, channel pi+ pi- e+ e-. The bump is eta -> pi+pi-pi0
// with an e+e- pair from the pi0 and the photon not used. The pair can come from
// the Dalitz decay pi0 -> e+e-gamma or from pi0 -> gamma gamma with one photon
// converted: both give the same M(e+e-gamma) and the same kinematic edge, and
// this study cannot tell them apart. The bump has no model-free yield, so three
// other measures are followed through every cut:
//
//   E  kinematic-edge step. For a fixed M(pi+pi-) the four-particle mass of a
//      bump event cannot exceed M4max(M(pipi)). In M(pipi) slices the spectrum
//      N - CB drops at this predicted place; the height h of the drop (rows
//      per 2.5 MeV, engine fit scan::FitEdgeStep) is fitted with the position
//      FIXED. E = sum of h over the slices 0.30-0.34, 0.34-0.37, 0.37-0.39,
//      0.39-0.40. Not in E: 0.40-0.41 (the bump is narrower than the fit range
//      there, no step can be separated) and the control slice 0.42-0.46.
//   T  ECAL-tagged eta. Rows with a photon whose M(e+e-gamma) lies in the pi0
//      window show an eta peak in the five-particle mass. The mass used is
//      m5c: M(pi+pi-e+e-gamma) with the photon energy rescaled so that
//      M(e+e-gamma) = m_pi0 (in the plain m5 the peak moves by 1 MeV per MeV
//      of M(e+e-gamma)). Efficiency of a cut for bump events = peak content
//      after the cut / before the cut.
//   W  N - CB and the CB share in the bump window of the engine (0.455-0.5225).
//
// Figure of merit: Q = E / sqrt(N), N = all rows in the bump window. gain of a
// cut = efficiency / sqrt(share of N kept), from E and from T. (Review: the
// first version divided by sqrt(share of N - CB kept); the noise is sqrt(N),
// and N does not depend on the CB estimate. scan_extra.C and scan_ribbons.C use
// the same definition.)
// W = N - CB: CB = 2 sqrt(N++ N--) has no k factor and the like-sign samples are
// very asymmetric (N--/N++ about 6): W carries an uncertainty of several per
// cent (the like-sign sum N++ + N-- is printed as the other extreme).
// E is a sum of step heights (rows per 2.5 MeV), NOT a yield.
// The kinematic gain of the candidates (about 10 %) rests on the tagged sample
// T; from E it lies inside the model spread. Its optimisation bias is measured
// with two folds in bump_candidates.txt (section I).
//
// Steps: 1 pi0 window and tagged sample, 2 reference edge fits (variations,
// free edge position, lower edge), 3 one-dimensional scans from "no cut" with
// two folds, 4 two-dimensional scans, 5 candidates (loose, nominal, pure;
// slides, old code), cuts varied around the nominal selection, isBest == 1,
// cor_ok fallback, 6 ROOT file, 7 figures.
//
// Output
//   tables   <prod>/scan/bump_pi0.txt, bump_reference.txt, bump_scan1d.txt,
//            bump_scan2d.txt, bump_candidates.txt
//   spectra  <prod>/scan/bump_scan.root
//   figures  plots/output/eta_scan/bump_*.png (.pdf for the main ones)
//
// Run (ACLiC, build directory outside the repository), about 35 s:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/aclic_bump", kTRUE)' 'plots/scan_bump.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "ScanEngine.h"

namespace bump {

using scan::Data;
using scan::Row;
using scan::Sel;
using scan::Spectra;
using scan::StepCfg;
using scan::StepFit;
using scan::PeakCfg;
using scan::PeakFit;
using scan::WinSum;

// ============================================================================
// A. Constants of this study
// ============================================================================

// ---- M(pi+pi-) slices of the edge fit -----------------------------------------
// 0..3 signal slices (their sum is E), 4 = 0.40-0.41 (reported, not in E),
// 5 = control slice above the kinematic limit m_eta - m_pi0 (edge set by hand)
const int kNSl = 6;
const int kNSig = 4;
const double kSlLo[kNSl] = {0.30, 0.34, 0.37, 0.39, 0.40, 0.42};
const double kSlHi[kNSl] = {0.34, 0.37, 0.39, 0.40, 0.41, 0.46};
const char* const kSlName[kNSl] = {"0.30-0.34", "0.34-0.37", "0.37-0.39", "0.39-0.40", "0.40-0.41", "0.42-0.46"};
const double kCtrlEdge = 0.5150;
// Fit range: edge -+ 45 MeV (engine default) for the first three slices.
// From 0.39 on the non-bump background is a threshold turn-on and the long
// range fails (chi2/ndf 68/38 for 0.39-0.40): there edge -+ 35 MeV, which ends
// below the eta peak (0.4625-0.5275, no eta term).
inline StepCfg SliceCfg(int i, double lo, double hi) {
    StepCfg c;
    c.mpLo = lo; c.mpHi = hi;
    if (i == 3 || i == 4) c.halfRange = 0.035;
    if (i == 5) { c.useEdge = true; c.edge = kCtrlEdge; }
    return c;
}

// ---- pi0 window and eta shape of the tagged sample (frozen; see bump_pi0.txt) --
// pi0 in M(e+e-gamma) of tagged eta events: position kPi0Pos, width kPi0Sig
// (Gaussian fit of the eta yield against M(e+e-gamma)); window = +-2 sigma.
const double kPi0Pos = 0.1459;
const double kPi0Sig = 0.0133;
const double kPi0Lo  = 0.1190;   // kPi0Pos - 2 kPi0Sig, rounded to 1 MeV
const double kPi0Hi  = 0.1730;   // kPi0Pos + 2 kPi0Sig
const double kMeeMax = 0.100;    // pairs from a pi0 (m_ee < m_pi0); keeps the photon rescaling factor finite
// eta peak in the constrained mass m5c for this window (free fit, pol2, 0.49-0.62)
const double kTagMean  = 0.5463;
const double kTagSigma = 0.0062;
const double kTagLo = 0.49, kTagHi = 0.62;
const int    kTagOrder = 2;
// m5c histogram: 2.5 MeV bins 0.40-0.80, CB in 20 MeV bins
const int    kTagNbin = 160;
const double kTagHLo = 0.40, kTagHHi = 0.80;
const int    kTagReb = 8;
// old production: plain m5, window 0.125-0.145, shape of the free fit there
const double kOldLo = 0.125, kOldHi = 0.145;
const double kOldMean = 0.5489, kOldSigma = 0.0088;

// five-particle mass with the photon energy rescaled so that M(e+e-gamma) = m_pi0:
// (P4 + f k)^2 = M4^2 + f (m5^2 - M4^2),  f = (m_pi0^2 - m_ee^2) / (M(eeg)^2 - m_ee^2)
inline double M5c(const Row& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
inline bool HasGam(const Row& r) { return r.g_idx > 0 && r.mee < kMeeMax; }
inline bool Tagged(const Row& r) { return HasGam(r) && r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi; }

// lower kinematic limit on the measured scale
inline double LowEdgeMeasured(double mpipi) {
    const double e = scan::M4min(mpipi / scan::kScalePiPi);
    return e > 0 ? scan::kScale * e : -1.0;
}

// ============================================================================
// B. Working point (a set of cuts)
// ============================================================================
struct WP {
    std::string name;
    double oaLo = -1;      // e+e- opening angle > oaLo [deg]          (-1 = no cut)
    double oaHi = -1;      // e+e- opening angle < oaHi [deg]
    double oaLab = -1;     // angle (pi+pi-),(e+e-) in the lab < oaLab [deg]
    double mpLo = -1;      // M(pi+pi-) > mpLo [GeV/c2]
    double mpHi = -1;      // M(pi+pi-) < mpHi
    double mm4 = -1;       // missing mass of the four particles < mm4 [GeV/c2]
    double etaRest = -1;   // angle (pi+pi-),(e+e-) in the frame forced to m_eta > etaRest [deg]
    int best = 0;          // 1: rows with isBest == 1 only
};
inline bool Pass(const WP& w, const Row& r) {
    if (w.best && r.isBest != 1) return false;
    if (w.mpLo >= 0 && !(r.mpipi > w.mpLo)) return false;
    if (w.mpHi >= 0 && !(r.mpipi < w.mpHi)) return false;
    if (w.oaLab >= 0 && !(r.oa_lab < w.oaLab)) return false;
    if (w.oaLo >= 0 && !(r.oa > w.oaLo)) return false;
    if (w.oaHi >= 0 && !(r.oa < w.oaHi)) return false;
    if (w.mm4 >= 0 && !(r.mm4 < w.mm4)) return false;
    if (w.etaRest >= 0 && !(r.oa_eta > w.etaRest)) return false;
    return true;
}
inline WP None() { WP w; w.name = "none"; return w; }
inline WP Slides() { WP w; w.name = "slides"; w.oaLab = 50; w.mpHi = 0.46; w.etaRest = 140; w.oaLo = 2; return w; }
inline WP OldCode() { WP w; w.name = "old code"; w.oaLab = 50; w.mpHi = 0.42; w.etaRest = 140; return w; }
inline WP OldCodeOA() { WP w = OldCode(); w.name = "old code + OA 2-6"; w.oaLo = 2; w.oaHi = 6; return w; }
// candidates of the bump selection (chosen from the scans, see bump_candidates.txt)
inline WP Loose()   { WP w; w.name = "loose";   w.oaLab = 50; w.mpHi = 0.42; return w; }
inline WP Nominal() { WP w; w.name = "nominal"; w.oaLab = 40; w.mpLo = 0.30; w.mpHi = 0.41; return w; }
inline WP Pure()    { WP w; w.name = "pure";    w.oaLab = 35; w.mpLo = 0.30; w.mpHi = 0.40; w.oaHi = 20; return w; }
std::string Describe(const WP& w) {
    std::string s;
    auto add = [&](const std::string& t) { if (!s.empty()) s += ", "; s += t; };
    if (w.oaLab >= 0) add(Form("oa_lab<%.0f", w.oaLab));
    if (w.mpLo >= 0) add(Form("Mpipi>%.2f", w.mpLo));
    if (w.mpHi >= 0) add(Form("Mpipi<%.2f", w.mpHi));
    if (w.oaLo >= 0) add(Form("OA>%.1f", w.oaLo));
    if (w.oaHi >= 0) add(Form("OA<%.0f", w.oaHi));
    if (w.mm4 >= 0) add(Form("mm4<%.2f", w.mm4));
    if (w.etaRest >= 0) add(Form("eta_rest>%.0f", w.etaRest));
    if (w.best) add("isBest==1");
    if (s.empty()) s = "no cut";
    return s;
}

// the seven scanned cuts
enum CutId { kOaLo = 0, kOaHi, kOaLab, kMpHi, kMpLo, kMm4, kEtaRest, kNCut };
const char* const kCutName[kNCut] = {"OA lower", "OA upper", "oa_lab upper", "M(pipi) upper", "M(pipi) lower",
                                     "mm4 upper", "eta-rest lower"};
const char* const kCutTag[kNCut] = {"oalo", "oahi", "oalab", "mphi", "mplo", "mm4", "etarest"};
const char* const kCutAxis[kNCut] = {"lower cut on OA(e^{+}e^{-}) [deg]", "upper cut on OA(e^{+}e^{-}) [deg]",
                                     "upper cut on oa_lab((#pi#pi),(ee)) [deg]",
                                     "upper cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]",
                                     "lower cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]",
                                     "upper cut on the missing mass of #pi#piee [GeV/c^{2}]",
                                     "lower cut on the #eta-rest angle [deg] (flag only)"};
// thresholds; -1 = no cut
inline std::vector<double> Thresholds(int c) {
    std::vector<double> v;
    switch (c) {
        case kOaLo:    for (int i = 0; i <= 8; ++i) v.push_back(0.5 * i); break;
        case kOaHi:    v = {4, 6, 8, 10, 15, 20, 30, -1}; break;
        case kOaLab:   for (int i = 25; i <= 70; i += 5) v.push_back(i); v.push_back(-1); break;
        case kMpHi:    for (int i = 38; i <= 43; ++i) v.push_back(0.01 * i); v.push_back(-1); break;
        case kMpLo:    v = {-1, 0.28, 0.30, 0.32, 0.34}; break;
        case kMm4:     for (int i = 0; i <= 7; ++i) v.push_back(2.60 + 0.05 * i); v.push_back(-1); break;
        case kEtaRest: v = {-1, 100, 120, 140}; break;
    }
    return v;
}
inline double& CutRef(WP& w, int c) {
    switch (c) {
        case kOaLo: return w.oaLo;
        case kOaHi: return w.oaHi;
        case kOaLab: return w.oaLab;
        case kMpHi: return w.mpHi;
        case kMpLo: return w.mpLo;
        case kMm4: return w.mm4;
        default: return w.etaRest;
    }
}
inline std::string ThrName(int c, double t) {
    if (t < 0 || (c == kOaLo && t == 0)) return "none";
    if (c == kMpHi || c == kMpLo || c == kMm4) return Form("%.2f", t);
    if (c == kOaLo) return Form("%.1f", t);
    return Form("%.0f", t);
}

// ============================================================================
// C. Measurement of a list of working points
// ============================================================================
struct Job {
    WP w;                 // the selection
    WP base;              // reference of the efficiencies: pass = w, fail = base and not w
    int fold = -1;        // -1 all event groups, 0 / 1 one fold (group parity)
    bool doT = true;      // tagged-eta fits
    bool doFail = true;   // edge fits of the failing rows
};
struct Res {
    // edge step
    bool okP[kNSl], okF[kNSig];
    StepFit fp[kNSl], ff[kNSig];
    double E = 0, eE = 0;           // sum of the step heights of the signal slices, selection
    double Ef = 0, eEf = 0;         // the same for the rows of the base that fail
    double E0 = 0, eE0 = 0;         // E of the base selection
    double effE = 0, eEffE = 0;     // E / E0 (the error takes the common rows into account)
    double closeE = 0;              // (E + Ef) / E0: 1 if the fits of the two parts add up
    double chi2 = 0; int ndf = 0;   // summed over the signal slices
    double bkg = 0;                 // summed polynomial at the edges
    // bump window
    WinSum win, win0;               // selection, base
    double effW = 0;                // (N - CB) / (N - CB of the base)
    double effN = 0;                // N / N of the base: rows kept in the bump window
    // tagged eta
    PeakFit tp, tf, t0;             // passing rows, failing rows, base
    double effT = 0, eEffT = 0;     // T_pass / T_base (the error takes the common rows into account)
    double closeT = 0;              // (T_pass + T_fail) / T_base: 1 if the fits of the two parts add up
    // figures of merit
    double Q = 0, eQ = 0;           // E / sqrt(N)
    double gainE = 0, eGainE = 0;   // effE / sqrt(effN)
    double gainT = 0, eGainT = 0;   // effT / sqrt(effN)
};

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
// share of two disjoint samples: p / (p + f)
inline void Ratio(double p, double ep, double f, double ef, double& eff, double& err) {
    const double t = p + f;
    if (!(t > 0)) { eff = 0; err = 0; return; }
    eff = p / t;
    err = std::sqrt(f * f * ep * ep + p * p * ef * ef) / (t * t);
}
// efficiency of a subsample: p / b, where the rows of p are part of b.
// With b = p + f (independent parts): var = [(1 - 2 eff) ep^2 + eff^2 eb^2] / b^2.
inline void SubRatio(double p, double ep, double b, double eb, double& eff, double& err) {
    if (!(b > 0)) { eff = 0; err = 0; return; }
    eff = p / b;
    err = std::sqrt(std::max(0.0, (1.0 - 2.0 * eff) * ep * ep + eff * eff * eb * eb)) / b;
}

// slice bounds of the selection (intersection with its M(pipi) cuts)
inline void PassSlice(const WP& w, int i, double& lo, double& hi) {
    lo = kSlLo[i]; hi = kSlHi[i];
    if (w.mpLo >= 0) lo = std::max(lo, w.mpLo);
    if (w.mpHi >= 0) hi = std::min(hi, w.mpHi);
}

// dE: rows for the M4 spectra; dT: rows with a photon (may be null if no job has doT)
std::vector<Res> Measure(const Data& dE, const Data* dT, const std::vector<Job>& jobs) {
    const int nj = (int)jobs.size();
    std::vector<Res> out(nj);
    // ---- the distinct base selections: E and the tagged peak of each ----
    std::vector<std::string> bkey;
    std::vector<int> bjob, bof(nj, 0);
    for (int j = 0; j < nj; ++j) {
        const std::string k = Describe(jobs[j].base) + Form(" fold %d", jobs[j].fold);
        size_t q = 0;
        while (q < bkey.size() && bkey[q] != k) ++q;
        if (q == bkey.size()) { bkey.push_back(k); bjob.push_back(j); }
        bof[j] = (int)q;
    }
    std::vector<double> bE(bkey.size(), 0.), bEe(bkey.size(), 0.);
    std::vector<PeakFit> bT(bkey.size());
    {
        std::vector<Sel> sels, ts;
        for (size_t q = 0; q < bkey.size(); ++q) {
            const WP b = jobs[bjob[q]].base;
            const int fold = jobs[bjob[q]].fold;
            for (int i = 0; i < kNSig; ++i) {
                const float lo = kSlLo[i], hi = kSlHi[i];
                sels.push_back([=](const Row& r) {
                    return r.mpipi >= lo && r.mpipi < hi && (fold < 0 || r.fold == fold) && Pass(b, r); });
            }
            ts.push_back([=](const Row& r) { return Tagged(r) && (fold < 0 || r.fold == fold) && Pass(b, r); });
        }
        const std::vector<Spectra> sp = scan::FillMany(dE, sels);
        for (size_t q = 0; q < bkey.size(); ++q) {
            const WP& b = jobs[bjob[q]].base;
            double v = 0;
            for (int i = 0; i < kNSig; ++i) {
                double lo, hi; PassSlice(b, i, lo, hi);
                if (hi - lo < 0.0045 || sp[q * kNSig + i].nsel[0] < 50) continue;
                const StepFit f = scan::FitEdgeStep(sp[q * kNSig + i], SliceCfg(i, lo, hi));
                if (!f.ok) continue;
                bE[q] += f.h; v += f.eh * f.eh;
            }
            bEe[q] = std::sqrt(v);
        }
        if (dT) {
            const std::vector<Spectra> st = scan::FillMany(*dT, ts, TagOpt());
            for (size_t q = 0; q < bkey.size(); ++q) bT[q] = scan::FitPeakMain(st[q], TagCfg());
        }
    }
    const int perJob = kNSl + kNSig + 2;
    const int batch = 24;
    for (int j0 = 0; j0 < nj; j0 += batch) {
        const int j1 = std::min(nj, j0 + batch);
        // ---- M4 spectra ----
        std::vector<Sel> sels;
        for (int j = j0; j < j1; ++j) {
            const WP w = jobs[j].w, b = jobs[j].base;
            const int fold = jobs[j].fold;
            for (int i = 0; i < kNSl; ++i) {
                const float lo = kSlLo[i], hi = kSlHi[i];
                sels.push_back([=](const Row& r) {
                    return r.mpipi >= lo && r.mpipi < hi && (fold < 0 || r.fold == fold) && Pass(w, r); });
            }
            for (int i = 0; i < kNSig; ++i) {
                const float lo = kSlLo[i], hi = kSlHi[i];
                const bool doFail = jobs[j].doFail;
                sels.push_back([=](const Row& r) {
                    return doFail && r.mpipi >= lo && r.mpipi < hi && (fold < 0 || r.fold == fold) &&
                           Pass(b, r) && !Pass(w, r); });
            }
            sels.push_back([=](const Row& r) { return (fold < 0 || r.fold == fold) && Pass(w, r); });
            sels.push_back([=](const Row& r) { return (fold < 0 || r.fold == fold) && Pass(b, r); });
        }
        const std::vector<Spectra> sp = scan::FillMany(dE, sels);
        for (int j = j0; j < j1; ++j) {
            Res& R = out[j];
            const WP& w = jobs[j].w;
            const Spectra* s = &sp[(j - j0) * perJob];
            double v = 0, vf = 0;
            for (int i = 0; i < kNSl; ++i) {
                double lo, hi; PassSlice(w, i, lo, hi);
                R.okP[i] = false;
                if (hi - lo < 0.0045 || s[i].nsel[0] < 50) continue;
                R.fp[i] = scan::FitEdgeStep(s[i], SliceCfg(i, lo, hi));
                R.okP[i] = R.fp[i].ok;
                if (!R.okP[i] || i >= kNSig) continue;
                R.E += R.fp[i].h; v += R.fp[i].eh * R.fp[i].eh;
                R.chi2 += R.fp[i].chi2; R.ndf += R.fp[i].ndf; R.bkg += R.fp[i].bkgAtEdge;
            }
            R.eE = std::sqrt(v);
            for (int i = 0; i < kNSig; ++i) {
                R.okF[i] = false;
                if (!jobs[j].doFail || s[kNSl + i].nsel[0] < 50) continue;
                // template of the failing rows: the part of the slice outside the M(pipi) cut
                double lo, hi; PassSlice(w, i, lo, hi);
                double flo = kSlLo[i], fhi = kSlHi[i];
                if (hi - lo >= 0.0045) {
                    if (hi < kSlHi[i] - 1e-9 && kSlHi[i] - hi >= 0.0045) flo = hi;
                    else if (lo > kSlLo[i] + 1e-9 && lo - kSlLo[i] >= 0.0045) fhi = lo;
                }
                R.ff[i] = scan::FitEdgeStep(s[kNSl + i], SliceCfg(i, flo, fhi));
                R.okF[i] = R.ff[i].ok;
                if (!R.okF[i]) continue;
                R.Ef += R.ff[i].h; vf += R.ff[i].eh * R.ff[i].eh;
            }
            R.eEf = std::sqrt(vf);
            R.E0 = bE[bof[j]]; R.eE0 = bEe[bof[j]];
            SubRatio(R.E, R.eE, R.E0, R.eE0, R.effE, R.eEffE);
            R.closeE = R.E0 != 0 ? (R.E + R.Ef) / R.E0 : 0;
            R.win = s[kNSl + kNSig].win[1];
            R.win0 = s[kNSl + kNSig + 1].win[1];
            R.effW = R.win0.sig() > 0 ? R.win.sig() / R.win0.sig() : 0;
            R.effN = R.win0.n[0] > 0 ? R.win.n[0] / R.win0.n[0] : 0;
            if (R.win.n[0] > 0) { R.Q = R.E / std::sqrt(R.win.n[0]); R.eQ = R.eE / std::sqrt(R.win.n[0]); }
            if (R.effN > 0) { R.gainE = R.effE / std::sqrt(R.effN); R.eGainE = R.eEffE / std::sqrt(R.effN); }
        }
        // ---- tagged eta ----
        if (!dT) continue;
        std::vector<Sel> ts;
        std::vector<int> owner;
        for (int j = j0; j < j1; ++j) {
            if (!jobs[j].doT) continue;
            const WP w = jobs[j].w, b = jobs[j].base;
            const int fold = jobs[j].fold;
            ts.push_back([=](const Row& r) { return Tagged(r) && (fold < 0 || r.fold == fold) && Pass(w, r); });
            ts.push_back([=](const Row& r) {
                return Tagged(r) && (fold < 0 || r.fold == fold) && Pass(b, r) && !Pass(w, r); });
            owner.push_back(j);
        }
        if (ts.empty()) continue;
        const std::vector<Spectra> st = scan::FillMany(*dT, ts, TagOpt());
        const PeakCfg tc = TagCfg();
        for (size_t q = 0; q < owner.size(); ++q) {
            Res& R = out[owner[q]];
            R.tp = scan::FitPeakMain(st[2 * q], tc);
            R.tf = scan::FitPeakMain(st[2 * q + 1], tc);
            const double sf = (st[2 * q + 1].nsel[0] > 0 && R.tf.ok) ? R.tf.S : 0.0;
            const double ef = (st[2 * q + 1].nsel[0] > 0 && R.tf.ok) ? R.tf.eS : 0.0;
            R.t0 = bT[bof[owner[q]]];
            if (R.t0.ok && R.t0.S > 0) {
                R.closeT = (R.tp.S + sf) / R.t0.S;
                SubRatio(R.tp.S, R.tp.eS, R.t0.S, R.t0.eS, R.effT, R.eEffT);
            }
            if (R.effN > 0) { R.gainT = R.effT / std::sqrt(R.effN); R.eGainT = R.eEffT / std::sqrt(R.effN); }
        }
    }
    return out;
}

// ---- model variations of the edge fit ------------------------------------------
// order 2, 3; half range 35, 45, 55 MeV; smearing 4.4, 5.4, 6.4 MeV (18 fits per slice).
// Per slice: half of (max - min) of h over the fits with p > 0.05 (the nominal fit
// always counts). For E the slice values are added linearly (the variations act
// in the same direction in all slices).
struct ESpread { double half[kNSig]; int kept[kNSig]; double sum = 0; };
ESpread EdgeSpread(const std::vector<Spectra>& sl, const WP& w) {
    ESpread s;
    for (int i = 0; i < kNSig; ++i) {
        s.half[i] = 0; s.kept[i] = 0;
        double lo, hi; PassSlice(w, i, lo, hi);
        if (hi - lo < 0.0045) continue;
        const StepFit nom = scan::FitEdgeStep(sl[i], SliceCfg(i, lo, hi));
        if (!nom.ok) continue;
        double mn = nom.h, mx = nom.h;
        for (int ord = 2; ord <= 3; ++ord)
            for (int ih = 0; ih < 3; ++ih)
                for (int is = -1; is <= 1; ++is) {
                    StepCfg c = SliceCfg(i, lo, hi);
                    c.order = ord; c.halfRange = 0.035 + 0.010 * ih; c.s = scan::kEdgeSigma + 0.001 * is;
                    const StepFit f = scan::FitEdgeStep(sl[i], c);
                    if (!f.ok || !(TMath::Prob(f.chi2, f.ndf) > 0.05)) continue;
                    s.kept[i]++;
                    mn = std::min(mn, f.h); mx = std::max(mx, f.h);
                }
        s.half[i] = 0.5 * (mx - mn);
        s.sum += s.half[i];
    }
    return s;
}

// ---- printing --------------------------------------------------------------------
FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list a, b;
    va_start(a, fmt);
    va_copy(b, a);
    vprintf(fmt, a);
    if (gOut) vfprintf(gOut, fmt, b);
    va_end(a);
    va_end(b);
}
void Open(const char* name) {
    if (gOut) fclose(gOut);
    gOut = fopen(Form("%s/%s", scan::kScanDir, name), "w");
}
void Close() { if (gOut) fclose(gOut); gOut = nullptr; }

void Filter(const Data& in, Data& out, const std::function<bool(const Row&)>& keep) {
    for (int k = 0; k < 3; ++k) {
        out.s[k].clear();
        for (const auto& r : in.s[k]) if (keep(r)) out.s[k].push_back(r);
        out.s[k].shrink_to_fit();
    }
    out.fallback = in.fallback;
}

// header and one line of a scan table
void ScanHead() {
    P("   %-8s | %6s %5s %8s | %6s %5s %5s | %6s %5s | %6s %5s %6s %5s %11s | %6s %5s %5s | %7s %6s %7s %5s %6s %6s %6s | %5s %5s | %6s %5s %6s %5s | %s\n",
      "thresh", "E", "err", "chi2/ndf", "eff_E", "err", "close", "E_fail", "err", "T_pass", "err", "T_fail", "err", "chi2 p / f",
      "eff_T", "err", "close", "N", "CB", "W=N-CB", "err", "CB/N", "eff_W", "eff_N", "Q", "err", "gain_E", "err", "gain_T", "err",
      "h per slice (0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 | 0.40-0.41 | control)");
}
void ScanLine(const char* thr, const Res& R) {
    P("   %-8s | %6.0f %5.0f %5.0f/%-3d | %6.3f %5.3f %5.3f | %6.0f %5.0f | %6.0f %5.0f %6.0f %5.0f %5.0f/%-5.0f | %6.3f %5.3f %5.3f | %7.0f %6.0f %7.0f %5.0f %6.3f %6.3f %6.3f | %5.2f %5.2f | %6.3f %5.3f %6.3f %5.3f |",
      thr, R.E, R.eE, R.chi2, R.ndf, R.effE, R.eEffE, R.closeE, R.Ef, R.eEf, R.tp.S, R.tp.eS, R.tf.S, R.tf.eS, R.tp.chi2, R.tf.chi2,
      R.effT, R.eEffT, R.closeT,
      R.win.n[0], R.win.cb, R.win.sig(), R.win.sigErr(), R.win.cbShare(), R.effW, R.effN, R.Q, R.eQ,
      R.gainE, R.eGainE, R.gainT, R.eGainT);
    for (int i = 0; i < kNSl; ++i) {
        if (i == kNSig || i == kNSig + 1) P(" |");
        if (R.okP[i]) P(" %4.0f+-%-3.0f", R.fp[i].h, R.fp[i].eh); else P(" %8s", "-");
    }
    P("\n");
}

// ---- free edge position: profile of chi2 over the mass scale of the edge ----------
// The fit range and the model are those of the nominal fit; only the factor between
// the true limit and the measured edge is varied (nominal scan::kScale).
struct EdgeFree {
    bool ok = false;
    double scale = 0, eLo = 0, eHi = 0;    // best factor and its -+ errors (delta chi2 = 1)
    double chi2min = 0, chi2nom = 0; int ndf = 0;
    double edge = 0, pred = 0;             // edge at the slice centre: fitted, predicted
    double h = 0, eh = 0;                  // step height at the best position
    std::vector<double> xs, cs;            // scanned edge positions (slice centre) and chi2
};
EdgeFree FreeEdge(const Spectra& sp, int i) {
    EdgeFree r;
    const StepCfg c0 = SliceCfg(i, kSlLo[i], kSlHi[i]);
    const StepFit nom = scan::FitEdgeStep(sp, c0);
    if (!nom.ok) return r;
    const double ec = scan::M4max(0.5 * (kSlLo[i] + kSlHi[i]) / scan::kScalePiPi);
    r.pred = scan::kScale * ec; r.chi2nom = nom.chi2; r.ndf = nom.ndf;
    std::vector<double> sc, hh, eh;
    int ib = -1;
    for (double s = 0.9700; s < 1.01501; s += 0.00025) {
        StepCfg v = c0;
        v.scale = s; v.lo = nom.lo; v.hi = nom.hi; v.withEta = nom.withEta ? 1 : 0;
        const StepFit f = scan::FitEdgeStep(sp, v);
        if (!f.ok) continue;
        sc.push_back(s); r.cs.push_back(f.chi2); r.xs.push_back(s * ec); hh.push_back(f.h); eh.push_back(f.eh);
        if (ib < 0 || f.chi2 < r.cs[ib]) ib = (int)sc.size() - 1;
    }
    if (ib <= 0 || ib >= (int)sc.size() - 1) return r;      // minimum at the border: no free position
    r.scale = sc[ib]; r.chi2min = r.cs[ib]; r.h = hh[ib]; r.eh = eh[ib];
    int a = ib, b = ib;
    while (a > 0 && r.cs[a] < r.chi2min + 1.0) --a;
    while (b < (int)sc.size() - 1 && r.cs[b] < r.chi2min + 1.0) ++b;
    r.eLo = r.scale - sc[a]; r.eHi = sc[b] - r.scale;
    r.edge = r.scale * ec;
    r.ok = true;
    return r;
}

// mean of N - CB per fine bin between a and b, with its error
void Level(const Spectra& sp, double a, double b, double& mean, double& err) {
    double s = 0, v = 0; int n = 0;
    for (int i = 0; i < sp.nbin; ++i) {
        if (sp.x(i) < a || sp.x(i) >= b) continue;
        s += sp.sig[i]; v += sp.sigVar[i]; ++n;
    }
    mean = n ? s / n : 0; err = n ? std::sqrt(v) / n : 0;
}

// ---- 1D scan of one cut on top of a base selection ---------------------------------
struct Scan1D {
    int cut = 0;
    WP base;
    std::vector<double> thr;
    std::vector<Res> all, f0, f1;
};
Scan1D RunScan(const Data& dE, const Data& dG, const WP& base, int cut, bool folds) {
    Scan1D s; s.cut = cut; s.base = base; s.thr = Thresholds(cut);
    std::vector<Job> jobs;
    for (int f = -1; f <= (folds ? 1 : -1); ++f)
        for (double t : s.thr) {
            Job j; j.base = base; j.w = base; CutRef(j.w, cut) = t; j.fold = f; j.doT = (f < 0);
            jobs.push_back(j);
        }
    const std::vector<Res> r = Measure(dE, &dG, jobs);
    const size_t n = s.thr.size();
    s.all.assign(r.begin(), r.begin() + n);
    if (folds) { s.f0.assign(r.begin() + n, r.begin() + 2 * n); s.f1.assign(r.begin() + 2 * n, r.begin() + 3 * n); }
    return s;
}
void PrintScan(const Scan1D& s) {
    P("\n%s  (base: %s)\n", kCutName[s.cut], Describe(s.base).c_str());
    ScanHead();
    for (size_t i = 0; i < s.thr.size(); ++i) ScanLine(ThrName(s.cut, s.thr[i]).c_str(), s.all[i]);
}
// plateau of a figure of merit: thresholds within one statistical error of the best value
void Plateau(const Scan1D& s, bool useT) {
    int ib = 0;
    for (size_t i = 0; i < s.thr.size(); ++i) {
        const double g = useT ? s.all[i].gainT : s.all[i].gainE, gb = useT ? s.all[ib].gainT : s.all[ib].gainE;
        if (g > gb) ib = (int)i;
    }
    const double gb = useT ? s.all[ib].gainT : s.all[ib].gainE;
    const double eb = useT ? s.all[ib].eGainT : s.all[ib].eGainE;
    std::string in;
    for (size_t i = 0; i < s.thr.size(); ++i) {
        const double g = useT ? s.all[i].gainT : s.all[i].gainE;
        if (g >= gb - eb) { if (!in.empty()) in += ", "; in += ThrName(s.cut, s.thr[i]); }
    }
    P("   %-15s gain_%s: best %.3f +- %.3f at %s; within one error of the best: %s\n", kCutName[s.cut],
      useT ? "T" : "E", gb, eb, ThrName(s.cut, s.thr[ib]).c_str(), in.c_str());
}

// ---- 2D scan --------------------------------------------------------------------------
struct Scan2D {
    int cx = 0, cy = 0;
    WP base;
    std::vector<double> tx, ty;
    std::vector<Res> r;          // index ix * ny + iy
    const Res& at(int ix, int iy) const { return r[ix * ty.size() + iy]; }
};
Scan2D RunScan2D(const Data& dE, const Data& dG, const WP& base, int cx, int cy) {
    Scan2D s; s.cx = cx; s.cy = cy; s.base = base; s.tx = Thresholds(cx); s.ty = Thresholds(cy);
    std::vector<Job> jobs;
    for (double a : s.tx)
        for (double b : s.ty) {
            Job j; j.base = base; j.w = base; CutRef(j.w, cx) = a; CutRef(j.w, cy) = b;
            jobs.push_back(j);
        }
    s.r = Measure(dE, &dG, jobs);
    return s;
}
void PrintScan2D(const Scan2D& s) {
    const char* qn[8] = {"E (sum of step heights, rows / 2.5 MeV)", "statistical error of E", "eff_E = E / (E + E_fail)",
                         "eff_T (tagged eta)", "error of eff_T", "eff_N (rows kept in the bump window)",
                         "gain_E = eff_E / sqrt(eff_N)", "gain_T = eff_T / sqrt(eff_N)"};
    P("\n2D: %s (rows) against %s (columns); base: %s\n", kCutName[s.cx], kCutName[s.cy], Describe(s.base).c_str());
    for (int q = 0; q < 8; ++q) {
        P("  %s\n  %10s", qn[q], "");
        for (double b : s.ty) P(" %7s", ThrName(s.cy, b).c_str());
        P("\n");
        for (size_t ix = 0; ix < s.tx.size(); ++ix) {
            P("  %10s", ThrName(s.cx, s.tx[ix]).c_str());
            for (size_t iy = 0; iy < s.ty.size(); ++iy) {
                const Res& R = s.at(ix, iy);
                const double v[8] = {R.E, R.eE, R.effE, R.effT, R.eEffT, R.effN, R.gainE, R.gainT};
                P(q < 2 ? " %7.0f" : " %7.3f", v[q]);
            }
            P("\n");
        }
    }
}

// ---- model spread of the candidate numbers ------------------------------------------------
// spreadE   per-slice half-ranges of h (EdgeSpread) added linearly
// effE      eff_E = E / E(base) recomputed in 28 variations of the edge fit applied to the selection
//           and to the base alike: half range -10, 0, +10 MeV (first three slices), smearing
//           -1, 0, +1 MeV, mass scale of the edge -0.0009, 0, +0.0009 (order 2), and order 3 with
//           the half range 55 MeV; half of (max - min)
// effT      eff_T = T / T(base) in 18 variations of the peak fit (order 2, 3; ranges 0.49-0.62,
//           0.48-0.66, 0.50-0.60; width -0.3, 0, +0.3 MeV); half of (max - min)
struct CandVar { double spreadE = 0, effE = 0, effELo = 0, effEHi = 0, effT = 0, effTLo = 0, effTHi = 0; };
CandVar CandSpread(const Data& dE, const Data& dG, const WP& w, const WP& base) {
    CandVar cv;
    std::vector<Sel> sels;
    for (int v = 0; v < 2; ++v)
        for (int i = 0; i < kNSig; ++i) {
            const float lo = kSlLo[i], hi = kSlHi[i];
            const WP x = v ? base : w;
            sels.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && Pass(x, r); });
        }
    const std::vector<Spectra> sp = scan::FillMany(dE, sels);
    const std::vector<Spectra> pass(sp.begin(), sp.begin() + kNSig);
    cv.spreadE = EdgeSpread(pass, w).sum;
    bool first = true;
    for (int q = 0; q < 28; ++q) {
        const int ih = q < 27 ? q % 3 - 1 : 0, is = q < 27 ? (q / 3) % 3 - 1 : 0, ic = q < 27 ? q / 9 - 1 : 0;
        double e[2] = {0, 0}; bool good = true;
        for (int v = 0; v < 2 && good; ++v)
            for (int i = 0; i < kNSig; ++i) {
                double lo, hi; PassSlice(v ? base : w, i, lo, hi);
                if (hi - lo < 0.0045 || sp[v * kNSig + i].nsel[0] < 50) continue;
                StepCfg c = SliceCfg(i, lo, hi);
                if (q == 27) { c.order = 3; if (i < 3) c.halfRange = 0.055; }
                else { if (i < 3) c.halfRange += 0.010 * ih; c.s += 0.001 * is; c.scale += 0.0009 * ic; }
                const StepFit f = scan::FitEdgeStep(sp[v * kNSig + i], c);
                if (!f.ok) { good = false; break; }
                e[v] += f.h;
            }
        if (!good || !(e[1] > 0)) continue;
        const double x = e[0] / e[1];
        if (first) { cv.effELo = cv.effEHi = x; first = false; }
        cv.effELo = std::min(cv.effELo, x); cv.effEHi = std::max(cv.effEHi, x);
    }
    cv.effE = 0.5 * (cv.effEHi - cv.effELo);
    std::vector<Sel> ts;
    ts.push_back([=](const Row& r) { return Tagged(r) && Pass(w, r); });
    ts.push_back([=](const Row& r) { return Tagged(r) && Pass(base, r); });
    const std::vector<Spectra> st = scan::FillMany(dG, ts, TagOpt());
    const double rl[3] = {0.49, 0.48, 0.50}, rh[3] = {0.62, 0.66, 0.60};
    first = true;
    for (int ord = 2; ord <= 3; ++ord)
        for (int ir = 0; ir < 3; ++ir)
            for (int iw = -1; iw <= 1; ++iw) {
                PeakCfg c = TagCfg(); c.order = ord; c.lo = rl[ir]; c.hi = rh[ir]; c.sigma = kTagSigma + 0.0003 * iw;
                const PeakFit a = scan::FitPeakMain(st[0], c), b = scan::FitPeakMain(st[1], c);
                if (!a.ok || !b.ok || !(b.S > 0)) continue;
                const double x = a.S / b.S;
                if (first) { cv.effTLo = cv.effTHi = x; first = false; }
                cv.effTLo = std::min(cv.effTLo, x); cv.effTHi = std::max(cv.effTHi, x);
            }
    cv.effT = 0.5 * (cv.effTHi - cv.effTLo);
    return cv;
}

// ============================================================================
// D. Drawing
// ============================================================================
const char* const kXt = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
const char* const kYt = "rows / 2.5 MeV/c^{2}";
const char* const kHdr = "p+p 4.53 GeV, all combinations (isBest #geq 0), corrected momenta";

void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.09);
    gStyle->SetLabelSize(0.05, "xy");
    gStyle->SetTitleSize(0.052, "xy");
    gStyle->SetTitleOffset(1.55, "y");
    gStyle->SetTitleOffset(1.15, "x");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetPaintTextFormat(".2f");
}
void Deco(TH1* h, int col, int mk, double ms = 0.8) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->GetXaxis()->SetNdivisions(506);
}
void Txt(double x, double y, const char* t, double size = 0.045, int col = kBlack) {
    TLatex l; l.SetNDC(); l.SetTextSize(size); l.SetTextFont(42); l.SetTextColor(col);
    l.DrawLatex(x, y, t);
}
void VLine(double x, double y1, double y2, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x, y1, x, y2);
    l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid);
    l->Draw();
}
void HLine(double x1, double x2, double y, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x1, y, x2, y);
    l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid);
    l->Draw();
}
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs(Form("%s/%s.png", scan::kFigDir, name));
    if (pdf) c->SaveAs(Form("%s/%s.pdf", scan::kFigDir, name));
}
TGraphErrors* Gr(int col, int mk, double ms = 1.3) {
    TGraphErrors* g = new TGraphErrors();
    g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineWidth(2);
    return g;
}
void Add(TGraphErrors* g, double x, double y, double e) {
    const int n = g->GetN();
    g->SetPoint(n, x, y); g->SetPointError(n, 0, e);
}
// frame with one labelled bin per threshold
TH1F* Frame(const char* name, int cut, const std::vector<double>& thr, double y1, double y2, const char* yt) {
    const int n = (int)thr.size();
    TH1F* h = new TH1F(name, "", n, 0, n);
    h->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) h->GetXaxis()->SetBinLabel(i + 1, ThrName(cut, thr[i]).c_str());
    h->GetXaxis()->SetLabelSize(n > 9 ? 0.052 : 0.065); h->GetXaxis()->SetTitle(kCutAxis[cut]); h->GetXaxis()->SetTitleSize(0.046);
    h->GetXaxis()->SetTitleOffset(1.35); h->GetXaxis()->LabelsOption("h");
    h->GetYaxis()->SetTitle(yt);
    h->SetMinimum(y1); h->SetMaximum(y2);
    return h;
}
// one pad with N - CB of a slice and its edge fit; kind 0 signal slice, 1 not used in E, 2 control
void DrawSlice(const Spectra& sp, const StepFit& f, const char* name, const char* title, int kind) {
    TH1D* h = sp.Hist("sig", name);
    Deco(h, kBlue + 1, 20, 0.9);
    h->GetXaxis()->SetTitle(kXt); h->GetYaxis()->SetTitle(kYt);
    const double x1 = f.lo - 0.015, x2 = f.hi + 0.015;
    h->GetXaxis()->SetRangeUser(x1 + 1e-4, x2 - 1e-4);
    double mx = 0, mn = 1e30;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double x = h->GetBinCenter(i);
        if (x < x1 || x > x2) continue;
        mx = std::max(mx, h->GetBinContent(i) + h->GetBinError(i));
        mn = std::min(mn, h->GetBinContent(i) - h->GetBinError(i));
    }
    const double lo = mn > 0 ? std::max(0.0, mn - 0.35 * (mx - mn)) : mn - 0.05 * (mx - mn);
    h->SetMinimum(lo); h->SetMaximum(mx + 0.60 * (mx - lo));
    h->Draw("e");
    if (f.ok) {
        TGraph* g = new TGraph(); TGraph* g0 = new TGraph();
        for (size_t i = 0; i < f.tx.size(); ++i) { g->SetPoint((int)i, f.tx[i], f.fy[i]); g0->SetPoint((int)i, f.tx[i], f.fy0[i]); }
        g->SetLineWidth(2); g0->SetLineWidth(2); g0->SetLineStyle(2);
        g0->Draw("l same"); g->Draw("l same");
        VLine(f.edge, lo, lo + 0.58 * (h->GetMaximum() - lo), kGreen + 2, 7);
        Txt(0.21, 0.84, Form("step h = %.0f #pm %.0f rows / 2.5 MeV", f.h, f.eh), 0.044);
        Txt(0.21, 0.785, Form("#chi^{2}/ndf = %.1f/%d,   #Delta#chi^{2}(step) = %.1f", f.chi2, f.ndf, f.dchi2), 0.044);
        Txt(0.21, 0.73, kind == 2 ? Form("green: edge set by hand at %.4f (none expected)", f.edge)
                                  : Form("green: predicted edge %.4f", f.edge), 0.040, kGreen + 2);
        Txt(0.21, 0.68, "full line: fit; dashed: the fit without its step", 0.036);
        if (kind == 1) Txt(0.58, 0.20, "not used in E", 0.044, kRed + 1);
    }
    Txt(0.17, 0.935, title, 0.038);
}

} // namespace bump

// ============================================================================
void scan_bump() {
    using namespace bump;
    TStopwatch swAll;
    gSystem->mkdir(scan::kFigDir, kTRUE);

    // ------------------------------------------------------------------------
    // data: dE = rows with 0.30 <= M4 < 0.62 (edge fits, bump window),
    //       dG = rows with a photon and m_ee < 0.1 (pi0 and tagged eta)
    // ------------------------------------------------------------------------
    Data dE, dG;
    {
        Data d;
        if (!scan::Load(d)) return;
        Filter(d, dE, [](const Row& r) { return r.m4_fb >= 0.30f && r.m4_fb < 0.62f; });
        Filter(d, dG, [](const Row& r) { return HasGam(r) && r.m_eeg >= 0.05f && r.m_eeg < 0.25f; });
    }
    printf("rows: M4 sample %zu / %zu / %zu, photon sample %zu / %zu / %zu\n", dE.s[0].size(), dE.s[1].size(),
           dE.s[2].size(), dG.s[0].size(), dG.s[1].size(), dG.s[2].size());

    // ========================================================================
    // 1. pi0 window and the tagged eta (bump_pi0.txt)
    // ========================================================================
    Open("bump_pi0.txt");
    P("==========================================================================================\n");
    P("bump_pi0.txt  --  written by plots/scan_bump.C\n");
    P("ECAL-tagged eta -> pi+pi-pi0 with an e+e- pair from the pi0 (Dalitz decay pi0 -> e+e-gamma, or pi0 -> gamma gamma with one photon\n");
    P("converted: the same M(e+e-gamma), not told apart here). Rows = stored combinations (isBest >= 0), cor_ok == 1,\n");
    P("photon = the good photon (of the three most energetic) with M(e+e-gamma) closest to 0.140, m_ee < %.2f.\n", kMeeMax);
    P("==========================================================================================\n\n");

    // ---- (A) the eta peak moves with M(e+e-gamma) in the plain m5, not in m5c ----
    {
        const int nb = 5;
        const double bl[nb] = {0.115, 0.130, 0.140, 0.150, 0.160}, bh[nb] = {0.130, 0.140, 0.150, 0.160, 0.175};
        std::vector<Sel> sb;
        for (int i = 0; i < nb; ++i) {
            const float a = bl[i], b = bh[i];
            sb.push_back([=](const Row& r) { return HasGam(r) && r.m_eeg >= a && r.m_eeg < b; });
        }
        scan::FillOpt o5 = TagOpt();
        o5.x = [](const Row& r) { return (double)r.m5; };
        const std::vector<Spectra> s5 = scan::FillMany(dG, sb, o5);
        const std::vector<Spectra> s5c = scan::FillMany(dG, sb, TagOpt());
        P("A. Position of the eta peak (free fit, pol2) in bands of M(e+e-gamma): plain m5 against the constrained m5c.\n");
        P("   m5c = M(pi+pi-e+e-gamma) with the photon energy rescaled so that M(e+e-gamma) = m_pi0.\n");
        P("   %-13s | %-42s | %-42s\n", "M(eeg) band", "plain m5: position, width [MeV], S, chi2/ndf", "m5c: position, width [MeV], S, chi2/ndf");
        for (int i = 0; i < nb; ++i) {
            PeakCfg c5 = TagCfg();
            c5.mean = 0.5475 + 0.96 * (0.5 * (bl[i] + bh[i]) - 0.135); c5.sigma = 0.008; c5.lo = c5.mean - 0.065; c5.hi = c5.mean + 0.085;
            const scan::FreeFit f5 = scan::FitPeakFree(s5[i], c5, false);
            const scan::FreeFit fc = scan::FitPeakFree(s5c[i], TagCfg(), false);
            P("   %.3f-%.3f   | %d %.4f+-%.4f %5.1f+-%-4.1f %5.0f+-%-4.0f %5.1f/%-2d | %d %.4f+-%.4f %5.1f+-%-4.1f %5.0f+-%-4.0f %5.1f/%-2d\n",
              bl[i], bh[i], (int)f5.ok, f5.mean, f5.eMean, 1e3 * f5.sigma, 1e3 * f5.eSigma, f5.S, f5.eS, f5.chi2, f5.ndf,
              (int)fc.ok, fc.mean, fc.eMean, 1e3 * fc.sigma, 1e3 * fc.eSigma, fc.S, fc.eS, fc.chi2, fc.ndf);
        }
    }

    // ---- (B) eta yield against M(e+e-gamma): the pi0 of tagged eta events ----
    const int nPi = 32; const double piLo = 0.060, piBw = 0.005;
    std::vector<double> piS[2], piE[2];       // [0] all photons, [1] rows with exactly one photon
    double pi0Par[2][3], pi0Err[2][3], pi0Chi2[2]; int pi0Ndf[2];
    Spectra spMeeg;                           // M(e+e-gamma) of all rows with a photon
    {
        for (int v = 0; v < 2; ++v) {
            std::vector<Sel> sb;
            for (int i = 0; i < nPi; ++i) {
                const float a = piLo + i * piBw, b = a + piBw;
                sb.push_back([=](const Row& r) { return HasGam(r) && (v == 0 || r.n_gam == 1) && r.m_eeg >= a && r.m_eeg < b; });
            }
            const std::vector<Spectra> s = scan::FillMany(dG, sb, TagOpt());
            TGraphErrors g;
            for (int i = 0; i < nPi; ++i) {
                const PeakFit f = scan::FitPeakMain(s[i], TagCfg());
                piS[v].push_back(f.S); piE[v].push_back(f.eS);
                g.SetPoint(i, piLo + (i + 0.5) * piBw, f.S); g.SetPointError(i, 0, f.eS);
            }
            TF1 fg("fg", "[0]*0.005/(sqrt(2*pi)*[2])*exp(-0.5*((x-[1])/[2])^2)", 0.095, 0.200);
            fg.SetParameters(3500, 0.146, 0.012);
            TFitResultPtr fr = g.Fit(&fg, "RQSN");
            for (int k = 0; k < 3; ++k) { pi0Par[v][k] = fr->Parameter(k); pi0Err[v][k] = fr->ParError(k); }
            pi0Chi2[v] = fr->Chi2(); pi0Ndf[v] = fr->Ndf();
        }
        P("\nB. Eta yield (m5c fit, frozen shape %.4f / %.1f MeV, pol%d, %.2f-%.2f) in 5 MeV bins of M(e+e-gamma).\n",
          kTagMean, 1e3 * kTagSigma, kTagOrder, kTagLo, kTagHi);
        P("   %-13s | %-16s | %-16s\n", "M(eeg) bin", "all photons", "exactly one photon");
        for (int i = 0; i < nPi; ++i)
            P("   %.3f-%.3f   | %7.1f +- %5.1f | %7.1f +- %5.1f\n", piLo + i * piBw, piLo + (i + 1) * piBw, piS[0][i], piE[0][i], piS[1][i], piE[1][i]);
        for (int v = 0; v < 2; ++v)
            P("   Gaussian fit 0.095-0.200 (%s): position %.4f +- %.4f, width %.2f +- %.2f MeV, content %.0f +- %.0f, chi2/ndf %.1f/%d\n",
              v ? "one photon" : "all photons", pi0Par[v][1], pi0Err[v][1], 1e3 * pi0Par[v][2], 1e3 * pi0Err[v][2], pi0Par[v][0], pi0Err[v][0], pi0Chi2[v], pi0Ndf[v]);
        P("   frozen in this macro: position %.4f, width %.1f MeV, window +-2 sigma = %.3f-%.3f\n", kPi0Pos, 1e3 * kPi0Sig, kPi0Lo, kPi0Hi);

        // M(e+e-gamma) spectrum of all rows with a photon
        scan::FillOpt og; og.nbin = 100; og.lo = 0.05; og.hi = 0.25; og.reb = 1;
        og.x = [](const Row& r) { return (double)r.m_eeg; };
        spMeeg = scan::Fill(dG, [](const Row& r) { return HasGam(r); }, og);
    }

    // ---- (C) windows: content, share of the old window, S/B ----
    struct TagRef { double lo, hi; PeakFit m, p; scan::FreeFit ff; scan::Spread sprd; double rows[3]; };
    std::vector<TagRef> tagRef;
    Spectra spTag, spTagBest, spOld, spOldBest;
    {
        const int nw = 8;
        const double wl[nw] = {kPi0Lo, 0.125, 0.130, 0.122, 0.115, 0.110, 0.100, 0.080};
        const double wh[nw] = {kPi0Hi, 0.145, 0.160, 0.158, 0.165, 0.180, 0.190, 0.220};
        std::vector<Sel> sb;
        for (int i = 0; i < nw; ++i) {
            const float a = wl[i], b = wh[i];
            sb.push_back([=](const Row& r) { return HasGam(r) && r.m_eeg >= a && r.m_eeg < b; });
        }
        // complements of the old window inside the widest one and inside the frozen one
        sb.push_back([](const Row& r) { return HasGam(r) && r.m_eeg >= 0.080f && r.m_eeg < 0.220f && !(r.m_eeg >= 0.125f && r.m_eeg < 0.145f); });
        sb.push_back([](const Row& r) { return HasGam(r) && r.m_eeg >= 0.080f && r.m_eeg < 0.220f && !Tagged(r); });
        sb.push_back([](const Row& r) { return Tagged(r) && r.isBest == 1; });
        const std::vector<Spectra> s = scan::FillMany(dG, sb, TagOpt());
        P("\nC. Tagged eta in m5c for several M(e+e-gamma) windows. Frozen shape as in B; S/B = S in +-1.9 sigma / (polynomial + CB).\n");
        P("   free = position and width free. spread = half of (max - min) of S over orders 2, 3, ranges 0.49-0.62, 0.48-0.66, 0.50-0.60,\n");
        P("   width -+0.3 MeV (18 fits, those with p > 0.05).\n");
        P("   %-13s | %8s %6s %7s | %6s %5s %9s %6s %6s | %6s %5s | %-34s | %6s\n", "window", "rows +-", "++", "--", "S", "err",
          "chi2/ndf", "S/B", "S/err", "S_pois", "err", "free: position, width [MeV], S", "spread");
        for (int i = 0; i < nw; ++i) {
            TagRef t; t.lo = wl[i]; t.hi = wh[i];
            for (int k = 0; k < 3; ++k) t.rows[k] = s[i].nsel[k];
            t.m = scan::FitPeakMain(s[i], TagCfg()); t.p = scan::FitPeakPois(s[i], TagCfg());
            t.ff = scan::FitPeakFree(s[i], TagCfg(), false);
            // model spread
            double mn = t.m.S, mx = t.m.S;
            const double rl[3] = {0.49, 0.48, 0.50}, rh[3] = {0.62, 0.66, 0.60};
            for (int ord = 2; ord <= 3; ++ord) for (int ir = 0; ir < 3; ++ir) for (int iw = -1; iw <= 1; ++iw) {
                PeakCfg c = TagCfg(); c.order = ord; c.lo = rl[ir]; c.hi = rh[ir]; c.sigma = kTagSigma + 0.0003 * iw;
                const PeakFit f = scan::FitPeakMain(s[i], c);
                if (!f.ok || !(f.prob > 0.05)) continue;
                mn = std::min(mn, f.S); mx = std::max(mx, f.S);
            }
            t.sprd.half = 0.5 * (mx - mn);
            tagRef.push_back(t);
            P("   %.3f-%.3f%s | %8.0f %6.0f %7.0f | %6.0f %5.0f %6.1f/%-2d %6.3f %6.1f | %6.0f %5.0f | %d %.4f+-%.4f %4.1f+-%-3.1f %5.0f+-%-4.0f | %6.0f\n",
              wl[i], wh[i], i == 0 ? "*" : " ", t.rows[0], t.rows[1], t.rows[2], t.m.S, t.m.eS, t.m.chi2, t.m.ndf, t.m.SoB, t.m.F,
              t.p.S, t.p.eS, (int)t.ff.ok, t.ff.mean, t.ff.eMean, 1e3 * t.ff.sigma, 1e3 * t.ff.eSigma, t.ff.S, t.ff.eS, t.sprd.half);
        }
        P("   (* = frozen window of the tagged sample T)\n");
        P("   tagged sample T: %.0f unlike-sign rows (%.0f e+e+, %.0f e-e-); in the peak window %.4f-%.4f: N = %.0f, CB = %.0f, S inside = %.0f\n",
          tagRef[0].rows[0], tagRef[0].rows[1], tagRef[0].rows[2], TagCfg().winLo, TagCfg().winHi, tagRef[0].m.N, tagRef[0].m.CB, tagRef[0].m.Swin);
        const PeakFit in = tagRef[1].m, out = scan::FitPeakMain(s[nw], TagCfg());
        double f, ef; Ratio(in.S, in.eS, out.S, out.eS, f, ef);
        P("   share of the old window 0.125-0.145 in 0.080-0.220: S_in %.0f +- %.0f, S_out %.0f +- %.0f, share %.3f +- %.3f\n", in.S, in.eS, out.S, out.eS, f, ef);
        const PeakFit in2 = tagRef[0].m, out2 = scan::FitPeakMain(s[nw + 1], TagCfg());
        Ratio(in2.S, in2.eS, out2.S, out2.eS, f, ef);
        P("   share of the frozen window %.3f-%.3f in 0.080-0.220: S_in %.0f +- %.0f, S_out %.0f +- %.0f, share %.3f +- %.3f\n", kPi0Lo, kPi0Hi, in2.S, in2.eS, out2.S, out2.eS, f, ef);
        spTag = s[0]; spTagBest = s[nw + 2];
        const PeakFit tb = scan::FitPeakMain(spTagBest, TagCfg());
        P("   frozen window, isBest == 1 only: rows %.0f, S = %.0f +- %.0f, chi2/ndf %.1f/%d, S/B %.3f\n", spTagBest.nsel[0], tb.S, tb.eS, tb.chi2, tb.ndf, tb.SoB);

        // plain m5 in the old window (comparison with the old production and the engine's spot check)
        scan::FillOpt o5 = TagOpt(); o5.nbin = 80; o5.reb = 4; o5.windows.clear();
        o5.x = [](const Row& r) { return (double)r.m5; };
        PeakCfg c5 = TagCfg(); c5.mean = kOldMean; c5.sigma = kOldSigma; c5.lo = 0.48; c5.hi = 0.64;
        c5.winLo = c5.mean - 1.9 * c5.sigma; c5.winHi = c5.mean + 1.9 * c5.sigma;
        std::vector<Sel> so;
        so.push_back([](const Row& r) { return HasGam(r) && r.m_eeg >= (float)kOldLo && r.m_eeg < (float)kOldHi; });
        so.push_back([](const Row& r) { return HasGam(r) && r.m_eeg >= (float)kOldLo && r.m_eeg < (float)kOldHi && r.isBest == 1; });
        so.push_back([](const Row& r) { return Tagged(r); });
        const std::vector<Spectra> s5 = scan::FillMany(dG, so, o5);
        spOld = s5[0]; spOldBest = s5[1];
        P("\nD. Plain m5 = M(pi+pi-e+e-gamma), 5 MeV bins, pol2, 0.48-0.64 (the old way).\n");
        for (int i = 0; i < 2; ++i) {
            const scan::FreeFit ff = scan::FitPeakFree(s5[i], c5, false);
            const PeakFit m = scan::FitPeakMain(s5[i], c5), p = scan::FitPeakPois(s5[i], c5);
            P("   window 0.125-0.145, %s: free position %.4f +- %.4f, width %.1f +- %.1f MeV, S %.0f +- %.0f; frozen %.4f / %.1f MeV: S = %.0f +- %.0f (chi2/ndf %.1f/%d, S/B %.3f), Poisson fit %.0f +- %.0f\n",
              i ? "isBest == 1" : "isBest >= 0", ff.mean, ff.eMean, 1e3 * ff.sigma, 1e3 * ff.eSigma, ff.S, ff.eS, kOldMean, 1e3 * kOldSigma, m.S, m.eS, m.chi2, m.ndf, m.SoB, p.S, p.eS);
        }
        {
            PeakCfg cw = c5; cw.mean = 0.556; cw.sigma = 0.012; cw.lo = 0.48; cw.hi = 0.66;
            const scan::FreeFit ff = scan::FitPeakFree(s5[2], cw, false);
            if (ff.ok) P("   frozen window %.3f-%.3f in the plain m5: free position %.4f, width %.1f +- %.1f MeV, S %.0f +- %.0f, chi2/ndf %.1f/%d\n",
                         kPi0Lo, kPi0Hi, ff.mean, 1e3 * ff.sigma, 1e3 * ff.eSigma, ff.S, ff.eS, ff.chi2, ff.ndf);
            else P("   frozen window %.3f-%.3f in the plain m5: the free fit does not converge (the peak is wider than the engine's limit of 12 MeV).\n", kPi0Lo, kPi0Hi);
        }
    }
    Close();

    // ========================================================================
    // 2. edge step: reference fits, variations, edge position, lower edge (bump_reference.txt)
    // ========================================================================
    Open("bump_reference.txt");
    P("==========================================================================================\n");
    P("bump_reference.txt  --  written by plots/scan_bump.C\n");
    P("Edge-step fit E_step of the engine in M(pi+pi-) slices. Rows = stored combinations (isBest >= 0), cor_ok == 1.\n");
    P("Model: N - CB = polynomial (order 2) [+ eta Gaussian] + h * smeared step at the FIXED predicted edge.\n");
    P("h = step height in rows per 2.5 MeV. E = sum of h over the first four slices.\n");
    P("==========================================================================================\n");
    const int nRef = 2;
    WP refWP[nRef] = {None(), None()};
    refWP[1].name = "oa_lab<50"; refWP[1].oaLab = 50;
    std::vector<Spectra> refSl[nRef];
    std::vector<StepFit> refFit[nRef];
    std::vector<EdgeFree> refFree;
    for (int b = 0; b < nRef; ++b) {
        std::vector<Sel> sels;
        const WP w = refWP[b];
        for (int i = 0; i < kNSl; ++i) {
            const float lo = kSlLo[i], hi = kSlHi[i];
            sels.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && Pass(w, r); });
        }
        refSl[b] = scan::FillMany(dE, sels);
        P("\nA%d. Selection: %s\n", b + 1, Describe(w).c_str());
        P("   %-10s %7s %15s %15s %3s | %6s %5s %6s %6s | %9s %5s %7s %6s | %6s %5s | %6s %4s\n", "M(pipi)", "edge", "edge range",
          "fit range", "eta", "h", "err", "bkg", "h/bkg", "chi2/ndf", "p", "no step", "dchi2", "S_eta", "err", "spread", "kept");
        const ESpread es = EdgeSpread(refSl[b], w);
        double E = 0, v = 0;
        for (int i = 0; i < kNSl; ++i) {
            const StepFit f = scan::FitEdgeStep(refSl[b][i], SliceCfg(i, kSlLo[i], kSlHi[i]));
            refFit[b].push_back(f);
            P("   %-10s %7.4f %7.4f-%7.4f %7.4f-%7.4f %3d | %6.1f %5.1f %6.0f %6.3f | %6.1f/%-2d %5.3f %7.1f %6.1f | %6.0f %5.0f |",
              kSlName[i], f.edge, f.edgeLo, f.edgeHi, f.lo, f.hi, (int)f.withEta, f.h, f.eh, f.bkgAtEdge,
              f.bkgAtEdge > 0 ? f.h / f.bkgAtEdge : 0., f.chi2, f.ndf, TMath::Prob(f.chi2, f.ndf), f.chi2NoStep, f.dchi2, f.Seta, f.eSeta);
            if (i < kNSig) { P(" %6.1f %4d\n", es.half[i], es.kept[i]); E += f.h; v += f.eh * f.eh; }
            else P(" %s\n", i == 4 ? "(not in E)" : "(control, edge set by hand)");
        }
        P("   E = %.0f +- %.0f (stat) +- %.0f (model spread, slices added linearly)\n", E, std::sqrt(v), es.sum);
    }
    // ---- free edge position ----
    P("\nB. Edge position fitted free (selection: no cut). The factor between the true limit M4max(M(pipi)/%.4f) and the measured\n", scan::kScalePiPi);
    P("   edge is scanned in steps of 0.00025; fit range and model of the nominal fit. Nominal factor %.4f +- 0.0009 (eta peak).\n", scan::kScale);
    P("   %-10s | %7s %7s %7s | %9s %9s %7s | %8s %8s %7s | %6s %5s\n", "M(pipi)", "factor", "-err", "+err", "edge fit", "predicted",
      "diff", "chi2 min", "chi2 nom", "ndf", "h", "err");
    for (int i = 0; i < kNSig; ++i) {
        const EdgeFree f = FreeEdge(refSl[0][i], i);
        refFree.push_back(f);
        if (!f.ok) { P("   %-10s | no minimum inside the scanned interval 0.970-1.015\n", kSlName[i]); continue; }
        P("   %-10s | %7.4f %7.4f %7.4f | %9.4f %9.4f %+7.4f | %8.1f %8.1f %7d | %6.1f %5.1f\n", kSlName[i], f.scale, f.eLo, f.eHi,
          f.edge, f.pred, f.edge - f.pred, f.chi2min, f.chi2nom, f.ndf, f.h, f.eh);
    }
    {
        double sw = 0, sx = 0;
        for (int i = 0; i < kNSig; ++i) {
            if (!refFree[i].ok) continue;
            const double e = 0.5 * (refFree[i].eLo + refFree[i].eHi), w = e > 0 ? 1.0 / (e * e) : 0;
            sw += w; sx += w * refFree[i].scale;
        }
        if (sw > 0) P("   weighted mean factor %.4f +- %.4f (eta peak: %.4f +- 0.0009)\n", sx / sw, 1.0 / std::sqrt(sw), scan::kScale);
    }
    // ---- lower edge ----
    std::vector<Spectra> lowSp, upSp;
    {
        std::vector<Sel> sels;
        for (int i = 0; i < kNSig; ++i) {
            const float lo = kSlLo[i], hi = kSlHi[i];
            sels.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi; });
        }
        scan::FillOpt ol; ol.nbin = 120; ol.lo = 0.05; ol.hi = 0.35; ol.reb = 4;
        ol.x = [](const Row& r) { const double e = LowEdgeMeasured(r.mpipi); return e > 0 ? r.m4 - e + 0.2 : -1.0; };
        scan::FillOpt ou = ol;
        ou.x = [](const Row& r) { const double e = scan::EdgeMeasured(r.mpipi); return e > 0 ? r.m4 - e + 0.2 : -1.0; };
        lowSp = scan::FillMany(dE, sels, ol);
        upSp = scan::FillMany(dE, sels, ou);
        P("\nC. Lower edge M4min^2 = m_eta^2 - m_eta (E0 + p0) (selection: no cut). Rows are aligned row by row at the predicted\n");
        P("   lower (upper) edge: x = M4 - edge(M(pipi)). level = mean N - CB per 2.5 MeV. A box-shaped bump would give the same\n");
        P("   level just above the lower edge as just below the upper edge.\n");
        P("   %-10s | %13s | %17s %17s %17s | %15s | %-28s\n", "M(pipi)", "M4min (centre)", "level -10..+10 MeV", "level +10..+30", "below upper edge",
          "ratio at M4min", "step fit at M4min: rise, chi2/ndf, dchi2");
        for (int i = 0; i < kNSig; ++i) {
            double l0, e0, l1, e1, lu, eu;
            Level(lowSp[i], 0.19, 0.21, l0, e0); Level(lowSp[i], 0.21, 0.23, l1, e1); Level(upSp[i], 0.175, 0.195, lu, eu);
            StepCfg c; c.useEdge = true; c.edge = 0.2; c.halfRange = 0.04; c.withEta = 0;
            const StepFit f = scan::FitEdgeStep(lowSp[i], c);
            P("   %-10s | %13.4f | %8.1f +- %-5.1f %8.1f +- %-5.1f %8.1f +- %-5.1f | %6.3f +- %5.3f | %7.1f +- %-5.1f %6.1f/%-2d %6.1f\n", kSlName[i],
              LowEdgeMeasured(0.5 * (kSlLo[i] + kSlHi[i])), l0, e0, l1, e1, lu, eu, lu > 0 ? l0 / lu : 0., lu > 0 ? e0 / lu : 0.,
              -f.h, f.eh, f.chi2, f.ndf, f.dchi2);
        }
    }
    Close();

    // ========================================================================
    // 3. one-dimensional scans from "no cut" (bump_scan1d.txt)
    // ========================================================================
    Open("bump_scan1d.txt");
    P("==========================================================================================\n");
    P("bump_scan1d.txt  --  written by plots/scan_bump.C\n");
    P("One cut at a time on top of the base selection. Rows = stored combinations (isBest >= 0), cor_ok == 1.\n");
    P("  E       sum of the edge-step heights (rows / 2.5 MeV) of the slices 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40; statistical error\n");
    P("  eff_E   E / E(base); its error takes into account that the rows are common. E_fail = the same fit on the rows of the base that FAIL\n");
    P("          the cut; close = (E + E_fail) / E(base), 1 when the two parts add up.\n");
    P("  T_pass, T_fail  eta peak content in m5c of the tagged rows (pi0 window %.3f-%.3f) passing / failing the cut, with the chi2 of the two\n", kPi0Lo, kPi0Hi);
    P("          fits (ndf 48); eff_T = T_pass / T(base); close = (T_pass + T_fail) / T(base).\n");
    P("  N, CB, W = N - CB in the bump window %.4f-%.4f; eff_W = W / W(base); eff_N = N / N(base)\n", scan::kBumpLo, scan::kBumpHi);
    P("  Q = E / sqrt(N); gain_E = eff_E / sqrt(eff_N) = Q / Q(base); gain_T = eff_T / sqrt(eff_N)   (noise = sqrt of ALL rows; no CB estimate enters)\n");
    P("  CB = 2 sqrt(N++ N--) without a k factor, N--/N++ about 6: CB, CB/N and W are estimates (W uncertain by several per cent).\n");
    P("  E is a sum of step heights, not a yield.\n");
    P("  A cut helps the bump when gain > 1 by more than its error.\n");
    P("  Not valid: E with the eta-rest cut at 120 and 140 (chi2 of the edge fits); T_fail for upper M(pipi) cuts (the failing rows start at\n");
    P("  M(pipi) + m_pi0, inside the fit range; this is why eff_T uses T(base)). M(pipi) cuts outside 0.30-0.40 cannot change E by construction.\n");
    P("==========================================================================================\n");
    std::vector<Scan1D> sc1;
    for (int c = 0; c < kNCut; ++c) {
        sc1.push_back(RunScan(dE, dG, None(), c, true));
        PrintScan(sc1.back());
    }
    P("\nPlateaus (thresholds whose gain lies within one statistical error of the best gain):\n");
    for (int c = 0; c < kNCut; ++c) { Plateau(sc1[c], false); Plateau(sc1[c], true); }
    P("\nTwo folds (event groups of even / odd index), E and eff_E only:\n");
    for (int c = 0; c < kNCut; ++c) {
        const Scan1D& s = sc1[c];
        P("\n%s\n   %-8s | %6s %5s %6s %5s | %6s %5s %6s %5s | %7s %5s | %6s %5s %6s %5s\n", kCutName[c], "thresh", "E(A)", "err", "eff_E", "err",
          "E(B)", "err", "eff_E", "err", "A - B", "err", "gain_E A", "err", "gain_E B", "err");
        for (size_t i = 0; i < s.thr.size(); ++i) {
            const Res& a = s.f0[i]; const Res& b = s.f1[i];
            P("   %-8s | %6.0f %5.0f %6.3f %5.3f | %6.0f %5.0f %6.3f %5.3f | %7.0f %5.0f | %6.3f %5.3f %6.3f %5.3f\n", ThrName(c, s.thr[i]).c_str(),
              a.E, a.eE, a.effE, a.eEffE, b.E, b.eE, b.effE, b.eEffE, a.E - b.E, std::sqrt(a.eE * a.eE + b.eE * b.eE),
              a.gainE, a.eGainE, b.gainE, b.eGainE);
        }
    }
    Close();

    // ========================================================================
    // 4. two-dimensional checks from "no cut" (bump_scan2d.txt)
    // ========================================================================
    Open("bump_scan2d.txt");
    P("==========================================================================================\n");
    P("bump_scan2d.txt  --  written by plots/scan_bump.C   (quantities as in bump_scan1d.txt)\n");
    P("==========================================================================================\n");
    const Scan2D s2a = RunScan2D(dE, dG, None(), kOaLab, kMpHi);
    PrintScan2D(s2a);
    const Scan2D s2b = RunScan2D(dE, dG, None(), kOaLo, kOaHi);
    PrintScan2D(s2b);
    Close();

    // ========================================================================
    // 5. candidates of the bump selection (bump_candidates.txt)
    // ========================================================================
    Open("bump_candidates.txt");
    P("==========================================================================================\n");
    P("bump_candidates.txt  --  written by plots/scan_bump.C   (quantities as in bump_scan1d.txt; base of the efficiencies: no cut)\n");
    P("==========================================================================================\n");
    const int nCand = 7;
    const WP cand[nCand] = {None(), Loose(), Nominal(), Pure(), Slides(), OldCode(), OldCodeOA()};
    std::vector<Res> cr;
    std::vector<CandVar> cvar;
    {
        std::vector<Job> cj;
        for (int f = -1; f <= 1; ++f)
            for (int c = 0; c < nCand; ++c) { Job j; j.w = cand[c]; j.base = None(); j.fold = f; cj.push_back(j); }
        cr = Measure(dE, &dG, cj);
        for (int c = 0; c < nCand; ++c) cvar.push_back(CandSpread(dE, dG, cand[c], None()));
    }
    P("\nA. Selections\n");
    for (int c = 0; c < nCand; ++c) P("   %-18s %s\n", cand[c].name.c_str(), Describe(cand[c]).c_str());
    P("\nB. Edge step. E: stat = fit error; spread = model variations (order 2, 3; half range 35, 45, 55 MeV; smearing 4.4, 5.4, 6.4 MeV;\n");
    P("   fits with p > 0.05), per-slice half-ranges added linearly. eff_E = E / E(no cut): err = statistical, spread = half of (max - min)\n");
    P("   over 28 variations applied to both selections (half range -+10 MeV, smearing -+1 MeV, edge scale -+0.0009; order 3 with 55 MeV).\n");
    P("   %-18s | %6s %5s %6s %8s %6s | %6s %5s %6s %5s | %6s %5s | %s\n", "selection", "E", "stat", "spread", "chi2/ndf", "h/bkg",
      "eff_E", "err", "spread", "close", "E_fail", "err", "h per slice (0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 | 0.40-0.41 | control)");
    for (int c = 0; c < nCand; ++c) {
        const Res& R = cr[c];
        P("   %-18s | %6.0f %5.0f %6.0f %4.0f/%-3d %6.3f | %6.3f %5.3f %6.3f %5.3f | %6.0f %5.0f |", cand[c].name.c_str(), R.E, R.eE,
          cvar[c].spreadE, R.chi2, R.ndf, R.bkg > 0 ? R.E / R.bkg : 0., R.effE, R.eEffE, cvar[c].effE, R.closeE, R.Ef, R.eEf);
        for (int i = 0; i < kNSl; ++i) {
            if (i == kNSig || i == kNSig + 1) P(" |");
            if (R.okP[i]) P(" %4.0f+-%-3.0f", R.fp[i].h, R.fp[i].eh); else P(" %8s", "-");
        }
        P("\n");
    }
    P("   E is NOT valid for the three selections with the eta-rest cut (chi2/ndf of the edge fits).\n");
    P("\nC. Bump window %.4f-%.4f and the figure of merit Q = E / sqrt(N). gain_E = eff_E / sqrt(eff_N) = Q / Q(no cut).\n", scan::kBumpLo, scan::kBumpHi);
    P("   CB = 2 sqrt(N++ N--) (no k factor). N++, N-- = like-sign rows; LSsum = N++ + N-- is the other extreme of a like-sign background estimate;\n");
    P("   dW = change of W = N - CB if LSsum were used: the size of the uncertainty of W that the statistical error does not show.\n");
    P("   %-18s | %7s %6s %7s %5s %6s %6s %6s | %5s %5s | %6s %5s %6s | %6s %6s %7s %6s %6s\n", "selection", "N", "CB", "W=N-CB", "err", "CB/N", "eff_W", "eff_N", "Q", "err", "gain_E", "err", "spread",
      "N++", "N--", "N--/N++", "LSsum", "dW [%]");
    for (int c = 0; c < nCand; ++c) {
        const Res& R = cr[c];
        P("   %-18s | %7.0f %6.0f %7.0f %5.0f %6.3f %6.3f %6.3f | %5.2f %5.2f | %6.3f %5.3f %6.3f | %6.0f %6.0f %7.1f %6.0f %6.1f\n", cand[c].name.c_str(), R.win.n[0], R.win.cb,
          R.win.sig(), R.win.sigErr(), R.win.cbShare(), R.effW, R.effN, R.Q, R.eQ, R.gainE, R.eGainE, R.effN > 0 ? cvar[c].effE / std::sqrt(R.effN) : 0.,
          R.win.n[1], R.win.n[2], R.win.n[1] > 0 ? R.win.n[2] / R.win.n[1] : 0., R.win.n[1] + R.win.n[2],
          R.win.sig() > 0 ? -100 * (R.win.n[1] + R.win.n[2] - R.win.cb) / R.win.sig() : 0.);
    }
    P("\nD. Tagged eta (m5c, pi0 window %.3f-%.3f). eff_T = T_pass / T(no cut): err = statistical, spread = half of (max - min) over 18\n", kPi0Lo, kPi0Hi);
    P("   variations of the peak fit applied to both selections. close = (T_pass + T_fail) / T(no cut). gain_T = eff_T / sqrt(eff_N).\n");
    P("   %-18s | %6s %5s %8s %6s | %6s %5s %6s %5s | %6s %5s %8s | %6s %5s %6s\n", "selection", "T_pass", "err", "chi2/ndf", "S/B",
      "eff_T", "err", "spread", "close", "T_fail", "err", "chi2/ndf", "gain_T", "err", "spread");
    for (int c = 0; c < nCand; ++c) {
        const Res& R = cr[c];
        P("   %-18s | %6.0f %5.0f %5.1f/%-2d %6.3f | %6.3f %5.3f %6.3f %5.3f | %6.0f %5.0f %5.1f/%-2d | %6.3f %5.3f %6.3f\n", cand[c].name.c_str(),
          R.tp.S, R.tp.eS, R.tp.chi2, R.tp.ndf, R.tp.SoB, R.effT, R.eEffT, cvar[c].effT, R.closeT, R.tf.S, R.tf.eS, R.tf.chi2, R.tf.ndf,
          R.gainT, R.eGainT, R.effN > 0 ? cvar[c].effT / std::sqrt(R.effN) : 0.);
    }
    P("\nE. Two folds (event groups of even / odd index)\n");
    P("   %-18s | %6s %5s %6s %5s %6s %5s | %6s %5s %6s %5s %6s %5s | %6s %5s %6s %5s\n", "selection", "E(A)", "err", "eff_E", "err", "gain_E", "err",
      "E(B)", "err", "eff_E", "err", "gain_E", "err", "gT(A)", "err", "gT(B)", "err");
    for (int c = 0; c < nCand; ++c) {
        const Res& a = cr[nCand + c]; const Res& b = cr[2 * nCand + c];
        P("   %-18s | %6.0f %5.0f %6.3f %5.3f %6.3f %5.3f | %6.0f %5.0f %6.3f %5.3f %6.3f %5.3f | %6.3f %5.3f %6.3f %5.3f\n", cand[c].name.c_str(),
          a.E, a.eE, a.effE, a.eEffE, a.gainE, a.eGainE, b.E, b.eE, b.effE, b.eEffE, b.gainE, b.eGainE, a.gainT, a.eGainT, b.gainT, b.eGainT);
    }
    P("   no cut, step height per slice (0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40): fold A");
    for (int i = 0; i < kNSig; ++i) P(" %4.0f+-%-3.0f", cr[nCand].fp[i].h, cr[nCand].fp[i].eh);
    P("; fold B");
    for (int i = 0; i < kNSig; ++i) P(" %4.0f+-%-3.0f", cr[2 * nCand].fp[i].h, cr[2 * nCand].fp[i].eh);
    P("\n   E(A) - E(B) without cut: %.0f +- %.0f\n", cr[nCand].E - cr[2 * nCand].E, std::sqrt(cr[nCand].eE * cr[nCand].eE + cr[2 * nCand].eE * cr[2 * nCand].eE));

    // ---- each cut varied around the nominal selection ----
    P("\nF. Each cut varied with the other cuts of the nominal selection in place (base = nominal without this cut).\n");
    std::vector<Scan1D> scN;
    for (int c = 0; c < kNCut; ++c) {
        WP base = Nominal(); base.name = "nominal without the scanned cut";
        CutRef(base, c) = -1;
        scN.push_back(RunScan(dE, dG, base, c, false));
        PrintScan(scN.back());
    }
    P("\nPlateaus around the nominal selection:\n");
    for (int c = 0; c < kNCut; ++c) { Plateau(scN[c], false); Plateau(scN[c], true); }

    // ---- isBest == 1 against all combinations ----
    std::vector<Res> rBest;
    {
        std::vector<Job> bj;
        const WP bs[4] = {None(), Loose(), Nominal(), Pure()};
        for (int i = 0; i < 4; ++i) { Job j; j.base = bs[i]; j.w = bs[i]; j.w.best = 1; bj.push_back(j); }
        rBest = Measure(dE, &dG, bj);
        P("\nG. isBest == 1 against all combinations (isBest >= 0). share = measure with isBest == 1 rows only / measure with all rows;\n");
        P("   'extra' = the same fit on the rows with isBest == 0. Q ratio = Q(isBest == 1) / Q(all): from E, and from the tagged eta (share_T / sqrt(share_N)).\n");
        P("   %-8s | %-13s %-13s %-13s %-15s | %-13s %-13s %-13s %-15s | %7s %7s %7s | %-13s %-13s | %-15s %-15s\n", "selection", "E all", "E isBest==1", "E extra", "share",
          "T all", "T isBest==1", "T extra", "share", "W all", "W best", "share", "Q all", "Q isBest==1", "Q ratio (E)", "Q ratio (T)");
        for (int i = 0; i < 4; ++i) {
            const Res& R = rBest[i];
            const double q0 = R.win0.n[0] > 0 ? R.E0 / std::sqrt(R.win0.n[0]) : 0, eq0 = R.win0.n[0] > 0 ? R.eE0 / std::sqrt(R.win0.n[0]) : 0;
            P("   %-8s | %5.0f +- %-4.0f %5.0f +- %-4.0f %5.0f +- %-4.0f %5.3f +- %5.3f | %5.0f +- %-4.0f %5.0f +- %-4.0f %5.0f +- %-4.0f %5.3f +- %5.3f | %7.0f %7.0f %7.3f | %4.2f +- %4.2f  %4.2f +- %4.2f  | %5.3f +- %5.3f %5.3f +- %5.3f\n",
              bs[i].name.c_str(), R.E0, R.eE0, R.E, R.eE, R.Ef, R.eEf, R.effE, R.eEffE, R.t0.S, R.t0.eS, R.tp.S, R.tp.eS, R.tf.S, R.tf.eS, R.effT, R.eEffT,
              R.win0.sig(), R.win.sig(), R.effW, q0, eq0, R.Q, R.eQ, R.gainE, R.eGainE, R.gainT, R.eGainT);
        }
        P("   CB share in the bump window, all / isBest == 1:");
        for (int i = 0; i < 4; ++i) P("  %s %.3f / %.3f", bs[i].name.c_str(), rBest[i].win0.cbShare(), rBest[i].win.cbShare());
        P("\n");
    }

    // ---- cor_ok fallback variant ----
    std::vector<Res> rFb;
    {
        std::vector<Job> fj;
        { Job j; j.w = None(); j.base = None(); fj.push_back(j); }
        { Job j; j.w = Nominal(); j.base = None(); fj.push_back(j); }
        dE.SetFallback(true); dG.SetFallback(true);
        rFb = Measure(dE, &dG, fj);
        dE.SetFallback(false); dG.SetFallback(false);
        P("\nH. Fallback variant: rows with cor_ok == 0 enter with the reconstructed pi- momentum (main analysis: cor_ok == 1 only).\n");
        P("   %-10s %-9s | %6s %5s %8s | %6s %5s %6s | %7s %6s %7s %6s | %5s %5s | %6s %5s %6s %5s\n", "selection", "variant", "E", "err", "chi2/ndf",
          "T", "err", "S/B", "N", "CB", "W=N-CB", "CB/N", "Q", "err", "eff_E", "err", "eff_T", "err");
        for (int i = 0; i < 2; ++i) {
            const Res& a = cr[i == 0 ? 0 : 2]; const Res& b = rFb[i];
            const Res* rr[2] = {&a, &b};
            for (int v = 0; v < 2; ++v) {
                const Res& R = *rr[v];
                P("   %-10s %-9s | %6.0f %5.0f %4.0f/%-3d | %6.0f %5.0f %6.3f | %7.0f %6.0f %7.0f %6.3f | %5.2f %5.2f | %6.3f %5.3f %6.3f %5.3f\n", i ? "nominal" : "none",
                  v ? "fallback" : "cor_ok", R.E, R.eE, R.chi2, R.ndf, R.tp.S, R.tp.eS, R.tp.SoB, R.win.n[0], R.win.cb, R.win.sig(), R.win.cbShare(), R.Q, R.eQ,
                  R.effE, R.eEffE, R.effT, R.eEffT);
            }
        }
    }
    // ---- bias control of the kinematic candidate: fixed procedure on one fold, measured on the other ----
    {
        P("\nI. Bias control of the kinematic gain (two folds = event groups of even / odd index).\n");
        P("   Fixed procedure applied to ONE fold: 1D scans from 'no cut' of oa_lab (upper), M(pipi) (upper) and M(pipi) (lower); for each cut the plateau =\n");
        P("   thresholds around the best gain_T whose gain_T lies within one statistical error of the best; chosen value = centre of the plateau (the looser\n");
        P("   value for an even number). The three values together are then measured on the fold that chose them and on the other fold.\n");
        const int cuts[3] = {kOaLab, kMpHi, kMpLo};
        std::vector<Job> jobs;
        for (int f = 0; f < 2; ++f)
            for (int q = 0; q < 3; ++q)
                for (double t : Thresholds(cuts[q])) {
                    Job j; j.base = None(); j.w = None(); CutRef(j.w, cuts[q]) = t; j.fold = f; j.doFail = false;
                    jobs.push_back(j);
                }
        const std::vector<Res> r = Measure(dE, &dG, jobs);
        WP ch[2];
        size_t k = 0;
        for (int f = 0; f < 2; ++f) {
            ch[f] = None(); ch[f].name = Form("chosen on fold %c", f ? 'B' : 'A');
            for (int q = 0; q < 3; ++q) {
                const std::vector<double> thr = Thresholds(cuts[q]);
                const int n = (int)thr.size();
                int ib = 0;
                for (int i = 0; i < n; ++i) if (r[k + i].gainT > r[k + ib].gainT) ib = i;
                const double lim = r[k + ib].gainT - r[k + ib].eGainT;
                int lo = ib, hi = ib;
                while (lo > 0 && r[k + lo - 1].gainT >= lim) --lo;
                while (hi < n - 1 && r[k + hi + 1].gainT >= lim) ++hi;
                const bool looseLast = thr.back() < 0;            // "no cut" is the last threshold (upper cuts) or the first
                const int s2 = lo + hi;
                const int ic = (s2 % 2 == 0) ? s2 / 2 : (looseLast ? (s2 + 1) / 2 : (s2 - 1) / 2);
                CutRef(ch[f], cuts[q]) = thr[ic];
                P("   fold %c, %-14s: best gain_T %.3f +- %.3f at %s; plateau %s .. %s; chosen %s\n", f ? 'B' : 'A', kCutName[cuts[q]], r[k + ib].gainT,
                  r[k + ib].eGainT, ThrName(cuts[q], thr[ib]).c_str(), ThrName(cuts[q], thr[lo]).c_str(), ThrName(cuts[q], thr[hi]).c_str(),
                  ThrName(cuts[q], thr[ic]).c_str());
                k += n;
            }
        }
        std::vector<Job> ej;
        const WP sel3[3] = {ch[0], ch[1], Nominal()};
        for (int s3 = 0; s3 < 3; ++s3)
            for (int f = 0; f < 2; ++f) { Job j; j.base = None(); j.w = sel3[s3]; j.fold = f; ej.push_back(j); }
        const std::vector<Res> e = Measure(dE, &dG, ej);
        P("   %-22s %-34s | %-15s %-15s | %-15s %-15s\n", "selection", "cuts", "gain_T fold A", "gain_T fold B", "gain_E fold A", "gain_E fold B");
        for (int s3 = 0; s3 < 3; ++s3)
            P("   %-22s %-34s | %5.3f +- %5.3f  %5.3f +- %5.3f  | %5.3f +- %5.3f  %5.3f +- %5.3f\n", s3 < 2 ? sel3[s3].name.c_str() : "nominal (all rows)",
              Describe(sel3[s3]).c_str(), e[2 * s3].gainT, e[2 * s3].eGainT, e[2 * s3 + 1].gainT, e[2 * s3 + 1].eGainT, e[2 * s3].gainE, e[2 * s3].eGainE,
              e[2 * s3 + 1].gainE, e[2 * s3 + 1].eGainE);
        const double bT = 0.5 * ((e[0].gainT - e[1].gainT) + (e[3].gainT - e[2].gainT));
        const double bE = 0.5 * ((e[0].gainE - e[1].gainE) + (e[3].gainE - e[2].gainE));
        const double uT = 0.5 * (e[1].gainT + e[2].gainT), uE = 0.5 * (e[1].gainE + e[2].gainE);
        const double euT = 0.5 * std::sqrt(e[1].eGainT * e[1].eGainT + e[2].eGainT * e[2].eGainT);
        const double euE = 0.5 * std::sqrt(e[1].eGainE * e[1].eGainE + e[2].eGainE * e[2].eGainE);
        P("   optimisation bias b = 1/2 [(g_A(c_A) - g_B(c_A)) + (g_B(c_B) - g_A(c_B))] (the difference of the two folds cancels): gain_T %+.3f, gain_E %+.3f\n", bT, bE);
        P("   gain measured on the fold that did NOT choose (mean of the two): gain_T = %.3f +- %.3f, gain_E = %.3f +- %.3f\n", uT, euT, uE, euE);
        P("   READING: the kinematic gain is supported by the tagged sample (gain_T). From the edge step (gain_E) it is not separated from 1: gain_E also\n");
        P("   has a model spread of 2-5 %% (section C) and cannot see M(pipi) cuts outside 0.30-0.40.\n");
    }
    Close();

    // ========================================================================
    // 6. spectra for the figures, ROOT file
    // ========================================================================
    // M4 spectra of the candidates; slices of the nominal selection (all and isBest == 1)
    std::vector<Spectra> candSp;
    std::vector<Spectra> nomSl, nomSlBest;
    std::vector<StepFit> nomFit, nomFitBest;
    {
        std::vector<Sel> sels;
        for (int c = 0; c < nCand; ++c) { const WP w = cand[c]; sels.push_back([=](const Row& r) { return Pass(w, r); }); }
        candSp = scan::FillMany(dE, sels);
        sels.clear();
        for (int v = 0; v < 2; ++v)
            for (int i = 0; i < kNSl; ++i) {
                WP w = Nominal(); w.best = v;
                if (i == kNSl - 1) w.mpHi = -1;      // control slice: the nominal cuts without the upper M(pipi) cut
                const float lo = kSlLo[i], hi = kSlHi[i];
                sels.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && Pass(w, r); });
            }
        const std::vector<Spectra> s = scan::FillMany(dE, sels);
        for (int i = 0; i < kNSl; ++i) {
            nomSl.push_back(s[i]); nomSlBest.push_back(s[kNSl + i]);
            nomFit.push_back(scan::FitEdgeStep(s[i], SliceCfg(i, kSlLo[i], kSlHi[i])));
            nomFitBest.push_back(scan::FitEdgeStep(s[kNSl + i], SliceCfg(i, kSlLo[i], kSlHi[i])));
        }
    }
    // tagged eta: pass / fail spectra of the main cuts; M4 of the tagged eta events
    const int nPF = 7;
    WP pfWP[nPF];
    const char* pfName[nPF] = {"oa_lab < 40#circ", "M(#pi#pi) > 0.30", "M(#pi#pi) < 0.41", "OA > 2#circ", "OA < 20#circ", "mm4 < 2.80", "#eta-rest > 140#circ"};
    pfWP[0].oaLab = 40; pfWP[1].mpLo = 0.30; pfWP[2].mpHi = 0.41; pfWP[3].oaLo = 2; pfWP[4].oaHi = 20; pfWP[5].mm4 = 2.80; pfWP[6].etaRest = 140;
    std::vector<Spectra> pfSp;
    Spectra m4Peak, m4Side;
    {
        std::vector<Sel> sels;
        for (int i = 0; i < nPF; ++i) {
            const WP w = pfWP[i];
            sels.push_back([=](const Row& r) { return Tagged(r) && Pass(w, r); });
            sels.push_back([=](const Row& r) { return Tagged(r) && !Pass(w, r); });
        }
        pfSp = scan::FillMany(dG, sels, TagOpt());
        const double a1 = kTagMean - 1.9 * kTagSigma, a2 = kTagMean + 1.9 * kTagSigma;
        const double b1 = kTagMean - 5.7 * kTagSigma, b2 = kTagMean - 3.8 * kTagSigma, b3 = kTagMean + 3.8 * kTagSigma, b4 = kTagMean + 5.7 * kTagSigma;
        std::vector<Sel> s2;
        s2.push_back([=](const Row& r) { const double m = M5c(r); return Tagged(r) && m >= a1 && m < a2; });
        s2.push_back([=](const Row& r) { const double m = M5c(r); return Tagged(r) && ((m >= b1 && m < b2) || (m >= b3 && m < b4)); });
        scan::FillOpt o4; o4.nbin = 60; o4.lo = 0.35; o4.hi = 0.65; o4.reb = 1;
        const std::vector<Spectra> s4 = scan::FillMany(dG, s2, o4);
        m4Peak = s4[0]; m4Side = s4[1];
    }
    {
        TFile fo(Form("%s/bump_scan.root", scan::kScanDir), "RECREATE");
        for (int c = 0; c < nCand; ++c)
            for (const char* w : {"us", "cb", "sig"}) {
                TH1D* h = candSp[c].Hist(w, Form("m4_%s_cand%d", w, c));
                h->SetTitle(Form("%s, %s;%s;%s", w, cand[c].name.c_str(), kXt, kYt)); h->Write(); delete h;
            }
        for (int i = 0; i < kNSl; ++i) {
            TH1D* h = refSl[0][i].Hist("sig", Form("slice%d_none", i)); h->SetTitle(Form("N-CB, no cut, M(pipi) %s;%s;%s", kSlName[i], kXt, kYt)); h->Write(); delete h;
            h = nomSl[i].Hist("sig", Form("slice%d_nominal", i)); h->SetTitle(Form("N-CB, nominal, M(pipi) %s;%s;%s", kSlName[i], kXt, kYt)); h->Write(); delete h;
        }
        for (const char* w : {"us", "cb", "sig"}) {
            TH1D* h = spTag.Hist(w, Form("m5c_%s_tagged", w)); h->SetTitle(Form("%s, tagged rows;m5c [GeV/c^{2}];rows / 2.5 MeV/c^{2}", w)); h->Write(); delete h;
            h = spMeeg.Hist(w, Form("meeg_%s", w)); h->SetTitle(Form("%s, rows with a photon;M(e^{+}e^{-}#gamma) [GeV/c^{2}];rows / 2 MeV/c^{2}", w)); h->Write(); delete h;
        }
        for (int c = 0; c < kNCut; ++c) {
            TGraphErrors gE, gT, gW;
            for (size_t i = 0; i < sc1[c].thr.size(); ++i) {
                const Res& R = sc1[c].all[i];
                gE.SetPoint((int)i, sc1[c].thr[i], R.E); gE.SetPointError((int)i, 0, R.eE);
                gT.SetPoint((int)i, sc1[c].thr[i], R.effT); gT.SetPointError((int)i, 0, R.eEffT);
                gW.SetPoint((int)i, sc1[c].thr[i], R.win.sig()); gW.SetPointError((int)i, 0, R.win.sigErr());
            }
            gE.Write(Form("E_vs_%s", kCutTag[c])); gT.Write(Form("effT_vs_%s", kCutTag[c])); gW.Write(Form("W_vs_%s", kCutTag[c]));
        }
        fo.Close();
    }

    // ========================================================================
    // 7. figures
    // ========================================================================
    Style();

    // ---- 7.1 slices with the fitted steps ----
    for (int v = 0; v < 2; ++v) {
        TCanvas* c = new TCanvas(Form("c_sl%d", v), "", 2400, 1500);
        c->Divide(3, 2, 0.002, 0.002);
        for (int i = 0; i < kNSl; ++i) {
            c->cd(i + 1);
            const Spectra& sp = v ? nomSl[i] : refSl[0][i];
            const StepFit& f = v ? nomFit[i] : refFit[0][i];
            const char* sel = v ? (i == kNSl - 1 ? "oa_lab < 40#circ" : "nominal") : "no cut";
            DrawSlice(sp, f, Form("hsl_%d_%d", v, i), Form("N - CB, %s, %s%.2f < M(#pi^{+}#pi^{-}) < %.2f GeV/c^{2}", sel,
                      i == kNSl - 1 ? "control: " : "", kSlLo[i], kSlHi[i]), i == 4 ? 1 : (i == 5 ? 2 : 0));
        }
        Save(c, v ? "bump_edge_slices_nominal" : "bump_edge_slices_none", true);
    }

    // ---- 7.2 - 7.4 measures against the thresholds ----
    for (int v = 0; v < 2; ++v) {
        const std::vector<Scan1D>& sc = v ? scN : sc1;
        const char* tag = v ? "nominal" : "none";
        const char* what = v ? "other cuts: nominal selection" : "one cut at a time, no other cut";
        // E
        {
            TCanvas* c = new TCanvas(Form("c_E%d", v), "", 2600, 1300);
            c->Divide(4, 2, 0.002, 0.002);
            for (int k = 0; k < kNCut; ++k) {
                c->cd(k + 1);
                const Scan1D& s = sc[k];
                const int n = (int)s.thr.size();
                const double e0 = s.all[0].E0, ee0 = s.all[0].eE0;
                TH1F* fr = Frame(Form("frE_%d_%d", v, k), k, s.thr, 0, 1.75 * e0, "step-height sum E [rows / 2.5 MeV] (not a yield)");
                fr->Draw();
                TBox* bx = new TBox(0, e0 - ee0, n, e0 + ee0); bx->SetFillColorAlpha(kBlue - 9, 0.45); bx->SetLineColor(0); bx->Draw();
                HLine(0, n, e0, kBlue - 7, 2, 1);
                TGraphErrors* g = Gr(kBlue + 1, 20, 1.5); TGraphErrors* ga = Gr(kBlue - 4, 26, 1.2); TGraphErrors* gb = Gr(kBlue - 4, 32, 1.2);
                for (int i = 0; i < n; ++i) {
                    Add(g, i + 0.5, s.all[i].E, s.all[i].eE);
                    if (!s.f0.empty()) { Add(ga, i + 0.32, 2 * s.f0[i].E, 2 * s.f0[i].eE); Add(gb, i + 0.68, 2 * s.f1[i].E, 2 * s.f1[i].eE); }
                }
                if (!s.f0.empty()) { ga->Draw("p same"); gb->Draw("p same"); }
                g->Draw("p same");
                if (k == kEtaRest) Txt(0.21, 0.20, "120#circ, 140#circ: edge fits fail (#chi^{2}), E not valid", 0.040, kRed + 1);
                if (k == kMpHi || k == kMpLo) Txt(0.21, 0.20, "cuts outside 0.30-0.40 cannot change E", 0.040, kRed + 1);
                Txt(0.17, 0.935, Form("edge-step sum E, scanned cut: %s", kCutName[k]), 0.045);
                gPad->RedrawAxis();
            }
            c->cd(8);
            Txt(0.05, 0.90, "p+p 4.53 GeV, all combinations (isBest #geq 0),", 0.042);
            Txt(0.05, 0.84, "corrected momenta", 0.042);
            Txt(0.05, 0.76, what, 0.045);
            Txt(0.05, 0.68, "E = sum of the step heights (not a yield) at the predicted", 0.042);
            Txt(0.05, 0.62, "kinematic edge in the M(#pi^{+}#pi^{-}) slices", 0.045);
            Txt(0.05, 0.56, "0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40", 0.045);
            TLegend* lg = new TLegend(0.05, 0.22, 0.95, 0.50); lg->SetTextSize(0.045);
            TGraphErrors* l1 = Gr(kBlue + 1, 20, 1.5); TGraphErrors* l2 = Gr(kBlue - 4, 26, 1.2); TGraphErrors* l3 = Gr(kBlue - 4, 32, 1.2);
            TBox* lb = new TBox(); lb->SetFillColorAlpha(kBlue - 9, 0.45); lb->SetLineColor(0);
            lg->AddEntry(l1, "all event groups", "pe");
            if (!sc[0].f0.empty()) { lg->AddEntry(l2, "2 #times fold A (even groups)", "pe"); lg->AddEntry(l3, "2 #times fold B (odd groups)", "pe"); }
            lg->AddEntry(lb, "without this cut, #pm stat. error", "f");
            lg->Draw();
            Save(c, Form("bump_E_vs_cuts_%s", tag), true);
        }
        // W
        {
            TCanvas* c = new TCanvas(Form("c_W%d", v), "", 2600, 1500);
            c->Divide(4, 2, 0.002, 0.002);
            for (int k = 0; k < kNCut; ++k) {
                c->cd(k + 1);
                TPad* up = new TPad(Form("upW_%d_%d", v, k), "", 0, 0.42, 1, 1); TPad* dn = new TPad(Form("dnW_%d_%d", v, k), "", 0, 0, 1, 0.42);
                up->SetBottomMargin(0.02); up->SetTopMargin(0.11); dn->SetTopMargin(0.03); dn->SetBottomMargin(0.36);
                up->Draw(); dn->Draw();
                const Scan1D& s = sc[k];
                const int n = (int)s.thr.size();
                up->cd(); up->SetLogy();
                TH1F* fr = Frame(Form("frW_%d_%d", v, k), k, s.thr, 300, 1.5e6, "rows in the bump window");
                fr->GetXaxis()->SetLabelSize(0); fr->GetXaxis()->SetTitleSize(0); fr->GetYaxis()->SetTitleSize(0.07); fr->GetYaxis()->SetLabelSize(0.065);
                fr->GetYaxis()->SetTitleOffset(1.0);
                fr->Draw();
                TGraphErrors* gn = Gr(kBlack, 20, 1.3); TGraphErrors* gc = Gr(kRed + 1, 20, 1.3); TGraphErrors* gs = Gr(kBlue + 1, 24, 1.3);
                TGraphErrors* gf = Gr(kRed + 1, 20, 1.3);
                for (int i = 0; i < n; ++i) {
                    const WinSum& w = s.all[i].win;
                    Add(gn, i + 0.5, w.n[0], std::sqrt(w.var[0])); Add(gc, i + 0.5, w.cb, std::sqrt(w.cbVar)); Add(gs, i + 0.5, w.sig(), w.sigErr());
                    Add(gf, i + 0.5, 100 * w.cbShare(), w.n[0] > 0 ? 100 * std::sqrt(w.cbVar) / w.n[0] : 0);
                }
                gn->Draw("pl same"); gc->Draw("pl same"); gs->Draw("pl same");
                Txt(0.17, 0.92, Form("bump window %.3f-%.4f, scanned cut: %s", scan::kBumpLo, scan::kBumpHi, kCutName[k]), 0.062);
                dn->cd();
                TH1F* f2 = Frame(Form("frW2_%d_%d", v, k), k, s.thr, 0, 34, "CB / N [%]");
                f2->GetXaxis()->SetLabelSize(0.13); f2->GetXaxis()->SetTitleSize(0.085); f2->GetXaxis()->SetTitleOffset(1.6);
                f2->GetYaxis()->SetTitleSize(0.095); f2->GetYaxis()->SetLabelSize(0.09); f2->GetYaxis()->SetTitleOffset(0.72); f2->GetYaxis()->SetNdivisions(505);
                f2->Draw();
                gf->Draw("pl same");
            }
            c->cd(8);
            Txt(0.05, 0.90, "p+p 4.53 GeV, all combinations (isBest #geq 0),", 0.042);
            Txt(0.05, 0.84, "corrected momenta", 0.042);
            Txt(0.05, 0.76, what, 0.045);
            Txt(0.05, 0.69, Form("rows with %.3f < M(#pi^{+}#pi^{-}e^{+}e^{-}) < %.4f GeV/c^{2}", scan::kBumpLo, scan::kBumpHi), 0.042);
            TLegend* lg = new TLegend(0.05, 0.30, 0.95, 0.62); lg->SetTextSize(0.045);
            TGraphErrors* l1 = Gr(kBlack, 20, 1.3); TGraphErrors* l2 = Gr(kRed + 1, 20, 1.3); TGraphErrors* l3 = Gr(kBlue + 1, 24, 1.3);
            lg->AddEntry(l1, "all pairs e^{+}e^{-} (N)", "pl"); lg->AddEntry(l2, "CB = 2#sqrt{N_{++}N_{--}}", "pl");
            lg->AddEntry(l3, "W = N - CB", "pl"); lg->AddEntry(l2, "lower panels: CB share of N (no k factor)", "pl");
            lg->Draw();
            Save(c, Form("bump_W_vs_cuts_%s", tag), true);
        }
        // efficiencies and gains
        {
            TCanvas* c = new TCanvas(Form("c_F%d", v), "", 2600, 1600);
            c->Divide(4, 2, 0.002, 0.002);
            for (int k = 0; k < kNCut; ++k) {
                c->cd(k + 1);
                TPad* up = new TPad(Form("upF_%d_%d", v, k), "", 0, 0.45, 1, 1); TPad* dn = new TPad(Form("dnF_%d_%d", v, k), "", 0, 0, 1, 0.45);
                up->SetBottomMargin(0.02); up->SetTopMargin(0.11); dn->SetTopMargin(0.03); dn->SetBottomMargin(0.34);
                up->Draw(); dn->Draw();
                const Scan1D& s = sc[k];
                const int n = (int)s.thr.size();
                up->cd();
                TH1F* fr = Frame(Form("frF_%d_%d", v, k), k, s.thr, 0, 1.3, "fraction kept by the cut");
                fr->GetXaxis()->SetLabelSize(0); fr->GetXaxis()->SetTitleSize(0); fr->GetYaxis()->SetTitleSize(0.07); fr->GetYaxis()->SetLabelSize(0.065);
                fr->GetYaxis()->SetTitleOffset(1.0);
                fr->Draw();
                HLine(0, n, 1, kGray + 1, 2, 1);
                TGraphErrors* gn = Gr(kBlack, 1); TGraphErrors* gc = Gr(kRed + 1, 1); TGraphErrors* gw = Gr(kBlue + 1, 1);
                gn->SetLineStyle(2); gc->SetLineStyle(2); gw->SetLineStyle(1);
                TGraphErrors* gt = Gr(kBlue + 1, 20, 1.4); TGraphErrors* ge = Gr(kBlue + 1, 25, 1.4);
                TGraphErrors* qt = Gr(kBlue + 1, 20, 1.4); TGraphErrors* qe = Gr(kBlue + 1, 25, 1.4);
                for (int i = 0; i < n; ++i) {
                    const Res& R = s.all[i];
                    Add(gn, i + 0.5, R.win0.n[0] > 0 ? R.win.n[0] / R.win0.n[0] : 0, 0);
                    Add(gc, i + 0.5, R.win0.cb > 0 ? R.win.cb / R.win0.cb : 0, 0);
                    Add(gw, i + 0.5, R.effW, 0);
                    Add(gt, i + 0.40, R.effT, R.eEffT); Add(ge, i + 0.60, R.effE, R.eEffE);
                    Add(qt, i + 0.40, R.gainT, R.eGainT); Add(qe, i + 0.60, R.gainE, R.eGainE);
                }
                gn->Draw("l same"); gc->Draw("l same"); gw->Draw("l same"); ge->Draw("p same"); gt->Draw("p same");
                Txt(0.17, 0.92, Form("scanned cut: %s", kCutName[k]), 0.07);
                if (k == kEtaRest) Txt(0.21, 0.10, "squares not valid at 120#circ, 140#circ (fits fail)", 0.055, kRed + 1);
                if (k == kMpHi || k == kMpLo) Txt(0.21, 0.10, "cuts outside 0.30-0.40 cannot change E", 0.055, kRed + 1);
                dn->cd();
                TH1F* f2 = Frame(Form("frF2_%d_%d", v, k), k, s.thr, 0.55, 1.30, "gain");
                f2->GetXaxis()->SetLabelSize(0.12); f2->GetXaxis()->SetTitleSize(0.08); f2->GetXaxis()->SetTitleOffset(1.6);
                f2->GetYaxis()->SetTitleSize(0.08); f2->GetYaxis()->SetLabelSize(0.08); f2->GetYaxis()->SetTitleOffset(0.85); f2->GetYaxis()->SetNdivisions(505);
                f2->Draw();
                HLine(0, n, 1, kGray + 1, 2, 1);
                qe->Draw("p same"); qt->Draw("p same");
            }
            c->cd(8);
            Txt(0.05, 0.93, "p+p 4.53 GeV, all combinations (isBest #geq 0),", 0.042);
            Txt(0.05, 0.88, "corrected momenta", 0.042);
            Txt(0.05, 0.82, what, 0.045);
            TLegend* lg = new TLegend(0.05, 0.36, 0.98, 0.78); lg->SetTextSize(0.040);
            TGraphErrors* l1 = Gr(kBlue + 1, 20, 1.4); TGraphErrors* l2 = Gr(kBlue + 1, 25, 1.4);
            TGraphErrors* l3 = Gr(kBlue + 1, 1); TGraphErrors* l4 = Gr(kBlack, 1); TGraphErrors* l5 = Gr(kRed + 1, 1);
            l4->SetLineStyle(2); l5->SetLineStyle(2);
            lg->AddEntry(l1, "tagged #eta peak (ECAL photon): share kept", "pe");
            lg->AddEntry(l2, "edge-step sum E: share kept", "pe");
            lg->AddEntry(l3, "W = N - CB in the bump window", "l");
            lg->AddEntry(l4, "N (all pairs) in the bump window", "l");
            lg->AddEntry(l5, "CB in the bump window", "l");
            lg->Draw();
            Txt(0.05, 0.29, "lower panels: gain = share kept / #sqrt{share of N kept},", 0.040);
            Txt(0.05, 0.23, "for E (squares) and for the tagged #eta (circles).", 0.040);
            Txt(0.05, 0.17, "A cut helps when the gain is above 1 by more", 0.040);
            Txt(0.05, 0.11, "than its error.", 0.040);
            Save(c, Form("bump_eff_vs_cuts_%s", tag), true);
        }
    }

    // ---- 7.5 tagged eta peak ----
    {
        const PeakCfg tc = TagCfg();
        const PeakFit f = scan::FitPeakMain(spTag, tc);
        TCanvas* c = new TCanvas("c_tag", "", 2400, 1700);
        c->Divide(2, 2, 0.002, 0.002);
        c->cd(1);
        TH1D* hn = spTag.Hist("us", "tg_us"); TH1D* hc = spTag.Hist("cb", "tg_cb"); TH1D* hs = spTag.Hist("sig", "tg_sig");
        Deco(hn, kBlack, 20); Deco(hc, kRed + 1, 20); Deco(hs, kBlue + 1, 20);
        hn->GetXaxis()->SetTitle("m5c = M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma), M(e^{+}e^{-}#gamma) fixed to m_{#pi^{0}} [GeV/c^{2}]");
        hn->GetXaxis()->SetTitleSize(0.045); hn->GetYaxis()->SetTitle(kYt);
        hn->GetXaxis()->SetRangeUser(0.4501, 0.6999); hn->SetMinimum(0); hn->SetMaximum(1.25 * hn->GetBinContent(hn->FindBin(0.699)));
        hn->Draw("e"); hc->Draw("e same"); hs->Draw("e same");
        TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
        for (int i = 0; i <= 300; ++i) { const double m = tc.lo + (tc.hi - tc.lo) * i / 300; gt->SetPoint(i, m, f.Total(m)); gb->SetPoint(i, m, f.Bkg(m)); }
        gt->SetLineWidth(2); gb->SetLineWidth(2); gb->SetLineStyle(2); gb->Draw("l same"); gt->Draw("l same");
        Txt(0.17, 0.935, Form("tagged rows: %.3f < M(e^{+}e^{-}#gamma) < %.3f GeV/c^{2}, no other cut", kPi0Lo, kPi0Hi), 0.040);
        TLegend* lg = new TLegend(0.20, 0.62, 0.64, 0.88); lg->SetTextSize(0.042);
        lg->AddEntry(hn, "all pairs e^{+}e^{-} (N)", "pe"); lg->AddEntry(hc, "CB = 2#sqrt{N_{++}N_{--}}", "pe"); lg->AddEntry(hs, "N - CB", "pe");
        lg->AddEntry(gt, "fit: Gaussian + pol2", "l"); lg->Draw();
        c->cd(2);
        TH1D* hr = (TH1D*)hs->Clone("tg_res");
        for (int i = 1; i <= hr->GetNbinsX(); ++i) hr->SetBinContent(i, hs->GetBinContent(i) - f.Bkg(hr->GetBinCenter(i)));
        hr->GetXaxis()->SetTitle("m5c [GeV/c^{2}]"); hr->GetYaxis()->SetTitle(kYt);
        hr->GetXaxis()->SetRangeUser(tc.lo + 1e-4, tc.hi - 1e-4);
        hr->SetMinimum(-0.35 * f.par[0] * 0.16); hr->SetMaximum(1.6 * f.par[0] * 0.16);
        hr->Draw("e");
        TGraph* gp = new TGraph();
        for (int i = 0; i <= 300; ++i) { const double m = tc.lo + (tc.hi - tc.lo) * i / 300; gp->SetPoint(i, m, f.Total(m) - f.Bkg(m)); }
        gp->SetLineWidth(2); gp->Draw("l same"); HLine(tc.lo, tc.hi, 0, kGray + 1, 2, 1);
        Txt(0.17, 0.935, "the same N - CB after subtracting the fitted polynomial: the tagged #eta peak", 0.038);
        Txt(0.21, 0.84, Form("S = %.0f #pm %.0f (stat) #pm %.0f (model)", f.S, f.eS, tagRef[0].sprd.half), 0.045);
        Txt(0.21, 0.78, Form("frozen shape %.4f, #sigma = %.1f MeV", kTagMean, 1e3 * kTagSigma), 0.045);
        Txt(0.21, 0.72, Form("#chi^{2}/ndf = %.1f/%d,  S/B = %.2f", f.chi2, f.ndf, f.SoB), 0.045);
        c->cd(3);
        {
            PeakCfg c5 = TagCfg(); c5.mean = kOldMean; c5.sigma = kOldSigma; c5.lo = 0.48; c5.hi = 0.64;
            c5.winLo = c5.mean - 1.9 * c5.sigma; c5.winHi = c5.mean + 1.9 * c5.sigma;
            const PeakFit f5 = scan::FitPeakMain(spOld, c5);
            TH1D* h5 = spOld.Hist("sig", "tg_old"); Deco(h5, kBlue + 1, 20);
            h5->GetXaxis()->SetTitle("plain m5 = M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) [GeV/c^{2}]"); h5->GetYaxis()->SetTitle("rows / 5 MeV/c^{2}");
            h5->GetXaxis()->SetRangeUser(0.4501, 0.6999); h5->SetMinimum(0); h5->SetMaximum(1.25 * h5->GetBinContent(h5->FindBin(0.699)));
            h5->Draw("e");
            TGraph* g5 = new TGraph(); TGraph* b5 = new TGraph();
            for (int i = 0; i <= 300; ++i) { const double m = c5.lo + (c5.hi - c5.lo) * i / 300; g5->SetPoint(i, m, f5.Total(m)); b5->SetPoint(i, m, f5.Bkg(m)); }
            g5->SetLineWidth(2); b5->SetLineWidth(2); b5->SetLineStyle(2); b5->Draw("l same"); g5->Draw("l same");
            Txt(0.17, 0.935, Form("the old way: N - CB in the plain m5, %.3f < M(e^{+}e^{-}#gamma) < %.3f GeV/c^{2}", kOldLo, kOldHi), 0.038);
            Txt(0.21, 0.84, Form("S = %.0f #pm %.0f,  #sigma = %.1f MeV", f5.S, f5.eS, 1e3 * kOldSigma), 0.045);
            Txt(0.21, 0.78, Form("#chi^{2}/ndf = %.1f/%d,  S/B = %.2f", f5.chi2, f5.ndf, f5.SoB), 0.045);
        }
        c->cd(4);
        {
            TH1D* hp = m4Peak.Hist("us", "tg_m4"); TH1D* hsb = m4Side.Hist("us", "tg_m4s");
            hp->Add(hsb, -1.0);
            Deco(hp, kBlue + 1, 20);
            hp->GetXaxis()->SetTitle(kXt); hp->GetYaxis()->SetTitle("rows / 5 MeV/c^{2}");
            hp->GetXaxis()->SetRangeUser(0.3501, 0.5999);
            const double ym = 1.55 * hp->GetMaximum();
            hp->SetMaximum(ym); hp->SetMinimum(-0.12 * ym);
            hp->Draw("e");
            HLine(0.35, 0.60, 0, kGray + 1, 2, 1);
            VLine(scan::kBumpLo, 0, 0.62 * ym, kGreen + 2, 7); VLine(scan::kBumpHi, 0, 0.62 * ym, kGreen + 2, 7);
            Txt(0.17, 0.935, "four-particle mass of the tagged #eta events (e^{+}e^{-} pairs, CB not subtracted)", 0.038);
            Txt(0.21, 0.84, "rows in the m5c peak window (#pm1.9#sigma)", 0.040, kBlue + 1);
            Txt(0.21, 0.79, "minus the two side bands (3.8-5.7#sigma);", 0.040, kBlue + 1);
            Txt(0.21, 0.74, "not exact above 0.51 (m5c > M4)", 0.040, kBlue + 1);
            Txt(0.21, 0.69, Form("green: bump window %.3f-%.4f", scan::kBumpLo, scan::kBumpHi), 0.040, kGreen + 2);
        }
        Save(c, "bump_tagged_eta", true);
    }

    // ---- 7.6 pass / fail fits of the tagged peak ----
    {
        const PeakCfg tc = TagCfg();
        TCanvas* c = new TCanvas("c_pf", "", 3500, 1150);
        c->Divide(nPF, 2, 0.001, 0.002);
        for (int i = 0; i < nPF; ++i)
            for (int v = 0; v < 2; ++v) {
                c->cd(v * nPF + i + 1);
                gPad->SetLeftMargin(0.19);
                const Spectra& sp = pfSp[2 * i + v];
                const PeakFit f = scan::FitPeakMain(sp, tc);
                TH1D* h = sp.Hist("sig", Form("pf_%d_%d", i, v)); Deco(h, kBlue + 1, 20, 0.8);
                h->GetXaxis()->SetTitle("m5c [GeV/c^{2}]"); h->GetYaxis()->SetTitle("N - CB [rows / 2.5 MeV/c^{2}]");
                h->GetYaxis()->SetTitleOffset(1.6);
                h->GetXaxis()->SetRangeUser(tc.lo - 0.02 + 1e-4, tc.hi + 0.02 - 1e-4);
                double mx = 0;
                for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinCenter(b) > tc.lo - 0.02 && h->GetBinCenter(b) < tc.hi + 0.02) mx = std::max(mx, h->GetBinContent(b) + h->GetBinError(b));
                h->SetMinimum(-0.08 * mx); h->SetMaximum(1.45 * mx);
                h->Draw("e");
                if (f.ok) {
                    TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
                    for (int q = 0; q <= 300; ++q) { const double m = tc.lo + (tc.hi - tc.lo) * q / 300; gt->SetPoint(q, m, f.Total(m)); gb->SetPoint(q, m, f.Bkg(m)); }
                    gt->SetLineWidth(2); gb->SetLineWidth(2); gb->SetLineStyle(2); gb->Draw("l same"); gt->Draw("l same");
                    Txt(0.23, 0.83, Form("S = %.0f #pm %.0f", f.S, f.eS), 0.055);
                    Txt(0.23, 0.76, Form("#chi^{2}/ndf = %.0f/%d", f.chi2, f.ndf), 0.055);
                    if (i == 2 && v == 1) Txt(0.23, 0.69, "fit not valid", 0.055, kRed + 1);
                }
                Txt(0.19, 0.935, Form("tagged rows %s %s", v ? "FAILING" : "passing", pfName[i]), 0.050);
            }
        Save(c, "bump_tagged_passfail", true);
    }

    // ---- 7.7 M(e+e-gamma) and the pi0 of the tagged eta ----
    {
        TCanvas* c = new TCanvas("c_pi0", "", 2400, 1000);
        c->Divide(2, 1, 0.002, 0.002);
        c->cd(1);
        TH1D* hn = spMeeg.Hist("us", "pi_us"); TH1D* hc = spMeeg.Hist("cb", "pi_cb"); TH1D* hs = spMeeg.Hist("sig", "pi_sig");
        Deco(hn, kBlack, 20); Deco(hc, kRed + 1, 20); Deco(hs, kBlue + 1, 24);
        hn->GetXaxis()->SetTitle("M(e^{+}e^{-}#gamma) [GeV/c^{2}]"); hn->GetYaxis()->SetTitle("rows / 2 MeV/c^{2}");
        gPad->SetLeftMargin(0.20); hn->GetYaxis()->SetTitleOffset(1.9);
        hn->SetMinimum(0); hn->SetMaximum(1.3 * hn->GetMaximum());
        hn->Draw("e"); hc->Draw("e same"); hs->Draw("e same");
        const double ym = hn->GetMaximum();
        VLine(kPi0Lo, 0, 0.7 * ym, kGreen + 2, 1); VLine(kPi0Hi, 0, 0.7 * ym, kGreen + 2, 1);
        VLine(kOldLo, 0, 0.5 * ym, kGray + 2, 2); VLine(kOldHi, 0, 0.5 * ym, kGray + 2, 2);
        Txt(0.20, 0.935, "all rows with a photon (m_{ee} < 0.1 GeV/c^{2}), no other cut", 0.040);
        TLegend* lg = new TLegend(0.67, 0.48, 0.97, 0.88); lg->SetTextSize(0.031); lg->SetFillStyle(0);
        lg->AddEntry(hn, "all pairs e^{+}e^{-} (N)", "pe"); lg->AddEntry(hc, "CB = 2#sqrt{N_{++}N_{--}}", "pe"); lg->AddEntry(hs, "N - CB", "pe");
        TLine* la = new TLine(); la->SetLineColor(kGreen + 2); la->SetLineWidth(2); TLine* lb = new TLine(); lb->SetLineColor(kGray + 2); lb->SetLineStyle(2); lb->SetLineWidth(2);
        lg->AddEntry(la, Form("new: %.3f-%.3f", kPi0Lo, kPi0Hi), "l"); lg->AddEntry(lb, Form("old: %.3f-%.3f", kOldLo, kOldHi), "l");
        lg->Draw();
        Txt(0.23, 0.84, "photon: the one", 0.034); Txt(0.23, 0.80, "closest to 0.140", 0.034); Txt(0.23, 0.76, "of up to three", 0.034);
        c->cd(2);
        TGraphErrors* g0 = Gr(kBlue + 1, 20, 1.3); TGraphErrors* g1 = Gr(kBlue + 1, 24, 1.3);
        double mx = 0;
        for (int i = 0; i < nPi; ++i) {
            Add(g0, piLo + (i + 0.5) * piBw, piS[0][i], piE[0][i]); Add(g1, piLo + (i + 0.5) * piBw + 0.0008, piS[1][i], piE[1][i]);
            mx = std::max(mx, piS[0][i]);
        }
        TH1F* fr = new TH1F("pi_fr", "", 10, 0.06, 0.22); fr->SetDirectory(nullptr);
        fr->GetXaxis()->SetTitle("M(e^{+}e^{-}#gamma) [GeV/c^{2}]"); fr->GetYaxis()->SetTitle("#eta peak content / 5 MeV/c^{2}");
        fr->SetMinimum(-0.12 * mx); fr->SetMaximum(1.5 * mx); fr->Draw();
        HLine(0.06, 0.22, 0, kGray + 1, 2, 1);
        TF1* fg = new TF1("pi_fg", "[0]*0.005/(sqrt(2*pi)*[2])*exp(-0.5*((x-[1])/[2])^2)", 0.095, 0.200);
        fg->SetParameters(pi0Par[0][0], pi0Par[0][1], pi0Par[0][2]); fg->SetLineColor(kBlack); fg->SetLineWidth(2); fg->Draw("same");
        VLine(kPi0Lo, 0, 1.05 * mx, kGreen + 2, 1); VLine(kPi0Hi, 0, 1.05 * mx, kGreen + 2, 1);
        VLine(kOldLo, 0, 0.8 * mx, kGray + 2, 2); VLine(kOldHi, 0, 0.8 * mx, kGray + 2, 2);
        g1->Draw("p same"); g0->Draw("p same");
        Txt(0.17, 0.935, "the #pi^{0} of tagged #eta events: #eta peak content against M(e^{+}e^{-}#gamma)", 0.040);
        Txt(0.20, 0.85, Form("Gaussian fit: %.4f #pm %.4f, #sigma = %.1f #pm %.1f MeV", pi0Par[0][1], pi0Err[0][1], 1e3 * pi0Par[0][2], 1e3 * pi0Err[0][2]), 0.036);
        TLegend* l2 = new TLegend(0.20, 0.70, 0.55, 0.82); l2->SetTextSize(0.034);
        l2->AddEntry(g0, "any number of photons", "pe"); l2->AddEntry(g1, "exactly one photon", "pe"); l2->Draw();
        Txt(0.66, 0.79, "green: new window", 0.034, kGreen + 2); Txt(0.66, 0.75, "(#pm2#sigma)", 0.034, kGreen + 2); Txt(0.66, 0.70, "grey: old window", 0.034, kGray + 2);
        Save(c, "bump_pi0_meeg", true);
    }

    // ---- 7.8 two-dimensional maps ----
    for (int v = 0; v < 2; ++v) {
        const Scan2D& s = v ? s2b : s2a;
        TCanvas* c = new TCanvas(Form("c_2d%d", v), "", 3000, 950);
        c->Divide(3, 1, 0.002, 0.002);
        const char* qn[3] = {"gain in tagged #eta / #sqrt{N}", "gain in E / #sqrt{N}", "share of the tagged #eta peak kept"};
        for (int q = 0; q < 3; ++q) {
            c->cd(q + 1); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.16);
            const int nx = (int)s.tx.size(), ny = (int)s.ty.size();
            TH2D* h = new TH2D(Form("h2_%d_%d", v, q), "", nx, 0, nx, ny, 0, ny); h->SetDirectory(nullptr);
            for (int ix = 0; ix < nx; ++ix) {
                h->GetXaxis()->SetBinLabel(ix + 1, ThrName(s.cx, s.tx[ix]).c_str());
                for (int iy = 0; iy < ny; ++iy) {
                    if (ix == 0) h->GetYaxis()->SetBinLabel(iy + 1, ThrName(s.cy, s.ty[iy]).c_str());
                    const Res& R = s.at(ix, iy);
                    h->SetBinContent(ix + 1, iy + 1, q == 0 ? R.gainT : (q == 1 ? R.gainE : R.effT));
                }
            }
            h->GetXaxis()->SetTitle(kCutAxis[s.cx]); h->GetYaxis()->SetTitle(kCutAxis[s.cy]);
            h->GetXaxis()->SetLabelSize(0.06); h->GetYaxis()->SetLabelSize(0.06); h->GetXaxis()->SetTitleSize(0.042); h->GetYaxis()->SetTitleSize(0.042);
            h->GetYaxis()->SetTitleOffset(1.75); h->GetXaxis()->SetTitleOffset(1.4); h->GetXaxis()->LabelsOption("h");
            h->SetMinimum(q == 2 ? 0.0 : 0.6); h->SetMaximum(q == 2 ? 1.05 : 1.2);
            h->SetMarkerSize(1.5);
            h->Draw("colz text");
            Txt(0.16, 0.935, Form("%s (no other cut)", qn[q]), 0.045);
        }
        Save(c, v ? "bump_2d_oa" : "bump_2d_oalab_mpipi", false);
    }

    // ---- 7.9 edge position and the lower edge ----
    {
        TCanvas* c = new TCanvas("c_pos", "", 2400, 1500);
        c->Divide(4, 2, 0.002, 0.002);
        for (int i = 0; i < kNSig; ++i) {
            c->cd(i + 1);
            const EdgeFree& f = refFree[i];
            TGraph* g = new TGraph();
            double cmin = 1e30;
            for (size_t q = 0; q < f.cs.size(); ++q) cmin = std::min(cmin, f.cs[q]);
            for (size_t q = 0; q < f.cs.size(); ++q) g->SetPoint((int)q, f.xs[q], f.cs[q] - cmin);
            TH1F* fr = new TH1F(Form("pos_fr%d", i), "", 10, f.pred - 0.012, f.pred + 0.012); fr->SetDirectory(nullptr);
            fr->GetXaxis()->SetTitle("edge position at the slice centre [GeV/c^{2}]"); fr->GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{min}"); fr->GetYaxis()->SetTitleOffset(1.2);
            fr->GetXaxis()->SetNdivisions(505); fr->GetXaxis()->SetTitleSize(0.048);
            fr->SetMinimum(0); fr->SetMaximum(30); fr->Draw();
            TBox* bx = new TBox(f.pred * (1 - 0.0009 / scan::kScale), 0, f.pred * (1 + 0.0009 / scan::kScale), 30); bx->SetFillColorAlpha(kGreen - 9, 0.6); bx->SetLineColor(0); bx->Draw();
            VLine(f.pred, 0, 30, kGreen + 2, 1); HLine(f.pred - 0.012, f.pred + 0.012, 1, kGray + 1, 2, 1);
            g->SetLineColor(kBlue + 1); g->SetLineWidth(3); g->Draw("l same");
            Txt(0.17, 0.935, Form("edge position free, no cut, M(#pi^{+}#pi^{-}) %s GeV/c^{2}", kSlName[i]), 0.042);
            if (f.ok) Txt(0.36, 0.84, Form("fit %.4f ^{+%.4f}_{-%.4f}", f.edge, f.eHi * f.pred / scan::kScale, f.eLo * f.pred / scan::kScale), 0.046, kBlue + 1);
            Txt(0.36, 0.77, Form("predicted %.4f", f.pred), 0.042, kGreen + 2);
            Txt(0.36, 0.72, "(band: mass scale error)", 0.036, kGreen + 2);
            gPad->RedrawAxis();
        }
        for (int i = 0; i < kNSig; ++i) {
            c->cd(kNSig + i + 1);
            TH1D* h = new TH1D(Form("low_%d", i), "", 80, -0.05, 0.15); h->SetDirectory(nullptr);
            for (int b = 1; b <= 80; ++b) {
                const int src = lowSp[i].bin(h->GetBinCenter(b) + 0.2);
                if (src < 0) continue;
                h->SetBinContent(b, lowSp[i].sig[src]); h->SetBinError(b, std::sqrt(lowSp[i].sigVar[src]));
            }
            Deco(h, kBlue + 1, 20, 0.8);
            h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) - predicted lower edge [GeV/c^{2}]"); h->GetXaxis()->SetTitleSize(0.042);
            h->GetYaxis()->SetTitle(kYt);
            h->SetMinimum(-0.03 * h->GetMaximum()); h->SetMaximum(1.3 * h->GetMaximum());
            h->Draw("e");
            VLine(0, 0, 0.5 * h->GetMaximum(), kMagenta + 1, 1);
            const double up = scan::EdgeMeasured(0.5 * (kSlLo[i] + kSlHi[i])) - LowEdgeMeasured(0.5 * (kSlLo[i] + kSlHi[i]));
            if (up < 0.15) VLine(up, 0, 0.5 * h->GetMaximum(), kGreen + 2, 7);
            Txt(0.17, 0.935, Form("N - CB, no cut, M(#pi^{+}#pi^{-}) %s GeV/c^{2}", kSlName[i]), 0.042);
            Txt(0.21, 0.84, "magenta: predicted lower edge", 0.040, kMagenta + 1);
            Txt(0.21, 0.79, up < 0.15 ? "green: upper edge (slice centre)" : "upper edge: right of the range", 0.040, kGreen + 2);
            Txt(0.21, 0.74, "no rows at the lower edge:", 0.040); Txt(0.21, 0.69, "the spectrum only starts there", 0.040);
        }
        Save(c, "bump_edge_position", false);
    }

    // ---- 7.10 candidates: spectra ----
    {
        TCanvas* c = new TCanvas("c_cand", "", 2700, 1650);
        c->Divide(3, 2, 0.002, 0.002);
        for (int k = 0; k < 6; ++k) {
            c->cd(k + 1); gPad->SetLeftMargin(0.19);
            const Spectra& sp = candSp[k];
            TH1D* hn = sp.Hist("us", Form("cd_us%d", k)); TH1D* hc = sp.Hist("cb", Form("cd_cb%d", k)); TH1D* hs = sp.Hist("sig", Form("cd_sig%d", k));
            Deco(hn, kBlack, 20); Deco(hc, kRed + 1, 20); Deco(hs, kBlue + 1, 20);
            hn->GetXaxis()->SetTitle(kXt); hn->GetYaxis()->SetTitle(kYt); hn->GetYaxis()->SetTitleOffset(1.75);
            hn->GetXaxis()->SetRangeUser(0.4001, 0.5999);
            double mx = 0;
            for (int b = hn->FindBin(0.401); b <= hn->FindBin(0.599); ++b) mx = std::max(mx, hn->GetBinContent(b));
            hn->SetMinimum(0); hn->SetMaximum(2.0 * mx);
            hn->Draw("e"); hc->Draw("e same"); hs->Draw("e same");
            VLine(scan::kBumpLo, 0, 0.75 * mx, kGreen + 2, 7); VLine(scan::kBumpHi, 0, 0.75 * mx, kGreen + 2, 7);
            const Res& R = cr[k];
            Txt(0.19, 0.935, Form("%s: %s", cand[k].name.c_str(), k == 0 ? "no cut" : Describe(cand[k]).c_str()), 0.038);
            if (k < 4) {
                Txt(0.22, 0.84, Form("step-height sum E = %.0f #pm %.0f (not a yield)", R.E, R.eE), 0.038);
                Txt(0.22, 0.785, Form("N - CB = %.0f,  CB/N = %.1f %% (no k factor)", R.win.sig(), 100 * R.win.cbShare()), 0.038);
            } else {
                Txt(0.22, 0.84, "E not valid (#eta-rest cut distorts the edge fits)", 0.036, kRed + 1);
                Txt(0.22, 0.785, Form("N - CB = %.0f,  CB/N = %.1f %% (no k factor)", R.win.sig(), 100 * R.win.cbShare()), 0.038);
            }
            Txt(0.22, 0.73, Form("tagged #eta kept: %.2f #pm %.2f,  gain %.2f #pm %.2f", R.effT, R.eEffT, R.gainT, R.eGainT), 0.038);
            if (cand[k].mpHi > 0) {
                Txt(0.22, 0.675, Form("CAUTION: the cut M(#pi#pi) < %.2f bends the spectrum", cand[k].mpHi), 0.036, kRed + 1);
                Txt(0.22, 0.625, "down above 0.52: the hump is NOT the bump shape", 0.036, kRed + 1);
            }
            if (k == 0) {
                TLegend* lg = new TLegend(0.22, 0.46, 0.52, 0.62); lg->SetTextSize(0.036); lg->SetFillStyle(0);
                lg->AddEntry(hn, "all pairs (N)", "pe"); lg->AddEntry(hc, "CB", "pe"); lg->AddEntry(hs, "N - CB", "pe"); lg->Draw();
                Txt(0.22, 0.675, "green: bump window;  gain = tagged #eta kept / #sqrt{N kept}", 0.034, kGreen + 2);
            }
        }
        Save(c, "bump_candidates", true);
    }

    // ---- 7.11 isBest == 1 against all combinations (nominal selection) ----
    {
        TCanvas* c = new TCanvas("c_best", "", 2800, 800);
        c->Divide(4, 1, 0.002, 0.002);
        for (int i = 0; i < kNSig; ++i) {
            c->cd(i + 1); gPad->SetLeftMargin(0.20);
            TH1D* ha = nomSl[i].Hist("sig", Form("bs_a%d", i)); TH1D* hb = nomSlBest[i].Hist("sig", Form("bs_b%d", i));
            Deco(ha, kBlue + 1, 20, 0.9); Deco(hb, kBlue + 1, 24, 0.9);
            const StepFit& fa = nomFit[i]; const StepFit& fb = nomFitBest[i];
            ha->GetXaxis()->SetTitle(kXt); ha->GetYaxis()->SetTitle(kYt); ha->GetYaxis()->SetTitleOffset(1.9);
            ha->GetXaxis()->SetRangeUser(fa.lo - 0.015 + 1e-4, fa.hi + 0.015 - 1e-4);
            double mx = 0;
            for (int b = 1; b <= ha->GetNbinsX(); ++b) if (ha->GetBinCenter(b) > fa.lo - 0.015 && ha->GetBinCenter(b) < fa.hi + 0.015) mx = std::max(mx, ha->GetBinContent(b));
            ha->SetMinimum(0); ha->SetMaximum(1.5 * mx);
            ha->Draw("e"); hb->Draw("e same");
            for (int v = 0; v < 2; ++v) {
                const StepFit& f = v ? fb : fa;
                TGraph* g = new TGraph();
                for (size_t q = 0; q < f.tx.size(); ++q) g->SetPoint((int)q, f.tx[q], f.fy[q]);
                g->SetLineWidth(2); g->SetLineStyle(v ? 2 : 1); g->Draw("l same");
            }
            VLine(fa.edge, 0, 0.6 * mx, kGreen + 2, 7);
            Txt(0.20, 0.935, Form("N - CB, nominal, M(#pi^{+}#pi^{-}) %s GeV/c^{2}", kSlName[i]), 0.040);
            Txt(0.24, 0.84, Form("all combinations (full): h = %.0f #pm %.0f", fa.h, fa.eh), 0.040);
            Txt(0.24, 0.78, Form("isBest == 1 (open, dashed): h = %.0f #pm %.0f", fb.h, fb.eh), 0.040);
            Txt(0.24, 0.72, "h = step at the green edge, rows / 2.5 MeV", 0.034);
        }
        Save(c, "bump_isbest_slices", false);
    }

    printf("total run time %.1f s\n", swAll.RealTime());
}
