// ============================================================================
// scan_ribbons.C  --  ribbon cutflow, N-1 companion, final spectra, summary table
// ============================================================================
// HADES p+p, T = 4.53 GeV, channel pi+ pi- e+ e-. Summary figures of the eta
// cut scan in M4 = M(pi+pi-e+e-), built on plots/ScanEngine.h and on the
// results of three studies (plots/scan_peak.C, scan_bump.C, scan_extra.C).
// All rows = stored particle combinations (isBest >= 0), corrected momenta,
// cor_ok == 1.
//
// Two targets, two selections:
//   PEAK  eta -> pi+pi-e+e- (whole peak content at the eta mass)
//         nominal of scan_peak.C : oa_lab < 55, M(pipi) > 0.30
//         + PROPOSED vertex cuts : vz > -160 mm, both pions used in the vertex fit
//         optional, NOT applied  : no pi0 photon (0.110 <= M(eeg) < 0.170); it has
//                                  no control sample (shown as a shaded column)
//   BUMP  eta -> pi+pi-pi0 with an e+e- pair from the pi0 (Dalitz decay or
//         converted photon; not told apart), photon not used
//         nominal of scan_bump.C : oa_lab < 40, 0.30 < M(pipi) < 0.41
//         + PROPOSED vertex cuts : vz > -160 mm, both pions used in the vertex fit
//
// Corrections after the review:
//   - the gain in S/sigma_S comes from the two vertex cuts, not from the kinematic
//     cuts. vz > -160 tightens the binding preselection (z > -500) and needs an
//     explicit decision; the efficiency of isUsedVertex must come from simulation.
//   - S/sigma_S (F_stat) is statistical only with a frozen peak shape: it compares
//     selections, it is not a significance. Yields are S +- stat +- model.
//   - CB = 2 sqrt(N++ N--) has no k factor: CB shares and N - CB are estimates.
//   - bump gain = (kept fraction of E or T) / sqrt(rows kept in the bump window):
//     one definition in scan_bump.C, scan_extra.C and here.
//
// Measures
//   peak window (engine)  N, CB, N - CB, B = fitted polynomial under the peak,
//                         S = peak content (P_main, whole Gaussian), model
//                         spread (48 variations), S/B, F_stat = S/sigma_S
//   bump window (engine)  N, CB, W = N - CB (W is an estimate: no k factor)
//   E   edge-step sum of scan_bump.C: step heights of N - CB at the predicted
//       kinematic edge in the M(pipi) slices 0.30-0.34, 0.34-0.37, 0.37-0.39,
//       0.39-0.40 (rows per 2.5 MeV). NOT a yield.
//   T   ECAL-tagged eta of scan_bump.C: eta peak in m5c = M(pi+pi-e+e-gamma)
//       with the photon energy rescaled to M(e+e-gamma) = m_pi0, pi0 window
//       0.119-0.173. Efficiency of a selection = T / T(preselection).
//
// Ribbon cutflow: x = successive cuts; line = cuts 1..j at the nominal values;
// band at step j = range of the quantity when cut j alone is moved by one
// a-priori step up and down. A cut that the selection does not use keeps its
// column ("not applied"): the line stays flat and a hatched box shows what the
// cut would do at its customary value (and one step around it).
//
// Output
//   tables   <prod>/scan/ribbons_selections.txt, ribbons_peak.txt,
//            ribbons_bump.txt, ribbons_nm1.txt, ribbons_table.txt
//   spectra  <prod>/scan/ribbons_spectra.root
//   figures  plots/output/eta_scan/ribbons_cutflow_{peak,bump}, ribbons_nm1_{peak,bump},
//            ribbons_spectra_{peak,bump}, ribbons_summary (.png and .pdf)
//
// Run (ACLiC, build directory outside the repository), about 1 min:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/aclic_ribbons", kTRUE)' 'plots/scan_ribbons.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "ScanEngine.h"

namespace rib {

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
// A. Cuts and working points
// ============================================================================
const double kNo = -9999.;      // "no cut"
inline bool On(double t) { return t > kNo + 1.; }

enum CutId { cOaLab = 0, cMpHi, cMpLo, cOaLo, cOaHi, cVz, cUsed, cPi0, cEtaRest, cNCut };
// pi0-photon veto of scan_extra.C: the stored photon (M(e+e-gamma) closest to 0.140)
// lies in kVetoLo..kVetoHi. The threshold of the cut cPi0 is the widening d of the
// window on both sides [GeV]: veto = [kVetoLo - d, kVetoHi + d).
const float kVetoLo = 0.110f, kVetoHi = 0.170f;

inline bool PassCut(int c, double t, const Row& r) {
    const float x = (float)t;
    switch (c) {
        case cOaLab:   return r.oa_lab < x;     // angle (pi+pi-),(e+e-) in the lab [deg]
        case cMpHi:    return r.mpipi < x;      // M(pi+pi-) [GeV/c2]
        case cMpLo:    return r.mpipi > x;
        case cOaLo:    return r.oa > x;         // e+e- opening angle [deg]
        case cOaHi:    return r.oa < x;
        case cVz:      return r.vz > x;         // event vertex z [mm]
        case cUsed: {                           // tracks used in the vertex fit: 1 = at least one pion,
            const int np = (r.pip_isUsedVertex != 0) + (r.pim_isUsedVertex != 0);   // 2 = both pions, 3 = all four
            if (t < 1.5) return np >= 1;
            if (t < 2.5) return np == 2;
            return np == 2 && r.ep_isUsedVertex != 0 && r.em_isUsedVertex != 0;
        }
        case cPi0:     return !(r.m_eeg >= kVetoLo - x && r.m_eeg < kVetoHi + x);
        case cEtaRest: return r.oa_eta > x;     // flag of the two reference chains only
    }
    return true;
}

struct WP {
    std::string name;
    double v[cNCut];
    int best = 0;       // 1: isBest == 1 only (explicit comparison)
    int fold = -1;      // 0 / 1: one half of the event groups
    WP() { for (int c = 0; c < cNCut; ++c) v[c] = kNo; }
};
inline bool Pass(const WP& w, const Row& r) {
    if (w.best && r.isBest != 1) return false;
    if (w.fold >= 0 && r.fold != w.fold) return false;
    for (int c = 0; c < cNCut; ++c)
        if (On(w.v[c]) && !PassCut(c, w.v[c], r)) return false;
    return true;
}
inline WP With(WP w, int c, double t, const char* name = nullptr) { w.v[c] = t; if (name) w.name = name; return w; }

// plain text of one cut (tables)
std::string CutText(int c, double t) {
    switch (c) {
        case cOaLab:   return Form("oa_lab<%.0f", t);
        case cMpHi:    return Form("Mpipi<%.2f", t);
        case cMpLo:    return Form("Mpipi>%.2f", t);
        case cOaLo:    return Form("OA>%.1f", t);
        case cOaHi:    return Form("OA<%.0f", t);
        case cVz:      return Form("vz>%.0f", t);
        case cUsed:    return t < 1.5 ? "vtxfit>=1pi" : (t < 2.5 ? "vtxfit 2pi" : "vtxfit all4");
        case cPi0:     return Form("no pi0 gamma %.3f-%.3f", kVetoLo - t, kVetoHi + t);
        case cEtaRest: return Form("eta_rest>%.0f", t);
    }
    return "?";
}
// TLatex text of one cut (figures)
std::string CutTex(int c, double t) {
    switch (c) {
        case cOaLab:   return Form("oa_{lab}(#pi#pi,ee) < %.0f#circ", t);
        case cMpHi:    return Form("M(#pi#pi) < %.2f", t);
        case cMpLo:    return Form("M(#pi#pi) > %.2f", t);
        case cOaLo:    return Form("OA(ee) > %.4g#circ", t);
        case cOaHi:    return Form("OA(ee) < %.0f#circ", t);
        case cVz:      return Form("z_{vtx} > %.0f mm", t);
        case cUsed:    return t < 1.5 ? "#geq 1 #pi in vertex fit" : (t < 2.5 ? "both #pi in vertex fit" : "all 4 tracks in vertex fit");
        case cPi0:     return "no #pi^{0} photon";
        case cEtaRest: return Form("#eta-rest angle > %.0f#circ", t);
    }
    return "?";
}
std::string Describe(const WP& w) {
    std::string s;
    for (int c = 0; c < cNCut; ++c) {
        if (!On(w.v[c])) continue;
        if (!s.empty()) s += ", ";
        s += CutText(c, w.v[c]);
    }
    if (w.best) s += s.empty() ? "isBest==1" : ", isBest==1";
    if (w.fold >= 0) s += Form("%sfold %d", s.empty() ? "" : ", ", w.fold);
    if (s.empty()) s = "no cut";
    return s;
}
// TLatex description in at most two lines
void DescribeTex(const WP& w, std::string& l1, std::string& l2) {
    std::vector<std::string> p;
    for (int c = 0; c < cNCut; ++c) if (On(w.v[c])) p.push_back(CutTex(c, w.v[c]));
    l1.clear(); l2.clear();
    if (p.empty()) { l1 = "no cut beyond the preselection"; return; }
    const size_t half = p.size() > 3 ? (p.size() + 1) / 2 : p.size();
    for (size_t i = 0; i < p.size(); ++i) {
        std::string& l = i < half ? l1 : l2;
        if (!l.empty()) l += ", ";
        l += p[i];
    }
}
std::string Key(const WP& w) {
    std::string k;
    for (int c = 0; c < cNCut; ++c) k += Form("%.4f|", w.v[c]);
    k += Form("b%d|f%d", w.best, w.fold);
    return k;
}

// ---- the selections ----------------------------------------------------------
inline WP None()     { WP w; w.name = "none"; return w; }
inline WP Slides()   { WP w; w.name = "slides"; w.v[cOaLab] = 50; w.v[cMpHi] = 0.46; w.v[cEtaRest] = 140; w.v[cOaLo] = 2; return w; }
inline WP OldCode()  { WP w; w.name = "old code"; w.v[cOaLab] = 50; w.v[cMpHi] = 0.42; w.v[cEtaRest] = 140; return w; }
// nominal selections of the two studies (kinematic cuts only)
inline WP PeakKin()  { WP w; w.name = "peak nominal (kin.)"; w.v[cOaLab] = 55; w.v[cMpLo] = 0.30; return w; }
inline WP BumpKin()  { WP w; w.name = "bump nominal (kin.)"; w.v[cOaLab] = 40; w.v[cMpLo] = 0.30; w.v[cMpHi] = 0.41; return w; }
// the same with the two PROPOSED vertex cuts of scan_extra.C (vz > -160 needs a decision: the binding
// preselection is z > -500). The pi0-photon veto is optional and unverified: it is NOT part of the selection.
const double kVzCut = -160.;
inline WP PeakFull() { WP w = PeakKin(); w.name = "peak nominal + vertex"; w.v[cVz] = kVzCut; w.v[cUsed] = 2; return w; }
inline WP PeakVeto() { WP w = PeakFull(); w.name = "peak nom. + vtx + veto"; w.v[cPi0] = 0; return w; }
inline WP BumpFull() { WP w = BumpKin(); w.name = "bump nominal + vertex"; w.v[cVz] = kVzCut; w.v[cUsed] = 2; return w; }

// ============================================================================
// B. Bump measures (definitions copied from plots/scan_bump.C)
// ============================================================================
const int kNSl = 6;        // 0..3 signal slices (their sum is E), 4 = 0.40-0.41 (not in E), 5 = control
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
inline void PassSlice(const WP& w, int i, double& lo, double& hi) {
    lo = kSlLo[i]; hi = kSlHi[i];
    if (On(w.v[cMpLo])) lo = std::max(lo, w.v[cMpLo]);
    if (On(w.v[cMpHi])) hi = std::min(hi, w.v[cMpHi]);
}
// tagged eta
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730, kMeeMax = 0.100;
const double kTagMean = 0.5463, kTagSigma = 0.0062, kTagLo = 0.49, kTagHi = 0.62;
const int    kTagOrder = 2, kTagNbin = 160, kTagReb = 8;
const double kTagHLo = 0.40, kTagHHi = 0.80;
inline double M5c(const Row& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
inline bool Tagged(const Row& r) {
    return r.g_idx > 0 && r.mee < kMeeMax && r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi;
}
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
// efficiency of a subsample p of b (b = p + independent rest)
inline void SubRatio(double p, double ep, double b, double eb, double& eff, double& err) {
    if (!(b > 0)) { eff = 0; err = 0; return; }
    eff = p / b;
    err = std::sqrt(std::max(0.0, (1.0 - 2.0 * eff) * ep * ep + eff * eff * eb * eb)) / b;
}
const int kNEvar = 28;     // variations of the edge fit for ratios (as CandSpread of scan_bump.C)
const int kNTvar = 18;     // variations of the tagged-peak fit

// ============================================================================
// C. Measurement of a list of working points (cached)
// ============================================================================
struct Res {
    WP w;
    // peak window
    PeakFit pm, pp, l1, l2;          // P_main, P_pois, null fits on e+e+ and e-e-
    scan::Spread spr;                // model spread of S
    scan::FreeFit fr;                // position and width free (check of the frozen shape)
    // bump window
    double bN = 0, bCB = 0, W = 0, eW = 0;
    // edge step
    bool okf[kNSl];
    StepFit sf[kNSl];
    double E = 0, eE = 0, sprE = 0, Echi2 = 0; int Endf = 0;
    double Evar[kNEvar]; bool EvarOk[kNEvar];
    // tagged eta
    PeakFit tg;
    double tagRows = 0;
    double Tvar[kNTvar]; bool TvarOk[kNTvar];
    double sprT = 0;

    double Sstat() const { return pm.eS; }
    double Stot() const { return std::sqrt(pm.eS * pm.eS + spr.half * spr.half); }
    double D() const { return pm.N - pm.CB; }
    double Q() const { return bN > 0 ? E / std::sqrt(bN) : 0; }       // E / sqrt(all rows in the bump window)
    double eQ() const { return bN > 0 ? eE / std::sqrt(bN) : 0; }
    double Eprob() const { return Endf > 0 ? TMath::Prob(Echi2, Endf) : 0; }
    bool Evalid() const { return Eprob() > 1e-3; }
};
std::map<std::string, Res> gRes;
const Res& R(const WP& w) {
    const auto it = gRes.find(Key(w));
    if (it == gRes.end()) { printf("R(): selection '%s' was not measured\n", Describe(w).c_str()); static Res dummy; return dummy; }
    return it->second;
}

// d: all rows; dT: tagged rows only (rib::Tagged)
void Measure(const Data& d, const Data& dT, const std::vector<WP>& list) {
    std::vector<WP> todo;
    for (const WP& w : list) {
        const std::string k = Key(w);
        if (gRes.count(k)) continue;
        bool dup = false;
        for (const WP& t : todo) if (Key(t) == k) { dup = true; break; }
        if (!dup) todo.push_back(w);
    }
    const size_t batch = 16;
    const int per = 1 + kNSl;
    for (size_t j0 = 0; j0 < todo.size(); j0 += batch) {
        const size_t j1 = std::min(todo.size(), j0 + batch);
        std::vector<Sel> sels, ts;
        for (size_t j = j0; j < j1; ++j) {
            const WP w = todo[j];
            sels.push_back([w](const Row& r) { return Pass(w, r); });
            for (int i = 0; i < kNSl; ++i) {
                const float lo = kSlLo[i], hi = kSlHi[i];
                sels.push_back([w, lo, hi](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && Pass(w, r); });
            }
            ts.push_back([w](const Row& r) { return Pass(w, r); });
        }
        const std::vector<Spectra> sp = scan::FillMany(d, sels);
        const std::vector<Spectra> st = scan::FillMany(dT, ts, TagOpt());
        for (size_t j = j0; j < j1; ++j) {
            Res X;
            X.w = todo[j];
            const WP& w = X.w;
            const Spectra* s = &sp[(j - j0) * per];
            // ---- peak ----
            X.pm = scan::FitPeakMain(s[0]);
            X.pp = scan::FitPeakPois(s[0]);
            X.l1 = scan::FitPeakLikeSign(s[0], 1);
            X.l2 = scan::FitPeakLikeSign(s[0], 2);
            X.spr = scan::PeakModelSpread(s[0]);
            X.fr = scan::FitPeakFree(s[0]);
            // ---- bump window ----
            const WinSum& wb = s[0].win[1];
            X.bN = wb.n[0]; X.bCB = wb.cb; X.W = wb.sig(); X.eW = wb.sigErr();
            // ---- edge step ----
            double v = 0;
            for (int i = 0; i < kNSl; ++i) {
                X.okf[i] = false;
                double lo, hi; PassSlice(w, i, lo, hi);
                if (hi - lo < 0.0045 || s[1 + i].nsel[0] < 50) continue;
                X.sf[i] = scan::FitEdgeStep(s[1 + i], SliceCfg(i, lo, hi));
                X.okf[i] = X.sf[i].ok;
                if (!X.okf[i] || i >= kNSig) continue;
                X.E += X.sf[i].h; v += X.sf[i].eh * X.sf[i].eh;
                X.Echi2 += X.sf[i].chi2; X.Endf += X.sf[i].ndf;
            }
            X.eE = std::sqrt(v);
            // model spread of E: order 2, 3; half range 35, 45, 55 MeV; smearing 4.4, 5.4, 6.4 MeV;
            // per slice half of (max - min) over the fits with p > 0.05; slices added linearly
            for (int i = 0; i < kNSig; ++i) {
                if (!X.okf[i]) continue;
                double lo, hi; PassSlice(w, i, lo, hi);
                double mn = X.sf[i].h, mx = X.sf[i].h;
                for (int ord = 2; ord <= 3; ++ord)
                    for (int ih = 0; ih < 3; ++ih)
                        for (int is = -1; is <= 1; ++is) {
                            StepCfg c = SliceCfg(i, lo, hi);
                            c.order = ord; c.halfRange = 0.035 + 0.010 * ih; c.s = scan::kEdgeSigma + 0.001 * is;
                            const StepFit f = scan::FitEdgeStep(s[1 + i], c);
                            if (!f.ok || !(TMath::Prob(f.chi2, f.ndf) > 0.05)) continue;
                            mn = std::min(mn, f.h); mx = std::max(mx, f.h);
                        }
                X.sprE += 0.5 * (mx - mn);
            }
            // variations for ratios of E between two selections
            for (int q = 0; q < kNEvar; ++q) {
                const int ih = q < 27 ? q % 3 - 1 : 0, is = q < 27 ? (q / 3) % 3 - 1 : 0, ic = q < 27 ? q / 9 - 1 : 0;
                X.Evar[q] = 0; X.EvarOk[q] = true;
                for (int i = 0; i < kNSig; ++i) {
                    double lo, hi; PassSlice(w, i, lo, hi);
                    if (hi - lo < 0.0045 || s[1 + i].nsel[0] < 50) continue;
                    StepCfg c = SliceCfg(i, lo, hi);
                    if (q == 27) { c.order = 3; if (i < 3) c.halfRange = 0.055; }
                    else { if (i < 3) c.halfRange += 0.010 * ih; c.s += 0.001 * is; c.scale += 0.0009 * ic; }
                    const StepFit f = scan::FitEdgeStep(s[1 + i], c);
                    if (!f.ok) { X.EvarOk[q] = false; break; }
                    X.Evar[q] += f.h;
                }
            }
            // ---- tagged eta ----
            const Spectra& t = st[j - j0];
            X.tagRows = t.nsel[0];
            X.tg = scan::FitPeakMain(t, TagCfg());
            {
                const double rl[3] = {0.49, 0.48, 0.50}, rh[3] = {0.62, 0.66, 0.60};
                int q = 0; double mn = 0, mx = 0; bool first = true;
                for (int ord = 2; ord <= 3; ++ord)
                    for (int ir = 0; ir < 3; ++ir)
                        for (int iw = -1; iw <= 1; ++iw, ++q) {
                            PeakCfg c = TagCfg(); c.order = ord; c.lo = rl[ir]; c.hi = rh[ir]; c.sigma = kTagSigma + 0.0003 * iw;
                            const PeakFit a = scan::FitPeakMain(t, c);
                            X.TvarOk[q] = a.ok; X.Tvar[q] = a.S;
                            if (!a.ok) continue;
                            if (first) { mn = mx = a.S; first = false; }
                            mn = std::min(mn, a.S); mx = std::max(mx, a.S);
                        }
                X.sprT = 0.5 * (mx - mn);
            }
            gRes[Key(w)] = X;
        }
    }
}

// efficiency of the tagged eta and of the edge step relative to a base selection
struct Eff { double v = 0, e = 0, spr = 0; };
Eff EffT(const Res& a, const Res& b) {
    Eff x;
    if (!a.tg.ok || !b.tg.ok || !(b.tg.S > 0)) return x;
    SubRatio(a.tg.S, a.tg.eS, b.tg.S, b.tg.eS, x.v, x.e);
    double mn = 0, mx = 0; bool first = true;
    for (int q = 0; q < kNTvar; ++q) {
        if (!a.TvarOk[q] || !b.TvarOk[q] || !(b.Tvar[q] > 0)) continue;
        const double r = a.Tvar[q] / b.Tvar[q];
        if (first) { mn = mx = r; first = false; }
        mn = std::min(mn, r); mx = std::max(mx, r);
    }
    x.spr = 0.5 * (mx - mn);
    return x;
}
Eff EffE(const Res& a, const Res& b) {
    Eff x;
    if (!(b.E > 0)) return x;
    SubRatio(a.E, a.eE, b.E, b.eE, x.v, x.e);
    double mn = 0, mx = 0; bool first = true;
    for (int q = 0; q < kNEvar; ++q) {
        if (!a.EvarOk[q] || !b.EvarOk[q] || !(b.Evar[q] > 0)) continue;
        const double r = a.Evar[q] / b.Evar[q];
        if (first) { mn = mx = r; first = false; }
        mn = std::min(mn, r); mx = std::max(mx, r);
    }
    x.spr = 0.5 * (mx - mn);
    return x;
}

// ============================================================================
// D. Cut flow
// ============================================================================
struct Step {
    int cut = -1;                  // -1 = preselection
    bool applied = true;
    double nom = kNo;              // nominal threshold; for a cut that is not applied: its customary value
    std::vector<double> alt;       // applied: one step down and up; not applied: customary value and one step around it
    std::string lab, band, note;   // TLatex: cut, values of the band, optional third line
    std::string why;               // plain remark for the table
};
struct Flow {
    std::string tag, target;
    std::vector<Step> st;                  // st[0] = preselection
    std::vector<WP> cum;                   // cuts 1..j at their nominal values
    std::vector<std::vector<WP>> alt;      // step j: cut j alone moved
    std::vector<WP> nm1;                   // nm1[0] = full selection; j: cut j removed (applied) or added (not applied)
    int n() const { return (int)st.size(); }
    const WP& full() const { return cum.back(); }
};
Step MkStep(int cut, bool applied, double nom, std::vector<double> alt, const char* band, const char* why, const char* note = "") {
    Step s; s.cut = cut; s.applied = applied; s.nom = nom; s.alt = alt; s.band = band; s.why = why; s.note = note;
    s.lab = CutTex(cut, nom);
    return s;
}
void Build(Flow& f) {
    f.cum.clear(); f.alt.clear(); f.nm1.clear();
    WP w; w.name = "preselection";
    f.cum.push_back(w);
    f.alt.push_back(std::vector<WP>());
    for (int j = 1; j < f.n(); ++j) {
        const Step& s = f.st[j];
        std::vector<WP> a;
        for (double t : s.alt) a.push_back(With(w, s.cut, t));
        if (s.applied) w = With(w, s.cut, s.nom);
        f.cum.push_back(w);
        f.alt.push_back(a);
    }
    f.nm1.push_back(w);
    for (int j = 1; j < f.n(); ++j) {
        const Step& s = f.st[j];
        f.nm1.push_back(With(w, s.cut, s.applied ? kNo : s.nom));
    }
}
Flow PeakFlow() {
    Flow f; f.tag = "peak"; f.target = "eta peak";
    Step p; p.lab = "preselection"; p.band = "isBest #geq 0, vertex, PT3"; p.why = "all stored combinations, eVertReco_z > -500, PT3, cor_ok == 1";
    f.st.push_back(p);
    f.st.push_back(MkStep(cOaLab, true, 55, {50, 60}, "band: 50#circ, 60#circ", "nominal of scan_peak.C (plateau 35 - none)"));
    f.st.push_back(MkStep(cMpHi, false, 0.46, {0.45, 0.46, 0.47}, "if applied: 0.45, 0.46, 0.47", "not applied: every upper M(pipi) cut is flagged as sculpting by scan_peak.C; 0.46 = slides / start point"));
    f.st.push_back(MkStep(cMpLo, true, 0.30, {0.28, 0.32}, "band: 0.28, 0.32", "nominal of scan_peak.C (plateau none - 0.34)"));
    f.st.push_back(MkStep(cOaLo, false, 2, {1.5, 2, 2.5}, "if applied: 1.5, 2, 2.5#circ", "not applied: only loses (scan_peak.C); 2 deg = slides"));
    f.st.push_back(MkStep(cOaHi, false, 20, {18, 20, 22}, "if applied: 18, 20, 22#circ", "not applied: no gain in F (scan_peak.C); 20 deg removes most of the like-sign estimate of the CB"));
    f.st.push_back(MkStep(cVz, true, kVzCut, {-170, -150}, "band: -170, -150 mm", "PROPOSED, needs a decision (binding cut: z > -500). scan_extra.C: tagged efficiency 0.998 +- 0.007, peak G = 1.06 +- 0.03; no signal below -160 mm"));
    f.st.push_back(MkStep(cUsed, true, 2, {1, 3}, "band: #geq1 #pi; all four", "PROPOSED. scan_extra.C: tagged efficiency 0.950 +- 0.013 for 0.67 of the background (G_tag = 1.16), peak G = 1.12 +- 0.05; reconstruction flag: efficiency from simulation; alternative |z_pi - vz| < 40 mm"));
    f.st.push_back(MkStep(cPi0, false, 0, {-0.010, 0, 0.010}, "if applied: edges #pm10 MeV", "OPTIONAL, not applied: no control sample; it raises S/B, a gain in S/sigma_S is not shown (alone G = 1.04 +- 0.03)",
                          "M(ee#gamma): 0.110 - 0.170"));
    Build(f);
    return f;
}
Flow BumpFlow() {
    Flow f; f.tag = "bump"; f.target = "bump";
    Step p; p.lab = "preselection"; p.band = "isBest #geq 0, vertex, PT3"; p.why = "all stored combinations, eVertReco_z > -500, PT3, cor_ok == 1";
    f.st.push_back(p);
    f.st.push_back(MkStep(cOaLab, true, 40, {35, 45}, "band: 35#circ, 45#circ", "nominal of scan_bump.C (plateau 30 - 45)"));
    f.st.push_back(MkStep(cMpHi, true, 0.41, {0.40, 0.42}, "band: 0.40, 0.42", "nominal of scan_bump.C; kinematic limit m_eta - m_pi0 = 0.413"));
    f.st.push_back(MkStep(cMpLo, true, 0.30, {0.28, 0.32}, "band: 0.28, 0.32", "nominal of scan_bump.C"));
    f.st.push_back(MkStep(cOaLo, false, 2, {1.5, 2, 2.5}, "if applied: 1.5, 2, 2.5#circ", "not applied: every value loses (scan_bump.C); 2 deg = slides"));
    f.st.push_back(MkStep(cOaHi, false, 20, {18, 20, 22}, "if applied: 18, 20, 22#circ", "not applied: gain 1.02-1.04 from the tagged eta with the sqrt(N) definition (scan_bump.C), small; it removes rows with a large like-sign share"));
    f.st.push_back(MkStep(cVz, true, kVzCut, {-170, -150}, "band: -170, -150 mm", "PROPOSED, needs a decision (binding cut: z > -500). scan_extra.C: tagged efficiency 0.998 +- 0.007, Gb = 1.052 +- 0.007; common value with the peak"));
    f.st.push_back(MkStep(cUsed, true, 2, {1, 3}, "band: #geq1 #pi; all four", "PROPOSED. scan_extra.C: Gb = 1.135 +- 0.015 from the tagged eta, 1.07 +- 0.04 from the edge step; efficiency from simulation"));
    Build(f);
    return f;
}
void Collect(const Flow& f, std::vector<WP>& v) {
    for (const WP& w : f.cum) v.push_back(w);
    for (const auto& a : f.alt) for (const WP& w : a) v.push_back(w);
    for (const WP& w : f.nm1) v.push_back(w);
}

// one quantity along the flow: central value, band of the cut variation, error of the central value
typedef std::function<double(const Res&)> Qf;
struct Ser { std::vector<double> c, lo, hi, e; };
Ser Series(const Flow& f, const Qf& q, const Qf& err = Qf()) {
    Ser s;
    for (int j = 0; j < f.n(); ++j) {
        const double c = q(R(f.cum[j]));
        double lo = c, hi = c;
        for (const WP& a : f.alt[j]) { const double v = q(R(a)); lo = std::min(lo, v); hi = std::max(hi, v); }
        s.c.push_back(c); s.lo.push_back(lo); s.hi.push_back(hi);
        s.e.push_back(err ? err(R(f.cum[j])) : 0.);
    }
    return s;
}
// the same for the N-1 list (no band)
Ser SeriesNm1(const Flow& f, const Qf& q, const Qf& err = Qf()) {
    Ser s;
    for (int j = 0; j < f.n(); ++j) {
        const double c = q(R(f.nm1[j]));
        s.c.push_back(c); s.lo.push_back(c); s.hi.push_back(c);
        s.e.push_back(err ? err(R(f.nm1[j])) : 0.);
    }
    return s;
}

// ============================================================================
// E. Printing
// ============================================================================
FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    char buf[4096];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof(buf), fmt, ap); va_end(ap);
    fputs(buf, stdout);
    if (gOut) fputs(buf, gOut);
}
void Open(const char* name) {
    gOut = fopen(Form("%s/%s", scan::kScanDir, name), "w");
    if (!gOut) printf("cannot write %s/%s\n", scan::kScanDir, name);
    P("# %s  --  written by plots/scan_ribbons.C\n", name);
    P("# data: %s ; rows = stored particle combinations (isBest >= 0), corrected momenta, cor_ok == 1\n", scan::kCompact);
}
void Close() { if (gOut) fclose(gOut); gOut = nullptr; }

// ============================================================================
// F. Drawing
// ============================================================================
const char* const kXt = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
const char* const kYt = "rows / 2.5 MeV/c^{2}";
const char* const kHdr = "p+p 4.53 GeV, all combinations (isBest #geq 0), corrected momenta";
const int kColN = kBlack, kColCB = kRed + 1, kColS = kBlue + 1, kColB = kGray + 2, kColT = kGreen + 2, kColE = kBlue + 1;
int ShadeCol() { static const int c = TColor::GetColor(255, 240, 214); return c; }

void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetHatchesLineWidth(1);
    gStyle->SetHatchesSpacing(0.7);
    gStyle->SetLineScalePS(1.5);
}
void Save(TCanvas* c, const char* name) {
    c->SaveAs(Form("%s/%s.png", scan::kFigDir, name));
    c->SaveAs(Form("%s/%s.pdf", scan::kFigDir, name));
}
// text in pixels (font 43), NDC coordinates of the current pad
void Txt(double x, double y, const char* t, double px = 24, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(43); l.SetTextSize(px); l.SetTextColor(col); l.SetTextAlign(align);
    l.DrawLatex(x, y, t);
}
// text in user coordinates
void TxtU(double x, double y, const char* t, double px = 20, int col = kBlack, int align = 21) {
    TLatex l; l.SetTextFont(43); l.SetTextSize(px); l.SetTextColor(col); l.SetTextAlign(align);
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
void AxPx(TAxis* a, double lab, double tit, double off) {
    a->SetLabelFont(43); a->SetLabelSize(lab); a->SetTitleFont(43); a->SetTitleSize(tit); a->SetTitleOffset(off);
}
void DecoH(TH1* h, int col, int mk, double ms = 0.7) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
}

// ---- the cut-flow canvas: header, three pads with a common x axis ----------------
const double kLm = 0.085, kRm = 0.065;
const double kHalfW = 0.30;                 // half width of a band box in x (step = 1)
struct Pads { TCanvas* c = nullptr; TPad* head = nullptr; TPad* p[3] = {nullptr, nullptr, nullptr}; };
Pads MakePads(const char* name) {
    Pads q;
    q.c = new TCanvas(name, name, 2400, 1650);
    const double y[5] = {1.0, 0.89, 0.47, 0.30, 0.0};
    q.head = new TPad(Form("%s_h", name), "", 0, y[1], 1, y[0]);
    q.head->SetMargin(0, 0, 0, 0); q.head->Draw();
    for (int i = 0; i < 3; ++i) {
        q.c->cd();
        q.p[i] = new TPad(Form("%s_p%d", name, i), "", 0, y[i + 2], 1, y[i + 1]);
        q.p[i]->SetLeftMargin(kLm); q.p[i]->SetRightMargin(kRm);
        q.p[i]->SetTopMargin(i == 0 ? 0.02 : 0.05);
        q.p[i]->SetBottomMargin(i == 2 ? 0.36 : 0.03);
        q.p[i]->Draw();
    }
    return q;
}
// frame of one pad; columns of cuts that are not applied get a light background
TH1F* Frame(const Flow& f, const char* name, double y1, double y2, const char* yt, bool logy) {
    const int n = f.n();
    TH1F* h = new TH1F(name, "", n, 0, n);
    h->SetDirectory(nullptr);
    h->SetMinimum(y1); h->SetMaximum(y2);
    h->GetXaxis()->SetLabelSize(0); h->GetXaxis()->SetTickLength(0);
    h->GetXaxis()->SetNdivisions(n, 0, 0, kFALSE);
    h->GetYaxis()->SetTitle(yt);
    AxPx(h->GetYaxis(), 30, 32, 1.75);
    h->GetYaxis()->CenterTitle(true);
    h->GetYaxis()->SetTickLength(0.006);
    h->Draw("axis");
    gPad->SetLogy(logy);
    for (int j = 1; j < n; ++j) {
        if (!f.st[j].applied) {
            TBox* b = new TBox(j, y1, j + 1, y2);
            b->SetFillColor(ShadeCol()); b->SetFillStyle(1001); b->SetLineWidth(0);
            b->Draw();
        }
    }
    for (int j = 1; j < n; ++j) VLine(j, y1, y2, kGray + 1, 3, 1);
    return h;
}
// how one quantity is drawn
struct Opt {
    int col = kBlack, mk = 20, lsty = 1, lwid = 2;
    double ms = 1.6;
    bool band = true;        // ribbon of the cut variation
    bool hatch = false;      // hatched box = error of the central value
    bool bars = false;       // error bars
    double dx = 0;           // shift in x (two quantities in one pad)
    double alpha = 0.32;
};
Opt MkOpt(int col, int mk, double ms = 1.6, bool bars = false, bool band = true, bool hatch = false, double dx = 0, int lsty = 1, int lwid = 2) {
    Opt o; o.col = col; o.mk = mk; o.ms = ms; o.bars = bars; o.band = band; o.hatch = hatch; o.dx = dx; o.lsty = lsty; o.lwid = lwid;
    return o;
}
// ribbon of one quantity: band joined from step to step, dashed box for cuts not applied, line, markers
void DrawSer(const Flow& f, const Ser& s, const Opt& o) {
    const int n = f.n();
    // error of the central value: hatched box (not repeated in the columns of cuts that are not applied)
    if (o.hatch) {
        for (int j = 0; j < n; ++j) {
            if (j > 0 && !f.st[j].applied) continue;
            TBox* b = new TBox(j + 0.5 + o.dx - kHalfW, s.c[j] - s.e[j], j + 0.5 + o.dx + kHalfW, s.c[j] + s.e[j]);
            b->SetFillColor(o.col); b->SetFillStyle(3345); b->SetLineColor(o.col); b->SetLineWidth(1);
            b->Draw();
        }
    }
    if (o.band) {
        TGraph* g = new TGraph();
        int k = 0;
        for (int j = 0; j < n; ++j) {
            const double hi = f.st[j].applied ? s.hi[j] : s.c[j];
            g->SetPoint(k++, j + 0.5 + o.dx - kHalfW, hi); g->SetPoint(k++, j + 0.5 + o.dx + kHalfW, hi);
        }
        for (int j = n - 1; j >= 0; --j) {
            const double lo = f.st[j].applied ? s.lo[j] : s.c[j];
            g->SetPoint(k++, j + 0.5 + o.dx + kHalfW, lo); g->SetPoint(k++, j + 0.5 + o.dx - kHalfW, lo);
        }
        g->SetFillColorAlpha(o.col, o.alpha); g->SetFillStyle(1001);
        g->Draw("f same");
        // cuts not applied: dashed box from the smallest to the largest value with the cut
        for (int j = 1; j < n; ++j) {
            if (f.st[j].applied) continue;
            TBox* b2 = new TBox(j + 0.5 + o.dx - kHalfW, s.lo[j], j + 0.5 + o.dx + kHalfW, s.hi[j]);
            b2->SetFillColorAlpha(o.col, 0.13); b2->SetFillStyle(1001); b2->SetLineWidth(0);
            b2->Draw();
            TBox* b = new TBox(j + 0.5 + o.dx - kHalfW, s.lo[j], j + 0.5 + o.dx + kHalfW, s.hi[j]);
            b->SetFillStyle(0); b->SetLineColor(o.col); b->SetLineStyle(2); b->SetLineWidth(2);
            b->Draw("l");
        }
    }
    TGraphErrors* gl = new TGraphErrors();
    for (int j = 0; j < n; ++j) { gl->SetPoint(j, j + 0.5 + o.dx, s.c[j]); gl->SetPointError(j, 0, o.bars ? s.e[j] : 0); }
    gl->SetLineColor(o.col); gl->SetLineStyle(o.lsty); gl->SetLineWidth(o.lwid);
    gl->SetMarkerColor(o.col); gl->SetMarkerStyle(o.mk); gl->SetMarkerSize(o.ms);
    if (o.mk > 0) {
        // the line without error bars, then the markers with them
        TGraph* l = new TGraph();
        for (int j = 0; j < n; ++j) l->SetPoint(j, j + 0.5 + o.dx, s.c[j]);
        l->SetLineColor(o.col); l->SetLineStyle(o.lsty); l->SetLineWidth(o.lwid);
        l->Draw("l same");
        gl->Draw(o.bars ? "pz same" : "px same");
    } else {
        gl->Draw("lx same");
    }
}
// points of the N-1 figure: horizontal line at the full selection, markers with errors
// (open markers for cuts that are added); stat: errors of the bars when the hatched box shows the total
void DrawNm1(const Flow& f, const Ser& s, const Opt& o, const Ser* stat = nullptr) {
    const int n = f.n();
    HLine(0, n, s.c[0], o.col, 2, 1);
    if (o.hatch) {
        for (int j = 0; j < n; ++j) {
            TBox* b = new TBox(j + 0.5 + o.dx - kHalfW, s.c[j] - s.e[j], j + 0.5 + o.dx + kHalfW, s.c[j] + s.e[j]);
            b->SetFillColor(o.col); b->SetFillStyle(3345); b->SetLineColor(o.col); b->SetLineWidth(1);
            b->Draw();
        }
    }
    if (o.mk <= 0) {
        for (int j = 0; j < n; ++j) HLine(j + 0.5 - kHalfW, j + 0.5 + kHalfW, s.c[j], o.col, o.lsty, 3);
        return;
    }
    TGraphErrors* g = new TGraphErrors(); TGraphErrors* go = new TGraphErrors();
    for (int j = 0; j < n; ++j) {
        const double e = stat ? stat->e[j] : (o.bars ? s.e[j] : 0);
        TGraphErrors* t = (j == 0 || f.st[j].applied) ? g : go;
        const int k = t->GetN();
        t->SetPoint(k, j + 0.5 + o.dx, s.c[j]); t->SetPointError(k, 0, e);
    }
    const int open = o.mk == 20 ? 24 : (o.mk == 21 ? 25 : (o.mk == 22 ? 26 : (o.mk == 33 ? 27 : o.mk)));
    g->SetLineColor(o.col); g->SetMarkerColor(o.col); g->SetMarkerStyle(o.mk); g->SetMarkerSize(o.ms); g->SetLineWidth(2);
    go->SetLineColor(o.col); go->SetMarkerColor(o.col); go->SetMarkerStyle(open); go->SetMarkerSize(o.ms); go->SetLineWidth(2);
    if (g->GetN()) g->Draw("pz same");
    if (go->GetN()) go->Draw("pz same");
}
// column labels under the lowest pad (three lines)
void StepLabels(const Flow& f, bool nm1) {
    const int n = f.n();
    const double bm = gPad->GetBottomMargin();
    for (int j = 0; j < n; ++j) {
        const Step& s = f.st[j];
        const double x = kLm + (j + 0.5) / n * (1.0 - kLm - kRm);
        std::string l1 = s.lab, l2 = s.band, l3 = s.note;
        int col = kBlack;
        if (!nm1) {
            if (j > 0 && !s.applied) { l2 = "NOT APPLIED"; l3 = s.band; col = kOrange + 4; }
        } else {
            l3.clear();
            if (j == 0) { l1 = "all cuts"; l2 = "of the selection"; }
            else if (s.applied) { l2 = "REMOVED"; }
            else { l2 = "ADDED"; col = kOrange + 4; }
        }
        Txt(x, bm - 0.035, l1.c_str(), 26, col, 23);
        Txt(x, bm - 0.130, l2.c_str(), 20, col, 23);
        if (!l3.empty()) Txt(x, bm - 0.205, l3.c_str(), 20, col, 23);
    }
}
void Range(const std::vector<Ser>& v, double& lo, double& hi, bool withErr) {
    lo = 1e30; hi = -1e30;
    for (const Ser& s : v)
        for (size_t j = 0; j < s.c.size(); ++j) {
            const double e = withErr ? s.e[j] : 0;
            lo = std::min(lo, std::min(s.lo[j], s.c[j] - e));
            hi = std::max(hi, std::max(s.hi[j], s.c[j] + e));
        }
}
// number above (up = true) or below the point of a series, clear of its band and error bar
void Num(const Ser& s, int j, const char* txt, int col, bool up, double gap, double dx = 0, bool band = true, double px = 20) {
    const double hi = band ? s.hi[j] : s.c[j], lo = band ? s.lo[j] : s.c[j];
    const double y = up ? std::max(hi, s.c[j] + s.e[j]) + gap : std::min(lo, s.c[j] - s.e[j]) - gap;
    TxtU(j + 0.5 + dx, y, txt, px, col, up ? 21 : 23);
}
TGraph* LegMark(const Opt& o, bool fill, int fsty = 1001) {
    TGraph* g = new TGraph();
    g->SetLineColor(o.col); g->SetMarkerColor(o.col); g->SetMarkerStyle(o.mk > 0 ? o.mk : 1); g->SetMarkerSize(o.ms);
    g->SetLineStyle(o.lsty); g->SetLineWidth(o.lwid);
    if (fill) { if (fsty == 1001) g->SetFillColorAlpha(o.col, o.alpha); else g->SetFillColor(o.col); g->SetFillStyle(fsty); }
    return g;
}

} // namespace rib

// ============================================================================
void scan_ribbons() {
    using namespace rib;
    TStopwatch swAll;
    Style();
    TH1::AddDirectory(kFALSE);
    gSystem->mkdir(scan::kFigDir, kTRUE);

    Data d, dT;
    if (!scan::Load(d)) return;
    for (int k = 0; k < 3; ++k) {
        for (const auto& r : d.s[k]) if (Tagged(r)) dT.s[k].push_back(r);
        dT.s[k].shrink_to_fit();
    }
    printf("tagged rows (pi0 window %.3f-%.3f, m_ee < %.2f): %zu / %zu / %zu\n", kPi0Lo, kPi0Hi, kMeeMax,
           dT.s[0].size(), dT.s[1].size(), dT.s[2].size());

    // ------------------------------------------------------------------------
    // 1. selections, flows, measurement
    // ------------------------------------------------------------------------
    const Flow fp = PeakFlow(), fb = BumpFlow();
    const WP wNone = None(), wSlides = Slides(), wOld = OldCode();
    const WP wPk = PeakKin(), wPf = PeakFull(), wBk = BumpKin(), wBf = BumpFull();
    std::vector<WP> all = {wNone, wSlides, wOld, wPk, wPf, wBk, wBf};
    Collect(fp, all); Collect(fb, all);
    // checks of the final selections: two halves of the event groups, isBest == 1
    std::vector<WP> chk;
    for (const WP& w : {wNone, wPf, wBf}) {
        for (int fo = 0; fo < 2; ++fo) { WP x = w; x.fold = fo; chk.push_back(x); }
        WP x = w; x.best = 1; chk.push_back(x);
    }
    // the optional pi0-photon veto on top of the peak selection; the two vertex cuts without kinematic cuts
    const WP wPv = PeakVeto();
    WP wVtx; wVtx.name = "vertex cuts only"; wVtx.v[cVz] = kVzCut; wVtx.v[cUsed] = 2;
    for (const WP& w : chk) all.push_back(w);
    all.push_back(wPv); all.push_back(wVtx);
    {
        TStopwatch sw;
        Measure(d, dT, all);
        printf("measured %zu selections in %.1f s\n", gRes.size(), sw.RealTime());
    }
    const Res& R0 = R(wNone);

    // ------------------------------------------------------------------------
    // 2. tables
    // ------------------------------------------------------------------------
    // ---- selections and reproduction of the study numbers ----
    Open("ribbons_selections.txt");
    P("#\n# SELECTIONS\n");
    for (const WP& w : {wNone, wSlides, wOld, wPk, wPf, wBk, wBf})
        P("   %-22s %s\n", w.name.c_str(), Describe(w).c_str());
    P("#\n# 'no pi0 gamma a-b': the stored photon (M(e+e-gamma) closest to 0.140) is NOT in a <= M(e+e-gamma) < b.\n");
    P("# 'vtxfit 2pi': both pions have isUsedVertex == 1. vz = eVertReco_z [mm].\n");
    for (const Flow* f : {&fp, &fb}) {
        P("#\n# CUT ORDER AND A-PRIORI STEPS, %s\n", f->target.c_str());
        P("   %-2s %-28s %-11s %-30s %s\n", "j", "cut", "status", "values of the band", "basis");
        for (int j = 0; j < f->n(); ++j) {
            const Step& s = f->st[j];
            std::string a;
            for (double t : s.alt) { if (!a.empty()) a += " ; "; a += s.cut == cPi0 ? Form("%.3f-%.3f", kVetoLo - t, kVetoHi + t) : CutText(s.cut, t); }
            P("   %-2d %-28s %-11s %-30s %s\n", j, j == 0 ? "preselection" : CutText(s.cut, s.nom).c_str(),
              j == 0 ? "-" : (s.applied ? "applied" : "not applied"), a.c_str(), s.why.c_str());
        }
    }
    P("#\n# REPRODUCTION OF THE STUDY NUMBERS (this macro / study report)\n");
    {
        const Res& a = R(wNone); const Res& b = R(wPk); const Res& c = R(wPf); const Res& e = R(wBk); const Res& s = R(wPv); const Res& v = R(wVtx);
        const Res& g = R(wBf);
        P("   (study numbers after the review corrections: scan_peak.C, scan_bump.C, scan_extra.C)\n");
        P("   peak, no cut           N %.0f / 86011   S %.0f +- %.0f +- %.0f / 4657 +- 390 +- 1233   F %.2f / 11.93\n", a.pm.N, a.pm.S, a.pm.eS, a.spr.half, a.pm.F);
        P("   peak nominal (kin.)    N %.0f / 68431   S %.0f +- %.0f +- %.0f / 4319 +- 348 +- 995    F %.2f / 12.41\n", b.pm.N, b.pm.S, b.pm.eS, b.spr.half, b.pm.F);
        P("   peak nominal + vertex  N %.0f / 45433   S %.0f +- %.0f +- %.0f / 3949 +- 276 +- 936    F %.2f / 14.31\n", c.pm.N, c.pm.S, c.pm.eS, c.spr.half, c.pm.F);
        P("   the same + veto (opt.) N %.0f / 37974   S %.0f +- %.0f +- %.0f / 3826 +- 252 +- 879    F %.2f / 15.16\n", s.pm.N, s.pm.S, s.pm.eS, s.spr.half, s.pm.F);
        P("   vertex cuts, no kin.   N %.0f / 55947   S %.0f +- %.0f +- %.0f / 4388 +- 307 +- 1078   F %.2f / 14.31\n", v.pm.N, v.pm.S, v.pm.eS, v.spr.half, v.pm.F);
        P("   bump nominal + vertex  N(bump) %.0f / 77735   W %.0f / 73770\n", g.bN, g.W);
        P("   bump, no cut           W %.0f +- %.0f / 148091 +- 469   E %.0f +- %.0f +- %.0f / 1241 +- 76 +- 319   T %.0f +- %.0f +- %.0f / 3301 +- 108 +- 125\n",
          a.W, a.eW, a.E, a.eE, a.sprE, a.tg.S, a.tg.eS, a.sprT);
        const Eff t = EffT(e, a), x = EffE(e, a);
        P("   bump nominal (kin.)    W %.0f / 101602   E %.0f +- %.0f +- %.0f / 1123 +- 64 +- 356   eff_T %.3f +- %.3f +- %.3f / 0.911 +- 0.014 +- 0.011   eff_E %.3f +- %.3f +- %.3f / 0.905 +- 0.030 +- 0.044\n",
          e.W, e.E, e.eE, e.sprE, t.v, t.e, t.spr, x.v, x.e, x.spr);
    }
    Close();

    // ---- cut flow tables ----
    auto peakLine = [&](const char* role, const char* cut, const Res& x) {
        P("   %-9s %-34s | %6.0f %6.0f %6.3f %7.0f %7.0f | %6.0f %5.0f %5.0f | %6.4f %6.2f | %5.1f/%-2d | %6.0f %5.0f\n",
          role, cut, x.pm.N, x.pm.CB, x.pm.N > 0 ? x.pm.CB / x.pm.N : 0, x.D(), x.pm.B, x.pm.S, x.pm.eS, x.spr.half,
          x.pm.SoB, x.pm.F, x.pm.chi2, x.pm.ndf, x.pp.ok ? x.pp.S : 0., x.pp.ok ? x.pp.eS : 0.);
    };
    auto bumpLine = [&](const char* role, const char* cut, const Res& x) {
        const Eff t = EffT(x, R0), e = EffE(x, R0);
        const double ew = R0.W > 0 ? x.W / R0.W : 0, en = R0.bN > 0 ? x.bN / R0.bN : 0;
        P("   %-9s %-34s | %7.0f %6.0f %6.3f %7.0f %5.3f %5.3f | %5.0f %4.0f %4.0f %6.1f/%-3d | %5.3f %5.3f | %5.0f %4.0f | %5.3f %5.3f | %5.2f %4.2f | %5.3f %5.3f\n",
          role, cut, x.bN, x.bCB, x.bN > 0 ? x.bCB / x.bN : 0, x.W, ew, en, x.E, x.eE, x.sprE, x.Echi2, x.Endf, e.v, e.e,
          x.tg.S, x.tg.eS, t.v, t.e, x.Q(), x.eQ(), en > 0 ? e.v / std::sqrt(en) : 0, en > 0 ? t.v / std::sqrt(en) : 0);
    };
    const char* peakHead = "   role      cut of this line                   |      N     CB   CB/N    N-CB       B |      S  stat model |    S/B  F_stat | chi2/ndf | S_pois   err\n";
    const char* bumpHead = "   role      cut of this line                   |       N     CB   CB/N       W  W/W0  N/N0 |     E stat modl   chi2/ndf | eff_E   err |     T stat | eff_T   err |     Q  err | gainE gainT\n";

    Open("ribbons_peak.txt");
    P("#\n# RIBBON CUT FLOW, ETA PEAK. Peak window %.5f-%.5f. N, CB: rows in the window; B = fitted pol2 under the peak (N - CB - B = S in the window);\n", scan::kPeakLo, scan::kPeakHi);
    P("# S = whole Gaussian (P_main, frozen 0.5439 / 5.4 MeV, 0.5225-0.6100), stat = fit error, model = half of (max - min) over the kept model fits.\n");
    P("# F_stat = S/stat is statistical only with a frozen peak shape: it compares selections, it is NOT a significance. CB = 2 sqrt(N++ N--), no k factor: an estimate.\n");
    P("# 'central' = cuts 1..j at the nominal values; 'band' = cut j alone moved; 'if-appl.' = a cut that is NOT applied, put on top of the cuts so far.\n#\n");
    P("%s", peakHead);
    for (int j = 0; j < fp.n(); ++j) {
        const Step& s = fp.st[j];
        peakLine("central", j == 0 ? "preselection" : (s.applied ? CutText(s.cut, s.nom).c_str() : Form("(%s not applied)", CutText(s.cut, s.nom).c_str())), R(fp.cum[j]));
        for (size_t a = 0; a < s.alt.size(); ++a) peakLine(s.applied ? "  band" : "  if-appl.", CutText(s.cut, s.alt[a]).c_str(), R(fp.alt[j][a]));
    }
    P("#\n# bands (smallest - largest value of the quantity over the lines of a step), relative to the central value [%%]\n");
    P("   %-34s | %-15s %-15s %-15s %-15s %-15s %-15s\n", "step", "N", "CB", "B", "S", "S/B", "F_stat");
    {
        const std::vector<Qf> q = {[](const Res& x) { return x.pm.N; }, [](const Res& x) { return x.pm.CB; }, [](const Res& x) { return x.pm.B; },
                                   [](const Res& x) { return x.pm.S; }, [](const Res& x) { return x.pm.SoB; }, [](const Res& x) { return x.pm.F; }};
        std::vector<Ser> ss;
        for (const Qf& f : q) ss.push_back(Series(fp, f));
        for (int j = 1; j < fp.n(); ++j) {
            P("   %-34s |", (std::string(fp.st[j].applied ? "" : "(n.a.) ") + CutText(fp.st[j].cut, fp.st[j].nom)).c_str());
            for (const Ser& s : ss) P(" %+6.1f .. %+5.1f", 100 * (s.lo[j] / s.c[j] - 1), 100 * (s.hi[j] / s.c[j] - 1));
            P("\n");
        }
    }
    Close();

    Open("ribbons_bump.txt");
    P("#\n# RIBBON CUT FLOW, BUMP. Bump window %.4f-%.4f: N, CB, W = N - CB. E = sum of the edge-step heights in four M(pipi) slices [rows / 2.5 MeV]\n", scan::kBumpLo, scan::kBumpHi);
    P("# (stat = fit error, modl = model spread, chi2/ndf summed; E is not a yield). T = tagged eta peak content in m5c (pi0 window %.3f-%.3f).\n", kPi0Lo, kPi0Hi);
    P("# eff_E, eff_T = E, T relative to the preselection (error with the common rows taken into account). Q = E/sqrt(N). gainE = eff_E/sqrt(N/N0), gainT = eff_T/sqrt(N/N0)\n");
    P("# (N = all rows in the bump window: the noise is sqrt(N); the same definition as in scan_bump.C and scan_extra.C). W = N - CB is an estimate (no k factor).\n");
    P("# 'central' = cuts 1..j at the nominal values; 'band' = cut j alone moved; 'if-appl.' = a cut that is NOT applied, put on top of the cuts so far.\n#\n");
    P("%s", bumpHead);
    for (int j = 0; j < fb.n(); ++j) {
        const Step& s = fb.st[j];
        bumpLine("central", j == 0 ? "preselection" : (s.applied ? CutText(s.cut, s.nom).c_str() : Form("(%s not applied)", CutText(s.cut, s.nom).c_str())), R(fb.cum[j]));
        for (size_t a = 0; a < s.alt.size(); ++a) bumpLine(s.applied ? "  band" : "  if-appl.", CutText(s.cut, s.alt[a]).c_str(), R(fb.alt[j][a]));
    }
    Close();

    Open("ribbons_nm1.txt");
    P("#\n# N-1: each cut of the full selection removed in turn ('removed'); cuts that the selection does not use are added at their customary value ('added').\n");
    P("#\n# ETA PEAK, full selection: %s\n", Describe(fp.full()).c_str());
    P("%s", peakHead);
    for (int j = 0; j < fp.n(); ++j)
        peakLine(j == 0 ? "full" : (fp.st[j].applied ? "removed" : "added"), j == 0 ? "-" : CutText(fp.st[j].cut, fp.st[j].nom).c_str(), R(fp.nm1[j]));
    P("#\n# BUMP, full selection: %s\n", Describe(fb.full()).c_str());
    P("%s", bumpHead);
    for (int j = 0; j < fb.n(); ++j)
        bumpLine(j == 0 ? "full" : (fb.st[j].applied ? "removed" : "added"), j == 0 ? "-" : CutText(fb.st[j].cut, fb.st[j].nom).c_str(), R(fb.nm1[j]));
    Close();

    // ---- summary table ----
    Open("ribbons_table.txt");
    P("#\n# SUMMARY TABLE. Peak window %.5f-%.5f, bump window %.4f-%.4f.\n", scan::kPeakLo, scan::kPeakHi, scan::kBumpLo, scan::kBumpHi);
    P("# S = peak content (whole Gaussian) +- statistical error +- model spread (half of max - min over the kept model fits). F_stat = S/stat.\n");
    P("# S/B = S in the window / (fitted polynomial + CB). E = edge-step sum [rows / 2.5 MeV] +- stat +- model spread: a visibility measure, NOT a yield;\n");
    P("# 'not valid' = the edge fits do not describe the spectra (summed chi2 probability below 0.001). eff_T = tagged-eta content relative to 'none'\n");
    P("# +- stat +- spread over 18 fit variations; n/a where the selection vetoes the tagging photon.\n");
    P("# READING: F_stat is statistical only with a frozen peak shape -- it compares selections and is NOT a significance (S/sqrt(stat^2 + model^2) is 3.6-4.8).\n");
    P("# CB = 2 sqrt(N++ N--) without a k factor; N--/N++ is 5-6 and like-sign pairs at small opening angle are almost absent: CB, CB/N and N - CB are estimates\n");
    P("# (N - CB of the bump window moves by -6 %% if the like-sign sum is used). '+ vertex' = vz > %.0f mm and both pions used in the vertex fit: PROPOSED;\n", kVzCut);
    P("# vz > %.0f tightens the binding preselection (z > -500) and needs a decision; the efficiency of isUsedVertex must come from simulation.\n#\n", kVzCut);
    P("   %-22s | %6s %6s %5s | %-19s | %6s %6s | %7s | %-28s | %-21s\n", "selection", "N", "CB", "CB/N", "S +- stat +- model", "S/B", "F_stat", "N - CB", "E +- stat +- model", "eff_T +- stat +- sprd");
    for (const WP& w : {wNone, wSlides, wOld, wPk, wPf, wBk, wBf}) {
        const Res& x = R(w);
        const Eff t = EffT(x, R0);
        std::string es = Form("%4.0f +- %2.0f +- %3.0f", x.E, x.eE, x.sprE);
        if (!x.Evalid()) es += Form(" not valid");
        std::string ts = On(w.v[cPi0]) ? "n/a (photon veto)" : Form("%5.3f +- %5.3f +- %5.3f", t.v, t.e, t.spr);
        P("   %-22s | %6.0f %6.0f %5.3f | %4.0f +- %3.0f +- %4.0f | %6.4f %6.2f | %7.0f | %-28s | %-21s\n", w.name.c_str(),
          x.pm.N, x.pm.CB, x.pm.CB / x.pm.N, x.pm.S, x.pm.eS, x.spr.half, x.pm.SoB, x.pm.F, x.W, es.c_str(), ts.c_str());
    }
    P("#\n# for information (NOT proposed): the two vertex cuts without any kinematic cut, and the optional pi0-photon veto on top of the peak selection\n");
    for (const WP& w : {wVtx, wPv}) {
        const Res& x = R(w);
        const Eff t = EffT(x, R0);
        std::string es = Form("%4.0f +- %2.0f +- %3.0f", x.E, x.eE, x.sprE);
        if (!x.Evalid()) es += Form(" not valid");
        std::string ts = On(w.v[cPi0]) ? "n/a (photon veto)" : Form("%5.3f +- %5.3f +- %5.3f", t.v, t.e, t.spr);
        P("   %-22s | %6.0f %6.0f %5.3f | %4.0f +- %3.0f +- %4.0f | %6.4f %6.2f | %7.0f | %-28s | %-21s\n", w.name.c_str(),
          x.pm.N, x.pm.CB, x.pm.CB / x.pm.N, x.pm.S, x.pm.eS, x.spr.half, x.pm.SoB, x.pm.F, x.W, es.c_str(), ts.c_str());
    }
    P("#\n# more columns for the same selections\n");
    P("   %-22s | %-9s %5s | %6s %5s | %-11s %-11s | %7s %5s | %10s %6s | %5s %5s | %6s %5s %6s\n", "selection", "chi2/ndf", "kept", "S_pois", "err",
      "null ++", "null --", "N(bump)", "CB/N", "E chi2/ndf", "p", "Q", "err", "T", "stat", "spread");
    for (const WP& w : {wNone, wSlides, wOld, wPk, wPf, wBk, wBf}) {
        const Res& x = R(w);
        P("   %-22s | %5.1f/%-3d %5d | %6.0f %5.0f | %4.0f+-%-5.0f %4.0f+-%-5.0f | %7.0f %5.3f | %6.1f/%-3d %6.4f | %5.2f %5.2f | %6.0f %5.0f %6.0f\n", w.name.c_str(),
          x.pm.chi2, x.pm.ndf, x.spr.npass, x.pp.ok ? x.pp.S : 0., x.pp.ok ? x.pp.eS : 0., x.l1.ok ? x.l1.S : 0., x.l1.ok ? x.l1.eS : 0.,
          x.l2.ok ? x.l2.S : 0., x.l2.ok ? x.l2.eS : 0., x.bN, x.bN > 0 ? x.bCB / x.bN : 0, x.Echi2, x.Endf, x.Eprob(), x.Q(), x.eQ(),
          x.tg.S, x.tg.eS, x.sprT);
    }
    P("#\n# check of the frozen peak shape (0.5439, 5.4 MeV): P_main with position and width free\n");
    P("   %-22s | %-19s %-15s | %6s %5s | %8s\n", "selection", "position", "width [MeV]", "S", "err", "chi2/ndf");
    for (const WP& w : {wNone, wSlides, wOld, wPk, wPf, wBk, wBf}) {
        const Res& x = R(w);
        if (!x.fr.ok) { P("   %-22s | free fit did not converge\n", w.name.c_str()); continue; }
        P("   %-22s | %.5f +- %.5f  %5.2f +- %4.2f   | %6.0f %5.0f | %5.1f/%-2d\n", w.name.c_str(), x.fr.mean, x.fr.eMean, 1e3 * x.fr.sigma, 1e3 * x.fr.eSigma,
          x.fr.S, x.fr.eS, x.fr.chi2, x.fr.ndf);
    }
    P("#\n# gains of the two final selections over 'none'\n");
    {
        const Res& x = R(wPf);
        const double r = x.pm.eS / R0.pm.eS;
        P("   peak nominal + vertex: F_stat %.2f -> %.2f, difference %.2f +- %.2f (error of the difference of nested samples, sqrt(2 (1 - sigma_pass/sigma_all)));\n",
          R0.pm.F, x.pm.F, x.pm.F - R0.pm.F, std::sqrt(std::max(0., 2 * (1 - r))));
        {
            const Res& k = R(wPk); const Res& v = R(wVtx);
            P("                          the gain comes from the vertex cuts, not from the kinematic cuts: kinematic cuts alone %.2f -> %.2f (%+.2f +- %.2f),\n",
              R0.pm.F, k.pm.F, k.pm.F - R0.pm.F, std::sqrt(std::max(0., 2 * (1 - k.pm.eS / R0.pm.eS))));
            P("                          vertex cuts alone %.2f -> %.2f (%+.2f +- %.2f)\n", R0.pm.F, v.pm.F, v.pm.F - R0.pm.F,
              std::sqrt(std::max(0., 2 * (1 - v.pm.eS / R0.pm.eS))));
        }
        P("                          S/B %.4f -> %.4f, CB share %.3f -> %.3f, S kept %.3f, rows in the window kept %.3f\n",
          R0.pm.SoB, x.pm.SoB, R0.pm.CB / R0.pm.N, x.pm.CB / x.pm.N, x.pm.S / R0.pm.S, x.pm.N / R0.pm.N);
        P("                          F_tot = S / sqrt(stat^2 + model^2): %.2f -> %.2f\n", R0.pm.S / R0.Stot(), x.pm.S / x.Stot());
        {
            // the two vertex cuts of the peak selection measured on the tagged eta (the photon veto cannot be measured there)
            const Res& a = R(fp.cum[5]); const Res& b = R(fp.cum[7]);
            const Eff t = EffT(b, a);
            double es = 0, ees = 0; SubRatio(b.pm.S, b.pm.eS, a.pm.S, a.pm.eS, es, ees);
            P("                          vz > -160 and both pions in the vertex fit, on top of the kinematic cuts: S kept %.3f +- %.3f, tagged eta kept %.3f +- %.3f +- %.3f,\n", es, ees, t.v, t.e, t.spr);
            const Res& pv = R(wPv);
            P("                          rows in the peak window kept %.3f; the optional photon veto (NOT applied) would keep S %.3f, rows %.3f, F_stat %.2f -> %.2f (+- %.2f)\n",
              b.pm.N / a.pm.N, pv.pm.S / b.pm.S, pv.pm.N / b.pm.N, b.pm.F, pv.pm.F, std::sqrt(std::max(0., 2 * (1 - pv.pm.eS / b.pm.eS))));
        }
        for (const WP& w : {wBk, wBf}) {
            const Res& b = R(w);
            const Eff t = EffT(b, R0), e = EffE(b, R0);
            const double ew = b.W / R0.W, en = b.bN / R0.bN;
            P("   %-22s: rows kept %.3f (N - CB kept %.3f), CB share %.3f -> %.3f; eff_E %.3f +- %.3f +- %.3f, eff_T %.3f +- %.3f +- %.3f;\n", w.name.c_str(), en, ew,
              R0.bCB / R0.bN, b.bCB / b.bN, e.v, e.e, e.spr, t.v, t.e, t.spr);
            P("                          gain = kept / sqrt(rows kept): from E %.3f +- %.3f +- %.3f, from T %.3f +- %.3f +- %.3f; Q %.2f +- %.2f -> %.2f +- %.2f; tagged S/B %.2f -> %.2f\n",
              e.v / std::sqrt(en), e.e / std::sqrt(en), e.spr / std::sqrt(en), t.v / std::sqrt(en), t.e / std::sqrt(en), t.spr / std::sqrt(en),
              R0.Q(), R0.eQ(), b.Q(), b.eQ(), R0.tg.SoB, b.tg.SoB);
        }
    }
    P("#\n# checks of the final selections: two halves of the event groups (fold 0 / 1) and the comparison with isBest == 1\n");
    P("   %-34s | %6s %5s %6s %6s | %7s %5s %4s %5s %4s\n", "selection", "S", "stat", "F_stat", "S/B", "W", "E", "stat", "T", "stat");
    for (const WP& w : chk) {
        const Res& x = R(w);
        const std::string nm = w.name + (w.best ? ", isBest == 1" : Form(", fold %d", w.fold));
        P("   %-34s | %6.0f %5.0f %6.2f %6.4f | %7.0f %5.0f %4.0f %5.0f %4.0f\n", nm.c_str(), x.pm.S, x.pm.eS, x.pm.F, x.pm.SoB, x.W, x.E, x.eE, x.tg.S, x.tg.eS);
    }
    Close();

    // ------------------------------------------------------------------------
    // 3. ribbon cut flow and N-1 figures
    // ------------------------------------------------------------------------
    for (int nm = 0; nm < 2; ++nm) {
        const char* howto = nm ? "Dashed lines: value with all cuts. Full markers: that cut removed. Open markers in shaded columns: a cut that the selection does not use, added at its customary value."
                               : "Line: cuts 1..j at their nominal values. Band at step j: cut j alone moved by one step down and up. Shaded columns: cut not used, line stays flat; dashed box = range if it were applied.";
        // ================= PEAK =================
        {
            const Flow& f = fp;
            const int n = f.n();
            auto ser = [&](const Qf& q, const Qf& e = Qf()) { return nm ? SeriesNm1(f, q, e) : Series(f, q, e); };
            const Ser sN = ser([](const Res& x) { return x.pm.N; }, [](const Res& x) { return x.pm.eN; });
            const Ser sCB = ser([](const Res& x) { return x.pm.CB; }, [](const Res& x) { return x.pm.eCB; });
            const Ser sD = ser([](const Res& x) { return x.D(); });
            const Ser sB = ser([](const Res& x) { return x.pm.B; });
            const Ser sS = ser([](const Res& x) { return x.pm.S; }, [](const Res& x) { return x.Stot(); });
            const Ser sSstat = ser([](const Res& x) { return x.pm.S; }, [](const Res& x) { return x.pm.eS; });
            const Ser sSoB = ser([](const Res& x) { return x.pm.SoB; }, [](const Res& x) { return x.pm.S > 0 ? x.pm.SoB * x.pm.eS / x.pm.S : 0.; });
            // error of F: error of the difference to the first column (nested samples)
            const double eRef = nm ? R(f.nm1[0]).pm.eS : R0.pm.eS;
            const Ser sF = ser([](const Res& x) { return x.pm.F; }, [eRef](const Res& x) {
                const double a = std::min(x.pm.eS, eRef), b = std::max(x.pm.eS, eRef);
                return std::sqrt(std::max(0., 2 * (1 - a / b))); });
            const Opt oN = MkOpt(kColN, 20), oCB = MkOpt(kColCB, 21, 1.5), oB = MkOpt(kColB, 22, 1.7);
            const Opt oD = MkOpt(kColS, 0, 0, false, false, false, 0, 7, 3);
            Opt oS = MkOpt(kColS, 20, 1.6, false, true, true);
            const Opt oR = MkOpt(kColS, 20, 1.6, true);

            Pads pd = MakePads(nm ? "ribbons_nm1_peak" : "ribbons_cutflow_peak");
            // ---- header ----
            pd.head->cd();
            Txt(0.012, 0.80, nm ? "#eta peak in M(#pi^{+}#pi^{-}e^{+}e^{-}): each cut removed in turn from the full peak selection (N-1)"
                                : "#eta peak in M(#pi^{+}#pi^{-}e^{+}e^{-}): ribbon cut flow of the peak selection", 34);
            Txt(0.012, 0.56, Form("%s.  Rows in the peak window %.4f - %.4f GeV/c^{2}.", kHdr, scan::kPeakLo, scan::kPeakHi), 23);
            Txt(0.012, 0.36, howto, 21);
            {
                TLegend* lg = new TLegend(0.012, 0.02, 0.99, 0.28);
                lg->SetNColumns(5); lg->SetTextFont(43); lg->SetTextSize(22); lg->SetFillStyle(0); lg->SetMargin(0.14);
                lg->SetColumnSeparation(0.01);
                lg->AddEntry(LegMark(oN, !nm), "all pairs N", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oCB, !nm), "CB = 2#sqrt{N_{++}N_{#minus#minus}} (no k factor)", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oD, false), "N - CB", "l");
                lg->AddEntry(LegMark(oB, !nm), "fitted background under the peak (pol2)", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oS, true, 3345), "peak content S (hatched: stat #oplus model spread)", "pf");
                lg->Draw();
            }
            // ---- main pad ----
            pd.p[0]->cd();
            double lo, hi; Range({sN, sCB, sB, sS}, lo, hi, true);
            const double y1 = 0.50 * lo, y2 = 1.9 * hi;
            Frame(f, Form("fr_pk0_%d", nm), y1, y2, "rows in the peak window", true);
            if (!nm) {
                DrawSer(f, sS, oS); DrawSer(f, sCB, oCB); DrawSer(f, sB, oB); DrawSer(f, sD, oD); DrawSer(f, sN, oN);
            } else {
                DrawNm1(f, sS, oS, &sSstat); DrawNm1(f, sCB, oCB); DrawNm1(f, sB, oB); DrawNm1(f, sD, oD); DrawNm1(f, sN, oN);
            }
            for (int j = 0; j < n; ++j) {
                const Res& x = R(nm ? f.nm1[j] : f.cum[j]);
                const bool app = j == 0 || f.st[j].applied;
                TxtU(j + 0.5, (app ? std::max(sN.hi[j], sN.c[j]) : sN.c[j]) * 1.10, Form("%.0f", sN.c[j]), 20, kColN, 21);
                // CB: above its band when it lies above S, below otherwise; next to the marker when the band is wide
                const bool cbAbove = sCB.c[j] > sS.c[j];
                const bool wide = app && sCB.hi[j] > 1.5 * sCB.lo[j];
                if (wide) TxtU(j + 0.5 + 0.05, sCB.c[j] * (cbAbove ? 1.09 : 1 / 1.09), Form("%.0f", sCB.c[j]), 20, kColCB, cbAbove ? 11 : 13);
                else TxtU(j + 0.5, cbAbove ? (app ? std::max(sCB.hi[j], sCB.c[j]) : sCB.c[j]) * 1.10 : (app ? std::min(sCB.lo[j], sCB.c[j]) : sCB.c[j]) / 1.10,
                          Form("%.0f", sCB.c[j]), 20, kColCB, cbAbove ? 21 : 23);
                TxtU(j + 0.5, y1 * 1.10, Form("S = %.0f #pm %.0f #pm %.0f", x.pm.S, x.pm.eS, x.spr.half), 18, kColS, 21);
            }
            gPad->RedrawAxis();
            // ---- S/B ----
            pd.p[1]->cd();
            Range({sSoB}, lo, hi, true);
            double pad = 0.30 * (hi - lo);
            Frame(f, Form("fr_pk1_%d", nm), lo - pad, hi + 2.3 * pad, "S/B", false);
            if (!nm) DrawSer(f, sSoB, oR); else DrawNm1(f, sSoB, oR);
            for (int j = 0; j < n; ++j) Num(sSoB, j, Form("%.3f", sSoB.c[j]), kColS, true, 0.30 * pad, 0, j == 0 || f.st[j].applied);
            Txt(kLm + 0.008, 0.855, "S/B = S in the window / (fitted background + CB); bars: statistical error of S", 19, kGray + 3);
            gPad->RedrawAxis();
            // ---- F_stat ----
            pd.p[2]->cd();
            Range({sF}, lo, hi, true);
            pad = 0.30 * (hi - lo);
            Frame(f, Form("fr_pk2_%d", nm), lo - pad, hi + 2.3 * pad, "S/#sigma_{S} (stat. only)", false);
            if (!nm) DrawSer(f, sF, oR); else DrawNm1(f, sF, oR);
            for (int j = 0; j < n; ++j) Num(sF, j, Form("%.1f", sF.c[j]), kColS, true, 0.30 * pad, 0, j == 0 || f.st[j].applied);
            Txt(kLm + 0.008, 0.905, nm ? "statistical only, frozen peak shape: compares selections, NOT a significance.  Bars: error of the difference to the full selection (nested samples)"
                                       : "statistical only, frozen peak shape: compares selections, NOT a significance.  Bars: error of the difference to the preselection (nested samples)", 19, kGray + 3);
            StepLabels(f, nm);
            gPad->RedrawAxis();
            Save(pd.c, nm ? "ribbons_nm1_peak" : "ribbons_cutflow_peak");
        }
        // ================= BUMP =================
        {
            const Flow& f = fb;
            const int n = f.n();
            const Res& base = R0;                         // efficiencies are relative to the preselection
            auto ser = [&](const Qf& q, const Qf& e = Qf()) { return nm ? SeriesNm1(f, q, e) : Series(f, q, e); };
            const Ser sN = ser([](const Res& x) { return x.bN; });
            const Ser sCB = ser([](const Res& x) { return x.bCB; });
            const Ser sW = ser([](const Res& x) { return x.W; }, [](const Res& x) { return x.eW; });
            const Ser sEE = ser([&base](const Res& x) { return EffE(x, base).v; }, [&base](const Res& x) { return &x == &base ? 0. : EffE(x, base).e; });
            const Ser sET = ser([&base](const Res& x) { return EffT(x, base).v; }, [&base](const Res& x) { return &x == &base ? 0. : EffT(x, base).e; });
            const Ser sRW = ser([&base](const Res& x) { return std::sqrt(x.bN / base.bN); });
            const double q0 = base.Q();
            const Ser sQ = ser([](const Res& x) { return x.Q(); }, [&base](const Res& x) {
                if (&x == &base) return 0.;
                const double ew = x.bN / base.bN; return base.Q() * EffE(x, base).e / std::sqrt(ew); });
            const Ser sQT = ser([&base](const Res& x) { const double ew = x.bN / base.bN; return base.Q() * EffT(x, base).v / std::sqrt(ew); },
                                [&base](const Res& x) { if (&x == &base) return 0.; const double ew = x.bN / base.bN; return base.Q() * EffT(x, base).e / std::sqrt(ew); });
            const double dxE = -0.07, dxT = 0.07;
            const Opt oN = MkOpt(kColN, 20), oCB = MkOpt(kColCB, 21, 1.5), oW = MkOpt(kColS, 20);
            const Opt oE = MkOpt(kColE, 21, 1.5, true, true, false, dxE), oT = MkOpt(kColT, 20, 1.6, true, true, false, dxT);
            const Opt oBE = MkOpt(kGray + 2, 0, 0, false, false, false, 0, 2, 3);

            Pads pd = MakePads(nm ? "ribbons_nm1_bump" : "ribbons_cutflow_bump");
            pd.head->cd();
            Txt(0.012, 0.80, nm ? "Bump (#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}, e^{+}e^{-} from the #pi^{0}: Dalitz decay or converted #gamma) in M(#pi^{+}#pi^{-}e^{+}e^{-}): each cut removed in turn from the bump selection (N-1)"
                                : "Bump (#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}, e^{+}e^{-} from the #pi^{0}: Dalitz decay or converted #gamma) in M(#pi^{+}#pi^{-}e^{+}e^{-}): ribbon cut flow of the bump selection", 32);
            Txt(0.012, 0.56, Form("%s.  Rows in the bump window %.3f - %.4f GeV/c^{2}.", kHdr, scan::kBumpLo, scan::kBumpHi), 23);
            Txt(0.012, 0.36, howto, 21);
            {
                TLegend* lg = new TLegend(0.012, 0.02, 0.99, 0.28);
                lg->SetNColumns(6); lg->SetTextFont(43); lg->SetTextSize(21); lg->SetFillStyle(0); lg->SetMargin(0.15);
                lg->SetColumnSeparation(0.01);
                lg->AddEntry(LegMark(oN, !nm), "all pairs N", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oCB, !nm), "CB = 2#sqrt{N_{++}N_{#minus#minus}} (no k)", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oW, !nm), "W = N - CB", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oE, !nm), "E: edge-step sum (not a yield)", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oT, !nm), "T: ECAL-tagged #eta peak", nm ? "p" : "lpf");
                lg->AddEntry(LegMark(oBE, false), "#sqrt{N/N_{0}} (break-even)", "l");
                lg->Draw();
            }
            // ---- main pad ----
            pd.p[0]->cd();
            double lo, hi; Range({sN, sCB, sW}, lo, hi, false);
            Frame(f, Form("fr_bp0_%d", nm), 0.62 * lo, 1.9 * hi, "rows in the bump window", true);
            if (!nm) { DrawSer(f, sCB, oCB); DrawSer(f, sW, oW); DrawSer(f, sN, oN); }
            else { DrawNm1(f, sCB, oCB); DrawNm1(f, sW, oW); DrawNm1(f, sN, oN); }
            for (int j = 0; j < n; ++j) {
                const bool app = j == 0 || f.st[j].applied;
                TxtU(j + 0.5, (app ? std::max(sN.hi[j], sN.c[j]) : sN.c[j]) * 1.12, Form("%.0f", sN.c[j]), 20, kColN, 21);
                TxtU(j + 0.5, (app ? std::min(sW.lo[j], sW.c[j]) : sW.c[j]) / 1.14, Form("%.0f", sW.c[j]), 20, kColS, 23);
                TxtU(j + 0.5, (app ? std::max(sCB.hi[j], sCB.c[j]) : sCB.c[j]) * 1.10, Form("%.0f", sCB.c[j]), 20, kColCB, 21);
            }
            gPad->RedrawAxis();
            // ---- E and T relative to the preselection ----
            pd.p[1]->cd();
            Range({sEE, sET, sRW}, lo, hi, true);
            double pad = 0.12 * (hi - lo);
            Frame(f, Form("fr_bp1_%d", nm), lo - pad, hi + 3.4 * pad, "kept fraction", false);
            if (!nm) { DrawSer(f, sRW, oBE); DrawSer(f, sEE, oE); DrawSer(f, sET, oT); }
            else {
                for (int j = 0; j < n; ++j) HLine(j + 0.5 - kHalfW, j + 0.5 + kHalfW, sRW.c[j], kGray + 2, 2, 3);
                DrawNm1(f, sEE, oE); DrawNm1(f, sET, oT);
            }
            Txt(kLm + 0.008, 0.855, "E/E_{0} (blue) and T/T_{0} (green), relative to the preselection. A measure above the dashed line #sqrt{N/N_{0}} (rows kept) means a gain.  E does not see cuts on M(#pi#pi) outside 0.30 - 0.40.", 19, kGray + 3);
            gPad->RedrawAxis();
            // ---- Q = E / sqrt(W) ----
            pd.p[2]->cd();
            Range({sQ, sQT}, lo, hi, true);
            pad = 0.25 * (hi - lo);
            const double ya = lo - 1.3 * pad, yb = hi + 2.4 * pad;
            Frame(f, Form("fr_bp2_%d", nm), ya, yb, "E/#sqrt{N}", false);
            HLine(0, n, q0, kGray + 2, 3, 1);
            if (!nm) { DrawSer(f, sQ, oE); DrawSer(f, sQT, oT); }
            else { DrawNm1(f, sQ, oE); DrawNm1(f, sQT, oT); }
            for (int j = 0; j < n; ++j) {
                const bool eTop = sQ.c[j] > sQT.c[j];
                const bool app = j == 0 || f.st[j].applied;
                Num(sQ, j, Form("%.2f", sQ.c[j] / q0), kColE, eTop, 0.25 * pad, dxE, app);
                Num(sQT, j, Form("%.2f", sQT.c[j] / q0), kColT, !eTop, 0.25 * pad, dxT, app);
            }
            {
                TGaxis* ax = new TGaxis(n, ya, n, yb, ya / q0, yb / q0, 506, "+L");
                ax->SetLabelFont(43); ax->SetLabelSize(28); ax->SetTitleFont(43); ax->SetTitleSize(30); ax->SetTitleOffset(1.55);
                ax->SetTickSize(0.004); ax->CenterTitle(true);
                ax->SetTitle("gain");
                ax->Draw();
            }
            Txt(kLm + 0.008, 0.905, Form("blue: E/#sqrt{N} (preselection: %.2f #pm %.2f);  green: the same from the tagged #eta, %.2f #times (T/T_{0}) / #sqrt{N/N_{0}};  numbers and right axis: gain over the preselection;  bars: error relative to the preselection",
                                         q0, base.eQ(), q0), 19, kGray + 3);
            StepLabels(f, nm);
            gPad->RedrawAxis();
            Save(pd.c, nm ? "ribbons_nm1_bump" : "ribbons_cutflow_bump");
        }
    }

    // ------------------------------------------------------------------------
    // 4. final spectra: none, slides, final selection
    // ------------------------------------------------------------------------
    TFile* fout = TFile::Open(Form("%s/ribbons_spectra.root", scan::kScanDir), "RECREATE");
    for (int tg = 0; tg < 2; ++tg) {
        const WP sel[3] = {wNone, wSlides, tg == 0 ? wPf : wBf};
        const char* ttl[3] = {"none", "slides", tg == 0 ? "peak selection" : "bump selection"};
        const char* key[3] = {"none", "slides", tg == 0 ? "peak" : "bump"};
        std::vector<Sel> sels, ssl, sst;
        for (int i = 0; i < 3; ++i) {
            const WP w = sel[i];
            sels.push_back([w](const Row& r) { return Pass(w, r); });
            const float lo = kSlLo[2], hi = kSlHi[2];
            ssl.push_back([w, lo, hi](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && Pass(w, r); });
            sst.push_back([w](const Row& r) { return Pass(w, r); });
        }
        const std::vector<Spectra> sp = scan::FillMany(d, sels);
        const std::vector<Spectra> sl = scan::FillMany(d, ssl);
        const std::vector<Spectra> st = scan::FillMany(dT, sst, TagOpt());
        const int nrow = tg == 0 ? 3 : 4;
        TCanvas* c = new TCanvas(tg == 0 ? "ribbons_spectra_peak" : "ribbons_spectra_bump", "", 2100, tg == 0 ? 2050 : 2650);
        const double top = tg == 0 ? 0.055 : 0.043;
        TPad* head = new TPad("head", "", 0, 1 - top, 1, 1); head->SetMargin(0, 0, 0, 0); head->Draw(); head->cd();
        Txt(0.012, 0.52, tg == 0 ? "M(#pi^{+}#pi^{-}e^{+}e^{-}) for the #eta-peak selection, next to 'none' and 'slides'"
                                 : "M(#pi^{+}#pi^{-}e^{+}e^{-}) for the bump selection, next to 'none' and 'slides'", 40);
        Txt(0.012, 0.10, Form("%s.  Black: all pairs, red: CB = 2#sqrt{N_{++}N_{#minus#minus}}, blue: all - CB.", kHdr), 27);
        TPad* pp[4][3];
        for (int r = 0; r < nrow; ++r)
            for (int i = 0; i < 3; ++i) {
                c->cd();
                const double h = (1 - top) / nrow;
                pp[r][i] = new TPad(Form("p%d%d", r, i), "", i / 3., std::max(0., (1 - top) - (r + 1) * h), (i + 1) / 3., (1 - top) - r * h);
                pp[r][i]->SetMargin(0.19, 0.045, 0.15, 0.13);
                pp[r][i]->Draw();
            }
        auto ax = [](TH1* h, const char* xt, const char* yt) {
            h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
            AxPx(h->GetXaxis(), 28, 30, 1.25); AxPx(h->GetYaxis(), 26, 30, 2.2);
            h->GetXaxis()->SetNdivisions(506);
        };
        auto title = [&](int i, const char* what) {
            std::string l1, l2; DescribeTex(sel[i], l1, l2);
            Txt(0.19, 0.955, Form("#bf{%s}  %s", ttl[i], what), 23);
            Txt(0.19, 0.915, l2.empty() ? l1.c_str() : Form("%s,", l1.c_str()), 19, kGray + 3);
            if (!l2.empty()) Txt(0.19, 0.882, l2.c_str(), 19, kGray + 3);
        };
        // ---- row 0: log, full range, common y ----
        TH1D* hu[3]; TH1D* hc[3]; TH1D* hs[3];
        double ymax = 0, ymin = 1e30;
        for (int i = 0; i < 3; ++i) {
            hu[i] = sp[i].Hist("us", Form("h_%s_%s_all", tg ? "bump" : "peak", key[i]));
            hc[i] = sp[i].Hist("cb", Form("h_%s_%s_cb", tg ? "bump" : "peak", key[i]));
            hs[i] = sp[i].Hist("sig", Form("h_%s_%s_sig", tg ? "bump" : "peak", key[i]));
            DecoH(hu[i], kColN, 20, 0.6); DecoH(hc[i], kColCB, 21, 0.5); DecoH(hs[i], kColS, 24, 0.6);
            ymax = std::max(ymax, hu[i]->GetMaximum());
            for (int b = 1; b <= hc[i]->GetNbinsX(); ++b) if (hc[i]->GetBinContent(b) > 0) ymin = std::min(ymin, hc[i]->GetBinContent(b));
            fout->cd();
            hu[i]->SetTitle(Form("%s, all pairs;M(pi+pi-e+e-) [GeV/c2];rows / 2.5 MeV", Describe(sel[i]).c_str())); hu[i]->Write();
            hc[i]->SetTitle(Form("%s, CB;M(pi+pi-e+e-) [GeV/c2];rows / 2.5 MeV", Describe(sel[i]).c_str())); hc[i]->Write();
            hs[i]->SetTitle(Form("%s, all - CB;M(pi+pi-e+e-) [GeV/c2];rows / 2.5 MeV", Describe(sel[i]).c_str())); hs[i]->Write();
        }
        ymin = std::max(ymin, 3.0);
        for (int i = 0; i < 3; ++i) {
            pp[0][i]->cd(); gPad->SetLogy();
            TH1D* h = (TH1D*)hu[i]->Clone(Form("fr0_%d_%d", tg, i));
            ax(h, kXt, kYt);
            h->SetMinimum(0.5 * ymin); h->SetMaximum(15 * ymax);
            h->Draw("e");
            hc[i]->Draw("e same"); hs[i]->Draw("e same"); h->Draw("e same");
            const double a1 = tg ? scan::kBumpLo : scan::kPeakLo, a2 = tg ? scan::kBumpHi : scan::kPeakHi;
            VLine(a1, 0.5 * ymin, 0.9 * ymax, kGray + 2, 2, 2); VLine(a2, 0.5 * ymin, 0.9 * ymax, kGray + 2, 2, 2);
            title(i, "- full range, log scale");
            if (i == 0) {
                TLegend* lg = new TLegend(0.52, 0.19, 0.93, 0.40);
                lg->SetTextFont(43); lg->SetTextSize(23); lg->SetFillStyle(0);
                lg->AddEntry(hu[i], "all pairs", "pe"); lg->AddEntry(hc[i], "CB", "pe"); lg->AddEntry(hs[i], "all - CB", "pe");
                lg->Draw();
            }
            Txt(0.23, 0.82, tg ? "dashed: bump window 0.455 - 0.5225" : "dashed: peak window 0.5336 - 0.5542", 20, kGray + 2);
            Txt(0.23, 0.78, "same y range in the three panels", 20, kGray + 2);
        }
        if (tg == 0) {
            // ---- row 1: linear zoom of N - CB with the fit (own y range, zero suppressed) ----
            for (int i = 0; i < 3; ++i) {
                pp[1][i]->cd();
                const Res& x = R(sel[i]);
                TH1D* h = (TH1D*)hs[i]->Clone(Form("fr1_%d", i));
                DecoH(h, kColS, 20, 0.9);
                ax(h, kXt, "all - CB [rows / 2.5 MeV/c^{2}]");
                h->GetXaxis()->SetRangeUser(0.5001, 0.6199);
                double ylo = 1e30, yhi = -1e30;
                for (int b = 1; b <= h->GetNbinsX(); ++b) {
                    const double m = h->GetBinCenter(b);
                    if (m < 0.50 || m > 0.62) continue;
                    ylo = std::min(ylo, h->GetBinContent(b) - h->GetBinError(b));
                    yhi = std::max(yhi, h->GetBinContent(b) + h->GetBinError(b));
                }
                const double y1 = ylo - 0.12 * (yhi - ylo), y2 = yhi + 0.55 * (yhi - ylo);
                h->SetMinimum(y1); h->SetMaximum(y2);
                h->Draw("e");
                TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
                for (int q = 0; q <= 350; ++q) {
                    const double m = scan::kFitLo + (scan::kFitHi - scan::kFitLo) * q / 350.;
                    gt->SetPoint(q, m, x.pm.Total(m)); gb->SetPoint(q, m, x.pm.Bkg(m));
                }
                gt->SetLineColor(kColS); gt->SetLineWidth(3); gb->SetLineColor(kColS); gb->SetLineWidth(2); gb->SetLineStyle(2);
                gb->Draw("l same"); gt->Draw("l same");
                VLine(scan::kPeakLo, y1, y1 + 0.50 * (y2 - y1), kGray + 2, 2, 1); VLine(scan::kPeakHi, y1, y1 + 0.50 * (y2 - y1), kGray + 2, 2, 1);
                title(i, "- linear zoom with the fit");
                Txt(0.22, 0.82, Form("S = %.0f #pm %.0f (stat) #pm %.0f (model)", x.pm.S, x.pm.eS, x.spr.half), 23, kColS);
                Txt(0.22, 0.775, Form("S/#sigma_{S} = %.1f (stat. only, not a significance)", x.pm.F), 23, kColS);
                Txt(0.22, 0.73, Form("S/B = %.3f,  CB share %.3f (no k),  #chi^{2}/ndf = %.1f/%d", x.pm.SoB, x.pm.CB / x.pm.N, x.pm.chi2, x.pm.ndf), 21, kColS);
                Txt(0.52, 0.19, "y axis: zero suppressed, own range", 20, kGray + 2);
            }
            // ---- row 2: fitted polynomial subtracted, common y ----
            double smin = 1e30, smax = -1e30;
            TH1D* hr[3];
            for (int i = 0; i < 3; ++i) {
                const Res& x = R(sel[i]);
                hr[i] = (TH1D*)hs[i]->Clone(Form("h_peak_%s_sub", key[i]));
                for (int b = 1; b <= hr[i]->GetNbinsX(); ++b) {
                    const double m = hr[i]->GetBinCenter(b);
                    hr[i]->SetBinContent(b, hs[i]->GetBinContent(b) - x.pm.Bkg(m));
                    if (m < scan::kFitLo || m > scan::kFitHi) continue;
                    smin = std::min(smin, hr[i]->GetBinContent(b) - hr[i]->GetBinError(b));
                    smax = std::max(smax, hr[i]->GetBinContent(b) + hr[i]->GetBinError(b));
                }
            }
            for (int i = 0; i < 3; ++i) {
                pp[2][i]->cd();
                const Res& x = R(sel[i]);
                TH1D* h = hr[i];
                DecoH(h, kColS, 20, 0.9);
                ax(h, kXt, "all - CB - polynomial");
                h->GetXaxis()->SetRangeUser(scan::kFitLo + 1e-4, scan::kFitHi - 1e-4);
                h->SetMinimum(smin - 0.08 * (smax - smin)); h->SetMaximum(smax + 0.35 * (smax - smin));
                h->Draw("e");
                HLine(scan::kFitLo, scan::kFitHi, 0, kGray + 2, 2, 1);
                TGraph* gt = new TGraph();
                for (int q = 0; q <= 350; ++q) {
                    const double m = scan::kFitLo + (scan::kFitHi - scan::kFitLo) * q / 350.;
                    gt->SetPoint(q, m, x.pm.Total(m) - x.pm.Bkg(m));
                }
                gt->SetLineColor(kColS); gt->SetLineWidth(3); gt->Draw("l same");
                title(i, "- fitted polynomial subtracted");
                Txt(0.50, 0.82, "same y range in the three panels", 20, kGray + 2);
                Txt(0.50, 0.78, Form("peak height %.0f rows / 2.5 MeV", x.pm.S * scan::GaussFrac(scan::kPeakMean - 0.5 * scan::kBinW, scan::kPeakMean + 0.5 * scan::kBinW, scan::kPeakMean, scan::kPeakSigma)), 20, kColS);
            }
        } else {
            // ---- row 1: linear zoom, common y ----
            for (int i = 0; i < 3; ++i) {
                pp[1][i]->cd();
                const Res& x = R(sel[i]);
                TH1D* h = (TH1D*)hu[i]->Clone(Form("frb1_%d", i));
                DecoH(h, kColN, 20, 0.8);
                ax(h, kXt, kYt);
                h->GetXaxis()->SetRangeUser(0.4001, 0.6199);
                double m = 0;
                for (int b = 1; b <= hu[0]->GetNbinsX(); ++b) if (hu[0]->GetBinCenter(b) > 0.40 && hu[0]->GetBinCenter(b) < 0.62) m = std::max(m, hu[0]->GetBinContent(b));
                const double y2 = 1.45 * m;
                h->SetMinimum(0); h->SetMaximum(y2);
                h->Draw("e");
                TBox* bx = new TBox(scan::kBumpLo, 0, scan::kBumpHi, 0.70 * y2);
                bx->SetFillColor(ShadeCol()); bx->SetFillStyle(1001); bx->Draw();
                TH1D* c1 = (TH1D*)hc[i]->Clone(Form("frb1c_%d", i)); DecoH(c1, kColCB, 21, 0.7); c1->Draw("e same");
                TH1D* s1 = (TH1D*)hs[i]->Clone(Form("frb1s_%d", i)); DecoH(s1, kColS, 20, 0.8); s1->Draw("e same");
                h->Draw("e same");
                VLine(scan::kPeakMean, 0, 0.70 * y2, kGray + 2, 3, 2);
                gPad->RedrawAxis();
                title(i, "- linear zoom");
                Txt(0.22, 0.82, Form("bump window (shaded): N = %.0f, CB/N = %.3f", x.bN, x.bCB / x.bN), 22);
                Txt(0.22, 0.775, Form("N - CB = %.0f #pm %.0f  (not a bump yield)", x.W, x.eW), 22, kColS);
                Txt(0.22, 0.73, "dotted: #eta peak position;  same y range in the three panels", 20, kGray + 2);
                if (On(sel[i].v[cMpHi])) Txt(0.22, 0.685, Form("fall above 0.52: made by M(#pi#pi) < %.2f, NOT the bump shape", sel[i].v[cMpHi]), 20, kRed + 1);
            }
            // ---- row 2: edge step in the slice 0.37-0.39, common y ----
            double emax = 0;
            TH1D* he[3];
            for (int i = 0; i < 3; ++i) {
                he[i] = sl[i].Hist("sig", Form("h_bump_%s_slice037_039_sig", key[i]));
                for (int b = 1; b <= he[i]->GetNbinsX(); ++b) {
                    const double m = he[i]->GetBinCenter(b);
                    if (m > 0.445 && m < 0.575) emax = std::max(emax, he[i]->GetBinContent(b) + he[i]->GetBinError(b));
                }
                fout->cd(); he[i]->SetTitle(Form("%s, 0.37 <= M(pipi) < 0.39, all - CB;M(pi+pi-e+e-) [GeV/c2];rows / 2.5 MeV", Describe(sel[i]).c_str())); he[i]->Write();
            }
            for (int i = 0; i < 3; ++i) {
                pp[2][i]->cd();
                const Res& x = R(sel[i]);
                const StepFit& sf = x.sf[2];
                TH1D* h = he[i];
                DecoH(h, kColS, 20, 0.9);
                ax(h, kXt, "all - CB [rows / 2.5 MeV/c^{2}]");
                h->GetXaxis()->SetRangeUser(0.4451, 0.5749);
                h->SetMinimum(0); h->SetMaximum(1.75 * emax);
                h->Draw("e");
                if (x.okf[2]) {
                    TGraph* g = new TGraph(); TGraph* g0 = new TGraph();
                    for (size_t q = 0; q < sf.tx.size(); ++q) { g->SetPoint((int)q, sf.tx[q], sf.fy[q]); g0->SetPoint((int)q, sf.tx[q], sf.fy0[q]); }
                    g->SetLineWidth(3); g0->SetLineWidth(2); g0->SetLineStyle(2);
                    g0->Draw("l same"); g->Draw("l same");
                    VLine(sf.edge, 0, 0.95 * emax, kGreen + 2, 7, 2);
                    Txt(0.22, 0.82, Form("0.37 #leq M(#pi#pi) < 0.39:  step h = %.0f #pm %.0f", sf.h, sf.eh), 22);
                    Txt(0.22, 0.775, Form("#chi^{2}/ndf = %.1f/%d,  #Delta#chi^{2}(step) = %.1f%s", sf.chi2, sf.ndf, sf.dchi2,
                                         TMath::Prob(sf.chi2, sf.ndf) < 1e-3 ? "  #rightarrow fit NOT valid" : ""), 22,
                        TMath::Prob(sf.chi2, sf.ndf) < 1e-3 ? kRed + 1 : kBlack);
                    Txt(0.22, 0.73, Form("green: predicted edge %.4f; dashed: fit without its step", sf.edge), 20, kGreen + 2);
                    Txt(0.22, 0.685, Form("step-height sum of four slices (not a yield): E = %.0f #pm %.0f #pm %.0f%s", x.E, x.eE, x.sprE, x.Evalid() ? "" : "  (not valid)"), 19, kGray + 3);
                }
                title(i, "- kinematic edge, one M(#pi#pi) slice");
            }
            // ---- row 3: tagged eta in m5c, common y ----
            double tmax = 0, tmin = 0;
            TH1D* ht[3];
            for (int i = 0; i < 3; ++i) {
                ht[i] = st[i].Hist("sig", Form("h_bump_%s_tagged_m5c_sig", key[i]));
                for (int b = 1; b <= ht[i]->GetNbinsX(); ++b) {
                    const double m = ht[i]->GetBinCenter(b);
                    if (m > 0.47 && m < 0.64) { tmax = std::max(tmax, ht[i]->GetBinContent(b) + ht[i]->GetBinError(b)); tmin = std::min(tmin, ht[i]->GetBinContent(b) - ht[i]->GetBinError(b)); }
                }
                fout->cd(); ht[i]->SetTitle(Form("%s, tagged rows, all - CB;m5c [GeV/c2];rows / 2.5 MeV", Describe(sel[i]).c_str())); ht[i]->Write();
            }
            for (int i = 0; i < 3; ++i) {
                pp[3][i]->cd();
                const Res& x = R(sel[i]);
                const Eff e = EffT(x, R0);
                TH1D* h = ht[i];
                DecoH(h, kColS, 20, 0.9);
                ax(h, "M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma), M(e^{+}e^{-}#gamma) fixed to m_{#pi^{0}} [GeV/c^{2}]", "all - CB [rows / 2.5 MeV/c^{2}]");
                h->GetXaxis()->SetRangeUser(0.4701, 0.6399);
                h->SetMinimum(tmin < 0 ? 1.2 * tmin : 0); h->SetMaximum(1.55 * tmax);
                h->Draw("e");
                if (x.tg.ok) {
                    TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
                    for (int q = 0; q <= 300; ++q) {
                        const double m = kTagLo + (kTagHi - kTagLo) * q / 300.;
                        gt->SetPoint(q, m, x.tg.Total(m)); gb->SetPoint(q, m, x.tg.Bkg(m));
                    }
                    gt->SetLineColor(kColS); gt->SetLineWidth(3); gb->SetLineColor(kColS); gb->SetLineWidth(2); gb->SetLineStyle(2);
                    gb->Draw("l same"); gt->Draw("l same");
                    Txt(0.22, 0.82, Form("tagged #eta: T = %.0f #pm %.0f (stat) #pm %.0f (model)", x.tg.S, x.tg.eS, x.sprT), 22, kColS);
                    if (i > 0) Txt(0.22, 0.775, Form("kept: %.3f #pm %.3f #pm %.3f of 'none'", e.v, e.e, e.spr), 22, kColS);
                    Txt(0.22, i > 0 ? 0.73 : 0.775, Form("S/B = %.2f,  #chi^{2}/ndf = %.1f/%d", x.tg.SoB, x.tg.chi2, x.tg.ndf), 22, kColS);
                }
                title(i, "- ECAL-tagged #eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
                Txt(0.22, i > 0 ? 0.685 : 0.73, Form("#pi^{0} window %.3f - %.3f; same y range in the three panels", kPi0Lo, kPi0Hi), 20, kGray + 2);
            }
        }
        Save(c, tg == 0 ? "ribbons_spectra_peak" : "ribbons_spectra_bump");
    }

    // ------------------------------------------------------------------------
    // 5. the summary table as a figure
    // ------------------------------------------------------------------------
    {
        const std::vector<WP> sel = {wNone, wSlides, wOld, wPk, wPf, wBk, wBf};
        const char* lab[7] = {"none", "slides", "old code", "peak", "peak+", "bump", "bump+"};
        const int n = (int)sel.size();
        TCanvas* c = new TCanvas("ribbons_summary", "", 2400, 1500);
        TPad* head = new TPad("sum_head", "", 0, 0.90, 1, 1); head->SetMargin(0, 0, 0, 0); head->Draw(); head->cd();
        Txt(0.012, 0.66, "Selections compared: #eta peak (upper row) and bump (lower row) in M(#pi^{+}#pi^{-}e^{+}e^{-})", 36);
        Txt(0.012, 0.36, Form("%s.  peak = oa_{lab} < 55#circ, M(#pi#pi) > 0.30;  bump = oa_{lab} < 40#circ, 0.30 < M(#pi#pi) < 0.41;", kHdr), 23);
        Txt(0.012, 0.10, "\"+\" = with the two proposed vertex cuts: z_{vtx} > -160 mm (new: needs a decision) and both #pi in the vertex fit; no #pi^{0}-photon veto.  slides, old code: reference chains (with the #eta-rest angle cut).", 22);
        TPad* pq[6];
        for (int i = 0; i < 6; ++i) {
            c->cd();
            pq[i] = new TPad(Form("sum_p%d", i), "", (i % 3) / 3., i < 3 ? 0.45 : 0.0, (i % 3 + 1) / 3., i < 3 ? 0.90 : 0.45);
            pq[i]->SetMargin(0.17, 0.03, 0.14, 0.12);
            pq[i]->Draw();
        }
        auto frame = [&](int ip, double y1, double y2, const char* yt, const char* ttl) {
            pq[ip]->cd();
            TH1F* h = new TH1F(Form("sum_fr%d", ip), "", n, 0, n);
            h->SetDirectory(nullptr);
            for (int i = 0; i < n; ++i) h->GetXaxis()->SetBinLabel(i + 1, lab[i]);
            h->GetXaxis()->LabelsOption("h");
            h->GetXaxis()->SetLabelFont(43); h->GetXaxis()->SetLabelSize(22); h->GetXaxis()->SetLabelOffset(0.012);
            h->GetXaxis()->SetTickLength(0);
            h->GetYaxis()->SetTitle(yt); AxPx(h->GetYaxis(), 28, 30, 2.0); h->GetYaxis()->CenterTitle(true);
            h->SetMinimum(y1); h->SetMaximum(y2);
            h->Draw("axis");
            for (int i = 3; i < n; i += 2) VLine(i, y1, y2, kGray + 1, 3, 1);
            Txt(0.17, 0.915, ttl, 25);
        };
        auto gr = [](int col, int mk, double ms) {
            TGraphErrors* g = new TGraphErrors();
            g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineWidth(2);
            return g;
        };
        auto add = [](TGraphErrors* g, double x, double y, double e) { const int k = g->GetN(); g->SetPoint(k, x, y); g->SetPointError(k, 0, e); };
        // (0) S
        {
            double hi = 0;
            for (const WP& w : sel) hi = std::max(hi, R(w).pm.S + R(w).Stot());
            frame(0, 0, 1.25 * hi, "peak content S [rows]", "Peak content S: bar = statistical, hatched = stat #oplus model spread");
            TGraphErrors* g = gr(kColS, 20, 1.8);
            for (int i = 0; i < n; ++i) {
                const Res& x = R(sel[i]);
                TBox* b = new TBox(i + 0.2, x.pm.S - x.Stot(), i + 0.8, x.pm.S + x.Stot());
                b->SetFillColor(kColS); b->SetFillStyle(3345); b->SetLineColor(kColS); b->Draw();
                add(g, i + 0.5, x.pm.S, x.pm.eS);
                TxtU(i + 0.5, x.pm.S + x.Stot() + 0.03 * hi, Form("%.0f", x.pm.S), 22, kColS, 21);
            }
            g->Draw("pz same");
            gPad->RedrawAxis();
        }
        // (1) F_stat
        {
            frame(1, 8, 18, "S/#sigma_{S} (stat. only)", "Peak: S/#sigma_{S}, stat. only: NOT a significance");
            TGraphErrors* g = gr(kColS, 20, 1.8);
            HLine(0, n, R0.pm.F, kColS, 2, 1);
            for (int i = 0; i < n; ++i) {
                const Res& x = R(sel[i]);
                const double e = std::sqrt(std::max(0., 2 * (1 - x.pm.eS / R0.pm.eS)));
                add(g, i + 0.5, x.pm.F, e);
                TxtU(i + 0.5, x.pm.F + e + 0.35, Form("%.1f", x.pm.F), 22, kColS, 21);
            }
            g->Draw("pz same");
            Txt(0.20, 0.19, "bars: error of the difference to 'none' (nested samples)", 20, kGray + 3);
            gPad->RedrawAxis();
        }
        // (2) S/B and CB share, peak window
        {
            frame(2, 0, 0.21, "ratio", "Peak window: S/B (blue) and CB share CB/N (red; no k factor)");
            TGraphErrors* g = gr(kColS, 20, 1.8); TGraphErrors* gc = gr(kColCB, 21, 1.6);
            for (int i = 0; i < n; ++i) {
                const Res& x = R(sel[i]);
                add(g, i + 0.5, x.pm.SoB, x.pm.SoB * x.pm.eS / x.pm.S);
                add(gc, i + 0.5, x.pm.CB / x.pm.N, 0);
                const bool sTop = x.pm.SoB > x.pm.CB / x.pm.N;
                TxtU(i + 0.5, x.pm.SoB + (sTop ? 0.014 : -0.014), Form("%.3f", x.pm.SoB), 21, kColS, sTop ? 21 : 23);
                TxtU(i + 0.5, x.pm.CB / x.pm.N + (sTop ? -0.008 : 0.008), Form("%.3f", x.pm.CB / x.pm.N), 21, kColCB, sTop ? 23 : 21);
            }
            gc->Draw("pz same"); g->Draw("pz same");
            gPad->RedrawAxis();
        }
        // (3) bump window: W kept and CB share
        {
            frame(3, 0, 1.25, "fraction", "Bump window: rows N kept (blue), CB share (red; no k factor)");
            TGraphErrors* g = gr(kColS, 20, 1.8); TGraphErrors* gc = gr(kColCB, 21, 1.6);
            for (int i = 0; i < n; ++i) {
                const Res& x = R(sel[i]);
                add(g, i + 0.5, x.bN / R0.bN, 0); add(gc, i + 0.5, x.bCB / x.bN, 0);
                TxtU(i + 0.5, x.bN / R0.bN + 0.05, Form("%.2f", x.bN / R0.bN), 21, kColS, 21);
                TxtU(i + 0.5, x.bCB / x.bN + 0.05, Form("%.3f", x.bCB / x.bN), 21, kColCB, 21);
            }
            gc->Draw("pz same"); g->Draw("pz same");
            gPad->RedrawAxis();
        }
        // (4) kept fractions of E and T
        {
            frame(4, 0.4, 1.25, "kept fraction", "Kept: step sum E (blue; not a yield), tagged #eta T (green)");
            TGraphErrors* g = gr(kColE, 21, 1.6); TGraphErrors* go = gr(kGray + 1, 25, 1.6); TGraphErrors* gt = gr(kColT, 20, 1.8);
            HLine(0, n, 1, kGray + 2, 3, 1);
            for (int i = 0; i < n; ++i) {
                const Res& x = R(sel[i]);
                const Eff e = EffE(x, R0), t = EffT(x, R0);
                add(x.Evalid() ? g : go, i + 0.5 - 0.12, e.v, i == 0 ? 0 : e.e);
                if (!On(sel[i].v[cPi0])) add(gt, i + 0.5 + 0.12, t.v, i == 0 ? 0 : t.e);
                else TxtU(i + 0.5 + 0.14, 0.66, "T n/a", 20, kColT, 21);
            }
            go->Draw("pz same"); g->Draw("pz same"); gt->Draw("pz same");
            Txt(0.20, 0.19, "open grey: edge fits not valid (#chi^{2});  T n/a: the tagging photon is vetoed", 20, kGray + 3);
            gPad->RedrawAxis();
        }
        // (5) gain
        {
            frame(5, 0.55, 1.45, "gain over 'none'", "Bump visibility gain: (kept fraction) / #sqrt{rows kept}");
            TGraphErrors* g = gr(kColE, 21, 1.6); TGraphErrors* go = gr(kGray + 1, 25, 1.6); TGraphErrors* gt = gr(kColT, 20, 1.8);
            HLine(0, n, 1, kGray + 2, 3, 1);
            for (int i = 0; i < n; ++i) {
                const Res& x = R(sel[i]);
                const Eff e = EffE(x, R0), t = EffT(x, R0);
                const double rw = std::sqrt(x.bN / R0.bN);
                add(x.Evalid() ? g : go, i + 0.5 - 0.12, e.v / rw, i == 0 ? 0 : e.e / rw);
                if (x.Evalid()) TxtU(i + 0.5 - 0.12, e.v / rw - (i == 0 ? 0 : e.e / rw) - 0.03, Form("%.2f", e.v / rw), 21, kColE, 23);
                if (!On(sel[i].v[cPi0])) {
                    add(gt, i + 0.5 + 0.12, t.v / rw, i == 0 ? 0 : t.e / rw);
                    TxtU(i + 0.5 + 0.12, t.v / rw + (i == 0 ? 0 : t.e / rw) + 0.03, Form("%.2f", t.v / rw), 21, kColT, 21);
                } else TxtU(i + 0.5 + 0.14, 1.22, "T n/a", 20, kColT, 21);
            }
            go->Draw("pz same"); g->Draw("pz same"); gt->Draw("pz same");
            Txt(0.20, 0.20, "blue: from E; green: from T; bars: statistical, relative to 'none'", 20, kGray + 3);
            gPad->RedrawAxis();
        }
        Save(c, "ribbons_summary");
    }
    fout->Close();
    delete fout;
    printf("scan_ribbons: done in %.1f s\n", swAll.RealTime());
}
