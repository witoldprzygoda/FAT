// ============================================================================
// ScanEngine.h  --  shared basis of the eta cut-scan study in M(pi+pi-e+e-)
// ============================================================================
// HADES p+p, T = 4.53 GeV, channel pi+ pi- e+ e-.
// Header-only; used by plots/scan_prepare.C (writes the compact file) and by
// all later plots/scan_*.C macros (read it, fill spectra, fit).
//
// Contents
//   A. frozen constants (binning, windows, fit range, peak shape, mass scale)
//   B. Row (one stored particle combination) and the loader of
//      scan_compact.root (trees us, pp, mm)
//   C. Spectra: M4 histograms of the three samples for an arbitrary selection,
//      CB = 2 sqrt(N++ N--) in 10 MeV bins interpolated to 2.5 MeV,
//      window sums with the per-event-group variance
//   D. fits: FitPeakMain (P_main), FitPeakPois (P_pois), PeakModelSpread,
//      FitEdgeStep (E_step)
//
// Conventions
//   - masses in GeV/c2, momenta in MeV/c, angles in deg, lengths in mm
//   - all rows are particle combinations (isBest >= 0); one event group can
//     give several rows (grp = running index of the group inside its sample)
//   - main mass scale = energy-loss corrected momenta; rows with cor_ok == 0
//     have NaN in every variable containing the pi- (use the *_fb columns or
//     Data::SetFallback(true) for the variant with the reconstructed pi-)
//   - like-sign samples: "ep" is the first and "em" the second lepton of the
//     same charge; use lepton quantities only symmetrically
//
// Compile macros that include this header with ACLiC and keep the build
// directory outside the repository, e.g.
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/.../scan/aclic_x", kTRUE)' \
//        'plots/scan_reference.C+'
// ============================================================================
#ifndef SCAN_ENGINE_H
#define SCAN_ENGINE_H

#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

namespace scan {

// ============================================================================
// A. Frozen constants
// ============================================================================

// ---- particle masses [GeV/c2] (same values as src/pparticle.h) -------------
const double kMeta  = 0.547862;
const double kMpi0  = 0.1349768;
const double kMpi   = 0.13957;      // 139.56995 MeV in FAT
const double kMK0   = 0.497611;
const double kMe    = 0.00051099895;
const double kMp    = 0.93827231;
const double kTbeam = 4530.0;       // beam kinetic energy [MeV]

// ---- files -----------------------------------------------------------------
const char* const kProdDir    = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const char* const kScanDir    = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan";
const char* const kCompact    = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/scan_compact.root";
const char* const kFigDir     = "/home/przygoda/HADES/FAT/FAT/plots/output/eta_scan";

// ---- M4 binning: 2.5 MeV from 0.35 to 0.95; CB in 10 MeV bins ---------------
const int    kNbin  = 240;
const double kMlo   = 0.35;
const double kMhi   = 0.95;
const double kBinW  = 0.0025;
const int    kCBreb = 4;            // 4 x 2.5 MeV = 10 MeV bins for the CB

// ---- frozen peak shape -------------------------------------------------------
// Free fit (P_main, position and width free, pol2, range kFitLo..kFitHi) on the
// reference selection oa_lab < 50, cor_ok == 1, isBest >= 0:
//   position 0.54388 +- 0.00049 GeV/c2, width 5.44 +- 0.60 MeV
// (Poisson fit of the unsubtracted spectrum: 0.54361 +- 0.00049, 5.50 +- 0.59).
// The width depends on the lower fit bound: 5.0 MeV (0.5150), 5.2 (0.5200),
// 5.4 (0.5225), 5.6 (0.5250); see reference_numbers.txt.
const double kPeakMean  = 0.5439;   // GeV/c2
const double kPeakSigma = 0.0054;   // GeV/c2

// ---- frozen fit range and background order of P_main / P_pois ---------------
// Lower bound: the highest kinematic edge of the bump (M(pipi) = 2 m_pi) is at
// 0.5241 true = 0.5202 on the measured scale; chi2 of the fit stops falling
// when the bound reaches 0.5225 (oa_lab < 50: 46.3/36, 31.1/34, 25.4/32,
// 20.9/31, 20.6/30 for 0.510, 0.515, 0.520, 0.5225, 0.525). The peak window
// starts at 0.5336, so the lower side band is only 11 MeV (2 sigma) wide.
// Consequence: S and the free width depend on this bound (S 4114, 4205, 4319,
// 4359 and width 5.0, 5.2, 5.4, 5.6 MeV for 0.515, 0.520, 0.5225, 0.525).
const double kFitLo    = 0.5225;
const double kFitHi    = 0.6100;
const int    kFitOrder = 2;

// ---- windows ----------------------------------------------------------------
const double kWinNsig  = 1.9;                                  // peak window half-width in sigma
const double kPeakLo   = kPeakMean - kWinNsig * kPeakSigma;    // 0.53364
const double kPeakHi   = kPeakMean + kWinNsig * kPeakSigma;    // 0.55416
// bump window: from the old first-look value to the start of the peak fit range
// (= highest edge 0.5202 rounded up to the bin grid)
const double kBumpLo   = 0.4550;
const double kBumpHi   = 0.5225;

// ---- mass scale of the corrected momenta: measured = scale * true ------------
// Four-particle mass, from the eta peak: 0.5439 / 0.547862 = 0.9928 (+- 0.0009).
// Pion pair, from K0s -> pi+pi- with a displaced pion crossing point
// (dz_pipi > 20 mm, doca_pipi < 10 mm): 0.49612 +- 0.00020, i.e. 0.9970 +- 0.0004.
// The four-particle mass is 0.4 % lower than the pion pair alone (leptons).
const double kScale     = kPeakMean / kMeta;
const double kScalePiPi = 0.9970;
// resolution used for the smeared kinematic edge (= peak width)
const double kEdgeSigma = kPeakSigma;

// ---- model variations for the systematic spread of S ------------------------
const int    kNVarRange = 4;
const double kVarLo[kNVarRange] = {0.520, 0.515, 0.510, 0.505};
const double kVarHi[kNVarRange] = {0.590, 0.610, 0.620, 0.650};
const int    kVarOrdMin = 1;
const int    kVarOrdMax = 4;
const double kVarDSigma = 0.0005;   // width +- 0.5 MeV
const double kVarPmin   = 0.05;     // keep fits with p > 0.05

// ---- kinematics of eta -> pi+ pi- pi0 ---------------------------------------
// kinetic energy of the pi0 in the eta rest frame for a given M(pi+pi-)
inline double TPi0(double mpipi) {
    return (kMeta * kMeta + kMpi0 * kMpi0 - mpipi * mpipi) / (2.0 * kMeta) - kMpi0;
}
// upper limit of M(pi+pi-e+e-) when the photon of pi0 -> e+e-gamma is lost
// (true masses, m_ee -> 0); -1 above the kinematic limit M(pipi) = m_eta - m_pi0
inline double M4max(double mpipi) {
    const double E0 = (kMeta * kMeta + kMpi0 * kMpi0 - mpipi * mpipi) / (2.0 * kMeta);
    const double p2 = E0 * E0 - kMpi0 * kMpi0;
    if (!(p2 >= 0.0)) return -1.0;
    const double m2 = kMeta * kMeta - kMeta * (E0 - std::sqrt(p2));
    return m2 > 0.0 ? std::sqrt(m2) : -1.0;
}
// lower limit (photon takes the largest possible momentum along the pi0)
inline double M4min(double mpipi) {
    const double E0 = (kMeta * kMeta + kMpi0 * kMpi0 - mpipi * mpipi) / (2.0 * kMeta);
    const double p2 = E0 * E0 - kMpi0 * kMpi0;
    if (!(p2 >= 0.0)) return -1.0;
    const double m2 = kMeta * kMeta - kMeta * (E0 + std::sqrt(p2));
    return m2 > 0.0 ? std::sqrt(m2) : -1.0;
}
// edge on the MEASURED (corrected) mass scale for a MEASURED M(pi+pi-):
// M(pi+pi-) is brought to the true scale with scalePiPi, the limit back to the
// measured four-particle scale with scale
inline double EdgeMeasured(double mpipi_meas, double scale = kScale, double scalePiPi = kScalePiPi) {
    const double e = M4max(mpipi_meas / scalePiPi);
    return e > 0.0 ? scale * e : -1.0;
}

// ============================================================================
// B. Row and loader
// ============================================================================
// The two lists define the columns of scan_compact.root AND the members of Row.

#define SCAN_INT_FIELDS(X)                                                        \
    X(file_id) X(pat_entry) X(run) X(grp) X(fold) X(isBest) X(nrow)               \
    X(n_pip) X(n_pim) X(n_lep1) X(n_lep2) X(n_lep) X(cor_ok)                      \
    X(ep_system) X(em_system) X(pip_system) X(pim_system)                         \
    X(ep_isUsedVertex) X(em_isUsedVertex) X(pip_isUsedVertex) X(pim_isUsedVertex) \
    X(ep_kIsLepton) X(em_kIsLepton) X(shared_ring)                                \
    X(pip_isring) X(pim_isring) X(evtStart)                                       \
    X(n_gam) X(g_idx) X(neutr_mult) X(fwdet_mult)

#define SCAN_FLT_FIELDS(X)                                                        \
    X(m4) X(mee) X(mpipi) X(oa) X(oa_lab) X(oa_eta) X(mm4) X(mm_pipi) X(pt4) X(y4) \
    X(m4_fb) X(mee_fb) X(mpipi_fb) X(oa_fb) X(oa_lab_fb) X(oa_eta_fb)             \
    X(mm4_fb) X(mm_pipi_fb) X(pt4_fb) X(y4_fb)                                    \
    X(mm_ee)                                                                      \
    X(vx) X(vy) X(vz) X(vchi2)                                                    \
    X(ep_p) X(em_p) X(pip_p) X(pim_p)                                             \
    X(ep_p_rec) X(em_p_rec) X(pip_p_rec) X(pim_p_rec)                             \
    X(ep_theta) X(em_theta) X(pip_theta) X(pim_theta)                             \
    X(ep_phi) X(em_phi) X(pip_phi) X(pim_phi)                                     \
    X(ep_rkchi2) X(em_rkchi2) X(pip_rkchi2) X(pim_rkchi2)                         \
    X(ep_richmatchquality) X(em_richmatchquality) X(ep_ringchi2) X(em_ringchi2)   \
    X(ep_oa_lept) X(em_oa_lept) X(ep_oa_hadr) X(em_oa_hadr)                       \
    X(pip_mass2) X(pim_mass2) X(pip_dedx_mdc) X(pim_dedx_mdc)                     \
    X(hyp_chi2)                                                                   \
    X(dz_ee) X(r_ee) X(rv_ee) X(doca_ee)                                          \
    X(dz_pipi) X(r_pipi) X(rv_pipi) X(doca_pipi)                                  \
    X(g_e) X(g_theta) X(g_phi) X(m_eeg) X(m5) X(mm5) X(m5_fb) X(mm5_fb)           \
    X(T_pi0) X(m4max)                                                             \
    X(fw1_p) X(fw1_theta)

// One stored particle combination.
//   keys      file_id, pat_entry (row in the PAT file), run, grp (event group),
//             fold (grp parity), isBest, nrow (rows in the group), n_pip, n_pim,
//             n_lep1, n_lep2 (distinct tracks in slot pi+, pi-, lepton 1,
//             lepton 2 of the group), n_lep (distinct lepton tracks, both slots)
//   masses    m4 = M(pi+pi-e+e-), mee, mpipi, mm4 = missing mass of the four,
//             mm_pipi, mm_ee (missing mass of e+e-); oa = e+e- opening angle,
//             oa_lab = angle between (pi+pi-) and (e+e-) in the lab, oa_eta =
//             the same in the rest frame forced to m_eta; pt4 [MeV/c], y4 (cms)
//   *_fb      the same with the reconstructed pi- momentum when cor_ok == 0
//   tracks    X_p (corrected; NaN for the pi- when cor_ok == 0), X_p_rec,
//             X_theta, X_phi, X_system, X_isUsedVertex, X_rkchi2
//   crossing  point of closest approach of the two lepton (pion) straight
//             lines: dz = z of the point - vz, r = distance from the z axis,
//             rv = transverse distance from (vx, vy), doca = distance of the
//             two lines; -9999 / -1 / -1 when the lines are parallel
//   photon    good photon (of g1..g3) with M(e+e-gamma) closest to 0.140:
//             g_idx (1..3), g_e [MeV], g_theta, g_phi, m_eeg, m5, mm5; -1 if none
//   derived   T_pi0 = TPi0(mpipi), m4max = M4max(mpipi) (true-mass formula
//             evaluated at the measured mpipi; -1 above the limit)
struct Row {
#define X(n) Int_t n;
    SCAN_INT_FIELDS(X)
#undef X
#define X(n) Float_t n;
    SCAN_FLT_FIELDS(X)
#undef X
};

// create (write = true) or connect (write = false) the branches of a compact tree
inline void BindRow(TTree* t, Row& r, bool write) {
#define X(n) if (write) t->Branch(#n, &r.n, #n "/I"); else t->SetBranchAddress(#n, &r.n);
    SCAN_INT_FIELDS(X)
#undef X
#define X(n) if (write) t->Branch(#n, &r.n, #n "/F"); else t->SetBranchAddress(#n, &r.n);
    SCAN_FLT_FIELDS(X)
#undef X
}

enum Sample { kUS = 0, kPP = 1, kMM = 2 };
inline const char* SampleName(int s) { return s == 0 ? "us" : (s == 1 ? "pp" : "mm"); }

// The three samples in memory.
struct Data {
    std::vector<Row> s[3];          // us, pp, mm
    bool fallback = false;          // true: rows with cor_ok == 0 carry the *_fb values
    std::vector<Row>& us() { return s[0]; }
    std::vector<Row>& pp() { return s[1]; }
    std::vector<Row>& mm() { return s[2]; }

    // Variant with the reconstructed pi- momentum for cor_ok == 0 rows:
    // on = true copies the *_fb values into m4, mee, mpipi, ... of these rows,
    // on = false restores NaN. cor_ok itself is never changed.
    void SetFallback(bool on) {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        for (int k = 0; k < 3; ++k) {
            for (auto& r : s[k]) {
                if (r.cor_ok) continue;
                if (on) {
                    r.m4 = r.m4_fb; r.mpipi = r.mpipi_fb; r.oa_lab = r.oa_lab_fb;
                    r.oa_eta = r.oa_eta_fb; r.mm4 = r.mm4_fb; r.mm_pipi = r.mm_pipi_fb;
                    r.pt4 = r.pt4_fb; r.y4 = r.y4_fb; r.m5 = r.m5_fb; r.mm5 = r.mm5_fb;
                    r.T_pi0 = (float)TPi0(r.mpipi); r.m4max = (float)M4max(r.mpipi);
                } else {
                    r.m4 = r.mpipi = r.oa_lab = r.oa_eta = r.mm4 = r.mm_pipi = nan;
                    r.pt4 = r.y4 = r.T_pi0 = r.m4max = nan;
                    r.m5 = (r.g_idx > 0) ? nan : -1.f;
                    r.mm5 = (r.g_idx > 0) ? nan : -1.f;
                }
            }
        }
        fallback = on;
    }
};

// Load the compact file into memory (about 2.6 GB, 15-25 s).
//   mlo, mhi : keep only rows with mlo <= m4_fb < mhi (default: everything)
//   samples  : bit mask, 1 = us, 2 = pp, 4 = mm
inline bool Load(Data& d, const char* path = kCompact,
                 double mlo = -1e9, double mhi = 1e9, int samples = 7) {
    TStopwatch sw;
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { printf("scan::Load: cannot open %s\n", path); return false; }
    for (int k = 0; k < 3; ++k) {
        d.s[k].clear();
        if (!(samples & (1 << k))) continue;
        TTree* t = (TTree*)f->Get(SampleName(k));
        if (!t) { printf("scan::Load: tree %s missing\n", SampleName(k)); return false; }
        Row r;
        BindRow(t, r, false);
        const Long64_t n = t->GetEntries();
        d.s[k].reserve(n);
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            if (!(r.m4_fb >= mlo && r.m4_fb < mhi)) continue;
            d.s[k].push_back(r);
        }
        t->ResetBranchAddresses();
    }
    f->Close();
    delete f;
    d.fallback = false;
    printf("scan::Load: us %zu, pp %zu, mm %zu rows in %.1f s\n",
           d.s[0].size(), d.s[1].size(), d.s[2].size(), sw.RealTime());
    return true;
}

// ============================================================================
// C. Spectra for an arbitrary selection
// ============================================================================

typedef std::function<bool(const Row&)> Sel;
typedef std::function<double(const Row&)> Var;

// ---- the reference selections ------------------------------------------------
inline Sel SelNone()      { return [](const Row&) { return true; }; }
inline Sel SelOaLab50()   { return [](const Row& r) { return r.oa_lab < 50.f; }; }
// "slides": oa_lab < 50, M(pipi) < 0.46, eta_rest > 140, OA > 2
inline Sel SelSlides()    { return [](const Row& r) {
    return r.oa_lab < 50.f && r.mpipi < 0.46f && r.oa_eta > 140.f && r.oa > 2.f; }; }
// "old code": oa_lab < 50, M(pipi) < 0.42, eta_rest > 140
inline Sel SelOldCode()   { return [](const Row& r) {
    return r.oa_lab < 50.f && r.mpipi < 0.42f && r.oa_eta > 140.f; }; }
// "old code + OA 2-6": the same with 2 < OA < 6 (flags pass_oa_min, pass_oa_max)
inline Sel SelOldCodeOA() { return [](const Row& r) {
    return r.oa_lab < 50.f && r.mpipi < 0.42f && r.oa_eta > 140.f && r.oa > 2.f && r.oa < 6.f; }; }
// add "isBest == 1" to a selection (only for explicit comparisons)
inline Sel SelBest(const Sel& s) { return [s](const Row& r) { return r.isBest == 1 && s(r); }; }
inline Sel SelAnd(const Sel& a, const Sel& b) { return [a, b](const Row& r) { return a(r) && b(r); }; }

// Sum in a mass window [lo, hi). n = rows, var = sum over event groups of k^2
// (k = selected rows of the group inside the window); index 0 us, 1 pp, 2 mm.
// cb = integral of the interpolated CB over the window, cbVar = its variance
// (linear propagation of the like-sign counts, per event group).
struct WinSum {
    double lo = 0, hi = 0;
    double n[3] = {0, 0, 0};
    double var[3] = {0, 0, 0};
    double cb = 0, cbVar = 0;
    double sig() const { return n[0] - cb; }                       // N - CB
    double sigErr() const { return std::sqrt(var[0] + cbVar); }
    double cbShare() const { return n[0] > 0 ? cb / n[0] : 0.; }
};

struct FillOpt {
    int    nbin = kNbin;            // fine bins
    double lo = kMlo, hi = kMhi;
    int    reb = kCBreb;            // fine bins per CB bin
    Var    x;                       // histogrammed variable; empty = Row::m4
    std::vector<std::pair<double, double>> windows;   // extra windows (after peak and bump)
};

// Histograms of one selection. Fine bins of width (hi-lo)/nbin.
struct Spectra {
    int nbin = 0, reb = 1, nw = 0;
    double lo = 0, hi = 0, bw = 0;
    std::vector<double> n[3];            // rows per fine bin: us, pp, mm
    std::vector<double> w[3];            // rows per CB (wide) bin
    std::vector<double> cbW, cbWVar;     // CB = 2 sqrt(N++ N--) per wide bin, var = N++ + N--
    std::vector<double> cb, cbVar;       // interpolated to fine bins
    std::vector<double> sig, sigVar;     // N - CB and its variance N + cbVar
    std::vector<WinSum> win;             // win[0] = peak window, win[1] = bump window, then extras
    double nsel[3] = {0, 0, 0};          // selected rows at any mass (cor_ok or fallback)
    std::vector<std::pair<int, int>> ls[2];   // like-sign rows: (grp, wide bin), for the CB variance

    double x(int i) const { return lo + (i + 0.5) * bw; }           // centre of fine bin i
    int bin(double m) const { return (m >= lo && m < hi) ? (int)((m - lo) / bw) : -1; }
    // interpolation of the wide-bin CB to the centre of fine bin i: indices and weight
    void Interp(int i, int& j0, int& j1, double& f) const {
        const double t = (x(i) - (lo + 0.5 * reb * bw)) / (reb * bw);
        if (t <= 0) { j0 = j1 = 0; f = 0; return; }
        if (t >= nw - 1) { j0 = j1 = nw - 1; f = 0; return; }
        j0 = (int)t; j1 = j0 + 1; f = t - j0;
    }
    // covariance of the interpolated CB in the fine bins i and k (they share wide bins)
    double CbCov(int i, int k) const {
        int a0, a1, b0, b1; double fa, fb;
        Interp(i, a0, a1, fa); Interp(k, b0, b1, fb);
        const int ja[2] = {a0, a1}, jb[2] = {b0, b1};
        const double wa[2] = {(1 - fa) / reb, fa / reb}, wb[2] = {(1 - fb) / reb, fb / reb};
        double s = 0;
        for (int p = 0; p < 2; ++p)
            for (int q = 0; q < 2; ++q)
                if (ja[p] == jb[q]) s += wa[p] * wb[q] * cbWVar[ja[p]];
        return s;
    }
    // fraction of fine bin i inside [a, b)
    double Frac(int i, double a, double b) const {
        const double x1 = lo + i * bw, x2 = x1 + bw;
        const double o = std::min(x2, b) - std::max(x1, a);
        return o > 0 ? o / bw : 0.;
    }
    // histogram of one component: "us", "pp", "mm", "cb", "sig"
    TH1D* Hist(const char* what, const char* name) const {
        TH1D* h = new TH1D(name, "", nbin, lo, hi);
        h->SetDirectory(nullptr);
        const std::string k = what;
        for (int i = 0; i < nbin; ++i) {
            double v = 0, e = 0;
            if (k == "us") { v = n[0][i]; e = std::sqrt(v); }
            else if (k == "pp") { v = n[1][i]; e = std::sqrt(v); }
            else if (k == "mm") { v = n[2][i]; e = std::sqrt(v); }
            else if (k == "cb") { v = cb[i]; e = std::sqrt(cbVar[i]); }
            else { v = sig[i]; e = std::sqrt(sigVar[i]); }
            h->SetBinContent(i + 1, v);
            h->SetBinError(i + 1, e);
        }
        return h;
    }
    // sum over a window that was not requested at fill time: fractional fine
    // bins, Poisson variance (no event-group correction)
    WinSum WindowApprox(double a, double b) const {
        WinSum s; s.lo = a; s.hi = b;
        std::vector<double> c(nw, 0.);
        for (int i = 0; i < nbin; ++i) {
            const double u = Frac(i, a, b);
            if (u <= 0) continue;
            for (int k = 0; k < 3; ++k) { s.n[k] += u * n[k][i]; s.var[k] += u * n[k][i]; }
            int j0, j1; double f; Interp(i, j0, j1, f);
            c[j0] += u * (1 - f) / reb; c[j1] += u * f / reb;
        }
        for (int j = 0; j < nw; ++j) { s.cb += c[j] * cbW[j]; s.cbVar += c[j] * c[j] * cbWVar[j]; }
        return s;
    }
};

// Fill the spectra of several selections in ONE pass over the rows.
// A row is used when cor_ok == 1 (or always after Data::SetFallback(true)),
// the selection accepts it and the mass is a number. Windows: 0 = peak
// (kPeakLo..kPeakHi), 1 = bump (kBumpLo..kBumpHi), then opt.windows.
inline std::vector<Spectra> FillMany(const Data& d, const std::vector<Sel>& sels,
                                     const FillOpt& opt = FillOpt()) {
    const int ns = (int)sels.size();
    std::vector<std::pair<double, double>> wl;
    wl.push_back(std::make_pair(kPeakLo, kPeakHi));
    wl.push_back(std::make_pair(kBumpLo, kBumpHi));
    for (const auto& p : opt.windows) wl.push_back(p);
    const int nwin = (int)wl.size();

    std::vector<Spectra> out(ns);
    for (auto& sp : out) {
        sp.nbin = opt.nbin; sp.reb = opt.reb; sp.nw = opt.nbin / opt.reb;
        sp.lo = opt.lo; sp.hi = opt.hi; sp.bw = (opt.hi - opt.lo) / opt.nbin;
        for (int k = 0; k < 3; ++k) { sp.n[k].assign(sp.nbin, 0.); sp.w[k].assign(sp.nw, 0.); }
        sp.win.resize(nwin);
        for (int w = 0; w < nwin; ++w) { sp.win[w].lo = wl[w].first; sp.win[w].hi = wl[w].second; }
    }
    const double lo = opt.lo, hi = opt.hi, bw = (opt.hi - opt.lo) / opt.nbin;
    const bool useX = (bool)opt.x;

    std::vector<int> kw(ns * nwin, 0);
    for (int k = 0; k < 3; ++k) {
        const std::vector<Row>& v = d.s[k];
        int cur = -1;
        bool dirty = false;
        for (size_t i = 0; i < v.size(); ++i) {
            const Row& r = v[i];
            if (!r.cor_ok && !d.fallback) continue;
            if (r.grp != cur) {
                if (dirty) {
                    for (int s = 0; s < ns; ++s)
                        for (int w = 0; w < nwin; ++w) {
                            int& c = kw[s * nwin + w];
                            if (c) { out[s].win[w].var[k] += (double)c * c; c = 0; }
                        }
                    dirty = false;
                }
                cur = r.grp;
            }
            const double m = useX ? opt.x(r) : (double)r.m4;
            if (!(m == m)) continue;
            const int ib = (m >= lo && m < hi) ? (int)((m - lo) / bw) : -1;
            for (int s = 0; s < ns; ++s) {
                if (!sels[s](r)) continue;
                Spectra& sp = out[s];
                sp.nsel[k] += 1;
                if (ib >= 0 && ib < sp.nbin) {
                    sp.n[k][ib] += 1;
                    if (k > 0) sp.ls[k - 1].push_back(std::make_pair(r.grp, ib / sp.reb));
                }
                for (int w = 0; w < nwin; ++w) {
                    if (m >= wl[w].first && m < wl[w].second) {
                        sp.win[w].n[k] += 1;
                        kw[s * nwin + w]++;
                        dirty = true;
                    }
                }
            }
        }
        if (dirty) {
            for (int s = 0; s < ns; ++s)
                for (int w = 0; w < nwin; ++w) {
                    int& c = kw[s * nwin + w];
                    if (c) { out[s].win[w].var[k] += (double)c * c; c = 0; }
                }
        }
    }

    // ---- combinatorial background ----
    for (auto& sp : out) {
        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < sp.nbin; ++i) sp.w[k][i / sp.reb] += sp.n[k][i];
        sp.cbW.assign(sp.nw, 0.); sp.cbWVar.assign(sp.nw, 0.);
        for (int j = 0; j < sp.nw; ++j) {
            sp.cbW[j] = 2.0 * std::sqrt(sp.w[1][j] * sp.w[2][j]);
            sp.cbWVar[j] = sp.w[1][j] + sp.w[2][j];
        }
        sp.cb.assign(sp.nbin, 0.); sp.cbVar.assign(sp.nbin, 0.);
        sp.sig.assign(sp.nbin, 0.); sp.sigVar.assign(sp.nbin, 0.);
        for (int i = 0; i < sp.nbin; ++i) {
            int j0, j1; double f; sp.Interp(i, j0, j1, f);
            const double a0 = (1 - f) / sp.reb, a1 = f / sp.reb;
            sp.cb[i] = a0 * sp.cbW[j0] + a1 * sp.cbW[j1];
            sp.cbVar[i] = a0 * a0 * sp.cbWVar[j0] + a1 * a1 * sp.cbWVar[j1];
            sp.sig[i] = sp.n[0][i] - sp.cb[i];
            sp.sigVar[i] = sp.n[0][i] + sp.cbVar[i];
        }
        // CB in the windows: sum_j c_j CB_j; variance per event group with the
        // weights dCB_j/dN++_j = sqrt(N--_j/N++_j) and dCB_j/dN--_j = sqrt(N++_j/N--_j)
        for (auto& ws : sp.win) {
            std::vector<double> c(sp.nw, 0.);
            for (int i = 0; i < sp.nbin; ++i) {
                const double u = sp.Frac(i, ws.lo, ws.hi);
                if (u <= 0) continue;
                int j0, j1; double f; sp.Interp(i, j0, j1, f);
                c[j0] += u * (1 - f) / sp.reb; c[j1] += u * f / sp.reb;
            }
            ws.cb = 0; ws.cbVar = 0;
            for (int j = 0; j < sp.nw; ++j) ws.cb += c[j] * sp.cbW[j];
            for (int q = 0; q < 2; ++q) {            // q = 0: pp rows, 1: mm rows
                std::vector<double> wt(sp.nw, 0.);
                for (int j = 0; j < sp.nw; ++j) {
                    const double own = sp.w[1 + q][j], oth = sp.w[2 - q][j];
                    wt[j] = c[j] * ((own > 0 && oth > 0) ? std::sqrt(oth / own) : 1.0);
                }
                int cur = -1; double acc = 0;
                for (const auto& e : sp.ls[q]) {
                    if (e.first != cur) { ws.cbVar += acc * acc; acc = 0; cur = e.first; }
                    acc += wt[e.second];
                }
                ws.cbVar += acc * acc;
            }
        }
    }
    return out;
}

// Spectra of one selection (see FillMany).
inline Spectra Fill(const Data& d, const Sel& sel, const FillOpt& opt = FillOpt()) {
    std::vector<Sel> v(1, sel);
    return FillMany(d, v, opt)[0];
}

// ============================================================================
// D. Fits
// ============================================================================

// ---- linear least squares, up to 10 parameters -------------------------------
const int kMaxPar = 10;
struct LinRes {
    bool ok = false;
    int np = 0, npts = 0;
    double par[kMaxPar];
    double cov[kMaxPar][kMaxPar];
    double chi2 = 0;
};
// invert a symmetric positive definite matrix in place (Cholesky)
inline bool CholInvert(int n, double a[kMaxPar][kMaxPar]) {
    double l[kMaxPar][kMaxPar];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = a[i][j];
            for (int k = 0; k < j; ++k) s -= l[i][k] * l[j][k];
            if (i == j) { if (!(s > 0)) return false; l[i][i] = std::sqrt(s); }
            else l[i][j] = s / l[j][j];
        }
    double li[kMaxPar][kMaxPar];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) li[i][j] = 0;
    for (int i = 0; i < n; ++i) {
        li[i][i] = 1.0 / l[i][i];
        for (int j = 0; j < i; ++j) {
            double s = 0;
            for (int k = j; k < i; ++k) s -= l[i][k] * li[k][j];
            li[i][j] = s / l[i][i];
        }
    }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = std::max(i, j); k < n; ++k) s += li[k][i] * li[k][j];
            a[i][j] = s;
        }
    return true;
}
// minimise sum (y_i - sum_k p_k F[i*np+k])^2 / var_i
inline LinRes LinFit(int npts, int np, const double* y, const double* var, const double* F) {
    LinRes r; r.np = np; r.npts = npts;
    for (int k = 0; k < kMaxPar; ++k) { r.par[k] = 0; for (int l = 0; l < kMaxPar; ++l) r.cov[k][l] = 0; }
    if (np > kMaxPar || npts <= np) return r;
    double b[kMaxPar];
    // columns are scaled to unit norm before the inversion (conditioning)
    double sc[kMaxPar];
    for (int k = 0; k < np; ++k) {
        double s = 0;
        for (int i = 0; i < npts; ++i) s += F[i * np + k] * F[i * np + k] / var[i];
        sc[k] = s > 0 ? 1.0 / std::sqrt(s) : 1.0;
    }
    for (int k = 0; k < np; ++k) {
        b[k] = 0;
        for (int l = 0; l < np; ++l) r.cov[k][l] = 0;
    }
    for (int i = 0; i < npts; ++i) {
        const double w = 1.0 / var[i];
        for (int k = 0; k < np; ++k) {
            const double fk = F[i * np + k] * sc[k] * w;
            b[k] += fk * y[i];
            for (int l = 0; l <= k; ++l) r.cov[k][l] += fk * F[i * np + l] * sc[l];
        }
    }
    for (int k = 0; k < np; ++k) for (int l = k + 1; l < np; ++l) r.cov[k][l] = r.cov[l][k];
    if (!CholInvert(np, r.cov)) return r;
    for (int k = 0; k < np; ++k) {
        double s = 0;
        for (int l = 0; l < np; ++l) s += r.cov[k][l] * b[l];
        r.par[k] = s * sc[k];
    }
    for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) r.cov[k][l] *= sc[k] * sc[l];
    r.chi2 = 0;
    for (int i = 0; i < npts; ++i) {
        double f = 0;
        for (int k = 0; k < np; ++k) f += r.par[k] * F[i * np + k];
        r.chi2 += (y[i] - f) * (y[i] - f) / var[i];
    }
    r.ok = true;
    return r;
}

// The same with a full covariance matrix V (npts*npts, row-major; destroyed):
// minimise (y - F p)^T V^-1 (y - F p). V = L L^T, then ordinary least squares
// on L^-1 y and L^-1 F.
inline LinRes LinFitCov(int npts, int np, const double* y, std::vector<double>& V, const double* F) {
    LinRes bad; bad.np = np; bad.npts = npts;
    for (int i = 0; i < npts; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = V[i * npts + j];
            for (int k = 0; k < j; ++k) s -= V[i * npts + k] * V[j * npts + k];
            if (i == j) { if (!(s > 0)) return bad; V[i * npts + i] = std::sqrt(s); }
            else V[i * npts + j] = s / V[j * npts + j];
        }
    std::vector<double> z(npts), G(npts * np), one(npts, 1.0);
    for (int i = 0; i < npts; ++i) {
        double s = y[i];
        for (int k = 0; k < i; ++k) s -= V[i * npts + k] * z[k];
        z[i] = s / V[i * npts + i];
        for (int c = 0; c < np; ++c) {
            double g = F[i * np + c];
            for (int k = 0; k < i; ++k) g -= V[i * npts + k] * G[k * np + c];
            G[i * np + c] = g / V[i * npts + i];
        }
    }
    return LinFit(npts, np, z.data(), one.data(), G.data());
}

// fraction of a Gaussian inside [x1, x2]
inline double GaussFrac(double x1, double x2, double mean, double sigma) {
    const double s = 1.0 / (std::sqrt(2.0) * sigma);
    return 0.5 * (std::erf((x2 - mean) * s) - std::erf((x1 - mean) * s));
}

// ---- eta peak: Gaussian with frozen position and width + polynomial -----------
struct PeakCfg {
    double mean = kPeakMean, sigma = kPeakSigma;   // frozen shape
    double lo = kFitLo, hi = kFitHi;               // fit range (bin centres inside)
    int order = kFitOrder;                         // polynomial order
    double winLo = kPeakLo, winHi = kPeakHi;       // window of N, CB, B
};
struct PeakFit {
    bool ok = false;
    int method = 0;              // 0 = P_main (chi2 on N - CB), 1 = P_pois (likelihood on N)
    PeakCfg cfg;
    double S = 0, eS = 0;        // peak content (whole Gaussian) and its statistical error
    double chi2 = 0; int ndf = 0; double prob = 0;
    double N = 0, eN = 0;        // unlike-sign rows in the window, error from the event groups
    double CB = 0, eCB = 0;      // CB in the window
    double B = 0;                // fitted non-peak, non-CB background in the window
    double Swin = 0;             // part of S inside the window
    double SoB = 0;              // Swin / (B + CB)
    double F = 0;                // S / eS
    int np = 0;
    double par[kMaxPar];         // par[0] = S, par[1+j] = coefficient of t^j, t = (m - x0)/dx
    double x0 = 0, dx = 1, bw = kBinW;
    // fitted background and total, in rows per fine bin, at mass m
    double Bkg(double m) const {
        const double t = (m - x0) / dx; double p = 0, tp = 1;
        for (int j = 0; j <= cfg.order; ++j) { p += par[1 + j] * tp; tp *= t; }
        return p;
    }
    double Total(double m) const {
        return Bkg(m) + par[0] * GaussFrac(m - 0.5 * bw, m + 0.5 * bw, cfg.mean, cfg.sigma);
    }
};

// Core of the peak fits on arrays over the fine bins of a spectrum.
//   y    fitted spectrum (N - CB for P_main, N for P_pois)
//   add  expected rows not described by the model (CB for P_main), may be null
//   cov  spectrum whose CB covariance is added to the Poisson term, may be null
// chi2 mode: covariance = diag(expected rows = model + add, at least 1)
//            + covariance of the interpolated CB (neighbouring fine bins share
//            the 10 MeV CB bins); iterated (first pass with the observed rows);
//            linear in all parameters (generalised least squares).
// poisson  : binned Poisson likelihood for y (Newton steps), chi2 = likelihood-ratio chi2.
inline PeakFit FitPeakArrays(const double* y, const double* add, const Spectra* cov,
                             int nbin, double lo, double bw, const PeakCfg& c, bool poisson) {
    PeakFit r; r.cfg = c; r.method = poisson ? 1 : 0; r.bw = bw;
    for (int k = 0; k < kMaxPar; ++k) r.par[k] = 0;
    const int np = c.order + 2;
    r.np = np;
    r.x0 = 0.5 * (c.lo + c.hi); r.dx = 0.5 * (c.hi - c.lo);
    std::vector<int> idx;
    for (int i = 0; i < nbin; ++i) {
        const double xc = lo + (i + 0.5) * bw;
        if (xc >= c.lo && xc < c.hi) idx.push_back(i);
    }
    const int npts = (int)idx.size();
    if (np > kMaxPar || npts <= np) return r;
    std::vector<double> F(npts * np), yy(npts), var(npts), ad(npts, 0.), V0, V;
    for (int q = 0; q < npts; ++q) {
        const int i = idx[q];
        const double x1 = lo + i * bw, xc = x1 + 0.5 * bw;
        F[q * np] = GaussFrac(x1, x1 + bw, c.mean, c.sigma);
        const double t = (xc - r.x0) / r.dx; double tp = 1;
        for (int j = 0; j <= c.order; ++j) { F[q * np + 1 + j] = tp; tp *= t; }
        yy[q] = y[i];
        if (add) ad[q] = add[i];
        var[q] = std::max(yy[q] + ad[q], 1.0);
    }
    if (cov && !poisson) {
        V0.assign(npts * npts, 0.);
        for (int q = 0; q < npts; ++q)
            for (int p = 0; p <= q; ++p) V0[q * npts + p] = V0[p * npts + q] = cov->CbCov(idx[q], idx[p]);
    }
    // one weighted linear fit with the current Poisson variances
    auto solve = [&]() -> LinRes {
        if (V0.empty()) return LinFit(npts, np, yy.data(), var.data(), F.data());
        V = V0;
        for (int q = 0; q < npts; ++q) V[q * npts + q] += var[q];
        return LinFitCov(npts, np, yy.data(), V, F.data());
    };
    LinRes lr = solve();
    if (!lr.ok) return r;
    if (!poisson) {
        for (int it = 0; it < 6; ++it) {
            const double sOld = lr.par[0];
            for (int q = 0; q < npts; ++q) {
                double f = 0;
                for (int k = 0; k < np; ++k) f += lr.par[k] * F[q * np + k];
                var[q] = std::max(f + ad[q], 1.0);
            }
            lr = solve();
            if (!lr.ok) return r;
            if (std::fabs(lr.par[0] - sOld) < 1e-4 * std::sqrt(lr.cov[0][0])) break;
        }
        r.chi2 = lr.chi2;
    } else {
        double p[kMaxPar], H[kMaxPar][kMaxPar];
        for (int k = 0; k < np; ++k) p[k] = lr.par[k];
        bool conv = false;
        for (int it = 0; it < 50 && !conv; ++it) {
            double g[kMaxPar];
            for (int k = 0; k < np; ++k) { g[k] = 0; for (int l = 0; l < np; ++l) H[k][l] = 0; }
            for (int q = 0; q < npts; ++q) {
                double mu = 0;
                for (int k = 0; k < np; ++k) mu += p[k] * F[q * np + k];
                if (mu < 1e-3) mu = 1e-3;
                const double a = 1.0 - yy[q] / mu;
                const double h = std::max(yy[q], 1e-3) / (mu * mu);
                for (int k = 0; k < np; ++k) {
                    g[k] += a * F[q * np + k];
                    for (int l = 0; l <= k; ++l) H[k][l] += h * F[q * np + k] * F[q * np + l];
                }
            }
            for (int k = 0; k < np; ++k) for (int l = k + 1; l < np; ++l) H[k][l] = H[l][k];
            // scale to unit diagonal before the inversion
            double sc[kMaxPar];
            for (int k = 0; k < np; ++k) sc[k] = H[k][k] > 0 ? 1.0 / std::sqrt(H[k][k]) : 1.0;
            for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) H[k][l] *= sc[k] * sc[l];
            if (!CholInvert(np, H)) return r;
            for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) H[k][l] *= sc[k] * sc[l];
            double step2 = 0;
            for (int k = 0; k < np; ++k) {
                double dlt = 0;
                for (int l = 0; l < np; ++l) dlt -= H[k][l] * g[l];
                p[k] += dlt;
                step2 += dlt * dlt / H[k][k];
            }
            if (step2 < 1e-10) conv = true;
        }
        if (!conv) return r;
        for (int k = 0; k < np; ++k) { lr.par[k] = p[k]; for (int l = 0; l < np; ++l) lr.cov[k][l] = H[k][l]; }
        double c2 = 0;
        for (int q = 0; q < npts; ++q) {
            double mu = 0;
            for (int k = 0; k < np; ++k) mu += p[k] * F[q * np + k];
            if (mu < 1e-3) mu = 1e-3;
            c2 += 2.0 * (mu - yy[q] + (yy[q] > 0 ? yy[q] * std::log(yy[q] / mu) : 0.0));
        }
        r.chi2 = c2;
    }
    for (int k = 0; k < np; ++k) r.par[k] = lr.par[k];
    r.S = lr.par[0];
    r.eS = std::sqrt(lr.cov[0][0]);
    r.ndf = npts - np;
    r.prob = TMath::Prob(r.chi2, r.ndf);
    r.F = r.eS > 0 ? r.S / r.eS : 0;
    r.Swin = r.S * GaussFrac(c.winLo, c.winHi, c.mean, c.sigma);
    // polynomial integrated over the window, in rows (midpoint rule)
    {
        const int nst = 400; const double h = (c.winHi - c.winLo) / nst; double s = 0;
        for (int q = 0; q < nst; ++q) s += r.Bkg(c.winLo + (q + 0.5) * h);
        r.B = s * h / bw;
    }
    r.ok = true;
    return r;
}

// window sum matching the fit configuration (exact if requested at fill time)
inline WinSum PeakWindow(const Spectra& sp, const PeakCfg& c) {
    for (const auto& w : sp.win)
        if (std::fabs(w.lo - c.winLo) < 1e-9 && std::fabs(w.hi - c.winHi) < 1e-9) return w;
    return sp.WindowApprox(c.winLo, c.winHi);
}

// P_main: chi2 fit of N - CB, Gaussian (frozen position and width) + polynomial
// in the frozen range. B = polynomial in the window, CB separately.
inline PeakFit FitPeakMain(const Spectra& sp, const PeakCfg& c = PeakCfg()) {
    PeakFit r = FitPeakArrays(sp.sig.data(), sp.cb.data(), &sp, sp.nbin, sp.lo, sp.bw, c, false);
    const WinSum w = PeakWindow(sp, c);
    r.N = w.n[0]; r.eN = std::sqrt(w.var[0]); r.CB = w.cb; r.eCB = std::sqrt(w.cbVar);
    r.SoB = (r.B + r.CB) > 0 ? r.Swin / (r.B + r.CB) : 0;
    return r;
}

// P_pois: Poisson likelihood fit of the unsubtracted N with the same shape;
// the CB stays inside the polynomial. B = polynomial in the window - CB.
inline PeakFit FitPeakPois(const Spectra& sp, const PeakCfg& c = PeakCfg()) {
    PeakFit r = FitPeakArrays(sp.n[0].data(), nullptr, nullptr, sp.nbin, sp.lo, sp.bw, c, true);
    const WinSum w = PeakWindow(sp, c);
    r.N = w.n[0]; r.eN = std::sqrt(w.var[0]); r.CB = w.cb; r.eCB = std::sqrt(w.cbVar);
    const double btot = r.B;
    r.B = btot - r.CB;
    r.SoB = btot > 0 ? r.Swin / btot : 0;
    return r;
}

// Null test: the same Poisson fit applied to a like-sign spectrum
// (k = 1: e+e+, k = 2: e-e-). N = like-sign rows in the window; B = polynomial.
inline PeakFit FitPeakLikeSign(const Spectra& sp, int k, const PeakCfg& c = PeakCfg()) {
    PeakFit r = FitPeakArrays(sp.n[k].data(), nullptr, nullptr, sp.nbin, sp.lo, sp.bw, c, true);
    const WinSum w = PeakWindow(sp, c);
    r.N = w.n[k]; r.eN = std::sqrt(w.var[k]); r.CB = 0; r.eCB = 0;
    r.SoB = r.B > 0 ? r.Swin / r.B : 0;
    return r;
}

// ---- free peak shape (used once, to determine the frozen numbers) -------------
// Position and width are varied; for each pair the remaining parameters are
// linear (P_main or P_pois). Newton steps on a 3x3 stencil of chi2(mean, sigma).
// Errors from the curvature of chi2 (delta chi2 = 1); eS includes the shape.
struct FreeFit {
    bool ok = false;
    double mean = 0, eMean = 0, sigma = 0, eSigma = 0, rho = 0;
    double S = 0, eS = 0, eSfixed = 0, chi2 = 0; int ndf = 0; double prob = 0;
};
inline FreeFit FitPeakFree(const Spectra& sp, PeakCfg c = PeakCfg(), bool poisson = false) {
    FreeFit r;
    const double hm = 0.0001, hs = 0.0001;
    auto fit = [&](double m, double s) -> PeakFit {
        PeakCfg q = c; q.mean = m; q.sigma = s;
        return poisson ? FitPeakPois(sp, q) : FitPeakMain(sp, q);
    };
    double m = c.mean, s = c.sigma;
    double cov[2][2] = {{0, 0}, {0, 0}}, dS[2] = {0, 0};
    for (int it = 0; it < 60; ++it) {
        double f[3][3]; double S[3][3];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                const PeakFit p = fit(m + (a - 1) * hm, s + (b - 1) * hs);
                if (!p.ok) return r;
                f[a][b] = p.chi2; S[a][b] = p.S;
            }
        const double gm = (f[2][1] - f[0][1]) / (2 * hm), gs = (f[1][2] - f[1][0]) / (2 * hs);
        const double hmm = (f[2][1] - 2 * f[1][1] + f[0][1]) / (hm * hm);
        const double hss = (f[1][2] - 2 * f[1][1] + f[1][0]) / (hs * hs);
        const double hms = (f[2][2] - f[2][0] - f[0][2] + f[0][0]) / (4 * hm * hs);
        const double det = hmm * hss - hms * hms;
        if (!(det > 0) || !(hmm > 0)) return r;
        double dm = -(hss * gm - hms * gs) / det, ds = -(-hms * gm + hmm * gs) / det;
        dm = std::max(-0.001, std::min(0.001, dm));
        ds = std::max(-0.001, std::min(0.001, ds));
        m += dm; s += ds;
        if (s < 0.002) s = 0.002;
        if (s > 0.012) s = 0.012;
        cov[0][0] = 2 * hss / det; cov[1][1] = 2 * hmm / det; cov[0][1] = cov[1][0] = -2 * hms / det;
        dS[0] = (S[2][1] - S[0][1]) / (2 * hm); dS[1] = (S[1][2] - S[1][0]) / (2 * hs);
        if (std::fabs(dm) < 2e-7 && std::fabs(ds) < 2e-7) { r.ok = true; break; }
    }
    if (!r.ok) return r;
    const PeakFit p = fit(m, s);
    r.mean = m; r.sigma = s;
    r.eMean = std::sqrt(cov[0][0]); r.eSigma = std::sqrt(cov[1][1]);
    r.rho = cov[0][1] / (r.eMean * r.eSigma);
    r.S = p.S; r.eSfixed = p.eS;
    r.eS = std::sqrt(p.eS * p.eS + dS[0] * dS[0] * cov[0][0] + dS[1] * dS[1] * cov[1][1] +
                     2 * dS[0] * dS[1] * cov[0][1]);
    r.chi2 = p.chi2; r.ndf = p.ndf - 2; r.prob = TMath::Prob(r.chi2, r.ndf);
    return r;
}

// ---- model variations: spread of S -------------------------------------------
// P_main repeated for polynomial order 1..4, the four ranges kVarLo/kVarHi and
// width -0.5, 0, +0.5 MeV (48 fits). Kept: fits with p > 0.05.
struct Spread {
    int nfit = 0, npass = 0;
    double smin = 0, smax = 0;   // over the kept fits
    double half = 0;             // (smax - smin) / 2  = the quoted model spread
    double rms = 0, mean = 0;
    double halfBkg = 0;          // the same for the 16 fits with the nominal width only
    double halfWidth = 0;        // |S(width + 0.5) - S(width - 0.5)| / 2 in the nominal range and order
    int npassBkg = 0;
    std::vector<PeakFit> fits;   // all 48 (when keep = true)
};
inline Spread PeakModelSpread(const Spectra& sp, const PeakCfg& nominal = PeakCfg(), bool keep = false) {
    Spread s;
    double sum = 0, sum2 = 0, bmin = 0, bmax = 0;
    {
        PeakCfg a = nominal, b = nominal;
        a.sigma = nominal.sigma - kVarDSigma; b.sigma = nominal.sigma + kVarDSigma;
        s.halfWidth = 0.5 * std::fabs(FitPeakMain(sp, b).S - FitPeakMain(sp, a).S);
    }
    for (int ir = 0; ir < kNVarRange; ++ir)
        for (int ord = kVarOrdMin; ord <= kVarOrdMax; ++ord)
            for (int iw = -1; iw <= 1; ++iw) {
                PeakCfg c = nominal;
                c.lo = kVarLo[ir]; c.hi = kVarHi[ir]; c.order = ord;
                c.sigma = nominal.sigma + iw * kVarDSigma;
                const PeakFit f = FitPeakMain(sp, c);
                if (keep) s.fits.push_back(f);
                s.nfit++;
                if (!f.ok || !(f.prob > kVarPmin)) continue;
                if (s.npass == 0) { s.smin = s.smax = f.S; }
                s.smin = std::min(s.smin, f.S); s.smax = std::max(s.smax, f.S);
                sum += f.S; sum2 += f.S * f.S;
                s.npass++;
                if (iw == 0) {
                    if (s.npassBkg == 0) { bmin = bmax = f.S; }
                    bmin = std::min(bmin, f.S); bmax = std::max(bmax, f.S);
                    s.npassBkg++;
                }
            }
    s.halfBkg = 0.5 * (bmax - bmin);
    if (s.npass > 0) {
        s.mean = sum / s.npass;
        s.rms = std::sqrt(std::max(0., sum2 / s.npass - s.mean * s.mean));
        s.half = 0.5 * (s.smax - s.smin);
    }
    return s;
}

// ---- kinematic edge of the bump in one M(pi+pi-) slice ------------------------
// N - CB = polynomial + [eta Gaussian] + h * step(m), where
// step(m) = 0.5 erfc((m - m_edge)/(sqrt(2) s)), averaged over the slice
// mpLo..mpHi (m_edge = EdgeMeasured(mpipi), uniform in mpipi), s = kEdgeSigma.
// Position and smearing are FIXED; h is the step height in rows per fine bin.
struct StepCfg {
    double mpLo = 0, mpHi = 0;      // M(pi+pi-) slice, measured scale
    bool   useEdge = false;         // true: one explicit edge position "edge" (control
    double edge = 0;                //       slices, or another variable such as M4 - edge)
    double scale = kScale;          // measured / true four-particle mass
    double scalePiPi = kScalePiPi;  // measured / true M(pi+pi-)
    double s = kEdgeSigma;          // mass resolution at the edge
    double halfRange = 0.045;       // fit range = edge -+ halfRange ...
    double lo = -1, hi = -1;        // ... unless given explicitly
    int order = 2;                  // polynomial order
    int withEta = -1;               // eta Gaussian in the model: 1 yes, 0 no, -1 auto:
                                    //   range ends below mean - sigma -> cut at mean - 3 sigma, no eta
                                    //   otherwise -> extended to at least mean + 3 sigma, with eta
    double etaMean = kPeakMean, etaSigma = kPeakSigma;
};
struct StepFit {
    bool ok = false;
    StepCfg cfg;
    double edge = 0, edgeLo = 0, edgeHi = 0;   // predicted edge: slice centre and the two ends
    double lo = 0, hi = 0;                     // fit range used
    bool withEta = false;
    double h = 0, eh = 0;                      // step height per fine bin and its error
    double chi2 = 0, chi2NoStep = 0; int ndf = 0;
    double dchi2 = 0;                          // chi2(no step) - chi2(step)
    double bkgAtEdge = 0;                      // polynomial at the edge (rows per fine bin)
    double Seta = 0, eSeta = 0;                // eta content in the slice (if in the model)
    int np = 0; double par[kMaxPar];           // polynomial (order+1), [eta], step
    double x0 = 0, dx = 1;
    std::vector<double> tx, ty;                // step template (bin centres, value) in the range
    std::vector<double> fy, fy0;               // fitted model at tx: total, and without the step term
    double Poly(double m) const {
        const double t = (m - x0) / dx; double p = 0, tp = 1;
        for (int j = 0; j <= cfg.order; ++j) { p += par[j] * tp; tp *= t; }
        return p;
    }
};
inline double StepTemplate(double m, const StepCfg& c) {
    const double q = 1.0 / (std::sqrt(2.0) * c.s);
    if (c.useEdge) return 0.5 * std::erfc((m - c.edge) * q);
    const int n = 21; double s = 0; int used = 0;
    for (int i = 0; i < n; ++i) {
        const double mp = c.mpLo + (i + 0.5) * (c.mpHi - c.mpLo) / n;
        const double e = EdgeMeasured(mp, c.scale, c.scalePiPi);
        if (e <= 0) continue;
        s += 0.5 * std::erfc((m - e) * q); used++;
    }
    return used > 0 ? s / used : 0.;
}
inline StepFit FitEdgeStep(const Spectra& sp, const StepCfg& c) {
    StepFit r; r.cfg = c;
    for (int k = 0; k < kMaxPar; ++k) r.par[k] = 0;
    if (c.useEdge) { r.edge = r.edgeLo = r.edgeHi = c.edge; }
    else {
        r.edge = EdgeMeasured(0.5 * (c.mpLo + c.mpHi), c.scale, c.scalePiPi);
        r.edgeHi = EdgeMeasured(c.mpLo, c.scale, c.scalePiPi);
        r.edgeLo = EdgeMeasured(c.mpHi, c.scale, c.scalePiPi);
        if (r.edge <= 0) return r;
    }
    r.lo = c.lo > 0 ? c.lo : r.edge - c.halfRange;
    r.hi = c.hi > 0 ? c.hi : r.edge + c.halfRange;
    if (c.withEta >= 0) r.withEta = (c.withEta == 1);
    else if (r.hi <= c.etaMean - 3.0 * c.etaSigma) r.withEta = false;
    else if (c.hi > 0) r.withEta = true;                      // explicit range: keep it
    else if (r.hi < c.etaMean - c.etaSigma) { r.hi = c.etaMean - 3.0 * c.etaSigma; r.withEta = false; }
    else { r.hi = std::max(r.hi, c.etaMean + 3.0 * c.etaSigma); r.withEta = true; }
    // range ends on the bin grid
    r.lo = sp.lo + sp.bw * std::floor((r.lo - sp.lo) / sp.bw + 0.5);
    r.hi = sp.lo + sp.bw * std::floor((r.hi - sp.lo) / sp.bw + 0.5);
    r.x0 = 0.5 * (r.lo + r.hi); r.dx = 0.5 * (r.hi - r.lo);
    std::vector<int> idx;
    for (int i = 0; i < sp.nbin; ++i) if (sp.x(i) >= r.lo && sp.x(i) < r.hi) idx.push_back(i);
    const int npts = (int)idx.size();
    const int npol = c.order + 1;
    const int np = npol + (r.withEta ? 1 : 0) + 1;
    r.np = np;
    if (np > kMaxPar || npts <= np) return r;
    std::vector<double> F(npts * np), F0(npts * (np - 1)), y(npts), var(npts), ad(npts), V0(npts * npts), V;
    for (int q = 0; q < npts; ++q) {
        const int i = idx[q];
        const double xc = sp.x(i), t = (xc - r.x0) / r.dx;
        double tp = 1; int k = 0;
        for (int j = 0; j < npol; ++j) { F[q * np + k++] = tp; tp *= t; }
        if (r.withEta) F[q * np + k++] = GaussFrac(xc - 0.5 * sp.bw, xc + 0.5 * sp.bw, c.etaMean, c.etaSigma);
        const double st = StepTemplate(xc, c);
        F[q * np + k++] = st;
        for (int j = 0; j < np - 1; ++j) F0[q * (np - 1) + j] = F[q * np + j];
        r.tx.push_back(xc); r.ty.push_back(st);
        y[q] = sp.sig[i]; ad[q] = sp.cb[i];
        var[q] = std::max(y[q] + ad[q], 1.0);
    }
    for (int q = 0; q < npts; ++q)
        for (int p = 0; p <= q; ++p) V0[q * npts + p] = V0[p * npts + q] = sp.CbCov(idx[q], idx[p]);
    auto solve = [&](int npar, const double* FF) -> LinRes {
        V = V0;
        for (int q = 0; q < npts; ++q) V[q * npts + q] += var[q];
        return LinFitCov(npts, npar, y.data(), V, FF);
    };
    LinRes a = solve(np, F.data());
    if (!a.ok) return r;
    for (int it = 0; it < 6; ++it) {
        const double hOld = a.par[np - 1];
        for (int q = 0; q < npts; ++q) {
            double f = 0;
            for (int k = 0; k < np; ++k) f += a.par[k] * F[q * np + k];
            var[q] = std::max(f + ad[q], 1.0);
        }
        a = solve(np, F.data());
        if (!a.ok) return r;
        if (std::fabs(a.par[np - 1] - hOld) < 1e-4 * std::sqrt(a.cov[np - 1][np - 1])) break;
    }
    // the same covariance for the fit without the step: the two chi2 are comparable
    const LinRes b = solve(np - 1, F0.data());
    if (!b.ok) return r;
    for (int k = 0; k < np; ++k) r.par[k] = a.par[k];
    for (int q = 0; q < npts; ++q) {
        double f = 0;
        for (int k = 0; k < np; ++k) f += a.par[k] * F[q * np + k];
        r.fy.push_back(f);
        r.fy0.push_back(f - a.par[np - 1] * F[q * np + np - 1]);
    }
    r.h = a.par[np - 1]; r.eh = std::sqrt(a.cov[np - 1][np - 1]);
    r.chi2 = a.chi2; r.chi2NoStep = b.chi2; r.ndf = npts - np;
    r.dchi2 = b.chi2 - a.chi2;
    r.bkgAtEdge = r.Poly(r.edge);
    if (r.withEta) { r.Seta = a.par[npol]; r.eSeta = std::sqrt(a.cov[npol][npol]); }
    r.ok = true;
    return r;
}

} // namespace scan

#endif // SCAN_ENGINE_H
