// ============================================================================
// sim_eta_peak.C  --  the eta PEAK in simulation: direct decay + conversion
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SMASH std + lep +
// eta owned rows (std and lep in partition mode A), weight w (PDG BR).
// The peak at the eta mass has two simulated components:
//   direct     eta -> pi+pi-e+e- (four-particle class 5, eta set,
//              W = 2 x 2.68e-4 per row)
//   conversion eta -> pi+pi-gamma, photon converted (class 4, std set)
//   total      direct + conversion (CompPeakTotal)
// X17 (class 6) only as information, BR(eta -> pi+pi- X17, X17 -> e+e-) = 1e-5
// (illustrative, of the order of the BESIII limit 1.5e-5), never in a sum.
// No combinatorial background in simulation. Absolute scale: sigma = sum w /
// L_sim (221.05 nb^-1); expected data rows = sum w * L_data / L_sim with
// L_data = 6080.5 nb^-1 (UNCONFIRMED). Data: rows of scan_compact.root,
// selection presets of SimEngine.h (vz > -160 mm, cor_ok, isBest >= 0, PT3);
// data peak content = P_main fit of N - CB (ScanEngine: Gaussian frozen at
// 0.5439 / 5.4 MeV, pol2, 0.5225-0.61), the method of the data study.
//
// Part 1  composition versus selection (pre, pre+OA>2, slides, old code,
//         peak nominal) and in classes of OA(e+e-) and M(ee); data peak
//         contents per class (own fits, and the numbers quoted by the data
//         study); fractions, absolute rows, a two-component fit of the data
//         class contents (conversion and direct scaled separately); data
//         N - CB against the simulated stack near the peak.
// Part 2  position, width and mass scale of the simulated peak (direct,
//         conversion, total; reconstructed corrected M4 and true M4) with
//         several fit choices; the data peak with the same selection.
// Part 3  efficiency of every cut for direct, conversion and total, against
//         the data-derived peak efficiencies (data study tables, as parsed in
//         tables/cuts_compare_data.txt, and own pass/fail fits at base "pre").
// Part 4  shapes: M(pipi), M(ee), OA, p(pi+pi-), p(e+e-) of direct and
//         conversion; X17 (BR 1e-5) in M(ee) and M4.
//
// Output: figs/eta_peak_*.png (.pdf for the main ones), tables/eta_peak_*.txt
// Run from /home/przygoda/HADES/FAT/FAT (part 0 = all; parts 1-4 can run in
// parallel after one compilation with part -1):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta_peak", kTRUE)' 'plots/eta_sim/sim_eta_peak.C+(0)'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"
#include "sim_ghostvar.h"

namespace sep {
using namespace simeng;

// ---------------------------------------------------------------- output helpers
FILE* gF = nullptr;
void P(const char* fmt, ...) {
    char buf[16384];
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
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs((std::string(kFigDir) + "/" + name + ".png").c_str());
    if (pdf) c->SaveAs((std::string(kFigDir) + "/" + name + ".pdf").c_str());
}

// ---------------------------------------------------------------- constants
const double kMetaG = 0.547862;                 // PDG eta mass [GeV/c2]
const double kK0sSD = 1.0048, kK0sSDe = 0.0007; // K0S M(pipi) scale sim / data (review finding 8)
const double kDataMeanStudy = 0.5439, kDataSigmaStudy = 0.0054;
const double kToRows = kLdata / kLsim;           // sum w -> expected data rows
const int kNP = 5;
const int kPresets[kNP] = {kPre, kPreOA2, kSlides, kOldCode, kPeakNom};
const char* kPresetShort[kNP] = {"pre", "pre+OA>2", "slides", "old code", "peak nominal"};

// ---------------------------------------------------------------- components
enum { cDir = 0, cConv, cTot, cX17, cDirNG, cConvNG, cDirKW, cNC };
std::vector<Comp> PeakComps() {
    return {CompPeakDirect(), CompPeakConv(), CompPeakTotal(), CompX17(1e-5), CompSigCls(5, true), CompSigCls(4, true)};
}
const int kColDir = kColPeakDirect, kColConv = kColPeakConv, kColTot = kColPeakTotal;

// keep only the rows of the peak classes (and X17) in memory
void Reduce(SimData& d, bool keepX17) {
    size_t n0 = 0, n1 = 0;
    for (int k = 0; k < kNTreeAll; ++k) {
        std::vector<Row> v;
        for (const Row& r : d.t[k]) {
            ++n0;
            if (r.tmpl == 1 && (r.cls == 4 || r.cls == 5 || (keepX17 && r.cls == 6))) v.push_back(r);
        }
        n1 += v.size();
        d.t[k].swap(v);
        std::vector<Row>().swap(v);
        d.t[k].shrink_to_fit();
    }
    printf("sep::Reduce: %zu -> %zu rows (peak classes%s)\n", n0, n1, keepX17 ? " + X17" : "");
}

// ---------------------------------------------------------------- cuts (simulation and data rows)
const double kNo = -9999.;
struct Cut {
    double vzLo = -160, oaLo = kNo, oaHi = kNo, labHi = kNo, mppLo = kNo, mppHi = kNo, etaLo = kNo, mm4Hi = kNo,
           pt4Lo = kNo, dzHi = kNo;
    int vtx = 0;   // 0 none, 1 >= 1 pion, 2 both pions, 3 all four in the vertex fit
};
inline bool On(double v) { return v != kNo; }
inline float DzPi(const Row& r) { return std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz)); }
inline float DzPi(const scan::Row&) { return std::numeric_limits<float>::quiet_NaN(); }
template <class R> inline bool PassBase(const R& r) { return r.isBest >= 0 && r.vz > -500.f && r.cor_ok == 1 && PT3(r); }
template <class R> inline bool PassCut(const R& r, const Cut& c) {
    if (c.vzLo > -1e8 && !(r.vz > c.vzLo)) return false;
    if (On(c.oaLo) && c.oaLo > 0 && !(r.oa > c.oaLo)) return false;
    if (On(c.oaHi) && !(r.oa < c.oaHi)) return false;
    if (On(c.labHi) && !(r.oa_lab < c.labHi)) return false;
    if (On(c.mppLo) && !(r.mpipi > c.mppLo)) return false;
    if (On(c.mppHi) && !(r.mpipi < c.mppHi)) return false;
    if (On(c.etaLo) && !(r.oa_eta > c.etaLo)) return false;
    if (On(c.mm4Hi) && !(r.mm4 < c.mm4Hi)) return false;
    if (On(c.pt4Lo) && c.pt4Lo > 0 && !(r.pt4 > c.pt4Lo)) return false;
    if (On(c.dzHi) && !(DzPi(r) < c.dzHi)) return false;
    const int nv = (r.pip_isUsedVertex == 1) + (r.pim_isUsedVertex == 1);
    if (c.vtx == 1 && nv < 1) return false;
    if (c.vtx == 2 && nv < 2) return false;
    if (c.vtx == 3 && !(nv == 2 && r.ep_isUsedVertex == 1 && r.em_isUsedVertex == 1)) return false;
    return true;
}
Cut CPre() { return Cut(); }
Cut CPre0() { Cut c; c.vzLo = -1e9; return c; }
Cut CPeakStart() { Cut c = CPre0(); c.labHi = 55; c.mppLo = 0.30; return c; }     // scan_extra.C (txt) peak start
Cut CPeakStartMD() { Cut c = CPre0(); c.labHi = 50; c.mppHi = 0.46; return c; }   // scan_extra.md peak start
Sel SimSel(const Cut& c) { return [c](const Row& r) { return PassBase(r) && PassCut(r, c); }; }
scan::Sel DatSel(const Cut& c) { return [c](const scan::Row& r) { return PassBase(r) && PassCut(r, c); }; }

enum Kind { kVz, kOaLo, kOaHi, kLab, kMppLo, kMppHi, kEta, kMm4, kPt4, kDz, kVtx, kNKind };
void SetCut(Cut& c, int k, double v) {
    switch (k) {
        case kVz: c.vzLo = (v == kNo) ? -1e9 : v; break;
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
std::string ValTxt(int k, double v) {   // the same text as sim_cuts.C (keys of cuts_compare_data.txt)
    if (k == kVtx) {
        static const char* n[4] = {"none", ">=1 pion", "both pions", "all four"};
        const int i = (int)v;
        return (i >= 0 && i < 4) ? n[i] : "?";
    }
    if (v == kNo || ((k == kOaLo || k == kPt4) && v == 0)) return "none";
    if (k == kMppLo || k == kMppHi || k == kMm4) return Form("%.2f", v);
    return Form("%g", v);
}
std::vector<double> Steps(double a, double b, double d) {
    std::vector<double> v;
    const int n = (int)std::lround((b - a) / d);
    for (int i = 0; i <= n; ++i) v.push_back(std::round((a + i * d) * 1000.) / 1000.);
    return v;
}
struct Fam {
    int kind; std::string key, title; std::vector<double> v; bool pre0, start; double xNone;
};
std::vector<Fam> Families() {
    std::vector<Fam> f;
    auto add = [&](int k, const char* key, const char* t, std::vector<double> v, bool p0, bool st, double xn) {
        Fam x; x.kind = k; x.key = key; x.title = t; x.v = v; x.pre0 = p0; x.start = st; x.xNone = xn; f.push_back(x);
    };
    add(kVz, "vz", "v_{z} > [mm]", {kNo, -200, -190, -180, -170, -165, -160, -155, -150, -140}, true, true, -210);
    add(kOaLo, "oa_lo", "OA > [deg]", Steps(0, 4, 0.5), false, false, 0);
    add(kOaHi, "oa_hi", "OA < [deg]", {4, 6, 8, 10, 15, 20, 30, kNo}, false, false, 35);
    std::vector<double> lab = Steps(25, 70, 5); lab.push_back(kNo);
    add(kLab, "oa_lab", "#alpha_{lab} < [deg]", lab, false, false, 75);
    std::vector<double> mlo = {kNo}; for (double x : Steps(0.28, 0.38, 0.02)) mlo.push_back(x);
    add(kMppLo, "mpipi_lo", "M_{#pi#pi} > [GeV/c^{2}]", mlo, false, false, 0.26);
    std::vector<double> mhi = Steps(0.38, 0.50, 0.01); mhi.push_back(kNo);
    add(kMppHi, "mpipi_hi", "M_{#pi#pi} < [GeV/c^{2}]", mhi, false, false, 0.52);
    add(kEta, "eta_rest", "#alpha_{#eta} > [deg]", {kNo, 100, 110, 120, 130, 140, 150, 160}, false, false, 90);
    std::vector<double> mm = Steps(2.00, 2.95, 0.05); mm.push_back(kNo);
    add(kMm4, "mm4_hi", "MM_{4} < [GeV/c^{2}]", mm, false, false, 3.0);
    add(kPt4, "pt4_lo", "p_{T,4} > [MeV/c]", Steps(0, 450, 50), false, true, 0);
    add(kDz, "dz_pi", "|z_{#pi} - v_{z}| < [mm], both #pi", {kNo, 150, 100, 60, 40, 30, 20}, false, true, 170);
    add(kVtx, "vtx", "tracks in the vertex fit", {0, 1, 2, 3}, false, true, 0);
    return f;
}
double XPlot(const Fam& f, double v) {
    if (f.kind == kVtx) return v;
    if (v == kNo) return f.xNone;
    return v;
}

// ---------------------------------------------------------------- data peak fits
// review finding 4: every data peak content carries the fit-model spread of the data study
// (scan::PeakModelSpread: P_main repeated for pol 1-4, four fit ranges 0.505-0.520 .. 0.59-0.65 and width
// 5.4 -+ 0.5 MeV = 48 fits; half of the range of S over the fits with p > 0.05). eM = that half range;
// eT() = stat (+) model. The S of every variant is kept so that pass/fail efficiencies can be recomputed
// variant by variant (the model error of an efficiency is the half range of the per-variant efficiencies).
struct DFit {
    bool ok = false; double S = 0, eS = 0, chi2 = 0, N = 0; int ndf = 0;
    double eM = 0; int nVarPass = 0;
    std::vector<double> vS; std::vector<char> vOk;
    double eT() const { return std::sqrt(eS * eS + eM * eM); }
};
DFit DoFit(const scan::Spectra& sp) {
    const scan::PeakFit f = scan::FitPeakMain(sp);
    DFit o; o.ok = f.ok; o.S = f.S; o.eS = f.eS; o.chi2 = f.chi2; o.ndf = f.ndf; o.N = f.N;
    const scan::Spread s = scan::PeakModelSpread(sp, scan::PeakCfg(), true);
    o.nVarPass = s.npass;
    for (const scan::PeakFit& v : s.fits) { o.vS.push_back(v.S); o.vOk.push_back(v.ok && v.prob > scan::kVarPmin); }
    if (s.npass >= 2) o.eM = s.half;
    else {   // fewer than two acceptable variants: looser p > 1e-3, nominal included; otherwise undetermined (eM = 0, flagged)
        double lo = o.S, hi = o.S; int n = 0;
        for (const scan::PeakFit& v : s.fits) if (v.ok && v.prob > 1e-3) { lo = std::min(lo, v.S); hi = std::max(hi, v.S); ++n; }
        o.eM = n >= 2 ? 0.5 * (hi - lo) : 0.;
        o.nVarPass = n >= 2 ? -n : 0;   // negative: loosened criterion
    }
    return o;
}
std::vector<DFit> FitDataSels(const scan::Data& dd, const std::vector<scan::Sel>& sels) {
    std::vector<DFit> out;
    for (size_t i0 = 0; i0 < sels.size(); i0 += 40) {
        const size_t i1 = std::min(sels.size(), i0 + 40);
        std::vector<scan::Sel> b(sels.begin() + i0, sels.begin() + i1);
        const std::vector<scan::Spectra> sp = scan::FillMany(dd, b);
        for (const auto& s : sp) out.push_back(DoFit(s));
    }
    return out;
}
// err = stat (+) model (used in every pull since the review); estat = statistical only; esys = half range
// of the efficiency over the model variants where both the pass and the fail fit are acceptable (p > 0.05)
// mode: 0 = variants with both fits acceptable (p > 0.05); 1 = fewer than two such variants, the pass fit acceptable;
// 2 = the fail fit acceptable; 3 = no acceptable variant (model error undetermined, esys = 0). The nominal value is
// always inside the range, so the model error covers the distance from the nominal to the variants.
struct Eff { double e = 0, err = 0, estat = 0, esys = 0; int nv = 0, mode = 0; bool ok = false; };
Eff PassFail(const DFit& p, const DFit& f) {
    Eff o;
    const double t = p.S + f.S;
    if (!(t > 0)) return o;
    o.e = p.S / t;
    o.estat = std::sqrt(std::pow((1 - o.e) * p.eS, 2) + std::pow(o.e * f.eS, 2)) / t;
    for (int mode = 0; mode < 3; ++mode) {
        double lo = o.e, hi = o.e; int n = 0;
        for (size_t i = 0; i < p.vS.size() && i < f.vS.size(); ++i) {
            const bool use = mode == 0 ? (p.vOk[i] && f.vOk[i]) : (mode == 1 ? (bool)p.vOk[i] : (bool)f.vOk[i]);
            const double tv = p.vS[i] + f.vS[i];
            if (!use || !(tv > 0)) continue;
            const double ev = p.vS[i] / tv;
            lo = std::min(lo, ev); hi = std::max(hi, ev); ++n;
        }
        if (n >= 2) { o.esys = 0.5 * (hi - lo); o.nv = n; o.mode = mode; break; }
        o.mode = 3;
    }
    o.err = std::hypot(o.estat, o.esys);
    o.ok = true;
    return o;
}

// ---------------------------------------------------------------- classes of OA and M(ee)
struct ClassSet {
    std::string tag, title; int var;   // var 0: OA [deg], 1: M(ee) [GeV/c2]
    std::vector<double> lo, hi; std::vector<std::string> name;
    std::vector<double> given, givenE; std::string givenSrc;   // data contents quoted elsewhere (rows), NaN = none
};
template <class R> inline double CVar(int var, const R& r) { return var == 0 ? (double)r.oa : (double)r.mee; }
const double kNaN = std::numeric_limits<double>::quiet_NaN();
std::vector<ClassSet> ClassSets() {
    std::vector<ClassSet> v;
    {
        ClassSet c; c.tag = "OA5"; c.title = "OA(e^{+}e^{-}) classes"; c.var = 0;
        c.lo = {0, 2, 4, 6, 10}; c.hi = {2, 4, 6, 10, 1e9}; c.name = {"<2", "2-4", "4-6", "6-10", ">10"};
        v.push_back(c);
    }
    {
        ClassSet c; c.tag = "MEE5"; c.title = "M(e^{+}e^{-}) classes"; c.var = 1;
        c.lo = {0, 0.008, 0.020, 0.050, 0.135}; c.hi = {0.008, 0.020, 0.050, 0.135, 1e9};
        c.name = {"<8", "8-20", "20-50", "50-135", ">135"};
        c.given = {1541, 1450, 688, 550, 83}; c.givenE = {194, 205, 130, 141, 68};
        c.givenSrc = "scan_peak.md sec. 5 (nominal, NO vz cut)";
        v.push_back(c);
    }
    {
        ClassSet c; c.tag = "OA8t"; c.title = "OA(e^{+}e^{-}) fine classes"; c.var = 0;
        c.lo = {0, 1, 2, 3, 4, 6, 10, 20}; c.hi = {1, 2, 3, 4, 6, 10, 20, 1e9};
        c.name = {"0-1", "1-2", "2-3", "3-4", "4-6", "6-10", "10-20", ">20"};
        c.given = {241, 798, 855, 851, 667, 365, 553, 374}; c.givenE = {78, 143, 159, 137, 131, 109, 112, 149};
        c.givenSrc = "assistant summary of today (vz > -160)";
        v.push_back(c);
    }
    {
        ClassSet c; c.tag = "MEE5t"; c.title = "M(e^{+}e^{-}) classes (today)"; c.var = 1;
        c.lo = {0, 0.010, 0.020, 0.050, 0.135}; c.hi = {0.010, 0.020, 0.050, 0.135, 0.270};
        c.name = {"0-10", "10-20", "20-50", "50-135", "135-270"};
        c.given = {2309, 1070, 751, 465, 99}; c.givenE = {kNaN, 175, 137, 139, 64};
        c.givenSrc = "assistant summary of today (vz > -160)";
        v.push_back(c);
    }
    {
        ClassSet c; c.tag = "OA4sp"; c.title = "OA(e^{+}e^{-}) classes (data study)"; c.var = 0;
        c.lo = {0, 2, 6, 15}; c.hi = {2, 6, 15, 1e9}; c.name = {"<2", "2-6", "6-15", ">15"};
        c.given = {944, 2131, 577, 667}; c.givenE = {154, 230, 130, 165};
        c.givenSrc = "scan_peak.md sec. 5 (nominal, NO vz cut)";
        v.push_back(c);
    }
    return v;
}
ClassSet FindCS(const std::string& tag) { for (const ClassSet& c : ClassSets()) if (c.tag == tag) return c; return ClassSet(); }

// two-component fit of data class contents: D_i = a C_i + b X_i (conversion C, direct X; expected rows)
struct MixRes { bool ok = false; double a = 0, ea = 0, b = 0, eb = 0, rho = 0, chi2 = 0; int ndf = 0;
                double a1 = 0, ea1 = 0, chi21 = 0; int ndf1 = 0; double fdir = 0, efdir = 0, fdirSim = 0; };
MixRes MixFit(const std::vector<double>& D, const std::vector<double>& eD, const std::vector<double>& C,
              const std::vector<double>& eC, const std::vector<double>& X, const std::vector<double>& eX) {
    MixRes m;
    const int n = (int)D.size();
    double a = 0.3, b = 0.3, Saa = 0, Sab = 0, Sbb = 0;
    for (int it = 0; it < 50; ++it) {
        double Sda = 0, Sdb = 0; Saa = Sab = Sbb = 0;
        for (int i = 0; i < n; ++i) {
            if (!(eD[i] > 0)) continue;
            const double w = 1. / (eD[i] * eD[i] + a * a * eC[i] * eC[i] + b * b * eX[i] * eX[i]);
            Saa += w * C[i] * C[i]; Sab += w * C[i] * X[i]; Sbb += w * X[i] * X[i];
            Sda += w * D[i] * C[i]; Sdb += w * D[i] * X[i];
        }
        const double det = Saa * Sbb - Sab * Sab;
        if (!(det > 0)) return m;
        a = (Sda * Sbb - Sdb * Sab) / det; b = (Saa * Sdb - Sab * Sda) / det;
    }
    const double det = Saa * Sbb - Sab * Sab;
    m.a = a; m.b = b; m.ea = std::sqrt(Sbb / det); m.eb = std::sqrt(Saa / det); m.rho = -Sab / det / (m.ea * m.eb);
    int np = 0;
    for (int i = 0; i < n; ++i) {
        if (!(eD[i] > 0)) continue;
        ++np;
        m.chi2 += std::pow(D[i] - a * C[i] - b * X[i], 2) / (eD[i] * eD[i] + a * a * eC[i] * eC[i] + b * b * eX[i] * eX[i]);
    }
    m.ndf = np - 2;
    // one factor for the simulated mixture
    double a1 = 0.3, S11 = 0;
    for (int it = 0; it < 50; ++it) {
        double Sd = 0; S11 = 0;
        for (int i = 0; i < n; ++i) {
            if (!(eD[i] > 0)) continue;
            const double T = C[i] + X[i], w = 1. / (eD[i] * eD[i] + a1 * a1 * (eC[i] * eC[i] + eX[i] * eX[i]));
            S11 += w * T * T; Sd += w * D[i] * T;
        }
        a1 = Sd / S11;
    }
    m.a1 = a1; m.ea1 = 1. / std::sqrt(S11);
    for (int i = 0; i < n; ++i) {
        if (!(eD[i] > 0)) continue;
        m.chi21 += std::pow(D[i] - a1 * (C[i] + X[i]), 2) / (eD[i] * eD[i] + a1 * a1 * (eC[i] * eC[i] + eX[i] * eX[i]));
    }
    m.ndf1 = np - 1;
    // direct fraction of the data peak implied by the fit, and in the simulation
    double sC = 0, sX = 0;
    for (int i = 0; i < n; ++i) { sC += C[i]; sX += X[i]; }
    const double T = a * sC + b * sX;
    m.fdir = b * sX / T;
    const double dfa = -b * sX * sC / (T * T), dfb = a * sC * sX / (T * T);
    const double cab = m.rho * m.ea * m.eb;
    m.efdir = std::sqrt(dfa * dfa * m.ea * m.ea + dfb * dfb * m.eb * m.eb + 2 * dfa * dfb * cab);
    m.fdirSim = sX / (sC + sX);
    m.ok = true;
    return m;
}

// ---------------------------------------------------------------- drawing helpers
void StyleH(TH1* h, int col, int lw, int ls, int fill = 0, int fstyle = 1001) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls);
    h->SetMarkerColor(col); h->SetMarkerSize(0);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(fstyle); } else h->SetFillStyle(0);
}
void StyleG(TGraph* g, int col, int mk, double ms, int lw = 2) {
    g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineWidth(lw);
}
void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy = false, double lm = 0.14,
            double rm = 0.03, double tm = 0.10, double bm = 0.14) {
    TPad* p = new TPad(Nm(name).c_str(), "", x1, y1, x2, y2);
    p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
    p->SetLogy(logy);
    p->Draw();
    p->cd();
    return p;
}
TH1D* Frame(double xlo, double xhi, double ylo, double yhi, const char* xt, const char* yt, double ts = 0.05) {
    TH1D* h = new TH1D(Nm("fr").c_str(), "", 100, xlo, xhi);
    h->SetDirectory(nullptr);
    h->SetMinimum(ylo); h->SetMaximum(yhi);
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 0.85); h->GetYaxis()->SetLabelSize(ts * 0.85);
    h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetMaxDigits(4);
    h->Draw("AXIS");
    return h;
}
TH1D* CatFrame(const std::vector<std::string>& lab, double ylo, double yhi, const char* xt, const char* yt, double ts = 0.05) {
    const int n = (int)lab.size();
    TH1D* h = new TH1D(Nm("cf").c_str(), "", n, 0, n);
    h->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) h->GetXaxis()->SetBinLabel(i + 1, lab[i].c_str());
    h->SetMinimum(ylo); h->SetMaximum(yhi);
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 1.05); h->GetYaxis()->SetLabelSize(ts * 0.85);
    h->GetXaxis()->SetTitleOffset(1.15); h->GetYaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetMaxDigits(4);
    h->Draw("AXIS");
    return h;
}
TH1D* CatHist(const std::vector<double>& y, const std::vector<double>& e) {
    const int n = (int)y.size();
    TH1D* h = new TH1D(Nm("ch").c_str(), "", n, 0, n);
    h->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) { h->SetBinContent(i + 1, y[i]); h->SetBinError(i + 1, e[i]); }
    return h;
}
TGraphErrors* CatGraph(const std::vector<double>& y, const std::vector<double>& e, double dx) {
    TGraphErrors* g = new TGraphErrors();
    for (size_t i = 0; i < y.size(); ++i) {
        if (!(y[i] == y[i])) continue;
        const int k = g->GetN();
        g->SetPoint(k, i + 0.5 + dx, y[i]);
        g->SetPointError(k, 0, (e[i] == e[i]) ? e[i] : 0.);
    }
    return g;
}
TH1D* Dummy(int lc, int lw, int ls, int fc = 0, int mk = 0, double ms = 0) {
    TH1D* h = new TH1D(Nm("dm").c_str(), "", 1, 0, 1); h->SetDirectory(nullptr);
    h->SetLineColor(lc); h->SetLineWidth(lw); h->SetLineStyle(ls);
    if (fc) { h->SetFillColor(fc); h->SetFillStyle(1001); } else h->SetFillStyle(0);
    h->SetMarkerColor(lc); h->SetMarkerStyle(mk ? mk : 1); h->SetMarkerSize(ms);
    return h;
}
const char* kSimNote = "SMASH std + lep + eta, owned rows (mode A), w (PDG); #sigma = #Sigma w / L_{sim}, L_{sim} = 221.05 nb^{-1}; no efficiency correction";
const char* kDataNote = "data: peak content = P_{main} fit of N - CB (Gaussian frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61); L_{data} = 6080.5 nb^{-1}, UNCONFIRMED";

// weighted quantile
double WQuant(std::vector<std::pair<double, double>> v, double q) {
    if (v.empty()) return kNaN;
    std::sort(v.begin(), v.end());
    double t = 0; for (const auto& p : v) t += p.second;
    double a = 0;
    for (const auto& p : v) { a += p.second; if (a >= q * t) return p.first; }
    return v.back().first;
}

// ---------------------------------------------------------------- M(ee) model check of the direct decay
// Kroll-Wada factor times a simple pi pi gamma* rate: dGamma/dq / Gamma(pi pi gamma) = (2 alpha / 3 pi q)
// (1 + 2 m^2/q^2) sqrt(1 - 4 m^2/q^2) F(q), F(q) = int dm (p_pi^3/m) k(q,m)^3 / int dm (p_pi^3/m) k(0,m)^3
// (P-wave pion pair, photon-momentum cubed; no form factor, no box anomaly, no rho). Integrated it gives
// 6.5e-3 against the PDG ratio 2.68e-4 / 4.22e-2 = 6.35e-3. Used only as a reference SHAPE.
// KW reweighting variant of the simulated direct decay: w' = w x model(q)/accepted(q) per bin of the TRUE
// M(ee) q, normalised so that the direct content in "pre" is unchanged. It assumes that the acceptance does
// not depend on q (UNVERIFIED): it shows what the class fractions would be if the generator followed QED.
double KWLam(double a, double b, double c) { return a * a + b * b + c * c - 2 * a * b - 2 * a * c - 2 * b * c; }
double KWF(double q) {
    const double me = kMetaG, mp = 0.13957;
    double s0 = 0, s1 = 0;
    for (int i = 0; i < 1000; ++i) {
        const double m = 2 * mp + (me - 2 * mp) * (i + 0.5) / 1000.;
        const double p = std::sqrt(std::max(0., m * m / 4 - mp * mp));
        const double k0 = (me * me - m * m) / (2 * me);
        const double l = KWLam(me * me, m * m, q * q);
        const double k = l > 0 ? std::sqrt(l) / (2 * me) : 0.;
        const double w = std::pow(p, 3) / m;
        s0 += w * std::pow(k0, 3);
        if (m + q < me) s1 += w * std::pow(k, 3);
    }
    return s1 / s0;
}
double KWDens(double q) {   // dGamma/dq / Gamma(pi pi gamma) [1/GeV]
    const double a = 1 / 137.036, m = 0.00051099895;
    if (q <= 2 * m) return 0.;
    return (2 * a / (3 * M_PI * q)) * (1 + 2 * m * m / (q * q)) * std::sqrt(1 - 4 * m * m / (q * q)) * KWF(q);
}
struct KWTab {
    std::vector<double> e, model, acc, fac;   // edges [GeV], model and accepted fractions, factor model/acc
    double norm = 1, ratioBR = 0;
    int bin(double q) const {
        const int n = (int)e.size() - 1;
        if (q < e[0]) return 0;
        for (int k = 0; k < n; ++k) if (q < e[k + 1]) return k;
        return n - 1;
    }
    double W(double q) const { return fac[bin(q)] * norm; }
};
KWTab MakeKW(const SimData& d) {
    KWTab t;
    const int nb = 14;
    for (int i = 0; i <= nb; ++i) t.e.push_back(0.001022 * std::pow(0.2685 / 0.001022, (double)i / nb));
    t.model.assign(nb, 0.); t.acc.assign(nb, 0.); t.fac.assign(nb, 1.);
    double tm = 0;
    for (int k = 0; k < nb; ++k) {
        const int ns = 100;
        for (int j = 0; j < ns; ++j) {
            const double q = t.e[k] + (t.e[k + 1] - t.e[k]) * (j + 0.5) / ns;
            t.model[k] += KWDens(q) * (t.e[k + 1] - t.e[k]) / ns;
        }
        tm += t.model[k];
    }
    t.ratioBR = tm;
    for (double& x : t.model) x /= tm;
    double ta = 0;
    for (int set = 1; set <= 3; ++set)
        for (const Row& r : d.tree(set, kUS)) {
            if (!(r.tmpl == 1 && r.cls == 5) || !Owned(r, SimOpt()) || !PassPreset(kPre, r) || !(r.mee_true > 0)) continue;
            t.acc[t.bin(r.mee_true)] += r.w; ta += r.w;
        }
    for (double& x : t.acc) x /= ta;
    for (int k = 0; k < nb; ++k) t.fac[k] = t.acc[k] > 0 ? t.model[k] / t.acc[k] : 1.;
    double s0 = 0, s1 = 0;
    for (int set = 1; set <= 3; ++set)
        for (const Row& r : d.tree(set, kUS)) {
            if (!(r.tmpl == 1 && r.cls == 5) || !Owned(r, SimOpt()) || !PassPreset(kPre, r)) continue;
            s0 += r.w; s1 += r.w * t.fac[t.bin(r.mee_true)];
        }
    t.norm = s0 / s1;
    return t;
}
Comp CompDirectKW(const KWTab& t) {
    Comp c{"peak direct, KW-reweighted M(ee)", [](const Row& r) { return r.tmpl == 1 && r.cls == 5; }};
    c.wf = [t](const Row& r) { return r.w * t.W(r.mee_true); };
    return c;
}
Comp CompTotalKW(const KWTab& t) {
    Comp c{"peak total, KW direct", [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5); }};
    c.wf = [t](const Row& r) { return r.cls == 5 ? r.w * t.W(r.mee_true) : r.w; };
    return c;
}

// =====================================================================================================
// Part 1: composition, classes, data against simulation near the peak
// =====================================================================================================
void Part1(SimData& d, const scan::Data& dd) {
    TStopwatch sw;
    // ---------------- 1a. data N - CB against the simulated stack near the peak (full simulation)
    {
        Binning b; b.nbin = 72; b.lo = 0.46; b.hi = 0.64;   // 2.5 MeV
        const std::vector<Comp> cc = {
            CompPhysBkg(),
            {"bump", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); }},
            CompPeakConv(), CompPeakDirect(), CompSim(), CompSimVar(1),
            {"sim without direct", [](const Row& r) { return ((r.tmpl == 1 && !X17Row(r)) || r.tmpl == 2) && r.cls != 5; }},
            CompX17(1e-5)};
        enum { kB = 0, kBu, kPC, kPD, kS, kSa, kSnd, kX, kNCC };
        TH1D* hs[kNP][kNCC]; TH1D* hd[kNP]; TH1D* hr[kNP]; TH1D* hra[kNP];
        OpenTab("eta_peak_m4_datasim.txt");
        P("# eta_peak_m4_datasim: data N - CB against the simulated M4 spectrum near the peak, absolute units\n");
        P("# simulation: std + lep + eta owned rows, mode A, w (PDG); T1 + T2 (default) and variant (a) T1 + T2 + ghost-flagged true pairs\n");
        P("# data: (N - CB) / L_data, CB = 2 sqrt(N++ N--) in 10 MeV bins; L_data = %.1f nb^-1 (UNCONFIRMED); 2.5 MeV bins, sums over bins whose centre is inside\n", kLdata);
        P("%-14s %-16s %20s %20s %20s %10s %10s %18s\n", "selection", "M4 range", "data N-CB [nb]", "sim T1+T2 [nb]", "sim (a) [nb]", "data/sim", "data/(a)", "sim peak (d+c) [nb]");
        for (int ip = 0; ip < kNP; ++ip) {
            const int p = kPresets[ip];
            const std::vector<WSpec> s = FillMany(d, SelPreset(p), cc, SimOpt(), b);
            for (int c = 0; c < kNCC; ++c) hs[ip][c] = s[c].Hist(Nm("hs").c_str(), 2);
            const DSpec ds = DataSpectrumPreset(dd, p, b.nbin, b.lo, b.hi);
            hd[ip] = ds.Hist("sig", Nm("hd").c_str(), 2);
            hr[ip] = (TH1D*)hd[ip]->Clone(Nm("hr").c_str()); hr[ip]->Divide(hs[ip][kS]);
            hra[ip] = (TH1D*)hd[ip]->Clone(Nm("hra").c_str()); hra[ip]->Divide(hs[ip][kSa]);
            const double rg[4][2] = {{0.46, scan::kBumpHi}, {scan::kBumpHi, scan::kPeakLo}, {scan::kPeakLo, scan::kPeakHi}, {0.56, 0.64}};
            const char* rn[4] = {"0.46-0.5225", "0.5225-0.5336", "peak window", "0.56-0.64"};
            for (int k = 0; k < 4; ++k) {
                double D = 0, eD = 0, S = 0, eS = 0, Sa = 0, eSa = 0, Pk = 0, ePk = 0;
                for (int i = 1; i <= b.nbin; ++i) {
                    const double x = hd[ip]->GetXaxis()->GetBinCenter(i);
                    if (!(x >= rg[k][0] && x < rg[k][1])) continue;
                    D += hd[ip]->GetBinContent(i) * b.bw(); eD += std::pow(hd[ip]->GetBinError(i) * b.bw(), 2);
                    S += hs[ip][kS]->GetBinContent(i) * b.bw(); eS += std::pow(hs[ip][kS]->GetBinError(i) * b.bw(), 2);
                    Sa += hs[ip][kSa]->GetBinContent(i) * b.bw(); eSa += std::pow(hs[ip][kSa]->GetBinError(i) * b.bw(), 2);
                    Pk += (hs[ip][kPC]->GetBinContent(i) + hs[ip][kPD]->GetBinContent(i)) * b.bw();
                    ePk += (std::pow(hs[ip][kPC]->GetBinError(i), 2) + std::pow(hs[ip][kPD]->GetBinError(i), 2)) * b.bw() * b.bw();
                }
                const double R = D / S, eR = R * std::sqrt(eD / (D * D) + eS / (S * S));
                const double Ra = D / Sa, eRa = Ra * std::sqrt(eD / (D * D) + eSa / (Sa * Sa));
                P("%-14s %-16s %9.4f +- %-8.4f %9.4f +- %-8.4f %9.4f +- %-8.4f %5.3f+-%.3f %5.3f+-%.3f %8.4f +- %-7.4f\n", kPresetShort[ip], rn[k],
                  D, std::sqrt(eD), S, std::sqrt(eS), Sa, std::sqrt(eSa), R, eR, Ra, eRa, Pk, std::sqrt(ePk));
            }
        }
        CloseTab();
        // figure
        gStyle->SetOptStat(0); gStyle->SetEndErrorSize(0);
        TCanvas* cv = new TCanvas(Nm("cM4").c_str(), "", 3000, 1500);
        double yMax = 0;
        for (int ip = 0; ip < kNP; ++ip)
            for (int i = 1; i <= hd[ip]->GetNbinsX(); ++i) {
                yMax = std::max(yMax, hd[ip]->GetBinContent(i) + hd[ip]->GetBinError(i));
                yMax = std::max(yMax, hs[ip][kSa]->GetBinContent(i));
            }
        cv->cd();
        TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.90, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
        Txt(0.01, 0.70, "d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} near the #eta peak: data N - CB against the simulated components (stacked), absolute units, M_{4} < 1 GeV/c^{2} region", 0.26);
        Txt(0.01, 0.40, Form("%s; data (N - CB)/L_{data}/bin width, L_{data} = 6080.5 nb^{-1} UNCONFIRMED; X17 not in the sum", kSimNote), 0.19);
        Txt(0.01, 0.12, "every selection: isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta; 2.5 MeV bins; dashed black: variant (a) with ghost-flagged true pairs; dotted: simulation without the direct decay", 0.19);
        for (int ip = 0; ip < kNP; ++ip) {
            cv->cd();
            const double x1 = ip * 0.2, x2 = x1 + 0.2;
            MkPad("pt", x1, 0.32, x2, 0.90, false, 0.17, 0.02, 0.09, 0.02);
            TH1D* fr = Frame(0.46, 0.64, 0, yMax * 1.25, "", ip == 0 ? "d#sigma/dM [nb/(GeV/c^{2})]" : "", 0.055);
            fr->GetXaxis()->SetLabelSize(0); fr->GetXaxis()->SetNdivisions(505); fr->GetYaxis()->SetTitleOffset(1.5);
            TH1D* st[4];
            st[0] = (TH1D*)hs[ip][kB]->Clone(Nm("st").c_str());
            st[1] = (TH1D*)st[0]->Clone(Nm("st").c_str()); st[1]->Add(hs[ip][kBu]);
            st[2] = (TH1D*)st[1]->Clone(Nm("st").c_str()); st[2]->Add(hs[ip][kPC]);
            st[3] = (TH1D*)st[2]->Clone(Nm("st").c_str()); st[3]->Add(hs[ip][kPD]);
            StyleH(st[3], kColDir, 1, 1, kColDir); StyleH(st[2], kColConv, 1, 1, kColConv);
            StyleH(st[1], kColBumpDalitz, 1, 1, kColBumpDalitz); StyleH(st[0], kColPhysBkg, 1, 1, kColPhysFill);
            for (int k = 3; k >= 0; --k) st[k]->Draw("HIST SAME");
            TH1D* hsa = hs[ip][kSa]; StyleH(hsa, kBlack, 2, 2); hsa->Draw("HIST SAME");
            TH1D* hnd = hs[ip][kSnd]; StyleH(hnd, kBlack, 2, 3); hnd->Draw("HIST SAME");
            TH1D* hx = hs[ip][kX]; StyleH(hx, kColX17, 2, 2); hx->Draw("HIST SAME");
            hd[ip]->SetMarkerStyle(20); hd[ip]->SetMarkerSize(1.0); hd[ip]->SetMarkerColor(kBlue); hd[ip]->SetLineColor(kBlue);
            hd[ip]->Draw("E1 SAME");
            for (double x : {scan::kPeakLo, scan::kPeakHi}) { TLine* l = new TLine(x, 0, x, yMax * 1.25); l->SetLineStyle(3); l->SetLineColor(kGray + 1); l->Draw(); }
            gPad->RedrawAxis();
            Txt(0.20, 0.93, Form("#bf{%s}", kPresetShort[ip]), 0.065);
            if (ip == 0) {
                TLegend* lg = new TLegend(0.19, 0.50, 0.97, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
                lg->AddEntry(hd[ip], "data N - CB", "pe");
                lg->AddEntry(st[3], "peak: #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} direct", "f");
                lg->AddEntry(st[2], "peak: #pi^{+}#pi^{-}#gamma, conversion", "f");
                lg->AddEntry(st[1], "bump (#pi^{0} Dalitz + conv.)", "f");
                lg->AddEntry(st[0], "physical background", "f");
                lg->AddEntry(hsa, "sum, variant (a)", "l");
                lg->AddEntry(hnd, "sum without direct", "l");
                lg->AddEntry(hx, "X17, BR 10^{-5} (info)", "l");
                lg->Draw();
            }
            cv->cd();
            MkPad("pr", x1, 0.0, x2, 0.32, false, 0.17, 0.02, 0.02, 0.26);
            TH1D* fr2 = Frame(0.46, 0.64, 0.3, 1.9, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", ip == 0 ? "data / sim" : "", 0.10);
            fr2->GetYaxis()->SetNdivisions(505); fr2->GetYaxis()->SetTitleOffset(0.7); fr2->GetXaxis()->SetNdivisions(505);
            TLine* l1 = new TLine(0.46, 1, 0.64, 1); l1->SetLineStyle(2); l1->Draw();
            hr[ip]->SetMarkerStyle(20); hr[ip]->SetMarkerSize(0.8); hr[ip]->SetMarkerColor(kBlue); hr[ip]->SetLineColor(kBlue); hr[ip]->Draw("E1 SAME");
            hra[ip]->SetMarkerStyle(24); hra[ip]->SetMarkerSize(0.8); hra[ip]->SetMarkerColor(kGray + 2); hra[ip]->SetLineColor(kGray + 2); hra[ip]->Draw("E1 SAME");
            if (ip == 0) Txt(0.20, 0.85, "#color[4]{#bullet} / T1+T2,  #color[922]{#circ} / variant (a)", 0.085);
        }
        Save(cv, "eta_peak_m4_datasim", true);
        printf("Part1a done %.0f s\n", sw.RealTime()); sw.Continue();
    }

    Reduce(d, true);
    const KWTab kw = MakeKW(d);
    std::vector<Comp> comps = PeakComps();
    comps.push_back(CompDirectKW(kw));

    // ---------------- data fits: totals per preset, classes, OA > 2 pass/fail
    struct Combo { int preset; std::string cs; };   // preset -1: nominal without the vz cut (data only; sim uses peak nominal)
    const std::vector<Combo> combos = {{kPre, "OA5"}, {kPre, "MEE5"}, {kPre, "OA8t"}, {kPre, "MEE5t"}, {kPre, "OA4sp"},
                                       {kPreOA2, "MEE5"}, {kPeakNom, "OA4sp"}, {kPeakNom, "MEE5"}, {-1, "OA4sp"}, {-1, "MEE5"}};
    auto dataPresetSel = [](int p) -> scan::Sel {
        if (p == -1) return [](const scan::Row& r) { return r.isBest >= 0 && r.vz > -500.f && r.cor_ok == 1 && r.oa_lab < 55.f && r.mpipi > 0.30f; };
        return DataSelPreset(p);
    };
    std::vector<scan::Sel> dsel;
    for (int ip = 0; ip < kNP; ++ip) dsel.push_back(DataSelPreset(kPresets[ip]));          // 0..4 totals
    dsel.push_back(dataPresetSel(-1));                                                      // 5 nominal no vz
    dsel.push_back([](const scan::Row& r) { return PassPreset(kPre, r) && r.oa <= 2.f; });  // 6 fail of OA > 2
    std::vector<std::vector<int>> cidx(combos.size());
    for (size_t k = 0; k < combos.size(); ++k) {
        const ClassSet cs = FindCS(combos[k].cs);
        const scan::Sel base = dataPresetSel(combos[k].preset);
        for (size_t i = 0; i < cs.lo.size(); ++i) {
            const double lo = cs.lo[i], hi = cs.hi[i]; const int var = cs.var;
            cidx[k].push_back((int)dsel.size());
            dsel.push_back([base, lo, hi, var](const scan::Row& r) { if (!base(r)) return false; const double x = CVar(var, r); return x >= lo && x < hi; });
        }
    }
    const std::vector<DFit> df = FitDataSels(dd, dsel);
    printf("Part1: %zu data fits done %.0f s\n", df.size(), sw.RealTime()); sw.Continue();

    // ---------------- 1b. composition versus selection
    OpenTab("eta_peak_composition.txt");
    P("# eta_peak_composition: the simulated eta peak (direct eta -> pi+pi-e+e- + conversion eta -> pi+pi-gamma) per selection\n");
    P("# %s\n", "simulation: std + lep + eta owned rows, mode A, w (PDG); sigma = sum w / L_sim (221.05 nb^-1), within the acceptance, no efficiency correction");
    P("# expected data rows = sum w x L_data / L_sim, L_data = %.1f nb^-1 (UNCONFIRMED). Errors: event groups. All M4 unless stated.\n", kLdata);
    P("# data: P_main fit of N - CB (Gaussian frozen 0.5439 / 5.4 MeV = whole Gaussian, pol2, 0.5225-0.61), same preset (vz > -160)\n");
    P("# dir KW: ILLUSTRATION ONLY - direct decay reweighted to the QED (Kroll-Wada) M(ee) shape with flat acceptance in the true M(ee) and fixed content in 'pre'\n");
    P("#   ASSUMED (both arbitrary; the OA > 2 pass share rises from 0.52 to 1.00 with the true M(ee), eta_peak_mee_true.txt): not a result\n");
    P("# CAUTION (review finding 2): the simulated direct decay has a NON-QED generator M(ee) (accepted density largest at threshold, close to 1/M^2):\n");
    P("#   every simulated direct-decay number (yield, f_dir, class shapes, OA efficiencies) is labelled 'generator M(ee) model not QED' and is not for physics\n");
    P("# CAUTION (review finding 1): data/sim here = METHOD A only (data P_main, whole Gaussian, against the simulated content over ALL M4, 23 %% of it\n");
    P("#   outside the window); the peak ratio depends on the method (0.3 to about 1): see tables/eta_summary_peak_methods.txt; not a headline number\n");
    P("# data S: +- stat +- model (fit-model spread of the data study: pol 1-4, 4 ranges, width -+ 0.5 MeV, half range over fits with p > 0.05)\n");
    P("# data/sim: +- stat +- model; NOT included: absolute scale incl. L_data (about +-20 %%, reconcile.md) and the ghost treatment (column ghost;\n");
    P("#   'data/sim (b)' = against the simulated peak WITHOUT rows with a ghost-flagged lepton)\n");
    P("# f_dir = direct / (direct + conversion). ghost: share of rows with a ghost-flagged lepton track (kept in the signal by the template rule)\n\n");
    P("%-13s %16s %16s %16s %13s %16s %16s %21s %16s %21s %14s %15s %15s %12s\n", "selection", "sig dir [nb]", "sig conv [nb]", "sig tot [nb]", "f_dir",
      "sig tot win [nb]", "exp. rows tot", "data S +-st +-mod", "data S/L [nb]", "data/sim +-st +-mod", "data/sim (b)", "ghost dir/conv", "X17 1e-5 [nb]", "dir KW [nb]");
    std::vector<WSum> comp0(cNC), compWin(cNC);
    double simTot[kNP], simTotE[kNP], simDir[kNP], simDirE[kNP], simConv[kNP], simConvE[kNP], fdir[kNP], fdirE[kNP], fdirW[kNP], fdirWE[kNP];
    double datS[kNP], datSE[kNP], datSM[kNP];
    for (int ip = 0; ip < kNP; ++ip) {
        const int p = kPresets[ip];
        const std::vector<WSum> s = SumMany(d, SelPreset(p), comps);
        const std::vector<WSum> w = SumMany(d, SelAnd(SelPreset(p), [](const Row& r) { return r.m4 >= scan::kPeakLo && r.m4 < scan::kPeakHi; }), comps);
        if (ip == 0) { comp0 = s; compWin = w; }
        const double D = s[cDir].S, C = s[cConv].S, eD = std::sqrt(s[cDir].V), eC = std::sqrt(s[cConv].V);
        const double f = D / (D + C), ef = std::sqrt(std::pow(C * eD, 2) + std::pow(D * eC, 2)) / std::pow(D + C, 2);
        const double Dw = w[cDir].S, Cw = w[cConv].S, eDw = std::sqrt(w[cDir].V), eCw = std::sqrt(w[cConv].V);
        const double fw = Dw / (Dw + Cw), efw = std::sqrt(std::pow(Cw * eDw, 2) + std::pow(Dw * eCw, 2)) / std::pow(Dw + Cw, 2);
        simTot[ip] = s[cTot].S * kToRows; simTotE[ip] = std::sqrt(s[cTot].V) * kToRows;
        simDir[ip] = D * kToRows; simDirE[ip] = eD * kToRows; simConv[ip] = C * kToRows; simConvE[ip] = eC * kToRows;
        fdir[ip] = f; fdirE[ip] = ef; fdirW[ip] = fw; fdirWE[ip] = efw;
        datS[ip] = df[ip].S; datSE[ip] = df[ip].eS; datSM[ip] = df[ip].eM;
        const double R = datS[ip] / simTot[ip], eR = R * std::sqrt(std::pow(datSE[ip] / datS[ip], 2) + std::pow(simTotE[ip] / simTot[ip], 2));
        const double eRm = R * datSM[ip] / datS[ip];
        const double gD = 1. - s[cDirNG].S / s[cDir].S, gC = 1. - s[cConvNG].S / s[cConv].S;
        const double simNG = (s[cDirNG].S + s[cConvNG].S) * kToRows, simNGE = std::sqrt(s[cDirNG].V + s[cConvNG].V) * kToRows;
        const double Rb = datS[ip] / simNG, eRb = Rb * std::sqrt(std::pow(datSE[ip] / datS[ip], 2) + std::pow(simNGE / simNG, 2));
        P("%-13s %7.4f +- %-6.4f %7.4f +- %-6.4f %7.4f +- %-6.4f %5.3f+-%.3f %7.4f +- %-6.4f %7.0f +- %-5.0f %5.0f +- %-4.0f +- %-5.0f %7.4f +- %-6.4f %5.3f+-%.3f+-%.3f %5.3f+-%.3f %6.3f / %-6.3f %7.5f %12.4f\n",
          kPresetShort[ip], D / kLsim, eD / kLsim, C / kLsim, eC / kLsim, s[cTot].S / kLsim, std::sqrt(s[cTot].V) / kLsim, f, ef,
          w[cTot].S / kLsim, std::sqrt(w[cTot].V) / kLsim, simTot[ip], simTotE[ip], datS[ip], datSE[ip], datSM[ip], datS[ip] / kLdata,
          datSE[ip] / kLdata, R, eR, eRm, Rb, eRb, gD, gC, s[cX17].S / kLsim, s[cDirKW].S / kLsim);
    }
    P("\n# f_dir in the peak window 0.5336-0.5542:");
    for (int ip = 0; ip < kNP; ++ip) P("  %s %.3f+-%.3f", kPresetShort[ip], fdirW[ip], fdirWE[ip]);
    P("\n# share of the simulated peak content inside the window: direct %.3f, conversion %.3f, total %.3f (pre)\n",
      compWin[cDir].S / comp0[cDir].S, compWin[cConv].S / comp0[cConv].S, compWin[cTot].S / comp0[cTot].S);
    // efficiency of OA > 2 relative to pre
    {
        const EffRes eD = Efficiency(d, SelPreset(kPre), SelPreset(kPreOA2), CompPeakDirect());
        const EffRes eC = Efficiency(d, SelPreset(kPre), SelPreset(kPreOA2), CompPeakConv());
        const EffRes eT = Efficiency(d, SelPreset(kPre), SelPreset(kPreOA2), CompPeakTotal());
        const Eff ed = PassFail(df[1], df[6]);
        P("# efficiency of OA > 2 relative to pre: sim direct %.3f+-%.3f (generator M(ee) model not QED), conversion %.3f+-%.3f, total %.3f+-%.3f; data (pass/fail fits, this macro) %.3f +- %.3f (stat) +- %.3f (model, %d variants);"
          " data (assistant, today) 0.778\n", eD.eps, eD.err, eC.eps, eC.err, eT.eps, eT.err, ed.e, ed.estat, ed.esys, ed.nv);
        P("# data pass/fail for OA > 2: S_pass %.0f +- %.0f +- %.0f, S_fail %.0f +- %.0f +- %.0f; S(pre) %.0f +- %.0f +- %.0f (stat, model)\n", df[1].S, df[1].eS, df[1].eM,
          df[6].S, df[6].eS, df[6].eM, df[0].S, df[0].eS, df[0].eM);
    }
    P("# data, nominal WITHOUT vz cut (reproduction of scan_peak.md 4319 +- 348): %.0f +- %.0f (stat) +- %.0f (model; data study: 19-26 %% of S)\n", df[5].S, df[5].eS, df[5].eM);
    P("# model spread relative to S:");
    for (int ip = 0; ip < kNP; ++ip) P("  %s %.2f", kPresetShort[ip], datSM[ip] / datS[ip]);
    P("\n");
    CloseTab();

    // ---------------- 1c. classes
    OpenTab("eta_peak_classes.txt");
    P("# eta_peak_classes: peak content per class of OA(e+e-) [deg] and M(ee) [GeV/c2] (reconstructed), simulation against data\n");
    P("# simulation: expected data rows = sum w x L_data / L_sim (L_data = %.1f nb^-1, UNCONFIRMED), all M4, std + lep + eta owned rows, mode A, w (PDG)\n", kLdata);
    P("# data (this macro): P_main fit of N - CB per class (frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61), same preset and class\n");
    P("# data errors (review finding 4): total = stat (+) model; model = fit-model spread per class (pol 1-4, 4 ranges, width -+ 0.5 MeV, half range, p > 0.05)\n");
    P("#   every fit below uses the TOTAL error; the stat-only fit (as before the review) is printed for comparison only\n");
    P("# sim direct: generator M(ee) model not QED (review finding 2); 'dir KW' / 'fr KW' / KW fits: ILLUSTRATION ONLY (flat acceptance and fixed content assumed)\n");
    P("# frac = class content / sum over the classes; share_dir = direct / (direct + conversion) in the class; R = data / sim total\n");
    P("# 'given' = numbers quoted by the data study or today's summary (other selection where stated); X17 = BR 1e-5 illustration (rows), not in any sum\n");
    P("# mixture fit: data_i = a conv_i + b dir_i (a, b = data/SMASH scale of each component); one-factor fit: data_i = a1 (conv_i + dir_i)\n");
    struct ClassOut { std::vector<double> sd, sde, sc, sce, st, ste, sx, sk, ske, dS, dSe, dSs, dSm, gv, gve; MixRes mix, mixKW, mixStat; };
    std::vector<ClassOut> cout_(combos.size());
    for (size_t k = 0; k < combos.size(); ++k) {
        const ClassSet cs = FindCS(combos[k].cs);
        const int p = combos[k].preset == -1 ? kPeakNom : combos[k].preset;
        ClassOut& o = cout_[k];
        const int n = (int)cs.lo.size();
        for (int i = 0; i < n; ++i) {
            const double lo = cs.lo[i], hi = cs.hi[i]; const int var = cs.var;
            const std::vector<WSum> s = SumMany(d, SelAnd(SelPreset(p), [lo, hi, var](const Row& r) { const double x = CVar(var, r); return x >= lo && x < hi; }), comps);
            o.sd.push_back(s[cDir].S * kToRows); o.sde.push_back(std::sqrt(s[cDir].V) * kToRows);
            o.sc.push_back(s[cConv].S * kToRows); o.sce.push_back(std::sqrt(s[cConv].V) * kToRows);
            o.st.push_back(s[cTot].S * kToRows); o.ste.push_back(std::sqrt(s[cTot].V) * kToRows);
            o.sx.push_back(s[cX17].S * kToRows);
            o.sk.push_back(s[cDirKW].S * kToRows); o.ske.push_back(std::sqrt(s[cDirKW].V) * kToRows);
            const DFit& f = df[cidx[k][i]];
            o.dS.push_back(f.S); o.dSe.push_back(f.eT()); o.dSs.push_back(f.eS); o.dSm.push_back(f.eM);
            o.gv.push_back(i < (int)cs.given.size() ? cs.given[i] : kNaN);
            o.gve.push_back(i < (int)cs.givenE.size() ? cs.givenE[i] : kNaN);
        }
        o.mix = MixFit(o.dS, o.dSe, o.sc, o.sce, o.sd, o.sde);
        o.mixKW = MixFit(o.dS, o.dSe, o.sc, o.sce, o.sk, o.ske);
        o.mixStat = MixFit(o.dS, o.dSs, o.sc, o.sce, o.sd, o.sde);
        double tD = 0, tC = 0, tT = 0, tS = 0, tG = 0, vS = 0, tK = 0;
        for (int i = 0; i < n; ++i) { tD += o.sd[i]; tC += o.sc[i]; tT += o.st[i]; tS += o.dS[i]; vS += o.dSe[i] * o.dSe[i]; tK += o.sk[i]; if (o.gv[i] == o.gv[i]) tG += o.gv[i]; }
        P("\n## %s, %s%s\n", combos[k].preset == -1 ? "peak nominal WITHOUT vz cut (data); simulation: peak nominal" : PresetName(p), cs.tag.c_str(),
          cs.var == 0 ? " (OA [deg])" : " (M(ee) [MeV/c2])");
        P("%-8s | %9s %9s %16s | %9s | %6s %6s %6s | %16s %13s | %13s | %17s %6s | %7s | %8s %6s | %13s\n", "class", "sim dir", "sim conv", "sim total", "share_dir",
          "fr dir", "fr cnv", "fr tot", "data S (tot)", "data frac", "R data/sim", "given", "g frac", "X17", "dir KW", "fr KW", "stat / model");
        for (int i = 0; i < n; ++i) {
            const double fd = o.sd[i] / tD, fc = o.sc[i] / tC, ft = o.st[i] / tT, fs = o.dS[i] / tS;
            const double efs = std::sqrt(std::pow((1 - fs) * o.dSe[i], 2) + fs * fs * (vS - o.dSe[i] * o.dSe[i])) / tS;
            const double R = o.dS[i] / o.st[i], eR = std::fabs(R) * std::sqrt(std::pow(o.dSe[i] / o.dS[i], 2) + std::pow(o.ste[i] / o.st[i], 2));
            const double sh = o.sd[i] / o.st[i];
            const std::string g = (o.gv[i] == o.gv[i]) ? (o.gve[i] == o.gve[i] ? Form("%5.0f +- %-4.0f", o.gv[i], o.gve[i]) : Form("%5.0f (no err)", o.gv[i])) : std::string("-");
            P("%-8s | %9.0f %9.0f %7.0f +- %-5.0f | %9.3f | %6.3f %6.3f %6.3f | %6.0f +- %-5.0f %6.3f+-%.3f | %5.3f+-%.3f | %17s %6.3f | %7.1f | %8.0f %6.3f | %5.0f / %-5.0f\n",
              cs.name[i].c_str(), o.sd[i], o.sc[i], o.st[i], o.ste[i], sh, fd, fc, ft, o.dS[i], o.dSe[i], fs, efs, R, eR, g.c_str(),
              (o.gv[i] == o.gv[i] && tG > 0) ? o.gv[i] / tG : kNaN, o.sx[i], o.sk[i], o.sk[i] / tK, o.dSs[i], o.dSm[i]);
        }
        const int it = combos[k].preset == -1 ? 5 : (int)(std::find(kPresets, kPresets + kNP, p) - kPresets);
        P("%-8s | %9.0f %9.0f %16.0f | %9.3f | %20s | %6.0f +- %-5.0f (whole selection fit %.0f +- %.0f) | given sum %.0f (%s)\n", "sum", tD, tC, tT, tD / tT, "",
          tS, std::sqrt(vS), df[it].S, df[it].eS, tG, cs.givenSrc.empty() ? "-" : cs.givenSrc.c_str());
        {
            std::string fl;
            for (int i = 0; i < n; ++i) {
                const DFit& f = df[cidx[k][i]];
                if (f.nVarPass <= 0) fl += Form(" %s(%s)", cs.name[i].c_str(), f.nVarPass < 0 ? "p > 1e-3 used" : "UNDETERMINED");
            }
            if (!fl.empty()) P("# model spread: fewer than two variants with p > 0.05 in:%s\n", fl.c_str());
        }
        const MixRes& m = o.mix;
        if (m.ok)
            P("mixture fit: a(conv) = %.3f +- %.3f, b(direct) = %.3f +- %.3f, rho %.2f, chi2/ndf %.1f/%d (p %.3f); one factor a1 = %.3f +- %.3f, chi2/ndf %.1f/%d (p %.3f);"
              " direct fraction of the data peak %.3f +- %.3f (simulation %.3f)\n",
              m.a, m.ea, m.b, m.eb, m.rho, m.chi2, m.ndf, TMath::Prob(m.chi2, std::max(1, m.ndf)), m.a1, m.ea1, m.chi21, m.ndf1,
              TMath::Prob(m.chi21, std::max(1, m.ndf1)), m.fdir, m.efdir, m.fdirSim);
        if (m.ok && (m.a < 0 || m.b < 0)) P("   -> a negative factor (total errors): the two simulated shapes do not combine physically into the data classes\n");
        const MixRes& ms = o.mixStat;
        if (ms.ok)
            P("   stat-only errors (as before the review, for comparison): a(conv) = %.3f +- %.3f, b(direct) = %.3f +- %.3f, chi2/ndf %.1f/%d; one factor %.3f +- %.3f, chi2/ndf %.1f/%d\n",
              ms.a, ms.ea, ms.b, ms.eb, ms.chi2, ms.ndf, ms.a1, ms.ea1, ms.chi21, ms.ndf1);
        const MixRes& q = o.mixKW;
        if (q.ok)
            P("ILLUSTRATION ONLY (flat acceptance and fixed content assumed, not a result); mixture fit with the KW-reweighted direct decay: a(conv) = %.3f +- %.3f, b(direct KW) = %.3f +- %.3f, rho %.2f, chi2/ndf %.1f/%d (p %.3f); one factor %.3f +- %.3f, chi2/ndf %.1f/%d (p %.3f);"
              " direct fraction of the data peak %.3f +- %.3f (simulation %.3f)\n",
              q.a, q.ea, q.b, q.eb, q.rho, q.chi2, q.ndf, TMath::Prob(q.chi2, std::max(1, q.ndf)), q.a1, q.ea1, q.chi21, q.ndf1,
              TMath::Prob(q.chi21, std::max(1, q.ndf1)), q.fdir, q.efdir, q.fdirSim);
    }
    P("\n# KW model check (pre, class 5): bins of the TRUE M(ee) [MeV], model fraction, accepted fraction, factor model/accepted; norm %.4f; model Gamma(pipiee)/Gamma(pipigamma) = %.3e (PDG ratio 6.35e-3)\n",
      kw.norm, kw.ratioBR);
    for (size_t k = 0; k + 1 < kw.e.size(); ++k)
        P("#   %7.2f - %-7.2f  model %.4f  accepted %.4f  factor %.3f\n", kw.e[k] * 1e3, kw.e[k + 1] * 1e3, kw.model[k], kw.acc[k], kw.fac[k]);
    CloseTab();
    printf("Part1 tables done %.0f s\n", sw.RealTime()); sw.Continue();

    // ---------------- figures: composition
    gStyle->SetOptStat(0); gStyle->SetEndErrorSize(4);
    {
        const std::vector<std::string> lab = {"pre", "+OA>2", "slides", "old code", "nominal"};
        TCanvas* cv = new TCanvas(Nm("cComp").c_str(), "", 2700, 1000);
        cv->cd();
        TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.88, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
        Txt(0.01, 0.66, "#eta peak per selection: simulated direct decay (#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, eta set) and conversion (#eta#rightarrow#pi^{+}#pi^{-}#gamma, std), against the data peak content", 0.24);
        Txt(0.01, 0.38, Form("%s; all M_{4}; direct decay: generator M_{ee} model not QED (yield and share not for physics)", kSimNote), 0.18);
        Txt(0.01, 0.10, kDataNote, 0.18);
        std::vector<double> yD(kNP), eD(kNP), yC(kNP), eC(kNP), yT(kNP), eT(kNP), yS(kNP), eS(kNP), yF(kNP), eF(kNP), yFw(kNP), eFw(kNP), yR(kNP), eR(kNP);
        for (int ip = 0; ip < kNP; ++ip) {
            yD[ip] = simDir[ip]; eD[ip] = simDirE[ip]; yC[ip] = simConv[ip]; eC[ip] = simConvE[ip]; yT[ip] = simTot[ip]; eT[ip] = simTotE[ip];
            const double eTotD = std::hypot(datSE[ip], datSM[ip]);   // stat (+) fit-model spread (review finding 4)
            yS[ip] = datS[ip]; eS[ip] = eTotD; yF[ip] = fdir[ip]; eF[ip] = fdirE[ip]; yFw[ip] = fdirW[ip]; eFw[ip] = fdirWE[ip];
            yR[ip] = datS[ip] / simTot[ip]; eR[ip] = yR[ip] * std::sqrt(std::pow(eTotD / datS[ip], 2) + std::pow(simTotE[ip] / simTot[ip], 2));
        }
        // (a) expected rows
        cv->cd(); MkPad("pa", 0.0, 0.0, 0.36, 0.88, false, 0.15, 0.02, 0.06, 0.12);
        CatFrame(lab, 0, *std::max_element(yT.begin(), yT.end()) * 1.35, "selection", "peak content [rows at L_{data}]  (#sigma = rows / 6080.5 nb)", 0.045);
        TGraphErrors* gD = CatGraph(yD, eD, -0.24); StyleG(gD, kColDir, 21, 1.6); gD->Draw("P SAME");
        TGraphErrors* gC = CatGraph(yC, eC, -0.08); StyleG(gC, kColConv, 22, 1.8); gC->Draw("P SAME");
        TGraphErrors* gT = CatGraph(yT, eT, 0.08); StyleG(gT, kColTot, 33, 2.4); gT->Draw("P SAME");
        TGraphErrors* gS = CatGraph(yS, eS, 0.24); StyleG(gS, kBlack, 20, 1.6); gS->Draw("P SAME");
        TLegend* lg = new TLegend(0.17, 0.72, 0.97, 0.93); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036); lg->SetNColumns(2);
        lg->AddEntry(gD, "sim: direct #pi^{+}#pi^{-}e^{+}e^{-}", "pe"); lg->AddEntry(gC, "sim: conversion #pi^{+}#pi^{-}#gamma", "pe");
        lg->AddEntry(gT, "sim: total (direct + conv.)", "pe"); lg->AddEntry(gS, "data: P_{main} (stat #oplus model)", "pe");
        lg->Draw();
        Txt(0.17, 0.66, "(a) absolute: expected data rows (L_{data} UNCONFIRMED)", 0.030);
        // (b) direct fraction
        cv->cd(); MkPad("pb", 0.36, 0.0, 0.68, 0.88, false, 0.15, 0.02, 0.06, 0.12);
        CatFrame(lab, 0.3, 0.8, "selection", "direct / (direct + conversion), simulation", 0.045);
        TGraphErrors* gF = CatGraph(yF, eF, -0.1); StyleG(gF, kColTot, 20, 1.6); gF->Draw("P SAME");
        TGraphErrors* gFw = CatGraph(yFw, eFw, 0.1); StyleG(gFw, kColDir, 24, 1.6); gFw->Draw("P SAME");
        TLegend* lb = new TLegend(0.17, 0.80, 0.97, 0.93); lb->SetBorderSize(0); lb->SetFillStyle(0); lb->SetTextSize(0.038);
        lb->AddEntry(gF, "all M_{4}", "pe"); lb->AddEntry(gFw, "peak window 0.5336-0.5542 GeV/c^{2}", "pe"); lb->Draw();
        Txt(0.17, 0.16, "(b) BR(#pi^{+}#pi^{-}e^{+}e^{-}) = 2.68#times10^{-4} (PDG); conversion: std", 0.032);
        // (c) data / sim
        cv->cd(); MkPad("pc", 0.68, 0.0, 1.0, 0.88, false, 0.15, 0.02, 0.06, 0.12);
        CatFrame(lab, 0.0, 1.2, "selection", "data P_{main} / simulated content (all M_{4})", 0.045);
        { TLine* l1 = new TLine(0, 1, kNP, 1); l1->SetLineStyle(3); l1->Draw(); }
        TGraphErrors* gR = CatGraph(yR, eR, 0); StyleG(gR, kBlack, 20, 1.8); gR->Draw("P SAME");
        Txt(0.17, 0.88, "(c) METHOD A only: data P_{main} / sim, all M_{4}", 0.034);
        Txt(0.17, 0.78, "#color[2]{method dependent: see eta_summary_f_peak_methods}", 0.030);
        Txt(0.17, 0.73, "errors: stat #oplus fit model (no abs. scale #pm20%)", 0.030);
        Txt(0.17, 0.83, "L_{data} = 6080.5 nb^{-1} UNCONFIRMED", 0.032);
        Save(cv, "eta_peak_composition", true);
    }
    // ---------------- figures: classes (pre, OA5 and MEE5)
    {
        TCanvas* cv = new TCanvas(Nm("cCls").c_str(), "", 2700, 1800);
        cv->cd();
        TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.93, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
        Txt(0.01, 0.62, "#eta peak content per class of OA(e^{+}e^{-}) (top) and M(e^{+}e^{-}) (bottom), selection pre (isBest #geq 0, PT3, v_{z} > -160 mm): simulation against data", 0.30);
        Txt(0.01, 0.18, Form("%s; %s; data errors stat #oplus fit-model spread; sim direct decay: generator M_{ee} model not QED", "sim: std + lep + eta owned rows, w (PDG), all M_{4}, rows at L_{data}", kDataNote), 0.20);
        const int rows[2] = {0, 1};   // combos: pre OA5, pre MEE5
        for (int rr = 0; rr < 2; ++rr) {
            const int k = rows[rr];
            const ClassSet cs = FindCS(combos[k].cs);
            const ClassOut& o = cout_[k];
            const int n = (int)cs.lo.size();
            double tD = 0, tC = 0, tT = 0, tS = 0, vS = 0, tK = 0;
            for (int i = 0; i < n; ++i) { tD += o.sd[i]; tC += o.sc[i]; tT += o.st[i]; tS += o.dS[i]; vS += o.dSe[i] * o.dSe[i]; tK += o.sk[i]; }
            std::vector<double> fD(n), fDe(n), fC(n), fCe(n), fT(n), fTe(n), fS(n), fSe(n), R(n), Re(n), sC(n), sT(n), sTe(n), X(n), fK(n), fKe(n);
            for (int i = 0; i < n; ++i) {
                fD[i] = o.sd[i] / tD; fDe[i] = o.sde[i] / tD; fC[i] = o.sc[i] / tC; fCe[i] = o.sce[i] / tC; fT[i] = o.st[i] / tT; fTe[i] = o.ste[i] / tT;
                fK[i] = o.sk[i] / tK; fKe[i] = o.ske[i] / tK;
                fS[i] = o.dS[i] / tS; fSe[i] = std::sqrt(std::pow((1 - fS[i]) * o.dSe[i], 2) + fS[i] * fS[i] * (vS - o.dSe[i] * o.dSe[i])) / tS;
                R[i] = o.dS[i] / o.st[i]; Re[i] = std::fabs(R[i]) * std::sqrt(std::pow(o.dSe[i] / o.dS[i], 2) + std::pow(o.ste[i] / o.st[i], 2));
                sC[i] = o.sc[i]; sT[i] = o.st[i]; sTe[i] = o.ste[i]; X[i] = o.sx[i];
            }
            const char* xt = cs.var == 0 ? "OA(e^{+}e^{-}) class [deg]" : "M(e^{+}e^{-}) class [MeV/c^{2}]";
            const double y1 = rr == 0 ? 0.465 : 0.0, y2 = rr == 0 ? 0.93 : 0.465;
            // (1) fractions
            cv->cd(); MkPad("cfa", 0.0, y1, 0.34, y2, false, 0.15, 0.02, 0.09, 0.15);
            double ym = 0; for (int i = 0; i < n; ++i) ym = std::max({ym, fD[i] + fDe[i], fC[i] + fCe[i], fS[i] + fSe[i], fK[i]});
            CatFrame(cs.name, 0, ym * 1.45, xt, "fraction of the peak content", 0.05);
            TH1D* hD = CatHist(fD, fDe); StyleH(hD, kColDir, 3, 1); hD->Draw("HIST SAME"); hD->Draw("E1 SAME");
            TH1D* hC = CatHist(fC, fCe); StyleH(hC, kColConv, 3, 1); hC->Draw("HIST SAME"); hC->Draw("E1 SAME");
            TH1D* hT = CatHist(fT, fTe); StyleH(hT, kColTot, 3, 2); hT->Draw("HIST SAME");
            TH1D* hK = CatHist(fK, fKe); StyleH(hK, kColDir, 3, 7); hK->Draw("HIST SAME");
            TGraphErrors* gS = CatGraph(fS, fSe, 0); StyleG(gS, kBlack, 20, 1.6); gS->Draw("P SAME");
            TLegend* lg = new TLegend(0.40, 0.56, 0.97, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
            lg->AddEntry(hD, "sim direct (gen. M_{ee} not QED)", "l"); lg->AddEntry(hC, "sim conversion", "l"); lg->AddEntry(hT, "sim total", "l");
            lg->AddEntry(gS, "data (fits, stat #oplus model)", "pe");
            lg->AddEntry(hK, "direct, QED M_{ee}: illustration (flat acc., fixed content)", "l");
            lg->Draw();
            Txt(0.17, 0.93, "(1) shapes: each sample normalised to 1", 0.045);
            // (2) absolute rows
            cv->cd(); MkPad("cfb", 0.34, y1, 0.67, y2, false, 0.17, 0.02, 0.09, 0.15);
            double ya = 0; for (int i = 0; i < n; ++i) ya = std::max({ya, sT[i] + sTe[i], o.dS[i] + o.dSe[i]});
            TH1D* fb = CatFrame(cs.name, 0, ya * 1.6, xt, "peak content [rows at L_{data}]", 0.05);
            fb->GetYaxis()->SetMaxDigits(5); fb->GetYaxis()->SetTitleOffset(1.6);
            TH1D* hTa = CatHist(sT, sTe); StyleH(hTa, kColDir, 1, 1, kColDir); hTa->Draw("HIST SAME");
            std::vector<double> z(n, 0.);
            TH1D* hCa = CatHist(sC, z); StyleH(hCa, kColConv, 1, 1, kColConv); hCa->Draw("HIST SAME");
            TH1D* hTe = CatHist(sT, sTe); StyleH(hTe, kBlack, 1, 1); hTe->SetMarkerSize(0); hTe->Draw("E1 SAME");
            TH1D* hX = CatHist(X, z); StyleH(hX, kColX17, 3, 2); hX->Draw("HIST SAME");
            TGraphErrors* gA = CatGraph(o.dS, o.dSe, 0); StyleG(gA, kBlack, 20, 1.6); gA->Draw("P SAME");
            TLegend* la = new TLegend(0.50, 0.62, 0.97, 0.89); la->SetBorderSize(0); la->SetFillStyle(0); la->SetTextSize(0.040);
            la->AddEntry(hTa, "sim direct (on top)", "f"); la->AddEntry(hCa, "sim conversion", "f"); la->AddEntry(gA, "data (fits)", "pe");
            la->AddEntry(hX, "X17, BR 10^{-5} (info)", "l"); la->Draw();
            Txt(0.17, 0.93, "(2) absolute, L_{data} = 6080.5 nb^{-1} UNCONFIRMED", 0.042);
            // (3) ratio
            cv->cd(); MkPad("cfc", 0.67, y1, 1.0, y2, false, 0.15, 0.02, 0.09, 0.15);
            double yr = 0; for (int i = 0; i < n; ++i) yr = std::max(yr, R[i] + Re[i]);
            CatFrame(cs.name, 0, std::max(1.0, yr * 1.6), xt, "data / sim total", 0.05);
            TGraphErrors* gR = CatGraph(R, Re, 0); StyleG(gR, kBlack, 20, 1.6); gR->Draw("P SAME");
            const MixRes& m = o.mix;
            TLine* l = new TLine(0, m.a1, n, m.a1); l->SetLineColor(kColTot); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
            TLine* l1 = new TLine(0, 1, n, 1); l1->SetLineStyle(3); l1->Draw();
            Txt(0.17, 0.93, "(3) data / simulated total per class", 0.045);
            Txt(0.17, 0.85, Form("one factor: %.3f #pm %.3f, #chi^{2}/ndf = %.1f/%d (stat only: %.0f/%d)", m.a1, m.ea1, m.chi21, m.ndf1, o.mixStat.chi21, o.mixStat.ndf1), 0.034, kColTot);
            Txt(0.17, 0.80, Form("free factors conv. / direct: %.2f #pm %.2f / %.2f #pm %.2f%s", m.a, m.ea, m.b, m.eb, m.a < 0 ? " (unphysical)" : ""), 0.034);
            Txt(0.17, 0.75, "data errors: stat #oplus fit-model spread; qualitative: data have far more peak at OA > 10#circ / M_{ee} > 20 MeV", 0.030);
        }
        Save(cv, "eta_peak_classes", true);
    }
    // ---------------- figure: finer classes and the quoted data numbers
    {
        TCanvas* cv = new TCanvas(Nm("cClsF").c_str(), "", 2400, 1800);
        cv->cd();
        TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.94, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
        Txt(0.01, 0.60, "#eta peak content per class, fractions: simulation (lines) against own data fits (full circles) and the numbers quoted by the data study (open circles)", 0.30);
        Txt(0.01, 0.15, "top: pre (v_{z} > -160 mm), quoted = today's summary (selection inferred: pre); bottom: peak nominal (#alpha_{lab} < 55#circ, M_{#pi#pi} > 0.30), quoted = scan_peak.md (no v_{z} cut); own fits there also without v_{z} cut", 0.22);
        const int ks[4] = {2, 3, 8, 9};   // pre OA8t, pre MEE5t, nominal(no vz) OA4sp, MEE5
        for (int q = 0; q < 4; ++q) {
            const int k = ks[q];
            const ClassSet cs = FindCS(combos[k].cs);
            const ClassOut& o = cout_[k];
            const int n = (int)cs.lo.size();
            double tD = 0, tC = 0, tT = 0, tS = 0, vS = 0, tG = 0, vG = 0, tK = 0;
            for (int i = 0; i < n; ++i) {
                tD += o.sd[i]; tC += o.sc[i]; tT += o.st[i]; tS += o.dS[i]; vS += o.dSe[i] * o.dSe[i]; tK += o.sk[i];
                tG += o.gv[i]; vG += (o.gve[i] == o.gve[i]) ? o.gve[i] * o.gve[i] : o.dSe[i] * o.dSe[i];
            }
            std::vector<double> fD(n), fDe(n), fC(n), fCe(n), fT(n), fTe(n), fS(n), fSe(n), fG(n), fGe(n), fK(n), fKe(n);
            for (int i = 0; i < n; ++i) {
                fD[i] = o.sd[i] / tD; fDe[i] = o.sde[i] / tD; fC[i] = o.sc[i] / tC; fCe[i] = o.sce[i] / tC; fT[i] = o.st[i] / tT; fTe[i] = o.ste[i] / tT;
                fS[i] = o.dS[i] / tS; fSe[i] = std::sqrt(std::pow((1 - fS[i]) * o.dSe[i], 2) + fS[i] * fS[i] * (vS - o.dSe[i] * o.dSe[i])) / tS;
                fK[i] = o.sk[i] / tK; fKe[i] = o.ske[i] / tK;
                const double ge = (o.gve[i] == o.gve[i]) ? o.gve[i] : kNaN;
                fG[i] = o.gv[i] / tG;
                fGe[i] = (ge == ge) ? std::sqrt(std::pow((1 - fG[i]) * ge, 2) + fG[i] * fG[i] * (vG - ge * ge)) / tG : 0.;
            }
            const char* xt = cs.var == 0 ? "OA(e^{+}e^{-}) class [deg]" : "M(e^{+}e^{-}) class [MeV/c^{2}]";
            cv->cd(); MkPad("cq", (q % 2) * 0.5, q < 2 ? 0.47 : 0.0, (q % 2) * 0.5 + 0.5, q < 2 ? 0.94 : 0.47, false, 0.12, 0.02, 0.09, 0.15);
            double ym = 0; for (int i = 0; i < n; ++i) ym = std::max({ym, fD[i] + fDe[i], fC[i] + fCe[i], fS[i] + fSe[i], fG[i] + fGe[i]});
            CatFrame(cs.name, 0, ym * 1.45, xt, "fraction of the peak content", 0.05);
            TH1D* hD = CatHist(fD, fDe); StyleH(hD, kColDir, 3, 1); hD->Draw("HIST SAME"); hD->Draw("E1 SAME");
            TH1D* hC = CatHist(fC, fCe); StyleH(hC, kColConv, 3, 1); hC->Draw("HIST SAME"); hC->Draw("E1 SAME");
            TH1D* hT = CatHist(fT, fTe); StyleH(hT, kColTot, 3, 2); hT->Draw("HIST SAME");
            TH1D* hK = CatHist(fK, fKe); StyleH(hK, kColDir, 3, 7); hK->Draw("HIST SAME");
            TGraphErrors* gS = CatGraph(fS, fSe, -0.1); StyleG(gS, kBlack, 20, 1.5); gS->Draw("P SAME");
            TGraphErrors* gG = CatGraph(fG, fGe, 0.1); StyleG(gG, kBlack, 24, 1.5); gG->Draw("P SAME");
            TLegend* lg = new TLegend(0.14, 0.70, 0.97, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040); lg->SetNColumns(3);
            lg->AddEntry(hD, "sim direct", "l"); lg->AddEntry(hC, "sim conversion", "l"); lg->AddEntry(hT, "sim total", "l");
            lg->AddEntry(gS, "data, own fits (stat #oplus model)", "pe"); lg->AddEntry(gG, "data, quoted", "pe"); lg->AddEntry(hK, "direct, QED M_{ee} (illustration)", "l"); lg->Draw();
            Txt(0.14, 0.93, Form("#bf{%s}, %s", combos[k].preset == -1 ? "peak nominal (data: no v_{z} cut)" : PresetName(combos[k].preset), cs.tag.c_str()), 0.045);
        }
        Save(cv, "eta_peak_classes_fine", false);
    }
    printf("Part1 done %.0f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 2: position, width, mass scale
// =====================================================================================================
struct GFit { double m = 0, em = 0, s = 0, es = 0, chi2 = 0; int ndf = 0; bool ok = false; };
GFit GaussIter(TH1D* h, double nsig, const char* opt = "QNR0", double m0 = 0.546, double s0 = 0.0065) {
    GFit o;
    TF1* fg = new TF1(Nm("fg").c_str(), "gaus", m0 - 2 * s0, m0 + 2 * s0);
    double mx = 0; for (int i = 1; i <= h->GetNbinsX(); ++i) mx = std::max(mx, h->GetBinContent(i));
    fg->SetParameters(mx, m0, s0);
    h->Fit(fg, opt);
    for (int it = 0; it < 8; ++it) {
        const double m = fg->GetParameter(1), s = std::fabs(fg->GetParameter(2));
        if (!(s > 0.0005 && s < 0.05)) break;
        fg->SetRange(m - nsig * s, m + nsig * s);
        h->Fit(fg, opt);
    }
    o.m = fg->GetParameter(1); o.em = fg->GetParError(1); o.s = std::fabs(fg->GetParameter(2)); o.es = fg->GetParError(2);
    o.chi2 = fg->GetChisquare(); o.ndf = fg->GetNDF(); o.ok = o.s > 0 && o.em > 0;
    return o;
}
void Part2(SimData& d, const scan::Data& dd) {
    TStopwatch sw;
    Reduce(d, false);
    OpenTab("eta_peak_position.txt");
    P("# eta_peak_position: position, width and mass scale of the simulated eta peak (direct, conversion, total) and of the data peak\n");
    P("# simulation: reconstructed corrected M4 of the class rows (std + lep + eta owned, mode A, w (PDG)), weighted histograms with event-group errors\n");
    P("# fits: Gaussian, chi2 with the bin errors (empty bins skipped) unless 'WL' (weighted likelihood); range mean +- n sigma iterated 8 times\n");
    P("# data: FitPeakFree of ScanEngine (Gaussian with free position and width + pol2 on N - CB, chi2 scan), same preset\n");
    P("# PDG m_eta = %.6f GeV/c2; K0S scale sim/data = %.4f +- %.4f (review: K0S M(pipi), sim 1.0018 / data 0.9970)\n\n", kMetaG, kK0sSD, kK0sSDe);
    const char* cn[3] = {"direct", "conversion", "total"};
    const std::vector<Comp> cc = {CompPeakDirect(), CompPeakConv(), CompPeakTotal()};
    const int col[3] = {kColDir, kColConv, kColTot};
    struct R { int c, ib; std::string meth; GFit g; };
    std::vector<R> res;
    const double bws[3] = {0.0005, 0.001, 0.002};
    const double nsg[4] = {1.5, 2.0, 2.5, 3.0};
    TH1D* hShow[3] = {nullptr, nullptr, nullptr};
    P("%-11s %-8s %-14s %22s %20s %12s\n", "component", "bin", "method", "mean [GeV/c2]", "sigma [MeV]", "chi2/ndf");
    for (int ib = 0; ib < 3; ++ib) {
        Binning b; b.lo = 0.49; b.hi = 0.60; b.nbin = (int)std::lround((b.hi - b.lo) / bws[ib]);
        const std::vector<WSpec> w = FillMany(d, SelPreset(kPre), cc, SimOpt(), b);
        for (int c = 0; c < 3; ++c) {
            TH1D* h = w[c].Hist(Nm("pk").c_str(), 0);
            if (ib == 1) hShow[c] = h;
            for (int ir = 0; ir < 4; ++ir) {
                const GFit g = GaussIter(h, nsg[ir]);
                res.push_back({c, ib, Form("gaus +-%.1f", nsg[ir]), g});
                P("%-11s %-8s %-14s %10.5f +- %-8.5f %8.2f +- %-7.2f %6.1f/%-4d\n", cn[c], Form("%.1f MeV", bws[ib] * 1e3), Form("gaus +-%.1f", nsg[ir]),
                  g.m, g.em, g.s * 1e3, g.es * 1e3, g.chi2, g.ndf);
            }
            const GFit gl = GaussIter(h, 2.0, "QNR0WL");
            res.push_back({c, ib, "gaus WL +-2", gl});
            P("%-11s %-8s %-14s %10.5f +- %-8.5f %8.2f +- %-7.2f %6s\n", cn[c], Form("%.1f MeV", bws[ib] * 1e3), "gaus WL +-2", gl.m, gl.em, gl.s * 1e3, gl.es * 1e3, "-");
            if (ib >= 1) {   // Gaussian + linear term for tails, 0.52-0.575
                TF1* fp = new TF1(Nm("fp").c_str(), "gaus(0)+pol1(3)", 0.52, 0.575);
                double mx = 0; for (int i = 1; i <= h->GetNbinsX(); ++i) mx = std::max(mx, h->GetBinContent(i));
                fp->SetParameters(mx, 0.546, 0.0065, 0, 0);
                h->Fit(fp, "QNR0");
                GFit g; g.m = fp->GetParameter(1); g.em = fp->GetParError(1); g.s = std::fabs(fp->GetParameter(2)); g.es = fp->GetParError(2);
                g.chi2 = fp->GetChisquare(); g.ndf = fp->GetNDF(); g.ok = true;
                res.push_back({c, ib, "gaus+pol1", g});
                P("%-11s %-8s %-14s %10.5f +- %-8.5f %8.2f +- %-7.2f %6.1f/%-4d\n", cn[c], Form("%.1f MeV", bws[ib] * 1e3), "gaus+pol1 .52-.575",
                  g.m, g.em, g.s * 1e3, g.es * 1e3, g.chi2, g.ndf);
            }
        }
    }
    // non-fit estimators: weighted median, half of the 16-84 % interval, truncated mean / rms (+-12 MeV around the median)
    P("\n# estimators without a fit (pre, rows with 0.49 <= M4 < 0.60): median, (q84 - q16)/2, mean and rms within +-12 MeV of the median; error of the mean = rms / sqrt(N_eff)\n");
    P("%-11s %12s %12s %20s %12s %10s\n", "component", "median", "(q84-q16)/2", "trunc. mean", "trunc. rms", "N_eff");
    double medTot = 0;
    std::vector<std::pair<double, double>> vals[3], dlt[3], rat[3], tru[3];
    std::vector<std::pair<double, double>> trkR[3][4], trkRrec[3][4];
    for (int set = 1; set <= 3; ++set)
        for (const Row& r : d.tree(set, kUS)) {
            if (!Owned(r, SimOpt()) || !PassPreset(kPre, r)) continue;
            const int c = r.cls == 5 ? 0 : (r.cls == 4 ? 1 : -1);
            if (c < 0) continue;
            for (int cc2 : {c, 2}) {
                if (r.m4 >= 0.49f && r.m4 < 0.60f) vals[cc2].push_back({r.m4, r.w});
                if (r.m4_true > 0) {
                    dlt[cc2].push_back({(r.m4 - r.m4_true) * 1e3, r.w});
                    rat[cc2].push_back({r.m4 / r.m4_true, r.w});
                    tru[cc2].push_back({r.m4_true, r.w});
                }
                const float pc[4] = {r.ep_p, r.em_p, r.pip_p, r.pim_p}, pr[4] = {r.ep_p_rec, r.em_p_rec, r.pip_p_rec, r.pim_p_rec};
                const float pt[4] = {r.ep_sim_p, r.em_sim_p, r.pip_sim_p, r.pim_sim_p};
                for (int t = 0; t < 4; ++t)
                    if (pt[t] > 0) { trkR[cc2][t].push_back({pc[t] / pt[t], r.w}); trkRrec[cc2][t].push_back({pr[t] / pt[t], r.w}); }
            }
        }
    double estM[3], estMe[3], estS[3];
    for (int c = 0; c < 3; ++c) {
        const double med = WQuant(vals[c], 0.5), hw = 0.5 * (WQuant(vals[c], 0.84) - WQuant(vals[c], 0.16));
        double sw0 = 0, sw1 = 0, sw2 = 0, sww = 0;
        for (const auto& p : vals[c]) if (std::fabs(p.first - med) < 0.012) { sw0 += p.second; sw1 += p.second * p.first; sw2 += p.second * p.first * p.first; sww += p.second * p.second; }
        const double mean = sw1 / sw0, rms = std::sqrt(std::max(0., sw2 / sw0 - mean * mean)), neff = sw0 * sw0 / sww;
        estM[c] = mean; estMe[c] = rms / std::sqrt(neff); estS[c] = rms;
        if (c == 2) medTot = med;
        P("%-11s %12.5f %10.2f MeV %10.5f +- %-7.5f %8.2f MeV %10.0f\n", cn[c], med, hw * 1e3, mean, rms / std::sqrt(neff), rms * 1e3, neff);
    }
    (void)medTot;
    // reco vs true
    P("\n# reconstructed (corrected) against true M4 (pre, all rows of the class): Delta = M4 - M4_true; ratio = M4 / M4_true\n");
    P("%-11s %14s %14s %14s %22s %22s %16s\n", "component", "median Delta", "median ratio", "(q84-q16)/2", "gaus core Delta [MeV]", "gaus sigma [MeV]", "median M4_true");
    double medRatio[3] = {0, 0, 0};
    TH1D* hDel[3];
    for (int c = 0; c < 3; ++c) {
        const double md = WQuant(dlt[c], 0.5), mr = WQuant(rat[c], 0.5), hw = 0.5 * (WQuant(dlt[c], 0.84) - WQuant(dlt[c], 0.16));
        medRatio[c] = mr;
        hDel[c] = new TH1D(Nm("hdel").c_str(), "", 160, -40, 40); hDel[c]->SetDirectory(nullptr); hDel[c]->Sumw2();
        for (const auto& p : dlt[c]) hDel[c]->Fill(p.first, p.second);
        const GFit g = GaussIter(hDel[c], 2.0, "QNR0", 1.5, 6.0);
        P("%-11s %10.2f MeV %14.5f %10.2f MeV %10.2f +- %-8.2f %10.2f +- %-8.2f %16.5f\n", cn[c], md, mr, hw, g.m, g.em, g.s, g.es, WQuant(tru[c], 0.5));
    }
    P("\n# per-track momentum: median of p / p_true (corrected p) and of p_rec / p_true (reconstructed, before the energy-loss correction), pre\n");
    P("%-11s %-6s %10s %10s\n", "component", "track", "corr/true", "rec/true");
    const char* tn[4] = {"e+", "e-", "pi+", "pi-"};
    for (int c = 0; c < 3; ++c)
        for (int t = 0; t < 4; ++t) P("%-11s %-6s %10.5f %10.5f\n", cn[c], tn[t], WQuant(trkR[c][t], 0.5), WQuant(trkRrec[c][t], 0.5));
    printf("Part2 sim done %.0f s\n", sw.RealTime()); sw.Continue();

    // other selections (total and components, 1 MeV, +-2 sigma)
    P("\n# selection dependence: Gaussian +-2 sigma, 1 MeV bins, and data FitPeakFree (P_main default), same preset\n");
    P("%-13s %-11s %22s %20s\n", "selection", "sample", "mean [GeV/c2]", "sigma [MeV]");
    struct DPos { double m, em, s, es, S, eS, chi2; int ndf; bool ok; };
    std::vector<DPos> dpos;
    std::vector<DPos> spos;   // sim total per preset
    for (int ip = 0; ip < kNP; ++ip) {
        Binning b; b.lo = 0.49; b.hi = 0.60; b.nbin = 110;
        const std::vector<WSpec> w = FillMany(d, SelPreset(kPresets[ip]), cc, SimOpt(), b);
        for (int c = 0; c < 3; ++c) {
            TH1D* h = w[c].Hist(Nm("pkp").c_str(), 0);
            const GFit g = GaussIter(h, 2.0);
            P("%-13s %-11s %10.5f +- %-8.5f %8.2f +- %-7.2f\n", kPresetShort[ip], cn[c], g.m, g.em, g.s * 1e3, g.es * 1e3);
            if (c == 2) spos.push_back({g.m, g.em, g.s, g.es, 0, 0, 0, 0, g.ok});
        }
        const scan::Spectra sp = scan::Fill(dd, DataSelPreset(kPresets[ip]));
        const scan::FreeFit f = scan::FitPeakFree(sp);
        const scan::FreeFit fpo = scan::FitPeakFree(sp, scan::PeakCfg(), true);
        dpos.push_back({f.mean, f.eMean, f.sigma, f.eSigma, f.S, f.eS, f.chi2, f.ndf, f.ok});
        P("%-13s %-11s %10.5f +- %-8.5f %8.2f +- %-7.2f  (S %.0f +- %.0f, chi2/ndf %.1f/%d; Poisson: %.5f +- %.5f, %.2f +- %.2f MeV)\n", kPresetShort[ip], "DATA",
          f.mean, f.eMean, f.sigma * 1e3, f.eSigma * 1e3, f.S, f.eS, f.chi2, f.ndf, fpo.mean, fpo.eMean, fpo.sigma * 1e3, fpo.eSigma * 1e3);
    }
    // data variants (pre): ranges x orders
    const scan::Spectra spPre = scan::Fill(dd, DataSelPreset(kPre));
    std::vector<double> dm, ds;
    P("\n# data peak (pre), FitPeakFree variants: fit range x polynomial order\n");
    for (int ir = 0; ir < scan::kNVarRange; ++ir)
        for (int ord = 1; ord <= 3; ++ord) {
            scan::PeakCfg c; c.lo = scan::kVarLo[ir]; c.hi = scan::kVarHi[ir]; c.order = ord;
            const scan::FreeFit f = scan::FitPeakFree(spPre, c);
            if (!f.ok) { P("range %.3f-%.3f order %d: no minimum\n", c.lo, c.hi, ord); continue; }
            dm.push_back(f.mean); ds.push_back(f.sigma);
            P("range %.3f-%.3f order %d: mean %.5f +- %.5f, sigma %.2f +- %.2f MeV, S %.0f +- %.0f, chi2/ndf %.1f/%d\n", c.lo, c.hi, ord, f.mean, f.eMean,
              f.sigma * 1e3, f.eSigma * 1e3, f.S, f.eS, f.chi2, f.ndf);
        }
    auto spread = [](std::vector<double> v) { if (v.empty()) return 0.; return 0.5 * (*std::max_element(v.begin(), v.end()) - *std::min_element(v.begin(), v.end())); };
    // summary over the standard fit choices (1 and 2 MeV bins, gaus +-1.5..3 sigma, WL, gaus+pol1)
    P("\n# SUMMARY (pre). Simulation: median over the fits with 1 and 2 MeV bins (gaus +-1.5..3 sigma, WL, gaus+pol1); 'fit' = half range over them\n");
    double sMed[3], sMedE[3], sSpr[3], wMed[3], wMedE[3], wSpr[3];
    for (int c = 0; c < 3; ++c) {
        std::vector<double> m, s, em, es;
        for (const R& r : res) if (r.c == c && r.ib >= 1 && r.g.ok) { m.push_back(r.g.m); s.push_back(r.g.s); em.push_back(r.g.em); es.push_back(r.g.es); }
        auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); const size_t n = v.size(); return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]); };
        sMed[c] = med(m); sMedE[c] = med(em); sSpr[c] = spread(m); wMed[c] = med(s); wMedE[c] = med(es); wSpr[c] = spread(s);
        P("sim %-11s position %.5f +- %.5f (stat) +- %.5f (fit) GeV/c2 = m_eta x %.5f; width %.2f +- %.2f (stat) +- %.2f (fit) MeV\n", cn[c], sMed[c], sMedE[c],
          sSpr[c], sMed[c] / kMetaG, wMed[c] * 1e3, wMedE[c] * 1e3, wSpr[c] * 1e3);
    }
    const DPos& D0 = dpos[0];
    P("data (pre)        position %.5f +- %.5f (stat) +- %.5f (fit) GeV/c2 = m_eta x %.5f; width %.2f +- %.2f (stat) +- %.2f (fit) MeV\n", D0.m, D0.em, spread(dm),
      D0.m / kMetaG, D0.s * 1e3, D0.es * 1e3, spread(ds) * 1e3);
    P("data study (oa_lab < 50, no vz cut): 0.5439, 5.4 MeV\n");
    for (int c = 0; c < 3; ++c) {
        const double r = sMed[c] / D0.m, er = r * std::sqrt(std::pow(sMedE[c] / sMed[c], 2) + std::pow(D0.em / D0.m, 2));
        const double rf = r * std::sqrt(std::pow(sSpr[c] / sMed[c], 2) + std::pow(spread(dm) / D0.m, 2));
        P("scale sim %-10s / data: %.5f +- %.5f (stat) +- %.5f (fit); offset %+.2f +- %.2f MeV; K0S sim/data %.4f +- %.4f -> expected offset %+.2f MeV; after K0S scale: %+.2f +- %.2f MeV\n",
          cn[c], r, er, rf, (sMed[c] - D0.m) * 1e3, std::sqrt(sMedE[c] * sMedE[c] + D0.em * D0.em) * 1e3, kK0sSD, kK0sSDe, (kK0sSD - 1) * D0.m * 1e3,
          (sMed[c] / kK0sSD - D0.m) * 1e3, std::sqrt(sMedE[c] * sMedE[c] + D0.em * D0.em) * 1e3);
    }
    P("reco/true (median of M4/M4_true): direct %.5f, conversion %.5f, total %.5f; data position / m_eta %.5f\n", medRatio[0], medRatio[1], medRatio[2], D0.m / kMetaG);
    P("width: sim total %.2f +- %.2f +- %.2f MeV against data %.2f +- %.2f +- %.2f MeV; difference %.2f +- %.2f (stat) MeV\n", wMed[2] * 1e3, wMedE[2] * 1e3,
      wSpr[2] * 1e3, D0.s * 1e3, D0.es * 1e3, spread(ds) * 1e3, (wMed[2] - D0.s) * 1e3, std::sqrt(wMedE[2] * wMedE[2] + D0.es * D0.es) * 1e3);
    CloseTab();

    // ---------------- figure
    gStyle->SetOptStat(0); gStyle->SetEndErrorSize(3);
    TCanvas* cv = new TCanvas(Nm("cPos").c_str(), "", 2600, 1900);
    cv->cd();
    TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.93, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
    Txt(0.01, 0.62, "#eta peak position, width and mass scale: simulation (direct, conversion, total) against the data peak, selection pre", 0.30);
    Txt(0.01, 0.18, "simulation: reconstructed corrected M_{4}, std + lep + eta owned rows, w (PDG); data: N - CB minus the fitted pol2 (FitPeakFree); K0S: M(#pi#pi) scale sim/data = 1.0048 #pm 0.0007", 0.22);
    // (a) shapes: sim total, direct, conversion (density), data background-subtracted
    {
        cv->cd(); MkPad("pa", 0.0, 0.465, 0.5, 0.93, false, 0.13, 0.02, 0.08, 0.13);
        const double lo = scan::kFitLo, hi = scan::kFitHi; const int nb = 35;   // 2.5 MeV
        Binning b; b.lo = lo; b.hi = hi; b.nbin = nb;
        const std::vector<WSpec> w = FillMany(d, SelPreset(kPre), cc, SimOpt(), b);
        Binning bs = b; bs.x = [](const Row& r) { return (double)r.m4 / kK0sSD; };
        const WSpec wsh = FillMany(d, SelPreset(kPre), {CompPeakTotal()}, SimOpt(), bs)[0];
        TH1D* hs[4];
        for (int c = 0; c < 3; ++c) { hs[c] = w[c].Hist(Nm("sh").c_str(), 0); hs[c]->Scale(1. / (hs[c]->Integral() * b.bw())); }
        hs[3] = wsh.Hist(Nm("sh").c_str(), 0); hs[3]->Scale(1. / (hs[3]->Integral() * b.bw()));
        // data
        scan::PeakCfg pc; pc.mean = D0.m; pc.sigma = D0.s;
        const scan::PeakFit pf = scan::FitPeakMain(spPre, pc);
        TH1D* hd = new TH1D(Nm("hdsub").c_str(), "", nb, lo, hi); hd->SetDirectory(nullptr);
        double sum = 0;
        for (int i = 0; i < spPre.nbin; ++i) {
            const double x = spPre.x(i);
            if (x < lo || x >= hi) continue;
            const int j = hd->FindBin(x);
            hd->SetBinContent(j, spPre.sig[i] - pf.Bkg(x)); hd->SetBinError(j, std::sqrt(spPre.sigVar[i]));
            sum += spPre.sig[i] - pf.Bkg(x);
        }
        hd->Scale(1. / (sum * b.bw()));
        double ym = 0; for (int i = 1; i <= nb; ++i) ym = std::max({ym, hd->GetBinContent(i) + hd->GetBinError(i), hs[0]->GetBinContent(i), hs[1]->GetBinContent(i)});
        Frame(lo, hi, -0.15 * ym, ym * 1.3, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "1/N dN/dM [(GeV/c^{2})^{-1}]", 0.045);
        StyleH(hs[0], kColDir, 3, 1); hs[0]->Draw("HIST SAME");
        StyleH(hs[1], kColConv, 3, 1); hs[1]->Draw("HIST SAME");
        StyleH(hs[2], kColTot, 4, 1); hs[2]->Draw("HIST SAME"); hs[2]->Draw("E1 SAME");
        StyleH(hs[3], kColTot, 3, 2); hs[3]->Draw("HIST SAME");
        hd->SetMarkerStyle(20); hd->SetMarkerColor(kBlue); hd->SetLineColor(kBlue); hd->SetMarkerSize(1.3); hd->Draw("E1 SAME");
        TLine* l0 = new TLine(lo, 0, hi, 0); l0->SetLineStyle(3); l0->Draw();
        TLine* lm = new TLine(kMetaG, -0.15 * ym, kMetaG, ym * 1.3); lm->SetLineStyle(3); lm->SetLineColor(kGray + 2); lm->Draw();
        TLegend* lg = new TLegend(0.55, 0.55, 0.97, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
        lg->AddEntry(hs[0], "sim direct", "l"); lg->AddEntry(hs[1], "sim conversion", "l"); lg->AddEntry(hs[2], "sim total", "l");
        lg->AddEntry(hs[3], "sim total, M_{4} / 1.0048 (K0S scale)", "l"); lg->AddEntry(hd, "data N - CB - pol2 (pre)", "pe");
        lg->AddEntry(lm, "PDG m_{#eta}", "l"); lg->Draw();
        Txt(0.15, 0.93, "(a) shapes, each normalised to 1 in 0.5225-0.61 GeV/c^{2}, 2.5 MeV bins", 0.040);
    }
    // (b) position per fit choice
    {
        cv->cd(); MkPad("pb", 0.5, 0.465, 1.0, 0.93, false, 0.13, 0.02, 0.08, 0.13);
        const int nfit = 6 * 3;   // per bin width: 4 gaus + WL + gaus+pol1
        Frame(-0.5, nfit - 0.5, 0.5415, 0.5505, "fit choice (bins 0.5 / 1 / 2 MeV; gaus #pm1.5, 2, 2.5, 3 #sigma, WL, gaus+pol1)", "peak position [GeV/c^{2}]", 0.040);
        TBox* bx = new TBox(-0.5, D0.m - D0.em, nfit - 0.5, D0.m + D0.em); bx->SetFillColorAlpha(kBlue, 0.15); bx->Draw();
        TLine* ld = new TLine(-0.5, D0.m, nfit - 0.5, D0.m); ld->SetLineColor(kBlue); ld->SetLineWidth(2); ld->Draw();
        TLine* lk = new TLine(-0.5, D0.m * kK0sSD, nfit - 0.5, D0.m * kK0sSD); lk->SetLineColor(kBlue); lk->SetLineWidth(2); lk->SetLineStyle(2); lk->Draw();
        TLine* le = new TLine(-0.5, kMetaG, nfit - 0.5, kMetaG); le->SetLineColor(kGray + 2); le->SetLineStyle(3); le->Draw();
        TGraphErrors* g[3];
        for (int c = 0; c < 3; ++c) {
            g[c] = new TGraphErrors();
            int k = 0;
            for (int ib = 0; ib < 3; ++ib) {
                int j = 0;
                for (const R& r : res) {
                    if (r.c != c || r.ib != ib) continue;
                    const double x = ib * 6 + j + (c - 1) * 0.22;
                    g[c]->SetPoint(k, x, r.g.m); g[c]->SetPointError(k, 0, r.g.em); ++k; ++j;
                }
            }
            StyleG(g[c], col[c], c == 0 ? 21 : (c == 1 ? 22 : 33), c == 2 ? 2.0 : 1.4); g[c]->Draw("P SAME");
        }
        TLegend* lg = new TLegend(0.15, 0.76, 0.97, 0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034); lg->SetNColumns(3);
        lg->AddEntry(g[0], "sim direct", "pe"); lg->AddEntry(g[1], "sim conversion", "pe"); lg->AddEntry(g[2], "sim total", "pe");
        lg->AddEntry(ld, "data (pre) #pm stat", "l"); lg->AddEntry(lk, "data #times 1.0048 (K0S)", "l"); lg->AddEntry(le, "PDG m_{#eta}", "l"); lg->Draw();
        Txt(0.15, 0.93, "(b) position for each fit choice", 0.040);
    }
    // (c) width per fit choice
    {
        cv->cd(); MkPad("pc", 0.0, 0.0, 0.5, 0.465, false, 0.13, 0.02, 0.08, 0.13);
        const int nfit = 18;
        Frame(-0.5, nfit - 0.5, 2.0, 11.0, "fit choice (bins 0.5 / 1 / 2 MeV; gaus #pm1.5, 2, 2.5, 3 #sigma, WL, gaus+pol1)", "peak width #sigma [MeV/c^{2}]", 0.040);
        TBox* bx = new TBox(-0.5, (D0.s - D0.es) * 1e3, nfit - 0.5, (D0.s + D0.es) * 1e3); bx->SetFillColorAlpha(kBlue, 0.15); bx->Draw();
        TLine* ld = new TLine(-0.5, D0.s * 1e3, nfit - 0.5, D0.s * 1e3); ld->SetLineColor(kBlue); ld->SetLineWidth(2); ld->Draw();
        TGraphErrors* g[3];
        for (int c = 0; c < 3; ++c) {
            g[c] = new TGraphErrors();
            int k = 0;
            for (int ib = 0; ib < 3; ++ib) {
                int j = 0;
                for (const R& r : res) {
                    if (r.c != c || r.ib != ib) continue;
                    g[c]->SetPoint(k, ib * 6 + j + (c - 1) * 0.22, r.g.s * 1e3); g[c]->SetPointError(k, 0, r.g.es * 1e3); ++k; ++j;
                }
            }
            StyleG(g[c], col[c], c == 0 ? 21 : (c == 1 ? 22 : 33), c == 2 ? 2.0 : 1.4); g[c]->Draw("P SAME");
        }
        Txt(0.15, 0.93, "(c) width for each fit choice; blue band: data (pre) #pm stat", 0.040);
    }
    // (d) Delta = M4 - M4_true
    {
        cv->cd(); MkPad("pd", 0.5, 0.0, 1.0, 0.465, true, 0.13, 0.02, 0.08, 0.13);
        double ym = 0;
        hDel[1]->Rebin(4);
        for (int c = 0; c < 2; ++c) { hDel[c]->Scale(1. / (hDel[c]->Integral() * hDel[c]->GetBinWidth(1))); ym = std::max(ym, hDel[c]->GetMaximum()); }
        Frame(-40, 40, ym * 2e-4, ym * 4, "M_{4} - M_{4}^{true} [MeV/c^{2}]", "1/N dN/d#Delta [(MeV/c^{2})^{-1}]", 0.040);
        StyleH(hDel[0], kColDir, 3, 1); hDel[0]->Draw("HIST SAME");
        StyleH(hDel[1], kColConv, 3, 1); hDel[1]->Draw("HIST SAME");
        TLine* l0 = new TLine(0, ym * 2e-4, 0, ym * 4); l0->SetLineStyle(3); l0->Draw();
        TLegend* lg = new TLegend(0.15, 0.74, 0.60, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
        lg->AddEntry(hDel[0], Form("direct (0.5 MeV bins), median M_{4}/M_{4}^{true} %.4f", medRatio[0]), "l");
        lg->AddEntry(hDel[1], Form("conversion (2 MeV bins), median %.4f", medRatio[1]), "l"); lg->Draw();
        Txt(0.15, 0.93, "(d) reconstructed minus true M_{4} (simulation, pre)", 0.040);
    }
    Save(cv, "eta_peak_position", true);
    printf("Part2 done %.0f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 3: cut efficiencies
// =====================================================================================================
struct OldD { double e, err; };
std::map<std::string, OldD> ParseOld() {
    std::map<std::string, OldD> m;
    std::ifstream in((std::string(kTabDir) + "/cuts_compare_data.txt").c_str());
    std::string line;
    auto trim = [](std::string s) { const size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' '); return a == std::string::npos ? std::string() : s.substr(a, b - a + 1); };
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.size() < 40) continue;
        const size_t bar = line.find('|');
        if (bar == std::string::npos || bar < 30) continue;
        const std::string fam = trim(line.substr(0, 9)), val = trim(line.substr(10, 11)), base = trim(line.substr(22, bar - 22));
        const std::string rest = line.substr(bar + 1);
        double e = 0, err = 0;
        if (sscanf(rest.c_str(), " %lf +- %lf", &e, &err) == 2) m[fam + "|" + val + "|" + base] = {e, err};
    }
    // scan_extra.md start (oa_lab < 50, M(pipi) < 0.46), quoted in the same table
    m["vtx|both pions|md"] = {0.912, 0.045};
    m["vtx|all four|md"] = {0.771, 0.049};
    m["vz|-160|md"] = {1.008, 0.029};
    m["dz_pi|40|md"] = {0.974, 0.038};
    return m;
}
void Part3(SimData& d, const scan::Data& dd) {
    TStopwatch sw;
    Reduce(d, false);
    const KWTab kw = MakeKW(d);
    const std::vector<Fam> fams = Families();
    const std::map<std::string, OldD> old = ParseOld();
    printf("Part3: %zu data efficiencies parsed\n", old.size());
    const std::vector<Comp> cc = {CompPeakDirect(), CompPeakConv(), CompPeakTotal(), CompDirectKW(kw), CompTotalKW(kw)};
    // data pass/fail at base pre (pre0 for vz), all families except dz_pi (no track z in the data compact file)
    struct Item { int f, i; int ip = -1, ifl = -1; };
    std::vector<Item> items;
    std::vector<scan::Sel> dsel;
    for (size_t f = 0; f < fams.size(); ++f) {
        if (fams[f].kind == kDz) continue;
        const Cut base = fams[f].pre0 ? CPre0() : CPre();
        for (size_t i = 0; i < fams[f].v.size(); ++i) {
            Cut cp = base; SetCut(cp, fams[f].kind, fams[f].v[i]);
            Item it; it.f = (int)f; it.i = (int)i;
            it.ip = (int)dsel.size(); dsel.push_back(DatSel(cp));
            it.ifl = (int)dsel.size(); dsel.push_back([base, cp](const scan::Row& r) { return PassBase(r) && PassCut(r, base) && !PassCut(r, cp); });
            items.push_back(it);
        }
    }
    const std::vector<DFit> df = FitDataSels(dd, dsel);
    std::map<std::pair<int, int>, Eff> dOwn;
    for (const Item& it : items) dOwn[{it.f, it.i}] = PassFail(df[it.ip], df[it.ifl]);
    printf("Part3: %zu data fits %.0f s\n", df.size(), sw.RealTime()); sw.Continue();

    OpenTab("eta_peak_eff.txt");
    P("# eta_peak_eff: efficiency of each cut alone for the simulated eta peak: direct (eta -> pi+pi-e+e-), conversion (eta -> pi+pi-gamma), total\n");
    P("# simulation: std + lep + eta owned rows, mode A, w (PDG), all M4; eps = S_pass / S_base, event-group errors (pass and fail rows of a group correlated)\n");
    P("# base 'pre' = isBest >= 0, PT3, vz > -160, z > -500, cor_ok (vz family: 'pre0' = without vz > -160); 'start' = pre0 + oa_lab < 55 + M(pipi) > 0.30\n");
    P("#   (scan_extra.C peak start); 'md' = pre0 + oa_lab < 50 + M(pipi) < 0.46 (scan_extra.md)\n");
    P("# data(old): data-study pass/fail peak fits as parsed in tables/cuts_compare_data.txt (1D scans: base 'no cut' = NO vz cut; start families: start base)\n");
    P("# data(this): pass/fail P_main fits of N - CB (frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61) at the SAME base as the simulation (pre / pre0), this macro\n");
    P("#   CAUTION: the data study flagged M(pipi) < 0.40-0.50 and MM4 < 2.60-2.90 as SCULPTING (the pass fraction turns on inside the fit range): data values biased there\n");
    P("# pulls: (sim - data) / sqrt(err_sim^2 + err_data^2), data(this) unless only data(old) exists\n\n");
    P("# dir KW / tot KW / p_tKW: ILLUSTRATION ONLY - direct decay reweighted to the QED (Kroll-Wada) M(ee) shape, flat acceptance and fixed content assumed; not a result\n");
    P("# sim direct: generator M(ee) model not QED (review finding 2); its OA efficiencies are not for physics\n");
    P("# data(this) error = stat (+) model since the review (model = half range of the efficiency over the fit-model variants where pass and fail fits\n");
    P("#   both have p > 0.05); the pulls use it; data(old) errors are those of the data-study table\n");
    P("%-9s %-11s %-6s | %15s %15s %15s | %15s %15s | %15s %15s | %6s %6s %6s %6s\n", "family", "value", "base", "sim direct", "sim conv", "sim total",
      "dir KW", "tot KW", "data(old)", "data(this)", "p_dir", "p_conv", "p_tot", "p_tKW");
    struct Pt { double x, sd, sde, sc, sce, st, ste, sk, ske, tk, tke, dold, dolde, dthis, dthise; bool hasOld, hasThis; };
    std::vector<std::vector<Pt>> pts(fams.size()), ptsStart(fams.size());
    std::vector<double> pullTot, pullConv, pullDir, pullTK;
    std::map<std::string, std::vector<double>> famChi;   // family -> sum pull^2 (dir, conv, tot, totKW), points
    std::vector<std::string> effFlags;                     // data points whose model error used a fallback (review finding 4)
    for (size_t f = 0; f < fams.size(); ++f) {
        const Fam& F = fams[f];
        for (int pass = 0; pass < 2; ++pass) {
            if (pass == 1 && !F.start) continue;
            const Cut base = pass == 1 ? CPeakStart() : (F.pre0 ? CPre0() : CPre());
            const std::string bname = pass == 1 ? "start" : (F.pre0 ? "pre0" : "pre");
            const std::string oldBase = pass == 1 ? "start" : (F.kind == kVz ? "start" : "no cut");
            for (size_t i = 0; i < F.v.size(); ++i) {
                Cut cp = base; SetCut(cp, F.kind, F.v[i]);
                EffRes e[5];
                for (int c = 0; c < 5; ++c) e[c] = Efficiency(d, SimSel(base), SimSel(cp), cc[c]);
                Pt p; p.x = XPlot(F, F.v[i]);
                p.sd = e[0].eps; p.sde = e[0].err; p.sc = e[1].eps; p.sce = e[1].err; p.st = e[2].eps; p.ste = e[2].err;
                p.sk = e[3].eps; p.ske = e[3].err; p.tk = e[4].eps; p.tke = e[4].err;
                const std::string key = F.key + "|" + ValTxt(F.kind, F.v[i]) + "|" + oldBase;
                auto io = old.find(key);
                p.hasOld = (io != old.end()) && !(F.kind == kVz && pass == 0 && false);
                // the vz family of the old table is at the start base: compare it with the simulation at the start base only
                if (F.kind == kVz && pass == 0) p.hasOld = false;
                if (p.hasOld) { p.dold = io->second.e; p.dolde = io->second.err; } else { p.dold = p.dolde = 0; }
                p.hasThis = false; p.dthis = p.dthise = 0;
                if (pass == 0) {
                    auto it = dOwn.find({(int)f, (int)i});
                    if (it != dOwn.end() && it->second.ok) {
                        p.hasThis = true; p.dthis = it->second.e; p.dthise = it->second.err;
                        if (it->second.mode != 0) effFlags.push_back(Form("%s %s: model error from %s (%d variants)", F.key.c_str(), ValTxt(F.kind, F.v[i]).c_str(),
                                                                         it->second.mode == 1 ? "variants with an acceptable PASS fit only" :
                                                                         it->second.mode == 2 ? "variants with an acceptable FAIL fit only" : "no acceptable variant: UNDETERMINED", it->second.nv));
                    }
                }
                const bool none = (F.kind == kVtx ? F.v[i] == 0 : (F.v[i] == kNo || ((F.kind == kOaLo || F.kind == kPt4) && F.v[i] == 0)));
                double dv = 0, de = 0; bool hd = false;
                if (p.hasThis && !none) { dv = p.dthis; de = p.dthise; hd = true; }
                else if (p.hasOld && !none) { dv = p.dold; de = p.dolde; hd = true; }
                auto pull = [&](double s, double se) { return hd ? (s - dv) / std::sqrt(se * se + de * de) : 0.; };
                const double pd = pull(p.sd, p.sde), pc = pull(p.sc, p.sce), pt = pull(p.st, p.ste), ptk = pull(p.tk, p.tke);
                if (hd && !(F.kind == kMppHi || F.kind == kMm4)) { pullTot.push_back(pt); pullConv.push_back(pc); pullDir.push_back(pd); pullTK.push_back(ptk); }
                if (hd) {
                    std::vector<double>& fc = famChi[F.key + " (" + bname + ")"];
                    if (fc.empty()) fc.assign(5, 0.);
                    fc[0] += pd * pd; fc[1] += pc * pc; fc[2] += pt * pt; fc[3] += ptk * ptk; fc[4] += 1;
                }
                const std::string so = p.hasOld ? Form("%6.3f+-%-6.3f", p.dold, p.dolde) : "-";
                const std::string st = p.hasThis ? Form("%6.3f+-%-6.3f", p.dthis, p.dthise) : "-";
                P("%-9s %-11s %-6s | %6.3f+-%-6.3f %6.3f+-%-6.3f %6.3f+-%-6.3f | %6.3f+-%-6.3f %6.3f+-%-6.3f | %15s %15s | %+6.1f %+6.1f %+6.1f %+6.1f\n", F.key.c_str(),
                  ValTxt(F.kind, F.v[i]).c_str(), bname.c_str(), p.sd, p.sde, p.sc, p.sce, p.st, p.ste, p.sk, p.ske, p.tk, p.tke, so.c_str(), st.c_str(), pd, pc, pt, ptk);
                (pass == 1 ? ptsStart[f] : pts[f]).push_back(p);
            }
        }
    }
    // scan_extra.md start
    P("\n# scan_extra.md start (pre0 + oa_lab < 50 + M(pipi) < 0.46)\n");
    {
        struct Q { const char* name; int kind; double v; const char* key; };
        const Q qs[] = {{"both pions in vertex fit", kVtx, 2, "vtx|both pions|md"}, {"all four in vertex fit", kVtx, 3, "vtx|all four|md"},
                        {"vz > -160", kVz, -160, "vz|-160|md"}, {"|z_pi - vz| < 40", kDz, 40, "dz_pi|40|md"}};
        for (const Q& q : qs) {
            const Cut base = CPeakStartMD(); Cut cp = base; SetCut(cp, q.kind, q.v);
            EffRes e[5];
            for (int c = 0; c < 5; ++c) e[c] = Efficiency(d, SimSel(base), SimSel(cp), cc[c]);
            const OldD o = old.at(q.key);
            P("%-26s | sim direct %.3f+-%.3f  conv %.3f+-%.3f  total %.3f+-%.3f | data %.3f+-%.3f | pulls %+5.1f %+5.1f %+5.1f\n", q.name, e[0].eps, e[0].err,
              e[1].eps, e[1].err, e[2].eps, e[2].err, o.e, o.err, (e[0].eps - o.e) / std::hypot(e[0].err, o.err),
              (e[1].eps - o.e) / std::hypot(e[1].err, o.err), (e[2].eps - o.e) / std::hypot(e[2].err, o.err));
        }
    }
    auto chi = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return s; };
    P("\n# data(this) points whose model error could not use variants with BOTH fits acceptable (nominal value always inside the range):\n");
    for (const std::string& x : effFlags) P("#   %s\n", x.c_str());
    P("\n# sum of pull^2 over the cut values with data (correlated points of one family; M(pipi) upper and MM4 upper excluded because of sculpting):\n");
    P("#   direct %.1f, conversion %.1f, total %.1f, total with KW direct %.1f for %zu points\n", chi(pullDir), chi(pullConv), chi(pullTot), chi(pullTK), pullTot.size());
    P("# per family (all points with data, incl. sculpted families): sum pull^2 direct / conversion / total / total KW, points\n");
    for (const auto& kv : famChi) P("#   %-22s %7.1f %7.1f %7.1f %7.1f  (%.0f)\n", kv.first.c_str(), kv.second[0], kv.second[1], kv.second[2], kv.second[3], kv.second[4]);
    CloseTab();
    printf("Part3 table %.0f s\n", sw.RealTime()); sw.Continue();

    // ---------------- figure
    gStyle->SetOptStat(0); gStyle->SetEndErrorSize(2);
    TCanvas* cv = new TCanvas(Nm("cEff").c_str(), "", 3000, 2300);
    cv->cd();
    TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.94, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
    Txt(0.01, 0.62, "Efficiency of each cut alone for the #eta peak: simulated direct decay, conversion and total, against the data-derived peak efficiencies (pass/fail fits)", 0.30);
    Txt(0.01, 0.18, "sim: std + lep + eta owned rows, w (PDG), all M_{4}, base pre (v_{z} family: no v_{z} cut); dashed lines and open squares: start base (#alpha_{lab} < 55#circ, M_{#pi#pi} > 0.30, no v_{z} cut); data(this): same base as sim; M_{#pi#pi} upper and MM_{4} upper: data sculpted", 0.21);
    const int ncol = 4, nrow = 3;
    for (size_t f = 0; f < fams.size(); ++f) {
        const Fam& F = fams[f];
        const int ic = f % ncol, ir = f / ncol;
        cv->cd();
        MkPad("pe", ic * 0.25, 0.94 - (ir + 1) * 0.313, (ic + 1) * 0.25, 0.94 - ir * 0.313, false, 0.14, 0.02, 0.07, 0.15);
        double xlo = 1e9, xhi = -1e9;
        for (const Pt& p : pts[f]) { xlo = std::min(xlo, p.x); xhi = std::max(xhi, p.x); }
        const double dx = (xhi - xlo) * 0.06 + 1e-6;
        TH1D* fr = Frame(xlo - dx, xhi + dx, 0, 1.15, F.title.c_str(), "efficiency", 0.055);
        if (F.kind == kVtx) {
            fr->GetXaxis()->SetNdivisions(4);
        }
        TLine* l1 = new TLine(xlo - dx, 1, xhi + dx, 1); l1->SetLineStyle(3); l1->Draw();
        auto mkG = [&](const std::vector<Pt>& v, int which) {
            TGraphErrors* g = new TGraphErrors();
            for (const Pt& p : v) {
                double y = 0, e = 0; bool ok = true;
                if (which == 0) { y = p.sd; e = p.sde; } else if (which == 1) { y = p.sc; e = p.sce; } else if (which == 2) { y = p.st; e = p.ste; }
                else if (which == 5) { y = p.tk; e = p.tke; }
                else if (which == 3) { ok = p.hasOld; y = p.dold; e = p.dolde; } else { ok = p.hasThis; y = p.dthis; e = p.dthise; }
                if (!ok) continue;
                const int k = g->GetN(); g->SetPoint(k, p.x, y); g->SetPointError(k, 0, e);
            }
            g->Sort();
            return g;
        };
        TGraphErrors* gD = mkG(pts[f], 0); StyleG(gD, kColDir, 1, 0, 3); gD->SetFillColorAlpha(kColDir, 0.25); gD->Draw("L3 SAME");
        TGraphErrors* gC = mkG(pts[f], 1); StyleG(gC, kColConv, 1, 0, 3); gC->SetFillColorAlpha(kColConv, 0.25); gC->Draw("L3 SAME");
        TGraphErrors* gT = mkG(pts[f], 2); StyleG(gT, kColTot, 1, 0, 4); gT->SetLineStyle(1); gT->Draw("L SAME");
        TGraphErrors* gK = mkG(pts[f], 5); StyleG(gK, kColTot, 1, 0, 3); gK->SetLineStyle(7); gK->SetLineColor(kMagenta + 2); gK->Draw("L SAME");
        if (!ptsStart[f].empty()) {
            TGraphErrors* sT = mkG(ptsStart[f], 2); StyleG(sT, kColTot, 1, 0, 3); sT->SetLineStyle(3); sT->Draw("L SAME");
            TGraphErrors* sO = mkG(ptsStart[f], 3); StyleG(sO, kBlack, 25, 1.3, 2); if (sO->GetN()) sO->Draw("P SAME");
        }
        TGraphErrors* gO = mkG(pts[f], 3); StyleG(gO, kGray + 2, 24, 1.3, 2); if (gO->GetN()) gO->Draw("P SAME");
        TGraphErrors* gH = mkG(pts[f], 4); StyleG(gH, kBlack, 20, 1.2, 2); if (gH->GetN()) gH->Draw("P SAME");
        if (F.kind == kMppHi || F.kind == kMm4) Txt(0.17, 0.20, "data sculpted (pass fraction turns on in the fit range)", 0.038, kRed + 1);
        if (F.kind == kDz) Txt(0.17, 0.20, "data(this) n/a: no track z in the data compact file", 0.038, kGray + 2);
        Txt(0.16, 0.95, Form("#bf{%s}", F.key.c_str()), 0.06);
        if (f == 0) {
            TLegend* lg = new TLegend(0.16, 0.14, 0.97, 0.66); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.046);
            lg->AddEntry(gD, "sim direct #pi^{+}#pi^{-}e^{+}e^{-}", "lf"); lg->AddEntry(gC, "sim conversion #pi^{+}#pi^{-}#gamma", "lf");
            lg->AddEntry(gT, "sim total", "l"); lg->AddEntry(gK, "sim total, QED-M_{ee} direct (illustration)", "l"); lg->AddEntry(gH, "data (this), stat #oplus model", "pe");
            lg->AddEntry(gO, "data (old table), no v_{z} cut", "pe");
            lg->AddEntry(Dummy(kColTot, 3, 3), "sim total, start base", "l");
            lg->AddEntry(Dummy(kBlack, 2, 1, 0, 25, 1.3), "data (old), start base", "p");
            lg->Draw();
        }
    }
    // last panel: text
    cv->cd();
    MkPad("pt", 3 * 0.25, 0.0, 1.0, 0.313, false, 0.05, 0.02, 0.07, 0.15);
    Txt(0.05, 0.86, "#bf{reading}", 0.065);
    Txt(0.05, 0.76, "efficiency relative to the base, all M_{4}", 0.048);
    Txt(0.05, 0.67, "data: S_{pass}/(S_{pass}+S_{fail}), P_{main} fits of N - CB", 0.048);
    Txt(0.05, 0.58, "(Gaussian frozen 0.5439 / 5.4 MeV)", 0.048);
    Txt(0.05, 0.49, Form("#Sigma pull^{2}, %zu points (no M_{#pi#pi}^{hi}, MM_{4}):", pullTot.size()), 0.048);
    Txt(0.05, 0.40, Form("direct %.0f, conversion %.0f, total %.0f", chi(pullDir), chi(pullConv), chi(pullTot)), 0.048);
    Txt(0.05, 0.31, Form("QED-M_{ee} direct (illustration only): %.0f", chi(pullTK)), 0.048);
    Txt(0.05, 0.22, "(points correlated; data errors stat #oplus fit model)", 0.048);
    Txt(0.05, 0.13, "sim direct: generator M_{ee} model not QED", 0.048);
    Save(cv, "eta_peak_eff", true);
    printf("Part3 done %.0f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 4: shapes
// =====================================================================================================
void Part4(SimData& d) {
    TStopwatch sw;
    Reduce(d, true);
    struct V { const char* key; const char* title; int nb; double lo, hi; bool logy; Var x; };
    const std::vector<V> vars = {
        {"mpipi", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", 30, 0.27, 0.42, false, [](const Row& r) { return (double)r.mpipi; }},
        {"mee", "M_{e^{+}e^{-}} [GeV/c^{2}]", 70, 0.0, 0.28, true, [](const Row& r) { return (double)r.mee; }},
        {"oa", "OA(e^{+}e^{-}) [deg]", 40, 0.0, 40.0, true, [](const Row& r) { return (double)r.oa; }},
        {"p_pipi", "p_{#pi^{+}#pi^{-}} [GeV/c]", 25, 0.0, 2.5, false, [](const Row& r) { return 1e-3 * r.p_pipi; }},
        {"p_ee", "p_{e^{+}e^{-}} [GeV/c]", 30, 0.0, 1.5, false, [](const Row& r) { return 1e-3 * r.p_ee; }},
        {"oa_lab", "#alpha_{lab}(#pi#pi, ee) [deg]", 30, 0.0, 90.0, false, [](const Row& r) { return (double)r.oa_lab; }}};
    const std::vector<Comp> cc = {CompPeakDirect(), CompPeakConv(), CompX17(1e-5)};
    OpenTab("eta_peak_shapes.txt");
    P("# eta_peak_shapes: distributions of the simulated eta peak components, selection pre (isBest >= 0, PT3, vz > -160), all M4\n");
    P("# std + lep + eta owned rows, mode A, w (PDG); X17 = BR(eta -> pi+pi- X17, X17 -> e+e-) = 1e-5, ILLUSTRATION only, not in any sum\n");
    P("# per variable: weighted mean, median, 10 %% and 90 %% quantiles (corrected momenta, reconstructed quantities)\n\n");
    P("%-8s %-12s %10s %10s %10s %10s\n", "variable", "component", "mean", "median", "q10", "q90");
    std::vector<std::vector<TH1D*>> h(vars.size());
    const char* cn[3] = {"direct", "conversion", "X17 1e-5"};
    for (size_t k = 0; k < vars.size(); ++k) {
        Binning b; b.nbin = vars[k].nb; b.lo = vars[k].lo; b.hi = vars[k].hi; b.x = vars[k].x;
        const std::vector<WSpec> w = FillMany(d, SelPreset(kPre), cc, SimOpt(), b);
        for (int c = 0; c < 3; ++c) h[k].push_back(w[c].Hist(Nm("sh").c_str(), 2));
        // quantiles
        std::vector<std::pair<double, double>> q[3];
        for (int set = 1; set <= 3; ++set)
            for (const Row& r : d.tree(set, kUS)) {
                if (!Owned(r, SimOpt()) || !PassPreset(kPre, r)) continue;
                const int c = r.cls == 5 ? 0 : r.cls == 4 ? 1 : r.cls == 6 ? 2 : -1;
                if (c < 0) continue;
                const double x = vars[k].x(r);
                if (x == x) q[c].push_back({x, c == 2 ? 2e-5 : r.w});
            }
        for (int c = 0; c < 3; ++c) {
            double s0 = 0, s1 = 0; for (const auto& p : q[c]) { s0 += p.second; s1 += p.second * p.first; }
            P("%-8s %-12s %10.4f %10.4f %10.4f %10.4f\n", vars[k].key, cn[c], s1 / s0, WQuant(q[c], 0.5), WQuant(q[c], 0.1), WQuant(q[c], 0.9));
        }
    }
    // fractions in the cut-relevant regions
    P("\n# fraction of each component (pre, all M4) passing: OA > 2, M(ee) > 20 MeV, M(pipi) > 0.30, M(pipi) < 0.46, oa_lab < 50, oa_eta > 140\n");
    {
        const std::vector<std::pair<const char*, Sel>> cuts = {
            {"OA > 2", [](const Row& r) { return r.oa > 2.f; }}, {"M(ee) > 20 MeV", [](const Row& r) { return r.mee > 0.020f; }},
            {"M(pipi) > 0.30", [](const Row& r) { return r.mpipi > 0.30f; }}, {"M(pipi) < 0.46", [](const Row& r) { return r.mpipi < 0.46f; }},
            {"oa_lab < 50", [](const Row& r) { return r.oa_lab < 50.f; }}, {"oa_eta > 140", [](const Row& r) { return r.oa_eta > 140.f; }}};
        for (const auto& ct : cuts) {
            P("%-16s", ct.first);
            for (int c = 0; c < 3; ++c) {
                const EffRes e = Efficiency(d, SelPreset(kPre), SelAnd(SelPreset(kPre), ct.second), cc[c]);
                P("  %s %.3f+-%.3f", cn[c], e.eps, e.err);
            }
            P("\n");
        }
    }
    // X17 in M4
    Binning bm; bm.nbin = 60; bm.lo = 0.50; bm.hi = 0.60;
    const std::vector<WSpec> wm = FillMany(d, SelPreset(kPre), cc, SimOpt(), bm);
    TH1D* hm[3]; for (int c = 0; c < 3; ++c) hm[c] = wm[c].Hist(Nm("hm4").c_str(), 2);
    {
        const std::vector<WSum> s = SumMany(d, SelPreset(kPre), cc);
        P("\n# absolute (pre, all M4): direct %.4f nb, conversion %.4f nb, X17 (BR 1e-5) %.5f nb = %.0f expected data rows (L_data %.1f nb^-1 UNCONFIRMED); X17/direct per unit BR = %.3f\n",
          s[0].S / kLsim, s[1].S / kLsim, s[2].S / kLsim, s[2].S * kToRows, kLdata, (s[2].S / 1e-5) / (s[0].S / 2.68e-4));
    }
    CloseTab();

    gStyle->SetOptStat(0); gStyle->SetEndErrorSize(0);
    TCanvas* cv = new TCanvas(Nm("cSh").c_str(), "", 3000, 1650);
    cv->cd();
    TPad* ph = new TPad(Nm("ph").c_str(), "", 0, 0.92, 1, 1); ph->SetMargin(0, 0, 0, 0); ph->Draw(); ph->cd();
    Txt(0.01, 0.62, "Simulated #eta peak components: direct decay #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (eta set) against conversion #eta#rightarrow#pi^{+}#pi^{-}#gamma (std), and the X17 illustration", 0.28);
    Txt(0.01, 0.18, "selection pre (isBest #geq 0, PT3, v_{z} > -160 mm), all M_{4}, std + lep + eta owned rows, w (PDG); shapes normalised to 1; last two panels absolute d#sigma/dx [nb/unit]; X17: BR = 10^{-5} illustration, not in any sum", 0.21);
    for (size_t k = 0; k < vars.size() + 2; ++k) {
        const int ic = k % 4, ir = k / 4;
        cv->cd();
        MkPad("ps", ic * 0.25, ir == 0 ? 0.46 : 0.0, (ic + 1) * 0.25, ir == 0 ? 0.92 : 0.46, k == 1 || k == 2 || k == 6, 0.18, 0.02, 0.08, 0.15);
        if (k < vars.size()) {
            TH1D* a = (TH1D*)h[k][0]->Clone(Nm("a").c_str()); a->Scale(1. / (a->Integral() * a->GetBinWidth(1)));
            TH1D* b = (TH1D*)h[k][1]->Clone(Nm("b").c_str()); b->Scale(1. / (b->Integral() * b->GetBinWidth(1)));
            TH1D* x = (TH1D*)h[k][2]->Clone(Nm("x").c_str()); x->Scale(1. / (x->Integral() * x->GetBinWidth(1)));
            const double ym = std::max({a->GetMaximum(), b->GetMaximum()});
            const bool lg = vars[k].logy;
            Frame(vars[k].lo, vars[k].hi, lg ? ym * 1e-4 : 0, lg ? ym * 5 : ym * 1.3, vars[k].title, "1/N dN/dx", 0.05);
            StyleH(a, kColDir, 3, 1); a->Draw("HIST SAME"); a->Draw("E1 SAME");
            StyleH(b, kColConv, 3, 1); b->Draw("HIST SAME"); b->Draw("E1 SAME");
            StyleH(x, kColX17, 2, 2);
            if (k == 1 || k == 2) x->Draw("HIST SAME");
            if (k == 0) {
                TLegend* l = new TLegend(0.20, 0.70, 0.75, 0.90); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.040);
                l->AddEntry(a, "direct #pi^{+}#pi^{-}e^{+}e^{-}", "l"); l->AddEntry(b, "conversion #pi^{+}#pi^{-}#gamma", "l");
                l->AddEntry(x, "X17 shape (M_{ee}, OA panels)", "l"); l->Draw();
            }
            Txt(0.17, 0.93, Form("shape: %s", vars[k].key), 0.05);
        } else if (k == vars.size()) {
            // X17 in M(ee), absolute
            TH1D* a = h[1][0]; TH1D* b = h[1][1]; TH1D* x = h[1][2];
            const double ym = std::max({a->GetMaximum(), b->GetMaximum(), x->GetMaximum()});
            Frame(0, 0.10, ym * 1e-4, ym * 5, "M_{e^{+}e^{-}} [GeV/c^{2}]", "d#sigma/dM_{ee} [nb/(GeV/c^{2})]", 0.05);
            gPad->SetLogy(1);
            StyleH(a, kColDir, 3, 1); a->Draw("HIST SAME");
            StyleH(b, kColConv, 3, 1); b->Draw("HIST SAME");
            StyleH(x, kColX17, 3, 2); x->Draw("HIST SAME");
            Txt(0.17, 0.93, "absolute: M_{ee} (X17 at BR 10^{-5}, info)", 0.05);
        } else {
            const double ym = std::max({hm[0]->GetMaximum(), hm[1]->GetMaximum(), hm[2]->GetMaximum()});
            Frame(0.50, 0.60, ym * 1e-3, ym * 5, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "d#sigma/dM [nb/(GeV/c^{2})]", 0.05);
            gPad->SetLogy(1);
            StyleH(hm[0], kColDir, 3, 1); hm[0]->Draw("HIST SAME");
            StyleH(hm[1], kColConv, 3, 1); hm[1]->Draw("HIST SAME");
            StyleH(hm[2], kColX17, 3, 2); hm[2]->Draw("HIST SAME");
            Txt(0.17, 0.93, "absolute: M_{4} (X17 at BR 10^{-5}, info)", 0.05);
        }
    }
    Save(cv, "eta_peak_shapes", true);

    // ---------------- true M(ee) of the accepted direct decay against the QED (Kroll-Wada) shape
    {
        const KWTab kw = MakeKW(d);
        OpenTab("eta_peak_mee_true.txt");
        P("# eta_peak_mee_true: TRUE M(ee) of the accepted simulated direct decay (class 5, pre, eta set) against a QED reference shape\n");
        P("# reference: Kroll-Wada factor x pi pi gamma* rate (P-wave pions, k^3), no form factor; integrated %.3e x Gamma(pi pi gamma) (PDG ratio 6.35e-3)\n", kw.ratioBR);
        P("# 'accepted' = fraction of the accepted rows (generator x acceptance); a ratio accepted/model far from a smooth function of order 1\n");
        P("# points to the generator's M(ee) distribution (acceptance assumed to vary slowly; UNVERIFIED: no generator-level rows available)\n");
        P("%-20s %10s %10s %14s\n", "true M(ee) [MeV]", "model", "accepted", "accepted/model");
        double fm = 0, fa = 0;
        for (size_t k = 0; k + 1 < kw.e.size(); ++k) {
            P("%8.2f - %-9.2f %10.4f %10.4f %14.3f\n", kw.e[k] * 1e3, kw.e[k + 1] * 1e3, kw.model[k], kw.acc[k], kw.acc[k] / kw.model[k]);
            if (kw.e[k] >= 0.0199) { fm += kw.model[k]; fa += kw.acc[k]; }
        }
        P("# fraction above about 20 MeV: model %.3f, accepted %.3f\n", fm, fa);
        // ---- review findings 2 and 3: threshold behaviour, OA > 2 pass share against the true M(ee), X17 acceptance
        {
            // fine bins near threshold: accepted density per MeV (QED vanishes at 2 m_e through beta = sqrt(1 - 4 m^2/q^2))
            std::vector<double> fe;
            for (int i = 0; i <= 40; ++i) fe.push_back(0.001022 * std::pow(0.270 / 0.001022, i / 40.));
            const int nf = (int)fe.size() - 1;
            std::vector<double> wAll(nf, 0.), wOA(nf, 0.), nAll(nf, 0.);
            double tot = 0;
            auto fbin = [&](double q) { for (int k = 0; k < nf; ++k) if (q >= fe[k] && q < fe[k + 1]) return k; return -1; };
            for (int set = 1; set <= 3; ++set)
                for (const Row& r : d.tree(set, kUS)) {
                    if (!(r.tmpl == 1 && r.cls == 5) || !Owned(r, SimOpt()) || !PassPreset(kPre, r)) continue;
                    tot += r.w;
                    const int k = fbin(r.mee_true);
                    if (k < 0) continue;
                    wAll[k] += r.w; nAll[k] += 1; if (r.oa > 2.f) wOA[k] += r.w;
                }
            P("\n## review findings 2/3: accepted direct decay (class 5, pre) in fine bins of the TRUE M(ee)\n");
            P("# density = accepted fraction per MeV; QED/KW density (generator level, no acceptance) for comparison; OA>2 = share of the bin passing OA > 2\n");
            P("%-20s %12s %12s %12s %10s %8s\n", "true M(ee) [MeV]", "acc frac", "acc /MeV", "QED /MeV", "OA>2 share", "rows");
            double kwN = 0; for (int j = 0; j < 4000; ++j) { const double q = 0.001022 + (0.27 - 0.001022) * (j + 0.5) / 4000; kwN += KWDens(q) * (0.27 - 0.001022) / 4000; }
            for (int k = 0; k < nf; ++k) {
                const double dq = (fe[k + 1] - fe[k]) * 1e3;
                double qk = 0; for (int j = 0; j < 50; ++j) { const double q = fe[k] + (fe[k + 1] - fe[k]) * (j + 0.5) / 50; qk += KWDens(q) * (fe[k + 1] - fe[k]) / 50; }
                P("%8.3f - %-9.3f %12.4f %12.4f %12.4f %10.3f %8.0f\n", fe[k] * 1e3, fe[k + 1] * 1e3, wAll[k] / tot, wAll[k] / tot / dq, qk / kwN / dq,
                  wAll[k] > 0 ? wOA[k] / wAll[k] : kNaN, nAll[k]);
            }
            // largest accepted density per MeV: in which bin
            int kmax = 0; for (int k = 0; k < nf; ++k) if (wAll[k] / (fe[k + 1] - fe[k]) > wAll[kmax] / (fe[kmax + 1] - fe[kmax])) kmax = k;
            P("# largest accepted density per MeV in bin %d (%.3f-%.3f MeV): %s\n", kmax, fe[kmax] * 1e3, fe[kmax + 1] * 1e3,
              kmax == 0 ? "the FIRST bin above threshold -> not QED (QED density -> 0 at 2 m_e); acceptance cannot create this" : "not the first bin");
            // power-law check: slope of ln(density) against ln(M) between 1.5 and 10 MeV
            double sx = 0, sy = 0, sxx = 0, sxy = 0; int nn = 0;
            for (int k = 0; k < nf; ++k) {
                const double qc = std::sqrt(fe[k] * fe[k + 1]);
                if (qc < 0.0015 || qc > 0.010 || !(wAll[k] > 0)) continue;
                const double x = std::log(qc), y = std::log(wAll[k] / (fe[k + 1] - fe[k]));
                sx += x; sy += y; sxx += x * x; sxy += x * y; ++nn;
            }
            const double slope = (nn * sxy - sx * sy) / (nn * sxx - sx * sx);
            P("# power-law slope of the accepted density dN/dM between 1.5 and 10 MeV: %.2f (1/M^2 law: -2; QED: about -1 well above threshold)\n", slope);
            // X17 against the direct decay: both forced channels have 50 %% per eta (mail_weights.md), so the row ratio = acceptance ratio
            double n5 = 0, n6 = 0, n5p = 0, n6p = 0;
            for (const Row& r : d.tree(3, kUS)) {
                if (r.tmpl != 1 || !Owned(r, SimOpt())) continue;
                const bool pre = PassPreset(kPre, r);
                if (r.cls == 5) { n5 += 1; if (pre) n5p += 1; }
                if (r.cls == 6) { n6 += 1; if (pre) n6p += 1; }
            }
            P("# accepted rows per generated forced decay, X17 (M(ee) = 17 MeV) / direct pi+pi-e+e- (generator M(ee)): all loaded rows %.0f / %.0f = %.3f; pre %.0f / %.0f = %.3f\n",
              n6, n5, n6 / n5, n6p, n5p, n6p / n5p);
            P("#   -> a QED-shaped direct decay (more pairs at 10-100 MeV, larger OA) would probably be accepted MORE often than the generator one\n");
            P("#      (UNVERIFIED): the KW variant's 'fixed content' assumption is arbitrary; it is an illustration of shapes only\n");
        }
        CloseTab();
        TCanvas* c2 = new TCanvas(Nm("cMt").c_str(), "", 2000, 900);
        c2->cd();
        TPad* ph2 = new TPad(Nm("ph").c_str(), "", 0, 0.88, 1, 1); ph2->SetMargin(0, 0, 0, 0); ph2->Draw(); ph2->cd();
        Txt(0.01, 0.62, "Simulated #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (eta set): true M(e^{+}e^{-}) of the accepted rows against the QED (Kroll-Wada) shape", 0.32);
        Txt(0.01, 0.18, "selection pre, class 5 rows; reference: Kroll-Wada #times #pi#pi#gamma* rate (P-wave #pi#pi, k^{3}), no form factor; both normalised to 1; acceptance NOT unfolded (no generator rows)", 0.22);
        const int nb = (int)kw.e.size() - 1;
        TH1D* hm = new TH1D(Nm("hmod").c_str(), "", nb, kw.e.data()); hm->SetDirectory(nullptr);
        TH1D* ha = new TH1D(Nm("hacc").c_str(), "", nb, kw.e.data()); ha->SetDirectory(nullptr);
        TH1D* hr = new TH1D(Nm("hrat").c_str(), "", nb, kw.e.data()); hr->SetDirectory(nullptr);
        for (int k = 0; k < nb; ++k) {
            const double lw = std::log(kw.e[k + 1] / kw.e[k]);
            hm->SetBinContent(k + 1, kw.model[k] / lw); ha->SetBinContent(k + 1, kw.acc[k] / lw);
            hr->SetBinContent(k + 1, kw.acc[k] / kw.model[k]);
        }
        c2->cd(); MkPad("pm1", 0.0, 0.0, 0.5, 0.88, true, 0.13, 0.03, 0.06, 0.13); gPad->SetLogx(1);
        Frame(0.001, 0.27, 0.005, 2.0, "true M_{e^{+}e^{-}} [GeV/c^{2}]", "dN/d ln M_{ee} (normalised to 1)", 0.045);
        StyleH(hm, kBlack, 3, 2); hm->Draw("HIST SAME");
        StyleH(ha, kColDir, 3, 1); ha->Draw("HIST SAME");
        TLegend* l1 = new TLegend(0.17, 0.16, 0.80, 0.30); l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.040);
        l1->AddEntry(ha, "simulation: accepted rows (gen. #times acc.)", "l"); l1->AddEntry(hm, "QED reference shape (generator level)", "l"); l1->Draw();
        c2->cd(); MkPad("pm2", 0.5, 0.0, 1.0, 0.88, true, 0.13, 0.03, 0.06, 0.13); gPad->SetLogx(1);
        Frame(0.001, 0.27, 0.05, 20, "true M_{e^{+}e^{-}} [GeV/c^{2}]", "accepted / QED reference", 0.045);
        StyleH(hr, kColDir, 3, 1); hr->Draw("HIST SAME");
        TLine* l0 = new TLine(0.001, 1, 0.27, 1); l0->SetLineStyle(3); l0->Draw();
        Txt(0.17, 0.88, Form("fraction above 20 MeV: QED %.2f, accepted simulation %.2f", fm, fa), 0.040);
        Save(c2, "eta_peak_mee_true", false);
    }
    printf("Part4 done %.0f s\n", sw.RealTime());
}

}  // namespace sep

void sim_eta_peak(int part = 0) {
    using namespace simeng;
    using namespace sep;
    if (part < 0) { printf("sim_eta_peak: compiled only\n"); return; }
    TStopwatch sw;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    const bool needData = part == 0 || part == 1 || part == 2 || part == 3;
    scan::Data dd;
    if (needData && !scan::Load(dd, scan::kCompact, 0.30, 1.10)) return;
    auto run = [&](int p) {
        SimData d;
        if (!Load(d, kCompact, kMaskUS)) return;
        if (p == 1) Part1(d, dd);
        if (p == 2) Part2(d, dd);
        if (p == 3) Part3(d, dd);
        if (p == 4) Part4(d);
    };
    if (part == 0) { for (int p = 1; p <= 4; ++p) run(p); }
    else run(part);
    printf("sim_eta_peak(%d): total %.0f s\n", part, sw.RealTime());
}
