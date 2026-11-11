// ============================================================================
// scan_peak.C  --  cut scan for the eta PEAK in M(pi+pi-e+e-)
// ============================================================================
// HADES p+p, T = 4.53 GeV. Target: the whole content of the peak at the eta
// mass (eta -> pi+pi-e+e- and eta -> pi+pi-gamma with a converted photon).
// Uses scan_compact.root through plots/ScanEngine.h (frozen peak shape, fit
// range, P_main / P_pois, model variations). All rows = stored particle
// combinations (isBest >= 0), corrected momenta, cor_ok == 1.
//
//   1. single-cut scans (each cut alone) and N-1 scans around the start point
//      (oa_lab < 50; the former start cut M(pipi) < 0.46 sculpts the fit range
//      and was removed after the review): S, sigma_S, F_stat = S/sigma_S, S/B,
//      CB share, eps_S from pass/fail fits, eps_B; folds 0 and 1 separately
//
// Reading of the numbers (review):
//   - F_stat = S/sigma_S is statistical only, with the peak position and width
//     frozen. It compares selections; it is NOT a significance. Every yield is
//     quoted as S +- stat +- model spread.
//   - CB = 2 sqrt(N++ N--) without a k factor. The like-sign samples are very
//     asymmetric (N--/N++ about 5) and almost empty at small opening angle, so
//     CB shares are estimates and lower limits at small opening angle.
//   - The classes of section 6 (OA, M(ee), shared ring) are signs of what the
//     peak holds; they cannot be turned into decay-mode fractions without a
//     simulation.
//   2. guard against sculpting: pass fraction against M4 on like-sign,
//      unlike-sign and mixed rows; fake peak content from the mixed rows
//   3. two-dimensional maps of F_stat and S
//   4. bias control: working points defined on one half of the events,
//      evaluated on the other half (fold 0 / 1 and random halves)
//   5. candidates loose / nominal / pure with model spread, Poisson fit,
//      like-sign null test, displaced-peak scan; none / slides / old code
//   6. content of the peak by class (width free per class)
//   7. cor_ok fallback variant
//
// Output: <prod>/scan/peak_*.txt, <prod>/scan/peak_wp.root,
//         plots/output/eta_scan/peak_*.png (.pdf for the main ones)
//
// Run (from the repository root):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/aclic_peak", kTRUE)' 'plots/scan_peak.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "ScanEngine.h"

namespace spk {

// ============================================================================
// printing
// ============================================================================
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
void OpenTable(const char* name) {
    if (gOut) fclose(gOut);
    gOut = fopen(Form("%s/%s", scan::kScanDir, name), "w");
    P("# %s  --  written by plots/scan_peak.C\n", name);
    P("# rows = stored particle combinations (isBest >= 0), corrected momenta, cor_ok == 1 unless stated\n");
    P("# peak shape frozen: position %.4f, width %.1f MeV; P_main = chi2 fit of N - CB, pol%d, range %.4f-%.4f\n",
      scan::kPeakMean, 1e3 * scan::kPeakSigma, scan::kFitOrder, scan::kFitLo, scan::kFitHi);
    P("# window %.5f-%.5f; S = whole Gaussian; F = S/sigma_S; S/B = S in the window / (polynomial + CB)\n",
      scan::kPeakLo, scan::kPeakHi);
    P("# F (F_stat) is statistical only, with a frozen peak shape: it compares selections and is NOT a significance (model spread of S: 19-26 %%).\n");
    P("# CB = 2 sqrt(N++ N--), no k factor; N--/N++ is about 5, like-sign pairs at small opening angle are almost absent:\n");
    P("# CB and CB/N are estimates, and lower limits at small opening angle. S, sigma_S and S/B do not depend on the CB (Poisson fit of N agrees).\n\n");
}
void CloseTable() { if (gOut) fclose(gOut); gOut = nullptr; }

// ============================================================================
// cut grid
// ============================================================================
const double kNo = 1e9;
enum { kOaLo = 0, kOaHi, kOaLab, kMpLo, kMpHi, kMm4, kNV };
const int kMaxIdx = 16;

struct VarDef {
    std::string key, name, title;
    std::vector<double> val;        // value used in the comparison (+-kNo = no cut)
    std::vector<std::string> lab;   // label of the value
    int loose = 1;                  // +1: "no cut" is the last index, -1: the first
    int start = 0;                  // index at the start point
    int none = 0;                   // index of "no cut"
    int n() const { return (int)val.size(); }
};
std::vector<VarDef> gVar;
std::vector<std::vector<int>> gFlag;   // [var][idx]: 1 = sculpting flag, not used

void BuildGrid() {
    gVar.assign(kNV, VarDef());
    {   // OA lower: 0 (= no cut) ... 4 deg, step 0.5
        VarDef& v = gVar[kOaLo]; v.key = "oa_lo"; v.name = "OA >"; v.title = "lower cut on OA(e^{+}e^{-}) [deg]";
        for (int i = 0; i <= 8; ++i) { v.val.push_back(i == 0 ? -kNo : 0.5 * i); v.lab.push_back(i == 0 ? "0" : Form("%.1f", 0.5 * i)); }
        v.loose = -1; v.none = 0; v.start = 0;
    }
    {   // OA upper: 4, 6, 8, 10, 15, 20, 30, none
        VarDef& v = gVar[kOaHi]; v.key = "oa_hi"; v.name = "OA <"; v.title = "upper cut on OA(e^{+}e^{-}) [deg]";
        for (double x : {4., 6., 8., 10., 15., 20., 30.}) { v.val.push_back(x); v.lab.push_back(Form("%.0f", x)); }
        v.val.push_back(kNo); v.lab.push_back("none");
        v.loose = 1; v.none = v.n() - 1; v.start = v.none;
    }
    {   // oa_lab: 25 ... 70 step 5, none
        VarDef& v = gVar[kOaLab]; v.key = "oa_lab"; v.name = "oa_lab <"; v.title = "upper cut on oa_lab: lab angle (#pi^{+}#pi^{-}),(e^{+}e^{-}) [deg]";
        for (int i = 0; i < 10; ++i) { v.val.push_back(25. + 5. * i); v.lab.push_back(Form("%.0f", 25. + 5. * i)); }
        v.val.push_back(kNo); v.lab.push_back("none");
        v.loose = 1; v.none = v.n() - 1; v.start = 5;
    }
    {   // M(pipi) lower: none, 0.28 ... 0.38 step 0.02
        VarDef& v = gVar[kMpLo]; v.key = "mpipi_lo"; v.name = "M(pipi) >"; v.title = "lower cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]";
        v.val.push_back(-kNo); v.lab.push_back("none");
        for (int i = 0; i < 6; ++i) { v.val.push_back(0.28 + 0.02 * i); v.lab.push_back(Form("%.2f", 0.28 + 0.02 * i)); }
        v.loose = -1; v.none = 0; v.start = 0;
    }
    {   // M(pipi) upper: 0.40 ... 0.50 step 0.01, none
        VarDef& v = gVar[kMpHi]; v.key = "mpipi_hi"; v.name = "M(pipi) <"; v.title = "upper cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]";
        for (int i = 0; i <= 10; ++i) { v.val.push_back(0.40 + 0.01 * i); v.lab.push_back(Form("%.2f", 0.40 + 0.01 * i)); }
        v.val.push_back(kNo); v.lab.push_back("none");
        // start: open. The former start value 0.46 (index 6) is flagged by the sculpting guard
        // (S_fake = -203 +- 11 rows) and is kept only as the labelled comparison point kOldStartMpHi.
        v.loose = 1; v.none = v.n() - 1; v.start = v.none;
    }
    {   // missing mass of the four particles, upper: 2.60 ... 2.95 step 0.05, none
        VarDef& v = gVar[kMm4]; v.key = "mm4_hi"; v.name = "MM4 <"; v.title = "upper cut on MM4: missing mass of #pi#pi ee [GeV/c^{2}]";
        for (int i = 0; i < 8; ++i) { v.val.push_back(2.60 + 0.05 * i); v.lab.push_back(Form("%.2f", 2.60 + 0.05 * i)); }
        v.val.push_back(kNo); v.lab.push_back("none");
        v.loose = 1; v.none = v.n() - 1; v.start = v.none;
    }
    gFlag.assign(kNV, std::vector<int>());
    for (int v = 0; v < kNV; ++v) gFlag[v].assign(gVar[v].n(), 0);
}

// a working point: one grid index per variable
struct Pt {
    int i[kNV];
    bool operator==(const Pt& o) const { for (int v = 0; v < kNV; ++v) if (i[v] != o.i[v]) return false; return true; }
    bool operator<(const Pt& o) const { for (int v = 0; v < kNV; ++v) if (i[v] != o.i[v]) return i[v] < o.i[v]; return false; }
};
Pt PtNone() { Pt p; for (int v = 0; v < kNV; ++v) p.i[v] = gVar[v].none; return p; }
Pt PtStart() { Pt p; for (int v = 0; v < kNV; ++v) p.i[v] = gVar[v].start; return p; }
const int kOldStartMpHi = 6;            // M(pipi) < 0.46: start cut of the first version, flagged
Pt PtOldStart() { Pt p = PtStart(); p.i[kMpHi] = kOldStartMpHi; return p; }
Pt With(Pt p, int v, int idx) { p.i[v] = idx; return p; }
std::string PtText(const Pt& p) {
    std::string s;
    for (int v = 0; v < kNV; ++v) {
        if (p.i[v] == gVar[v].none) continue;
        if (!s.empty()) s += ", ";
        s += gVar[v].name + " " + gVar[v].lab[p.i[v]];
    }
    return s.empty() ? "no cut" : s;
}
bool PtFlagged(const Pt& p) { for (int v = 0; v < kNV; ++v) if (gFlag[v][p.i[v]]) return true; return false; }

struct Cut { float x[kNV]; };
Cut ToCut(const Pt& p) { Cut c; for (int v = 0; v < kNV; ++v) c.x[v] = (float)gVar[v].val[p.i[v]]; return c; }
inline bool PassVar(const scan::Row& r, int v, float x) {
    switch (v) {
        case kOaLo:  return r.oa > x;
        case kOaHi:  return r.oa < x;
        case kOaLab: return r.oa_lab < x;
        case kMpLo:  return r.mpipi > x;
        case kMpHi:  return r.mpipi < x;
        default:     return r.mm4 < x;
    }
}
inline bool PassCut(const scan::Row& r, const Cut& c) {
    return r.oa > c.x[kOaLo] && r.oa < c.x[kOaHi] && r.oa_lab < c.x[kOaLab] &&
           r.mpipi > c.x[kMpLo] && r.mpipi < c.x[kMpHi] && r.mm4 < c.x[kMm4];
}

// ---- sample: all rows, fold (parity of the event group), or a random half ----
struct Samp { int mode = 0; unsigned seed = 0; int half = 0; };
inline bool InSamp(const scan::Row& r, const Samp& s) {
    if (s.mode == 0) return true;
    if (s.mode == 1) return r.fold == s.half;
    unsigned x = (unsigned)r.grp * 2654435761u + s.seed * 0x9E3779B9u;
    x ^= x >> 16; x *= 0x85ebca6bu; x ^= x >> 13; x *= 0xc2b2ae35u; x ^= x >> 16;
    return (int)(x & 1u) == s.half;
}
Samp SampAll() { return Samp(); }
Samp SampFold(int f) { Samp s; s.mode = 1; s.half = f; return s; }
Samp SampRnd(unsigned seed, int h) { Samp s; s.mode = 2; s.seed = seed; s.half = h; return s; }
Samp Other(Samp s) { s.half = 1 - s.half; return s; }

// selection of a working point; failVar >= 0: rows passing the point with this
// variable open, but FAILING the cut value failVal on it
scan::Sel MkSel(const Pt& p, const Samp& s = Samp(), int failVar = -1, int failIdx = 0) {
    Cut c = ToCut(failVar >= 0 ? With(p, failVar, gVar[failVar].none) : p);
    const float fv = failVar >= 0 ? (float)gVar[failVar].val[failIdx] : 0.f;
    return [c, s, failVar, fv](const scan::Row& r) {
        if (!InSamp(r, s)) return false;
        if (!PassCut(r, c)) return false;
        if (failVar >= 0 && PassVar(r, failVar, fv)) return false;
        return true;
    };
}

// ============================================================================
// evaluation of one spectrum
// ============================================================================
struct Ev {
    bool ok = false;
    double N = 0, CB = 0, B = 0, S = 0, eS = 0, F = 0, SoB = 0, cbs = 0, chi2 = 0, prob = 0, Swin = 0;
    double npp = 0, nmm = 0;            // like-sign rows (e+e+, e-e-) in the peak window
    int ndf = 0;
};
Ev Eval(const scan::Spectra& sp, const scan::PeakCfg& c = scan::PeakCfg()) {
    const scan::PeakFit f = scan::FitPeakMain(sp, c);
    Ev e; e.ok = f.ok;
    e.N = f.N; e.CB = f.CB; e.B = f.B; e.S = f.S; e.eS = f.eS; e.F = f.F; e.SoB = f.SoB;
    e.cbs = f.N > 0 ? f.CB / f.N : 0; e.chi2 = f.chi2; e.ndf = f.ndf; e.prob = f.prob; e.Swin = f.Swin;
    const scan::WinSum w = scan::PeakWindow(sp, c);
    e.npp = w.n[1]; e.nmm = w.n[2];
    return e;
}

// one point of a one-dimensional scan
struct ScanPt {
    Ev pass, fail;
    double epsS = 1, eEpsS = 0, epsB = 1, sqB = 1, z = 0, clos = 0;
};
void PassFail(ScanPt& q, const Ev& all) {
    const double sp = q.pass.S, sf = q.fail.S, tot = sp + sf;
    if (std::fabs(tot) > 1e-9) {
        q.epsS = sp / tot;
        q.eEpsS = std::sqrt(sf * sf * q.pass.eS * q.pass.eS + sp * sp * q.fail.eS * q.fail.eS) / (tot * tot);
    }
    const double bp = q.pass.B + q.pass.CB, bf = q.fail.B + q.fail.CB;
    q.epsB = (bp + bf) > 0 ? bp / (bp + bf) : 1;
    q.sqB = std::sqrt(std::max(q.epsB, 0.));
    q.z = q.eEpsS > 0 ? (q.epsS - q.sqB) / q.eEpsS : 0;
    q.clos = tot - all.S;
}

// N-1 scans of all six variables around the point 'base' for several samples
// in one pass: out[sample][var][idx]
typedef std::vector<std::vector<ScanPt>> Scan6;
std::vector<Scan6> Scan1D(const scan::Data& d, const Pt& base, const std::vector<Samp>& samp, bool withFail) {
    std::vector<scan::Sel> sels;
    for (const Samp& s : samp)
        for (int v = 0; v < kNV; ++v)
            for (int i = 0; i < gVar[v].n(); ++i) {
                sels.push_back(MkSel(With(base, v, i), s));
                if (withFail) sels.push_back(MkSel(base, s, v, i));
            }
    const std::vector<scan::Spectra> sp = scan::FillMany(d, sels);
    std::vector<Scan6> out(samp.size(), Scan6(kNV));
    size_t k = 0;
    for (size_t is = 0; is < samp.size(); ++is)
        for (int v = 0; v < kNV; ++v) {
            out[is][v].resize(gVar[v].n());
            for (int i = 0; i < gVar[v].n(); ++i) {
                out[is][v][i].pass = Eval(sp[k++]);
                if (withFail) out[is][v][i].fail = Eval(sp[k++]);
            }
            if (withFail)
                for (int i = 0; i < gVar[v].n(); ++i) PassFail(out[is][v][i], out[is][v][gVar[v].none].pass);
        }
    return out;
}

// ---- plateau of a one-dimensional scan --------------------------------------
// best = allowed value with the largest F; plateau = contiguous allowed values
// around it with F >= F_best - 1, widened to two grid steps where possible;
// centre = middle of the plateau (a half-integer middle goes to the looser value)
struct Plat { int lo = 0, hi = 0, ctr = 0, best = 0; double fmax = 0; };
Plat Plateau1D(const std::vector<double>& F, const std::vector<int>& flag, int looseDir) {
    const int n = (int)F.size();
    Plat p; p.best = -1;
    for (int i = 0; i < n; ++i) if (!flag[i] && (p.best < 0 || F[i] > F[p.best])) p.best = i;
    if (p.best < 0) { p.best = 0; return p; }
    p.fmax = F[p.best]; p.lo = p.hi = p.best;
    while (p.lo > 0 && !flag[p.lo - 1] && F[p.lo - 1] >= p.fmax - 1.0) p.lo--;
    while (p.hi < n - 1 && !flag[p.hi + 1] && F[p.hi + 1] >= p.fmax - 1.0) p.hi++;
    while (p.hi - p.lo < 2) {
        const bool cl = p.lo > 0 && !flag[p.lo - 1], ch = p.hi < n - 1 && !flag[p.hi + 1];
        if (!cl && !ch) break;
        if (cl && (!ch || F[p.lo - 1] >= F[p.hi + 1])) p.lo--; else p.hi++;
    }
    const int s = p.lo + p.hi;
    p.ctr = (s % 2 == 0) ? s / 2 : (looseDir > 0 ? (s + 1) / 2 : (s - 1) / 2);
    return p;
}


// ============================================================================
// guard against sculpting
// ============================================================================
// Pass fraction of every cut value against M4 in the fine bins of the engine,
// for four row sets: 0 = unlike-sign, 1 = like-sign (++ and -- added),
// 2, 3 = mixed rows (two independent halves of the mixing offsets).
// Context 0: the cut alone; context 1: the other variables at the start point.
const int kNSrc = 4;
struct GuardHist {
    std::vector<double> pas[kNSrc][2][kNV][kMaxIdx];
    std::vector<double> tot[kNSrc][2][kNV];
    GuardHist() {
        for (int s = 0; s < kNSrc; ++s)
            for (int c = 0; c < 2; ++c)
                for (int v = 0; v < kNV; ++v) {
                    tot[s][c][v].assign(scan::kNbin, 0.);
                    for (int i = 0; i < kMaxIdx; ++i) pas[s][c][v][i].assign(scan::kNbin, 0.);
                }
    }
};
inline void GuardFill(GuardHist& g, int src, const scan::Row& r, double m4) {
    if (!(m4 >= scan::kMlo && m4 < scan::kMhi)) return;
    const int ib = (int)((m4 - scan::kMlo) / scan::kBinW);
    if (ib < 0 || ib >= scan::kNbin) return;
    bool ps[kNV];
    for (int v = 0; v < kNV; ++v) ps[v] = PassVar(r, v, (float)gVar[v].val[gVar[v].start]);
    for (int v = 0; v < kNV; ++v) {
        bool oth = true;
        for (int w = 0; w < kNV; ++w) if (w != v && !ps[w]) oth = false;
        g.tot[src][0][v][ib] += 1;
        if (oth) g.tot[src][1][v][ib] += 1;
        for (int i = 0; i < gVar[v].n(); ++i)
            if (PassVar(r, v, (float)gVar[v].val[i])) {
                g.pas[src][0][v][i][ib] += 1;
                if (oth) g.pas[src][1][v][i][ib] += 1;
            }
    }
}
std::vector<double> Add(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> c(a.size());
    for (size_t i = 0; i < a.size(); ++i) c[i] = a[i] + b[i];
    return c;
}

// ---- mixed rows ---------------------------------------------------------------
// (pi+pi-) of one unlike-sign row with (e+e-) of a row of another event group,
// 20 fixed offsets in the row order. M4, the two pair angles and the missing
// mass are recomputed from the corrected momenta; M(pipi), OA and M(ee) stay.
// A mixed row has no eta and no correlation between the pairs, but the
// kinematic limits (M4 > M(pipi) + M(ee), MM4 < sqrt(s) - M4) hold.
const int kNOff = 20;
const int kOff[kNOff] = {1009, 2003, 3001, 4001, 5003, 6007, 7001, 8009, 9001, 10007,
                         11003, 12007, 13001, 14009, 15013, 16001, 17011, 18013, 19001, 20011};
struct MixPool {
    const std::vector<scan::Row>* rows = nullptr;
    std::vector<int> idx;           // unlike-sign rows with cor_ok == 1
    std::vector<float> pp, ee;      // px, py, pz, E [MeV] of the pion pair and of the lepton pair
    double pb = 0, eIni = 0;        // beam momentum, total energy
};
struct MixKin { double m4 = 0, oal = 0, mm = 0, oae = 0; };
void BuildPool(const scan::Data& d, MixPool& mp) {
    const double mPi = 1e3 * scan::kMpi, mE = 1e3 * scan::kMe, mP = 1e3 * scan::kMp, d2r = TMath::DegToRad();
    const double eb = scan::kTbeam + mP;
    mp.pb = std::sqrt(eb * eb - mP * mP); mp.eIni = eb + mP;
    mp.rows = &d.s[0];
    const std::vector<scan::Row>& v = d.s[0];
    for (size_t i = 0; i < v.size(); ++i) if (v[i].cor_ok) mp.idx.push_back((int)i);
    const size_t n = mp.idx.size();
    mp.pp.resize(4 * n); mp.ee.resize(4 * n);
    auto add = [&](float* o, double p, double th, double ph, double m) {
        const double st = std::sin(th * d2r);
        o[0] += (float)(p * st * std::cos(ph * d2r)); o[1] += (float)(p * st * std::sin(ph * d2r));
        o[2] += (float)(p * std::cos(th * d2r)); o[3] += (float)std::sqrt(p * p + m * m);
    };
    for (size_t k = 0; k < n; ++k) {
        const scan::Row& r = v[mp.idx[k]];
        float* a = &mp.pp[4 * k]; float* b = &mp.ee[4 * k];
        a[0] = a[1] = a[2] = a[3] = 0; b[0] = b[1] = b[2] = b[3] = 0;
        add(a, r.pip_p, r.pip_theta, r.pip_phi, mPi); add(a, r.pim_p, r.pim_theta, r.pim_phi, mPi);
        add(b, r.ep_p, r.ep_theta, r.ep_phi, mE); add(b, r.em_p, r.em_theta, r.em_phi, mE);
    }
}
inline void Kin(const MixPool& mp, size_t k, size_t j, MixKin& q, bool withEta) {
    const float* a = &mp.pp[4 * k]; const float* b = &mp.ee[4 * j];
    const double x = (double)a[0] + b[0], y = (double)a[1] + b[1], z = (double)a[2] + b[2], e = (double)a[3] + b[3];
    const double p2 = x * x + y * y + z * z, m2 = e * e - p2;
    q.m4 = 1e-3 * (m2 > 0 ? std::sqrt(m2) : -std::sqrt(-m2));
    const double pa = std::sqrt((double)a[0] * a[0] + (double)a[1] * a[1] + (double)a[2] * a[2]);
    const double pc = std::sqrt((double)b[0] * b[0] + (double)b[1] * b[1] + (double)b[2] * b[2]);
    double cs = ((double)a[0] * b[0] + (double)a[1] * b[1] + (double)a[2] * b[2]) / (pa * pc);
    q.oal = std::acos(std::max(-1., std::min(1., cs))) * TMath::RadToDeg();
    const double q2 = (mp.eIni - e) * (mp.eIni - e) - x * x - y * y - (mp.pb - z) * (mp.pb - z);
    q.mm = 1e-3 * (q2 > 0 ? std::sqrt(q2) : -std::sqrt(-q2));
    if (withEta) {     // angle of the two pairs in the frame of (p_total, m_eta), as in the stored oa_eta
        const double mEta = 1e3 * scan::kMeta, eh = std::sqrt(p2 + mEta * mEta), g = eh / mEta;
        const double bx = x / eh, by = y / eh, bz = z / eh, b2 = bx * bx + by * by + bz * bz;
        double u[2][3];
        const float* w[2] = {a, b};
        for (int t = 0; t < 2; ++t) {
            const double bp = bx * w[t][0] + by * w[t][1] + bz * w[t][2];
            const double f = (b2 > 0 ? (g - 1) * bp / b2 : 0.) - g * w[t][3];
            u[t][0] = w[t][0] + f * bx; u[t][1] = w[t][1] + f * by; u[t][2] = w[t][2] + f * bz;
        }
        const double n0 = std::sqrt(u[0][0] * u[0][0] + u[0][1] * u[0][1] + u[0][2] * u[0][2]);
        const double n1 = std::sqrt(u[1][0] * u[1][0] + u[1][1] * u[1][1] + u[1][2] * u[1][2]);
        const double c2 = (u[0][0] * u[1][0] + u[0][1] * u[1][1] + u[0][2] * u[1][2]) / (n0 * n1);
        q.oae = std::acos(std::max(-1., std::min(1., c2))) * TMath::RadToDeg();
    }
}
// closure of the recomputation on unmixed rows, then the pass counts of the cut grid
struct MixCheck { double dM4 = 0, dOaLab = 0, dMm4 = 0, dOaEta = 0; long nmix = 0; };
MixCheck FillMixedGrid(const MixPool& mp, GuardHist& g) {
    const std::vector<scan::Row>& v = *mp.rows;
    const size_t n = mp.idx.size();
    MixCheck mc; MixKin q;
    for (size_t k = 0; k < n; k += 7) {
        const scan::Row& r = v[mp.idx[k]];
        Kin(mp, k, k, q, true);
        mc.dM4 = std::max(mc.dM4, std::fabs(q.m4 - r.m4));
        mc.dOaLab = std::max(mc.dOaLab, std::fabs(q.oal - r.oa_lab));
        mc.dOaEta = std::max(mc.dOaEta, std::fabs(q.oae - r.oa_eta));
        if (r.mm4 > 0.5) mc.dMm4 = std::max(mc.dMm4, std::fabs(q.mm - r.mm4));
    }
    for (size_t k = 0; k < n; ++k) {
        const scan::Row& a = v[mp.idx[k]];
        for (int o = 0; o < kNOff; ++o) {
            const size_t j = (k + kOff[o]) % n;
            const scan::Row& b = v[mp.idx[j]];
            if (a.grp == b.grp) continue;
            Kin(mp, k, j, q, false);
            if (!(q.m4 >= scan::kMlo && q.m4 < scan::kMhi)) continue;
            scan::Row m = a;
            m.oa = b.oa; m.mee = b.mee; m.m4 = (float)q.m4; m.oa_lab = (float)q.oal; m.mm4 = (float)q.mm;
            GuardFill(g, 2 + (o & 1), m, q.m4);
            mc.nmix++;
        }
    }
    return mc;
}
// pass counts of arbitrary selections on the mixed rows, fit range only
// (pas[half][selection][fine bin], tot[half][fine bin])
struct MixSpec { std::vector<std::vector<double>> pas[2]; std::vector<double> tot[2]; };
MixSpec MixedFill(const MixPool& mp, const std::vector<scan::Sel>& sels) {
    MixSpec ms;
    for (int h = 0; h < 2; ++h) {
        ms.tot[h].assign(scan::kNbin, 0.);
        ms.pas[h].assign(sels.size(), std::vector<double>(scan::kNbin, 0.));
    }
    const std::vector<scan::Row>& v = *mp.rows;
    const size_t n = mp.idx.size();
    MixKin q;
    for (size_t k = 0; k < n; ++k) {
        const scan::Row& a = v[mp.idx[k]];
        for (int o = 0; o < kNOff; ++o) {
            const size_t j = (k + kOff[o]) % n;
            const scan::Row& b = v[mp.idx[j]];
            if (a.grp == b.grp) continue;
            Kin(mp, k, j, q, false);
            if (!(q.m4 >= scan::kFitLo && q.m4 < scan::kFitHi)) continue;
            Kin(mp, k, j, q, true);
            const int ib = (int)((q.m4 - scan::kMlo) / scan::kBinW);
            scan::Row m = a;
            m.oa = b.oa; m.mee = b.mee; m.m4 = (float)q.m4; m.oa_lab = (float)q.oal; m.mm4 = (float)q.mm;
            m.oa_eta = (float)q.oae;
            const int h = o & 1;
            ms.tot[h][ib] += 1;
            for (size_t s = 0; s < sels.size(); ++s) if (sels[s](m)) ms.pas[h][s][ib] += 1;
        }
    }
    return ms;
}

// ---- turn-on interval of a pass fraction ---------------------------------------
// f in bins of 'reb' fine bins between 0.40 and 0.80; fmin, fmax = extremes;
// mA, mB = lowest and highest mass where f crosses fmin + 0.1 (fmax - fmin) or
// fmin + 0.9 (fmax - fmin); dfit = change of f inside the fit range.
// edgeIn: mA or mB inside the fit range and dfit >= 0.05.
struct Turn { bool ok = false; double fmin = 0, fmax = 0, delta = 0, mA = -1, mB = -1, fLo = 0, fHi = 0, dfit = 0; bool edgeIn = false; };
Turn TurnOn(const std::vector<double>& pas, const std::vector<double>& tot, int reb, double minRows) {
    Turn t;
    std::vector<double> x, f;
    for (int i = 0; i + reb <= scan::kNbin; i += reb) {
        const double xc = scan::kMlo + (i + 0.5 * reb) * scan::kBinW;
        if (xc < 0.40 || xc > 0.80) continue;
        double a = 0, b = 0;
        for (int j = 0; j < reb; ++j) { a += pas[i + j]; b += tot[i + j]; }
        if (b < minRows) continue;
        x.push_back(xc); f.push_back(a / b);
    }
    const int n = (int)x.size();
    if (n < 5) return t;
    t.ok = true;
    t.fmin = *std::min_element(f.begin(), f.end()); t.fmax = *std::max_element(f.begin(), f.end());
    t.delta = t.fmax - t.fmin;
    double fl = 1e9, fh = -1e9; bool first = true;
    for (int i = 0; i < n; ++i) {
        if (x[i] < scan::kFitLo || x[i] > scan::kFitHi) continue;
        if (first) { t.fLo = f[i]; first = false; }
        t.fHi = f[i];
        fl = std::min(fl, f[i]); fh = std::max(fh, f[i]);
    }
    t.dfit = fh > fl ? fh - fl : 0;
    if (t.delta < 1e-6) return t;
    double mA = 1e9, mB = -1e9;
    for (int i = 0; i + 1 < n; ++i)
        for (double L : {0.1, 0.9}) {
            const double g0 = (f[i] - t.fmin) / t.delta - L, g1 = (f[i + 1] - t.fmin) / t.delta - L;
            if (g0 * g1 > 0 || g0 == g1) continue;
            const double m = x[i] + (x[i + 1] - x[i]) * g0 / (g0 - g1);
            mA = std::min(mA, m); mB = std::max(mB, m);
        }
    if (mA < 1e8) { t.mA = mA; t.mB = mB; }
    const bool aIn = t.mA >= scan::kFitLo && t.mA <= scan::kFitHi, bIn = t.mB >= scan::kFitLo && t.mB <= scan::kFitHi;
    t.edgeIn = (aIn || bIn) && t.dfit >= 0.05;
    return t;
}

// ---- fake peak content from a pass fraction -----------------------------------
// Gaussian + pol2 fitted to an arbitrary spectrum y in the fit range with the
// variances w (linear least squares); returns the Gaussian content. coef =
// derivative of the result with respect to y_i (the estimator is linear).
double GaussCoef(const std::vector<double>& y, const std::vector<double>& w, std::vector<double>* coef = nullptr) {
    std::vector<double> F, yy, var;
    const double x0 = 0.5 * (scan::kFitLo + scan::kFitHi), dx = 0.5 * (scan::kFitHi - scan::kFitLo);
    int npts = 0;
    for (int i = 0; i < scan::kNbin; ++i) {
        const double x1 = scan::kMlo + i * scan::kBinW, xc = x1 + 0.5 * scan::kBinW;
        if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
        const double t = (xc - x0) / dx;
        F.push_back(scan::GaussFrac(x1, x1 + scan::kBinW, scan::kPeakMean, scan::kPeakSigma));
        F.push_back(1.); F.push_back(t); F.push_back(t * t);
        yy.push_back(y[i]); var.push_back(std::max(w[i], 1e-9));
        npts++;
    }
    const scan::LinRes r = scan::LinFit(npts, 4, yy.data(), var.data(), F.data());
    if (!r.ok) return 0;
    if (coef) {
        coef->assign(scan::kNbin, 0.);
        int q = 0;
        for (int i = 0; i < scan::kNbin; ++i) {
            const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
            if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
            double a = 0;
            for (int l = 0; l < 4; ++l) a += r.cov[0][l] * F[q * 4 + l] / var[q];
            (*coef)[i] = a; q++;
        }
    }
    return r.par[0];
}
// y = Bref(m) * f(m): Bref = smooth background of the base selection (the pol2
// of the Poisson fit of the unsubtracted spectrum, so the CB is inside), f =
// pass fraction of the mixed rows; fitted with Gaussian + pol2 and the weights
// of a data fit. A pol2 background stays a pol2 only if f is a constant.
double FakeOne(const scan::PeakFit& smooth, const std::vector<double>& pas, const std::vector<double>& tot, double* ebin = nullptr) {
    std::vector<double> y(scan::kNbin, 0.), w(scan::kNbin, 1.), ef(scan::kNbin, 0.), coef;
    for (int i = 0; i < scan::kNbin; ++i) {
        const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
        if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
        const double b = smooth.Bkg(xc);
        const double f = tot[i] > 0 ? pas[i] / tot[i] : 0;
        y[i] = b * f; w[i] = std::max(y[i], 1.);
        ef[i] = tot[i] > 0 ? b * std::sqrt(std::max(f * (1 - f), 0.) / tot[i]) : 0;
    }
    const double s = GaussCoef(y, w, &coef);
    if (ebin) {
        double e2 = 0;
        for (int i = 0; i < scan::kNbin; ++i) e2 += coef[i] * coef[i] * ef[i] * ef[i];
        *ebin = std::sqrt(e2);
    }
    return s;
}
// both mixed halves together; error = the larger of the binomial error and of
// half the difference of the two halves
struct Fake { double s = 0, e = 0; };
Fake FakeS(const scan::PeakFit& smooth, const std::vector<double>& p0, const std::vector<double>& t0,
           const std::vector<double>& p1, const std::vector<double>& t1) {
    Fake f; double eb = 0;
    f.s = FakeOne(smooth, Add(p0, p1), Add(t0, t1), &eb);
    const double a = FakeOne(smooth, p0, t0), b = FakeOne(smooth, p1, t1);
    f.e = std::max(eb, 0.5 * std::fabs(a - b));
    return f;
}

// ============================================================================
// definition of the candidates on a sample
// ============================================================================
scan::Data Split(const scan::Data& src, const Samp& s) {
    scan::Data o;
    for (int k = 0; k < 3; ++k)
        for (const auto& r : src.s[k]) if (InSamp(r, s)) o.s[k].push_back(r);
    return o;
}
// start point with the flagged values opened
Pt PtBase0() { Pt p = PtStart(); for (int v = 0; v < kNV; ++v) if (gFlag[v][p.i[v]]) p.i[v] = gVar[v].none; return p; }

// grid of all combinations of allowed (not flagged) values of a fixed sub-list
std::vector<Pt> MediumGrid() {
    const std::vector<std::vector<int>> want = {
        {0, 1, 2, 3, 4},                 // OA lower: 0 ... 2.0
        {3, 4, 5, 6, 7},                 // OA upper: 10, 15, 20, 30, none
        {3, 4, 5, 6, 7, 8, 9, 10},       // oa_lab: 40 ... 70, none
        {0, 1, 2, 3, 4, 5},              // M(pipi) lower: none, 0.28 ... 0.36
        {0, 2, 4, 6, 8, 11},             // M(pipi) upper: 0.40, 0.42, 0.44, 0.46, 0.48, none
        {0, 2, 4, 5, 8}};                // MM4 upper: 2.60, 2.70, 2.80, 2.85, none
    std::vector<std::vector<int>> use(kNV);
    for (int v = 0; v < kNV; ++v) {
        for (int i : want[v]) if (!gFlag[v][i]) use[v].push_back(i);
        if (use[v].empty()) use[v].push_back(gVar[v].none);
    }
    std::vector<Pt> g;
    Pt p;
    for (int a : use[0]) for (int b : use[1]) for (int c : use[2]) for (int e : use[3]) for (int f : use[4]) for (int h : use[5]) {
        p.i[0] = a; p.i[1] = b; p.i[2] = c; p.i[3] = e; p.i[4] = f; p.i[5] = h;
        g.push_back(p);
    }
    return g;
}
std::vector<Ev> EvalPts(const scan::Data& d, const std::vector<Pt>& pts, const Samp& s = Samp()) {
    std::vector<scan::Sel> sels;
    for (const Pt& p : pts) sels.push_back(MkSel(p, s));
    const std::vector<scan::Spectra> sp = scan::FillMany(d, sels);
    std::vector<Ev> ev;
    for (const auto& x : sp) ev.push_back(Eval(x));
    return ev;
}

// nominal = self-consistent plateau centre: N-1 scans around the current
// point; each variable moves to the centre of its plateau (allowed values
// with F >= F_best - 1), or stays open when "no cut" has the largest F;
// repeated until the point no longer moves.
// loose  = grid point with the most rows in the window among F >= F_best - 1
// pure   = grid point with the largest S/B among F >= F_best - 1
// best   = grid point with the largest F (the biased choice, for comparison)
struct Defn {
    Pt nominal, loose, pure, best;
    Scan6 nm1;
    std::vector<Plat> plat;
    std::vector<Pt> grid; std::vector<Ev> gev;
    Ev evNominal;
    double fbest = 0; int nPlat = 0, nIter = 0; bool converged = false;
    std::vector<Pt> path;
};
Defn Define(const scan::Data& dd, bool verbose) {
    Defn D;
    Pt p = PtBase0();
    D.path.push_back(p);
    for (int it = 0; it < 6; ++it) {
        D.nm1 = Scan1D(dd, p, {Samp()}, false)[0];
        D.plat.assign(kNV, Plat());
        Pt q = p;
        for (int v = 0; v < kNV; ++v) {
            std::vector<double> F;
            for (const auto& x : D.nm1[v]) F.push_back(x.pass.F);
            D.plat[v] = Plateau1D(F, gFlag[v], gVar[v].loose);
            q.i[v] = D.plat[v].best == gVar[v].none ? gVar[v].none : D.plat[v].ctr;
        }
        D.nIter = it + 1;
        if (verbose) printf("   iteration %d: %s  ->  %s\n", it + 1, PtText(p).c_str(), PtText(q).c_str());
        if (q == p) { D.converged = true; break; }
        bool seen = false;
        for (const Pt& o : D.path) if (o == q) seen = true;
        p = q; D.path.push_back(p);
        if (seen) {                       // a cycle: keep this point, scan once more around it
            D.nm1 = Scan1D(dd, p, {Samp()}, false)[0];
            for (int v = 0; v < kNV; ++v) {
                std::vector<double> F;
                for (const auto& x : D.nm1[v]) F.push_back(x.pass.F);
                D.plat[v] = Plateau1D(F, gFlag[v], gVar[v].loose);
            }
            break;
        }
    }
    D.nominal = p;
    D.evNominal = D.nm1[0][p.i[0]].pass;
    D.grid = MediumGrid();
    D.gev = EvalPts(dd, D.grid);
    int ib = 0;
    for (size_t k = 0; k < D.grid.size(); ++k) if (D.gev[k].F > D.gev[ib].F) ib = (int)k;
    D.best = D.grid[ib]; D.fbest = std::max(D.gev[ib].F, D.evNominal.F);
    int il = -1, ip = -1;
    for (size_t k = 0; k < D.grid.size(); ++k) {
        if (D.gev[k].F < D.fbest - 1.0) continue;
        D.nPlat++;
        if (il < 0 || D.gev[k].N > D.gev[il].N) il = (int)k;
        if (ip < 0 || D.gev[k].SoB > D.gev[ip].SoB) ip = (int)k;
    }
    if (il < 0) il = ib;
    if (ip < 0) ip = ib;
    D.loose = D.grid[il]; D.pure = D.grid[ip];
    return D;
}

// ============================================================================
// full evaluation of one selection (candidates and reference chains)
// ============================================================================
const int kNDisp = 13;                 // displaced-peak positions 0.575 ... 0.695, step 10 MeV
struct Full {
    std::string name, text;
    scan::Sel sel;
    Pt pt; bool hasPt = false;         // grid point, when the selection is one
    scan::Spectra sp;
    scan::PeakFit pm, pp, lsP, lsM;
    scan::Spread sd;
    double ftot = 0;
    double dPos[kNDisp], dS[kNDisp], dE[kNDisp], dChi2[kNDisp];
    double pullRms = 0, pullMax = 0, pullMean = 0;   // all 13 positions
    double pullRms8 = 0, pullMean8 = 0;              // the 8 positions whose fit range ends below 0.72
    Fake fake;
    Ev fold[2];
    Ev fail;                           // rows of "no cut" failing this selection
    double dF = 0, edF = 0;            // gain in F over "no cut" and its error
    double epsS = 1, eEpsS = 0, epsB = 1;
};
void FullFit(Full& f) {
    f.pm = scan::FitPeakMain(f.sp);
    f.pp = scan::FitPeakPois(f.sp);
    f.lsP = scan::FitPeakLikeSign(f.sp, 1);
    f.lsM = scan::FitPeakLikeSign(f.sp, 2);
    f.sd = scan::PeakModelSpread(f.sp, scan::PeakCfg(), true);
    f.ftot = f.pm.S / std::sqrt(f.pm.eS * f.pm.eS + f.sd.half * f.sd.half);
    double s1 = 0, s2 = 0, t1 = 0, t2 = 0;
    for (int i = 0; i < kNDisp; ++i) {
        const double pos = 0.575 + 0.010 * i, sh = pos - scan::kPeakMean;
        scan::PeakCfg c; c.mean = pos; c.lo = scan::kFitLo + sh; c.hi = scan::kFitHi + sh;
        c.winLo = scan::kPeakLo + sh; c.winHi = scan::kPeakHi + sh;
        const scan::PeakFit q = scan::FitPeakMain(f.sp, c);
        f.dPos[i] = pos; f.dS[i] = q.S; f.dE[i] = q.eS; f.dChi2[i] = q.chi2;
        const double pull = q.eS > 0 ? q.S / q.eS : 0;
        s1 += pull; s2 += pull * pull; f.pullMax = std::max(f.pullMax, std::fabs(pull));
        if (i < 8) { t1 += pull; t2 += pull * pull; }
    }
    f.pullMean = s1 / kNDisp; f.pullRms = std::sqrt(s2 / kNDisp);
    f.pullMean8 = t1 / 8; f.pullRms8 = std::sqrt(t2 / 8);
}

// peak fit with a free width (position frozen): chi2 of P_main scanned in the
// width, minimum by a parabola, error from delta chi2 = 1; the error of S
// includes the width. ok = false when chi2 has no minimum inside 2.5-11 MeV.
struct WFit { bool ok = false; double sigma = 0, eSigma = 0, S = 0, eS = 0, eSfix = 0, chi2 = 0, dchi2 = 0; int ndf = 0; };
WFit WidthFit(const scan::Spectra& sp) {
    WFit w;
    const double s0 = 0.0025, h = 0.0001; const int n = 86;
    std::vector<double> c2(n), S(n), eS(n);
    int im = 0;
    for (int i = 0; i < n; ++i) {
        scan::PeakCfg c; c.sigma = s0 + h * i;
        const scan::PeakFit f = scan::FitPeakMain(sp, c);
        c2[i] = f.chi2; S[i] = f.S; eS[i] = f.eS;
        if (c2[i] < c2[im]) im = i;
    }
    w.ndf = scan::FitPeakMain(sp).ndf - 1;
    if (im == 0 || im == n - 1) return w;
    const double a = (c2[im - 1] - 2 * c2[im] + c2[im + 1]) / (h * h), b = (c2[im + 1] - c2[im - 1]) / (2 * h);
    if (!(a > 0)) return w;
    w.sigma = s0 + h * im - b / a;
    w.chi2 = c2[im] - 0.5 * b * b / a;
    // delta chi2 = 1 on both sides (linear interpolation of the scan)
    double lo = -1, hi = -1;
    for (int i = im; i > 0; --i) if (c2[i - 1] - w.chi2 >= 1.0) { lo = s0 + h * (i - (1.0 - (c2[i] - w.chi2)) / (c2[i - 1] - c2[i])); break; }
    for (int i = im; i < n - 1; ++i) if (c2[i + 1] - w.chi2 >= 1.0) { hi = s0 + h * (i + (1.0 - (c2[i] - w.chi2)) / (c2[i + 1] - c2[i])); break; }
    if (lo < 0 && hi < 0) return w;
    w.eSigma = (lo > 0 && hi > 0) ? 0.5 * (hi - lo) : (lo > 0 ? w.sigma - lo : hi - w.sigma);
    scan::PeakCfg c; c.sigma = w.sigma;
    const scan::PeakFit f = scan::FitPeakMain(sp, c);
    w.S = f.S; w.eSfix = f.eS;
    const double dSds = (S[im + 1] - S[im - 1]) / (2 * h);
    w.eS = std::sqrt(f.eS * f.eS + dSds * dSds * w.eSigma * w.eSigma);
    w.dchi2 = scan::FitPeakMain(sp).chi2 - w.chi2;
    w.ok = true;
    return w;
}

// ============================================================================
// figures
// ============================================================================
const char* kXt = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
const char* kYt = "rows / 2.5 MeV/c^{2}";
const char* kHdr = "p+p 4.53 GeV, all combinations (isBest #geq 0), corrected momenta";
const int kColS = kBlue + 1, kColF0 = kGray + 2, kColF1 = kOrange + 7;

void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.15); gStyle->SetPadTopMargin(0.09);
    gStyle->SetLabelSize(0.05, "xyz"); gStyle->SetTitleSize(0.055, "xyz");
    gStyle->SetTitleOffset(1.3, "y"); gStyle->SetTitleOffset(1.15, "x");
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0); gStyle->SetNumberContours(60);
    gStyle->SetPalette(kBird);
}
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs(Form("%s/%s.png", scan::kFigDir, name));
    if (pdf) c->SaveAs(Form("%s/%s.pdf", scan::kFigDir, name));
}
void Txt(double x, double y, const char* t, double size = 0.045, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align);
    l.DrawLatex(x, y, t);
}
void HLine(double x1, double x2, double y, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid); l->Draw();
}
void VLine(double x, double y1, double y2, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x, y1, x, y2); l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid); l->Draw();
}
void Box(double x1, double y1, double x2, double y2, int col, int fill = 1001, double alpha = 1.) {
    TBox* b = new TBox(x1, y1, x2, y2);
    b->SetFillStyle(fill); b->SetLineWidth(0);
    if (alpha < 1.) b->SetFillColorAlpha(col, alpha); else b->SetFillColor(col);
    b->Draw();
}
TLegend* Leg(double x1, double y1, double x2, double y2, double size, int ncol = 1) {
    TLegend* l = new TLegend(x1, y1, x2, y2);
    l->SetTextSize(size); l->SetNColumns(ncol); l->SetFillStyle(0); l->SetBorderSize(0);
    return l;
}
int gFrame = 0;
// frame with one bin per cut value
TH1D* CatFrame(const VarDef& v, double y0, double y1, const char* yt, double ts = 0.06) {
    TH1D* h = new TH1D(Form("fr%d", gFrame++), "", v.n(), 0, v.n());
    h->SetDirectory(nullptr);
    for (int i = 0; i < v.n(); ++i) h->GetXaxis()->SetBinLabel(i + 1, v.lab[i].c_str());
    h->GetXaxis()->SetTitle(v.title.c_str()); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetLabelSize((v.n() > 9 ? 1.0 : 1.25) * ts); h->GetXaxis()->SetTitleSize(ts); h->GetXaxis()->SetTitleOffset(1.15);
    h->GetYaxis()->SetLabelSize(ts); h->GetYaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleOffset(1.15);
    h->SetMinimum(y0); h->SetMaximum(y1);
    h->Draw("axis");
    return h;
}
void Deco(TGraph* g, int col, int mk, double ms = 1.3, int lsty = 1) {
    g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineStyle(lsty); g->SetLineWidth(2);
}
void DecoH(TH1* h, int col, int mk, double ms = 0.9) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->GetXaxis()->SetTitle(kXt); h->GetYaxis()->SetTitle(kYt); h->GetXaxis()->SetNdivisions(506);
}

// ---- one-dimensional scans: F_stat and S for the six variables -------------------
// sc[0] = all rows, sc[1], sc[2] = folds (S x 2, F x sqrt 2 to compare).
void DrawScans(const char* fname, const char* head, const std::vector<Scan6>& sc, const Pt& base, bool pdf) {
    TCanvas* c = new TCanvas(Form("c_%s", fname), "", 3000, 2700);
    TPad* top = new TPad(Form("top_%s", fname), "", 0, 0.94, 1, 1); top->Draw();
    TPad* body = new TPad(Form("body_%s", fname), "", 0, 0, 1, 0.94); body->Draw();
    top->cd();
    Txt(0.01, 0.74, head, 0.34, kBlack, 12);
    Txt(0.01, 0.42, "p+p 4.53 GeV, isBest #geq 0, corrected momenta, fit P_main (frozen peak shape).  Blue: all rows; grey, orange open: the two halves of the events (S #times 2, F #times #sqrt{2})", 0.25, kBlack, 12);
    Txt(0.01, 0.13, "green: plateau (allowed values with F #geq F_{best} #minus 1) and its centre; hatched: flagged by the sculpting guard, not used; black tick: value of the base point", 0.25, kBlack, 12);
    body->Divide(3, 4, 0.001, 0.001);
    for (int v = 0; v < kNV; ++v) {
        const VarDef& V = gVar[v];
        const int n = V.n();
        std::vector<double> F(n);
        double f0 = 1e9, f1 = -1e9, s0 = 1e9, s1 = -1e9;
        for (int i = 0; i < n; ++i) {
            F[i] = sc[0][v][i].pass.F;
            for (size_t k = 0; k < sc.size(); ++k) {
                const double sf = k == 0 ? 1. : std::sqrt(2.), ss = k == 0 ? 1. : 2.;
                const Ev& e = sc[k][v][i].pass;
                f0 = std::min(f0, sf * e.F); f1 = std::max(f1, sf * e.F);
                s0 = std::min(s0, ss * (e.S - e.eS)); s1 = std::max(s1, ss * (e.S + e.eS));
            }
        }
        const Plat pl = Plateau1D(F, gFlag[v], V.loose);
        for (int row = 0; row < 2; ++row) {
            body->cd((v / 3) * 6 + row * 3 + v % 3 + 1);
            gPad->SetLeftMargin(0.155); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.07); gPad->SetBottomMargin(0.19);
            const double y0 = row == 0 ? f0 - 0.4 : s0 - 0.05 * (s1 - s0), y1 = row == 0 ? f1 + 0.9 : s1 + 0.22 * (s1 - s0);
            TH1D* fr = CatFrame(V, y0, y1, row == 0 ? "F_{stat} = S / #sigma_{S} (stat. only)" : "peak content S [rows]", 0.062);
            fr->GetYaxis()->SetTitleOffset(1.2);
            for (int i = 0; i < n; ++i) if (gFlag[v][i]) Box(i + 0.05, y0, i + 0.95, y1, kGray + 1, 3354);
            if (row == 0) {
                Box(pl.lo + 0.05, pl.fmax - 1.0, pl.hi + 0.95, pl.fmax, kGreen - 9);
                HLine(0, n, pl.fmax - 1.0, kGreen + 2, 2);
                VLine(pl.ctr + 0.5, y0, pl.fmax, kGreen + 2, 1, 3);
            }
            VLine(base.i[v] + 0.5, y0, y0 + 0.12 * (y1 - y0), kBlack, 1, 4);
            for (int k = (int)sc.size() - 1; k >= 0; --k) {
                TGraphErrors* g = new TGraphErrors();
                const double sf = k == 0 ? 1. : std::sqrt(2.), ss = k == 0 ? 1. : 2.;
                for (int i = 0; i < n; ++i) {
                    const Ev& e = sc[k][v][i].pass;
                    const double x = i + 0.5 + (k == 0 ? 0. : (k == 1 ? -0.2 : 0.2));
                    g->SetPoint(i, x, row == 0 ? sf * e.F : ss * e.S);
                    g->SetPointError(i, 0, row == 0 ? 0. : ss * e.eS);
                }
                Deco(g, k == 0 ? kColS : (k == 1 ? kColF0 : kColF1), k == 0 ? 20 : (k == 1 ? 24 : 25), k == 0 ? 1.8 : 1.3);
                g->Draw(k == 0 ? "lp same" : "p same");
            }
            Txt(0.185, 0.86, Form("%s", row == 0 ? "F_{stat}" : "S #pm #sigma_{S}"), 0.065, kColS);
            if (row == 0) Txt(0.97, 0.86, Form("plateau %s to %s, centre %s", V.lab[pl.lo].c_str(), V.lab[pl.hi].c_str(), V.lab[pl.ctr].c_str()), 0.058, kGreen + 3, 31);
            gPad->RedrawAxis();
        }
    }
    Save(c, fname, pdf);
}

// ---- efficiencies of the single cuts ---------------------------------------------
void DrawEff(const char* fname, const char* head, const Scan6& sc, bool pdf) {
    TCanvas* c = new TCanvas(Form("c_%s", fname), "", 3000, 2700);
    TPad* top = new TPad(Form("top_%s", fname), "", 0, 0.94, 1, 1); top->Draw();
    TPad* body = new TPad(Form("body_%s", fname), "", 0, 0, 1, 0.94); body->Draw();
    top->cd(); Txt(0.01, 0.74, head, 0.34, kBlack, 12);
    Txt(0.01, 0.42, "rows 1, 3: signal efficiency #varepsilon_{S} = S_{pass}/(S_{pass}+S_{fail}) from pass/fail fits (blue) against #sqrt{#varepsilon_{B}} (black line): a cut helps the statistical precision only above the line", 0.25, kBlack, 12);
    Txt(0.01, 0.13, "rows 2, 4: S/B (blue) and CB share of the rows in the peak window (red); hatched: flagged by the sculpting guard (there #varepsilon_{S} is not reliable)", 0.25, kBlack, 12);
    body->Divide(3, 4, 0.001, 0.001);
    for (int v = 0; v < kNV; ++v) {
        const VarDef& V = gVar[v]; const int n = V.n();
        body->cd((v / 3) * 6 + v % 3 + 1);
        gPad->SetLeftMargin(0.155); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.07); gPad->SetBottomMargin(0.19);
        CatFrame(V, 0.2, 1.15, "#varepsilon_{S}  and  #sqrt{#varepsilon_{B}}", 0.062)->GetYaxis()->SetTitleOffset(1.2);
        for (int i = 0; i < n; ++i) if (gFlag[v][i]) Box(i + 0.05, 0.2, i + 0.95, 1.15, kGray + 1, 3354);
        TGraphErrors* ge = new TGraphErrors(); TGraph* gb = new TGraph();
        for (int i = 0; i < n; ++i) {
            ge->SetPoint(i, i + 0.5, sc[v][i].epsS); ge->SetPointError(i, 0, sc[v][i].eEpsS);
            gb->SetPoint(i, i + 0.5, sc[v][i].sqB);
        }
        Deco(ge, kColS, 20, 1.8); Deco(gb, kBlack, 1, 1., 1);
        gb->Draw("l same"); ge->Draw("p same");
        gPad->RedrawAxis();
        body->cd((v / 3) * 6 + 3 + v % 3 + 1);
        gPad->SetLeftMargin(0.155); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.07); gPad->SetBottomMargin(0.19);
        double ym = 0;
        for (int i = 0; i < n; ++i) ym = std::max(ym, std::max(sc[v][i].pass.SoB, sc[v][i].pass.cbs));
        CatFrame(V, 0, 1.25 * ym, "S/B  and  CB/N (no k factor)", 0.062)->GetYaxis()->SetTitleOffset(1.2);
        for (int i = 0; i < n; ++i) if (gFlag[v][i]) Box(i + 0.05, 0, i + 0.95, 1.25 * ym, kGray + 1, 3354);
        TGraph* g1 = new TGraph(); TGraph* g2 = new TGraph();
        for (int i = 0; i < n; ++i) { g1->SetPoint(i, i + 0.5, sc[v][i].pass.SoB); g2->SetPoint(i, i + 0.5, sc[v][i].pass.cbs); }
        Deco(g1, kColS, 20, 1.8); Deco(g2, kRed + 1, 21, 1.6);
        g2->Draw("lp same"); g1->Draw("lp same");
        TLegend* lg = Leg(0.19, 0.80, 0.80, 0.90, 0.060, 2);
        lg->AddEntry(g1, "S/B", "pl"); lg->AddEntry(g2, "CB/N", "pl"); lg->Draw();
        gPad->RedrawAxis();
    }
    Save(c, fname, pdf);
}

// ---- two-dimensional maps -----------------------------------------------------------
struct Map2D { int vx = 0, vy = 0; Pt base; std::vector<std::vector<Ev>> ev; double fmax = 0; int bx = 0, by = 0; };
void DrawMaps(const char* fname, const char* head, std::vector<Map2D>& maps, bool drawS, bool pdf) {
    gStyle->SetPaintTextFormat(drawS ? "4.0f" : "4.1f");
    TCanvas* c = new TCanvas(Form("c_%s", fname), "", 3000, 2500);
    TPad* top = new TPad(Form("top_%s", fname), "", 0, 0.95, 1, 1); top->Draw();
    TPad* body = new TPad(Form("body_%s", fname), "", 0, 0, 1, 0.95); body->Draw();
    top->cd(); Txt(0.01, 0.5, head, 0.50, kBlack, 12);
    body->Divide(2, 2, 0.001, 0.001);
    for (size_t m = 0; m < maps.size(); ++m) {
        Map2D& M = maps[m];
        const VarDef& X = gVar[M.vx]; const VarDef& Y = gVar[M.vy];
        body->cd((int)m + 1);
        gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13); gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.14);
        TH2D* h = new TH2D(Form("%s_%zu", fname, m), "", X.n(), 0, X.n(), Y.n(), 0, Y.n());
        TH2D* hf = new TH2D(Form("%s_f%zu", fname, m), "", X.n(), 0, X.n(), Y.n(), 0, Y.n());
        h->SetDirectory(nullptr); hf->SetDirectory(nullptr);
        M.fmax = -1e9;
        for (int i = 0; i < X.n(); ++i)
            for (int j = 0; j < Y.n(); ++j) {
                const bool fl = gFlag[M.vx][i] || gFlag[M.vy][j];
                if (!fl && M.ev[i][j].F > M.fmax) { M.fmax = M.ev[i][j].F; M.bx = i; M.by = j; }
            }
        double zmin = 1e9, zmax = -1e9;
        for (int i = 0; i < X.n(); ++i)
            for (int j = 0; j < Y.n(); ++j) {
                const Ev& e = M.ev[i][j];
                const bool fl = gFlag[M.vx][i] || gFlag[M.vy][j];
                h->SetBinContent(i + 1, j + 1, drawS ? e.S : e.F);
                hf->SetBinContent(i + 1, j + 1, fl ? M.fmax - 3.0 : e.F);
                if (!fl) { zmin = std::min(zmin, drawS ? e.S : e.F); zmax = std::max(zmax, drawS ? e.S : e.F); }
            }
        for (int i = 0; i < X.n(); ++i) h->GetXaxis()->SetBinLabel(i + 1, X.lab[i].c_str());
        for (int j = 0; j < Y.n(); ++j) h->GetYaxis()->SetBinLabel(j + 1, Y.lab[j].c_str());
        h->GetXaxis()->SetTitle(X.title.c_str()); h->GetYaxis()->SetTitle(Y.title.c_str());
        h->GetXaxis()->SetLabelSize(0.055); h->GetYaxis()->SetLabelSize(0.055);
        h->GetXaxis()->SetTitleSize(0.040); h->GetYaxis()->SetTitleSize(0.040);
        h->GetXaxis()->SetTitleOffset(1.45); h->GetYaxis()->SetTitleOffset(1.5);
        h->GetZaxis()->SetLabelSize(0.04);
        h->SetMinimum(drawS ? zmin - 0.02 * (zmax - zmin) : std::max(zmin, M.fmax - 3.0)); h->SetMaximum(zmax + 0.02 * (zmax - zmin));
        h->SetMarkerSize(1.25);
        h->Draw("colz text");
        for (int i = 0; i < X.n(); ++i)
            for (int j = 0; j < Y.n(); ++j)
                if (gFlag[M.vx][i] || gFlag[M.vy][j]) Box(i, j, i + 1, j + 1, kGray + 2, 3354);
        // contour F = F_best - 1: edges between cells inside and outside the plateau
        for (int i = 0; i < X.n(); ++i)
            for (int j = 0; j < Y.n(); ++j) {
                const bool in = hf->GetBinContent(i + 1, j + 1) >= M.fmax - 1.0;
                if (!in) continue;
                const bool l = i > 0 && hf->GetBinContent(i, j + 1) >= M.fmax - 1.0;
                const bool r = i < X.n() - 1 && hf->GetBinContent(i + 2, j + 1) >= M.fmax - 1.0;
                const bool dn = j > 0 && hf->GetBinContent(i + 1, j) >= M.fmax - 1.0;
                const bool up = j < Y.n() - 1 && hf->GetBinContent(i + 1, j + 2) >= M.fmax - 1.0;
                if (!l) VLine(i, j, j + 1, kBlack, 1, 5);
                if (!r) VLine(i + 1, j, j + 1, kBlack, 1, 5);
                if (!dn) HLine(i, i + 1, j, kBlack, 1, 5);
                if (!up) HLine(i, i + 1, j + 1, kBlack, 1, 5);
            }
        std::string oth;
        for (int v = 0; v < kNV; ++v) {
            if (v == M.vx || v == M.vy || M.base.i[v] == gVar[v].none) continue;
            if (!oth.empty()) oth += ", ";
            oth += gVar[v].name + " " + gVar[v].lab[M.base.i[v]];
        }
        Txt(0.13, 0.955, Form("%s;  other cuts: %s", drawS ? "peak content S [rows]" : "F_{stat} = S/#sigma_{S} (stat. only, not a significance)", oth.empty() ? "none" : oth.c_str()), 0.034);
        Txt(0.13, 0.915, Form("black line: region F_{stat} #geq F_{best} #minus 1 (F_{best} = %.2f); hatched: flagged by the sculpting guard", M.fmax), 0.030);
        gPad->RedrawAxis();
    }
    Save(c, fname, pdf);
}

// ---- fitted spectra -------------------------------------------------------------------
TGraph* Curve(const scan::PeakFit& f, int what, int col, int sty) {      // 0 total, 1 background, 2 Gaussian
    TGraph* g = new TGraph();
    for (int i = 0; i <= 400; ++i) {
        const double m = f.cfg.lo + (f.cfg.hi - f.cfg.lo) * i / 400;
        g->SetPoint(i, m, what == 0 ? f.Total(m) : (what == 1 ? f.Bkg(m) : f.Total(m) - f.Bkg(m)));
    }
    g->SetLineColor(col); g->SetLineStyle(sty); g->SetLineWidth(3);
    return g;
}
void DrawFits(const char* fname, const char* head, const std::vector<const Full*>& fs, bool pdf) {
    const int n = (int)fs.size();
    TCanvas* c = new TCanvas(Form("c_%s", fname), "", 1000 * n, 2700);
    TPad* top = new TPad(Form("top_%s", fname), "", 0, 0.96, 1, 1); top->Draw();
    TPad* body = new TPad(Form("body_%s", fname), "", 0, 0, 1, 0.96); body->Draw();
    top->cd(); Txt(0.01, 0.5, head, 0.50, kBlack, 12);
    body->Divide(n, 3, 0.001, 0.001);
    for (int k = 0; k < n; ++k) {
        const Full& f = *fs[k];
        TH1D* hn = f.sp.Hist("us", Form("%s_us%d", fname, k)); TH1D* hc = f.sp.Hist("cb", Form("%s_cb%d", fname, k));
        TH1D* hs = f.sp.Hist("sig", Form("%s_sig%d", fname, k));
        DecoH(hn, kBlack, 20); DecoH(hc, kRed + 1, 20); DecoH(hs, kBlue + 1, 20);
        // row 1: all pairs, CB, N - CB
        body->cd(k + 1);
        gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.12);
        hn->GetXaxis()->SetRangeUser(0.4501, 0.6499);
        double top1 = 0;
        for (int i = hn->FindBin(0.4501); i <= hn->FindBin(0.6499); ++i) top1 = std::max(top1, hn->GetBinContent(i));
        hn->SetMinimum(0); hn->SetMaximum(1.55 * top1);
        hn->GetYaxis()->SetTitleOffset(1.75);
        hn->Draw("e"); hc->Draw("e same"); hs->Draw("e same");
        Txt(0.19, 0.955, Form("#bf{%s}", f.name.c_str()), 0.055);
        Txt(0.19, 0.895, f.text.c_str(), 0.040);
        TLegend* lg = Leg(0.22, 0.66, 0.77, 0.86, 0.045);
        lg->AddEntry(hn, "all pairs e^{+}e^{-} (N)", "pe"); lg->AddEntry(hc, "CB = 2#sqrt{N_{++}N_{--}}", "pe");
        lg->AddEntry(hs, "N #minus CB", "pe"); lg->Draw();
        VLine(scan::kPeakLo, 0, 0.55 * hn->GetMaximum(), kGray + 2, 2); VLine(scan::kPeakHi, 0, 0.55 * hn->GetMaximum(), kGray + 2, 2);
        // row 2: zoom on N - CB with the fit
        body->cd(n + k + 1);
        gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04);
        TH1D* hz = (TH1D*)hs->Clone(Form("%s_z%d", fname, k));
        hz->GetXaxis()->SetRangeUser(0.5001, 0.6299);
        double a = 1e12, b = -1e12;
        for (int i = hz->FindBin(0.5001); i <= hz->FindBin(0.6299); ++i) { a = std::min(a, hz->GetBinContent(i)); b = std::max(b, hz->GetBinContent(i)); }
        hz->SetMinimum(a - 0.75 * (b - a)); hz->SetMaximum(b + 0.22 * (b - a));
        hz->GetYaxis()->SetTitle("N #minus CB  [rows / 2.5 MeV/c^{2}]"); hz->GetYaxis()->SetTitleOffset(1.75);
        hz->SetMarkerSize(1.2);
        hz->Draw("e");
        Curve(f.pm, 1, kBlack, 2)->Draw("l same"); Curve(f.pm, 0, kBlack, 1)->Draw("l same");
        hz->Draw("e same");
        Txt(0.23, 0.90, "linear zoom: N #minus CB, fit (solid), pol2 (dashed)", 0.043);
        Txt(0.40, 0.40, Form("S = %.0f #pm %.0f (stat) #pm %.0f (model)", f.pm.S, f.pm.eS, f.sd.half), 0.040);
        Txt(0.40, 0.345, Form("#chi^{2}/ndf = %.1f/%d;  Poisson: %.0f #pm %.0f", f.pm.chi2, f.pm.ndf, f.pp.S, f.pp.eS), 0.040);
        Txt(0.40, 0.29, Form("window: N = %.0f, CB/N = %.3f", f.pm.N, f.pm.CB / f.pm.N), 0.040);
        Txt(0.40, 0.235, Form("S/B = %.3f,  F_{stat} = %.1f,  F_{tot} = %.1f", f.pm.SoB, f.pm.F, f.ftot), 0.040);
        Txt(0.40, 0.18, Form("fake content from sculpting: %.0f #pm %.0f", f.fake.s, f.fake.e), 0.040);
        // row 3: after subtracting the polynomial
        body->cd(2 * n + k + 1);
        gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04);
        TH1D* hr = new TH1D(Form("%s_r%d", fname, k), "", f.sp.nbin, f.sp.lo, f.sp.hi); hr->SetDirectory(nullptr);
        for (int i = 0; i < f.sp.nbin; ++i) { hr->SetBinContent(i + 1, f.sp.sig[i] - f.pm.Bkg(f.sp.x(i))); hr->SetBinError(i + 1, std::sqrt(f.sp.sigVar[i])); }
        DecoH(hr, kBlue + 1, 20, 1.2);
        hr->GetXaxis()->SetRangeUser(0.5001, 0.6299);
        const double pk = std::max(f.pm.Total(scan::kPeakMean) - f.pm.Bkg(scan::kPeakMean), 200.);
        hr->SetMinimum(-0.9 * pk); hr->SetMaximum(1.9 * pk);
        hr->GetYaxis()->SetTitle("N #minus CB #minus pol2  [rows / 2.5 MeV/c^{2}]"); hr->GetYaxis()->SetTitleOffset(1.75);
        hr->Draw("e");
        HLine(0.50, 0.63, 0, kBlack, 3, 1);
        Curve(f.pm, 2, kBlack, 1)->Draw("l same");
        VLine(scan::kPeakLo, -0.9 * pk, 1.3 * pk, kGray + 2, 2); VLine(scan::kPeakHi, -0.9 * pk, 1.3 * pk, kGray + 2, 2);
        VLine(scan::kFitLo, -0.9 * pk, 0.5 * pk, kGray + 2, 3); VLine(scan::kFitHi, -0.9 * pk, 0.5 * pk, kGray + 2, 3);
        Txt(0.23, 0.90, "pol2 subtracted; dashed: window, dotted: fit range", 0.040);
        Txt(0.23, 0.845, Form("null test on like-sign pairs: %.0f #pm %.0f (++), %.0f #pm %.0f (#minus#minus)", f.lsP.S, f.lsP.eS, f.lsM.S, f.lsM.eS), 0.037);
    }
    Save(c, fname, pdf);
}

// ---- table of a one-dimensional scan ----------------------------------------------
void ScanTable(const char* file, const char* title, const std::vector<Scan6>& sc) {
    OpenTable(file);
    P("%s\n", title);
    P("pass = rows passing the cut, fail = rows of the same base failing it. eps_S = S_pass/(S_pass+S_fail); eps_B from polynomial + CB in the window.\n");
    P("z = (eps_S - sqrt(eps_B)) / error of eps_S: the cut helps the statistical precision only if z > 0. clos = S_pass + S_fail - S_base.\n");
    P("eps_S is not reliable where the fail fit is bad (chi2_fail large: the removed rows have a threshold inside the fit range).\n");
    P("f0 / f1 = the two halves of the event groups (parity of the group index). flag = 1: value flagged by the sculpting guard (peak_guard.txt), not used.\n\n");
    for (int v = 0; v < kNV; ++v) {
        P("%-9s %5s %4s | %7s %6s | %6s %5s %6s %6s %8s | %6s %5s %8s | %6s %6s %6s %6s %6s | %5s | %6s %5s %6s | %6s %5s %6s\n",
          gVar[v].key.c_str(), "value", "flag", "N", "CB/N", "S", "err", "F", "S/B", "chi2/ndf", "S_fail", "err", "chi2/ndf",
          "eps_S", "err", "eps_B", "sqrt", "z", "clos", "S f0", "err", "F f0", "S f1", "err", "F f1");
        for (int i = 0; i < gVar[v].n(); ++i) {
            const ScanPt& q = sc[0][v][i];
            P("%-9s %5s %4d | %7.0f %6.3f | %6.0f %5.0f %6.2f %6.4f %5.1f/%-2d | %6.0f %5.0f %5.1f/%-2d | %6.3f %6.3f %6.3f %6.3f %6.2f | %5.0f | %6.0f %5.0f %6.2f | %6.0f %5.0f %6.2f\n",
              gVar[v].key.c_str(), gVar[v].lab[i].c_str(), gFlag[v][i], q.pass.N, q.pass.cbs, q.pass.S, q.pass.eS, q.pass.F, q.pass.SoB,
              q.pass.chi2, q.pass.ndf, q.fail.S, q.fail.eS, q.fail.chi2, q.fail.ndf, q.epsS, q.eEpsS, q.epsB, q.sqB, q.z, q.clos,
              sc[1][v][i].pass.S, sc[1][v][i].pass.eS, sc[1][v][i].pass.F, sc[2][v][i].pass.S, sc[2][v][i].pass.eS, sc[2][v][i].pass.F);
        }
        std::vector<double> F;
        for (const auto& x : sc[0][v]) F.push_back(x.pass.F);
        const Plat pl = Plateau1D(F, gFlag[v], gVar[v].loose);
        P("   allowed values: best F = %.2f at %s; plateau (F >= best - 1): %s to %s; centre %s\n\n", pl.fmax, gVar[v].lab[pl.best].c_str(),
          gVar[v].lab[pl.lo].c_str(), gVar[v].lab[pl.hi].c_str(), gVar[v].lab[pl.ctr].c_str());
    }
    CloseTable();
}

// ---- output tree of all evaluated points --------------------------------------------
struct TreeOut {
    TFile* f = nullptr; TTree* t = nullptr;
    char set[40], name[80];
    Int_t sample, idx[kNV], flagged, ndf;
    Float_t val[kNV], N, CB, B, S, eS, F, SoB, cbShare, chi2, epsS, eEpsS, epsB, spread, Ftot, Spois, eSpois, Sfake, eSfake;
    void Open(const char* path) {
        f = TFile::Open(path, "RECREATE");
        t = new TTree("wp", "eta peak: evaluated selections (P_main, frozen shape)");
        t->Branch("set", set, "set/C"); t->Branch("name", name, "name/C");
        t->Branch("sample", &sample, "sample/I");
        t->Branch("idx", idx, Form("idx[%d]/I", kNV)); t->Branch("val", val, Form("val[%d]/F", kNV));
        t->Branch("flagged", &flagged, "flagged/I");
        t->Branch("N", &N, "N/F"); t->Branch("CB", &CB, "CB/F"); t->Branch("B", &B, "B/F");
        t->Branch("S", &S, "S/F"); t->Branch("eS", &eS, "eS/F"); t->Branch("F", &F, "F/F");
        t->Branch("SoB", &SoB, "SoB/F"); t->Branch("cbShare", &cbShare, "cbShare/F");
        t->Branch("chi2", &chi2, "chi2/F"); t->Branch("ndf", &ndf, "ndf/I");
        t->Branch("epsS", &epsS, "epsS/F"); t->Branch("eEpsS", &eEpsS, "eEpsS/F"); t->Branch("epsB", &epsB, "epsB/F");
        t->Branch("spread", &spread, "spread/F"); t->Branch("Ftot", &Ftot, "Ftot/F");
        t->Branch("Spois", &Spois, "Spois/F"); t->Branch("eSpois", &eSpois, "eSpois/F");
        t->Branch("Sfake", &Sfake, "Sfake/F"); t->Branch("eSfake", &eSfake, "eSfake/F");
        t->SetTitle("eta peak selections. val[] = oa_lo, oa_hi, oa_lab, mpipi_lo, mpipi_hi, mm4_hi (-999 = no cut); sample 0 = all, 1 = fold 0, 2 = fold 1; -1 = not filled");
    }
    void Fill(const char* s, const char* n, int smp, const Pt* p, const Ev& e, const ScanPt* q = nullptr, const Full* fu = nullptr) {
        snprintf(set, sizeof(set), "%s", s); snprintf(name, sizeof(name), "%s", n);
        sample = smp; flagged = 0;
        for (int v = 0; v < kNV; ++v) {
            idx[v] = p ? p->i[v] : -1;
            val[v] = (p && p->i[v] != gVar[v].none) ? (float)gVar[v].val[p->i[v]] : -999.f;
            if (p && gFlag[v][p->i[v]]) flagged = 1;
        }
        N = e.N; CB = e.CB; B = e.B; S = e.S; eS = e.eS; F = e.F; SoB = e.SoB; cbShare = e.cbs; chi2 = e.chi2; ndf = e.ndf;
        epsS = q ? q->epsS : -1; eEpsS = q ? q->eEpsS : -1; epsB = q ? q->epsB : -1;
        spread = fu ? fu->sd.half : -1; Ftot = fu ? fu->ftot : -1; Spois = fu ? fu->pp.S : -1; eSpois = fu ? fu->pp.eS : -1;
        Sfake = fu ? fu->fake.s : -1; eSfake = fu ? fu->fake.e : -1;
        t->Fill();
    }
    void Close() { f->cd(); t->Write(); f->Close(); }
};

} // namespace spk

// nSplit = number of random halvings of the events used for the bias measurement
void scan_peak(int nSplit = 40) {
    using namespace spk;
    TStopwatch swAll, sw;
    Style();
    gSystem->mkdir(scan::kFigDir, kTRUE);
    BuildGrid();

    scan::Data d;
    if (!scan::Load(d)) return;
    // rows around the fit range only: enough for P_main in the frozen range
    // (the CB interpolation needs 0.51-0.62) and about 7 times faster to scan
    scan::Data dz;
    for (int k = 0; k < 3; ++k)
        for (const auto& r : d.s[k]) if (r.m4_fb >= 0.50f && r.m4_fb < 0.63f) dz.s[k].push_back(r);
    {
        const scan::PeakFit a = scan::FitPeakMain(scan::Fill(d, MkSel(PtStart())));
        const scan::PeakFit b = scan::FitPeakMain(scan::Fill(dz, MkSel(PtStart())));
        printf("rows with 0.50 <= M4 < 0.63: us %zu, pp %zu, mm %zu; start point on all rows: S = %.3f +- %.3f, on these rows: S = %.3f +- %.3f\n",
               dz.s[0].size(), dz.s[1].size(), dz.s[2].size(), a.S, a.eS, b.S, b.eS);
    }
    TreeOut to; to.Open(Form("%s/peak_wp.root", scan::kScanDir));

    // ------------------------------------------------------------------------
    // 1. one-dimensional scans: each cut alone; N-1 around the start point
    // ------------------------------------------------------------------------
    const std::vector<Samp> samp3 = {SampAll(), SampFold(0), SampFold(1)};
    const std::vector<Scan6> scSingle = Scan1D(dz, PtNone(), samp3, true);
    const std::vector<Scan6> scNm1 = Scan1D(dz, PtStart(), samp3, true);
    const Ev evNone = scSingle[0][0][gVar[0].none].pass;
    printf("one-dimensional scans: %.1f s\n", sw.RealTime());

    // ------------------------------------------------------------------------
    // 2. guard against sculpting
    // ------------------------------------------------------------------------
    sw.Start();
    GuardHist* gh = new GuardHist();
    for (int k = 0; k < 3; ++k)
        for (const auto& r : d.s[k]) if (r.cor_ok) GuardFill(*gh, k == 0 ? 0 : 1, r, r.m4);
    MixPool pool; BuildPool(d, pool);
    const MixCheck mc = FillMixedGrid(pool, *gh);
    printf("guard: %.1f s; mixed rows in 0.35-0.95: %ld; unmixed rows recomputed: max |dM4| %.1e GeV, |d oa_lab| %.1e deg, |d oa_eta| %.1e deg, |d MM4| %.1e GeV\n",
           sw.RealTime(), mc.nmix, mc.dM4, mc.dOaLab, mc.dOaEta, mc.dMm4);
    const scan::Spectra spNoneZ = scan::Fill(dz, MkSel(PtNone()));
    const scan::PeakFit smoothNone = scan::FitPeakPois(spNoneZ);
    const double fakeThr = 0.25 * evNone.eS;
    struct GRes { Turn tl, tu, tm; Fake fk; int fG = 0, fM = 0, fI = 0; };
    static GRes gres[2][kNV][kMaxIdx];
    {
        scan::PeakFit smooth[2][kNV];
        {
            std::vector<scan::Sel> sels;
            for (int v = 0; v < kNV; ++v) sels.push_back(MkSel(With(PtStart(), v, gVar[v].none)));
            const std::vector<scan::Spectra> sp = scan::FillMany(dz, sels);
            for (int v = 0; v < kNV; ++v) { smooth[0][v] = smoothNone; smooth[1][v] = scan::FitPeakPois(sp[v]); }
        }
        for (int c = 0; c < 2; ++c)
            for (int v = 0; v < kNV; ++v)
                for (int i = 0; i < gVar[v].n(); ++i) {
                    GRes& g = gres[c][v][i];
                    g.tl = TurnOn(gh->pas[1][c][v][i], gh->tot[1][c][v], 8, 300);
                    g.tu = TurnOn(gh->pas[0][c][v][i], gh->tot[0][c][v], 4, 300);
                    g.tm = TurnOn(Add(gh->pas[2][c][v][i], gh->pas[3][c][v][i]), Add(gh->tot[2][c][v], gh->tot[3][c][v]), 4, 300);
                    g.fk = FakeS(smooth[c][v], gh->pas[2][c][v][i], gh->tot[2][c][v], gh->pas[3][c][v][i], gh->tot[3][c][v]);
                    if (i == gVar[v].none) continue;
                    g.fG = g.tl.edgeIn ? 1 : 0;
                    g.fM = std::fabs(g.fk.s) > fakeThr ? 1 : 0;
                }
        // the decision uses context 0 (the cut alone); an unflagged value between
        // two flagged ones is not used either (no robust safe value)
        for (int v = 0; v < kNV; ++v) {
            for (int i = 0; i < gVar[v].n(); ++i) gFlag[v][i] = (gres[0][v][i].fG || gres[0][v][i].fM) ? 1 : 0;
            for (int i = 0; i < gVar[v].n(); ++i) {
                if (gFlag[v][i] || i == gVar[v].none) continue;
                bool lo = false, hi = false;
                for (int j = 0; j < i; ++j) if (gres[0][v][j].fG || gres[0][v][j].fM) lo = true;
                for (int j = i + 1; j < gVar[v].n(); ++j) if (gres[0][v][j].fG || gres[0][v][j].fM) hi = true;
                // a value identical to "no cut" inside the fit range stays allowed
                if (lo && hi && gres[0][v][i].tu.dfit > 1e-4) { gres[0][v][i].fI = 1; gFlag[v][i] = 1; }
            }
        }
        // the mixed spectrum itself fitted with Gaussian + pol2 (is it a pol2 without any cut?)
        const std::vector<double> mt = Add(gh->tot[2][0][0], gh->tot[3][0][0]);
        std::vector<double> coef;
        const double sMix = GaussCoef(mt, mt, &coef);
        double e2 = 0, nMixWin = 0, bWin = 0;
        for (int i = 0; i < scan::kNbin; ++i) {
            e2 += coef[i] * coef[i] * mt[i];
            const double fr = spNoneZ.Frac(i, scan::kPeakLo, scan::kPeakHi);
            nMixWin += fr * mt[i]; bWin += fr * smoothNone.Bkg(spNoneZ.x(i));
        }
        const double kap = bWin / nMixWin;

        OpenTable("peak_guard.txt");
        P("GUARD AGAINST SCULPTING: pass fraction f of each cut value against M4, and the fake peak content a cut makes on a smooth background.\n");
        P("row sets: LS = like-sign rows (++ and -- added), US = unlike-sign rows, MIX = mixed rows: pion pair of one event with the lepton pair of\n");
        P("          another (%d fixed offsets in the row order, %ld rows in 0.35-0.95; M4, oa_lab, missing mass recomputed).\n", kNOff, mc.nmix);
        P("          Recomputation checked on unmixed rows: max |dM4| = %.1e GeV, |d oa_lab| = %.1e deg, |d MM4| = %.1e GeV.\n", mc.dM4, mc.dOaLab, mc.dMm4);
        P("f lo, f hi = pass fraction at the two ends of the fit range (%.4f, %.4f); dfit = its change inside the fit range.\n", scan::kFitLo, scan::kFitHi);
        P("M10-M90 = lowest and highest mass where f crosses 10 %% or 90 %% of its total change between 0.40 and 0.80 (LS: 20 MeV bins; US, MIX: 10 MeV bins).\n");
        P("S_fake  = Gaussian content when (smooth background of the base selection, pol2 from the Poisson fit) x f_MIX is fitted with Gaussian + pol2.\n");
        P("          A pol2 background stays a pol2 only if f is flat or a straight line in the fit range. err = larger of the binomial error and of half\n");
        P("          the difference between the two halves of the mixing offsets.\n");
        P("FLAGS (decision from context 0, the cut alone):\n");
        P("   G = 1: M10 or M90 of the like-sign pass fraction lies inside the fit range and dfit(LS) >= 0.05  (turn-on edge inside the fit range)\n");
        P("   M = 1: |S_fake| > 0.25 sigma_S of the uncut sample = %.0f rows\n", fakeThr);
        P("   I = 1: unflagged value between two flagged ones (no robust safe value); flag = G or M or I: value NOT USED for plateaus and candidates.\n");
        P("check of the method: the mixed spectrum without any cut fitted with Gaussian + pol2 gives %.0f +- %.0f mixed rows = %.0f +- %.0f rows at the\n",
          sMix, std::sqrt(e2), kap * sMix, kap * std::sqrt(e2));
        P("   data normalisation (mixed rows in the window %.0f, smooth data background %.0f).\n", nMixWin, bWin);
        P("check against the data: cuts that remove (almost) no row from the window but change S: M(pipi) < 0.50: S_fake %.0f, observed S(cut) - S(no cut) = %.0f;\n",
          gres[0][kMpHi][10].fk.s, scSingle[0][kMpHi][10].pass.S - evNone.S);
        P("   M(pipi) < 0.48: %.0f, %.0f; MM4 < 2.90: %.0f, %.0f.\n\n", gres[0][kMpHi][8].fk.s, scSingle[0][kMpHi][8].pass.S - evNone.S,
          gres[0][kMm4][6].fk.s, scSingle[0][kMm4][6].pass.S - evNone.S);
        for (int c = 0; c < 2; ++c) {
            const std::vector<Scan6>& sc = c == 0 ? scSingle : scNm1;
            P("%s\n", c == 0 ? "CONTEXT 0: the cut alone" : Form("CONTEXT 1: the other cuts at the start point (%s); S_fake relative to the base of this context; for information", PtText(PtStart()).c_str()));
            P("%-9s %5s | %6s %6s %6s %13s | %6s %6s %6s %13s | %6s %6s %6s %13s | %6s %4s | %8s %5s | %s\n", "variable", "value",
              "LS lo", "f hi", "dfit", "M10-M90", "US lo", "f hi", "dfit", "M10-M90", "MIX lo", "f hi", "dfit", "M10-M90",
              "S_fake", "err", "chi2_fail", "clos", c == 0 ? "G M I flag" : "G M");
            for (int v = 0; v < kNV; ++v) {
                for (int i = 0; i < gVar[v].n(); ++i) {
                    const GRes& g = gres[c][v][i];
                    const ScanPt& q = sc[0][v][i];
                    P("%-9s %5s | %6.3f %6.3f %6.3f %6.3f-%6.3f | %6.3f %6.3f %6.3f %6.3f-%6.3f | %6.3f %6.3f %6.3f %6.3f-%6.3f | %6.0f %4.0f | %5.1f/%-2d %5.0f | ",
                      gVar[v].key.c_str(), gVar[v].lab[i].c_str(), g.tl.fLo, g.tl.fHi, g.tl.dfit, g.tl.mA, g.tl.mB, g.tu.fLo, g.tu.fHi, g.tu.dfit, g.tu.mA, g.tu.mB,
                      g.tm.fLo, g.tm.fHi, g.tm.dfit, g.tm.mA, g.tm.mB, g.fk.s, g.fk.e, q.fail.chi2, q.fail.ndf, q.clos);
                    if (c == 0) P("%d %d %d  %d\n", g.fG, g.fM, g.fI, gFlag[v][i]); else P("%d %d\n", g.fG, g.fM);
                }
                P("\n");
            }
        }
        P("VALUES NOT USED:\n");
        for (int v = 0; v < kNV; ++v) {
            P("   %-10s", gVar[v].name.c_str());
            int nf = 0;
            for (int i = 0; i < gVar[v].n(); ++i) if (gFlag[v][i]) { P(" %s", gVar[v].lab[i].c_str()); nf++; }
            P("%s\n", nf ? "" : " (none flagged)");
        }
        CloseTable();
    }

    // ---- figures of the guard ----
    {
        const int rep[kNV][3] = {{2, 4, 8}, {1, 4, 6}, {2, 5, 9}, {2, 4, 6}, {2, 6, 10}, {2, 4, 6}};
        const int col[3] = {kAzure + 2, kOrange + 7, kGreen + 3};
        TCanvas* c = new TCanvas("c_guard", "", 3000, 2100);
        TPad* top = new TPad("g_top", "", 0, 0.94, 1, 1); top->Draw();
        TPad* body = new TPad("g_body", "", 0, 0, 1, 0.94); body->Draw();
        top->cd();
        Txt(0.01, 0.60, "Sculpting guard: fraction of rows passing a single cut against M(#pi^{+}#pi^{-}e^{+}e^{-}); p+p 4.53 GeV, isBest #geq 0, corrected momenta", 0.40, kBlack, 12);
        Txt(0.01, 0.16, "markers: like-sign rows; solid line: mixed rows; dashed line: unlike-sign rows; yellow: fit range of the peak; a flat or straight line inside the band is harmless; * = flagged, not used", 0.30, kBlack, 12);
        body->Divide(3, 2, 0.001, 0.001);
        for (int v = 0; v < kNV; ++v) {
            body->cd(v + 1);
            gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08);
            TH1D* fr = new TH1D(Form("gfr%d", v), "", 10, 0.40, 0.80); fr->SetDirectory(nullptr);
            fr->SetMinimum(0); fr->SetMaximum(1.25);
            fr->GetXaxis()->SetTitle(kXt); fr->GetYaxis()->SetTitle("fraction of rows passing the cut");
            fr->GetYaxis()->SetTitleOffset(1.0); fr->GetXaxis()->SetNdivisions(508);
            fr->Draw("axis");
            Box(scan::kFitLo, 0, scan::kFitHi, 1.25, kYellow - 9);
            VLine(scan::kPeakLo, 0, 1.05, kGray + 2, 2); VLine(scan::kPeakHi, 0, 1.05, kGray + 2, 2);
            TLegend* lg = Leg(0.13, 0.835, 0.975, 0.915, 0.035, 3);
            for (int q = 0; q < 3; ++q) {
                const int i = rep[v][q];
                TGraphErrors* gl = new TGraphErrors(); TGraph* gm = new TGraph(); TGraph* gu = new TGraph();
                for (int b = 0; b + 8 <= scan::kNbin; b += 8) {
                    const double xc = scan::kMlo + (b + 4) * scan::kBinW;
                    if (xc < 0.40 || xc > 0.80) continue;
                    double a = 0, t = 0;
                    for (int j = 0; j < 8; ++j) { a += gh->pas[1][0][v][i][b + j]; t += gh->tot[1][0][v][b + j]; }
                    if (t < 100) continue;
                    const int np = gl->GetN();
                    gl->SetPoint(np, xc, a / t); gl->SetPointError(np, 0, std::sqrt(std::max(a / t * (1 - a / t), 0.) / t));
                }
                for (int b = 0; b + 4 <= scan::kNbin; b += 4) {
                    const double xc = scan::kMlo + (b + 2) * scan::kBinW;
                    if (xc < 0.40 || xc > 0.80) continue;
                    double a = 0, t = 0, au = 0, tu = 0;
                    for (int j = 0; j < 4; ++j) {
                        a += gh->pas[2][0][v][i][b + j] + gh->pas[3][0][v][i][b + j]; t += gh->tot[2][0][v][b + j] + gh->tot[3][0][v][b + j];
                        au += gh->pas[0][0][v][i][b + j]; tu += gh->tot[0][0][v][b + j];
                    }
                    if (t > 0) gm->SetPoint(gm->GetN(), xc, a / t);
                    if (tu > 0) gu->SetPoint(gu->GetN(), xc, au / tu);
                }
                Deco(gl, col[q], 20, 1.3); Deco(gm, col[q], 1, 1., 1); Deco(gu, col[q], 1, 1., 7);
                gm->Draw("l same"); gu->Draw("l same"); gl->Draw("p same");
                lg->AddEntry(gl, Form("%s %s%s", gVar[v].name.c_str(), gVar[v].lab[i].c_str(), gFlag[v][i] ? " *" : ""), "pl");
            }
            lg->Draw();
            Txt(0.12, 0.94, gVar[v].title.c_str(), 0.040);
            gPad->RedrawAxis();
        }
        Save(c, "peak_guard_passfrac", true);

        TCanvas* c2 = new TCanvas("c_fake", "", 3000, 1900);
        TPad* top2 = new TPad("f_top", "", 0, 0.93, 1, 1); top2->Draw();
        TPad* body2 = new TPad("f_body", "", 0, 0, 1, 0.93); body2->Draw();
        top2->cd();
        Txt(0.01, 0.60, "Sculpting guard: fake peak content S_{fake} = Gaussian part of (smooth pol2 background) #times (pass fraction of mixed rows), fitted like the data", 0.34, kBlack, 12);
        Txt(0.01, 0.16, Form("blue: the cut alone (error from the mixed sample); green band: #pm 0.25 #sigma_{S} = #pm %.0f rows; hatched: flagged; open grey (lower row, right): S(cut) #minus S(no cut) in the data, where it also holds the true signal loss", fakeThr), 0.26, kBlack, 12);
        body2->Divide(3, 2, 0.001, 0.001);
        for (int v = 0; v < kNV; ++v) {
            body2->cd(v + 1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.17);
            const int n = gVar[v].n();
            CatFrame(gVar[v], -520, 420, "S_{fake} [rows]", 0.055)->GetYaxis()->SetTitleOffset(1.35);
            Box(0, -fakeThr, n, fakeThr, kGreen - 9);
            for (int i = 0; i < n; ++i) if (gFlag[v][i]) Box(i + 0.05, -520, i + 0.95, 420, kGray + 1, 3354);
            HLine(0, n, 0, kBlack, 3, 1);
            TGraphErrors* g = new TGraphErrors(); TGraph* go = new TGraph();
            for (int i = 0; i < n; ++i) {
                g->SetPoint(i, i + 0.5, gres[0][v][i].fk.s); g->SetPointError(i, 0, gres[0][v][i].fk.e);
                go->SetPoint(i, i + 0.5, scSingle[0][v][i].pass.S - evNone.S);
            }
            Deco(g, kColS, 20, 1.6); Deco(go, kGray + 2, 24, 1.4);
            if (v == kMpHi || v == kMm4) go->Draw("p same");
            g->Draw("p same");
            gPad->RedrawAxis();
        }
        Save(c2, "peak_guard_fake", false);
    }

    // ---- tables and figures of the scans (now with the flags) ----
    ScanTable("peak_scan_single.txt", "SINGLE-CUT SCANS: each cut alone, applied to the preselection (no other cut)", scSingle);
    ScanTable("peak_scan_nm1.txt", Form("N-1 SCANS around the start point (%s; M(pipi) < 0.46 removed: flagged as sculpting): one cut varied, the others at the start point", PtText(PtStart()).c_str()), scNm1);
    DrawScans("peak_scan_single", "#eta peak, single-cut scans: each cut alone on the preselection (p+p 4.53 GeV, isBest #geq 0, corrected momenta, P_main with frozen shape)", scSingle, PtNone(), true);
    DrawScans("peak_scan_nm1", "#eta peak, N#minus1 scans around the start point (oa_lab < 50#circ; no M(#pi#pi) cut): one cut varied, the others fixed", scNm1, PtStart(), true);
    DrawEff("peak_scan_eff", "#eta peak, single-cut scans (each cut alone on the preselection): efficiencies, S/B and CB share; p+p 4.53 GeV, isBest #geq 0", scSingle[0], false);
    for (int t = 0; t < 2; ++t)
        for (int s = 0; s < 3; ++s)
            for (int v = 0; v < kNV; ++v)
                for (int i = 0; i < gVar[v].n(); ++i) {
                    const Pt p = With(t == 0 ? PtNone() : PtStart(), v, i);
                    const ScanPt& q = (t == 0 ? scSingle : scNm1)[s][v][i];
                    to.Fill(t == 0 ? "single" : "nm1_start", gVar[v].key.c_str(), s, &p, q.pass, &q);
                }

    // ------------------------------------------------------------------------
    // 3. two-dimensional maps (other cuts: start point without flagged values)
    // ------------------------------------------------------------------------
    const Pt base0 = PtBase0();
    std::vector<Map2D> maps(4);
    {
        const int mv[4][2] = {{kOaLab, kMpLo}, {kOaLo, kOaHi}, {kMm4, kOaLab}, {kMpLo, kMpHi}};
        std::vector<Pt> pts;
        for (int m = 0; m < 4; ++m) {
            maps[m].vx = mv[m][0]; maps[m].vy = mv[m][1]; maps[m].base = base0;
            for (int i = 0; i < gVar[mv[m][0]].n(); ++i)
                for (int j = 0; j < gVar[mv[m][1]].n(); ++j) pts.push_back(With(With(base0, mv[m][0], i), mv[m][1], j));
        }
        const std::vector<Ev> ev = EvalPts(dz, pts);
        size_t k = 0;
        OpenTable("peak_maps.txt");
        P("TWO-DIMENSIONAL MAPS of F_stat and S. Other cuts: the start point without flagged values = %s.\n", PtText(base0).c_str());
        P("Rows and columns of flagged values are marked with *.\n\n");
        for (int m = 0; m < 4; ++m) {
            const VarDef& X = gVar[maps[m].vx]; const VarDef& Y = gVar[maps[m].vy];
            maps[m].ev.assign(X.n(), std::vector<Ev>(Y.n()));
            for (int i = 0; i < X.n(); ++i)
                for (int j = 0; j < Y.n(); ++j) {
                    maps[m].ev[i][j] = ev[k];
                    to.Fill(Form("map_%s_%s", X.key.c_str(), Y.key.c_str()), "", 0, &pts[k], ev[k]);
                    k++;
                }
            for (int w = 0; w < 2; ++w) {
                P("%s: columns %s, rows %s\n%10s", w == 0 ? "F_stat" : "S", X.name.c_str(), Y.name.c_str(), "");
                for (int i = 0; i < X.n(); ++i) P(" %6s%s", X.lab[i].c_str(), gFlag[maps[m].vx][i] ? "*" : " ");
                P("\n");
                for (int j = 0; j < Y.n(); ++j) {
                    P("%9s%s", Y.lab[j].c_str(), gFlag[maps[m].vy][j] ? "*" : " ");
                    for (int i = 0; i < X.n(); ++i) P(w == 0 ? " %6.2f " : " %6.0f ", w == 0 ? maps[m].ev[i][j].F : maps[m].ev[i][j].S);
                    P("\n");
                }
                P("\n");
            }
        }
        DrawMaps("peak_map_F", "#eta peak: maps of F_{stat} = S/#sigma_{S} (P_main, frozen shape; p+p 4.53 GeV, isBest #geq 0, corrected momenta)", maps, false, true);
        DrawMaps("peak_map_S", "#eta peak: maps of the peak content S [rows] (P_main, frozen shape; p+p 4.53 GeV, isBest #geq 0, corrected momenta)", maps, true, true);
        for (int m = 0; m < 4; ++m)
            P("map %s x %s: best allowed F = %.2f at %s %s, %s %s\n", gVar[maps[m].vx].key.c_str(), gVar[maps[m].vy].key.c_str(), maps[m].fmax,
              gVar[maps[m].vx].name.c_str(), gVar[maps[m].vx].lab[maps[m].bx].c_str(), gVar[maps[m].vy].name.c_str(), gVar[maps[m].vy].lab[maps[m].by].c_str());
        CloseTable();
    }

    // ------------------------------------------------------------------------
    // 4. candidates on all rows, bias control on halves
    // ------------------------------------------------------------------------
    sw.Start();
    printf("definition of the candidates on all rows:\n");
    const Defn full = Define(dz, true);
    for (size_t k = 0; k < full.grid.size(); ++k) to.Fill("grid", "", 0, &full.grid[k], full.gev[k]);
    const std::vector<Scan6> scNom = Scan1D(dz, full.nominal, samp3, true);
    ScanTable("peak_scan_nm1_nominal.txt", Form("N-1 SCANS around the NOMINAL candidate (%s): one cut varied, the others at the nominal values", PtText(full.nominal).c_str()), scNom);
    DrawScans("peak_scan_nm1_nominal", Form("#eta peak, N#minus1 scans around the nominal candidate (%s)", PtText(full.nominal).c_str()), scNom, full.nominal, true);
    for (int s = 0; s < 3; ++s)
        for (int v = 0; v < kNV; ++v)
            for (int i = 0; i < gVar[v].n(); ++i) {
                const Pt p = With(full.nominal, v, i);
                to.Fill("nm1_nominal", gVar[v].key.c_str(), s, &p, scNom[s][v][i].pass, &scNom[s][v][i]);
            }

    const int kNType = 4;
    const char* typeName[kNType] = {"nominal", "loose", "pure", "best-F"};
    auto typePt = [](const Defn& D, int t) -> Pt { return t == 0 ? D.nominal : (t == 1 ? D.loose : (t == 2 ? D.pure : D.best)); };
    // bias of one halving: selecting half minus independent half, both directions averaged
    struct Bias { double dF[kNType], dS[kNType], dSsig[kNType]; bool same[kNType]; };
    auto biasOf = [&](const scan::Data& A, const scan::Data& B, const Defn& DA, const Defn& DB) {
        std::vector<Pt> pts;
        for (int t = 0; t < kNType; ++t) { pts.push_back(typePt(DA, t)); pts.push_back(typePt(DB, t)); }
        const std::vector<Ev> ea = EvalPts(A, pts), eb = EvalPts(B, pts);
        Bias b;
        for (int t = 0; t < kNType; ++t) {
            const Ev& aA = ea[2 * t]; const Ev& aB = eb[2 * t];         // point of A on A, on B
            const Ev& bA = ea[2 * t + 1]; const Ev& bB = eb[2 * t + 1]; // point of B on A, on B
            b.dF[t] = 0.5 * ((aA.F - aB.F) + (bB.F - bA.F));
            b.dS[t] = 0.5 * ((aA.S - aB.S) + (bB.S - bA.S));
            b.dSsig[t] = b.dS[t] / (0.25 * (aA.eS + aB.eS + bA.eS + bB.eS));
            b.same[t] = typePt(DA, t) == typePt(DB, t);
        }
        return b;
    };

    OpenTable("peak_bias.txt");
    P("BIAS CONTROL. The candidates are defined by a fixed procedure applied to one half of the events and evaluated on the other half.\n");
    P("procedure: nominal = self-consistent plateau centre (N-1 scans around the current point; each variable goes to the centre of the allowed\n");
    P("   values with F >= F_best - 1, or stays open when 'no cut' has the largest F; repeated until stable); start = %s.\n", PtText(base0).c_str());
    P("   grid = %zu combinations of allowed values; loose = most rows in the window among F >= F_best - 1; pure = largest S/B among F >= F_best - 1;\n", full.grid.size());
    P("   best-F = largest F of the grid (the biased choice, shown for comparison only).\n\n");
    Defn Dh[2];
    scan::Data dzF[2];
    std::vector<std::string> wpName[2]; std::vector<Ev> wpSel[2], wpInd[2];
    for (int h = 0; h < 2; ++h) { dzF[h] = Split(dz, SampFold(h)); printf("definition on fold %d:\n", h); Dh[h] = Define(dzF[h], true); }
    int nWp = 0;
    for (int h = 0; h < 2; ++h) {
        std::vector<std::string> nm; std::vector<Pt> pts;
        nm.push_back("fixed: no cut"); pts.push_back(PtNone());
        nm.push_back("fixed: oa_lab<50"); pts.push_back(With(PtNone(), kOaLab, 5));
        nm.push_back("fixed: +M(pipi)<0.46"); pts.push_back(PtOldStart());
        for (int t = 0; t < kNType; ++t) { nm.push_back(typeName[t]); pts.push_back(typePt(Dh[h], t)); }
        for (int v = 0; v < kNV; ++v)
            for (int e = 0; e < 2; ++e) {
                const int i = e == 0 ? Dh[h].plat[v].lo : Dh[h].plat[v].hi;
                if (i == Dh[h].nominal.i[v]) continue;
                if (std::fabs(Dh[h].nm1[v][i].pass.N - Dh[h].evNominal.N) < 0.5) continue;      // the same rows as the nominal point
                nm.push_back(Form("nominal, %s -> %s", gVar[v].key.c_str(), gVar[v].lab[i].c_str())); pts.push_back(With(Dh[h].nominal, v, i));
            }
        std::vector<Ev> es = EvalPts(dzF[h], pts), ei = EvalPts(dzF[1 - h], pts);
        // the two reference chains with the eta-rest angle
        for (int r = 0; r < 2; ++r) {
            nm.push_back(r == 0 ? "fixed: slides" : "fixed: old code");
            const scan::Sel s = r == 0 ? scan::SelSlides() : scan::SelOldCode();
            es.push_back(Eval(scan::Fill(dzF[h], s))); ei.push_back(Eval(scan::Fill(dzF[1 - h], s)));
        }
        P("WORKING POINTS DEFINED ON FOLD %d (%zu points; iterations %d, %s; plateau of the grid: %d of %zu points, F_best = %.2f)\n", h, nm.size(),
          Dh[h].nIter, Dh[h].converged ? "converged" : "stopped on a cycle", Dh[h].nPlat, Dh[h].grid.size(), Dh[h].fbest);
        P("%-28s | %6s %5s %6s | %6s %5s %6s | %7s | %s\n", "point", "S sel", "err", "F sel", "S ind", "err", "F ind", "dF", "cuts");
        for (size_t k = 0; k < nm.size(); ++k) {
            P("%-28s | %6.0f %5.0f %6.2f | %6.0f %5.0f %6.2f | %7.2f | %s\n", nm[k].c_str(), es[k].S, es[k].eS, es[k].F, ei[k].S, ei[k].eS, ei[k].F,
              es[k].F - ei[k].F, k < pts.size() ? PtText(pts[k]).c_str() : (k == pts.size() ? "oa_lab<50, M(pipi)<0.46, eta_rest>140, OA>2" : "oa_lab<50, M(pipi)<0.42, eta_rest>140"));
            to.Fill(Form("wp_def_fold%d", h), nm[k].c_str(), 1 + h, k < pts.size() ? &pts[k] : nullptr, es[k]);
            to.Fill(Form("wp_def_fold%d", h), nm[k].c_str(), 2 - h, k < pts.size() ? &pts[k] : nullptr, ei[k]);
        }
        P("   (sel = on the half that defined the point, ind = on the other half; dF = F sel - F ind also contains the plain difference of the two halves,\n");
        P("    see the fixed points)\n\n");
        wpName[h] = nm; wpSel[h] = es; wpInd[h] = ei; nWp += (int)nm.size();
    }
    const Bias b2 = biasOf(dzF[0], dzF[1], Dh[0], Dh[1]);
    P("MEASURED OPTIMISATION BIAS, two folds: b = 1/2 [ (X_0(c_0) - X_1(c_0)) + (X_1(c_1) - X_0(c_1)) ], c_h = point chosen on fold h, X_h = value on fold h.\n");
    P("   The difference of the two halves cancels; b = 0 when both folds choose the same point. sigma_S is that of one half.\n");
    P("%-10s | %7s %8s %9s | %s\n", "candidate", "b(F)", "b(S)", "b(S)/sig", "same point in both folds");
    for (int t = 0; t < kNType; ++t) P("%-10s | %7.2f %8.0f %9.2f | %s\n", typeName[t], b2.dF[t], b2.dS[t], b2.dSsig[t], b2.same[t] ? "yes" : "no");
    P("\n");
    // random halvings
    double m1[kNType][3], m2[kNType][3]; int nSame[kNType];
    for (int t = 0; t < kNType; ++t) { nSame[t] = 0; for (int j = 0; j < 3; ++j) m1[t][j] = m2[t][j] = 0; }
    std::map<Pt, int> nomCount;
    std::vector<Bias> allB;
    for (int k = 0; k < nSplit; ++k) {
        const scan::Data A = Split(dz, SampRnd(1000 + k, 0)), B = Split(dz, SampRnd(1000 + k, 1));
        const Defn DA = Define(A, false), DB = Define(B, false);
        const Bias b = biasOf(A, B, DA, DB);
        allB.push_back(b);
        nomCount[DA.nominal]++; nomCount[DB.nominal]++;
        for (int t = 0; t < kNType; ++t) {
            const double x[3] = {b.dF[t], b.dS[t], b.dSsig[t]};
            for (int j = 0; j < 3; ++j) { m1[t][j] += x[j]; m2[t][j] += x[j] * x[j]; }
            if (b.same[t]) nSame[t]++;
        }
        if (k % 5 == 4) printf("   random halving %d of %d done (%.0f s)\n", k + 1, nSplit, sw.RealTime()), sw.Continue();
    }
    double bMean[kNType][3], bErr[kNType][3];
    for (int t = 0; t < kNType; ++t)
        for (int j = 0; j < 3; ++j) {
            bMean[t][j] = nSplit > 0 ? m1[t][j] / nSplit : 0;
            const double var = nSplit > 1 ? (m2[t][j] / nSplit - bMean[t][j] * bMean[t][j]) * nSplit / (nSplit - 1) : 0;
            bErr[t][j] = nSplit > 1 ? std::sqrt(std::max(var, 0.) / nSplit) : 0;
        }
    P("MEASURED OPTIMISATION BIAS, %d random halvings of the event groups (same estimator; mean +- rms/sqrt(n); the halvings reuse the same data,\n", nSplit);
    P("   so the errors are too small by an unknown factor below sqrt(n)):\n");
    P("%-10s | %15s | %15s | %17s | %s\n", "candidate", "b(F)", "b(S) rows/half", "b(S)/sigma_S", "same point in both halves");
    for (int t = 0; t < kNType; ++t)
        P("%-10s | %6.2f +- %5.2f | %6.0f +- %5.0f | %7.2f +- %6.2f | %d of %d\n", typeName[t], bMean[t][0], bErr[t][0], bMean[t][1], bErr[t][1],
          bMean[t][2], bErr[t][2], nSame[t], nSplit);
    P("nominal points chosen on the %d random halves:\n", 2 * nSplit);
    for (const auto& kv : nomCount) P("   %3d x  %s\n", kv.second, PtText(kv.first).c_str());
    P("\n");
    printf("bias control: %.1f s\n", sw.RealTime());

    // ---- figure of the bias control ----
    {
        TCanvas* c = new TCanvas("c_bias", "", 3000, 1350);
        TPad* top = new TPad("b_top", "", 0, 0.93, 1, 1); top->Draw();
        TPad* body = new TPad("b_body", "", 0, 0, 1, 0.93); body->Draw();
        top->cd();
        Txt(0.01, 0.62, "#eta peak, bias control: working points defined on one half of the events (filled) and evaluated on the other half (open); P_main, frozen shape", 0.36, kBlack, 12);
        Txt(0.01, 0.17, "left, middle: F_{stat} of each point on the defining half and on the independent half.  Right: measured optimisation bias (selecting minus independent, half-to-half difference cancelled)", 0.29, kBlack, 12);
        body->Divide(3, 1, 0.001, 0.001);
        for (int h = 0; h < 2; ++h) {
            body->cd(h + 1);
            gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.02); gPad->SetBottomMargin(0.40); gPad->SetTopMargin(0.07);
            const int n = (int)wpName[h].size();
            TH1D* fr = new TH1D(Form("bfr%d", h), "", n, 0, n); fr->SetDirectory(nullptr);
            double y0 = 1e9, y1 = -1e9;
            for (int k = 0; k < n; ++k) {
                fr->GetXaxis()->SetBinLabel(k + 1, wpName[h][k].c_str());
                y0 = std::min(y0, std::min(wpSel[h][k].F, wpInd[h][k].F)); y1 = std::max(y1, std::max(wpSel[h][k].F, wpInd[h][k].F));
            }
            fr->LabelsOption("v"); fr->GetXaxis()->SetLabelSize(0.042);
            fr->GetYaxis()->SetTitle("F_{stat} = S/#sigma_{S} (stat. only) on one half"); fr->GetYaxis()->SetTitleOffset(1.15);
            fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045);
            fr->SetMinimum(y0 - 0.4); fr->SetMaximum(y1 + 1.0);
            fr->Draw("axis");
            TGraph* gs = new TGraph(); TGraph* gi = new TGraph();
            for (int k = 0; k < n; ++k) { gs->SetPoint(k, k + 0.5, wpSel[h][k].F); gi->SetPoint(k, k + 0.5, wpInd[h][k].F); }
            Deco(gs, kColS, 20, 1.6); Deco(gi, kColS, 24, 1.6);
            gs->Draw("p same"); gi->Draw("p same");
            TLegend* lg = Leg(0.15, 0.85, 0.97, 0.92, 0.038, 2);
            lg->AddEntry(gs, Form("on fold %d (defining)", h), "p"); lg->AddEntry(gi, Form("on fold %d (independent)", 1 - h), "p");
            lg->Draw();
            Txt(0.11, 0.945, Form("points defined on fold %d", h), 0.042);
        }
        body->cd(3);
        gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.02); gPad->SetBottomMargin(0.40); gPad->SetTopMargin(0.07);
        TH1D* fr = new TH1D("bfr2", "", kNType, 0, kNType); fr->SetDirectory(nullptr);
        for (int t = 0; t < kNType; ++t) fr->GetXaxis()->SetBinLabel(t + 1, typeName[t]);
        fr->GetXaxis()->SetLabelSize(0.06); fr->GetYaxis()->SetTitle("bias of S in units of #sigma_{S} (one half)"); fr->GetYaxis()->SetTitleOffset(1.35);
        fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045);
        fr->SetMinimum(-0.6); fr->SetMaximum(2.2);
        fr->Draw("axis");
        HLine(0, kNType, 0, kBlack, 3, 1);
        TGraphErrors* gr = new TGraphErrors(); TGraph* g2 = new TGraph();
        for (int t = 0; t < kNType; ++t) {
            gr->SetPoint(t, t + 0.4, bMean[t][2]); gr->SetPointError(t, 0, bErr[t][2]);
            g2->SetPoint(t, t + 0.6, b2.dSsig[t]);
        }
        Deco(gr, kColS, 20, 1.8); Deco(g2, kGray + 2, 25, 1.6);
        gr->Draw("p same"); g2->Draw("p same");
        TLegend* lg = Leg(0.18, 0.80, 0.95, 0.92, 0.038);
        lg->AddEntry(gr, Form("mean of %d random halvings", nSplit), "pe"); lg->AddEntry(g2, "folds 0 / 1 (one halving)", "p");
        lg->Draw();
        Txt(0.14, 0.945, "measured optimisation bias", 0.042);
        Save(c, "peak_bias", false);
    }

    // ------------------------------------------------------------------------
    // 5. candidates and reference chains: full fits on all rows
    // ------------------------------------------------------------------------
    sw.Start();
    std::vector<Full> fu;
    auto addPt = [&](const char* name, const Pt& p) {
        Full f; f.name = name; f.text = PtText(p); f.sel = MkSel(p); f.pt = p; f.hasPt = true; fu.push_back(f);
    };
    auto addSel = [&](const char* name, const char* text, const scan::Sel& s) { Full f; f.name = name; f.text = text; f.sel = s; fu.push_back(f); };
    addPt("no cut", PtNone());                                           // 0
    addPt("oa_lab < 50", With(PtNone(), kOaLab, 5));                     // 1
    addPt("+Mpipi<0.46 (flag)", PtOldStart());                           // 2: start point of the first version, flagged as sculpting
    addSel("slides", "oa_lab < 50, M(pipi) < 0.46, eta_rest > 140, OA > 2", scan::SelSlides());      // 3
    addSel("old code", "oa_lab < 50, M(pipi) < 0.42, eta_rest > 140", scan::SelOldCode());           // 4
    addPt("loose", full.loose);                                          // 5
    addPt("nominal", full.nominal);                                      // 6
    addPt("pure", full.pure);                                            // 7
    addPt("best-F grid point", full.best);                               // 8
    const int nFu = (int)fu.size(), iNone = 0, iLoose = 5, iNom = 6, iPure = 7;
    Fake diag[4];
    {
        std::vector<scan::Sel> sels;
        for (int k = 0; k < nFu; ++k) {
            const scan::Sel s = fu[k].sel;
            sels.push_back(s);
            for (int h = 0; h < 2; ++h) sels.push_back([s, h](const scan::Row& r) { return r.fold == h && s(r); });
            sels.push_back([s](const scan::Row& r) { return !s(r); });
        }
        const std::vector<scan::Spectra> sp = scan::FillMany(d, sels);
        std::vector<scan::Sel> msel;
        for (int k = 0; k < nFu; ++k) msel.push_back(fu[k].sel);
        msel.push_back([](const scan::Row& r) { return r.oa_eta > 140.f; });
        msel.push_back([](const scan::Row& r) { return r.oa_lab < 50.f && r.mpipi < 0.42f; });
        msel.push_back([](const scan::Row& r) { return r.oa_lab < 50.f && r.oa > 2.f; });
        msel.push_back([](const scan::Row& r) { return r.oa_lab < 50.f && r.oa_eta > 140.f; });
        const MixSpec ms = MixedFill(pool, msel);
        for (int k = 0; k < 4; ++k) diag[k] = FakeS(smoothNone, ms.pas[0][nFu + k], ms.tot[0], ms.pas[1][nFu + k], ms.tot[1]);
        for (int k = 0; k < nFu; ++k) {
            Full& f = fu[k];
            f.sp = sp[4 * k]; FullFit(f);
            f.fold[0] = Eval(sp[4 * k + 1]); f.fold[1] = Eval(sp[4 * k + 2]); f.fail = Eval(sp[4 * k + 3]);
            f.fake = FakeS(smoothNone, ms.pas[0][k], ms.tot[0], ms.pas[1][k], ms.tot[1]);
            const double sa = fu[iNone].pm.eS, sp_ = f.pm.eS, sf = f.fail.eS;
            f.dF = f.pm.F - fu[iNone].pm.F;
            f.edF = k == iNone ? 0 : std::sqrt((1 - sp_ / sa) * (1 - sp_ / sa) + (sf / sa) * (sf / sa));
            const double tot = f.pm.S + f.fail.S;
            f.epsS = tot != 0 ? f.pm.S / tot : 1;
            f.eEpsS = tot != 0 ? std::sqrt(f.fail.S * f.fail.S * f.pm.eS * f.pm.eS + f.pm.S * f.pm.S * f.fail.eS * f.fail.eS) / (tot * tot) : 0;
            const double bp = f.pm.B + f.pm.CB, bf = f.fail.B + f.fail.CB;
            f.epsB = bp / (bp + bf);
        }
    }
    printf("full fits: %.1f s\n", sw.RealTime());

    OpenTable("peak_candidates.txt");
    P("CANDIDATES AND REFERENCE CHAINS, all rows (cor_ok == 1), full fits.\n");
    P("definition of the candidates (all rows): %d iterations, %s; grid of %zu allowed combinations, %d of them with F >= F_best - 1 (F_best = %.2f).\n",
      full.nIter, full.converged ? "converged" : "stopped on a cycle", full.grid.size(), full.nPlat, full.fbest);
    for (int k = 0; k < nFu; ++k) P("   %-18s %s\n", fu[k].name.c_str(), fu[k].text.c_str());
    P("\nS = whole Gaussian, P_main. spread = half of (max - min) of S over the kept model fits (orders 1-4, four ranges, width +-0.5 MeV; p > 0.05).\n");
    P("F_tot = S / sqrt(stat^2 + spread^2). S_pois = Poisson fit of the unsubtracted spectrum. S(++), S(--) = the same fit on the like-sign spectra.\n");
    P("S_fake = fake peak content from the pass fraction of mixed rows on a smooth background (see peak_guard.txt): a bias estimate, not subtracted.\n\n");
    P("%-18s | %7s %6s | %6s %5s %6s %5s %4s | %6s %6s %6s %8s | %6s %5s | %5s %4s | %5s %4s | %6s %4s\n", "selection", "N", "CB/N",
      "S", "stat", "spread", "rms", "kept", "S/B", "F_stat", "F_tot", "chi2/ndf", "S_pois", "err", "S(++)", "err", "S(--)", "err", "S_fake", "err");
    for (int k = 0; k < nFu; ++k) {
        const Full& f = fu[k];
        P("%-18s | %7.0f %6.3f | %6.0f %5.0f %6.0f %5.0f %4d | %6.4f %6.2f %6.2f %5.1f/%-2d | %6.0f %5.0f | %5.0f %4.0f | %5.0f %4.0f | %6.0f %4.0f\n",
          f.name.c_str(), f.pm.N, f.pm.CB / f.pm.N, f.pm.S, f.pm.eS, f.sd.half, f.sd.rms, f.sd.npass, f.pm.SoB, f.pm.F, f.ftot, f.pm.chi2, f.pm.ndf,
          f.pp.S, f.pp.eS, f.lsP.S, f.lsP.eS, f.lsM.S, f.lsM.eS, f.fake.s, f.fake.e);
        to.Fill("candidate", f.name.c_str(), 0, f.hasPt ? &f.pt : nullptr, Eval(f.sp), nullptr, &f);
    }
    P("S_fake of parts of the reference chains: eta_rest > 140 alone %.0f +- %.0f; oa_lab < 50 and M(pipi) < 0.42: %.0f +- %.0f;\n", diag[0].s, diag[0].e, diag[1].s, diag[1].e);
    P("   oa_lab < 50 and OA > 2: %.0f +- %.0f; oa_lab < 50 and eta_rest > 140: %.0f +- %.0f\n", diag[2].s, diag[2].e, diag[3].s, diag[3].e);
    P("\nmodel spread in parts: background model only (nominal width, kept of 16) and width +-0.5 MeV only; S range of the kept fits\n");
    for (int k = 0; k < nFu; ++k)
        P("%-18s | background %5.0f (kept %2d) | width %5.0f | S from %5.0f to %5.0f, mean %5.0f\n", fu[k].name.c_str(), fu[k].sd.halfBkg, fu[k].sd.npassBkg,
          fu[k].sd.halfWidth, fu[k].sd.smin, fu[k].sd.smax, fu[k].sd.mean);
    P("\nGAIN OVER 'no cut'. dF = F_stat - F_stat(no cut); its error from the independent rows removed by the selection:\n");
    P("   sigma(dF)^2 = (1 - sigma_pass/sigma_all)^2 + (sigma_fail/sigma_all)^2. eps_S = S/(S + S_fail), fail = all rows not in the selection.\n");
    P("   two halves: S and F on fold 0 and fold 1.\n");
    P("%-18s | %6s %5s %7s | %6s %5s %6s %6s | %6s %5s | %6s %5s %6s | %6s %5s %6s\n", "selection", "dF", "err", "dF/err", "eps_S", "err", "eps_B", "sqrt",
      "S_fail", "err", "S f0", "err", "F f0", "S f1", "err", "F f1");
    for (int k = 0; k < nFu; ++k) {
        const Full& f = fu[k];
        P("%-18s | %6.2f %5.2f %7.2f | %6.3f %5.3f %6.3f %6.3f | %6.0f %5.0f | %6.0f %5.0f %6.2f | %6.0f %5.0f %6.2f\n", f.name.c_str(), f.dF, f.edF,
          f.edF > 0 ? f.dF / f.edF : 0., f.epsS, f.eEpsS, f.epsB, std::sqrt(f.epsB), f.fail.S, f.fail.eS, f.fold[0].S, f.fold[0].eS, f.fold[0].F,
          f.fold[1].S, f.fold[1].eS, f.fold[1].F);
    }
    P("\nGAIN AGAINST BIAS. Measured optimisation bias on half samples (%d random halvings), in units of F: nominal %.2f +- %.2f, loose %.2f +- %.2f,\n",
      nSplit, bMean[0][0], bErr[0][0], bMean[1][0], bErr[1][0]);
    P("   pure %.2f +- %.2f, best-F %.2f +- %.2f (errors too small: the halvings reuse the data). The same procedures on all rows give:\n", bMean[2][0], bErr[2][0], bMean[3][0], bErr[3][0]);
    P("   nominal: dF = %.2f +- %.2f, minus bias %.2f = %.2f;  best-F grid point: dF = %.2f +- %.2f, minus bias %.2f = %.2f.\n",
      fu[iNom].dF, fu[iNom].edF, bMean[0][0], fu[iNom].dF - bMean[0][0], fu[8].dF, fu[8].edF, bMean[3][0], fu[8].dF - bMean[3][0]);
    P("   (the bias on all rows is at most that of the half samples: true differences in F grow by sqrt(2), their noise does not)\n");
    P("\nDISPLACED-PEAK SCAN: the same fit (frozen width, pol2, range and window moved with the peak) at 13 positions without a known peak.\n");
    P("pull = S / sigma_S. Expected for a correct background model: mean 0, rms 1.\n");
    P("%-18s |", "position");
    for (int i = 0; i < kNDisp; ++i) P(" %6.3f", 0.575 + 0.010 * i);
    P(" | %5s %5s %5s | %6s %5s\n", "mean", "rms", "max", "mean8", "rms8");
    for (int k = 0; k < nFu; ++k) {
        P("%-18s |", fu[k].name.c_str());
        for (int i = 0; i < kNDisp; ++i) P(" %6.2f", fu[k].dS[i] / fu[k].dE[i]);
        P(" | %5.2f %5.2f %5.2f | %6.2f %5.2f\n", fu[k].pullMean, fu[k].pullRms, fu[k].pullMax, fu[k].pullMean8, fu[k].pullRms8);
    }
    P("mean8, rms8 = the 8 positions 0.575-0.645, whose fit range ends below 0.72 (above it the spectrum bends: omega region).\n");
    P("The fit ranges of neighbouring positions overlap (87.5 MeV wide, 10 MeV steps): the pulls are strongly correlated.\n");
    P("A mean away from 0 (it is negative for the candidates and for 'no cut') means that pol2 does not describe the background at these positions;\n");
    P("the same model error can sit under the eta peak. It is covered by the model spread of S, not by the statistical error:\n");
    P("F_stat must not be read as a significance.\n");
    P("%-18s |", "S at 'nominal'");
    for (int i = 0; i < kNDisp; ++i) P(" %6.0f", fu[iNom].dS[i]);
    P("\n%-18s |", "chi2 at 'nominal'");
    for (int i = 0; i < kNDisp; ++i) P(" %6.1f", fu[iNom].dChi2[i]);
    P("\n");
    CloseTable();

    // ---- figures: fitted spectra, displaced peak, null test, summary ----
    DrawFits("peak_fits_candidates", Form("#eta peak, candidates: spectra and fit P_main (%s)", kHdr), {&fu[iLoose], &fu[iNom], &fu[iPure]}, true);
    DrawFits("peak_fits_references", Form("#eta peak, reference chains: spectra and fit P_main (%s)", kHdr), {&fu[0], &fu[3], &fu[4]}, true);
    {
        TCanvas* c = new TCanvas("c_disp", "", 2400, 1100);
        c->Divide(2, 1, 0.001, 0.001);
        const int set[2][3] = {{iLoose, iNom, iPure}, {0, 3, 4}};
        const int col[3] = {kAzure + 2, kBlue + 2, kViolet - 3}; const int mk[3] = {24, 20, 21};
        for (int p = 0; p < 2; ++p) {
            c->cd(p + 1); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.12);
            TH1D* fr = new TH1D(Form("dfr%d", p), "", 10, 0.565, 0.705); fr->SetDirectory(nullptr);
            fr->SetMinimum(-4.5); fr->SetMaximum(6.5);
            fr->GetXaxis()->SetTitle("assumed peak position in M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); fr->GetYaxis()->SetTitle("pull = S / #sigma_{S}");
            fr->GetYaxis()->SetTitleOffset(0.95); fr->GetXaxis()->SetNdivisions(507);
            fr->Draw("axis");
            Box(0.565, -1, 0.705, 1, kGreen - 9);
            Box(0.650, -4.5, 0.705, 3.6, kGray + 1, 3354);
            HLine(0.565, 0.705, 0, kBlack, 1, 1); HLine(0.565, 0.705, 2, kGray + 1, 2, 1); HLine(0.565, 0.705, -2, kGray + 1, 2, 1);
            TLegend* lg = Leg(0.14, 0.70, 0.95, 0.87, 0.026);
            for (int q = 0; q < 3; ++q) {
                const Full& f = fu[set[p][q]];
                TGraph* g = new TGraph();
                for (int i = 0; i < kNDisp; ++i) g->SetPoint(i, f.dPos[i] + 0.0012 * (q - 1), f.dS[i] / f.dE[i]);
                Deco(g, col[q], mk[q], 1.5);
                g->Draw("lp same");
                lg->AddEntry(g, Form("%s: mean %+.2f, rms %.2f (13 pos.); %+.2f, %.2f (first 8)", f.name.c_str(), f.pullMean, f.pullRms, f.pullMean8, f.pullRms8), "pl");
            }
            lg->Draw();
            Txt(0.12, 0.945, p == 0 ? "displaced-peak scan, candidates (peak 5.4 MeV + pol2)" : "displaced-peak scan, reference chains", 0.036);
            Txt(0.12, 0.895, "green: |pull| < 1; hatched: range ends above 0.72; no peak there: expect mean 0, rms 1", 0.027);
            gPad->RedrawAxis();
        }
        Save(c, "peak_displaced", false);
    }
    {
        const Full& f = fu[iNom];
        TCanvas* c = new TCanvas("c_null", "", 2400, 1100);
        c->Divide(2, 1, 0.001, 0.001);
        for (int k = 1; k <= 2; ++k) {
            c->cd(k); gPad->SetLeftMargin(0.16); gPad->SetTopMargin(0.12);
            TH1D* h = f.sp.Hist(k == 1 ? "pp" : "mm", Form("null_%d", k)); DecoH(h, kRed + 1, 20, 1.2);
            h->GetXaxis()->SetRangeUser(0.5001, 0.6299);
            h->SetMinimum(0.6 * h->GetBinContent(h->FindBin(0.523))); h->SetMaximum(1.40 * h->GetBinContent(h->FindBin(0.609)));
            h->GetYaxis()->SetTitleOffset(1.45);
            h->Draw("e");
            const scan::PeakFit& q = k == 1 ? f.lsP : f.lsM;
            Curve(q, 1, kBlack, 2)->Draw("l same"); Curve(q, 0, kBlack, 1)->Draw("l same");
            VLine(scan::kPeakLo, h->GetMinimum(), h->GetMaximum() * 0.8, kGray + 2, 2); VLine(scan::kPeakHi, h->GetMinimum(), h->GetMaximum() * 0.8, kGray + 2, 2);
            Txt(0.16, 0.945, Form("null test, nominal selection: like-sign pairs %s", k == 1 ? "e^{+}e^{+}" : "e^{#minus}e^{#minus}"), 0.042);
            Txt(0.16, 0.895, f.text.c_str(), 0.034);
            Txt(0.20, 0.82, Form("peak fit (Poisson): S = %.0f #pm %.0f, #chi^{2}/ndf = %.1f/%d", q.S, q.eS, q.chi2, q.ndf), 0.038);
            Txt(0.20, 0.765, "solid: Gaussian + pol2, dashed: pol2, vertical: window", 0.034);
        }
        Save(c, "peak_null_likesign", false);
    }
    {
        const int order[8] = {0, 1, 2, 3, 4, iLoose, iNom, iPure};
        TCanvas* c = new TCanvas("c_sum", "", 2000, 2300);
        TPad* top = new TPad("s_top", "", 0, 0.955, 1, 1); top->Draw();
        TPad* body = new TPad("s_body", "", 0, 0, 1, 0.955); body->Draw();
        top->cd(); Txt(0.01, 0.5, "#eta peak: reference chains and candidates (p+p 4.53 GeV, isBest #geq 0, corrected momenta, fit P_main)", 0.36, kBlack, 12);
        body->Divide(1, 3, 0.001, 0.001);
        const char* yt[3] = {"peak content S [rows]", "S / #sigma  (not a significance)", "S/B  and  CB/N"};
        for (int p = 0; p < 3; ++p) {
            body->cd(p + 1); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.16);
            TH1D* fr = new TH1D(Form("sfr%d", p), "", 8, 0, 8); fr->SetDirectory(nullptr);
            for (int k = 0; k < 8; ++k) fr->GetXaxis()->SetBinLabel(k + 1, fu[order[k]].name.c_str());
            fr->GetXaxis()->SetLabelSize(0.075); fr->GetYaxis()->SetTitle(yt[p]); fr->GetYaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetLabelSize(0.055);
            fr->GetYaxis()->SetTitleOffset(0.75);
            fr->SetMinimum(p == 0 ? 1500 : 0); fr->SetMaximum(p == 0 ? 7200 : (p == 1 ? 18 : 0.24));
            fr->Draw("axis");
            VLine(5, fr->GetMinimum(), fr->GetMaximum(), kGray + 1, 2, 1);
            TGraphErrors* g1 = new TGraphErrors(); TGraphErrors* g2 = new TGraphErrors();
            for (int k = 0; k < 8; ++k) {
                const Full& f = fu[order[k]];
                if (p == 0) {
                    g1->SetPoint(k, k + 0.5, f.pm.S); g1->SetPointError(k, 0, f.pm.eS);
                    g2->SetPoint(k, k + 0.5, f.pm.S); g2->SetPointError(k, 0, std::sqrt(f.pm.eS * f.pm.eS + f.sd.half * f.sd.half));
                } else if (p == 1) { g1->SetPoint(k, k + 0.5, f.pm.F); g2->SetPoint(k, k + 0.5, f.ftot); }
                else { g1->SetPoint(k, k + 0.5, f.pm.SoB); g2->SetPoint(k, k + 0.5, f.pm.CB / f.pm.N); }
            }
            Deco(g1, kColS, 20, 2.0); g1->SetLineWidth(5);
            Deco(g2, p == 2 ? kRed + 1 : kColS, p == 2 ? 21 : 24, 2.0); g2->SetLineWidth(1);
            if (p == 0) { g2->Draw("p same"); g1->Draw("p same"); } else { g2->Draw("p same"); g1->Draw("p same"); }
            TLegend* lg = Leg(0.12, 0.79, 0.75, 0.91, 0.055, 2);
            lg->AddEntry(g1, p == 0 ? "S #pm stat (thick)" : (p == 1 ? "F_{stat} = S/#sigma_{S}: statistical only, frozen shape" : "S/B"), p == 0 ? "pe" : "p");
            lg->AddEntry(g2, p == 0 ? "#pm stat #oplus model spread (thin)" : (p == 1 ? "F_{tot} = S/#sqrt{stat^{2} + model^{2}}" : "CB/N in the window (no k factor: an estimate)"), p == 0 ? "pe" : "p");
            lg->Draw();
            gPad->RedrawAxis();
        }
        Save(c, "peak_summary", true);
    }

    // ------------------------------------------------------------------------
    // 6. content of the peak by class (nominal candidate; width free per class)
    // ------------------------------------------------------------------------
    struct Cls { const char* group; const char* name; const char* tex; std::function<bool(const scan::Row&)> f; };
    const std::vector<Cls> cls = {
        {"OA", "OA < 2", "< 2", [](const scan::Row& r) { return r.oa < 2.f; }},
        {"OA", "OA 2-6", "2#minus6", [](const scan::Row& r) { return r.oa >= 2.f && r.oa < 6.f; }},
        {"OA", "OA 6-15", "6#minus15", [](const scan::Row& r) { return r.oa >= 6.f && r.oa < 15.f; }},
        {"OA", "OA > 15", "> 15", [](const scan::Row& r) { return r.oa >= 15.f; }},
        {"M(ee)", "M(ee) < 8 MeV", "< 8", [](const scan::Row& r) { return r.mee < 0.008f; }},
        {"M(ee)", "M(ee) 8-20", "8#minus20", [](const scan::Row& r) { return r.mee >= 0.008f && r.mee < 0.020f; }},
        {"M(ee)", "M(ee) 20-50", "20#minus50", [](const scan::Row& r) { return r.mee >= 0.020f && r.mee < 0.050f; }},
        {"M(ee)", "M(ee) 50-135", "50#minus135", [](const scan::Row& r) { return r.mee >= 0.050f && r.mee < 0.135f; }},
        {"M(ee)", "M(ee) > 135", "> 135", [](const scan::Row& r) { return r.mee >= 0.135f; }},
        {"ring", "shared RICH ring", "shared", [](const scan::Row& r) { return r.shared_ring == 1; }},
        {"ring", "two RICH rings", "two rings", [](const scan::Row& r) { return r.shared_ring != 1; }},
        {"rows", "event with 1 row", "1 row", [](const scan::Row& r) { return r.nrow == 1; }},
        {"rows", "event with > 1 row", "> 1 row", [](const scan::Row& r) { return r.nrow > 1; }},
        {"isBest", "isBest == 1", "= 1", [](const scan::Row& r) { return r.isBest == 1; }},
        {"isBest", "isBest == 0", "= 0", [](const scan::Row& r) { return r.isBest != 1; }}};
    const int nCls = (int)cls.size();
    struct ClsRes { Ev ev; WFit wf; };
    std::vector<ClsRes> cres[2]; ClsRes ctot[2];
    {
        const scan::Sel base[2] = {fu[iNom].sel, fu[iNone].sel};
        std::vector<scan::Sel> sels;
        for (int b = 0; b < 2; ++b) {
            sels.push_back(base[b]);
            for (int k = 0; k < nCls; ++k) { const scan::Sel s = base[b]; const auto g = cls[k].f; sels.push_back([s, g](const scan::Row& r) { return s(r) && g(r); }); }
        }
        const std::vector<scan::Spectra> sp = scan::FillMany(d, sels);
        OpenTable("peak_classes.txt");
        P("CONTENT OF THE PEAK BY CLASS. Each class is fitted separately (P_main, pol2, frozen position %.4f).\n", scan::kPeakMean);
        P("S(5.4) = width frozen at 5.4 MeV; S(free) = width free (chi2 scanned from 2.5 to 11 MeV; error of S includes the width); n/a = no minimum.\n");
        P("share = S(5.4) of the class / sum of S(5.4) over the classes of the group. dchi2 = chi2(5.4) - chi2(free width).\n");
        P("The classes of one group are disjoint and add up to the selection ('sum' against 'all' checks the additivity of the fits).\n");
        P("N++, N-- = like-sign rows in the peak window; CB/N = 2 sqrt(N++ N--)/N without a k factor. Where N++ is a handful of rows (small OA, small M(ee),\n");
        P("shared ring) the CB share is only a lower limit. arith = (N++ + N--)/N, the share if the like-sign sum were used instead.\n");
        P("READING: the classes are SIGNS of what the peak holds (eta -> pi+pi-e+e- and eta -> pi+pi-gamma with a converted photon both give small OA and\n");
        P("small M(ee)); they cannot be turned into decay-mode fractions without a simulation.\n");
        for (int b = 0; b < 2; ++b) {
            P("\nSELECTION: %s (%s)\n", fu[b == 0 ? iNom : iNone].name.c_str(), fu[b == 0 ? iNom : iNone].text.c_str());
            P("%-7s %-20s | %7s %6s %6s | %6s %5s %6s %6s %8s | %6s %5s %6s %5s %6s | %6s %6s %7s %6s\n", "group", "class", "N", "CB/N", "S/B", "S(5.4)", "err", "share", "err",
              "chi2/ndf", "S(free)", "err", "width", "err", "dchi2", "N++", "N--", "N--/N++", "arith");
            ctot[b].ev = Eval(sp[b * (nCls + 1)]); ctot[b].wf = WidthFit(sp[b * (nCls + 1)]);
            cres[b].resize(nCls);
            for (int k = 0; k < nCls; ++k) { cres[b][k].ev = Eval(sp[b * (nCls + 1) + 1 + k]); cres[b][k].wf = WidthFit(sp[b * (nCls + 1) + 1 + k]); }
            auto line = [&](const char* g, const char* n, const ClsRes& r, double share, double eshare) {
                P("%-7s %-20s | %7.0f %6.3f %6.3f | %6.0f %5.0f %6.3f %6.3f %5.1f/%-2d | ", g, n, r.ev.N, r.ev.cbs, r.ev.SoB, r.ev.S, r.ev.eS, share, eshare, r.ev.chi2, r.ev.ndf);
                if (r.wf.ok) P("%6.0f %5.0f %6.2f %5.2f %6.1f", r.wf.S, r.wf.eS, 1e3 * r.wf.sigma, 1e3 * r.wf.eSigma, r.wf.dchi2);
                else P("%6s %5s %6s %5s %6s", "n/a", "", "n/a", "", "");
                P(" | %6.0f %6.0f %7.1f %6.3f\n", r.ev.npp, r.ev.nmm, r.ev.npp > 0 ? r.ev.nmm / r.ev.npp : 0., r.ev.N > 0 ? (r.ev.npp + r.ev.nmm) / r.ev.N : 0.);
            };
            line("all", "whole selection", ctot[b], 1., 0.);
            for (int k = 0; k < nCls; ++k) {
                double sum = 0, sum2 = 0, sumF = 0, sumF2 = 0; bool allOk = true;
                for (int j = 0; j < nCls; ++j) if (std::string(cls[j].group) == cls[k].group) {
                    sum += cres[b][j].ev.S; sum2 += cres[b][j].ev.eS * cres[b][j].ev.eS;
                    if (cres[b][j].wf.ok) { sumF += cres[b][j].wf.S; sumF2 += cres[b][j].wf.eS * cres[b][j].wf.eS; } else allOk = false;
                }
                const double s = cres[b][k].ev.S, e = cres[b][k].ev.eS, rest = sum - s, erest2 = sum2 - e * e;
                const double eshare = sum != 0 ? std::sqrt(rest * rest * e * e + s * s * erest2) / (sum * sum) : 0;
                line(cls[k].group, cls[k].name, cres[b][k], sum != 0 ? s / sum : 0, eshare);
                if (k == nCls - 1 || std::string(cls[k + 1].group) != cls[k].group) {
                    P("%-7s %-20s | %7s %6s %6s | %6.0f %5.0f %6s %6s %8s | ", cls[k].group, "sum of the classes", "", "", "", sum, std::sqrt(sum2), "", "", "");
                    if (allOk) P("%6.0f %5.0f\n", sumF, std::sqrt(sumF2)); else P("%6s\n", "n/a");
                    if (std::string(cls[k].group) == "M(ee)") {
                        // cumulative shares of the content above 20 and 50 MeV (classes in rising order of M(ee))
                        int first = -1;
                        for (int j = 0; j < nCls; ++j) if (first < 0 && std::string(cls[j].group) == "M(ee)") first = j;
                        for (int cutAt = 2; cutAt <= 3; ++cutAt) {
                            double a = 0, a2 = 0, b0 = 0, b2 = 0;
                            for (int j = first; j <= k; ++j) {
                                if (j - first >= cutAt) { a += cres[b][j].ev.S; a2 += cres[b][j].ev.eS * cres[b][j].ev.eS; }
                                else { b0 += cres[b][j].ev.S; b2 += cres[b][j].ev.eS * cres[b][j].ev.eS; }
                            }
                            const double t = a + b0;
                            P("%-7s %-20s | share of the content %.3f +- %.3f (S = %.0f +- %.0f)\n", "M(ee)", cutAt == 2 ? "M(ee) > 20 MeV" : "M(ee) > 50 MeV",
                              t != 0 ? a / t : 0., t != 0 ? std::sqrt(b0 * b0 * a2 + a * a * b2) / (t * t) : 0., a, std::sqrt(a2));
                        }
                    }
                }
            }
        }
        CloseTable();
    }
    {
        const char* grp[5] = {"OA", "M(ee)", "ring", "rows", "isBest"};
        const char* gt[5] = {"OA(e^{+}e^{#minus}) [deg]", "M(e^{+}e^{#minus}) [MeV/c^{2}]", "RICH rings of the pair", "rows of the event", "isBest of the row"};
        TCanvas* c = new TCanvas("c_cls", "", 3200, 1700);
        TPad* top = new TPad("c_top", "", 0, 0.93, 1, 1); top->Draw();
        TPad* body = new TPad("c_body", "", 0, 0, 1, 0.93); body->Draw();
        top->cd();
        Txt(0.01, 0.55, Form("#eta peak content by class, nominal selection (%s); each class fitted separately, P_main, pol2", fu[iNom].text.c_str()), 0.40, kBlack, 12);
        Txt(0.01, 0.12, "upper row: peak content S; filled: width frozen at 5.4 MeV, open: width free.  Lower row: fitted width; line: 5.4 MeV.  Numbers: S/B and CB share (no k factor; a lower limit at small OA).  Classes are signs, not decay-mode fractions", 0.30, kBlack, 12);
        body->Divide(5, 2, 0.001, 0.001);
        for (int g = 0; g < 5; ++g) {
            std::vector<int> id;
            for (int k = 0; k < nCls; ++k) if (std::string(cls[k].group) == grp[g]) id.push_back(k);
            const int n = (int)id.size();
            double ym = 0;
            for (int k : id) ym = std::max(ym, std::max(cres[0][k].ev.S + cres[0][k].ev.eS, cres[0][k].wf.ok ? cres[0][k].wf.S + cres[0][k].wf.eS : 0.));
            for (int row = 0; row < 2; ++row) {
                body->cd(row * 5 + g + 1);
                gPad->SetLeftMargin(0.22); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.20);
                TH1D* fr = new TH1D(Form("cfr%d_%d", g, row), "", n, 0, n); fr->SetDirectory(nullptr);
                for (int q = 0; q < n; ++q) fr->GetXaxis()->SetBinLabel(q + 1, cls[id[q]].tex);
                fr->GetXaxis()->SetLabelSize(n > 4 ? 0.058 : 0.08); fr->GetXaxis()->SetTitle(gt[g]); fr->GetXaxis()->SetTitleSize(0.06); fr->GetXaxis()->SetTitleOffset(1.4);
                fr->GetYaxis()->SetTitle(row == 0 ? "peak content S [rows]" : "peak width #sigma [MeV/c^{2}]"); fr->GetYaxis()->SetTitleSize(0.06);
                fr->GetYaxis()->SetLabelSize(0.055); fr->GetYaxis()->SetTitleOffset(1.85);
                fr->SetMinimum(row == 0 ? -300 : 0); fr->SetMaximum(row == 0 ? 1.45 * ym : 13);
                fr->Draw("axis");
                HLine(0, n, row == 0 ? 0 : 5.4, kBlack, row == 0 ? 3 : 2, 1);
                TGraphErrors* g1 = new TGraphErrors(); TGraphErrors* g2 = new TGraphErrors();
                for (int q = 0; q < n; ++q) {
                    const ClsRes& r = cres[0][id[q]];
                    if (row == 0) {
                        g1->SetPoint(q, q + 0.38, r.ev.S); g1->SetPointError(q, 0, r.ev.eS);
                        if (r.wf.ok) { const int np = g2->GetN(); g2->SetPoint(np, q + 0.62, r.wf.S); g2->SetPointError(np, 0, r.wf.eS); }
                        TLatex l; l.SetTextFont(42); l.SetTextSize(n > 4 ? 0.036 : 0.042); l.SetTextAlign(21);
                        l.DrawLatex(q + 0.5, 1.34 * ym, Form("%.3f", r.ev.SoB));
                        l.DrawLatex(q + 0.5, 1.24 * ym, Form("%.0f%%", 100 * r.ev.cbs));
                        if (q == 0) { l.SetTextAlign(11); l.SetTextSize(0.036); l.DrawLatex(0.03, 1.14 * ym, "above: S/B, CB share (no k)"); }
                    } else if (r.wf.ok) { const int np = g2->GetN(); g2->SetPoint(np, q + 0.5, 1e3 * r.wf.sigma); g2->SetPointError(np, 0, 1e3 * r.wf.eSigma); }
                }
                Deco(g1, kColS, 20, 1.8); Deco(g2, kColS, 24, 1.8);
                if (row == 0) g1->Draw("p same");
                g2->Draw("p same");
                gPad->RedrawAxis();
            }
        }
        Save(c, "peak_classes", true);
    }

    // ------------------------------------------------------------------------
    // 7. cor_ok fallback variant at the nominal candidate
    // ------------------------------------------------------------------------
    {
        const scan::Sel sn = fu[iNom].sel;
        d.SetFallback(true);
        const std::vector<scan::Spectra> sp = scan::FillMany(d, {sn, [sn](const scan::Row& r) { return !r.cor_ok && sn(r); }, fu[iNone].sel,
                                                                 [](const scan::Row& r) { return !r.cor_ok; }});
        d.SetFallback(false);
        OpenTable("peak_fallback.txt");
        P("COR_OK FALLBACK VARIANT: rows with cor_ok == 0 (no corrected pi- momentum in PAT) enter with the reconstructed momentum of that track.\n");
        P("added rows = the cor_ok == 0 rows alone (Poisson fit, independent of the main sample).\n\n");
        P("%-34s | %7s %6s | %6s %5s %6s %6s %8s | %6s %5s\n", "selection", "N", "CB/N", "S", "err", "F", "S/B", "chi2/ndf", "S_pois", "err");
        for (int b = 0; b < 2; ++b) {
            const Full& f = fu[b == 0 ? iNom : iNone];
            const scan::PeakFit a = scan::FitPeakMain(sp[2 * b]), ap = scan::FitPeakPois(sp[2 * b]), x = scan::FitPeakPois(sp[2 * b + 1]);
            P("%-34s | %7.0f %6.3f | %6.0f %5.0f %6.2f %6.4f %5.1f/%-2d | %6.0f %5.0f\n", Form("%s, cor_ok == 1 (main)", f.name.c_str()), f.pm.N, f.pm.CB / f.pm.N,
              f.pm.S, f.pm.eS, f.pm.F, f.pm.SoB, f.pm.chi2, f.pm.ndf, f.pp.S, f.pp.eS);
            P("%-34s | %7.0f %6.3f | %6.0f %5.0f %6.2f %6.4f %5.1f/%-2d | %6.0f %5.0f\n", Form("%s, with cor_ok == 0 rows", f.name.c_str()), a.N, a.CB / a.N,
              a.S, a.eS, a.F, a.SoB, a.chi2, a.ndf, ap.S, ap.eS);
            P("%-34s | %7.0f %6s | %6s %5s %6s %6s %8s | %6.0f %5.0f\n", Form("%s, added rows alone", f.name.c_str()), x.N, "", "", "", "", "", "", x.S, x.eS);
            P("   change of S: %+.0f rows (%.1f %%), of sigma_S: %+.0f rows, of F: %+.2f; rows in the window: +%.1f %%\n\n", a.S - f.pm.S, 100 * (a.S / f.pm.S - 1),
              a.eS - f.pm.eS, a.F - f.pm.F, 100 * (a.N / f.pm.N - 1));
            to.Fill("fallback", f.name.c_str(), 0, nullptr, Eval(sp[2 * b]));
        }
        CloseTable();
    }

    // ---- spectra of the candidates into the output file ----
    to.f->cd();
    for (int k = 0; k < nFu; ++k) {
        std::string tag = fu[k].name;
        for (auto& ch : tag) if (!isalnum(ch)) ch = '_';
        for (const char* w : {"us", "pp", "mm", "cb", "sig"}) {
            TH1D* h = fu[k].sp.Hist(w, Form("%s_%s", w, tag.c_str()));
            h->SetTitle(Form("%s, %s: %s;%s;%s", w, fu[k].name.c_str(), fu[k].text.c_str(), kXt, kYt));
            h->Write(); delete h;
        }
    }
    to.Close();
    printf("total run time %.1f s\n", swAll.RealTime());
}
