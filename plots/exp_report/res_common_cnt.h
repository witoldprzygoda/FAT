// res_common_cnt.h -- COPY of plots/exp_report/res_common.h made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// res_common.h  --  shared code of the chapter "Results without simulation" (prefix res)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Header used by plots/exp_report/res_peak.C, res_bump.C and res_mee.C (it is not
// a macro of its own). It holds
//   A. paths, luminosity and trigger constants
//   B. the selections of the chapter (defined exactly as in the task text)
//   C. output helpers (text tables, figure style, vector PDF + PNG)
//   D. the M4 side-band subtraction of N - CB (used for the M(ee) / OA content
//      of the eta peak and for the X17 window)
//   E. the ECAL-tagged eta (copy of the frozen definition of plots/scan_bump.C)
//   F. luminosity per day / period from the per-run table
//
// Conventions of all res_* macros
//   - rows = stored particle combinations; main variant all rows (isBest >= 0),
//     comparison isBest == 1; corrected momenta, rows with cor_ok == 1 only
//     (the engine skips cor_ok == 0 rows; no fallback)
//   - in the FAT production already: event vertex z > -500 mm, PT3 (trigbit == 8192)
//   - CB = 2 sqrt(N++ N--) per 10 MeV M4 bin, interpolated to 2.5 MeV (ScanEngine), k = 1
//   - data windows: peak 0.53364-0.55416, bump 0.4550-0.5225 GeV/c2
//
// Compile the macros with ACLiC and a build directory outside the repository:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_res", kTRUE)' 'plots/exp_report/res_peak.C+'
// ============================================================================
#ifndef RES_COMMON_CNT_H
#define RES_COMMON_CNT_H

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../ScanEngine.h"

namespace res {

using scan::Data;
using scan::PeakCfg;
using scan::PeakFit;
using scan::Row;
using scan::Sel;
using scan::Spectra;

// ============================================================================
// A. paths and constants
// ============================================================================
const char* const kFigOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp";
const char* const kTabOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res/tab";
const char* const kWork   = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res/work";
const char* const kCntScratch = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res/figs";
const char* const kCntScope = " res_bump_slices res_bump_tagged res_peak_fits res_peak_overview res_tagged_period res_x17_window fix_spread ";
// per-run luminosity table with the coverage of the FAT input (copy of
// scratchpad final/remote/lumi_coverage_perrun.csv; columns dst_stem, runid, day,
// n_dst_list, n_dst_merged, L6_SF, L5_SF, no_hf)
const char* const kLumiCsv = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res/lumi_coverage_perrun.csv";

// START-free luminosity (HADES elastic-scattering analysis, LGADfree_luminosity_talk_draft_v2)
const double kLall   = 6080.5;     // nb^-1, all 17 615 runs, six sectors (recommended)
const double kLdata  = 6079.2;     // nb^-1, runs present in our data
const double kLrelSyst  = 0.0556;  // HADES systematic (+-338 nb^-1)
const double kLrelFit   = 0.0103;  // Argonne fit
const double kLrelNorm  = 0.040;   // normalisation
const double kLrelTot   = 0.069;   // total
const double kLlive     = 0.932;   // PT3 live fraction relative to PT2: one-sided, L x 0.932 at most (-6.8 %)
const double kPT3       = 64.0 / 63.0;   // trigger word 8192 holds 63/64 of the PT3 events

// data side bands of the peak window (same width as in plots/scan_extra.C, 11.14 MeV each)
const double kSb1Lo = 0.5225, kSb1Hi = scan::kPeakLo;   // 0.5225 - 0.53364
const double kSb2Lo = scan::kPeakHi, kSb2Hi = 0.56530;  // 0.55416 - 0.56530

// periods of the beam time (file_id = day)
const int kNPer = 4;
const char* const kPerName[kNPer] = {"I", "IIA", "IIB", "III"};
inline int Period(int day) {
    if (day >= 34 && day <= 46) return 0;
    if (day >= 47 && day <= 48) return 1;
    if (day >= 49 && day <= 53) return 2;
    if (day >= 59 && day <= 67) return 3;
    return -1;
}

// ============================================================================
// B. selections (cor_ok == 1 is applied by the engine; vz = eVertReco_z)
// ============================================================================
inline bool Pre(const Row& r) { return r.vz > -160.f; }                         // + isBest >= 0 (all rows)
inline bool PiVtx(const Row& r) { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; }
inline bool PeakNomKin(const Row& r) { return r.oa_lab < 55.f && r.mpipi > 0.30f; }
inline bool BumpNomKin(const Row& r) { return r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; }
inline bool X17Mee(const Row& r) { return r.mee >= 0.0155f && r.mee < 0.0245f; }

struct NSel {
    std::string key;    // short name for tables
    std::string tex;    // LaTeX name
    std::string def;    // definition in words
    Sel f;
};
inline NSel SelPre() { return {"pre", "pre", "cor_ok && isBest >= 0 && vz > -160 mm", [](const Row& r) { return Pre(r); }}; }
inline NSel SelPreVtx() { return {"pre+vertex", "pre + vertex", "pre + both pions in the vertex fit",
                                  [](const Row& r) { return Pre(r) && PiVtx(r); }}; }
inline NSel SelPeakNom() { return {"peak nominal", "peak nominal", "pre + oa_lab < 55 + M(pipi) > 0.30",
                                   [](const Row& r) { return Pre(r) && PeakNomKin(r); }}; }
inline NSel SelPeakNomVtx() { return {"peak nominal+vertex", "peak nominal + vertex",
                                      "pre + oa_lab < 55 + M(pipi) > 0.30 + both pions in the vertex fit",
                                      [](const Row& r) { return Pre(r) && PeakNomKin(r) && PiVtx(r); }}; }
inline NSel SelBumpNom() { return {"bump nominal", "bump nominal", "pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41",
                                   [](const Row& r) { return Pre(r) && BumpNomKin(r); }}; }
inline NSel SelBumpNomVtx() { return {"bump nominal+vertex", "bump nominal + vertex",
                                      "pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41 + both pions in the vertex fit",
                                      [](const Row& r) { return Pre(r) && BumpNomKin(r) && PiVtx(r); }}; }
// frozen selections of the simulation study (FAT_simwt plots/eta_sim/tables/s5_opt_eta_frozen.txt)
inline NSel SelPEAK() { return {"PEAK", "PEAK", "pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2",
                                [](const Row& r) { return Pre(r) && r.dz_pi < 40.f && r.mm4 < 2.85f; }}; }
inline NSel SelDIRECT() { return {"DIRECT", "DIRECT", "pre + dz_pi < 40 mm",
                                  [](const Row& r) { return Pre(r) && r.dz_pi < 40.f; }}; }
// X17 selection without its M4 window (the M4 spectrum is fitted / the window applied later)
inline NSel SelX17() { return {"X17 region", "X17 region", "pre + 15.5 <= M(ee) < 24.5 MeV/c2 + both pions in the vertex fit",
                               [](const Row& r) { return Pre(r) && X17Mee(r) && PiVtx(r); }}; }

// variants
inline Sel Best1(const Sel& s) { return [s](const Row& r) { return r.isBest == 1 && s(r); }; }
inline Sel OA2(const Sel& s) { return [s](const Row& r) { return r.oa > 2.f && s(r); }; }
inline Sel And(const Sel& a, const Sel& b) { return [a, b](const Row& r) { return a(r) && b(r); }; }
// variant index v: 0 = isBest >= 0, 1 = isBest == 1, 2 = isBest >= 0 + OA > 2, 3 = isBest == 1 + OA > 2
inline Sel Variant(const Sel& s, int v) {
    Sel o = s;
    if (v & 2) o = OA2(o);
    if (v & 1) o = Best1(o);
    return o;
}
inline const char* VariantName(int v) {
    static const char* n[4] = {"isBest>=0", "isBest==1", "isBest>=0, OA>2", "isBest==1, OA>2"};
    return n[v & 3];
}
inline const char* VariantTex(int v) {
    static const char* n[4] = {"$\\geq 0$", "$=1$", "$\\geq 0$, OA$>2^\\circ$", "$=1$, OA$>2^\\circ$"};
    return n[v & 3];
}

// ============================================================================
// C. output helpers
// ============================================================================
FILE* gOut = nullptr;
inline void P(const char* fmt, ...) {
    va_list a, b;
    va_start(a, fmt);
    va_copy(b, a);
    vprintf(fmt, a);
    if (gOut) vfprintf(gOut, fmt, b);
    va_end(a);
    va_end(b);
}
inline void OpenTxt(const char* name) {
    if (gOut) fclose(gOut);
    gOut = fopen(Form("%s/%s", kTabOut, name), "w");
}
inline void CloseTxt() { if (gOut) fclose(gOut); gOut = nullptr; }

// LaTeX snippet: plain text written as given
struct Tex {
    FILE* f = nullptr;
    explicit Tex(const char* name) { f = fopen(Form("%s/%s", kTabOut, name), "w"); }
    ~Tex() { if (f) fclose(f); }
    void operator()(const char* fmt, ...) {
        va_list a;
        va_start(a, fmt);
        if (f) vfprintf(f, fmt, a);
        va_end(a);
    }
};

inline void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelSize(0.050, "xyz"); gStyle->SetTitleSize(0.056, "xyz");
    gStyle->SetTitleOffset(1.05, "x"); gStyle->SetTitleOffset(1.30, "y");
    gStyle->SetPadLeftMargin(0.16); gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.09); gStyle->SetPadBottomMargin(0.14);
    gStyle->SetEndErrorSize(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
}
// vector PDF cropped to the canvas (EPS -> ps2pdf -dEPSCrop, as in the other chapters) + PNG preview
inline void Save(TCanvas* c, const char* name) {
    const std::string n = name;
    const bool inScope = std::string(kCntScope).find(" " + n + " ") != std::string::npos;
    const std::string eb = std::string(kCntScratch) + "/" + n;                      // EPS always in the scratch folder
    const std::string b = inScope ? std::string(kFigOut) + "/" + n + "_cnt" : eb;   // new name next to the original
    c->SaveAs((eb + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf", eb.c_str(), b.c_str()));
    c->SaveAs((b + ".png").c_str());
}

inline void Filter(const Data& in, Data& out, const std::function<bool(const Row&)>& keep) {
    for (int k = 0; k < 3; ++k) {
        out.s[k].clear();
        for (const auto& r : in.s[k]) if (keep(r)) out.s[k].push_back(r);
        out.s[k].shrink_to_fit();
    }
    out.fallback = in.fallback;
}

// copy of a spectrum with the CB multiplied by k (sensitivity to the CB normalisation)
inline Spectra ScaleCB(const Spectra& sp, double k) {
    Spectra o = sp;
    for (auto& v : o.cbW) v *= k;
    for (auto& v : o.cbWVar) v *= k * k;
    for (int i = 0; i < o.nbin; ++i) {
        o.cb[i] = k * sp.cb[i];
        o.cbVar[i] = k * k * sp.cbVar[i];
        o.sig[i] = o.n[0][i] - o.cb[i];
        o.sigVar[i] = o.n[0][i] + o.cbVar[i];
    }
    for (auto& w : o.win) { w.cb *= k; w.cbVar *= k * k; }
    return o;
}

// visible yield per luminosity [nb] = S x 64/63 / L
inline double Yield(double s) { return s * kPT3 / kLdata; }

// ============================================================================
// D. M4 side-band subtraction of N - CB
// ============================================================================
// For a selection (e.g. the eta-peak base selection AND one M(ee) bin):
//   D = (N - CB)_W - r (N - CB)_SB
//   W  = data peak window 0.53364-0.55416, SB = 0.5225-0.53364 and 0.55416-0.56530
//   r  = integral of the fitted P_main polynomial (base selection, all M(ee)) over W
//        divided by its integral over SB: the non-CB background is assumed to have
//        the same M4 shape in every bin of M(ee) or OA
//   S  = D / (fW - r fSB): whole-Gaussian content (fW, fSB = shares of the frozen
//        Gaussian 0.5439 / 5.4 MeV inside W and inside SB; the tails in the side
//        bands are subtracted with the background)
// Errors: unlike-sign rows per event group with the weights +1 (W) and -r (SB),
// so rows of one group in W and SB are correlated correctly; CB = 2 sqrt(N++ N--)
// in 10 MeV bins, interpolated to 2.5 MeV and integrated over W and SB exactly as
// in ScanEngine; its variance per like-sign event group with dCB/dN++ =
// sqrt(N--/N++) (and vice versa).
struct SBRes {
    double r = 0, fW = 0, fSB = 0;
    double nW = 0, nSB = 0;              // unlike-sign rows
    double cbW = 0, cbSB = 0;            // CB integrated over W and SB
    double ppW = 0, mmW = 0;             // like-sign rows in W
    double DW = 0, eDW = 0;              // (N - CB)_W
    double DSB = 0, eDSB = 0;            // (N - CB)_SB
    double D = 0, eD = 0;                // (N - CB)_W - r (N - CB)_SB
    double S = 0, eS = 0;                // D / (fW - r fSB)
    double bkgW() const { return r * DSB; }   // non-CB background in W
};
// integral of a fitted P_main polynomial over [a, b) (rows / fine bin, midpoint rule)
inline double PolyInt(const PeakFit& f, double a, double b) {
    const int n = 400; const double h = (b - a) / n; double s = 0;
    for (int q = 0; q < n; ++q) s += f.Bkg(a + (q + 0.5) * h);
    return s * h;
}
inline double RatioWSB(const PeakFit& f) {
    const double w = PolyInt(f, scan::kPeakLo, scan::kPeakHi);
    const double b = PolyInt(f, kSb1Lo, kSb1Hi) + PolyInt(f, kSb2Lo, kSb2Hi);
    return b > 0 ? w / b : 0;
}
// several selections in one pass over the rows
inline std::vector<SBRes> SideBandMany(const Data& d, const std::vector<Sel>& sels, double r,
                                       double mean = scan::kPeakMean, double sigma = scan::kPeakSigma) {
    const int ns = (int)sels.size();
    std::vector<SBRes> out(ns);
    // geometry of the CB interpolation: dummy spectrum with the engine binning
    Spectra g; g.nbin = scan::kNbin; g.lo = scan::kMlo; g.hi = scan::kMhi; g.bw = scan::kBinW;
    g.reb = scan::kCBreb; g.nw = scan::kNbin / scan::kCBreb;
    const int nw = g.nw;
    std::vector<double> cW(nw, 0.), cSB(nw, 0.);
    for (int i = 0; i < g.nbin; ++i) {
        const double uW = g.Frac(i, scan::kPeakLo, scan::kPeakHi);
        const double uS = g.Frac(i, kSb1Lo, kSb1Hi) + g.Frac(i, kSb2Lo, kSb2Hi);
        if (uW <= 0 && uS <= 0) continue;
        int j0, j1; double f; g.Interp(i, j0, j1, f);
        cW[j0] += uW * (1 - f) / g.reb; cW[j1] += uW * f / g.reb;
        cSB[j0] += uS * (1 - f) / g.reb; cSB[j1] += uS * f / g.reb;
    }
    const double fW = scan::GaussFrac(scan::kPeakLo, scan::kPeakHi, mean, sigma);
    const double fSB = scan::GaussFrac(kSb1Lo, kSb1Hi, mean, sigma) + scan::GaussFrac(kSb2Lo, kSb2Hi, mean, sigma);
    // unlike-sign rows
    std::vector<double> aD(ns, 0.), aW(ns, 0.), aS(ns, 0.), vD(ns, 0.), vW(ns, 0.), vS(ns, 0.);
    // like-sign rows: wide-bin counts and (grp, bin) lists
    std::vector<std::vector<double>> wl[2];
    std::vector<std::vector<std::pair<int, int>>> lsr[2];
    for (int q = 0; q < 2; ++q) { wl[q].assign(ns, std::vector<double>(nw, 0.)); lsr[q].assign(ns, {}); }
    for (int k = 0; k < 3; ++k) {
        const std::vector<Row>& v = d.s[k];
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const Row& row = v[i];
            if (!row.cor_ok && !d.fallback) continue;
            if (k == 0 && row.grp != cur) {
                for (int s = 0; s < ns; ++s) {
                    vD[s] += aD[s] * aD[s]; vW[s] += aW[s] * aW[s]; vS[s] += aS[s] * aS[s];
                    aD[s] = aW[s] = aS[s] = 0;
                }
                cur = row.grp;
            }
            const double m = row.m4;
            if (!(m >= 0.50 && m < 0.58)) continue;
            const bool inW = m >= scan::kPeakLo && m < scan::kPeakHi;
            const bool inS = (m >= kSb1Lo && m < kSb1Hi) || (m >= kSb2Lo && m < kSb2Hi);
            const int jb = (int)((m - scan::kMlo) / (g.reb * g.bw));
            for (int s = 0; s < ns; ++s) {
                if (!sels[s](row)) continue;
                if (k == 0) {
                    if (inW) { out[s].nW += 1; aD[s] += 1; aW[s] += 1; }
                    else if (inS) { out[s].nSB += 1; aD[s] -= r; aS[s] += 1; }
                } else {
                    if (jb >= 0 && jb < nw) { wl[k - 1][s][jb] += 1; lsr[k - 1][s].push_back(std::make_pair(row.grp, jb)); }
                    if (inW) { if (k == 1) out[s].ppW += 1; else out[s].mmW += 1; }
                }
            }
        }
        if (k == 0)
            for (int s = 0; s < ns; ++s) {
                vD[s] += aD[s] * aD[s]; vW[s] += aW[s] * aW[s]; vS[s] += aS[s] * aS[s];
                aD[s] = aW[s] = aS[s] = 0;
            }
    }
    for (int s = 0; s < ns; ++s) {
        SBRes& o = out[s];
        o.r = r; o.fW = fW; o.fSB = fSB;
        std::vector<double> cb(nw, 0.);
        for (int j = 0; j < nw; ++j) cb[j] = 2.0 * std::sqrt(wl[0][s][j] * wl[1][s][j]);
        for (int j = 0; j < nw; ++j) { o.cbW += cW[j] * cb[j]; o.cbSB += cSB[j] * cb[j]; }
        // CB variance of a linear combination c of the wide bins
        auto cbVar = [&](const std::vector<double>& c) {
            double var = 0;
            for (int q = 0; q < 2; ++q) {
                std::vector<double> wt(nw, 0.);
                for (int j = 0; j < nw; ++j) {
                    const double own = wl[q][s][j], oth = wl[1 - q][s][j];
                    wt[j] = c[j] * ((own > 0 && oth > 0) ? std::sqrt(oth / own) : 1.0);
                }
                int cur = -1; double acc = 0;
                for (const auto& e : lsr[q][s]) {
                    if (e.first != cur) { var += acc * acc; acc = 0; cur = e.first; }
                    acc += wt[e.second];
                }
                var += acc * acc;
            }
            return var;
        };
        std::vector<double> cD(nw);
        for (int j = 0; j < nw; ++j) cD[j] = cW[j] - r * cSB[j];
        o.DW = o.nW - o.cbW;   o.eDW = std::sqrt(vW[s] + cbVar(cW));
        o.DSB = o.nSB - o.cbSB; o.eDSB = std::sqrt(vS[s] + cbVar(cSB));
        o.D = o.DW - r * o.DSB; o.eD = std::sqrt(vD[s] + cbVar(cD));
        const double den = fW - r * fSB;
        o.S = o.D / den; o.eS = o.eD / den;
    }
    return out;
}

// share a / (a + b) of two independent parts
inline void Share(double a, double ea, double b, double eb, double& f, double& ef) {
    const double t = a + b;
    if (!(std::fabs(t) > 0)) { f = 0; ef = 0; return; }
    f = a / t;
    ef = std::sqrt(b * b * ea * ea + a * a * eb * eb) / (t * t);
}
// ratio of two independent numbers
inline void Ratio(double a, double ea, double b, double eb, double& f, double& ef) {
    if (!(std::fabs(b) > 0)) { f = 0; ef = 0; return; }
    f = a / b;
    ef = std::fabs(f) * std::sqrt((a != 0 ? ea * ea / (a * a) : 0) + eb * eb / (b * b));
}

// ============================================================================
// E. ECAL-tagged eta (frozen definition of plots/scan_bump.C)
// ============================================================================
// Rows with a good ECAL photon (g_idx > 0) and m_ee < 0.10 GeV/c2 whose
// M(e+e-gamma) lies in the pi0 window 0.119-0.173 (pi0 at 0.1459, 13.3 MeV, +-2 sigma;
// M(e+e-gamma) peaks 8 % above m_pi0 (ECAL not calibrated: the photon energy is about 17 % too high,
// median rescaling factor 0.856 in the eta window, verify_res_numbers.C). m5c =
// M(pi+pi-e+e-gamma) with the photon energy rescaled so that M(e+e-gamma) = m_pi0.
// eta peak in m5c: frozen Gaussian 0.5463 / 6.2 MeV + pol2 in 0.49-0.62, 2.5 MeV
// bins, CB in 20 MeV bins (P_main of the engine); T = Gaussian content.
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730, kMeeMax = 0.100;
const double kTagMean = 0.5463, kTagSigma = 0.0062, kTagLo = 0.49, kTagHi = 0.62;
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
    c.mean = kTagMean; c.sigma = kTagSigma; c.lo = kTagLo; c.hi = kTagHi; c.order = 2;
    c.winLo = c.mean - scan::kWinNsig * c.sigma; c.winHi = c.mean + scan::kWinNsig * c.sigma;
    return c;
}
inline scan::FillOpt TagOpt() {
    scan::FillOpt o;
    o.nbin = 160; o.lo = 0.40; o.hi = 0.80; o.reb = 8;
    o.x = [](const Row& r) { return Tagged(r) ? M5c(r) : -1.0; };
    const PeakCfg c = TagCfg();
    o.windows.push_back(std::make_pair(c.winLo, c.winHi));
    return o;
}
// model spread of T: orders 2, 3; ranges 0.49-0.62, 0.48-0.66, 0.50-0.60; width -+0.3 MeV (as scan_bump.C)
inline double TagSpread(const Spectra& s, int* kept = nullptr) {
    const PeakFit nom = scan::FitPeakMain(s, TagCfg());
    double mn = nom.S, mx = nom.S; int n = 0;
    const double rl[3] = {0.49, 0.48, 0.50}, rh[3] = {0.62, 0.66, 0.60};
    for (int ord = 2; ord <= 3; ++ord)
        for (int ir = 0; ir < 3; ++ir)
            for (int iw = -1; iw <= 1; ++iw) {
                PeakCfg c = TagCfg(); c.order = ord; c.lo = rl[ir]; c.hi = rh[ir]; c.sigma = kTagSigma + 0.0003 * iw;
                const PeakFit f = scan::FitPeakMain(s, c);
                if (!f.ok || !(f.prob > 0.05)) continue;
                mn = std::min(mn, f.S); mx = std::max(mx, f.S); n++;
            }
    if (kept) *kept = n;
    return 0.5 * (mx - mn);
}

// ============================================================================
// F. luminosity per day and per period
// ============================================================================
// L6_SF of each run, scaled by the merged share of its DSTs (two runs are partial;
// 58 runs without any merged DST are absent from the data). Sum = 6079.17 nb^-1.
inline bool LumiPerDay(std::map<int, double>& L) {
    std::ifstream in(kLumiCsv);
    if (!in) { printf("res::LumiPerDay: cannot open %s\n", kLumiCsv); return false; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.compare(0, 8, "dst_stem") == 0) continue;
        std::stringstream ss(line);
        std::string c[8];
        for (int k = 0; k < 8; ++k) std::getline(ss, c[k], ',');
        const int day = std::atoi(c[2].c_str());
        const double nl = std::atof(c[3].c_str()), nm = std::atof(c[4].c_str()), l6 = std::atof(c[5].c_str());
        if (nl > 0) L[day] += l6 * nm / nl;
    }
    return true;
}

} // namespace res

#endif // RES_COMMON_H
