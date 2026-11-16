// cb_common_cnt.h -- COPY of plots/exp_report/cb_common.h made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/cb
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// cb_common.h  --  shared code of the combinatorial-background chapter (prefix cb)
// ============================================================================
// Used by plots/exp_report/cb_checks.C and plots/exp_report/cb_kmix.C.
// Experimental data only. Contents:
//   - output directories of the chapter
//   - the named selections (on scan::Row; cor_ok == 1 is required by
//     scan::FillMany; vertex z > -500 mm and PT3 are applied in FAT)
//   - CB variants on an existing scan::Spectra (scale, per-bin factor k,
//     step / cubic interpolation, like-sign sum shape)
//   - like-sign spectra in the place of N - CB (null tests)
//   - the bump edge-step sum E over the four M(pi+pi-) slices of scan_bump.C
//   - plotting helpers (style, EPS -> PDF)
// ============================================================================
#ifndef CB_COMMON_CNT_H
#define CB_COMMON_CNT_H

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../ScanEngine.h"

namespace cbc {

using scan::Data;
using scan::PeakCfg;
using scan::PeakFit;
using scan::Row;
using scan::Sel;
using scan::Spectra;
using scan::StepCfg;
using scan::StepFit;
using scan::WinSum;

const char* const kFigOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp";
const char* const kTabOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/cb/tab";
const char* const kWork   = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/cb/work";
const char* const kCntScratch = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/cb/figs";
const char* const kCntScope = " cb_null cb_spectra ";
const char* const kKmixFile = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb/cb_kmix.root";

// ---- selection pieces --------------------------------------------------------------
inline bool Pre(const Row& r)   { return r.isBest >= 0 && r.vz > -160.f; }
inline bool Vtx(const Row& r)   { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; }
inline bool Oa2(const Row& r)   { return r.oa > 2.f; }
inline bool PeakK(const Row& r) { return r.oa_lab < 55.f && r.mpipi > 0.30f; }
inline bool BumpK(const Row& r) { return r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; }
inline bool DzPi(const Row& r)  { return r.dz_pi < 40.f; }
inline bool Mm4(const Row& r)   { return r.mm4 < 2.85f; }
inline bool MeeX(const Row& r)  { return r.mee >= 0.0155f && r.mee < 0.0245f; }

struct SD { std::string key, name; Sel sel; bool best = false; };
inline std::vector<SD> Selections() {
    std::vector<SD> v;
    auto add = [&](const char* k, const char* n, Sel s) { SD d; d.key = k; d.name = n; d.sel = s; v.push_back(d); };
    add("fat",      "FAT only (vz>-500)",       [](const Row& r) { return r.isBest >= 0; });
    add("pre",      "pre",                      [](const Row& r) { return Pre(r); });
    add("preoa",    "pre + OA>2",               [](const Row& r) { return Pre(r) && Oa2(r); });
    add("prevtx",   "pre + vertex",             [](const Row& r) { return Pre(r) && Vtx(r); });
    add("direct",   "DIRECT (pre+dz_pi<40)",    [](const Row& r) { return Pre(r) && DzPi(r); });
    add("peaksim",  "PEAK (DIRECT+MM4<2.85)",   [](const Row& r) { return Pre(r) && DzPi(r) && Mm4(r); });
    add("peaknom",  "peak nominal",             [](const Row& r) { return Pre(r) && PeakK(r); });
    add("peakvtx",  "peak nominal + vertex",    [](const Row& r) { return Pre(r) && PeakK(r) && Vtx(r); });
    add("peakoa",   "peak nominal + OA>2",      [](const Row& r) { return Pre(r) && PeakK(r) && Oa2(r); });
    add("bumpnom",  "bump nominal",             [](const Row& r) { return Pre(r) && BumpK(r); });
    add("bumpvtx",  "bump nominal + vertex",    [](const Row& r) { return Pre(r) && BumpK(r) && Vtx(r); });
    add("bumpoa",   "bump nominal + OA>2",      [](const Row& r) { return Pre(r) && BumpK(r) && Oa2(r); });
    add("x17",      "X17 (Mee 15.5-24.5+vtx)",  [](const Row& r) { return Pre(r) && MeeX(r) && Vtx(r); });
    const int nAll = (int)v.size();
    const char* bk[] = {"pre", "preoa", "prevtx", "direct", "peaknom", "peakvtx", "bumpnom", "bumpvtx"};
    for (const char* k : bk)
        for (int i = 0; i < nAll; ++i)
            if (v[i].key == k) {
                SD d = v[i]; d.key = std::string(k) + "_b1"; d.name = v[i].name + " [isBest==1]";
                d.sel = scan::SelBest(v[i].sel); d.best = true; v.push_back(d);
            }
    return v;
}
inline int Find(const std::vector<SD>& v, const std::string& k) {
    for (size_t i = 0; i < v.size(); ++i) if (v[i].key == k) return (int)i;
    printf("cb: no selection %s\n", k.c_str());
    return 0;
}

// ---- M(pi+pi-) slices of the bump edge-step fit (as plots/scan_bump.C) --------------
const int kNSig = 4;
const double kSlLo[kNSig] = {0.30, 0.34, 0.37, 0.39};
const double kSlHi[kNSig] = {0.34, 0.37, 0.39, 0.40};
inline StepCfg SliceCfg(int i) {
    StepCfg c; c.mpLo = kSlLo[i]; c.mpHi = kSlHi[i];
    if (i == 3) c.halfRange = 0.035;
    return c;
}
inline int SliceOf(double mpipi) {
    for (int i = 0; i < kNSig; ++i) if (mpipi >= kSlLo[i] && mpipi < kSlHi[i]) return i;
    return -1;
}

// ---- output ------------------------------------------------------------------------
inline FILE*& OutFile() { static FILE* f = nullptr; return f; }
inline void P(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    if (OutFile()) { va_start(ap, fmt); vfprintf(OutFile(), fmt, ap); va_end(ap); }
}
inline std::string TexName(std::string nm) {
    for (size_t p = 0; (p = nm.find('_', p)) != std::string::npos; p += 2) nm.replace(p, 1, "\\_");
    for (size_t p = 0; (p = nm.find('>', p)) != std::string::npos; p += 3) nm.replace(p, 1, "$>$");
    for (size_t p = 0; (p = nm.find('<', p)) != std::string::npos; p += 3) nm.replace(p, 1, "$<$");
    for (size_t p = 0; (p = nm.find("==", p)) != std::string::npos; p += 1) nm.replace(p, 2, "=");
    return nm;
}

// ============================================================================
// CB variants on an existing spectrum
// ============================================================================
// mode 0 linear interpolation (engine), 1 step (constant inside the 10 MeV bin),
// 2 cubic (Catmull-Rom through the wide-bin values), 3 like-sign sum shape
// (N++ + N--) x R (R = integrated 2 sqrt(N++ N--) / (N++ + N--)), 4 plain sum
// N++ + N-- (R = 1). f = global factor, kw = optional factor per wide bin.
// The fit covariance of the CB (Spectra::CbCov) keeps the linear-interpolation
// weights in every mode (approximation for modes 1 and 2); the window CB
// variance is scaled with the window CB.
inline Spectra Modify(const Spectra& in, double f, int mode = 0, const std::vector<double>* kw = nullptr) {
    Spectra s = in;
    double s1 = 0, s2 = 0;
    for (int j = 0; j < s.nw; ++j) { s1 += s.w[1][j]; s2 += s.w[2][j]; }
    const double R = (mode == 3 && s1 + s2 > 0) ? 2.0 * std::sqrt(s1 * s2) / (s1 + s2) : 1.0;
    for (int j = 0; j < s.nw; ++j) {
        const double k = f * (kw ? (*kw)[j] : 1.0);
        double base = 2.0 * std::sqrt(s.w[1][j] * s.w[2][j]), var = s.w[1][j] + s.w[2][j];
        if (mode >= 3) { base = (s.w[1][j] + s.w[2][j]) * R; var = (s.w[1][j] + s.w[2][j]) * R * R; }
        s.cbW[j] = k * base; s.cbWVar[j] = k * k * var;
    }
    for (int i = 0; i < s.nbin; ++i) {
        int j0, j1; double u; s.Interp(i, j0, j1, u);
        const double rb = s.reb;
        if (mode == 1) {
            const int j = std::min(i / s.reb, s.nw - 1);
            s.cb[i] = s.cbW[j] / rb; s.cbVar[i] = s.cbWVar[j] / (rb * rb);
        } else if (mode == 2 && j0 != j1 && j0 >= 1 && j1 <= s.nw - 2) {
            const int jj[4] = {j0 - 1, j0, j1, j1 + 1};
            const double u2 = u * u, u3 = u2 * u;
            const double w[4] = {0.5 * (-u3 + 2 * u2 - u), 0.5 * (3 * u3 - 5 * u2 + 2),
                                 0.5 * (-3 * u3 + 4 * u2 + u), 0.5 * (u3 - u2)};
            double c = 0, v = 0;
            for (int q = 0; q < 4; ++q) { c += w[q] * s.cbW[jj[q]] / rb; v += w[q] * w[q] * s.cbWVar[jj[q]] / (rb * rb); }
            s.cb[i] = c; s.cbVar[i] = v;
        } else {
            const double a0 = (1 - u) / rb, a1 = u / rb;
            s.cb[i] = a0 * s.cbW[j0] + a1 * s.cbW[j1];
            s.cbVar[i] = a0 * a0 * s.cbWVar[j0] + a1 * a1 * s.cbWVar[j1];
        }
        s.sig[i] = s.n[0][i] - s.cb[i];
        s.sigVar[i] = s.n[0][i] + s.cbVar[i];
    }
    for (size_t w = 0; w < s.win.size(); ++w) {
        WinSum& ws = s.win[w];
        double c = 0;
        for (int i = 0; i < s.nbin; ++i) { const double fr = s.Frac(i, ws.lo, ws.hi); if (fr > 0) c += fr * s.cb[i]; }
        const double old = in.win[w].cb;
        ws.cbVar = old > 0 ? in.win[w].cbVar * (c / old) * (c / old) : 0;
        ws.cb = c;
    }
    return s;
}

// like-sign spectrum k (1 = ++, 2 = --, 3 = sum) placed in "sig" (and in n[0]) with
// no CB, so that FitEdgeStep / FitPeakMain run on the like-sign rows (null tests).
// For k = 3 also n[1] holds the sum (for scan::FitPeakLikeSign(s, 1)).
inline Spectra LikeSignAsSignal(const Spectra& in, int k) {
    Spectra s = in;
    for (int i = 0; i < s.nbin; ++i) {
        const double y = (k == 3) ? in.n[1][i] + in.n[2][i] : in.n[k][i];
        s.n[0][i] = y; s.sig[i] = y; s.sigVar[i] = y; s.cb[i] = 0; s.cbVar[i] = 0;
    }
    for (int j = 0; j < s.nw; ++j) { s.cbW[j] = 0; s.cbWVar[j] = 0; }
    if (k == 3) for (int i = 0; i < s.nbin; ++i) s.n[1][i] = in.n[1][i] + in.n[2][i];
    return s;
}

// ---- E = sum of the step heights over the four slices -------------------------------
struct EVal { double E = 0, eE = 0; double h[kNSig], eh[kNSig]; bool ok = false; };
inline EVal EdgeSum(const std::vector<Spectra>& sl) {
    EVal e; double v = 0; int nok = 0;
    for (int i = 0; i < kNSig; ++i) {
        const StepFit f = scan::FitEdgeStep(sl[i], SliceCfg(i));
        e.h[i] = f.ok ? f.h : 0; e.eh[i] = f.ok ? f.eh : 0;
        if (!f.ok) continue;
        e.E += f.h; v += f.eh * f.eh; ++nok;
    }
    e.eE = std::sqrt(v); e.ok = (nok == kNSig);
    return e;
}
inline std::vector<Spectra> ModifyAll(const std::vector<Spectra>& sl, double f, int mode,
                                      const std::vector<std::vector<double>>* kw = nullptr) {
    std::vector<Spectra> o;
    for (size_t i = 0; i < sl.size(); ++i) o.push_back(Modify(sl[i], f, mode, kw ? &(*kw)[i] : nullptr));
    return o;
}

// ============================================================================
// plotting helpers
// ============================================================================
inline void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelSize(0.055, "xyz"); gStyle->SetTitleSize(0.060, "xyz");
    gStyle->SetTitleOffset(1.05, "x"); gStyle->SetTitleOffset(1.20, "y");
    gStyle->SetPadLeftMargin(0.16); gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.08); gStyle->SetPadBottomMargin(0.15);
    gStyle->SetEndErrorSize(0);
}
inline void Deco(TH1* h, int col, int mk, double ms = 0.7) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms); h->SetLineWidth(1);
}
inline void Deco(TGraph* g, int col, int mk, double ms = 0.9) {
    g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineWidth(2);
}
// vector PDF cropped to the canvas (EPS -> ps2pdf -dEPSCrop) and a PNG preview
inline void Save(TCanvas* c, const char* name) {
    const std::string n = name;
    const bool inScope = std::string(kCntScope).find(" " + n + " ") != std::string::npos;
    const std::string eb = std::string(kCntScratch) + "/" + n;                       // EPS always in the scratch folder
    const std::string fb = inScope ? std::string(kFigOut) + "/" + n + "_cnt" : eb;   // new name next to the original
    c->SaveAs((eb + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf", eb.c_str(), fb.c_str()));
    c->SaveAs((fb + ".png").c_str());
}
inline void Txt(double x, double y, const char* t, double size = 0.050, int col = kBlack) {
    static TLatex* tx = nullptr;
    if (!tx) { tx = new TLatex(); tx->SetNDC(); tx->SetTextFont(42); }
    tx->SetTextSize(size); tx->SetTextColor(col); tx->DrawLatex(x, y, t);
}

} // namespace cbc

#endif // CB_COMMON_H
