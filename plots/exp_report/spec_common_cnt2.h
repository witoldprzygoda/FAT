// spec_common_cnt2.h -- COPY of plots/exp_report/spec_common.h (via spec_common_cnt.h of figure agent A) made by the final checker (integ4, 2026-11-15): figure strings
// 'rows' -> 'counts' (analyst's terminology); figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png};
// tables, ROOT files and EPS -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/check/figs/spec. Physics code unchanged.
// spec_common_cnt.h -- COPY of plots/exp_report/spec_common.h made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/spec
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// spec_common.h  --  shared definitions of the "spectra and kinematics" chapter
//                    (prefix spec) of the experimental-data report
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22, 29 days), channel pi+ pi- e+ e-. DATA ONLY.
// Header-only; included by plots/exp_report/spec_m4.C, spec_vars.C and
// spec_vz500.C (and nothing else). Builds on ../ScanEngine.h (rows, CB, fits).
//
// Selections (rows = stored particle combinations; the FAT production already
// applies event vertex z > -500 mm and the PT3 trigger, trigbit == 8192):
//   pre          = cor_ok == 1 && isBest >= 0 && vz > -160 mm
//   +vertex      = both pions in the event-vertex fit
//                  (pip_isUsedVertex != 0 && pim_isUsedVertex != 0)
//   peak nominal = pre + oa_lab < 55 deg + M(pi+pi-) > 0.30 GeV/c2
//   bump nominal = pre + oa_lab < 40 deg + 0.30 < M(pi+pi-) < 0.41 GeV/c2
//   frozen simulation selections (defined blind in the simulation study):
//     PEAK   = pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2
//     DIRECT = pre + dz_pi < 40 mm
//     BUMP   = pre (no further cut)
//   variants: isBest == 1 (comparison), OA(e+e-) > 2 deg (not a fixed cut)
// Mass windows: data peak window 0.53364-0.55416, data bump window
// 0.4550-0.5225 GeV/c2 (ScanEngine kPeakLo/kPeakHi, kBumpLo/kBumpHi); the
// simulation-centred windows of the frozen selections (0.5333-0.5587,
// 0.4576-0.5255) are only marked for reference.
// CB = 2 sqrt(N++ N--) of the like-sign hypotheses, k = 1.
// Luminosity: L = 6079.2 nb^-1 (START-free elastic-scattering luminosity of the
// runs present in the data; all 17 615 runs: 6080.5 nb^-1); PT3 rows x 64/63.
// ============================================================================
#ifndef SPEC_COMMON_CNT2_H
#define SPEC_COMMON_CNT2_H

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../ScanEngine.h"

namespace spec {

using scan::Row;
using scan::Sel;

// ---- output locations ----------------------------------------------------------
const char* const kFigOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp";
const char* const kTabOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/check/figs/spec/tab";
const char* const kWork   = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/check/figs/spec/work";
const char* const kCntScratch = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/check/figs/spec/figs";
const char* const kCntScope = " spec_m4_overview spec_m4_isbest spec_vertexcut spec_vz500 spec_var_mee spec_var_oa spec_var_mpipi spec_var_oalab spec_var_oaeta spec_var_mm4 spec_var_pt4 spec_var_y4 spec_var_vz spec_var_dzpi ";

// ---- normalisation -------------------------------------------------------------
const double kLumi     = 6079.2;        // nb^-1, runs present in the data (START-free, six sectors)
const double kLumiAll  = 6080.5;        // nb^-1, all 17 615 runs (L_SF)
const double kLumiRel  = 0.069;         // total relative uncertainty of L_SF (5.56 % HADES syst, 1.03 % fit, 4.0 % norm.)
const double kPT3      = 64.0 / 63.0;   // trigger word 8192 holds 63/64 of the PT3 events
const double kCBsyst   = 0.03;          // +-3 % CB sensitivity quoted with every N - CB

// ---- windows [GeV/c2] ----------------------------------------------------------
const double kSimPeakLo = 0.5333, kSimPeakHi = 0.5587;   // frozen simulation windows (sim-centred)
const double kSimBumpLo = 0.4576, kSimBumpHi = 0.5255;

// ---- cut values -------------------------------------------------------------------
const float kVz     = -160.f;   // mm
const float kDzPi   = 40.f;     // mm
const float kMM4    = 2.85f;    // GeV/c2
const float kOAmin  = 2.f;      // deg (variant only)

inline bool PreNoVz(const Row& r) { return r.cor_ok == 1 && r.isBest >= 0; }
inline bool Pre(const Row& r)     { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > kVz; }
inline bool PiVtx(const Row& r)   { return r.pip_isUsedVertex != 0 && r.pim_isUsedVertex != 0; }
inline bool PeakKin(const Row& r) { return r.oa_lab < 55.f && r.mpipi > 0.30f; }
inline bool BumpKin(const Row& r) { return r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; }

// one named selection
struct NSel {
    std::string tag;     // short tag for files and tables
    std::string name;    // plain-text name (tables)
    std::string root;    // TLatex name (figures)
    Sel sel;
};

// the base selections of the chapter (isBest >= 0, no OA cut)
inline std::vector<NSel> BaseSelections() {
    std::vector<NSel> v;
    v.push_back({"pre", "pre (= frozen BUMP)", "pre (= frozen BUMP)",
                 [](const Row& r) { return Pre(r); }});
    v.push_back({"prevtx", "pre + vertex", "pre + #pi#pi in vertex fit",
                 [](const Row& r) { return Pre(r) && PiVtx(r); }});
    v.push_back({"peaknom", "peak nominal", "peak nominal",
                 [](const Row& r) { return Pre(r) && PeakKin(r); }});
    v.push_back({"peaknomvtx", "peak nominal + vertex", "peak nominal + vertex",
                 [](const Row& r) { return Pre(r) && PeakKin(r) && PiVtx(r); }});
    v.push_back({"bumpnom", "bump nominal", "bump nominal",
                 [](const Row& r) { return Pre(r) && BumpKin(r); }});
    v.push_back({"bumpnomvtx", "bump nominal + vertex", "bump nominal + vertex",
                 [](const Row& r) { return Pre(r) && BumpKin(r) && PiVtx(r); }});
    v.push_back({"PEAK", "frozen PEAK", "frozen PEAK (#Deltaz_{#pi} < 40 mm, MM_{4} < 2.85)",
                 [](const Row& r) { return Pre(r) && r.dz_pi < kDzPi && r.mm4 < kMM4; }});
    v.push_back({"DIRECT", "frozen DIRECT", "frozen DIRECT (#Deltaz_{#pi} < 40 mm)",
                 [](const Row& r) { return Pre(r) && r.dz_pi < kDzPi; }});
    return v;
}

// variant: OA(e+e-) > 2 deg added
inline Sel WithOA(const Sel& s) { return [s](const Row& r) { return r.oa > kOAmin && s(r); }; }

// ---- plot style ----------------------------------------------------------------
inline void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gStyle->SetLabelSize(0.050, "xyz"); gStyle->SetTitleSize(0.058, "xyz");
    gStyle->SetTitleOffset(1.05, "x"); gStyle->SetTitleOffset(1.30, "y");
    gStyle->SetPadLeftMargin(0.16); gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.08); gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0);
}
inline void Deco(TH1* h, int col, int mk, double ms = 0.7) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->SetLineWidth(1);
}
// TLatex in PostScript/PDF: a blank in front of a Symbol-font command (#minus,
// #geq, #pm, #pi, ...) is drawn about three blanks wide; replace it by a
// narrow blank (#scale[0.4]{ }). Text commands (#font, #it, ...) are left alone.
inline std::string Fx(const std::string& s) {
    static const char* const textCmd[] = {"#font", "#color", "#splitline", "#scale", "#kern", "#lower", "#it", "#bf", "#frac", "#bar", "#hat"};
    std::string o;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == ' ' && i + 1 < s.size() && s[i + 1] == '#') {
            bool text = false;
            for (const char* t : textCmd) if (s.compare(i + 1, std::strlen(t), t) == 0) text = true;
            if (!text) { o += "#scale[0.4]{ }"; continue; }
        }
        o += s[i];
    }
    return o;
}
// N-- in TLatex (two #minus would merge into one long dash)
const char* const kNmm = "N_{#minus#scale[0.3]{ }#minus}";
inline void Txt(double x, double y, const char* s, double size = 0.050, int align = 11, int col = kBlack) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextAlign(align); l.SetTextColor(col);
    l.DrawLatex(x, y, Fx(s).c_str());
}
inline TLegend* Leg(double x1, double y1, double x2, double y2, double size) {
    TLegend* l = new TLegend(x1, y1, x2, y2);
    l->SetTextSize(size); l->SetFillStyle(0); l->SetBorderSize(0); l->SetTextFont(42);
    return l;
}
inline void LE(TLegend* l, TObject* o, const char* s, const char* opt) { l->AddEntry(o, Fx(s).c_str(), opt); }
// light solid colours of the mass windows
inline int ColBump() { return TColor::GetColor("#d9f2c9"); }
inline int ColPeak() { return TColor::GetColor("#cfe0fb"); }
// solid light band (EPS has no transparency: draw it before the data and redraw
// the axes afterwards)
inline TBox* Band(double x1, double x2, double y1, double y2, int col) {
    TBox* b = new TBox(x1, y1, x2, y2);
    b->SetFillColor(col); b->SetFillStyle(1001); b->SetLineWidth(0);
    b->Draw();
    return b;
}
inline TLine* VLine(double x, double y1, double y2, int col, int sty, int w = 2) {
    TLine* l = new TLine(x, y1, x, y2);
    l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(w);
    l->Draw();
    return l;
}
// vector PDF cropped to the canvas (EPS -> ps2pdf -dEPSCrop, as in the simulation
// report; keeps transparency-free boxes vector) + PNG preview
inline void Save(TCanvas* c, const std::string& name) {
    const bool inScope = std::string(kCntScope).find(" " + name + " ") != std::string::npos;
    const std::string eb = std::string(kCntScratch) + "/" + name;                          // EPS always in the scratch folder
    const std::string base = inScope ? std::string(kFigOut) + "/" + name + "_cnt" : eb;   // new name next to the original
    c->SaveAs((base + ".png").c_str());
    c->SaveAs((eb + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf", eb.c_str(), base.c_str()));
}

// ---- merged display histogram of a Spectra component ---------------------------
// what = "us", "pp", "mm", "cb", "sig"; f = number of fine bins merged.
// CB and N - CB variances include the covariance of the interpolated CB
// (neighbouring fine bins share the wide CB bins).
inline TH1D* Merge(const scan::Spectra& sp, const char* what, int f, const char* name, double xscale = 1.0) {
    const int nb = sp.nbin / f;
    TH1D* h = new TH1D(name, "", nb, sp.lo * xscale, (sp.lo + nb * f * sp.bw) * xscale);
    h->SetDirectory(nullptr);
    const std::string k = what;
    for (int b = 0; b < nb; ++b) {
        double v = 0, var = 0;
        for (int i = b * f; i < (b + 1) * f; ++i) {
            if (k == "us") { v += sp.n[0][i]; var += sp.n[0][i]; }
            else if (k == "pp") { v += sp.n[1][i]; var += sp.n[1][i]; }
            else if (k == "mm") { v += sp.n[2][i]; var += sp.n[2][i]; }
            else if (k == "cb") { v += sp.cb[i]; }
            else { v += sp.sig[i]; var += sp.n[0][i]; }
        }
        if (k == "cb" || k == "sig")
            for (int i = b * f; i < (b + 1) * f; ++i)
                for (int j = b * f; j < (b + 1) * f; ++j) var += sp.CbCov(i, j);
        h->SetBinContent(b + 1, v);
        h->SetBinError(b + 1, std::sqrt(std::max(var, 0.)));
    }
    return h;
}

} // namespace spec

#endif
