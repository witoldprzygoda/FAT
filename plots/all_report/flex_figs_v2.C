// ============================================================================
// flex_figs_v2.C  --  package flex of report_all: figures of the flexible SMASH description, version 2 (figure fix)
// ============================================================================
// Copy of flex_figs.C (unchanged, kept) with the fixes asked for by the flex verification (flex_results.md Sect. 6.3):
//   1. Part III palette (report_all common colours, fix_wp1.C / fix_figE2.C): data black (blue in ratio panels after the fit),
//      CB red, G_D light grey (G_D omega/rho0 pions middle grey), G_C dark grey, G_eta violet, C2 orange, C3 magenta,
//      C4 green, C5 cyan. C3 and C4 are drawn separately (no longer one cyan "C3 + C4").
//   2. Negative components (s_GC = -2.28 in the 0.40-1.00 fit, s_C3 = -0.09 in F8 and F7bj) are no longer hidden inside the
//      cumulative stack: on linear plots the positive terms are stacked from zero and every negative term is drawn as a
//      hatched band BELOW ZERO; the fitted model mu (= stack top + negative band) is drawn as a solid line. In the log-scale
//      M(ee) projections the negative C3 is a hatched band at the stack top whose lower edge is the fitted model.
//   3. flex_mee: legends built from the drawn objects (ratio markers blue = after the fit, open black = SMASH as is), the
//      "SMASH as is" line thicker.
//   4. TLatex gaps ("isBest   >= 0", "other   eta", "sigma   x"): a blank before a symbol-font command is replaced by
//      #kern[0.5]{#cmd} (helper Lx); the last x tick label is no longer cut (wider right margin).
//   5. Panel chi2 = Pearson sum (N - mu)^2 / V with the fit variance V, so the panels add up to the global Pearson value.
//   6. NEW flex_scales_v2: every scale of the nominal fit with stat and stat (+) MC, the fit-method spread of the verification
//      (template-noise weight lambda 0.5-1.5, vflex_lamscan.json) as a box, SMASH (1) and the common eta scale of C2, C4, C5
//      (Wald, vflex_ties.json). The old comparison figure (variants and Part III) is kept as flex_scales_cmp_v2.
// Every output name gets the suffix _v2; the original figures are not touched.
// Reads /hdd3/.../flex/flex_fit.root, flex_results.txt, flex_bias_F8.txt and hists_v1.root (B_mee_R0/R1/R5).
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_figs_v2", kTRUE)' 'plots/all_report/flex_figs_v2.C+O'
// ============================================================================
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "defs.h"
#include "hists_io.h"

namespace flexg2 {
using namespace alldefs;
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/flex";
const char* const kVfl = "/hdd3/przygoda/HADES/FAT_out/all_report/vflex";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";

// ---------------------------------------------------------------- Part III palette
const int cCB = kRed - 7, cCBdark = kRed + 1;
const int cGDnom = kGray, cGDom = kGray + 1, cGC = kGray + 2, cGC2 = kGray + 3;
const int cGeta = kViolet - 9, cC2 = kOrange + 7, cC3 = kMagenta + 1, cC4 = kGreen + 2, cC5 = kCyan + 1;
const int cT3 = kYellow - 6;          // fake-lepton template: not in the Part III palette, not drawn in these figures
const int cFit = kBlue + 1;           // "after the fit" in the M(ee) ratio panels

// ---------------------------------------------------------------- TLatex: no gap before symbol-font characters
// ROOT 6.32 widens a blank that precedes a symbol-font character (Greek letters, >=, x, +-, (+), ...), both in PNG and PDF.
// " #geq" -> "#kern[0.5]{#geq}" gives a normal word space (checked in the PDF and the PNG output).
std::string Lx(const std::string& s) {
    static const std::set<std::string> sym = {
        "alpha", "beta", "gamma", "delta", "epsilon", "varepsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi",
        "pi", "rho", "sigma", "tau", "upsilon", "phi", "varphi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi",
        "Sigma", "Upsilon", "Phi", "Psi", "Omega", "geq", "leq", "times", "pm", "mp", "oplus", "otimes", "approx", "sim", "rightarrow",
        "leftarrow", "Rightarrow", "bullet", "circ", "infty", "propto", "neq", "equiv", "cdot", "diamond", "partial", "ll", "gg", "perp",
        "nabla", "downarrow", "uparrow", "leftrightarrow", "dagger"};
    std::string o;
    o.reserve(s.size() + 64);
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == ' ' && i + 1 < s.size() && s[i + 1] == '#') {
            size_t j = i + 2;
            while (j < s.size() && std::isalpha((unsigned char)s[j])) ++j;
            const std::string run = s.substr(i + 2, j - i - 2);
            std::string best;
            for (size_t L = run.size(); L > 0; --L)
                if (sym.count(run.substr(0, L))) { best = run.substr(0, L); break; }
            if (!best.empty()) { o += "#kern[0.5]{#" + best + "}"; i = i + 1 + best.size(); continue; }
        }
        o += s[i];
    }
    return o;
}
void Txt(TLatex& t, double x, double y, const std::string& s) { t.DrawLatex(x, y, Lx(s).c_str()); }

void Save(TCanvas* c, const std::string& base) {
    c->SaveAs((std::string(kFig) + "/" + base + "_v2.pdf").c_str());
    c->SaveAs((std::string(kFig) + "/" + base + "_v2.png").c_str());
}

// filled step polygon between two step functions (low = nullptr: constant base); bins with centre outside [xlo, xhi] skipped
TGraph* Poly(const TH1* up, const TH1* low, double base, double ymin, double xlo, double xhi) {
    std::vector<int> bins;
    for (int b = 1; b <= up->GetNbinsX(); ++b) {
        const double xc = up->GetBinCenter(b);
        if (xc > xlo && xc < xhi) bins.push_back(b);
    }
    TGraph* g = new TGraph();
    int k = 0;
    for (int b : bins) {
        const double y = std::max(ymin, up->GetBinContent(b));
        g->SetPoint(k++, up->GetBinLowEdge(b), y); g->SetPoint(k++, up->GetBinLowEdge(b + 1), y);
    }
    for (int q = (int)bins.size() - 1; q >= 0; --q) {
        const int b = bins[q];
        const double y = std::max(ymin, low ? low->GetBinContent(b) : base);
        g->SetPoint(k++, up->GetBinLowEdge(b + 1), y); g->SetPoint(k++, up->GetBinLowEdge(b), y);
    }
    return g;
}
// legend proxy with a fill (and optional hatch)
TH1D* LegBox(int col, int style = 1001, int lcol = -1) {
    static int n = 0;
    TH1D* h = new TH1D(Form("legbox%d", n++), "", 1, 0, 1);
    h->SetDirectory(nullptr);
    h->SetFillColor(col); h->SetFillStyle(style); h->SetLineColor(lcol >= 0 ? lcol : col); h->SetLineWidth(lcol >= 0 ? 2 : 1);
    return h;
}

// ---------------------------------------------------------------- results parser (as flex_figs.C)
struct FitInfo { std::vector<std::string> names; std::vector<double> v, es, em, et; std::vector<std::vector<double>> cor; double chi2P = 0, chi2H = 0; int ndf = 0; };
std::map<std::string, FitInfo> gFit;
std::map<std::string, double> gP3, gP3e;
std::string gNominal = "F8", gAlt = "F7bj";
void ReadResults() {
    std::ifstream in(Form("%s/flex_results.txt", kOut));
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream is(line);
        std::string key, tag; is >> key;
        if (key == "nominal") { is >> gNominal; continue; }
        if (key == "alternative") { is >> gAlt; continue; }
        if (key == "p3equiv") { std::string n; double v; while (is >> n >> v) gP3[n] = v; continue; }
        if (key == "p3equiv_err") { std::string n, w; double v, s, m, t; is >> n >> v >> w >> s >> w >> m >> w >> t; gP3e[n] = t; continue; }
        if (key == "names" || key == "values" || key == "err_stat" || key == "err_mc" || key == "err_total") {
            is >> tag;
            FitInfo& f = gFit[tag];
            if (key == "names") { f.names.clear(); std::string n; while (is >> n) f.names.push_back(n); continue; }
            std::vector<double> x; double d; while (is >> d) x.push_back(d);
            if (key == "values") f.v = x; else if (key == "err_stat") f.es = x; else if (key == "err_mc") f.em = x; else f.et = x;
            continue;
        }
        if (key == "gof") {
            is >> tag;
            std::string w; FitInfo& f = gFit[tag];
            while (is >> w) {
                if (w == "chi2P") is >> f.chi2P; else if (w == "chi2GLShyb") is >> f.chi2H; else if (w == "ndf") is >> f.ndf;
            }
            continue;
        }
        if (key == "cor_total") {
            int n; is >> tag >> n;
            FitInfo& f = gFit[tag];
            f.cor.assign(n, std::vector<double>(n, 0.));
            for (int a = 0; a < n; ++a) { std::getline(in, line); std::istringstream r(line); for (int b = 0; b < n; ++b) r >> f.cor[a][b]; }
        }
    }
    gP3e["C3"] = gP3e["C3C4"]; gP3e["C4"] = gP3e["C3C4"];
}
int Idx(const FitInfo& f, const std::string& n) { for (size_t k = 0; k < f.names.size(); ++k) if (f.names[k] == n) return (int)k; return -1; }
double Val(const std::string& tag, const std::string& n, double def = -999.) {
    if (!gFit.count(tag)) return def;
    const FitInfo& f = gFit[tag];
    const int k = Idx(f, n);
    return (k >= 0 && k < (int)f.v.size()) ? f.v[k] : def;
}
std::string ScTxt(double s) { return s < 0 ? Form("#times (%.2f)", s) : Form("#times %.2f", s); }

// ---------------------------------------------------------------- spectra (per slice and summed)
struct Part { std::string name, label, sname; int col; };
std::vector<Part> PartsFor(TDirectory* d, const std::string& h) {
    // bottom -> top in the Part III order: G_D, G_C, (T3), G_eta, C2, C3, C4, C5
    std::vector<Part> v;
    auto has = [&](const std::string& n) { return d->Get((n + "_" + h).c_str()) != nullptr; };
    if (has("GD")) v.push_back({"GD", "G_{D} non-#eta Dalitz pair", "s_GD", cGDnom});
    if (has("GD_nom")) v.push_back({"GD_nom", "G_{D} other pions", "s_GD_nom", cGDnom});
    if (has("GD_om")) v.push_back({"GD_om", "G_{D} #omega/#rho^{0} pions", "s_GD_om", cGDom});
    if (has("GC")) v.push_back({"GC", "G_{C} non-#eta conv. pair", "s_GC", cGC});
    if (has("GC_nom")) v.push_back({"GC_nom", "G_{C} other pions", "s_GC_nom", cGC});
    if (has("GC_om")) v.push_back({"GC_om", "G_{C} #omega/#rho^{0} pions", "s_GC_om", cGC2});
    if (has("T3")) v.push_back({"T3", "T3 fake lepton", "s_T3", cT3});
    if (has("Geta")) v.push_back({"Geta", "G_{#eta} other #eta decays", "s_Geta", cGeta});
    v.push_back({"C2", "C2 bump, #pi^{0} Dalitz", "s_C2", cC2});
    v.push_back({"C3", "C3 bump, conversion", "s_C3", cC3});
    v.push_back({"C4", "C4 #eta#rightarrow#pi^{+}#pi^{-}#gamma, conv.", "s_C4", cC4});
    v.push_back({"C5", "C5 #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}", "s_C5", cC5});
    return v;
}
TH1D* Get(TDirectory* d, const std::string& n) {
    TH1D* h = (TH1D*)d->Get(n.c_str());
    if (!h) return nullptr;
    TH1D* c = (TH1D*)h->Clone((n + Form("_c%p", (void*)d)).c_str());
    c->SetDirectory(nullptr);
    return c;
}
// per MeV/c2 and restricted to [lo, hi)
TH1D* PerMeV(TH1D* h, double lo, double hi) {
    std::vector<double> e;
    for (int b = 1; b <= h->GetNbinsX() + 1; ++b) { const double x = h->GetBinLowEdge(b); if (x >= lo - 1e-9 && x <= hi + 1e-9) e.push_back(x); }
    static int n = 0;
    TH1D* o = new TH1D(Form("%s_pm%d", h->GetName(), n++), "", (int)e.size() - 1, e.data());
    o->SetDirectory(nullptr);
    for (int b = 1; b <= o->GetNbinsX(); ++b) {
        const int bb = h->FindBin(o->GetBinCenter(b));
        const double w = 1000. * h->GetBinWidth(bb);
        o->SetBinContent(b, h->GetBinContent(bb) / w); o->SetBinError(b, h->GetBinError(bb) / w);
    }
    return o;
}

// isFit: the directory is a fit (legend gives the fitted scales); notes: lines of the footnote
void DrawSpectra(TFile* f, const std::string& dir, const std::string& cls, double lo, double hi, const std::string& base, const std::string& title,
                 const std::vector<std::string>& notes, bool sumOnly, bool isFit) {
    TDirectory* d = (TDirectory*)f->Get(dir.c_str());
    if (!d) { printf("flex_figs_v2: %s missing\n", dir.c_str()); return; }
    const int ncol = sumOnly ? 1 : 6;
    const double xLeg = sumOnly ? 0.665 : 0.835;   // panels left of xLeg, legend right of it
    const double yNote = sumOnly ? 0.085 : 0.075;  // footnote strip below
    const double yTop = sumOnly ? 0.94 : 0.945;
    const double ySplit = yNote + (yTop - yNote) * 0.335;
    TCanvas* c = new TCanvas(("c" + base).c_str(), "", sumOnly ? 1350 : 2400, sumOnly ? 860 : 1000);
    c->cd();
    TPad* top = new TPad(("top" + base).c_str(), "", 0, yTop, 1, 1); top->Draw();
    top->cd(); { TLatex tt; tt.SetNDC(); tt.SetTextSize(sumOnly ? 0.50 : 0.50); tt.SetTextFont(42); Txt(tt, 0.01, 0.28, title); }
    std::vector<std::string> hs;
    for (int j = 1; j <= 5; ++j) hs.push_back(Form("S%d%s", j, cls.c_str()));
    const char* slTxt[6] = {"S1: M(#pi^{+}#pi^{-}) < 0.34", "S2: 0.34-0.37", "S3: 0.37-0.39", "S4: 0.39-0.43", "S5: M(#pi^{+}#pi^{-}) #geq 0.43", "sum of the slices"};
    std::vector<Part> parts = PartsFor(d, hs[0]);
    // legend proxies
    TH1D* legN = nullptr; TH1D* legMu = nullptr;
    std::vector<bool> isNeg(parts.size(), false);
    const double kE = isFit ? Val(dir, "k_E2", 1.0) : 1.0;
    const double fsz = sumOnly ? 1.0 : 1.0;
    for (int col = 0; col < ncol; ++col) {
        const int jcol = sumOnly ? 5 : col;
        c->cd();
        const double x0 = xLeg * col / ncol, x1 = xLeg * (col + 1) / ncol;
        TPad* p1 = new TPad(Form("p1_%s_%d", base.c_str(), col), "", x0, ySplit, x1, yTop);
        TPad* p2 = new TPad(Form("p2_%s_%d", base.c_str(), col), "", x0, yNote, x1, ySplit);
        const double lm = sumOnly ? 0.15 : 0.235, rm = sumOnly ? 0.035 : 0.075;
        p1->SetBottomMargin(0.015); p1->SetTopMargin(0.075); p1->SetLeftMargin(lm); p1->SetRightMargin(rm);
        p2->SetTopMargin(0.03); p2->SetBottomMargin(0.30); p2->SetLeftMargin(lm); p2->SetRightMargin(rm);
        p1->Draw(); p2->Draw();
        std::vector<std::string> use;
        if (jcol < 5) use.push_back(hs[jcol]); else use = hs;
        TH1D *hN = nullptr, *hCB = nullptr, *hMu = nullptr, *hEG = nullptr;
        std::vector<TH1D*> hp(parts.size(), nullptr);
        double ch = 0; int nb = 0;
        for (const auto& h : use) {
            auto acc = [&](TH1D*& a, TH1D* b) { if (!b) return; if (!a) a = b; else { a->Add(b); delete b; } };
            // Pearson per slice with the fit variance V (= mu error^2), as the global chi2P
            TH1D* n1 = Get(d, "N_" + h); TH1D* m1 = Get(d, "mu_" + h);
            if (n1 && m1)
                for (int b = 1; b <= n1->GetNbinsX(); ++b) {
                    const double x = n1->GetBinCenter(b);
                    const double V = m1->GetBinError(b) * m1->GetBinError(b);
                    if (x < lo || x > hi || V <= 0) continue;
                    ch += std::pow(n1->GetBinContent(b) - m1->GetBinContent(b), 2) / V; ++nb;
                }
            acc(hN, n1); acc(hMu, m1); acc(hCB, Get(d, "kCB_" + h)); acc(hEG, Get(d, "vEG_" + h));
            for (size_t q = 0; q < parts.size(); ++q) acc(hp[q], Get(d, parts[q].name + "_" + h));
        }
        if (!hN) continue;
        for (int b = 1; b <= hN->GetNbinsX(); ++b) hN->SetBinError(b, std::sqrt(std::max(0., hEG->GetBinContent(b))));
        TH1D* N = PerMeV(hN, lo, hi);
        TH1D* CB = PerMeV(hCB, lo, hi);
        TH1D* MU = PerMeV(hMu, lo, hi);
        // positive terms stacked from zero, negative terms below zero
        std::vector<std::pair<TH1D*, int>> layers;   // cumulative top, colour
        TH1D* cum = (TH1D*)CB->Clone(Form("cum_%s_%d", base.c_str(), col)); cum->SetDirectory(nullptr);
        layers.push_back({(TH1D*)cum->Clone(Form("lay_%s_%d_cb", base.c_str(), col)), cCB});
        struct NegBand { TH1D *up, *low; int col; };
        std::vector<NegBand> negs;
        TH1D* ncum = (TH1D*)CB->Clone(Form("ncum_%s_%d", base.c_str(), col)); ncum->Reset(); ncum->SetDirectory(nullptr);
        std::vector<TH1D*> xs(parts.size(), nullptr);
        for (size_t q = 0; q < parts.size(); ++q) {
            if (!hp[q]) continue;
            TH1D* x = PerMeV(hp[q], lo, hi);
            xs[q] = x;
            if (x->Integral() < 0) { isNeg[q] = true; continue; }
            cum->Add(x);
            layers.push_back({(TH1D*)cum->Clone(Form("lay_%s_%d_%zu", base.c_str(), col, q)), parts[q].col});
        }
        // negative terms below zero, mirrored order (the term nearest the top of the list sits next to zero)
        for (int q = (int)parts.size() - 1; q >= 0; --q) {
            if (!xs[q] || !isNeg[q]) continue;
            TH1D* up = (TH1D*)ncum->Clone(Form("nu_%s_%d_%d", base.c_str(), col, q));
            ncum->Add(xs[q]);
            TH1D* low = (TH1D*)ncum->Clone(Form("nl_%s_%d_%d", base.c_str(), col, q));
            negs.push_back({up, low, parts[q].col});
        }
        p1->cd();
        double ymax = 0, ymin = 0, sN = 0, sMu = 0;
        for (int b = 1; b <= N->GetNbinsX(); ++b) {
            ymax = std::max(ymax, std::max(N->GetBinContent(b) + N->GetBinError(b), std::max(MU->GetBinContent(b), cum->GetBinContent(b))));
            ymin = std::min(ymin, ncum->GetBinContent(b));
            const double w = 1000. * N->GetBinWidth(b);
            sN += N->GetBinContent(b) * w; sMu += MU->GetBinContent(b) * w;
        }
        TH1D* fr = (TH1D*)N->Clone(Form("fr_%s_%d", base.c_str(), col)); fr->Reset();
        fr->SetMaximum(1.12 * ymax); fr->SetMinimum(ymin < 0 ? 1.25 * ymin : 0.);
        fr->GetYaxis()->SetTitle("rows per MeV/c^{2}");
        fr->GetYaxis()->SetTitleSize((sumOnly ? 0.050 : 0.064) * fsz); fr->GetYaxis()->SetLabelSize((sumOnly ? 0.046 : 0.060) * fsz);
        fr->GetYaxis()->SetTitleOffset(sumOnly ? 1.45 : 1.62);
        fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetMaxDigits(4); fr->GetXaxis()->SetNdivisions(505);
        fr->Draw("AXIS");
        for (int q = (int)layers.size() - 1; q >= 0; --q) {
            TGraph* g = Poly(layers[q].first, nullptr, 0., -1e30, lo, hi);
            g->SetFillColor(layers[q].second); g->SetFillStyle(1001); g->SetLineColor(layers[q].second); g->Draw("F");
        }
        for (auto& nbd : negs) {
            TGraph* g = Poly(nbd.up, nbd.low, 0., -1e30, lo, hi);
            g->SetFillColor(nbd.col); g->SetFillStyle(3354); g->SetLineColor(nbd.col); g->Draw("F");
            TH1D* edge = (TH1D*)nbd.low->Clone(); edge->SetLineColor(nbd.col); edge->SetLineWidth(1); edge->SetFillStyle(0); edge->Draw("HIST SAME");
        }
        if (!negs.empty()) { TLine* z = new TLine(lo, 0, hi, 0); z->SetLineColor(kBlack); z->SetLineWidth(1); z->Draw(); }
        if (hi > 0.7001) {   // end of the nominal fit range
            TLine* l7 = new TLine(0.70, fr->GetMinimum(), 0.70, fr->GetMaximum()); l7->SetLineStyle(3); l7->SetLineColor(kGray + 3); l7->SetLineWidth(2); l7->Draw();
            TLatex t7; t7.SetTextSize(0.034); t7.SetTextFont(42); t7.SetTextAlign(31); t7.SetTextColor(kGray + 3);
            Txt(t7, 0.695, fr->GetMinimum() + 0.04 * (fr->GetMaximum() - fr->GetMinimum()), "nominal fit range ends at 0.70");
        }
        TH1D* MUl = (TH1D*)MU->Clone(Form("mul_%s_%d", base.c_str(), col)); MUl->SetLineColor(kBlack); MUl->SetLineWidth(1); MUl->SetFillStyle(0); MUl->Draw("HIST SAME");
        N->SetMarkerStyle(20); N->SetMarkerSize(sumOnly ? 0.9 : 0.8); N->SetMarkerColor(kBlack); N->SetLineColor(kBlack); N->Draw("E1 SAME");
        fr->Draw("AXIS SAME");
        TLatex lt; lt.SetNDC(); lt.SetTextSize((sumOnly ? 0.050 : 0.064) * fsz); lt.SetTextFont(42); Txt(lt, lm + 0.02, 0.945, slTxt[jcol]);
        {
            TLatex ls; ls.SetNDC(); ls.SetTextSize(sumOnly ? 0.042 : 0.054); ls.SetTextFont(42);
            Txt(ls, lm + 0.03, 0.86, Form("#chi^{2}_{P} %.0f / %d bins", ch, nb));
            if (jcol == 5) { if (!sumOnly) ls.SetTextSize(0.050); Txt(ls, lm + 0.03, sumOnly ? 0.80 : 0.79, Form("model/data, summed: %.3f", sMu / sN)); }
        }
        legN = N; legMu = MUl;
        // ratio
        p2->cd();
        TH1D* R = (TH1D*)N->Clone(Form("R_%s_%d", base.c_str(), col));
        R->Divide(MU);
        for (int b = 1; b <= R->GetNbinsX(); ++b) R->SetBinError(b, MU->GetBinContent(b) > 0 ? N->GetBinError(b) / MU->GetBinContent(b) : 0.);
        R->SetMinimum(0.55); R->SetMaximum(1.45);
        const double rs = sumOnly ? 0.085 : 0.094;
        R->GetYaxis()->SetTitle("data / model"); R->GetYaxis()->SetTitleSize(rs * 1.05); R->GetYaxis()->SetLabelSize(rs);
        R->GetYaxis()->SetTitleOffset(sumOnly ? 0.85 : 1.12); R->GetYaxis()->SetNdivisions(505);
        R->GetXaxis()->SetTitle(Lx("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]").c_str()); R->GetXaxis()->SetTitleSize(rs * 1.05); R->GetXaxis()->SetLabelSize(rs);
        R->GetXaxis()->SetTitleOffset(1.12); R->GetXaxis()->SetNdivisions(505); R->GetXaxis()->SetLabelOffset(0.02);
        R->SetMarkerStyle(20); R->SetMarkerSize(sumOnly ? 0.7 : 0.6);
        R->Draw("E1");
        TLine* l1 = new TLine(lo, 1, hi, 1); l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
        if (hi > 0.7001) { TLine* l7 = new TLine(0.70, 0.55, 0.70, 1.45); l7->SetLineStyle(3); l7->SetLineColor(kGray + 3); l7->SetLineWidth(2); l7->Draw(); }
    }
    // legend column
    c->cd();
    TPad* pl = new TPad(("pleg" + base).c_str(), "", xLeg, ySplit - 0.02, 1, yTop); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.02, 0.03, 0.99, 0.97); lg->SetTextSize(sumOnly ? 0.044 : 0.050); lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetMargin(0.17);
    if (legN) {
        lg->AddEntry(legN, Lx("data N (S0, isBest #geq 0)").c_str(), "PE");
        lg->AddEntry(legMu, isFit ? "fitted model (total)" : "SMASH as is + CB (total)", "L");
        for (int q = (int)parts.size() - 1; q >= 0; --q) {
            std::string lab = parts[q].label;
            const double s = Val(dir, parts[q].sname, 1.0);
            if (isFit) lab += " " + ScTxt(s);
            if (isNeg[q]) { lab += " < 0"; lg->AddEntry(LegBox(parts[q].col, 3354, parts[q].col), Lx(lab).c_str(), "F"); }
            else lg->AddEntry(LegBox(parts[q].col), Lx(lab).c_str(), "F");
        }
        lg->AddEntry(LegBox(cCB), Lx(isFit ? Form("k CB, k = %.2f", kE) : "CB = 2#sqrt{N_{++}N_{--}}, k = 1").c_str(), "F");
    }
    lg->Draw();
    c->cd();
    TLatex ln; ln.SetNDC(); ln.SetTextSize(sumOnly ? 0.021 : 0.0185); ln.SetTextFont(42);
    for (size_t k = 0; k < notes.size(); ++k) Txt(ln, 0.01, (sumOnly ? 0.050 : 0.043) - (sumOnly ? 0.030 : 0.026) * k, notes[k]);
    Save(c, base);
}

// ---------------------------------------------------------------- M(ee) projections (hists_v1 BEFORE layouts)
void DrawMee(TFile* fh) {
    const FitInfo& F = gFit[gNominal];
    auto sc = [&](const std::string& n) { const int k = Idx(F, n); return k >= 0 ? F.v[k] : -999.; };
    struct C { std::string name, label; std::vector<int> parts; int col; double s; };
    std::vector<C> comps;   // bottom -> top (Part III order)
    const double sGD = sc("s_GD"), sGDom = sc("s_GD_om"), sGDnom = sc("s_GD_nom");
    if (sGD > -900) comps.push_back({"GD", "G_{D}", {kGD}, cGDnom, sGD});
    else {
        comps.push_back({"GD_nom", "G_{D} other pions", {kGpi0Dnom, kGvecnom, kGothnom}, cGDnom, sGDnom});
        comps.push_back({"GD_om", "G_{D} #omega/#rho^{0} pions", {kGpi0Dom, kGvecom, kGothom}, cGDom, sGDom});
    }
    comps.push_back({"GC", "G_{C} conv. pair", {kGconv}, cGC, sc("s_GC")});
    comps.push_back({"Geta", "G_{#eta} other #eta", {kGeta}, cGeta, sc("s_Geta")});
    comps.push_back({"C2", "C2 bump, #pi^{0} Dalitz", {kC2}, cC2, sc("s_C2")});
    comps.push_back({"C3", "C3 bump, conv.", {kC3}, cC3, sc("s_C3")});
    comps.push_back({"C4", "C4 #eta#rightarrow#pi^{+}#pi^{-}#gamma", {kC4}, cC4, sc("s_C4")});
    comps.push_back({"C5", "C5 #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}", {kC5}, cC5, sc("s_C5")});
    const double kE2 = sc("k_E2");
    const int regs[3] = {0, 1, 5};
    const char* rtxt[3] = {"M4 < 0.455 GeV/c^{2}", "bump window B: 0.455-0.5225", "CR-U: 0.565-0.700"};
    TCanvas* c = new TCanvas("cmee2", "", 2300, 1050);
    c->cd();
    TPad* top = new TPad("tmee2", "", 0, 0.935, 1, 1); top->Draw(); top->cd();
    TLatex tt; tt.SetNDC(); tt.SetTextSize(0.42); tt.SetTextFont(42);
    Txt(tt, 0.01, 0.28, Form("M(e^{+}e^{-}) projections, S0, all slices: data vs SMASH as is (dashed) and after the flexible fit %s (stacked); the grey region below 10 MeV is not fitted", gNominal.c_str()));
    const double xmax = 0.16;
    const double yNote = 0.06, ySplit = 0.36;
    TH1D *legNd = nullptr, *legR1 = nullptr, *legR0 = nullptr, *legFit = nullptr, *legAsis = nullptr;
    double negFracHi = 0, negFracLo = 0;   // largest |negative term| / fitted model above and below 10 MeV
    for (int q = 0; q < 3; ++q) {
        const int r = regs[q];
        c->cd();
        TPad* p1 = new TPad(Form("m1_%d", q), "", 0.80 * q / 3., ySplit, 0.80 * (q + 1) / 3., 0.935);
        TPad* p2 = new TPad(Form("m2_%d", q), "", 0.80 * q / 3., yNote, 0.80 * (q + 1) / 3., ySplit);
        p1->SetBottomMargin(0.015); p1->SetLeftMargin(0.17); p1->SetRightMargin(0.04); p1->SetTopMargin(0.075); p1->SetLogy();
        p2->SetTopMargin(0.03); p2->SetBottomMargin(0.28); p2->SetLeftMargin(0.17); p2->SetRightMargin(0.04);
        p1->Draw(); p2->Draw();
        const std::string lay = BVarName(kBmee, r);
        TH1* N = allio::Data(fh, "S0", "all", "N", lay);
        TH1* CB = allio::Data(fh, "S0", "all", "CB", lay);
        TH1D* asis = (TH1D*)N->Clone(Form("asis2_%d", q)); asis->Reset();
        TH1D* fit = (TH1D*)N->Clone(Form("fit2_%d", q)); fit->Reset();
        TH1D* cbk = (TH1D*)N->Clone(Form("cbk2_%d", q)); cbk->Reset();
        for (int b = 1; b <= N->GetNbinsX(); ++b) {
            const double k = N->GetBinCenter(b) >= 0.010 ? kE2 : 1.0;
            cbk->SetBinContent(b, k * CB->GetBinContent(b));
            asis->SetBinContent(b, CB->GetBinContent(b));
            fit->SetBinContent(b, k * CB->GetBinContent(b));
        }
        std::vector<std::pair<TH1D*, int>> layers;
        TH1D* cum = (TH1D*)cbk->Clone(Form("mcum_%d", q));
        layers.push_back({(TH1D*)cum->Clone(Form("mlay_%d_cb", q)), cCB});
        std::vector<std::pair<TH1D*, int>> negs;   // negative terms (absolute values)
        for (size_t k = 0; k < comps.size(); ++k) {
            const auto& cc = comps[k];
            TH1D* h = (TH1D*)N->Clone(Form("mc2_%s_%d", cc.name.c_str(), q)); h->Reset();
            for (int p : cc.parts) { TH1* x = allio::Sim(fh, "S0", "nom", p, lay); h->Add(x, kNbToRaw); delete x; }
            asis->Add(h);
            h->Scale(cc.s);
            fit->Add(h);
            if (cc.s >= 0) { cum->Add(h); layers.push_back({(TH1D*)cum->Clone(Form("mlay_%d_%zu", q, k)), cc.col}); }
            else { TH1D* a = (TH1D*)h->Clone(Form("mneg_%d_%zu", q, k)); a->Scale(-1.); negs.push_back({a, cc.col}); }
        }
        for (auto& ng : negs)
            for (int b = 1; b <= N->GetNbinsX(); ++b) {
                if (N->GetBinCenter(b) > xmax || fit->GetBinContent(b) <= 0) continue;
                const double fr0 = ng.first->GetBinContent(b) / fit->GetBinContent(b);
                if (N->GetBinCenter(b) >= 0.010) negFracHi = std::max(negFracHi, fr0); else negFracLo = std::max(negFracLo, fr0);
            }
        p1->cd();
        TH1D* fr = (TH1D*)N->Clone(Form("mfr2_%d", q)); fr->Reset();
        fr->GetXaxis()->SetRangeUser(0, xmax);
        double ymax = 0, ymin = 1e9;
        for (int b = 1; b <= N->GetNbinsX(); ++b) if (N->GetBinCenter(b) < xmax) { ymax = std::max(ymax, N->GetBinContent(b)); if (N->GetBinContent(b) > 0) ymin = std::min(ymin, N->GetBinContent(b)); }
        fr->SetMaximum(3 * ymax); fr->SetMinimum(std::max(0.5, 0.3 * ymin));
        fr->GetYaxis()->SetTitle("rows per 2 MeV/c^{2}"); fr->GetYaxis()->SetTitleSize(0.060); fr->GetYaxis()->SetLabelSize(0.054); fr->GetYaxis()->SetTitleOffset(1.32); fr->GetXaxis()->SetLabelSize(0);
        fr->Draw("AXIS");
        const double fmin = fr->GetMinimum();
        for (int l = (int)layers.size() - 1; l >= 0; --l) {
            TGraph* g = Poly(layers[l].first, nullptr, fmin, fmin, 0., xmax);
            g->SetFillColor(layers[l].second); g->SetFillStyle(1001); g->Draw("F");
        }
        // negative terms: hatched band from the fitted model up to the top of the positive stack
        {
            TH1D* lowc = (TH1D*)fit->Clone(Form("mnl_%d", q));
            for (auto& ng : negs) {
                TH1D* upc = (TH1D*)lowc->Clone(Form("mnu_%d_%d", q, ng.second)); upc->Add(ng.first);
                TGraph* g = Poly(upc, lowc, fmin, fmin, 0., xmax);
                g->SetFillColor(ng.second); g->SetFillStyle(3354); g->Draw("F");
                lowc = upc;
            }
        }
        TH1D* fitl = (TH1D*)fit->Clone(Form("fitl_%d", q)); fitl->SetLineColor(cFit); fitl->SetLineWidth(2); fitl->SetFillStyle(0);
        fitl->GetXaxis()->SetRangeUser(0, xmax); fitl->Draw("HIST SAME");
        asis->SetLineColor(kBlack); asis->SetLineStyle(7); asis->SetLineWidth(3); asis->SetFillStyle(0);
        asis->GetXaxis()->SetRangeUser(0, xmax); asis->Draw("HIST SAME");
        TH1D* Nd = (TH1D*)N->Clone(Form("Nd2_%d", q));
        Nd->SetMarkerStyle(20); Nd->SetMarkerSize(0.75); Nd->SetMarkerColor(kBlack); Nd->SetLineColor(kBlack);
        Nd->GetXaxis()->SetRangeUser(0, xmax); Nd->Draw("E1 SAME");
        TBox* bx = new TBox(0, fr->GetMinimum(), 0.010, fr->GetMaximum()); bx->SetFillColorAlpha(kGray, 0.40); bx->Draw();
        fr->Draw("AXIS SAME");
        TLatex lt; lt.SetNDC(); lt.SetTextSize(0.060); lt.SetTextFont(42); Txt(lt, 0.18, 0.94, rtxt[q]);
        legNd = Nd; legFit = fitl; legAsis = asis;
        p2->cd();
        TH1D* R1 = (TH1D*)N->Clone(Form("mr1_%d", q)); R1->Divide(fit);
        TH1D* R0 = (TH1D*)N->Clone(Form("mr0_%d", q)); R0->Divide(asis);
        for (int b = 1; b <= N->GetNbinsX(); ++b) {
            R1->SetBinError(b, fit->GetBinContent(b) > 0 ? N->GetBinError(b) / fit->GetBinContent(b) : 0);
            R0->SetBinError(b, asis->GetBinContent(b) > 0 ? N->GetBinError(b) / asis->GetBinContent(b) : 0);
        }
        R1->GetXaxis()->SetRangeUser(0, xmax);
        R1->SetMinimum(0.3); R1->SetMaximum(1.9);
        R1->GetYaxis()->SetTitle("data / model"); R1->GetYaxis()->SetTitleSize(0.098); R1->GetYaxis()->SetLabelSize(0.088); R1->GetYaxis()->SetTitleOffset(0.82); R1->GetYaxis()->SetNdivisions(505);
        R1->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); R1->GetXaxis()->SetTitleSize(0.098); R1->GetXaxis()->SetLabelSize(0.088); R1->GetXaxis()->SetTitleOffset(1.18); R1->GetXaxis()->SetNdivisions(505);
        R1->SetMarkerStyle(20); R1->SetMarkerSize(0.65); R1->SetMarkerColor(cFit); R1->SetLineColor(cFit);
        R0->SetMarkerStyle(24); R0->SetMarkerSize(0.65); R0->SetMarkerColor(kBlack); R0->SetLineColor(kBlack);
        R1->Draw("E1"); R0->Draw("E1 SAME"); R1->Draw("E1 SAME");
        TLine* l1 = new TLine(0, 1, xmax, 1); l1->SetLineStyle(2); l1->Draw();
        TBox* bx2 = new TBox(0, 0.3, 0.010, 1.9); bx2->SetFillColorAlpha(kGray, 0.40); bx2->Draw();
        legR1 = R1; legR0 = R0;
    }
    c->cd();
    TPad* pl = new TPad("mleg2", "", 0.80, ySplit - 0.05, 1.0, 0.935); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.03, 0.02, 0.99, 0.98); lg->SetTextSize(0.049); lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetMargin(0.17);
    lg->AddEntry(legNd, "data N", "PE");
    lg->AddEntry(legFit, "fitted model (total)", "L");
    lg->AddEntry(legAsis, "SMASH as is + CB", "L");
    for (int k = (int)comps.size() - 1; k >= 0; --k) {
        std::string lab = comps[k].label + " " + ScTxt(comps[k].s);
        if (comps[k].s < 0) { lab += ", hatched"; lg->AddEntry(LegBox(comps[k].col, 3354, comps[k].col), Lx(lab).c_str(), "F"); }
        else lg->AddEntry(LegBox(comps[k].col), Lx(lab).c_str(), "F");
    }
    lg->AddEntry(LegBox(cCB), Lx(Form("k CB, k = %.2f (#geq 10 MeV)", kE2)).c_str(), "F");
    lg->AddEntry(LegBox(kGray, 1001), "M(ee) < 10 MeV/c^{2}: not fitted", "F");
    lg->Draw();
    c->cd();
    TPad* pr = new TPad("mleg2r", "", 0.80, yNote, 1.0, ySplit - 0.05); pr->Draw(); pr->cd();
    TLegend* lr = new TLegend(0.03, 0.30, 0.99, 0.95); lr->SetHeader("ratio panels:"); lr->SetTextSize(0.105); lr->SetFillStyle(0); lr->SetBorderSize(0); lr->SetMargin(0.17);
    lr->AddEntry(legR1, "data / fitted model", "PE");
    lr->AddEntry(legR0, "data / SMASH as is", "PE");
    lr->Draw();
    c->cd();
    TLatex ln; ln.SetNDC(); ln.SetTextSize(0.019); ln.SetTextFont(42);
    Txt(ln, 0.01, 0.034, Form("hists_v1.root B_mee_R0, R1, R5 (unblinded regions); simulation #times %.1f raw rows per nb (L = 6079.2 nb^{-1} #times 63/64); scales of %s (no tilt); below 10 MeV/c^{2} the CB is taken at k = 1.",
                              kNbToRaw, gNominal.c_str()));
    Txt(ln, 0.01, 0.010, Form("C3 enters with s = %.2f: it is not stacked; the hatched band (too thin to see: at most %.1f %% of the model above 10 MeV, %.1f %% below) is the part it removes; blue line = fitted model. R0 includes M4 0.30-0.40.",
                              sc("s_C3"), 100 * negFracHi, 100 * negFracLo));
    printf("flex_mee_v2: largest negative C3 fraction of the model: %.4f (>= 10 MeV), %.4f (< 10 MeV)\n", negFracHi, negFracLo);
    Save(c, "flex_mee");
}

// ---------------------------------------------------------------- bias test results (scenario 0)
struct Verd { std::map<std::string, std::string> v; std::map<std::string, double> bias; };
Verd ReadBias(const std::string& tag) {
    Verd o;
    std::ifstream in(Form("%s/flex_bias_%s.txt", kOut, tag.c_str()));
    std::string line; int scen = -1;
    while (std::getline(in, line)) {
        if (line.rfind("# scenario 0", 0) == 0) scen = 0;
        else if (line.rfind("# scenario 1", 0) == 0) scen = 1;
        if (scen != 0 || line.rfind("s_", 0) != 0) continue;
        std::istringstream is(line);
        std::string n; double tr, ss, sm, mcs, rf; is >> n >> tr >> ss >> sm >> mcs >> rf;
        std::string w, last; while (is >> w) last = w;
        o.v[n] = last; o.bias[n] = rf;
    }
    return o;
}

// ---------------------------------------------------------------- NEW: compact scale summary (flex_scales_v2)
// lambda scan of the verification (vflex_lamscan.json): the weight of the BB-lite template variance in the fit weights
struct LamScan { std::vector<double> lam; std::vector<std::vector<double>> th; std::vector<std::string> names; };
LamScan ReadLamScan() {
    LamScan L;
    std::ifstream in(Form("%s/vflex_lamscan.json", kVfl));
    std::stringstream ss; ss << in.rdbuf();
    const std::string s = ss.str();
    auto numsAfter = [&](size_t p, size_t stop) {
        std::vector<double> v;
        while (p < stop) {
            while (p < stop && !(std::isdigit((unsigned char)s[p]) || s[p] == '-' || s[p] == '.')) ++p;
            if (p >= stop) break;
            size_t e = p; while (e < stop && (std::isdigit((unsigned char)s[e]) || s[e] == '-' || s[e] == '.' || s[e] == 'e' || s[e] == 'E' || s[e] == '+')) ++e;
            v.push_back(std::atof(s.substr(p, e - p).c_str())); p = e;
        }
        return v;
    };
    size_t a = s.find("\"lams\"");
    if (a != std::string::npos) { size_t b = s.find('[', a), e = s.find(']', b); L.lam = numsAfter(b + 1, e); }
    a = s.find("\"th\"");
    if (a != std::string::npos) {
        size_t p = s.find('[', a) + 1;
        for (size_t k = 0; k < L.lam.size(); ++k) { size_t b = s.find('[', p), e = s.find(']', b); L.th.push_back(numsAfter(b + 1, e)); p = e + 1; }
    }
    a = s.find("\"names\"");
    if (a != std::string::npos) {
        size_t b = s.find('[', a), e = s.find(']', b), p = b;
        while (true) { size_t q1 = s.find('"', p + 1); if (q1 == std::string::npos || q1 > e) break; size_t q2 = s.find('"', q1 + 1); L.names.push_back(s.substr(q1 + 1, q2 - q1 - 1)); p = q2; }
    }
    return L;
}
double LamVal(const LamScan& L, const std::string& n, double lam) {
    int k = -1, j = -1;
    for (size_t i = 0; i < L.names.size(); ++i) if (L.names[i] == n) k = (int)i;
    for (size_t i = 0; i < L.lam.size(); ++i) if (std::fabs(L.lam[i] - lam) < 1e-6) j = (int)i;
    if (k < 0 || j < 0 || k >= (int)L.th[j].size()) return -999.;
    return L.th[j][k];
}

void DrawScalesSummary() {
    const FitInfo& F = gFit[gNominal];
    const LamScan L = ReadLamScan();
    // common eta scale of C2, C4, C5 (Wald combination with the F8 stat + MC covariance; vflex_ties.json)
    const double aEta = 0.4563875, eaEta = 0.0578898, pEta = 0.206;            // "wald C2,C4,C5"
    const double aEtc = 0.5706941, eaEtc = 0.0713602, pEtc = 0.507;            // "wald C2,C4,C5 (att.-corr.)"
    struct Row { std::string nm, lab, info; int col, cline; int grp; };
    // Fisher-information shares: flex_info.txt (reproduced by the verification to 0.1 %)
    const std::vector<Row> rows = {
        {"s_C2", "C2: bump, #pi^{0} Dalitz", "bump B 82 %, M4 < 0.455: 17 %", cC2, kOrange + 7, 0},
        {"s_C4", "C4: #eta#rightarrow#pi^{+}#pi^{-}#gamma, conversion", "peak P 75 %, 0.5225-0.534: 20 %", cC4, kGreen + 2, 0},
        {"s_C5", "C5: #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, direct", "peak P 96 %", cC5, kCyan + 2, 0},
        {"s_C3", "C3: bump, conversion", "B 74 %; needs M(ee) < 10 MeV", cC3, kMagenta + 1, 1},
        {"s_Geta", "G_{#eta}: other #eta decays", "CR-U 76 % (S1: 39 %)", cGeta, kViolet - 5, 2},
        {"s_GD_nom", "G_{D}: Dalitz pair, other pions", "CR-U 77 % (S5: 49 %)", cGDnom, kGray + 2, 2},
        {"s_GD_om", "G_{D}: Dalitz pair, #omega/#rho^{0} pions", "CR-U 98.5 % (S5: 76 %)", cGDom, kGray + 3, 2},
        {"s_GC", "G_{C}: conversion pair", "CR-U 85 % (S5: 45 %)", cGC, kBlack, 2},
        {"k_E2", "k_{E2}: CB scale (nuisance)", "constraint 1 #pm 0.05 (CB replicas)", cCB, cCBdark, 3}};
    const int n = (int)rows.size();
    const double xlo = -0.5, xhi = 3.2;
    TCanvas* c = new TCanvas("csc2", "", 1750, 1280);
    c->SetLeftMargin(0.215); c->SetRightMargin(0.285); c->SetBottomMargin(0.135); c->SetTopMargin(0.170);
    TH2D* fr = new TH2D("scfr2", "", 10, xlo, xhi, n, 0, n);
    for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, Lx(rows[i].lab).c_str());
    fr->GetYaxis()->SetLabelSize(0.024); fr->GetXaxis()->SetLabelSize(0.028); fr->GetXaxis()->SetTitleSize(0.027); fr->GetXaxis()->SetTitleOffset(1.10);
    fr->GetXaxis()->SetTitle(Lx("scale relative to SMASH as is (#sigma #times BR #times V3)").c_str()); fr->SetStats(0);
    fr->GetYaxis()->SetTickLength(0);
    fr->Draw();
    // common eta scale: spans the rows it was determined from (C2, C4, C5)
    const double yEtaLo = n - 3, yEtaHi = n;
    TBox* be = new TBox(aEta - eaEta, yEtaLo, aEta + eaEta, yEtaHi); be->SetFillColorAlpha(kYellow - 9, 0.85); be->Draw();
    TBox* bc = new TBox(aEtc - eaEtc, yEtaLo, aEtc + eaEtc, yEtaHi); bc->SetFillColor(kOrange - 3); bc->SetFillStyle(3354); bc->Draw();
    TLine* le = new TLine(aEta, yEtaLo, aEta, yEtaHi); le->SetLineColor(kOrange - 3); le->SetLineWidth(2); le->Draw();
    TLine* one = new TLine(1, 0, 1, n); one->SetLineStyle(2); one->SetLineColor(kGray + 3); one->SetLineWidth(2); one->Draw();
    TLine* zero = new TLine(0, 0, 0, n); zero->SetLineStyle(3); zero->SetLineColor(kGray + 2); zero->Draw();
    // group separators and headers
    const char* ghead[4] = {"#eta signal classes: one common #eta scale fits (p = 0.21)", "not determined with M(ee) #geq 10 MeV/c^{2}",
                            "CR-U background: 80-90 % correlated, split not determined", "CB normalisation"};
    int lastg = -1;
    for (int i = 0; i < n; ++i) {
        if (rows[i].grp != lastg) {
            const double y = n - i;
            if (i > 0) { TLine* l = new TLine(xlo, y, xhi, y); l->SetLineColor(kGray + 2); l->SetLineWidth(1); l->Draw(); }
            TLatex t; t.SetTextSize(0.0185); t.SetTextFont(52); t.SetTextAlign(33); t.SetTextColor(kGray + 3);
            Txt(t, xhi - 0.03, y - 0.04, ghead[rows[i].grp]);
            lastg = rows[i].grp;
        }
    }
    for (int i = 0; i < n; ++i) { TLine* l = new TLine(xlo, n - i, xhi, n - i); l->SetLineColor(kGray); l->SetLineStyle(3); if (i > 0 && rows[i].grp == rows[i - 1].grp) l->Draw(); }
    for (int i = 0; i < n; ++i) {
        const int k = Idx(F, rows[i].nm);
        if (k < 0) continue;
        const double v = F.v[k], es = F.es[k], et = F.et[k], em = F.em[k];
        const double y = n - i - 0.55;
        // fit-method box (lambda 0.5 - 1.5) and the full lambda 0 - 2 range (thin dotted)
        const double l05 = LamVal(L, rows[i].nm, 0.5), l15 = LamVal(L, rows[i].nm, 1.5), l0 = LamVal(L, rows[i].nm, 0.0), l2 = LamVal(L, rows[i].nm, 2.0);
        if (l05 > -900 && l15 > -900) {
            const double a = std::max(xlo, std::min(l05, l15)), b = std::min(xhi, std::max(l05, l15));
            TBox* bm = new TBox(a, y - 0.24, std::max(b, a + 0.012), y + 0.24); bm->SetFillColorAlpha(rows[i].col, 0.55); bm->SetLineColor(rows[i].cline); bm->Draw();
            TBox* bo = new TBox(a, y - 0.24, std::max(b, a + 0.012), y + 0.24); bo->SetFillStyle(0); bo->SetLineColor(rows[i].cline); bo->SetLineWidth(1); bo->Draw("L");
        }
        if (l0 > -900 && l2 > -900) {
            const double lo2 = std::min(l0, l2), hi2 = std::max(l0, l2);
            const double a = std::max(xlo, lo2), b = std::min(xhi, hi2), yd = y - 0.34;
            TLine* ld = new TLine(a, yd, b, yd); ld->SetLineColor(rows[i].cline); ld->SetLineStyle(3); ld->SetLineWidth(2); ld->Draw();
            if (lo2 >= xlo) { TLine* la = new TLine(a, yd - 0.06, a, yd + 0.06); la->SetLineColor(rows[i].cline); la->SetLineWidth(2); la->Draw(); }
            if (hi2 <= xhi) { TLine* lb = new TLine(b, yd - 0.06, b, yd + 0.06); lb->SetLineColor(rows[i].cline); lb->SetLineWidth(2); lb->Draw(); }
            else { TArrow* ar = new TArrow(xhi - 0.10, yd, xhi - 0.005, yd, 0.010, "|>"); ar->SetLineColor(rows[i].cline); ar->SetFillColor(rows[i].cline); ar->SetLineWidth(2); ar->Draw(); }
        }
        TLine* lt = new TLine(std::max(xlo, v - et), y, std::min(xhi, v + et), y); lt->SetLineColor(kBlack); lt->SetLineWidth(2); lt->Draw();
        for (double xe : {v - et, v + et}) if (xe > xlo && xe < xhi) { TLine* tk = new TLine(xe, y - 0.11, xe, y + 0.11); tk->SetLineWidth(2); tk->Draw(); }
        TLine* ls = new TLine(v - es, y, v + es, y); ls->SetLineColor(kBlack); ls->SetLineWidth(9); ls->Draw();
        TGraph* g = new TGraph(1); g->SetPoint(0, v, y); g->SetMarkerStyle(20); g->SetMarkerColor(rows[i].col); g->SetMarkerSize(2.3); g->Draw("P");
        TGraph* go = new TGraph(1); go->SetPoint(0, v, y); go->SetMarkerStyle(24); go->SetMarkerColor(kBlack); go->SetMarkerSize(2.3); go->Draw("P");
        // right column
        TLatex t; t.SetTextSize(0.0178); t.SetTextFont(42); t.SetTextAlign(12);
        Txt(t, xhi + 0.07, y + 0.27, Form("%.2f #pm %.2f (stat) #pm %.2f (MC)", v, es, em));
        if (l05 > -900) {
            const double l2v = LamVal(L, rows[i].nm, 2.0);
            if (l2v > xhi) Txt(t, xhi + 0.07, y + 0.0, Form("fit method #lambda 0.5-1.5: %.2f to %.2f (2: %.2f)", l05, l15, l2v));
            else Txt(t, xhi + 0.07, y + 0.0, Form("fit method #lambda 0.5-1.5: %.2f to %.2f", l05, l15));
        }
        Txt(t, xhi + 0.07, y - 0.27, "information: " + rows[i].info);
    }
    // legend
    TLegend* lg = new TLegend(0.005, 0.835, 0.995, 0.995); lg->SetNColumns(2); lg->SetTextSize(0.0185); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetMargin(0.12);
    TGraph* gp = new TGraph(1); gp->SetMarkerStyle(20); gp->SetMarkerSize(1.8); gp->SetMarkerColor(kGray + 1); gp->SetLineColor(kBlack); gp->SetLineWidth(9);
    TGraph* gt = new TGraph(1); gt->SetLineColor(kBlack); gt->SetLineWidth(2);
    lg->AddEntry(gp, Lx(Form("nominal fit %s (colour of the channel); thick bar: stat", gNominal.c_str())).c_str(), "PL");
    lg->AddEntry(LegBox(kYellow - 9), Lx(Form("common #eta scale of C2, C4, C5: %.2f #pm %.2f (Wald, p = %.2f)", aEta, eaEta, pEta)).c_str(), "F");
    lg->AddEntry(gt, Lx("thin bar: stat #oplus MC (MC = template noise of the simulation)").c_str(), "L");
    lg->AddEntry(LegBox(kOrange - 3, 3354, kOrange - 3), Lx(Form("the same after the template-noise correction: %.2f #pm %.2f (p = %.2f)", aEtc, eaEtc, pEtc)).c_str(), "F");
    lg->AddEntry(LegBox(kGray + 1, 1001, kBlack), Lx("box: fit method, template-noise weight #lambda = 0.5-1.5 (verification)").c_str(), "F");
    TGraph* gl = new TGraph(1); gl->SetLineColor(kGray + 3); gl->SetLineStyle(2); gl->SetLineWidth(2);
    lg->AddEntry(gl, Lx("1 = SMASH as is").c_str(), "L");
    TGraph* gd = new TGraph(1); gd->SetLineColor(kBlack); gd->SetLineStyle(3); gd->SetLineWidth(2);
    lg->AddEntry(gd, Lx("dotted, below the bar: full #lambda range 0-2 (#lambda = 0: Poisson weights)").c_str(), "L");
    lg->Draw();
    TLatex tn; tn.SetNDC(); tn.SetTextSize(0.0172); tn.SetTextFont(42);
    Txt(tn, 0.01, 0.040, Form("M(e^{+}e^{-}) #geq 10 MeV/c^{2}, S0, isBest #geq 0, M4 0.40-0.70 GeV/c^{2} (350 bins), model %s. The scales describe the data in this range only: they are a documented", gNominal.c_str()));
    Txt(tn, 0.01, 0.017, "description, not measurements (D36), and cannot be carried above M4 = 0.70. The flexible model is still incomplete (p < 0.001). Information = share of the Fisher information.");
    Save(c, "flex_scales");
}

// ---------------------------------------------------------------- old comparison figure (variants, Part III), gaps fixed
void DrawScalesCmp() {
    const FitInfo& F = gFit[gNominal];
    const FitInfo& A = gFit[gAlt];
    const FitInfo& E = gFit["V_E1k"];
    const FitInfo& Es = gFit["V_E2split"];
    Verd vb = ReadBias(gNominal);
    std::vector<std::string> rows = {"s_C2", "s_C3", "s_C4", "s_C5", "s_Geta", "s_GD_om", "s_GD_nom", "s_GC", "k_E2"};
    std::map<std::string, std::string> lab = {{"s_C2", "C2: #pi^{0} Dalitz (bump)"}, {"s_C3", "C3: conversion (bump)"}, {"s_C4", "C4: #eta#rightarrow#pi^{+}#pi^{-}#gamma conv."},
                                              {"s_C5", "C5: #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}"}, {"s_Geta", "G_{#eta}: other #eta"}, {"s_GD_om", "G_{D}: #omega/#rho^{0} pions"},
                                              {"s_GD_nom", "G_{D}: other pions"}, {"s_GC", "G_{C}: conversion pair"}, {"k_E2", "k_{E2}: CB scale"}};
    std::map<std::string, std::string> p3 = {{"s_C2", "C2"}, {"s_C3", "C3"}, {"s_C4", "C4"}, {"s_C5", "C5"}, {"s_Geta", "Geta"}, {"s_GD_om", "GD"}, {"s_GD_nom", "GD"}, {"s_GC", "GC"}, {"k_E2", "k_E2"}};
    const int cA = kBlue + 1, cS = kGreen + 3, cE = kOrange + 3, cP = kRed + 1;
    TCanvas* c = new TCanvas("csc_cmp2", "", 1700, 1150);
    c->SetLeftMargin(0.22); c->SetRightMargin(0.27); c->SetBottomMargin(0.12); c->SetTopMargin(0.155);
    const int n = (int)rows.size();
    TH2D* fr = new TH2D("scfr_cmp2", "", 10, -0.6, 3.4, n, 0, n);
    for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, Lx(lab[rows[i]]).c_str());
    fr->GetYaxis()->SetLabelSize(0.027); fr->GetXaxis()->SetLabelSize(0.027); fr->GetXaxis()->SetTitleSize(0.028); fr->GetXaxis()->SetTitleOffset(1.25);
    fr->GetXaxis()->SetTitle(Lx("scale relative to SMASH (#sigma #times PDG BR #times V3)").c_str()); fr->SetStats(0);
    fr->Draw();
    TLine* one = new TLine(1, 0, 1, n); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->SetLineWidth(2); one->Draw();
    TLine* zero = new TLine(0, 0, 0, n); zero->SetLineStyle(3); zero->SetLineColor(kGray + 1); zero->Draw();
    auto pt = [&](const FitInfo& f, const std::string& nm, double y, int col, int mk, bool thick, const std::string& alt = "") {
        int k = Idx(f, nm);
        if (k < 0 && !alt.empty()) k = Idx(f, alt);
        if (k < 0) return;
        const double v = f.v[k], et = f.et[k], es = f.es[k];
        TLine* lt = new TLine(std::max(-0.6, v - et), y, std::min(3.4, v + et), y); lt->SetLineColor(col); lt->SetLineWidth(2); lt->Draw();
        if (thick) { TLine* ls = new TLine(v - es, y, v + es, y); ls->SetLineColor(col); ls->SetLineWidth(7); ls->Draw(); }
        TGraph* g = new TGraph(1); g->SetPoint(0, std::max(-0.58, std::min(3.38, v)), y); g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetMarkerSize(1.6); g->Draw("P");
    };
    for (int i = 0; i < n; ++i) {
        const double y0 = n - i - 0.5;
        pt(F, rows[i], y0 + 0.27, kBlack, 20, true);
        pt(A, rows[i], y0 + 0.09, cA, 24, false, rows[i] == "s_GD_om" || rows[i] == "s_GD_nom" ? "s_GD" : "");
        pt(Es, rows[i], y0 - 0.09, cS, 26, false);
        pt(E, rows[i], y0 - 0.27, cE, 32, false);
        if (gP3.count(p3[rows[i]])) {
            const double v = gP3[p3[rows[i]]], e = gP3e.count(p3[rows[i]]) ? gP3e[p3[rows[i]]] : 0.;
            const double y = y0 - 0.40;
            TLine* l = new TLine(v - e, y, v + e, y); l->SetLineColor(cP); l->SetLineWidth(2); l->Draw();
            TGraph* g = new TGraph(1); g->SetPoint(0, v, y); g->SetMarkerStyle(25); g->SetMarkerColor(cP); g->SetMarkerSize(1.4); g->Draw("P");
        }
        const int k = Idx(F, rows[i]);
        TLatex t; t.SetTextSize(0.0185); t.SetTextFont(42); t.SetTextAlign(12);
        if (k >= 0) {
            const double rel = F.v[k] != 0 ? F.et[k] / std::fabs(F.v[k]) : 99;
            std::string vv = vb.v.count(rows[i]) ? vb.v[rows[i]] : (rows[i] == "k_E2" ? "nuisance" : "");
            Txt(t, 3.47, y0 + 0.27, Form("%.2f #pm %.2f #pm %.2f", F.v[k], F.es[k], F.em[k]));
            Txt(t, 3.47, y0 + 0.00, Form("%s (#sigma/s = %.0f %%)", rel < 0.3 ? "constrained" : rel < 0.5 ? "weakly constrained" : "not constrained", 100 * rel));
            if (vb.bias.count(rows[i])) {
                if (std::fabs(F.v[k]) > F.et[k]) Txt(t, 3.47, y0 - 0.27, Form("bias test %s, %+.0f %%", vv.c_str(), vb.bias[rows[i]]));
                else Txt(t, 3.47, y0 - 0.27, Form("bias test %s (truth #approx 0)", vv.c_str()));
            }
        }
    }
    for (int i = 1; i < n; ++i) { TLine* l = new TLine(-0.6, i, 3.4, i); l->SetLineColor(kGray); l->SetLineStyle(3); l->Draw(); }
    TLegend* lg = new TLegend(0.01, 0.855, 0.99, 0.995); lg->SetNColumns(2); lg->SetTextSize(0.022); lg->SetBorderSize(0); lg->SetFillStyle(0);
    auto mkg = [](int mk, int col) { TGraph* g = new TGraph(1); g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2); return g; };
    lg->AddEntry(mkg(20, kBlack), Lx(Form("nominal %s: M(ee) #geq 10 MeV (thick: stat, thin: stat #oplus MC)", gNominal.c_str())).c_str(), "PL");
    lg->AddEntry(mkg(26, cS), Lx("same model, M(ee) 10-20 and #geq 20 MeV fitted separately").c_str(), "PL");
    lg->AddEntry(mkg(24, cA), Lx(Form("alternative %s: per-slice tilts instead of the #omega/#rho^{0} split", gAlt.c_str())).c_str(), "PL");
    lg->AddEntry(mkg(32, cE), Lx("same model + E1 (M(ee) < 10 MeV) with a free CB-like k_{E1}").c_str(), "PL");
    lg->AddEntry(mkg(25, cP), Lx("Part III F-main: a_{#eta}, a#kappa, a#rho_{d}, ag; c_{j}(1+bx), c_{j}#kappa_{B} (weighted)").c_str(), "PL");
    lg->Draw();
    TLatex tn; tn.SetNDC(); tn.SetTextSize(0.019); tn.SetTextFont(42);
    Txt(tn, 0.01, 0.025, "S0, isBest #geq 0, M4 0.40-0.70, 5 M(#pi^{+}#pi^{-}) slices; 1 = SMASH as is. Right: nominal value #pm stat #pm MC; bias test = pseudo-data from the full simulation");
    Txt(tn, 0.01, 0.005, "at the fitted scales, fitted with bootstrap-replica templates (realistic template noise); PASS = |bias| < max(0.5 #sigma_{stat}, 10 %) and 68 % coverage in 0.60-0.76.");
    Save(c, "flex_scales_cmp");
}

// ---------------------------------------------------------------- correlation matrix
std::string Nice(const std::string& n) {
    static const std::map<std::string, std::string> m = {{"s_C2", "C2"}, {"s_C3", "C3"}, {"s_C4", "C4"}, {"s_C5", "C5"}, {"s_Geta", "G_{#eta}"},
                                                         {"s_GD", "G_{D}"}, {"s_GD_om", "G_{D} #omega/#rho"}, {"s_GD_nom", "G_{D} other"}, {"s_GC", "G_{C}"},
                                                         {"k_E2", "k_{E2}"}, {"k_E1", "k_{E1}"}, {"s_T3", "T3"}};
    auto it = m.find(n);
    if (it != m.end()) return it->second;
    if (n.rfind("b_", 0) == 0) return "b_{" + n.substr(2) + "}";
    return n;
}
void DrawCorr(const std::string& tag, const std::string& base) {
    const FitInfo& F = gFit[tag];
    const int n = (int)F.names.size();
    TCanvas* c = new TCanvas(("cc2" + tag).c_str(), "", 1100, 950);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.13); c->SetTopMargin(0.08);
    TH2D* h = new TH2D(("cor2" + tag).c_str(), "", n, 0, n, n, 0, n);
    for (int a = 0; a < n; ++a) {
        h->GetXaxis()->SetBinLabel(a + 1, Lx(Nice(F.names[a])).c_str());
        h->GetYaxis()->SetBinLabel(n - a, Lx(Nice(F.names[a])).c_str());
        for (int b = 0; b < n; ++b) h->SetBinContent(a + 1, n - b, F.cor[a][b]);
    }
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetPalette(kTemperatureMap);
    h->SetMinimum(-1); h->SetMaximum(1); h->SetStats(0);
    h->GetXaxis()->SetLabelSize(0.038); h->GetYaxis()->SetLabelSize(0.038); h->SetMarkerSize(1.3);
    h->GetXaxis()->LabelsOption("v");
    h->Draw("COLZ TEXT");
    TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextFont(42);
    Txt(t, 0.15, 0.94, Form("Correlations (stat #oplus MC) of the flexible fit %s, F10 E2, M4 0.40-0.70", tag.c_str()));
    Save(c, base);
}

// ---------------------------------------------------------------- bias test
void DrawBias(const std::string& tag) {
    std::ifstream in(Form("%s/flex_bias_%s.txt", kOut, tag.c_str()));
    std::string line; int scen = -1;
    std::vector<std::string> names[2];
    std::vector<double> bsig[2], cov[2], brel[2], msig[2];
    while (std::getline(in, line)) {
        if (line.rfind("# scenario 0", 0) == 0) scen = 0;
        else if (line.rfind("# scenario 1", 0) == 0) scen = 1;
        if (scen < 0 || line.rfind("s_", 0) != 0) continue;
        std::istringstream is(line);
        std::string n; double tr, ss, sm, mcs, rf, rfs, pm, pw, cv, mr, ms, mw;
        is >> n >> tr >> ss >> sm >> mcs >> rf >> rfs >> pm >> pw >> cv >> mr >> ms >> mw;
        names[scen].push_back(n); bsig[scen].push_back(rfs); cov[scen].push_back(cv); brel[scen].push_back(rf); msig[scen].push_back(ms);
    }
    const int n = (int)names[0].size();
    if (n == 0) return;
    TCanvas* c = new TCanvas(("cb2" + tag).c_str(), "", 1700, 760);
    c->Divide(2, 1);
    TGraph* gl0 = nullptr; TGraph* gl1 = nullptr;
    for (int pnl = 0; pnl < 2; ++pnl) {
        c->cd(pnl + 1);
        gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.15); gPad->SetRightMargin(0.04);
        TH2D* fr = new TH2D(Form("bfr2%d", pnl), "", 10, pnl == 0 ? -80 : 0, pnl == 0 ? 40 : 1, n, 0, n);
        for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, Lx(Nice(names[0][i])).c_str());
        fr->GetXaxis()->SetTitle(Lx(pnl == 0 ? "bias of the realistic fit [% of the truth]" : "68 % coverage of stat #oplus MC").c_str());
        fr->GetXaxis()->SetTitleSize(0.050); fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.060); fr->SetStats(0);
        fr->Draw();
        TBox* bx = pnl == 0 ? new TBox(-10, 0, 10, n) : new TBox(0.60, 0, 0.76, n);
        bx->SetFillColorAlpha(kGreen - 9, 0.5); bx->Draw();
        fr->Draw("AXIS SAME");
        for (int s = 0; s < 2; ++s)
            for (int i = 0; i < (int)names[s].size(); ++i) {
                double x = pnl == 0 ? brel[s][i] : cov[s][i];
                if (pnl == 0 && names[s][i] == "s_C3" && s == 0) x = -999;   // truth -0.09: relative bias meaningless (in sigma, right)
                const double y = n - i - 0.5 + (s == 0 ? 0.15 : -0.15);
                TGraph* g = new TGraph(1);
                g->SetPoint(0, pnl == 0 ? std::max(-79., std::min(39., x)) : x, y);
                g->SetMarkerStyle(s == 0 ? 20 : 24); g->SetMarkerSize(1.6); g->SetMarkerColor(s == 0 ? kBlack : kBlue + 1);
                if (x > -900) g->Draw("P");
                if (s == 0) gl0 = g; else gl1 = g;
                if (pnl == 0 && s == 0) { TLatex t; t.SetTextSize(0.034); t.SetTextFont(42); t.SetTextAlign(12); Txt(t, 16, y, Form("%+.1f #sigma_{tot}", bsig[s][i])); }
            }
        if (pnl == 1) { TLatex t; t.SetNDC(); t.SetTextSize(0.040); t.SetTextFont(42); Txt(t, 0.17, 0.90, Form("500 pseudo-experiments per truth, model %s", tag.c_str())); }
    }
    c->cd(0);
    TLegend* lg = new TLegend(0.08, 0.90, 0.50, 0.995); lg->SetNColumns(2); lg->SetTextSize(0.034); lg->SetBorderSize(0); lg->SetFillStyle(0);
    TGraph* p0 = new TGraph(1); p0->SetMarkerStyle(20); p0->SetMarkerSize(1.6); p0->SetMarkerColor(kBlack);
    TGraph* p1 = new TGraph(1); p1->SetMarkerStyle(24); p1->SetMarkerSize(1.6); p1->SetMarkerColor(kBlue + 1);
    lg->AddEntry(p0, "truth = flexible fit of the data", "P");
    lg->AddEntry(p1, "truth = SMASH as is", "P");
    lg->Draw();
    (void)gl0; (void)gl1;
    c->cd(0);
    Save(c, "flex_bias");
}

// ---------------------------------------------------------------- variants (stability)
void DrawVariants() {
    const std::vector<std::pair<std::string, std::string>> V = {
        {gNominal, "nominal " + gNominal}, {gAlt, "alt. " + gAlt + " (tilts)"}, {"F7", "F7 (7 scales)"}, {"F9", "F9 (G_{C} split)"},
        {"V_E2split", "E2 split 10-20/#geq20"}, {"V_E1k", "+E1, k_{E1} free"}, {"V_E1T3", "+E1, T3 free"}, {"V_E1k3", "+E1, E2 split"},
        {"V_r42", "M4 0.42-0.70"}, {"V_r45", "M4 0.45-0.70"}, {"V_rawBg", "raw G_{D}, G_{C}"}, {"V_T3all", "+T3 in E2"}, {"V_r100", "M4 0.40-1.00"}, {"V_r100raw", "0.40-1.00, raw bkg"}};
    const std::vector<std::string> sc = {"s_C2", "s_C3", "s_C4", "s_C5", "s_Geta", "s_GD_om", "s_GD_nom", "s_GC", "k_E2"};
    const std::vector<std::string> ttl = {"C2 bump, #pi^{0} Dalitz", "C3 bump, conversion", "C4 #eta#rightarrow#pi^{+}#pi^{-}#gamma", "C5 #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}",
                                          "G_{#eta} other #eta", "G_{D} #omega/#rho^{0} pions", "G_{D} other pions", "G_{C} conversion pair", "k_{E2} CB scale"};
    const int tcol[9] = {kOrange + 7, kMagenta + 1, kGreen + 2, kCyan + 2, kViolet - 5, kGray + 2, kGray + 2, kGray + 3, kRed + 1};
    const double xlo[9] = {0, -0.6, -0.8, 0, 0, -1.0, 0.6, -3, 0.3}, xhi[9] = {0.8, 0.8, 1.4, 1.4, 6.5, 3.0, 1.8, 3.5, 1.6};
    TCanvas* c = new TCanvas("cvar2", "", 2000, 1250);
    c->Divide(5, 2, 0.002, 0.004);
    const int nv = (int)V.size();
    for (int q = 0; q < (int)sc.size(); ++q) {
        c->cd(q + 1);
        gPad->SetLeftMargin(q % 5 == 0 ? 0.42 : 0.06); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.10);
        TH2D* fr = new TH2D(Form("vfr2%d", q), "", 10, xlo[q], xhi[q], nv, 0, nv);
        for (int i = 0; i < nv; ++i) fr->GetYaxis()->SetBinLabel(nv - i, q % 5 == 0 ? Lx(V[i].second).c_str() : "");
        fr->GetYaxis()->SetLabelSize(0.075); fr->GetXaxis()->SetLabelSize(0.06); fr->SetStats(0); fr->GetXaxis()->SetNdivisions(505);
        fr->Draw();
        TLine* one = new TLine(1, 0, 1, nv); one->SetLineStyle(2); one->SetLineColor(kGray + 2); if (1 > xlo[q] && 1 < xhi[q]) one->Draw();
        const FitInfo& N = gFit[gNominal];
        const int kn = Idx(N, sc[q]);
        if (kn >= 0) { TBox* b = new TBox(std::max(xlo[q], N.v[kn] - N.et[kn]), 0, std::min(xhi[q], N.v[kn] + N.et[kn]), nv); b->SetFillColorAlpha(kGray, 0.4); b->Draw(); }
        TBox* bad = new TBox(xlo[q], 0, xhi[q], 2); bad->SetFillColorAlpha(kGray + 1, 0.18); bad->Draw();   // 0.40-1.00 rows
        for (int i = 0; i < nv; ++i) {
            if (!gFit.count(V[i].first)) continue;
            const FitInfo& f = gFit[V[i].first];
            int k = Idx(f, sc[q]);
            if (k < 0 && (sc[q] == "s_GD_om" || sc[q] == "s_GD_nom")) k = Idx(f, "s_GD");
            if (k < 0) continue;
            const double y = nv - i - 0.5;
            const double v = f.v[k], et = f.et[k];
            const bool r100 = V[i].first.find("r100") != std::string::npos;
            const int col = i == 0 ? kBlack : r100 ? kGray + 2 : kBlue + 1;
            TLine* l = new TLine(std::max(xlo[q], v - et), y, std::min(xhi[q], v + et), y); l->SetLineWidth(2); l->SetLineColor(col); l->Draw();
            TGraph* g = new TGraph(1); g->SetPoint(0, std::max(xlo[q], std::min(xhi[q], v)), y); g->SetMarkerStyle(i == 0 ? 20 : r100 ? 25 : 21); g->SetMarkerSize(1.2); g->SetMarkerColor(col); g->Draw("P");
        }
        TLatex t; t.SetNDC(); t.SetTextSize(0.072); t.SetTextFont(42); t.SetTextColor(tcol[q]);
        Txt(t, q % 5 == 0 ? 0.43 : 0.08, 0.935, ttl[q]);
    }
    c->cd(10);
    TLatex t; t.SetNDC(); t.SetTextSize(0.062); t.SetTextFont(42);
    Txt(t, 0.04, 0.84, "Every scale in every variant");
    Txt(t, 0.04, 0.75, "(value #pm stat #oplus MC).");
    Txt(t, 0.04, 0.64, "Grey band: nominal fit.");
    Txt(t, 0.04, 0.55, "Dashed: SMASH (1).");
    Txt(t, 0.04, 0.44, "Variants use the nominal model");
    Txt(t, 0.04, 0.36, "(no tilt, G_{D} split), except");
    Txt(t, 0.04, 0.28, "alt., F7, F9.");
    Txt(t, 0.04, 0.17, "Grey open squares, shaded rows:");
    Txt(t, 0.04, 0.09, "M4 0.40-1.00, not a usable fit.");
    Save(c, "flex_variants");
}

}  // namespace flexg2

void flex_figs_v2() {
    using namespace flexg2;
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetHatchesLineWidth(1); gStyle->SetHatchesSpacing(0.9);
    // PDF line widths are absolute while these canvases are 1300-2400 px wide: with the default scale (3) a width-2 line
    // of the fitted model became a 3 pt band that hid the data. Scale 1.5 keeps frames and error bars at about 0.75 pt.
    gStyle->SetLineScalePS(1.5);
    ReadResults();
    TFile* f = TFile::Open(Form("%s/flex_fit.root", kOut));
    TFile* fh = TFile::Open(kHistFile);
    const FitInfo& A = gFit["ASIS"];
    const FitInfo& N = gFit[gNominal];
    const FitInfo& Al = gFit[gAlt];
    const std::string toy = "pseudo-data from a perfect model give 356 #pm 33 (largest of 2000: 461), so p < 0.001.";
    DrawSpectra(f, "ASIS", "E2", 0.40, 0.70, "flex_spectra_asis",
                "M(e^{+}e^{-}) #geq 10 MeV/c^{2}: data vs SMASH as is (every component at the SMASH cross section, PDG branching ratios and V3; CB with k = 1)",
                {Form("S0, isBest #geq 0, L = 6079.2 nb^{-1}, raw rows; data errors from the event groups. No parameter is fitted. Hybrid GLS #chi^{2} %.0f / 350 bins, Pearson %.0f; panel #chi^{2}_{P} = #Sigma(N-#mu)^{2}/V.", A.chi2H, A.chi2P)},
                false, false);
    DrawSpectra(f, gNominal, "E2", 0.40, 0.70, "flex_spectra_fit",
                Form("M(e^{+}e^{-}) #geq 10 MeV/c^{2}: data vs SMASH after the flexible fit %s (one free scale per component, no tilt)", gNominal.c_str()),
                {Form("S0, isBest #geq 0, L = 6079.2 nb^{-1}, raw rows; data errors from the event groups. Hybrid GLS #chi^{2} %.0f / %d (Pearson %.0f); %s Panel #chi^{2}_{P} = #Sigma(N-#mu)^{2}/V.", N.chi2H, N.ndf, N.chi2P, toy.c_str()),
                 Form("The negative C3 term (s_{C3} = %.2f) is drawn below zero (hatched); the black line is the fitted model = top of the stack minus that band.", Val(gNominal, "s_C3", 0.))},
                false, true);
    DrawSpectra(f, gAlt, "E2", 0.40, 0.70, "flex_spectra_alt",
                Form("M(e^{+}e^{-}) #geq 10 MeV/c^{2}: data vs SMASH after the alternative flexible fit %s (7 scales plus per-slice tilts of the Dalitz and conversion backgrounds)", gAlt.c_str()),
                {Form("S0, isBest #geq 0, L = 6079.2 nb^{-1}, raw rows; data errors from the event groups. Hybrid GLS #chi^{2} %.0f / %d (Pearson %.0f). Panel #chi^{2}_{P} = #Sigma(N-#mu)^{2}/V.", Al.chi2H, Al.ndf, Al.chi2P),
                 Form("The negative C3 term (s_{C3} = %.2f) is drawn below zero (hatched); the black line is the fitted model = top of the stack minus that band.", Val(gAlt, "s_C3", 0.))},
                false, true);
    DrawMee(fh);
    DrawScalesSummary();
    DrawScalesCmp();
    DrawCorr(gNominal, "flex_corr");
    DrawCorr(gAlt, "flex_corr_alt");
    DrawBias(gNominal);
    DrawVariants();
    {
        const FitInfo& R = gFit["V_r100"];
        const FitInfo& Rr = gFit["V_r100raw"];
        DrawSpectra(f, "V_r100_asis", "E2", 0.40, 1.00, "flex_r100_asis", "M(e^{+}e^{-}) #geq 10 MeV/c^{2}, M4 0.40-1.00: SMASH as is (sum of the slices)",
                    {"S0, isBest #geq 0, raw rows; every component at #sigma_{SMASH} #times PDG BR #times V3, CB k = 1; background templates smoothed per piece (0.40-0.70, 0.70-1.00)."}, true, false);
        DrawSpectra(f, "V_r100raw", "E2", 0.40, 1.00, "flex_r100raw",
                    Form("M4 0.40-1.00 after the flexible fit with raw background templates (model %s, sum of the slices)", gNominal.c_str()),
                    {Form("Pearson #chi^{2}/ndf %.0f/%d, hybrid GLS %.0f; G_{C} %.2f, G_{#eta} %.1f, k_{E2} %.2f (far from its constraint 1 #pm 0.05).", Rr.chi2P, Rr.ndf, Rr.chi2H,
                          Val("V_r100raw", "s_GC", 0.), Val("V_r100raw", "s_Geta", 0.), Val("V_r100raw", "k_E2", 0.)),
                     "Scaled SMASH components cannot describe 0.70-1.00; the scales of 0.40-0.70 must not be carried above 0.70."},
                    true, true);
        DrawSpectra(f, "V_r100", "E2", 0.40, 1.00, "flex_r100",
                    Form("M4 0.40-1.00 after the flexible fit (model %s, background smoothed per piece, sum of the slices)", gNominal.c_str()),
                    {Form("Pearson #chi^{2}/ndf %.0f/%d, hybrid GLS %.0f; s_{GC} = %.2f (unphysical) is drawn below zero (hatched), the black line is the fitted model", R.chi2P, R.ndf, R.chi2H, Val("V_r100", "s_GC", 0.)),
                     "(top of the stack minus that band). Step at 0.70: the two smoothing pieces meet. The components cannot describe 0.70-1.00."},
                    true, true);
    }
    printf("flex_figs_v2: done\n");
}
