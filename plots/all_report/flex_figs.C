// ============================================================================
// flex_figs.C  --  package flex of report_all: figures of the flexible SMASH description
// ============================================================================
// Reads /hdd3/.../flex/flex_fit.root, flex_results.txt, flex_bias_<tag>.txt and hists_v1.root (M(ee) projections, BEFORE
// layouts B_mee_R<k>), writes docs/latex/figs_all/flex_*.pdf/.png:
//   flex_spectra_asis   data vs SMASH as is (all scales 1), per M(pi+pi-) slice and summed, M(ee) >= 10 MeV
//   flex_spectra_fit    the same after the flexible fit (nominal model)
//   flex_spectra_alt    the same for the alternative model (per-slice tilts instead of the omega/rho0 split)
//   flex_mee            M(ee) projections in M4 < 0.455, the bump window B and CR-U: data vs SMASH as is and after the fit
//   flex_scales         all scales with errors against 1 (SMASH) and against the Part III values
//   flex_corr           correlation matrix (stat + MC) of the nominal fit
//   flex_bias           bias test on pseudo-data (realistic template noise): bias and 68 % coverage of every scale
//   flex_variants       every scale in every variant (stability)
//   flex_r100           the M4 range extended to 1.00: data vs SMASH as is and after the fit (summed over slices)
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_flex", kTRUE)' 'plots/all_report/flex_figs.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "defs.h"
#include "hists_io.h"

namespace flexg {
using namespace alldefs;
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/flex";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";

void Save(TCanvas* c, const std::string& base) {
    c->SaveAs((std::string(kFig) + "/" + base + ".pdf").c_str());
    c->SaveAs((std::string(kFig) + "/" + base + ".png").c_str());
}

// cumulative stack drawn from the top down (THStack with "SAME" crashes in ROOT 6.32; same helper as wp5_extract.C)
struct Stack {
    std::vector<TH1D*> h;
    void Add(TH1D* x) { h.push_back(x); }
    void Draw(const std::string& nm) {
        std::vector<TH1D*> cum;
        for (size_t q = 0; q < h.size(); ++q) {
            TH1D* c = (TH1D*)h[q]->Clone(Form("%s_cum%zu", nm.c_str(), q)); c->SetDirectory(nullptr);
            if (q > 0) c->Add(cum.back());
            c->SetFillColor(h[q]->GetFillColor()); c->SetLineColor(h[q]->GetLineColor()); c->SetLineWidth(1); c->SetFillStyle(1001);
            cum.push_back(c);
        }
        for (int q = (int)cum.size() - 1; q >= 0; --q) cum[q]->Draw("HIST SAME");
    }
};

// ---------------------------------------------------------------- results parser
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
    // Part III errors of GD, GC from p3equiv_err (C3C4 for C3 and C4)
    gP3e["C3"] = gP3e["C3C4"]; gP3e["C4"] = gP3e["C3C4"];
}
int Idx(const FitInfo& f, const std::string& n) { for (size_t k = 0; k < f.names.size(); ++k) if (f.names[k] == n) return (int)k; return -1; }

// ---------------------------------------------------------------- spectra (per slice and summed)
struct Part { std::string name, label; int col; };
std::vector<Part> PartsFor(TDirectory* d, const std::string& h) {
    std::vector<Part> v;
    auto has = [&](const std::string& n) { return d->Get((n + "_" + h).c_str()) != nullptr; };
    if (has("GC")) v.push_back({"GC", "G_{C} (non-#eta, conversion pair)", kViolet - 9});
    if (has("GC_nom")) v.push_back({"GC_nom", "G_{C}, other pions", kViolet - 9});
    if (has("GC_om")) v.push_back({"GC_om", "G_{C}, #omega/#rho^{0} pions", kViolet - 4});
    if (has("GD")) v.push_back({"GD", "G_{D} (non-#eta, Dalitz pair)", kGray});
    if (has("GD_nom")) v.push_back({"GD_nom", "G_{D}, other pions", kGray});
    if (has("GD_om")) v.push_back({"GD_om", "G_{D}, #omega/#rho^{0} pions", kGray + 2});
    if (has("T3")) v.push_back({"T3", "T3 (fake lepton)", kMagenta - 9});
    if (has("Geta")) v.push_back({"Geta", "G_{#eta} (other #eta)", kOrange - 3});
    v.push_back({"C34", "C3 + C4 (conversions)", kCyan - 7});
    v.push_back({"C2", "C2 (#pi^{0} Dalitz, bump)", kAzure + 1});
    v.push_back({"C5", "C5 (#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-})", kGreen + 2});
    return v;
}
TH1D* Get(TDirectory* d, const std::string& n) {
    TH1D* h = (TH1D*)d->Get(n.c_str());
    if (!h) return nullptr;
    TH1D* c = (TH1D*)h->Clone((n + Form("_c%p", (void*)d)).c_str());
    c->SetDirectory(nullptr);
    return c;
}
TH1D* GetPart(TDirectory* d, const std::string& p, const std::string& h) {
    if (p == "C34") { TH1D* a = Get(d, "C3_" + h); TH1D* b = Get(d, "C4_" + h); if (a && b) a->Add(b); return a; }
    return Get(d, p + "_" + h);
}
// per MeV/c2 and restricted to [lo, hi)
TH1D* PerMeV(TH1D* h, double lo, double hi) {
    std::vector<double> e;
    for (int b = 1; b <= h->GetNbinsX() + 1; ++b) { const double x = h->GetBinLowEdge(b); if (x >= lo - 1e-9 && x <= hi + 1e-9) e.push_back(x); }
    TH1D* o = new TH1D(Form("%s_pm", h->GetName()), "", (int)e.size() - 1, e.data());
    o->SetDirectory(nullptr);
    for (int b = 1; b <= o->GetNbinsX(); ++b) {
        const int bb = h->FindBin(o->GetBinCenter(b));
        const double w = 1000. * h->GetBinWidth(bb);
        o->SetBinContent(b, h->GetBinContent(bb) / w); o->SetBinError(b, h->GetBinError(bb) / w);
    }
    return o;
}
void DrawSpectra(TFile* f, const std::string& dir, const std::string& cls, double lo, double hi, const std::string& base, const std::string& title,
                 const std::string& note, bool sumOnly = false) {
    TDirectory* d = (TDirectory*)f->Get(dir.c_str());
    if (!d) { printf("flex_figs: %s missing\n", dir.c_str()); return; }
    const int ncol = sumOnly ? 1 : 6;
    const double xLeg = sumOnly ? 0.70 : 0.865;   // panels left of xLeg, legend right of it
    TCanvas* c = new TCanvas(("c" + base).c_str(), "", sumOnly ? 1300 : 2400, sumOnly ? 820 : 900);
    c->cd();
    TPad* top = new TPad("top", "", 0, 0.945, 1, 1); top->Draw();
    top->cd(); TLatex tt; tt.SetNDC(); tt.SetTextSize(0.45); tt.SetTextFont(42); tt.DrawLatex(0.01, 0.30, title.c_str());
    std::vector<std::string> hs;
    for (int j = 1; j <= 5; ++j) hs.push_back(Form("S%d%s", j, cls.c_str()));
    const char* slTxt[6] = {"S1: M(#pi^{+}#pi^{-}) < 0.34", "S2: 0.34-0.37", "S3: 0.37-0.39", "S4: 0.39-0.43", "S5: #geq 0.43", "sum of the slices"};
    TH1D* legN = nullptr; TH1D* legCB = nullptr;
    std::vector<TH1D*> legP; std::vector<Part> legParts;
    for (int col = 0; col < ncol; ++col) {
        const int jcol = sumOnly ? 5 : col;
        c->cd();
        const double x0 = xLeg * col / ncol, x1 = xLeg * (col + 1) / ncol;
        TPad* p1 = new TPad(Form("p1_%d", col), "", x0, 0.36, x1, 0.945);
        TPad* p2 = new TPad(Form("p2_%d", col), "", x0, 0.045, x1, 0.36);
        p1->SetBottomMargin(0.015); p1->SetTopMargin(0.075); p1->SetLeftMargin(0.21); p1->SetRightMargin(0.03);
        p2->SetTopMargin(0.03); p2->SetBottomMargin(0.27); p2->SetLeftMargin(0.21); p2->SetRightMargin(0.03);
        p1->Draw(); p2->Draw();
        std::vector<std::string> use;
        if (jcol < 5) use.push_back(hs[jcol]); else use = hs;
        TH1D *hN = nullptr, *hCB = nullptr, *hMu = nullptr, *hEG = nullptr;
        std::vector<Part> parts = PartsFor(d, hs[0]);
        std::vector<TH1D*> hp(parts.size(), nullptr);
        for (const auto& h : use) {
            auto acc = [&](TH1D*& a, TH1D* b) { if (!b) return; if (!a) a = b; else { a->Add(b); delete b; } };
            acc(hN, Get(d, "N_" + h)); acc(hCB, Get(d, "kCB_" + h)); acc(hMu, Get(d, "mu_" + h)); acc(hEG, Get(d, "vEG_" + h));
            for (size_t q = 0; q < parts.size(); ++q) acc(hp[q], GetPart(d, parts[q].name, h));
        }
        if (!hN) continue;
        for (int b = 1; b <= hN->GetNbinsX(); ++b) hN->SetBinError(b, std::sqrt(std::max(0., hEG->GetBinContent(b))));
        TH1D* N = PerMeV(hN, lo, hi);
        TH1D* CB = PerMeV(hCB, lo, hi);
        TH1D* MU = PerMeV(hMu, lo, hi);
        Stack* st = new Stack();
        CB->SetFillColor(kRed - 7); CB->SetLineColor(kRed - 7); st->Add(CB);
        std::vector<TH1D*> pp;
        for (size_t q = 0; q < parts.size(); ++q) {
            if (!hp[q]) { pp.push_back(nullptr); continue; }
            TH1D* x = PerMeV(hp[q], lo, hi);
            x->SetFillColor(parts[q].col); x->SetLineColor(parts[q].col); x->SetLineWidth(1);
            st->Add(x);
            pp.push_back(x);
        }
        p1->cd();
        double ymax = 0;
        for (int b = 1; b <= N->GetNbinsX(); ++b) ymax = std::max(ymax, std::max(N->GetBinContent(b) + N->GetBinError(b), MU->GetBinContent(b)));
        TH1D* fr = (TH1D*)N->Clone(Form("fr%d", col)); fr->Reset(); fr->SetMaximum(1.12 * ymax); fr->SetMinimum(0);
        fr->GetYaxis()->SetTitle("rows per MeV/c^{2}"); fr->GetYaxis()->SetTitleSize(sumOnly ? 0.050 : 0.062); fr->GetYaxis()->SetLabelSize(sumOnly ? 0.045 : 0.056);
        fr->GetYaxis()->SetTitleOffset(sumOnly ? 1.9 : 1.6);
        fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetMaxDigits(4);
        fr->Draw("AXIS");
        st->Draw(Form("%s_%d", base.c_str(), col));
        N->SetMarkerStyle(20); N->SetMarkerSize(sumOnly ? 0.9 : 0.75); N->SetMarkerColor(kBlack); N->SetLineColor(kBlack); N->Draw("E1 SAME");
        fr->Draw("AXIS SAME");
        TLatex lt; lt.SetNDC(); lt.SetTextSize(sumOnly ? 0.050 : 0.062); lt.SetTextFont(42); lt.DrawLatex(0.23, 0.945, slTxt[jcol]);
        legN = N; legCB = CB; legP = pp; legParts = parts;
        // ratio
        p2->cd();
        TH1D* R = (TH1D*)N->Clone(Form("R%d", col));
        R->Divide(MU);
        for (int b = 1; b <= R->GetNbinsX(); ++b) R->SetBinError(b, MU->GetBinContent(b) > 0 ? N->GetBinError(b) / MU->GetBinContent(b) : 0.);
        R->SetMinimum(0.55); R->SetMaximum(1.45);
        const double rs = sumOnly ? 0.075 : 0.088;
        R->GetYaxis()->SetTitle("data / model"); R->GetYaxis()->SetTitleSize(rs * 1.1); R->GetYaxis()->SetLabelSize(rs); R->GetYaxis()->SetTitleOffset(sumOnly ? 1.25 : 1.05); R->GetYaxis()->SetNdivisions(505);
        R->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); R->GetXaxis()->SetTitleSize(rs * 1.05); R->GetXaxis()->SetLabelSize(rs); R->GetXaxis()->SetTitleOffset(1.15);
        R->GetXaxis()->SetNdivisions(505);
        R->SetMarkerStyle(20); R->SetMarkerSize(0.6);
        R->Draw("E1");
        TLine* l1 = new TLine(lo, 1, hi, 1); l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
        // chi2 of the panel (Pearson with the model variance and the event-group excess: information)
        double ch = 0; int nb = 0;
        for (int b = 1; b <= hN->GetNbinsX(); ++b) {
            const double x = hN->GetBinCenter(b);
            if (x < lo || x > hi || hMu->GetBinContent(b) <= 0) continue;
            const double v = hMu->GetBinError(b) * hMu->GetBinError(b) + std::max(0., hEG->GetBinContent(b) - hN->GetBinContent(b));
            ch += std::pow(hN->GetBinContent(b) - hMu->GetBinContent(b), 2) / v; ++nb;
        }
        TLatex l2; l2.SetNDC(); l2.SetTextSize(rs); l2.SetTextFont(42); l2.DrawLatex(0.25, 0.86, Form("#chi^{2}_{P} %.0f / %d bins", ch, nb));
    }
    // legend column
    c->cd();
    TPad* pl = new TPad("pleg", "", xLeg, 0.36, 1, 0.945); pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.02, 0.10, 0.99, 0.92); lg->SetTextSize(sumOnly ? 0.050 : 0.062); lg->SetFillStyle(0); lg->SetBorderSize(0);
    if (legN) {
        lg->AddEntry(legN, "data N (S0, isBest #geq 0)", "PE");
        for (int q = (int)legParts.size() - 1; q >= 0; --q) if (legP[q]) lg->AddEntry(legP[q], legParts[q].label.c_str(), "F");
        lg->AddEntry(legCB, "k CB, CB = 2#sqrt{N_{++}N_{--}}", "F");
    }
    lg->Draw();
    c->cd();
    if (!note.empty()) { TLatex ln; ln.SetNDC(); ln.SetTextSize(sumOnly ? 0.022 : 0.020); ln.SetTextFont(42); ln.DrawLatex(0.01, 0.012, note.c_str()); }
    Save(c, base);
}

// ---------------------------------------------------------------- M(ee) projections (hists_v1 BEFORE layouts)
void DrawMee(TFile* fh) {
    const FitInfo& F = gFit[gNominal];
    auto sc = [&](const std::string& n) { const int k = Idx(F, n); return k >= 0 ? F.v[k] : -999.; };
    struct C { std::string name, label; std::vector<int> parts; int col; double s; };
    std::vector<C> comps;
    const double sGD = sc("s_GD"), sGDom = sc("s_GD_om"), sGDnom = sc("s_GD_nom");
    comps.push_back({"GC", "G_{C}", {kGconv}, kViolet - 9, sc("s_GC")});
    if (sGD > -900) comps.push_back({"GD", "G_{D}", {kGD}, kGray, sGD});
    else {
        comps.push_back({"GD_nom", "G_{D}, other pions", {kGpi0Dnom, kGvecnom, kGothnom}, kGray, sGDnom});
        comps.push_back({"GD_om", "G_{D}, #omega/#rho^{0} pions", {kGpi0Dom, kGvecom, kGothom}, kGray + 2, sGDom});
    }
    comps.push_back({"Geta", "G_{#eta}", {kGeta}, kOrange - 3, sc("s_Geta")});
    comps.push_back({"C3", "C3", {kC3}, kCyan - 7, sc("s_C3")});
    comps.push_back({"C4", "C4", {kC4}, kCyan + 2, sc("s_C4")});
    comps.push_back({"C2", "C2", {kC2}, kAzure + 1, sc("s_C2")});
    comps.push_back({"C5", "C5", {kC5}, kGreen + 2, sc("s_C5")});
    const double kE2 = sc("k_E2");
    const int regs[3] = {0, 1, 5};
    const char* rtxt[3] = {"M4 < 0.455 GeV/c^{2}", "bump window B: 0.455-0.5225", "CR-U: 0.565-0.700"};
    TCanvas* c = new TCanvas("cmee", "", 2300, 1000);
    c->cd();
    TPad* top = new TPad("tmee", "", 0, 0.93, 1, 1); top->Draw(); top->cd();
    TLatex tt; tt.SetNDC(); tt.SetTextSize(0.40); tt.SetTextFont(42);
    tt.DrawLatex(0.01, 0.30, Form("M(e^{+}e^{-}) projections (S0, all slices): data vs SMASH as is (dashed) and after the flexible fit %s (stacked); M(ee) < 10 MeV/c^{2} is not fitted", gNominal.c_str()));
    const double xmax = 0.16;
    for (int q = 0; q < 3; ++q) {
        const int r = regs[q];
        c->cd();
        TPad* p1 = new TPad(Form("m1_%d", q), "", 0.80 * q / 3., 0.36, 0.80 * (q + 1) / 3., 0.93);
        TPad* p2 = new TPad(Form("m2_%d", q), "", 0.80 * q / 3., 0.045, 0.80 * (q + 1) / 3., 0.36);
        p1->SetBottomMargin(0.015); p1->SetLeftMargin(0.17); p1->SetRightMargin(0.03); p1->SetTopMargin(0.07); p1->SetLogy();
        p2->SetTopMargin(0.03); p2->SetBottomMargin(0.27); p2->SetLeftMargin(0.17); p2->SetRightMargin(0.03);
        p1->Draw(); p2->Draw();
        const std::string lay = BVarName(kBmee, r);
        TH1* N = allio::Data(fh, "S0", "all", "N", lay);
        TH1* CB = allio::Data(fh, "S0", "all", "CB", lay);
        // model: CB (k_E2 above 10 MeV, 1 below) + sum s_c comp; as is: CB + sum comp
        TH1D* asis = (TH1D*)N->Clone(Form("asis%d", q)); asis->Reset();
        TH1D* fit = (TH1D*)N->Clone(Form("fit%d", q)); fit->Reset();
        TH1D* cbk = (TH1D*)N->Clone(Form("cbk%d", q)); cbk->Reset();
        for (int b = 1; b <= N->GetNbinsX(); ++b) {
            const double k = N->GetBinCenter(b) >= 0.010 ? kE2 : 1.0;
            cbk->SetBinContent(b, k * CB->GetBinContent(b));
            asis->SetBinContent(b, CB->GetBinContent(b));
            fit->SetBinContent(b, k * CB->GetBinContent(b));
        }
        Stack* st = new Stack();
        cbk->SetFillColor(kRed - 7); cbk->SetLineColor(kRed - 7); st->Add(cbk);
        std::vector<TH1D*> hc;
        for (const auto& cc : comps) {
            TH1D* h = (TH1D*)N->Clone(Form("mc_%s_%d", cc.name.c_str(), q)); h->Reset();
            for (int p : cc.parts) { TH1* x = allio::Sim(fh, "S0", "nom", p, lay); h->Add(x, kNbToRaw); delete x; }
            asis->Add(h);
            h->Scale(cc.s);
            fit->Add(h);
            h->SetFillColor(cc.col); h->SetLineColor(cc.col);
            st->Add(h);
            hc.push_back(h);
        }
        p1->cd();
        TH1D* fr = (TH1D*)N->Clone(Form("mfr%d", q)); fr->Reset();
        fr->GetXaxis()->SetRangeUser(0, xmax);
        double ymax = 0, ymin = 1e9;
        for (int b = 1; b <= N->GetNbinsX(); ++b) if (N->GetBinCenter(b) < xmax) { ymax = std::max(ymax, N->GetBinContent(b)); if (N->GetBinContent(b) > 0) ymin = std::min(ymin, N->GetBinContent(b)); }
        fr->SetMaximum(3 * ymax); fr->SetMinimum(std::max(0.5, 0.3 * ymin));
        fr->GetYaxis()->SetTitle("rows per 2 MeV/c^{2}"); fr->GetYaxis()->SetTitleSize(0.058); fr->GetYaxis()->SetLabelSize(0.052); fr->GetYaxis()->SetTitleOffset(1.35); fr->GetXaxis()->SetLabelSize(0);
        fr->Draw("AXIS");
        st->Draw(Form("mee%d", q));
        asis->SetLineColor(kBlack); asis->SetLineStyle(2); asis->SetLineWidth(2); asis->Draw("HIST SAME");
        TH1D* Nd = (TH1D*)N->Clone(Form("Nd%d", q));
        Nd->SetMarkerStyle(20); Nd->SetMarkerSize(0.7); Nd->SetMarkerColor(kBlack); Nd->SetLineColor(kBlack); Nd->Draw("E1 SAME");
        TBox* bx = new TBox(0, fr->GetMinimum(), 0.010, fr->GetMaximum()); bx->SetFillColorAlpha(kGray, 0.35); bx->Draw();
        fr->Draw("AXIS SAME");
        TLatex lt; lt.SetNDC(); lt.SetTextSize(0.058); lt.SetTextFont(42); lt.DrawLatex(0.18, 0.94, rtxt[q]);
        if (q == 2) {
            c->cd();
            TPad* pl = new TPad("mleg", "", 0.80, 0.36, 1.0, 0.93); pl->Draw(); pl->cd();
            TLegend* lg = new TLegend(0.03, 0.08, 0.99, 0.95); lg->SetTextSize(0.062); lg->SetFillStyle(0); lg->SetBorderSize(0);
            lg->AddEntry(Nd, "data N", "PE");
            for (int k = (int)comps.size() - 1; k >= 0; --k) lg->AddEntry(hc[k], Form("%s #times %.2f", comps[k].label.c_str(), comps[k].s), "F");
            lg->AddEntry(cbk, Form("k CB, k = %.2f (#geq 10 MeV)", kE2), "F");
            lg->AddEntry(asis, "SMASH as is + CB", "L");
            lg->Draw();
            p1->cd();
        }
        p2->cd();
        TH1D* R1 = (TH1D*)N->Clone(Form("mr1%d", q)); R1->Divide(fit);
        TH1D* R0 = (TH1D*)N->Clone(Form("mr0%d", q)); R0->Divide(asis);
        for (int b = 1; b <= N->GetNbinsX(); ++b) {
            R1->SetBinError(b, fit->GetBinContent(b) > 0 ? N->GetBinError(b) / fit->GetBinContent(b) : 0);
            R0->SetBinError(b, asis->GetBinContent(b) > 0 ? N->GetBinError(b) / asis->GetBinContent(b) : 0);
        }
        R1->GetXaxis()->SetRangeUser(0, xmax);
        R1->SetMinimum(0.3); R1->SetMaximum(1.9);
        R1->GetYaxis()->SetTitle("data / model"); R1->GetYaxis()->SetTitleSize(0.095); R1->GetYaxis()->SetLabelSize(0.085); R1->GetYaxis()->SetTitleOffset(0.85); R1->GetYaxis()->SetNdivisions(505);
        R1->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); R1->GetXaxis()->SetTitleSize(0.095); R1->GetXaxis()->SetLabelSize(0.085); R1->GetXaxis()->SetTitleOffset(1.2); R1->GetXaxis()->SetNdivisions(505);
        R1->SetMarkerStyle(20); R1->SetMarkerSize(0.6); R1->SetMarkerColor(kBlue + 1); R1->SetLineColor(kBlue + 1);
        R0->SetMarkerStyle(24); R0->SetMarkerSize(0.6); R0->SetMarkerColor(kGray + 2); R0->SetLineColor(kGray + 2);
        R1->Draw("E1"); R0->Draw("E1 SAME");
        TLine* l1 = new TLine(0, 1, xmax, 1); l1->SetLineStyle(2); l1->Draw();
        if (q == 0) { TLatex l2; l2.SetNDC(); l2.SetTextSize(0.085); l2.SetTextFont(42); l2.DrawLatex(0.30, 0.86, "#bullet after the fit   #circ SMASH as is"); }
    }
    c->cd();
    TLatex ln; ln.SetNDC(); ln.SetTextSize(0.019); ln.SetTextFont(42);
    ln.DrawLatex(0.01, 0.012, Form("hists_v1.root B_mee_R0, R1, R5 (unblinded regions); simulation x %.1f raw rows per nb (L = 6079.2 nb^{-1} x 63/64); scales of %s (no tilt); grey band: M(ee) < 10 MeV/c^{2}, CB at k = 1, not fitted",
                                   kNbToRaw, gNominal.c_str()));
    Save(c, "flex_mee");
}

// ---------------------------------------------------------------- scales summary
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
void DrawScales() {
    const FitInfo& F = gFit[gNominal];
    const FitInfo& A = gFit[gAlt];
    const FitInfo& E = gFit["V_E1k"];
    const FitInfo& Es = gFit["V_E2split"];
    Verd vb = ReadBias(gNominal);
    std::vector<std::string> rows = {"s_C2", "s_C3", "s_C4", "s_C5", "s_Geta", "s_GD_om", "s_GD_nom", "s_GC", "k_E2"};
    std::map<std::string, std::string> lab = {{"s_C2", "C2  #pi^{0} Dalitz (bump)"}, {"s_C3", "C3  conversion (bump)"}, {"s_C4", "C4  #eta#rightarrow#pi^{+}#pi^{-}#gamma conv."},
                                              {"s_C5", "C5  #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}"}, {"s_Geta", "G_{#eta}  other #eta"}, {"s_GD_om", "G_{D}  #omega/#rho^{0} pions"},
                                              {"s_GD_nom", "G_{D}  other pions"}, {"s_GC", "G_{C}  conversion pair"}, {"k_E2", "k_{E2}  CB scale"}};
    std::map<std::string, std::string> p3 = {{"s_C2", "C2"}, {"s_C3", "C3"}, {"s_C4", "C4"}, {"s_C5", "C5"}, {"s_Geta", "Geta"}, {"s_GD_om", "GD"}, {"s_GD_nom", "GD"}, {"s_GC", "GC"}, {"k_E2", "k_E2"}};
    TCanvas* c = new TCanvas("csc", "", 1700, 1150);
    c->SetLeftMargin(0.22); c->SetRightMargin(0.27); c->SetBottomMargin(0.12); c->SetTopMargin(0.155);
    const int n = (int)rows.size();
    TH2D* fr = new TH2D("scfr", "", 10, -0.6, 3.4, n, 0, n);
    for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, lab[rows[i]].c_str());
    fr->GetYaxis()->SetLabelSize(0.027); fr->GetXaxis()->SetLabelSize(0.027); fr->GetXaxis()->SetTitleSize(0.028); fr->GetXaxis()->SetTitleOffset(1.25);
    fr->GetXaxis()->SetTitle("scale relative to SMASH (#sigma #times PDG BR #times V3)"); fr->SetStats(0);
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
        pt(A, rows[i], y0 + 0.09, kBlue + 1, 24, false, rows[i] == "s_GD_om" || rows[i] == "s_GD_nom" ? "s_GD" : "");
        pt(Es, rows[i], y0 - 0.09, kGreen + 2, 26, false);
        pt(E, rows[i], y0 - 0.27, kOrange + 7, 32, false);
        if (gP3.count(p3[rows[i]])) {
            const double v = gP3[p3[rows[i]]], e = gP3e.count(p3[rows[i]]) ? gP3e[p3[rows[i]]] : 0.;
            const double y = y0 - 0.40;
            TLine* l = new TLine(v - e, y, v + e, y); l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->Draw();
            TGraph* g = new TGraph(1); g->SetPoint(0, v, y); g->SetMarkerStyle(25); g->SetMarkerColor(kRed + 1); g->SetMarkerSize(1.4); g->Draw("P");
        }
        // right margin: nominal value and the bias-test verdict
        const int k = Idx(F, rows[i]);
        TLatex t; t.SetTextSize(0.0185); t.SetTextFont(42); t.SetTextAlign(12);
        if (k >= 0) {
            const double rel = F.v[k] != 0 ? F.et[k] / std::fabs(F.v[k]) : 99;
            std::string vv = vb.v.count(rows[i]) ? vb.v[rows[i]] : (rows[i] == "k_E2" ? "nuisance" : "");
            t.DrawLatex(3.47, y0 + 0.27, Form("%.2f #pm %.2f #pm %.2f", F.v[k], F.es[k], F.em[k]));
            t.DrawLatex(3.47, y0 + 0.00, Form("%s (#sigma/s = %.0f %%)", rel < 0.3 ? "constrained" : rel < 0.5 ? "weakly constrained" : "not constrained", 100 * rel));
            if (vb.bias.count(rows[i])) {
                if (std::fabs(F.v[k]) > F.et[k]) t.DrawLatex(3.47, y0 - 0.27, Form("bias test %s, %+.0f %%", vv.c_str(), vb.bias[rows[i]]));
                else t.DrawLatex(3.47, y0 - 0.27, Form("bias test %s (truth #approx 0)", vv.c_str()));
            }
        }
    }
    for (int i = 1; i < n; ++i) { TLine* l = new TLine(-0.6, i, 3.4, i); l->SetLineColor(kGray); l->SetLineStyle(3); l->Draw(); }
    TLegend* lg = new TLegend(0.01, 0.855, 0.99, 0.995); lg->SetNColumns(2); lg->SetTextSize(0.022); lg->SetBorderSize(0); lg->SetFillStyle(0);
    TGraph* g1 = new TGraph(1); g1->SetMarkerStyle(20); g1->SetMarkerColor(kBlack); g1->SetLineColor(kBlack); g1->SetLineWidth(2);
    TGraph* g2 = new TGraph(1); g2->SetMarkerStyle(24); g2->SetMarkerColor(kBlue + 1); g2->SetLineColor(kBlue + 1); g2->SetLineWidth(2);
    TGraph* g3 = new TGraph(1); g3->SetMarkerStyle(26); g3->SetMarkerColor(kGreen + 2); g3->SetLineColor(kGreen + 2); g3->SetLineWidth(2);
    TGraph* g4 = new TGraph(1); g4->SetMarkerStyle(32); g4->SetMarkerColor(kOrange + 7); g4->SetLineColor(kOrange + 7); g4->SetLineWidth(2);
    TGraph* g5 = new TGraph(1); g5->SetMarkerStyle(25); g5->SetMarkerColor(kRed + 1); g5->SetLineColor(kRed + 1); g5->SetLineWidth(2);
    lg->AddEntry(g1, Form("nominal %s: M(ee) #geq 10 MeV (thick: stat, thin: stat #oplus MC)", gNominal.c_str()), "PL");
    lg->AddEntry(g3, "same model, M(ee) 10-20 and #geq 20 MeV fitted separately", "PL");
    lg->AddEntry(g2, Form("alternative %s: per-slice tilts instead of the #omega/#rho^{0} split", gAlt.c_str()), "PL");
    lg->AddEntry(g4, "same model + E1 (M(ee) < 10 MeV) with a free CB-like k_{E1}", "PL");
    lg->AddEntry(g5, "Part III F-main: a_{#eta}, a#kappa, a#rho_{d}, ag; c_{j}(1+bx), c_{j}#kappa_{B} (weighted)", "PL");
    lg->Draw();
    TLatex tn; tn.SetNDC(); tn.SetTextSize(0.019); tn.SetTextFont(42);
    tn.DrawLatex(0.01, 0.025, "S0, isBest #geq 0, M4 0.40-0.70, 5 M(#pi^{+}#pi^{-}) slices; 1 = SMASH as is. Right: nominal value #pm stat #pm MC; bias test = pseudo-data from the full simulation");
    tn.DrawLatex(0.01, 0.005, "at the fitted scales, fitted with bootstrap-replica templates (realistic template noise); PASS = |bias| < max(0.5 #sigma_{stat}, 10 %) and 68 % coverage in 0.60-0.76.");
    Save(c, "flex_scales");
}

// ---------------------------------------------------------------- correlation matrix
void DrawCorr(const std::string& tag, const std::string& base) {
    const FitInfo& F = gFit[tag];
    const int n = (int)F.names.size();
    TCanvas* c = new TCanvas(("cc" + tag).c_str(), "", 1100, 950);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.13); c->SetTopMargin(0.08);
    TH2D* h = new TH2D(("cor" + tag).c_str(), "", n, 0, n, n, 0, n);
    for (int a = 0; a < n; ++a) {
        h->GetXaxis()->SetBinLabel(a + 1, F.names[a].c_str());
        h->GetYaxis()->SetBinLabel(n - a, F.names[a].c_str());
        for (int b = 0; b < n; ++b) h->SetBinContent(a + 1, n - b, F.cor[a][b]);
    }
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetPalette(kTemperatureMap);
    h->SetMinimum(-1); h->SetMaximum(1); h->SetStats(0);
    h->GetXaxis()->SetLabelSize(0.040); h->GetYaxis()->SetLabelSize(0.040); h->SetMarkerSize(1.3);
    h->Draw("COLZ TEXT");
    TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextFont(42);
    t.DrawLatex(0.15, 0.94, Form("Correlations (stat #oplus MC) of the flexible fit %s, F10 E2, M4 0.40-0.70", tag.c_str()));
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
    TCanvas* c = new TCanvas(("cb" + tag).c_str(), "", 1700, 720);
    c->Divide(2, 1);
    for (int pnl = 0; pnl < 2; ++pnl) {
        c->cd(pnl + 1);
        gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.12); gPad->SetRightMargin(0.04);
        TH2D* fr = new TH2D(Form("bfr%d", pnl), "", 10, pnl == 0 ? -80 : 0, pnl == 0 ? 40 : 1, n, 0, n);
        for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, names[0][i].c_str());
        fr->GetXaxis()->SetTitle(pnl == 0 ? "bias of the realistic fit [% of the truth]" : "68 % coverage of stat #oplus MC");
        fr->GetXaxis()->SetTitleSize(0.050); fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.060); fr->SetStats(0);
        fr->Draw();
        TBox* bx = pnl == 0 ? new TBox(-10, 0, 10, n) : new TBox(0.60, 0, 0.76, n);
        bx->SetFillColorAlpha(kGreen - 9, 0.5); bx->Draw();
        fr->Draw("AXIS SAME");
        for (int s = 0; s < 2; ++s)
            for (int i = 0; i < (int)names[s].size(); ++i) {
                double x = pnl == 0 ? brel[s][i] : cov[s][i];
                if (pnl == 0 && names[s][i] == "s_C3" && s == 0) x = -999;   // truth -0.09: relative bias meaningless (in sigma below)
                const double y = n - i - 0.5 + (s == 0 ? 0.15 : -0.15);
                TGraph* g = new TGraph(1);
                g->SetPoint(0, pnl == 0 ? std::max(-79., std::min(39., x)) : x, y);
                g->SetMarkerStyle(s == 0 ? 20 : 24); g->SetMarkerSize(1.6); g->SetMarkerColor(s == 0 ? kBlack : kBlue + 1);
                if (x > -900) g->Draw("P");
                if (pnl == 0) {
                    TLatex t; t.SetTextSize(0.034); t.SetTextFont(42); t.SetTextAlign(12);
                    if (s == 0) t.DrawLatex(16, y, Form("%+.1f #sigma_{tot}", bsig[s][i]));
                }
            }
        TLatex t; t.SetNDC(); t.SetTextSize(0.040); t.SetTextFont(42);
        if (pnl == 0) t.DrawLatex(0.17, 0.94, "#bullet truth = flexible fit of the data   #circ truth = SMASH as is");
        else t.DrawLatex(0.17, 0.94, Form("500 pseudo-experiments per truth, model %s", tag.c_str()));
    }
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
    const double xlo[9] = {0, -0.6, -0.8, 0, 0, -1.0, 0.6, -3, 0.3}, xhi[9] = {0.8, 0.8, 1.4, 1.4, 6.5, 3.0, 1.8, 3.5, 1.6};
    TCanvas* c = new TCanvas("cvar", "", 2000, 1250);
    c->Divide(5, 2, 0.002, 0.004);
    const int nv = (int)V.size();
    for (int q = 0; q < (int)sc.size(); ++q) {
        c->cd(q + 1 + (q >= 5 ? 0 : 0));
        gPad->SetLeftMargin(q % 5 == 0 ? 0.42 : 0.06); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.10);
        TH2D* fr = new TH2D(Form("vfr%d", q), "", 10, xlo[q], xhi[q], nv, 0, nv);
        for (int i = 0; i < nv; ++i) fr->GetYaxis()->SetBinLabel(nv - i, q % 5 == 0 ? V[i].second.c_str() : "");
        fr->GetYaxis()->SetLabelSize(0.075); fr->GetXaxis()->SetLabelSize(0.06); fr->SetStats(0); fr->GetXaxis()->SetNdivisions(505);
        fr->Draw();
        TLine* one = new TLine(1, 0, 1, nv); one->SetLineStyle(2); one->SetLineColor(kGray + 2); if (1 > xlo[q] && 1 < xhi[q]) one->Draw();
        // nominal band
        const FitInfo& N = gFit[gNominal];
        const int kn = Idx(N, sc[q]);
        if (kn >= 0) { TBox* b = new TBox(std::max(xlo[q], N.v[kn] - N.et[kn]), 0, std::min(xhi[q], N.v[kn] + N.et[kn]), nv); b->SetFillColorAlpha(kGray, 0.4); b->Draw(); }
        for (int i = 0; i < nv; ++i) {
            if (!gFit.count(V[i].first)) continue;
            const FitInfo& f = gFit[V[i].first];
            int k = Idx(f, sc[q]);
            if (k < 0 && (sc[q] == "s_GD_om" || sc[q] == "s_GD_nom")) k = Idx(f, "s_GD");
            if (k < 0) continue;
            const double y = nv - i - 0.5;
            const double v = f.v[k], et = f.et[k];
            TLine* l = new TLine(std::max(xlo[q], v - et), y, std::min(xhi[q], v + et), y); l->SetLineWidth(2); l->SetLineColor(i == 0 ? kBlack : kBlue + 1); l->Draw();
            TGraph* g = new TGraph(1); g->SetPoint(0, std::max(xlo[q], std::min(xhi[q], v)), y); g->SetMarkerStyle(i == 0 ? 20 : 21); g->SetMarkerSize(1.2); g->SetMarkerColor(i == 0 ? kBlack : kBlue + 1); g->Draw("P");
        }
        TLatex t; t.SetNDC(); t.SetTextSize(0.075); t.SetTextFont(42); t.DrawLatex(q % 5 == 0 ? 0.45 : 0.10, 0.93, sc[q].c_str());
    }
    c->cd(10);
    TLatex t; t.SetNDC(); t.SetTextSize(0.062); t.SetTextFont(42);
    t.DrawLatex(0.04, 0.80, "Every scale in every variant");
    t.DrawLatex(0.04, 0.70, "(value #pm stat #oplus MC).");
    t.DrawLatex(0.04, 0.58, "Grey band: nominal fit.");
    t.DrawLatex(0.04, 0.48, "Dashed: SMASH (1).");
    t.DrawLatex(0.04, 0.36, "Variants use the nominal model");
    t.DrawLatex(0.04, 0.28, "(no tilt, G_{D} split), except");
    t.DrawLatex(0.04, 0.20, "alt., F7, F9.");
    Save(c, "flex_variants");
}

}  // namespace flexg

void flex_figs() {
    using namespace flexg;
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ");
    ReadResults();
    TFile* f = TFile::Open(Form("%s/flex_fit.root", kOut));
    TFile* fh = TFile::Open(kHistFile);
    const FitInfo& A = gFit["ASIS"];
    const FitInfo& N = gFit[gNominal];
    const FitInfo& Al = gFit[gAlt];
    DrawSpectra(f, "ASIS", "E2", 0.40, 0.70, "flex_spectra_asis",
                "Data vs SMASH as is (every component at #sigma_{SMASH} #times PDG BR #times V3, CB k = 1): M(e^{+}e^{-}) #geq 10 MeV/c^{2}",
                Form("S0, isBest #geq 0, L = 6079.2 nb^{-1}, raw rows. Hybrid GLS #chi^{2} %.0f / %d bins (Pearson %.0f). No parameter is fitted.", A.chi2H, 350, A.chi2P));
    DrawSpectra(f, gNominal, "E2", 0.40, 0.70, "flex_spectra_fit",
                Form("Data vs SMASH after the flexible fit %s (one free scale per component, no tilt): M(e^{+}e^{-}) #geq 10 MeV/c^{2}", gNominal.c_str()),
                Form("S0, isBest #geq 0, L = 6079.2 nb^{-1}, raw rows. Hybrid GLS #chi^{2} %.0f / %d (Pearson %.0f); expected from template noise 363 #pm 42 (pseudo-data).", N.chi2H, N.ndf, N.chi2P));
    DrawSpectra(f, gAlt, "E2", 0.40, 0.70, "flex_spectra_alt",
                Form("Data vs SMASH after the alternative flexible fit %s (7 scales + per-slice tilts of G_{D}, G_{C}): M(e^{+}e^{-}) #geq 10 MeV/c^{2}", gAlt.c_str()),
                Form("S0, isBest #geq 0, L = 6079.2 nb^{-1}, raw rows. Hybrid GLS #chi^{2} %.0f / %d (Pearson %.0f).", Al.chi2H, Al.ndf, Al.chi2P));
    DrawMee(fh);
    DrawScales();
    DrawCorr(gNominal, "flex_corr");
    DrawCorr(gAlt, "flex_corr_alt");
    DrawBias(gNominal);
    DrawVariants();
    // range to 1.00: summed over the slices, as is and after the fit
    {
        const FitInfo& R = gFit["V_r100"];
        DrawSpectra(f, "V_r100_asis", "E2", 0.40, 1.00, "flex_r100_asis", "M4 0.40-1.00, SMASH as is (sum of the slices, M(ee) #geq 10 MeV/c^{2})", "", true);
        const FitInfo& Rr = gFit["V_r100raw"];
        DrawSpectra(f, "V_r100raw", "E2", 0.40, 1.00, "flex_r100raw",
                    Form("M4 0.40-1.00 after the flexible fit with raw background templates (model %s, sum of the slices)", gNominal.c_str()),
                    Form("Pearson #chi^{2}/ndf %.0f/%d, hybrid GLS %.0f", Rr.chi2P, Rr.ndf, Rr.chi2H), true);
        DrawSpectra(f, "V_r100", "E2", 0.40, 1.00, "flex_r100",
                    Form("M4 0.40-1.00 after the flexible fit (model %s, sum of the slices)", gNominal.c_str()),
                    Form("Pearson #chi^{2}/ndf %.0f/%d, hybrid GLS %.0f; s_{GC} = %.2f (unphysical): the components cannot describe 0.70-1.00", R.chi2P, R.ndf, R.chi2H, R.v.size() > 7 ? R.v[7] : 0.), true);
    }
    printf("flex_figs: done\n");
}
