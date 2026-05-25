// mass_ee_ratio_const_fit_sim.C — three ratio panels for SMASH simulation
// with FIXED 20-bin binning over [0, 1.4] GeV/c², round-dot markers,
// vertical error bars only, and a constant fit y = a drawn on each panel.
//
// Analog of mass_ee_ratio_const_fit_exp.C but:
//   * no factor 63 (no PT2 downscale in sim)
//   * no CB subtraction (sim has pure signal after MC purity gate)
//   * three histogram variants per name: H (no trigger), H_pt3, H_pt2
//
// Three panels (all on weighted signal):
//   1) N_PT2 / N_PT3        — correction factor for the PT3-selected sample
//   2) N_PT3 / N_PT2        — efficiency PT3 conditional on PT2 (= 1/panel-1)
//   3) N_PT3 / N_none       — ABSOLUTE PT3 efficiency vs the unbiased sample
//                             (sim-specific — exp data has no unbiased ref)
//
// On each panel: red horizontal line = fit  y = a , with fitted a ± error
// and χ²/ndf annotated in red just above the line.
//
// Arguments (all optional):
//   1) oa_cut_deg   — OA threshold, default 0 ⇒ no cut
//   2) use_weights  — true → weighted χ², false → unweighted (each bin = 1)
//                     default true
//   3) fit_xmin     — lower bound of fit range [GeV/c²], default -1 ⇒ full
//   4) fit_xmax     — upper bound of fit range [GeV/c²], default -1 ⇒ full
//
// Every event is weighted by `sim_genweight`. Filename gets "_oa{N}",
// "_unw", "_fit{lo}-{hi}" suffixes only when the option differs from
// default so successive runs do not overwrite each other.
//
// Usage:
//   root -l -b -q plots/mass_ee_ratio_const_fit_sim.C
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_sim.C(4)'
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_sim.C(0, false)'
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_sim.C(0, true, 0.1, 0.8)'
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_sim.C(4, true, 0.1, 0.8)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";
constexpr const char* kNtName    = "dilepton_nt";

// Detect whether sim_genweight is filled (lepton mode) or empty (std mode)
// and return the appropriate TTree::Draw weight expression. Macro-level
// decision — NOT per-event — so stray zero-weight events in lepton data
// don't get over-weighted by a fallback ternary.
std::string detectWeightExpr(TTree* nt) {
    TH1D* h = new TH1D("h_wsum_detect", "", 1, 0, 2);
    nt->Draw("1>>h_wsum_detect", "sim_genweight", "goff");
    const double sum = h->Integral();
    delete h;
    if (sum > 0.0) {
        std::cout << "[sim] sim_genweight integral = " << sum
                  << "  → using weighted fills (sim_genweight)\n";
        return "sim_genweight";
    }
    std::cout << "[sim] sim_genweight integral = 0  → unweighted (1.0)\n";
    return "1.0";
}

TH1D* drawFromNt(TTree* t, const std::string& cut_w, const std::string& name,
                 int nb, double xmin, double xmax) {
    TH1D* h = new TH1D(name.c_str(), "", nb, xmin, xmax);
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut_w.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d != 0.0 && std::isfinite(d) && std::isfinite(n)) {
            const double val = n / d;
            const double rel = std::sqrt(
                (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
                std::pow(ed / d, 2));
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

void styleDot(TH1D* h, Color_t color = kBlue + 1) {
    h->SetMarkerStyle(20);   // full circle
    h->SetMarkerSize(0.9);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

FitResult fitConst(TH1D* h, double fit_xmin, double fit_xmax,
                   bool use_weights, bool draw, Color_t color) {
    const std::string fname = std::string(h->GetName()) +
                              (use_weights ? "_fconst_w" : "_fconst_u");
    TF1* f = new TF1(fname.c_str(), "[0]", fit_xmin, fit_xmax);
    f->SetParName(0, "a");
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    const std::string opt =
        std::string("RQ") + (use_weights ? "" : "W") + (draw ? "" : "N0");
    h->Fit(f, opt.c_str());
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

}  // anonymous namespace

void mass_ee_ratio_const_fit_sim(double oa_cut_deg  = 0.0,
                                 bool   use_weights = true,
                                 double fit_xmin    = -1.0,
                                 double fit_xmax    = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);   // no horizontal error bars

    const int    nb   = 20;
    const double xmin = 0.0;
    const double xmax = 1.4;

    // Resolve fit range — sentinel -1 ⇒ full histogram range.
    const bool   custom_fit_range =
        (fit_xmin >= 0.0 && fit_xmin < xmax) ||
        (fit_xmax >  0.0 && fit_xmax <= xmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < xmax) ? fit_xmin : xmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= xmax) ? fit_xmax : xmax;

    TFile* f = TFile::Open(kInputFile, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << kInputFile << "\n"; return; }
    TTree* nt = dynamic_cast<TTree*>(f->Get(kNtName));
    if (!nt) { std::cerr << kNtName << " missing\n"; return; }

    // Detect once whether the sample carries per-event weights.
    const std::string w_expr = detectWeightExpr(nt);

    // --- Compose per-event weight expressions. The OA cut is "oa > oa_cut_deg";
    //     0 (default) ⇒ no filtering. The weight is sim_genweight (lepton)
    //     or 1.0 (std), decided above by detectWeightExpr.
    const std::string oa_clause =
        (oa_cut_deg > 0.0) ? std::string(" * (oa>") + std::to_string(oa_cut_deg) + ")"
                           : std::string{};
    const std::string w_none = w_expr + oa_clause;
    const std::string w_pt3  = std::string("(pt3==1) * ") + w_expr + oa_clause;
    const std::string w_pt2  = std::string("(pt2==1) * ") + w_expr + oa_clause;
    std::cout << "Weight expressions:\n"
              << "  none = " << w_none << "\n"
              << "  PT3  = " << w_pt3  << "\n"
              << "  PT2  = " << w_pt2  << "\n"
              << "Fit mode (drawn): "
              << (use_weights ? "WEIGHTED (uses bin errors as 1/sigma^2)"
                              : "UNWEIGHTED (all bins weight 1, simple mean)")
              << "\nFit range:        [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_fit_range ? "  (custom)" : "  (full)") << "\n";

    TH1D* h_none = drawFromNt(nt, w_none, "h_none_fit", nb, xmin, xmax);
    TH1D* h_pt3  = drawFromNt(nt, w_pt3,  "h_pt3_fit",  nb, xmin, xmax);
    TH1D* h_pt2  = drawFromNt(nt, w_pt2,  "h_pt2_fit",  nb, xmin, xmax);

    TH1D* r_pt2_pt3  = makeRatio(h_pt2, h_pt3, "r_pt2_pt3");
    TH1D* r_pt3_pt2  = makeRatio(h_pt3, h_pt2, "r_pt3_pt2");
    TH1D* r_pt3_none = makeRatio(h_pt3, h_none, "r_pt3_none");

    styleDot(r_pt2_pt3); styleDot(r_pt3_pt2); styleDot(r_pt3_none);

    const std::string oa_tag =
        (oa_cut_deg > 0.0)
            ? Form(", OA > %g#circ", oa_cut_deg)
            : std::string{};
    r_pt2_pt3->SetTitle((std::string(
        "N_{PT2} / N_{PT3} (signal, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}").c_str());
    r_pt3_pt2->SetTitle((std::string(
        "N_{PT3} / N_{PT2} (signal, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT3}/N_{PT2}").c_str());
    r_pt3_none->SetTitle((std::string(
        "N_{PT3} / N_{none} (signal, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT3}/N_{none}").c_str());

    TCanvas* c = new TCanvas("c_mass_ee_ratio_const_fit_sim",
                             "SMASH M_{ee} ratio constant fit (20 fixed bins)",
                             1800, 600);
    c->Divide(3, 1, 0.001, 0.001);

    auto drawPad = [&](int idx, TH1D* h, double ref_y,
                       const std::string& /*ref_label*/,
                       int ref_style, int ref_color, double y_max_floor) {
        c->cd(idx);
        const double pad_lm = 0.14, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        // Auto Y range, +25% margin, never below y_max_floor.
        double ymax = 0.0;
        for (int b = 1; b <= nb; ++b) {
            const double v = h->GetBinContent(b) + h->GetBinError(b);
            if (std::isfinite(v) && v > ymax) ymax = v;
        }
        const double y_hi = std::max(ymax * 1.25, y_max_floor);
        h->GetYaxis()->SetRangeUser(0.0, y_hi);

        h->Draw("E1");

        TLine* lref = new TLine(xmin, ref_y, xmax, ref_y);
        lref->SetLineStyle(ref_style);
        lref->SetLineColor(ref_color);
        if (ref_style == 2) lref->SetLineWidth(2);
        lref->Draw();

        // Both fits computed; drawn one respects use_weights.
        const FitResult fit  = fitConst(h, fxlo, fxhi, use_weights,
                                        /*draw=*/true,  kRed + 1);
        const FitResult fit2 = fitConst(h, fxlo, fxhi, !use_weights,
                                        /*draw=*/false, kRed + 1);

        // Annotate fit value ± error just above the fit line itself.
        const double y_frac = (y_hi > 0) ? (fit.a / y_hi) : 0.0;
        const double y_line_ndc =
            pad_bm + (1.0 - pad_bm - pad_tm) * y_frac;
        const double y_text =
            (y_line_ndc < 0.70) ? y_line_ndc + 0.045 : y_line_ndc - 0.075;

        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.050);
        tex.SetTextColor(kRed + 1);
        tex.DrawLatex(0.52, y_text,
            Form("a = %.4f #pm %.4f", fit.a, fit.e));
        tex.SetTextSize(0.035);
        tex.SetTextColor(kGray + 3);
        tex.DrawLatex(0.52, y_text - 0.045,
            Form("#chi^{2}/ndf = %.1f / %d",
                 fit.chi2_ndf * fit.ndf, fit.ndf));

        const char* tag_drawn = use_weights ? "weighted  " : "unweighted";
        const char* tag_other = use_weights ? "unweighted" : "weighted  ";
        std::cout << "  pad " << idx
                  << "  [" << tag_drawn << " drawn] a = " << fit.a
                  << " ± " << fit.e
                  << "   chi2/ndf = " << (fit.chi2_ndf * fit.ndf)
                  << "/" << fit.ndf << "\n"
                  << "          [" << tag_other << "]       a = " << fit2.a
                  << " ± " << fit2.e
                  << "   chi2/ndf = " << (fit2.chi2_ndf * fit2.ndf)
                  << "/" << fit2.ndf << "\n";
    };

    std::cout << "Constant fit y = a:\n";
    drawPad(1, r_pt2_pt3,  1.0, "1.0 (PT3 = PT2)",       3, kGray + 2, 2.0);
    drawPad(2, r_pt3_pt2,  1.0, "1.0 (PT3 = PT2)",       3, kGray + 2, 1.2);
    drawPad(3, r_pt3_none, 1.0, "1.0 (PT3 = none)",      3, kGray + 2, 1.2);

    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_fit_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_ratio_const_fit_sim" + oa_suffix + wt_suffix + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "Saved: " << base << ".{pdf,png}\n";

    f->Close();
}
