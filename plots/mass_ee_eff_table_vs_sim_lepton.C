// mass_ee_eff_table_vs_sim_lepton.C — overlay of IZA trigger-efficiency data
// (provided externally as a HARDCODED table of (M, y) pairs) with SMASH
// LEPTON simulation, plus a (sim·scale)/IZA shape-ratio panel.
//
// Two panels:
//
//   Pad 1 — Trigger efficiency overlay:
//     IZA: 28 fixed bins on [0, 1.12] GeV/c² with values taken directly
//          from the kIzaTable below (PT3/(63·PT2) per bin from external
//          analysis).
//     sim: SAME 28 bins on [0, 1.12] from output_epem_sim.root, weighted
//          by sim_genweight (lepton mode). Computed as N_PT3 / N_PT2
//          (no factor 63 in sim — no PT2 downscale).
//     Both fitted with constant y = a (separate fits), values annotated.
//
//   Pad 2 — (sim · scale) / IZA shape ratio:
//     scale = a_IZA / a_sim from the two fits in pad 1.
//     A third constant fit is overlaid on the ratio; values near 1.0 mean
//     sim and IZA have the same mass shape, deviations expose residual
//     differences after normalisation.
//
// Bin layout: 28 bins of width 40 MeV, edges 0.00, 0.04, ..., 1.12.
// Bin i (1-based) center = 0.02 + (i-1)·0.04, matching kIzaTable order.
//
// IZA errors are NOT provided in the table; we set them to a small
// symbolic value (0.01 absolute) so both weighted and unweighted fits
// produce meaningful results. Tune `kIzaBinErr` below if you have real
// per-bin uncertainties to plug in.
//
// Arguments (all optional, semantics identical to the other macros):
//   1) oa_cut_deg   — OA threshold applied to SIM only (table is fixed);
//                     default 0 ⇒ no cut
//   2) use_weights  — true → weighted χ² fit; false → unweighted; default true
//   3) fit_xmin     — fit range lower bound [GeV/c²]; -1 ⇒ full
//   4) fit_xmax     — fit range upper bound [GeV/c²]; -1 ⇒ full
//
// Usage:
//   root -l -b -q plots/mass_ee_eff_table_vs_sim_lepton.C
//   root -l -b -q 'plots/mass_ee_eff_table_vs_sim_lepton.C(4)'
//   root -l -b -q 'plots/mass_ee_eff_table_vs_sim_lepton.C(0, true, 0.1, 0.8)'
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
#include <limits>
#include <string>
#include <vector>

namespace {

// -------- Binning (table-driven) --------
constexpr int    kNb   = 28;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.12;     // 28 bins × 0.04 GeV
constexpr double kBinW = 0.04;

// -------- IZA table — bin centers and PT3/(63·PT2) values --------
struct IzaPoint { double m; double y; };
const std::vector<IzaPoint> kIzaTable = {
    {0.02, 0.678}, {0.06, 0.752}, {0.10, 0.788}, {0.14, 0.787},
    {0.18, 0.806}, {0.22, 0.823}, {0.26, 0.832}, {0.30, 0.832},
    {0.34, 0.823}, {0.38, 0.858}, {0.42, 0.839}, {0.46, 0.849},
    {0.50, 0.810}, {0.54, 0.792}, {0.58, 0.812}, {0.62, 0.802},
    {0.66, 0.803}, {0.70, 0.793}, {0.74, 0.789}, {0.78, 0.782},
    {0.82, 0.772}, {0.86, 0.760}, {0.90, 0.722}, {0.94, 0.722},
    {0.98, 0.697}, {1.02, 0.712}, {1.06, 0.581}, {1.10, 0.730},
};
constexpr double kIzaBinErr = 0.01;   // symbolic, no real errors supplied

// -------- Colors --------
constexpr Color_t kColIza    = kBlue   + 1;
constexpr Color_t kColIzaFit = kBlue   + 2;
constexpr Color_t kColSim    = kRed    + 1;
constexpr Color_t kColSimFit = kRed    + 2;
constexpr Color_t kColRatio  = kViolet + 1;

// -------- Helpers --------

TH1D* buildIzaTable(const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->SetDirectory(nullptr);
    h->Sumw2();
    for (size_t i = 0; i < kIzaTable.size(); ++i) {
        const int b = static_cast<int>(i) + 1;
        if (b > kNb) break;
        h->SetBinContent(b, kIzaTable[i].y);
        h->SetBinError  (b, kIzaBinErr);
    }
    return h;
}

TH1D* drawFromNt(TTree* t, const std::string& cut_or_weight,
                 const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut_or_weight.c_str(), "goff");
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

TH1D* makeScaledSimOverExp(TH1D* sim, TH1D* exp_h, double scale,
                           double scale_err, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(sim->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = sim->GetNbinsX();
    const double rel_scale =
        (scale != 0.0) ? std::abs(scale_err / scale) : 0.0;
    for (int b = 1; b <= nb; ++b) {
        const double s  = sim->GetBinContent(b);
        const double es = sim->GetBinError(b);
        const double e  = exp_h->GetBinContent(b);
        const double ee = exp_h->GetBinError(b);
        if (e != 0.0 && std::isfinite(e) && std::isfinite(s)) {
            const double val = scale * s / e;
            const double rel = std::sqrt(
                (s != 0.0 ? std::pow(es / s, 2) : 0.0) +
                std::pow(ee / e, 2) +
                rel_scale * rel_scale);
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

void styleDot(TH1D* h, Color_t color, Style_t marker) {
    h->SetMarkerStyle(marker);
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
                   bool use_weights, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fconst";
    TF1* f = new TF1(fname.c_str(), "[0]", fit_xmin, fit_xmax);
    f->SetParName(0, "a");
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    const std::string opt = std::string("RQ") + (use_weights ? "" : "W");
    h->Fit(f, opt.c_str());
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

double safeMaxAbs(TH1D* h) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = std::abs(h->GetBinContent(b)) + h->GetBinError(b);
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

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

}  // anonymous namespace

void mass_ee_eff_table_vs_sim_lepton(double oa_cut_deg  = 0.0,
                                     bool   use_weights = true,
                                     double fit_xmin    = -1.0,
                                     double fit_xmax    = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    // Resolve fit range: sentinel −1 ⇒ use the histogram's full range.
    const bool   custom_range =
        (fit_xmin >= 0.0 && fit_xmin < kXmax) ||
        (fit_xmax >  0.0 && fit_xmax <= kXmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kXmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kXmax;

    // -------- IZA from table --------
    TH1D* h_iza = buildIzaTable("h_iza_eff_table");

    // -------- SIM lepton from file --------
    TFile* fs = TFile::Open("output_epem_sim.root", "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open output_epem_sim.root\n";
        return;
    }
    TTree* nt_sim = dynamic_cast<TTree*>(fs->Get("dilepton_nt"));
    if (!nt_sim) {
        std::cerr << "dilepton_nt missing in output_epem_sim.root\n";
        return;
    }
    const std::string w_expr = detectWeightExpr(nt_sim);
    const std::string oa_clause =
        (oa_cut_deg > 0.0)
            ? std::string(" * (oa>") + std::to_string(oa_cut_deg) + ")"
            : std::string{};
    const std::string w_pt3 = "(pt3==1) * " + w_expr + oa_clause;
    const std::string w_pt2 = "(pt2==1) * " + w_expr + oa_clause;

    TH1D* h_sim_pt3 = drawFromNt(nt_sim, w_pt3, "h_sim_pt3_table");
    TH1D* h_sim_pt2 = drawFromNt(nt_sim, w_pt2, "h_sim_pt2_table");
    TH1D* h_sim     = makeRatio(h_sim_pt3, h_sim_pt2, "h_sim_eff_table");

    std::cout << "Fit mode: "
              << (use_weights ? "WEIGHTED" : "UNWEIGHTED")
              << "\nFit range: [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_range ? "  (custom)" : "  (full)") << "\n";

    // -------- Styling and titles --------
    styleDot(h_iza, kColIza, 20);
    styleDot(h_sim, kColSim, 21);

    const std::string oa_tag =
        (oa_cut_deg > 0.0)
            ? Form(", SIM OA > %g#circ", oa_cut_deg)
            : std::string{};
    h_iza->SetTitle((std::string(
        "Trigger efficiency  IZA(table) vs sim(lepton) — 28 fixed bins") +
        oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "PT3 / (63 #upoint PT2) (IZA)   PT3/PT2 (sim)").c_str());

    // -------- Fits (BEFORE drawing so canvas reuses results) --------
    const FitResult f_iza = fitConst(h_iza, fxlo, fxhi, use_weights, kColIzaFit);
    const FitResult f_sim = fitConst(h_sim, fxlo, fxhi, use_weights, kColSimFit);
    std::cout << "\nPad 1 fits (constant y = a):\n"
              << "  IZA(table)   a = " << f_iza.a << " ± " << f_iza.e
              << "   chi2/ndf = " << (f_iza.chi2_ndf * f_iza.ndf)
              << "/" << f_iza.ndf << "\n"
              << "  sim(lepton)  a = " << f_sim.a << " ± " << f_sim.e
              << "   chi2/ndf = " << (f_sim.chi2_ndf * f_sim.ndf)
              << "/" << f_sim.ndf << "\n";

    // -------- Scale + scaled-sim/IZA ratio --------
    const double scale =
        (f_sim.a != 0.0) ? f_iza.a / f_sim.a : 1.0;
    const double scale_err = std::abs(scale) * std::sqrt(
        (f_iza.a != 0.0 ? std::pow(f_iza.e / f_iza.a, 2) : 0.0) +
        (f_sim.a != 0.0 ? std::pow(f_sim.e / f_sim.a, 2) : 0.0));

    TH1D* h_ratio = makeScaledSimOverExp(
        h_sim, h_iza, scale, scale_err, "h_ratio_eff_table");
    styleDot(h_ratio, kColRatio, 20);
    h_ratio->SetTitle((std::string(
        "(sim #upoint scale) / IZA — shape ratio (28 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "(PT3/PT2)_{sim} #upoint scale / (PT3/(63#upointPT2))_{IZA}").c_str());

    const FitResult f_ratio =
        fitConst(h_ratio, fxlo, fxhi, use_weights, kColRatio);
    std::cout << "\nPad 2 fit (constant y = a on the ratio):\n"
              << "  scale = a_IZA / a_sim = " << scale << " ± " << scale_err << "\n"
              << "  ratio fit a = " << f_ratio.a << " ± " << f_ratio.e
              << "   chi2/ndf = " << (f_ratio.chi2_ndf * f_ratio.ndf)
              << "/" << f_ratio.ndf << "\n";

    // -------- Canvas: 2 panels --------
    TCanvas* c = new TCanvas(
        "c_mass_ee_eff_table_vs_sim_lepton",
        "Trigger eff (IZA table vs sim lepton) + scaled-sim/IZA ratio",
        1500, 600);
    c->Divide(2, 1, 0.001, 0.001);

    // Pad 1 — Efficiency overlay
    c->cd(1);
    {
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        const double ymax = std::max(safeMaxAbs(h_iza), safeMaxAbs(h_sim));
        const double y_hi = std::max(ymax * 1.30, 1.2);
        h_iza->GetYaxis()->SetRangeUser(0.0, y_hi);
        h_iza->Draw("E1");          // fit attached
        h_sim->Draw("E1 SAME");

        TLine* l1 = new TLine(kXmin, 1.0, kXmax, 1.0);
        l1->SetLineStyle(3); l1->SetLineColor(kGray + 2); l1->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColIzaFit);
        tex.DrawLatex(0.16, 0.85,
            Form("IZA:  a = %.4f #pm %.4f", f_iza.a, f_iza.e));
        tex.SetTextSize(0.030);
        tex.DrawLatex(0.16, 0.81,
            Form("       #chi^{2}/ndf = %.1f / %d",
                 f_iza.chi2_ndf * f_iza.ndf, f_iza.ndf));

        tex.SetTextAlign(31);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColSimFit);
        tex.DrawLatex(0.95, 0.85,
            Form("sim:  a = %.4f #pm %.4f", f_sim.a, f_sim.e));
        tex.SetTextSize(0.030);
        tex.DrawLatex(0.95, 0.81,
            Form("#chi^{2}/ndf = %.1f / %d       ",
                 f_sim.chi2_ndf * f_sim.ndf, f_sim.ndf));

        TLegend* leg = new TLegend(0.40, 0.74, 0.65, 0.83);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_iza, "IZA (table)",   "lpe");
        leg->AddEntry(h_sim, "sim (lepton)",  "lpe");
        leg->Draw();
    }

    // Pad 2 — (sim·scale)/IZA ratio with a constant fit
    c->cd(2);
    {
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        // Y range: centred on 1, robust ignoring big-error bins.
        double ymax = 0.0, ymin = std::numeric_limits<double>::infinity();
        for (int b = 1; b <= kNb; ++b) {
            const double v = h_ratio->GetBinContent(b);
            const double e = h_ratio->GetBinError(b);
            if (!std::isfinite(v) || v == 0.0) continue;
            if (e > 0.0 && std::abs(e / v) > 0.5) continue;
            ymax = std::max(ymax, v + e);
            ymin = std::min(ymin, std::max(0.0, v - e));
        }
        if (!std::isfinite(ymin)) ymin = 0.0;
        if (ymax <= 0.0)          ymax = 2.0;
        const double y_lo = std::max(0.0, std::min(ymin, 0.5) - 0.15);
        const double y_hi = std::max(ymax * 1.20, 1.5);
        h_ratio->GetYaxis()->SetRangeUser(y_lo, y_hi);

        h_ratio->Draw("E1");

        // Reference line at 1.0 (perfect agreement after rescaling).
        TLine* l1 = new TLine(kXmin, 1.0, kXmax, 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kGray + 3); l1->SetLineWidth(2);
        l1->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColRatio);
        tex.DrawLatex(0.16, 0.85,
            Form("scale = a_{IZA}/a_{sim} = %.4f #pm %.4f",
                 scale, scale_err));
        tex.SetTextSize(0.035);
        tex.DrawLatex(0.16, 0.80,
            Form("ratio fit  a = %.4f #pm %.4f", f_ratio.a, f_ratio.e));
        tex.SetTextSize(0.030);
        tex.SetTextColor(kGray + 3);
        tex.DrawLatex(0.16, 0.76,
            Form("#chi^{2}/ndf = %.1f / %d",
                 f_ratio.chi2_ndf * f_ratio.ndf, f_ratio.ndf));
        tex.DrawLatex(0.16, 0.72,
            "fit = 1 #Rightarrow sim and IZA identical mass shape");
    }

    // -------- Save --------
    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_eff_table_vs_sim_lepton" +
        oa_suffix + wt_suffix + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "Saved: " << base << ".{pdf,png}\n";

    fs->Close();
}
