// research/plots/compare_mu_corr.C — overlay μ vs OA before/after ECAL
// energy-scale correction.
//
// Reads:
//   fit_results_<fl>.root           — uncorrected (μ_orig)
//   fit_results_<fl>_ecalcor.root   — corrected   (μ_corr)
// for the π⁰ pipeline, and the analogous fit_results_eta_*.root files
// for the η pipeline. Produces:
//   plots/output/compare_mu_pi0_<fl>_vs_oa.{pdf,png}
//   plots/output/compare_mu_eta_<fl>_vs_oa.{pdf,png}
//
// Each panel overlays two TGraphErrors on the same axes:
//   - black: μ original
//   - red:   μ ecalcor
// plus a horizontal PDG reference line.
//
// Usage (from research/):
//   root -l -b -q plots/compare_mu_corr.C            # REC, both pi0 + eta
//   root -l -b -q 'plots/compare_mu_corr.C("cor")'   # COR (if available)

#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace {
    constexpr double kPi0PDG = 0.13498;
    constexpr double kEtaPDG = 0.5478;
}

// Read μ vs OA from a fit_results TTree, sorted by oa_lo. Skips panel_idx == 0
// (the integrated full panel).
struct OAGraph {
    std::vector<double> X, EX, Y, EY;
};

OAGraph readMuVsOA(const std::string& fpath) {
    OAGraph g;
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fpath << "\n"; return g;
    }
    auto* t = (TTree*)f->Get("fit_results");
    if (!t) {
        std::cerr << "TTree 'fit_results' not in " << fpath << "\n"; return g;
    }
    int   panel_idx = 0;
    float oa_lo = 0, oa_hi = 0, mu = 0, sigma = 0;
    t->SetBranchAddress("panel_idx", &panel_idx);
    t->SetBranchAddress("oa_lo",     &oa_lo);
    t->SetBranchAddress("oa_hi",     &oa_hi);
    t->SetBranchAddress("mu",        &mu);
    t->SetBranchAddress("sigma",     &sigma);

    struct Row { double x, ex, y, ey; };
    std::vector<Row> rows;
    for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
        t->GetEntry(ev);
        if (panel_idx == 0) continue;
        rows.push_back({0.5 * (oa_lo + oa_hi),
                        0.5 * (oa_hi - oa_lo),
                        mu, sigma});
    }
    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) { return a.x < b.x; });
    for (const auto& r : rows) {
        g.X .push_back(r.x);
        g.EX.push_back(r.ex);
        g.Y .push_back(r.y);
        g.EY.push_back(r.ey);
    }
    f->Close();
    return g;
}

void plotPanel(const std::string& title, const std::string& xtitle,
               const std::string& ytitle,
               const OAGraph& g_orig, const OAGraph& g_corr,
               double pdg_ref, const std::string& out_base)
{
    if (g_orig.X.empty() && g_corr.X.empty()) return;

    auto* g1 = new TGraphErrors((int)g_orig.X.size(),
        g_orig.X.data(), g_orig.Y.data(),
        g_orig.EX.data(), g_orig.EY.data());
    g1->SetMarkerColor(kBlack); g1->SetLineColor(kBlack);
    g1->SetMarkerStyle(20); g1->SetMarkerSize(0.9);

    auto* g2 = new TGraphErrors((int)g_corr.X.size(),
        g_corr.X.data(), g_corr.Y.data(),
        g_corr.EX.data(), g_corr.EY.data());
    g2->SetMarkerColor(kRed); g2->SetLineColor(kRed);
    g2->SetMarkerStyle(21); g2->SetMarkerSize(0.9);

    const TString ctitle = TString::Format("%s;%s;%s",
                                           title.c_str(),
                                           xtitle.c_str(), ytitle.c_str());

    auto* c = new TCanvas(("c_cmp_" + out_base).c_str(),
                          out_base.c_str(), 1100, 700);
    c->SetMargin(0.13, 0.05, 0.12, 0.08);
    c->SetGrid();
    g1->SetTitle(ctitle);
    g1->Draw("AP");
    g1->GetXaxis()->SetLimits(0.0, 15.0);
    g2->Draw("P SAME");

    auto* lref = new TLine(0.0, pdg_ref, 15.0, pdg_ref);
    lref->SetLineColor(kBlue);
    lref->SetLineStyle(2);
    lref->SetLineWidth(2);
    lref->Draw();

    auto* leg = new TLegend(0.55, 0.72, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.034);
    leg->AddEntry(g1,   "#mu original (no correction)",        "lpe");
    leg->AddEntry(g2,   "#mu ECAL-corrected",                  "lpe");
    leg->AddEntry(lref, TString::Format("PDG = %.5f", pdg_ref).Data(), "l");
    leg->Draw();

    c->Update();
    const std::string out_pdf = "plots/output/" + out_base + ".pdf";
    const std::string out_png = "plots/output/" + out_base + ".png";
    c->SaveAs(out_pdf.c_str());
    c->SaveAs(out_png.c_str());
    std::cout << "  wrote " << out_pdf << "\n";
    std::cout << "  wrote " << out_png << "\n";
}

void compare_mu_corr(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    // SIM file naming: fit_results_<fl>_sim.root and
    // fit_results_<fl>_ecalcor_sim.root.

    // -- π⁰ ---------------------------------------------------------------
    OAGraph pi0_orig = readMuVsOA("fit_results_" + fl + "_sim.root");
    OAGraph pi0_corr = readMuVsOA("fit_results_" + fl + "_ecalcor_sim.root");
    plotPanel(
        TString::Format("#pi^{0} #mu vs OA (%s, SIM) — before/after ECAL correction", fl.c_str()).Data(),
        "OA(e^{+}e^{-}) [deg]",
        "#mu_{CB} [GeV/c^{2}]  (vertical bars = #pm#sigma_{CB})",
        pi0_orig, pi0_corr, kPi0PDG,
        "compare_mu_pi0_" + fl + "_sim_vs_oa");

    // -- η ----------------------------------------------------------------
    OAGraph eta_orig = readMuVsOA("fit_results_eta_" + fl + "_sim.root");
    OAGraph eta_corr = readMuVsOA("fit_results_eta_" + fl + "_ecalcor_sim.root");
    plotPanel(
        TString::Format("#eta #mu vs OA (%s, SIM) — before/after ECAL correction", fl.c_str()).Data(),
        "OA(e^{+}e^{-}) [deg]",
        "#mu_{CB} [GeV/c^{2}]  (vertical bars = #pm#sigma_{CB})",
        eta_orig, eta_corr, kEtaPDG,
        "compare_mu_eta_" + fl + "_sim_vs_oa");
}
