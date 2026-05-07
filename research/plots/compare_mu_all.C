// research/plots/compare_mu_all.C — overlay μ vs OA for exp + sim,
// before and after ECAL energy-scale correction.
//
// Reads four fit_results TTrees per meson:
//   fit_results[_eta]_<fl>.root              (exp original)
//   fit_results[_eta]_<fl>_ecalcor.root      (exp corrected)
//   fit_results[_eta]_<fl>_sim.root          (sim original)
//   fit_results[_eta]_<fl>_ecalcor_sim.root  (sim corrected)
//
// Produces side-by-side overlays:
//   plots/output/compare_mu_pi0_<fl>_all.{pdf,png}
//   plots/output/compare_mu_eta_<fl>_all.{pdf,png}
//
// Usage (from research/):
//   root -l -b -q plots/compare_mu_all.C            # REC, both pi0 + eta

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

void plotPanelAll(const std::string& title, const std::string& xtitle,
                  const std::string& ytitle, double pdg_ref,
                  const OAGraph& exp_orig, const OAGraph& exp_corr,
                  const OAGraph& sim_orig, const OAGraph& sim_corr,
                  const std::string& out_base)
{
    auto build = [](const OAGraph& g, int color, int marker) {
        auto* gg = new TGraphErrors((int)g.X.size(),
            g.X.data(), g.Y.data(),
            g.EX.data(), g.EY.data());
        gg->SetMarkerColor(color); gg->SetLineColor(color);
        gg->SetMarkerStyle(marker); gg->SetMarkerSize(0.9);
        return gg;
    };

    auto* g_exp_o = build(exp_orig, kBlack,    20);   // ●
    auto* g_exp_c = build(exp_corr, kRed,      21);   // ■
    auto* g_sim_o = build(sim_orig, kBlue + 2, 24);   // ○
    auto* g_sim_c = build(sim_corr, kGreen + 2,25);   // □

    const TString ctitle = TString::Format("%s;%s;%s",
        title.c_str(), xtitle.c_str(), ytitle.c_str());

    auto* c = new TCanvas(("c_all_" + out_base).c_str(),
                          out_base.c_str(), 1200, 750);
    c->SetMargin(0.13, 0.05, 0.12, 0.08);
    c->SetGrid();
    g_exp_o->SetTitle(ctitle);
    g_exp_o->Draw("AP");
    g_exp_o->GetXaxis()->SetLimits(0.0, 15.0);
    g_exp_c->Draw("P SAME");
    g_sim_o->Draw("P SAME");
    g_sim_c->Draw("P SAME");

    auto* lref = new TLine(0.0, pdg_ref, 15.0, pdg_ref);
    lref->SetLineColor(kMagenta + 2);
    lref->SetLineStyle(2);
    lref->SetLineWidth(2);
    lref->Draw();

    auto* leg = new TLegend(0.55, 0.62, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(g_exp_o, "exp original",     "lpe");
    leg->AddEntry(g_exp_c, "exp ECAL-corrected","lpe");
    leg->AddEntry(g_sim_o, "sim original",     "lpe");
    leg->AddEntry(g_sim_c, "sim ECAL-corrected","lpe");
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

void compare_mu_all(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    // -- π⁰ ---------------------------------------------------------------
    OAGraph pi0_exp_o = readMuVsOA("fit_results_" + fl + ".root");
    OAGraph pi0_exp_c = readMuVsOA("fit_results_" + fl + "_ecalcor.root");
    OAGraph pi0_sim_o = readMuVsOA("fit_results_" + fl + "_sim.root");
    OAGraph pi0_sim_c = readMuVsOA("fit_results_" + fl + "_ecalcor_sim.root");
    plotPanelAll(
        TString::Format("#pi^{0} #mu vs OA (%s) — exp + sim, before/after ECAL correction", fl.c_str()).Data(),
        "OA(e^{+}e^{-}) [deg]",
        "#mu_{CB} [GeV/c^{2}]  (vertical bars = #pm#sigma_{CB})",
        kPi0PDG, pi0_exp_o, pi0_exp_c, pi0_sim_o, pi0_sim_c,
        "compare_mu_pi0_" + fl + "_all_vs_oa");

    // -- η ----------------------------------------------------------------
    OAGraph eta_exp_o = readMuVsOA("fit_results_eta_" + fl + ".root");
    OAGraph eta_exp_c = readMuVsOA("fit_results_eta_" + fl + "_ecalcor.root");
    OAGraph eta_sim_o = readMuVsOA("fit_results_eta_" + fl + "_sim.root");
    OAGraph eta_sim_c = readMuVsOA("fit_results_eta_" + fl + "_ecalcor_sim.root");
    plotPanelAll(
        TString::Format("#eta #mu vs OA (%s) — exp + sim, before/after ECAL correction", fl.c_str()).Data(),
        "OA(e^{+}e^{-}) [deg]",
        "#mu_{CB} [GeV/c^{2}]  (vertical bars = #pm#sigma_{CB})",
        kEtaPDG, eta_exp_o, eta_exp_c, eta_sim_o, eta_sim_c,
        "compare_mu_eta_" + fl + "_all_vs_oa");
}
