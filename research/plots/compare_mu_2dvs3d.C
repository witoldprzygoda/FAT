// research/plots/compare_mu_2dvs3d.C — compare original / 2D-corrected /
// 3D-corrected μ_π⁰, μ_η vs OA on a single plot.
//
// Reads three sets of fit_results TTrees:
//   fit_results[_eta]_<fl>.root              (original, no correction)
//   fit_results[_eta]_<fl>_ecalcor.root      (corrected with 2D map)
//   fit_results[_eta]_<fl>_ecalcor_3dmap.root (corrected with 3D map)
//
// Usage (from research/):
//   root -l -b -q plots/compare_mu_2dvs3d.C            # REC, both pi0+eta

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
    if (!t) { f->Close(); return g; }
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

void plot3way(const std::string& title, double pdg_ref,
              const OAGraph& orig, const OAGraph& corr2d, const OAGraph& corr3d,
              const std::string& out_base)
{
    auto build = [](const OAGraph& g, int color, int marker) {
        auto* gg = new TGraphErrors((int)g.X.size(),
            g.X.data(), g.Y.data(), g.EX.data(), g.EY.data());
        gg->SetMarkerColor(color); gg->SetLineColor(color);
        gg->SetMarkerStyle(marker); gg->SetMarkerSize(0.9);
        return gg;
    };

    auto* g_o  = build(orig,   kBlack,    20);
    auto* g_2d = build(corr2d, kRed,      21);
    auto* g_3d = build(corr3d, kGreen + 2,22);

    auto* c = new TCanvas(("c_23_" + out_base).c_str(),
                          out_base.c_str(), 1100, 700);
    c->SetMargin(0.13, 0.05, 0.12, 0.08);
    c->SetGrid();
    g_o->SetTitle(TString::Format("%s;OA(e^{+}e^{-}) [deg];#mu_{CB} [GeV/c^{2}]  (bars = #pm#sigma_{CB})",
                                  title.c_str()));
    g_o ->Draw("AP");
    g_o ->GetXaxis()->SetLimits(0.0, 15.0);
    g_2d->Draw("P SAME");
    g_3d->Draw("P SAME");

    auto* lref = new TLine(0.0, pdg_ref, 15.0, pdg_ref);
    lref->SetLineColor(kBlue);
    lref->SetLineStyle(2);
    lref->SetLineWidth(2);
    lref->Draw();

    auto* leg = new TLegend(0.55, 0.70, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    leg->AddEntry(g_o,  "original (no correction)", "lpe");
    leg->AddEntry(g_2d, "corrected with 2D map (E,#theta)",   "lpe");
    leg->AddEntry(g_3d, "corrected with 3D map (E,#theta,#phi_{loc})", "lpe");
    leg->AddEntry(lref, TString::Format("PDG = %.5f", pdg_ref).Data(), "l");
    leg->Draw();

    c->Update();
    c->SaveAs(("plots/output/" + out_base + ".pdf").c_str());
    c->SaveAs(("plots/output/" + out_base + ".png").c_str());
    std::cout << "  wrote plots/output/" << out_base << ".{pdf,png}\n";
}

void compare_mu_2dvs3d(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    OAGraph pi0_o  = readMuVsOA("fit_results_" + fl + ".root");
    OAGraph pi0_2d = readMuVsOA("fit_results_" + fl + "_ecalcor.root");
    OAGraph pi0_3d = readMuVsOA("fit_results_" + fl + "_ecalcor_3dmap.root");
    plot3way("#pi^{0} #mu vs OA — 2D vs 3D ECAL map",
             kPi0PDG, pi0_o, pi0_2d, pi0_3d,
             "compare_mu_pi0_" + fl + "_2dvs3d_vs_oa");

    OAGraph eta_o  = readMuVsOA("fit_results_eta_" + fl + ".root");
    OAGraph eta_2d = readMuVsOA("fit_results_eta_" + fl + "_ecalcor.root");
    OAGraph eta_3d = readMuVsOA("fit_results_eta_" + fl + "_ecalcor_3dmap.root");
    plot3way("#eta #mu vs OA — 2D vs 3D ECAL map",
             kEtaPDG, eta_o, eta_2d, eta_3d,
             "compare_mu_eta_" + fl + "_2dvs3d_vs_oa");
}
