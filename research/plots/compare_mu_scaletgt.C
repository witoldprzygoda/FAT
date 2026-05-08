// research/plots/compare_mu_scaletgt.C — compare original / 2D / 3D-map /
// scale_target map approaches on exp pi0+eta.

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

namespace { constexpr double kPi0PDG = 0.13498; constexpr double kEtaPDG = 0.5478; }

struct OAGraph { std::vector<double> X, EX, Y, EY; };

OAGraph readMuVsOA(const std::string& fpath) {
    OAGraph g;
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "missing " << fpath << "\n"; return g; }
    auto* t = (TTree*)f->Get("fit_results");
    if (!t) { f->Close(); return g; }
    int   panel_idx = 0;
    float oa_lo = 0, oa_hi = 0, mu = 0, sigma = 0;
    t->SetBranchAddress("panel_idx", &panel_idx);
    t->SetBranchAddress("oa_lo", &oa_lo);
    t->SetBranchAddress("oa_hi", &oa_hi);
    t->SetBranchAddress("mu", &mu);
    t->SetBranchAddress("sigma", &sigma);
    struct Row { double x, ex, y, ey; };
    std::vector<Row> rows;
    for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
        t->GetEntry(ev);
        if (panel_idx == 0) continue;
        rows.push_back({0.5*(oa_lo+oa_hi), 0.5*(oa_hi-oa_lo), mu, sigma});
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){ return a.x<b.x; });
    for (const auto& r : rows) { g.X.push_back(r.x); g.EX.push_back(r.ex); g.Y.push_back(r.y); g.EY.push_back(r.ey); }
    f->Close();
    return g;
}

void plot4(const std::string& title, double pdg,
           const OAGraph& orig, const OAGraph& c2d, const OAGraph& c3d, const OAGraph& cst,
           const std::string& out_base)
{
    auto build = [](const OAGraph& g, int color, int marker) {
        auto* gg = new TGraphErrors((int)g.X.size(),
            g.X.data(), g.Y.data(), g.EX.data(), g.EY.data());
        gg->SetMarkerColor(color); gg->SetLineColor(color);
        gg->SetMarkerStyle(marker); gg->SetMarkerSize(0.9);
        return gg;
    };
    auto* g_o   = build(orig, kBlack, 20);
    auto* g_2d  = build(c2d,  kRed,   21);
    auto* g_3d  = build(c3d,  kGreen+2, 22);
    auto* g_st  = build(cst,  kBlue,  23);

    auto* c = new TCanvas(("c_st_"+out_base).c_str(), out_base.c_str(), 1200, 700);
    c->SetMargin(0.13, 0.05, 0.12, 0.08);
    c->SetGrid();
    g_o->SetTitle(TString::Format("%s;OA(e^{+}e^{-}) [deg];#mu_{CB} [GeV/c^{2}]", title.c_str()));
    g_o ->Draw("AP");
    g_o ->GetXaxis()->SetLimits(0, 15);
    g_2d->Draw("P SAME");
    g_3d->Draw("P SAME");
    g_st->Draw("P SAME");
    auto* lref = new TLine(0, pdg, 15, pdg);
    lref->SetLineColor(kMagenta+2); lref->SetLineStyle(2); lref->SetLineWidth(2);
    lref->Draw();
    auto* leg = new TLegend(0.55, 0.66, 0.94, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
    leg->AddEntry(g_o,  "original",                              "lpe");
    leg->AddEntry(g_2d, "2D map (peak fit, E,#theta)",           "lpe");
    leg->AddEntry(g_3d, "3D map (peak fit, E,#theta,#phi)",      "lpe");
    leg->AddEntry(g_st, "scaletgt 3D (per-event #pi^{0}+#eta)",  "lpe");
    leg->AddEntry(lref, TString::Format("PDG = %.5f", pdg).Data(), "l");
    leg->Draw();
    c->Update();
    c->SaveAs(("plots/output/"+out_base+".pdf").c_str());
    c->SaveAs(("plots/output/"+out_base+".png").c_str());
    std::cout << "  wrote plots/output/" << out_base << ".{pdf,png}\n";
}

void compare_mu_scaletgt(const char* flavour = "rec") {
    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    // π⁰
    plot4("#pi^{0} #mu vs OA — 4 correction approaches", kPi0PDG,
          readMuVsOA("fit_results_" + fl + ".root"),
          readMuVsOA("fit_results_" + fl + "_ecalcor.root"),
          readMuVsOA("fit_results_" + fl + "_ecalcor_3dmap.root"),
          readMuVsOA("fit_results_" + fl + "_ecalcor_scaletgt.root"),
          "compare_mu_pi0_" + fl + "_scaletgt");

    // η
    plot4("#eta #mu vs OA — 4 correction approaches", kEtaPDG,
          readMuVsOA("fit_results_eta_" + fl + ".root"),
          readMuVsOA("fit_results_eta_" + fl + "_ecalcor.root"),
          readMuVsOA("fit_results_eta_" + fl + "_ecalcor_3dmap.root"),
          readMuVsOA("fit_results_eta_" + fl + "_ecalcor_scaletgt.root"),
          "compare_mu_eta_" + fl + "_scaletgt");
}
