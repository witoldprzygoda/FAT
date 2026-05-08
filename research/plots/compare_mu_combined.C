// Compare orig vs 2D vs 3D-pi0 vs combined-pi0+eta, on exp.

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
#include <vector>
#include <algorithm>

namespace { constexpr double kPi0PDG = 0.13498; constexpr double kEtaPDG = 0.5478; }

struct OAGraph { std::vector<double> X, EX, Y, EY; };

OAGraph readMuVsOA(const std::string& fpath) {
    OAGraph g;
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr<<"missing "<<fpath<<"\n"; return g; }
    auto* t = (TTree*)f->Get("fit_results");
    if (!t) { f->Close(); return g; }
    int   panel_idx = 0;
    float oa_lo=0, oa_hi=0, mu=0, sigma=0;
    t->SetBranchAddress("panel_idx", &panel_idx);
    t->SetBranchAddress("oa_lo", &oa_lo); t->SetBranchAddress("oa_hi", &oa_hi);
    t->SetBranchAddress("mu", &mu); t->SetBranchAddress("sigma", &sigma);
    struct Row { double x,ex,y,ey; };
    std::vector<Row> rows;
    for (Long64_t ev=0; ev<t->GetEntries(); ++ev) {
        t->GetEntry(ev);
        if (panel_idx == 0) continue;
        rows.push_back({0.5*(oa_lo+oa_hi), 0.5*(oa_hi-oa_lo), mu, sigma});
    }
    std::sort(rows.begin(), rows.end(), [](const Row&a, const Row&b){return a.x<b.x;});
    for (const auto& r : rows) { g.X.push_back(r.x); g.EX.push_back(r.ex); g.Y.push_back(r.y); g.EY.push_back(r.ey); }
    f->Close();
    return g;
}

void plot4(const std::string& title, double pdg,
           const OAGraph& a, const OAGraph& b, const OAGraph& c, const OAGraph& d,
           const std::string& la, const std::string& lb, const std::string& lc, const std::string& ld,
           const std::string& out)
{
    auto build = [](const OAGraph& g, int color, int marker) {
        auto* gg = new TGraphErrors((int)g.X.size(),
            g.X.data(), g.Y.data(), g.EX.data(), g.EY.data());
        gg->SetMarkerColor(color); gg->SetLineColor(color);
        gg->SetMarkerStyle(marker); gg->SetMarkerSize(0.9);
        return gg;
    };
    auto* ga = build(a, kBlack,    20);
    auto* gb = build(b, kRed,      21);
    auto* gc = build(c, kGreen+2,  22);
    auto* gd = build(d, kBlue,     23);

    auto* cv = new TCanvas(("c_"+out).c_str(), out.c_str(), 1200, 700);
    cv->SetMargin(0.13, 0.05, 0.12, 0.08); cv->SetGrid();
    ga->SetTitle(TString::Format("%s;OA(e^{+}e^{-}) [deg];#mu_{CB} [GeV/c^{2}]", title.c_str()));
    ga->Draw("AP");
    ga->GetXaxis()->SetLimits(0, 15);
    gb->Draw("P SAME"); gc->Draw("P SAME"); gd->Draw("P SAME");

    auto* lref = new TLine(0, pdg, 15, pdg);
    lref->SetLineColor(kMagenta+2); lref->SetLineStyle(2); lref->SetLineWidth(2);
    lref->Draw();

    auto* leg = new TLegend(0.55, 0.66, 0.94, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
    leg->AddEntry(ga, la.c_str(), "lpe");
    leg->AddEntry(gb, lb.c_str(), "lpe");
    leg->AddEntry(gc, lc.c_str(), "lpe");
    leg->AddEntry(gd, ld.c_str(), "lpe");
    leg->AddEntry(lref, TString::Format("PDG = %.5f", pdg).Data(), "l");
    leg->Draw();

    cv->Update();
    cv->SaveAs(("plots/output/"+out+".pdf").c_str());
    cv->SaveAs(("plots/output/"+out+".png").c_str());
    std::cout << "  wrote " << out << ".{pdf,png}\n";
}

void compare_mu_combined(const char* flavour = "rec") {
    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    plot4("#pi^{0} #mu vs OA — 4 approaches (SIM)", kPi0PDG,
          readMuVsOA("fit_results_" + fl + "_sim.root"),
          readMuVsOA("fit_results_" + fl + "_ecalcor_sim.root"),
          readMuVsOA("fit_results_" + fl + "_ecalcor_3dmap_sim.root"),
          readMuVsOA("fit_results_" + fl + "_ecalcor_combined_sim.root"),
          "sim original", "sim 2D peak-fit (E,#theta)",
          "sim 3D peak-fit #pi^{0} (E,#theta,#phi)",
          "sim 3D combined #pi^{0}+#eta",
          "compare_mu_pi0_" + fl + "_combined_sim");

    plot4("#eta #mu vs OA — 4 approaches (SIM)", kEtaPDG,
          readMuVsOA("fit_results_eta_" + fl + "_sim.root"),
          readMuVsOA("fit_results_eta_" + fl + "_ecalcor_sim.root"),
          readMuVsOA("fit_results_eta_" + fl + "_ecalcor_3dmap_sim.root"),
          readMuVsOA("fit_results_eta_" + fl + "_ecalcor_combined_sim.root"),
          "sim original", "sim 2D peak-fit (E,#theta)",
          "sim 3D peak-fit #pi^{0} (E,#theta,#phi)",
          "sim 3D combined #pi^{0}+#eta",
          "compare_mu_eta_" + fl + "_combined_sim");
}
