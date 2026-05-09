// Classic μ vs OA plots after the combined π⁰+η 3D calibration.
// Four plots: π⁰/η × exp/sim, OA range 0–15°.

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

static OAGraph readMuVsOA(const std::string& fpath) {
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

static void plot1(const std::string& title, double pdg,
                  const OAGraph& g,
                  double yLo, double yHi,
                  const std::string& out)
{
    auto* gg = new TGraphErrors((int)g.X.size(),
        g.X.data(), g.Y.data(), g.EX.data(), g.EY.data());
    gg->SetMarkerColor(kBlack); gg->SetLineColor(kBlack);
    gg->SetMarkerStyle(20); gg->SetMarkerSize(0.9);
    gg->SetTitle(TString::Format(
        "%s;OA(e^{+}e^{-}) [deg];#mu_{CB} [GeV/c^{2}]  (bars = #pm#sigma_{CB})",
        title.c_str()));

    auto* cv = new TCanvas(("c_"+out).c_str(), out.c_str(), 1100, 700);
    cv->SetMargin(0.13, 0.05, 0.12, 0.08); cv->SetGrid();
    gg->Draw("AP");
    gg->GetXaxis()->SetLimits(0.0, 15.0);
    gg->GetYaxis()->SetRangeUser(yLo, yHi);

    auto* lref = new TLine(0.0, pdg, 15.0, pdg);
    lref->SetLineColor(kMagenta+2); lref->SetLineStyle(2); lref->SetLineWidth(2);
    lref->Draw();

    auto* leg = new TLegend(0.55, 0.78, 0.94, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
    leg->AddEntry(gg,   "fit #mu (bars = #pm#sigma_{CB})", "lpe");
    leg->AddEntry(lref, TString::Format("PDG = %.5f", pdg).Data(), "l");
    leg->Draw();

    cv->Update();
    cv->SaveAs(("plots/output/"+out+".pdf").c_str());
    cv->SaveAs(("plots/output/"+out+".png").c_str());
    std::cout << "  wrote " << out << ".{pdf,png}\n";
}

void mu_combined_4plots() {
    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    plot1("#pi^{0} #mu vs OA (EXP, 3D combined #pi^{0}+#eta)", kPi0PDG,
          readMuVsOA("fit_results_rec_ecalcor_combined.root"),
          0.115, 0.165,
          "mu_pi0_rec_combined_exp_vs_oa");

    plot1("#pi^{0} #mu vs OA (SIM, 3D combined #pi^{0}+#eta)", kPi0PDG,
          readMuVsOA("fit_results_rec_ecalcor_combined_sim.root"),
          0.115, 0.165,
          "mu_pi0_rec_combined_sim_vs_oa");

    plot1("#eta #mu vs OA (EXP, 3D combined #pi^{0}+#eta)", kEtaPDG,
          readMuVsOA("fit_results_eta_rec_ecalcor_combined.root"),
          0.49, 0.62,
          "mu_eta_rec_combined_exp_vs_oa");

    plot1("#eta #mu vs OA (SIM, 3D combined #pi^{0}+#eta)", kEtaPDG,
          readMuVsOA("fit_results_eta_rec_ecalcor_combined_sim.root"),
          0.49, 0.62,
          "mu_eta_rec_combined_sim_vs_oa");
}
