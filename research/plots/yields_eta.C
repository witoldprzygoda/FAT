// research/plots/yields_eta.C — derived plots from fit_results_eta_<fl>.root.
//
// Reads research/fit_results_eta_<flavour>.root produced by plots/fit_eta.C
// and makes:
//
//   1) yield_eta_<fl>_vs_oa.{pdf,png}
//      η signal counts vs opening angle, three series (μ ± 1σ, ±2σ, ±3σ).
//      X — slice center (deg);    X-bars = slice half-width = 0.1°.
//      Y — yield(signal − bg);    Y-bars = IntegralAndError statistical error.
//
//   2) mu_eta_<fl>_vs_oa.{pdf,png}
//      Crystal Ball peak position μ vs opening angle.
//      X — slice center (deg);    X-bars = slice half-width.
//      Y — μ from CB fit;         Y-bars = ±σ_CB (peak width).
//
// Usage (from research/):
//   root -l -b -q plots/yields_eta.C            # REC (default)
//   root -l -b -q 'plots/yields_eta.C("cor")'   # COR

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <vector>

void yields_eta(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    if (fl != "rec" && fl != "cor") {
        std::cerr << "Unknown flavour '" << flavour
                  << "' (expected 'rec' or 'cor')\n";
        return;
    }

    const std::string fpath = "fit_results_eta_" + fl + ".root";
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fpath << "\n"; return;
    }
    auto* t = (TTree*)f->Get("fit_results");
    if (!t) {
        std::cerr << "TTree 'fit_results' not in " << fpath << "\n"; return;
    }

    int    panel_idx = 0;
    float  oa_lo = 0, oa_hi = 0;
    float  mu = 0, sigma = 0;
    float  y1 = 0, e1 = 0, y2 = 0, e2 = 0, y3 = 0, e3 = 0;
    float  yf = 0, ef = 0;
    float  yam = 0, eam = 0;   // asymmetric [μ−5σ, μ+3σ]

    t->SetBranchAddress("panel_idx",           &panel_idx);
    t->SetBranchAddress("oa_lo",               &oa_lo);
    t->SetBranchAddress("oa_hi",               &oa_hi);
    t->SetBranchAddress("mu",                  &mu);
    t->SetBranchAddress("sigma",               &sigma);
    t->SetBranchAddress("yield_data_1sig",     &y1);
    t->SetBranchAddress("yield_data_1sig_err", &e1);
    t->SetBranchAddress("yield_data_2sig",     &y2);
    t->SetBranchAddress("yield_data_2sig_err", &e2);
    t->SetBranchAddress("yield_data_3sig",     &y3);
    t->SetBranchAddress("yield_data_3sig_err", &e3);
    t->SetBranchAddress("yield_data_full",     &yf);
    t->SetBranchAddress("yield_data_full_err", &ef);
    t->SetBranchAddress("yield_data_m5p3sig",     &yam);
    t->SetBranchAddress("yield_data_m5p3sig_err", &eam);

    std::vector<double> X, EX, Y1, EY1, Y2, EY2, Y3, EY3, YF, EYF,
                        YAM, EYAM, Mu, Sig;

    for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
        t->GetEntry(ev);
        if (panel_idx == 0) continue;

        X .push_back(0.5 * (oa_lo + oa_hi));
        EX.push_back(0.5 * (oa_hi - oa_lo));
        Y1.push_back(y1);   EY1.push_back(e1);
        Y2.push_back(y2);   EY2.push_back(e2);
        Y3.push_back(y3);   EY3.push_back(e3);
        YF.push_back(yf);   EYF.push_back(ef);
        YAM.push_back(yam); EYAM.push_back(eam);
        Mu.push_back(mu);   Sig.push_back(sigma);
    }

    const int N = (int)X.size();
    if (N == 0) { std::cerr << "No slice rows found in TTree.\n"; return; }

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    // --- Plot 1: yields vs OA (3 series) ------------------------------------
    auto* g1  = new TGraphErrors(N, X.data(), Y1.data(),  EX.data(), EY1.data());
    auto* g2  = new TGraphErrors(N, X.data(), Y2.data(),  EX.data(), EY2.data());
    auto* g3  = new TGraphErrors(N, X.data(), Y3.data(),  EX.data(), EY3.data());
    auto* gf  = new TGraphErrors(N, X.data(), YF.data(),  EX.data(), EYF.data());
    auto* gam = new TGraphErrors(N, X.data(), YAM.data(), EX.data(), EYAM.data());

    g1->SetMarkerColor(kBlue);      g1->SetLineColor(kBlue);
    g1->SetMarkerStyle(20);         g1->SetMarkerSize(0.9);
    g2->SetMarkerColor(kGreen + 2); g2->SetLineColor(kGreen + 2);
    g2->SetMarkerStyle(21);         g2->SetMarkerSize(0.9);
    g3->SetMarkerColor(kRed);       g3->SetLineColor(kRed);
    g3->SetMarkerStyle(22);         g3->SetMarkerSize(1.0);
    gf->SetMarkerColor(kBlack);     gf->SetLineColor(kBlack);
    gf->SetMarkerStyle(33);         gf->SetMarkerSize(1.2);
    gam->SetMarkerColor(kMagenta + 1); gam->SetLineColor(kMagenta + 1);
    gam->SetMarkerStyle(29);           gam->SetMarkerSize(1.2);

    auto* mg = new TMultiGraph();
    mg->Add(g1,  "P");
    mg->Add(g2,  "P");
    mg->Add(g3,  "P");
    mg->Add(gam, "P");
    mg->Add(gf,  "P");
    mg->SetTitle(TString::Format(
        "#eta signal yield vs opening angle (%s);"
        "OA(e^{+}e^{-}) [deg];yield (signal #minus bg)",
        fl.c_str()));

    auto* c1 = new TCanvas(("c_yields_eta_" + fl).c_str(), "yields vs OA", 1000, 700);
    c1->SetMargin(0.13, 0.05, 0.12, 0.08);
    c1->SetGrid();
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(0.0, 15.0);

    auto* leg1 = new TLegend(0.55, 0.62, 0.94, 0.90);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(g1,  "yield(#mu #pm 1#sigma)",                "lpe");
    leg1->AddEntry(g2,  "yield(#mu #pm 2#sigma)",                "lpe");
    leg1->AddEntry(g3,  "yield(#mu #pm 3#sigma)",                "lpe");
    leg1->AddEntry(gam, "yield([#mu#minus5#sigma, #mu+3#sigma])", "lpe");
    leg1->AddEntry(gf,  "yield(full fit range)",                 "lpe");
    leg1->Draw();

    c1->Update();
    const std::string out1 = "plots/output/yield_eta_" + fl + "_vs_oa";
    c1->SaveAs((out1 + ".pdf").c_str());
    c1->SaveAs((out1 + ".png").c_str());

    // --- Plot 2: μ vs OA (vertical bars = ±σ_CB) -----------------------------
    auto* gm = new TGraphErrors(N, X.data(), Mu.data(), EX.data(), Sig.data());
    gm->SetMarkerColor(kBlack);  gm->SetLineColor(kBlack);
    gm->SetMarkerStyle(20);      gm->SetMarkerSize(0.9);
    gm->SetTitle(TString::Format(
        "#eta peak position vs opening angle (%s);"
        "OA(e^{+}e^{-}) [deg];"
        "#mu_{CB} [GeV/c^{2}]  (vertical bars = #pm#sigma_{CB})",
        fl.c_str()));

    auto* c2 = new TCanvas(("c_mu_eta_" + fl).c_str(), "mu vs OA", 1000, 700);
    c2->SetMargin(0.13, 0.05, 0.12, 0.08);
    c2->SetGrid();
    gm->Draw("AP");
    gm->GetXaxis()->SetLimits(0.0, 15.0);

    constexpr double kEtaPDG = 0.5478;
    auto* lref = new TLine(0.0, kEtaPDG, 15.0, kEtaPDG);
    lref->SetLineColor(kRed);
    lref->SetLineStyle(2);
    lref->SetLineWidth(2);
    lref->Draw();

    auto* leg2 = new TLegend(0.62, 0.78, 0.94, 0.90);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(gm,   "fit #mu (bars = #pm#sigma_{CB})", "lpe");
    leg2->AddEntry(lref, TString::Format("PDG #eta = %.4f", kEtaPDG).Data(), "l");
    leg2->Draw();

    c2->Update();
    const std::string out2 = "plots/output/mu_eta_" + fl + "_vs_oa";
    c2->SaveAs((out2 + ".pdf").c_str());
    c2->SaveAs((out2 + ".png").c_str());

    f->Close();

    std::cout << "Wrote:\n"
              << "  " << out1 << ".{pdf,png}\n"
              << "  " << out2 << ".{pdf,png}\n"
              << "  ("<< N <<" OA slices)\n";
}
