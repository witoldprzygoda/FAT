// k_factor_mee.C
//
// The like-sign acceptance correction k versus the e+e- invariant mass, from
// existing outputs (all pair combinations, no isBest, no opening-angle cut):
//   k_mix  (exp data): N+-_mix / 2 sqrt(N++_mix N--_mix), event mixing without
//                      multiplicity classes (h_mixnc_* of mixing_betanew.root);
//   k_true (January SMASH sim, MC truth): B_truth / 2 sqrt(N++ N--) from
//                      cb_closure_sim_jan26.root, for true-lepton pairs only and
//                      for all pairs including fake legs.
// The smooth curve is a spline through k_mix on coarse knots, drawn to guide the
// eye. It is NOT the final k = k_mix,data * k_true / k_mix,sim: k_mix,sim needs
// a simulation production with like-sign trees and enough statistics.
//
// Run: root -l -b -q 'plots/k_factor_mee.C+'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSpline.h>
#include <TStyle.h>
#include <cmath>
#include <cstdio>
#include <vector>

namespace kfm {
// k on the bins [e[i], e[i+1]) from three spectra on a finer common grid
TGraphErrors* makeKGraph(TH1* n, TH1* pp, TH1* mm, const std::vector<double>& e, double minN) {
    TGraphErrors* g = new TGraphErrors();
    for (size_t i = 0; i + 1 < e.size(); ++i) {
        const int b1 = n->FindBin(e[i] + 1e-6), b2 = n->FindBin(e[i + 1] - 1e-6);
        const double N = n->Integral(b1, b2), P = pp->Integral(b1, b2), M = mm->Integral(b1, b2);
        if (N < minN || P <= 0 || M <= 0) continue;
        const double k = N / (2 * std::sqrt(P * M));
        const int ip = g->GetN();
        g->SetPoint(ip, 0.5 * (e[i] + e[i + 1]), k);
        g->SetPointError(ip, 0.5 * (e[i + 1] - e[i]), k * std::sqrt(1 / N + 0.25 / P + 0.25 / M));
    }
    return g;
}
void frameStyle(TH1D* f) {
    f->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]"); f->GetYaxis()->SetTitle("k = N_{+-} / 2#sqrt{N_{++}N_{--}}");
    f->GetXaxis()->SetTitleSize(0.055); f->GetYaxis()->SetTitleSize(0.055);
    f->GetXaxis()->SetLabelSize(0.05); f->GetYaxis()->SetLabelSize(0.05);
    f->GetYaxis()->SetTitleOffset(1.15); f->GetXaxis()->SetTitleOffset(1.1);
}
}

void k_factor_mee(const char* mix_file = "plots/output/mixing_betanew.root",
                  const char* sim_file = "plots/output/cb_closure_sim_jan26.root",
                  const char* out = "plots/output/k_factor_mee")
{
    using namespace kfm;
    gStyle->SetOptStat(0);
    TFile fm(mix_file), fs(sim_file);
    TH1* mpm = (TH1*)fm.Get("h_mixnc_pm_PT3"); TH1* mpp = (TH1*)fm.Get("h_mixnc_pp_PT3"); TH1* mmm = (TH1*)fm.Get("h_mixnc_mm_PT3");
    TH1* bAll = (TH1*)fs.Get("Btruth_all"); TH1* nPP = (TH1*)fs.Get("Npp_all"); TH1* nMM = (TH1*)fs.Get("Nmm_all");
    TH1* bLep = (TH1*)fs.Get("h_all_pm_c3"); TH1* pLep = (TH1*)fs.Get("h_all_pp_c3"); TH1* mLep = (TH1*)fs.Get("h_all_mm_c3");
    if (!mpm || !mpp || !mmm || !bAll || !nPP || !nMM || !bLep || !pLep || !mLep) { printf("missing histograms\n"); return; }

    // analysis grid for the data, coarser for the simulation (statistics)
    std::vector<double> gridD; for (int i = 0; i <= 20; ++i) gridD.push_back(i * 0.0025);
    for (double x : {0.08, 0.11, 0.14}) gridD.push_back(x);
    for (int i = 1; i <= 33; ++i) gridD.push_back(0.14 + i * 0.02);
    gridD.push_back(0.95); gridD.push_back(1.40);
    const std::vector<double> gridS = {0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.11, 0.14, 0.18, 0.22, 0.30, 0.40, 0.50, 0.60, 0.80};
    const std::vector<double> knots = {0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.11, 0.14, 0.18, 0.22, 0.30, 0.40, 0.50, 0.60, 0.80, 0.95, 1.40};

    TGraphErrors* gMix = makeKGraph(mpm, mpp, mmm, gridD, 100);
    TGraphErrors* gLep = makeKGraph(bLep, pLep, mLep, gridS, 20);
    TGraphErrors* gAll = makeKGraph(bAll, nPP, nMM, gridS, 20);
    TGraphErrors* gKn  = makeKGraph(mpm, mpp, mmm, knots, 100);
    TSpline3* sp = new TSpline3("sp", gKn);

    auto dump = [](const char* t, TGraphErrors* g) { printf("%s:", t); for (int i = 0; i < g->GetN(); ++i) printf(" %.3f:%.3f(%.3f)", g->GetX()[i], g->GetY()[i], g->GetEY()[i]); printf("\n"); };
    dump("k_mix knots", gKn); dump("k_true leptons", gLep); dump("k_true all", gAll);

    gMix->SetMarkerStyle(20); gMix->SetMarkerSize(0.7); gMix->SetMarkerColor(kBlack); gMix->SetLineColor(kBlack);
    gLep->SetMarkerStyle(21); gLep->SetMarkerSize(1.0); gLep->SetMarkerColor(kBlue + 1); gLep->SetLineColor(kBlue + 1);
    gAll->SetMarkerStyle(25); gAll->SetMarkerSize(1.0); gAll->SetMarkerColor(kRed + 1); gAll->SetLineColor(kRed + 1);
    sp->SetLineColor(kGray + 2); sp->SetLineWidth(2);

    TCanvas* c = new TCanvas("c", "", 1500, 620);
    c->Divide(2, 1);
    const double xmax[2] = {1.0, 0.20};
    const double ymin[2] = {0.6, 0.6}, ymax[2] = {1.8, 1.8};
    for (int p = 0; p < 2; ++p) {
        c->cd(p + 1);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetGridy();
        TH1D* f = new TH1D(Form("f%d", p), "", 10, 0., xmax[p]);
        f->SetMinimum(ymin[p]); f->SetMaximum(ymax[p]); frameStyle(f); f->Draw();
        TLine* l = new TLine(0, 1, xmax[p], 1); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
        sp->Draw("same");
        gAll->Draw("P same"); gLep->Draw("P same"); gMix->Draw("P same");
        if (p == 0) {
            TLegend* lg = new TLegend(0.36, 0.68, 0.97, 0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
            lg->AddEntry(gMix, "k_{mix}: event mixing, data", "pe");
            lg->AddEntry(gLep, "k_{true}: simulation, true leptons", "pe");
            lg->AddEntry(gAll, "k_{true}: simulation, with fake legs", "pe");
            lg->AddEntry(sp, "spline through k_{mix} (guide)", "l");
            lg->Draw();
        }
        TLatex t; t.SetNDC(); t.SetTextSize(0.045);
        t.DrawLatex(0.13, 0.935, p == 0 ? "all pair combinations, no OA cut" : "low-mass zoom");
    }
    c->SaveAs(Form("%s.pdf", out));
    c->SaveAs(Form("%s.png", out));
}
