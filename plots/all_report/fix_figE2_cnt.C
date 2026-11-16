// fix_figE2_cnt.C -- COPY of fix_figE2.C (figure agent B, 2026-11-15): y title in "counts" instead of "rows"; output fix_E2_cnt.{pdf,png}.
// Run (from /home/przygoda/HADES/FAT/FAT): root -l -b -q 'plots/all_report/fix_figE2_cnt.C'
// ============================================================================
// fix_figE2.C -- review fix (prefix fix): v5_figE2.C (verifier v5) redrawn in the common palette of report_all
// (C2 orange, C3+C4 magenta (conversions), C5 cyan, G_eta violet, as Fig. B1); only colours, names and output paths changed.
// ============================================================================
// Same content as docs/latex/figs_all/wp5_E2.pdf (WP5 nominal F-main, E2): per slice, data N - k CB - fitted non-eta
// background (blue, errors = event-group variance of N + CB replica variance), the fitted eta classes stacked, and all eta
// classes at SMASH x PDG BR (dashed). The ONLY change: in WP5's version the white-filled legend in the S5 panel covered the
// S5 data in M4 0.40-0.57 (where the residuals are negative, cf. the S5 pulls of Fig. E1); here the legend sits in the
// empty upper-right part of the "all slices" panel and is transparent.
// Input: /hdd3/przygoda/HADES/FAT_out/all_report/v5/v5_E2_input.txt (v5_figE2_prep.py, from WP5 results.txt + hists_v1.root)
// Output: docs/latex/figs_all/v5_wp5_E2_corrected.{pdf,png}
// Run (from /home/przygoda/HADES/FAT/FAT): root -l -b -q 'plots/all_report/v5_figE2.C'
// ============================================================================
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

void fix_figE2_cnt() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.050, "XYZ"); gStyle->SetTitleSize(0.055, "XYZ"); gStyle->SetLegendBorderSize(0);
    gStyle->SetEndErrorSize(0); gStyle->SetFrameLineWidth(2); gStyle->SetHistLineWidth(2);
    std::vector<double> e;
    for (int i = 0; i <= 24; ++i) e.push_back(0.400 + 0.005 * i);
    for (int i = 1; i <= 20; ++i) e.push_back(0.520 + 0.0025 * i);
    for (int i = 1; i <= 26; ++i) e.push_back(0.570 + 0.005 * i);
    const int colEta = kViolet - 9, colC2 = kOrange + 7, colC34 = kMagenta + 1, colC5 = kCyan + 1;
    // histograms per panel (0..4 slices, 5 = sum)
    TH1D *hD[6], *h2[6], *h34[6], *h5[6], *hE[6], *hS[6];
    for (int p = 0; p < 6; ++p) {
        hD[p] = new TH1D(Form("D%d", p), "", 70, e.data()); h2[p] = new TH1D(Form("c2%d", p), "", 70, e.data());
        h34[p] = new TH1D(Form("c34%d", p), "", 70, e.data()); h5[p] = new TH1D(Form("c5%d", p), "", 70, e.data());
        hE[p] = new TH1D(Form("ge%d", p), "", 70, e.data()); hS[p] = new TH1D(Form("sm%d", p), "", 70, e.data());
    }
    std::vector<double> var[6];
    for (int p = 0; p < 6; ++p) var[p].assign(71, 0.);
    std::ifstream in("/hdd3/przygoda/HADES/FAT_out/all_report/v5/v5_E2_input.txt");
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int j, x; double lo, hi, D, vD, c2, c34, c5, ge, sm;
        ss >> j >> x >> lo >> hi >> D >> vD >> c2 >> c34 >> c5 >> ge >> sm;
        const double w = (hi - lo) * 1000.;
        for (int p : {j, 5}) {
            hD[p]->AddBinContent(x, D / w); var[p][x] += vD / (w * w);
            h2[p]->AddBinContent(x, c2 / w); h34[p]->AddBinContent(x, c34 / w); h5[p]->AddBinContent(x, c5 / w);
            hE[p]->AddBinContent(x, ge / w); hS[p]->AddBinContent(x, sm / w);
        }
    }
    TCanvas* c = new TCanvas("cE2fix", "", 2100, 1300);
    c->Divide(3, 2, 0.002, 0.002);
    for (int p = 0; p < 6; ++p) {
        c->cd(p + 1);
        gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.07);
        for (int x = 1; x <= 70; ++x) hD[p]->SetBinError(x, std::sqrt(var[p][x]));
        double ymax = 0, ymin = 0;
        for (int x = 1; x <= 70; ++x) {
            ymax = std::max(ymax, hD[p]->GetBinContent(x) + hD[p]->GetBinError(x)); ymin = std::min(ymin, hD[p]->GetBinContent(x) - hD[p]->GetBinError(x));
            ymax = std::max(ymax, hS[p]->GetBinContent(x));
        }
        TH1D* fr = new TH1D(Form("fr%d", p), "", 70, e.data());
        fr->SetMaximum(1.15 * ymax); fr->SetMinimum(1.15 * ymin);
        fr->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); fr->GetYaxis()->SetTitle("N - k CB - background [counts/(MeV/c^{2})]");
        fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleOffset(1.6); fr->GetXaxis()->SetTitleSize(0.052);
        fr->Draw("AXIS");
        // manual stack (cumulative), drawn from the top down
        TH1D* comp[4] = {hE[p], h34[p], h2[p], h5[p]};
        const int col[4] = {colEta, colC34, colC2, colC5};
        TH1D* cum[4];
        for (int q = 0; q < 4; ++q) {
            cum[q] = (TH1D*)comp[q]->Clone(Form("cum%d_%d", p, q));
            if (q > 0) cum[q]->Add(cum[q - 1]);
            cum[q]->SetFillColor(col[q]); cum[q]->SetLineColor(col[q]); cum[q]->SetLineWidth(1); cum[q]->SetFillStyle(1001);
        }
        for (int q = 3; q >= 0; --q) cum[q]->Draw("HIST SAME");
        hS[p]->SetLineColor(kBlack); hS[p]->SetLineStyle(2); hS[p]->SetLineWidth(2); hS[p]->Draw("HIST SAME");
        hD[p]->SetMarkerStyle(20); hD[p]->SetMarkerSize(0.8); hD[p]->SetMarkerColor(kBlue + 1); hD[p]->SetLineColor(kBlue + 1); hD[p]->Draw("E1 SAME");
        TLine* l0 = new TLine(0.40, 0, 0.70, 0); l0->SetLineStyle(3); l0->Draw();
        fr->Draw("AXIS SAME");
        TLatex lt; lt.SetNDC(); lt.SetTextSize(0.05);
        lt.DrawLatex(0.21, 0.88, p < 5 ? Form("S%d, M_{ee} #geq 10 MeV/c^{2}", p + 1) : "all slices, M_{ee} #geq 10 MeV/c^{2}");
        if (p == 4) {
            lt.SetTextSize(0.040);
            lt.DrawLatex(0.21, 0.80, "a_{#eta} = 0.382, #kappa = 0.83, #rho_{d} = 1.32, g = 1.72");
            lt.DrawLatex(0.21, 0.74, "CROSS-CHECK (F-main failed WP4)");
        }
        if (p == 5) {
            TLegend* lg = new TLegend(0.50, 0.58, 0.96, 0.87); lg->SetTextSize(0.034); lg->SetFillStyle(0);
            lg->AddEntry(hD[p], "data N - k CB - fitted bkg.", "pe");
            lg->AddEntry(cum[3], "C5 (a_{#eta} #rho_{d})", "f");
            lg->AddEntry(cum[2], "C2 (a_{#eta})", "f");
            lg->AddEntry(cum[1], "C3 + C4, conversions (a_{#eta} #kappa)", "f");
            lg->AddEntry(cum[0], "G_{#eta} (a_{#eta} g)", "f");
            lg->AddEntry(hS[p], "all #eta, SMASH #times PDG BR", "l");
            lg->Draw();
        }
    }
    c->SaveAs("/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/fix_E2_cnt.pdf");
    c->SaveAs("/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/fix_E2_cnt.png");
}
