// x17p26_proj.C -- X17 PDG 2026 (ruling D25), prefix x17p26, 2026-11-15: the pre-data projection of x17ab (x17ab_proj_mass.pdf) redrawn
// in PDG 2026 units (x 2.76/2.68; a pure conversion of r -> BR, no refit), with the drawing code of x17ab_model.C block (c) copied
// verbatim except the factor, one text line and the output name x17ab_proj_mass_pdg26. Input: docs/latex/tables_all/x17ab_proj_mass.txt.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17p26", kTRUE)' 'plots/all_report/x17p26_proj.C+O'
#include <cstdio>
#include <vector>
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
void x17p26_proj() {
    gStyle->SetOptStat(0);
    const double FBR = 2.76 / 2.68, kAproj = 0.40;
    const char* kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
    std::vector<double> mgrid, ulm[5];
    FILE* fp = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/x17ab_proj_mass.txt", "r");
    if (!fp) { printf("x17p26_proj: input missing\n"); return; }
    char line[512];
    while (fgets(line, sizeof line, fp)) {
        double m, rx, a[5];
        if (line[0] != '#' && sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &m, &rx, &a[0], &a[1], &a[2], &a[3], &a[4]) == 7) {
            mgrid.push_back(m);
            for (int q = 0; q < 5; ++q) ulm[q].push_back(FBR * a[q]);
        }
    }
    fclose(fp);
    TCanvas c("cp", "", 900, 700);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12); gPad->SetLogy();
    const int n = (int)mgrid.size();
    TGraphAsymmErrors b2(n), b1(n);
    TGraph med(n);
    for (int k = 0; k < n; ++k) {
        b2.SetPoint(k, mgrid[k], ulm[2][k]); b2.SetPointError(k, 0.5, 0.5, ulm[2][k] - ulm[0][k], ulm[4][k] - ulm[2][k]);
        b1.SetPoint(k, mgrid[k], ulm[2][k]); b1.SetPointError(k, 0.5, 0.5, ulm[2][k] - ulm[1][k], ulm[3][k] - ulm[2][k]);
        med.SetPoint(k, mgrid[k], ulm[2][k]);
        if (mgrid[k] > 16.95 && mgrid[k] < 17.05) printf("x17p26_proj: 17 MeV/c2: %.4e [%.4e, %.4e] [%.4e, %.4e]\n", ulm[2][k], ulm[1][k], ulm[3][k], ulm[0][k], ulm[4][k]);
    }
    b2.SetFillColor(kOrange); b1.SetFillColor(kGreen + 1);
    b2.SetTitle(";m_{X} [MeV/c^{2}];90% CL_{s} on BR(#eta#rightarrow#pi^{+}#pi^{-}X17) BR(X17#rightarrowe^{+}e^{-})");
    b2.GetXaxis()->SetTitleSize(0.045); b2.GetYaxis()->SetTitleSize(0.04);
    b2.SetMinimum(3e-6); b2.SetMaximum(5e-4);
    b2.Draw("A3"); b1.Draw("3 SAME");
    med.SetLineStyle(2); med.SetLineWidth(2); med.Draw("L SAME");
    TLegend lg(0.42, 0.74, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.033); lg.SetFillStyle(0);
    lg.AddEntry(&med, "expected (simulation Asimov, pre-data)", "l"); lg.AddEntry(&b1, "#pm 1#sigma", "f"); lg.AddEntry(&b2, "#pm 2#sigma", "f");
    lg.Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
    tx.DrawLatex(0.16, 0.17, Form("S0, X2D; a_{#eta} = %.2f, background 1.35 #times SMASH, data CB; L = 6079.2 nb^{-1}", kAproj));
    tx.DrawLatex(0.16, 0.22, "BR(#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}) = 2.76#times10^{-4} (PDG 2026)");
    c.SaveAs(Form("%s/x17ab_proj_mass_pdg26.pdf", kFig)); c.SaveAs(Form("%s/x17ab_proj_mass_pdg26.png", kFig));
}
