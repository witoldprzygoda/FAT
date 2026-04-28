// pippim_masses_cor.C — Three CORRECTED-kinematics mass spectra:
//   M(gamma gamma)        — mass_gg_cor
//   M(pi+ pi-)            — mass_pippim_cor
//   M(pi+ pi- gamma gamma) — mass_pippimgg_cor
//
// Hadronic channel — no combinatorial background, single output file.
//
// Output: plots/output/{m_gg_cor,m_pippim_cor,m_pippimgg_cor}.{pdf,png}
//
// Usage:
//   root -l -b -q plots/pippim_masses_cor.C                     # default file
//   root -l -b -q 'plots/pippim_masses_cor.C("/path/file.root")'

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
#include <iostream>
#include <string>

namespace {
    void applyGlobalStyle() {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetTitleSize(0.055, "XYZ");
        gStyle->SetTitleOffset(1.05, "X");
        gStyle->SetTitleOffset(1.20, "Y");
        gStyle->SetLabelSize(0.045, "XYZ");
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetPadLeftMargin(0.135);
        gStyle->SetPadRightMargin(0.04);
        gStyle->SetPadTopMargin(0.06);
        gStyle->SetPadBottomMargin(0.13);
        // Subtle vertical grid: thin dotted gray instead of default solid black
        gStyle->SetGridColor(kGray + 1);
        gStyle->SetGridStyle(3);   // 3 = dotted
        gStyle->SetGridWidth(1);
    }

    // Compose nice axis title from a particle-mass label (uses LaTeX).
    void styleHist(TH1* h, const char* xtitle, const char* ytitle = "Counts") {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.6);
        h->SetMarkerColor(kBlack);
        // Thin blue line connecting bin centres (drawn via the "L" suffix in
        // the Draw option below). LineColor also colours the error bars.
        h->SetLineColor(kBlue);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(xtitle);
        h->GetYaxis()->SetTitle(ytitle);
        h->GetYaxis()->SetMaxDigits(3);   // avoid 1e6 in axis labels
    }

    // Draw a TLatex annotation in pad-relative NDC coords (top-right by default).
    void drawLabel(const char* text, double x = 0.94, double y = 0.94, double size = 0.04) {
        TLatex* l = new TLatex();
        l->SetNDC(true);
        l->SetTextAlign(33);   // top-right
        l->SetTextSize(size);
        l->DrawLatex(x, y, text);
    }

    void saveAs(TCanvas* c, const std::string& base) {
        gSystem->mkdir("plots/output", kTRUE);
        c->SaveAs(("plots/output/" + base + ".pdf").c_str());
        c->SaveAs(("plots/output/" + base + ".png").c_str());
        std::cout << "Saved plots/output/" << base << ".{pdf,png}\n";
    }
}

void pippim_masses_cor(const char* file = "/home/przygoda/HADES/FAT/FAT/output_pippim.root") {

    applyGlobalStyle();

    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << file << std::endl;
        return;
    }

    // ------------------------------------------------------------------
    // 1. M(gamma gamma) — sharp pi0 peak expected at ~0.135 GeV/c^2
    // ------------------------------------------------------------------
    {
        auto h = (TH1*)f->Get("compound_cor/mass_gg_cor");
        if (!h) { std::cerr << "missing compound_cor/mass_gg_cor\n"; }
        else {
            styleHist(h,
                      "M_{#gamma#gamma} [GeV/c^{2}]",
                      "Counts / 5 MeV/c^{2}");
            auto* c = new TCanvas("c_m_gg_cor", "M(gg)", 800, 600);
            h->SetMinimum(0);
            h->GetXaxis()->SetRangeUser(0.0, 0.6);
            h->Draw("HIST L");          // thin blue polyline through bin centres
            h->Draw("PE SAME");         // overlay markers + error bars
            drawLabel("HADES p+p @ 4.5 GeV   (CORRECTED)", 0.94, 0.88, 0.034);
            drawLabel("ECAL: ecal_mult == 2, ecal_quality", 0.94, 0.83, 0.030);
            c->Update();
            saveAs(c, "m_gg_cor");
        }
    }

    // ------------------------------------------------------------------
    // 2. M(pi+ pi-) — wide spectrum, structure from rho/omega/etc.
    // ------------------------------------------------------------------
    {
        auto h = (TH1*)f->Get("compound_cor/mass_pippim_cor");
        if (!h) { std::cerr << "missing compound_cor/mass_pippim_cor\n"; }
        else {
            styleHist(h,
                      "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
                      "Counts / 10 MeV/c^{2}");
            auto* c = new TCanvas("c_m_pippim_cor", "M(pippim)", 800, 600);
            c->SetGridx();
            h->SetMinimum(0);
            h->GetXaxis()->SetRangeUser(0.2, 1.2);
            h->Draw("HIST L");          // thin blue polyline through bin centres
            h->Draw("PE SAME");         // overlay markers + error bars
            drawLabel("HADES p+p @ 4.5 GeV   (CORRECTED)", 0.94, 0.88, 0.034);
            drawLabel("Tree: PipPim_ID, isBest+vertex_z", 0.94, 0.83, 0.030);
            c->Update();
            saveAs(c, "m_pippim_cor");
        }
    }

    // ------------------------------------------------------------------
    // 3. M(pi+ pi- gamma gamma) — eta candidate (~0.547 GeV/c^2)
    //    after pi0 mass window on M(gg)
    // ------------------------------------------------------------------
    {
        auto h = (TH1*)f->Get("compound_cor/mass_pippimgg_cor");
        if (!h) { std::cerr << "missing compound_cor/mass_pippimgg_cor\n"; }
        else {
            styleHist(h,
                      "M_{#pi^{+}#pi^{-}#gamma#gamma} [GeV/c^{2}]",
                      "Counts / 6 MeV/c^{2}");
            auto* c = new TCanvas("c_m_pippimgg_cor", "M(pippimgg)", 800, 600);
            c->SetGridx();
            h->SetMinimum(0);
            h->GetXaxis()->SetRangeUser(0.4, 1.4);
            h->Draw("HIST L");          // thin blue polyline through bin centres
            h->Draw("PE SAME");         // overlay markers + error bars
            drawLabel("HADES p+p @ 4.5 GeV   (CORRECTED)", 0.94, 0.88, 0.034);
            drawLabel("ecal_mult==2 + #pi^{0} window on M(#gamma#gamma)",
                      0.94, 0.83, 0.030);
            c->Update();
            saveAs(c, "m_pippimgg_cor");
        }
    }

    f->Close();
    std::cout << "\nDone. Three plots saved under plots/output/\n";
}
