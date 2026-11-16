// ============================================================================
// fix_wp0defs.C -- review fix (prefix fix): Fig. 1 of wp0_check.C (WP0, "definitions on the simulated bump") redrawn with the
// title "C2 + C3" written without TLatex subscripts (the subscripts produced a stray space). Block copied unchanged otherwise.
// Reads hists_v1.root only. Output: docs/latex/figs_all/fix_definitions.{pdf,png}
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_fix", kTRUE)' 'plots/all_report/fix_wp0defs.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "defs.h"
#include "hists_io.h"
#include "s5_blind_ana_fit.h"
void fix_wp0defs() {
    using namespace alldefs;
    bla::Style();
    gStyle->SetPalette(kGreyScale);
    TFile* f = TFile::Open(kHistFile, "READ");
    if (!f || f->IsZombie()) { printf("fix_wp0defs: cannot open the factory file\n"); return; }
    // ---------------------------------------------------------------- figure 1: definitions on the simulated bump
    {
        TCanvas* c = new TCanvas("c1", "", 1000, 800);
        c->SetRightMargin(0.15);
        c->SetLeftMargin(0.13);
        TH2* m2 = (TH2*)allio::Sim(f, "S0", "nom", kC2, "B_map");
        TH2* m3 = (TH2*)allio::Sim(f, "S0", "nom", kC3, "B_map");
        m2->Add(m3);
        m2->Scale(1. / 0.01 / 0.01 * 1e-3);   // nb per bin -> ub / (GeV/c2)^2
        m2->GetXaxis()->SetRangeUser(0.40, 0.70);
        m2->GetYaxis()->SetRangeUser(0.27, 0.55);
        m2->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
        m2->GetYaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]");
        m2->GetZaxis()->SetTitle("d^{2}#sigma/dM_{4}dM_{#pi#pi} [#mub/(GeV/c^{2})^{2}]");
        m2->GetZaxis()->SetTitleOffset(1.3);
        m2->Draw("COLZ");
        // windows B and P (vertical bands), CR-U, CR-K, CR-L boxes
        auto box = [&](double x1, double y1, double x2, double y2, int col, int style) {
            TBox* b = new TBox(x1, y1, x2, y2);
            b->SetFillStyle(0); b->SetLineColor(col); b->SetLineWidth(3); b->SetLineStyle(style); b->Draw();
        };
        box(kWin[kWB].lo, 0.27, kWin[kWB].hi, 0.55, kBlue + 1, 1);
        box(kWin[kWP].lo, 0.27, kWin[kWP].hi, 0.55, kRed + 1, 1);
        box(kWin[kWCRU].lo, 0.27, kWin[kWCRU].hi, 0.55, kGreen + 2, 1);
        box(kWin[kWCRK].lo, 0.43, kWin[kWCRK].hi, 0.55, kMagenta + 1, 1);
        box(kWin[kWCRL].lo, 0.27, kWin[kWCRL].hi, 0.55, kOrange + 7, 2);
        for (int s = 0; s < kNSlice - 1; ++s) {
            TLine* l = new TLine(0.40, kSliceEdge[s], 0.70, kSliceEdge[s]);
            l->SetLineStyle(3); l->SetLineWidth(2); l->Draw();
        }
        // edge E(m, q) for q = 0 and 60 MeV
        for (int k = 0; k < 2; ++k) {
            TGraph* g = new TGraph();
            for (int i = 0; i <= 200; ++i) {
                const double m = 0.27 + (kEdgeLimit - 0.27) * i / 200.;
                const double e = Edge(m, k == 0 ? 0. : 0.060);
                if (e > 0) g->SetPoint(g->GetN(), e, m);
            }
            g->SetLineColor(k == 0 ? kBlack : kGray + 2); g->SetLineWidth(3); g->SetLineStyle(k == 0 ? 1 : 7);
            g->Draw("L");
        }
        TLatex tx; tx.SetTextSize(0.033);
        tx.DrawLatex(0.480, 0.535, "#color[602]{B}");
        tx.DrawLatex(0.537, 0.535, "#color[634]{P}");
        tx.DrawLatex(0.60, 0.535, "#color[418]{CR-U}");
        tx.DrawLatex(0.44, 0.45, "#color[617]{CR-K (S5)}");
        for (int s = 0; s < kNSlice; ++s) {
            const double y = s == 0 ? 0.30 : s == kNSlice - 1 ? 0.505 : 0.5 * (kSliceEdge[s - 1] + kSliceEdge[s]);
            tx.DrawLatex(0.680, y - 0.004, Form("S%d", s + 1));
        }
        tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.13, 0.95, "simulation: bump C2 + C3, S0, w #times V3;  M_{4}' = 0.99615 M_{4},  M_{#pi#pi}' = 0.9952 M_{#pi#pi}");
        TLegend* lg = new TLegend(0.54, 0.66, 0.80, 0.78);
        lg->SetFillStyle(1001);
        TGraph* g0 = new TGraph(); g0->SetLineWidth(3);
        TGraph* g1 = new TGraph(); g1->SetLineWidth(3); g1->SetLineStyle(7); g1->SetLineColor(kGray + 2);
        lg->AddEntry(g0, "edge E(m, 0)", "l");
        lg->AddEntry(g1, "edge E(m, 60 MeV/c^{2})", "l");
        lg->SetTextSize(0.032);
        lg->Draw();
        bla::SaveFig(c, "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/fix_definitions");
        delete c;
    }
    f->Close();
}
