// rep_pimslot.C -- the pi- slot of the PipPimEpEm hypothesis before the PID cut, split by the true particle
// (figs/banana_pimslot.pdf). Input: banana_lep.root (banana_fill.C on the lepton PID sample, see rep_banana.C).
// Curves: expected beta_new of an electron / a muon measured with the pion mass, 1/b = (1/b_true + 1/b_pi)/2.
// Run: root -l -b -q 'docs/latex/macros/rep_pimslot.C+'
#include <TFile.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <cmath>

void rep_pimslot()
{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetPalette(kBird);
   gStyle->SetLineScalePS(1.6);
   TFile f("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/report/banana_lep.root");
   TFile fc("/hdd3/przygoda/HADES/files/sim/pp45_pid_cuts_sim_betanew_v6.root");
   TCutG* c = (TCutG*)fc.Get("pim_cut_3sig");
   TGraph* g = new TGraph(c->GetN());
   for (int i = 0; i < c->GetN(); ++i) g->SetPoint(i, -c->GetX()[i] / 1000., c->GetY()[i]);
   g->SetLineColor(kRed + 1);
   g->SetLineWidth(3);
   auto half = [](double m) {
      TF1* fn = new TF1(Form("half%g", m), [m](double* x, double*) {
         double p = -x[0], mpi = 0.13957;
         double bt = p / std::sqrt(p * p + m * m), bp = p / std::sqrt(p * p + mpi * mpi);
         return 1. / (0.5 * (1. / bt + 1. / bp)); }, -1.0, -0.06, 0);
      fn->SetNpx(400);
      return fn; };
   TF1* fe = half(0.000511);
   fe->SetLineColor(kBlack); fe->SetLineStyle(2); fe->SetLineWidth(3);
   TF1* fm = half(0.10566);
   fm->SetLineColor(kMagenta + 1); fm->SetLineStyle(7); fm->SetLineWidth(3);
   const char* sp[3] = {"pim", "em", "other"};
   const char* tt[3] = {"true #pi^{-}", "true e^{-}", "other (#mu^{-}, K^{-}, ...)"};
   TCanvas cv("cv", "", 1600, 1350);
   cv.Divide(2, 2, 0.002, 0.002);
   for (int k = 0; k < 3; ++k) {
      TVirtualPad* p = cv.cd(k + 1);
      p->SetLogz(); p->SetLeftMargin(0.14); p->SetRightMargin(0.14); p->SetBottomMargin(0.14); p->SetTopMargin(0.08);
      TH2D* s = (TH2D*)f.Get(Form("hSlotTrue_pim_%s", sp[k]));
      TH2D* h = new TH2D(Form("h%d", k), "", s->GetNbinsX(), s->GetXaxis()->GetXmin() / 1000., s->GetXaxis()->GetXmax() / 1000.,
                         s->GetNbinsY(), s->GetYaxis()->GetXmin(), s->GetYaxis()->GetXmax());
      for (int i = 1; i <= s->GetNbinsX(); ++i) for (int j = 1; j <= s->GetNbinsY(); ++j) h->SetBinContent(i, j, s->GetBinContent(i, j));
      h->GetXaxis()->SetRangeUser(-1.0, 0.0);
      h->GetYaxis()->SetRangeUser(0.5, 1.12);
      h->GetXaxis()->SetTitle("charge #times momentum  [GeV/c]");
      h->GetYaxis()->SetTitle("#beta_{new}");
      for (TAxis* a : {h->GetXaxis(), h->GetYaxis()}) { a->SetTitleSize(0.055); a->SetLabelSize(0.048); }
      h->GetYaxis()->SetTitleOffset(1.2);
      h->GetXaxis()->SetNdivisions(505);
      h->GetZaxis()->SetLabelSize(0.042);
      h->SetMinimum(1);
      h->Draw("COLZ");
      g->Draw("L same"); fe->Draw("same"); fm->Draw("same");
      TLatex t; t.SetNDC(); t.SetTextSize(0.06);
      t.DrawLatex(0.17, 0.935, Form("#pi^{-} slot: %s  (%.0f)", tt[k], s->GetEntries()));
   }
   cv.cd(4);
   TLegend* lg = new TLegend(0.08, 0.30, 0.95, 0.75);
   lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.060);
   lg->AddEntry(g, "#pi^{-} cut, simulation v6, 3#sigma", "l");
   lg->AddEntry(fe, "e^{-} measured as #pi", "l");
   lg->AddEntry(fm, "#mu^{-} measured as #pi", "l");
   lg->Draw();
   TLatex fx; fx.SetNDC(); fx.SetTextSize(0.058);
   fx.DrawLatex(0.10, 0.20, "1/#beta_{new} = (1/#beta_{true} + 1/#beta_{#pi}) / 2");
   const char* base = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs/banana_pimslot";
   cv.SaveAs(Form("%s.eps", base));
   gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", base, base, base));
   cv.SaveAs(Form("%s.png", base));
}
