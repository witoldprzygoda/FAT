// rep_banana.C -- beta_new versus q*p for the report (figs/banana_*.pdf, tables/banana_inside.tex)
//
// Inputs (filled on vae25 by banana_fill.C from simulation PID samples written WITHOUT PID cut;
// isBest==1 && eVertReco_z>-500, no trigger condition):
//   banana_lep.root : lepton PID sample (pp45_sim_leptons, LEP_merged, 101 files), tree PipPimEpEm,
//                     slots pip, ep (+p) and pim, em (-p)
//   banana_hh.root  : hadron PID sample (pp45_sim_pidsamples, HH_elossfix, first 150 files), tree PPipPim,
//                     slots p, pip (+p) and pim (-p)
// Cuts: simulation pp45_pid_cuts_sim_betanew_v6.root (as used by the simulation PAT: p, pi+-: 3 sigma;
// e+-: 5 sigma) solid; experiment pp45_pid_cuts_exp_betanew_v1.root (same levels) dashed.
//
// Figures:
//   banana_all.pdf     every track in a slot whose mass hypothesis is its own kind: HH slots p, pip, pim
//                      + LEP slots ep, em (the "banana plot" with all cuts)
//   banana_leptons.pdf lepton region, all four slots of PipPimEpEm
//   banana_ppee.pdf    all four slots of PipPimEpEm (what the PID cut acts on in this analysis)
//   banana_species.pdf the tracks of banana_ppee split by the GEANT identity
//
// Run: root -l -b -q 'docs/latex/macros/rep_banana.C+'
#include <TFile.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TString.h>
#include <cstdio>
#include <string>

namespace rb {
void savePdf(TCanvas& c, const char* base);
const char* kDir  = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/report/";
const char* kSim  = "/hdd3/przygoda/HADES/files/sim/pp45_pid_cuts_sim_betanew_v6.root";
const char* kOut  = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs/";
const char* kTab  = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables/";

// species: p, pi+, pi-, e+, e-
const char* kSp[5]    = {"p", "pip", "pim", "ep", "em"};
const char* kLab[5]   = {"p", "#pi^{+}", "#pi^{-}", "e^{+}", "e^{-}"};
const char* kTex[6]   = {"p", "$\\pi^+$", "$\\pi^-$", "$e^+$", "$e^-$", "other"};
const char* kLevel[5] = {"3sig", "3sig", "3sig", "5sig", "5sig"};
const int   kSign[5]  = {+1, +1, -1, +1, -1};
const int   kCol[5]   = {kOrange + 7, kBlue + 1, kCyan + 2, kRed + 1, kMagenta + 1};
const char* kTrue[6]  = {"p", "pip", "pim", "ep", "em", "other"};

TH2D* get(TFile& f, const char* n)
{
   TH2D* h = (TH2D*)f.Get(n);
   if (!h) { printf("missing %s in %s\n", n, f.GetName()); return nullptr; }
   return h;
}

// histogram in MeV/c -> GeV/c, optional rebin
TH2D* toGeV(TH2D* h, const char* name, int rx, int ry)
{
   TH2D* c = (TH2D*)h->Clone(Form("%s_tmp", name));
   c->Rebin2D(rx, ry);
   TH2D* g = new TH2D(name, "", c->GetNbinsX(), c->GetXaxis()->GetXmin() / 1000., c->GetXaxis()->GetXmax() / 1000.,
                      c->GetNbinsY(), c->GetYaxis()->GetXmin(), c->GetYaxis()->GetXmax());
   for (int i = 1; i <= c->GetNbinsX(); ++i)
      for (int j = 1; j <= c->GetNbinsY(); ++j) g->SetBinContent(i, j, c->GetBinContent(i, j));
   g->SetEntries(c->GetEntries());
   delete c;
   return g;
}

// polygon (p, beta) in MeV/c -> graph (q*p, beta) in GeV/c
TGraph* toGraph(TCutG* cut, int sign, int col, int style, int width)
{
   TGraph* g = new TGraph(cut->GetN());
   for (int i = 0; i < cut->GetN(); ++i) g->SetPoint(i, sign * cut->GetX()[i] / 1000., cut->GetY()[i]);
   g->SetLineColor(col);
   g->SetLineStyle(style);
   g->SetLineWidth(width);
   g->SetFillStyle(0);
   return g;
}

void axes(TH2D* h, double sz)
{
   h->GetXaxis()->SetTitle("charge #times momentum  [GeV/c]");
   h->GetYaxis()->SetTitle("#beta_{new}");
   h->GetXaxis()->SetTitleSize(sz);
   h->GetYaxis()->SetTitleSize(sz);
   h->GetXaxis()->SetLabelSize(sz * 0.85);
   h->GetYaxis()->SetLabelSize(sz * 0.85);
   h->GetZaxis()->SetLabelSize(sz * 0.75);
   h->GetXaxis()->SetTitleOffset(1.0);
   h->GetYaxis()->SetTitleOffset(0.75);
   h->SetStats(0);
}

// fraction of the entries of h (MeV/c) inside the polygon (bin centres, 10 MeV/c x 0.0025)
double fracInside(TH2D* h, TCutG* cut, int sign)
{
   double in = 0, all = 0;
   for (int i = 1; i <= h->GetNbinsX(); ++i)
      for (int j = 1; j <= h->GetNbinsY(); ++j) {
         double w = h->GetBinContent(i, j);
         if (w <= 0) continue;
         double x = h->GetXaxis()->GetBinCenter(i), y = h->GetYaxis()->GetBinCenter(j);
         all += w;
         if (x * sign > 0 && cut->IsInside(sign * x, y)) in += w;
      }
   return all > 0 ? in / all : 0;
}

void drawCuts(TCutG** cS, TCutG** cE)
{
   for (int k = 0; k < 5; ++k) {
      if (cE) toGraph(cE[k], kSign[k], kCol[k], 2, 2)->Draw("L same");
      toGraph(cS[k], kSign[k], kCol[k], 1, 3)->Draw("L same");
   }
}

TLegend* cutLegend(double x0, double y0, double x1, double y1)
{
   TLegend* lg = new TLegend(x0, y0, x1, y1);
   lg->SetBorderSize(0);
   lg->SetFillStyle(1001);
   lg->SetFillColor(kWhite);
   lg->SetTextSize(0.030);
   TGraph* gs = new TGraph(1);
   TGraph* ge = new TGraph(1);
   gs->SetLineWidth(3); gs->SetLineColor(kBlack);
   ge->SetLineWidth(2); ge->SetLineStyle(2); ge->SetLineColor(kBlack);
   lg->AddEntry(gs, "simulation cuts (v6)", "l");
   lg->AddEntry(ge, "experiment cuts (v1)", "l");
   lg->AddEntry((TObject*)nullptr, "p, #pi^{#pm}: 3#sigma;   e^{#pm}: 5#sigma", "");
   return lg;
}

void labels(double sz)
{
   TLatex t;
   t.SetTextSize(sz);
   t.SetTextFont(62);
   t.SetTextColor(kCol[0]); t.DrawLatex(1.55, 0.60, kLab[0]);
   t.SetTextColor(kCol[1]); t.DrawLatex(0.50, 0.86, kLab[1]);
   t.SetTextColor(kCol[2]); t.DrawLatex(-0.80, 0.86, kLab[2]);
   t.SetTextColor(kCol[3]); t.DrawLatex(1.95, 1.08, kLab[3]);
   t.SetTextColor(kCol[4]); t.DrawLatex(-1.95, 1.08, kLab[4]);
}

void bigPlot(TH2D* h, const char* file, const char* head, double xlo, double xhi, double ylo, double yhi,
             TCutG** cS, TCutG** cE, bool withLabels, bool withLegend)
{
   gStyle->SetPalette(kGreyScale);
   TColor::InvertPalette();
   gStyle->SetNumberContours(64);
   TCanvas c(Form("c_%s", file), "", 1800, 1150);
   c.SetLeftMargin(0.085);
   c.SetRightMargin(0.105);
   c.SetBottomMargin(0.11);
   c.SetTopMargin(0.03);
   c.SetLogz();
   axes(h, 0.05);
   h->GetXaxis()->SetRangeUser(xlo, xhi);
   h->GetYaxis()->SetRangeUser(ylo, yhi);
   h->SetMinimum(1);
   h->Draw("COLZ");
   drawCuts(cS, cE);
   if (xlo < 0 && xhi > 0) {
      TLine* l0 = new TLine(0, ylo, 0, yhi);
      l0->SetLineStyle(3);
      l0->Draw();
   }
   if (withLabels) labels(0.055);
   if (withLegend) cutLegend(0.60, 0.15, 0.885, 0.33)->Draw();
   (void)head; // the description is in the figure caption of the report
   savePdf(c, Form("%s%s", kOut, file));
   c.SaveAs(Form("%s%s.png", kOut, file));
   TColor::InvertPalette();
}

// vector PDF with the exact canvas size: EPS (tight bounding box) -> ps2pdf -dEPSCrop
void savePdf(TCanvas& c, const char* base)
{
   TString eps = Form("%s.eps", base), pdf = Form("%s.pdf", base);
   c.SaveAs(eps);
   gSystem->Exec(Form("ps2pdf -dEPSCrop %s %s && rm -f %s", eps.Data(), pdf.Data(), eps.Data()));
}
} // namespace rb

void rep_banana()
{
   using namespace rb;
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetTextFont(42);
   gStyle->SetLabelFont(42, "XYZ");
   gStyle->SetTitleFont(42, "XYZ");
   gStyle->SetLineScalePS(1.6);

   TFile fl(Form("%sbanana_lep.root", kDir)), fhh(Form("%sbanana_hh.root", kDir));
   TFile fs(kSim), fe(Form("%spp45_pid_cuts_exp_betanew_v1.root", kDir));
   TCutG *cS[5], *cE[5];
   for (int k = 0; k < 5; ++k) {
      cS[k] = (TCutG*)fs.Get(Form("%s_cut_%s", kSp[k], kLevel[k]));
      cE[k] = (TCutG*)fe.Get(Form("%s_cut_%s", kSp[k], kLevel[k]));
   }

   // ---- figure 1: every track in its own mass hypothesis
   TH2D* hComb = (TH2D*)get(fhh, "hSlot_p")->Clone("hComb");
   hComb->Add(get(fhh, "hSlot_pip"));
   hComb->Add(get(fhh, "hSlot_pim"));
   hComb->Add(get(fl, "hSlot_ep"));
   hComb->Add(get(fl, "hSlot_em"));
   bigPlot(toGeV(hComb, "hCombG", 2, 2), "banana_all",
           "SMASH p+p 4.5 GeV, simulation without PID cut: p, #pi^{#pm} from p#pi^{+}#pi^{-}, e^{#pm} from #pi^{+}#pi^{-}e^{+}e^{-} (isBest = 1)",
           -2.0, 3.0, 0.0, 1.2, cS, cE, true, true);

   // ---- figure 2: lepton region, PipPimEpEm all slots
   TH2D* hAll = get(fl, "hAll");
   {
      TH2D* h = toGeV(hAll, "hAllZ", 1, 1);
      bigPlot(h, "banana_leptons",
              "Hypothesis #pi^{+}#pi^{-}e^{+}e^{-} without PID cut, all four slots: lepton region",
              -1.5, 1.5, 0.80, 1.12, cS, cE, false, false);
   }
   // ---- figure 3: PipPimEpEm all slots
   bigPlot(toGeV(hAll, "hAllG", 2, 2), "banana_ppee",
           "Hypothesis #pi^{+}#pi^{-}e^{+}e^{-} without PID cut, all four slots (isBest = 1)",
           -2.0, 3.0, 0.0, 1.2, cS, cE, true, true);

   // ---- figure 4: true species in PipPimEpEm
   {
      gStyle->SetPalette(kBird);
      TCanvas c("c_species", "", 1800, 1150);
      c.Divide(3, 2, 0.002, 0.002);
      const char* ttl[6] = {"true p", "true #pi^{+}", "true #pi^{-}", "true e^{+}", "true e^{-}", "other (#mu, K, d, ...)"};
      for (int k = 0; k < 6; ++k) {
         TVirtualPad* p = c.cd(k + 1);
         p->SetLogz();
         p->SetLeftMargin(0.13);
         p->SetRightMargin(0.13);
         p->SetBottomMargin(0.13);
         p->SetTopMargin(0.08);
         TH2D* src = get(fl, Form("hTrue_%s", kTrue[k]));
         TH2D* h = toGeV(src, Form("hSpG%d", k), 2, 2);
         axes(h, 0.055);
         h->GetYaxis()->SetTitleOffset(1.05);
         h->GetXaxis()->SetRangeUser(-2.0, 3.0);
         h->GetYaxis()->SetRangeUser(0.0, 1.2);
         h->SetMinimum(1);
         h->Draw("COLZ");
         for (int m = 0; m < 5; ++m) toGraph(cS[m], kSign[m], (m == k ? kCol[m] : kGray + 1), 1, (m == k ? 3 : 1))->Draw("L same");
         TLatex s;
         s.SetNDC();
         s.SetTextSize(0.065);
         s.DrawLatex(0.15, 0.935, Form("%s   (%.0f tracks)", ttl[k], src->GetEntries()));
      }
      savePdf(c, Form("%sbanana_species", kOut));
      c.SaveAs(Form("%sbanana_species.png", kOut));
   }

   // ---- table: composition of each slot and fraction inside the cut of the slot
   struct Slot { TFile* f; const char* slot; int cut; const char* tex; };
   Slot sl[7] = {{&fhh, "p", 0, "p$\\pi^+\\pi^-$: p"},     {&fhh, "pip", 1, "p$\\pi^+\\pi^-$: $\\pi^+$"},
                 {&fhh, "pim", 2, "p$\\pi^+\\pi^-$: $\\pi^-$"},
                 {&fl, "pip", 1, "$\\pi^+\\pi^-e^+e^-$: $\\pi^+$"}, {&fl, "ep", 3, "$\\pi^+\\pi^-e^+e^-$: $e^+$"},
                 {&fl, "pim", 2, "$\\pi^+\\pi^-e^+e^-$: $\\pi^-$"}, {&fl, "em", 4, "$\\pi^+\\pi^-e^+e^-$: $e^-$"}};
   FILE* f = fopen(Form("%sbanana_inside.tex", kTab), "w");
   FILE* g = fopen(Form("%sbanana_inside.txt", kTab), "w");
   fprintf(g, "# slot: true species -> tracks, fraction inside the simulation cut v6 of the slot (bin centres 10 MeV/c x 0.0025)\n");
   fprintf(f, "%% generated by macros/rep_banana.C\n");
   fprintf(f, "\\begin{table}[H]\n\\centering\n");
   fprintf(f, "\\caption{Composition of the slots before the PID cut and the fraction of each kind of track that lies "
              "inside the simulation cut of its slot (v6; p, $\\pi^\\pm$ $3\\sigma$, $e^\\pm$ $5\\sigma$). "
              "Each cell: share of the slot's tracks / fraction inside the cut. Simulation PID samples, "
              "$\\mathrm{isBest}=1$, vertex $z>-500$\\,mm, no trigger condition. Shares below 0.1\\,\\%% are omitted.}\n");
   fprintf(f, "\\label{tab:inside}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}l r l l l l l l@{}}\n\\toprule\n");
   fprintf(f, "Hypothesis: slot & tracks");
   for (int k = 0; k < 6; ++k) fprintf(f, " & %s", kTex[k]);
   fprintf(f, " \\\\\n\\midrule\n");
   for (int i = 0; i < 7; ++i) {
      double tot = get(*sl[i].f, Form("hSlot_%s", sl[i].slot))->GetEntries();
      fprintf(f, "%s & %.0f", sl[i].tex, tot);
      fprintf(g, "%s %-4s %10.0f", sl[i].f == &fhh ? "HH " : "LEP", sl[i].slot, tot);
      for (int k = 0; k < 6; ++k) {
         TH2D* h = get(*sl[i].f, Form("hSlotTrue_%s_%s", sl[i].slot, kTrue[k]));
         double n = h->GetEntries(), share = tot > 0 ? n / tot : 0;
         double in = fracInside(h, cS[sl[i].cut], kSign[sl[i].cut]);
         fprintf(g, "  %s %9.0f (%.4f) in %.4f", kTrue[k], n, share, in);
         if (share < 0.001) fprintf(f, " & --");
         else fprintf(f, " & %.3f\\,/\\,%.3f", share, in);
      }
      fprintf(f, " \\\\\n");
      fprintf(g, "\n");
      if (i == 2) fprintf(f, "\\midrule\n");
   }
   fprintf(f, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
   fclose(f);
   fclose(g);
   printf("written %sbanana_*.pdf/png, %sbanana_inside.{tex,txt}\n", kOut, kTab);
}
