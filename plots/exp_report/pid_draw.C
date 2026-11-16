// pid_draw.C -- PID chapter of the experimental-data report: figures and tables.
//
// What it does
//   Draws beta_new versus charge x momentum for the data PID samples and for the final analysed sample,
//   with the PID polygons of the data production (betanew_v1; p, pi+-: 3 sigma solid, e+-: 5 sigma solid,
//   1 sigma dashed), the overlay of the data (v1) and simulation (v6) polygons, the START-based DST beta
//   of the leptons of the final sample, and writes the tables of the chapter.
//
// Inputs (all read-only)
//   IN = /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/inputs/
//     banana_exp_hh_PipPim_ID.root      June 2026 hadron production, tree PipPim_ID (both tracks with the pion
//                                       mass), days 034, 045, 049, 062, first 40 M rows per file, ALL rows
//                                       (isBest -1/0/1 = no PID cut), vz > -500 mm, no trigger condition
//     banana_exp_hh_PP_ID.root          the same files, tree PP_ID (both tracks with the proton mass)
//     banana_exp_june_EpEm_ID.root      June lepton production (no PID cut), tree EpEm_ID, days 036/045/050/062
//     banana_exp_june_PipPimEpEm_ID.root June lepton production (no PID cut), tree PipPimEpEm_ID, 29 days
//     banana_exp_EpEm_ID_betanew.root   betanew production (PID applied, e+- 5 sigma), EpEm_ID, 4 days
//     pp45_pid_cuts_exp_betanew_v1.root data polygons (md5 b8f0cc2f00c3c0e6d50627ea0bfc8c22)
//     pp45_pid_cuts_sim_betanew_v6.root simulation polygons
//     (the banana files were filled with banana_fill_exp.C, copy in IN; x = sign x <slot>_p, uncorrected
//      momentum, y = <slot>_beta_new; 10 MeV/c x 0.0025 bins)
//   FINAL = /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/pid_final.root
//     (pid_final_fill.C on the 29 betanew day files, merged by pid_run_final.sh)
// Outputs
//   /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/pid_*.pdf and .png:
//     pid_banana_hadrons, pid_banana_leptons, pid_banana_leptons_zoom, pid_banana_all, pid_banana_final,
//     pid_banana_final_slots, pid_polygons_v1_v6, pid_dstbeta
//   /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/pid_samples.{txt,tex}, pid_polygons.{txt,tex},
//     pid_dstbeta.{txt,tex}
// Run:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid", kTRUE)' \
//        '/home/przygoda/HADES/FAT/FAT/plots/exp_report/pid_draw.C+'
// Changes by the verifier (2026-11-14, figures only, tables unchanged): pid_dstbeta tof panel y maximum x3000 (the legend
//   covered the peak), pid_polygons_v1_v6 last panel y range 0-0.030 and a shorter axis title (title overlapped the top label),
//   wider left margin for the e+- panels.
//   fracInside (Table pid_samples, banana samples): denominator now includes under/overflow (it was the in-range sum).
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TF1.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

namespace pd {
const char* kIn    = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/inputs/";
const char* kFinal = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/pid_final.root";
const char* kFig   = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/";
const char* kTab   = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/";

// species: p, pi+, pi-, e+, e-
const char* kSp[5]   = {"p", "pip", "pim", "ep", "em"};
const char* kLab[5]  = {"p", "#pi^{+}", "#pi^{-}", "e^{+}", "e^{-}"};
const char* kPatLv[5] = {"3sig", "3sig", "3sig", "5sig", "5sig"};
const int   kSign[5] = {+1, +1, -1, +1, -1};
const int   kCol[5]  = {kOrange + 7, kBlue + 1, kCyan + 2, kRed + 1, kMagenta + 1};
const char* kLev[6]  = {"1sig", "2sig", "25sig", "3sig", "35sig", "5sig"};

TCutG* gV1[5][6];
TCutG* gV6[5][6];

// TLatex in PostScript widens a plain space in front of a Symbol-font character (pi, times, beta, ...):
// replace " #pi" by "#kern[0.8]{#pi}" etc. (checked in the PDF output)
TString L(const char* s)
{
   TString r(s);
   const char* sym[] = {"#pi", "#times", "#beta", "#geq", "#leq", "#sigma", "#approx"};
   for (const char* x : sym) r.ReplaceAll(TString(" ") + x, Form("#kern[0.8]{%s}", x));
   return r;
}

void savePdf(TCanvas& c, const char* base)
{
   TString eps = Form("%s%s.eps", kFig, base), pdf = Form("%s%s.pdf", kFig, base);
   c.SaveAs(eps);
   gSystem->Exec(Form("ps2pdf -dEPSCrop %s %s && rm -f %s", eps.Data(), pdf.Data(), eps.Data()));
   c.SaveAs(Form("%s%s.png", kFig, base));
}

TH2D* get2(TFile& f, const char* n)
{
   TH2D* h = (TH2D*)f.Get(n);
   if (!h) { printf("missing %s in %s\n", n, f.GetName()); return nullptr; }
   return h;
}

// histogram in MeV/c -> GeV/c, optional rebin
TH2D* toGeV(const TH2D* h, const char* name, int rx, int ry)
{
   TH2D* c = (TH2D*)h->Clone(Form("%s_tmp", name));
   c->Rebin2D(rx, ry);
   TH2D* g = new TH2D(name, "", c->GetNbinsX(), c->GetXaxis()->GetXmin() / 1000., c->GetXaxis()->GetXmax() / 1000.,
                      c->GetNbinsY(), c->GetYaxis()->GetXmin(), c->GetYaxis()->GetXmax());
   for (int i = 1; i <= c->GetNbinsX(); ++i)
      for (int j = 1; j <= c->GetNbinsY(); ++j) g->SetBinContent(i, j, c->GetBinContent(i, j));
   g->SetEntries(c->GetEntries());
   g->SetDirectory(nullptr);
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

// all data polygons: PAT level solid, 1 sigma dashed; hl >= 0: only that species in colour, others grey
void drawCuts(int hl = -1)
{
   for (int k = 0; k < 5; ++k) {
      const int lv = (k < 3) ? 3 : 5;
      const int col = (hl < 0 || hl == k) ? kCol[k] : kGray + 1;
      const int w = (hl < 0 || hl == k) ? 3 : 2;
      toGraph(gV1[k][0], kSign[k], col, 2, w - 1)->Draw("L same");
      toGraph(gV1[k][lv], kSign[k], col, 1, w)->Draw("L same");
   }
}

void axes2(TH2D* h, double sz, const char* ytitle = "#beta_{new}")
{
   h->GetXaxis()->SetTitle(L("charge #times momentum  [GeV/c]"));
   h->GetYaxis()->SetTitle(L(ytitle));
   h->GetXaxis()->SetTitleSize(sz);
   h->GetYaxis()->SetTitleSize(sz);
   h->GetXaxis()->SetLabelSize(sz * 0.85);
   h->GetYaxis()->SetLabelSize(sz * 0.85);
   h->GetZaxis()->SetLabelSize(sz * 0.75);
   h->GetXaxis()->SetTitleOffset(1.0);
   h->GetYaxis()->SetTitleOffset(1.0);
   h->SetStats(0);
}

void padStyle(TVirtualPad* p)
{
   p->SetLeftMargin(0.11);
   p->SetRightMargin(0.12);
   p->SetBottomMargin(0.13);
   p->SetTopMargin(0.09);
   p->SetLogz();
}

void header(const char* txt, double sz = 0.048)
{
   TLatex t;
   t.SetNDC();
   t.SetTextSize(sz);
   t.SetTextFont(42);
   t.DrawLatex(0.11, 0.935, L(txt));
}

void speciesLabels(double sz)
{
   TLatex t;
   t.SetTextSize(sz);
   t.SetTextFont(62);
   t.SetTextColor(kCol[0]); t.DrawLatex(1.55, 0.58, kLab[0]);
   t.SetTextColor(kCol[1]); t.DrawLatex(0.45, 0.86, kLab[1]);
   t.SetTextColor(kCol[2]); t.DrawLatex(-0.75, 0.86, kLab[2]);
   t.SetTextColor(kCol[3]); t.DrawLatex(1.95, 1.09, kLab[3]);
   t.SetTextColor(kCol[4]); t.DrawLatex(-1.95, 1.09, kLab[4]);
}

void lineLegend(double x0, double y0, double x1, double y1, double sz)
{
   TLegend* lg = new TLegend(x0, y0, x1, y1);
   lg->SetBorderSize(0);
   lg->SetFillStyle(0);
   lg->SetTextSize(sz);
   TGraph* gs = new TGraph(1); gs->SetLineWidth(3); gs->SetLineColor(kBlack);
   TGraph* gd = new TGraph(1); gd->SetLineWidth(2); gd->SetLineStyle(2); gd->SetLineColor(kBlack);
   lg->AddEntry(gs, "data PID cuts (v1)", "l");
   lg->AddEntry((TObject*)nullptr, L("p, #pi^{#pm}: 3#sigma; e^{#pm}: 5#sigma"), "");
   lg->AddEntry(gd, "1#sigma", "l");
   lg->Draw();
}

// one panel: histogram (GeV/c) with polygons
void panel(TVirtualPad* p, TH2D* h, double xlo, double xhi, double ylo, double yhi, double zmax, const char* head,
           int hl = -1, bool labels = false, const char* ytitle = "#beta_{new}", bool cuts = true)
{
   p->cd();
   padStyle(p);
   axes2(h, 0.050, ytitle);
   h->GetXaxis()->SetRangeUser(xlo, xhi);
   h->GetYaxis()->SetRangeUser(ylo, yhi);
   h->SetMinimum(0.9);
   if (zmax > 0) h->SetMaximum(zmax);
   h->Draw("COLZ");
   if (cuts) drawCuts(hl);
   if (xlo < 0 && xhi > 0) {
      TLine* l0 = new TLine(0, ylo, 0, yhi);
      l0->SetLineStyle(3);
      l0->Draw();
   }
   if (labels) speciesLabels(0.050);
   header(head);
   p->RedrawAxis();
}

double maxIn(TH2D* h, double xlo, double xhi, double ylo, double yhi)
{
   double m = 0;
   for (int i = h->GetXaxis()->FindBin(xlo); i <= h->GetXaxis()->FindBin(xhi - 1e-9); ++i)
      for (int j = h->GetYaxis()->FindBin(ylo); j <= h->GetYaxis()->FindBin(yhi - 1e-9); ++j) m = std::max(m, h->GetBinContent(i, j));
   return m;
}

// fraction of the entries of h (MeV/c) inside the polygon (bin centres, 10 MeV/c x 0.0025)
// Denominator = ALL entries, including under- and overflow (beta_new outside 0-1.3, q p outside -2.5..3.5 GeV/c): such
// tracks are outside every polygon except protons above 3.5 GeV/c (0.55 % of the PP_ID slots), which are counted as outside.
// (Verifier fix 2026-11-14: the denominator used to be the in-range sum only, which inflated the fractions of the
//  reference slots, e.g. e- of PipPimEpEm_ID 63.4 % instead of 54.7 % and pi- of PipPim_ID 57.4 % instead of 47.6 %.)
double fracInside(const TH2D* h, TCutG* cut, int sign)
{
   double in = 0;
   const double all = h->Integral(0, h->GetNbinsX() + 1, 0, h->GetNbinsY() + 1);
   for (int i = 1; i <= h->GetNbinsX(); ++i) {
      const double x = h->GetXaxis()->GetBinCenter(i);
      for (int j = 1; j <= h->GetNbinsY(); ++j) {
         const double w = h->GetBinContent(i, j);
         if (w <= 0) continue;
         const double y = h->GetYaxis()->GetBinCenter(j);
         if (x * sign > 0 && cut->IsInside(sign * x, y)) in += w;
      }
   }
   return all > 0 ? in / all : 0;
}

// lower and upper edge of a polygon at momentum x (MeV/c): crossings of the vertical line
bool yRange(TCutG* c, double x, double& lo, double& hi)
{
   lo = 1e9; hi = -1e9;
   const int n = c->GetN();
   for (int i = 0; i < n; ++i) {
      const int j = (i + 1) % n;
      const double x1 = c->GetX()[i], x2 = c->GetX()[j], y1 = c->GetY()[i], y2 = c->GetY()[j];
      if ((x1 - x) * (x2 - x) > 0 || x1 == x2) continue;
      const double y = y1 + (y2 - y1) * (x - x1) / (x2 - x1);
      lo = std::min(lo, y); hi = std::max(hi, y);
   }
   return hi > lo;
}

TString pct(double f, int prec = 1) { return Form("%.*f", prec, 100. * f); }
} // namespace pd

void pid_draw()
{
   using namespace pd;
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetTextFont(42);
   gStyle->SetLabelFont(42, "XYZ");
   gStyle->SetTitleFont(42, "XYZ");
   gStyle->SetLineScalePS(1.6);
   gStyle->SetPalette(kGreyScale);
   TColor::InvertPalette();
   gStyle->SetNumberContours(64);

   TFile f1(Form("%spp45_pid_cuts_exp_betanew_v1.root", kIn)), f6(Form("%spp45_pid_cuts_sim_betanew_v6.root", kIn));
   for (int k = 0; k < 5; ++k)
      for (int l = 0; l < 6; ++l) {
         gV1[k][l] = (TCutG*)f1.Get(Form("%s_cut_%s", kSp[k], kLev[l]));
         gV6[k][l] = (TCutG*)f6.Get(Form("%s_cut_%s", kSp[k], kLev[l]));
         if (!gV1[k][l] || !gV6[k][l]) { printf("missing polygon %s %s\n", kSp[k], kLev[l]); return; }
      }
   TFile fHH(Form("%sbanana_exp_hh_PipPim_ID.root", kIn)), fPP(Form("%sbanana_exp_hh_PP_ID.root", kIn));
   TFile fEE(Form("%sbanana_exp_june_EpEm_ID.root", kIn)), fPE(Form("%sbanana_exp_june_PipPimEpEm_ID.root", kIn));
   TFile fBN(Form("%sbanana_exp_EpEm_ID_betanew.root", kIn));
   TFile fF(kFinal);

   const double X0 = -2.0, X1 = 3.0, Y0 = 0.0, Y1 = 1.2;

   // ===== (a) hadrons without ring: PipPim_ID (pion mass) and PP_ID (proton mass)
   TH2D* hHH = toGeV(get2(fHH, "hAll_noring"), "hHH", 2, 2);
   TH2D* hPP = toGeV(get2(fPP, "hAll_noring"), "hPP", 2, 2);
   {
      const double zm = std::max(maxIn(hHH, X0, X1, Y0, Y1), maxIn(hPP, X0, X1, Y0, Y1));
      TCanvas c("c_had", "", 2200, 950);
      c.Divide(2, 1, 0.002, 0.002);
      panel(c.cd(1), hHH, X0, X1, Y0, Y1, zm, "hypothesis #pi^{+}#pi^{-}: both tracks with the pion mass", -1, true);
      lineLegend(0.115, 0.14, 0.42, 0.255, 0.030);
      panel(c.cd(2), hPP, X0, X1, Y0, Y1, zm, "hypothesis pp: both tracks with the proton mass", -1, true);
      savePdf(c, "pid_banana_hadrons");
   }
   // ===== (b) lepton candidates with a ring: June EpEm_ID and the lepton slots of June PipPimEpEm_ID
   TH2D* hEEm = (TH2D*)get2(fEE, "hAll")->Clone("hEEm");
   TH2D* hPEm = (TH2D*)get2(fPE, "hSlot_ep")->Clone("hPEm");
   hPEm->Add(get2(fPE, "hSlot_em"));
   {
      TH2D* a = toGeV(hEEm, "hEE2", 2, 2);
      TH2D* b = toGeV(hPEm, "hPE2", 2, 2);
      const double zm = std::max(maxIn(a, X0, X1, Y0, Y1), maxIn(b, X0, X1, Y0, Y1));
      TCanvas c("c_lep", "", 2200, 950);
      c.Divide(2, 1, 0.002, 0.002);
      panel(c.cd(1), a, X0, X1, Y0, Y1, zm, "hypothesis e^{+}e^{-}: slots e^{+}, e^{-} (RICH ring), no PID cut", -1, true);
      lineLegend(0.115, 0.14, 0.42, 0.255, 0.030);
      panel(c.cd(2), b, X0, X1, Y0, Y1, zm, "hypothesis #pi^{+}#pi^{-}e^{+}e^{-}: slots e^{+}, e^{-}, no PID cut", -1, true);
      savePdf(c, "pid_banana_leptons");
      // zoom on the lepton band
      TH2D* az = toGeV(hEEm, "hEEz", 1, 1);
      TH2D* bz = toGeV(hPEm, "hPEz", 1, 1);
      const double zx0 = -2.0, zx1 = 2.0, zy0 = 0.80, zy1 = 1.12;
      const double zz = std::max(maxIn(az, zx0, zx1, zy0, zy1), maxIn(bz, zx0, zx1, zy0, zy1));
      TCanvas cz("c_lepz", "", 2200, 950);
      cz.Divide(2, 1, 0.002, 0.002);
      panel(cz.cd(1), az, zx0, zx1, zy0, zy1, zz, "e^{+}e^{-} hypothesis, lepton band (no PID cut)");
      panel(cz.cd(2), bz, zx0, zx1, zy0, zy1, zz, "#pi^{+}#pi^{-}e^{+}e^{-} hypothesis, slots e^{#pm} (no PID cut)");
      savePdf(cz, "pid_banana_leptons_zoom");
   }
   // ===== (c) all species, each in a slot of its own mass hypothesis
   {
      TH2D* s = (TH2D*)get2(fPP, "hAll_noring")->Clone("hSum");
      s->Add(get2(fHH, "hAll_noring"));
      s->Add(hEEm);
      TH2D* g = toGeV(s, "hSumG", 2, 2);
      TCanvas c("c_all", "", 1800, 1150);
      TVirtualPad* p = c.cd();
      panel(p, g, X0, X1, Y0, Y1, -1, "no PID cut: p from pp, #pi^{#pm} from #pi^{+}#pi^{-} (no ring), e^{#pm} from e^{+}e^{-} (ring)", -1, true);
      lineLegend(0.115, 0.14, 0.40, 0.25, 0.028);
      savePdf(c, "pid_banana_all");
   }
   // ===== (d) final analysed sample: four slots, beta_new and the DST (START) beta of the same tracks
   const char* slot4[4] = {"pip", "pim", "ep", "em"};
   {
      TH2D* bn = (TH2D*)get2(fF, "hBN_fat_pip")->Clone("hBNsum");
      TH2D* bd = (TH2D*)get2(fF, "hBD_fat_pip")->Clone("hBDsum");
      for (int k = 1; k < 4; ++k) { bn->Add(get2(fF, Form("hBN_fat_%s", slot4[k]))); bd->Add(get2(fF, Form("hBD_fat_%s", slot4[k]))); }
      TH2D* a = toGeV(bn, "hFinBN", 2, 2);
      TH2D* b = toGeV(bd, "hFinBD", 2, 2);
      const double zm = std::max(maxIn(a, X0, X1, Y0, 1.3), maxIn(b, X0, X1, Y0, 1.3));
      TCanvas c("c_fin", "", 2200, 950);
      c.Divide(2, 1, 0.002, 0.002);
      panel(c.cd(1), a, X0, X1, Y0, 1.3, zm, "final sample, slots #pi^{+}, #pi^{-}, e^{+}, e^{-}: #beta_{new}", -1, true);
      lineLegend(0.115, 0.14, 0.42, 0.255, 0.030);
      panel(c.cd(2), b, X0, X1, Y0, 1.3, zm, "the same tracks: DST #beta (START time)", -1, false, "#beta (DST, START)", false);
      {
         // the sorter's lepton pass needs DST beta >= 0.9
         TLine* l = new TLine(X0, 0.9, X1, 0.9);
         l->SetLineColor(kRed + 1); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
         TLatex t; t.SetTextSize(0.040); t.SetTextColor(kRed + 1); t.DrawLatex(-1.9, 0.86, "#beta = 0.9");
      }
      savePdf(c, "pid_banana_final");

      // per-slot zoom with the slot's own polygons highlighted
      TCanvas cs("c_slots", "", 2000, 1500);
      cs.Divide(2, 2, 0.002, 0.002);
      const int own[4] = {1, 2, 3, 4};
      const double sx0[4] = {0.0, -2.5, 0.0, -1.9}, sx1[4] = {2.5, 0.0, 1.9, 0.0};
      const double sy0[4] = {0.35, 0.35, 0.88, 0.88}, sy1[4] = {1.10, 1.10, 1.10, 1.10};
      const char* ttl[4] = {"final sample, slot #pi^{+}", "final sample, slot #pi^{-}", "final sample, slot e^{+}", "final sample, slot e^{-}"};
      const int pos[4] = {1, 2, 3, 4};
      for (int k = 0; k < 4; ++k) {
         TH2D* h = toGeV(get2(fF, Form("hBN_fat_%s", slot4[k])), Form("hSl%d", k), k < 2 ? 2 : 1, k < 2 ? 2 : 1);
         TVirtualPad* p = cs.cd(pos[k]);
         panel(p, h, sx0[k], sx1[k], sy0[k], sy1[k], -1, Form("%s  (%.2f M rows)", ttl[k], h->GetEntries() / 1e6), own[k]);
         // also the 2 sigma level of the own species, dotted
         toGraph(gV1[own[k]][1], kSign[own[k]], kCol[own[k]], 3, 2)->Draw("L same");
      }
      savePdf(cs, "pid_banana_final_slots");
   }

   // ===== (e) data v1 versus simulation v6 polygons, and the half width of the 1 sigma polygon
   FILE* tp = fopen(Form("%spid_polygons.txt", kTab), "w");
   {
      TCanvas c("c_poly", "", 2100, 1300);
      c.Divide(3, 2, 0.002, 0.002);
      const double px0[5] = {0.0, 0.0, 0.0, 0.0, 0.0}, px1[5] = {4.0, 2.5, 2.5, 1.95, 1.95};
      const double py0[5] = {0.0, 0.33, 0.33, 0.89, 0.89}, py1[5] = {1.06, 1.07, 1.07, 1.09, 1.09};
      for (int k = 0; k < 5; ++k) {
         TVirtualPad* p = c.cd(k + 1);
         p->SetLeftMargin(k < 3 ? 0.14 : 0.17); p->SetRightMargin(0.03); p->SetBottomMargin(0.13); p->SetTopMargin(0.09);   // wider labels for e (verifier fix)
         TH1F* fr = p->DrawFrame(px0[k], py0[k], px1[k], py1[k]);
         fr->GetXaxis()->SetTitle("momentum  [GeV/c]");
         fr->GetYaxis()->SetTitle("#beta_{new}");
         fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.055);
         fr->GetXaxis()->SetLabelSize(0.047); fr->GetYaxis()->SetLabelSize(0.047);
         fr->GetYaxis()->SetTitleOffset(k < 3 ? 1.25 : 1.45);
         const int lv = (k < 3) ? 3 : 5;
         toGraph(gV6[k][lv], 1, kRed + 1, 1, 2)->Draw("L same");
         toGraph(gV6[k][0], 1, kRed + 1, 2, 2)->Draw("L same");
         toGraph(gV1[k][lv], 1, kBlack, 1, 2)->Draw("L same");
         toGraph(gV1[k][0], 1, kBlack, 2, 2)->Draw("L same");
         TLatex t; t.SetNDC(); t.SetTextSize(0.058);
         t.DrawLatex(0.14, 0.935, Form("%s:  %s and 1#sigma", kLab[k], k < 3 ? "3#sigma" : "5#sigma"));
         if (k == 0) {
            TLegend* lg = new TLegend(0.40, 0.17, 0.97, 0.40);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
            TGraph* a = new TGraph(1); a->SetLineColor(kBlack); a->SetLineWidth(2);
            TGraph* b = new TGraph(1); b->SetLineColor(kRed + 1); b->SetLineWidth(2);
            TGraph* d = new TGraph(1); d->SetLineColor(kGray + 2); d->SetLineWidth(2); d->SetLineStyle(2);
            lg->AddEntry(a, "data (betanew_v1)", "l");
            lg->AddEntry(b, "simulation (betanew_v6)", "l");
            lg->AddEntry(d, "dashed: 1#sigma", "l");
            lg->Draw();
         }
      }
      // half width of the 1 sigma polygon versus momentum
      TVirtualPad* p = c.cd(6);
      p->SetLeftMargin(0.17); p->SetRightMargin(0.03); p->SetBottomMargin(0.13); p->SetTopMargin(0.09);
      TH1F* fr = p->DrawFrame(0.0, 0.0, 2.0, 0.030);   // was 0.045: the axis title overlapped the top label (verifier fix)
      fr->GetXaxis()->SetTitle("momentum  [GeV/c]");
      fr->GetYaxis()->SetTitle(L("1#sigma half width in #beta_{new}"));
      fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.047); fr->GetYaxis()->SetLabelSize(0.047);
      fr->GetYaxis()->SetTitleOffset(1.55);
      TLegend* lg = new TLegend(0.55, 0.55, 0.96, 0.90);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
      const int sp3[4] = {0, 1, 2, 3};   // p, pi+, pi-, e (e- polygon = e+ polygon in both files)
      for (int q = 0; q < 4; ++q) {
         const int k = sp3[q];
         for (int set = 0; set < 2; ++set) {
            TCutG* cc = set == 0 ? gV1[k][0] : gV6[k][0];
            TGraph* g = new TGraph();
            const double pmin = (k == 0) ? 150 : (k < 3 ? 80 : 20);
            for (double x = pmin; x <= 2000; x += 10) { double lo, hi; if (yRange(cc, x, lo, hi)) g->SetPoint(g->GetN(), x / 1000., 0.5 * (hi - lo)); }
            g->SetLineColor(kCol[k]); g->SetLineWidth(2); g->SetLineStyle(set == 0 ? 1 : 2);
            g->Draw("L same");
            if (set == 0) lg->AddEntry(g, Form("%s data", k == 3 ? "e^{#pm}" : kLab[k]), "l");
         }
      }
      TGraph* gdash = new TGraph(1); gdash->SetLineColor(kGray + 2); gdash->SetLineStyle(2); gdash->SetLineWidth(2);
      lg->AddEntry(gdash, "dashed: simulation", "l");
      lg->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.058); t.DrawLatex(0.17, 0.935, "1#sigma half width: data vs simulation");
      savePdf(c, "pid_polygons_v1_v6");

      // table: levels, ranges, centre and 1 sigma half width at fixed momenta
      fprintf(tp, "# PID polygons, data betanew_v1 (PAT of the analysed production) and simulation betanew_v6\n");
      fprintf(tp, "# plane (p [MeV/c] uncorrected, beta_new); levels 1/2/2.5/3/3.5/5 sigma in both files; PAT uses p, pi+-: 3 sigma, e+-: 5 sigma\n");
      fprintf(tp, "# per species: number of points, p range; at p = 200/500/1000/1500 MeV/c: centre of the 1 sigma polygon and its half width (data | simulation), half width of the PAT-level polygon\n");
      const double pq[4] = {200, 500, 1000, 1500};
      for (int k = 0; k < 5; ++k) {
         const int lv = (k < 3) ? 3 : 5;
         double xmin = 1e9, xmax = -1e9;
         for (int i = 0; i < gV1[k][0]->GetN(); ++i) { xmin = std::min(xmin, gV1[k][0]->GetX()[i]); xmax = std::max(xmax, gV1[k][0]->GetX()[i]); }
         fprintf(tp, "%-3s PAT level %s  points v1 %d v6 %d  p range %.0f-%.0f MeV/c\n", kSp[k], kPatLv[k], gV1[k][0]->GetN(), gV6[k][0]->GetN(), xmin, xmax);
         for (int q = 0; q < 4; ++q) {
            double l1, h1, l6, h6, L1, H1, L6, H6;
            const bool a = yRange(gV1[k][0], pq[q], l1, h1), b = yRange(gV6[k][0], pq[q], l6, h6);
            const bool A = yRange(gV1[k][lv], pq[q], L1, H1), B = yRange(gV6[k][lv], pq[q], L6, H6);
            if (a && b && A && B)
               fprintf(tp, "    p=%4.0f  centre %.4f | %.4f   1sig half width %.4f | %.4f   %s half width %.4f | %.4f\n", pq[q], 0.5 * (l1 + h1), 0.5 * (l6 + h6),
                       0.5 * (h1 - l1), 0.5 * (h6 - l6), kPatLv[k], 0.5 * (H1 - L1), 0.5 * (H6 - L6));
         }
      }
   }
   fclose(tp);
   {
      // LaTeX version of the polygon table
      FILE* f = fopen(Form("%spid_polygons.tex", kTab), "w");
      fprintf(f, "%% generated by plots/exp_report/pid_draw.C\n\\begin{table}[H]\n\\centering\n");
      fprintf(f, "\\caption{PID polygons in the plane (momentum $p$, $\\beta_{\\mathrm{new}}$): data cut file \\code{betanew\\_v1} "
                 "(used by PAT for the analysed production) and simulation cut file \\code{betanew\\_v6}. Both files hold the levels "
                 "1, 2, 2.5, 3, 3.5 and $5\\sigma$; PAT uses $3\\sigma$ for p and $\\pi^\\pm$ and $5\\sigma$ for $e^\\pm$. "
                 "$p$ is the reconstructed momentum without energy-loss correction. Centre and half width of the $1\\sigma$ polygon "
                 "and half width of the polygon used in PAT, data\\,/\\,simulation. In both files the $e^-$ polygon is identical to the $e^+$ one.}\n");
      fprintf(f, "\\label{tab:pid_polygons}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}l l c c c c c@{}}\n\\toprule\n");
      fprintf(f, "Species & PAT level & $p$ range [GeV/$c$] & $p$ [GeV/$c$] & centre ($1\\sigma$) & half width $1\\sigma$ & half width PAT level \\\\\n\\midrule\n");
      const double pq[4] = {200, 500, 1000, 1500};
      const char* tex[5] = {"p", "$\\pi^+$", "$\\pi^-$", "$e^+$", "$e^-$"};
      for (int k = 0; k < 5; ++k) {
         const int lv = (k < 3) ? 3 : 5;
         double xmin = 1e9, xmax = -1e9;
         for (int i = 0; i < gV1[k][0]->GetN(); ++i) { xmin = std::min(xmin, gV1[k][0]->GetX()[i]); xmax = std::max(xmax, gV1[k][0]->GetX()[i]); }
         bool first = true;
         for (int q = 0; q < 4; ++q) {
            double l1, h1, l6, h6, L1, H1, L6, H6;
            if (!(yRange(gV1[k][0], pq[q], l1, h1) && yRange(gV6[k][0], pq[q], l6, h6) && yRange(gV1[k][lv], pq[q], L1, H1) && yRange(gV6[k][lv], pq[q], L6, H6))) continue;
            fprintf(f, "%s & %s & %s & %.1f & %.4f\\,/\\,%.4f & %.4f\\,/\\,%.4f & %.4f\\,/\\,%.4f \\\\\n",
                    first ? tex[k] : "", first ? (k < 3 ? "$3\\sigma$" : "$5\\sigma$") : "", first ? Form("%.2f--%.1f", xmin / 1000., xmax / 1000.) : "",
                    pq[q] / 1000., 0.5 * (l1 + h1), 0.5 * (l6 + h6), 0.5 * (h1 - l1), 0.5 * (h6 - l6), 0.5 * (H1 - L1), 0.5 * (H6 - L6));
            first = false;
         }
         if (k < 4) fprintf(f, "\\midrule\n");
      }
      fprintf(f, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
      fclose(f);
   }

   // ===== (f) DST (START) beta of the leptons of the final sample
   {
      TCanvas c("c_dst", "", 2000, 1500);
      c.Divide(2, 2, 0.002, 0.002);
      const char* ls[2] = {"ep", "em"};
      const char* ll[2] = {"e^{+}", "e^{-}"};
      double ymax = 0;
      for (int q = 0; q < 2; ++q) ymax = std::max(ymax, ((TH1D*)fF.Get(Form("hDst_fat_%s", ls[q])))->GetMaximum());
      for (int q = 0; q < 2; ++q) {
         TVirtualPad* p = c.cd(q + 1);
         p->SetLogy(); p->SetLeftMargin(0.13); p->SetRightMargin(0.03); p->SetBottomMargin(0.13); p->SetTopMargin(0.09);
         TH1D* a = (TH1D*)fF.Get(Form("hDst_fat_%s", ls[q]));
         TH1D* l = (TH1D*)fF.Get(Form("hDstL_fat_%s", ls[q]));
         TH1D* h = (TH1D*)fF.Get(Form("hDstH_fat_%s", ls[q]));
         a->SetTitle(Form(";%s;%s tracks / 0.005", L("#beta (DST, measured against the START time)").Data(), ll[q]));
         a->GetXaxis()->SetRange(a->GetXaxis()->FindBin(0.0), a->GetNbinsX());   // 0 .. 1.6, no overflow
         a->SetMinimum(0.5); a->SetMaximum(ymax * 60);
         a->GetXaxis()->SetTitleSize(0.052); a->GetYaxis()->SetTitleSize(0.052);
         a->GetXaxis()->SetLabelSize(0.045); a->GetYaxis()->SetLabelSize(0.045);
         a->SetMarkerStyle(20); a->SetMarkerSize(0.6); a->SetMarkerColor(kBlack); a->SetLineColor(kBlack);
         l->SetLineColor(kGreen + 2); l->SetLineWidth(2);
         h->SetLineColor(kOrange + 7); h->SetLineWidth(2);
         a->Draw("P");
         l->Draw("HIST same");
         h->Draw("HIST same");
         a->Draw("P same");
         TLine* ln = new TLine(0.9, 0.5, 0.9, ymax * 2);
         ln->SetLineStyle(2); ln->SetLineWidth(2); ln->Draw();
         TLegend* lg = new TLegend(0.15, 0.68, 0.55, 0.89);
         lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
         lg->AddEntry(a, Form("slot %s, all", ll[q]), "p");
         lg->AddEntry(l, "lepton pass (kIsLepton = 1)", "l");
         lg->AddEntry(h, "hadron pass (kIsLepton = 0)", "l");
         lg->AddEntry(ln, L("lepton pass needs #beta #geq 0.9"), "l");
         lg->Draw();
         TH1D* hc = (TH1D*)fF.Get(Form("hCnt_fat_%s", ls[q]));
         TLatex t; t.SetNDC(); t.SetTextSize(0.048);
         t.DrawLatex(0.13, 0.935, L(Form("final sample, slot %s:  DST #beta < 0.9 for %.1f %% of the tracks", ll[q], 100. * hc->GetBinContent(20) / hc->GetBinContent(1))));
      }
      // per-day share
      {
         TVirtualPad* p = c.cd(3);
         p->SetLeftMargin(0.13); p->SetRightMargin(0.03); p->SetBottomMargin(0.13); p->SetTopMargin(0.09);
         TH1F* fr = p->DrawFrame(33, 0.0, 68, 0.20);
         fr->GetXaxis()->SetTitle("day of 2022 (beam time feb22)");
         fr->GetYaxis()->SetTitle(L("share of tracks with DST #beta < 0.9"));
         fr->GetXaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetTitleSize(0.052);
         fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
         const int mk[2] = {20, 24};
         const int cl[2] = {kBlack, kGray + 2};
         TLegend* lg = new TLegend(0.16, 0.74, 0.60, 0.88);
         lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
         for (int q = 0; q < 2; ++q) {
            TH1D* n = (TH1D*)fF.Get(Form("hDay_fat_%s", ls[q]));
            TH1D* lo = (TH1D*)fF.Get(Form("hDayLow_fat_%s", ls[q]));
            TGraphErrors* g = new TGraphErrors();
            for (int b = 1; b <= n->GetNbinsX(); ++b) {
               const double N = n->GetBinContent(b);
               if (N <= 0) continue;
               const double r = lo->GetBinContent(b) / N;
               const int i = g->GetN();
               g->SetPoint(i, n->GetXaxis()->GetBinCenter(b) + (q == 0 ? -0.12 : 0.12), r);
               g->SetPointError(i, 0, std::sqrt(r * (1 - r) / N));
            }
            g->SetMarkerStyle(mk[q]); g->SetMarkerColor(cl[q]); g->SetLineColor(cl[q]); g->SetMarkerSize(1.1);
            g->Draw("P same");
            lg->AddEntry(g, Form("slot %s", ll[q]), "p");
         }
         lg->Draw();
         TLatex t; t.SetNDC(); t.SetTextSize(0.048);
         t.DrawLatex(0.13, 0.935, "final sample: share per day (days 034-053, 059-067)");
      }
      // DST time of flight of the leptons near zero
      {
         TVirtualPad* p = c.cd(4);
         p->SetLogy(); p->SetLeftMargin(0.13); p->SetRightMargin(0.03); p->SetBottomMargin(0.13); p->SetTopMargin(0.09);
         TTree* t = (TTree*)fF.Get("tag_us");
         TH1D* ha = new TH1D("hTofA", ";t_{DST} of the e^{#pm} tracks  [ns];tracks / 0.1 ns", 160, -2, 14);
         TH1D* hh = new TH1D("hTofH", "", 160, -2, 14);
         t->Draw("tof1>>hTofA", "", "goff"); t->Draw("tof2>>+hTofA", "", "goff");
         t->Draw("tof1>>hTofH", "k1==0", "goff"); t->Draw("tof2>>+hTofH", "k2==0", "goff");
         ha->SetMinimum(0.5); ha->SetMaximum(ha->GetMaximum() * 3000);   // headroom: the legend must not cover the peak (verifier fix)
         ha->GetXaxis()->SetTitleSize(0.052); ha->GetYaxis()->SetTitleSize(0.052);
         ha->GetXaxis()->SetLabelSize(0.045); ha->GetYaxis()->SetLabelSize(0.045);
         ha->SetMarkerStyle(20); ha->SetMarkerSize(0.6); ha->SetMarkerColor(kBlack); ha->SetLineColor(kBlack);
         hh->SetLineColor(kOrange + 7); hh->SetLineWidth(2);
         ha->Draw("P");
         hh->Draw("HIST same");
         TLine* ln = new TLine(0, 0.5, 0, ha->GetMaximum() / 30.);
         ln->SetLineStyle(2); ln->SetLineWidth(2); ln->SetLineColor(kRed + 1); ln->Draw();
         TLegend* lg = new TLegend(0.15, 0.70, 0.62, 0.89);
         lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
         lg->AddEntry(ha, "e^{+} and e^{-} slots, all", "p");
         lg->AddEntry(hh, "hadron pass (kIsLepton = 0)", "l");
         lg->AddEntry(ln, L("hadron pass needs t_{DST} #geq 0 (#beta #geq 0)"), "l");
         lg->Draw();
         TLatex tx; tx.SetNDC(); tx.SetTextSize(0.048);
         tx.DrawLatex(0.13, 0.935, "final sample: DST time of flight of the leptons");
      }
      savePdf(c, "pid_dstbeta");
   }

   // ===== tables: samples and fractions inside the polygons
   {
      struct S { TFile* f; std::vector<std::string> h; std::vector<std::string> u; int sign; int own; std::string sample, slot, mass; };
      std::vector<S> rows = {
         {&fHH, {"hSlot_pip_noring"}, {"hUniq_pip"}, +1, 1, "H: hadrons, no ring", "$\\pi^+$", "$\\pi$"},
         {&fHH, {"hSlot_pim_noring"}, {"hUniq_pim"}, -1, 2, "", "$\\pi^-$", "$\\pi$"},
         {&fPP, {"hSlot_p1_noring", "hSlot_p2_noring"}, {"hUniq_p1", "hUniq_p2"}, +1, 0, "", "p (both)", "p"},
         {&fEE, {"hSlot_ep"}, {"hUniq_ep"}, +1, 3, "L: lepton slots", "$e^+$", "$e$"},
         {&fEE, {"hSlot_em"}, {"hUniq_em"}, -1, 4, "", "$e^-$", "$e$"},
         {&fPE, {"hSlot_pip"}, {"hUniq_pip"}, +1, 1, "L: all slots", "$\\pi^+$", "$\\pi$"},
         {&fPE, {"hSlot_pim"}, {"hUniq_pim"}, -1, 2, "", "$\\pi^-$", "$\\pi$"},
         {&fPE, {"hSlot_ep"}, {"hUniq_ep"}, +1, 3, "", "$e^+$", "$e$"},
         {&fPE, {"hSlot_em"}, {"hUniq_em"}, -1, 4, "", "$e^-$", "$e$"},
         {&fBN, {"hSlot_ep"}, {"hUniq_ep"}, +1, 3, "betanew, PID applied", "$e^+$", "$e$"},
         {&fBN, {"hSlot_em"}, {"hUniq_em"}, -1, 4, "", "$e^-$", "$e$"},
      };
      const char* hypo[11] = {"$\\pi^+\\pi^-$", "", "pp", "$e^+e^-$", "", "$\\pi^+\\pi^-e^+e^-$", "", "", "", "$e^+e^-$", ""};
      FILE* ft = fopen(Form("%spid_samples.txt", kTab), "w");
      FILE* fx = fopen(Form("%spid_samples.tex", kTab), "w");
      fprintf(ft, "# Samples and fractions inside the data PID polygons (betanew_v1). Columns: sample | hypothesis | slot | track entries | unique tracks |\n");
      fprintf(ft, "#   inside p 3sig | inside pi(same charge) 3sig | inside e(same charge) 5sig | own species 1sig | own 2sig | own 3sig\n");
      fprintf(ft, "# banana samples: bin centres (10 MeV/c x 0.0025), denominator = all entries incl. under/overflow (fixed by the verifier); final sample: exact test on each track; '-' = not applicable\n");
      fprintf(fx, "%% generated by plots/exp_report/pid_draw.C\n\\begin{table}[H]\n\\centering\n");
      fprintf(fx, "\\caption{Track samples in the plane ($p$, $\\beta_{\\mathrm{new}}$) and the fraction of the tracks of each slot inside the data "
                  "PID polygons (\\code{betanew\\_v1}; p and $\\pi$ at $3\\sigma$, $e$ at $5\\sigma$, polygon of the same charge as the slot) and inside "
                  "the $1\\sigma$ and $2\\sigma$ polygons of the slot's own species. Track entries count one track per row (a track used in several "
                  "combinations is counted several times; its $\\beta_{\\mathrm{new}}$ differs between rows); unique tracks count each track once per "
                  "event (for H they include the 0.3\\,\\%% of tracks with a ring). All samples: event vertex $z>-500$\\,mm. "
                  "H: June 2026 hadron production without PID cut, all rows (all isBest values, no trigger condition), days 034, 045, 049, 062, "
                  "first 40\\,M rows per file, tracks without RICH ring. L: June 2026 lepton production without PID cut, all rows; $e^+e^-$ "
                  "hypothesis days 036, 045, 050, 062; $\\pi^+\\pi^-e^+e^-$ hypothesis 29 days. betanew: the analysed production (PID applied), "
                  "$e^+e^-$ hypothesis, days 036, 045, 050, 062. For these samples the fraction is computed at the bin centres "
                  "(10\\,MeV/$c$ $\\times$ 0.0025); tracks outside the histogram range ($\\beta_{\\mathrm{new}}$ outside 0--1.3, which happens mainly for the reference track of the hypothesis: 17\\,\\%% of the $\\pi^-$ in $\\pi^+\\pi^-$, 13\\,\\%% of the $e^-$ in $\\pi^+\\pi^-e^+e^-$) count as outside. Final sample: betanew production, $\\pi^+\\pi^-e^+e^-$ rows kept by FAT (isBest $\\geq 0$, PT3), "
                  "29 days; exact test per track.}\n");
      fprintf(fx, "\\label{tab:pid_samples}\n\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{@{}l l l r r r r r r r@{}}\n\\toprule\n");
      fprintf(fx, "Sample & hyp. & slot & entries & unique & in p & in $\\pi$ & in $e$ & own $1\\sigma$ & own $2\\sigma$ \\\\\n\\midrule\n");
      for (size_t i = 0; i < rows.size(); ++i) {
         S& r = rows[i];
         TH2D* h = (TH2D*)get2(*r.f, r.h[0].c_str())->Clone(Form("hTab%zu", i));
         for (size_t j = 1; j < r.h.size(); ++j) h->Add(get2(*r.f, r.h[j].c_str()));
         double nu = 0;
         for (auto& u : r.u) nu += get2(*r.f, u.c_str())->GetEntries();
         const double n = h->GetEntries();
         double fin[3] = {-1, -1, -1};   // p, pi, e of the same charge
         if (r.sign > 0) {
            fin[0] = fracInside(h, gV1[0][3], +1);
            fin[1] = fracInside(h, gV1[1][3], +1);
            fin[2] = fracInside(h, gV1[3][5], +1);
         } else {
            fin[1] = fracInside(h, gV1[2][3], -1);
            fin[2] = fracInside(h, gV1[4][5], -1);
         }
         const int sg = kSign[r.own];
         const double o1 = fracInside(h, gV1[r.own][0], sg), o2 = fracInside(h, gV1[r.own][1], sg), o3 = fracInside(h, gV1[r.own][3], sg);
         fprintf(ft, "%-45s %-22s %-10s %11.0f %11.0f  %6s %6s %6s  %6s %6s %6s\n", r.sample.empty() ? "\"" : r.sample.c_str(), hypo[i], r.slot.c_str(), n, nu,
                 fin[0] < 0 ? "-" : pct(fin[0]).Data(), pct(fin[1]).Data(), pct(fin[2]).Data(), pct(o1).Data(), pct(o2).Data(), pct(o3).Data());
         fprintf(fx, "%s & %s & %s & %.0f & %.0f & %s & %s & %s & %s & %s \\\\\n", r.sample.c_str(), hypo[i], r.slot.c_str(), n, nu,
                 fin[0] < 0 ? "--" : pct(fin[0]).Data(), pct(fin[1]).Data(), pct(fin[2]).Data(), pct(o1).Data(), pct(o2).Data());
         if (i == 2 || i == 4 || i == 8) fprintf(fx, "\\midrule\n");
      }
      fprintf(fx, "\\midrule\n");
      // final sample: exact counters (bins: 1 n, 2 n_unique, 3..8 in own 1..5 sigma, 15..19 in p, pip, pim, ep, em)
      const char* fsl[4] = {"pip", "pim", "ep", "em"};
      const char* ftx[4] = {"$\\pi^+$", "$\\pi^-$", "$e^+$", "$e^-$"};
      const int fsg[4] = {+1, -1, +1, -1};
      for (int k = 0; k < 4; ++k) {
         TH1D* c = (TH1D*)fF.Get(Form("hCnt_fat_%s", fsl[k]));
         const double n = c->GetBinContent(1), nu = c->GetBinContent(2);
         const double inP = c->GetBinContent(15) / n, inPi = c->GetBinContent(fsg[k] > 0 ? 16 : 17) / n, inE = c->GetBinContent(fsg[k] > 0 ? 18 : 19) / n;
         const double o1 = c->GetBinContent(3) / n, o2 = c->GetBinContent(4) / n, o3 = c->GetBinContent(6) / n;
         // bin-centre cross-check on the filled histogram
         const int own = k + 1;
         const double o1b = fracInside(get2(fF, Form("hBN_fat_%s", fsl[k])), gV1[own][0], kSign[own]);
         fprintf(ft, "%-45s %-22s %-10s %11.0f %11.0f  %6s %6s %6s  %6s %6s %6s   (1sig at bin centres %s)\n", k == 0 ? "final sample, betanew prod. after PID" : "\"",
                 k == 0 ? "$\\pi^+\\pi^-e^+e^-$" : "", ftx[k], n, nu, fsg[k] > 0 ? pct(inP, 2).Data() : "-", pct(inPi, 2).Data(), pct(inE, 2).Data(),
                 pct(o1).Data(), pct(o2).Data(), pct(o3, 2).Data(), pct(o1b).Data());
         fprintf(fx, "%s & %s & %s & %.0f & %.0f & %s & %s & %s & %s & %s \\\\\n", k == 0 ? "final sample" : "", k == 0 ? "$\\pi^+\\pi^-e^+e^-$" : "",
                 ftx[k], n, nu, fsg[k] > 0 ? pct(inP).Data() : "--", pct(inPi).Data(), pct(inE).Data(), pct(o1).Data(), pct(o2).Data());
      }
      fprintf(fx, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
      fclose(ft);
      fclose(fx);
   }

   // ===== table: DST beta < 0.9 of the lepton slots per selection
   {
      const char* sel[6] = {"fat", "pre", "pre_best", "pre_vtx", "peak", "bump"};
      const char* stx[6] = {"FAT sample (isBest $\\geq 0$, $z>-500$, PT3)", "pre ($+$ cor\\_ok, $z>-160$)", "pre, isBest $=1$",
                            "pre $+$ both pions in vertex fit", "peak nominal, M4 in 0.5336--0.5542", "bump nominal, M4 in 0.4550--0.5225"};
      FILE* ft = fopen(Form("%spid_dstbeta.txt", kTab), "w");
      FILE* fx = fopen(Form("%spid_dstbeta.tex", kTab), "w");
      fprintf(ft, "# Lepton slots of the final sample: share with DST (START-based) beta < 0.9 = tracks that entered through the hadron pass of the track sorter\n");
      fprintf(ft, "# columns: selection slot entries | share beta<0.9 (entries) | share (unique tracks) | share kIsLepton==0 | beta<0.9 with kIsLepton==1 | share beta>1 | rows evtStart==0 | share beta<0.9 in evtStart==0 rows | beta<0 | tof>60\n");
      fprintf(fx, "%% generated by plots/exp_report/pid_draw.C\n\\begin{table}[H]\n\\centering\n");
      fprintf(fx, "\\caption{Leptons of the analysed $\\pi^+\\pi^-e^+e^-$ rows whose START-based DST velocity is below the lepton-pass "
                  "limit of the track sorter ($\\beta<0.9$): they entered through the hadron pass (sorter flag kIsLepton $=0$). Shares of "
                  "the track entries (one per row) and of the unique tracks (each track once per event). M4 from energy-loss corrected "
                  "momenta; peak nominal: $\\mathrm{oa_{lab}}<55^\\circ$, $M(\\pi\\pi)>0.30$\\,GeV/$c^2$; bump nominal: $\\mathrm{oa_{lab}}<40^\\circ$, "
                  "$0.30<M(\\pi\\pi)<0.41$\\,GeV/$c^2$; data windows in GeV/$c^2$. isBest $\\geq 0$ unless stated.}\n");
      fprintf(fx, "\\label{tab:pid_dstbeta}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}l r r r r r@{}}\n\\toprule\n");
      fprintf(fx, "Selection & rows & $e^+$: $\\beta<0.9$ [\\%%] & $e^-$: $\\beta<0.9$ [\\%%] & $e^+$ unique [\\%%] & $e^-$ unique [\\%%] \\\\\n\\midrule\n");
      for (int s = 0; s < 6; ++s) {
         double sh[2], shu[2], err[2];
         double nrows = 0;
         for (int q = 0; q < 2; ++q) {
            TH1D* c = (TH1D*)fF.Get(Form("hCnt_%s_%s", sel[s], q == 0 ? "ep" : "em"));
            const double n = c->GetBinContent(1), nu = c->GetBinContent(2);
            nrows = n;
            sh[q] = c->GetBinContent(20) / n; shu[q] = c->GetBinContent(21) / nu; err[q] = std::sqrt(sh[q] * (1 - sh[q]) / n);
            fprintf(ft, "%-9s %-2s %9.0f | %.4f +- %.4f | %.4f | %.4f | %.0f | %.4f | %.4f | %.4f | %.0f | %.0f\n", sel[s], q == 0 ? "ep" : "em", n, sh[q], err[q], shu[q],
                    c->GetBinContent(22) / n, c->GetBinContent(24), c->GetBinContent(25) / n, c->GetBinContent(26) / n,
                    c->GetBinContent(26) > 0 ? c->GetBinContent(27) / c->GetBinContent(26) : 0., c->GetBinContent(29), c->GetBinContent(30));
         }
         fprintf(fx, "%s & %.0f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f & %.2f \\\\\n", stx[s], nrows, 100 * sh[0], 100 * err[0], 100 * sh[1], 100 * err[1], 100 * shu[0], 100 * shu[1]);
      }
      // rows with at least one / both leptons below 0.9 (unlike-sign and like-sign FAT rows)
      const char* tg[3] = {"tag_us", "tag_pp", "tag_mm"};
      for (int i = 0; i < 3; ++i) {
         TTree* t = (TTree*)fF.Get(tg[i]);
         const double n = t->GetEntries();
         const double any = t->GetEntries("b1<0.9||b2<0.9"), both = t->GetEntries("b1<0.9&&b2<0.9");
         const double l1 = t->GetEntries("b1<0.9"), l2 = t->GetEntries("b2<0.9");
         const double s0 = t->GetEntries("evtStart==0"), s0any = t->GetEntries("evtStart==0&&(b1<0.9||b2<0.9)");
         fprintf(ft, "%s rows %.0f : lepton share beta<0.9 %.4f (slot1 %.4f slot2 %.4f); rows with >=1 lepton below %.4f, both %.4f; evtStart==0 rows %.4f, of them >=1 below %.4f\n",
                 tg[i], n, (l1 + l2) / (2 * n), l1 / n, l2 / n, any / n, both / n, s0 / n, s0 > 0 ? s0any / s0 : 0.);
      }
      fprintf(fx, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
      fclose(ft);
      fclose(fx);
   }
   // ===== estimate of the tracks removed by the hadron-pass condition (DST tof >= 0, i.e. beta >= 0):
   // lepton DST tof near 0 fitted with an exponential (and a straight line) in 0..1.5 ns and integrated below 0
   {
      TTree* t = (TTree*)fF.Get("tag_us");
      TH1D h("hTofLow", "", 30, 0., 3.);
      t->Draw("tof1>>hTofLow", "", "goff"); t->Draw("tof2>>+hTofLow", "", "goff");
      const double nlep = 2. * t->GetEntries();
      TF1 fe("fe", "expo", 0., 1.5), fl("fl", "pol1", 0., 1.5);
      h.Fit(&fe, "QR0L"); h.Fit(&fl, "QR0L");
      const double bw = h.GetBinWidth(1);
      const double lostExp = std::exp(fe.GetParameter(0)) / fe.GetParameter(1) / bw;            // integral of A exp(b t) from -inf to 0
      const double a0 = fl.GetParameter(0), a1 = fl.GetParameter(1);
      const double lostLin = (a1 > 0) ? 0.5 * a0 * a0 / a1 / bw : -1;                             // triangle down to the zero crossing
      TH1D hh("hTofHigh", "", 70, 0., 70.);
      t->Draw("tof1>>hTofHigh", "", "goff"); t->Draw("tof2>>+hTofHigh", "", "goff");
      const double n4060 = hh.Integral(hh.FindBin(40.5), hh.FindBin(59.5));
      FILE* f = fopen(Form("%spid_sorter_loss.txt", kTab), "w");
      fprintf(f, "# Tracks removed by the hadron pass of the DST track sorter (beta >= 0 && tof <= 60 ns), estimated from the leptons of the final sample\n");
      fprintf(f, "# (FAT rows of PipPimEpEm_ID: isBest >= 0, vz > -500, PT3; %.0f lepton entries = 2 per row)\n", nlep);
      fprintf(f, "lepton entries with tof in [0,0.1) ns: %.0f ; [0,1) ns: %.0f ; [0,3) ns: %.0f\n", h.GetBinContent(1), h.Integral(1, 10), h.Integral());
      fprintf(f, "exponential fit 0-1.5 ns: slope %.3f /ns -> extrapolated entries below tof = 0: %.0f = %.5f of the lepton entries\n", fe.GetParameter(1), lostExp, lostExp / nlep);
      fprintf(f, "straight-line fit 0-1.5 ns: zero at %.3f ns -> extrapolated entries below tof = 0: %.0f = %.5f of the lepton entries\n", -a0 / a1, lostLin, lostLin / nlep);
      fprintf(f, "tof in 40-60 ns: %.0f entries (%.1f per ns) -> a flat continuation beyond 60 ns adds about %.1f entries per ns of the START window\n", n4060, n4060 / 20., n4060 / 20.);
      fprintf(f, "rows: the START shift is common to all tracks of an event, the lepton (fastest track) leaves first; the row loss is of the order of the lepton estimate\n");
      fclose(f);
      printf("sorter loss estimate: exp %.0f (%.5f), lin %.0f (%.5f) of %.0f lepton entries\n", lostExp, lostExp / nlep, lostLin, lostLin / nlep, nlep);
   }

   // ===== protons in the pi+ slot of the final sample, from the DST (START) beta:
   // per momentum bin, the share of pi+ slot tracks below the midpoint between the pion and proton velocity,
   // minus the same share of the pi- slot (pion reference without antiprotons), divided by (1 - pi- share)
   {
      const char* sel[4] = {"fat", "pre", "peak", "bump"};
      const double edges[6] = {500, 1000, 1500, 1900, 2200, 2500};   // the pi polygons end at 2500 MeV/c
      const double mpi = 139.57, mp = 938.27;
      FILE* f = fopen(Form("%spid_pip_protons.txt", kTab), "w");
      FILE* x = fopen(Form("%spid_pip_protons.tex", kTab), "w");
      fprintf(f, "# Protons in the pi+ slot (they pass the pi+ polygon at high momentum, where beta_new of a proton measured with the pion mass lies\n");
      fprintf(f, "# inside the pi+ band). Estimate from the DST beta: f(pi+) = share below the pion/proton midpoint, f(pi-) the same for pi- (reference);\n");
      fprintf(f, "# proton share = (f(pi+) - f(pi-)) / (1 - f(pi-)). Columns: selection, p bin [MeV/c], pi+ tracks, share of all pi+, f(pi+), f(pi-), proton share, protons, protons / all pi+\n");
      fprintf(x, "%% generated by plots/exp_report/pid_draw.C\n\\begin{table}[H]\n\\centering\n");
      fprintf(x, "\\caption{Protons in the $\\pi^+$ slot of the analysed rows (isBest $\\geq 0$, $z>-500$\\,mm, PT3), estimated from the START-based DST velocity: "
                 "share of $\\pi^+$-slot tracks with $\\beta_{\\mathrm{DST}}$ below the midpoint of the pion and proton velocities, minus the same share for "
                 "the $\\pi^-$ slot (no antiprotons), divided by one minus the $\\pi^-$ share. Momentum $p$ without energy-loss correction.}\n");
      fprintf(x, "\\label{tab:pid_pip_protons}\n\\footnotesize\n\\begin{tabular}{@{}l r r r r r@{}}\n\\toprule\n");
      fprintf(x, "$p$ [GeV/$c$] & $\\pi^+$ tracks & share of slot & $f(\\pi^+)$ & $f(\\pi^-)$ & protons [\\%%] \\\\\n\\midrule\n");
      for (int s = 0; s < 4; ++s) {
         TH2D* hp = get2(fF, Form("hBD_%s_pip", sel[s]));
         TH2D* hm = get2(fF, Form("hBD_%s_pim", sel[s]));
         const double tot = hp->GetEntries();
         double sumProt = 0;
         for (int i = 0; i < 5; ++i) {
            const double a = edges[i], b = edges[i + 1], pc = 0.5 * (a + b);
            const double mid = 0.5 * (pc / std::sqrt(pc * pc + mpi * mpi) + pc / std::sqrt(pc * pc + mp * mp));
            const int ya = hp->GetYaxis()->FindBin(mid);
            const int xa = hp->GetXaxis()->FindBin(a + 1), xb = hp->GetXaxis()->FindBin(b - 1);
            const int ma = hm->GetXaxis()->FindBin(-b + 1), mb = hm->GetXaxis()->FindBin(-a - 1);
            const double Np = hp->Integral(xa, xb, 1, hp->GetNbinsY()), Lp = hp->Integral(xa, xb, 1, ya - 1);
            const double Nm = hm->Integral(ma, mb, 1, hm->GetNbinsY()), Lm = hm->Integral(ma, mb, 1, ya - 1);
            if (Np <= 0) continue;
            const double fp = Lp / Np, fm = Nm > 0 ? Lm / Nm : 0;
            const double ex = (Nm > 20) ? (fp - fm) / (1 - fm) : fp;
            sumProt += ex * Np;
            fprintf(f, "%-5s %4.0f-%4.0f  %9.0f  %.4f  f(pi+) %.4f  f(pi-) %.4f (%6.0f pi-)  proton share %.3f  protons %8.0f  of all pi+ %.5f\n", sel[s], a, b, Np, Np / tot, fp, fm, Nm, ex, ex * Np, ex * Np / tot);
            if (s == 0) fprintf(x, "%.1f--%.1f & %.0f & %.4f & %.3f & %.3f & %.1f \\\\\n", a / 1000., b / 1000., Np, Np / tot, fp, fm, 100. * ex);
         }
         fprintf(f, "%-5s total protons %.0f of %.0f pi+ slot tracks = %.4f\n", sel[s], sumProt, tot, sumProt / tot);
         if (s == 0) fprintf(x, "\\midrule\nall & %.0f & & & & %.2f \\\\\n", tot, 100. * sumProt / tot);
      }
      fprintf(x, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
      fclose(f);
      fclose(x);
   }

   TColor::InvertPalette();
   printf("written %spid_*.pdf/png and %spid_*.{txt,tex}\n", kFig, kTab);
}
