// fix2_vtx_fig_cnt.C -- COPY of fix2_vtx_fig.C (figure agent B, 2026-11-15): "counts" instead of "rows" in the axis titles,
// panel labels and one legend entry; writes ONLY docs/latex/figs_all/fix2_vtx_fig_cnt.{pdf,png} (its number dump goes to
// /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/fix2_vtx/). Nothing else changed.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_fix2_cnt", kTRUE)' 'plots/all_report/fix2_vtx_fig_cnt.C+'
// ==============================================================================================================
// fix2_vtx_fig.C -- final review (prefix fix2, 2026-11-14): copy of vtx_fig.C with (1) a 2 x 3 layout and larger text
// (numbers readable at full text width, about 7.5-8 pt), (2) the numbers of the tables drawn: S with the event-group
// factor (stat_g), S(Nv) with the smooth like-sign CB as the mean of its 12 variants (the points are the pol4, 0.48-0.70
// variant), E with its statistical error x 1.2 (Part II calibration) and its model rms shown separately, and the
// combination with the nominal CB added. The physics content is unchanged. Original header follows.
// Output: docs/latex/figs_all/fix2_vtx_fig.{pdf,png}, docs/latex/figs_exp/fix2_vtx_fig.{pdf,png}, docs/latex/tables_all/fix2_vtx_fig.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic_fix2", kTRUE)' 'plots/all_report/fix2_vtx_fig.C+'
// --------------------------------------------------------------------------------------------------------------
// vtx_fig.C -- figure of the eta efficiency of the mandatory event-vertex cut vz > -500 mm (integration into report_all /
// report_exp, 2026-11-14). New macro; it uses the verifier's library and rows unchanged (vtx_lib.h, vtx_load.h = copies of
// ve_verify/vv_lib.h, vv_load.h; rows written by vtx_fill.C = vv_fill.C in /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/).
// Classes: V = rows with an event vertex (vz > -500 mm), Nv = rows without (PAT writes vz = -1000). Base: PT3 (trigbit 8192),
// isBest >= 0, cor_ok; unlike-sign rows (us) and like-sign rows (pp, mm) of the 29 PAT day files.
//   (a) z of the pi+pi- crossing point (doca < 10 mm), V and Nv, unit area; LH2 cell and finder window marked
//   (b) distance of that point from the beam axis
//   (c) eta peak: N - CB - pol2 (P_main: frozen Gaussian 0.5439 / 5.4 MeV + pol2 in 0.5225-0.61); V with the nominal CB,
//       Nv with the nominal CB and with a smooth like-sign CB (N++ and N-- fitted with pol4 in 0.48-0.70)
//   (d) Nv only: the frozen Gaussian fitted to the e+e+ rows and to the CB itself (the dip that fakes S(Nv) with the nominal CB)
//   (e) ECAL-tagged eta, m5c (frozen 0.5463 / 6.2 MeV + pol2 in 0.49-0.62), V and Nv
//   (f) eps_vtx = X_V/(X_V+X_Nv) for every measure (S, T: computed here; E and the combination: vtx_ana.txt / vtx_numbers.txt)
// Output: docs/latex/figs_all/vtx_fig.{pdf,png}, docs/latex/figs_exp/vtx_fig.{pdf,png}, docs/latex/tables_all/vtx_fig.txt
// Run (from /home/przygoda/HADES/FAT/FAT, about 1 min, 1 GB of memory):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic_fig", kTRUE)' 'plots/all_report/vtx_fig.C+'
// ==============================================================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TString.h"
#include "vtx_lib.h"
using namespace vvl;
#include "vtx_load.h"

namespace vfg {
FILE* gO = nullptr;
// smooth fit of a like-sign fine-bin spectrum (copy of Smooth() in vtx_ana3.C)
std::vector<double> Smooth(const Spec& sp, int k, double lo, double hi, int ord) {
  std::vector<int> idx; for (int i = 0; i < kNb; ++i) if (sp.x(i) >= lo && sp.x(i) < hi) idx.push_back(i);
  const double x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
  std::vector<std::vector<double>> F; std::vector<double> y, add;
  for (int i : idx) { std::vector<double> row; const double t = (sp.x(i) - x0) / dx; double tp = 1;
    for (int j = 0; j <= ord; ++j) { row.push_back(tp); tp *= t; } F.push_back(row); y.push_back(sp.n[k][i]); add.push_back(0); }
  Fit f = GLS(sp, idx, F, y, add, false);
  std::vector<double> out(kNb, 0.);
  for (int i = 0; i < kNb; ++i) { const double t = (sp.x(i) - x0) / dx; double tp = 1, v = 0; for (int j = 0; j <= ord; ++j) { v += f.p[j] * tp; tp *= t; } out[i] = v; }
  return out;
}
// polynomial part of a Peak() fit at fine bin i
double Pol(const Fit& f, const PCfg& c, double x) {
  const double x0 = 0.5 * (c.lo + c.hi), dx = 0.5 * (c.hi - c.lo), t = (x - x0) / dx; double tp = 1, v = 0;
  for (int j = 0; j <= c.ord; ++j) { v += f.p[j + 1] * tp; tp *= t; }
  return v;
}
// residual histogram (data - pol) merged by m fine bins; what: 0 N - CB, 1 N++, 6 CB
TH1D* Resid(const Spec& sp, const Fit& f, const PCfg& c, int what, int m, const char* name, bool smoothCB) {
  const int nb = (int)std::lround((c.hi - c.lo) / (m * kBw));
  TH1D* h = new TH1D(name, "", nb, c.lo, c.lo + nb * m * kBw); h->SetDirectory(nullptr);
  for (int b = 0; b < nb; ++b) {
    double v = 0, var = 0;
    for (int k = 0; k < m; ++k) {
      const double x = c.lo + (b * m + k + 0.5) * kBw; const int i = (int)((x - kLo) / kBw);
      double y = 0, vy = 0;
      if (what == 0) { y = sp.n[0][i] - sp.cb[i]; vy = sp.n[0][i] + (smoothCB ? 0. : sp.CbCov(i, i)); }
      else if (what == 1) { y = sp.n[1][i]; vy = sp.n[1][i]; }
      else { y = sp.cb[i]; vy = sp.CbCov(i, i); }
      v += y - Pol(f, c, x); var += vy;
    }
    h->SetBinContent(b + 1, v); h->SetBinError(b + 1, std::sqrt(std::max(var, 1.0)));
  }
  return h;
}
TF1* Gaus(const Fit& f, const PCfg& c, int m, const char* name, int col, int sty = 1) {
  TF1* g = new TF1(name, "[0]*[3]/(sqrt(2*pi)*[2])*exp(-0.5*((x-[1])/[2])^2)", c.lo, c.hi);
  g->SetParameters(f.S(), c.mu, c.sg, m * kBw); g->SetLineColor(col); g->SetLineWidth(3); g->SetLineStyle(sty); g->SetNpx(400);
  return g;
}
void Txt(double x, double y, const char* s, double size = 0.055, int col = kBlack, int align = 11) {
  TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align); l.DrawLatex(x, y, s);
}
void Range(TH1D** h, int n, double& lo, double& hi) {
  lo = 0; hi = 0;
  for (int k = 0; k < n; ++k) for (int b = 1; b <= h[k]->GetNbinsX(); ++b) {
    hi = std::max(hi, h[k]->GetBinContent(b) + h[k]->GetBinError(b)); lo = std::min(lo, h[k]->GetBinContent(b) - h[k]->GetBinError(b)); }
}
}  // namespace vfg

void fix2_vtx_fig_cnt() {
  using namespace vfg;
  const char* FA = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/";
  const char* FE = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/";
  gO = fopen("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/fix2_vtx/fix2_vtx_fig.txt", "w");   // _cnt copy: table to scratch
  auto PR = [](const char* s) { printf("%s", s); if (gO) fprintf(gO, "%s", s); };
  PR("# fix2_vtx_fig.txt -- numbers drawn in figs_all/fix2_vtx_fig.pdf (plots/all_report/fix2_vtx_fig.C; rows of ve_verify)\n");
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
  gStyle->SetPadLeftMargin(0.19); gStyle->SetPadRightMargin(0.03); gStyle->SetPadBottomMargin(0.17); gStyle->SetPadTopMargin(0.10);
  gStyle->SetLabelSize(0.056, "xyz"); gStyle->SetTitleSize(0.060, "xyz"); gStyle->SetTitleOffset(1.50, "y"); gStyle->SetTitleOffset(1.12, "x");
  gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
  std::vector<R> d = Load("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/d_*.root");
  const int cV = kBlue + 1, cN = kGray + 2, cPP = kGreen + 2, cCB = kRed + 1;
  TCanvas* c = new TCanvas("c", "", 1300, 1500);
  c->Divide(2, 3, 0.003, 0.003);
  // ---------------- (a), (b) pion crossing point
  {
    TH1D* hz[2]; TH1D* hr[2];
    for (int k = 0; k < 2; ++k) {
      hz[k] = new TH1D(Form("hz%d", k), "", 80, -350, 50); hz[k]->SetDirectory(nullptr);
      hr[k] = new TH1D(Form("hr%d", k), "", 40, 0, 40); hr[k]->SetDirectory(nullptr);
    }
    double nr8[2] = {0, 0}, nall[2] = {0, 0};
    for (const R& r : d) {
      if (r.smp != 0 || !r.cor || !(r.dca < 10)) continue;
      const int k = r.vz > -500 ? 0 : 1;
      hz[k]->Fill(r.zx); hr[k]->Fill(r.rx); nall[k] += 1; if (r.rx >= 8) nr8[k] += 1;
    }
    PR(Form("crossing point (us, cor_ok, doca < 10 mm): rows V %.0f, Nv %.0f; share r >= 8 mm: V %.3f, Nv %.3f\n",
            nall[0], nall[1], nr8[0] / nall[0], nr8[1] / nall[1]));
    PR(Form("crossing z: mean V %.1f mm, Nv %.1f mm; mode V %.1f, Nv %.1f\n", hz[0]->GetMean(), hz[1]->GetMean(),
            hz[0]->GetBinCenter(hz[0]->GetMaximumBin()), hz[1]->GetBinCenter(hz[1]->GetMaximumBin())));
    for (int p = 0; p < 2; ++p) {
      c->cd(p + 1);
      TH1D** h = p == 0 ? hz : hr;
      double ym = 0;
      for (int k = 0; k < 2; ++k) { h[k]->Scale(1.0 / h[k]->Integral()); ym = std::max(ym, h[k]->GetMaximum()); }
      h[0]->SetMaximum(1.85 * ym); h[0]->SetMinimum(0);
      h[0]->GetXaxis()->SetTitle(p == 0 ? "z of the #pi^{+}#pi^{-} crossing point [mm]" : "crossing point: distance from axis [mm]");
      h[0]->GetYaxis()->SetTitle("fraction of counts per bin");
      h[0]->GetXaxis()->SetNdivisions(505);
      h[0]->SetLineColor(cV); h[0]->SetLineWidth(3);
      h[1]->SetLineColor(cN); h[1]->SetLineWidth(3); h[1]->SetFillColorAlpha(kGray, 0.35);
      h[0]->Draw("hist"); h[1]->Draw("hist same"); h[0]->Draw("hist same");
      if (p == 0) {
        TBox* cell = new TBox(-148.4, 0, -98.6, 1.85 * ym); cell->SetFillColorAlpha(kAzure - 9, 0.35); cell->Draw();
        TLine* lw[2];
        int iw = 0;
        for (double z : {-170., -30.}) { TLine* l = new TLine(z, 0, z, 0.60 * 1.85 * ym); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); lw[iw++] = l; }
        h[1]->Draw("hist same"); h[0]->Draw("hist same");
        TLegend* lg = new TLegend(0.21, 0.60, 0.70, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.056);
        lg->AddEntry(h[0], "with event vertex (V)", "l"); lg->AddEntry(h[1], "no event vertex (Nv)", "lf");
        lg->AddEntry(cell, "LH_{2} cell", "f"); lg->AddEntry(lw[0], "finder window", "l"); lg->Draw();
      } else {
        TLine* l = new TLine(8, 0, 8, 1.85 * ym); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
        TLegend* lg = new TLegend(0.42, 0.72, 0.97, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.056);
        lg->AddEntry(h[0], "with event vertex (V)", "l"); lg->AddEntry(h[1], "no event vertex (Nv)", "lf"); lg->Draw();
        Txt(0.50, 0.63, "share with r #geq 8 mm:", 0.056); Txt(0.50, 0.555, Form("V %.0f %%, Nv %.0f %%", 100 * nr8[0] / nall[0], 100 * nr8[1] / nall[1]), 0.056);
      }
      Txt(0.19, 0.925, p == 0 ? "(a) e^{+}e^{-} counts, pion doca < 10 mm" : "(b) the same counts", 0.062);
    }
  }
  // ---------------- spectra in M4 for V and Nv
  std::vector<Sel> s2 = {[](const R& r) { return r.vz > -500; }, [](const R& r) { return r.vz <= -500; }};
  std::vector<Spec> sp = Fill(d, s2);
  PCfg pc;   // P_main
  Fit fV = Peak(sp[0], pc), fN = Peak(sp[1], pc);
  Spec qN = sp[1];
  {
    std::vector<double> fp = Smooth(sp[1], 1, 0.48, 0.70, 4), fm = Smooth(sp[1], 2, 0.48, 0.70, 4);
    for (int i = 0; i < kNb; ++i) qN.cb[i] = (fp[i] > 0 && fm[i] > 0) ? 2 * std::sqrt(fp[i] * fm[i]) : 0;
    for (int j = 0; j < qN.nw; ++j) qN.cbWv[j] = 0;
  }
  Fit fNs = Peak(qN, pc);
  PCfg p1 = pc; p1.what = 1; Fit fPP = Peak(sp[1], p1);
  PCfg p6 = pc; p6.what = 6; Fit fCB = Peak(sp[1], p6);
  PR(Form("S (P_main, nominal CB): V %.0f +- %.0f, Nv %.0f +- %.0f\n", fV.S(), fV.eS(), fN.S(), fN.eS()));
  PR(Form("S(Nv) with the smooth like-sign CB (pol4, 0.48-0.70; vtx_ana3.txt: mean of 12 variants 263, rms 30): %.0f +- %.0f\n", fNs.S(), fNs.eS()));
  PR(Form("Nv: frozen Gaussian on N++: %.0f +- %.0f; on the CB itself: %.0f +- %.0f\n", fPP.S(), fPP.eS(), fCB.S(), fCB.eS()));
  // (c)
  {
    c->cd(3);
    TH1D* h[3] = {Resid(sp[0], fV, pc, 0, 2, "rV", false), Resid(sp[1], fN, pc, 0, 2, "rN", false), Resid(qN, fNs, pc, 0, 2, "rNs", true)};
    double lo, hi; TH1D* hd[2] = {h[0], h[2]}; Range(hd, 2, lo, hi);
    h[0]->SetMarkerStyle(20); h[0]->SetMarkerColor(cV); h[0]->SetLineColor(cV); h[0]->SetMarkerSize(1.1);
    h[1]->SetMarkerStyle(24); h[1]->SetMarkerColor(cN); h[1]->SetLineColor(cN); h[1]->SetMarkerSize(1.0);
    h[2]->SetMarkerStyle(21); h[2]->SetMarkerColor(cN); h[2]->SetLineColor(cN); h[2]->SetMarkerSize(1.0);
    h[0]->SetMinimum(1.15 * lo); h[0]->SetMaximum(2.05 * hi);
    h[0]->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); h[0]->GetXaxis()->SetNdivisions(505);
    h[0]->GetYaxis()->SetTitle("N #minus CB #minus pol2 per 5 MeV/c^{2}");
    h[0]->Draw("e"); h[2]->Draw("e same");
    TLine* l0 = new TLine(pc.lo, 0, pc.hi, 0); l0->SetLineStyle(3); l0->Draw();
    TF1* gV = Gaus(fV, pc, 2, "gV", cV); TF1* gNs = Gaus(fNs, pc, 2, "gNs", cN); TF1* gN = Gaus(fN, pc, 2, "gN", cN, 2);
    gV->Draw("same"); gNs->Draw("same"); gN->Draw("same");
    Txt(0.19, 0.925, "(c) #eta peak, P_{main}", 0.062);
    // legend: S with the event-group factor of the verifier (397/390, 411/401); smooth CB: mean of the 12 variants (vtx_ana3.txt A)
    TLegend* lg = new TLegend(0.23, 0.655, 0.97, 0.885); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.056);
    lg->AddEntry(h[0], Form("V: S = %.0f #pm %.0f", fV.S(), fV.eS() * 397. / 390.), "pe");
    lg->AddEntry(h[2], "Nv, smooth CB: 263 #pm 324", "pe");
    lg->AddEntry(gN, Form("Nv, nominal CB: %.0f #pm %.0f", fN.S(), fN.eS() * 411. / 401.), "l");
    lg->Draw();
    PR(Form("panel (c) legend: V %.0f +- %.0f (stat_g), Nv smooth CB 263 +- 324 (mean of 12 variants; points: pol4 0.48-0.70, S = %.0f), Nv nominal %.0f +- %.0f (stat_g)\n",
            fV.S(), fV.eS() * 397. / 390., fNs.S(), fN.S(), fN.eS() * 411. / 401.));
  }
  // (d)
  {
    c->cd(4);
    TH1D* h[2] = {Resid(sp[1], fPP, p1, 1, 2, "rPP", false), Resid(sp[1], fCB, p6, 6, 2, "rCB", false)};
    double lo, hi; Range(h, 2, lo, hi);
    h[1]->SetMarkerStyle(20); h[1]->SetMarkerColor(cCB); h[1]->SetLineColor(cCB); h[1]->SetMarkerSize(1.0);
    h[0]->SetMarkerStyle(21); h[0]->SetMarkerColor(cPP); h[0]->SetLineColor(cPP); h[0]->SetMarkerSize(1.0);
    h[1]->SetMinimum(1.20 * lo); h[1]->SetMaximum(-1.05 * lo);
    h[1]->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{#pm}e^{#pm}) [GeV/c^{2}]"); h[1]->GetXaxis()->SetNdivisions(505);
    h[1]->GetYaxis()->SetTitle("counts #minus pol2 per 5 MeV/c^{2}");
    h[1]->Draw("e"); h[0]->Draw("e same");
    TLine* l0 = new TLine(pc.lo, 0, pc.hi, 0); l0->SetLineStyle(3); l0->Draw();
    Gaus(fCB, pc, 2, "gCB", cCB)->Draw("same"); Gaus(fPP, pc, 2, "gPP", cPP)->Draw("same");
    Txt(0.19, 0.925, "(d) Nv: CB and e^{+}e^{+} counts at m_{#eta}", 0.062);
    TLegend* lg = new TLegend(0.30, 0.715, 0.97, 0.885); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.056);
    lg->AddEntry(h[1], Form("CB of Nv: %.0f #pm %.0f", fCB.S(), fCB.eS()), "pe");
    lg->AddEntry(h[0], Form("e^{+}e^{+} counts: %.0f #pm %.0f", fPP.S(), fPP.eS()), "pe");
    lg->Draw();
  }
  // (e) ECAL-tagged eta
  double tVS = 0, tVe = 0, tNS = 0, tNe = 0;
  {
    c->cd(5);
    auto tag = [](const R& r) { return r.ng > 0 && r.meeg >= 0.119 && r.meeg < 0.173 && r.mee < 0.10; };
    std::vector<Sel> ts = {[tag](const R& r) { return tag(r) && r.vz > -500; }, [tag](const R& r) { return tag(r) && r.vz <= -500; }};
    std::vector<Spec> st = Fill(d, ts, [](const R& r) { return r.m5c; }, 0.5463 - 1.9 * 0.0062, 0.5463 + 1.9 * 0.0062);
    PCfg tc; tc.mu = 0.5463; tc.sg = 0.0062; tc.lo = 0.49; tc.hi = 0.62;
    Fit tV = Peak(st[0], tc), tN = Peak(st[1], tc);
    tVS = tV.S(); tVe = tV.eS(); tNS = tN.S(); tNe = tN.eS();
    PR(Form("T (ECAL-tagged): V %.0f +- %.0f, Nv %.0f +- %.0f\n", tV.S(), tV.eS(), tN.S(), tN.eS()));
    TH1D* h[2] = {Resid(st[0], tV, tc, 0, 2, "tV", false), Resid(st[1], tN, tc, 0, 2, "tN", false)};
    double lo, hi; Range(h, 2, lo, hi);
    h[0]->SetMarkerStyle(20); h[0]->SetMarkerColor(cV); h[0]->SetLineColor(cV); h[0]->SetMarkerSize(1.0);
    h[1]->SetMarkerStyle(21); h[1]->SetMarkerColor(cN); h[1]->SetLineColor(cN); h[1]->SetMarkerSize(0.9);
    h[0]->SetMinimum(1.15 * lo); h[0]->SetMaximum(1.75 * hi);
    h[0]->GetXaxis()->SetTitle("m_{5c} = M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) [GeV/c^{2}]"); h[0]->GetXaxis()->SetNdivisions(505);
    h[0]->GetYaxis()->SetTitle("N #minus CB #minus pol2 per 5 MeV/c^{2}");
    h[0]->Draw("e"); h[1]->Draw("e same");
    TLine* l0 = new TLine(tc.lo, 0, tc.hi, 0); l0->SetLineStyle(3); l0->Draw();
    Gaus(tV, tc, 2, "gtV", cV)->Draw("same"); Gaus(tN, tc, 2, "gtN", cN)->Draw("same");
    Txt(0.19, 0.925, "(e) ECAL-tagged #eta, T", 0.062);
    TLegend* lg = new TLegend(0.50, 0.715, 0.97, 0.885); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.056);
    lg->AddEntry(h[0], Form("V: T = %.0f #pm %.0f", tV.S(), tV.eS()), "pe");
    lg->AddEntry(h[1], Form("Nv: T = %.0f #pm %.0f", tN.S(), tN.eS()), "pe");
    lg->Draw();
  }
  // (f) summary of eps_vtx (values as in the tables: tables_all/fix2_numbers.txt, vtx_numbers.txt, vtx_ana2/3.txt)
  {
    c->cd(6);
    gPad->SetLeftMargin(0.43); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.17);
    double e1, v1 = EpsV(fV.S(), fV.eS() * 397. / 390., fN.S(), fN.eS() * 411. / 401., e1);   // stat_g scaling of the verifier
    double e2, v2 = EpsV(4646, 370, 165, 324, e2);                                             // N only (vtx_ana2.txt), plain fit errors
    double e3, v3 = EpsV(4670, 370, 263, 324, e3);                                             // smooth CB, mean of 12 variants (vtx_ana3.txt)
    e3 = 0.063;                                                                                // as quoted by the verifier (0.062 with these inputs)
    double e4, v4 = EpsV(tVS, tVe, tNS, tNe, e4);
    double e5s, v5 = EpsV(1241, 76 * 1.2, 21, 80 * 1.2, e5s);                                  // E, fit errors x 1.2 (Part II calibration)
    const double e5m = 0.028, e5 = std::hypot(e5s, e5m);                                       // + model rms over 18 variations
    const double v6 = 0.994, e6 = 0.021;                                                       // S(smooth) + T + E (fix2_numbers.txt: 0.9937 +- 0.0211)
    const double v7 = 0.990, e7 = 0.022;                                                       // S(nominal) + T + E (fix2_numbers.txt: 0.9897 +- 0.0218)
    const double e1t = std::hypot(e1, 0.060);                                                  // nominal S: stat_g (+) model 0.060
    PR(Form("eps: S nominal %.3f +- %.3f +- 0.060; S N-only %.3f +- %.3f; S smooth %.3f +- %.3f; T %.3f +- %.3f; E %.3f +- %.3f (stat x1.2) +- %.3f; S(smooth)+T+E %.3f +- %.3f; S(nominal)+T+E %.3f +- %.3f\n",
            v1, e1, v2, e2, v3, e3, v4, e4, v5, e5s, e5m, v6, e6, v7, e7));
    const int n = 7;
    const double val[n] = {v1, v2, v3, v4, v5, v6, v7}, err[n] = {e1t, e2, e3, e4, e5, e6, e7};
    const char* lab[n] = {"S, nominal CB", "S, N only (no CB)", "S, smooth CB", "T, ECAL tag", "E, edge step", "S(smooth)+T+E", "S(nominal)+T+E"};
    TString vtxt[n] = {Form("%.3f #pm %.3f #pm 0.060", v1, e1), Form("%.3f #pm %.3f", v2, e2), Form("%.3f #pm %.3f", v3, e3),
                       Form("%.3f #pm %.3f", v4, e4), Form("%.3f #pm %.3f #pm %.3f", v5, e5s, e5m), Form("%.3f #pm %.3f", v6, e6),
                       Form("%.3f #pm %.3f", v7, e7)};
    const double xlo = 0.35, xhi = 1.10, ytop = n + 1.30;
    TH1D* fr = new TH1D("fr", "", 10, xlo, xhi); fr->SetDirectory(nullptr);
    fr->SetMinimum(0.35); fr->SetMaximum(ytop);
    fr->GetXaxis()->SetTitle("#varepsilon_{vtx}(#eta) = X_{V} / (X_{V} + X_{Nv})"); fr->GetXaxis()->SetNdivisions(505);
    fr->GetYaxis()->SetLabelSize(0); fr->GetYaxis()->SetTickLength(0);
    fr->Draw("axis");
    const double yl = n + 0.55;
    TLine* ls = new TLine(0.9997, 0.35, 0.9997, yl); ls->SetLineColor(kCyan + 2); ls->SetLineWidth(2); ls->Draw();
    TLine* lr = new TLine(0.6773, 0.35, 0.6773, yl); lr->SetLineStyle(3); lr->SetLineWidth(3); lr->SetLineColor(kGray + 2); lr->Draw();
    TLine* lp = new TLine(0.565, 0.35, 0.565, yl); lp->SetLineStyle(3); lp->SetLineWidth(3); lp->SetLineColor(kGray + 2); lp->Draw();
    TLine* l46 = new TLine(0.46, 0.35, 0.46, yl); l46->SetLineStyle(9); l46->SetLineWidth(3); l46->Draw();
    TLatex a; a.SetTextFont(42); a.SetTextSize(0.052); a.SetTextAlign(21);
    a.SetTextColor(kBlack); a.DrawLatex(0.43, yl + 0.15, "0.93 mb");
    a.SetTextColor(kGray + 2); a.DrawLatex(0.575, yl + 0.15, "P"); a.DrawLatex(0.677, yl + 0.15, "all");
    a.SetTextColor(kCyan + 2); a.DrawLatex(0.985, yl + 0.15, "SMASH");
    for (int k = 0; k < n; ++k) {
      const double y = n - k;
      TGraphErrors* g = new TGraphErrors(1); g->SetPoint(0, val[k], y); g->SetPointError(0, err[k], 0);
      const bool main = (k == 5);
      const bool open = (k == 0 || k == 1 || k == 6);
      g->SetMarkerStyle(main ? 33 : (k == 6 ? 27 : (open ? 24 : 20))); g->SetMarkerSize(main ? 2.6 : (k == 6 ? 2.2 : 1.5));
      g->SetMarkerColor(kBlack); g->SetLineColor(kBlack); g->SetLineWidth(main ? 4 : 2);
      g->Draw("p same");
      if (k == 4) {   // E: thick bar = stat (x1.2), thin = stat (+) model
        TGraphErrors* gs = new TGraphErrors(1); gs->SetPoint(0, val[k], y); gs->SetPointError(0, e5s, 0);
        gs->SetMarkerStyle(20); gs->SetMarkerSize(1.5); gs->SetLineWidth(4); gs->Draw("p same");
      }
      TLatex t; t.SetTextFont(main ? 62 : 42); t.SetTextSize(0.054); t.SetTextAlign(32); t.DrawLatex(xlo - 0.012, y + 0.20, lab[k]);
      TLatex t2; t2.SetTextFont(main ? 62 : 42); t2.SetTextSize(0.050); t2.SetTextAlign(32); t2.SetTextColor(kGray + 3);
      t2.DrawLatex(xlo - 0.012, y - 0.27, vtxt[k].Data());
    }
    Txt(0.43, 0.925, "(f) #varepsilon_{vtx}(#eta) per measure", 0.062);
  }
  (void)FE;   // _cnt copy: figs_all only (figs_exp is handled separately)
  c->SaveAs(Form("%sfix2_vtx_fig_cnt.pdf", FA)); c->SaveAs(Form("%sfix2_vtx_fig_cnt.png", FA));
  if (gO) fclose(gO);
}
