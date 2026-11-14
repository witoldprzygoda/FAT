// ============================================================================
// s5_blind_ana_summary.C  --  tables and figures of the blind analysis
// (label blind_ana). Reads ONLY the results written by s5_blind_ana_eta.C and
// s5_blind_ana_x17.C, the half-A templates and (for the side-band figure) the
// X17 pseudo-data of three toys. No sealed file.
// Output: tables/s5_blind_ana_eta_summary.txt, tables/s5_blind_ana_x17_summary.txt,
//         figs/s5_blind_ana_*.pdf/.png
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_ana_build", kTRUE)' plots/eta_sim/s5_blind_ana_summary.C+
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TTree.h"

#include "s5_blind_ana_fit.h"
#include "s5_blind_ana_tmpl.h"

using namespace bla;

namespace bsum {

const char* const kRes = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/figs";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables";
const char* kSelN[4] = {"PEAK", "DIRECT", "BUMP", "X17SEL"};
const char* kQtyN[4] = {"peak", "bump", "direct", "conv"};
const char* kMethN[11] = {"", "M1", "M2", "M2b", "M2d", "M3", "M3c", "B1", "B3", "B4", "D5"};

struct Rec { double Y, raw, cal, eS, eF, eT, eP, eTot, sp, chi2; int ndf; };
typedef std::tuple<int, int, int> Key;   // sel, qty, meth

void Txt(double x, double y, const char* s, double size = 0.050, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align); l.DrawLatex(x, y, s);
}

void Run() {
    Style();
    // ------------------------------------------------------------- eta results
    TFile fe(Form("%s/eta_results.root", kRes));
    TTree* t = (TTree*)fe.Get("eta");
    Int_t toy, sel, meth, qty, ndf; Double_t Y, raw, cal, eS, eF, eT, eP, eTot, sp, c2;
    t->SetBranchAddress("toy", &toy); t->SetBranchAddress("sel", &sel); t->SetBranchAddress("meth", &meth); t->SetBranchAddress("qty", &qty);
    t->SetBranchAddress("Y", &Y); t->SetBranchAddress("raw", &raw); t->SetBranchAddress("cal", &cal); t->SetBranchAddress("eStat", &eS);
    t->SetBranchAddress("eFit", &eF); t->SetBranchAddress("eTmpl", &eT); t->SetBranchAddress("eParent", &eP); t->SetBranchAddress("eTot", &eTot);
    t->SetBranchAddress("spread", &sp); t->SetBranchAddress("chi2", &c2); t->SetBranchAddress("ndf", &ndf);
    std::map<Key, std::vector<Rec>> M;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); M[Key(sel, qty, meth)].push_back(Rec{Y, raw, cal, eS, eF, eT, eP, eTot, sp, c2, ndf}); }
    // half-A expectations (reference, NOT truth)
    TFile ftp(kTemplFile);
    double ref[4][4], refE[4][4];
    for (int s = 0; s < 4; ++s) {
        M4Templ tt; LoadM4(ftp, s, tt);
        ref[s][0] = 0.5 * (tt.regPeak[kDirect5] + tt.regPeak[kConv4]); refE[s][0] = 0.5 * std::hypot(tt.regPeakE[kDirect5], tt.regPeakE[kConv4]);
        ref[s][1] = 0.5 * (tt.regBump[kBump2] + tt.regBump[kBump3]); refE[s][1] = 0.5 * std::hypot(tt.regBumpE[kBump2], tt.regBumpE[kBump3]);
        ref[s][2] = 0.5 * tt.regPeak[kDirect5]; refE[s][2] = 0.5 * tt.regPeakE[kDirect5];
        ref[s][3] = 0.5 * tt.regPeak[kConv4]; refE[s][3] = 0.5 * tt.regPeakE[kConv4];
    }
    FILE* fo = fopen(Form("%s/s5_blind_ana_eta_summary.txt", kTab), "w");
    fprintf(fo, "# s5_blind_ana_eta_summary.txt -- blind analysis of the 200 eta toys (label blind_ana), written by s5_blind_ana_summary.C AFTER the frozen run\n");
    fprintf(fo, "# rows at L_toy = 110.525 nb^-1 (simulation only, no CB); procedure: tables/s5_blind_ana_procedure.txt; no sealed file read\n");
    fprintf(fo, "# mean, rms = over the 200 toys; <eStat> event-group error, <eTmpl> half-A template MC, <ePar> expected half-B parent MC noise (common offset), mean eTot = sqrt(eStat^2+eTmpl^2)\n");
    fprintf(fo, "# rms/<eStat> checks the per-toy error; <chi2/ndf> of the fit (Pearson, Poisson variance); half-A = template expectation (x0.5), a REFERENCE, not the truth\n");
    fprintf(fo, "%-7s %-6s %-4s %8s %7s %7s %7s %7s %7s %7s %6s %7s %7s %14s\n", "sel", "qty", "meth", "mean", "rms", "<eStat>", "<eFit>", "<eTmpl>", "<ePar>", "<eTot>", "rms/eS", "chi2/nd", "<sprd>", "half-A ref");
    for (auto& kv : M) {
        const auto& v = kv.second; const int n = v.size();
        double s = 0, s2 = 0, a[10] = {0};
        for (const auto& r : v) { s += r.Y; s2 += r.Y * r.Y; a[0] += r.eS; a[1] += r.eF; a[2] += r.eT; a[3] += r.eP; a[4] += r.eTot; a[5] += r.ndf > 0 ? r.chi2 / r.ndf : 0; a[6] += r.sp; }
        const double mean = s / n, rms = std::sqrt(std::max(0., s2 / n - mean * mean));
        const int ss = std::get<0>(kv.first), qq = std::get<1>(kv.first), mm = std::get<2>(kv.first);
        fprintf(fo, "%-7s %-6s %-4s %8.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %6.2f %7.2f %7.1f %7.1f+-%-5.1f\n", kSelN[ss], kQtyN[qq], kMethN[mm], mean, rms,
                a[0] / n, a[1] / n, a[2] / n, a[3] / n, a[4] / n, rms / (a[0] / n), a[5] / n, mm == 1 ? a[6] / n : 0., ref[ss][qq], refE[ss][qq]);
    }
    // method differences per toy
    fprintf(fo, "\n# per-toy differences between methods (same toy): mean +- error of the mean, rms of the difference\n");
    auto diff = [&](int s, int q, int m1, int m2) {
        const auto& a = M[Key(s, q, m1)]; const auto& b = M[Key(s, q, m2)];
        if (a.empty() || b.empty()) return;
        double d1 = 0, d2 = 0; const int n = std::min(a.size(), b.size());
        for (int i = 0; i < n; ++i) { const double d = a[i].Y - b[i].Y; d1 += d; d2 += d * d; }
        const double m = d1 / n, r = std::sqrt(std::max(0., d2 / n - m * m));
        fprintf(fo, "%-7s %-6s %-4s - %-4s : %8.1f +- %5.1f   rms %6.1f\n", kSelN[s], kQtyN[q], kMethN[m1], kMethN[m2], m, r / std::sqrt(n), r);
    };
    for (int s = 0; s < 4; ++s) {
        diff(s, 0, 1, 2); diff(s, 0, 3, 2); diff(s, 0, 4, 2); diff(s, 0, 5, 2); diff(s, 0, 6, 2);
        diff(s, 1, 7, 2); diff(s, 1, 8, 2); diff(s, 1, 9, 2); diff(s, 1, 3, 2);
    }
    diff(1, 2, 10, 4); diff(1, 3, 10, 4);
    fclose(fo);

    // ------------------------------------------------------------- figure 1: toy 0, PEAK selection, full spectrum + zoom
    for (int s : {0, 2}) {
        TH1D* h = (TH1D*)fe.Get(Form("toy0_m4_%s", kSelN[s]));
        TH1D* m1 = (TH1D*)fe.Get(Form("toy0_M1_%s", kSelN[s]));
        TH1D* m1b = (TH1D*)fe.Get(Form("toy0_M1bkg_%s", kSelN[s]));
        TH1D* cP = (TH1D*)fe.Get(Form("toy0_M2_peak45_%s", kSelN[s]));
        TH1D* cB = (TH1D*)fe.Get(Form("toy0_M2_bump23_%s", kSelN[s]));
        TH1D* cK = (TH1D*)fe.Get(Form("toy0_M2_bkgT2_%s", kSelN[s]));
        TH1D* b1 = (TH1D*)fe.Get(Form("toy0_B1_%s", kSelN[s]));
        TH1D* b1b = (TH1D*)fe.Get(Form("toy0_B1bump_%s", kSelN[s]));
        TCanvas c("c", "", 1400, 650);
        c.Divide(2, 1, 0.002, 0.002);
        // left: full spectrum, log y, M2 components
        c.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.16);
        TH1D* fr = (TH1D*)h->Clone("fr"); fr->Reset(); fr->UseCurrentStyle();
        fr->GetXaxis()->SetRangeUser(0.35, 1.0); fr->SetMinimum(0.5); fr->SetMaximum(30 * h->GetMaximum());
        fr->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); fr->GetYaxis()->SetTitle("rows / 2 MeV/c^{2}");
        fr->Draw("AXIS");
        TH1D* sumc = (TH1D*)cP->Clone("sumc"); sumc->Add(cB); sumc->Add(cK);
        cK->SetLineColor(kAzure + 2); cK->SetLineStyle(2); cB->SetLineColor(kOrange + 1); cP->SetLineColor(kRed + 1); sumc->SetLineColor(kGray + 2);
        for (TH1D* q : {cK, cB, cP, sumc}) { q->SetLineWidth(3); q->GetXaxis()->SetRangeUser(0.400, 0.700); q->Draw("HIST SAME ][]"); }
        h->SetMarkerStyle(20); h->SetMarkerSize(0.7); h->SetLineColor(kBlack); h->Draw("E0 SAME");
        TLegend lg(0.62, 0.16, 0.97, 0.50); lg.SetTextSize(0.042);
        lg.AddEntry(h, Form("toy 0, %s", kSelN[s]), "pe"); lg.AddEntry(sumc, "M2 sum", "l");
        lg.AddEntry(cP, "M2 #eta peak (4+5)", "l"); lg.AddEntry(cB, "M2 bump (2+3)", "l"); lg.AddEntry(cK, "M2 phys. bkg", "l");
        lg.Draw();
        Txt(0.19, 0.86, "simulation (pseudo-data), L = 110.5 nb^{-1}", 0.045);
        // right: zoom with M1 and B1
        c.cd(2); gPad->SetLeftMargin(0.16);
        TH1D* fz = (TH1D*)h->Clone("fz"); fz->Reset(); fz->UseCurrentStyle();
        fz->GetXaxis()->SetRangeUser(0.40, 0.62); double mx = 0; for (int b = h->FindBin(0.40); b <= h->FindBin(0.62); ++b) mx = std::max(mx, h->GetBinContent(b));
        fz->SetMinimum(0); fz->SetMaximum(1.75 * mx);
        fz->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); fz->GetYaxis()->SetTitle("rows / 2 MeV/c^{2}");
        fz->GetXaxis()->SetNdivisions(505);
        fz->Draw("AXIS");
        TBox bw(kBumpLo, 0, kBumpHi, 1.12 * mx); bw.SetFillColorAlpha(kOrange, 0.15); bw.Draw();
        TBox pw(kPeakLo, 0, kPeakHi, 1.12 * mx); pw.SetFillColorAlpha(kRed, 0.15); pw.Draw();
        b1->SetLineColor(kGray + 2); b1->SetLineWidth(3); b1->Draw("HIST SAME ][]");
        b1b->SetLineColor(kOrange + 1); b1b->SetLineWidth(3); b1b->SetLineStyle(2); b1b->Draw("HIST SAME ][]");
        m1->SetLineColor(kRed + 1); m1->SetLineWidth(3); m1->Draw("HIST SAME ][]");
        m1b->SetLineColor(kAzure + 2); m1b->SetLineWidth(3); m1b->SetLineStyle(2); m1b->Draw("HIST SAME ][]");
        h->Draw("E0 SAME");
        fz->Draw("AXIS SAME");
        TLegend lz(0.18, 0.70, 0.66, 0.93); lz.SetTextSize(0.040);
        lz.AddEntry(m1, "M1: Gaussian + pol2", "l"); lz.AddEntry(m1b, "M1: pol2 part", "l");
        lz.AddEntry(b1, "B1: bump + #eta + pol2", "l"); lz.AddEntry(b1b, "B1: bump part", "l");
        lz.Draw();
        TLegend lw(0.67, 0.76, 0.97, 0.93); lw.SetTextSize(0.040);
        lw.AddEntry(&bw, "bump window", "f"); lw.AddEntry(&pw, "peak window", "f");
        lw.Draw();
        SaveFig(&c, Form("%s/s5_blind_ana_toy0_%s", kFig, kSelN[s]));
    }

    // ------------------------------------------------------------- figure 2: yields per method (mean +- rms over toys, mean error)
    {
        struct Item { int s, q, m; };
        std::vector<Item> peakItems = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}, {0, 0, 4}, {0, 0, 5}, {0, 0, 6}, {1, 0, 1}, {1, 0, 2}, {1, 0, 5}, {2, 0, 1}, {2, 0, 2}, {2, 0, 5},
                                       {1, 2, 4}, {1, 2, 10}};
        std::vector<Item> bumpItems = {{0, 1, 2}, {0, 1, 3}, {0, 1, 7}, {0, 1, 8}, {0, 1, 9}, {2, 1, 2}, {2, 1, 3}, {2, 1, 7}, {2, 1, 8}, {2, 1, 9}};
        for (int which = 0; which < 2; ++which) {
            const auto& items = which == 0 ? peakItems : bumpItems;
            const int n = items.size();
            TCanvas c("c2", "", 1400, 850);
            gPad->SetBottomMargin(0.30); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.17);
            TH1D fr("fr2", "", n, 0, n);
            double ymin = 1e9, ymax = -1e9;
            TGraphErrors gR(n), gE(n), gRef(n);
            for (int i = 0; i < n; ++i) {
                const auto& v = M[Key(items[i].s, items[i].q, items[i].m)];
                double s = 0, s2 = 0, e = 0; for (const auto& r : v) { s += r.Y; s2 += r.Y * r.Y; e += r.eTot; }
                const double mean = s / v.size(), rms = std::sqrt(std::max(0., s2 / v.size() - mean * mean));
                gR.SetPoint(i, i + 0.4, mean); gR.SetPointError(i, 0, rms);
                gE.SetPoint(i, i + 0.6, mean); gE.SetPointError(i, 0, e / v.size());
                gRef.SetPoint(i, i + 0.5, ref[items[i].s][items[i].q]); gRef.SetPointError(i, 0.45, refE[items[i].s][items[i].q]);
                ymin = std::min(ymin, mean - 1.3 * rms); ymax = std::max(ymax, mean + 1.6 * rms);
                fr.GetXaxis()->SetBinLabel(i + 1, Form("#splitline{%s %s}{%s}", kMethN[items[i].m], kQtyN[items[i].q], kSelN[items[i].s]));
            }
            fr.SetMinimum(std::min(0., ymin)); fr.SetMaximum(ymax);
            fr.GetXaxis()->SetLabelSize(0.055); fr.GetXaxis()->LabelsOption("v");
            fr.GetYaxis()->SetTitle(which == 0 ? "yield in the peak window [rows]" : "yield in the bump window [rows]");
            fr.GetYaxis()->SetTitleOffset(0.95);
            fr.Draw("AXIS");
            gRef.SetFillColorAlpha(kGray, 0.5); gRef.SetLineColor(kGray + 2); gRef.Draw("2 SAME");
            gR.SetMarkerStyle(20); gR.SetMarkerSize(1.4); gR.SetMarkerColor(which == 0 ? kRed + 1 : kOrange + 1); gR.SetLineColor(which == 0 ? kRed + 1 : kOrange + 1); gR.SetLineWidth(3); gR.Draw("P SAME");
            gE.SetMarkerStyle(25); gE.SetMarkerSize(1.4); gE.SetLineWidth(3); gE.SetLineColor(kAzure + 2); gE.SetMarkerColor(kAzure + 2); gE.Draw("P SAME");
            TLegend lg(0.12, 0.835, 0.98, 0.995); lg.SetNColumns(2); lg.SetTextSize(0.045);
            lg.AddEntry(&gR, "mean #pm rms over 200 toys", "pe"); lg.AddEntry(&gE, "mean #pm mean per-toy error", "pe"); lg.AddEntry(&gRef, "half-A template (reference, not truth)", "f");
            lg.Draw();
            SaveFig(&c, Form("%s/s5_blind_ana_%s_methods", kFig, which == 0 ? "peak" : "bump"));
        }
    }
    // ------------------------------------------------------------- figure 3: per-toy scatter M1 vs M2 (peak) and B4 vs M2 (bump), PEAK selection
    {
        TCanvas c("c3", "", 1300, 650); c.Divide(2, 1, 0.002, 0.002);
        const int pairs[2][3] = {{0, 1, 2}, {1, 9, 2}};   // qty, method x, method y=M2
        TGraph* g[2];
        for (int k = 0; k < 2; ++k) {
            c.cd(k + 1); gPad->SetLeftMargin(0.20);
            const auto& a = M[Key(0, pairs[k][0], pairs[k][1])]; const auto& b = M[Key(0, pairs[k][0], pairs[k][2])];
            g[k] = new TGraph(a.size());
            double lo = 1e9, hi = -1e9;
            for (size_t i = 0; i < a.size(); ++i) { g[k]->SetPoint(i, a[i].Y, b[i].Y); lo = std::min({lo, a[i].Y, b[i].Y}); hi = std::max({hi, a[i].Y, b[i].Y}); }
            TH2D fr(Form("fr3%d", k), "", 10, lo - 10, hi + 10, 10, lo - 10, hi + 10);
            fr.GetXaxis()->SetTitle(Form("%s %s yield [rows]", kMethN[pairs[k][1]], kQtyN[pairs[k][0]]));
            fr.GetYaxis()->SetTitle(Form("M2 %s yield [rows]", kQtyN[pairs[k][0]]));
            fr.GetXaxis()->SetNdivisions(505); fr.GetYaxis()->SetNdivisions(505); fr.GetYaxis()->SetTitleOffset(1.55);
            fr.DrawCopy("AXIS");
            TLine l(lo - 10, lo - 10, hi + 10, hi + 10); l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.DrawClone();
            g[k]->SetMarkerStyle(20); g[k]->SetMarkerSize(0.8); g[k]->SetMarkerColor(k == 0 ? kRed + 1 : kOrange + 1); g[k]->Draw("P SAME");
            Txt(0.20, 0.86, Form("PEAK selection, 200 toys, r = %.2f", g[k]->GetCorrelationFactor()), 0.045);
        }
        SaveFig(&c, Form("%s/s5_blind_ana_scatter", kFig));
    }

    // ------------------------------------------------------------- X17
    TFile fx(Form("%s/x17_results.root", kRes));
    TTree* x = (TTree*)fx.Get("x17");
    const char* xv[] = {"BR", "eBRstat", "eBRtot", "S", "S_sideband", "UL90_CLs", "UL90_CLs_stat", "UL90_PL", "aD", "aO", "eBRparent", "Nwin",
                        "X2_BR", "X2_eBR", "X2_UL90_CLs", "X3_BR", "X3_eBRstat", "X3_S", "X3_eBRtot", "X3_UL90_CLs",
                        "X4_BR", "X4_eBRstat", "X4_S", "X4_eBRtot", "X4_UL90_CLs", "X4_aD", "X5_BR", "X5_eBRstat", "X5_eBRtot", "X5_UL90_CLs"};
    FILE* fx2 = fopen(Form("%s/s5_blind_ana_x17_summary.txt", kTab), "w");
    fprintf(fx2, "# s5_blind_ana_x17_summary.txt -- blind analysis of the 1500 X17 toys (label blind_ana), written AFTER the frozen run\n");
    fprintf(fx2, "# L_data = 6080.5 nb^-1 (UNCONFIRMED); X1 = frozen primary; X2, X3 frozen variants; X4, X5 post-freeze variants (declared in the procedure file)\n");
    fprintf(fx2, "# distribution over the 1500 toys (the injected BR differs between toys and is unknown here)\n");
    fprintf(fx2, "%-14s %11s %11s %11s %11s %11s %11s %11s\n", "quantity", "mean", "rms", "min", "q16", "median", "q84", "max");
    std::map<std::string, std::vector<double>> X;
    for (const char* v : xv) {
        double val; x->SetBranchAddress(v, &val);
        std::vector<double> all;
        for (Long64_t i = 0; i < x->GetEntries(); ++i) { x->GetEntry(i); all.push_back(val); }
        x->ResetBranchAddresses();
        X[v] = all;
        std::vector<double> so = all; std::sort(so.begin(), so.end());
        double s = 0, s2 = 0; for (double a : all) { s += a; s2 += a * a; }
        const int n = all.size();
        fprintf(fx2, "%-14s %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n", v, s / n, std::sqrt(std::max(0., s2 / n - s * s / n / n)), so[0], so[n * 16 / 100], so[n / 2], so[n * 84 / 100], so[n - 1]);
    }
    // correlations between estimators
    auto corr = [&](const char* a, const char* b) {
        const auto& u = X[a]; const auto& w = X[b]; const int n = u.size();
        double su = 0, sw = 0, suu = 0, sww = 0, suw = 0;
        for (int i = 0; i < n; ++i) { su += u[i]; sw += w[i]; suu += u[i] * u[i]; sww += w[i] * w[i]; suw += u[i] * w[i]; }
        const double cu = suu / n - su * su / n / n, cw = sww / n - sw * sw / n / n, cuw = suw / n - su * sw / n / n;
        double d1 = 0, d2 = 0; for (int i = 0; i < n; ++i) { d1 += u[i] - w[i]; d2 += (u[i] - w[i]) * (u[i] - w[i]); }
        fprintf(fx2, "corr(%s, %s) = %.3f; slope d%s/d%s = %.3f; mean diff %.3e, rms diff %.3e\n", a, b, cuw / std::sqrt(cu * cw), b, a, cuw / cu, d1 / n, std::sqrt(d2 / n - d1 * d1 / n / n));
    };
    corr("X5_BR", "BR"); corr("X5_BR", "X3_BR"); corr("X5_BR", "X4_BR"); corr("X5_BR", "X2_BR");
    fclose(fx2);

    // figure 4: X17 example fit (toy index 0)
    {
        TH1D* hd = (TH1D*)fx.Get("ex0_data"); TH1D* hm = (TH1D*)fx.Get("ex0_model"); TH1D* hb = (TH1D*)fx.Get("ex0_bkg"); TH1D* hxx = (TH1D*)fx.Get("ex0_x17");
        TH1D* tX = (TH1D*)fx.Get("tpl_X17_perBR");
        TCanvas c("c4", "", 1000, 750);
        gPad->SetLeftMargin(0.15);
        TH1D* fr = (TH1D*)hd->Clone("fr4"); fr->Reset(); fr->UseCurrentStyle(); fr->GetXaxis()->SetRangeUser(5, 45);
        double mx = 0; for (int b = hd->FindBin(5.1); b <= hd->FindBin(44.9); ++b) mx = std::max(mx, hd->GetBinContent(b));
        fr->SetMinimum(0); fr->SetMaximum(1.3 * mx);
        fr->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]"); fr->GetYaxis()->SetTitle("rows / 0.5 MeV/c^{2}");
        fr->Draw("AXIS");
        TBox w(13.0, 0, 26.5, 1.3 * mx); w.SetFillColorAlpha(kOrange, 0.12); w.Draw();
        hd->SetMarkerStyle(20); hd->SetMarkerSize(0.8); hd->Draw("E0 SAME");
        hb->SetLineColor(kAzure + 2); hb->SetLineWidth(3); hb->SetLineStyle(2); hb->Draw("HIST SAME ][]");
        hm->SetLineColor(kGray + 2); hm->SetLineWidth(3); hm->Draw("HIST SAME ][]");
        TH1D* ref = (TH1D*)tX->Clone("ref1e5"); ref->Scale(2e-5); ref->SetLineColor(kRed + 1); ref->SetLineWidth(3); ref->GetXaxis()->SetRangeUser(5, 45); ref->Draw("HIST SAME ][]");
        TLegend lg(0.47, 0.66, 0.95, 0.92); lg.SetTextSize(0.040);
        lg.AddEntry(hd, "X17 toy (index 0)", "pe"); lg.AddEntry(hm, "X1 fit in the window", "l"); lg.AddEntry(hb, "X1 background (direct + other)", "l");
        lg.AddEntry(ref, "X17 line for BR = 2#times10^{-5}", "l");
        lg.Draw();
        Txt(0.52, 0.58, "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)", 0.040);
        Txt(0.52, 0.53, "simulation (pseudo-data)", 0.040);
        SaveFig(&c, Form("%s/s5_blind_ana_x17_example", kFig));
    }
    // figure 5: side band of three toys vs half-A template: the frozen MC pattern
    {
        TFile fp("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/pseudo/x17_toys.root");
        TTree* tt = (TTree*)fp.Get("toys"); Int_t id; tt->SetBranchAddress("id", &id);
        TCanvas c("c5", "", 1100, 700); gPad->SetLeftMargin(0.14);
        TH1D* oraw = (TH1D*)fx.Get("tpl_other_raw"); TH1D* dd = (TH1D*)fx.Get("tpl_direct5");
        TH1D* tsum = (TH1D*)oraw->Clone("tsum"); tsum->Add(dd);
        TH1D* fr = (TH1D*)tsum->Clone("fr5"); fr->Reset(); fr->UseCurrentStyle(); fr->GetXaxis()->SetRangeUser(26.5, 45); fr->SetMinimum(0); fr->SetMaximum(650);
        fr->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]"); fr->GetYaxis()->SetTitle("rows / 0.5 MeV/c^{2}");
        fr->Draw("AXIS");
        tsum->SetLineColor(kGray + 1); tsum->SetFillColorAlpha(kGray, 0.35); tsum->SetMarkerSize(0); tsum->GetXaxis()->SetRangeUser(26.5, 45); tsum->Draw("E2 SAME");
        const int col[3] = {kRed + 1, kOrange + 1, kMagenta + 1};
        TLegend lg(0.52, 0.62, 0.95, 0.92); lg.SetTextSize(0.045);
        for (int k = 0; k < 3; ++k) {
            tt->GetEntry(k);
            TH1D* h = (TH1D*)fp.Get(Form("toy_%d_mee", id)); h->SetDirectory(nullptr);
            h->SetLineColor(col[k]); h->SetMarkerColor(col[k]); h->SetMarkerStyle(20 + k); h->SetMarkerSize(0.9); h->GetXaxis()->SetRangeUser(26.5, 45);
            h->Draw("P SAME");
            lg.AddEntry(h, Form("toy index %d (side band)", k), "p");
        }
        lg.AddEntry(tsum, "half-A template, direct + other", "f");
        lg.Draw();
        Txt(0.17, 0.86, "L_{data} = 6080.5 nb^{-1}", 0.045);
        Txt(0.17, 0.80, "(UNCONFIRMED)", 0.045);
        Txt(0.17, 0.74, "side band above the window", 0.045);
        SaveFig(&c, Form("%s/s5_blind_ana_x17_sideband", kFig));
    }
    // figure 6: distributions of BR-hat (X1, X3, X4, X5) and UL90 (X1)
    {
        TCanvas c("c6", "", 1400, 650); c.Divide(2, 1, 0.002, 0.002);
        c.cd(1); gPad->SetLeftMargin(0.15);
        const char* nm[4] = {"BR", "X3_BR", "X4_BR", "X5_BR"};
        const char* lab[4] = {"X1 (frozen primary)", "X3 (frozen, 8-40 MeV/c^{2})", "X4 (post-freeze, a_{D} constrained)", "X5 (post-freeze, counting)"};
        const int col[4] = {kRed + 1, kOrange + 1, kAzure + 2, kGray + 2};
        TH1D* hh[4]; double mx = 0;
        for (int k = 0; k < 4; ++k) {
            hh[k] = new TH1D(Form("hbr%d", k), "", 70, -110, 70);
            for (double v : X[nm[k]]) hh[k]->Fill(v * 1e6);
            mx = std::max(mx, hh[k]->GetMaximum());
        }
        hh[0]->SetMaximum(1.5 * mx);
        hh[0]->GetXaxis()->SetTitle("BR estimate [10^{-6}]"); hh[0]->GetYaxis()->SetTitle("toys / 2.6#times10^{-6}");
        for (int k = 0; k < 4; ++k) { hh[k]->SetLineColor(col[k]); hh[k]->SetLineWidth(3); hh[k]->Draw(k ? "HIST SAME" : "HIST"); }
        TLegend lg(0.17, 0.66, 0.80, 0.92); lg.SetTextSize(0.040);
        for (int k = 0; k < 4; ++k) lg.AddEntry(hh[k], lab[k], "l");
        lg.Draw();
        c.cd(2); gPad->SetLeftMargin(0.18);
        TH1D hu("hu", "", 60, 0, 120), hs("hs", "", 60, 0, 120), h5("h5", "", 60, 0, 120);
        for (double v : X["UL90_CLs"]) hu.Fill(v * 1e6);
        for (double v : X["UL90_CLs_stat"]) hs.Fill(v * 1e6);
        for (double v : X["X5_UL90_CLs"]) h5.Fill(v * 1e6);
        hs.SetMaximum(1.2 * std::max({hs.GetMaximum(), hu.GetMaximum(), h5.GetMaximum()}));
        hs.GetXaxis()->SetTitle("UL_{90} on BR [10^{-6}]"); hs.GetYaxis()->SetTitle("toys / 2#times10^{-6}"); hs.GetYaxis()->SetTitleOffset(1.45);
        hs.SetLineColor(kOrange + 1); hs.SetLineWidth(3); hs.Draw("HIST");
        hu.SetLineColor(kRed + 1); hu.SetLineWidth(3); hu.Draw("HIST SAME");
        h5.SetLineColor(kGray + 2); h5.SetLineWidth(3); h5.SetLineStyle(2); h5.Draw("HIST SAME");
        TLegend l2(0.35, 0.70, 0.95, 0.92); l2.SetTextSize(0.040);
        l2.AddEntry(&hu, "X1, CLs with #sigma_{tot} (primary)", "l"); l2.AddEntry(&hs, "X1, CLs with #sigma_{stat} only", "l"); l2.AddEntry(&h5, "X5 counting (post-freeze)", "l");
        l2.Draw();
        Txt(0.40, 0.64, "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)", 0.040);
        SaveFig(&c, Form("%s/s5_blind_ana_x17_br", kFig));
    }
    printf("summary written\n");
}
}
void s5_blind_ana_summary() { bsum::Run(); }
