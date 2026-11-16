// ============================================================================
// x17c_summary.C  --  WP9 X-c (prefix x17c): toy summary (expected bands from 2000 background-only toys, injection
//                     tests) and the expected-limit figures. Inputs: x17c_state.root (x17c_freeze.C), toys_*.root
//                     (x17c_toys.C). Output: tables_all/x17c_expected.txt, x17c_inject.txt; figs_all/x17c_expected,
//                     x17c_toys_ul, x17c_inject (.pdf/.png).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c", kTRUE)' \
//        'plots/all_report/x17c_summary.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "x17c_lib.h"

namespace x17cs {
using namespace x17cl;
struct ToySet { double m = 0, inj = 0; std::vector<double> ul[3], rhat[3], rerr[3]; int nTot = 0, nOk[3] = {0, 0, 0}; };
bool ReadToys(double m, double inj, const std::vector<int>& firsts, ToySet& T) {
    T.m = m; T.inj = inj;
    for (int fst : firsts) {
        TFile* f = TFile::Open(Form("%s/toys_m%04d_i%04d_%04d.root", kOutDir, (int)std::lround(100 * m), (int)std::lround(inj * 1e7), fst), "READ");
        if (!f || f->IsZombie()) { printf("x17c_summary: toys m %.2f inj %.1e first %d missing\n", m, inj, fst); return false; }
        TTree* t = (TTree*)f->Get("toys");
        int ok[3]; double ul[3], rh[3], re[3];
        t->SetBranchAddress("ok", ok); t->SetBranchAddress("ul", ul); t->SetBranchAddress("rhat", rh); t->SetBranchAddress("rerr", re);
        for (Long64_t k = 0; k < t->GetEntries(); ++k) {
            t->GetEntry(k);
            ++T.nTot;
            for (int a = 0; a < 3; ++a) {
                if (!ok[a] || !std::isfinite(ul[a])) continue;
                ++T.nOk[a];
                T.ul[a].push_back(ul[a]); T.rhat[a].push_back(rh[a]); T.rerr[a].push_back(re[a]);
            }
        }
        f->Close();
    }
    return true;
}
double Quant(std::vector<double> v, double q) {
    if (v.empty()) return (double)NAN;
    std::sort(v.begin(), v.end());
    const double x = q * (v.size() - 1);
    const size_t i = (size_t)std::floor(x);
    const double w = x - i;
    return i + 1 < v.size() ? (1 - w) * v[i] + w * v[i + 1] : v.back();
}
const double kQ[5] = {0.02275, 0.15866, 0.5, 0.84134, 0.97725};
}  // namespace x17cs

void x17c_summary() {
    using namespace x17cs;
    gStyle->SetOptStat(0);
    Frozen z;
    if (!ReadConfig(z)) return;
    // Asimov expectation per mass
    TFile* fs = TFile::Open(Form("%s/x17c_state.root", kOutDir), "READ");
    TTree* te = (TTree*)fs->Get("exp");
    double m, RX, e1[5], e2[5], e3[5], nX, s1x;
    te->SetBranchAddress("m", &m); te->SetBranchAddress("RX", &RX); te->SetBranchAddress("e1", e1); te->SetBranchAddress("e2", e2); te->SetBranchAddress("e3", e3);
    te->SetBranchAddress("nX", &nX); te->SetBranchAddress("s1x", &s1x);
    std::vector<double> am, ae[3][5], anX;
    for (Long64_t k = 0; k < te->GetEntries(); ++k) {
        te->GetEntry(k);
        am.push_back(m); anX.push_back(nX);
        for (int q = 0; q < 5; ++q) { ae[0][q].push_back(e1[q]); ae[1][q].push_back(e2[q]); ae[2][q].push_back(e3[q]); }
    }
    auto asimovAt = [&](double mm, int norm, int q) -> double { for (size_t k = 0; k < am.size(); ++k) if (std::fabs(am[k] - mm) < 1e-6) return ae[norm][q][k]; return (double)NAN; };
    // toys
    const std::vector<int> f4 = {0, 500, 1000, 1500}, f1 = {0};
    std::vector<ToySet> bo;
    for (double mm : {16.7, 16.85, 17.0}) { ToySet T; if (ReadToys(mm, 0., f4, T)) bo.push_back(T); }
    for (double mm : {10.0, 25.0, 40.0}) { ToySet T; if (ReadToys(mm, 0., f1, T)) bo.push_back(T); }
    std::vector<ToySet> inj;
    for (double b : {1e-5, 3e-5, 6e-5}) { ToySet T; if (ReadToys(17.0, b, f1, T)) inj.push_back(T); }
    FILE* tb = fopen(Form("%s/x17c_expected.txt", kTab), "w");
    fprintf(tb, "# x17c_expected.txt -- expected 90 %% CLs upper limits on BR(eta -> pi+pi-X17) BR(X17 -> e+e-) of the FROZEN procedure (x17_frozen.txt).\n");
    fprintf(tb, "# Data side bands (M(ee) 2-12 and 28-62 MeV/c2) fitted with r = 0, masked rows 12-28 MeV/c2 predicted -> background-only truth.\n");
    fprintf(tb, "# Asimov: asymptotic bands of the Asimov set. Toys: quantiles (2.3, 15.9, 50, 84.1, 97.7 %%) of the observed-style UL of background-only toys\n");
    fprintf(tb, "# (over-dispersed by qScale = %.4f, constraint centres drawn), each fitted on all cells. S0, isBest >= 0, X2D; L = 6079.2 nb^-1.\n", z.qScale);
    fprintf(tb, "# norm (1) r x 2.68e-4, N5 fitted (PRIMARY); (2) eta yield from a_eta(WP5) = %.3f, sX %.3f; (3) from the peak a_P = %.3f, sX %.3f.\n", z.a0, z.xSig2, z.aP,
            z.xSig3);
    fprintf(tb, "%-6s %-4s %-7s %6s %10s %10s %10s %10s %10s   %s\n", "m", "norm", "method", "ntoys", "-2s", "-1s", "median", "+1s", "+2s", "toy/Asimov median");
    struct Q { double m; int norm; double t[5], a[5]; int n; };
    std::vector<Q> qs;
    for (const auto& T : bo) {
        for (int a = 0; a < 3; ++a) {
            Q q; q.m = T.m; q.norm = a + 1; q.n = (int)T.ul[a].size();
            for (int k = 0; k < 5; ++k) { q.t[k] = Quant(T.ul[a], kQ[k]); q.a[k] = asimovAt(T.m, a, k); }
            qs.push_back(q);
            fprintf(tb, "%6.2f %4d %-7s %6s %10.3e %10.3e %10.3e %10.3e %10.3e\n", T.m, a + 1, "Asimov", "-", q.a[0], q.a[1], q.a[2], q.a[3], q.a[4]);
            fprintf(tb, "%6.2f %4d %-7s %6d %10.3e %10.3e %10.3e %10.3e %10.3e   %.3f (converged %d/%d)\n", T.m, a + 1, "toys", q.n, q.t[0], q.t[1], q.t[2], q.t[3], q.t[4],
                    q.t[2] / q.a[2], T.nOk[a], T.nTot);
            printf("x17c_summary: m %.2f norm %d: toys median %.3e [%.3e %.3e | %.3e %.3e] (%d of %d), Asimov median %.3e [%.3e %.3e | %.3e %.3e]\n", T.m, a + 1, q.t[2],
                   q.t[0], q.t[1], q.t[3], q.t[4], q.n, T.nTot, q.a[2], q.a[0], q.a[1], q.a[3], q.a[4]);
        }
        // r_hat pulls under background only (normalisation 1)
        double s = 0, s2 = 0; int nn = 0;
        for (size_t k = 0; k < T.rhat[0].size(); ++k) if (T.rerr[0][k] > 0) { const double p = T.rhat[0][k] / T.rerr[0][k]; s += p; s2 += p * p; ++nn; }
        const double mp = s / std::max(1, nn), rp = std::sqrt(std::max(0., s2 / std::max(1, nn) - mp * mp));
        fprintf(tb, "%6.2f   background-only r_hat/sigma (norm 1): mean %+.3f, RMS %.3f (%d toys)\n", T.m, mp, rp, nn);
        printf("x17c_summary: m %.2f: r_hat pull mean %+.3f RMS %.3f\n", T.m, mp, rp);
    }
    fclose(tb);
    // injection tests (17 MeV/c2, normalisation (1) and (2)/(3) fitted on the same toys)
    tb = fopen(Form("%s/x17c_inject.txt", kTab), "w");
    fprintf(tb, "# x17c_inject.txt -- X17 injected at 17 MeV/c2 into the background-only truth of the frozen procedure (normalisation (1): rows = BR/2.68e-4 x N5 R_X),\n");
    fprintf(tb, "# 500 toys per BR; fitted with the frozen model in the three normalisations. cov = fraction of toys with UL >= injected BR (90 %% CLs: >= 0.90 expected).\n");
    fprintf(tb, "%-9s %4s %10s %10s %10s %8s %8s %8s\n", "BR_inj", "norm", "mean_rhat", "med_rhat", "RMS", "pull_m", "pull_rms", "cov");
    std::vector<double> ix, iy[3], iyE[3];
    for (const auto& T : inj) {
        for (int a = 0; a < 3; ++a) {
            const auto& r = T.rhat[a];
            double mn = 0; for (double v : r) mn += v; mn /= std::max<size_t>(1, r.size());
            double rms = 0; for (double v : r) rms += (v - mn) * (v - mn); rms = std::sqrt(rms / std::max<size_t>(1, r.size() - 1));
            double ps = 0, ps2 = 0; int np = 0, nc = 0;
            for (size_t k = 0; k < r.size(); ++k) {
                if (T.rerr[a][k] > 0) { const double p = (r[k] - T.inj) / T.rerr[a][k]; ps += p; ps2 += p * p; ++np; }
                if (T.ul[a][k] >= T.inj) ++nc;
            }
            const double pm = ps / std::max(1, np), pr = std::sqrt(std::max(0., ps2 / std::max(1, np) - pm * pm));
            const double cov = (double)nc / std::max<size_t>(1, r.size());
            fprintf(tb, "%9.2e %4d %10.3e %10.3e %10.3e %+8.3f %8.3f %8.3f\n", T.inj, a + 1, mn, Quant(r, 0.5), rms, pm, pr, cov);
            printf("x17c_summary: inj %.1e norm %d: mean rhat %.3e (ratio %.3f), RMS %.3e, pull %+.2f/%.2f, coverage %.3f\n", T.inj, a + 1, mn, mn / T.inj, rms, pm, pr, cov);
            if (a == 0) ix.push_back(T.inj);
            iy[a].push_back(mn); iyE[a].push_back(rms / std::sqrt((double)std::max<size_t>(1, r.size())));
        }
    }
    // slope through the origin and the background-only point (17 MeV, mean r_hat)
    for (int a = 0; a < 3; ++a) {
        double sxy = 0, sxx = 0, b0 = 0;
        for (const auto& T : bo) if (std::fabs(T.m - 17.0) < 1e-6) { double s = 0; for (double v : T.rhat[a]) s += v; b0 = s / std::max<size_t>(1, T.rhat[a].size()); }
        for (size_t k = 0; k < ix.size(); ++k) { sxy += ix[k] * (iy[a][k] - b0); sxx += ix[k] * ix[k]; }
        const double slope = sxx > 0 ? sxy / sxx : NAN;
        fprintf(tb, "norm %d: offset (background-only mean r_hat) %.3e, slope of (mean r_hat - offset) vs BR_inj = %.4f\n", a + 1, b0, slope);
        printf("x17c_summary: norm %d: offset %.3e, slope %.4f\n", a + 1, b0, slope);
    }
    fclose(tb);

    // ---------------------------------------------------------------- figure: expected limit vs m_X
    {
        TCanvas c("ce", "", 950, 750);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.12); gPad->SetLogy(); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04);
        const int n = (int)am.size();
        TGraphAsymmErrors b2(n), b1(n);
        TGraph med(n), med2(n), med3(n), proj;
        for (int k = 0; k < n; ++k) {
            b2.SetPoint(k, am[k], ae[0][2][k]); b2.SetPointError(k, 0, 0, ae[0][2][k] - ae[0][0][k], ae[0][4][k] - ae[0][2][k]);
            b1.SetPoint(k, am[k], ae[0][2][k]); b1.SetPointError(k, 0, 0, ae[0][2][k] - ae[0][1][k], ae[0][3][k] - ae[0][2][k]);
            med.SetPoint(k, am[k], ae[0][2][k]); med2.SetPoint(k, am[k], ae[1][2][k]); med3.SetPoint(k, am[k], ae[2][2][k]);
        }
        FILE* fp = fopen(Form("%s/x17ab_proj_mass.txt", kTab), "r");
        if (fp) {
            char line[512];
            while (fgets(line, sizeof line, fp)) { double mm, rx, a0, a1, a2; if (line[0] != '#' && sscanf(line, "%lf %lf %lf %lf %lf", &mm, &rx, &a0, &a1, &a2) == 5) proj.SetPoint(proj.GetN(), mm, a2); }
            fclose(fp);
        }
        b2.SetFillColor(kOrange); b1.SetFillColor(kGreen + 1); b2.SetLineColor(kOrange); b1.SetLineColor(kGreen + 1);
        b2.SetTitle(";m_{X} [MeV/c^{2}];90% CL_{s} UL on BR(#eta#rightarrow#pi^{+}#pi^{-}X) BR(X#rightarrowe^{+}e^{-})");
        b2.GetXaxis()->SetTitleSize(0.045); b2.GetYaxis()->SetTitleSize(0.042); b2.GetXaxis()->SetLabelSize(0.04); b2.GetYaxis()->SetLabelSize(0.04);
        b2.GetYaxis()->SetTitleOffset(1.6);
        b2.SetMinimum(3e-6); b2.SetMaximum(2e-3); b2.GetXaxis()->SetLimits(9.5, 40.5);
        b2.Draw("A3"); b1.Draw("3 SAME");
        med.SetLineStyle(2); med.SetLineWidth(3); med.Draw("L SAME");
        med2.SetLineColor(kBlue + 1); med2.SetLineWidth(2); med2.SetLineStyle(7); med2.Draw("L SAME");
        med3.SetLineColor(kMagenta + 2); med3.SetLineWidth(2); med3.SetLineStyle(3); med3.Draw("L SAME");
        if (proj.GetN() > 0) { proj.SetLineColor(kGray + 2); proj.SetLineWidth(2); proj.SetLineStyle(4); proj.Draw("L SAME"); }
        // toy quantiles (normalisation 1)
        TGraphAsymmErrors tq;
        for (const auto& q : qs) if (q.norm == 1) { const int k = tq.GetN(); tq.SetPoint(k, q.m, q.t[2]); tq.SetPointError(k, 0, 0, q.t[2] - q.t[1], q.t[3] - q.t[2]); }
        tq.SetMarkerStyle(20); tq.SetMarkerColor(kBlack); tq.SetLineColor(kBlack); tq.SetLineWidth(2); tq.Draw("P SAME");
        TLegend lg(0.30, 0.68, 0.96, 0.95); lg.SetBorderSize(0); lg.SetTextSize(0.027); lg.SetFillStyle(0); lg.SetNColumns(1);
        lg.AddEntry(&med, "(1) N5 fitted: Asimov median (PRIMARY)", "l");
        lg.AddEntry(&b1, "(1) #pm 1#sigma", "f"); lg.AddEntry(&b2, "(1) #pm 2#sigma", "f");
        lg.AddEntry(&tq, "(1) toys: median, 16-84 % (2000 / 500 per mass)", "pe");
        lg.AddEntry(&med2, Form("(2) #eta yield from a_{#eta} = %.2f #pm %.2f: median", z.a0, z.aSig), "l");
        lg.AddEntry(&med3, Form("(3) #eta yield from the peak, a_{P} = %.2f #pm %.2f: median", z.aP, z.aPSig), "l");
        if (proj.GetN() > 0) lg.AddEntry(&proj, "pre-data projection (x17ab, simulation Asimov)", "l");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
        tx.DrawLatex(0.17, 0.16, "EXPECTED (no data in M(e^{+}e^{-}) 12-28 MeV/c^{2} used); S0, isBest #geq 0, X2D");
        tx.DrawLatex(0.17, 0.205, "background from the data side bands; L = 6079.2 nb^{-1}; 2D likelihood, q~_{r}");
        c.SaveAs(Form("%s/x17c_expected.pdf", kFig)); c.SaveAs(Form("%s/x17c_expected.png", kFig));
    }
    // ---------------------------------------------------------------- figure: toy UL distribution at 17 MeV/c2
    for (const auto& T : bo) {
        if (std::fabs(T.m - 17.0) > 1e-6) continue;
        TCanvas c("cu", "", 950, 700);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogx();
        const int nb = 60;
        std::vector<double> be(nb + 1);
        for (int k = 0; k <= nb; ++k) be[k] = 5e-6 * std::pow(40., (double)k / nb);
        TH1D h1("hu1", ";90% CL_{s} UL on BR_{X} BR_{ee} (background-only toys, m_{X} = 17 MeV/c^{2});toys per bin", nb, be.data());
        TH1D h2("hu2", "", nb, be.data()), h3("hu3", "", nb, be.data());
        for (double v : T.ul[0]) h1.Fill(v);
        for (double v : T.ul[1]) h2.Fill(v);
        for (double v : T.ul[2]) h3.Fill(v);
        h1.SetLineColor(kBlack); h1.SetLineWidth(3); h2.SetLineColor(kBlue + 1); h2.SetLineWidth(2); h2.SetLineStyle(7);
        h3.SetLineColor(kMagenta + 2); h3.SetLineWidth(2); h3.SetLineStyle(3);
        h1.SetMaximum(1.35 * std::max(h1.GetMaximum(), std::max(h2.GetMaximum(), h3.GetMaximum())));
        h1.GetXaxis()->SetTitleSize(0.042); h1.GetYaxis()->SetTitleSize(0.045);
        h1.Draw("HIST"); h2.Draw("HIST SAME"); h3.Draw("HIST SAME");
        std::vector<TLine*> ls;
        for (int k = 0; k < 5; ++k) {
            const double a = asimovAt(17.0, 0, k);
            TLine* l = new TLine(a, 0, a, h1.GetMaximum() * 0.8);
            l->SetLineColor(k == 2 ? kRed + 1 : kOrange + 7); l->SetLineStyle(k == 2 ? 1 : 2); l->SetLineWidth(2); l->Draw(); ls.push_back(l);
        }
        gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.04);
        TLegend lg(0.56, 0.66, 0.95, 0.93); lg.SetBorderSize(0); lg.SetTextSize(0.03); lg.SetFillStyle(0);
        lg.AddEntry(&h1, Form("(1) N5 fitted, %zu toys", T.ul[0].size()), "l");
        lg.AddEntry(&h2, "(2) a_{#eta} normalisation", "l");
        lg.AddEntry(&h3, "(3) peak normalisation", "l");
        lg.AddEntry(ls[2], "(1) Asimov median", "l"); lg.AddEntry(ls[0], "(1) Asimov #pm1#sigma, #pm2#sigma", "l");
        lg.Draw();
        c.SaveAs(Form("%s/x17c_toys_ul.pdf", kFig)); c.SaveAs(Form("%s/x17c_toys_ul.png", kFig));
        for (auto* l : ls) delete l;
    }
    // ---------------------------------------------------------------- figure: injection
    if (!ix.empty()) {
        TCanvas c("ci", "", 900, 700);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TGraphErrors g[3];
        double b0[3] = {0, 0, 0};
        for (const auto& T : bo) if (std::fabs(T.m - 17.0) < 1e-6) for (int a = 0; a < 3; ++a) { double s = 0; for (double v : T.rhat[a]) s += v; b0[a] = s / std::max<size_t>(1, T.rhat[a].size()); }
        const int col[3] = {kBlack, kBlue + 1, kMagenta + 2};
        for (int a = 0; a < 3; ++a) {
            g[a].SetPoint(0, 0., b0[a] * 1e5);
            for (size_t k = 0; k < ix.size(); ++k) { g[a].SetPoint(k + 1, ix[k] * 1e5 * (1 + 0.01 * a), iy[a][k] * 1e5); g[a].SetPointError(k + 1, 0, iyE[a][k] * 1e5); }
            g[a].SetMarkerStyle(20 + a); g[a].SetMarkerColor(col[a]); g[a].SetLineColor(col[a]);
        }
        g[0].SetTitle(";injected BR_{X} BR_{ee} [10^{-5}];mean fitted BR_{X} BR_{ee} [10^{-5}]");
        g[0].GetXaxis()->SetTitleSize(0.045); g[0].GetYaxis()->SetTitleSize(0.045);
        g[0].GetXaxis()->SetLimits(-0.3, 6.6); g[0].SetMinimum(-0.5); g[0].SetMaximum(8);
        g[0].Draw("AP"); g[1].Draw("P SAME"); g[2].Draw("P SAME");
        TLine d(0, 0, 6.5, 6.5); d.SetLineStyle(2); d.Draw();
        TLegend lg(0.16, 0.70, 0.62, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.032);
        lg.AddEntry(&g[0], "(1) N5 fitted", "pe"); lg.AddEntry(&g[1], "(2) a_{#eta} normalisation", "pe"); lg.AddEntry(&g[2], "(3) peak normalisation", "pe");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
        tx.DrawLatex(0.30, 0.17, "m_{X} = 17 MeV/c^{2}; truth = data side-band fit (r = 0) + X17; 500 toys per point");
        tx.DrawLatex(0.30, 0.22, "dashed: fitted = injected; (2), (3): slope N5/(a N5_{ref}) expected");
        c.SaveAs(Form("%s/x17c_inject.pdf", kFig)); c.SaveAs(Form("%s/x17c_inject.png", kFig));
    }
    fs->Close();
    printf("x17c_summary: done\n");
}
