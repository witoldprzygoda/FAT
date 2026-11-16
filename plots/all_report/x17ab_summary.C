// ============================================================================
// x17ab_summary.C  --  WP9 X-b: summary of the X17 toy fits (x17ab_toys.C): linearity (slope), errors, coverage
// ============================================================================
// Reads /hdd3/.../x17ab/toys_cfg<k>_0.root; for the toys of the simulation study the injected BR is joined from the
// sealed key (blind/sealed/x17_key.root; that blind test was unblinded on 2026-11-13). BR estimates:
//   (1) BR1 = r_hat x 2.68e-4 (POI normalised to the fitted direct decays N5),
//   (2) BR2 = N_X / s1 = r_hat N5_hat R_X / (N5_true R_X / 2.68e-4) = r_hat x (N5_hat / N5_true) x 2.68e-4 (a_eta = 1 known in the
//       simulation; in the data s1 carries a_eta of WP5).
// Per injected BR: mean, RMS, mean error, pull, 68 % coverage of +-sigma, coverage of the 90 % CLs limit (UL >= BR), median
// limit. Linearity: weighted straight line through the five means (errors RMS/sqrt(n)).
// Output: docs/latex/tables_all/x17ab_toys.txt, figs_all/x17ab_toys_lin.pdf/.png, x17ab_toys_ul.pdf/.png
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' 'plots/all_report/x17ab_summary.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

namespace x17sum {
const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/all_report/x17ab";
const char* const kKey = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/sealed/x17_key.root";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const double kBRdir = 2.68e-4;
const double kBR[5] = {0., 1e-6, 3e-6, 1e-5, 3e-5};
struct Toy { int id; double br, rhat, rerr, ul, ulPL, N5, chi2; int ndf, ok; double kE1, kE2, a, g, thC; };
struct Stat {
    int n = 0;
    double m1 = 0, r1 = 0, m2 = 0, r2 = 0, err = 0, pm = 0, ps = 0, cov68 = 0, covUL = 0, covULc = 0, covUL2 = 0, covUL2c = 0, ulMed = 0, ul16 = 0, ul84 = 0, okf = 0, chi2 = 0, N5 = 0;
};
double Q(std::vector<double> v, double q) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const double x = q * (v.size() - 1);
    const size_t i = (size_t)x;
    return i + 1 < v.size() ? v[i] + (x - i) * (v[i + 1] - v[i]) : v.back();
}
}  // namespace x17sum

void x17ab_summary(const char* cfgList = "0,1,2,3,4,5,6,7,8,9,18,10,11,12,13,14,15,16,17") {
    using namespace x17sum;
    gStyle->SetOptStat(0);
    // key
    std::map<int, double> key;
    {
        TFile* f = TFile::Open(kKey, "READ");
        TTree* t = (TTree*)f->Get("key");
        int id; double br;
        t->SetBranchAddress("id", &id); t->SetBranchAddress("br", &br);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); key[id] = br; }
        f->Close(); delete f;
    }
    std::vector<int> cfgs;
    {
        std::string s(cfgList);
        size_t p = 0;
        while (p < s.size()) { size_t q = s.find(',', p); if (q == std::string::npos) q = s.size(); cfgs.push_back(atoi(s.substr(p, q - p).c_str())); p = q + 1; }
    }
    std::map<int, std::string> cfgText = {
        {0, "1500 sim toys + CB-like: NOMINAL fit (bkgMode 4)"}, {1, "1500 sim toys, no CB"},
        {2, "1500 sim toys, CB injected as per-row geometric mean, fit with the spread rule"},
        {3, "1500 sim toys, background pol2 (flat x (1 + b t + c t^2))"}, {4, "1500 sim toys, three background classes (bkgMode 3)"},
        {5, "1500 sim toys, background shape fixed"}, {6, "1500 sim toys, no conversion-tail nuisance"},
        {7, "1500 sim toys, X17 line from the reweighted C5 slice"}, {8, "1500 sim toys, CB input without like-sign fluctuation"},
        {9, "1500 sim toys, X17 yield nuisance 10 %"}, {18, "1500 sim toys, bkgMode 3 (no M(ee) dependence of the shape correction)"},
        {10, "smooth-parent toys + CB-like: NOMINAL fit"}, {11, "smooth-parent toys, CB geometric mean injected, spread fitted"},
        {12, "smooth-parent toys, truth and fit tilt-only (bkgMode 0)"},
        {13, "smooth-parent toys, truth with 3-class shapes, fit 2 classes bkgMode 3"},
        {14, "smooth-parent toys, truth with 3-class shapes, NOMINAL fit (2 classes, bkgMode 4)"},
        {15, "smooth-parent toys, 2-class truth, fit bkgMode 3"},
        {16, "smooth-parent toys, truth with 3-class shapes, fit 3 classes bkgMode 3"},
        {17, "smooth-parent toys, 2-class truth, fit 3 classes bkgMode 3"}};
    FILE* tb = fopen(Form("%s/x17ab_toys.txt", kTab), "w");
    fprintf(tb, "# x17ab_toys.txt -- WP9 X-b: X17 2D likelihood on toys (x17ab_toys.C, x17ab_summary.C). BR in 1e-6. BR1 = r_hat x 2.68e-4; BR2 = r_hat N5_hat/N5_true x 2.68e-4.\n");
    fprintf(tb, "# err = Hessian error of BR1; pull = (r_hat - r_true)/err; cov68 = |pull| <= 1; covUL = fraction UL(CLs 90%%, asymptotic) >= BR_true (criterion >= 0.88);\n");
    fprintf(tb, "# covULc = the same after subtracting the BR = 0 offset of BR1 from the limits (diagnostic of the parent-MC offset); ok = both fits converged.\n");
    fprintf(tb, "# covUL2 = coverage of the limit in normalisation (2), UL2 = UL x N5_hat/N5_true; covUL2c = after subtracting the BR = 0 offset of BR2.\n");
    fprintf(tb, "# N5_true = 7199.5 rows (half-A direct decays in the toy-grid M4 range, all M(ee), at L = 6080.5 nb^-1 of the simulation study).\n");
    std::map<int, std::vector<Stat>> all;
    std::map<int, double> slope1, slope1e, off1, slope2, slope2e, off2;
    for (int cfg : cfgs) {
        TFile* f = TFile::Open(Form("%s/toys_cfg%d_0.root", kDir, cfg), "READ");
        if (!f || f->IsZombie()) { printf("x17sum: cfg %d missing\n", cfg); continue; }
        TTree* t = (TTree*)f->Get("fits");
        Toy y;
        double brT = -1;
        t->SetBranchAddress("id", &y.id); t->SetBranchAddress("rhat", &y.rhat); t->SetBranchAddress("rhatErr", &y.rerr); t->SetBranchAddress("ul", &y.ul);
        t->SetBranchAddress("ulPL", &y.ulPL); t->SetBranchAddress("N5", &y.N5); t->SetBranchAddress("chi2", &y.chi2); t->SetBranchAddress("ndf", &y.ndf);
        t->SetBranchAddress("ok", &y.ok); t->SetBranchAddress("kE1", &y.kE1); t->SetBranchAddress("kE2", &y.kE2); t->SetBranchAddress("aEta", &y.a);
        t->SetBranchAddress("g", &y.g); t->SetBranchAddress("thC", &y.thC);
        if (t->GetBranch("brTrue")) t->SetBranchAddress("brTrue", &brT);
        // N5 truth: the half-A expectation (toy grid) -- printed by x17ab_toys.C; the same for every cfg
        const double N5true = 7199.5;
        std::vector<std::vector<Toy>> byBR(5);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            y.br = brT >= 0 ? brT : (key.count(y.id) ? key[y.id] : -1);
            int k = -1;
            for (int q = 0; q < 5; ++q) if (std::fabs(y.br - kBR[q]) < 1e-9) k = q;
            if (k < 0) continue;
            byBR[k].push_back(y);
        }
        f->Close(); delete f;
        std::vector<Stat> S(5);
        double off0 = 0, off20 = 0;
        for (int k = 0; k < 5; ++k) {
            Stat& s = S[k];
            std::vector<double> b1, b2, er, pl, ul, chi;
            int c68 = 0, cul = 0, okc = 0;
            for (const auto& z : byBR[k]) {
                const double rt = kBR[k] / kBRdir;
                b1.push_back(z.rhat * kBRdir);
                b2.push_back(z.rhat * z.N5 / N5true * kBRdir);
                er.push_back(z.rerr * kBRdir);
                const double p = z.rerr > 0 ? (z.rhat - rt) / z.rerr : 0;
                pl.push_back(p);
                if (std::fabs(p) <= 1) ++c68;
                ul.push_back(z.ul * kBRdir);
                if (z.ul >= rt) ++cul;
                okc += z.ok;
                if (z.ndf > 0) chi.push_back(z.chi2 / z.ndf);
            }
            s.n = (int)b1.size();
            if (!s.n) continue;
            auto mr = [&](const std::vector<double>& v, double& m, double& r) {
                m = 0; for (double x : v) m += x; m /= v.size();
                r = 0; for (double x : v) r += (x - m) * (x - m); r = std::sqrt(r / std::max<size_t>(1, v.size() - 1));
            };
            mr(b1, s.m1, s.r1); mr(b2, s.m2, s.r2);
            double tmp; mr(er, s.err, tmp); mr(pl, s.pm, s.ps);
            s.cov68 = (double)c68 / s.n; s.covUL = (double)cul / s.n; s.okf = (double)okc / s.n;
            s.ulMed = Q(ul, 0.5); s.ul16 = Q(ul, 0.16); s.ul84 = Q(ul, 0.84);
            s.chi2 = Q(chi, 0.5);
            double n5 = 0; for (const auto& z : byBR[k]) n5 += z.N5; s.N5 = n5 / s.n;
            if (k == 0) { off0 = s.m1; off20 = s.m2; }
        }
        for (int k = 0; k < 5; ++k) {
            int c = 0;
            int c2 = 0, c2c = 0;
            for (const auto& z : byBR[k]) {
                if (z.ul * kBRdir - off0 >= kBR[k]) ++c;
                const double ul2 = z.ul * z.N5 / N5true * kBRdir;   // limit in normalisation (2)
                if (ul2 >= kBR[k]) ++c2;
                if (ul2 - off20 >= kBR[k]) ++c2c;
            }
            S[k].covULc = S[k].n ? (double)c / S[k].n : 0;
            S[k].covUL2 = S[k].n ? (double)c2 / S[k].n : 0;
            S[k].covUL2c = S[k].n ? (double)c2c / S[k].n : 0;
        }
        // linearity
        auto lin = [&](int which, double& sl, double& sle, double& of) {
            TGraphErrors g;
            for (int k = 0; k < 5; ++k) {
                if (!S[k].n) continue;
                const double m = which == 1 ? S[k].m1 : S[k].m2, r = which == 1 ? S[k].r1 : S[k].r2;
                g.SetPoint(g.GetN(), kBR[k], m); g.SetPointError(g.GetN() - 1, 0, r / std::sqrt((double)S[k].n));
            }
            TF1 fl("fl", "pol1", -1e-6, 4e-5);
            g.Fit(&fl, "QN0");
            sl = fl.GetParameter(1); sle = fl.GetParError(1); of = fl.GetParameter(0);
        };
        lin(1, slope1[cfg], slope1e[cfg], off1[cfg]);
        lin(2, slope2[cfg], slope2e[cfg], off2[cfg]);
        all[cfg] = S;
        fprintf(tb, "\n## cfg %d: %s\n", cfg, cfgText.count(cfg) ? cfgText[cfg].c_str() : "");
        fprintf(tb, "BRtrue    n    <BR1>    rms1    <BR2>    rms2     err  pull_m pull_s  cov68  covUL covULc covUL2 covUL2c  ULmed  UL16  UL84    <N5>  chi2/ndf  ok\n");
        for (int k = 0; k < 5; ++k) {
            const Stat& s = S[k];
            fprintf(tb, "%5.0f %5d %8.2f %7.2f %8.2f %7.2f %7.2f %+6.2f %6.2f  %5.3f  %5.3f  %5.3f  %5.3f   %5.3f %6.2f %5.2f %5.2f %7.0f %8.2f %5.3f\n", kBR[k] * 1e6, s.n, s.m1 * 1e6,
                    s.r1 * 1e6, s.m2 * 1e6, s.r2 * 1e6, s.err * 1e6, s.pm, s.ps, s.cov68, s.covUL, s.covULc, s.covUL2, s.covUL2c, s.ulMed * 1e6, s.ul16 * 1e6, s.ul84 * 1e6, s.N5, s.chi2, s.okf);
        }
        fprintf(tb, "linearity BR1: slope %.4f +- %.4f, offset %+.2f e-6;  BR2: slope %.4f +- %.4f, offset %+.2f e-6\n", slope1[cfg], slope1e[cfg], off1[cfg] * 1e6,
                slope2[cfg], slope2e[cfg], off2[cfg] * 1e6);
        printf("x17sum: cfg %2d: slope1 %.4f +- %.4f off1 %+.2fe-6 | slope2 %.4f +- %.4f off2 %+.2fe-6 | covUL %.3f %.3f %.3f %.3f %.3f | covUL2c %.3f %.3f %.3f %.3f %.3f | ULmed(0) %.2fe-6 | N5(0) %.0f chi2 %.2f ok %.3f\n",
               cfg, slope1[cfg], slope1e[cfg], off1[cfg] * 1e6, slope2[cfg], slope2e[cfg], off2[cfg] * 1e6, S[0].covUL, S[1].covUL, S[2].covUL, S[3].covUL,
               S[4].covUL, S[0].covUL2c, S[1].covUL2c, S[2].covUL2c, S[3].covUL2c, S[4].covUL2c, S[0].ulMed * 1e6, S[0].N5, S[0].chi2, S[0].okf);
    }
    fclose(tb);
    // ---------------------------------------------------------------- figures (nominal: cfg 0 and cfg 10)
    auto drawLin = [&](int cfgA, int cfgB, const char* name) {
        if (!all.count(cfgA) || !all.count(cfgB)) return;
        TCanvas c(name, "", 1400, 600);
        c.Divide(2, 1);
        for (int pnl = 0; pnl < 2; ++pnl) {
            c.cd(pnl + 1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetGridy();
            const int cf = pnl == 0 ? cfgA : cfgB;
            const std::vector<Stat>& S = all[cf];
            TGraphErrors& g1 = *new TGraphErrors();
            TGraphErrors& g2 = *new TGraphErrors();
            for (int k = 0; k < 5; ++k) {
                g1.SetPoint(k, kBR[k] * 1e6 - 0.3, S[k].m1 * 1e6); g1.SetPointError(k, 0, S[k].r1 * 1e6);
                g2.SetPoint(k, kBR[k] * 1e6 + 0.3, S[k].m2 * 1e6); g2.SetPointError(k, 0, S[k].r2 * 1e6);
            }
            g1.SetMarkerStyle(20); g1.SetMarkerColor(kBlack); g1.SetLineColor(kBlack);
            g2.SetMarkerStyle(25); g2.SetMarkerColor(kBlue + 1); g2.SetLineColor(kBlue + 1);
            g1.SetTitle(Form("%s;injected BR_{X}BR_{ee} [10^{-6}];fitted BR_{X}BR_{ee} [10^{-6}]", pnl == 0 ? "simulation-study toys + CB-like" : "smooth-parent toys + CB-like"));
            g1.GetXaxis()->SetLimits(-2, 32);
            double lo = 1e9, hi = -1e9;
            for (int k = 0; k < 5; ++k) { lo = std::min({lo, (S[k].m1 - S[k].r1) * 1e6, (S[k].m2 - S[k].r2) * 1e6}); hi = std::max({hi, (S[k].m1 + S[k].r1) * 1e6, (S[k].m2 + S[k].r2) * 1e6}); }
            g1.SetMinimum(std::min(-5., lo - 3)); g1.SetMaximum(std::max(35., hi + 3));
            g1.GetXaxis()->SetTitleSize(0.05); g1.GetYaxis()->SetTitleSize(0.05);
            g1.Draw("AP"); g2.Draw("P SAME");
            TLine& l = *new TLine(-2, -2, 32, 32); l.SetLineStyle(2); l.SetLineColor(kRed + 1); l.Draw();
            TLegend& lg = *new TLegend(0.17, 0.70, 0.62, 0.88); lg.SetBorderSize(0); lg.SetTextSize(0.037);
            lg.AddEntry(&g1, Form("(1) r #times 2.68#times10^{-4}: slope %.3f #pm %.3f", slope1[cf], slope1e[cf]), "lp");
            lg.AddEntry(&g2, Form("(2) N_{X}/s_{1}: slope %.3f #pm %.3f", slope2[cf], slope2e[cf]), "lp");
            lg.AddEntry(&l, "fitted = injected", "l");
            lg.Draw();
        }
        c.SaveAs(Form("%s/%s.pdf", kFig, name)); c.SaveAs(Form("%s/%s.png", kFig, name));
    };
    drawLin(0, 10, "x17ab_toys_lin");
    // UL distributions at BR = 0 and pulls (cfg 0, 10)
    if (all.count(0) && all.count(10)) {
        TCanvas c("cu", "", 1400, 600);
        c.Divide(2, 1);
        const int cf[2] = {0, 10};
        const int col[2] = {kBlack, kBlue + 1};
        TH1D* hu[2]; TH1D* hp[2];
        for (int q = 0; q < 2; ++q) {
            hu[q] = new TH1D(Form("hu%d", q), ";90% CL_{s} upper limit on BR_{X}BR_{ee} at injected 0 [10^{-6}];toys", 40, 0, 40);
            hp[q] = new TH1D(Form("hp%d", q), ";pull (r_{fit} - r_{true})/#sigma_{r};toys", 40, -6, 6);
            hu[q]->SetDirectory(nullptr); hp[q]->SetDirectory(nullptr);
            TFile* f = TFile::Open(Form("%s/toys_cfg%d_0.root", kDir, cf[q]), "READ");
            TTree* t = (TTree*)f->Get("fits");
            int id; double rh, re, ul, brT = -1;
            t->SetBranchAddress("id", &id); t->SetBranchAddress("rhat", &rh); t->SetBranchAddress("rhatErr", &re); t->SetBranchAddress("ul", &ul);
            if (t->GetBranch("brTrue")) t->SetBranchAddress("brTrue", &brT);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
                t->GetEntry(i);
                const double br = brT >= 0 ? brT : (key.count(id) ? key[id] : -1);
                if (br < 0) continue;
                if (br == 0) hu[q]->Fill(ul * kBRdir * 1e6);
                if (re > 0) hp[q]->Fill((rh - br / kBRdir) / re);
            }
            f->Close(); delete f;
            hu[q]->SetLineColor(col[q]); hp[q]->SetLineColor(col[q]); hu[q]->SetLineWidth(2); hp[q]->SetLineWidth(2);
        }
        c.cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        hu[0]->SetMaximum(1.25 * std::max(hu[0]->GetMaximum(), hu[1]->GetMaximum()));
        hu[0]->GetXaxis()->SetTitleSize(0.045); hu[0]->Draw("HIST"); hu[1]->Draw("HIST SAME");
        TLegend lg(0.45, 0.72, 0.89, 0.88); lg.SetBorderSize(0); lg.SetTextSize(0.035);
        lg.AddEntry(hu[0], Form("simulation-study toys (median %.1f)", all[0][0].ulMed * 1e6), "l");
        lg.AddEntry(hu[1], Form("smooth-parent toys (median %.1f)", all[10][0].ulMed * 1e6), "l");
        lg.Draw();
        c.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        hp[0]->SetMaximum(1.25 * std::max(hp[0]->GetMaximum(), hp[1]->GetMaximum()));
        hp[0]->GetXaxis()->SetTitleSize(0.045); hp[0]->Draw("HIST"); hp[1]->Draw("HIST SAME");
        TLegend l2(0.15, 0.72, 0.60, 0.88); l2.SetBorderSize(0); l2.SetTextSize(0.035);
        l2.AddEntry(hp[0], Form("sim-study toys: mean %.2f, RMS %.2f", hp[0]->GetMean(), hp[0]->GetRMS()), "l");
        l2.AddEntry(hp[1], Form("smooth-parent: mean %.2f, RMS %.2f", hp[1]->GetMean(), hp[1]->GetRMS()), "l");
        l2.Draw();
        c.SaveAs(Form("%s/x17ab_toys_ul.pdf", kFig)); c.SaveAs(Form("%s/x17ab_toys_ul.png", kFig));
    }
    printf("x17sum: done\n");
}
