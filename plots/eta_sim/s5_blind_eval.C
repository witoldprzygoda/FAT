// s5_blind_eval.C -- unblinding of the blind test of the signal extraction (label blind_eval)
//
// SIMULATION ONLY (SMASH std + lep + eta, partition mode A, PDG weights x physics weights V3 main). No CB exists in the
// simulation and none is added. This macro does not touch any simulation row: it compares the frozen results of the blind
// analysis (blind/results/*_results.root, written by s5_blind_ana_*.C) with the sealed key (blind/sealed/*_key.root, written by
// s5_blind_prep.C) and with the half-A templates (blind/templates/templates.root, for the reference values only).
//
//   eta : per selection x method x quantity: bias against the half-B expectation mu (the parameter) and against the true rows
//         of every toy T; error of the bias from the toy spread and from the common MC terms (template, parent) quoted by the
//         analyst; pulls (mean, width) with the quoted errors; coverage of the 68 % intervals; Y against T (r, slope).
//   X17 : per estimator X1..X5 and injected BR: mean, spread, linearity, pulls, coverage of the 68 % interval, false positives
//         at BR = 0 (one-sided Z = BR-hat / sigma >= 3 primary, also 1.645, 2, 5), UL90 distribution and coverage, comparison
//         with the expected limit of the optimisation for the X17 selection of the blind test (the v1 selection;
//         tables/s5_opt_x17_result.txt, block "## v1 selection", rules v2 since 2026-11-13).
//
// Outputs: tables/s5_blind_eval_{eta,eta_ab,x17}.txt, figs/s5_blind_eval_*.{pdf,png}
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_eval_build", kTRUE)' \
//        plots/eta_sim/s5_blind_eval.C+ > /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_eval.log 2>&1
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace bev {

const char* const kBase = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/figs";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables";
const char* const kSelN[4] = {"PEAK", "DIRECT", "BUMP", "X17SEL"};
const char* const kMethN[11] = {"", "M1", "M2", "M2b", "M2d", "M3", "M3c", "B1", "B3", "B4", "D5"};
const char* const kQtyN[4] = {"peak", "bump", "direct", "conv"};
const char* const kQtyLong[4] = {"classes 4+5 in the peak window", "classes 2+3 in the bump window",
                                 "class 5 in the peak window", "class 4 in the peak window"};
const int kNToy = 200;
const double kSc = 6080.5 / 221.05;      // L_data (UNCONFIRMED) / L_sim
// expected 90 % CL limit of the optimisation for the X17 selection of the blind test (the v1 selection): median CLs and its
// 16-84 % band, read by ReadExpUL() from tables/s5_opt_x17_result.txt, block "## v1 selection" (s5_opt_x17.C). With the weight
// rules v2 (2026-11-13) 5.4805e-6, band 3.8217e-6 .. 7.9215e-6; the rules-v1 values (tables/s5_opt_x17_frozen_v1.txt) were
// 5.477e-6, 3.819e-6 .. 7.918e-6. NaN if the table cannot be read.
double kExpUL = std::numeric_limits<double>::quiet_NaN(), kExpUL16 = kExpUL, kExpUL84 = kExpUL;
bool ReadExpUL() {
    std::ifstream in(Form("%s/s5_opt_x17_result.txt", kTab));
    std::string l;
    bool blk = false;
    while (std::getline(in, l)) {
        if (l.compare(0, 3, "## ") == 0) { blk = l.compare(0, 16, "## v1 selection ") == 0; continue; }
        if (!blk || l.compare(0, 19, "expected 90% CL UL ") != 0) continue;
        const size_t a = l.find("median CLs (n = "), b = l.find("band n = ");
        if (a == std::string::npos || b == std::string::npos) break;
        const size_t a2 = l.find(") ", a), b2 = l.find(": ", b);
        double m = 0, lo = 0, hi = 0;
        if (a2 == std::string::npos || b2 == std::string::npos || sscanf(l.c_str() + a2 + 2, "%lf", &m) != 1 ||
            sscanf(l.c_str() + b2 + 2, "%lf..%lf", &lo, &hi) != 2)
            break;
        kExpUL = 1e-6 * m; kExpUL16 = 1e-6 * lo; kExpUL84 = 1e-6 * hi;
        return true;
    }
    printf("ERROR: expected X17 limit of the v1 selection not found in %s/s5_opt_x17_result.txt\n", kTab);
    return false;
}

// ----------------------------------------------------------------------------------------------------------------- helpers
double Mean(const std::vector<double>& v) { double s = 0; for (double x : v) s += x; return v.empty() ? 0 : s / v.size(); }
double Sd(const std::vector<double>& v) {
    if (v.size() < 2) return 0;
    const double m = Mean(v); double s = 0; for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (v.size() - 1));
}
double Quant(std::vector<double> v, double q) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const double p = q * (v.size() - 1); const int i = (int)std::floor(p); const double f = p - i;
    if (i + 1 >= (int)v.size()) return v.back();
    return v[i] * (1 - f) + v[i + 1] * f;
}
double Cov(const std::vector<double>& a, const std::vector<double>& b) {
    const double ma = Mean(a), mb = Mean(b); double s = 0;
    for (size_t i = 0; i < a.size(); ++i) s += (a[i] - ma) * (b[i] - mb);
    return a.size() > 1 ? s / (a.size() - 1) : 0;
}
double Frac(const std::vector<bool>& v) { int n = 0; for (bool b : v) n += b; return v.empty() ? 0 : double(n) / v.size(); }
double BinErr(double p, int n) { return n > 0 ? std::sqrt(std::max(p * (1 - p), 0.25 / n) / n) : 0; }

void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.055, "XYZ"); gStyle->SetTitleSize(0.060, "XYZ");
    gStyle->SetTitleOffset(1.05, "X"); gStyle->SetTitleOffset(1.20, "Y");
    gStyle->SetPadLeftMargin(0.16); gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.08);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0); gStyle->SetLegendTextSize(0.050);
    gStyle->SetEndErrorSize(0); gStyle->SetLineWidth(2); gStyle->SetFrameLineWidth(2); gStyle->SetHistLineWidth(3);
}
void SaveFig(TCanvas* c, const std::string& name) {
    const std::string base = std::string(kFig) + "/" + name;
    c->SaveAs((base + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", base.c_str(), base.c_str(), base.c_str()));
    c->SaveAs((base + ".png").c_str());
}
void Txt(double x, double y, const char* s, double size = 0.055, int col = kBlack, int align = 12) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.SetTextFont(42);
    t.DrawLatex(x, y, s);
}

// ----------------------------------------------------------------------------------------------------------------- eta
struct Rec { double Y = 0, eS = 0, eT = 0, eP = 0, eTot = 0; bool ok = false; };
struct EtaRes {
    int sel = 0, meth = 0, qty = 0, n = 0;
    double mu = 0, meanT = 0, refA = 0, refAe = 0;
    double meanY = 0, sdY = 0, bias = 0, eBiasStat = 0, common = 0, eBiasAll = 0, zBias = 0, rel = 0, eRelStat = 0, eRelAll = 0;
    double biasT = 0, eBiasT = 0, sdYT = 0;
    double mS = 0, mTm = 0, mP = 0, mTot = 0, spreadRatio = 0;
    double pMuTot_m = 0, pMuTot_s = 0, pMuAll_m = 0, pMuAll_s = 0, pMuStat_m = 0, pMuStat_s = 0, pTTot_m = 0, pTTot_s = 0, pTStat_m = 0, pTStat_s = 0;
    double covMuTot = 0, covMuAll = 0, covMuStat = 0, covTTot = 0, covTStat = 0;
    double r = 0, slope = 0, eSlope = 0;
    std::vector<double> Y, T, pullMuTot, pullTTot;
};
int Key(int s, int m, int q) { return s * 10000 + m * 10 + q; }

// ----------------------------------------------------------------------------------------------------------------- X17
const int kNEst = 5;
const char* const kEstN[kNEst] = {"X1", "X2", "X3", "X4", "X5"};
const char* const kEstLong[kNEst] = {"X1 frozen primary: window fit, a_D and a_O free, sigma x sqrt(S)",
                                     "X2 frozen: raw half-A 'other', template variance",
                                     "X3 frozen: fit 8-40 MeV/c2, sigma x sqrt(S)",
                                     "X4 post-freeze: a_D = 1 +- 0.10",
                                     "X5 post-freeze: counting, B from half A"};
const double kBR[5] = {0, 1e-6, 3e-6, 1e-5, 3e-5};
struct XToy { int id = 0, ibr = -1; double br = 0, mus = 0, mub = 0; int ns = 0, nb = 0, nobs = 0; double Nwin = -1;
              double BR[kNEst], sT[kNEst], sS[kNEst], UL[kNEst]; double ULstat1 = 0, ULpl1 = 0; };
struct XRes {
    int n = 0;
    double mean = 0, sd = 0, med = 0, off = 0, eOff = 0;
    double pTot_m = 0, pTot_s = 0, pStat_m = 0, pStat_s = 0, cov68 = 0, mST = 0, mSS = 0;
    double ulMed = 0, ul16 = 0, ul84 = 0, ul025 = 0, ul975 = 0, ulCov = 0, ulCovE = 0;
    double fz1645 = 0, fz2 = 0, fz3 = 0, fz5 = 0, fz3stat = 0, fzm3 = 0;
};

} // namespace bev

using namespace bev;

void s5_blind_eval() {
    Style();
    gSystem->mkdir(kFig, kTRUE); gSystem->mkdir(kTab, kTRUE);
    if (ReadExpUL()) printf("expected UL90 (v1 selection, s5_opt_x17_result.txt): median %.4e, band %.4e .. %.4e\n", kExpUL, kExpUL16, kExpUL84);

    // =================================================================================================== read eta results
    std::map<int, std::vector<Rec>> R;
    std::map<int, double> calSum;   // sum over the toys of the calibration Y / raw (second realisation below)
    std::map<int, int> calN;
    {
        TFile f(Form("%s/results/eta_results.root", kBase));
        TTree* t = (TTree*)f.Get("eta");
        Int_t toy, sel, meth, qty; Double_t Y, eS, eT, eP, eTot, raw;
        t->SetBranchAddress("toy", &toy); t->SetBranchAddress("sel", &sel); t->SetBranchAddress("meth", &meth); t->SetBranchAddress("qty", &qty);
        t->SetBranchAddress("Y", &Y); t->SetBranchAddress("eStat", &eS); t->SetBranchAddress("eTmpl", &eT); t->SetBranchAddress("eParent", &eP);
        t->SetBranchAddress("eTot", &eTot); t->SetBranchAddress("raw", &raw);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            auto& v = R[Key(sel, meth, qty)]; if (v.empty()) v.resize(kNToy);
            if (toy < 0 || toy >= kNToy) { printf("ERROR toy index %d\n", toy); continue; }
            v[toy] = Rec{Y, eS, eT, eP, eTot, true};
            if (raw != 0) { calSum[Key(sel, meth, qty)] += Y / raw; calN[Key(sel, meth, qty)]++; }
        }
        printf("eta results: %lld entries, %zu (sel, meth, qty) combinations\n", t->GetEntries(), R.size());
    }
    // =================================================================================================== read the sealed key
    // n[toy][sel][reg][grp]
    static int N[kNToy][4][8][9];
    double MU[4][8][9], MUT[4][8];
    {
        TFile f(Form("%s/sealed/eta_key.root", kBase));
        TTree* k = (TTree*)f.Get("key");
        Int_t toy, sel, reg, n[9], ntot;
        k->SetBranchAddress("toy", &toy); k->SetBranchAddress("sel", &sel); k->SetBranchAddress("reg", &reg);
        k->SetBranchAddress("n", n); k->SetBranchAddress("ntot", &ntot);
        int cnt = 0;
        for (Long64_t i = 0; i < k->GetEntries(); ++i) { k->GetEntry(i); for (int g = 0; g < 9; ++g) N[toy][sel][reg][g] = n[g]; ++cnt; }
        TTree* e = (TTree*)f.Get("expect");
        Double_t mu[9], mutot, vcoh, vind;
        e->SetBranchAddress("sel", &sel); e->SetBranchAddress("reg", &reg); e->SetBranchAddress("mu", mu); e->SetBranchAddress("mutot", &mutot);
        e->SetBranchAddress("vcoh", &vcoh); e->SetBranchAddress("vind", &vind);
        for (Long64_t i = 0; i < e->GetEntries(); ++i) { e->GetEntry(i); for (int g = 0; g < 9; ++g) MU[sel][reg][g] = mu[g]; MUT[sel][reg] = mutot; }
        printf("eta key: %d entries; expect: %lld entries\n", cnt, e->GetEntries());
    }
    auto Truth = [&](int toy, int s, int q) -> double {
        const int* n = N[toy][s][q == 1 ? 1 : 3];
        if (q == 0) return n[2] + n[3];
        if (q == 1) return n[0] + n[1];
        if (q == 2) return n[3];
        return n[2];
    };
    auto Expect = [&](int s, int q) -> double {
        const double* m = MU[s][q == 1 ? 1 : 3];
        if (q == 0) return m[2] + m[3];
        if (q == 1) return m[0] + m[1];
        if (q == 2) return m[3];
        return m[2];
    };
    // half-A reference (templates x 0.5 = rows at L_toy)
    double REF[4][4], REFE[4][4];
    double REFBKG[4][2], REFBKGE[4][2];   // half-A bkgT2 in the peak (0) and bump (1) window
    {
        TFile f(Form("%s/templates/templates.root", kBase));
        const char* comp[4] = {"peak45", "bump23", "direct5", "conv4"};
        for (int s = 0; s < 4; ++s) {
            for (int q = 0; q < 4; ++q) {
                TH1D* h = (TH1D*)f.Get(Form("reg_%s_%s", kSelN[s], comp[q]));
                const int b = q == 1 ? 2 : 4;
                REF[s][q] = h ? 0.5 * h->GetBinContent(b) : -1; REFE[s][q] = h ? 0.5 * h->GetBinError(b) : 0;
            }
            TH1D* h = (TH1D*)f.Get(Form("reg_%s_bkgT2", kSelN[s]));
            REFBKG[s][0] = 0.5 * h->GetBinContent(4); REFBKGE[s][0] = 0.5 * h->GetBinError(4);
            REFBKG[s][1] = 0.5 * h->GetBinContent(2); REFBKGE[s][1] = 0.5 * h->GetBinError(2);
        }
    }

    // =================================================================================================== eta statistics
    std::map<int, EtaRes> E;
    for (auto& kv : R) {
        const int key = kv.first, s = key / 10000, m = (key / 10) % 1000, q = key % 10;
        EtaRes r; r.sel = s; r.meth = m; r.qty = q; r.mu = Expect(s, q); r.refA = REF[s][q]; r.refAe = REFE[s][q];
        std::vector<double> Y, T, YT, eS, eTm, eP, eTot;
        std::vector<double> pMuTot, pMuAll, pMuStat, pTTot, pTStat;
        std::vector<bool> cMuTot, cMuAll, cMuStat, cTTot, cTStat;
        for (int toy = 0; toy < kNToy; ++toy) {
            const Rec& x = kv.second[toy];
            if (!x.ok || !std::isfinite(x.Y) || x.eTot <= 0) continue;
            const double t = Truth(toy, s, q);
            const double eAll = std::sqrt(x.eTot * x.eTot + x.eP * x.eP);
            Y.push_back(x.Y); T.push_back(t); YT.push_back(x.Y - t);
            eS.push_back(x.eS); eTm.push_back(x.eT); eP.push_back(x.eP); eTot.push_back(x.eTot);
            pMuTot.push_back((x.Y - r.mu) / x.eTot); pMuAll.push_back((x.Y - r.mu) / eAll); pMuStat.push_back((x.Y - r.mu) / x.eS);
            pTTot.push_back((x.Y - t) / x.eTot); pTStat.push_back((x.Y - t) / x.eS);
            cMuTot.push_back(std::fabs(x.Y - r.mu) <= x.eTot); cMuAll.push_back(std::fabs(x.Y - r.mu) <= eAll);
            cMuStat.push_back(std::fabs(x.Y - r.mu) <= x.eS); cTTot.push_back(std::fabs(x.Y - t) <= x.eTot);
            cTStat.push_back(std::fabs(x.Y - t) <= x.eS);
        }
        r.n = Y.size();
        r.meanY = Mean(Y); r.sdY = Sd(Y); r.meanT = Mean(T);
        r.bias = r.meanY - r.mu; r.eBiasStat = r.sdY / std::sqrt((double)r.n);
        r.mS = Mean(eS); r.mTm = Mean(eTm); r.mP = Mean(eP); r.mTot = Mean(eTot);
        r.common = std::sqrt(r.mTm * r.mTm + r.mP * r.mP);
        r.eBiasAll = std::sqrt(r.eBiasStat * r.eBiasStat + r.common * r.common);
        r.zBias = r.bias / r.eBiasAll;
        r.rel = r.bias / r.mu; r.eRelStat = r.eBiasStat / r.mu; r.eRelAll = r.eBiasAll / r.mu;
        r.biasT = Mean(YT); r.sdYT = Sd(YT); r.eBiasT = r.sdYT / std::sqrt((double)r.n);
        r.spreadRatio = r.sdY / r.mS;
        r.pMuTot_m = Mean(pMuTot); r.pMuTot_s = Sd(pMuTot); r.pMuAll_m = Mean(pMuAll); r.pMuAll_s = Sd(pMuAll);
        r.pMuStat_m = Mean(pMuStat); r.pMuStat_s = Sd(pMuStat); r.pTTot_m = Mean(pTTot); r.pTTot_s = Sd(pTTot);
        r.pTStat_m = Mean(pTStat); r.pTStat_s = Sd(pTStat);
        r.covMuTot = Frac(cMuTot); r.covMuAll = Frac(cMuAll); r.covMuStat = Frac(cMuStat); r.covTTot = Frac(cTTot); r.covTStat = Frac(cTStat);
        const double vT = Cov(T, T), vY = Cov(Y, Y), cYT = Cov(Y, T);
        r.r = (vT > 0 && vY > 0) ? cYT / std::sqrt(vT * vY) : 0;
        r.slope = vT > 0 ? cYT / vT : 0;
        r.eSlope = (vT > 0 && r.n > 2) ? std::sqrt(std::max(0.0, (vY / vT - r.slope * r.slope) / (r.n - 2))) : 0;
        r.Y = Y; r.T = T; r.pullMuTot = pMuTot; r.pullTTot = pTTot;
        E[key] = r;
    }

    // ------------------------------------------------------------------------------------------------ eta table
    {
        FILE* fo = fopen(Form("%s/s5_blind_eval_eta.txt", kTab), "w");
        fprintf(fo, "# s5_blind_eval_eta.txt -- unblinding of the eta toys (label blind_eval). SIMULATION ONLY, no CB.\n");
        fprintf(fo, "# 200 toys at L_toy = L_sim/2 = 110.525 nb^-1 (Poisson bootstrap of half B). Y = blind estimate (blind/results/eta_results.root),\n");
        fprintf(fo, "# T = true rows of the toy (sealed key), mu = half-B expectation per toy (sealed key; the parameter), refA = half-A template x 0.5.\n");
        fprintf(fo, "# bias = <Y> - mu +- eStat (toy rms/sqrt(n)) +- common (sqrt(<eTmpl>^2 + <eParent>^2), the MC terms common to all toys, as quoted by the analyst);\n");
        fprintf(fo, "# z = bias / sqrt(eStat^2 + common^2); rel = bias / mu. biasT = <Y - T> +- rms/sqrt(n).\n");
        fprintf(fo, "# pulls: (Y - mu)/eTot ['tot' = the quoted error sqrt(eStat^2+eTmpl^2)], (Y - mu)/eAll [eAll = sqrt(eTot^2 + eParent^2)],\n");
        fprintf(fo, "#        (Y - mu)/eStat, (Y - T)/eTot, (Y - T)/eStat: mean and rms. cov = fraction of toys inside +-1 error (nominal 0.683; binomial\n");
        fprintf(fo, "#        error 0.033, but all toys share one bias, so the coverage is conditional on this MC realisation).\n");
        fprintf(fo, "# rms/<eS> = toy spread over the mean per-toy event-group error. r, slope = correlation and slope of Y against T over the toys.\n");
        for (int q = 0; q < 4; ++q) {
            fprintf(fo, "\n## quantity %s: %s\n", kQtyN[q], kQtyLong[q]);
            fprintf(fo, "%-7s %-4s %8s %8s %8s | %8s %7s | %8s %6s %6s %6s %6s %7s | %8s %6s | %6s %6s %6s %6s | %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s | %5s %5s %5s %5s %5s | %5s %5s %5s\n",
                    "sel", "meth", "mu", "<T>", "refA", "<Y>", "rmsY", "bias", "eStat", "common", "z", "rel%", "eRelAll%",
                    "biasT", "eBT", "<eS>", "<eTm>", "<eP>", "rms/eS",
                    "pTm", "pTs", "pAm", "pAs", "pSm", "pSs", "qTm", "qTs", "qSm", "qSs",
                    "cTot", "cAll", "cStat", "cTtot", "cTst", "r", "slope", "eSl");
            for (int s = 0; s < 4; ++s)
                for (int m = 1; m <= 10; ++m) {
                    auto it = E.find(Key(s, m, q)); if (it == E.end()) continue;
                    const EtaRes& r = it->second;
                    fprintf(fo, "%-7s %-4s %8.2f %8.2f %8.2f | %8.2f %7.2f | %+8.2f %6.2f %6.2f %+6.2f %+6.1f %7.1f | %+8.2f %6.2f | %6.2f %6.2f %6.2f %6.3f | %+5.2f %5.2f %+5.2f %5.2f %+5.2f %5.2f %+5.2f %5.2f %+5.2f %5.2f | %5.3f %5.3f %5.3f %5.3f %5.3f | %5.3f %5.3f %5.3f\n",
                            kSelN[s], kMethN[m], r.mu, r.meanT, r.refA, r.meanY, r.sdY, r.bias, r.eBiasStat, r.common, r.zBias, 100 * r.rel, 100 * r.eRelAll,
                            r.biasT, r.eBiasT, r.mS, r.mTm, r.mP, r.spreadRatio,
                            r.pMuTot_m, r.pMuTot_s, r.pMuAll_m, r.pMuAll_s, r.pMuStat_m, r.pMuStat_s, r.pTTot_m, r.pTTot_s, r.pTStat_m, r.pTStat_s,
                            r.covMuTot, r.covMuAll, r.covMuStat, r.covTTot, r.covTStat, r.r, r.slope, r.eSlope);
                }
        }
        fprintf(fo, "\n# column key: pT = pull (Y-mu)/eTot, pA = (Y-mu)/eAll, pS = (Y-mu)/eStat, qT = (Y-T)/eTot, qS = (Y-T)/eStat; m = mean, s = rms.\n");
        fprintf(fo, "# cTot/cAll/cStat = coverage of mu by Y +- eTot / eAll / eStat; cTtot/cTst = coverage of T by Y +- eTot / eStat.\n");
        // second, independent realisation of the MC noise: the analyst's half-A Asimov closure (dev.log, before the toys) for the
        // methods that do not take their background from the templates. Raw values from
        // /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_ana_build/dev.log (raw = Gaussian or bifurcated-Gaussian
        // content in the window), times the mean calibration Y/raw of the toys (computed above from eta_results);
        // truth = half-A template x 0.5 (REF, from templates.root).
        // Weight rules v2 (2026-11-13): raw values of the v2 dev.log and of s5_write_m3exact.C. The rules-v1 values were M1 raw
        // 195.1 / 195.1 / 209.6 (calibration 1.01067 / 1.01298 / 1.02473, truth 220.31 / 230.23 / 241.00), M3 exact 130.09 / 129.07 /
        // 128.54, B1 raw 950.5 / 925.4 / 985.2 (calibration 1.07627 / 1.08416 / 1.08794, truth 957.93 / 1003.84 / 1064.39).
        fprintf(fo, "\n## second realisation: half-A Asimov closure (analyst dev.log) against the half-B toys (this table), relative bias in %%\n");
        fprintf(fo, "# sigma(avg) = <eParent>/mu/sqrt(2): the two realisations differ by independent parent MC noise; the average estimates the method bias.\n");
        fprintf(fo, "# M1, B1: analyst dev.log (binned); M3: exact unbinned half-A sums (s5_write_m3exact.C; dev.log binned values give -37.5/-41.2/-40.6 %%).\n");
        fprintf(fo, "# A cal = A raw x mean Y/raw of the toys (1 for M3); A truth = half-A template x 0.5.\n");
        fprintf(fo, "%-7s %-4s %9s %9s %9s | %8s %8s | %8s %8s\n", "sel", "meth", "A raw", "A cal", "A truth", "A rel%", "B rel%", "avg%", "sig(avg)");
        struct AA { int s, m, q; double raw; };
        const AA aa[] = {{0, 1, 0, 199.1}, {1, 1, 0, 199.4}, {2, 1, 0, 214.1},
                         // M3 rows: EXACT (unbinned) half-A sums from s5_write_m3exact.C (review finding L5); the binned dev.log
                         // values are 139.8 / 137.3 / 145.4 with truth 223.6 / 233.7 / 244.7 (rules v2)
                         {0, 5, 0, 134.15}, {1, 5, 0, 133.46}, {2, 5, 0, 133.19},
                         {0, 7, 1, 933.8}, {1, 7, 1, 908.3}, {2, 7, 1, 967.0}};
        for (const AA& a : aa) {
            const int key = Key(a.s, a.m, a.q);
            const EtaRes& r = E[key];
            const double cal = a.m == 5 ? 1.0 : (calN[key] > 0 ? calSum[key] / calN[key] : std::numeric_limits<double>::quiet_NaN());
            const double truth = REF[a.s][a.q];
            const double ac = a.raw * cal, ar = 100 * (ac / truth - 1), br = 100 * r.rel, sg = 100 * r.mP / r.mu / std::sqrt(2.);
            fprintf(fo, "%-7s %-4s %9.1f %9.1f %9.1f | %+8.1f %+8.1f | %+8.1f %8.1f\n", kSelN[a.s], kMethN[a.m], a.raw, ac, truth, ar, br, 0.5 * (ar + br), sg);
        }
        fclose(fo);
    }
    // ------------------------------------------------------------------------------------------------ half A vs half B table
    {
        FILE* fo = fopen(Form("%s/s5_blind_eval_eta_ab.txt", kTab), "w");
        fprintf(fo, "# s5_blind_eval_eta_ab.txt -- half A (templates, x 0.5) against half B (sealed expectation) and the toys, rows at L_toy = 110.525 nb^-1\n");
        fprintf(fo, "# SIMULATION ONLY. Targets in their windows (peak window 0.5333-0.5587, bump window 0.4576-0.5255 GeV/c2). T2 = physical background.\n");
        fprintf(fo, "# B - A / sqrt(eA^2 + eA^2): the half-B MC error is taken equal to the half-A error (same statistics).\n");
        fprintf(fo, "%-7s %-28s %10s %8s %10s %10s %8s %8s\n", "sel", "quantity", "half A", "eA", "half B mu", "toy <T>", "B-A", "pull");
        const char* nm[4] = {"peak (cl. 4+5)", "bump (cl. 2+3)", "direct (cl. 5)", "conv (cl. 4)"};
        for (int s = 0; s < 4; ++s) {
            for (int q = 0; q < 4; ++q) {
                double mt = 0; for (int t = 0; t < kNToy; ++t) mt += Truth(t, s, q); mt /= kNToy;
                const double b = Expect(s, q), a = REF[s][q], ea = REFE[s][q];
                fprintf(fo, "%-7s %-28s %10.2f %8.2f %10.2f %10.2f %+8.2f %+8.2f\n", kSelN[s], nm[q], a, ea, b, mt, b - a, ea > 0 ? (b - a) / (std::sqrt(2.) * ea) : 0.);
            }
            for (int w = 0; w < 2; ++w) {
                const int reg = w == 0 ? 3 : 1;
                double b = 0; for (int g = 4; g < 9; ++g) b += MU[s][reg][g];
                double mt = 0; for (int t = 0; t < kNToy; ++t) for (int g = 4; g < 9; ++g) mt += N[t][s][reg][g]; mt /= kNToy;
                const double a = REFBKG[s][w], ea = REFBKGE[s][w];
                fprintf(fo, "%-7s %-28s %10.2f %8.2f %10.2f %10.2f %+8.2f %+8.2f\n", kSelN[s], w == 0 ? "T2 bkg in peak window" : "T2 bkg in bump window",
                        a, ea, b, mt, b - a, (b - a) / (std::sqrt(2.) * ea));
            }
        }
        fclose(fo);
    }

    // =================================================================================================== X17
    std::vector<XToy> X;
    double XA = 0, BA = 0, BAe2 = 0;
    {
        TFile ft(Form("%s/templates/templates.root", kBase));
        TH1D* hx = (TH1D*)ft.Get("mee_X17_x17perBR"); TH1D* hb = (TH1D*)ft.Get("mee_X17_sumT12");
        for (int b = 27; b <= 53; ++b) { XA += kSc * hx->GetBinContent(b); BA += kSc * hb->GetBinContent(b); BAe2 += std::pow(kSc * hb->GetBinError(b), 2); }
    }
    {
        TFile fk(Form("%s/sealed/x17_key.root", kBase));
        TTree* k = (TTree*)fk.Get("key");
        Int_t id, index, ibr, ns, nb, nobs; Double_t br, mus, mub;
        k->SetBranchAddress("id", &id); k->SetBranchAddress("index", &index); k->SetBranchAddress("ibr", &ibr); k->SetBranchAddress("br", &br);
        k->SetBranchAddress("mu_s_win", &mus); k->SetBranchAddress("mu_b_win", &mub); k->SetBranchAddress("n_s_win", &ns);
        k->SetBranchAddress("n_b_win", &nb); k->SetBranchAddress("n_obs_win", &nobs);
        std::map<int, int> pos;
        for (Long64_t i = 0; i < k->GetEntries(); ++i) {
            k->GetEntry(i);
            XToy x; x.id = id; x.ibr = ibr; x.br = br; x.mus = mus; x.mub = mub; x.ns = ns; x.nb = nb; x.nobs = nobs;
            for (int e = 0; e < kNEst; ++e) { x.BR[e] = x.sT[e] = x.sS[e] = x.UL[e] = std::nan(""); }
            pos[id] = X.size(); X.push_back(x);
        }
        TFile fr(Form("%s/results/x17_results.root", kBase));
        TTree* t = (TTree*)fr.Get("x17");
        Int_t rid; Double_t v[24];
        const char* br_[24] = {"BR", "eBRtot", "eBRstat", "UL90_CLs", "X2_BR", "X2_eBR", "X2_eBR", "X2_UL90_CLs",
                               "X3_BR", "X3_eBRtot", "X3_eBRstat", "X3_UL90_CLs", "X4_BR", "X4_eBRtot", "X4_eBRstat", "X4_UL90_CLs",
                               "X5_BR", "X5_eBRtot", "X5_eBRstat", "X5_UL90_CLs", "UL90_CLs_stat", "UL90_PL", "Nwin", "S"};
        t->SetBranchAddress("id", &rid);
        Double_t x2e;
        for (int j = 0; j < 24; ++j) { if (j == 6) continue; t->SetBranchAddress(br_[j], &v[j]); }
        (void)x2e;
        int nmatch = 0, nmis = 0, nwinbad = 0;
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i); v[6] = v[5];
            auto it = pos.find(rid); if (it == pos.end()) { ++nmis; continue; }
            XToy& x = X[it->second]; ++nmatch;
            for (int e = 0; e < kNEst; ++e) { x.BR[e] = v[4 * e]; x.sT[e] = v[4 * e + 1]; x.sS[e] = v[4 * e + 2]; x.UL[e] = v[4 * e + 3]; }
            x.ULstat1 = v[20]; x.ULpl1 = v[21]; x.Nwin = v[22];
            if (std::fabs(x.Nwin - x.nobs) > 0.5) ++nwinbad;
        }
        printf("X17: key %zu toys, results %lld, matched %d, unmatched %d, window-count mismatches (results Nwin vs key n_obs_win) %d\n",
               X.size(), t->GetEntries(), nmatch, nmis, nwinbad);
        printf("X17: half-A X_win %.5g per unit BR, B_win %.2f +- %.2f rows (L_data UNCONFIRMED)\n", XA, BA, std::sqrt(BAe2));
    }
    double XB = 0, BB = 0; { int c = 0; for (auto& x : X) { if (x.br > 0) { XB += x.mus / x.br; ++c; } BB = x.mub; } XB /= std::max(c, 1); }
    printf("X17: half-B X_win %.5g per unit BR, B_win %.2f rows\n", XB, BB);

    XRes XR[kNEst][6];   // [est][ibr 0..4, 5 = all]
    double linA[kNEst], linAe[kNEst], linB[kNEst], linBe[kNEst], linChi2[kNEst], slopeReal[kNEst], slopeRealE[kNEst];
    for (int e = 0; e < kNEst; ++e) {
        for (int ib = 0; ib < 6; ++ib) {
            std::vector<double> b, pT, pS, ul, zz;
            std::vector<bool> c68, cul, z1645, z2, z3, z5, z3s, zm3;
            double sTs = 0, sSs = 0;
            for (auto& x : X) {
                if (ib < 5 && x.ibr != ib) continue;
                if (!std::isfinite(x.BR[e])) continue;
                b.push_back(x.BR[e]);
                pT.push_back((x.BR[e] - x.br) / x.sT[e]); pS.push_back((x.BR[e] - x.br) / x.sS[e]);
                c68.push_back(std::fabs(x.BR[e] - x.br) <= x.sT[e]);
                ul.push_back(x.UL[e]); cul.push_back(x.UL[e] >= x.br);
                const double z = x.BR[e] / x.sT[e], zs = x.BR[e] / x.sS[e];
                z1645.push_back(z >= 1.645); z2.push_back(z >= 2); z3.push_back(z >= 3); z5.push_back(z >= 5); z3s.push_back(zs >= 3); zm3.push_back(z <= -3);
                sTs += x.sT[e]; sSs += x.sS[e];
            }
            XRes& r = XR[e][ib]; r.n = b.size(); if (!r.n) continue;
            r.mean = Mean(b); r.sd = Sd(b); r.med = Quant(b, 0.5);
            r.off = ib < 5 ? r.mean - kBR[ib] : 0; r.eOff = r.sd / std::sqrt((double)r.n);
            r.pTot_m = Mean(pT); r.pTot_s = Sd(pT); r.pStat_m = Mean(pS); r.pStat_s = Sd(pS); r.cov68 = Frac(c68);
            r.mST = sTs / r.n; r.mSS = sSs / r.n;
            r.ulMed = Quant(ul, 0.5); r.ul16 = Quant(ul, 0.16); r.ul84 = Quant(ul, 0.84); r.ul025 = Quant(ul, 0.025); r.ul975 = Quant(ul, 0.975);
            r.ulCov = Frac(cul); r.ulCovE = BinErr(r.ulCov, r.n);
            r.fz1645 = Frac(z1645); r.fz2 = Frac(z2); r.fz3 = Frac(z3); r.fz5 = Frac(z5); r.fz3stat = Frac(z3s); r.fzm3 = Frac(zm3);
        }
        // linearity: weighted straight line through the 5 means (sigma = rms / sqrt(n))
        double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
        for (int ib = 0; ib < 5; ++ib) { const double w = 1.0 / std::pow(XR[e][ib].eOff, 2); S += w; Sx += w * kBR[ib]; Sy += w * XR[e][ib].mean; Sxx += w * kBR[ib] * kBR[ib]; Sxy += w * kBR[ib] * XR[e][ib].mean; }
        const double D = S * Sxx - Sx * Sx;
        linB[e] = (S * Sxy - Sx * Sy) / D; linA[e] = (Sxx * Sy - Sx * Sxy) / D; linBe[e] = std::sqrt(S / D); linAe[e] = std::sqrt(Sxx / D);
        linChi2[e] = 0; for (int ib = 0; ib < 5; ++ib) linChi2[e] += std::pow((XR[e][ib].mean - linA[e] - linB[e] * kBR[ib]) / XR[e][ib].eOff, 2);
        // per-toy slope against the realised signal (n_s_win / X_B), all toys
        std::vector<double> yr, xr; for (auto& x : X) { if (!std::isfinite(x.BR[e])) continue; yr.push_back(x.BR[e]); xr.push_back(x.ns / XB); }
        const double vx = Cov(xr, xr), cxy = Cov(xr, yr), vy = Cov(yr, yr);
        slopeReal[e] = cxy / vx; slopeRealE[e] = std::sqrt(std::max(0.0, (vy / vx - slopeReal[e] * slopeReal[e]) / (xr.size() - 2)));
    }
    // X1 variants of the limit at BR = 0
    std::vector<double> ulS0, ulP0; for (auto& x : X) if (x.ibr == 0) { ulS0.push_back(x.ULstat1); ulP0.push_back(x.ULpl1); }
    // X5 prediction: E[BR5] = (BR X_B + B_B - B_A)/X_A
    // ------------------------------------------------------------------------------------------------ X17 table
    {
        FILE* fo = fopen(Form("%s/s5_blind_eval_x17.txt", kTab), "w");
        fprintf(fo, "# s5_blind_eval_x17.txt -- unblinding of the X17 toys (label blind_eval). SIMULATION ONLY, no CB.\n");
        fprintf(fo, "# 1500 toys (300 per injected BR) at L_data = 6080.5 nb^-1 (UNCONFIRMED); estimates from blind/results/x17_results.root, truth from blind/sealed/x17_key.root.\n");
        fprintf(fo, "# BR values in 1e-6. off = <BR-hat> - BR_injected +- rms/sqrt(n). pulls (BR-hat - BR)/sigma_tot and /sigma_stat. cov68 = |BR-hat - BR| <= sigma_tot.\n");
        fprintf(fo, "# UL = the estimator's UL90 (asymptotic CLs with sigma_tot; X2: with its sigma); ULcov = fraction UL >= BR_injected (nominal >= 0.90).\n");
        fprintf(fo, "# Z = BR-hat / sigma_tot (one-sided): fraction of toys with Z >= 1.645 / 2 / 3 / 5; Z3s = BR-hat/sigma_stat >= 3; Zm3 = Z <= -3. At BR = 0 these are false positives.\n");
        fprintf(fo, "# expected limit of the optimisation (tables/s5_opt_x17_result.txt, v1 selection): median CLs %.3f, band %.3f-%.3f (counting, B known)\n", 1e6 * kExpUL, 1e6 * kExpUL16, 1e6 * kExpUL84);
        fprintf(fo, "# half-A window: X = %.5g rows per unit BR, B = %.1f +- %.1f rows; half-B (truth): X = %.5g, B = %.1f rows\n", XA, BA, std::sqrt(BAe2), XB, BB);
        for (int e = 0; e < kNEst; ++e) {
            fprintf(fo, "\n## %s\n", kEstLong[e]);
            fprintf(fo, "%-6s %4s %9s %8s %9s %9s %7s | %7s %7s | %6s %5s %6s %5s %5s | %7s %7s %7s %7s %7s %6s %5s | %5s %5s %5s %5s %5s %5s\n",
                    "BR", "n", "<BRhat>", "rms", "median", "off", "eOff", "<sTot>", "<sStat>", "pTm", "pTs", "pSm", "pSs", "cov68",
                    "ULmed", "UL16", "UL84", "UL2.5", "UL97.5", "ULcov", "eCov", "Z1.6", "Z2", "Z3", "Z5", "Z3s", "Zm3");
            for (int ib = 0; ib < 6; ++ib) {
                const XRes& r = XR[e][ib];
                fprintf(fo, "%-6s %4d %+9.2f %8.2f %+9.2f %+9.2f %7.2f | %7.2f %7.2f | %+6.2f %5.2f %+6.2f %5.2f %5.3f | %7.2f %7.2f %7.2f %7.2f %7.2f %6.3f %5.3f | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n",
                        ib < 5 ? Form("%g", 1e6 * kBR[ib]) : "all", r.n, 1e6 * r.mean, 1e6 * r.sd, 1e6 * r.med, ib < 5 ? 1e6 * r.off : 0., 1e6 * r.eOff,
                        1e6 * r.mST, 1e6 * r.mSS, r.pTot_m, r.pTot_s, r.pStat_m, r.pStat_s, r.cov68,
                        1e6 * r.ulMed, 1e6 * r.ul16, 1e6 * r.ul84, 1e6 * r.ul025, 1e6 * r.ul975, r.ulCov, r.ulCovE,
                        r.fz1645, r.fz2, r.fz3, r.fz5, r.fz3stat, r.fzm3);
            }
            fprintf(fo, "linearity (weighted line through the 5 means): BR-hat = (%+.2f +- %.2f)e-6 + (%.4f +- %.4f) x BR, chi2 = %.1f / 3\n",
                    1e6 * linA[e], 1e6 * linAe[e], linB[e], linBe[e], linChi2[e]);
            fprintf(fo, "slope against the realised signal n_s_win / X_B (all 1500 toys, OLS): %.4f +- %.4f\n", slopeReal[e], slopeRealE[e]);
            fprintf(fo, "median UL90 at BR = 0 / expected %.3fe-6 = %.2f\n", 1e6 * kExpUL, XR[e][0].ulMed / kExpUL);
        }
        fprintf(fo, "\n## X1, other limit definitions at BR = 0: CLs with sigma_stat: median %.2f (16-84 %%: %.2f-%.2f); profile likelihood without CLs, sigma_tot: median %.2f (%.2f-%.2f) [1e-6]\n",
                1e6 * Quant(ulS0, 0.5), 1e6 * Quant(ulS0, 0.16), 1e6 * Quant(ulS0, 0.84), 1e6 * Quant(ulP0, 0.5), 1e6 * Quant(ulP0, 0.16), 1e6 * Quant(ulP0, 0.84));
        {   // coverage of the X1 sigma_stat limit at every BR
            fprintf(fo, "## X1 CLs with sigma_stat: coverage per BR:");
            for (int ib = 0; ib < 5; ++ib) { int n = 0, c = 0; for (auto& x : X) if (x.ibr == ib) { ++n; c += x.ULstat1 >= x.br; } fprintf(fo, " %g: %.3f", 1e6 * kBR[ib], double(c) / n); }
            fprintf(fo, "\n");
        }
        fprintf(fo, "\n## X5 closure: predicted <BR-hat> = (BR X_B + B_B - B_A) / X_A against measured\n");
        for (int ib = 0; ib < 5; ++ib) {
            const double pred = (kBR[ib] * XB + BB - BA) / XA;
            fprintf(fo, "BR %5g: predicted %+8.3f, measured %+8.3f +- %.3f [1e-6], difference / error %+.2f\n", 1e6 * kBR[ib], 1e6 * pred, 1e6 * XR[4][ib].mean, 1e6 * XR[4][ib].eOff,
                    (XR[4][ib].mean - pred) / XR[4][ib].eOff);
        }
        fprintf(fo, "offset term (B_B - B_A)/X_A = %+.3fe-6; its expected MC sigma sqrt(eA^2 + eB^2)/X_A with eB = eA = %.1f rows: %.3fe-6\n",
                1e6 * (BB - BA) / XA, std::sqrt(BAe2), 1e6 * std::sqrt(2 * BAe2) / XA);
        // realised signal check
        {
            int n = 0; double sp = 0;
            for (auto& x : X) if (x.br > 0) { sp += (x.ns - x.mus) / std::sqrt(x.mus); ++n; }
            fprintf(fo, "sampled signal: mean of (n_s - mu_s)/sqrt(mu_s) over %d toys with BR > 0: %+.3f\n", n, sp / n);
        }
        fclose(fo);
    }

    // =================================================================================================== FIGURES
    const int selCol[4] = {kRed + 1, kOrange + 7, kAzure + 2, kGreen + 3};
    const int selMk[4] = {20, 21, 22, 33};
    struct Panel { int q; std::vector<int> meths; const char* title; };
    std::vector<Panel> panels = {{0, {1, 2, 3, 4, 5, 6}, "peak yield (classes 4+5)"},
                                 {1, {2, 3, 7, 8, 9}, "bump yield (classes 2+3)"},
                                 {2, {4, 10}, "direct yield (class 5)"}};
    // ---------------------------------------------------------------- Fig 1: bias and pulls per method
    {
        TCanvas c("c1", "", 2100, 1500);
        c.Divide(3, 2, 0.001, 0.001);
        for (int row = 0; row < 2; ++row)
            for (int ip = 0; ip < 3; ++ip) {
                c.cd(row * 3 + ip + 1);
                gPad->SetLeftMargin(0.20); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.12);
                gPad->SetGridy(0);
                const Panel& P = panels[ip];
                const int nm = P.meths.size();
                double lo = 1e9, hi = -1e9;
                for (int m : P.meths)
                    for (int s = 0; s < 3; ++s) {
                        auto it = E.find(Key(s, m, P.q)); if (it == E.end()) continue;
                        const EtaRes& r = it->second;
                        if (row == 0) { lo = std::min(lo, 100 * (r.rel - r.eRelAll)); hi = std::max(hi, 100 * (r.rel + r.eRelAll)); }
                        else { lo = std::min(lo, r.pMuTot_m - r.pMuTot_s); hi = std::max(hi, r.pMuTot_m + r.pMuTot_s); }
                    }
                if (row == 1) { lo = std::min(lo, -1.6); hi = std::max(hi, 1.6); }
                const double pad = 0.12 * (hi - lo);
                TH1D* fr = new TH1D(Form("fr1_%d_%d", row, ip), "", nm, 0, nm);
                for (int i = 0; i < nm; ++i) fr->GetXaxis()->SetBinLabel(i + 1, kMethN[P.meths[i]]);
                fr->SetMinimum(lo - pad); fr->SetMaximum(hi + (row == 0 && ip == 0 ? 4.5 : 2.2) * pad);
                fr->GetXaxis()->SetLabelSize(0.085); fr->GetYaxis()->SetLabelSize(0.058); fr->GetYaxis()->SetTitleSize(0.065);
                fr->GetYaxis()->SetTitleOffset(1.45);
                fr->GetYaxis()->SetTitle(row == 0 ? "bias (#LTY#GT - #mu)/#mu [%]" : "pull (Y - #mu)/#sigma_{tot}");
                fr->Draw("AXIS");
                TLine l0; l0.SetLineStyle(2); l0.SetLineWidth(2); l0.DrawLine(0, 0, nm, 0);
                if (row == 1) { l0.SetLineColor(kGray + 1); l0.DrawLine(0, 1, nm, 1); l0.DrawLine(0, -1, nm, -1); }
                for (int s = 0; s < 3; ++s) {
                    TGraphErrors* gi = new TGraphErrors(); TGraphErrors* go = new TGraphErrors();
                    for (int i = 0; i < nm; ++i) {
                        auto it = E.find(Key(s, P.meths[i], P.q)); if (it == E.end()) continue;
                        const EtaRes& r = it->second; const double x = i + 0.5 + (s - 1) * 0.22;
                        const int k = gi->GetN();
                        if (row == 0) {
                            gi->SetPoint(k, x, 100 * r.rel); gi->SetPointError(k, 0, 100 * r.eRelStat);
                            go->SetPoint(k, x, 100 * r.rel); go->SetPointError(k, 0, 100 * r.eRelAll);
                        } else {
                            gi->SetPoint(k, x, r.pMuTot_m); gi->SetPointError(k, 0, r.pMuTot_s);
                        }
                    }
                    if (row == 0) { go->SetLineColor(selCol[s]); go->SetLineWidth(2); go->SetMarkerSize(0); go->Draw("PZ SAME"); }
                    gi->SetLineColor(selCol[s]); gi->SetMarkerColor(selCol[s]); gi->SetMarkerStyle(selMk[s]); gi->SetMarkerSize(2.2);
                    gi->SetLineWidth(row == 0 ? 6 : 3); gi->Draw("PZ SAME");
                }
                Txt(0.19, 0.95, Form("%s", P.title), 0.065, kBlack, 12);
                if (row == 0 && ip == 0) {
                    TLegend* lg = new TLegend(0.23, 0.70, 0.65, 0.89); lg->SetTextSize(0.055); lg->SetFillStyle(0);
                    for (int s = 0; s < 3; ++s) {
                        TGraph* gg = new TGraph(); gg->SetMarkerStyle(selMk[s]); gg->SetMarkerColor(selCol[s]); gg->SetLineColor(selCol[s]); gg->SetMarkerSize(2.0);
                        lg->AddEntry(gg, Form("%s selection", kSelN[s]), "p");
                    }
                    lg->Draw();
                }
                if (row == 0 && ip == 2) {
                    Txt(0.24, 0.85, "thick: toy statistics", 0.052, kBlack, 12);
                    Txt(0.24, 0.79, "thin: + common MC terms", 0.052, kBlack, 12);
                    Txt(0.24, 0.73, "(template, parent)", 0.052, kBlack, 12);
                }
                if (row == 1 && ip == 2) {
                    Txt(0.24, 0.85, "marker: mean, bar: rms", 0.052, kBlack, 12);
                    Txt(0.24, 0.79, "#mu = half-B expectation", 0.052, kBlack, 12);
                    Txt(0.24, 0.73, "#sigma_{tot} = quoted error", 0.052, kBlack, 12);
                }
            }
        c.cd(0);
        SaveFig(&c, "s5_blind_eval_eta_bias");
    }
    // ---------------------------------------------------------------- Fig 2/3: pull distributions
    auto PullFig = [&](const char* name, const std::vector<std::pair<int, std::pair<int, int>>>& items, int nx, int ny, const char* note) {
        TCanvas c(name, "", 700 * nx, 620 * ny);
        c.Divide(nx, ny, 0.001, 0.001);
        int ipad = 0;
        for (auto& it : items) {
            const int s = it.first, m = it.second.first, q = it.second.second;
            c.cd(++ipad);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.16);
            auto f = E.find(Key(s, m, q)); if (f == E.end()) continue;
            const EtaRes& r = f->second;
            TH1D* h = new TH1D(Form("%s_h%d", name, ipad), "", 30, -6, 6);
            TH1D* g = new TH1D(Form("%s_g%d", name, ipad), "", 30, -6, 6);
            int nout = 0;
            for (size_t i = 0; i < r.pullMuTot.size(); ++i) {
                const double p = r.pullMuTot[i], pt = r.pullTTot[i];
                if (std::fabs(p) >= 6) ++nout;
                h->Fill(std::max(-5.99, std::min(5.99, p))); g->Fill(std::max(-5.99, std::min(5.99, pt)));
            }
            h->SetFillColorAlpha(kOrange + 1, 0.55); h->SetLineColor(kRed + 1); h->SetLineWidth(3);
            g->SetLineColor(kAzure + 2); g->SetLineWidth(3); g->SetFillStyle(0);
            const double mx = std::max(h->GetMaximum(), g->GetMaximum());
            h->SetMaximum(1.55 * std::max(mx, 0.4 * r.n / std::sqrt(2 * M_PI)));
            h->SetMinimum(0);
            h->GetXaxis()->SetTitle("pull (Y - true) / #sigma_{tot}");
            h->GetYaxis()->SetTitle("toys / 0.4");
            h->GetXaxis()->SetTitleSize(0.065); h->GetYaxis()->SetTitleSize(0.065); h->GetXaxis()->SetLabelSize(0.058); h->GetYaxis()->SetLabelSize(0.058);
            h->GetYaxis()->SetTitleOffset(1.05);
            h->Draw("HIST"); g->Draw("HIST SAME");
            TF1* u = new TF1(Form("%s_u%d", name, ipad), Form("%g*exp(-0.5*x*x)", r.n * 0.4 / std::sqrt(2 * M_PI)), -6, 6);
            u->SetLineColor(kBlack); u->SetLineStyle(7); u->SetLineWidth(3); u->SetNpx(300); u->Draw("SAME");
            Txt(0.17, 0.95, Form("%s, %s yield, %s sel.", kMethN[m], kQtyN[q], kSelN[s]), 0.062, kBlack, 12);
            Txt(0.18, 0.84, Form("vs #mu: %+.2f / %.2f / %.0f%%", r.pMuTot_m, r.pMuTot_s, 100 * r.covMuTot), 0.055, kRed + 1, 12);
            Txt(0.18, 0.77, Form("vs T: %+.2f / %.2f / %.0f%%", r.pTTot_m, r.pTTot_s, 100 * r.covTTot), 0.055, kAzure + 2, 12);
            if (nout) Txt(0.18, 0.70, Form("%d toys beyond #pm6 (in edge bins)", nout), 0.048, kBlack, 12);
        }
        if (ipad < nx * ny) {
            c.cd(++ipad);
            Txt(0.06, 0.88, "filled: (Y - #mu)/#sigma_{tot}, #mu = half-B expectation", 0.060, kRed + 1, 12);
            Txt(0.06, 0.78, "line: (Y - T)/#sigma_{tot}, T = true rows of the toy", 0.060, kAzure + 2, 12);
            Txt(0.06, 0.68, "dashed: unit Gaussian", 0.060, kBlack, 12);
            Txt(0.06, 0.56, "numbers: pull mean / rms / coverage", 0.060, kBlack, 12);
            Txt(0.06, 0.48, "of the #pm1#sigma_{tot} interval (68%)", 0.060, kBlack, 12);
            Txt(0.06, 0.36, note, 0.060, kBlack, 12);
            Txt(0.06, 0.24, "200 toys, L_{toy} = 110.5 nb^{-1}, simulation only", 0.055, kBlack, 12);
        }
        c.cd(0);
        SaveFig(&c, name);
    };
    PullFig("s5_blind_eval_eta_pulls_peak", {{0, {1, 0}}, {0, {2, 0}}, {0, {3, 0}}, {0, {4, 0}}, {0, {5, 0}}, {0, {6, 0}}}, 3, 2, "");
    PullFig("s5_blind_eval_eta_pulls_bump", {{2, {9, 1}}, {2, {2, 1}}, {2, {3, 1}}, {2, {7, 1}}, {2, {8, 1}}, {1, {4, 2}}, {1, {10, 2}}}, 4, 2,
            "M2d, D5: direct yield (class 5)");
    // ---------------------------------------------------------------- Fig 4: Y against T
    {
        std::vector<std::pair<int, std::pair<int, int>>> items = {{0, {1, 0}}, {0, {2, 0}}, {0, {5, 0}}, {2, {9, 1}},
                                                                  {2, {2, 1}}, {2, {7, 1}}, {1, {4, 2}}, {1, {10, 2}}};
        TCanvas c("c4", "", 2600, 1300);
        c.Divide(4, 2, 0.001, 0.001);
        int ipad = 0;
        for (auto& it : items) {
            const int s = it.first, m = it.second.first, q = it.second.second;
            c.cd(++ipad);
            gPad->SetLeftMargin(0.22); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.16);
            const EtaRes& r = E[Key(s, m, q)];
            double lo = 1e9, hi = -1e9;
            for (size_t i = 0; i < r.Y.size(); ++i) { lo = std::min({lo, r.Y[i], r.T[i]}); hi = std::max({hi, r.Y[i], r.T[i]}); }
            const double pd = 0.08 * (hi - lo); lo -= pd; hi += pd;
            TH1D* fr = new TH1D(Form("fr4_%d", ipad), "", 10, lo, hi);
            fr->SetMinimum(lo); fr->SetMaximum(hi + 0.25 * (hi - lo));
            fr->GetXaxis()->SetTitle("true rows T"); fr->GetYaxis()->SetTitle("extracted Y [rows]");
            fr->GetXaxis()->SetTitleSize(0.065); fr->GetYaxis()->SetTitleSize(0.065); fr->GetXaxis()->SetLabelSize(0.055); fr->GetYaxis()->SetLabelSize(0.055);
            fr->GetXaxis()->SetNdivisions(505); fr->GetYaxis()->SetNdivisions(505); fr->GetYaxis()->SetTitleOffset(1.60);
            fr->Draw("AXIS");
            TLine d; d.SetLineStyle(2); d.SetLineWidth(2); d.DrawLine(lo, lo, hi, hi);
            TGraph* g = new TGraph(r.Y.size());
            for (size_t i = 0; i < r.Y.size(); ++i) g->SetPoint(i, r.T[i], r.Y[i]);
            g->SetMarkerStyle(20); g->SetMarkerSize(0.9); g->SetMarkerColor(kRed + 1); g->Draw("P SAME");
            TGraph* gm = new TGraph(1); gm->SetPoint(0, r.mu, r.meanY); gm->SetMarkerStyle(34); gm->SetMarkerSize(3.0); gm->SetMarkerColor(kBlack); gm->Draw("P SAME");
            TLine v; v.SetLineColor(kGray + 2); v.SetLineStyle(3); v.SetLineWidth(2); v.DrawLine(r.mu, lo, r.mu, hi);
            Txt(0.22, 0.95, Form("%s, %s yield, %s sel.", kMethN[m], kQtyN[q], kSelN[s]), 0.062, kBlack, 12);
            Txt(0.25, 0.84, Form("r = %.2f, slope %.2f #pm %.2f", r.r, r.slope, r.eSlope), 0.055, kBlack, 12);
            Txt(0.25, 0.77, Form("#LTY-T#GT = %+.1f, rms %.1f", r.biasT, r.sdYT), 0.055, kBlack, 12);
        }
        c.cd(0);
        SaveFig(&c, "s5_blind_eval_eta_scatter");
    }
    // ---------------------------------------------------------------- Fig 5: coverage per method
    {
        TCanvas c("c5", "", 2100, 750);
        c.Divide(3, 1, 0.001, 0.001);
        for (int ip = 0; ip < 3; ++ip) {
            c.cd(ip + 1);
            gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.12);
            const Panel& P = panels[ip]; const int nm = P.meths.size();
            TH1D* fr = new TH1D(Form("fr5_%d", ip), "", nm, 0, nm);
            for (int i = 0; i < nm; ++i) fr->GetXaxis()->SetBinLabel(i + 1, kMethN[P.meths[i]]);
            fr->SetMinimum(0); fr->SetMaximum(1.25);
            fr->GetXaxis()->SetLabelSize(0.085); fr->GetYaxis()->SetLabelSize(0.055); fr->GetYaxis()->SetTitleSize(0.062);
            fr->GetYaxis()->SetTitle("coverage of the #pm1#sigma interval");
            fr->Draw("AXIS");
            TLine l; l.SetLineStyle(2); l.SetLineWidth(2); l.DrawLine(0, 0.683, nm, 0.683);
            for (int s = 0; s < 3; ++s) {
                TGraphErrors* g1 = new TGraphErrors(); TGraphErrors* g2 = new TGraphErrors();
                for (int i = 0; i < nm; ++i) {
                    auto it = E.find(Key(s, P.meths[i], P.q)); if (it == E.end()) continue;
                    const EtaRes& r = it->second; const double x = i + 0.5 + (s - 1) * 0.22;
                    int k = g1->GetN();
                    g1->SetPoint(k, x, r.covMuTot); g1->SetPointError(k, 0, BinErr(r.covMuTot, r.n));
                    g2->SetPoint(k, x, r.covMuAll); g2->SetPointError(k, 0, 0);
                }
                g1->SetMarkerStyle(selMk[s]); g1->SetMarkerColor(selCol[s]); g1->SetLineColor(selCol[s]); g1->SetMarkerSize(2.2); g1->SetLineWidth(3); g1->Draw("PZ SAME");
                g2->SetMarkerStyle(selMk[s] == 33 ? 27 : selMk[s] + 4); g2->SetMarkerColor(selCol[s]); g2->SetMarkerSize(2.2); g2->Draw("P SAME");
            }
            Txt(0.19, 0.95, panels[ip].title, 0.065, kBlack, 12);
            if (ip == 0) {
                TLegend* lg = new TLegend(0.20, 0.13, 0.60, 0.33); lg->SetTextSize(0.050); lg->SetFillStyle(0);
                for (int s = 0; s < 3; ++s) {
                    TGraph* gg = new TGraph(); gg->SetMarkerStyle(selMk[s]); gg->SetMarkerColor(selCol[s]); gg->SetMarkerSize(2.0);
                    lg->AddEntry(gg, Form("%s selection", kSelN[s]), "p");
                }
                lg->Draw();
            }
            if (ip == 2) {
                Txt(0.22, 0.30, "full: #sigma_{tot} (quoted)", 0.052, kBlack, 12);
                Txt(0.22, 0.23, "open: #sigma_{tot} #oplus #sigma_{parent}", 0.052, kBlack, 12);
                Txt(0.22, 0.16, "dashed: 68.3%", 0.052, kBlack, 12);
            }
        }
        c.cd(0);
        SaveFig(&c, "s5_blind_eval_eta_coverage");
    }
    // ---------------------------------------------------------------- X17 figures
    const int estCol[kNEst] = {kRed + 1, kMagenta + 1, kOrange + 7, kAzure + 2, kBlack};
    const int estMk[kNEst] = {20, 25, 21, 22, 33};
    const int estUse[5] = {0, 1, 2, 3, 4};
    const int nUse = 5;
    // Fig 6: linearity
    {
        TCanvas c("c6", "", 2000, 850);
        c.Divide(2, 1, 0.001, 0.001);
        c.cd(1); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14);
        TH1D* fr = new TH1D("fr6a", "", 10, -2, 33);
        fr->SetMinimum(-100); fr->SetMaximum(85);
        fr->GetXaxis()->SetTitle("injected BR [10^{-6}]"); fr->GetYaxis()->SetTitle("extracted BR: mean #pm rms [10^{-6}]");
        fr->Draw("AXIS");
        TLine d; d.SetLineStyle(2); d.SetLineWidth(2); d.DrawLine(-2, -2, 33, 33);
        TLine z; z.SetLineStyle(3); z.SetLineColor(kGray + 1); z.DrawLine(-2, 0, 33, 0);
        TLegend* lg = new TLegend(0.17, 0.62, 0.62, 0.90); lg->SetTextSize(0.042); lg->SetFillStyle(0);
        for (int k = 0; k < nUse; ++k) {
            const int e = estUse[k];
            TGraphErrors* g = new TGraphErrors();
            for (int ib = 0; ib < 5; ++ib) { g->SetPoint(ib, 1e6 * kBR[ib] + (k - 2.0) * 0.30, 1e6 * XR[e][ib].mean); g->SetPointError(ib, 0, 1e6 * XR[e][ib].sd); }
            g->SetMarkerStyle(estMk[e]); g->SetMarkerColor(estCol[e]); g->SetLineColor(estCol[e]); g->SetMarkerSize(2.0); g->SetLineWidth(3);
            g->Draw("PLZ SAME");
            lg->AddEntry(g, Form("%s: slope %.3f #pm %.3f", kEstN[e], linB[e], linBe[e]), "pl");
        }
        lg->Draw();
        Txt(0.17, 0.95, "X17 toys, L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)", 0.048, kBlack, 12);
        c.cd(2); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14);
        TH1D* fr2 = new TH1D("fr6b", "", 10, -2, 33);
        fr2->SetMinimum(-95); fr2->SetMaximum(45);
        fr2->GetXaxis()->SetTitle("injected BR [10^{-6}]"); fr2->GetYaxis()->SetTitle("#LTextracted#GT - injected [10^{-6}]");
        fr2->Draw("AXIS");
        z.DrawLine(-2, 0, 33, 0);
        TLegend* lg2 = new TLegend(0.40, 0.30, 0.95, 0.56); lg2->SetTextSize(0.042); lg2->SetFillColor(kWhite); lg2->SetFillStyle(1001);
        for (int k = 0; k < nUse; ++k) {
            const int e = estUse[k];
            TGraphErrors* g = new TGraphErrors(); TGraphErrors* gs = new TGraphErrors();
            for (int ib = 0; ib < 5; ++ib) {
                const double x = 1e6 * kBR[ib] + (k - 2.0) * 0.30;
                g->SetPoint(ib, x, 1e6 * XR[e][ib].off); g->SetPointError(ib, 0, 1e6 * XR[e][ib].eOff);
                gs->SetPoint(ib, x, 1e6 * XR[e][ib].off); gs->SetPointError(ib, 0, 1e6 * XR[e][ib].mST);
            }
            gs->SetLineColor(estCol[e]); gs->SetLineWidth(1); gs->SetLineStyle(1); gs->SetMarkerSize(0); gs->Draw("PZ SAME");
            g->SetMarkerStyle(estMk[e]); g->SetMarkerColor(estCol[e]); g->SetLineColor(estCol[e]); g->SetMarkerSize(2.0); g->SetLineWidth(4);
            g->Draw("PZ SAME");
            lg2->AddEntry(g, kEstN[e], "p");
        }
        // X5 prediction
        TGraph* gp = new TGraph();
        for (int ib = 0; ib < 5; ++ib) gp->SetPoint(ib, 1e6 * kBR[ib], 1e6 * ((kBR[ib] * XB + BB - BA) / XA - kBR[ib]));
        gp->SetLineColor(kGray + 2); gp->SetLineStyle(2); gp->SetLineWidth(3); gp->Draw("L SAME");
        lg2->AddEntry(gp, "X5 expected: (B_{B} - B_{A})/X", "l");
        lg2->Draw();
        Txt(0.17, 0.95, "thick: error of the mean; thin: mean quoted #sigma_{tot}", 0.045, kBlack, 12);
        c.cd(0);
        SaveFig(&c, "s5_blind_eval_x17_linearity");
    }
    // Fig 7: UL distributions at BR = 0 and UL coverage
    {
        TCanvas c("c7", "", 2000, 850);
        c.Divide(2, 1, 0.001, 0.001);
        c.cd(1); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14); gPad->SetLogx();
        const int nb = 50; double ed[nb + 1]; for (int i = 0; i <= nb; ++i) ed[i] = std::pow(10., 0.0 + 2.3 * i / nb);
        TH1D* hs[kNEst]; double mx = 0;
        for (int k = 0; k < nUse; ++k) {
            const int e = estUse[k];
            hs[e] = new TH1D(Form("ul0_%d", e), "", nb, ed);
            for (auto& x : X) if (x.ibr == 0) hs[e]->Fill(std::max(1.01, std::min(199., 1e6 * x.UL[e])));
            mx = std::max(mx, hs[e]->GetMaximum());
        }
        TH1D* fr = new TH1D("fr7a", "", nb, ed); fr->SetMinimum(0); fr->SetMaximum(2.0 * mx);
        fr->GetXaxis()->SetTitle("UL_{90} at BR = 0 [10^{-6}]"); fr->GetYaxis()->SetTitle("toys / bin");
        fr->GetXaxis()->SetMoreLogLabels(); fr->Draw("AXIS");
        TBox* band = new TBox(1e6 * kExpUL16, 0, 1e6 * kExpUL84, 2.0 * mx); band->SetFillColorAlpha(kGreen + 1, 0.25); band->Draw();
        TLine el; el.SetLineColor(kGreen + 3); el.SetLineWidth(4); el.DrawLine(1e6 * kExpUL, 0, 1e6 * kExpUL, 2.0 * mx);
        TLegend* lg = new TLegend(0.50, 0.58, 0.97, 0.90); lg->SetTextSize(0.042); lg->SetFillColor(kWhite); lg->SetFillStyle(1001);
        for (int k = 0; k < nUse; ++k) {
            const int e = estUse[k];
            hs[e]->SetLineColor(estCol[e]); hs[e]->SetLineWidth(4); hs[e]->Draw("HIST SAME");
            lg->AddEntry(hs[e], Form("%s: median %.1f", kEstN[e], 1e6 * XR[e][0].ulMed), "l");
        }
        lg->AddEntry(band, Form("expected, B known: %.1f", 1e6 * kExpUL), "f");
        lg->Draw();
        Txt(0.16, 0.95, "300 toys at BR = 0, L_{data} UNCONFIRMED", 0.048, kBlack, 12);
        c.cd(2); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14);
        TH1D* fr2 = new TH1D("fr7b", "", 5, 0, 5);
        const char* lab[5] = {"0", "1", "3", "10", "30"};
        for (int i = 0; i < 5; ++i) fr2->GetXaxis()->SetBinLabel(i + 1, lab[i]);
        fr2->GetXaxis()->SetLabelSize(0.065);
        fr2->SetMinimum(0); fr2->SetMaximum(1.18);
        fr2->GetXaxis()->SetTitle("injected BR [10^{-6}]"); fr2->GetYaxis()->SetTitle("fraction with UL_{90} #geq BR");
        fr2->Draw("AXIS");
        TLine l9; l9.SetLineStyle(2); l9.SetLineWidth(2); l9.DrawLine(0, 0.9, 5, 0.9);
        TLegend* lg2 = new TLegend(0.16, 0.16, 0.60, 0.46); lg2->SetTextSize(0.045); lg2->SetFillStyle(0);
        for (int k = 0; k < nUse; ++k) {
            const int e = estUse[k];
            TGraphErrors* g = new TGraphErrors();
            for (int ib = 0; ib < 5; ++ib) { g->SetPoint(ib, ib + 0.5 + (k - 2.0) * 0.10, XR[e][ib].ulCov); g->SetPointError(ib, 0, XR[e][ib].ulCovE); }
            g->SetMarkerStyle(estMk[e]); g->SetMarkerColor(estCol[e]); g->SetLineColor(estCol[e]); g->SetMarkerSize(2.0); g->SetLineWidth(3); g->Draw("PLZ SAME");
            lg2->AddEntry(g, kEstN[e], "pl");
        }
        {   // X1 with sigma_stat only
            TGraphErrors* g = new TGraphErrors();
            for (int ib = 0; ib < 5; ++ib) { int n = 0, cc = 0; for (auto& x : X) if (x.ibr == ib) { ++n; cc += x.ULstat1 >= x.br; } const double p = double(cc) / n; g->SetPoint(ib, ib + 0.5 + 0.30, p); g->SetPointError(ib, 0, BinErr(p, n)); }
            g->SetMarkerStyle(24); g->SetMarkerColor(kRed + 1); g->SetLineColor(kRed + 1); g->SetLineStyle(2); g->SetMarkerSize(2.0); g->SetLineWidth(3); g->Draw("PLZ SAME");
            lg2->AddEntry(g, "X1 with #sigma_{stat} only", "pl");
        }
        lg2->Draw();
        Txt(0.16, 0.95, "coverage of the 90% CL upper limits", 0.048, kBlack, 12);
        c.cd(0);
        SaveFig(&c, "s5_blind_eval_x17_limits");
    }
    // Fig 8: Z at BR = 0 and pulls
    {
        TCanvas c("c8", "", 2000, 850);
        c.Divide(2, 1, 0.001, 0.001);
        for (int pnl = 0; pnl < 2; ++pnl) {
            c.cd(pnl + 1); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14);
            TH1D* hh[kNEst]; double mx = 0;
            for (int k = 0; k < nUse; ++k) {
                const int e = estUse[k];
                hh[e] = new TH1D(Form("z%d_%d", pnl, e), "", 60, -9, 6);
                for (auto& x : X) {
                    if (pnl == 0 && x.ibr != 0) continue;
                    const double v = pnl == 0 ? x.BR[e] / x.sT[e] : (x.BR[e] - x.br) / x.sT[e];
                    hh[e]->Fill(std::max(-8.99, std::min(5.99, v)));
                }
                mx = std::max(mx, hh[e]->GetMaximum());
            }
            TH1D* fr = new TH1D(Form("fr8_%d", pnl), "", 60, -9, 6); fr->SetMinimum(0); fr->SetMaximum(2.0 * mx);
            fr->GetXaxis()->SetTitle(pnl == 0 ? "Z = BR-hat / #sigma_{tot} at BR = 0" : "pull (BR-hat - BR) / #sigma_{tot}, all BR");
            fr->GetYaxis()->SetTitle("toys / 0.25"); fr->Draw("AXIS");
            if (pnl == 0) { TLine t3; t3.SetLineColor(kGray + 2); t3.SetLineStyle(2); t3.SetLineWidth(3); t3.DrawLine(3, 0, 3, 2.0 * mx); Txt(0.80, 0.45, "Z = 3", 0.050, kGray + 2, 12); }
            TLegend* lg = new TLegend(0.16, 0.60, 0.62, 0.90); lg->SetTextSize(0.042); lg->SetFillStyle(0);
            for (int k = 0; k < nUse; ++k) {
                const int e = estUse[k];
                hh[e]->SetLineColor(estCol[e]); hh[e]->SetLineWidth(4); hh[e]->Draw("HIST SAME");
                if (pnl == 0) lg->AddEntry(hh[e], Form("%s: Z#geq3 in %.1f%%, Z#leq-3 in %.1f%%", kEstN[e], 100 * XR[e][0].fz3, 100 * XR[e][0].fzm3), "l");
                else lg->AddEntry(hh[e], Form("%s: mean %+.2f, rms %.2f", kEstN[e], XR[e][5].pTot_m, XR[e][5].pTot_s), "l");
            }
            lg->Draw();
            Txt(0.16, 0.95, pnl == 0 ? "false positives: 300 toys at BR = 0" : "1500 toys, L_{data} UNCONFIRMED", 0.048, kBlack, 12);
        }
        c.cd(0);
        SaveFig(&c, "s5_blind_eval_x17_pulls");
    }
    printf("done\n");
}
