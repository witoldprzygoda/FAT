// ============================================================================
// verify_smp_lumi_fit.C  --  independent check of the luminosity of the
//                            analysed runs, the per-day rates and the eta peak
//                            per period of the chapter "samples, luminosity
//                            and cut flow" (smp)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Verification macro (written by the verifier of the smp chapter); it does not
// use smp_common.h or ScanEngine.h.
// Inputs:
//   the ORIGINAL luminosity files (scratch copy of vae25, revision N):
//     .../scratchpad/final/remote/lumi_startfree_perrun.csv, lumi_coverage_perrun.csv,
//     lumi_startfree_perday.csv
//   outputs of verify_smp_counts.C (report_exp/verify_smp/counts.txt, runs_us.txt, hists.root)
// What it does:
//   1. L of all runs, of the runs with merged DSTs (n_dst_merged > 0), of the runs
//      with FAT e+e- rows; per period; L-weighted PT3/PT2 live fraction x 64/63
//   2. per-day rate R = rows(pre) x 64/63 / L(day), weighted mean, chi2/ndf,
//      rms, and the plain ratio sum(rows) x 64/63 / sum(L)
//   3. eta peak per period with a SIMPLE chi2 fit (not ScanEngine): N - CB in
//      2.5 MeV bins, CB = 2 sqrt(N++ N--) in 10 MeV bins linearly interpolated,
//      diagonal errors sqrt(N + var CB), Gaussian 0.5439 / 5.4 MeV integrated over
//      the bin + pol2, fit range 0.5225-0.6100 GeV/c2; s = S x 64/63 / L
// Output: report_exp/verify_smp/lumi_fit.txt
// Run (from /home/przygoda/HADES/FAT/FAT, after verify_smp_counts.C):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp_verify", kTRUE)' 'plots/exp_report/verify_smp_lumi_fit.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

namespace vlf {
const char* kRemote = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote";
const char* kOut = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/verify_smp";
std::vector<std::string> Split(const std::string& s) {
    std::vector<std::string> v; std::string c; bool q = false;
    for (char ch : s) { if (ch == '"') { q = !q; continue; } if (ch == ',' && !q) { v.push_back(c); c.clear(); } else if (ch != '\r') c += ch; }
    v.push_back(c); return v;
}
// rows (as maps column -> value) of a CSV with '#' comments
std::vector<std::map<std::string, std::string>> Csv(const std::string& p) {
    std::vector<std::map<std::string, std::string>> out;
    std::ifstream in(p.c_str()); std::string line; std::vector<std::string> head;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (head.empty()) { head = Split(line); continue; }
        const std::vector<std::string> v = Split(line);
        std::map<std::string, std::string> m;
        for (size_t i = 0; i < head.size() && i < v.size(); ++i) m[head[i]] = v[i];
        out.push_back(m);
    }
    return out;
}
int Period(int d) { if (d >= 34 && d <= 46) return 0; if (d >= 47 && d <= 48) return 1; if (d >= 49 && d <= 53) return 2; if (d >= 59 && d <= 67) return 3; return -1; }
}

void verify_smp_lumi_fit()
{
    using namespace vlf;
    FILE* o = fopen(Form("%s/lumi_fit.txt", kOut), "w");
    const double f6463 = 64.0 / 63.0;
    // ---- 1. luminosity ----
    auto pr = Csv(std::string(kRemote) + "/lumi_startfree_perrun.csv");
    auto cv = Csv(std::string(kRemote) + "/lumi_coverage_perrun.csv");
    auto pd = Csv(std::string(kRemote) + "/lumi_startfree_perday.csv");
    std::map<int, double> L6run; std::map<int, int> dayRun; std::map<int, int> noHF;
    for (auto& r : pr) { const int run = std::atoi(r["runnumber"].c_str()); L6run[run] = std::atof(r["L6_SF"].c_str()); dayRun[run] = std::atoi(r["day"].c_str()); noHF[run] = std::atoi(r["no_hf"].c_str()); }
    std::map<int, int> merged;
    for (auto& r : cv) { const long rid = std::atol(r["runid"].c_str()); merged[(int)(rid - 443000000L)] = std::atoi(r["n_dst_merged"].c_str()); }
    std::set<int> fatRuns; { std::ifstream in(Form("%s/runs_us.txt", kOut)); int r; while (in >> r) fatRuns.insert(r); }
    double Lall = 0, Lmer = 0, Lfat = 0, LperA[4] = {}, LperM[4] = {};
    std::map<int, double> LdayM;
    int nAll = 0, nMer = 0, nFatNotInTable = 0, nNoHF = 0, nMerNoHF = 0;
    for (auto& p : L6run) {
        const int run = p.first; const double L = p.second; const int d = dayRun[run];
        nAll++; Lall += L; LperA[Period(d)] += L; if (noHF[run]) nNoHF++;
        if (merged.count(run) && merged[run] > 0) { nMer++; Lmer += L; LperM[Period(d)] += L; LdayM[d] += L; if (noHF[run]) nMerNoHF++; }
        if (fatRuns.count(run)) Lfat += L;
    }
    for (int r : fatRuns) if (!L6run.count(r)) nFatNotInTable++;
    fprintf(o, "LUMI all runs %d L6 %.3f ; runs with merged DSTs %d (no_hf among them %d) L %.3f ; runs with FAT e+e- rows %zu (not in table %d) L %.3f ; no_hf runs %d\n",
            nAll, Lall, nMer, nMerNoHF, Lmer, fatRuns.size(), nFatNotInTable, Lfat, nNoHF);
    for (int q = 0; q < 4; ++q) fprintf(o, "LUMI period %d all %.3f merged %.3f\n", q, LperA[q], LperM[q]);
    // live fraction
    double sL = 0, sLf = 0, fmin = 9, fmax = 0;
    for (auto& r : pd) {
        const int d = std::atoi(r["day"].c_str());
        const double f = std::atof(r["LF3_over_LF2"].c_str()) * f6463;
        sL += LdayM[d]; sLf += LdayM[d] * f; fmin = std::min(fmin, f); fmax = std::max(fmax, f);
    }
    fprintf(o, "LIVE L-weighted LF3/LF2 x 64/63 = %.4f (per day %.4f-%.4f); L x f = %.1f; relative change %.2f %%\n", sLf / sL, fmin, fmax, sLf, 100 * (sLf / sL - 1));
    // ---- 2. rates ----
    std::map<std::string, std::map<int, std::vector<double>>> day;   // hyp -> day -> rows, var, best
    {
        std::ifstream in(Form("%s/counts.txt", kOut)); std::string line;
        while (std::getline(in, line)) {
            if (line.compare(0, 4, "DAY ") != 0) continue;
            std::istringstream is(line); std::string k, h; int d; double n, v, b; is >> k >> h >> d >> n >> v >> b;
            day[h][d] = {n, v, b};
        }
    }
    for (const char* h : {"us", "pp", "mm"}) {
        double sw = 0, swx = 0, sN = 0, sLd = 0; std::vector<double> R, E;
        double swa = 0, swxa = 0, swb = 0, swxb = 0, Na = 0, La = 0, Nb = 0, Lb = 0;
        for (auto& p : day[h]) {
            const int d = p.first; const double L = LdayM[d];
            const double r = p.second[0] * f6463 / L, e = std::sqrt(p.second[1]) * f6463 / L;
            R.push_back(r); E.push_back(e); sw += 1 / (e * e); swx += r / (e * e); sN += p.second[0]; sLd += L;
            if (d <= 44) { swa += 1 / (e * e); swxa += r / (e * e); Na += p.second[0]; La += L; } else { swb += 1 / (e * e); swxb += r / (e * e); Nb += p.second[0]; Lb += L; }
        }
        const double m = swx / sw; double chi2 = 0, s2 = 0;
        for (size_t i = 0; i < R.size(); ++i) { chi2 += (R[i] - m) * (R[i] - m) / (E[i] * E[i]); s2 += (R[i] / m - 1) * (R[i] / m - 1); }
        fprintf(o, "RATE %s weighted mean %.3f +- %.3f chi2/ndf %.0f/%zu rms %.4f ; plain sum ratio %.3f (sum rows %.0f, sum L %.3f) ; "
                   "034-044 weighted %.3f plain %.3f ; 045-067 weighted %.3f plain %.3f ; ratio weighted %.4f plain %.4f\n",
                h, m, 1 / std::sqrt(sw), chi2, R.size() - 1, std::sqrt(s2 / R.size()), sN * f6463 / sLd, sN, sLd,
                swxa / swa, Na * f6463 / La, swxb / swb, Nb * f6463 / Lb, (swxa / swa) / (swxb / swb), (Na / La) / (Nb / Lb));
    }
    // ---- 3. peak per period ----
    TFile* f = TFile::Open(Form("%s/hists.root", kOut), "READ");
    const double mu = 0.5439, sg = 0.0054, bw = 0.0025;
    TF1* fn = new TF1("fpk", [=](double* x, double* p) {
        const double a = (x[0] - 0.5 * bw - mu) / (std::sqrt(2.) * sg), b = (x[0] + 0.5 * bw - mu) / (std::sqrt(2.) * sg);
        const double u = x[0] - 0.56;
        return p[0] * 0.5 * (TMath::Erf(b) - TMath::Erf(a)) + p[1] + p[2] * u + p[3] * u * u; }, 0.5225, 0.6100, 4);
    const char* gname[5] = {"I", "IIA", "IIB", "III", "all"};
    for (const char* s : {"pre", "pnomv"}) {
        for (int g = 0; g < 5; ++g) {
            TH1D* hu = (TH1D*)f->Get(Form("m4_us_%s_%d", s, g));
            TH1D* hp = (TH1D*)f->Get(Form("m4_pp_%s_%d", s, g));
            TH1D* hm = (TH1D*)f->Get(Form("m4_mm_%s_%d", s, g));
            const int nb = hu->GetNbinsX(), nw = nb / 4;
            std::vector<double> cbW(nw), cbV(nw);
            for (int j = 0; j < nw; ++j) {
                double a = 0, b = 0; for (int i = 0; i < 4; ++i) { a += hp->GetBinContent(4 * j + i + 1); b += hm->GetBinContent(4 * j + i + 1); }
                cbW[j] = 2 * std::sqrt(a * b); cbV[j] = a + b;
            }
            TH1D* hs = new TH1D(Form("sig_%s_%d", s, g), "", nb, 0.35, 0.95);
            double Npk = 0, CBpk = 0;
            for (int i = 0; i < nb; ++i) {
                const double x = hu->GetBinCenter(i + 1);
                const double t = (x - (0.35 + 0.5 * 4 * bw)) / (4 * bw);
                int j0, j1; double fr;
                if (t <= 0) { j0 = j1 = 0; fr = 0; } else if (t >= nw - 1) { j0 = j1 = nw - 1; fr = 0; } else { j0 = (int)t; j1 = j0 + 1; fr = t - j0; }
                const double cb = ((1 - fr) * cbW[j0] + fr * cbW[j1]) / 4, cv = ((1 - fr) * (1 - fr) * cbV[j0] + fr * fr * cbV[j1]) / 16;
                const double n = hu->GetBinContent(i + 1);
                hs->SetBinContent(i + 1, n - cb); hs->SetBinError(i + 1, std::sqrt(std::max(n + cv, 1.0)));
                if (x > 0.53364 && x < 0.55416) { Npk += n; CBpk += cb; }
            }
            fn->SetParameters(300, 100, 0, 0);
            hs->Fit(fn, "Q0R");
            const double S = fn->GetParameter(0), eS = fn->GetParError(0);
            const double Lg = (g < 4) ? LperM[g] : Lmer;
            fprintf(o, "PEAK %-5s %-4s L %.2f S %.1f +- %.1f chi2/ndf %.1f/%d s %.4f +- %.4f ; CB/N (bin centres in data window) %.4f\n", s, gname[g], Lg, S, eS,
                    fn->GetChisquare(), fn->GetNDF(), S * f6463 / Lg, eS * f6463 / Lg, CBpk / Npk);
        }
    }
    fclose(o);
    printf("verify_smp_lumi_fit done\n");
}
