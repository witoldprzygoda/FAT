// ============================================================================
// wp4_eval.C  --  WP4 stress closure: bias, coverage, pulls, verdict (prefix wp4)
// ============================================================================
// Reads /hdd3/przygoda/HADES/FAT_out/all_report/wp4/toys/toys_<cfg>.root (wp4_toys.C) and writes
//   /hdd3/przygoda/HADES/FAT_out/all_report/wp4/verdict.txt          (machine-readable)
//   /hdd3/przygoda/HADES/FAT_out/all_report/wp4/wp4_summary.txt      (all numbers)
//   docs/latex/tables_all/wp4_V1.txt, wp4_V1.tex                     (Table V1)
//   docs/latex/figs_all/wp4_V1.pdf/.png, wp4_V1_pulls.pdf/.png       (Fig. V1)
// Criteria (SPEC WP4): |bias| < max(0.5 sigma_stat, f |true|), f = 5 % (a_eta, Y_bump), 10 % (Y_peak),
// 20 % (rho_d and the BR quantities); 68 % coverage 0.60-0.76; null scenario (g): |mean pull| < 0.5.
// sigma_stat = RMS of the toy estimates. Coverage = fraction of toys with |est - true| <= error.
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp4", kTRUE)' 'plots/all_report/wp4_eval.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TBox.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "wp4_lib.h"

namespace wp4ev {
using namespace wp4;

const int kNCfg = 6;
const char* const kCfg[kNCfg] = {"BA", "AA", "AB", "BB", "RF", "RB"};
const char* const kCfgText[kNCfg] = {
    "SPEC: truth half B (smoothed background, raw eta), templates half A",
    "method: truth and templates from half A (no MC mismatch)",
    "swap: truth half A, templates half B",
    "method: truth and templates from half B (no MC mismatch)",
    "realistic MC: truth full sample, templates = bootstrap replica (toy mod 100); error = fit (+) sigma_MC",
    "realistic MC + bootstrap bias correction: templates = parent replica, 20 children: est = 2 est_p - mean_c est_pc, error = sandwich (+) RMS_c"};

// gated quantities: estimator, quantity, relative threshold, used in null scenario
struct Gate { const char* est; int q; double f; bool ratio; };
const int kNGate = 11;
const Gate kGate[kNGate] = {
    {"F-main", qFMa, 0.05, false}, {"F-main", qFMYb, 0.05, false}, {"F-main", qFMYp, 0.10, false}, {"F-main", qFMrho, 0.20, true},
    {"F-free", qFFa2, 0.05, false}, {"F-free", qFFYb, 0.05, false}, {"F-free", qFFYp, 0.10, false}, {"F-free", qFFBRC, 0.20, true},
    {"F-gold", qFGBRG, 0.20, true},
    {"B4", qB4Yb, 0.05, false},
    {"peak SB", qPSBYp, 0.10, false}};
const int kNEst = 5;
const char* const kEst[kNEst] = {"F-main", "F-free", "F-gold", "B4", "peak SB"};
const int kGateScen[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};   // (a)-(h); (i) informational

struct Stat {
    bool have = false;
    int n = 0, nfail = 0;
    double tru = 0, asi = 0, mean = 0, med = 0, rms = 0, mErr = 0;
    double pullM = 0, pullR = 0, cov = 0, covD = 0, covT = 0, covS = 0, errMed = 0, errdMed = 0, errsMed = 0, mc = 0, mcMean = 0;
    double pullSM = 0, pullSR = 0, pullTM = 0, pullTR = 0;   // pulls with the sandwich error and with sandwich (+) sigma_MC
    double bias() const { return mean - tru; }
    double brel() const { return tru != 0 ? (mean - tru) / std::fabs(tru) : 0.; }
};
inline double Median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
// stats[cfg][scen][q]
Stat gS[kNCfg][kNScen][kNQ];
double gMC[kNCfg][kNScen][kNQ], gMCmean[kNCfg][kNScen][kNQ], gMCsum[kNCfg][kNScen][kNQ];
std::vector<double> gPullA[kNCfg][kNQ];   // pulls of scenario (a) with the role's error (filled after the loop)
bool gHaveCfg[kNCfg];

void Load(int c) {
    TChain* tu = new TChain("truth");
    TChain* t = new TChain("toys");
    const int nf = tu->Add(Form("%s/toys/toys_%s*.root", kOutDir, kCfg[c]));
    t->Add(Form("%s/toys/toys_%s*.root", kOutDir, kCfg[c]));
    gHaveCfg[c] = nf > 0 && t->GetEntries() > 0;
    if (!gHaveCfg[c]) { printf("wp4_eval: no toys for %s\n", kCfg[c]); return; }
    int scen; double val[kNQ], mc[kNQ], mcm[kNQ];
    tu->SetBranchAddress("scen", &scen); tu->SetBranchAddress("val", val);
    const bool hasMC = tu->GetBranch("mc") != nullptr;
    if (hasMC) { tu->SetBranchAddress("mc", mc); tu->SetBranchAddress("mcMean", mcm); }
    double tv[kNScen][kNQ];
    bool hs[kNScen] = {false};
    for (Long64_t i = 0; i < tu->GetEntries(); ++i) {
        tu->GetEntry(i);
        hs[scen] = true;
        for (int q = 0; q < kNQ; ++q) { tv[scen][q] = val[q]; gMC[c][scen][q] = hasMC ? mc[q] : 0.; gMCmean[c][scen][q] = hasMC ? mcm[q] : 0.; }
    }
    int toy, ok[kNQ];
    double est[kNQ], err[kNQ], errd[kNQ], errs[kNQ];
    t->SetBranchAddress("scen", &scen); t->SetBranchAddress("toy", &toy);
    t->SetBranchAddress("est", est); t->SetBranchAddress("err", err); t->SetBranchAddress("errd", errd); t->SetBranchAddress("ok", ok);
    t->SetBranchAddress("errs", errs);
    double mct[kNQ];
    const bool perToyMC = t->GetBranch("mct") != nullptr;
    if (perToyMC) t->SetBranchAddress("mct", mct);
    std::vector<double> E[kNScen][kNQ], P[kNScen][kNQ], ER[kNScen][kNQ], ED[kNScen][kNQ], ES[kNScen][kNQ], PS[kNScen][kNQ], PT[kNScen][kNQ];
    int nIn[kNScen][kNQ], nCov[kNScen][kNQ], nCovD[kNScen][kNQ], nCovT[kNScen][kNQ], nCovS[kNScen][kNQ], nFail[kNScen][kNQ];
    double asi[kNScen][kNQ];
    for (int s = 0; s < kNScen; ++s) for (int q = 0; q < kNQ; ++q) { nIn[s][q] = nCov[s][q] = nCovD[s][q] = nCovT[s][q] = nCovS[s][q] = nFail[s][q] = 0; asi[s][q] = 0; }
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        for (int q = 0; q < kNQ; ++q) {
            if (toy < 0) { asi[scen][q] = est[q]; continue; }
            if (!ok[q] || !(est[q] == est[q]) || !(err[q] > 0)) { ++nFail[scen][q]; continue; }
            const double d = est[q] - tv[scen][q];
            const double smc = perToyMC ? mct[q] : gMC[c][scen][q];
            const double et = std::sqrt(errs[q] * errs[q] + smc * smc);
            if (perToyMC) gMCsum[c][scen][q] += smc;
            E[scen][q].push_back(est[q]); P[scen][q].push_back(d / err[q]); ER[scen][q].push_back(err[q]); ED[scen][q].push_back(errd[q]);
            ES[scen][q].push_back(errs[q]); PS[scen][q].push_back(errs[q] > 0 ? d / errs[q] : 0.); PT[scen][q].push_back(et > 0 ? d / et : 0.);
            ++nIn[scen][q];
            if (std::fabs(d) <= err[q]) ++nCov[scen][q];
            if (std::fabs(d) <= errd[q]) ++nCovD[scen][q];
            if (std::fabs(d) <= errs[q]) ++nCovS[scen][q];
            if (std::fabs(d) <= et) ++nCovT[scen][q];
        }
    }
    for (int s = 0; s < kNScen; ++s) {
        if (!hs[s]) continue;
        for (int q = 0; q < kNQ; ++q) {
            Stat& X = gS[c][s][q];
            X.have = true; X.tru = tv[s][q]; X.asi = asi[s][q]; X.n = nIn[s][q]; X.nfail = nFail[s][q];
            X.mc = gMC[c][s][q]; X.mcMean = gMCmean[c][s][q];
            if (X.n < 2) continue;
            double m = 0, m2 = 0, pm = 0, p2 = 0;
            for (int i = 0; i < X.n; ++i) { m += E[s][q][i]; m2 += E[s][q][i] * E[s][q][i]; pm += P[s][q][i]; p2 += P[s][q][i] * P[s][q][i]; }
            m /= X.n; pm /= X.n;
            X.mean = m; X.rms = std::sqrt(std::max(0., m2 / X.n - m * m)); X.mErr = X.rms / std::sqrt((double)X.n);
            X.pullM = pm; X.pullR = std::sqrt(std::max(0., p2 / X.n - pm * pm));
            X.cov = (double)nCov[s][q] / X.n; X.covD = (double)nCovD[s][q] / X.n; X.covT = (double)nCovT[s][q] / X.n;
            X.covS = (double)nCovS[s][q] / X.n;
            X.med = Median(E[s][q]); X.errMed = Median(ER[s][q]); X.errdMed = Median(ED[s][q]); X.errsMed = Median(ES[s][q]);
            auto mr = [&](const std::vector<double>& v, double& m, double& r) {
                double a = 0, b = 0; for (double x : v) { a += x; b += x * x; }
                m = a / v.size(); r = std::sqrt(std::max(0., b / v.size() - m * m));
            };
            mr(PS[s][q], X.pullSM, X.pullSR); mr(PT[s][q], X.pullTM, X.pullTR);
        }
    }
    for (int q = 0; q < kNQ; ++q) {   // pulls of scenario (a) with the error of the role of this configuration
        const bool real = (c == 4 || c == 5), meth = (c == 1 || c == 3);
        gPullA[c][q] = real ? PT[0][q] : meth ? PS[0][q] : P[0][q];
    }
    if (perToyMC)   // per-toy sigma_MC: report the mean over the toys
        for (int s = 0; s < kNScen; ++s) for (int q = 0; q < kNQ; ++q) if (gS[c][s][q].n > 0) gS[c][s][q].mc = gMCsum[c][s][q] / gS[c][s][q].n;
    if (c == 4 || c == 5)   // 100 template replicas shared by 5 toys each: the MC part of the mean bias has 100 degrees of freedom
        for (int s = 0; s < kNScen; ++s)
            for (int q = 0; q < kNQ; ++q) {
                Stat& X = gS[c][s][q];
                if (X.n < 2) continue;
                const double vmc = std::max(0., X.rms * X.rms - X.errsMed * X.errsMed);
                X.mErr = std::sqrt(X.errsMed * X.errsMed / X.n + vmc / 100.);
            }
    delete t; delete tu;
}


// ---------------------------------------------------------------------------
// roles of the configurations
//   SPEC   BA         : sigma_stat = toy RMS; coverage with the fitter's error (BB-lite)          [SPEC literal]
//   METHOD AA and BB  : sigma_stat = toy RMS; coverage with the sandwich (stat) error              [no MC mismatch]
//   RF                : sigma_stat = median sandwich error; coverage with sandwich (+) sigma_MC   [realistic MC]
//   AB                : as SPEC, diagnostic only (swap)
// coverage is required in the scenarios whose truth the model contains: a, b-, b+, c-, c+, d-, d+, g, h;
// in (e) and (f) (deliberate mismodelling) only the bias criterion applies; (i) is informational
// ---------------------------------------------------------------------------
enum Role { kSpec = 0, kMethod, kReal, kDiag };
inline int RoleOf(int c) { return c == 0 ? kSpec : (c == 1 || c == 3) ? kMethod : (c == 4 || c == 5) ? kReal : kDiag; }
inline bool CovScen(int s) { return s != 7 && s != 8 && s != 11; }
inline double SigStat(const Stat& X, int c) { return RoleOf(c) == kReal ? X.errsMed : X.rms; }
inline double CovOf(const Stat& X, int c) { const int r = RoleOf(c); return r == kReal ? X.covT : r == kMethod ? X.covS : X.cov; }
inline double NullPull(const Stat& X, int c) { const int r = RoleOf(c); return r == kReal ? X.pullTM : r == kMethod ? X.pullSM : X.pullM; }

struct Verd { bool bias, cov, nul, fail; double thr; bool marg; };
// marg = every failed criterion is within 2 sigma of its limit (toy statistics: bias error of the mean; binomial error of
// the coverage sqrt(0.68 0.32 / N) = 0.021 for N = 500)
Verd Judge(const Stat& X, const Gate& g, int s, int c) {
    Verd v{true, true, true, true, 0, true};
    v.thr = std::max(0.5 * SigStat(X, c), g.f * std::fabs(X.tru));
    v.bias = std::fabs(X.bias()) < v.thr;
    if (!v.bias && std::fabs(X.bias()) - v.thr > 2 * X.mErr) v.marg = false;
    const double cv = CovOf(X, c);
    v.cov = !CovScen(s) || (cv >= 0.60 && cv <= 0.76);
    const double sb = std::sqrt(0.68 * 0.32 / std::max(1, X.n));
    if (!v.cov && (cv < 0.60 - 2 * sb || cv > 0.76 + 2 * sb)) v.marg = false;
    v.nul = (s != 9) || std::fabs(NullPull(X, c)) < 0.5;
    if (!v.nul) v.marg = false;
    v.fail = X.nfail <= 5;
    if (!v.fail) v.marg = false;
    return v;
}
inline bool Applies(const Gate& g, int s) { return !(g.ratio && s == 9); }   // ratios undefined at a_eta = 0
inline const char* ScenTag(int s) { return kScen[s].tag; }
inline int EstOf(const Gate& g) { int e = 0; while (std::string(kEst[e]) != g.est) ++e; return e; }

// relative bias for plotting/tables: scenario g relative to the true value of scenario a
inline double RelB(int c, int s, int q) {
    const Stat& X = gS[c][s][q];
    const double ref = s == 9 ? gS[c][0][q].tru : X.tru;
    return ref != 0 ? 100. * X.bias() / std::fabs(ref) : 0.;
}
inline double RelBerr(int c, int s, int q) {
    const Stat& X = gS[c][s][q];
    const double ref = s == 9 ? gS[c][0][q].tru : X.tru;
    return ref != 0 ? 100. * X.mErr / std::fabs(ref) : 0.;
}

}  // namespace wp4ev

void wp4_eval() {
    using namespace wp4ev;
    for (int c = 0; c < kNCfg; ++c) Load(c);
    const std::string tabDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all", figDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
    gSystem->mkdir(tabDir.c_str(), kTRUE); gSystem->mkdir(figDir.c_str(), kTRUE);
    // ------------------------------------------------------------------ summary (all numbers)
    FILE* sum = fopen(Form("%s/wp4_summary.txt", kOutDir), "w");
    fprintf(sum, "# wp4_summary.txt -- plots/all_report/wp4_eval.C; toys of wp4_toys.C (500 per scenario and configuration)\n");
    fprintf(sum, "# columns: cfg scen quantity N nfail true asimov mean median rms bias bias/|true| bias/rms errMed errdMed errsMed "
                 "pull(fit) mean rms, pull(sandwich) mean rms, pull(sandwich+MC) mean rms, cov(fit) cov(Fisher data) cov(sandwich) "
                 "cov(sandwich+MC) sigmaMC MCshift\n");
    for (int c = 0; c < kNCfg; ++c) {
        if (!gHaveCfg[c]) continue;
        fprintf(sum, "\n## cfg %s: %s\n", kCfg[c], kCfgText[c]);
        for (int s = 0; s < kNScen; ++s)
            for (int q = 0; q < kNQ; ++q) {
                const Stat& X = gS[c][s][q];
                if (!X.have) continue;
                fprintf(sum, "%s %-3s %-10s %3d %2d %12.6g %12.6g %12.6g %12.6g %11.5g %11.4g %8.4f %7.3f %10.4g %10.4g %10.4g "
                             "%6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %5.3f %5.3f %5.3f %5.3f %10.4g %10.4g\n",
                        kCfg[c], ScenTag(s), QName(q), X.n, X.nfail, X.tru, X.asi, X.mean, X.med, X.rms, X.bias(), X.brel(),
                        X.rms > 0 ? X.bias() / X.rms : 0., X.errMed, X.errdMed, X.errsMed, X.pullM, X.pullR, X.pullSM, X.pullSR, X.pullTM,
                        X.pullTR, X.cov, X.covD, X.covS, X.covT, X.mc, X.mcMean);
            }
    }
    fclose(sum);
    // ------------------------------------------------------------------ verdicts
    bool pass[kNCfg][kNEst], marg[kNCfg][kNEst];
    std::vector<std::string> why[kNCfg][kNEst];
    for (int c = 0; c < kNCfg; ++c) for (int e = 0; e < kNEst; ++e) { pass[c][e] = gHaveCfg[c]; marg[c][e] = true; }
    FILE* vf = fopen(Form("%s/verdict.txt", kOutDir), "w");
    fprintf(vf, "# WP4 stress closure verdict (SPEC report_all, prefix wp4) -- plots/all_report/wp4_eval.C, %s\n", "500 toys per scenario");
    fprintf(vf, "# criteria (SPEC WP4): |bias| < max(0.5 sigma_stat, f |true|), f = 0.05 (a_eta, Y_bump), 0.10 (Y_peak), 0.20 (rho_d, BR);\n");
    fprintf(vf, "#   68 %% coverage in [0.60, 0.76] (required in a b- b+ c- c+ d- d+ g h; not in the mismodelling scenarios e, f);\n");
    fprintf(vf, "#   null scenario g: |mean pull| < 0.5; at most 5 failed fits of 500. Scenario i is informational.\n");
    fprintf(vf, "# configurations: SPEC = BA (truth half B, templates half A; fitter error), METHOD = AA and BB (no template/truth MC\n");
    fprintf(vf, "#   mismatch; sandwich stat error), RF = truth full sample, templates bootstrap replicas (error = sandwich (+) sigma_MC,\n");
    fprintf(vf, "#   sigma_stat = median sandwich error), RB = RF + bootstrap bias correction 2 est - mean(20 child replicas) with error\n");
    fprintf(vf, "#   sandwich (+) child RMS, AB = swap (diagnostic). FAIL_MARGINAL = every failed criterion within 2 sigma of its limit.\n");
    fprintf(vf, "# ROW cfg scen estimator quantity true mean sigma_stat bias_rel bias/sigma_stat thr_rel coverage null_pull nfail pass_bias pass_cov pass_null\n");
    for (int c = 0; c < kNCfg; ++c) {
        if (!gHaveCfg[c]) continue;
        for (int gi = 0; gi < kNGate; ++gi) {
            const Gate& g = kGate[gi];
            const int e = EstOf(g);
            for (int s : kGateScen) {
                if (!Applies(g, s)) continue;
                const Stat& X = gS[c][s][g.q];
                if (!X.have || X.n < 10) { pass[c][e] = false; marg[c][e] = false; why[c][e].push_back(Form("%s:%s no toys", ScenTag(s), QName(g.q))); continue; }
                const Verd v = Judge(X, g, s, c);
                const double ss = SigStat(X, c);
                fprintf(vf, "ROW %s %-3s %-8s %-10s %12.6g %12.6g %11.5g %8.4f %7.3f %7.4f %5.3f %6.3f %3d %d %d %d\n", kCfg[c], ScenTag(s), g.est,
                        QName(g.q), X.tru, X.mean, ss, X.brel(), ss > 0 ? X.bias() / ss : 0., X.tru != 0 ? v.thr / std::fabs(X.tru) : 0.,
                        CovOf(X, c), NullPull(X, c), X.nfail, v.bias, v.cov, v.nul);
                if (!v.bias) why[c][e].push_back(Form("%s:%s bias %+.1f%% (thr %.1f%%)", ScenTag(s), QName(g.q), 100 * X.brel(),
                                                      X.tru != 0 ? 100 * v.thr / std::fabs(X.tru) : 0.));
                if (!v.cov) why[c][e].push_back(Form("%s:%s coverage %.2f", ScenTag(s), QName(g.q), CovOf(X, c)));
                if (!v.nul) why[c][e].push_back(Form("%s:%s null pull %.2f", ScenTag(s), QName(g.q), NullPull(X, c)));
                if (!v.fail) why[c][e].push_back(Form("%s:%s %d failed fits", ScenTag(s), QName(g.q), X.nfail));
                if (!(v.bias && v.cov && v.nul && v.fail)) { pass[c][e] = false; if (!v.marg) marg[c][e] = false; }
            }
        }
    }
    fprintf(vf, "\n# EST cfg estimator PASS|FAIL reasons (first 8)\n");
    for (int c = 0; c < kNCfg; ++c) {
        if (!gHaveCfg[c]) continue;
        for (int e = 0; e < kNEst; ++e) {
            std::string r;
            const auto& w = why[c][e];
            for (size_t i = 0; i < w.size() && i < 8; ++i) r += (i ? "; " : "") + w[i];
            if (w.size() > 8) r += Form("; ... (%zu failures in total)", w.size());
            fprintf(vf, "EST %s %-8s %s %s\n", kCfg[c], kEst[e], pass[c][e] ? "PASS" : (marg[c][e] ? "FAIL_MARGINAL" : "FAIL"), r.c_str());
        }
    }
    // MC limitation per gated quantity (RF, scenario a): sigma_MC / sigma_stat and the replica mean shift
    fprintf(vf, "\n# MC quantity sigma_stat_rel sigma_MC_rel MC_shift_rel(replica mean - true) sigma_MC/sigma_stat   [RF, scenario a]\n");
    bool mcLim[kNEst] = {false};
    for (int gi = 0; gi < kNGate; ++gi) {
        const Gate& g = kGate[gi];
        const Stat& X = gS[4][0][g.q];
        if (!X.have) continue;
        const double r = X.errsMed > 0 ? X.mc / X.errsMed : 0.;
        if (r > 1) mcLim[EstOf(g)] = true;
        fprintf(vf, "MC %-10s %8.4f %8.4f %+8.4f %6.2f\n", QName(g.q), X.errsMed / std::fabs(X.tru), X.mc / std::fabs(X.tru), X.mcMean / std::fabs(X.tru), r);
    }
    // S10 entries: |bias| of scenario (a)
    fprintf(vf, "\n# S10 quantity |bias(a)|_rel: RF (estimator as specified) RB (with bootstrap bias correction) METHOD (mean of AA, BB) SPEC (BA) AB\n");
    for (int gi = 0; gi < kNGate; ++gi) {
        const Gate& g = kGate[gi];
        auto rb = [&](int c) { return gHaveCfg[c] && gS[c][0][g.q].have ? std::fabs(gS[c][0][g.q].brel()) : -1.; };
        const double meth = std::fabs(0.5 * (gS[1][0][g.q].brel() + gS[3][0][g.q].brel()));
        fprintf(vf, "S10 %-10s %8.4f %8.4f %8.4f %8.4f %8.4f\n", QName(g.q), rb(4), rb(5), meth, rb(0), rb(2));
    }
    // decisions
    fprintf(vf, "\n# GATE estimator SPEC METHOD RF RB MC_LIMITED DECISION DECISION_WITH_BC\n");
    fprintf(vf, "#   DECISION = PASS if METHOD (AA and BB) and RF pass (estimator exactly as specified);\n");
    fprintf(vf, "#   DECISION_WITH_BC = PASS if METHOD and RB pass (estimator with the bootstrap bias correction 2 est - mean(replica est)\n");
    fprintf(vf, "#   and the replica RMS added to the error); the SPEC (BA) result is reported, see the results note for why it is not the gate;\n");
    fprintf(vf, "#   MC_LIMITED = sigma_MC > sigma_stat for a gated quantity (fall-back 'simulation-limited').\n");
    auto st = [&](int c, int e) { return !gHaveCfg[c] ? "NA" : pass[c][e] ? "PASS" : marg[c][e] ? "FAIL_MARGINAL" : "FAIL"; };
    auto comb = [&](int c1, int c2, int e) {   // both configurations
        if (pass[c1][e] && pass[c2][e]) return "PASS";
        const bool m1 = pass[c1][e] || marg[c1][e], m2 = pass[c2][e] || marg[c2][e];
        return (m1 && m2) ? "FAIL_MARGINAL" : "FAIL";
    };
    for (int e = 0; e < kNEst; ++e) {
        std::string meth = comb(1, 3, e);
        auto dec3 = [&](int cr) {
            const bool okM = (meth == "PASS"), okR = pass[cr][e];
            if (okM && okR) return std::string("PASS");
            const bool mM = meth != "FAIL", mR = pass[cr][e] || marg[cr][e];
            return std::string((mM && mR) ? "FAIL_MARGINAL" : "FAIL");
        };
        fprintf(vf, "GATE %-8s SPEC=%s METHOD=%s RF=%s RB=%s MC_LIMITED=%s DECISION=%s DECISION_WITH_BC=%s\n", kEst[e], st(0, e), meth.c_str(),
                st(4, e), st(5, e), mcLim[e] ? "yes" : "no", dec3(4).c_str(), gHaveCfg[5] ? dec3(5).c_str() : "NA");
    }
    fprintf(vf, "X17 NOT_TESTED (the X17 fit of WP9 (x17ab, x17_model.h) was still changing when WP4 ran; its toy validation is X-b)\n");
    fclose(vf);

    // ------------------------------------------------------------------ console digest
    for (int c = 0; c < kNCfg; ++c) {
        if (!gHaveCfg[c]) continue;
        printf("\n== %s (%s)\n%-10s", kCfg[c], kCfgText[c], "quantity");
        for (int s = 0; s < kNScen; ++s) printf(" %6s", ScenTag(s));
        printf("   [bias %% (g: of true(a)); then coverage of the role]\n");
        for (int gi = 0; gi < kNGate; ++gi) {
            const Gate& g = kGate[gi];
            printf("%-10s", QName(g.q));
            for (int s = 0; s < kNScen; ++s) { if (!gS[c][s][g.q].have || !Applies(g, s)) printf(" %6s", "-"); else printf(" %6.1f", RelB(c, s, g.q)); }
            printf("\n%-10s", "  cov");
            for (int s = 0; s < kNScen; ++s) { if (!gS[c][s][g.q].have || !Applies(g, s)) printf(" %6s", "-"); else printf(" %6.2f", CovOf(gS[c][s][g.q], c)); }
            printf("\n");
        }
        for (int e = 0; e < kNEst; ++e) printf("  %-8s %s\n", kEst[e], pass[c][e] ? "PASS" : "FAIL");
    }

    // ------------------------------------------------------------------ Table V1 (txt + tex)
    FILE* tt = fopen((tabDir + "/wp4_V1.txt").c_str(), "w");
    FILE* tx = fopen((tabDir + "/wp4_V1.tex").c_str(), "w");
    fprintf(tt, "# Table V1 -- WP4 stress closure (plots/all_report/wp4_eval.C); 500 toys per scenario; data S0, isBest >= 0, L = 6079.2 nb^-1\n");
    fprintf(tt, "# sigma_stat, sigma_MC, MC shift: RF configuration, scenario (a), relative to the true value\n");
    fprintf(tt, "# bias(a): relative bias in scenario (a) for SPEC (BA), swap (AB), METHOD (mean of AA and BB), RF\n");
    fprintf(tt, "# worst: largest |relative bias| over the gated scenarios a-h (g excluded), with the scenario\n");
    fprintf(tt, "# RB: RF with the bootstrap bias correction (est = 2 est - mean of 20 child replicas; error = sandwich (+) child RMS)\n");
    fprintf(tt, "# verdicts: pass / FAIL / fail(m) = every failed criterion within 2 sigma of its limit (toy statistics)\n");
    fprintf(tt, "%-8s %-10s %7s %7s %8s | %8s %8s %8s %8s %8s | %12s %12s %12s | %7s %7s %7s %7s\n", "est", "quantity", "s_stat", "s_MC", "MCshift",
            "BA(a)", "AB(a)", "METH(a)", "RF(a)", "RB(a)", "worstMETH", "worstRF", "worstRB", "SPEC", "METH", "RF", "RB");
    fprintf(tx, "%% Table V1 -- WP4 stress closure (plots/all_report/wp4_eval.C); numbers in per cent of the true value\n");
    fprintf(tx, "\\begin{tabular}{llrrr|rrrrr|rrr|cccc}\n\\toprule\n");
    fprintf(tx, " & & \\multicolumn{3}{c|}{realistic MC, scenario (a)} & \\multicolumn{5}{c|}{bias in scenario (a)} & "
                "\\multicolumn{3}{c|}{largest bias, (a)--(h)} & \\multicolumn{4}{c}{verdict} \\\\\n");
    fprintf(tx, "estimator & quantity & $\\sigma_\\mathrm{stat}$ & $\\sigma_\\mathrm{MC}$ & MC shift & BA & AB & method & RF & RF+BC & "
                "method & RF & RF+BC & SPEC & method & RF & RF+BC \\\\\n\\midrule\n");
    auto worst = [&](int c, const Gate& g, int& sw) {
        double w = 0; sw = 0;
        for (int s : kGateScen) {
            if (s == 9 || !Applies(g, s) || !gS[c][s][g.q].have) continue;
            if (std::fabs(gS[c][s][g.q].brel()) > std::fabs(w)) { w = gS[c][s][g.q].brel(); sw = s; }
        }
        return 100. * w;
    };
    const char* qTex[kNQ] = {"$a_\\eta$", "$\\kappa$", "$\\rho_d$", "$g$", "$Y_\\mathrm{bump}$", "$Y_\\mathrm{peak}$", "$Y_5$", "$Y_4$",
                             "$\\alpha_2$", "$\\alpha_3$", "$\\alpha_4$", "$\\alpha_5$", "$Y_\\mathrm{bump}$", "$Y_\\mathrm{peak}$", "$\\alpha_5/\\alpha_2$",
                             "$R_{34}$", "$\\mathrm{BR}_C$", "$\\alpha_2$", "$\\alpha_5$", "$\\mathrm{BR}_G$", "$Y_\\mathrm{bump}$", "$Y_\\mathrm{peak}$"};
    for (int gi = 0; gi < kNGate; ++gi) {
        const Gate& g = kGate[gi];
        const int e = EstOf(g);
        const Stat& R = gS[4][0][g.q];
        const double tr = std::fabs(R.tru);
        const double meth = 50. * (gS[1][0][g.q].brel() + gS[3][0][g.q].brel());
        int s1, s2, s3, s5;
        const double wA = worst(1, g, s1), wB = worst(3, g, s2), wR = worst(4, g, s3), wC = gHaveCfg[5] ? worst(5, g, s5) : 0.;
        if (!gHaveCfg[5]) s5 = 0;
        const double wM = std::fabs(wA) > std::fabs(wB) ? wA : wB;
        const int sM = std::fabs(wA) > std::fabs(wB) ? s1 : s2;
        // per-quantity verdict (all gated scenarios): 0 pass, 1 marginal fail, 2 fail
        auto qv = [&](int c) {
            if (!gHaveCfg[c]) return -1;
            int r = 0;
            for (int s : kGateScen) {
                if (!Applies(g, s) || !gS[c][s][g.q].have) continue;
                const Verd v = Judge(gS[c][s][g.q], g, s, c);
                if (!(v.bias && v.cov && v.nul && v.fail)) r = std::max(r, v.marg ? 1 : 2);
            }
            return r;
        };
        const int vS = qv(0), vM = std::max(qv(1), qv(3)), vR = qv(4), vB = qv(5);
        auto vt = [](int v) { return v < 0 ? "NA" : v == 0 ? "pass" : v == 1 ? "fail(m)" : "FAIL"; };
        auto vx = [](int v) { return v < 0 ? "--" : v == 0 ? "pass" : v == 1 ? "fail (m)" : "\\textbf{fail}"; };
        const double rbA = gHaveCfg[5] ? 100 * gS[5][0][g.q].brel() : 0.;
        fprintf(tt, "%-8s %-10s %7.1f %7.1f %+8.1f | %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f | %+8.1f (%-2s) %+8.1f (%-2s) %+8.1f (%-2s) | %7s %7s %7s %7s\n",
                g.est, QName(g.q), 100 * R.errsMed / tr, 100 * R.mc / tr, 100 * R.mcMean / tr, 100 * gS[0][0][g.q].brel(), 100 * gS[2][0][g.q].brel(),
                meth, 100 * R.brel(), rbA, wM, ScenTag(sM), wR, ScenTag(s3), wC, ScenTag(s5), vt(vS), vt(vM), vt(vR), vt(vB));
        fprintf(tx, "%s & %s & %.1f & %.1f & %+.1f & %+.1f & %+.1f & %+.1f & %+.1f & %+.1f & %+.1f (%s) & %+.1f (%s) & %+.1f (%s) & %s & %s & %s & %s \\\\\n",
                e == 4 ? "peak SB" : g.est, qTex[g.q], 100 * R.errsMed / tr, 100 * R.mc / tr, 100 * R.mcMean / tr, 100 * gS[0][0][g.q].brel(),
                100 * gS[2][0][g.q].brel(), meth, 100 * R.brel(), rbA, wM, ScenTag(sM), wR, ScenTag(s3), wC, ScenTag(s5), vx(vS), vx(vM), vx(vR), vx(vB));
    }
    fprintf(tx, "\\bottomrule\n\\end{tabular}\n");
    // per-scenario block (RF and METHOD)
    fprintf(tt, "\n# relative bias [%%] per scenario (g: in per cent of the true value of (a)); coverage of the role in brackets\n");
    for (int c : {1, 3, 4, 5, 0}) {
        if (!gHaveCfg[c]) continue;
        fprintf(tt, "## %s\n%-10s", kCfg[c], "quantity");
        for (int s = 0; s < kNScen; ++s) fprintf(tt, " %13s", ScenTag(s));
        fprintf(tt, "\n");
        for (int gi = 0; gi < kNGate; ++gi) {
            const Gate& g = kGate[gi];
            fprintf(tt, "%-10s", QName(g.q));
            for (int s = 0; s < kNScen; ++s) {
                if (!gS[c][s][g.q].have || !Applies(g, s)) { fprintf(tt, " %13s", "-"); continue; }
                fprintf(tt, " %+6.1f (%.2f)", RelB(c, s, g.q), CovOf(gS[c][s][g.q], c));
            }
            fprintf(tt, "\n");
        }
    }
    fclose(tt); fclose(tx);

    // ------------------------------------------------------------------ Fig. V1: relative bias per scenario
    bla::Style();
    gStyle->SetPadLeftMargin(0.17); gStyle->SetPadBottomMargin(0.16); gStyle->SetPadTopMargin(0.14); gStyle->SetPadRightMargin(0.03);
    const int nPan = kNGate;
    TCanvas* cv = new TCanvas("cV1", "", 1800, 1500);
    cv->Divide(4, 3, 0.002, 0.002);
    const int nPl = 6;
    const int cfgPlot[nPl] = {0, 2, 1, 3, 4, 5};
    const int col[nPl] = {kBlack, kGray + 1, kBlue + 1, kAzure + 7, kRed + 1, kMagenta + 1};
    const int mk[nPl] = {25, 26, 20, 24, 21, 33};
    const char* lab[nPl] = {"SPEC: truth B, templ. A", "swap: truth A, templ. B", "method: A/A", "method: B/B", "realistic MC (RF)",
                            "RF + bias correction"};
    const char* ttl[kNGate] = {"F-main: a_{#eta}", "F-main: Y_{bump}", "F-main: Y_{peak}", "F-main: #rho_{d}", "F-free: #alpha_{2}",
                               "F-free: Y_{bump}", "F-free: Y_{peak}", "F-free: BR_{C}", "F-gold: BR_{G}", "B4: Y_{bump}", "peak side bands: Y_{peak}"};
    std::vector<TObject*> keep;
    for (int p = 0; p < nPan; ++p) {
        const Gate& g = kGate[p];
        cv->cd(p + 1);
        const double ylim = g.f >= 0.2 ? 60. : 40.;
        TH1D* fr = new TH1D(Form("frV1_%d", p), "", 11, -0.5, 10.5);
        fr->SetMinimum(-ylim); fr->SetMaximum(ylim);
        for (int s = 0; s < 11; ++s) fr->GetXaxis()->SetBinLabel(s + 1, ScenTag(s));
        fr->GetXaxis()->SetLabelSize(0.075); fr->GetYaxis()->SetLabelSize(0.06);
        fr->GetYaxis()->SetTitle("bias / true value [%]"); fr->GetYaxis()->SetTitleSize(0.065); fr->GetYaxis()->SetTitleOffset(1.15);
        fr->GetXaxis()->SetTitle("scenario"); fr->GetXaxis()->SetTitleSize(0.065); fr->GetXaxis()->SetTitleOffset(1.0);
        fr->Draw("AXIS");
        keep.push_back(fr);
        TBox* bx = new TBox(-0.5, -100 * g.f, 10.5, 100 * g.f);
        bx->SetFillColorAlpha(kGreen + 1, 0.15); bx->SetLineWidth(0); bx->Draw();
        keep.push_back(bx);
        TLine* l0 = new TLine(-0.5, 0, 10.5, 0); l0->SetLineStyle(2); l0->Draw(); keep.push_back(l0);
        for (int k = 0; k < nPl; ++k) {
            const int c = cfgPlot[k];
            if (!gHaveCfg[c]) continue;
            TGraphErrors* gr = new TGraphErrors();
            TGraph* out = new TGraph();
            for (int s = 0; s < 11; ++s) {
                if (!Applies(g, s) || !gS[c][s][g.q].have) continue;
                const double x = s - 0.30 + 0.12 * k;
                double y = RelB(c, s, g.q);
                if (std::fabs(y) > ylim * 0.96) { out->SetPoint(out->GetN(), x, y > 0 ? ylim * 0.93 : -ylim * 0.93); continue; }
                const int n = gr->GetN();
                gr->SetPoint(n, x, y); gr->SetPointError(n, 0, RelBerr(c, s, g.q));
            }
            gr->SetMarkerStyle(mk[k]); gr->SetMarkerColor(col[k]); gr->SetLineColor(col[k]); gr->SetMarkerSize(1.3);
            gr->Draw("P"); keep.push_back(gr);
            out->SetMarkerStyle(k < 2 ? 32 : 23); out->SetMarkerColor(col[k]); out->SetMarkerSize(1.4);
            if (out->GetN()) out->Draw("P");
            keep.push_back(out);
        }
        TLatex* tl = new TLatex();
        tl->SetNDC(); tl->SetTextSize(0.075); tl->SetTextFont(42);
        tl->DrawLatex(0.19, 0.905, ttl[p]);
        keep.push_back(tl);
    }
    cv->cd(12);
    TLegend* lg = new TLegend(0.02, 0.10, 0.98, 0.95);
    lg->SetTextSize(0.07);
    for (int k = 0; k < nPl; ++k) {
        if (!gHaveCfg[cfgPlot[k]]) continue;
        TGraph* d = new TGraph(); d->SetMarkerStyle(mk[k]); d->SetMarkerColor(col[k]); d->SetMarkerSize(1.5); keep.push_back(d);
        lg->AddEntry(d, lab[k], "p");
    }
    TBox* bl = new TBox(); bl->SetFillColorAlpha(kGreen + 1, 0.15); keep.push_back(bl);
    lg->AddEntry(bl, "SPEC tolerance", "f");
    lg->Draw();
    bla::SaveFig(cv, figDir + "/wp4_V1");

    // ------------------------------------------------------------------ Fig. V1b: coverage per scenario
    TCanvas* cc = new TCanvas("cV1b", "", 1800, 1500);
    cc->Divide(4, 3, 0.002, 0.002);
    for (int p = 0; p < nPan; ++p) {
        const Gate& g = kGate[p];
        cc->cd(p + 1);
        TH1D* fr = new TH1D(Form("frV1b_%d", p), "", 11, -0.5, 10.5);
        fr->SetMinimum(0.0); fr->SetMaximum(1.05);
        for (int s = 0; s < 11; ++s) fr->GetXaxis()->SetBinLabel(s + 1, ScenTag(s));
        fr->GetXaxis()->SetLabelSize(0.075); fr->GetYaxis()->SetLabelSize(0.06);
        fr->GetYaxis()->SetTitle("68 % coverage"); fr->GetYaxis()->SetTitleSize(0.065); fr->GetYaxis()->SetTitleOffset(1.15);
        fr->GetXaxis()->SetTitle("scenario"); fr->GetXaxis()->SetTitleSize(0.065); fr->GetXaxis()->SetTitleOffset(1.0);
        fr->Draw("AXIS"); keep.push_back(fr);
        TBox* bx = new TBox(-0.5, 0.60, 10.5, 0.76); bx->SetFillColorAlpha(kGreen + 1, 0.2); bx->SetLineWidth(0); bx->Draw(); keep.push_back(bx);
        for (int k = 0; k < nPl; ++k) {
            const int c = cfgPlot[k];
            if (!gHaveCfg[c] || c == 2) continue;
            TGraph* gr = new TGraph();
            for (int s = 0; s < 11; ++s) {
                if (!Applies(g, s) || !gS[c][s][g.q].have) continue;
                gr->SetPoint(gr->GetN(), s - 0.30 + 0.12 * k, CovOf(gS[c][s][g.q], c));
            }
            gr->SetMarkerStyle(mk[k]); gr->SetMarkerColor(col[k]); gr->SetMarkerSize(1.3);
            gr->Draw("P"); keep.push_back(gr);
        }
        TLatex* tl = new TLatex(); tl->SetNDC(); tl->SetTextSize(0.075); tl->SetTextFont(42);
        tl->DrawLatex(0.19, 0.905, ttl[p]); keep.push_back(tl);
    }
    cc->cd(12);
    TLegend* lc = new TLegend(0.02, 0.10, 0.98, 0.95);
    lc->SetTextSize(0.07);
    const char* labc[nPl] = {"SPEC (fitter error)", "", "method A/A (sandwich)", "method B/B (sandwich)", "RF (sandw. #oplus #sigma_{MC})",
                             "RF + bias corr."};
    for (int k = 0; k < nPl; ++k) {
        if (k == 1 || !gHaveCfg[cfgPlot[k]]) continue;
        TGraph* d = new TGraph(); d->SetMarkerStyle(mk[k]); d->SetMarkerColor(col[k]); d->SetMarkerSize(1.5); keep.push_back(d);
        lc->AddEntry(d, labc[k], "p");
    }
    TBox* bl2 = new TBox(); bl2->SetFillColorAlpha(kGreen + 1, 0.2); keep.push_back(bl2);
    lc->AddEntry(bl2, "required 0.60-0.76", "f");
    lc->Draw();
    bla::SaveFig(cc, figDir + "/wp4_V1_coverage");
    // ------------------------------------------------------------------ Fig. V1c: pulls of scenario (a)
    {
        TCanvas* cp = new TCanvas("cV1c", "", 1800, 1500);
        cp->Divide(4, 3, 0.002, 0.002);
        const int cp3[3] = {1, 4, 5};
        const int colp[3] = {kBlue + 1, kRed + 1, kMagenta + 1};
        for (int p = 0; p < nPan; ++p) {
            const Gate& g = kGate[p];
            cp->cd(p + 1);
            double ymax = 0;
            std::vector<TH1D*> hs;
            for (int k = 0; k < 3; ++k) {
                const int c = cp3[k];
                TH1D* h = new TH1D(Form("pull_%d_%d", p, k), "", 40, -5, 5);
                h->SetDirectory(nullptr);
                if (gHaveCfg[c]) for (double x : gPullA[c][g.q]) h->Fill(std::max(-4.99, std::min(4.99, x)));
                h->SetLineColor(colp[k]); h->SetLineWidth(3);
                ymax = std::max(ymax, h->GetMaximum());
                hs.push_back(h); keep.push_back(h);
            }
            hs[0]->SetMaximum(1.35 * ymax);
            hs[0]->GetXaxis()->SetTitle("pull (estimate - true) / error"); hs[0]->GetXaxis()->SetTitleSize(0.06);
            hs[0]->GetYaxis()->SetTitle("toys / 0.25"); hs[0]->GetYaxis()->SetTitleSize(0.065);
            hs[0]->GetXaxis()->SetLabelSize(0.06); hs[0]->GetYaxis()->SetLabelSize(0.06);
            hs[0]->Draw("HIST");
            for (int k = 1; k < 3; ++k) if (gHaveCfg[cp3[k]]) hs[k]->Draw("HIST SAME");
            TF1* gs = new TF1(Form("gs_%d", p), "[0]*exp(-0.5*x*x)", -5, 5);
            gs->SetParameter(0, 500. * 0.25 / std::sqrt(2 * M_PI)); gs->SetLineColor(kBlack); gs->SetLineStyle(2); gs->SetLineWidth(2);
            gs->Draw("SAME"); keep.push_back(gs);
            TLatex* tl = new TLatex(); tl->SetNDC(); tl->SetTextSize(0.075); tl->SetTextFont(42);
            tl->DrawLatex(0.19, 0.905, ttl[p]); keep.push_back(tl);
        }
        cp->cd(12);
        TLegend* lp = new TLegend(0.02, 0.15, 0.98, 0.95);
        lp->SetTextSize(0.07);
        lp->SetHeader("scenario (a), 500 toys");
        const char* lbp[3] = {"method A/A (sandwich)", "RF (sandw. #oplus #sigma_{MC})", "RF + bias correction"};
        for (int k = 0; k < 3; ++k) { TH1D* d = new TH1D(Form("dl%d", k), "", 1, 0, 1); d->SetLineColor(colp[k]); d->SetLineWidth(3); keep.push_back(d); lp->AddEntry(d, lbp[k], "l"); }
        TF1* dg = new TF1("dg", "gaus", 0, 1); dg->SetLineColor(kBlack); dg->SetLineStyle(2); keep.push_back(dg);
        lp->AddEntry(dg, "N(0, 1)", "l");
        lp->Draw();
        bla::SaveFig(cp, figDir + "/wp4_V1_pulls");
    }
    printf("\nwp4_eval: wrote %s/verdict.txt, %s/wp4_V1.{txt,tex}, %s/wp4_V1{,_coverage}.pdf\n", kOutDir, tabDir.c_str(), figDir.c_str());
}
