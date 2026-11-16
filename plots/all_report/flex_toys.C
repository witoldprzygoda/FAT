// ============================================================================
// flex_toys.C  --  package flex of report_all: bias test of the flexible fit on pseudo-data (WP4 "RF" machinery)
// ============================================================================
// What it does (plain words):
//   Pseudo-data are made from the FULL simulation with KNOWN scales (truth):
//     truth mu_i = k CB_i (data CB, k = truth) + sum_t s_t T_t,i (full-sample templates: smoothed G_D, G_C, raw eta classes),
//   every toy: N_i ~ Poisson(mu_i); the measured CB = one of the 200 like-sign replicas of the data CB (CB measurement noise).
//   "Realistic" fit (RF, as WP4): the templates of toy t are the WP0 bootstrap replica t mod 100 (smoothing redone), so the
//   template-to-truth difference has exactly the MC noise of the real analysis; error = sandwich (Poisson) (+) sigma_MC,
//   sigma_MC = RMS of the 100 replica fits of the Asimov data (N = mu, CB = data CB).
//   "Method" fit: the same toys with the full-sample templates (no template noise): tests the estimator alone.
//   Truth scenarios: 0 = the flexible fit of the data (scales as fitted), 1 = SMASH as is (all scales 1, no tilt, k = 1).
// Pass rule (fixed before the run): a scale PASSES if |RF bias| < max(0.5 sigma_stat, 10 % |truth|) and the 68 % coverage of
// stat (+) MC is within 0.60-0.76; "marginal" = a failed criterion within 2 sigma of its limit (bias error = toy error of the
// mean (+) sigma_MC/10; coverage binomial). The truth-near-0 case (|truth| < sigma_tot) uses |bias| < 0.5 sigma_tot.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_flex", kTRUE)'
//   root -l -b -q -e "$B" -e '.L plots/all_report/flex_toys.C+O'          # compile once
//   for t in F8 F7bj; do for s in 0 1; do for p in 0 1 2 3 4; do echo "$t $s $p"; done; done; done | xargs -P 6 -L 1 sh -c \
//       'root -l -b -q -e "$B" "plots/all_report/flex_toys.C+O(\"$0\", $1, $2, 100)"'   # 5 s per 100 toys
//   root -l -b -q -e "$B" 'plots/all_report/flex_toys.C+O("F8", 0, -1)'    # evaluation (both scenarios), tables
// Output: /hdd3/.../flex/toys/toys_<tag>_s<scen>_p<part>.txt, flex_bias_<tag>.txt; docs/latex/tables_all/flex_bias_<tag>.{txt,tex}
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TRandom3.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "flex_lib.h"

namespace flext {
using namespace flex;

// model of a ladder tag of flex_fit.C: F7 / F8 / F9 (+ "b" common tilt, "bj" per-slice tilts), F7cj
inline bool ModelOf(const std::string& tag, ModelOpt& o) {
    if (tag.rfind("F7", 0) == 0) o.comps = {cC2, cC3, cC4, cC5, cGeta, cGD, cGC};
    else if (tag.rfind("F8", 0) == 0) o.comps = {cC2, cC3, cC4, cC5, cGeta, cGDom, cGDnom, cGC};
    else if (tag.rfind("F9", 0) == 0) o.comps = {cC2, cC3, cC4, cC5, cGeta, cGDom, cGDnom, cGCom, cGCnom};
    else return false;
    const std::string sfx = tag.substr(2);
    o.tilt = sfx == "b" ? 1 : sfx == "bj" ? 2 : sfx == "cj" ? 1 : 0;
    o.perSliceNorm = sfx == "cj";
    return true;
}

struct Setup {
    Layout L; Model M, Mt; Bins B; Templ T; BinT BT; RepSet RS; FitRes Fd; Vec truth; Vec muTrue; RepFits RA; Vec sigMC;
};
inline bool Prepare(const std::string& tag, int scen, Setup& s, FILE* log) {
    s.L.fc = {ClsE2()}; s.L.mlo = 0.40; s.L.mhi = 0.70;
    ModelOpt o;
    if (!ModelOf(tag, o)) { printf("flex_toys: unknown tag %s\n", tag.c_str()); return false; }
    s.M = MakeFlex(o, false, true, tag, tag);
    if (!MakeBins(s.L, s.B)) return false;
    SimSrc S; if (!S.Open()) return false;
    s.T = MakeTempl(S, s.L, s.M.Comps(), -1);
    s.BT = BinTempl(s.B, s.T, s.M);
    s.RS = MakeReplicaSet(S, s.L, s.B, s.M);
    s.Fd = Fit(s.B, s.BT, s.M, nullptr, true, &s.B.eg);
    if (!s.Fd.ok) { printf("flex_toys: data fit failed\n"); return false; }
    s.truth = s.Fd.th;
    if (scen == 1) { s.truth = s.M.start; }   // SMASH as is: scales 1, tilts 0, k = 1
    Vec V;
    Eval(s.B, s.BT, s.M, s.truth, s.muTrue, nullptr, &V, nullptr);
    // the Gaussian constraint k_E2 = 1 +- 0.05 is an auxiliary measurement: in the pseudo-experiments its centre is the truth
    // (Asimov) or the truth + 0.05 x gauss (toys)
    s.Mt = s.M;
    for (auto& c : s.Mt.cons) c.val = s.truth[c.k];
    // Asimov fits with every replica -> sigma_MC and the replica shift
    Bins BA = s.B;
    BA.y = s.muTrue;
    s.RA = RunReplicas(BA, s.Mt, s.RS, s.truth);
    s.sigMC.assign(s.M.np(), 0.);
    for (int k = 0; k < s.M.np(); ++k) s.sigMC[k] = std::sqrt(std::max(0., s.RA.cov[k][k]));
    // Asimov with the full templates (must return the truth)
    FitRes FA = Fit(BA, s.BT, s.Mt, nullptr, false);
    double dmax = 0;
    for (int k = 0; k < s.M.np(); ++k) dmax = std::max(dmax, std::fabs(FA.th[k] - s.truth[k]));
    if (log) {
        fprintf(log, "# tag %s scenario %d: bins %d, truth:", tag.c_str(), scen, s.B.n);
        for (int k = 0; k < s.M.np(); ++k) fprintf(log, " %s=%.5f", s.M.names[k].c_str(), s.truth[k]);
        fprintf(log, "\n# Asimov with the full templates returns the truth to %.3g; Asimov replica fits ok %d\n", dmax, s.RA.nok);
    }
    return true;
}

}  // namespace flext

void flex_toys(const char* tagc = "F7b", int scen = 0, int part = 0, int ntoy = 100) {
    using namespace flex;
    using namespace flext;
    const std::string tag = tagc;
    const std::string dir = std::string(kOut) + "/toys";
    gSystem->mkdir(dir.c_str(), kTRUE);
    TStopwatch sw;
    if (part >= 0) {
        FILE* out = fopen(Form("%s/toys_%s_s%d_p%d.txt", dir.c_str(), tag.c_str(), scen, part), "w");
        Setup s;
        if (!Prepare(tag, scen, s, out)) { fclose(out); return; }
        const int np = s.M.np();
        fprintf(out, "names"); for (auto& n : s.M.names) fprintf(out, " %s", n.c_str()); fprintf(out, "\n");
        fprintf(out, "truth"); for (double v : s.truth) fprintf(out, " %.8g", v); fprintf(out, "\n");
        fprintf(out, "sigmc"); for (double v : s.sigMC) fprintf(out, " %.8g", v); fprintf(out, "\n");
        fprintf(out, "asimov_rep_mean"); for (double v : s.RA.mean) fprintf(out, " %.8g", v); fprintf(out, "\n");
        TRandom3 rng(1000003u * (unsigned)(scen + 1) + 7919u * (unsigned)(part + 1));
        for (int t = part * ntoy; t < (part + 1) * ntoy; ++t) {
            Bins Bt = s.B;
            for (int i = 0; i < Bt.n; ++i) Bt.y[i] = rng.Poisson(s.muTrue[i]);
            // toy counts are independent Poisson: event-group covariance = diag(N)
            std::fill(Bt.eg.begin(), Bt.eg.end(), 0.);
            for (int i = 0; i < Bt.n; ++i) Bt.eg[(size_t)i * Bt.n + i] = Bt.y[i];
            Model Mt = s.Mt;
            for (auto& c : Mt.cons) c.val = s.truth[c.k] + c.sig * rng.Gaus();
            // measured CB = data CB replica (CB truth = data CB); the CB covariance of the fit stays the data one
            const Vec& cr = s.B.cbrep[t % kNRepC];
            for (int i = 0; i < Bt.n; ++i) Bt.cb[i] = cr[i];
            // the truth used the nominal CB with k_truth: the toy N contains k_truth CB_nominal; the fit sees the replica
            const int r = t % kNRepS;
            for (int mode = 0; mode < 2; ++mode) {   // 0 RF (replica templates), 1 method (full templates)
                const BinT& bt = mode == 0 ? s.RS.bt[r] : s.BT;
                FitRes F = Fit(Bt, bt, Mt, &s.truth, true, &Bt.eg);
                GLSRes G;
                if (F.ok && mode == 0) G = HybridGLS(Bt, bt, Mt, F, s.RS);
                fprintf(out, "toy %d mode %d ok %d nit %d chi2P %.3f chi2hyb %.3f ndf %d th", t, mode, F.ok ? 1 : 0, F.nit, F.chi2P, G.chi2, F.ndf);
                for (int k = 0; k < np; ++k) fprintf(out, " %.8g", F.ok ? F.th[k] : 0.);
                fprintf(out, " es");
                for (int k = 0; k < np; ++k) fprintf(out, " %.6g", F.ok ? std::sqrt(std::max(0., F.Cstat[k][k])) : 0.);
                fprintf(out, "\n");
            }
            if ((t + 1) % 20 == 0) { printf("flex_toys %s s%d p%d: %d toys (%.0f s)\n", tag.c_str(), scen, part, t + 1 - part * ntoy, sw.RealTime()); sw.Continue(); fflush(out); }
        }
        fclose(out);
        return;
    }
    // ------------------------------------------------------------------ evaluation of both scenarios
    FILE* tab = fopen(Form("%s/flex_bias_%s.txt", kOut, tag.c_str()), "w");
    FILE* tt = fopen(Form("%s/flex_bias_%s.txt", kTabDir, tag.c_str()), "w");
    FILE* tx = fopen(Form("%s/flex_bias_%s.tex", kTabDir, tag.c_str()), "w");
    fprintf(tt, "# flex_bias_%s.txt -- plots/all_report/flex_toys.C: bias test of the flexible fit %s on pseudo-data\n", tag.c_str(), tag.c_str());
    fprintf(tt, "# truth = full simulation x known scales + data CB (k truth); RF = templates from the 100 WP0 bootstrap replicas (the MC noise of\n");
    fprintf(tt, "# the analysis); method = full-sample templates. Toy N Poisson, measured CB = data CB replica. Error = sandwich (+) sigma_MC.\n");
    fprintf(tt, "# pass rule: |RF bias| < max(0.5 sigma_stat, 10%% |truth|) and coverage(68%%) in 0.60-0.76 (marginal within 2 sigma of a limit)\n");
    fprintf(tx, "%% flex_bias_%s.tex -- plots/all_report/flex_toys.C (bias test of the flexible fit %s)\n", tag.c_str(), tag.c_str());
    fprintf(tx, "\\begin{tabular}{llrrrrrrrrl}\n\\toprule\n");
    fprintf(tx, "scenario & scale & truth & $\\sigma_\\mathrm{stat}$ & $\\sigma_\\mathrm{MC}$ & MC shift [\\%%] & RF bias [\\%%] & bias/$\\sigma_\\mathrm{tot}$ & pull width & coverage & verdict \\\\\n\\midrule\n");
    for (int scen = 0; scen < 2; ++scen) {
        std::vector<std::string> names;
        Vec truth, sigmc, amean;
        std::vector<Vec> th[2], es[2];
        std::vector<double> c2p[2], c2h;
        int nfail[2] = {0, 0}, ndf = 0;
        for (int part = 0; part < 50; ++part) {
            std::ifstream in(Form("%s/toys_%s_s%d_p%d.txt", dir.c_str(), tag.c_str(), scen, part));
            if (!in) continue;
            std::string line;
            while (std::getline(in, line)) {
                std::istringstream is(line);
                std::string key; is >> key;
                if (key == "names") { names.clear(); std::string n; while (is >> n) names.push_back(n); }
                else if (key == "truth") { truth.clear(); double v; while (is >> v) truth.push_back(v); }
                else if (key == "sigmc") { sigmc.clear(); double v; while (is >> v) sigmc.push_back(v); }
                else if (key == "asimov_rep_mean") { amean.clear(); double v; while (is >> v) amean.push_back(v); }
                else if (key == "toy") {
                    int t, mode, ok, nit, nd; std::string w; double cp, ch;
                    is >> t >> w >> mode >> w >> ok >> w >> nit >> w >> cp >> w >> ch >> w >> nd >> w;
                    const int np = (int)names.size();
                    Vec a(np), e(np);
                    for (int k = 0; k < np; ++k) is >> a[k];
                    is >> w;
                    for (int k = 0; k < np; ++k) is >> e[k];
                    if (!ok) { nfail[mode]++; continue; }
                    th[mode].push_back(a); es[mode].push_back(e);
                    c2p[mode].push_back(cp); if (mode == 0) c2h.push_back(ch);
                    ndf = nd;
                }
            }
        }
        if (names.empty()) continue;
        const int np = (int)names.size();
        fprintf(tab, "\n# scenario %d (%s): toys RF %zu (failed %d), method %zu (failed %d)\n", scen, scen == 0 ? "truth = flexible fit of the data" : "truth = SMASH as is",
                th[0].size(), nfail[0], th[1].size(), nfail[1]);
        auto meanrms = [](const std::vector<double>& v, double& m, double& r) { m = 0; r = 0; for (double x : v) m += x; m /= std::max<size_t>(1, v.size()); for (double x : v) r += (x - m) * (x - m); r = std::sqrt(r / std::max<size_t>(1, v.size() - 1)); };
        double mP, rP, mH, rH, mPm, rPm;
        meanrms(c2p[0], mP, rP); meanrms(c2h, mH, rH); meanrms(c2p[1], mPm, rPm);
        fprintf(tab, "gof scenario %d ndf %d RF Pearson mean %.1f rms %.1f | RF hybrid GLS mean %.1f rms %.1f | method Pearson mean %.1f rms %.1f\n", scen, ndf, mP, rP, mH, rH, mPm, rPm);
        // quantiles of the RF hybrid chi2 (for the p-value of the data)
        std::vector<double> sh = c2h; std::sort(sh.begin(), sh.end());
        if (!sh.empty()) fprintf(tab, "gofq scenario %d hybrid GLS quantiles 5%% %.1f 50%% %.1f 95%% %.1f max %.1f\n", scen, sh[(size_t)(0.05 * sh.size())], sh[sh.size() / 2], sh[(size_t)(0.95 * sh.size())], sh.back());
        std::vector<double> sp = c2p[0]; std::sort(sp.begin(), sp.end());
        if (!sp.empty()) fprintf(tab, "gofq scenario %d Pearson quantiles 5%% %.1f 50%% %.1f 95%% %.1f max %.1f\n", scen, sp[(size_t)(0.05 * sp.size())], sp[sp.size() / 2], sp[(size_t)(0.95 * sp.size())], sp.back());
        fprintf(tab, "%-10s %9s %8s %8s %8s %9s %9s %8s %8s %8s %8s %9s %9s %s\n", "param", "truth", "s_stat", "s_MC", "MCshift%", "RFbias%", "RFbias/st", "pullmean", "pullwid", "cover68",
                "meth%", "meth/sst", "methwid", "verdict");
        for (int k = 0; k < np; ++k) {
            const double tr = truth[k];
            // RF
            std::vector<double> v, vs, pulls, vm, pm;
            int cov = 0;
            for (size_t q = 0; q < th[0].size(); ++q) {
                v.push_back(th[0][q][k]); vs.push_back(es[0][q][k]);
                const double et = std::sqrt(es[0][q][k] * es[0][q][k] + sigmc[k] * sigmc[k]);
                const double p = et > 0 ? (th[0][q][k] - tr) / et : 0.;
                pulls.push_back(p); if (std::fabs(p) < 1.) ++cov;
            }
            for (size_t q = 0; q < th[1].size(); ++q) { vm.push_back(th[1][q][k]); pm.push_back(es[1][q][k] > 0 ? (th[1][q][k] - tr) / es[1][q][k] : 0.); }
            double m, r, ms, rs, pmn, pr, mm, rm, pmm, prm;
            meanrms(v, m, r); meanrms(vs, ms, rs); meanrms(pulls, pmn, pr); meanrms(vm, mm, rm); meanrms(pm, pmm, prm);
            const bool fixedPar = r == 0 && sigmc[k] == 0;
            if (fixedPar) { fprintf(tab, "%-10s %9.4f fixed\n", names[k].c_str(), tr); continue; }
            const double sStat = ms, sTot = std::sqrt(ms * ms + sigmc[k] * sigmc[k]);
            const double bias = m - tr, biasErr = std::sqrt(r * r / std::max<size_t>(1, v.size()) + sigmc[k] * sigmc[k] / 100.);
            const double cv = (double)cov / std::max<size_t>(1, pulls.size());
            const double tol = std::fabs(tr) >= sTot ? std::max(0.5 * sStat, 0.10 * std::fabs(tr)) : 0.5 * sTot;
            const bool bOK = std::fabs(bias) < tol, bMarg = std::fabs(bias) < tol + 2 * biasErr;
            const double ce = 0.021 * std::sqrt(500. / std::max<size_t>(1, pulls.size()));
            const bool cOK = cv >= 0.60 && cv <= 0.76, cMarg = cv >= 0.60 - 2 * ce && cv <= 0.76 + 2 * ce;
            const char* verdict = (bOK && cOK) ? "PASS" : (bMarg && cMarg) ? "FAIL(marginal)" : "FAIL";
            const bool isK = names[k].rfind("k_", 0) == 0;
            const double rel = tr != 0 ? 100. * bias / tr : 0.;
            const double mcs = tr != 0 ? 100. * (amean[k] - tr) / tr : 0.;
            const double mrel = tr != 0 ? 100. * (mm - tr) / tr : 0.;
            fprintf(tab, "%-10s %9.4f %8.4f %8.4f %8.1f %9.1f %9.2f %8.2f %8.2f %8.3f %8.1f %9.2f %9.2f %s\n", names[k].c_str(), tr, sStat, sigmc[k], mcs, rel, bias / sTot,
                    pmn, pr, cv, mrel, sStat > 0 ? (mm - tr) / sStat : 0., prm, verdict);
            fprintf(tt, "scen %d %-10s truth %.4f s_stat %.4f s_MC %.4f MCshift %.1f%% RFbias %.1f%% (%.2f s_tot, +- %.4f) pull %.2f width %.2f cover %.3f | method %.1f%% width %.2f | %s\n",
                    scen, names[k].c_str(), tr, sStat, sigmc[k], mcs, rel, bias / sTot, biasErr, pmn, pr, cv, mrel, prm, isK ? "(nuisance)" : verdict);
            std::string nm = names[k];
            for (size_t p = 0; p < nm.size(); ++p) if (nm[p] == '_') { nm.insert(p, "\\"); ++p; }
            fprintf(tx, "%s & \\texttt{%s} & %.3f & %.3f & %.3f & %+.1f & %+.1f & %+.2f & %.2f & %.2f & %s \\\\\n", scen == 0 ? "data fit" : "SMASH", nm.c_str(), tr, sStat, sigmc[k],
                    mcs, rel, bias / sTot, pr, cv, isK ? "nuisance" : verdict);
        }
        if (scen == 0) fprintf(tx, "\\midrule\n");
    }
    fprintf(tx, "\\bottomrule\n\\end{tabular}\n");
    fclose(tab); fclose(tt); fclose(tx);
    printf("flex_toys: evaluation written (%s/flex_bias_%s.txt)\n", kOut, tag.c_str());
}
