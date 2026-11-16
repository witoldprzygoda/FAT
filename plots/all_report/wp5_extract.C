// ============================================================================
// wp5_extract.C  --  WP5 main extraction, stage B (SPEC report_all, prefix wp5)
// ============================================================================
// Fits of the data F10 / F5g histograms of hists_v1.root with the stage-B model (spec Sect. 2, D10-D13):
//   F-main (nominal after the fall-backs of WP3: M(ee) >= 10 MeV only, kappa = 1 +- 0.2), F-main as specified (E1 + E2),
//   F-free, F-gold, per-slice a_eta, folds, pre, data isBest == 1, periods (for the WP2 stability test), the S3/S4
//   model variants that the fall-back "stage B chi2 > 1.5" asks for, profile scans (a_eta, rho_d, kappa, s).
// Errors: sandwich (event groups of N, CB replicas, constraints) + MC (100 WP0 simulation replicas, refitted).
// WP4 verdict (verdict.txt): F-main, F-free and F-gold FAILED the stress closure -> every number here is a
// CROSS-CHECK (labelled); the main results are the side-band estimators of WP6. This macro provides the inputs that
// WP6 (other-eta terms), WP7 (variants), WP8 (ratios) and the X17 package (a_eta, background shapes) need.
//
// Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt (machine-readable), wp5_extract.log, wp5_fits.root;
//          docs/latex/tables_all/wp5_tableE1.{txt,tex}, wp5_tableE2.{txt,tex}, wp5_acceptance.txt;
//          docs/latex/figs_all/wp5_E1.pdf, wp5_E1_spec.pdf, wp5_E2.pdf, wp5_E3.pdf (+ .png).
// Run (from /home/przygoda/HADES/FAT/FAT; after wp5_rows.C; about 3 min, 3 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp5", kTRUE)' \
//        'plots/all_report/wp5_extract.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TBox.h"

#include "wp5_lib.h"

namespace wp5m {
using namespace wp5;

FILE* gLog = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); } while (0)

// one fit with all its errors
struct Block {
    std::string tag, text;
    Spec S;
    Bins B;
    FitRes F;
    Replicas Rp;
    bool hasRep = false;
    std::vector<DVal> dv;
    GLSRes G;
    bool hasGLS = false;
    const Templ* T = nullptr;
    const std::vector<Templ>* R = nullptr;
    const DataSet* D = nullptr;
    Mat CovMC() const { return hasRep ? Rp.cov : Zero(S.np(), S.np()); }
    Mat CovTot() const { return Add(F.Cstat, CovMC()); }
    double Err(int k, int what) const {   // 0 N, 1 CB, 2 stat, 3 MC, 4 total
        const Mat C = what == 0 ? F.Cn : what == 1 ? F.Ccb : what == 2 ? F.Cstat : what == 3 ? CovMC() : CovTot();
        return std::sqrt(std::max(0., C[k][k]));
    }
    const DVal* Get(const std::string& n) const { for (auto& d : dv) if (d.name == n) return &d; return nullptr; }
    double V(const std::string& n) const { const int k = S.Idx(n); return k >= 0 ? F.th[k] : NAN; }
};

Block Run(const std::string& tag, const std::string& text, const DataSet& D, const Templ& T, const std::vector<Templ>* R,
          const Spec& S, int lay, int useE, const std::vector<Deriv>& der, bool gls,
          const std::function<double(int, int)>& fixfn = nullptr,
          const std::function<void(Bins&, const Templ&, int)>& rebuild = nullptr, double mlo = kF10Lo, double mhi = kF10Hi) {
    Block b; b.tag = tag; b.text = text; b.S = S; b.T = &T; b.R = R; b.D = &D;
    MakeBins(D, T, lay, useE, b.B, mlo, mhi, fixfn);
    b.F = Fit(b.B, S, D);
    if (!b.F.ok) { P("  [%s] FIT FAILED (it %d)\n", tag.c_str(), b.F.nit); return b; }
    if (R) {
        if (rebuild) {
            // replica-dependent fixed terms (F-gold): rebuild the bins for every replica
            b.Rp.fit.resize(R->size());
            for (size_t r = 0; r < R->size(); ++r) {
                Bins Br = b.B;
                rebuild(Br, (*R)[r], (int)r);
                b.Rp.fit[r] = Fit(Br, S, D, &b.F.th, false);
            }
            const int np = S.np();
            b.Rp.mean.assign(np, 0.); b.Rp.cov = Zero(np, np); b.Rp.nok = 0;
            for (auto& f : b.Rp.fit) if (f.ok) { ++b.Rp.nok; for (int k = 0; k < np; ++k) b.Rp.mean[k] += f.th[k]; }
            for (int k = 0; k < np; ++k) b.Rp.mean[k] /= std::max(1, b.Rp.nok);
            for (auto& f : b.Rp.fit) if (f.ok) for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) b.Rp.cov[k][l] += (f.th[k] - b.Rp.mean[k]) * (f.th[l] - b.Rp.mean[l]);
            for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) b.Rp.cov[k][l] /= std::max(1, b.Rp.nok - 1);
        } else {
            b.Rp = RunReplicas(b.B, S, D, *R, b.F.th);
        }
        b.hasRep = true;
    }
    for (const auto& d : der) b.dv.push_back(Evaluate(d, b.F, T, b.hasRep ? &b.Rp : nullptr, R));
    if (gls && R) { b.G = GLS(b.B, S, D, b.F, *R); b.hasGLS = b.G.ok; }
    return b;
}

void PrintBlock(const Block& b) {
    P("\n--- %s: %s\n", b.tag.c_str(), b.text.c_str());
    if (!b.F.ok) { P("    failed\n"); return; }
    P("    bins %d, free %d, ndf %d, iterations %d; chi2: Pearson(V) %.1f (%.3f), Pearson(V x event-group factor) %.1f (%.3f), Baker-Cousins %.1f",
      b.F.n, b.F.nfree, b.F.ndf, b.F.nit, b.F.chi2P, b.F.chi2P / b.F.ndf, b.F.chi2E, b.F.chi2E / b.F.ndf, b.F.chi2BC);
    if (b.hasGLS) P("; GLS with MC covariance: at the fit %.1f (%.3f), diag data cov %.1f (%.3f), minimised %.1f (%.3f)",
                    b.G.chi2, b.G.chi2 / b.G.ndf, b.G.chi2diag, b.G.chi2diag / b.G.ndf, b.G.chi2min, b.G.chi2min / b.G.ndf);
    P("\n");
    if (b.hasRep) P("    replicas: %d of %zu converged\n", b.Rp.nok, b.Rp.fit.size());
    P("    %-10s %12s %10s %10s %10s %10s %10s %12s\n", "param", "value", "N", "CB", "stat", "MC", "total", "MC shift");
    for (int k = 0; k < b.S.np(); ++k)
        P("    %-10s %12.5g %10.4g %10.4g %10.4g %10.4g %10.4g %12.4g%s\n", b.S.names[k].c_str(), b.F.th[k], b.Err(k, 0), b.Err(k, 1),
          b.Err(k, 2), b.Err(k, 3), b.Err(k, 4), b.hasRep ? b.Rp.mean[k] - b.F.th[k] : 0., b.S.fixed[k] ? "  (fixed)" : "");
    if (!b.dv.empty()) {
        P("    %-12s %12s %10s %10s %10s %10s %10s %10s %12s\n", "derived", "value", "N", "CB", "constr", "stat", "MC", "total", "bias-corr.");
        for (const auto& d : b.dv) {
            double ec = 0; for (double x : d.eCons) ec += x * x;
            P("    %-12s %12.5g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %12.5g\n", d.name.c_str(), d.v, d.eN, d.eCB, std::sqrt(ec), d.eStat, d.eMC, d.eTot, d.bc);
        }
    }
}

void WriteMat(FILE* o, const char* key, const std::string& tag, const Mat& C) {
    fprintf(o, "%s %s %zu\n", key, tag.c_str(), C.size());
    for (const auto& r : C) { for (double x : r) fprintf(o, " %.8e", x); fprintf(o, "\n"); }
}
void WriteBlock(FILE* o, const Block& b) {
    fprintf(o, "\n# ---------------------------------------------------------------------------------------------\n");
    fprintf(o, "fit %s ok %d\n", b.tag.c_str(), b.F.ok ? 1 : 0);
    fprintf(o, "text %s %s\n", b.tag.c_str(), b.text.c_str());
    if (!b.F.ok) return;
    fprintf(o, "gof %s nbins %d nfree %d ndf %d chi2P %.3f chi2E %.3f chi2BC %.3f chi2GLS %.3f chi2GLSdiag %.3f chi2GLSmin %.3f nrep_ok %d\n",
            b.tag.c_str(), b.F.n, b.F.nfree, b.F.ndf, b.F.chi2P, b.F.chi2E, b.F.chi2BC, b.hasGLS ? b.G.chi2 : -1., b.hasGLS ? b.G.chi2diag : -1.,
            b.hasGLS ? b.G.chi2min : -1., b.hasRep ? b.Rp.nok : 0);
    fprintf(o, "names %s", b.tag.c_str()); for (auto& n : b.S.names) fprintf(o, " %s", n.c_str()); fprintf(o, "\n");
    fprintf(o, "fixed %s", b.tag.c_str()); for (int k : b.S.fixed) fprintf(o, " %d", k); fprintf(o, "\n");
    fprintf(o, "constraints %s", b.tag.c_str()); for (auto& c : b.S.cons) fprintf(o, " %s=%.3f+-%.3f", b.S.names[c.k].c_str(), c.val, c.sig); fprintf(o, "\n");
    fprintf(o, "values %s", b.tag.c_str()); for (double x : b.F.th) fprintf(o, " %.8g", x); fprintf(o, "\n");
    const char* ek[5] = {"err_N", "err_CB", "err_stat", "err_mc", "err_total"};
    for (int w = 0; w < 5; ++w) { fprintf(o, "%s %s", ek[w], b.tag.c_str()); for (int k = 0; k < b.S.np(); ++k) fprintf(o, " %.6g", b.Err(k, w)); fprintf(o, "\n"); }
    for (size_t q = 0; q < b.S.cons.size(); ++q) {
        fprintf(o, "err_cons_%s %s", b.S.names[b.S.cons[q].k].c_str(), b.tag.c_str());
        for (int k = 0; k < b.S.np(); ++k) fprintf(o, " %.6g", std::sqrt(std::max(0., b.F.Ccons[q][k][k])));
        fprintf(o, "\n");
    }
    if (b.hasRep) { fprintf(o, "mc_shift %s", b.tag.c_str()); for (int k = 0; k < b.S.np(); ++k) fprintf(o, " %.6g", b.Rp.mean[k] - b.F.th[k]); fprintf(o, "\n"); }
    WriteMat(o, "cov_stat", b.tag, b.F.Cstat);
    if (b.hasRep) WriteMat(o, "cov_mc", b.tag, b.Rp.cov);
    WriteMat(o, "cov_total", b.tag, b.CovTot());
    for (const auto& d : b.dv) {
        double ec = 0; for (double x : d.eCons) ec += x * x;
        fprintf(o, "derived %s %-12s %.8g N %.5g CB %.5g constr %.5g stat %.5g mc %.5g total %.5g mc_shift %.5g bias_corr %.8g mc68 %.5g mc_median %.8g\n", b.tag.c_str(),
                d.name.c_str(), d.v, d.eN, d.eCB, std::sqrt(ec), d.eStat, d.eMC, d.eTot, d.shift, d.bc, d.eMC68, d.median);
    }
}

// read "values" and "err_total" of stageA.txt
bool ReadStageA(std::vector<std::string>& names, Vec& val, Vec& err) {
    std::ifstream in("/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt");
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream s(line); std::string k; s >> k;
        if (k == "names") { std::string x; names.clear(); while (s >> x) names.push_back(x); }
        else if (k == "values") { double x; val.clear(); while (s >> x) val.push_back(x); }
        else if (k == "err_total") { double x; err.clear(); while (s >> x) err.push_back(x); }
    }
    return !val.empty();
}

// profile of parameter k: fix at the grid values, refit; dropCons removes a constraint on k
struct Prof { Vec x, chi2P, chi2BC; double min = 0, minBC = 0; };
Prof Profile(const Block& b, int k, const Vec& grid, bool dropCons) {
    Prof p;
    Spec S = b.S;
    if (dropCons) { std::vector<bla::Constraint> c; for (auto& q : S.cons) if (q.k != k) c.push_back(q); S.cons = c; }
    // fixed-weight chi2 (variance frozen at the nominal fit): its minimum is the fit itself
    FitRes F0 = Fit(b.B, S, *b.D, &b.F.th, false, 300, &b.F.V);
    p.min = F0.chi2P; p.minBC = F0.chi2BC;
    Spec Sf = S; Sf.fixed[k] = 1;
    Vec st = F0.th;
    for (double v : grid) {
        Sf.start = st; Sf.start[k] = v;
        FitRes F = Fit(b.B, Sf, *b.D, &Sf.start, false, 300, &b.F.V);
        if (!F.ok) continue;
        p.x.push_back(v); p.chi2P.push_back(F.chi2P); p.chi2BC.push_back(F.chi2BC);
        st = F.th;
    }
    return p;
}

// consistency of a set of values with a common mean (GLS): chi2, ndf, p
void Consistency(const Vec& a, const Mat& C, double& mean, double& emean, double& chi2, int& ndf, double& p) {
    const int n = (int)a.size();
    Mat Ci = C;
    bla::CholInv(Ci);
    double s1 = 0, sa = 0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { s1 += Ci[i][j]; sa += Ci[i][j] * a[j]; }
    mean = sa / s1; emean = 1. / std::sqrt(s1);
    chi2 = 0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) chi2 += (a[i] - mean) * Ci[i][j] * (a[j] - mean);
    ndf = n - 1; p = TMath::Prob(chi2, ndf);
}

}  // namespace wp5m

void wp5_extract(int nrep = 100) {
    using namespace wp5;
    using namespace wp5m;
    gSystem->mkdir(kOut, kTRUE);
    gLog = fopen(Form("%s/wp5_extract.log", kOut), "w");
    FILE* res = fopen(Form("%s/results.txt", kOut), "w");
    TFile* fh = TFile::Open(kHistFile);
    P("wp5_extract.C: stage-B extraction (SPEC report_all WP5), %d simulation replicas\n", nrep);

    // ---------------------------------------------------------------- inputs
    Templ T0;
    if (!LoadTempl(fh, "S0", "nom", 1.0, T0)) return;
    std::vector<Templ> R0;
    if (!LoadReplicas(fh, "S0", 1.0, T0, nrep, R0)) return;
    DataSet D0;
    if (!LoadData(fh, "S0", "all", D0, gLog)) return;
    P("inputs: S0 data (event-group covariance %s), S0__nom templates, %zu replicas\n", D0.egFull ? "full" : "diagonal", R0.size());
    std::vector<std::string> saN; Vec saV, saE;
    ReadStageA(saN, saV, saE);
    double bA = 0, bAerr = 0;
    for (size_t k = 0; k < saN.size(); ++k) if (saN[k] == "b") { bA = saV[k]; bAerr = saE[k]; }
    P("stage A (wp3/stageA.txt): b = %.4f +- %.4f (total)\n", bA, bAerr);

    // ---------------------------------------------------------------- 0. implementation check against wp4::NLFit
    Opt oSpec; oSpec.useE = 3; oSpec.kappaCons = false;
    {
        Bins Bs; MakeBins(D0, T0, 0, 3, Bs);
        Spec Ss = SpecMain(oSpec);
        FitRes Fs = Fit(Bs, Ss, D0, nullptr, false);
        wp4::FitData FD;
        for (int i = 0; i < Bs.n; ++i) {
            double s[wp4::kNSig], vs[wp4::kNSig];
            for (int t = 0; t < kNSig; ++t) { s[t] = Bs.s[t][i]; vs[t] = Bs.vs[t][i]; }
            FD.Push(Bs.y[i], Bs.cb[i], Bs.vcb[i], Bs.x[i], Bs.j[i], Bs.e[i], Bs.gd[i], Bs.gc[i], 0., s, vs);
        }
        const wp4::FitOut f4 = wp4::NLFit(FD, wp4::SpecFMain());
        P("\n0. check of the fitter against wp4::NLFit (F-main, E1 + E2, same bins): ");
        if (Fs.ok && f4.ok) {
            const double a4 = f4.th[0], k4 = f4.th[1] / f4.th[0], r4 = f4.th[2] / f4.th[0];
            P("a_eta %.6f vs %.6f, kappa %.6f vs %.6f, rho_d %.6f vs %.6f, k_E1 %.5f vs %.5f, chi2P %.3f vs %.3f\n",
              Fs.th[Ss.Idx("a_eta")], a4, Fs.th[Ss.Idx("kappa")], k4, Fs.th[Ss.Idx("rho_d")], r4, Fs.th[Ss.Idx("k_E1")], f4.th[11], Fs.chi2P, f4.chi2P);
        } else P("one of the fits failed (%d %d)\n", Fs.ok, f4.ok);
    }

    // ---------------------------------------------------------------- 1. F-main as specified (E1 + E2, k_E1 free)
    Block bSpec = Run("FMAIN_SPEC", "F-main as specified: F10 E1 + E2, k_E1 free, k_E2 = 1 +- 0.05, g = 1 +- 0.3, kappa free",
                      D0, T0, &R0, SpecMain(oSpec), 0, 3, DerivMain(SpecMain(oSpec)), false);
    PrintBlock(bSpec);
    const double kE1 = bSpec.V("k_E1");
    const bool fbE1 = !(kE1 >= 0.8 && kE1 <= 3.0);
    P("  fall-back check 'k_E1 outside 0.8-3': k_E1 = %.3f -> %s\n", kE1, fbE1 ? "FIRED (drop E1, kappa = 1 +- 0.2)" : "not fired");

    // ---------------------------------------------------------------- 2. nominal F-main
    Opt oNom; oNom.useE = fbE1 ? 2 : 3; oNom.kappaCons = fbE1;
    Spec SNom = SpecMain(oNom);
    Block bNom = Run("FMAIN", fbE1 ? "F-main NOMINAL after the fall-back: F10 E2 only (M(ee) >= 10 MeV), kappa = 1 +- 0.2, k_E2 = 1 +- 0.05, g = 1 +- 0.3"
                                   : "F-main NOMINAL: F10 E1 + E2", D0, T0, &R0, SNom, 0, oNom.useE, DerivMain(SNom), true);
    PrintBlock(bNom);

    // chi2 decomposition of the nominal fit (where the misfit sits)
    {
        P("\n  chi2 (Pearson, V x event-group factor) of FMAIN by slice and M4 region:\n");
        const double edges[6] = {0.400, 0.455, 0.5225, 0.565, 0.575, 0.700};
        const char* rn[5] = {"0.400-0.455", "B 0.455-0.5225", "0.5225-0.565", "0.565-0.575", "CR-U 0.575-0.700"};
        double tot[5] = {0}; int nb[5] = {0};
        for (int j = 0; j < kNSlice; ++j) {
            double c[5] = {0}; int m[5] = {0};
            for (int i = 0; i < bNom.B.n; ++i) {
                if (bNom.B.j[i] != j) continue;
                const double mc = 0.5 * (bNom.B.lo[i] + bNom.B.hi[i]);
                int r = 0; while (r < 4 && mc >= edges[r + 1]) ++r;
                const double se = bNom.B.y[i] > 0 ? D0.eg10[(size_t)bNom.B.cell[i] * 700 + bNom.B.cell[i]] / bNom.B.y[i] : 1.;
                const double v = bNom.F.V[i] + (se - 1.) * std::max(bNom.F.mu[i], 1.);
                c[r] += std::pow(bNom.B.y[i] - bNom.F.mu[i], 2) / v; ++m[r];
            }
            P("    S%d:", j + 1);
            for (int r = 0; r < 5; ++r) { P("  %s %.1f/%d", rn[r], c[r], m[r]); tot[r] += c[r]; nb[r] += m[r]; }
            P("\n");
        }
        P("    all:"); for (int r = 0; r < 5; ++r) P("  %s %.1f/%d", rn[r], tot[r], nb[r]); P("\n");
    }

    // ---------------------------------------------------------------- 3. model variants (fall-back "stage B chi2 > 1.5": S3/S4)
    std::vector<Block> vars;
    auto addVar = [&](const std::string& tag, const std::string& text, Opt o, const Templ* T = nullptr, const std::vector<Templ>* R = nullptr,
                      double mlo = kF10Lo, double mhi = kF10Hi) {
        Spec S = SpecMain(o);
        vars.push_back(Run(tag, text, D0, T ? *T : T0, R, S, 0, o.useE, DerivMain(S), R != nullptr, nullptr, nullptr, mlo, mhi));
        PrintBlock(vars.back());
    };
    { Opt o = oNom; o.fixB = true; o.bVal = bA; o.bSig = bAerr; addVar("FMAIN_bA", "nominal with b constrained to stage A (fall-back form for |rho(a_eta, b)| > 0.7)", o, nullptr, &R0); }
    { Opt o = oNom; o.quad = true; addVar("FMAIN_quad", "S3: nominal + quadratic term q x^2", o, nullptr, &R0); }
    { Opt o = oNom; o.perSliceTilt = true; addVar("FMAIN_tiltj", "S3: nominal with per-slice tilts b_j", o, nullptr, &R0); }
    { Opt o = oNom; o.kappaCons = false; addVar("FMAIN_kfree", "nominal with kappa free (no 1 +- 0.2)", o, nullptr, &R0); }
    { Opt o = oNom; o.gFree = true; addVar("FMAIN_gfree", "S8: nominal with g free", o, nullptr, &R0); }
    { Opt o = oNom; addVar("FMAIN_r42", "S5: nominal on M4 0.42-0.70", o, nullptr, &R0, 0.42, 0.70); }
    { Opt o = oNom; addVar("FMAIN_r45", "S5: nominal on M4 0.45-0.70", o, nullptr, &R0, 0.45, 0.70); }
    Templ Traw;
    LoadTempl(fh, "S0", "nom", 1.0, Traw, false); RawAsSmooth(Traw);
    std::vector<Templ> Rraw;
    LoadReplicas(fh, "S0", 1.0, T0, nrep, Rraw, -1, true);
    { Opt o = oNom; addVar("FMAIN_raw", "S4: nominal with the RAW non-eta templates (BB-lite on them too: their variance is not in V)", o, &Traw, &Rraw); }
    Templ T4;
    LoadTempl(fh, "S0", "nom", 1.0, T4, false); SmoothAll(T4, 2);
    std::vector<Templ> R4; LoadReplicas(fh, "S0", 1.0, T0, nrep, R4, 2);
    { Opt o = oNom; addVar("FMAIN_pol2", "S4: nominal with exp(pol2) smoothing of the non-eta templates", o, &T4, &R4); }

    // ---------------------------------------------------------------- 4. per-slice a_eta
    Spec SPs = SpecMain(oNom, true);
    std::vector<Deriv> dPs;
    for (int j = 0; j < kNSlice; ++j) { const int k = SPs.Idx(Form("a_eta_%d", j + 1)); dPs.push_back({Form("a_eta_%d", j + 1), "per-slice eta scale", [=](const Vec& th, const Templ&) { return th[k]; }}); }
    Block bPs = Run("FMAIN_SLICE", "nominal with a_eta per M(pipi) slice (kappa, rho_d, g common)", D0, T0, &R0, SPs, 0, oNom.useE, dPs, false);
    PrintBlock(bPs);
    double psMean = 0, psErr = 0, psChi2 = 0, psP = 0, psChi2t = 0, psPt = 0, psMeant = 0, psErrt = 0; int psNdf = 0;
    double ps4P = 0, ps4Chi2 = 0; int ps4Ndf = 0;
    {
        Vec a(kNSlice); Mat Cs = Zero(kNSlice, kNSlice), Ct = Zero(kNSlice, kNSlice);
        std::vector<int> ix; for (int j = 0; j < kNSlice; ++j) ix.push_back(SPs.Idx(Form("a_eta_%d", j + 1)));
        const Mat CT = bPs.CovTot();
        for (int j = 0; j < kNSlice; ++j) { a[j] = bPs.F.th[ix[j]]; for (int l = 0; l < kNSlice; ++l) { Cs[j][l] = bPs.F.Cstat[ix[j]][ix[l]]; Ct[j][l] = CT[ix[j]][ix[l]]; } }
        Consistency(a, Cs, psMean, psErr, psChi2, psNdf, psP);
        int nd; Consistency(a, Ct, psMeant, psErrt, psChi2t, nd, psPt);
        Vec a4(a.begin(), a.begin() + 4); Mat C4 = Zero(4, 4); for (int j = 0; j < 4; ++j) for (int l = 0; l < 4; ++l) C4[j][l] = Ct[j][l];
        double m4, e4; Consistency(a4, C4, m4, e4, ps4Chi2, ps4Ndf, ps4P);
        P("  per-slice a_eta: weighted mean %.4f +- %.4f (stat), chi2 %.2f/%d p = %.3g (stat only); with stat + MC: mean %.4f +- %.4f, chi2 %.2f/%d p = %.3g; "
          "S1-S4 only (stat + MC): chi2 %.2f/%d p = %.3g\n", psMean, psErr, psChi2, psNdf, psP, psMeant, psErrt, psChi2t, psNdf, psPt, ps4Chi2, ps4Ndf, ps4P);
    }

    // ---------------------------------------------------------------- 5. F-free (E2: nominal classes; E1 + E2: as specified)
    Spec SFr = SpecFree(oNom);
    Block bFree = Run("FFREE", "F-free on F10 E2 (nominal classes): alpha_2..alpha_5 free, G_eta x alpha_2 g", D0, T0, &R0, SFr, 0, oNom.useE, DerivFree(SFr), false);
    PrintBlock(bFree);
    Spec SFr12 = SpecFree(oSpec);
    Block bFree12 = Run("FFREE_SPEC", "F-free as specified on F10 E1 + E2 (k_E1 free)", D0, T0, &R0, SFr12, 0, 3, DerivFree(SFr12), false);
    PrintBlock(bFree12);

    // ---------------------------------------------------------------- 6. F-gold (F5g), C3 + C4 and kappa_B fixed from F-main
    const int nA = SNom.Idx("a_eta"), nK = SNom.Idx("kappa"), nKB = SNom.Idx("kappa_B");
    auto goldFix = [&](const Vec& th, const Templ& T) {
        return [&th, &T, nA, nK](int h, int x) { return th[nA] * th[nK] * (T.s[kS3][h][x] + T.s[kS4][h][x]); };
    };
    Spec SG = SpecGold(oNom, bNom.F.th[nKB]);
    Block bGold = Run("FGOLD", "F-gold on F5g (M(ee) >= 20 MeV): alpha_2, alpha_5, g, c_j, b, k_E2; kappa_B and a kappa (C3 + C4) fixed from F-main",
                      D0, T0, &R0, SG, 1, 2, DerivGold(SG), false, goldFix(bNom.F.th, T0),
                      [&](Bins& Br, const Templ& Tr, int r) {
                          const Vec& thr = bNom.Rp.fit[r].ok ? bNom.Rp.fit[r].th : bNom.F.th;
                          Spec Sr = SpecGold(oNom, thr[nKB]);
                          (void)Sr;
                          MakeBins(D0, Tr, 1, 2, Br, kF10Lo, kF10Hi, goldFix(thr, Tr));
                      });
    PrintBlock(bGold);
    // the replica fits above keep kappa_B of the nominal (Spec fixed); add the kappa_B input error by +-1 sigma (total) refits
    double goldKBerr = 0;
    {
        const double sKB = std::sqrt(bNom.CovTot()[nKB][nKB]);
        double brv[2];
        for (int s = 0; s < 2; ++s) {
            Spec Sx = SpecGold(oNom, bNom.F.th[nKB] + (s == 0 ? sKB : -sKB));
            Bins Bx; MakeBins(D0, T0, 1, 2, Bx, kF10Lo, kF10Hi, goldFix(bNom.F.th, T0));
            FitRes Fx = Fit(Bx, Sx, D0, &bGold.F.th, false);
            brv[s] = Fx.ok ? kBRdir * Fx.th[Sx.Idx("alpha_5")] / Fx.th[Sx.Idx("alpha_2")] : NAN;
        }
        goldKBerr = 0.5 * std::fabs(brv[0] - brv[1]);
        P("  F-gold: BR_G with kappa_B +- 1 sigma (%.3f): %.4g / %.4g -> input error %.3g\n", sKB, brv[0], brv[1], goldKBerr);
    }

    // ---------------------------------------------------------------- 7. mass scale: s scan (templates at s = 0.99521, 0.99525, 0.99615, 0.99705) and r
    struct SP { std::string var; double s, r; FitRes F; };
    std::vector<SP> sscan = {{"sK0", kSK0, 1, FitRes()}, {"sDn", kS - kSerr, 1, FitRes()}, {"nom", kS, 1, FitRes()}, {"sUp", kS + kSerr, 1, FitRes()},
                             {"r089", kS, kRlo, FitRes()}, {"r110", kS, kRhi, FitRes()}};
    P("\n--- mass-transfer scan (nominal spec, templates of the WP0 variants; chi2P = fixed-weight chi2 with the variance of FMAIN)\n");
    for (auto& x : sscan) {
        Templ Tv;
        if (x.var == "nom") { }
        else {
            LoadTempl(fh, "S0", x.var.c_str(), 1.0, Tv);
            Bins Bv; MakeBins(D0, Tv, 0, oNom.useE, Bv);
            x.F = Fit(Bv, SNom, D0, &bNom.F.th, false, 300, &bNom.F.V);
        }
        if (x.var == "nom") x.F = Fit(bNom.B, SNom, D0, &bNom.F.th, false, 300, &bNom.F.V);
        P("    %-5s s = %.5f r = %.2f: chi2P %.2f chi2BC %.2f a_eta %.4f rho_d %.4f kappa %.4f g %.4f\n", x.var.c_str(), x.s, x.r, x.F.chi2P, x.F.chi2BC,
          x.F.th[nA], x.F.th[SNom.Idx("rho_d")], x.F.th[nK], x.F.th[SNom.Idx("g")]);
    }
    // parabola through sDn, nom, sUp (chi2P), s_hat and sigma (delta chi2 = 1)
    double sHat = NAN, sSig = NAN, sSigScaled = NAN;
    {
        const double x0 = sscan[1].s, x1 = sscan[2].s, x2 = sscan[3].s;
        const double y0 = sscan[1].F.chi2P, y1 = sscan[2].F.chi2P, y2 = sscan[3].F.chi2P;
        const double h = x1 - x0;
        const double A = (y0 - 2 * y1 + y2) / (2 * h * h), Bc = (y2 - y0) / (2 * h);
        if (A > 0) { sHat = x1 - Bc / (2 * A); sSig = 1. / std::sqrt(A); }
        sSigScaled = sSig * std::sqrt(std::max(1., bNom.F.chi2E / bNom.F.ndf));
        P("    free s from the chi2 parabola (sDn, nom, sUp): s_hat = %.5f +- %.5f (delta chi2 = 1), +- %.5f (x sqrt(chi2/ndf)); nominal 0.99615: %.2f sigma (%.2f scaled)\n",
          sHat, sSig, sSigScaled, (sHat - kS) / sSig, (sHat - kS) / sSigScaled);
    }

    // ---------------------------------------------------------------- 8. folds
    std::vector<Block> folds;
    std::vector<Templ> Tf(2);
    std::vector<DataSet> Df(2);
    for (int q = 0; q < 2; ++q) {
        const char* sub = q == 0 ? "fold0" : "fold1";
        LoadData(fh, "S0", sub, Df[q], gLog);
        const double f = Df[q].Nwin[kYTot] / D0.Nwin[kYTot];
        LoadTempl(fh, "S0", "nom", f, Tf[q]);
        folds.push_back(Run(Form("FMAIN_%s", sub), Form("nominal on data %s (event-group parity), templates x %.5f (row fraction)", sub, f), Df[q], Tf[q], nullptr, SNom, 0, oNom.useE, DerivMain(SNom), false));
        PrintBlock(folds.back());
    }
    double foldChi2 = 0, foldP = 0; int foldNdf = 0;
    std::map<std::string, double> foldZ;
    {
        // joint test of the free parameters with the data parts of the two sandwich covariances (constraints are common priors)
        std::vector<int> fr; for (int k = 0; k < SNom.np(); ++k) if (!SNom.fixed[k]) fr.push_back(k);
        const int nf = (int)fr.size();
        Mat C = Zero(nf, nf); Vec d(nf);
        for (int a = 0; a < nf; ++a) {
            d[a] = folds[0].F.th[fr[a]] - folds[1].F.th[fr[a]];
            for (int b = 0; b < nf; ++b) C[a][b] = folds[0].F.Cn[fr[a]][fr[b]] + folds[0].F.Ccb[fr[a]][fr[b]] + folds[1].F.Cn[fr[a]][fr[b]] + folds[1].F.Ccb[fr[a]][fr[b]];
        }
        Mat Ci = C; bla::CholInv(Ci);
        for (int a = 0; a < nf; ++a) for (int b = 0; b < nf; ++b) foldChi2 += d[a] * Ci[a][b] * d[b];
        foldNdf = nf; foldP = TMath::Prob(foldChi2, foldNdf);
        P("  folds: joint chi2 of the parameter differences (data covariance) %.2f/%d, p = %.3g\n", foldChi2, foldNdf, foldP);
        for (const char* q : {"a_eta", "rho_d", "kappa", "g", "Y_bump", "Y_peak", "DR"}) {
            const DVal *u = folds[0].Get(q), *v = folds[1].Get(q);
            if (!u || !v) continue;
            const double e2 = u->eN * u->eN + u->eCB * u->eCB + v->eN * v->eN + v->eCB * v->eCB;
            foldZ[q] = (u->v - v->v) / std::sqrt(e2);
            P("    %-7s fold0 %.4g fold1 %.4g: z = %+.2f (p = %.3g)\n", q, u->v, v->v, foldZ[q], TMath::Prob(foldZ[q] * foldZ[q], 1));
        }
    }

    // ---------------------------------------------------------------- 9. pre (parallel column)
    Templ Tpre; LoadTempl(fh, "pre", "nom", 1.0, Tpre);
    std::vector<Templ> Rpre; LoadReplicas(fh, "pre", 1.0, Tpre, nrep, Rpre);
    DataSet Dpre; LoadData(fh, "pre", "all", Dpre, gLog);
    Block bPre = Run("FMAIN_pre", "nominal spec on selection pre (no pion-vertex requirement), pre templates and replicas", Dpre, Tpre, &Rpre, SNom, 0, oNom.useE, DerivMain(SNom), true);
    PrintBlock(bPre);

    // ---------------------------------------------------------------- 10. data isBest == 1 (simulation templates isBest >= 0, D2)
    DataSet Db; LoadData(fh, "S0_best", "all", Db, gLog);
    Block bBest = Run("FMAIN_best", "nominal spec on data S0 with isBest == 1 (templates: simulation isBest >= 0, D2); data-only yield ratios", Db, T0, &R0, SNom, 0, oNom.useE, DerivMain(SNom), false);
    PrintBlock(bBest);
    std::map<std::string, std::array<double, 3>> bestRatio;   // value, stat (independent-sample bound), MC (replica by replica)
    for (const char* q : {"a_eta", "rho_d", "Y_bump", "Y_peak", "Y2_B", "Y5_P", "DR"}) {
        const DVal *u = bBest.Get(q), *v = bNom.Get(q);
        if (!u || !v) continue;
        const double r = u->v / v->v;
        const double es = r * std::sqrt(std::pow(u->eStat / u->v, 2) + std::pow(v->eStat / v->v, 2));
        double emc = 0;
        if (u->rep.size() == v->rep.size() && u->rep.size() > 2) {
            double s = 0, s2 = 0; const int n = (int)u->rep.size();
            for (int k = 0; k < n; ++k) { const double x = u->rep[k] / v->rep[k]; s += x; s2 += x * x; }
            emc = std::sqrt(std::max(0., (s2 - s * s / n) / (n - 1)));
        }
        bestRatio[q] = {r, es, emc};
        P("  isBest==1 / isBest>=0: %-7s %.4f +- %.4f (stat, independent-sample bound) +- %.4f (MC)\n", q, r, es, emc);
    }

    // ---------------------------------------------------------------- 11. periods (input of the WP2 stability test): background shape fixed
    std::vector<Block> pers;
    double lsum = 0; for (int p = 0; p < kNPer; ++p) lsum += kLper[p];
    for (int p = 0; p < kNPer; ++p) {
        const char* sub = SubTag(kSubPerI + p);
        DataSet Dp; LoadData(fh, "S0", sub, Dp, gLog);
        const double f = kLper[p] / lsum;
        Templ Tp; LoadTempl(fh, "S0", "nom", f, Tp);
        // background of the nominal fit, scaled by f, baked into G_D; one free overall scale c_bg
        Spec S;
        const int ia = S.Add("a_eta", 0.4), ik = S.Add("kappa", bNom.F.th[nK]), ir = S.Add("rho_d", 1.3), ig = S.Add("g", bNom.F.th[SNom.Idx("g")]);
        S.fixed[ik] = 1; S.fixed[ig] = 1;
        S.terms = {{kS2, ia, -1, -1}, {kS3, ia, ik, -1}, {kS4, ia, ik, -1}, {kS5, ia, ir, -1}, {kSeta, ia, ig, -1}};
        const int ic = S.Add("c_bg", 1.0);
        for (int j = 0; j < kNSlice; ++j) S.iC[j] = ic;
        S.iB = S.Add("b", 0.); S.fixed[S.iB] = 1;
        S.iKB = S.Add("kappa_B", 0.); S.fixed[S.iKB] = 1;
        S.iK[1] = S.Add("k_E2", 1.0); S.cons.push_back({S.iK[1], 1.0, 0.05});
        Block b; b.tag = Form("PERIOD_%s", PeriodName(p)); b.text = Form("period %s (L share %.4f): a_eta, rho_d, overall background scale free; background shape, kappa, g fixed to FMAIN", PeriodName(p), f);
        b.S = S; b.T = &Tp; b.D = &Dp;
        MakeBins(Dp, Tp, 0, oNom.useE, b.B);
        const Vec& t = bNom.F.th;
        for (int i = 0; i < b.B.n; ++i) {
            const int j = b.B.j[i];
            b.B.gd[i] = t[SNom.iC[j]] * (1. + t[SNom.iB] * b.B.x[i]) * (b.B.gd[i] + t[nKB] * b.B.gc[i]);
            b.B.gc[i] = 0.;
        }
        b.F = Fit(b.B, S, Dp);
        for (const auto& d : DerivMain(S)) b.dv.push_back(Evaluate(d, b.F, Tp, nullptr, nullptr));
        PrintBlock(b);
        pers.push_back(b);
    }
    double perChi2 = 0, perP = 0, perMean = 0, perErr = 0; int perNdf = 0;
    double perChi2DR = 0, perPDR = 0, perMeanDR = 0, perErrDR = 0;
    {
        Vec a, dr; Mat C = Zero(kNPer, kNPer), Cd = Zero(kNPer, kNPer);
        for (int p = 0; p < kNPer; ++p) { const DVal* u = pers[p].Get("a_eta"); const DVal* w = pers[p].Get("DR"); a.push_back(u->v); dr.push_back(w->v); C[p][p] = u->eStat * u->eStat; Cd[p][p] = w->eStat * w->eStat; }
        Consistency(a, C, perMean, perErr, perChi2, perNdf, perP);
        int nd; Consistency(dr, Cd, perMeanDR, perErrDR, perChi2DR, nd, perPDR);
        P("  periods: a_eta chi2 %.2f/%d p = %.3g (mean %.4f); DR chi2 %.2f/%d p = %.3g (mean %.4f)\n", perChi2, perNdf, perP, perMean, perChi2DR, nd, perPDR, perMeanDR);
    }

    // ---------------------------------------------------------------- 12. local free s: peak region only (E2, M4 0.5225-0.565)
    const int nG = SNom.Idx("g"), nR = SNom.Idx("rho_d"), nB = SNom.iB;
    double sLoc = NAN, sLocSig = NAN, sLocScaled = NAN, locChi2Ndf = 1;
    std::vector<std::array<double, 3>> sLocPts;   // s, chi2P, a_eta
    {
        Opt o = oNom; o.fixB = true; o.bVal = bNom.F.th[nB]; o.bSig = 0; o.fixKB = true; o.kBval = bNom.F.th[nKB]; o.kappaCons = false;
        Spec SL = SpecMain(o);
        for (const char* q : {"a_eta", "kappa", "g"}) { const int k = SL.Idx(q); SL.fixed[k] = 1; SL.start[k] = bNom.F.th[SNom.Idx(q)]; }
        { std::vector<bla::Constraint> c; for (auto& q : SL.cons) if (!SL.fixed[q.k]) c.push_back(q); SL.cons = c; }
        Vec st = SL.start; for (int k = 0; k < SL.np(); ++k) { const int q = SNom.Idx(SL.names[k]); if (q >= 0) st[k] = bNom.F.th[q]; }
        Bins Bn; MakeBins(D0, T0, 0, oNom.useE, Bn, 0.5225, 0.565);
        const FitRes Fn = Fit(Bn, SL, D0, &st, false);
        P("\n--- local free s: peak region M4 0.5225-0.565, E2; a_eta, kappa, g, b, kappa_B fixed to FMAIN; rho_d, c_j, k_E2 free; fixed-weight chi2 (variance of the s = 0.99615 fit)\n");
        for (auto& x : sscan) {
            if (x.r != 1) continue;
            Templ Tv;
            if (x.var == "nom") Tv = T0; else LoadTempl(fh, "S0", x.var.c_str(), 1.0, Tv);
            Bins Bv; MakeBins(D0, Tv, 0, oNom.useE, Bv, 0.5225, 0.565);
            FitRes F = Fit(Bv, SL, D0, &Fn.th, false, 300, &Fn.V);
            P("    %-5s s = %.5f: chi2P %.2f / %d, a_eta %.4f rho_d %.4f\n", x.var.c_str(), x.s, F.chi2P, F.ndf, F.th[SL.Idx("a_eta")], F.th[SL.Idx("rho_d")]);
            sLocPts.push_back({x.s, F.chi2P, F.th[SL.Idx("a_eta")]});
            if (x.var == "nom") locChi2Ndf = Fn.chi2E / Fn.ndf;
        }
        // parabola through sDn, nom, sUp (entries 1, 2, 3)
        const double x0 = sLocPts[1][0], x1 = sLocPts[2][0], h = x1 - x0;
        const double y0 = sLocPts[1][1], y1 = sLocPts[2][1], y2 = sLocPts[3][1];
        const double A = (y0 - 2 * y1 + y2) / (2 * h * h), Bc = (y2 - y0) / (2 * h);
        if (A > 0) { sLoc = x1 - Bc / (2 * A); sLocSig = 1. / std::sqrt(A); }
        sLocScaled = sLocSig * std::sqrt(std::max(1., locChi2Ndf));
        P("    local s_hat = %.5f +- %.5f (delta chi2 = 1) +- %.5f (x sqrt(chi2/ndf) = %.2f); pull vs 0.99615: %.2f (%.2f scaled)\n", sLoc, sLocSig, sLocScaled,
          std::sqrt(std::max(1., locChi2Ndf)), (sLoc - kS) / sLocSig, (sLoc - kS) / sLocScaled);
    }

    // ---------------------------------------------------------------- 13. profiles (Fig. E3)
    Prof prA, prR, prK, prG;
    {
        const double ea = std::sqrt(bNom.F.C[nA][nA]), er = std::sqrt(bNom.F.C[nR][nR]);
        Vec ga, gr, gk, gg;
        for (int i = -24; i <= 24; ++i) { ga.push_back(bNom.F.th[nA] + i * ea / 6.); gr.push_back(bNom.F.th[nR] + i * er / 6.); }
        for (int i = 0; i <= 40; ++i) { gk.push_back(-0.6 + i * 0.06); gg.push_back(0.4 + i * 0.08); }
        prA = Profile(bNom, nA, ga, false);
        prR = Profile(bNom, nR, gr, false);
        prK = Profile(bNom, nK, gk, true);
        prG = Profile(bNom, nG, gg, true);
    }

    // ---------------------------------------------------------------- 14. acceptance (spec WP5)
    const Mat CT = bNom.CovTot(), CS = bNom.F.Cstat;
    const double rhoAbS = CS[nA][nB] / std::sqrt(CS[nA][nA] * CS[nB][nB]), rhoAbT = CT[nA][nB] / std::sqrt(CT[nA][nA] * CT[nB][nB]);
    // corr(rho_d, kappa)
    const double rhoRkS = CS[nR][nK] / std::sqrt(CS[nR][nR] * CS[nK][nK]), rhoRkT = CT[nR][nK] / std::sqrt(CT[nR][nR] * CT[nK][nK]);
    const double rhoRkSpec = bSpec.CovTot()[bSpec.S.Idx("rho_d")][bSpec.S.Idx("kappa")] /
                             std::sqrt(bSpec.CovTot()[bSpec.S.Idx("rho_d")][bSpec.S.Idx("rho_d")] * bSpec.CovTot()[bSpec.S.Idx("kappa")][bSpec.S.Idx("kappa")]);
    // pull runs in 0.52-0.56 (E2 histograms of the nominal fit)
    int maxRun = 0, maxRun1 = 0, nHistRun = 0, sumRun = 0, sumRun1 = 0;
    std::vector<std::array<int, 2>> runsJ(kNSlice);
    {
        std::vector<Vec> pj(kNSlice);
        Vec sumR(16, 0.), sumV(16, 0.);
        for (int i = 0; i < bNom.B.n; ++i) {
            const double c = 0.5 * (bNom.B.lo[i] + bNom.B.hi[i]);
            if (c < 0.52 || c > 0.56) continue;
            const double se = bNom.B.y[i] > 0 ? D0.eg10[(size_t)bNom.B.cell[i] * 700 + bNom.B.cell[i]] / bNom.B.y[i] : 1.;
            const double v = bNom.F.V[i] + (se - 1.) * std::max(bNom.F.mu[i], 1.);
            pj[bNom.B.j[i]].push_back((bNom.B.y[i] - bNom.F.mu[i]) / std::sqrt(v));
            const int k = (int)std::lround((c - 0.52125) / 0.0025);
            if (k >= 0 && k < 16) { sumR[k] += bNom.B.y[i] - bNom.F.mu[i]; sumV[k] += v; }
        }
        for (int j = 0; j < kNSlice; ++j) {
            int L, L1; Runs(pj[j], L, L1);
            runsJ[j] = {L, L1};
            maxRun = std::max(maxRun, L); maxRun1 = std::max(maxRun1, L1);
            if (L > 4) ++nHistRun;
        }
        Vec ps(16); for (int k = 0; k < 16; ++k) ps[k] = sumV[k] > 0 ? sumR[k] / std::sqrt(sumV[k]) : 0.;
        Runs(ps, sumRun, sumRun1);
        P("\n  pull runs in 0.52-0.56 (E2, 16 bins of 2.5 MeV): longest same-sign run per slice:");
        for (int j = 0; j < kNSlice; ++j) P(" S%d %d (|pull|>1: %d)", j + 1, runsJ[j][0], runsJ[j][1]);
        P("; summed over slices: %d (|pull|>1: %d)\n", sumRun, sumRun1);
        P("    pulls of the slice sum:"); for (double x : ps) P(" %+.2f", x); P("\n");
        const double p5 = ProbRunGE(16, 5);
        P("    random signs: P(longest run >= 5 in 16 bins) = %.3f per histogram, %.3f for at least one of 5\n", p5, 1 - std::pow(1 - p5, 5));
    }
    const double chi2gls = bNom.G.chi2 / bNom.G.ndf, chi2p = bNom.F.chi2E / bNom.F.ndf;
    const double errScale = std::sqrt(std::max(1., chi2gls));
    struct Acc { std::string crit, value, verdict; };
    std::vector<Acc> acc;
    acc.push_back({"chi2/ndf <= 1.2", Form("GLS (data + MC covariance) %.1f/%d = %.2f; minimised %.2f; Pearson %.1f/%d = %.2f", bNom.G.chi2, bNom.G.ndf, chi2gls,
                                           bNom.G.chi2min / bNom.G.ndf, bNom.F.chi2E, bNom.F.ndf, chi2p),
                   Form("FAIL -> total errors x sqrt(chi2_GLS/ndf) = x %.2f (scaled columns); fall-back 'stage B chi2/ndf > 1.5' FIRED", errScale)});
    acc.push_back({"|rho(a_eta, b)| < 0.7", Form("stat %.3f, stat+MC %.3f", rhoAbS, rhoAbT), (std::fabs(rhoAbT) < 0.7 && std::fabs(rhoAbS) < 0.7) ? "PASS" : "FAIL"});
    acc.push_back({"per-slice a_eta p > 0.01", Form("stat+MC p = %.3f (chi2 %.2f/%d); stat only p = %.3f (chi2 %.2f/%d)", psPt, psChi2t, psNdf, psP, psChi2, psNdf),
                   (psPt > 0.01 && psP > 0.01) ? "PASS" : (psPt > 0.01 ? "PASS (stat+MC), FAIL (stat only)" : "FAIL")});
    acc.push_back({"no pull run > 4 bins in 0.52-0.56", Form("longest same-sign run per slice %d %d %d %d %d, slice sum %d; with |pull| > 1: max %d, sum %d",
                                                             runsJ[0][0], runsJ[1][0], runsJ[2][0], runsJ[3][0], runsJ[4][0], sumRun, maxRun1, sumRun1),
                   (maxRun <= 4 && sumRun <= 4) ? "PASS" : Form("FAIL (%d of 5 slices%s)", nHistRun, sumRun > 4 ? ", and the slice sum" : "")});
    acc.push_back({"free s within 2 sigma of 0.99615", Form("peak region: s = %.5f +- %.5f (%.1f sigma; %.1f sigma with sqrt(chi2/ndf)); whole fit: %.5f +- %.5f (%.1f sigma)",
                                                             sLoc, sLocSig, (sLoc - kS) / sLocSig, (sLoc - kS) / sLocScaled, sHat, sSig, (sHat - kS) / sSig),
                   std::fabs((sLoc - kS) / sLocSig) < 2 ? "PASS (peak region)" : (std::fabs((sLoc - kS) / sLocScaled) < 2 ? "PASS only with the chi2-scaled error" : "FAIL")});
    acc.push_back({"folds p > 0.05", Form("joint chi2 %.2f/%d p = %.3f; a_eta z = %+.2f", foldChi2, foldNdf, foldP, foldZ["a_eta"]), foldP > 0.05 ? "PASS" : "FAIL"});
    P("\n=== WP5 acceptance (nominal FMAIN) ===\n");
    for (auto& a : acc) P("  %-36s | %s | %s\n", a.crit.c_str(), a.value.c_str(), a.verdict.c_str());
    // fall-backs
    std::vector<std::string> fbl;
    fbl.push_back(Form("k_E1 outside 0.8-3 (F-main E1+E2: k_E1 = %.2f +- %.2f): %s", kE1, bSpec.Err(bSpec.S.Idx("k_E1"), 4), fbE1 ? "FIRED -> E1 dropped, kappa = 1 +- 0.2 (nominal FMAIN)" : "not fired"));
    {
        double best = 1e9; std::string bt;
        for (auto& v : vars) if (v.hasGLS && v.G.chi2 / v.G.ndf < best) { best = v.G.chi2 / v.G.ndf; bt = v.tag; }
        fbl.push_back(Form("stage B chi2/ndf > 1.5 (GLS %.2f): FIRED -> S3/S4 variants run; best GLS chi2/ndf %.2f (%s); unresolved -> the side-band estimators are the main result (already required by WP4)", chi2gls, best, bt.c_str()));
    }
    fbl.push_back(Form("|rho(a_eta, b)| > 0.7: %s (%.2f)", std::fabs(rhoAbT) > 0.7 ? "FIRED" : "not fired", rhoAbT));
    fbl.push_back(Form("per-slice a_eta p < 0.01: %s (p = %.3f with MC, %.3f stat only)", (psPt < 0.01) ? "FIRED" : "not fired", psPt, psP));
    fbl.push_back(Form("corr(rho_d, kappa) > 0.9: not fired (nominal %.2f stat+MC, %.2f stat; F-main E1+E2 %.2f); direct-decay bias > 20 %%: WP4 RF +9.7 %% in (a), +22.6 %% only in the stress scenario (e) -> not fired", rhoRkT, rhoRkS, rhoRkSpec));
    {
        const DVal* r34 = bFree.Get("R_34"); const DVal* r34s = bFree12.Get("R_34");
        const double e34 = r34 ? std::sqrt(r34->eStat * r34->eStat + r34->eMC68 * r34->eMC68) : 0, e34s = r34s ? std::sqrt(r34s->eStat * r34s->eStat + r34s->eMC68 * r34s->eMC68) : 0;
        fbl.push_back(Form("R_34 off by > 2 sigma: F-free E2 R_34 = %.2f +- %.2f (alpha_3, alpha_4 not determined without E1); F-free E1+E2 R_34 = %.3f +- %.3f (%.1f sigma) -> BR_C WITHDRAWN in the nominal (E2) configuration; the E1+E2 value is a labelled cross-check",
                        r34 ? r34->v : 0, e34, r34s ? r34s->v : 0, e34s, r34s ? (r34s->v - 1) / e34s : 0));
    }
    fbl.push_back(Form("MC error > data statistical error: FIRED (a_eta MC/stat = %.2f, Y_bump %.2f, Y_peak %.2f) -> simulation-limited", bNom.Err(nA, 3) / bNom.Err(nA, 2),
                       bNom.Get("Y_bump")->eMC / bNom.Get("Y_bump")->eStat, bNom.Get("Y_peak")->eMC / bNom.Get("Y_peak")->eStat));
    fbl.push_back("WP4 gate: F-main, F-free, F-gold FAILED -> every number of this package is a CROSS-CHECK; B4 + peak side bands (WP6) are the main result");
    P("\n=== fall-backs ===\n"); for (auto& f : fbl) P("  - %s\n", f.c_str());

    // ---------------------------------------------------------------- 15. results.txt
    fclose(res);
    res = fopen(Form("%s/results.txt", kOut), "w");
    fprintf(res, "# =====================================================================================================\n");
    fprintf(res, "# WP5 stage-B extraction (SPEC report_all) -- plots/all_report/wp5_extract.C (+ wp5_lib.h, wp5_rows.C), 2026-11-14\n");
    fprintf(res, "# inputs: hists_v1.root (WP0), wp3/stageA.txt, wp4/verdict.txt, event-group covariances wp3_rows.root / wp5/wp5_egcov.root\n");
    fprintf(res, "# STATUS: WP4 gate FAILED for F-main, F-free, F-gold (template attenuation bias, simulation-limited) -> ALL fits here are\n");
    fprintf(res, "#         CROSS-CHECKS. Main results: B4 + peak side bands (WP6). The F-main parameters below are the inputs of the\n");
    fprintf(res, "#         'other eta' terms of WP6, of the X17 constraints and of WP7/WP8 ratios.\n");
    fprintf(res, "# Units: raw data rows (yields x 64/63 already in Y*); sim templates = nb x %.4f (L = %.1f nb^-1 x 63/64)\n", kNbToRaw, kL);
    fprintf(res, "# Errors: N = sandwich, event-group covariance of N; CB = sandwich, CB replica covariance (200); err_cons_<p> = part of the\n");
    fprintf(res, "#   Gaussian constraint on <p>; stat = N + CB + constraints; mc = covariance over %d WP0 simulation replicas (refitted,\n", nrep);
    fprintf(res, "#   smoothing redone); total = stat + mc. mc_shift = replica mean - nominal (template attenuation, WP4). bias_corr = 2 v - mean(rep)\n");
    fprintf(res, "#   (WP4 'RB', information only). mc68 = half of the central 68 %% replica interval (robust; use it for ratios).\n");
    fprintf(res, "# Fit model: mu = k_e CB + c_j (1 + b x)[G_D + kappa_B G_C](smoothed exp(pol3), WP4 rule) + a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta],\n");
    fprintf(res, "#   x = (M4 - 0.55)/0.1 at the F10 bin centre; F10 M4 0.400-0.700; slices S1..S5; E1 < 10 MeV, E2 >= 10 MeV\n");
    fprintf(res, "# =====================================================================================================\n");
    fprintf(res, "status FMAIN CROSS_CHECK (WP4 FAIL)\nstatus FFREE CROSS_CHECK (WP4 FAIL)\nstatus FGOLD CROSS_CHECK (WP4 FAIL; FAIL_MARGINAL with bias correction)\n");
    fprintf(res, "nominal FMAIN\nnominal_classes %s\n", fbE1 ? "E2" : "E1 E2");
    fprintf(res, "error_scale_chi2 %.4f  # sqrt(chi2_GLS/ndf) of FMAIN, to multiply the total errors (acceptance chi2/ndf <= 1.2 failed)\n", errScale);
    fprintf(res, "simulation_limited 1\nbump_label EXTRAPOLATED (stage A: tilt from CR-U only)\n");
    for (auto& f : fbl) fprintf(res, "fallback %s\n", f.c_str());
    for (auto& a : acc) fprintf(res, "acceptance %s | %s | %s\n", a.crit.c_str(), a.value.c_str(), a.verdict.c_str());
    fprintf(res, "rho_a_b stat %.4f total %.4f\nrho_rho_kappa stat %.4f total %.4f spec_total %.4f\n", rhoAbS, rhoAbT, rhoRkS, rhoRkT, rhoRkSpec);
    fprintf(res, "perslice_p stat %.4g total %.4g chi2_stat %.3f chi2_total %.3f ndf %d mean_stat %.5f %.5f mean_total %.5f %.5f s1s4_total_p %.4g\n", psP, psPt, psChi2, psChi2t, psNdf, psMean, psErr, psMeant, psErrt, ps4P);
    fprintf(res, "folds_p %.4g chi2 %.3f ndf %d\n", foldP, foldChi2, foldNdf);
    for (auto& z : foldZ) fprintf(res, "folds_z %s %.4f\n", z.first.c_str(), z.second);
    fprintf(res, "periods_aeta chi2 %.3f ndf %d p %.4g mean %.5f err %.5f\nperiods_DR chi2 %.3f ndf %d p %.4g mean %.5f err %.5f\n", perChi2, perNdf, perP, perMean, perErr,
            perChi2DR, perNdf, perPDR, perMeanDR, perErrDR);
    fprintf(res, "free_s_global s_hat %.6f sigma %.6f sigma_scaled %.6f  # parabola through sDn/nom/sUp (fixed-weight chi2)\n", sHat, sSig, sSigScaled);
    fprintf(res, "free_s_peak s_hat %.6f sigma %.6f sigma_scaled %.6f\n", sLoc, sLocSig, sLocScaled);
    for (auto& x : sscan) fprintf(res, "sscan %s s %.5f r %.2f chi2P %.3f chi2BC %.3f a_eta %.5f rho_d %.5f kappa %.5f g %.5f\n", x.var.c_str(), x.s, x.r, x.F.chi2P, x.F.chi2BC,
                                  x.F.th[nA], x.F.th[nR], x.F.th[nK], x.F.th[nG]);
    for (auto& br : bestRatio) fprintf(res, "best_ratio %s %.5f stat_bound %.5f mc %.5f\n", br.first.c_str(), br.second[0], br.second[1], br.second[2]);
    fprintf(res, "gold_kappaB_input_err BR_G %.4g\n", goldKBerr);
    WriteBlock(res, bNom);
    WriteBlock(res, bSpec);
    for (auto& v : vars) WriteBlock(res, v);
    WriteBlock(res, bPs);
    WriteBlock(res, bFree);
    WriteBlock(res, bFree12);
    WriteBlock(res, bGold);
    for (auto& f : folds) WriteBlock(res, f);
    WriteBlock(res, bPre);
    WriteBlock(res, bBest);
    for (auto& p : pers) WriteBlock(res, p);
    // per-bin model of FMAIN (all 10 F10 histograms; E1 = prediction with k_E1 = 1 when E1 is dropped)
    Bins Ball; MakeBins(D0, T0, 0, 3, Ball);
    Spec SAll = SNom;   // E1 bins use k = kFix = 1 when k_E1 is not a parameter
    Vec pcb, pbg, psig[kNSig];
    Parts(Ball, SAll, bNom.F.th, pcb, pbg, psig);
    Vec muAll; Eval(Ball, SAll, bNom.F.th, muAll, nullptr, nullptr, nullptr);
    // background errors per bin: stat (gradient) and MC (replicas)
    Vec bgStat(Ball.n, 0.), bgMC(Ball.n, 0.);
    {
        auto bgOf = [&](const Vec& th, int i, const Bins& B) {
            const int j = B.j[i];
            return th[SNom.iC[j]] * (1. + th[SNom.iB] * B.x[i]) * (B.gd[i] + th[SNom.iKB] * B.gc[i]);
        };
        for (int i = 0; i < Ball.n; ++i) {
            Vec g = NumGrad([&](const Vec& th) { return bgOf(th, i, Ball); }, bNom.F.th, bNom.F.C);
            bgStat[i] = std::sqrt(std::max(0., Quad(bNom.F.Cstat, g)));
        }
        std::vector<Vec> rv(Ball.n);
        for (size_t r = 0; r < R0.size(); ++r) {
            if (!bNom.Rp.fit[r].ok) continue;
            Bins Br = Ball; SwapTempl(Br, R0[r]);
            for (int i = 0; i < Ball.n; ++i) rv[i].push_back(bgOf(bNom.Rp.fit[r].th, i, Br));
        }
        for (int i = 0; i < Ball.n; ++i) {
            double s = 0, q = 0; const int n = (int)rv[i].size();
            for (double x : rv[i]) { s += x; q += x * x; }
            bgMC[i] = n > 1 ? std::sqrt(std::max(0., (q - s * s / n) / (n - 1))) : 0.;
        }
    }
    fprintf(res, "\n# per-bin FMAIN model (raw rows): hist bin M4lo M4hi N CB kCB bkg bkg_err_stat bkg_err_mc C2 C3 C4 C5 Geta mu in_fit\n");
    fprintf(res, "#   (C* = fitted amplitude x template; E1 rows: prediction from the E2 fit with k_E1 = 1 if E1 is dropped)\n");
    for (int i = 0; i < Ball.n; ++i)
        fprintf(res, "bin %s %2d %.5f %.5f %.0f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d\n", F10Name(Ball.j[i], Ball.e[i]).c_str(), Ball.xb[i], Ball.lo[i], Ball.hi[i],
                Ball.y[i], Ball.cb[i], pcb[i], pbg[i], bgStat[i], bgMC[i], psig[kS2][i], psig[kS3][i], psig[kS4][i], psig[kS5][i], psig[kSeta][i], muAll[i], ((oNom.useE >> Ball.e[i]) & 1));
    // window sums of the FMAIN background (exact integral of the smoothed densities with the tilt), per class
    fprintf(res, "\n# FMAIN background in windows (raw rows; c_j (1 + b x)[G_D + kappa_B G_C], exact integral), per M(ee) class: window E1 E2 (E1 = prediction)\n");
    {
        const double wl[4] = {kWin[kWB].lo, kWin[kWP].lo, kWin[kWSBlo].lo, kWin[kWCRU].lo}, wh[4] = {kWin[kWB].hi, kWin[kWP].hi, kWin[kWSBlo].hi, kWin[kWCRU].hi};
        const char* wn[5] = {"B", "P", "SBlo", "CRU", "SBhi"};
        auto bgW = [&](const Vec& th, const Templ& T, int e, double lo, double hi) {
            double s = 0;
            for (int j = 0; j < kNSlice; ++j) {
                const int h = 2 * j + e;
                double d0, d1, c0, c1; T.sgd[h].Int2(lo, hi, d0, d1); T.sgc[h].Int2(lo, hi, c0, c1);
                s += th[SNom.iC[j]] * ((d0 + th[SNom.iB] * d1) + th[SNom.iKB] * (c0 + th[SNom.iB] * c1));
            }
            return s;
        };
        for (int w = 0; w < 5; ++w) {
            const double lo = w < 4 ? wl[w] : kWin[kWSBhi].lo, hi = w < 4 ? wh[w] : kWin[kWSBhi].hi;
            fprintf(res, "bkgwin %s %.5f %.5f", wn[w], lo, hi);
            for (int e = 0; e < 2; ++e) {
                const double v = bgW(bNom.F.th, T0, e, lo, hi);
                Vec g = NumGrad([&](const Vec& th) { return bgW(th, T0, e, lo, hi); }, bNom.F.th, bNom.F.C);
                double s = 0, q = 0; int n = 0;
                for (size_t r = 0; r < R0.size(); ++r) if (bNom.Rp.fit[r].ok) { const double x = bgW(bNom.Rp.fit[r].th, R0[r], e, lo, hi); s += x; q += x * x; ++n; }
                fprintf(res, " E%d %.2f stat %.2f mc %.2f", e + 1, v, std::sqrt(std::max(0., Quad(bNom.F.Cstat, g))), n > 1 ? std::sqrt(std::max(0., (q - s * s / n) / (n - 1))) : 0.);
            }
            fprintf(res, "\n");
        }
        // X17: background shapes on the X2D M4 grid (0.5225-0.6100, 2.5 MeV), per class, summed over slices (raw rows of the F10 class)
        fprintf(res, "\n# X17 inputs: FMAIN background on the X2D M4 grid per M(ee) class (summed over slices; raw rows of the whole class; E1 = prediction)\n");
        fprintf(res, "x17_aeta %.6f stat %.6f total %.6f bias_corr %.6f\nx17_g %.5f total %.5f\nx17_kappaB %.5f total %.5f\nx17_b %.6f total %.6f  # per 0.1 GeV/c2\n",
                bNom.F.th[nA], bNom.Err(nA, 2), bNom.Err(nA, 4), bNom.Get("a_eta")->bc, bNom.F.th[nG], bNom.Err(nG, 4), bNom.F.th[nKB], bNom.Err(nKB, 4), bNom.F.th[nB], bNom.Err(nB, 4));
        for (int e = 0; e < 2; ++e)
            for (int i = 0; i < kX2Dnx; ++i) {
                const double lo = kX2Dxlo + i * 0.0025, hi = lo + 0.0025;
                double s = 0, q = 0; int n = 0;
                for (size_t r = 0; r < R0.size(); ++r) if (bNom.Rp.fit[r].ok) { const double x = bgW(bNom.Rp.fit[r].th, R0[r], e, lo, hi); s += x; q += x * x; ++n; }
                fprintf(res, "x17bkg E%d %2d %.5f %.5f %.3f mc %.3f\n", e + 1, i + 1, lo, hi, bgW(bNom.F.th, T0, e, lo, hi), n > 1 ? std::sqrt(std::max(0., (q - s * s / n) / (n - 1))) : 0.);
            }
    }
    fprintf(res, "\n# smoothing of the non-eta templates (wp4::SmoothFit): density [raw rows per GeV/c2] = exp(p0 + p1 x + p2 x^2 + p3 x^3), x = (M4 - 0.55)/0.1\n");
    for (int h = 0; h < kNH; ++h)
        fprintf(res, "smooth %s GD order %d %.8g %.8g %.8g %.8g neff %.1f scaled %d GC order %d %.8g %.8g %.8g %.8g neff %.1f scaled %d\n", HName(h).c_str(),
                T0.sgd[h].order, T0.sgd[h].p[0], T0.sgd[h].p[1], T0.sgd[h].p[2], T0.sgd[h].p[3], T0.sgd[h].neff, T0.sgd[h].scaled ? 1 : 0,
                T0.sgc[h].order, T0.sgc[h].p[0], T0.sgc[h].p[1], T0.sgc[h].p[2], T0.sgc[h].p[3], T0.sgc[h].neff, T0.sgc[h].scaled ? 1 : 0);
    fprintf(res, "\n# signal template window sums (raw rows at a_eta = kappa = rho_d = g = 1; x 64/63 not applied): comp B P SB CRU Fit tot\n");
    for (int t = 0; t < kNSig; ++t) { fprintf(res, "simwin %s", SigTag(t)); for (int w = 0; w < kNYW; ++w) fprintf(res, " %.3f", T0.W[t][w]); fprintf(res, "\n"); }
    fprintf(res, "datawin N"); for (int w = 0; w < kNYW; ++w) fprintf(res, " %.1f", D0.Nwin[w]); fprintf(res, "\n");
    fprintf(res, "datawin CB"); for (int w = 0; w < kNYW; ++w) fprintf(res, " %.2f", D0.CBwin[w]); fprintf(res, "\n");
    // profiles
    auto wprof = [&](const char* n, const Prof& p) { fprintf(res, "profile %s min %.4f", n, p.min); for (size_t i = 0; i < p.x.size(); ++i) fprintf(res, " %.6g:%.4f", p.x[i], p.chi2P[i] - p.min); fprintf(res, "\n"); };
    wprof("a_eta", prA); wprof("rho_d", prR); wprof("kappa_noconstraint", prK); wprof("g_noconstraint", prG);
    fclose(res);

    // ---------------------------------------------------------------- 16. tables
    gSystem->mkdir(kTabDir, kTRUE);
    auto ds = [](const DVal* d) { double c = 0; for (double x : d->eCons) c += x * x; return c; };
    // E1: parameters (nominal, spec E1+E2, pre, isBest == 1) and correlations of the nominal
    {
        FILE* t = fopen(Form("%s/wp5_tableE1.txt", kTabDir), "w");
        fprintf(t, "# Table E1 (WP5, CROSS-CHECK: F-main failed the WP4 closure): F-main parameters, data S0 (and pre, isBest == 1), L = %.1f nb^-1\n", kL);
        fprintf(t, "# value +- stat (sandwich: event groups + CB + constraints) +- MC (100 simulation replicas); [total], for the nominal column x %.2f = sqrt(chi2_GLS/ndf)\n", errScale);
        const Block* bl[4] = {&bNom, &bSpec, &bPre, &bBest};
        const char* bn[4] = {"FMAIN (E2, nominal)", "F-main E1+E2 (spec)", "pre (E2)", "isBest==1 (E2)"};
        fprintf(t, "%-10s", "param"); for (int q = 0; q < 4; ++q) fprintf(t, " | %-36s", bn[q]); fprintf(t, "\n");
        std::vector<std::string> pn = {"a_eta", "kappa", "rho_d", "g", "c_1", "c_2", "c_3", "c_4", "c_5", "b", "kappa_B", "k_E1", "k_E2"};
        for (auto& n : pn) {
            fprintf(t, "%-10s", n.c_str());
            for (int q = 0; q < 4; ++q) {
                const int k = bl[q]->S.Idx(n);
                if (k < 0) { fprintf(t, " | %-36s", "--"); continue; }
                fprintf(t, " | %8.4f +- %6.4f +- %6.4f [%6.4f]", bl[q]->F.th[k], bl[q]->Err(k, 2), bl[q]->Err(k, 3), bl[q]->Err(k, 4) * (q == 0 ? errScale : 1.));
            }
            fprintf(t, "\n");
        }
        fprintf(t, "%-10s", "chi2/ndf"); for (int q = 0; q < 4; ++q) fprintf(t, " | Pearson %6.1f/%3d GLS %-12s", bl[q]->F.chi2E, bl[q]->F.ndf, bl[q]->hasGLS ? Form("%.1f", bl[q]->G.chi2) : "--"); fprintf(t, "\n");
        fprintf(t, "\n# correlation matrix (stat + MC) of FMAIN\n%-8s", "");
        for (int k = 0; k < SNom.np(); ++k) fprintf(t, " %7s", SNom.names[k].c_str()); fprintf(t, "\n");
        for (int k = 0; k < SNom.np(); ++k) { fprintf(t, "%-8s", SNom.names[k].c_str()); for (int l = 0; l < SNom.np(); ++l) fprintf(t, " %7.3f", CT[k][l] / std::sqrt(CT[k][k] * CT[l][l])); fprintf(t, "\n"); }
        fclose(t);
        // tex
        FILE* x = fopen(Form("%s/wp5_tableE1.tex", kTabDir), "w");
        fprintf(x, "%% Table E1 (wp5_extract.C). CROSS-CHECK: F-main failed the WP4 closure.\n\\begin{tabular}{lrrrr}\n\\toprule\n");
        fprintf(x, "parameter & F-main (E2, nominal) & F-main E1+E2 (spec) & pre (E2) & isBest $=1$ (E2) \\\\\n\\midrule\n");
        std::map<std::string, std::string> tn = {{"a_eta", "$a_\\eta$"}, {"kappa", "$\\kappa$"}, {"rho_d", "$\\rho_d$"}, {"g", "$g$"}, {"c_1", "$c_1$"}, {"c_2", "$c_2$"}, {"c_3", "$c_3$"},
                                                 {"c_4", "$c_4$"}, {"c_5", "$c_5$"}, {"b", "$b$"}, {"kappa_B", "$\\kappa_B$"}, {"k_E1", "$k_{E1}$"}, {"k_E2", "$k_{E2}$"}};
        for (auto& n : pn) {
            fprintf(x, "%s", tn[n].c_str());
            for (int q = 0; q < 4; ++q) {
                const int k = bl[q]->S.Idx(n);
                if (k < 0) { fprintf(x, " & --"); continue; }
                const bool cons = [&]() { for (auto& c : bl[q]->S.cons) if (c.k == k) return true; return false; }();
                fprintf(x, " & $%.3f \\pm %.3f \\pm %.3f$%s", bl[q]->F.th[k], bl[q]->Err(k, 2), bl[q]->Err(k, 3), cons ? "$^{c}$" : "");
            }
            fprintf(x, " \\\\\n");
        }
        fprintf(x, "\\midrule\n$\\chi^2/\\mathrm{ndf}$ (Pearson)");
        for (int q = 0; q < 4; ++q) fprintf(x, " & %.0f/%d", bl[q]->F.chi2E, bl[q]->F.ndf);
        fprintf(x, " \\\\\n$\\chi^2/\\mathrm{ndf}$ (GLS, with MC cov.)");
        for (int q = 0; q < 4; ++q) fprintf(x, " & %s", bl[q]->hasGLS ? Form("%.0f/%d", bl[q]->G.chi2, bl[q]->G.ndf) : "--");
        fprintf(x, " \\\\\n\\bottomrule\n\\end{tabular}\n");
        fclose(x);
    }
    // E2: yields and ratios
    {
        struct Row { const Block* b; std::string q, label; double s10; };
        // S10 = |bias(a)| of WP4 RF (estimator as specified), relative
        std::vector<Row> rows = {
            {&bNom, "a_eta", "a_eta", 0.1006}, {&bNom, "Y2_B", "Y2 in B (C2)", 0}, {&bNom, "Y3_B", "Y3 in B (C3)", 0}, {&bNom, "Y_bump", "Y_bump = Y2 + Y3 in B", 0.1158},
            {&bNom, "Y4_P", "Y4 in P (C4)", 0}, {&bNom, "Y5_P", "Y5 in P (C5)", 0}, {&bNom, "Y_peak", "Y_peak = Y4 + Y5 in P", 0.0658},
            {&bNom, "Y2_tot", "Y2, all M4", 0}, {&bNom, "Y3_tot", "Y3, all M4", 0}, {&bNom, "Y4_tot", "Y4, all M4", 0}, {&bNom, "Y5_tot", "Y5, all M4", 0},
            {&bNom, "DR", "DR", 0}, {&bNom, "rho_d", "rho_d", 0.0969}, {&bNom, "BR_A", "BR_A = 2.68e-4 rho_d", 0.0969}, {&bGold, "BR_G", "BR_G (F-gold)", 0.1482},
            {&bFree12, "R_34", "R_34 (F-free E1+E2)", 0}, {&bFree12, "BR_C", "BR_C (F-free E1+E2; WITHDRAWN in E2)", 0.6687}, {&bFree, "R_34", "R_34 (F-free E2)", 0},
            {&bNom, "kappa", "kappa (1 +- 0.2 constrained)", 0}, {&bNom, "kappa_B", "kappa_B", 0}, {&bSpec, "kappa", "kappa (F-main E1+E2)", 0}, {&bSpec, "kappa_B", "kappa_B (F-main E1+E2)", 0}};
        FILE* t = fopen(Form("%s/wp5_tableE2.txt", kTabDir), "w");
        fprintf(t, "# Table E2 (WP5, CROSS-CHECK: F-main / F-free / F-gold failed the WP4 closure). Data S0, isBest >= 0, E2 fit (E1 dropped), L = %.1f nb^-1.\n", kL);
        fprintf(t, "# Yields = raw rows x 64/63 (not efficiency corrected). Columns: value, stat (event groups), CB (replicas + k_E2), nuisance (g, kappa constraints),\n");
        fprintf(t, "#   MC (simulation replicas; mc68 for ratios), total, total x sqrt(chi2_GLS/ndf) = x %.2f, S10 = WP4 RF bias of scenario (a) (relative),\n", errScale);
        fprintf(t, "#   bias-corrected value 2 v - mean(replicas) (information). Systematics: WP7.\n");
        fprintf(t, "%-40s %12s %10s %10s %10s %10s %10s %10s %7s %12s\n", "quantity", "value", "stat", "CB", "nuis", "MC", "total", "scaled", "S10", "bias-corr");
        FILE* x = fopen(Form("%s/wp5_tableE2.tex", kTabDir), "w");
        fprintf(x, "%% Table E2 (wp5_extract.C). CROSS-CHECK: the template fits failed the WP4 closure.\n\\begin{tabular}{lrrrrrrr}\n\\toprule\n");
        fprintf(x, "quantity & value & stat & CB & nuis. & MC & S10 & bias-corr. \\\\\n\\midrule\n");
        for (auto& r : rows) {
            const DVal* d = r.b->Get(r.q);
            if (!d) continue;
            // constraint parts: k_E2 -> CB, others -> nuisance
            double cbx = d->eCB * d->eCB, nu = 0;
            for (size_t q = 0; q < r.b->S.cons.size() && q < d->eCons.size(); ++q) {
                const std::string& pn = r.b->S.names[r.b->S.cons[q].k];
                if (pn == "k_E2") cbx += d->eCons[q] * d->eCons[q]; else nu += d->eCons[q] * d->eCons[q];
            }
            const bool ratio = r.q == "R_34" || r.q == "BR_C" || r.q == "BR_G" || r.q == "DR";
            const double mc = ratio ? d->eMC68 : d->eMC;
            const double tot = std::sqrt(d->eStat * d->eStat + mc * mc);
            fprintf(t, "%-40s %12.5g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %7s %12.5g\n", r.label.c_str(), d->v, d->eN, std::sqrt(cbx), std::sqrt(nu), mc, tot, tot * errScale,
                    r.s10 > 0 ? Form("%.1f%%", 100 * r.s10) : "--", d->bc);
            const bool br = r.q.find("BR") == 0;
            const double f = br ? 1e4 : 1.;
            const char* fmt = (std::fabs(d->v * f) > 100) ? "%.0f" : "%.3f";
            std::string lab = r.label; for (size_t p = 0; (p = lab.find('_', p)) != std::string::npos; p += 2) lab.replace(p, 1, "\\_");
            fprintf(x, "%s%s", lab.c_str(), br ? " [$10^{-4}$]" : "");
            for (double v : {d->v * f, d->eN * f, std::sqrt(cbx) * f, std::sqrt(nu) * f, mc * f}) fprintf(x, (std::string(" & ") + fmt).c_str(), v);
            fprintf(x, " & %s", r.s10 > 0 ? Form("%.0f\\,\\%%", 100 * r.s10) : "--");
            fprintf(x, (std::string(" & ") + fmt + " \\\\\n").c_str(), d->bc * f);
        }
        fprintf(x, "\\bottomrule\n\\end{tabular}\n");
        fclose(x);
        // isBest ratios and consistency tests
        fprintf(t, "\n# data-only ratios isBest == 1 / isBest >= 0 (same templates; stat = independent-sample upper bound, the samples are nested)\n");
        for (auto& br : bestRatio) fprintf(t, "%-10s %.4f +- %.4f (stat bound) +- %.4f (MC)\n", br.first.c_str(), br.second[0], br.second[1], br.second[2]);
        fprintf(t, "\n# F-gold: BR_G input error from kappa_B +- 1 sigma: %.3g\n", goldKBerr);
        fclose(t);
    }
    // acceptance file
    {
        FILE* t = fopen(Form("%s/wp5_acceptance.txt", kTabDir), "w");
        fprintf(t, "# WP5 acceptance criteria (SPEC report_all) for the nominal F-main (E2, cross-check) and the fall-backs\n");
        for (auto& a : acc) fprintf(t, "%-36s | %s | %s\n", a.crit.c_str(), a.value.c_str(), a.verdict.c_str());
        fprintf(t, "\n# fall-backs\n"); for (auto& f : fbl) fprintf(t, "- %s\n", f.c_str());
        fprintf(t, "\n# variants (S3/S4/S5/S8 members; GLS chi2 with MC covariance)\n");
        for (auto& v : vars) {
            const DVal* a = v.Get("a_eta"); const DVal* yb = v.Get("Y_bump"); const DVal* yp = v.Get("Y_peak"); const DVal* r = v.Get("rho_d");
            fprintf(t, "%-12s Pearson %.1f/%d GLS %s | a_eta %.4f (%+.1f%%) Y_bump %.0f (%+.1f%%) Y_peak %.0f (%+.1f%%) rho_d %.3f | %s\n", v.tag.c_str(), v.F.chi2E, v.F.ndf,
                    v.hasGLS ? Form("%.1f", v.G.chi2) : "--", a->v, 100 * (a->v / bNom.Get("a_eta")->v - 1), yb->v, 100 * (yb->v / bNom.Get("Y_bump")->v - 1), yp->v,
                    100 * (yp->v / bNom.Get("Y_peak")->v - 1), r->v, v.text.c_str());
        }
        fclose(t);
    }

    // ---------------------------------------------------------------- 17. figures
    bla::Style();
    gStyle->SetPadLeftMargin(0.16);
    const std::vector<double> e10 = F10Edges();
    auto mkH = [&](const char* n) { TH1D* h = new TH1D(n, "", kNB, e10.data()); h->SetDirectory(nullptr); return h; };
    // manual stack (THStack crashes in BuildAndPaint with this ROOT version): cumulative clones drawn from the top down
    auto drawStack = [&](const std::vector<TH1D*>& comp, const std::vector<int>& col, const std::string& nm) {
        std::vector<TH1D*> cum;
        for (size_t q = 0; q < comp.size(); ++q) {
            TH1D* h = (TH1D*)comp[q]->Clone(Form("%s_cum%zu", nm.c_str(), q)); h->SetDirectory(nullptr);
            if (q > 0) h->Add(cum.back());
            h->SetFillColor(col[q]); h->SetLineColor(col[q]); h->SetLineWidth(1); h->SetFillStyle(1001);
            cum.push_back(h);
        }
        for (int q = (int)cum.size() - 1; q >= 0; --q) cum[q]->Draw("HIST SAME");
        return cum.empty() ? 0. : cum.back()->GetMaximum();
    };
    const int colCB = kRed - 7, colBg = kGray, colEta = kOrange - 3, colC2 = kAzure + 1, colC34 = kCyan - 7, colC5 = kGreen + 2;
    // Fig E1: 10 panels + pulls for a fit block (E2 fitted; E1 fitted or predicted)
    auto figE1 = [&](const Block& b, const std::string& base, const std::string& title) {
        Bins Bx; MakeBins(*b.D, *b.T, 0, 3, Bx);
        Vec pc, pb, ps[kNSig], mu, V, Vd;
        Parts(Bx, b.S, b.F.th, pc, pb, ps);
        Eval(Bx, b.S, b.F.th, mu, nullptr, &V, &Vd);
        TCanvas* c = new TCanvas(Form("c_%s", base.c_str()), "", 2600, 1400);
        c->cd();
        TLatex tt; tt.SetNDC(); tt.SetTextSize(0.017); tt.DrawLatex(0.01, 0.975, title.c_str());
        for (int e = 1; e >= 0; --e)
            for (int j = 0; j < kNSlice; ++j) {
                const int row = 1 - e;   // row 0 = E2 (top), row 1 = E1
                const double x0 = 0.005 + j * 0.199, x1 = x0 + 0.199, ytop = 0.96 - row * 0.48, ybot = ytop - 0.48;
                const double ysplit = ybot + 0.30 * (ytop - ybot);
                TPad* pm = new TPad(Form("pm%d%d_%s", j, e, base.c_str()), "", x0, ysplit, x1, ytop);
                TPad* pp = new TPad(Form("pp%d%d_%s", j, e, base.c_str()), "", x0, ybot, x1, ysplit);
                pm->SetBottomMargin(0.012); pm->SetTopMargin(0.08); pm->SetLeftMargin(0.19); pm->SetRightMargin(0.03);
                pp->SetTopMargin(0.02); pp->SetBottomMargin(0.36); pp->SetLeftMargin(0.19); pp->SetRightMargin(0.03);
                c->cd(); pm->Draw(); pp->Draw();
                TH1D *hN = mkH(Form("hN%d%d_%s", j, e, base.c_str())), *hCB = mkH(Form("hCB%d%d%s", j, e, base.c_str())), *hBg = mkH(Form("hBg%d%d%s", j, e, base.c_str()));
                TH1D *hEt = mkH(Form("hEt%d%d%s", j, e, base.c_str())), *h2 = mkH(Form("h2%d%d%s", j, e, base.c_str())), *h34 = mkH(Form("h34%d%d%s", j, e, base.c_str()));
                TH1D *h5 = mkH(Form("h5%d%d%s", j, e, base.c_str())), *hP = mkH(Form("hP%d%d%s", j, e, base.c_str()));
                const bool fitted = (b.S.iK[e] >= 0);
                for (int i = 0; i < Bx.n; ++i) {
                    if (Bx.j[i] != j || Bx.e[i] != e) continue;
                    const int x = Bx.xb[i];
                    const double w = (Bx.hi[i] - Bx.lo[i]) * 1000.;
                    const double se = Bx.y[i] > 0 ? b.D->eg10[(size_t)Bx.cell[i] * 700 + Bx.cell[i]] / Bx.y[i] : 1.;
                    hN->SetBinContent(x, Bx.y[i] / w); hN->SetBinError(x, std::sqrt(std::max(0., Bx.y[i] * se)) / w);
                    hCB->SetBinContent(x, pc[i] / w); hBg->SetBinContent(x, pb[i] / w); hEt->SetBinContent(x, ps[kSeta][i] / w);
                    h2->SetBinContent(x, ps[kS2][i] / w); h34->SetBinContent(x, (ps[kS3][i] + ps[kS4][i]) / w); h5->SetBinContent(x, ps[kS5][i] / w);
                    const double vv = V[i] + (se - 1.) * std::max(mu[i], 1.);
                    hP->SetBinContent(x, fitted ? (Bx.y[i] - mu[i]) / std::sqrt(vv) : (mu[i] > 0 ? Bx.y[i] / mu[i] : 0.));
                    hP->SetBinError(x, fitted ? 0. : std::sqrt(vv) / std::max(mu[i], 1.));
                }
                struct L { TH1D* h; int c; const char* n; };
                L ls[6] = {{hCB, colCB, "CB = 2#sqrt{N_{++}N_{--}}"}, {hBg, colBg, "non-#eta background (G_{D} + #kappa_{B} G_{C})"}, {hEt, colEta, "G_{#eta}"},
                           {h34, colC34, "C3 + C4 (#gamma conversion)"}, {h2, colC2, "C2 (#pi^{0} Dalitz)"}, {h5, colC5, "C5 (#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-})"}};
                for (auto& l : ls) { l.h->SetFillColor(l.c); l.h->SetLineColor(l.c); l.h->SetLineWidth(1); }
                pm->cd();
                double stmax = 0; for (int x = 1; x <= kNB; ++x) { double t = 0; for (auto& l : ls) t += l.h->GetBinContent(x); stmax = std::max(stmax, t); }
                const double ymax = 1.25 * std::max(hN->GetMaximum(), stmax);
                TH1D* fr = mkH(Form("fr%d%d%s", j, e, base.c_str()));
                fr->SetMaximum(ymax); fr->SetMinimum(0.);
                fr->GetYaxis()->SetTitle("rows / (MeV/c^{2})"); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetLabelSize(0.065); fr->GetYaxis()->SetTitleOffset(1.25);
                fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetMaxDigits(3);
                fr->Draw("AXIS");
                drawStack({hCB, hBg, hEt, h34, h2, h5}, {colCB, colBg, colEta, colC34, colC2, colC5}, Form("stk%d%d%s", j, e, base.c_str()));
                hN->SetMarkerStyle(20); hN->SetMarkerSize(0.7); hN->SetMarkerColor(kBlack); hN->SetLineColor(kBlack); hN->Draw("E1 SAME");
                fr->Draw("AXIS SAME");
                TLatex lt; lt.SetNDC(); lt.SetTextSize(0.07);
                lt.DrawLatex(0.23, 0.84, Form("S%d, %s%s", j + 1, e == 0 ? "M_{ee} < 10" : "M_{ee} #geq 10", fitted ? " (fit)" : " (prediction)"));
                if (j == 4 && e == 1) {
                    TLegend* lg = new TLegend(0.21, 0.28, 0.66, 0.80); lg->SetTextSize(0.050); lg->SetFillStyle(1001); lg->SetFillColor(kWhite);
                    lg->AddEntry(hN, "data N (S0, isBest #geq 0)", "pe");
                    for (int q = 5; q >= 0; --q) lg->AddEntry(ls[q].h, ls[q].n, "f");
                    lg->Draw();
                }
                pp->cd();
                hP->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); hP->GetXaxis()->SetTitleSize(0.15); hP->GetXaxis()->SetLabelSize(0.13); hP->GetXaxis()->SetTitleOffset(1.0);
                hP->GetYaxis()->SetTitle(fitted ? "pull" : "data/model"); hP->GetYaxis()->SetTitleSize(0.15); hP->GetYaxis()->SetLabelSize(0.12); hP->GetYaxis()->SetTitleOffset(0.55);
                hP->GetYaxis()->SetNdivisions(505);
                if (fitted) { hP->SetMinimum(-7.5); hP->SetMaximum(7.5); hP->SetFillColor(kGray + 1); hP->SetLineColor(kGray + 2); hP->Draw("HIST"); }
                else { hP->SetMinimum(0.6); hP->SetMaximum(1.8); hP->SetMarkerStyle(20); hP->SetMarkerSize(0.5); hP->Draw("E1"); }
                TLine* l0 = new TLine(0.40, fitted ? 0. : 1., 0.70, fitted ? 0. : 1.); l0->SetLineStyle(2); l0->Draw();
                if (fitted) for (double y : {-2., 2.}) { TLine* l = new TLine(0.40, y, 0.70, y); l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw(); }
            }
        bla::SaveFig(c, Form("%s/%s", kFigDir, base.c_str()));
    };
    figE1(bNom, "wp5_E1", Form("F-main, CROSS-CHECK (WP4 fail). S0, isBest #geq 0, L = %.1f nb^{-1}. Top: fit, M_{ee} #geq 10 MeV/c^{2} (E1 dropped), pulls. "
                               "Bottom: M_{ee} < 10 MeV/c^{2}, prediction (k_{E1} = 1), data/model. #chi^{2}_{GLS}/ndf = %.2f, a_{#eta} = %.3f",
                               kL, chi2gls, bNom.F.th[nA]));
    figE1(bSpec, "wp5_E1_spec", Form("F-main as specified (E1 + E2, k_{E1} free = %.2f; CROSS-CHECK): S0, isBest #geq 0, L = %.1f nb^{-1}; pulls (N - #mu)/#sigma; a_{#eta} = %.3f",
                                     kE1, kL, bSpec.F.th[bSpec.S.Idx("a_eta")]));
    // Fig E2: N - kCB - background per slice (E2) with the eta templates
    {
        Vec pc, pb, ps[kNSig], mu, V, Vd;
        Parts(bNom.B, SNom, bNom.F.th, pc, pb, ps);
        Eval(bNom.B, SNom, bNom.F.th, mu, nullptr, &V, &Vd);
        // eta at SMASH (a = kappa = rho = g = 1)
        TCanvas* c = new TCanvas("cE2", "", 2100, 1300);
        c->Divide(3, 2, 0.002, 0.002);
        std::vector<TH1D*> hs[7];
        for (int p = 0; p < 6; ++p) {
            c->cd(p + 1);
            gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.03);
            TH1D *hD = mkH(Form("e2D%d", p)), *h2 = mkH(Form("e22%d", p)), *h34 = mkH(Form("e234%d", p)), *h5 = mkH(Form("e25%d", p)), *hE = mkH(Form("e2E%d", p)), *hS = mkH(Form("e2S%d", p));
            for (int i = 0; i < bNom.B.n; ++i) {
                if (p < 5 && bNom.B.j[i] != p) continue;
                const int x = bNom.B.xb[i];
                const double w = (bNom.B.hi[i] - bNom.B.lo[i]) * 1000.;
                const double se = bNom.B.y[i] > 0 ? D0.eg10[(size_t)bNom.B.cell[i] * 700 + bNom.B.cell[i]] / bNom.B.y[i] : 1.;
                const double vd = bNom.B.y[i] * se + std::pow(bNom.F.th[SNom.iK[1]], 2) * bNom.B.vcb[i];
                hD->SetBinContent(x, hD->GetBinContent(x) + (bNom.B.y[i] - pc[i] - pb[i]) / w);
                hD->SetBinError(x, std::sqrt(std::pow(hD->GetBinError(x), 2) + vd / (w * w)));
                h2->AddBinContent(x, ps[kS2][i] / w); h34->AddBinContent(x, (ps[kS3][i] + ps[kS4][i]) / w); h5->AddBinContent(x, ps[kS5][i] / w); hE->AddBinContent(x, ps[kSeta][i] / w);
                hS->AddBinContent(x, (bNom.B.s[kS2][i] + bNom.B.s[kS3][i] + bNom.B.s[kS4][i] + bNom.B.s[kS5][i] + bNom.B.s[kSeta][i]) / w);
            }
            struct L { TH1D* h; int c; const char* n; };
            L ls[4] = {{hE, colEta, "G_{#eta} (a_{#eta} g)"}, {h34, colC34, "C3 + C4 (a_{#eta} #kappa)"}, {h2, colC2, "C2 (a_{#eta})"}, {h5, colC5, "C5 (a_{#eta} #rho_{d})"}};
            for (auto& l : ls) { l.h->SetFillColor(l.c); l.h->SetLineColor(l.c); }
            double ymax = 0, ymin = 0;
            for (int x = 1; x <= kNB; ++x) { ymax = std::max(ymax, hD->GetBinContent(x) + hD->GetBinError(x)); ymin = std::min(ymin, hD->GetBinContent(x) - hD->GetBinError(x)); ymax = std::max(ymax, hS->GetBinContent(x)); }
            TH1D* fr = mkH(Form("e2fr%d", p));
            fr->SetMaximum(1.15 * ymax); fr->SetMinimum(1.15 * ymin);
            fr->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); fr->GetYaxis()->SetTitle("N - k CB - background [rows/(MeV/c^{2})]");
            fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleOffset(1.6); fr->GetXaxis()->SetTitleSize(0.052);
            fr->Draw("AXIS");
            drawStack({hE, h34, h2, h5}, {colEta, colC34, colC2, colC5}, Form("e2stk%d", p));
            hS->SetLineColor(kBlack); hS->SetLineStyle(2); hS->SetLineWidth(2); hS->Draw("HIST SAME");
            hD->SetMarkerStyle(20); hD->SetMarkerSize(0.8); hD->SetMarkerColor(kBlue + 1); hD->SetLineColor(kBlue + 1); hD->Draw("E1 SAME");
            TLine* l0 = new TLine(0.40, 0, 0.70, 0); l0->SetLineStyle(3); l0->Draw();
            fr->Draw("AXIS SAME");
            TLatex lt; lt.SetNDC(); lt.SetTextSize(0.05);
            lt.DrawLatex(0.21, 0.88, p < 5 ? Form("S%d, M_{ee} #geq 10 MeV/c^{2}", p + 1) : "all slices, M_{ee} #geq 10 MeV/c^{2}");
            // [v5 verification 2026-11-14] the legend was drawn white-filled in the S5 panel at (0.19, 0.20, 0.62, 0.56) and hid
            // the S5 data of M4 0.40-0.57 (negative residuals); it now sits, transparent, in the empty upper right of the
            // "all slices" panel. Corrected figure without rerunning WP5: docs/latex/figs_all/v5_wp5_E2_corrected.pdf (v5_figE2.C).
            if (p == 5) {
                TLegend* lg = new TLegend(0.50, 0.58, 0.96, 0.87); lg->SetTextSize(0.034); lg->SetFillStyle(0);
                lg->AddEntry(hD, "data N - k CB - fitted bkg.", "pe");
                for (int q = 3; q >= 0; --q) lg->AddEntry(ls[q].h, ls[q].n, "f");
                lg->AddEntry(hS, "all #eta, SMASH #times PDG BR", "l");
                lg->Draw();
            }
            if (p == 4) {
                lt.SetTextSize(0.040);
                lt.DrawLatex(0.21, 0.80, Form("a_{#eta} = %.3f, #kappa = %.2f, #rho_{d} = %.2f, g = %.2f", bNom.F.th[nA], bNom.F.th[nK], bNom.F.th[nR], bNom.F.th[nG]));
                lt.DrawLatex(0.21, 0.74, "CROSS-CHECK (F-main failed WP4)");
            }
        }
        bla::SaveFig(c, Form("%s/wp5_E2", kFigDir));
    }
    // Fig E3: profile scans
    {
        TCanvas* c = new TCanvas("cE3", "", 1800, 1350);
        c->Divide(3, 2, 0.002, 0.002);
        struct PP { const Prof* p; const char* xt; int k; bool cons; };
        PP pl[4] = {{&prA, "a_{#eta}", nA, false}, {&prR, "#rho_{d}", nR, false}, {&prK, "#kappa (constraint removed)", nK, true}, {&prG, "g (constraint removed)", nG, true}};
        for (int q = 0; q < 4; ++q) {
            c->cd(q + 1);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.15);
            const Prof& p = *pl[q].p;
            TGraph* g = new TGraph();
            for (size_t i = 0; i < p.x.size(); ++i) g->SetPoint(g->GetN(), p.x[i], p.chi2P[i] - p.min);
            g->SetLineWidth(3); g->SetLineColor(kBlack);
            double xmin = p.x.front(), xmax = p.x.back();
            TH1D* fr = new TH1D(Form("e3fr%d", q), "", 10, xmin, xmax); fr->SetDirectory(nullptr);
            fr->SetMinimum(0); fr->SetMaximum(9); fr->GetXaxis()->SetTitle(pl[q].xt); fr->GetYaxis()->SetTitle("#Delta#chi^{2}");
            fr->GetXaxis()->SetNdivisions(505);
            fr->Draw("AXIS");
            g->Draw("L SAME");
            // parabolas from the stat (sandwich) and total errors of the nominal fit
            int imin = 0; for (size_t i = 0; i < p.x.size(); ++i) if (p.chi2P[i] < p.chi2P[imin]) imin = (int)i;
            const double xh = pl[q].cons ? p.x[imin] : bNom.F.th[pl[q].k];
            const double es = bNom.Err(pl[q].k, 2), et = bNom.Err(pl[q].k, 4);
            TF1* f1 = new TF1(Form("e3s%d", q), Form("((x-%g)/%g)^2", xh, es), xmin, xmax); f1->SetLineColor(kBlue + 1); f1->SetLineStyle(2); f1->SetLineWidth(2);
            TF1* f2 = new TF1(Form("e3t%d", q), Form("((x-%g)/%g)^2", xh, et), xmin, xmax); f2->SetLineColor(kRed + 1); f2->SetLineStyle(7); f2->SetLineWidth(2);
            if (!pl[q].cons) { f1->Draw("SAME"); f2->Draw("SAME"); }
            for (double y : {1., 4.}) { TLine* l = new TLine(xmin, y, xmax, y); l->SetLineStyle(3); l->Draw(); }
            if (q == 0) {
                TLegend* lg = new TLegend(0.30, 0.66, 0.92, 0.90); lg->SetTextSize(0.040); lg->SetFillStyle(0);
                lg->AddEntry(g, "profile #Delta#chi^{2} (fixed fit weights)", "l");
                lg->AddEntry(f1, "stat (sandwich: event groups, CB)", "l");
                lg->AddEntry(f2, "stat #oplus MC (100 sim. replicas)", "l");
                lg->Draw();
            }
        }
        // s: global (parabola extrapolated) and peak region
        for (int q = 0; q < 2; ++q) {
            c->cd(5 + q);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.15);
            TGraph* g = new TGraph();
            double ymin = 1e30;
            if (q == 0) { for (int k = 0; k < 4; ++k) ymin = std::min(ymin, sscan[k].F.chi2P); for (int k = 0; k < 4; ++k) g->SetPoint(g->GetN(), sscan[k].s, sscan[k].F.chi2P - ymin); }
            else { for (auto& x : sLocPts) ymin = std::min(ymin, x[1]); for (auto& x : sLocPts) g->SetPoint(g->GetN(), x[0], x[1] - ymin); }
            g->SetMarkerStyle(20); g->SetMarkerSize(1.4);
            TH1D* fr = new TH1D(Form("e3s%d", q), "", 10, 0.9925, 0.9980); fr->SetDirectory(nullptr);
            double ym = 0; for (int i = 0; i < g->GetN(); ++i) ym = std::max(ym, g->GetY()[i]);
            fr->SetMinimum(0); fr->SetMaximum(std::max(9., 1.3 * ym)); fr->GetXaxis()->SetTitle("mass-transfer scale s"); fr->GetYaxis()->SetTitle("#Delta#chi^{2}");
            fr->GetXaxis()->SetNdivisions(505);
            fr->Draw("AXIS");
            const double sh = q == 0 ? sHat : sLoc, sg = q == 0 ? sSig : sLocSig;
            const double y1 = q == 0 ? sscan[2].F.chi2P - ymin : sLocPts[2][1] - ymin;
            TF1* f = new TF1(Form("e3sp%d", q), Form("%g + ((x-%g)/%g)^2 - ((%g-%g)/%g)^2", y1, sh, sg, kS, sh, sg), 0.9925, 0.9980);
            f->SetLineColor(kBlue + 1); f->SetLineStyle(2); f->Draw("SAME");
            g->Draw("P SAME");
            TLine* l = new TLine(kS, 0, kS, fr->GetMaximum()); l->SetLineColor(kRed + 1); l->SetLineStyle(7); l->Draw();
            TLatex lt; lt.SetNDC(); lt.SetTextSize(0.042);
            const double yt = q == 0 ? 0.45 : 0.85;
            lt.DrawLatex(0.20, yt, q == 0 ? "whole fit (0.40-0.70)" : "peak region 0.5225-0.565");
            lt.DrawLatex(0.20, yt - 0.06, Form("#hat{s} = %.5f #pm %.5f", sh, sg));
            lt.DrawLatex(0.20, yt - 0.12, "red: s = 0.99615 (nominal)");
        }
        bla::SaveFig(c, Form("%s/wp5_E3", kFigDir));
    }
    // Fig E4 (extra): a_eta and Y_peak/Y_bump of every fit (consistency summary)
    {
        struct Pt { std::string lab; double v, es, et; int col; };
        std::vector<Pt> pts;
        auto addp = [&](const Block& b, const std::string& lab, int col, const char* q = "a_eta") {
            const DVal* d = b.Get(q); if (!d) return;
            pts.push_back({lab, d->v, d->eStat, d->eTot, col});
        };
        addp(bNom, "F-main nominal (E2)", kBlack);
        addp(bSpec, "F-main E1+E2 (spec)", kGray + 2);
        std::map<std::string, std::string> vl = {{"FMAIN_bA", "b constrained to stage A"}, {"FMAIN_quad", "+ quadratic term (S3)"}, {"FMAIN_tiltj", "per-slice tilts (S3)"},
                                                 {"FMAIN_kfree", "#kappa free"}, {"FMAIN_gfree", "g free (S8)"}, {"FMAIN_r42", "range 0.42-0.70 (S5)"}, {"FMAIN_r45", "range 0.45-0.70 (S5)"},
                                                 {"FMAIN_raw", "raw non-#eta templates (S4)"}, {"FMAIN_pol2", "exp(pol2) smoothing (S4)"}};
        for (auto& v : vars) addp(v, vl.count(v.tag) ? vl[v.tag] : v.tag, kBlue + 1);
        for (int j = 0; j < kNSlice; ++j) { const int k = SPs.Idx(Form("a_eta_%d", j + 1)); pts.push_back({Form("slice S%d", j + 1), bPs.F.th[k], bPs.Err(k, 2), bPs.Err(k, 4), kGreen + 2}); }
        for (auto& f : folds) addp(f, f.tag == "FMAIN_fold0" ? "fold 0 (even groups)" : "fold 1 (odd groups)", kMagenta + 1);
        for (auto& p : pers) addp(p, "period " + p.tag.substr(7) + " (bkg fixed)", kOrange + 7);
        addp(bPre, "pre", kRed + 1);
        addp(bBest, "isBest == 1 (data)", kRed + 1);
        { const DVal* d = bGold.Get("alpha_2"); pts.push_back({"F-gold #alpha_{2}", d->v, d->eStat, d->eTot, kCyan + 2}); }
        { const DVal* d = bFree.Get("alpha_2"); pts.push_back({"F-free #alpha_{2} (E2)", d->v, d->eStat, d->eTot, kCyan + 2}); }
        TCanvas* c = new TCanvas("cE4", "", 1500, 1500);
        c->SetLeftMargin(0.36); c->SetBottomMargin(0.09); c->SetTopMargin(0.04);
        const int n = (int)pts.size();
        TH1D* fr = new TH1D("e4fr", "", 10, 0.15, 0.75); fr->SetDirectory(nullptr);
        fr->SetMinimum(0); fr->SetMaximum(n + 1); fr->GetYaxis()->SetLabelSize(0); fr->GetYaxis()->SetTickLength(0);
        fr->GetXaxis()->SetTitle("a_{#eta} (#eta scale relative to SMASH #times PDG BR)");
        fr->GetXaxis()->SetTitleSize(0.035); fr->GetXaxis()->SetLabelSize(0.033); fr->GetXaxis()->SetTitleOffset(1.1);
        fr->Draw("AXIS");
        TBox* band = new TBox(bNom.F.th[nA] - bNom.Err(nA, 2), 0, bNom.F.th[nA] + bNom.Err(nA, 2), n + 1); band->SetFillColorAlpha(kGray, 0.5); band->Draw();
        TLatex lt; lt.SetTextSize(0.024); lt.SetTextAlign(32);
        for (int i = 0; i < n; ++i) {
            const double y = n - i;
            TLine* lt2 = new TLine(pts[i].v - pts[i].et, y, pts[i].v + pts[i].et, y); lt2->SetLineColor(pts[i].col); lt2->SetLineWidth(1); lt2->Draw();
            TLine* ls = new TLine(pts[i].v - pts[i].es, y, pts[i].v + pts[i].es, y); ls->SetLineColor(pts[i].col); ls->SetLineWidth(5); ls->Draw();
            TGraph* g = new TGraph(1); g->SetPoint(0, pts[i].v, y); g->SetMarkerStyle(20); g->SetMarkerColor(pts[i].col); g->SetMarkerSize(1.3); g->Draw("P");
            lt.DrawLatex(0.14, y, pts[i].lab.c_str());
        }
        fr->Draw("AXIS SAME");
        TLatex l2; l2.SetNDC(); l2.SetTextSize(0.021); l2.DrawLatex(0.37, 0.968, "thick: stat; thin: stat #oplus MC; band: nominal (stat). Cross-checks");
        bla::SaveFig(c, Form("%s/wp5_E4", kFigDir));
    }
    TFile* fo = TFile::Open(Form("%s/wp5_fits.root", kOut), "RECREATE");
    fo->Close();
    P("\nwrote %s/results.txt, tables %s/wp5_table{E1,E2}.{txt,tex}, wp5_acceptance.txt, figures %s/wp5_E{1,1_spec,2,3,4}.pdf\n", kOut, kTabDir, kFigDir);
    fclose(gLog);
}

