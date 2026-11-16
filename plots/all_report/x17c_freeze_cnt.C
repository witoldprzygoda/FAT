// x17c_freeze_cnt.C -- COPY of x17c_freeze.C (figure agent B, 2026-11-15): "counts" instead of "rows" in the y titles of
// Figs. x17c_sb_mee and x17c_sb_m4 (-> docs/latex/figs_all/x17c_sb_mee_cnt, x17c_sb_m4_cnt .pdf/.png). Every other output (log,
// x17c_config.txt, x17c_state.root, the three tables, x17c_fill) goes to /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/;
// the frozen files are NOT touched. Reads the same inputs (hists_v1.root, wp5 results, x17c_egcov.root, x17_lines.root).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_x17c_cnt", kTRUE)' 'plots/all_report/x17c_freeze_cnt.C+O'
// ============================================================================
// x17c_freeze.C  --  WP9 X-c (prefix x17c): side-band fit (r = 0), nuisance sizes, variance factor, masked-row
//                    prediction and the Asimov expected limits that are frozen in x17_frozen.txt
// ============================================================================
// SPEC report_all WP9 X-c, D18, D19. Written 2026-11-14.
// BLINDING: only the X2D data cells with M(ee) 2-12 and 28-62 MeV/c2 are read (x17cl::LoadSideBands); the rows 12-28
// MeV/c2 (all M4 bins) never enter: they are masked in the likelihood, and their background is PREDICTED from the side
// bands (x17c::FillMaskedRows) before the background-only Asimov set is built.
// Steps:
//   1. Inputs: X2D (M4 0.5225-0.6100 x M(ee) 2-62 MeV/c2) of hists_v1.root, S0, data CB (like-sign, not blinded),
//      simulation templates S0__nom; background M4 shapes = WP5 x17bkg (E1, E2), X17 line = x17_lines.root.
//   2. Side-band fit with r = 0 (model of x17ab, bkgMode 4) and its goodness of fit (fall-back "Pearson chi2/ndf > 1.5").
//   3. Variants of the side-band fit (background shape, CB, peak/direct templates) -> N5-bias nuisance, the
//      direct-shape / peak-width / mass-scale terms of the X17 yield nuisance, and their effect on the expected limit.
//   4. Event-group variance factor phi of X17-like linear combinations (x17c_egcov.root of x17c_disp.C).
//   5. Masked-row prediction (n_j of 12-28 MeV/c2 from a pol2 in ln(n_j / simulated G_D + G_C)) and its variants.
//   6. Background-only Asimov set -> expected 90 % CLs limits (asymptotic, q~_r) on the mass grid, normalisations (1),
//      (2), (3); replica checks (CB, C4).
//   7. x17c_config.txt (frozen numbers for x17c_toys.C / x17c_summary.C / X-d), tables, side-band figures.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c", kTRUE)' \
//        'plots/all_report/x17c_freeze.C+O'
// ============================================================================
#include <algorithm>
#include <cstdarg>
#include <functional>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TMatrixDSym.h"
#include "TPad.h"
#include "THStack.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVectorD.h"

#include "x17c_lib.h"

namespace x17cf {
using namespace x17cl;
FILE* gLog = nullptr;
void P(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
void P(const char* fmt, ...) {
    va_list a;
    va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gLog) { va_start(a, fmt); vfprintf(gLog, fmt, a); va_end(a); }
}
// components of the model in cell (i, j): CB, other background, eta peak (D + C), eta-related terms, X17
void Comps(const x17c::Model& M, const std::vector<double>& p, int i, int j, double out[5]) {
    const double mu = M.Mu(p, i, j);
    std::vector<double> q = p;
    for (int jj = 0; jj < M.L.ny; ++jj) q[jj] = 0.;
    const double muB = mu - M.Mu(q, i, j);
    q = p; q[M.L.iN5] = 0; q[M.L.iN4] = 0;
    const double muP = mu - M.Mu(q, i, j);
    q = p; q[M.L.iA] = 0;
    const double muE = mu - M.Mu(q, i, j);
    q = p; q[M.L.iR] = 0;
    const double muX = mu - M.Mu(q, i, j);
    out[0] = p[M.L.iK + M.in->kcl[j]] * M.in->CB[M.in->G.idx(i, j)];
    out[1] = muB; out[2] = muP - (M.cfg.xNorm == 0 ? muX : 0.); out[3] = muE; out[4] = muX;
}
struct VarRes {
    std::string name, text;
    bool ok = false;
    double chi2 = 0, dev = 0; int ndf = 0;
    double N5 = 0, N4 = 0, a = 0, g = 0, RX = 0, N5ref = 0;
    double expMed = 0, expMed2 = 0, bX17c = 0;
    bool inN5spread = false;
};
}  // namespace x17cf

void x17c_freeze_cnt() {
    using namespace x17cf;
    gSystem->mkdir(kOutDir, kTRUE);
    gSystem->mkdir(kFrozenDir, kTRUE);
    gStyle->SetOptStat(0);
    gLog = fopen("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_freeze.log", "w");   // _cnt copy: scratch
    TFile* f = TFile::Open(alldefs::kHistFile, "READ");
    // ------------------------------------------------------------ 1. inputs
    Wp5 w = ReadWp5();
    if (!w.ok) return;
    // a_eta and g priors: WP5 total error x sqrt(chi2/ndf) (+) the largest shift of a WP5 variant that its chi2 does not
    // reject (S3 quadratic term: a_eta 0.5317, g 1.0056; wp5_results.md Sect. 2.4)
    const double aQuad = 0.53169366, gQuad = 1.0056172;
    Frozen z;
    z.a0 = w.a;
    z.aSig = std::hypot(w.aTot * w.chi2scale, aQuad - w.a);
    z.g0 = w.g;
    z.gSig = std::hypot(w.gTot * w.chi2scale, gQuad - w.g);
    // normalisation (3): peak-based eta scale a_P = Y_peak(data, WP5 F-main, rows x 64/63) / Y_peak(sim, SMASH, a = 1)
    {
        TH1* c4 = (TH1*)f->Get("sim/S0__nom/C4_win");
        TH1* c5 = (TH1*)f->Get("sim/S0__nom/C5_win");
        const double ysim = (c4->GetBinContent(alldefs::kWP + 1) + c5->GetBinContent(alldefs::kWP + 1)) * alldefs::kL;
        z.aP = w.Ypeak / ysim;
        z.aPSig = std::hypot(w.YpeakTot * w.chi2scale, w.YpeakQuad - w.Ypeak) / ysim;
        P("x17c: Y_peak sim (a = 1) %.1f rows x 64/63; data (WP5 F-main) %.1f +- %.1f (x %.2f) (+) quad %.1f -> a_P = %.4f +- %.4f\n", ysim, w.Ypeak,
          w.YpeakTot, w.chi2scale, w.YpeakQuad - w.Ypeak, z.aP, z.aPSig);
    }
    P("x17c: priors: a_eta %.4f +- %.4f (WP5 total %.4f x %.2f (+) quad shift %.4f); g %.3f +- %.3f\n", z.a0, z.aSig, w.aTot, w.chi2scale, aQuad - w.a, z.g0, z.gSig);
    BuildOpt bo;   // nominal: WP5 shapes, 2 classes, spread-rule CB
    x17c::Inputs I;
    if (!Build(f, I, w, bo, z.a0)) return;
    if (!SetLine(I, 17.0)) { P("x17c: line file missing\n"); return; }
    // conversion-tail nuisance from the 100 replicas of C4 (as x17ab_model.C)
    {
        TH3* h = (TH3*)f->Get("sim/S0__nom/C4_X2D_rep");
        TH2* h0 = (TH2*)f->Get("sim/S0__nom/C4_X2D");
        auto tail = [&](std::function<double(int, int)> get) {
            double t = 0, a = 0;
            for (int j = 0; j <= I.G.ny + 1; ++j) for (int i = 1; i <= I.G.nx; ++i) {
                const double v = get(i, j);
                a += v;
                if ((j >= 1 && j <= I.G.ny && I.Ctail[j - 1] > 0) || j == I.G.ny + 1) t += v;
            }
            return a > 0 ? t / a : 0.;
        };
        const double nom = tail([&](int i, int j) { return h0->GetBinContent(i, j); });
        double s = 0, s2 = 0; int nr = 0;
        for (int r = 1; r <= h->GetNbinsZ(); ++r) { const double v = tail([&](int i, int j) { return h->GetBinContent(i, j, r); }); s += v; s2 += v * v; ++nr; }
        const double m = s / nr;
        z.cTail = std::sqrt(std::max(0., s2 / nr - m * m)) / nom;
        P("x17c: conversion tail (M(ee) >= 12 MeV) fraction of C4 %.4f, replica RMS -> sC = %.3f\n", nom, z.cTail);
    }
    double xHist17 = 0, xLine17 = 0;
    {
        x17c::HistOpt o; x17c::Inputs Ih;
        x17c::BuildFromHists(f, Ih, o, alldefs::kNbToRaw);
        for (int j = 0; j < I.G.ny; ++j) { xHist17 += std::fabs(Ih.X[j] - I.X[j]); xLine17 = std::max(xLine17, std::fabs(Ih.X[j] - I.X[j])); }
        P("x17c: X17 line at 17 MeV: LineFor vs C6_X2D: sum|dX| %.2e, max %.2e; R_X %.4f vs %.4f; N5ref %.1f rows (a_eta = 1), s_X2D = %.4e rows per unit BR\n",
          xHist17, xLine17, I.RX, Ih.RX, I.N5ref, I.N5ref * I.RX / kBR);
    }
    // ------------------------------------------------------------ 2. data side bands and the nominal side-band fit
    std::vector<double> n;
    std::vector<char> use;
    if (!LoadSideBands(f, n, use)) return;
    double nSB = 0; for (double v : n) nSB += v;
    P("x17c: data side-band cells %d, rows %.0f (M(ee) 2-12 and 28-62 MeV/c2); rows 12-28 MeV/c2 masked\n", (int)std::count(use.begin(), use.end(), 1), nSB);
    // reference of the masked-row prediction: the simulated G_D + G_C row sums (raw), with their MC errors in the weights;
    // variant: locally smoothed simulation (SimRowSmooth)
    double smChi2D = 0, smChi2C = 0;
    const std::vector<double> simRow = SimRow(f, I);
    const std::vector<double> simRel = SimRowRelErr(f, I);
    const std::vector<double> simSm = SimRowSmooth(f, I, &smChi2D, &smChi2C);
    P("x17c: local smoothing of the simulated rows (variant): G_D exp(pol2) 8-62 MeV chi2/ndf %.2f, G_C exp(pol1) 8-28 MeV chi2/ndf %.2f\n", smChi2D, smChi2C);
    for (int j = 0; j < I.G.ny; ++j) P("x17c: sim row %2d: raw %9.1f +- %4.1f %%, local smooth %9.1f (raw/smooth %.3f)\n", j, simRow[j], 100 * simRel[j], simSm[j], simRow[j] / simSm[j]);
    x17c::Config C1 = MakeConfig(z, 1);
    x17c::Config Csb = C1; Csb.use = use;
    SideBand SB = FitSideBands(I, Csb, n, simRow, z.fillOrder, &simRel);
    x17c::Model Msb(I, Csb);
    const x17c::Layout& L = Msb.L;
    P("x17c: side-band fit (r = 0): ok %d, iter %d, %s; Pearson chi2 %.1f / ndf %d = %.3f (E1 %.1f/%d, E2 %.1f/%d); deviance %.1f\n", SB.fit.ok, SB.fit.iter,
      SB.fit.msg.c_str(), SB.chi2, SB.ndf, SB.chi2 / SB.ndf, SB.chi2c[0], SB.nc[0], SB.chi2c[1], SB.nc[1], SB.dev);
    const std::vector<double>& ps = SB.fit.p;
    const std::vector<double>& es = SB.fit.err;
    P("x17c: N5 %.1f +- %.1f (%.1f %%), N4 %.1f +- %.1f, k_E1 %.3f +- %.3f, k_E2 %.3f +- %.3f, a_eta %.4f +- %.4f, g %.3f +- %.3f, thC %.3f +- %.3f\n", ps[L.iN5],
      es[L.iN5], 100 * es[L.iN5] / ps[L.iN5], ps[L.iN4], es[L.iN4], ps[L.iK], es[L.iK], ps[L.iK + 1], es[L.iK + 1], ps[L.iA], es[L.iA], ps[L.iG], es[L.iG], ps[L.iTC],
      es[L.iTC]);
    for (int c = 0; c < I.nbcl; ++c)
        P("x17c: class %d: b %.4f +- %.4f, c %.4f +- %.4f, b' %.4f +- %.4f, c' %.4f +- %.4f\n", c, ps[L.iB + c], es[L.iB + c], ps[L.iC + c], es[L.iC + c], ps[L.iBy + c],
          es[L.iBy + c], ps[L.iCy + c], es[L.iCy + c]);
    for (int j = 0; j < I.G.ny; ++j) {
        const bool m = j >= kMaskLo && j <= kMaskHi;
        P("x17c: row %2d %4.0f-%4.0f MeV: %s n_j %10.1f +- %7.1f  sim %9.1f  ratio %.3f  chi2 %6.1f/35\n", j, 1e3 * I.G.ye[j], 1e3 * I.G.ye[j + 1], m ? "PRED" : "fit ",
          SB.pFill[j], m ? 0. : es[j], simRow[j], SB.pFill[j] / simRow[j], m ? 0. : SB.chi2row[j]);
    }
    double maxRowChi = 0; int maxRow = -1;
    for (int j = 0; j < I.G.ny; ++j) if (!(j >= kMaskLo && j <= kMaskHi) && SB.chi2row[j] > maxRowChi) { maxRowChi = SB.chi2row[j]; maxRow = j; }
    // residuals near m_eta (M4 0.5425-0.5525) per class
    double resE[2] = {0, 0}, varE[2] = {0, 0};
    for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) {
        const int b = I.G.idx(i, j);
        if (!use[b]) continue;
        if (!(I.G.xc(i) > 0.5425 && I.G.xc(i) < 0.5525)) continue;
        const double mu = Msb.Mu(ps, i, j);
        resE[I.kcl[j]] += n[b] - mu; varE[I.kcl[j]] += mu;
    }
    P("x17c: residual N - model at M4 0.5425-0.5525: E1 %.1f +- %.1f (%.2f sigma), E2 upper %.1f +- %.1f (%.2f sigma)\n", resE[0], std::sqrt(varE[0]),
      resE[0] / std::sqrt(varE[0]), resE[1], std::sqrt(varE[1]), resE[1] / std::sqrt(varE[1]));

    // ------------------------------------------------------------ 3. variance factor phi (event groups)
    double phiAll = 1, phiMax = 1, phiAvgUp = 0;
    std::vector<std::pair<double, double>> phiM;
    {
        TFile* fe = TFile::Open(Form("%s/x17c_egcov.root", kOutDir), "READ");
        TMatrixDSym* S = fe ? (TMatrixDSym*)fe->Get("egcov_sb") : nullptr;
        TVectorD* cnt = fe ? (TVectorD*)fe->Get("cnt_sb") : nullptr;
        if (!S || !cnt) { P("x17c: x17c_egcov.root missing: run x17c_disp.C first\n"); return; }
        const int nc = S->GetNrows();
        auto cellIJ = [&](int c, int& i, int& j) { const int jr = c / I.G.nx; i = c % I.G.nx; j = jr <= 4 ? jr : jr + 8; };
        // check: cnt equals the side-band data
        double md = 0;
        for (int c = 0; c < nc; ++c) { int i, j; cellIJ(c, i, j); md = std::max(md, std::fabs((*cnt)(c) - n[I.G.idx(i, j)])); }
        auto phiOf = [&](const std::vector<double>& wv) {
            double num = 0, den = 0;
            for (int a = 0; a < nc; ++a) {
                if (wv[a] == 0) continue;
                den += wv[a] * wv[a] * (*cnt)(a);
                for (int b = 0; b < nc; ++b) if (wv[b] != 0) num += wv[a] * wv[b] * (*S)(a, b);
            }
            return den > 0 ? num / den : 1.;
        };
        std::vector<double> w1(nc, 1.);
        phiAll = phiOf(w1);
        P("x17c: event groups: check max|cnt - N| = %.3g; phi(all side-band cells, w = 1) = %.4f\n", md, phiAll);
        for (int e = 0; e < 2; ++e) {
            std::vector<double> we(nc, 0.);
            for (int c = 0; c < nc; ++c) { int i, j; cellIJ(c, i, j); we[c] = I.kcl[j] == e ? 1. : 0.; }
            P("x17c: phi(class E%d, w = 1) = %.4f\n", e + 1, phiOf(we));
        }
        int nUp = 0;
        for (double m : {8., 10., 30., 32., 34., 36., 38., 40., 42., 44.}) {
            x17c::Inputs Iw = I;
            if (!SetLine(Iw, m)) continue;
            std::vector<double> wv(nc, 0.);
            for (int c = 0; c < nc; ++c) {
                int i, j; cellIJ(c, i, j);
                const double mu = Msb.Mu(ps, i, j);
                wv[c] = mu > 0 ? Iw.gX[i] * Iw.X[j] / mu : 0.;
            }
            const double ph = phiOf(wv);
            phiM.push_back({m, ph});
            phiMax = std::max(phiMax, ph);
            if (m >= 30) { phiAvgUp += ph; ++nUp; }
            P("x17c: phi(X17-like weights s/mu, line at %.0f MeV/c2, side-band cells only) = %.4f\n", m, ph);
        }
        phiAvgUp /= std::max(1, nUp);
        fe->Close();
    }
    // frozen variance factor: the largest phi of the X17-like combinations; the Pearson fall-back (> 1.5) overrides it
    const double chi2ndf = SB.chi2 / SB.ndf;
    const bool fbPearson = chi2ndf > 1.5;
    z.qScale = std::max(phiMax, fbPearson ? chi2ndf : 1.);
    P("x17c: variance factor: phi_max %.4f (upper-band mean %.4f), Pearson chi2/ndf %.3f -> fall-back %s; q~ divided by %.4f\n", phiMax, phiAvgUp, chi2ndf,
      fbPearson ? "FIRED" : "not fired", z.qScale);

    // ------------------------------------------------------------ 4. variants of the side-band fit
    // Every variant: its own side-band fit (r = 0) and masked-row prediction. Its expected median is computed after the
    // nuisance sizes are fixed (below), so that only the variant's choice differs from the nominal.
    std::vector<VarRes> vars;
    struct VarKeep { x17c::Inputs I; x17c::Config C; SideBand S; };
    std::vector<VarKeep> keep;
    double rxHistNom = 0;
    { x17c::Inputs It; Build(f, It, w, BuildOpt(), z.a0, &rxHistNom); }
    auto runVar = [&](const std::string& name, const std::string& text, const BuildOpt& b, const x17c::Config* Cmod, int fillOrder, int refMode, bool spread) {
        VarRes v; v.name = name; v.text = text; v.inN5spread = spread;
        x17c::Inputs Iv;
        double rxh = 0;
        if (!Build(f, Iv, w, b, z.a0, &rxh)) return;
        if (!SetLine(Iv, 17.0)) return;
        Iv.RX *= rxh / rxHistNom;   // R_X change of a simulation variant (C6/C5 of its directory) carried onto the line
        x17c::Config Cv = Cmod ? *Cmod : C1;
        Cv.use = use;
        std::vector<double> ref = refMode == 0 ? simRow : (refMode == 1 ? LsRow(f, Iv) : (refMode == 3 ? simSm : std::vector<double>(Iv.G.ny, 1.)));
        SideBand S = FitSideBands(Iv, Cv, n, ref, fillOrder, refMode == 0 ? &simRel : nullptr);
        x17c::Model Mv(Iv, Cv);
        v.ok = S.fit.ok; v.chi2 = S.chi2; v.dev = S.dev; v.ndf = S.ndf;
        v.N5 = S.fit.p[Mv.L.iN5]; v.N4 = S.fit.p[Mv.L.iN4]; v.a = S.fit.p[Mv.L.iA]; v.g = S.fit.p[Mv.L.iG]; v.RX = Iv.RX; v.N5ref = Iv.N5ref;
        vars.push_back(v);
        keep.push_back({Iv, Cv, S});
        P("x17c: variant %-10s ok %d chi2 %.1f/%d dev %.1f  N5 %.1f N4 %.1f a %.3f g %.3f R_X %.4f N5ref %.1f  [%s]\n", name.c_str(), v.ok, v.chi2, v.ndf, v.dev, v.N5,
          v.N4, v.a, v.g, v.RX, v.N5ref, text.c_str());
    };
    {
        BuildOpt b;
        runVar("nominal", "WP5 background shapes, bkgMode 4, 2 classes, spread-rule CB", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); b.bshape = 1; runVar("Bsim2", "background shape: simulated G_D + G_C, exp(pol2) (x17ab nominal)", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); b.bshape = 2; runVar("Bsim1", "background shape: simulated, exp(pol1)", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); b.bshape = 3; runVar("Bsim3", "background shape: simulated, exp(pol3)", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); b.bshape = 4; runVar("BkappaB", "background shape: simulated G_D + kappa_B(WP5) G_C, exp(pol2)", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); x17c::Config Cv = C1; Cv.bkgMode = 2; runVar("flatPol2", "background flat x (1 + b t + c t^2) per class (spec 'pol2' variant)", b, &Cv, z.fillOrder, 0, true);
        b = BuildOpt(); b.bshape = 1; b.nbcl = 3; runVar("3class", "three background classes (E1, 10-30, >= 30 MeV), simulated exp(pol2)", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); Cv = C1; Cv.bkgMode = 0; runVar("tiltOnly", "WP5 shapes x (1 + b t) per class", b, &Cv, z.fillOrder, 0, true);
        b = BuildOpt(); Cv = C1; Cv.bkgMode = 3; runVar("noMeeDep", "WP5 shapes x (1 + b t + c t^2), no M(ee) dependence", b, &Cv, z.fillOrder, 0, true);
        b = BuildOpt(); b.cbMode = 1; runVar("CBfact", "CB factorised (sum rule)", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); b.cbMode = 2; runVar("CBgeo", "CB per-row geometric mean", b, nullptr, z.fillOrder, 0, true);
        b = BuildOpt(); Cv = C1; Cv.kSig[0] = -1; runVar("kE1free", "k_E1 free", b, &Cv, z.fillOrder, 0, true);
        // fill variants (do not change the side-band fit)
        b = BuildOpt(); runVar("fill1", "masked rows: pol1 in ln(n_j / sim_j)", b, nullptr, 1, 0, false);
        runVar("fill3", "masked rows: pol3 in ln(n_j / sim_j)", b, nullptr, 3, 0, false);
        runVar("fillLS", "masked rows: reference = like-sign rows N++ + N--", b, nullptr, 2, 1, false);
        runVar("fillFlat", "masked rows: reference flat (pol2 in ln n_j)", b, nullptr, 2, 2, false);
        runVar("fillSmooth", "masked rows: reference = locally smoothed simulated G_D + G_C rows", b, nullptr, 2, 3, false);
        // peak / direct templates (yield factor of r)
        for (const char* var : {"V2", "GS", "r089", "r110", "sUp", "sDn"}) {
            b = BuildOpt(); b.var = var;
            runVar(std::string("sim_") + var, std::string("peak/direct templates of S0__") + var, b, nullptr, z.fillOrder, 0, false);
        }
    }
    const VarRes& v0 = vars[0];
    // N5-bias: half spread of N5 over the background-shape and CB variants whose side-band deviance is not worse than the
    // nominal by more than 2 sqrt(2 ndf) (rule fixed before the fits were run)
    double n5lo = 1e30, n5hi = -1e30;
    const double devCut = 2 * std::sqrt(2. * v0.ndf);
    for (auto& v : vars) {
        if (!v.inN5spread) continue;
        if (v.dev - v0.dev > devCut) { v.inN5spread = false; P("x17c: variant %s excluded from the N5 spread (deviance +%.1f > %.1f)\n", v.name.c_str(), v.dev - v0.dev, devCut); continue; }
        n5lo = std::min(n5lo, v.N5); n5hi = std::max(n5hi, v.N5);
    }
    const double sN5 = 0.5 * (n5hi - n5lo) / v0.N5;
    // yield-factor terms: direct shape (V2, GS), peak width (r089, r110), mass scale (sUp, sDn).
    // Normalisation (1): X17 rows per unit r = N5_hat R_X; (2)/(3): N5ref R_X (= C6 per unit BR x 2.68e-4 x L).
    auto rel = [&](const char* nm, int norm) {
        for (const auto& v : vars) if (v.name == nm) return norm == 1 ? (v.N5 * v.RX) / (v0.N5 * v0.RX) - 1. : (v.N5ref * v.RX) / (v0.N5ref * v0.RX) - 1.;
        return 0.;
    };
    double dDir[3] = {0, 0, 0}, dWid[3] = {0, 0, 0}, dScl[3] = {0, 0, 0};
    for (int nm = 1; nm <= 2; ++nm) {
        dDir[nm] = std::max(std::fabs(rel("sim_V2", nm)), std::fabs(rel("sim_GS", nm)));
        dWid[nm] = 0.5 * (std::fabs(rel("sim_r089", nm)) + std::fabs(rel("sim_r110", nm)));
        dScl[nm] = 0.5 * (std::fabs(rel("sim_sUp", nm)) + std::fabs(rel("sim_sDn", nm)));
    }
    // R_X kinematics over the mass grid (direct-like M(pipi), generator pion angle) and the slice/C6 validation (x17ab: 4.3 %)
    double dDL = 0, dNP = 0;
    for (double m : MassGrid()) {
        std::vector<double> X; double r0, r1, r3;
        x17c::LineFor(m, I.G, X, r0, 0); x17c::LineFor(m, I.G, X, r1, 1); x17c::LineFor(m, I.G, X, r3, 3);
        dDL = std::max(dDL, std::fabs(r1 / r0 - 1)); dNP = std::max(dNP, std::fabs(r3 / r0 - 1));
    }
    const double dVal = 0.043, dBRdir = 0.11 / 2.68;
    const double kin2 = dDL * dDL + dNP * dNP + dVal * dVal;
    const double sA = z.aSig / z.a0, sAP = z.aPSig / z.aP;
    z.xSig1 = std::sqrt(kin2 + dBRdir * dBRdir + sN5 * sN5 + dDir[1] * dDir[1] + dWid[1] * dWid[1] + dScl[1] * dScl[1]);
    z.xSig2 = std::sqrt(kin2 + dDir[2] * dDir[2] + dWid[2] * dWid[2] + dScl[2] * dScl[2] + sA * sA);
    z.xSig3 = std::sqrt(kin2 + dDir[2] * dDir[2] + dWid[2] * dWid[2] + dScl[2] * dScl[2] + sAP * sAP);
    P("x17c: X17 yield nuisance (1): DL %.3f NP %.3f valid %.3f BR_dir %.3f N5bias %.3f (N5 %.1f-%.1f) dir %.3f width %.3f scale %.3f -> sX1 = %.3f\n", dDL, dNP, dVal,
      dBRdir, sN5, n5lo, n5hi, dDir[1], dWid[1], dScl[1], z.xSig1);
    P("x17c: X17 yield nuisance (2): DL %.3f NP %.3f valid %.3f dir %.3f width %.3f scale %.3f a_eta %.3f -> sX2 = %.3f; (3): a_P %.3f -> sX3 = %.3f\n", dDL, dNP, dVal,
      dDir[2], dWid[2], dScl[2], sA, z.xSig2, sAP, z.xSig3);
    // fall-back: N5 error > 50 %
    const double n5rel = es[L.iN5] / ps[L.iN5];
    const bool fbN5 = n5rel > 0.5;
    P("x17c: fall-back N5 error > 50 %%: %s (%.1f %%)\n", fbN5 ? "FIRED -> normalisation (2) primary" : "not fired", 100 * n5rel);
    WriteConfig(z, "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_config.txt");   // _cnt copy: scratch (the frozen config is not touched)
    // expected medians of the variants (each from its own side-band fit and masked-row prediction), final nuisance sizes
    for (size_t k = 0; k < vars.size(); ++k) {
        for (int norm = 1; norm <= 2; ++norm) {
            x17c::Config Ca = MakeConfig(z, norm);
            Ca.bkgMode = keep[k].C.bkgMode; Ca.kSig[0] = keep[k].C.kSig[0]; Ca.kSig[1] = keep[k].C.kSig[1];
            x17c::Model Ma(keep[k].I, Ca);
            const std::vector<double> as = Ma.Expected(keep[k].S.pFill);
            if (norm == 1) vars[k].bX17c = WindowX17c(keep[k].I.G, as);
            Lim l = LimitBR(keep[k].I, Ca, as, true);
            (norm == 1 ? vars[k].expMed : vars[k].expMed2) = l.exp[2];
        }
        P("x17c: variant %-10s expected median at 17 MeV/c2: (1) %.3e (x %.3f)  (2) %.3e (x %.3f); background-only rows in X17c %.1f\n", vars[k].name.c_str(),
          vars[k].expMed, vars[k].expMed / vars[0].expMed, vars[k].expMed2, vars[k].expMed2 / vars[0].expMed2, vars[k].bX17c);
    }

    // ------------------------------------------------------------ 5./6. Asimov expected limits on the mass grid
    // nominal side-band fit again with the final configuration (only q scale / xSig changed: no effect on the r = 0 fit)
    C1 = MakeConfig(z, 1);
    Csb = C1; Csb.use = use;
    SB = FitSideBands(I, Csb, n, simRow, z.fillOrder, &simRel);
    x17c::Config Call = C1; Call.use.clear();
    x17c::Model Mall(I, Call);
    const std::vector<double> asimov = Mall.Expected(SB.pFill);
    double nAs = 0, nAsMask = 0;
    for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) { nAs += asimov[I.G.idx(i, j)]; if (j >= kMaskLo && j <= kMaskHi) nAsMask += asimov[I.G.idx(i, j)]; }
    P("x17c: background-only Asimov set: %.0f rows, of which %.0f in the masked rows (prediction); counting window X17c: %.1f rows (CB part %.1f)\n", nAs, nAsMask,
      WindowX17c(I.G, asimov), WindowX17c(I.G, I.CB));
    const std::vector<double> mg = MassGrid();
    struct MRow { double m, RX, e1[5], e2[5], e3[5], nX, N5, s1x; };
    std::vector<MRow> rows;
    for (double m : mg) {
        MRow r{}; r.m = m;
        x17c::Inputs Im = I;
        if (!SetLine(Im, m)) continue;
        r.RX = Im.RX;
        for (int norm = 1; norm <= 3; ++norm) {
            x17c::Config Cn = MakeConfig(z, norm);
            Lim l = LimitBR(Im, Cn, asimov, true);
            double* e = norm == 1 ? r.e1 : (norm == 2 ? r.e2 : r.e3);
            for (int k = 0; k < 5; ++k) e[k] = l.exp[k];
            if (norm == 1) { r.nX = l.nX; r.N5 = l.N5; }
        }
        r.s1x = Im.N5ref * Im.RX / kBR;
        rows.push_back(r);
        P("x17c: m %5.2f R_X %.4f  expected (1) %.3e [%.3e %.3e | %.3e %.3e]  (2) %.3e  (3) %.3e  rows at median(1) %.1f\n", m, r.RX, r.e1[2], r.e1[0], r.e1[1], r.e1[3],
          r.e1[4], r.e2[2], r.e3[2], r.nX);
    }
    // effect of single choices on the expected median at 17 MeV/c2 (Asimov from the nominal side-band fit)
    auto med17 = [&](const x17c::Inputs& Iv, const x17c::Config& Cv) { Lim l = LimitBR(Iv, Cv, asimov, true); return l.exp[2]; };
    x17c::Inputs I17 = I; SetLine(I17, 17.0);
    const double base17 = med17(I17, C1);
    struct Eff { std::string name; double v; };
    std::vector<Eff> eff;
    {
        x17c::Config Cv = C1; Cv.qScale = 1.; eff.push_back({"no variance factor (q~ not scaled)", med17(I17, Cv)});
        Cv = C1; Cv.xSig = 0.; eff.push_back({"no X17 yield nuisance", med17(I17, Cv)});
        Cv = C1; Cv.wSig = 0.; Cv.sSig = 0.; eff.push_back({"no line-shape nuisances", med17(I17, Cv)});
        Cv = C1; Cv.aSig = 0.; Cv.gSig = 0.; eff.push_back({"a_eta and g fixed", med17(I17, Cv)});
        Cv = C1; Cv.cTailSig = 0.; eff.push_back({"no conversion-tail nuisance", med17(I17, Cv)});
        Cv = C1; Cv.kSig[1] = 0.; eff.push_back({"k_E2 fixed", med17(I17, Cv)});
        x17c::Inputs Iv = I; SetLine(Iv, 17.0, 0, 1.10); eff.push_back({"line width x 1.10 (fixed)", med17(Iv, C1)});
        Iv = I; SetLine(Iv, 17.0, 0, 0.90); eff.push_back({"line width x 0.90 (fixed)", med17(Iv, C1)});
        Iv = I; SetLine(Iv, 17.0, 2); eff.push_back({"R_X with e+ decay angle 1 + cos^2 (quoted separately)", med17(Iv, C1)});
        Iv = I; SetLine(Iv, 17.0, 1); eff.push_back({"R_X with direct-like M(pipi)", med17(Iv, C1)});
    }
    for (const auto& e : eff) P("x17c: expected median at 17 MeV with '%s': %.3e (x %.3f)\n", e.name.c_str(), e.v, e.v / base17);
    // replica checks on the Asimov set: CB replicas and C4 replicas (r_hat spread, BR units)
    double cbRms = 0, c4Rms = 0, sigR = 0;
    {
        x17c::Model M17(I17, C1);
        x17c::Profile P0 = x17c::ProfileFits(M17, asimov, x17c::GlobOf(C1));
        sigR = P0.rhatErr * kBR;
        std::vector<double> rc, r4;
        for (int k = 0; k < 50; ++k) {
            BuildOpt b; b.cbRep = k;
            x17c::Inputs Iv; if (!Build(f, Iv, w, b, z.a0)) break; SetLine(Iv, 17.0);
            x17c::Model Mv(Iv, C1);
            rc.push_back(x17c::ProfileFits(Mv, asimov, x17c::GlobOf(C1)).rhat * kBR);
            b = BuildOpt(); b.simRep = k;
            x17c::Inputs Iw; if (!Build(f, Iw, w, b, z.a0)) break; SetLine(Iw, 17.0);
            x17c::Model Mw(Iw, C1);
            r4.push_back(x17c::ProfileFits(Mw, asimov, x17c::GlobOf(C1)).rhat * kBR);
        }
        auto rms = [](const std::vector<double>& v) { double m = 0, s = 0; for (double x : v) m += x; m /= std::max<size_t>(1, v.size()); for (double x : v) s += (x - m) * (x - m); return std::sqrt(s / std::max<size_t>(1, v.size() - 1)); };
        cbRms = rms(rc); c4Rms = rms(r4);
        P("x17c: Asimov at 17 MeV: sigma(r_hat) = %.3e (BR); CB replicas (50): r_hat RMS %.3e; C4 replicas (50): r_hat RMS %.3e\n", sigR, cbRms, c4Rms);
    }

    // ------------------------------------------------------------ 7. outputs
    // state file for x17c_toys.C / x17c_summary.C
    {
        TFile* fo = TFile::Open("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_state.root", "RECREATE");   // _cnt copy: scratch
        TH2D ha("asimov", "background-only Asimov set (side-band fit + predicted masked rows), raw rows", I.G.nx, I.G.xe.data(), I.G.ny, I.G.ye.data());
        TH2D hd("data_sb", "data side-band rows (masked rows = 0)", I.G.nx, I.G.xe.data(), I.G.ny, I.G.ye.data());
        for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) { ha.SetBinContent(i + 1, j + 1, asimov[I.G.idx(i, j)]); hd.SetBinContent(i + 1, j + 1, n[I.G.idx(i, j)]); }
        ha.Write(); hd.Write();
        TVectorD pv((int)SB.pFill.size());
        for (size_t k = 0; k < SB.pFill.size(); ++k) pv((int)k) = SB.pFill[k];
        pv.Write("pFill");
        TTree t("exp", "Asimov expected limits (BR_X BR_ee)");
        MRow r;
        t.Branch("m", &r.m); t.Branch("RX", &r.RX); t.Branch("e1", r.e1, "e1[5]/D"); t.Branch("e2", r.e2, "e2[5]/D"); t.Branch("e3", r.e3, "e3[5]/D");
        t.Branch("nX", &r.nX); t.Branch("s1x", &r.s1x);
        for (const auto& x : rows) { r = x; t.Fill(); }
        t.Write();
        fo->Close();
    }
    // tables
    FILE* tb = fopen("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_sidebands.txt", "w");   // _cnt copy: scratch
    fprintf(tb, "# x17c_sidebands.txt -- X17 side-band fit with r = 0 (x17c_freeze.C). Data S0, isBest >= 0, X2D: M4 0.5225-0.6100 GeV/c2 x M(ee)\n");
    fprintf(tb, "# 2-12 and 28-62 MeV/c2 (rows 12-28 MeV/c2 masked, D19). Raw rows; L = 6079.2 nb^-1. Model of x17ab (bkgMode 4) with the WP5 background shapes.\n");
    fprintf(tb, "fit ok %d; Pearson chi2 %.1f / ndf %d = %.3f (E1 %.1f / %d cells, E2 %.1f / %d cells); deviance %.1f; largest row chi2 %.1f/35 (row %.0f-%.0f MeV)\n", SB.fit.ok,
            SB.chi2, SB.ndf, SB.chi2 / SB.ndf, SB.chi2c[0], SB.nc[0], SB.chi2c[1], SB.nc[1], SB.dev, maxRowChi, 1e3 * I.G.ye[maxRow], 1e3 * I.G.ye[maxRow + 1]);
    fprintf(tb, "N5 %.1f +- %.1f; N4 %.1f +- %.1f; k_E1 %.3f +- %.3f; k_E2 %.3f +- %.3f; a_eta %.4f +- %.4f; g %.3f +- %.3f; thC %.3f +- %.3f\n", ps[L.iN5], es[L.iN5],
            ps[L.iN4], es[L.iN4], ps[L.iK], es[L.iK], ps[L.iK + 1], es[L.iK + 1], ps[L.iA], es[L.iA], ps[L.iG], es[L.iG], ps[L.iTC], es[L.iTC]);
    fprintf(tb, "residual N - model at M4 0.5425-0.5525: E1 %.1f +- %.1f, E2 %.1f +- %.1f rows\n", resE[0], std::sqrt(varE[0]), resE[1], std::sqrt(varE[1]));
    fprintf(tb, "variance factor: phi(all) %.4f; X17-like:", phiAll);
    for (auto& q : phiM) fprintf(tb, " %.0f MeV %.4f;", q.first, q.second);
    fprintf(tb, " q~ divided by %.4f\n", z.qScale);
    fprintf(tb, "# Mee_lo Mee_hi  status   n_j          err       sim_j     n_j/sim_j  chi2_row/35\n");
    for (int j = 0; j < I.G.ny; ++j) {
        const bool m = j >= kMaskLo && j <= kMaskHi;
        fprintf(tb, "%5.0f %5.0f  %s %11.1f %9.1f %11.1f %9.4f %9.1f\n", 1e3 * I.G.ye[j], 1e3 * I.G.ye[j + 1], m ? "predicted" : "fitted   ", SB.pFill[j], m ? 0. : es[j],
                simRow[j], SB.pFill[j] / simRow[j], m ? 0. : SB.chi2row[j]);
    }
    fclose(tb);
    tb = fopen("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_variants.txt", "w");   // _cnt copy: scratch
    fprintf(tb, "# x17c_variants.txt -- variants of the X17 side-band fit (r = 0; x17c_freeze.C). Each: side-band Pearson chi2/ndf, deviance, N5 (direct rows,\n");
    fprintf(tb, "# all M(ee)), N4, a_eta, g, R_X, and the Asimov expected median 90 %% CLs at 17 MeV/c2 built from the variant's own side-band fit.\n");
    fprintf(tb, "# 'N5spread' = used for the N5-bias nuisance (background-shape and CB variants with deviance <= nominal + 2 sqrt(2 ndf) = +%.1f).\n", devCut);
    fprintf(tb, "%-10s %10s %8s %9s %9s %7s %7s %7s %11s %11s %8s %s  %s\n", "variant", "chi2/ndf", "ddev", "N5", "N4", "a_eta", "g", "R_X", "exp(1)", "exp(2)", "B(X17c)",
            "N5spread", "description");
    for (const auto& v : vars)
        fprintf(tb, "%-10s %10.4f %+8.1f %9.1f %9.1f %7.4f %7.3f %7.4f %11.3e %11.3e %8.1f %8s  %s%s\n", v.name.c_str(), v.chi2 / v.ndf, v.dev - v0.dev, v.N5, v.N4, v.a,
                v.g, v.RX, v.expMed, v.expMed2, v.bX17c, v.inN5spread ? "yes" : "no", v.text.c_str(), v.ok ? "" : " [fit not converged]");
    fprintf(tb, "N5-bias = half spread %.1f-%.1f / %.1f = %.4f\n", n5lo, n5hi, v0.N5, sN5);
    fprintf(tb, "yield-factor terms (1): direct shape %.4f, peak width %.4f, mass scale %.4f; (2): %.4f, %.4f, %.4f\n", dDir[1], dWid[1], dScl[1], dDir[2], dWid[2], dScl[2]);
    fprintf(tb, "R_X kinematics (max over the mass grid): direct-like M(pipi) %.4f, generator pion angle %.4f; slice/C6 validation %.4f; BR_dir %.4f\n", dDL, dNP, dVal, dBRdir);
    fprintf(tb, "sX(1) = %.4f, sX(2) = %.4f (incl. a_eta %.4f), sX(3) = %.4f (incl. a_P %.4f); all log-normal\n", z.xSig1, z.xSig2, z.aSig / z.a0, z.xSig3, z.aPSig / z.aP);
    fprintf(tb, "# single choices on the expected median at 17 MeV/c2 (nominal %.4e):\n", base17);
    for (const auto& e : eff) fprintf(tb, "  %-55s %.4e  x %.3f\n", e.name.c_str(), e.v, e.v / base17);
    fprintf(tb, "Asimov at 17 MeV: sigma(r_hat) %.3e (BR); r_hat RMS over 50 CB replicas %.3e, over 50 C4 replicas %.3e (MC, quoted separately)\n", sigR, cbRms, c4Rms);
    fclose(tb);
    tb = fopen("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_expected_asimov.txt", "w");   // _cnt copy: scratch
    fprintf(tb, "# x17c_expected_asimov.txt -- expected 90 %% CLs upper limits on BR(eta -> pi+pi-X17) BR(X17 -> e+e-) (asymptotic, q~_r, Asimov bands),\n");
    fprintf(tb, "# from the background-only Asimov set of the DATA side-band fit (x17c_freeze.C). S0, isBest >= 0, X2D; L = 6079.2 nb^-1; q~ / %.4f.\n", z.qScale);
    fprintf(tb, "# (1) r x 2.68e-4 with N5 fitted (PRIMARY); (2) eta yield from a_eta = %.3f +- %.3f; (3) eta yield from the peak, a_P = %.3f +- %.3f (WP5, cross-check value).\n",
            z.a0, z.aSig, z.aP, z.aPSig);
    fprintf(tb, "# m[MeV]   R_X    (1): -2s        -1s        median     +1s        +2s     | (2) median  | (3) median | X17 rows at (1) median | s_X2D [rows/BR, a=1]\n");
    for (const auto& r : rows)
        fprintf(tb, "%6.2f %7.4f  %10.3e %10.3e %10.3e %10.3e %10.3e | %10.3e | %10.3e | %8.1f | %.4e\n", r.m, r.RX, r.e1[0], r.e1[1], r.e1[2], r.e1[3], r.e1[4], r.e2[2], r.e3[2],
                r.nX, r.s1x);
    fclose(tb);

    // ------------------------------------------------------------ figures: side-band fit (M(ee) and M4 projections)
    {
        const int ny = I.G.ny, nx = I.G.nx;
        x17c::Model Mf(I, Call);
        TH1D hN("hN", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hc[4];
        const char* cn[4] = {"hCB", "hBg", "hPk", "hEt"};
        for (int k = 0; k < 4; ++k) { hc[k] = new TH1D(cn[k], "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back()); hc[k]->SetDirectory(nullptr); }
        TH1D hTot("hTot", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        for (int j = 0; j < ny; ++j) {
            double sn = 0, c[5], s[4] = {0, 0, 0, 0};
            for (int i = 0; i < nx; ++i) { Comps(Mf, SB.pFill, i, j, c); for (int k = 0; k < 4; ++k) s[k] += c[k]; sn += n[I.G.idx(i, j)]; }
            for (int k = 0; k < 4; ++k) hc[k]->SetBinContent(j + 1, s[k] / 2.);
            hTot.SetBinContent(j + 1, (s[0] + s[1] + s[2] + s[3]) / 2.);
            if (!(j >= kMaskLo && j <= kMaskHi)) { hN.SetBinContent(j + 1, sn / 2.); hN.SetBinError(j + 1, std::sqrt(sn) / 2.); }
        }
        TCanvas c("csb", "", 900, 850);
        TPad p1("p1", "", 0, 0.3, 1, 1), p2("p2", "", 0, 0, 1, 0.3);
        p1.SetBottomMargin(0.02); p1.SetLeftMargin(0.13); p1.SetLogy(); p2.SetTopMargin(0.02); p2.SetBottomMargin(0.32); p2.SetLeftMargin(0.13);
        p1.Draw(); p2.Draw(); p1.cd();
        THStack st("st", ";M(e^{+}e^{-}) [MeV/c^{2}];counts per MeV/c^{2}");
        const int col[4] = {kRed + 1, kGray + 1, kCyan + 2, kViolet + 2};
        for (int k : {0, 3, 2, 1}) { hc[k]->SetFillColor(col[k]); hc[k]->SetLineColor(col[k]); st.Add(hc[k]); }
        st.SetMinimum(20); st.SetMaximum(1e5);
        st.Draw("HIST");
        st.GetYaxis()->SetTitleSize(0.05); st.GetYaxis()->SetLabelSize(0.045); st.GetXaxis()->SetLabelSize(0);
        TBox bx(12, 20, 28, 1e5); bx.SetFillColor(kGray); bx.SetFillStyle(3354); bx.Draw();
        hTot.SetLineColor(kBlue + 1); hTot.SetLineWidth(2); hTot.SetLineStyle(2); hTot.Draw("HIST SAME");
        hN.SetMarkerStyle(20); hN.SetMarkerColor(kBlack); hN.SetLineColor(kBlack);
        TH1D hNs = hN; for (int j = kMaskLo; j <= kMaskHi; ++j) { hNs.SetBinContent(j + 1, 0); hNs.SetBinError(j + 1, 0); }
        hNs.Draw("E1 SAME");
        TLegend lg(0.50, 0.55, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.036); lg.SetFillStyle(0);
        lg.AddEntry(&hNs, "data N (S0), side bands", "pe");
        lg.AddEntry(hc[1], "other background n_{j} B(M4)", "f");
        lg.AddEntry(hc[2], "#eta peak (D + C)", "f");
        lg.AddEntry(hc[3], "#eta-related (C2+C3, G_{#eta})", "f");
        lg.AddEntry(hc[0], "CB (like-sign data)", "f");
        lg.AddEntry(&hTot, "model; 12-28: predicted (no data)", "l");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
        tx.SetTextSize(0.03);
        tx.DrawLatex(0.15, 0.05, "S0, M4 0.5225-0.610 GeV/c^{2}; fit with r = 0; 12-28 MeV/c^{2} blinded (D19); L = 6079.2 nb^{-1}");
        p2.cd();
        TH1D hr("hr", ";M(e^{+}e^{-}) [MeV/c^{2}];data / model (n_{j} free)", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        for (int j = 0; j < ny; ++j) {
            if (j >= kMaskLo && j <= kMaskHi) continue;
            const double t = hTot.GetBinContent(j + 1);
            if (t > 0) { hr.SetBinContent(j + 1, hN.GetBinContent(j + 1) / t); hr.SetBinError(j + 1, hN.GetBinError(j + 1) / t); }
        }
        hr.SetMinimum(0.97); hr.SetMaximum(1.03);
        hr.GetXaxis()->SetTitleSize(0.12); hr.GetXaxis()->SetLabelSize(0.10); hr.GetYaxis()->SetTitleSize(0.10); hr.GetYaxis()->SetLabelSize(0.08);
        hr.GetYaxis()->SetTitleOffset(0.55); hr.GetYaxis()->SetNdivisions(505);
        hr.SetMarkerStyle(20); hr.Draw("E1");
        TLine l1(2, 1, 62, 1); l1.SetLineStyle(2); l1.Draw();
        TBox bx2(12, 0.97, 28, 1.03); bx2.SetFillColor(kGray); bx2.SetFillStyle(3354); bx2.Draw();
        c.SaveAs(Form("%s/x17c_sb_mee_cnt.pdf", kFig)); c.SaveAs(Form("%s/x17c_sb_mee_cnt.png", kFig));
        for (int k = 0; k < 4; ++k) delete hc[k];
    }
    {
        // M4 projections of the two side bands with pulls
        const int nx = I.G.nx;
        x17c::Model Mf(I, Call);
        TCanvas c("csb4", "", 1400, 800);
        for (int band = 0; band < 2; ++band) {
            c.cd();
            TPad* p1 = new TPad(Form("q1%d", band), "", 0.5 * band, 0.3, 0.5 * (band + 1), 1);
            TPad* p2 = new TPad(Form("q2%d", band), "", 0.5 * band, 0, 0.5 * (band + 1), 0.3);
            p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.15); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.15);
            p1->Draw(); p2->Draw(); p1->cd();
            TH1D* hN = new TH1D(Form("hN4%d", band), ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];counts per 2.5 MeV/c^{2}", nx, I.G.xe.data());
            TH1D* hc[4];
            for (int k = 0; k < 4; ++k) hc[k] = new TH1D(Form("hc4%d%d", band, k), "", nx, I.G.xe.data());
            for (int j = 0; j < I.G.ny; ++j) {
                const bool inb = band == 0 ? (j < kMaskLo) : (j > kMaskHi);
                if (!inb) continue;
                for (int i = 0; i < nx; ++i) {
                    double cc[5]; Comps(Mf, SB.pFill, i, j, cc);
                    for (int k = 0; k < 4; ++k) hc[k]->AddBinContent(i + 1, cc[k]);
                    hN->AddBinContent(i + 1, n[I.G.idx(i, j)]);
                }
            }
            for (int i = 1; i <= nx; ++i) hN->SetBinError(i, std::sqrt(hN->GetBinContent(i)));
            THStack* st = new THStack(Form("st4%d", band), ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];counts per 2.5 MeV/c^{2}");
            const int col[4] = {kRed + 1, kGray + 1, kCyan + 2, kViolet + 2};
            for (int k : {0, 3, 2, 1}) { hc[k]->SetFillColor(col[k]); hc[k]->SetLineColor(col[k]); st->Add(hc[k]); }
            double mx = hN->GetMaximum() * 1.25;
            if (band == 1) mx = hN->GetMaximum() * 1.75;
            st->SetMaximum(mx); st->SetMinimum(0);
            st->Draw("HIST");
            st->GetYaxis()->SetTitleSize(0.05); st->GetYaxis()->SetLabelSize(0.045); st->GetXaxis()->SetLabelSize(0); st->GetYaxis()->SetTitleOffset(1.4);
            hN->SetMarkerStyle(20); hN->Draw("E1 SAME");
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
            tx.DrawLatex(0.2, 0.85, band == 0 ? "M(e^{+}e^{-}) 2-12 MeV/c^{2}" : "M(e^{+}e^{-}) 28-62 MeV/c^{2}");
            if (band == 1) {
                TLegend* lg = new TLegend(0.18, 0.50, 0.55, 0.82); lg->SetBorderSize(0); lg->SetTextSize(0.038); lg->SetFillStyle(0);
                lg->AddEntry(hN, "data N (S0)", "pe"); lg->AddEntry(hc[1], "other background", "f"); lg->AddEntry(hc[2], "#eta peak (D + C)", "f");
                lg->AddEntry(hc[3], "#eta-related", "f"); lg->AddEntry(hc[0], "CB", "f"); lg->Draw();
            }
            p2->cd();
            TH1D* hp = new TH1D(Form("hp4%d", band), ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];(N - #mu)/#sqrt{#mu}", nx, I.G.xe.data());
            for (int i = 1; i <= nx; ++i) {
                double t = 0; for (int k = 0; k < 4; ++k) t += hc[k]->GetBinContent(i);
                if (t > 0) hp->SetBinContent(i, (hN->GetBinContent(i) - t) / std::sqrt(t));
            }
            hp->SetMinimum(-4); hp->SetMaximum(4); hp->SetFillColor(kGray + 2);
            hp->GetXaxis()->SetTitleSize(0.11); hp->GetXaxis()->SetLabelSize(0.09); hp->GetYaxis()->SetTitleSize(0.09); hp->GetYaxis()->SetLabelSize(0.08);
            hp->GetYaxis()->SetTitleOffset(0.7); hp->GetXaxis()->SetNdivisions(505);
            hp->Draw("HIST");
        }
        c.SaveAs(Form("%s/x17c_sb_m4_cnt.pdf", kFig)); c.SaveAs(Form("%s/x17c_sb_m4_cnt.png", kFig));
    }
    {
        // masked-row prediction: ln(n_j / sim_j) of the side-band rows and the fill curves
        TCanvas c("cfill", "", 900, 700);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);
        TGraphErrors gd;
        TGraph gp;
        for (int j = 0; j < I.G.ny; ++j) {
            const double y = 1e3 * I.G.yc(j);
            if (j >= kMaskLo && j <= kMaskHi) { gp.SetPoint(gp.GetN(), y, SB.pFill[j] / simRow[j]); continue; }
            gd.SetPoint(gd.GetN(), y, SB.pFill[j] / simRow[j]);
            gd.SetPointError(gd.GetN() - 1, 1., SB.pFill[j] / simRow[j] * std::hypot(SB.fit.err[j] / SB.pFill[j], simRel[j]));
        }
        gd.SetTitle(";M(e^{+}e^{-}) [MeV/c^{2}];n_{j} / simulated (G_{D} + G_{C}) row sum");
        gd.SetMarkerStyle(20); gd.Draw("AP");
        gd.GetXaxis()->SetLimits(0, 64); gd.GetXaxis()->SetTitleSize(0.045); gd.GetYaxis()->SetTitleSize(0.045);
        gp.SetMarkerStyle(24); gp.SetMarkerColor(kBlue + 1); gp.SetMarkerSize(1.3); gp.Draw("P SAME");
        TLegend lg(0.36, 0.74, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.028);
        lg.AddEntry(&gd, "side-band rows: fitted n_{j} (data); error: data (+) simulation", "pe");
        lg.AddEntry(&gp, "masked rows 12-28 MeV/c^{2}: prediction (pol2 in ln of the ratio)", "p");
        TGraph gr;
        for (int j = 0; j < I.G.ny; ++j) gr.SetPoint(gr.GetN(), 1e3 * I.G.yc(j), SB.pFill[j] / simSm[j]);
        gr.SetMarkerStyle(5); gr.SetMarkerColor(kGray + 2); gr.Draw("P SAME");
        lg.AddEntry(&gr, "n_{j} / locally smoothed simulation (variant reference)", "p");
        lg.Draw();
        c.SaveAs("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_fill.pdf"); c.SaveAs("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17c/x17c_fill.png");   // _cnt copy: scratch only (its "rows" are M(ee) rows of the X2D grid, unchanged)
    }
    f->Close();
    P("x17c_freeze: done\n");
    fclose(gLog);
}
