// ============================================================================
// x17d_unblind.C  --  WP9 X-d (prefix x17d): X17 UNBLINDING. The FROZEN procedure of x17_frozen.txt applied to the full
//                     data X2D (side bands of hists_v1.root + the masked rows 12-28 MeV/c2 of the sealed file).
// ============================================================================
// SPEC report_all WP9 X-d; x17_frozen.txt Sect. 10 steps 1, 2, 3 (local p0), 4, 5 (1 + cos^2), 6. Written 2026-11-14.
// Nothing of the frozen procedure is changed: model x17c_model.h, configuration x17c_config.txt (x17cl::ReadConfig),
// inputs x17cl::Build (nominal BuildOpt), line x17cl::SetLine, limit x17cl::LimitBR (asymptotic CLs, q~ / qScale).
// Steps:
//   0. Reproduction of the frozen chain from the side bands alone (B_pred(X17c) = 6311.7, Asimov median 2.759e-5 at
//      17 MeV/c2): proves that inputs, code and configuration are those of the freeze.
//   1. Full fit (all 1050 cells, r = 0 and r free at 17 MeV/c2): parameters, Pearson chi2/ndf (frozen fall-back:
//      qScale := max(qScale, chi2/ndf) if > 1.5), chi2 and peak-band residual per M(ee) row, fitted n_j of the formerly
//      masked rows against their frozen prediction.
//   2. Mass scan (38 masses of x17cl::MassGrid), normalisations (1), (2), (3): observed 90 % CLs UL, r_hat +- error,
//      Asimov bands of the OBSERVED r = 0 fit (labelled, in addition to the frozen expectation), local p0; the e+ decay
//      angle 1 + cos^2 limit (quoted separately).
//   3. Counting cross-check in X17c with the frozen B_pred and sigma_B.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/x17d/x17d_state.root, x17d_unblind.log;
//         docs/latex/tables_all/x17d_scan.txt, x17d_fit.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)' \
//        'plots/all_report/x17d_unblind.C+O'
// ============================================================================
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"

#include "x17d_lib.h"

namespace x17du {
FILE* gLog = nullptr;
void P(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
void P(const char* fmt, ...) {
    va_list a;
    va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gLog) { va_start(a, fmt); vfprintf(gLog, fmt, a); va_end(a); fflush(gLog); }
}
}  // namespace x17du

void x17d_unblind() {
    using namespace x17d;
    using x17du::P;
    TStopwatch sw;
    gSystem->mkdir(kOutD, kTRUE);
    x17du::gLog = fopen(Form("%s/x17d_unblind.log", kOutD), "w");
    Ctx c;
    if (!Init(c)) { P("x17d: init failed\n"); return; }
    x17c::Inputs& I = c.I;
    Frozen& z = c.z;
    const std::vector<double>& n = c.D.n;
    const int nx = I.G.nx, ny = I.G.ny;
    double nTot = 0; for (double v : n) nTot += v;
    P("x17d: frozen configuration: a0 %.6g aSig %.6g g0 %.6g gSig %.6g kE1 %.3g kE2 %.3g cTail %.5g sX %.6g/%.6g/%.6g qScale %.6g aP %.6g bkgMode %d bshape %d fill %d\n",
      z.a0, z.aSig, z.g0, z.gSig, z.kE1, z.kE2, z.cTail, z.xSig1, z.xSig2, z.xSig3, z.qScale, z.aP, z.bkgMode, z.bshape, z.fillOrder);
    P("x17d: DATA X2D (S0): %.0f rows in 1050 cells; side bands %.0f rows (hists_v1.root), formerly masked rows 12-28 MeV/c2 %.0f rows (sealed file);"
      " sealed side-band cells = hists_v1.root: max |diff| %.3g -> %s\n", nTot, c.D.nSB, c.D.nMasked, c.D.maxDiffSB, c.D.ok ? "IDENTICAL" : "DIFFERENT");
    P("x17d: N5ref %.1f (frozen 7458.8), R_X(17) %.4f (frozen 1.3455), s_X2D %.4e\n", I.N5ref, I.RX, I.N5ref * I.RX / kBR);

    // ------------------------------------------------------------ 0. reproduction of the frozen chain (side bands only)
    const std::vector<double> simRow = SimRow(c.f, I), simRel = SimRowRelErr(c.f, I);
    x17c::Config C1 = MakeConfig(z, 1);
    x17c::Config Csb = C1; Csb.use = c.D.useSB;
    SideBand SB = FitSideBands(I, Csb, c.D.nSBonly, simRow, z.fillOrder, &simRel);
    x17c::Model M1(I, C1);
    const x17c::Layout& L = M1.L;
    const std::vector<double> asimovFrozen = M1.Expected(SB.pFill);
    const double bPredRe = WindowX17c(I.G, asimovFrozen);
    Lim lFz = LimitBR(I, C1, asimovFrozen, true);
    P("x17d: REPRODUCTION of the freeze: side-band chi2 %.1f/%d, N5 %.1f; B_pred(X17c) %.1f (frozen %.1f); Asimov median 17 MeV/c2 (1) %.4e (frozen %.4e)\n",
      SB.chi2, SB.ndf, SB.fit.p[L.iN5], bPredRe, kBpredX17c, lFz.exp[2], kFrozenExpMed17);
    const bool reproOK = std::fabs(bPredRe - kBpredX17c) < 0.1 && std::fabs(lFz.exp[2] / kFrozenExpMed17 - 1) < 2e-3;
    P("x17d: reproduction %s\n", reproOK ? "OK" : "FAILED");

    // ------------------------------------------------------------ 1. full fit (all cells)
    const x17c::Glob g1 = x17c::GlobOf(C1);
    const x17c::Profile Pf = x17c::ProfileFits(M1, n, g1);
    Gof G0 = GoodnessOfFit(M1, Pf.zero, n), Gh = GoodnessOfFit(M1, Pf.hat, n);
    std::vector<char> useM(n.size(), 0);
    for (size_t b = 0; b < n.size(); ++b) useM[b] = !c.D.useSB[b];
    Gof G0m = GoodnessOfFit(M1, Pf.zero, n, &useM), G0s = GoodnessOfFit(M1, Pf.zero, n, &c.D.useSB);
    const int ndf0 = G0.ncell - G0.nfree, ndfh = Gh.ncell - Gh.nfree;
    const double chi2ndf = std::max(G0.chi2 / ndf0, Gh.chi2 / ndfh);
    const bool fbPearson = chi2ndf > 1.5;
    P("x17d: FULL FIT r = 0: ok %d iter %d; Pearson chi2 %.1f / ndf %d = %.4f (deviance %.1f); masked rows only %.1f / %d cells; side-band cells %.1f / %d cells\n",
      Pf.zero.ok, Pf.zero.iter, G0.chi2, ndf0, G0.chi2 / ndf0, G0.dev, G0m.chi2, G0m.ncell, G0s.chi2, G0s.ncell);
    P("x17d: FULL FIT r free (17 MeV/c2): ok %d; Pearson chi2 %.1f / ndf %d = %.4f; r_hat %.4e +- %.4e (BR)\n", Pf.hat.ok, Gh.chi2, ndfh, Gh.chi2 / ndfh, Pf.rhat * kBR,
      Pf.rhatErr * kBR);
    P("x17d: frozen fall-back 'Pearson chi2/ndf > 1.5 -> qScale := chi2/ndf': %s (%.4f); qScale stays %.4f\n", fbPearson ? "FIRED" : "not fired", chi2ndf,
      fbPearson ? chi2ndf : z.qScale);
    if (fbPearson) z.qScale = std::max(z.qScale, chi2ndf);
    const std::vector<double>& p0 = Pf.zero.p;
    const std::vector<double>& e0 = Pf.zero.err;
    P("x17d: r = 0 fit: N5 %.1f +- %.1f (%.1f %%), N4 %.1f +- %.1f, k_E1 %.3f +- %.3f, k_E2 %.4f +- %.4f, a_eta %.4f +- %.4f, g %.3f +- %.3f, thC %.3f +- %.3f\n", p0[L.iN5],
      e0[L.iN5], 100 * e0[L.iN5] / p0[L.iN5], p0[L.iN4], e0[L.iN4], p0[L.iK], e0[L.iK], p0[L.iK + 1], e0[L.iK + 1], p0[L.iA], e0[L.iA], p0[L.iG], e0[L.iG], p0[L.iTC],
      e0[L.iTC]);
    for (int cl = 0; cl < I.nbcl; ++cl)
        P("x17d: class %d: b %.4f +- %.4f, c %.4f +- %.4f, b' %.4f +- %.4f, c' %.4f +- %.4f\n", cl, p0[L.iB + cl], e0[L.iB + cl], p0[L.iC + cl], e0[L.iC + cl],
          p0[L.iBy + cl], e0[L.iBy + cl], p0[L.iCy + cl], e0[L.iCy + cl]);
    double predM = 0, fitM = 0;
    for (int j = 0; j < ny; ++j) {
        const bool msk = j >= kMaskLo && j <= kMaskHi;
        double nr = 0; for (int i = 0; i < nx; ++i) nr += n[I.G.idx(i, j)];
        const double sres = G0.varPeak[j] > 0 ? G0.resPeak[j] / std::sqrt(G0.varPeak[j]) : 0.;
        P("x17d: row %2d %4.0f-%4.0f MeV %s: N %8.0f  n_j %10.1f +- %6.1f  %s %10.1f  (n_j/pred %.4f)  chi2 %5.1f/35  peak-band N-mu %+7.1f (%+.2f sigma)\n", j,
          1e3 * I.G.ye[j], 1e3 * I.G.ye[j + 1], msk ? "UNBL" : "sb  ", nr, p0[j], e0[j], msk ? "pred" : "sbft", SB.pFill[j], p0[j] / SB.pFill[j], G0.chi2row[j],
          G0.resPeak[j], sres);
        if (msk) { predM += SB.pFill[j]; fitM += p0[j]; }
    }
    P("x17d: formerly masked rows: fitted sum n_j %.1f against the frozen prediction %.1f (ratio %.4f)\n", fitM, predM, fitM / predM);
    // peak-band residual summed over the unblinded rows (0.5325-0.5575)
    {
        double r = 0, v = 0, r2 = 0, v2 = 0;
        for (int j = kMaskLo; j <= kMaskHi; ++j) { r += G0.resPeak[j]; v += G0.varPeak[j]; }
        for (int j = 0; j < ny; ++j) if (I.G.ye[j] >= 0.014 - 1e-9 && I.G.ye[j + 1] <= 0.026 + 1e-9) { r2 += G0.resPeak[j]; v2 += G0.varPeak[j]; }
        P("x17d: peak band M4 0.5325-0.5575, r = 0 fit: N - mu in 12-28 MeV/c2 %+.1f +- %.1f (%+.2f sigma); in 14-26 MeV/c2 %+.1f +- %.1f (%+.2f sigma)\n", r, std::sqrt(v),
          r / std::sqrt(v), r2, std::sqrt(v2), r2 / std::sqrt(v2));
    }

    // ------------------------------------------------------------ 3. counting cross-check in X17c (frozen B_pred, sigma_B)
    double nX17cExact = 0, cbX17cExact = 0;
    {
        TH1* hw = (TH1*)c.f->Get("data/S0/all/N_win");
        TH1* hc = (TH1*)c.f->Get("data/S0/all/CB_win");
        if (hw) nX17cExact = hw->GetBinContent(alldefs::kWX17c + 1);
        if (hc) cbX17cExact = hc->GetBinContent(alldefs::kWX17c + 1);
    }
    const double nX17cOv = WindowX17c(I.G, n), cbX17cOv = WindowX17c(I.G, I.CB);
    auto countBR = [&](double N, double& ex, double& err, double& br, double& brErr) {
        ex = N - kBpredX17c;
        err = std::hypot(std::sqrt(N), kSigBX17c);
        br = ex / (z.a0 * kS1X17c);
        brErr = err / (z.a0 * kS1X17c);
    };
    double exO, erO, brO, bseO, exE, erE, brE, bseE;
    countBR(nX17cOv, exO, erO, brO, bseO);
    countBR(nX17cExact, exE, erE, brE, bseE);
    P("x17d: COUNTING X17c (M4 0.5312-0.5566, M(ee) 15.5-24.5): N (X2D overlap fractions) %.1f, N (exact row count, N_win) %.0f; CB %.1f / %.1f; B_pred %.1f +- %.1f (frozen)\n",
      nX17cOv, nX17cExact, cbX17cOv, cbX17cExact, kBpredX17c, kSigBX17c);
    P("x17d: COUNTING excess (overlap N, primary) %+.1f +- %.1f rows -> BR = excess / (a_eta s1) = %+.3e +- %.3e  [a_eta %.4f, s1 %.4e]; Gaussian 90 %% one-sided bound %.3e\n",
      exO, erO, brO, bseO, z.a0, kS1X17c, std::max(brO, 0.) + 1.2816 * bseO);
    P("x17d: COUNTING excess (exact N, information) %+.1f +- %.1f rows -> BR %+.3e +- %.3e\n", exE, erE, brE, bseE);

    // ------------------------------------------------------------ 2. mass scan
    const std::vector<double> mg = MassGrid();
    TFile* fo = TFile::Open(Form("%s/x17d_state.root", kOutD), "RECREATE");
    TTree t("obs", "X-d observed limits (BR_X BR_ee units) and local p0, frozen procedure");
    double m, RX, ul[3], rhat[3], rerr[3], ea[3][5], q0[3], Z[3], pv[3], nX[3], N5h, ulC2, s1x;
    int ok[3], okp[3];
    t.Branch("m", &m); t.Branch("RX", &RX); t.Branch("ul", ul, "ul[3]/D"); t.Branch("rhat", rhat, "rhat[3]/D"); t.Branch("rerr", rerr, "rerr[3]/D");
    t.Branch("expA", ea, "expA[3][5]/D"); t.Branch("q0", q0, "q0[3]/D"); t.Branch("Z", Z, "Z[3]/D"); t.Branch("p0", pv, "p0[3]/D"); t.Branch("nX", nX, "nX[3]/D");
    t.Branch("N5", &N5h); t.Branch("ulCos2", &ulC2); t.Branch("s1x", &s1x); t.Branch("ok", ok, "ok[3]/I"); t.Branch("okp", okp, "okp[3]/I");
    FILE* tb = fopen(Form("%s/x17d_scan.txt", kTabD), "w");
    fprintf(tb, "# x17d_scan.txt -- X17 UNBLINDED (x17d_unblind.C): OBSERVED 90 %% CLs upper limits on BR(eta -> pi+pi-X17) BR(X17 -> e+e-), asymptotic CLs\n");
    fprintf(tb, "# (q~_r / qScale %.4f), FROZEN procedure of x17_frozen.txt, all 1050 X2D cells (S0, isBest >= 0, M4 0.5225-0.6100 x M(ee) 2-62 MeV/c2);\n", z.qScale);
    fprintf(tb, "# L = 6079.2 nb^-1. (1) r x 2.68e-4, N5 fitted (PRIMARY); (2) eta yield a_eta = %.4f (WP5); (3) eta yield from the peak a_P = %.4f.\n", z.a0, z.aP);
    fprintf(tb, "# r_hat in BR units (norm 1) with its Hessian error; Z, p0 = local significance (asymptotic, q0 / qScale; Z = 0 for r_hat <= 0);\n");
    fprintf(tb, "# expA = Asimov median of the OBSERVED r = 0 fit (labelled; the frozen expectation is x17c_expected_asimov.txt); cos2 = (1) with e+ angle 1 + cos^2.\n");
    fprintf(tb, "# m[MeV]   R_X   UL(1)      UL(2)      UL(3)     | r_hat(1)    err      | Z(1)   p0(1)   | expA(1) med | N_X UL(1) rows | UL(1) cos2 | N5_hat\n");
    for (double mm : mg) {
        m = mm;
        x17c::Inputs Im = I;
        if (!SetLine(Im, mm)) continue;
        RX = Im.RX;
        s1x = Im.N5ref * Im.RX / kBR;
        for (int a = 0; a < 3; ++a) {
            x17c::Config Cn = MakeConfig(z, a + 1);
            Lim l = LimitBR(Im, Cn, n, true);
            ok[a] = l.ok; ul[a] = l.ul; rhat[a] = l.rhat; rerr[a] = l.rerr; nX[a] = l.nX;
            for (int k = 0; k < 5; ++k) ea[a][k] = l.exp[k];
            if (a == 0) N5h = l.N5;
            x17c::Model Mn(Im, Cn);
            P0 pz = LocalP0(Mn, n, x17c::GlobOf(Cn));
            okp[a] = pz.ok; q0[a] = pz.q0; Z[a] = pz.Z; pv[a] = pz.p0;
        }
        {
            x17c::Inputs Ic = I;
            SetLine(Ic, mm, 2);
            Lim l = LimitBR(Ic, C1, n, false);
            ulC2 = l.ul;
        }
        t.Fill();
        P("x17d: m %5.2f R_X %.4f  UL (1) %.3e (2) %.3e (3) %.3e | r_hat %+.3e +- %.3e | Z %.2f p0 %.3f | Asimov(obs) med %.3e | N_X %.1f | cos2 %.3e | N5 %.0f | ok %d%d%d %d%d%d\n",
          mm, RX, ul[0], ul[1], ul[2], rhat[0], rerr[0], Z[0], pv[0], ea[0][2], nX[0], ulC2, N5h, ok[0], ok[1], ok[2], okp[0], okp[1], okp[2]);
        fprintf(tb, "%6.2f %7.4f %10.3e %10.3e %10.3e | %+10.3e %9.3e | %5.2f %7.4f | %10.3e | %8.1f | %10.3e | %7.1f\n", mm, RX, ul[0], ul[1], ul[2], rhat[0], rerr[0], Z[0],
                pv[0], ea[0][2], nX[0], ulC2, N5h);
    }
    fclose(tb);
    t.Write();
    // parameters for the figures: r = 0 fit and free fit at 17 MeV/c2 (norm 1), side-band fit with the predicted rows
    {
        TVectorD a((int)p0.size()), b((int)Pf.hat.p.size()), s((int)SB.pFill.size()), e((int)e0.size());
        for (size_t k = 0; k < p0.size(); ++k) { a((int)k) = p0[k]; e((int)k) = e0[k]; }
        for (size_t k = 0; k < Pf.hat.p.size(); ++k) b((int)k) = Pf.hat.p[k];
        for (size_t k = 0; k < SB.pFill.size(); ++k) s((int)k) = SB.pFill[k];
        a.Write("pZero"); e.Write("eZero"); b.Write("pHat17"); s.Write("pSBfill");
        TH2D hd("data", "data N X2D (S0, all cells, unblinded)", nx, I.G.xe.data(), ny, I.G.ye.data());
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) hd.SetBinContent(i + 1, j + 1, n[I.G.idx(i, j)]);
        hd.Write();
        TVectorD misc(16);
        misc(0) = G0.chi2; misc(1) = ndf0; misc(2) = Gh.chi2; misc(3) = ndfh; misc(4) = z.qScale; misc(5) = nX17cOv; misc(6) = nX17cExact; misc(7) = cbX17cOv;
        misc(8) = exO; misc(9) = erO; misc(10) = brO; misc(11) = bseO; misc(12) = fitM; misc(13) = predM; misc(14) = G0m.chi2; misc(15) = G0m.ncell;
        misc.Write("misc");
    }
    fo->Close();
    // fit table
    tb = fopen(Form("%s/x17d_fit.txt", kTabD), "w");
    fprintf(tb, "# x17d_fit.txt -- X17 unblinded: full fit of the frozen model (normalisation (1), r = 0 and r free at 17 MeV/c2) to all 1050 X2D cells (x17d_unblind.C).\n");
    fprintf(tb, "# Data S0, isBest >= 0; raw rows; L = 6079.2 nb^-1. Rows 12-28 MeV/c2 were blinded until x17_frozen.txt; their n_j are compared with the frozen prediction.\n");
    fprintf(tb, "reproduction of the freeze (side bands only): B_pred(X17c) %.1f (frozen %.1f); Asimov median 17 MeV (1) %.4e (frozen %.4e): %s\n", bPredRe, kBpredX17c,
            lFz.exp[2], kFrozenExpMed17, reproOK ? "OK" : "FAILED");
    fprintf(tb, "data rows: all %.0f, side bands %.0f, formerly masked %.0f; sealed side-band cells = hists_v1.root (max |diff| %.3g)\n", nTot, c.D.nSB, c.D.nMasked, c.D.maxDiffSB);
    fprintf(tb, "r = 0 fit: Pearson chi2 %.1f / ndf %d = %.4f; deviance %.1f; formerly masked cells %.1f / %d; side-band cells %.1f / %d\n", G0.chi2, ndf0, G0.chi2 / ndf0, G0.dev,
            G0m.chi2, G0m.ncell, G0s.chi2, G0s.ncell);
    fprintf(tb, "r free (17 MeV): Pearson chi2 %.1f / ndf %d = %.4f; r_hat %.4e +- %.4e (BR); fall-back chi2/ndf > 1.5: %s; qScale %.4f\n", Gh.chi2, ndfh, Gh.chi2 / ndfh,
            Pf.rhat * kBR, Pf.rhatErr * kBR, fbPearson ? "FIRED" : "not fired", z.qScale);
    fprintf(tb, "N5 %.1f +- %.1f; N4 %.1f +- %.1f; k_E1 %.3f +- %.3f; k_E2 %.4f +- %.4f; a_eta %.4f +- %.4f; g %.3f +- %.3f; thC %.3f +- %.3f\n", p0[L.iN5], e0[L.iN5],
            p0[L.iN4], e0[L.iN4], p0[L.iK], e0[L.iK], p0[L.iK + 1], e0[L.iK + 1], p0[L.iA], e0[L.iA], p0[L.iG], e0[L.iG], p0[L.iTC], e0[L.iTC]);
    fprintf(tb, "masked rows: fitted sum n_j %.1f, frozen prediction %.1f, ratio %.4f\n", fitM, predM, fitM / predM);
    fprintf(tb, "counting X17c: N %.1f (overlap) / %.0f (exact); CB %.1f / %.1f; B_pred %.1f +- %.1f; excess %+.1f +- %.1f; BR %+.3e +- %.3e (exact N: %+.1f, %+.3e)\n",
            nX17cOv, nX17cExact, cbX17cOv, cbX17cExact, kBpredX17c, kSigBX17c, exO, erO, brO, bseO, exE, brE);
    fprintf(tb, "# Mee_lo Mee_hi status  N_row     n_j       err     prediction  n_j/pred  chi2_row/35  peakband(0.5325-0.5575) N-mu  sigma\n");
    for (int j = 0; j < ny; ++j) {
        const bool msk = j >= kMaskLo && j <= kMaskHi;
        double nr = 0; for (int i = 0; i < nx; ++i) nr += n[I.G.idx(i, j)];
        fprintf(tb, "%5.0f %5.0f  %s %8.0f %10.1f %8.1f %10.1f %9.4f %9.1f %+10.1f %+6.2f\n", 1e3 * I.G.ye[j], 1e3 * I.G.ye[j + 1], msk ? "unblinded" : "side-band", nr, p0[j],
                e0[j], SB.pFill[j], p0[j] / SB.pFill[j], G0.chi2row[j], G0.resPeak[j], G0.varPeak[j] > 0 ? G0.resPeak[j] / std::sqrt(G0.varPeak[j]) : 0.);
    }
    fclose(tb);
    c.f->Close();
    P("x17d_unblind: done, %.0f s\n", sw.RealTime());
    fclose(x17du::gLog);
}
